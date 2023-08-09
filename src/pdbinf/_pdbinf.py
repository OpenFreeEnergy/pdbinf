# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
"""pdbinf - assign PDB bond information from PDBx templates


PDBFile requirements
--------------------

If inputting a PDB file, this must be "PDB compliant".  In practical terms,
this means it should follow standard naming for atoms and residues.  One way to
enforce this is to first pass the file through `pdbfixer`
(https://github.com/openmm/pdbfixer)

RDKit Mol requirements
----------------------

RDKit molecules that are input are required to have:
- MonomerInfo() assigned, including
  - Name
  - ResidueName
  - ResidueNumber
  - InsertionCode
  - ChainId
- a Conformer (only the first will be used)

Template requirements
---------------------
Currently the following PDBx fields are used and required:

_chem_comp_atom
 - atom_id
   matches against rdkit atom MonomerInfo().GetName() property
 - pdbx_aromatic_flag
   used to set IsAromatic on rdkit Atoms

_chem_comp_bond
 - atom_id_1, atom_id_2
   must match against atom_id in _chem_comp_atom entries
 - pdbx_aromatic_flag
   if 'Y' assigns bonds to Chem.BondType.AROMATIC
 - value_order
   used to calculate the valence table, which is used for formal charge

"""
from collections import namedtuple
import itertools
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
import logging
import pathlib
import numpy as np
import os
from typing import Iterator, Union

from . import RESNAME_ALIASES, ATOMNAME_ALIASES

logger = logging.getLogger(__name__)

PT = Chem.GetPeriodicTable()
STANDARD_RESNAMES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}
MAX_AMIDE_LENGTH = 2.0
MAX_DISULPHIDE_LENGTH = 2.5

Residue = namedtuple('Residue', 'resname,resid,icode,chain')


def strip_bonds(m: Chem.Mol) -> Chem.Mol:
    em = AllChem.EditableMol(m)

    for b in m.GetBonds():
        em.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

    return em.GetMol()


def residue_spans(m: Chem.Mol) -> Iterator[tuple[int, int, Residue]]:
    """a generator of spans demarking different residues

    defines a residue as being a tuple of (resname, resnum) from the MI

    Returns
    -------
    generator of (int, int, resname, chainid) tuples
      indicating that the residue called *resname* is between the two ints
    """
    def reshash(atom) -> Residue:
        mi = atom.GetMonomerInfo()
        resname = mi.GetResidueName().strip()
        resnum = mi.GetResidueNumber()
        icode = mi.GetInsertionCode()
        chain = mi.GetChainId()

        return Residue(resname, resnum, icode, chain)

    begin = 0
    current = reshash(m.GetAtomWithIdx(0))

    for i in range(1, m.GetNumAtoms()):
        nxt = reshash(m.GetAtomWithIdx(i))

        if nxt != current:
            yield begin, i, current

            current = nxt
            begin = i

    yield begin, m.GetNumAtoms(), current


def assign_intra_props(mol, atom_span: range, reference_block):
    """for atoms within span, assign bonds and aromaticity based on NAMES

    """
    nm_2_idx = dict()

    # convert indices to names
    aliases = ATOMNAME_ALIASES.get(reference_block.name, {})
    for idx in atom_span:
        atom = mol.GetAtomWithIdx(idx)
        nm = atom.GetMonomerInfo().GetName().strip()
        nm = aliases.get(nm, nm)
        nm_2_idx[nm] = idx

    logger.debug(f'assigning intra props for {reference_block.name}')

    em = AllChem.EditableMol(mol)

    # we'll assign rdkit SINGLE, DOUBLE or AROMATIC bonds
    # but we'll also want to know the original *valence*
    valence = np.zeros(mol.GetNumAtoms(), dtype=int)

    # grab bond data from gemmi Block
    for row in reference_block.find(
            '_chem_comp_bond.',
            ['atom_id_1', 'atom_id_2', 'pdbx_aromatic_flag', 'value_order']
    ):
        nm1, nm2, arom, order = row.str(0), row.str(1), row[2], row[3]
        try:
            idx1, idx2 = nm_2_idx[nm1], nm_2_idx[nm2]
        except KeyError:
            continue

        v = 1 if order == 'SING' else 2  # 'DOUB'
        valence[idx1] += v
        valence[idx2] += v

        if arom == 'Y':
            order = Chem.BondType.AROMATIC
        elif v == 1:
            order = Chem.BondType.SINGLE
        else:  # order == 'DOUB'
            order = Chem.BondType.DOUBLE

        logger.debug(f"adding bond: {nm1}-{nm2} at {idx1} {idx2} {order}")

        em.AddBond(idx1, idx2, order=order)

    mol = em.GetMol()

    # find lone hydrogens, then attach to the nearest heavy atom
    additional_bonds = []
    heavy_atoms = []
    lone_H = []
    for idx in atom_span:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 1:
            heavy_atoms.append(idx)
            continue
        if atom.GetBonds():  # if any bonds, we're ok
            continue
        lone_H.append(idx)
    if lone_H:
        logger.debug(f"found lone hydrogens: {lone_H}")
        conf = mol.GetConformer()
        for idx in lone_H:
            pos = conf.GetAtomPosition(idx)
            minidx, mindist = -1, float('inf')
            for idx2 in heavy_atoms:
                pos2 = conf.GetAtomPosition(idx2)
                d = pos.Distance(pos2)
                if d > mindist:
                    continue
                minidx, mindist = idx2, d

            if mindist < MAX_AMIDE_LENGTH:
                logger.debug(f"attached hydrogen {idx} to heavy atom {minidx} d={mindist}")
                additional_bonds.append((idx, minidx))
    if additional_bonds:
        em = Chem.EditableMol(mol)
        for i, j in additional_bonds:
            em.AddBond(i, j, order=Chem.BondType.SINGLE)
        mol = em.GetMol()

    for nm, arom in reference_block.find('_chem_comp_atom.',
                                         ['atom_id', 'pdbx_aromatic_flag']):
        try:
            idx = nm_2_idx[nm]
        except KeyError:
            # todo: could check atom is marked as leaving atom
            continue

        atom = mol.GetAtomWithIdx(idx)
        atom.SetIsAromatic(arom == 'Y')

    return mol, valence


def load_pdb_file(pdb_path: Union[str, os.PathLike],
                  templates: list[gemmi.cif.Document]) -> Chem.Mol:
    """Load a PDB File and assign bonds (with order)

    Parameters
    ----------
    pdb_path : Union[str, pathlib.Path]
      path to the file to be loaded
    templates : list[gemmi.cif.Document]
    """
    m = Chem.MolFromPDBFile(str(pdb_path),  # rdkit doesn't like Path objects
                            flavor=1 & 8,
                            proximityBonding=False,
                            removeHs=False)
    for atom in m.GetAtoms():
        atom.SetNoImplicit(True)

    return assign_pdb_bonds(m, templates)


def gemmi_to_rdkit(block: gemmi.cif.Block):
    """Convert a gemmi Block to an RDKit mol

    used for PDBx loading while RDKit lacks support
    """
    pt = Chem.GetPeriodicTable()

    m = Chem.Mol()
    em = AllChem.EditableMol(m)
    pos = []

    for row in block.find('_atom_site.',
                      ['group_PDB', 'type_symbol', 'label_atom_id',
                       'label_comp_id',
                       'label_asym_id', 'label_seq_id', 'pdbx_PDB_ins_code',
                       'Cartn_x', 'Cartn_y', 'Cartn_z']):
        het_status, elem, name, resname, chainid, resnum, icode, posx, posy, posz = row

        if icode == '?':
            icode = ''

        a = Chem.Atom(pt.GetAtomicNumber(elem))

        mi = Chem.AtomPDBResidueInfo()
        mi.SetIsHeteroAtom(het_status.startswith('HETATM'))
        mi.SetName(name)
        mi.SetResidueName(resname)
        mi.SetChainId(chainid)
        mi.SetResidueNumber(int(resnum))
        mi.SetInsertionCode(icode)
        a.SetMonomerInfo(mi)

        em.AddAtom(a)
        pos.append((float(posx), float(posy), float(posz)))

    conf = Chem.Conformer(len(pos))
    conf.Set3D(True)
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, p)

    m = em.GetMol()
    m.AddConformer(conf)

    return m


def load_pdbx_file(pdbx_path: Union[str, os.PathLike],
                   templates: list[gemmi.cif.Document]) -> Chem.Mol:
    """Load a pdbx file and assign bond information from templates

    Parameters
    ----------
    pdbx_path: str or pathlib.Path
        path to the pdbx file
    templates : list of gemmi.cif.Document
        the templates to be applied

    Returns
    -------
    model : Chem.Mol

    Raises
    ------
    RuntimeError
      if there is more than one block in the file
    """
    # first load pdbx using gemmi
    d = gemmi.cif.read_file(str(pdbx_path))
    block = d.sole_block()

    # transfer gemmi object into rdkit
    m = gemmi_to_rdkit(block)

    return assign_pdb_bonds(m, templates)


def doc_contains(d: gemmi.cif.Document, label: str) -> bool:
    # hack while gemmi release not pushed out yet
    try:
        _ = d[label]
    except KeyError:
        return False
    else:
        return True


def get_atom_named(mol, i, j, resname, atomname) -> int:
    aliases = ATOMNAME_ALIASES[resname]

    for x in range(i, j):
        nm = mol.GetAtomWithIdx(x).GetMonomerInfo().GetName().strip()
        if aliases.get(nm, nm) == atomname:
            return x
    else:
        raise ValueError


def smallest_connection(mol, a, b, x, y, conf) -> tuple[int, int, float]:
    # todo: might this return a SS bond instead of the backbond bond expected?
    #       if so, return list of all possible bonds rather than single bond?
    mindist = float('inf')
    ret = -1, -1

    for i in range(a, b):
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 1:
            continue  # skip hydrogens
        posi = conf.GetAtomPosition(i)
        for j in range(x, y):
            if mol.GetAtomWithIdx(j).GetAtomicNum() == 1:
                continue

            d = posi.Distance(conf.GetAtomPosition(j))

            if d > mindist:
                continue
            mindist = d
            ret = i, j

    return *ret, mindist


def assign_inter_residue_bonds(mol):
    bonds = []
    conf = mol.GetConformer()

    logger.debug('trying to add inter residue bonds')

    valence = np.zeros(mol.GetNumAtoms(), dtype=int)

    res_gen = residue_spans(mol)
    res_A: Residue
    res_B: Residue
    i_A, j_A, res_A = next(res_gen)
    for i_B, j_B, res_B in res_gen:
        logger.debug(f'trying {res_A.resname}({i_A} {j_A}) '
                     f'{res_B.resname}({i_B} {j_B})')
        if res_A.chain != res_B.chain:
            logger.debug('skipping due to chain break '
                         f'(chain id {res_A.chain} vs {res_B.chain})')
            # new chain so skip adding bonds from this residue
            i_A, j_A, res_A = i_B, j_B, res_B
            continue

        canonical_resname_A = RESNAME_ALIASES.get(res_A.resname, res_A.resname)
        canonical_resname_B = RESNAME_ALIASES.get(res_B.resname, res_B.resname)
        if canonical_resname_A in STANDARD_RESNAMES and canonical_resname_B in STANDARD_RESNAMES:
            N_A = get_atom_named(mol, i_A, j_A, canonical_resname_A, 'N')
            C_A = get_atom_named(mol, i_A, j_A, canonical_resname_A, 'C')
            N_B = get_atom_named(mol, i_B, j_B, canonical_resname_B, 'N')
            C_B = get_atom_named(mol, i_B, j_B, canonical_resname_B, 'C')

            d1 = conf.GetAtomPosition(C_A).Distance(conf.GetAtomPosition(N_B))
            d2 = conf.GetAtomPosition(N_A).Distance(conf.GetAtomPosition(C_B))
            logger.debug(f"C-N distance {d1}")
            logger.debug(f"N-C distance {d2}")
            if d1 < MAX_AMIDE_LENGTH:
                logger.debug(f'accepting C-N between {C_A} {N_B}')
                bonds.append((C_A, N_B))
            elif d2 < MAX_AMIDE_LENGTH:
                logger.debug(f'accepting N-C between {N_A} {C_B}')
                bonds.append((N_A, C_B))
            else:
                logger.debug("Failed to find inter residue bond between "
                             f"{res_A} and {res_B}, assuming chain break")
        else:
            # find nonstandard residue bond
            i, j, d = smallest_connection(mol, i_A, j_A, i_B, j_B, conf)
            if d < MAX_AMIDE_LENGTH:
                logger.debug(f'accepting nonstandard bond {i}-{j} {d}')
                bonds.append((i, j))

        i_A, j_A, res_A = i_B, j_B, res_B

    logger.debug(f'found {len(bonds)} bonds')

    em = AllChem.EditableMol(mol)
    for i, j in bonds:
        logger.debug(f'adding inter bond {i} {j}')

        valence[i] += 1
        valence[j] += 1
        em.AddBond(i, j, order=Chem.BondType.SINGLE)

    return em.GetMol(), valence


def assign_sulphur_bonds(mol):
    sulphurs = [i for i in range(mol.GetNumAtoms())
                if mol.GetAtomWithIdx(i).GetAtomicNum() == 16]

    valence = np.zeros(mol.GetNumAtoms(), dtype=int)

    conf = mol.GetConformer()
    bonds = []
    for i, j in itertools.combinations(sulphurs, 2):
        di = conf.GetAtomPosition(i)
        if di.Distance(conf.GetAtomPosition(j)) < MAX_DISULPHIDE_LENGTH:
            bonds.append((i, j))

    em = AllChem.EditableMol(mol)

    for i, j in bonds:
        em.AddBond(i, j, order=Chem.BondType.SINGLE)
        valence[i] += 1
        valence[j] += 1

    return em.GetMol(), valence


def assign_charge(mol: Chem.Mol, valence) -> Chem.Mol:
    # hacked adaptation of rdkit's xyz2mol
    for atom, v in zip(mol.GetAtoms(), valence):
        atnum = atom.GetAtomicNum()
        if atnum == 1:
            continue

        if atnum == 5:
            chg = 3 - v
        elif atnum == 15 and v == 5:
            chg = 0
        elif atnum == 16 and v == 6:
            chg = 0
        else:
            chg = PT.GetNOuterElecs(atnum) - 8 + v

        atom.SetFormalCharge(int(chg))  # rdkit hates np.int64

    return mol


def histidine_check(mol, i, j) -> str:
    # returns either HIS (i.e. HIE) or HID, (HIP is covered by HIS)
    # these templates differ on where the double bond (and charge) is placed
    has_HD1 = False
    has_HE2 = False
    for x in range(i, j):
        nm = mol.GetAtomWithIdx(x).GetMonomerInfo().GetName().strip()

        has_HD1 |= nm == 'HD1'
        has_HE2 |= nm == 'HE2'

    if has_HD1 and not has_HE2:
        return 'HID'
    else:
        return 'HIS'


def assign_pdb_bonds(mol: Chem.Mol, templates: list[gemmi.cif.Document]) -> Chem.Mol:
    """Take an RDKit Molecule and assign aromaticity and bonds

    Requires that the bonds have the MonomerInfo property populated

    Parameters
    ----------
    mol : Chem.Mol
      the input molecule, this is modified in place and returned
    templates : list[gemmi.cif.Document]
      the template molecules from which to assign bonds

    Note
    ----
    Any bonds on the input molecule will be removed at the start of the process
    """
    mol = strip_bonds(mol)

    # 1) assign properties inside each Residue
    valence = np.zeros(mol.GetNumAtoms(), dtype=int)
    for i, j, res in residue_spans(mol):
        resname = RESNAME_ALIASES.get(res.resname, res.resname)
        # for ions, we can let this slide
        if (j - i) == 1:
            logger.debug(f"Skipping presumed ion: {res}")
            continue
        elif resname == 'HIS':
            resname = histidine_check(mol, i, j)

        for t in templates:
            # check if resname in template doc
            if not doc_contains(t, resname):
                continue
            mol, v = assign_intra_props(mol, range(i, j), t[resname])
            valence += v
            break  # avoid double assignment, ordering of templates is relevant
        else:
            raise ValueError(f"Failed to find template for residue: '{res}'")

    # 2) assign bonds between residues
    mol, v = assign_inter_residue_bonds(mol)
    valence += v

    # 3) sulphur bridges
    mol, v = assign_sulphur_bonds(mol)
    valence += v

    # 4) ch ch ch ch charges
    mol = assign_charge(mol, valence)

    return mol
