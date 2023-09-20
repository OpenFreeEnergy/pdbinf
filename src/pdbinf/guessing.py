# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
"""pdbinf.guessing

Tools for guessing residue and atom names from malformed inputs.
"""
from collections import defaultdict, deque
import gemmi
import itertools
from rdkit import Chem
from rdkit.Chem import rdMolHash, rdmolops
from typing import Optional


PT = Chem.GetPeriodicTable()


def copy_mol_subset(mol: Chem.Mol, idx: list[int]) -> Chem.Mol:
    """Extract a copy of a subset of a larger molecule"""
    sub = Chem.Mol()
    em = Chem.EditableMol(sub)

    old2new = dict()

    bonds = []

    for i, index in enumerate(idx):
        old2new[index] = i

        src = mol.GetAtomWithIdx(index)

        em.AddAtom(src)

        for b in src.GetBonds():
            other = b.GetOtherAtomIdx(index)

            if other < index:
                continue
            if other not in idx:
                continue

            bonds.append((index, other, b.GetBondType()))

    for i, j, o in bonds:
        i2, j2 = old2new[i], old2new[j]

        em.AddBond(i2, j2, o)

    return em.GetMol()


def block_to_mol(block: gemmi.cif.Block) -> Chem.Mol:
    """Convert a gemmi template Block into an rdkit Mol"""
    id_2_index = dict()

    m = Chem.Mol()
    em = Chem.EditableMol(m)

    for i, (elem, atomid) in enumerate(block.find('_chem_comp_atom.',
                                                  ['type_symbol', 'atom_id'])):
        if elem == 'H':
            continue

        id_2_index[atomid] = i

        a = Chem.Atom(PT.GetAtomicNumber(elem))
        em.AddAtom(a)

    for id1, id2, o, arom in block.find('_chem_comp_bond.',
                                        ['atom_id_1', 'atom_id_2',
                                         'value_order', 'pdbx_aromatic_flag']):
        try:
            index1, index2 = id_2_index[id1], id_2_index[id2]
        except KeyError:
            continue

        if arom == 'Y':
            rdko = Chem.BondType.AROMATIC
        else:
            rdko = {'SING': Chem.BondType.SINGLE,
                    'DOUB': Chem.BondType.DOUBLE}[o]

        em.AddBond(index1, index2, order=rdko)

    m = em.GetMol()
    m.UpdatePropertyCache()
    return m


def block_to_ix(block):
    """Generator of possible heavy atom indices representing combinations of leaving atoms

    Multiple molecules are returned representing different combinations of leaving molecules
    """
    m = block_to_mol(block)

    # find all leaving atoms
    heavy_count = 0
    leaving_atoms = []
    for elem, leave in block.find('_chem_comp_atom.', ['type_symbol', 'pdbx_leaving_atom_flag']):
        if elem == 'H':
            continue

        if leave == 'Y':
            leaving_atoms.append(heavy_count)
        heavy_count += 1

    static_atoms = [i for i, atom in enumerate(m.GetAtoms())
                    if i not in leaving_atoms]

    # iterate over combinations of leaving atoms true/false
    for includeX in itertools.product((True, False), repeat=len(leaving_atoms)):
        this_mol_ix = static_atoms.copy()
        for keep, ix in zip(includeX, leaving_atoms):
            if keep:
                this_mol_ix.append(ix)
        this_mol_ix.sort()

        yield this_mol_ix


def block_to_molhash(block: gemmi.cif.Block):
    """Convert a block to its possible element graph hashes

    Multiple possible hashes are returned based on the combination of possible leaving atoms
    """
    whole_mol = block_to_mol(block)
    for ix in block_to_ix(block):
        mol = copy_mol_subset(whole_mol, ix)

        yield rdMolHash.MolHash(mol, rdMolHash.HashFunction.ElementGraph)


def normalised_hash(mol: Chem.Mol):
    target_mol = Chem.RemoveAllHs(mol, sanitize=False)
    Chem.RemoveStereochemistry(target_mol)
    return rdMolHash.MolHash(target_mol, rdMolHash.HashFunction.ElementGraph)


def guess_residue_name(mol: Chem.Mol, templates: list[gemmi.cif.Block],
                       subset: Optional[list[int]] = None) -> str:
    """Guess the Residue name for *mol* from within *templates*"""
    if subset is not None:
        mol = copy_mol_subset(mol, subset)
    target = normalised_hash(mol)

    for t in templates:
        for h in block_to_molhash(t):
            if h == target:
                return t.name
    else:
        raise ValueError("Could not match any template against mol")


def block_to_query(block: gemmi.cif.Block, ix: list[int]) -> Chem.Mol:
    """Generate a query molecule from a gemmi Block"""
    id_2_index = dict()

    pieces = []

    heavy_atom_count = 0
    for i, (elem, atomid) in enumerate(block.find('_chem_comp_atom.',
                                                  ['type_symbol', 'atom_id'])):
        if elem == 'H':
            continue
        heavy_atom_count += 1
        if (heavy_atom_count - 1) not in ix:
            continue

        id_2_index[atomid] = i

        piece = f'[#{PT.GetAtomicNumber(elem)}:{i}]'
        pieces.append(piece)

    m = Chem.MolFromSmarts('.'.join(pieces))
    em = Chem.EditableMol(m)

    for id1, id2 in block.find('_chem_comp_bond.', ['atom_id_1', 'atom_id_2']):
        try:
            index1, index2 = id_2_index[id1], id_2_index[id2]
        except KeyError:
            continue

        em.AddBond(index1, index2, order=Chem.BondType.UNSPECIFIED)

    m = em.GetMol()
    m.UpdatePropertyCache()
    return m


def guess_atom_names(mol: Chem.Mol, template: gemmi.cif.Block,
                     subset: Optional[list[int]]) -> list[str]:
    """Get the canonical atom names for *mol* based on *template*

    Parameters
    ----------
    mol : Chem.Mol
      target molecule to guess atom names for, use subset to specify a subsection of this mol
    template : gemmi.cif.Block
      template to extract atom names from
    subset : optional list[int]
      atom indices that indicate a single monomer to guess for

    Returns
    -------
    names : list[str]
      the atom names for each atom
    """
    if subset is not None:
        mol = copy_mol_subset(mol, subset)
    target = normalised_hash(mol)
    whole_block_mol = block_to_mol(template)

    for ix in block_to_ix(template):
        # find appropriate template among possible leaving atoms
        subset_block_mol = copy_mol_subset(whole_block_mol, ix)

        if rdMolHash.MolHash(subset_block_mol, rdMolHash.HashFunction.ElementGraph) != target:
            continue

        q = block_to_query(template, ix)

        m = mol.GetSubstructMatch(q)

        if not m:
            raise ValueError

        # mapping of heavy atom names to molecule indices
        names = [r.str(0) for r in template.find('_chem_comp_atom.', ['atom_id'])]
        # rearrange to match input order and return
        heavy_atom_names = [names[i] for i in m]
        idx_to_name = {i: n for i, n in zip(m, heavy_atom_names)}

        # first identify hydrogen names
        hydrogen_names = {nm for elem, nm in template.find('_chem_comp_atom.', ['type_symbol', 'atom_id'])
                          if elem == 'H'}
        # create list of hydrogens per heavy atom
        hydrogens_per_heavy = defaultdict(deque)
        for i, j in template.find('_chem_comp_bond.', ['atom_id_1', 'atom_id_2']):
            if i in hydrogen_names:
                hydrogens_per_heavy[j].append(i)
            elif j in hydrogen_names:
                hydrogens_per_heavy[i].append(j)

        output_names = []
        for atom in mol.GetAtoms():
            try:
                nm = idx_to_name[atom.GetIdx()]
            except KeyError:
                # figure out which heavy atom we're attached to
                hvy_idx = atom.GetBonds()[0].GetOtherAtomIdx(atom.GetIdx())
                hvy_nm = idx_to_name[hvy_idx]

                # grab the next hydrogen name off that heavy atom
                nm = hydrogens_per_heavy[hvy_nm].popleft()

            output_names.append(nm)

        return output_names
    else:
        raise ValueError("Failed to find correct subset of atoms to match template to mol")
