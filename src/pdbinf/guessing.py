# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
"""pdbinf.guessing

Tools for guessing residue and atom names from malformed inputs.
"""
import gemmi
from rdkit import Chem
from rdkit.Chem import rdMolHash


PT = Chem.GetPeriodicTable()


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


def block_to_molhash(block: gemmi.cif.Block) -> str:
    """Convert a block to its element graph hash"""
    m = block_to_mol(block)

    return rdMolHash.MolHash(rdMolHash.HashFunction.ElementGraph)


def guess_residue_name(mol: Chem.Mol, templates: list[gemmi.cif.Block]) -> str:
    """Guess the Residue name for *mol* from within *templates*"""
    target = rdMolHash.MolHash(mol, rdMolHash.HashFunction.ElementGraph)

    for t in templates:
        h = block_to_molhash(t)

        if h == target:
            return t.name
    else:
        raise ValueError


def block_to_query(block: gemmi.cif.Block) -> Chem.Mol:
    """Generate a query molecule from a gemmi Block"""
    id_2_index = dict()

    pieces = []

    for i, (elem, atomid) in enumerate(block.find('_chem_comp_atom.',
                                                  ['type_symbol', 'atom_id'])):
        if elem == 'H':
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


def guess_atom_names(mol: Chem.Mol, template: gemmi.cif.Block) -> list[str]:
    """Get the canonical atom names for *mol* based on *template*"""
    q = block_to_query(template)

    m = mol.GetSubstructMatch(q)

    if not m:
        raise ValueError

    names = [r.str(0) for r in template.find('_chem_comp_atom.', ['atom_id'])]

    # rearrange to match input order and return
    return [names[i] for i in m]
