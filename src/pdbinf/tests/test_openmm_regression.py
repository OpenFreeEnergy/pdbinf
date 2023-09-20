# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
import pdbinf
import pytest

app = pytest.importorskip('openmm.app')


def assert_bonds_match(ref, other):
    # check that all bonds in ref (an openmm PDBFile) are in other (an rdkit mol)
    i = 0
    for b in ref.topology.bonds():
        ix1, ix2 = b.atom1.index, b.atom2.index

        x = other.GetBondBetweenAtoms(ix1, ix2)

        assert x is not None, (ix1, ix2)
        i += 1

    # checks that other cannot contain any additional bonds
    assert i == other.GetNumBonds()


def test_openmm_regression(plb_proteins, tpo_doc):
    ref = app.PDBFile(plb_proteins)

    if 'cdk2' in plb_proteins:
        t = [pdbinf.STANDARD_AA_DOC, tpo_doc]
    else:
        t = [pdbinf.STANDARD_AA_DOC]

    m = pdbinf.load_pdb_file(plb_proteins, templates=t)

    assert_bonds_match(ref, m)


def test_nucleic_regression(openmm_nucleic):
    ref = app.PDBFile(openmm_nucleic)

    m = pdbinf.load_pdb_file(openmm_nucleic, templates=[pdbinf.RNA_DOC, pdbinf.DNA_DOC])

    assert_bonds_match(ref, m)