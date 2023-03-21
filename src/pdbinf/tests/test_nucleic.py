# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf

"""
More involved testing of dna/rna support
"""

import pdbinf


def test_nucleic(openmm_nucleic):
    m = pdbinf.load_pdb_file(openmm_nucleic, templates=[pdbinf.DNA_DOC, pdbinf.RNA_DOC])

    assert m.GetNumAtoms() == 767
    # {'U': 3, 'T': 1, 'A': 5, 'C': 6, 'G': 9} -> 59 total rings
    assert m.GetNumBonds() == 827
