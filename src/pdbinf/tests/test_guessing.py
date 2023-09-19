import pytest
from rdkit import Chem

import pdbinf
from pdbinf import guessing as g


def test_mol_subset():
    m1 = Chem.MolFromSmiles('c1ccccc1CCO')

    m2 = g.copy_mol_subset(m1, [6, 7, 8])

    assert m2.GetNumAtoms() == 3
    assert m2.GetNumBonds() == 2
    assert Chem.MolToSmiles(m2) == 'CCO'
    # check ordering preserved
    assert m2.GetAtomWithIdx(2).GetAtomicNum() == 8
