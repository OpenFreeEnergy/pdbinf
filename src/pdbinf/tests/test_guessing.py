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


class TestCDK2Guessing:
    @pytest.fixture
    def mol_tpo_idx(self, cdk2_pdb):
        m = Chem.MolFromPDBFile(cdk2_pdb,
                                flavor=1 & 8,
                                proximityBonding=True,
                                removeHs=False,
                                )

        tpo_idx = [a.GetIdx() for a in m.GetAtoms()
                   if a.GetMonomerInfo().GetResidueName() == 'TPO']

        return m, tpo_idx

    def test_guess_residue(self, mol_tpo_idx, tpo_doc):
        m, tpo_idx = mol_tpo_idx

        assert g.guess_residue_name(m, tpo_doc, subset=tpo_idx) == 'TPO'

    def test_guess_atom_names(self, mol_tpo_idx, tpo_doc):
        m, tpo_idx = mol_tpo_idx

        names = g.guess_atom_names(m, tpo_doc['TPO'],
                                   subset=tpo_idx)

        assert len(names) == len(tpo_idx)
        assert names == ['N', 'CA', 'CB', 'CG2', 'OG1', 'P', 'O1P', 'O2P', 'O3P', 'C', 'O',
                         'H', 'HA', 'HB', 'HG21', 'HG22', 'HG23']
