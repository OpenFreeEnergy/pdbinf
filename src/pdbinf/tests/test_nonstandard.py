import gemmi
import pdbinf
from rdkit import Chem
import pytest


tpo = """\
REMARK   1 PDBFIXER FROM: tpo.pdb
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-19
HETATM    1  H1  ACE     1      -1.861   4.051   1.454  1.00  0.00           H  
HETATM    2  CH3 ACE     1      -1.252   3.158   1.245  1.00  0.00           C  
HETATM    3  H2  ACE     1      -1.461   2.423   2.061  1.00  0.00           H  
HETATM    4  H3  ACE     1      -0.170   3.423   1.346  1.00  0.00           H  
HETATM    5  C   ACE     1      -1.522   2.585  -0.102  1.00  0.00           C  
HETATM    6  O   ACE     1      -2.666   2.764  -0.584  1.00  0.00           O  
HETATM    7  N   TPO     2      -0.365   1.748  -0.928  1.00  0.00           N  
HETATM    8  CA  TPO     2      -0.159   0.351  -0.158  1.00  0.00           C  
HETATM    9  CB  TPO     2       1.315  -0.018   0.021  1.00  0.00           C  
HETATM   10  CG2 TPO     2       1.942   0.984   0.953  1.00  0.00           C  
HETATM   11  OG1 TPO     2       1.924  -0.131  -1.190  1.00  0.00           O  
HETATM   12  P   TPO     2       2.522  -1.636  -1.565  1.00  0.00           P  
HETATM   13  O1P TPO     2       2.285  -2.666  -0.492  1.00  0.00           O  
HETATM   14  O2P TPO     2       4.193  -1.597  -1.883  1.00  0.00           O  
HETATM   15  O3P TPO     2       1.769  -2.179  -2.981  1.00  0.00           O  
HETATM   16  C   TPO     2      -0.881  -0.709  -0.837  1.00  0.00           C  
HETATM   17  O   TPO     2      -1.555  -0.570  -1.856  1.00  0.00           O  
HETATM   18  H   TPO     2       0.320   2.242  -1.638  1.00  0.00           H  
HETATM   19  HA  TPO     2      -0.554   0.439   0.896  1.00  0.00           H  
HETATM   20  HB  TPO     2       1.344  -1.004   0.536  1.00  0.00           H  
HETATM   21 HG21 TPO     2       2.873   0.507   1.359  1.00  0.00           H  
HETATM   22 HG22 TPO     2       2.124   1.962   0.496  1.00  0.00           H  
HETATM   23 HG23 TPO     2       1.262   1.078   1.828  1.00  0.00           H  
HETATM   24  N   NME     3      -0.837  -2.343  -0.220  1.00  0.00           N  
HETATM   25  H   NME     3      -0.418  -3.247  -1.189  1.00  0.00           H  
HETATM   26  C   NME     3      -1.923  -2.794   0.887  1.00  0.00           C  
HETATM   27  H1  NME     3      -2.865  -1.683   1.480  1.00  0.00           H  
HETATM   28  H2  NME     3      -3.399  -3.145  -0.439  1.00  0.00           H  
HETATM   29  H3  NME     3      -1.984  -3.991   1.501  1.00  0.00           H  
TER      30      NME     3
END
"""


def test_tpo(tpo_doc):
    m = Chem.MolFromPDBBlock(tpo, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC, tpo_doc])

    assert m.GetNumAtoms() == 29
    assert m.GetNumBonds() == 28
    for i, atom in enumerate(m.GetAtoms()):
        if i in [13, 14]:  # this is P-[O-]
            assert atom.GetFormalCharge() == -1
        else:
            assert atom.GetFormalCharge() == 0

        if i == 12:  # this is the P=O
            b = atom.GetBonds()

            assert len(b) == 1
            assert b[0].GetBondType() == Chem.BondType.DOUBLE


def test_tpo_fail():
    m = Chem.MolFromPDBBlock(tpo, proximityBonding=False, removeHs=False)

    with pytest.raises(ValueError, match="Failed to find template"):
        m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

