import gemmi
from rdkit import Chem

import alibaby

ala_ala = """\
REMARK   1 PDBFIXER FROM: NTerminal_ALA_ALA.pdb
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-03
CRYST1   48.000   48.000   48.000  90.00  90.00  90.00 P 1           1 
ATOM      1  N   ALA     1      27.352  25.519  24.270  1.00  0.00           N  
ATOM      2  H   ALA     1      27.682  24.857  23.580  1.00  0.00           H  
ATOM      3  H2  ALA     1      26.346  25.594  24.227  1.00  0.00           H  
ATOM      4  H3  ALA     1      27.635  25.175  25.179  1.00  0.00           H  
ATOM      5  CA  ALA     1      27.988  26.835  24.037  1.00  0.00           C  
ATOM      6  HA  ALA     1      27.733  27.179  23.034  1.00  0.00           H  
ATOM      7  CB  ALA     1      27.480  27.872  25.050  1.00  0.00           C  
ATOM      8  HB1 ALA     1      27.750  27.571  26.064  1.00  0.00           H  
ATOM      9  HB2 ALA     1      27.924  28.846  24.840  1.00  0.00           H  
ATOM     10  HB3 ALA     1      26.396  27.963  24.978  1.00  0.00           H  
ATOM     11  C   ALA     1      29.507  26.678  24.103  1.00  0.00           C  
ATOM     12  O   ALA     1      29.961  25.613  24.504  1.00  0.00           O  
ATOM     13  N   ALA     2      30.269  27.698  23.708  1.00  0.00           N  
ATOM     14  H   ALA     2      29.842  28.573  23.431  1.00  0.00           H  
ATOM     15  CA  ALA     2      31.734  27.734  23.756  1.00  0.00           C  
ATOM     16  HA  ALA     2      32.086  27.184  24.631  1.00  0.00           H  
ATOM     17  CB  ALA     2      32.290  27.061  22.492  1.00  0.00           C  
ATOM     18  HB1 ALA     2      31.957  27.600  21.603  1.00  0.00           H  
ATOM     19  HB2 ALA     2      33.380  27.063  22.518  1.00  0.00           H  
ATOM     20  HB3 ALA     2      31.947  26.027  22.437  1.00  0.00           H  
ATOM     21  C   ALA     2      32.210  29.196  23.885  1.00  0.00           C  
ATOM     22  O   ALA     2      31.403  30.109  23.706  1.00  0.00           O  
HETATM   23  N   NME     3      33.495  29.410  24.192  1.00  0.00           N  
HETATM   24  H   NME     3      34.098  28.612  24.315  1.00  0.00           H  
HETATM   25  C   NME     3      34.087  30.739  24.339  1.00  0.00           C  
HETATM   26  H1  NME     3      33.505  31.332  25.048  1.00  0.00           H  
HETATM   27  H2  NME     3      35.114  30.659  24.701  1.00  0.00           H  
HETATM   28  H3  NME     3      34.088  31.254  23.376  1.00  0.00           H  
TER      29      NME     3
CONECT   21   23
CONECT   23   21   25   24
CONECT   24   23
CONECT   25   26   27   28   23
CONECT   26   25
CONECT   27   25
CONECT   28   25
END
"""


def test_ala_ala():
    m = Chem.MolFromPDBBlock(ala_ala, proximityBonding=False, removeHs=False)

    templates = [
        gemmi.cif.read('/home/richard/code/alibaby/data/ala.cif')
    ]

    m = alibaby.assign_pdb_bonds(m, templates=templates)

    assert m.GetNumAtoms() == 28
    assert m.GetNumBonds() == 27
    assert m.GetAtomWithIdx(0).GetFormalCharge() == 1
