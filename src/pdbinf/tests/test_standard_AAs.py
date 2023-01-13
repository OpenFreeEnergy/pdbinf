import gemmi
from rdkit import Chem

import pdbinf

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


arg = """\
REMARK   1 PDBFIXER FROM: MainChain_ARG.pdb                                                                                                 
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-13                                                                                              
CRYST1   48.000   48.000   48.000  90.00  90.00  90.00 P 1           1                                                                      
HETATM    1  H1  ACE     1      26.138  25.214  24.583  1.00  0.00           H                                                              
HETATM    2  CH3 ACE     1      26.124  26.209  24.139  1.00  0.00           C                                                              
HETATM    3  H2  ACE     1      25.413  26.835  24.676  1.00  0.00           H                                                              
HETATM    4  H3  ACE     1      25.837  26.135  23.092  1.00  0.00           H                                                              
HETATM    5  C   ACE     1      27.506  26.811  24.243  1.00  0.00           C                                                              
HETATM    6  O   ACE     1      28.418  26.182  24.765  1.00  0.00           O                                                              
ATOM      7  N   ARG     2      27.671  28.039  23.742  1.00  0.00           N                                                              
ATOM      8  H   ARG     2      26.871  28.490  23.320  1.00  0.00           H                                                              
ATOM      9  CA  ARG     2      28.942  28.780  23.768  1.00  0.00           C                                                              
ATOM     10  HA  ARG     2      29.476  28.515  24.682  1.00  0.00           H                                                              
ATOM     11  CB  ARG     2      29.787  28.342  22.555  1.00  0.00           C                                                              
ATOM     12  HB2 ARG     2      29.888  27.255  22.575  1.00  0.00           H                                                              
ATOM     13  HB3 ARG     2      29.263  28.617  21.637  1.00  0.00           H                                                              
ATOM     14  CG  ARG     2      31.194  28.957  22.536  1.00  0.00           C                                                              
ATOM     15  HG2 ARG     2      31.121  30.043  22.463  1.00  0.00           H                                                              
ATOM     16  HG3 ARG     2      31.716  28.692  23.456  1.00  0.00           H                                                              
ATOM     17  CD  ARG     2      31.981  28.424  21.331  1.00  0.00           C                                                              
ATOM     18  HD2 ARG     2      32.042  27.336  21.408  1.00  0.00           H                                                              
ATOM     19  HD3 ARG     2      31.438  28.683  20.419  1.00  0.00           H                                                              
ATOM     20  NE  ARG     2      33.344  28.991  21.272  1.00  0.00           N                                                              
ATOM     21  HE  ARG     2      33.600  29.629  22.006  1.00  0.00           H                                                              
ATOM     22  CZ  ARG     2      34.256  28.739  20.348  1.00  0.00           C                                                              
ATOM     23  NH1 ARG     2      34.029  27.928  19.353  1.00  0.00           N                                                              
ATOM     24 HH11 ARG     2      33.132  27.480  19.289  1.00  0.00           H                                                              
ATOM     25 HH12 ARG     2      34.732  27.748  18.660  1.00  0.00           H                                                              
ATOM     26  NH2 ARG     2      35.428  29.305  20.406  1.00  0.00           N                                                              
ATOM     27 HH21 ARG     2      35.650  29.938  21.155  1.00  0.00           H                                                              
ATOM     28 HH22 ARG     2      36.117  29.111  19.703  1.00  0.00           H                                                              
ATOM     29  C   ARG     2      28.683  30.288  23.768  1.00  0.00           C                                                              
ATOM     30  O   ARG     2      28.008  30.786  22.875  1.00  0.00           O                                                              
HETATM   31  N   NME     3      29.250  31.004  24.744  1.00  0.00           N                                                              
HETATM   32  H   NME     3      29.758  30.508  25.458  1.00  0.00           H                                                              
HETATM   33  C   NME     3      29.101  32.453  24.889  1.00  0.00           C                                                              
HETATM   34  H1  NME     3      28.055  32.699  25.084  1.00  0.00           H                                                              
HETATM   35  H2  NME     3      29.712  32.817  25.716  1.00  0.00           H                                                              
HETATM   36  H3  NME     3      29.407  32.954  23.968  1.00  0.00           H                                                              
TER      37      NME     3                                                                                                                  
CONECT    1    2                                                                                                                            
CONECT    2    5    1    3    4                                                                                                             
CONECT    3    2                                                                                                                            
CONECT    4    2                                                                                                                            
CONECT    5    2    6    7                                                                                                                  
CONECT    6    5                                                                                                                            
CONECT    7    5                                                                                                                            
CONECT   29   31                                                                                                                            
CONECT   31   29   33   32                                                                                                                  
CONECT   32   31                                                                                                                            
CONECT   33   34   35   36   31                                                                                                             
CONECT   34   33                                                                                                                            
CONECT   35   33                                                                                                                            
CONECT   36   33                                                                                                                            
END 
"""


def test_ala_ala():
    m = Chem.MolFromPDBBlock(ala_ala, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 28
    assert m.GetNumBonds() == 27
    assert m.GetAtomWithIdx(0).GetFormalCharge() == 1


def test_arg():
    m = Chem.MolFromPDBBlock(arg, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 36
    assert m.GetNumBonds() == 35
    # atom 25 is assigned the double bond, so should have FC +1
    atom = m.GetAtomWithIdx(25)
    assert atom.GetFormalCharge() == 1
    assert any(b.GetBondType() == Chem.BondType.DOUBLE
               for b in atom.GetBonds())