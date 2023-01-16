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
END 
"""


trp = """\
REMARK   1 PDBFIXER FROM: MainChain_TRP.pdb                                                                                                 
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-13                                                                                              
CRYST1   48.000   48.000   48.000  90.00  90.00  90.00 P 1           1                                                                      
HETATM    1  H1  ACE     1      25.987  25.343  25.068  1.00  0.00           H                                                              
HETATM    2  CH3 ACE     1      26.037  26.135  24.322  1.00  0.00           C                                                              
HETATM    3  H2  ACE     1      25.187  26.804  24.443  1.00  0.00           H                                                              
HETATM    4  H3  ACE     1      26.028  25.698  23.325  1.00  0.00           H                                                              
HETATM    5  C   ACE     1      27.318  26.912  24.518  1.00  0.00           C                                                              
HETATM    6  O   ACE     1      28.091  26.607  25.417  1.00  0.00           O                                                              
ATOM      7  N   TRP     2      27.538  27.920  23.677  1.00  0.00           N                                                              
ATOM      8  H   TRP     2      26.840  28.135  22.980  1.00  0.00           H                                                              
ATOM      9  CA  TRP     2      28.712  28.797  23.720  1.00  0.00           C                                                              
ATOM     10  HA  TRP     2      29.022  28.915  24.759  1.00  0.00           H                                                              
ATOM     11  CB  TRP     2      29.868  28.144  22.946  1.00  0.00           C                                                              
ATOM     12  HB2 TRP     2      30.130  27.203  23.433  1.00  0.00           H                                                              
ATOM     13  HB3 TRP     2      29.522  27.909  21.938  1.00  0.00           H                                                              
ATOM     14  CG  TRP     2      31.112  28.973  22.840  1.00  0.00           C                                                              
ATOM     15  CD1 TRP     2      31.621  29.485  21.697  1.00  0.00           C                                                              
ATOM     16  HD1 TRP     2      31.194  29.338  20.711  1.00  0.00           H                                                              
ATOM     17  NE1 TRP     2      32.746  30.233  21.984  1.00  0.00           N                                                              
ATOM     18  HE1 TRP     2      33.294  30.707  21.281  1.00  0.00           H                                                              
ATOM     19  CE2 TRP     2      33.010  30.254  23.337  1.00  0.00           C                                                              
ATOM     20  CZ2 TRP     2      33.997  30.881  24.111  1.00  0.00           C                                                              
ATOM     21  HZ2 TRP     2      34.755  31.491  23.642  1.00  0.00           H                                                              
ATOM     22  CH2 TRP     2      33.987  30.703  25.505  1.00  0.00           C                                                              
ATOM     23  HH2 TRP     2      34.741  31.178  26.118  1.00  0.00           H                                                              
ATOM     24  CZ3 TRP     2      32.998  29.900  26.102  1.00  0.00           C                                                              
ATOM     25  HZ3 TRP     2      32.999  29.759  27.174  1.00  0.00           H                                                              
ATOM     26  CE3 TRP     2      32.012  29.273  25.314  1.00  0.00           C                                                              
ATOM     27  HE3 TRP     2      31.266  28.648  25.784  1.00  0.00           H                                                              
ATOM     28  CD2 TRP     2      31.991  29.437  23.912  1.00  0.00           C                                                              
ATOM     29  C   TRP     2      28.381  30.195  23.175  1.00  0.00           C                                                              
ATOM     30  O   TRP     2      27.456  30.346  22.378  1.00  0.00           O                                                              
HETATM   31  N   NME     3      29.137  31.215  23.603  1.00  0.00           N                                                              
HETATM   32  H   NME     3      29.911  30.999  24.215  1.00  0.00           H                                                              
HETATM   33  C   NME     3      28.956  32.603  23.187  1.00  0.00           C                                                              
HETATM   34  H1  NME     3      27.915  32.898  23.332  1.00  0.00           H                                                              
HETATM   35  H2  NME     3      29.602  33.258  23.773  1.00  0.00           H                                                              
HETATM   36  H3  NME     3      29.204  32.703  22.129  1.00  0.00           H                                                              
TER      37      NME     3                                                                                                                                                                                                                                             
END                                                                                                                                         
"""


pro = """\
       REMARK   1 PDBFIXER FROM: MainChain_PRO.pdb                                                                                                 
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-14                                                                                              
CRYST1   48.000   48.000   48.000  90.00  90.00  90.00 P 1           1                                                                      
HETATM    1  H1  ACE     1      26.508  24.897  24.198  1.00  0.00           H                                                              
HETATM    2  CH3 ACE     1      26.278  25.954  24.075  1.00  0.00           C                                                              
HETATM    3  H2  ACE     1      25.662  26.287  24.908  1.00  0.00           H                                                              
HETATM    4  H3  ACE     1      25.754  26.105  23.134  1.00  0.00           H                                                              
HETATM    5  C   ACE     1      27.577  26.740  24.064  1.00  0.00           C                                                              
HETATM    6  O   ACE     1      28.637  26.123  24.044  1.00  0.00           O                                                              
ATOM      7  N   PRO     2      27.528  28.085  24.064  1.00  0.00           N                                                              
ATOM      8  CD  PRO     2      26.342  28.911  23.908  1.00  0.00           C                                                              
ATOM      9  HD2 PRO     2      25.907  29.109  24.889  1.00  0.00           H                                                              
ATOM     10  HD3 PRO     2      25.604  28.452  23.251  1.00  0.00           H                                                              
ATOM     11  CG  PRO     2      26.851  30.212  23.291  1.00  0.00           C                                                              
ATOM     12  HG2 PRO     2      26.195  31.053  23.521  1.00  0.00           H                                                              
ATOM     13  HG3 PRO     2      26.960  30.089  22.213  1.00  0.00           H                                                              
ATOM     14  CB  PRO     2      28.228  30.369  23.936  1.00  0.00           C                                                              
ATOM     15  HB2 PRO     2      28.110  30.844  24.912  1.00  0.00           H                                                              
ATOM     16  HB3 PRO     2      28.899  30.957  23.310  1.00  0.00           H                                                              
ATOM     17  CA  PRO     2      28.725  28.923  24.097  1.00  0.00           C                                                              
ATOM     18  HA  PRO     2      29.348  28.666  23.240  1.00  0.00           H                                                              
ATOM     19  C   PRO     2      29.531  28.734  25.393  1.00  0.00           C                                                              
ATOM     20  O   PRO     2      29.011  28.944  26.487  1.00  0.00           O                                                              
HETATM   21  N   NME     3      30.814  28.375  25.264  1.00  0.00           N                                                              
HETATM   22  H   NME     3      31.150  28.164  24.338  1.00  0.00           H                                                              
HETATM   23  C   NME     3      31.719  28.162  26.391  1.00  0.00           C                                                              
HETATM   24  H1  NME     3      31.372  27.311  26.981  1.00  0.00           H                                                              
HETATM   25  H2  NME     3      32.730  27.965  26.033  1.00  0.00           H                                                              
HETATM   26  H3  NME     3      31.727  29.047  27.031  1.00  0.00           H                                                              
TER      27      NME     3                                                                                                                                                                                                                                           
END                                                                                                                                                                                                                                                         
"""


cys = """\
REMARK   1 PDBFIXER FROM: MainChain_CYS.pdb                                                                                                 
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-14                                                                                              
CRYST1   48.000   48.000   48.000  90.00  90.00  90.00 P 1           1                                                                      
HETATM    1  H1  ACE     1      25.973  25.174  24.544  1.00  0.00           H                                                              
HETATM    2  CH3 ACE     1      25.997  26.195  24.165  1.00  0.00           C                                                              
HETATM    3  H2  ACE     1      25.393  26.837  24.805  1.00  0.00           H                                                              
HETATM    4  H3  ACE     1      25.611  26.217  23.148  1.00  0.00           H                                                              
HETATM    5  C   ACE     1      27.425  26.686  24.169  1.00  0.00           C                                                              
HETATM    6  O   ACE     1      28.326  25.956  24.559  1.00  0.00           O                                                              
ATOM      7  N   CYS     2      27.631  27.926  23.732  1.00  0.00           N                                                              
ATOM      8  H   CYS     2      26.835  28.482  23.454  1.00  0.00           H                                                              
ATOM      9  CA  CYS     2      28.930  28.601  23.716  1.00  0.00           C                                                              
ATOM     10  HA  CYS     2      29.492  28.302  24.603  1.00  0.00           H                                                              
ATOM     11  CB  CYS     2      29.705  28.139  22.471  1.00  0.00           C                                                              
ATOM     12  HB2 CYS     2      29.719  27.047  22.449  1.00  0.00           H                                                              
ATOM     13  HB3 CYS     2      29.201  28.502  21.573  1.00  0.00           H                                                              
ATOM     14  SG  CYS     2      31.419  28.741  22.503  1.00  0.00           S                                                              
ATOM     15  HG  CYS     2      31.840  28.100  21.406  1.00  0.00           H                                                              
ATOM     16  C   CYS     2      28.734  30.130  23.768  1.00  0.00           C                                                              
ATOM     17  O   CYS     2      27.658  30.623  23.428  1.00  0.00           O                                                              
HETATM   18  N   NME     3      29.759  30.876  24.196  1.00  0.00           N                                                              
HETATM   19  H   NME     3      30.627  30.400  24.395  1.00  0.00           H                                                              
HETATM   20  C   NME     3      29.728  32.333  24.310  1.00  0.00           C                                                              
HETATM   21  H1  NME     3      28.835  32.642  24.858  1.00  0.00           H                                                              
HETATM   22  H2  NME     3      30.614  32.689  24.837  1.00  0.00           H                                                              
HETATM   23  H3  NME     3      29.697  32.778  23.314  1.00  0.00           H                                                              
TER      24      NME     3                                                                                                                                                                                                                                        
END                
"""


cyx = """\
REMARK   1 PDBFIXER FROM: MainChain_CYX.pdb                                                                                                 
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-14                                                                                              
CRYST1   48.000   48.000   48.000  90.00  90.00  90.00 P 1           1                                                                      
HETATM    1  H1  ACE     1      28.573  24.032  26.099  1.00  0.00           H                                                              
HETATM    2  CH3 ACE     1      28.403  25.075  26.364  1.00  0.00           C                                                              
HETATM    3  H2  ACE     1      29.047  25.351  27.196  1.00  0.00           H                                                              
HETATM    4  H3  ACE     1      27.360  25.215  26.643  1.00  0.00           H                                                              
HETATM    5  C   ACE     1      28.726  25.945  25.173  1.00  0.00           C                                                              
HETATM    6  O   ACE     1      29.129  25.442  24.132  1.00  0.00           O                                                              
ATOM      7  N   CYS     2      28.533  27.252  25.319  1.00  0.00           N                                                              
ATOM      8  H   CYS     2      28.195  27.601  26.205  1.00  0.00           H                                                              
ATOM      9  CA  CYS     2      28.788  28.255  24.284  1.00  0.00           C                                                              
ATOM     10  HA  CYS     2      29.726  28.000  23.788  1.00  0.00           H                                                              
ATOM     11  CB  CYS     2      27.665  28.212  23.233  1.00  0.00           C                                                              
ATOM     12  HB2 CYS     2      27.542  27.183  22.893  1.00  0.00           H                                                              
ATOM     13  HB3 CYS     2      27.997  28.784  22.367  1.00  0.00           H                                                              
ATOM     14  SG  CYS     2      26.035  28.864  23.703  1.00  0.00           S                                                              
ATOM     15  C   CYS     2      28.948  29.658  24.899  1.00  0.00           C                                                              
ATOM     16  O   CYS     2      28.616  29.871  26.065  1.00  0.00           O                                                              
HETATM   17  N   NME     3      29.464  30.614  24.115  1.00  0.00           N                                                              
HETATM   18  H   NME     3      29.717  30.370  23.171  1.00  0.00           H                                                              
HETATM   19  C   NME     3      29.676  31.993  24.550  1.00  0.00           C                                                              
HETATM   20  H1  NME     3      28.717  32.449  24.805  1.00  0.00           H                                                              
HETATM   21  H2  NME     3      30.312  32.007  25.437  1.00  0.00           H                                                              
HETATM   22  H3  NME     3      30.151  32.572  23.757  1.00  0.00           H                                                              
TER      23      NME     3                                                                                                                  
HETATM   24  H1  ACE     4      23.123  21.295  23.934  1.00  0.00           H                                                              
HETATM   25  CH3 ACE     4      23.375  22.084  23.227  1.00  0.00           C                                                              
HETATM   26  H2  ACE     4      22.641  22.105  22.423  1.00  0.00           H                                                              
HETATM   27  H3  ACE     4      24.367  21.904  22.817  1.00  0.00           H                                                              
HETATM   28  C   ACE     4      23.365  23.411  23.949  1.00  0.00           C                                                              
HETATM   29  O   ACE     4      23.112  23.462  25.145  1.00  0.00           O                                                              
ATOM     30  N   CYS     5      23.650  24.486  23.219  1.00  0.00           N                                                              
ATOM     31  H   CYS     5      23.839  24.370  22.234  1.00  0.00           H                                                              
ATOM     32  CA  CYS     5      23.696  25.853  23.738  1.00  0.00           C                                                              
ATOM     33  HA  CYS     5      22.909  25.977  24.484  1.00  0.00           H                                                              
ATOM     34  CB  CYS     5      25.052  26.054  24.429  1.00  0.00           C                                                              
ATOM     35  HB2 CYS     5      25.136  25.300  25.213  1.00  0.00           H                                                              
ATOM     36  HB3 CYS     5      25.847  25.859  23.709  1.00  0.00           H                                                              
ATOM     37  SG  CYS     5      25.362  27.661  25.209  1.00  0.00           S                                                              
ATOM     38  C   CYS     5      23.449  26.867  22.607  1.00  0.00           C                                                              
ATOM     39  O   CYS     5      23.798  26.601  21.456  1.00  0.00           O                                                              
HETATM   40  N   NME     6      22.844  28.017  22.931  1.00  0.00           N                                                              
HETATM   41  H   NME     6      22.620  28.173  23.901  1.00  0.00           H                                                              
HETATM   42  C   NME     6      22.545  29.086  21.979  1.00  0.00           C                                                              
HETATM   43  H1  NME     6      23.464  29.613  21.717  1.00  0.00           H                                                              
HETATM   44  H2  NME     6      22.113  28.665  21.068  1.00  0.00           H                                                              
HETATM   45  H3  NME     6      21.838  29.792  22.415  1.00  0.00           H                                                              
TER      46      NME     6  
END
"""


lys = """\
REMARK   1 PDBFIXER FROM: lys.pdb                                                                                                           
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-15                                                                                              
HETATM    1  H1  ACE     1       3.851  -3.562   0.131  1.00  0.00           H                                                              
HETATM    2  CH3 ACE     1       3.389  -3.037  -0.713  1.00  0.00           C                                                              
HETATM    3  H2  ACE     1       3.132  -3.704  -1.541  1.00  0.00           H                                                              
HETATM    4  H3  ACE     1       4.006  -2.170  -1.057  1.00  0.00           H                                                              
HETATM    5  C   ACE     1       2.115  -2.449  -0.125  1.00  0.00           C                                                              
HETATM    6  O   ACE     1       1.100  -3.195  -0.165  1.00  0.00           O                                                              
ATOM      7  N   LYS     2       2.145  -1.000   0.467  1.00  0.00           N                                                              
ATOM      8  CA  LYS     2       0.684  -0.172   0.107  1.00  0.00           C                                                              
ATOM      9  C   LYS     2       1.161   1.161  -0.224  1.00  0.00           C                                                              
ATOM     10  O   LYS     2       0.958   1.659  -1.347  1.00  0.00           O                                                              
ATOM     11  CB  LYS     2      -0.253  -0.408   1.195  1.00  0.00           C                                                              
ATOM     12  CG  LYS     2      -1.599   0.233   1.121  1.00  0.00           C                                                              
ATOM     13  CD  LYS     2      -2.466  -0.159  -0.054  1.00  0.00           C                                                              
ATOM     14  CE  LYS     2      -3.754   0.617   0.150  1.00  0.00           C                                                              
ATOM     15  NZ  LYS     2      -4.785   0.362  -0.787  1.00  0.00           N                                                              
ATOM     16  H   LYS     2       2.562  -0.670   1.540  1.00  0.00           H                                                              
ATOM     17  HA  LYS     2       0.300  -0.689  -0.853  1.00  0.00           H                                                              
ATOM     18  HB2 LYS     2      -0.509  -1.551   1.129  1.00  0.00           H                                                              
ATOM     19  HB3 LYS     2       0.163  -0.267   2.211  1.00  0.00           H                                                              
ATOM     20  HG2 LYS     2      -1.613   1.336   1.159  1.00  0.00           H                                                              
ATOM     21  HG3 LYS     2      -2.164  -0.139   2.024  1.00  0.00           H                                                              
ATOM     22  HD2 LYS     2      -1.955   0.239  -0.949  1.00  0.00           H                                                              
ATOM     23  HD3 LYS     2      -2.615  -1.226  -0.207  1.00  0.00           H                                                              
ATOM     24  HE2 LYS     2      -4.119   0.453   1.180  1.00  0.00           H                                                              
ATOM     25  HE3 LYS     2      -3.466   1.712   0.126  1.00  0.00           H                                                              
ATOM     26  HZ1 LYS     2      -4.575   0.480  -1.785  1.00  0.00           H                                                              
ATOM     27  HZ2 LYS     2      -5.564   1.061  -0.557  1.00  0.00           H                                                              
ATOM     28  HZ3 LYS     2      -5.249  -0.563  -0.540  1.00  0.00           H                                                              
HETATM   29  N   NME     3       2.024   2.091   0.877  1.00  0.00           N                                                              
HETATM   30  H   NME     3       1.295   3.264   1.380  1.00  0.00           H                                                              
HETATM   31  C   NME     3       3.556   2.525   0.096  1.00  0.00           C                                                              
HETATM   32  H1  NME     3       4.519   3.117   1.252  1.00  0.00           H                                                              
HETATM   33  H2  NME     3       4.264   1.218  -0.292  1.00  0.00           H                                                              
HETATM   34  H3  NME     3       3.463   3.431  -1.036  1.00  0.00           H                                                              
TER      35      NME     3                                                                                                                                                                                                                                        
END                                                                                                                                         
"""


asp = """\
REMARK   1 PDBFIXER FROM: asp.pdb
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-16
HETATM    1  H1  ACE     1      -2.705  -3.217   0.456  1.00  0.00           H  
HETATM    2  CH3 ACE     1      -3.024  -2.157   0.693  1.00  0.00           C  
HETATM    3  H2  ACE     1      -3.650  -2.037   1.542  1.00  0.00           H  
HETATM    4  H3  ACE     1      -3.578  -1.864  -0.271  1.00  0.00           H  
HETATM    5  C   ACE     1      -1.784  -1.307   0.763  1.00  0.00           C  
HETATM    6  O   ACE     1      -1.321  -1.088   1.857  1.00  0.00           O  
ATOM      7  N   ASP     2      -1.103  -0.732  -0.532  1.00  0.00           N  
ATOM      8  CA  ASP     2      -0.051   0.501  -0.003  1.00  0.00           C  
ATOM      9  C   ASP     2       1.187   0.212  -0.739  1.00  0.00           C  
ATOM     10  O   ASP     2       1.613   0.828  -1.748  1.00  0.00           O  
ATOM     11  CB  ASP     2      -0.799   1.698  -0.478  1.00  0.00           C  
ATOM     12  CG  ASP     2      -0.320   3.035  -0.340  1.00  0.00           C  
ATOM     13  OD1 ASP     2       0.193   3.550  -1.430  1.00  0.00           O  
ATOM     14  OD2 ASP     2      -0.296   3.907   0.710  1.00  0.00           O  
ATOM     15  H   ASP     2      -0.664  -1.422  -1.301  1.00  0.00           H  
ATOM     16  HA  ASP     2       0.091   0.492   1.050  1.00  0.00           H  
ATOM     17  HB2 ASP     2      -1.138   1.514  -1.584  1.00  0.00           H  
ATOM     18  HB3 ASP     2      -1.827   1.550   0.014  1.00  0.00           H  
ATOM     19  HD2 ASP     2      -0.651   4.887   0.567  1.00  0.00           H  
HETATM   20  N   NME     3       2.060  -1.094  -0.134  1.00  0.00           N  
HETATM   21  H   NME     3       1.505  -2.428  -0.518  1.00  0.00           H  
HETATM   22  C   NME     3       3.545  -1.013   0.401  1.00  0.00           C  
HETATM   23  H1  NME     3       3.903  -2.282   1.333  1.00  0.00           H  
HETATM   24  H2  NME     3       4.537  -1.650  -1.016  1.00  0.00           H  
HETATM   25  H3  NME     3       4.277   0.116   0.704  1.00  0.00           H  
TER      26      NME     3
CONECT    1    2
CONECT    2    5    1    3    4
CONECT    3    2
CONECT    4    2
CONECT    5    2    6    7
CONECT    6    5
CONECT    7    5
CONECT    9   20
CONECT   20    9   22   21
CONECT   21   20
CONECT   22   23   24   25   20
CONECT   23   22
CONECT   24   22
CONECT   25   22
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


def test_trp():
    m = Chem.MolFromPDBBlock(trp, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 36
    assert m.GetNumBonds() == 37
    for at in m.GetAtoms():
        assert at.GetFormalCharge() == 0


def test_pro():
    m = Chem.MolFromPDBBlock(pro, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 26
    assert m.GetNumBonds() == 26
    for at in m.GetAtoms():
        assert at.GetFormalCharge() == 0


def test_cys():
    m = Chem.MolFromPDBBlock(cys, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 23
    assert m.GetNumBonds() == 22
    for at in m.GetAtoms():
        assert at.GetFormalCharge() == 0


def test_cyx():
    # has sulphur-sulphur bond
    m = Chem.MolFromPDBBlock(cyx, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 44
    assert m.GetNumBonds() == 43
    # S-S bond, will return None if this bond does not exist
    assert m.GetBondBetweenAtoms(13, 35) is not None


def test_lys():
    m = Chem.MolFromPDBBlock(lys, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 34
    assert m.GetNumBonds() == 33
    assert m.GetAtomWithIdx(14).GetFormalCharge() == 1


def test_asp():
    m = Chem.MolFromPDBBlock(asp, proximityBonding=False, removeHs=False)
    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 25
    assert m.GetNumBonds() == 24

