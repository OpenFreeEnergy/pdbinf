# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
import pdbinf
import pytest

DA_pdb = """\
ATOM      1  OP3  DA A   1       1.845  -1.282  -5.339  1.00 10.00           O
ATOM      2  P    DA A   1       0.934  -0.156  -4.636  1.00 10.00           P
ATOM      3  OP1  DA A   1       1.781   0.996  -4.255  1.00 10.00           O
ATOM      4  OP2  DA A   1      -0.204   0.331  -5.665  1.00 10.00           O
ATOM      5  O5'  DA A   1       0.241  -0.771  -3.320  1.00 10.00           O
ATOM      6  C5'  DA A   1      -0.549   0.270  -2.744  1.00 10.00           C
ATOM      7  C4'  DA A   1      -1.239  -0.251  -1.482  1.00 10.00           C
ATOM      8  O4'  DA A   1      -0.267  -0.564  -0.458  1.00 10.00           O
ATOM      9  C3'  DA A   1      -2.105   0.859  -0.835  1.00 10.00           C
ATOM     10  O3'  DA A   1      -3.409   0.895  -1.418  1.00 10.00           O
ATOM     11  C2'  DA A   1      -2.173   0.398   0.640  1.00 10.00           C
ATOM     12  C1'  DA A   1      -0.965  -0.545   0.797  1.00 10.00           C
ATOM     13  N9   DA A   1      -0.078  -0.047   1.852  1.00 10.00           N
ATOM     14  C8   DA A   1       0.962   0.817   1.689  1.00 10.00           C
ATOM     15  N7   DA A   1       1.535   1.044   2.835  1.00 10.00           N
ATOM     16  C5   DA A   1       0.897   0.346   3.805  1.00 10.00           C
ATOM     17  C6   DA A   1       1.069   0.196   5.191  1.00 10.00           C
ATOM     18  N6   DA A   1       2.079   0.869   5.856  1.00 10.00           N
ATOM     19  N1   DA A   1       0.236  -0.603   5.850  1.00 10.00           N
ATOM     20  C2   DA A   1      -0.729  -1.249   5.224  1.00 10.00           C
ATOM     21  N3   DA A   1      -0.925  -1.144   3.927  1.00 10.00           N
ATOM     22  C4   DA A   1      -0.142  -0.368   3.184  1.00 10.00           C
ATOM     23 HOP3  DA A   1       2.241  -0.873  -6.121  1.00 10.00           H
ATOM     24 HOP2  DA A   1      -0.732  -0.447  -5.887  1.00 10.00           H
ATOM     25  H5'  DA A   1      -1.302   0.594  -3.463  1.00 10.00           H
ATOM     26 H5''  DA A   1       0.092   1.112  -2.486  1.00 10.00           H
ATOM     27  H4'  DA A   1      -1.846  -1.126  -1.712  1.00 10.00           H
ATOM     28  H3'  DA A   1      -1.617   1.830  -0.918  1.00 10.00           H
ATOM     29 HO3'  DA A   1      -3.924   1.538  -0.913  1.00 10.00           H
ATOM     30  H2'  DA A   1      -3.103  -0.136   0.831  1.00 10.00           H
ATOM     31 H2''  DA A   1      -2.079   1.253   1.311  1.00 10.00           H
ATOM     32  H1'  DA A   1      -1.309  -1.549   1.046  1.00 10.00           H
ATOM     33  H8   DA A   1       1.266   1.250   0.748  1.00 10.00           H
ATOM     34  H61  DA A   1       2.185   0.761   6.814  1.00 10.00           H
ATOM     35  H62  DA A   1       2.683   1.447   5.363  1.00 10.00           H
ATOM     36  H2   DA A   1      -1.383  -1.889   5.798  1.00 10.00           H
END
"""

DC_pdb = """\
ATOM      1  OP3  DC A   1       1.941  -1.055  -4.672  1.00 10.00           O
ATOM      2  P    DC A   1       0.987  -0.017  -3.894  1.00 10.00           P
ATOM      3  OP1  DC A   1       1.802   1.099  -3.365  1.00 10.00           O
ATOM      4  OP2  DC A   1      -0.119   0.560  -4.910  1.00 10.00           O
ATOM      5  O5'  DC A   1       0.255  -0.772  -2.674  1.00 10.00           O
ATOM      6  C5'  DC A   1      -0.571   0.196  -2.027  1.00 10.00           C
ATOM      7  C4'  DC A   1      -1.300  -0.459  -0.852  1.00 10.00           C
ATOM      8  O4'  DC A   1      -0.363  -0.863   0.171  1.00 10.00           O
ATOM      9  C3'  DC A   1      -2.206   0.569  -0.129  1.00 10.00           C
ATOM     10  O3'  DC A   1      -3.488   0.649  -0.756  1.00 10.00           O
ATOM     11  C2'  DC A   1      -2.322  -0.040   1.288  1.00 10.00           C
ATOM     12  C1'  DC A   1      -1.106  -0.981   1.395  1.00 10.00           C
ATOM     13  N1   DC A   1      -0.267  -0.584   2.528  1.00 10.00           N
ATOM     14  C2   DC A   1       0.270   0.648   2.563  1.00 10.00           C
ATOM     15  O2   DC A   1       0.052   1.424   1.647  1.00 10.00           O
ATOM     16  N3   DC A   1       1.037   1.035   3.581  1.00 10.00           N
ATOM     17  C4   DC A   1       1.291   0.212   4.589  1.00 10.00           C
ATOM     18  N4   DC A   1       2.085   0.622   5.635  1.00 10.00           N
ATOM     19  C5   DC A   1       0.746  -1.088   4.580  1.00 10.00           C
ATOM     20  C6   DC A   1      -0.035  -1.465   3.541  1.00 10.00           C
ATOM     21 HOP3  DC A   1       2.359  -0.564  -5.392  1.00 10.00           H
ATOM     22 HOP2  DC A   1      -0.626  -0.198  -5.231  1.00 10.00           H
ATOM     23  H5'  DC A   1      -1.302   0.583  -2.737  1.00 10.00           H
ATOM     24 H5''  DC A   1       0.046   1.015  -1.659  1.00 10.00           H
ATOM     25  H4'  DC A   1      -1.885  -1.313  -1.193  1.00 10.00           H
ATOM     26  H3'  DC A   1      -1.731   1.549  -0.094  1.00 10.00           H
ATOM     27 HO3'  DC A   1      -4.031   1.232  -0.207  1.00 10.00           H
ATOM     28  H2'  DC A   1      -3.250  -0.602   1.387  1.00 10.00           H
ATOM     29 H2''  DC A   1      -2.266   0.742   2.046  1.00 10.00           H
ATOM     30  H1'  DC A   1      -1.444  -2.009   1.526  1.00 10.00           H
ATOM     31  H41  DC A   1       2.461   1.516   5.636  1.00 10.00           H
ATOM     32  H42  DC A   1       2.265   0.019   6.373  1.00 10.00           H
ATOM     33  H5   DC A   1       0.943  -1.771   5.394  1.00 10.00           H
ATOM     34  H6   DC A   1      -0.467  -2.454   3.514  1.00 10.00           H
END
"""

DG_pdb = """\
ATOM      1  OP3  DG A   1      -1.603  -1.547   5.624  1.00 10.00           O
ATOM      2  P    DG A   1      -0.818  -0.321   4.935  1.00 10.00           P
ATOM      3  OP1  DG A   1      -1.774   0.766   4.630  1.00 10.00           O
ATOM      4  OP2  DG A   1       0.312   0.224   5.941  1.00 10.00           O
ATOM      5  O5'  DG A   1      -0.126  -0.826   3.572  1.00 10.00           O
ATOM      6  C5'  DG A   1       0.550   0.300   3.011  1.00 10.00           C
ATOM      7  C4'  DG A   1       1.233  -0.113   1.706  1.00 10.00           C
ATOM      8  O4'  DG A   1       0.253  -0.471   0.705  1.00 10.00           O
ATOM      9  C3'  DG A   1       1.976   1.091   1.073  1.00 10.00           C
ATOM     10  O3'  DG A   1       3.294   1.218   1.612  1.00 10.00           O
ATOM     11  C2'  DG A   1       2.026   0.692  -0.421  1.00 10.00           C
ATOM     12  C1'  DG A   1       0.897  -0.345  -0.573  1.00 10.00           C
ATOM     13  N9   DG A   1      -0.068   0.111  -1.575  1.00 10.00           N
ATOM     14  C8   DG A   1      -1.172   0.877  -1.341  1.00 10.00           C
ATOM     15  N7   DG A   1      -1.804   1.094  -2.458  1.00 10.00           N
ATOM     16  C5   DG A   1      -1.145   0.482  -3.472  1.00 10.00           C
ATOM     17  C6   DG A   1      -1.361   0.377  -4.866  1.00 10.00           C
ATOM     18  O6   DG A   1      -2.321   0.914  -5.391  1.00 10.00           O
ATOM     19  N1   DG A   1      -0.473  -0.327  -5.601  1.00 10.00           N
ATOM     20  C2   DG A   1       0.593  -0.928  -5.003  1.00 10.00           C
ATOM     21  N2   DG A   1       1.474  -1.643  -5.774  1.00 10.00           N
ATOM     22  N3   DG A   1       0.804  -0.839  -3.709  1.00 10.00           N
ATOM     23  C4   DG A   1      -0.027  -0.152  -2.917  1.00 10.00           C
ATOM     24 HOP3  DG A   1      -2.002  -1.203   6.435  1.00 10.00           H
ATOM     25 HOP2  DG A   1       0.913  -0.513   6.114  1.00 10.00           H
ATOM     26  H5'  DG A   1       1.299   0.661   3.715  1.00 10.00           H
ATOM     27 H5''  DG A   1      -0.170   1.093   2.808  1.00 10.00           H
ATOM     28  H4'  DG A   1       1.921  -0.940   1.879  1.00 10.00           H
ATOM     29  H3'  DG A   1       1.411   2.013   1.211  1.00 10.00           H
ATOM     30 HO3'  DG A   1       3.732   1.921   1.114  1.00 10.00           H
ATOM     31  H2'  DG A   1       2.990   0.246  -0.665  1.00 10.00           H
ATOM     32 H2''  DG A   1       1.834   1.559  -1.053  1.00 10.00           H
ATOM     33  H1'  DG A   1       1.316  -1.306  -0.873  1.00 10.00           H
ATOM     34  H8   DG A   1      -1.477   1.248  -0.373  1.00 10.00           H
ATOM     35  H1   DG A   1      -0.601  -0.413  -6.559  1.00 10.00           H
ATOM     36  H21  DG A   1       2.240  -2.073  -5.363  1.00 10.00           H
ATOM     37  H22  DG A   1       1.329  -1.722  -6.730  1.00 10.00           H
END
"""

DT_pdb = """\
ATOM      1  OP3  DT A   1      -3.912  -2.311   1.636  1.00 10.00           O
ATOM      2  P    DT A   1      -3.968  -1.665   3.118  1.00 10.00           P
ATOM      3  OP1  DT A   1      -4.406  -2.599   4.208  1.00 10.00           O
ATOM      4  OP2  DT A   1      -4.901  -0.360   2.920  1.00 10.00           O
ATOM      5  O5'  DT A   1      -2.493  -1.028   3.315  1.00 10.00           O
ATOM      6  C5'  DT A   1      -2.005  -0.136   2.327  1.00 10.00           C
ATOM      7  C4'  DT A   1      -0.611   0.328   2.728  1.00 10.00           C
ATOM      8  O4'  DT A   1       0.247  -0.829   2.764  1.00 10.00           O
ATOM      9  C3'  DT A   1       0.008   1.286   1.720  1.00 10.00           C
ATOM     10  O3'  DT A   1       0.965   2.121   2.368  1.00 10.00           O
ATOM     11  C2'  DT A   1       0.710   0.360   0.754  1.00 10.00           C
ATOM     12  C1'  DT A   1       1.157  -0.778   1.657  1.00 10.00           C
ATOM     13  N1   DT A   1       1.164  -2.047   0.989  1.00 10.00           N
ATOM     14  C2   DT A   1       2.333  -2.544   0.374  1.00 10.00           C
ATOM     15  O2   DT A   1       3.410  -1.945   0.363  1.00 10.00           O
ATOM     16  N3   DT A   1       2.194  -3.793  -0.240  1.00 10.00           N
ATOM     17  C4   DT A   1       1.047  -4.570  -0.300  1.00 10.00           C
ATOM     18  O4   DT A   1       0.995  -5.663  -0.857  1.00 10.00           O
ATOM     19  C5   DT A   1      -0.143  -3.980   0.369  1.00 10.00           C
ATOM     20  C7   DT A   1      -1.420  -4.757   0.347  1.00 10.00           C
ATOM     21  C6   DT A   1      -0.013  -2.784   0.958  1.00 10.00           C
ATOM     22 HOP3  DT A   1      -4.684  -2.823   1.313  1.00 10.00           H
ATOM     23 HOP2  DT A   1      -5.874  -0.475   2.871  1.00 10.00           H
ATOM     24  H5'  DT A   1      -1.971  -0.654   1.366  1.00 10.00           H
ATOM     25 H5''  DT A   1      -2.683   0.718   2.254  1.00 10.00           H
ATOM     26  H4'  DT A   1      -0.644   0.753   3.736  1.00 10.00           H
ATOM     27  H3'  DT A   1      -0.722   1.941   1.236  1.00 10.00           H
ATOM     28 HO3'  DT A   1       0.743   3.042   2.146  1.00 10.00           H
ATOM     29  H2'  DT A   1      -0.003  -0.001   0.003  1.00 10.00           H
ATOM     30 H2''  DT A   1       1.556   0.818   0.232  1.00 10.00           H
ATOM     31  H1'  DT A   1       2.163  -0.604   2.054  1.00 10.00           H
ATOM     32  H3   DT A   1       3.023  -4.171  -0.689  1.00 10.00           H
ATOM     33  H71  DT A   1      -1.735  -4.916  -0.683  1.00 10.00           H
ATOM     34  H72  DT A   1      -2.195  -4.208   0.881  1.00 10.00           H
ATOM     35  H73  DT A   1      -1.272  -5.721   0.833  1.00 10.00           H
ATOM     36  H6   DT A   1      -0.857  -2.318   1.458  1.00 10.00           H
END
"""


@pytest.fixture
def DA(tmpdir):
    fn = tmpdir / 'DA.pdb'
    with open(fn, 'w') as out:
        out.write(DA_pdb)

    yield str(fn)


@pytest.fixture
def DC(tmpdir):
    fn = tmpdir / 'DC.pdb'
    with open(fn, 'w') as out:
        out.write(DC_pdb)

    yield str(fn)


@pytest.fixture
def DG(tmpdir):
    fn = tmpdir / 'DG.pdb'
    with open(fn, 'w') as out:
        out.write(DG_pdb)

    yield str(fn)


@pytest.fixture
def DT(tmpdir):
    fn = tmpdir / 'DT.pdb'
    with open(fn, 'w') as out:
        out.write(DT_pdb)

    yield str(fn)


def test_DA(DA):
    m = pdbinf.load_pdb_file(DA, templates=[pdbinf.DNA_DOC])

    assert m
    assert m.GetNumAtoms() == 36
    assert m.GetNumBonds() == 38


def test_DC(DC):
    m = pdbinf.load_pdb_file(DC, templates=[pdbinf.DNA_DOC])

    assert m
    assert m.GetNumAtoms() == 34
    assert m.GetNumBonds() == 35


def test_DG(DG):
    m = pdbinf.load_pdb_file(DG, templates=[pdbinf.DNA_DOC])

    assert m
    assert m.GetNumAtoms() == 37
    assert m.GetNumBonds() == 39


def test_DT(DT):
    m = pdbinf.load_pdb_file(DT, templates=[pdbinf.DNA_DOC])

    assert m
    assert m.GetNumAtoms() == 36
    assert m.GetNumBonds() == 37
