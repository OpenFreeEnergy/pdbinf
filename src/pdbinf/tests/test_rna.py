# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
import pdbinf
import pytest


A_pdb = """\
ATOM      1  OP3   A A   1       2.135  -1.141  -5.313  1.00 10.00           O
ATOM      2  P     A A   1       1.024  -0.137  -4.723  1.00 10.00           P
ATOM      3  OP1   A A   1       1.633   1.190  -4.488  1.00 10.00           O
ATOM      4  OP2   A A   1      -0.183   0.005  -5.778  1.00 10.00           O
ATOM      5  O5'   A A   1       0.456  -0.720  -3.334  1.00 10.00           O
ATOM      6  C5'   A A   1      -0.520   0.209  -2.863  1.00 10.00           C
ATOM      7  C4'   A A   1      -1.101  -0.287  -1.538  1.00 10.00           C
ATOM      8  O4'   A A   1      -0.064  -0.383  -0.538  1.00 10.00           O
ATOM      9  C3'   A A   1      -2.105   0.739  -0.969  1.00 10.00           C
ATOM     10  O3'   A A   1      -3.445   0.360  -1.287  1.00 10.00           O
ATOM     11  C2'   A A   1      -1.874   0.684   0.558  1.00 10.00           C
ATOM     12  O2'   A A   1      -3.065   0.271   1.231  1.00 10.00           O
ATOM     13  C1'   A A   1      -0.755  -0.367   0.729  1.00 10.00           C
ATOM     14  N9    A A   1       0.158   0.029   1.803  1.00 10.00           N
ATOM     15  C8    A A   1       1.265   0.813   1.672  1.00 10.00           C
ATOM     16  N7    A A   1       1.843   0.963   2.828  1.00 10.00           N
ATOM     17  C5    A A   1       1.143   0.292   3.773  1.00 10.00           C
ATOM     18  C6    A A   1       1.290   0.091   5.156  1.00 10.00           C
ATOM     19  N6    A A   1       2.344   0.664   5.846  1.00 10.00           N
ATOM     20  N1    A A   1       0.391  -0.656   5.787  1.00 10.00           N
ATOM     21  C2    A A   1      -0.617  -1.206   5.136  1.00 10.00           C
ATOM     22  N3    A A   1      -0.792  -1.051   3.841  1.00 10.00           N
ATOM     23  C4    A A   1       0.056  -0.320   3.126  1.00 10.00           C
ATOM     24 HOP3   A A   1       2.448  -0.755  -6.142  1.00 10.00           H
ATOM     25 HOP2   A A   1      -0.552  -0.879  -5.902  1.00 10.00           H
ATOM     26  H5'   A A   1      -1.319   0.301  -3.599  1.00 10.00           H
ATOM     27 H5''   A A   1      -0.052   1.182  -2.712  1.00 10.00           H
ATOM     28  H4'   A A   1      -1.586  -1.254  -1.677  1.00 10.00           H
ATOM     29  H3'   A A   1      -1.890   1.736  -1.353  1.00 10.00           H
ATOM     30 HO3'   A A   1      -4.024   1.035  -0.908  1.00 10.00           H
ATOM     31  H2'   A A   1      -1.543   1.654   0.930  1.00 10.00           H
ATOM     32 HO2'   A A   1      -3.740   0.936   1.037  1.00 10.00           H
ATOM     33  H1'   A A   1      -1.185  -1.346   0.940  1.00 10.00           H
ATOM     34  H8    A A   1       1.611   1.246   0.745  1.00 10.00           H
ATOM     35  H61   A A   1       2.432   0.522   6.801  1.00 10.00           H
ATOM     36  H62   A A   1       2.996   1.205   5.374  1.00 10.00           H
ATOM     37  H2    A A   1      -1.325  -1.807   5.688  1.00 10.00           H
END   
"""

C_pdb = """\
ATOM      1  OP3   C A   1       2.147  -1.021  -4.678  1.00 10.00           O
ATOM      2  P     C A   1       1.049  -0.039  -4.028  1.00 10.00           P
ATOM      3  OP1   C A   1       1.692   1.237  -3.646  1.00 10.00           O
ATOM      4  OP2   C A   1      -0.116   0.246  -5.102  1.00 10.00           O
ATOM      5  O5'   C A   1       0.415  -0.733  -2.721  1.00 10.00           O
ATOM      6  C5'   C A   1      -0.546   0.181  -2.193  1.00 10.00           C
ATOM      7  C4'   C A   1      -1.189  -0.419  -0.942  1.00 10.00           C
ATOM      8  O4'   C A   1      -0.190  -0.648   0.076  1.00 10.00           O
ATOM      9  C3'   C A   1      -2.178   0.583  -0.307  1.00 10.00           C
ATOM     10  O3'   C A   1      -3.518   0.283  -0.703  1.00 10.00           O
ATOM     11  C2'   C A   1      -2.001   0.373   1.215  1.00 10.00           C
ATOM     12  O2'   C A   1      -3.228  -0.059   1.806  1.00 10.00           O
ATOM     13  C1'   C A   1      -0.924  -0.729   1.317  1.00 10.00           C
ATOM     14  N1    C A   1      -0.036  -0.470   2.453  1.00 10.00           N
ATOM     15  C2    C A   1       0.652   0.683   2.514  1.00 10.00           C
ATOM     16  O2    C A   1       0.529   1.504   1.620  1.00 10.00           O
ATOM     17  N3    C A   1       1.467   0.945   3.535  1.00 10.00           N
ATOM     18  C4    C A   1       1.620   0.070   4.520  1.00 10.00           C
ATOM     19  N4    C A   1       2.464   0.350   5.569  1.00 10.00           N
ATOM     20  C5    C A   1       0.916  -1.151   4.483  1.00 10.00           C
ATOM     21  C6    C A   1       0.087  -1.399   3.442  1.00 10.00           C
ATOM     22 HOP3   C A   1       2.501  -0.569  -5.456  1.00 10.00           H
ATOM     23 HOP2   C A   1      -0.508  -0.608  -5.323  1.00 10.00           H
ATOM     24  H5'   C A   1      -1.315   0.371  -2.941  1.00 10.00           H
ATOM     25 H5''   C A   1      -0.052   1.118  -1.933  1.00 10.00           H
ATOM     26  H4'   C A   1      -1.699  -1.350  -1.188  1.00 10.00           H
ATOM     27  H3'   C A   1      -1.917   1.604  -0.586  1.00 10.00           H
ATOM     28 HO3'   C A   1      -4.088   0.939  -0.278  1.00 10.00           H
ATOM     29  H2'   C A   1      -1.653   1.290   1.689  1.00 10.00           H
ATOM     30 HO2'   C A   1      -3.874   0.644   1.656  1.00 10.00           H
ATOM     31  H1'   C A   1      -1.392  -1.708   1.418  1.00 10.00           H
ATOM     32  H41   C A   1       2.950   1.189   5.590  1.00 10.00           H
ATOM     33  H42   C A   1       2.571  -0.289   6.290  1.00 10.00           H
ATOM     34  H5    C A   1       1.030  -1.873   5.278  1.00 10.00           H
ATOM     35  H6    C A   1      -0.465  -2.326   3.393  1.00 10.00           H
END   
"""

G_pdb = """\
ATOM      1  OP3   G A   1      -1.945  -1.360   5.599  1.00 10.00           O
ATOM      2  P     G A   1      -0.911  -0.277   5.008  1.00 10.00           P
ATOM      3  OP1   G A   1      -1.598   1.022   4.844  1.00 10.00           O
ATOM      4  OP2   G A   1       0.325  -0.105   6.025  1.00 10.00           O
ATOM      5  O5'   G A   1      -0.365  -0.780   3.580  1.00 10.00           O
ATOM      6  C5'   G A   1       0.542   0.217   3.109  1.00 10.00           C
ATOM      7  C4'   G A   1       1.100  -0.200   1.748  1.00 10.00           C
ATOM      8  O4'   G A   1       0.033  -0.318   0.782  1.00 10.00           O
ATOM      9  C3'   G A   1       2.025   0.898   1.182  1.00 10.00           C
ATOM     10  O3'   G A   1       3.395   0.582   1.439  1.00 10.00           O
ATOM     11  C2'   G A   1       1.741   0.884  -0.338  1.00 10.00           C
ATOM     12  O2'   G A   1       2.927   0.560  -1.066  1.00 10.00           O
ATOM     13  C1'   G A   1       0.675  -0.220  -0.507  1.00 10.00           C
ATOM     14  N9    G A   1      -0.297   0.162  -1.534  1.00 10.00           N
ATOM     15  C8    G A   1      -1.440   0.880  -1.334  1.00 10.00           C
ATOM     16  N7    G A   1      -2.066   1.037  -2.464  1.00 10.00           N
ATOM     17  C5    G A   1      -1.364   0.431  -3.453  1.00 10.00           C
ATOM     18  C6    G A   1      -1.556   0.279  -4.846  1.00 10.00           C
ATOM     19  O6    G A   1      -2.534   0.755  -5.397  1.00 10.00           O
ATOM     20  N1    G A   1      -0.626  -0.401  -5.551  1.00 10.00           N
ATOM     21  C2    G A   1       0.459  -0.934  -4.923  1.00 10.00           C
ATOM     22  N2    G A   1       1.384  -1.626  -5.664  1.00 10.00           N
ATOM     23  N3    G A   1       0.649  -0.800  -3.630  1.00 10.00           N
ATOM     24  C4    G A   1      -0.226  -0.134  -2.868  1.00 10.00           C
ATOM     25 HOP3   G A   1      -2.247  -1.021   6.453  1.00 10.00           H
ATOM     26 HOP2   G A   1       0.745  -0.973   6.104  1.00 10.00           H
ATOM     27  H5'   G A   1       1.362   0.327   3.820  1.00 10.00           H
ATOM     28 H5''   G A   1       0.018   1.168   3.011  1.00 10.00           H
ATOM     29  H4'   G A   1       1.640  -1.144   1.833  1.00 10.00           H
ATOM     30  H3'   G A   1       1.772   1.868   1.610  1.00 10.00           H
ATOM     31 HO3'   G A   1       3.923   1.300   1.065  1.00 10.00           H
ATOM     32  H2'   G A   1       1.346   1.847  -0.662  1.00 10.00           H
ATOM     33 HO2'   G A   1       3.573   1.254  -0.871  1.00 10.00           H
ATOM     34  H1'   G A   1       1.148  -1.167  -0.769  1.00 10.00           H
ATOM     35  H8    G A   1      -1.776   1.261  -0.381  1.00 10.00           H
ATOM     36  H1    G A   1      -0.736  -0.518  -6.508  1.00 10.00           H
ATOM     37  H21   G A   1       2.165  -2.007  -5.232  1.00 10.00           H
ATOM     38  H22   G A   1       1.256  -1.736  -6.619  1.00 10.00           H
END   
"""

U_pdb = """\
ATOM      1  OP3   U A   1      -2.122   1.033  -4.690  1.00 10.00           O
ATOM      2  P     U A   1      -1.030   0.047  -4.037  1.00 10.00           P
ATOM      3  OP1   U A   1      -1.679  -1.228  -3.660  1.00 10.00           O
ATOM      4  OP2   U A   1       0.138  -0.241  -5.107  1.00 10.00           O
ATOM      5  O5'   U A   1      -0.399   0.736  -2.726  1.00 10.00           O
ATOM      6  C5'   U A   1       0.557  -0.182  -2.196  1.00 10.00           C
ATOM      7  C4'   U A   1       1.197   0.415  -0.942  1.00 10.00           C
ATOM      8  O4'   U A   1       0.194   0.645   0.074  1.00 10.00           O
ATOM      9  C3'   U A   1       2.181  -0.588  -0.301  1.00 10.00           C
ATOM     10  O3'   U A   1       3.524  -0.288  -0.686  1.00 10.00           O
ATOM     11  C2'   U A   1       1.995  -0.383   1.218  1.00 10.00           C
ATOM     12  O2'   U A   1       3.219   0.046   1.819  1.00 10.00           O
ATOM     13  C1'   U A   1       0.922   0.723   1.319  1.00 10.00           C
ATOM     14  N1    U A   1       0.028   0.464   2.451  1.00 10.00           N
ATOM     15  C2    U A   1      -0.690  -0.671   2.486  1.00 10.00           C
ATOM     16  O2    U A   1      -0.587  -1.474   1.580  1.00 10.00           O
ATOM     17  N3    U A   1      -1.515  -0.936   3.517  1.00 10.00           N
ATOM     18  C4    U A   1      -1.641  -0.055   4.530  1.00 10.00           C
ATOM     19  O4    U A   1      -2.391  -0.292   5.460  1.00 10.00           O
ATOM     20  C5    U A   1      -0.894   1.146   4.502  1.00 10.00           C
ATOM     21  C6    U A   1      -0.070   1.384   3.459  1.00 10.00           C
ATOM     22 HOP3   U A   1      -2.475   0.583  -5.470  1.00 10.00           H
ATOM     23 HOP2   U A   1       0.534   0.613  -5.325  1.00 10.00           H
ATOM     24  H5'   U A   1       1.329  -0.373  -2.942  1.00 10.00           H
ATOM     25 H5''   U A   1       0.060  -1.117  -1.940  1.00 10.00           H
ATOM     26  H4'   U A   1       1.712   1.345  -1.185  1.00 10.00           H
ATOM     27  H3'   U A   1       1.923  -1.609  -0.583  1.00 10.00           H
ATOM     28 HO3'   U A   1       4.094  -0.926  -0.234  1.00 10.00           H
ATOM     29  H2'   U A   1       1.643  -1.301   1.688  1.00 10.00           H
ATOM     30 HO2'   U A   1       3.865  -0.657   1.671  1.00 10.00           H
ATOM     31  H1'   U A   1       1.392   1.700   1.423  1.00 10.00           H
ATOM     32  H3    U A   1      -2.024  -1.762   3.528  1.00 10.00           H
ATOM     33  H5    U A   1      -0.982   1.863   5.305  1.00 10.00           H
ATOM     34  H6    U A   1       0.507   2.295   3.421  1.00 10.00           H
END   
"""


@pytest.fixture()
def A(tmpdir):
    fn = tmpdir / 'A.pdb'
    content = A_pdb

    with open(fn, 'w') as out:
        out.write(content)

    yield str(fn)


@pytest.fixture()
def C(tmpdir):
    fn = tmpdir / 'C.pdb'
    content = C_pdb

    with open(fn, 'w') as out:
        out.write(content)

    yield str(fn)


@pytest.fixture()
def G(tmpdir):
    fn = tmpdir / 'C.pdb'
    content = G_pdb

    with open(fn, 'w') as out:
        out.write(content)

    yield str(fn)


@pytest.fixture()
def U(tmpdir):
    fn = tmpdir / 'U.pdb'
    content = U_pdb

    with open(fn, 'w') as out:
        out.write(content)

    yield str(fn)


def test_A(A):
    m = pdbinf.load_pdb_file(A, templates=[pdbinf.RNA_DOC])

    assert m
    assert m.GetNumAtoms() == 37
    assert m.GetNumBonds() == 39  # 3 rings


def test_C(C):
    m = pdbinf.load_pdb_file(C, templates=[pdbinf.RNA_DOC])

    assert m
    assert m.GetNumAtoms() == 35
    assert m.GetNumBonds() == 36  # 2 rings


def test_G(G):
    m = pdbinf.load_pdb_file(G, templates=[pdbinf.RNA_DOC])

    assert m
    assert m.GetNumAtoms() == 38
    assert m.GetNumBonds() == 40  # 3 rings


def test_U(U):
    m = pdbinf.load_pdb_file(U, templates=[pdbinf.RNA_DOC])

    assert m
    assert m.GetNumAtoms() == 34
    assert m.GetNumBonds() == 35  # 2 rings
