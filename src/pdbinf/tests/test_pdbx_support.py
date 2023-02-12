# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf

import pdbinf
import pytest


pdbx_ALA_ALA = """\
data_DEPOSITED_XRAY
#

_cell.entry_id          DEPOSITED_XRAY
_cell.length_a          48.000
_cell.length_b          48.000
_cell.length_c          48.000
_cell.angle_alpha       90.00
_cell.angle_beta        90.00
_cell.angle_gamma       90.00
_cell.Z_PDB             1
_cell.pdbx_unique_axis  ?
#
_symmetry.entry_id                        DEPOSITED_XRAY
_symmetry.space_group_name_H-M            "P 1"
_symmetry.pdbx_full_space_group_name_H-M  ?
_symmetry.cell_setting                    ?
_symmetry.Int_Tables_number               ?
#
_atom_sites.entry_id                   DEPOSITED_XRAY
_atom_sites.fract_transf_matrix[1][1]  0.020833
_atom_sites.fract_transf_matrix[1][2]  0.000000
_atom_sites.fract_transf_matrix[1][3]  0.000000
_atom_sites.fract_transf_matrix[2][1]  0.000000
_atom_sites.fract_transf_matrix[2][2]  0.020833
_atom_sites.fract_transf_matrix[2][3]  0.000000
_atom_sites.fract_transf_matrix[3][1]  0.000000
_atom_sites.fract_transf_matrix[3][2]  0.000000
_atom_sites.fract_transf_matrix[3][3]  0.020833
_atom_sites.fract_transf_vector[1]     0.00000
_atom_sites.fract_transf_vector[2]     0.00000
_atom_sites.fract_transf_vector[3]     0.00000
#   #
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM   1  C  C     .  ACE  A  1  1  ?  27.376  26.711  24.184  1.00  0.00  ?  1  ACE  A  C     1  
HETATM   2  O  O     .  ACE  A  1  1  ?  28.219  26.023  24.745  1.00  0.00  ?  1  ACE  A  O     1  
HETATM   3  C  CH3   .  ACE  A  1  1  ?  25.965  26.197  24.023  1.00  0.00  ?  1  ACE  A  CH3   1  
HETATM   4  H  H1    .  ACE  A  1  1  ?  25.898  25.204  24.465  1.00  0.00  ?  1  ACE  A  H1    1  
HETATM   5  H  H2    .  ACE  A  1  1  ?  25.276  26.871  24.530  1.00  0.00  ?  1  ACE  A  H2    1  
HETATM   6  H  H3    .  ACE  A  1  1  ?  25.721  26.143  22.963  1.00  0.00  ?  1  ACE  A  H3    1  
ATOM     7  N  N     .  ALA  A  1  2  ?  27.627  27.923  23.696  1.00  0.00  ?  2  ALA  A  N     1  
ATOM     8  C  CA    .  ALA  A  1  2  ?  28.919  28.605  23.776  1.00  0.00  ?  2  ALA  A  CA    1  
ATOM     9  C  C     .  ALA  A  1  2  ?  28.737  30.134  23.791  1.00  0.00  ?  2  ALA  A  C     1  
ATOM    10  O  O     .  ALA  A  1  2  ?  27.696  30.645  23.377  1.00  0.00  ?  2  ALA  A  O     1  
ATOM    11  C  CB    .  ALA  A  1  2  ?  29.795  28.154  22.598  1.00  0.00  ?  2  ALA  A  CB    1  
ATOM    12  H  H     .  ALA  A  1  2  ?  26.870  28.444  23.279  1.00  0.00  ?  2  ALA  A  H     1  
ATOM    13  H  HA    .  ALA  A  1  2  ?  29.414  28.316  24.706  1.00  0.00  ?  2  ALA  A  HA    1  
ATOM    14  H  HB1   .  ALA  A  1  2  ?  29.321  28.429  21.654  1.00  0.00  ?  2  ALA  A  HB1   1  
ATOM    15  H  HB2   .  ALA  A  1  2  ?  30.774  28.629  22.659  1.00  0.00  ?  2  ALA  A  HB2   1  
ATOM    16  H  HB3   .  ALA  A  1  2  ?  29.928  27.071  22.630  1.00  0.00  ?  2  ALA  A  HB3   1  
ATOM    17  N  N     .  ALA  A  1  3  ?  29.757  30.863  24.249  1.00  0.00  ?  3  ALA  A  N     1  
ATOM    18  C  CA    .  ALA  A  1  3  ?  29.782  32.325  24.309  1.00  0.00  ?  3  ALA  A  CA    1  
ATOM    19  C  C     .  ALA  A  1  3  ?  31.224  32.854  24.202  1.00  0.00  ?  3  ALA  A  C     1  
ATOM    20  O  O     .  ALA  A  1  3  ?  32.178  32.113  24.438  1.00  0.00  ?  3  ALA  A  O     1  
ATOM    21  C  CB    .  ALA  A  1  3  ?  29.108  32.780  25.612  1.00  0.00  ?  3  ALA  A  CB    1  
ATOM    22  H  H     .  ALA  A  1  3  ?  30.604  30.392  24.534  1.00  0.00  ?  3  ALA  A  H     1  
ATOM    23  H  HA    .  ALA  A  1  3  ?  29.214  32.723  23.466  1.00  0.00  ?  3  ALA  A  HA    1  
ATOM    24  H  HB1   .  ALA  A  1  3  ?  29.658  32.391  26.470  1.00  0.00  ?  3  ALA  A  HB1   1  
ATOM    25  H  HB2   .  ALA  A  1  3  ?  29.091  33.869  25.662  1.00  0.00  ?  3  ALA  A  HB2   1  
ATOM    26  H  HB3   .  ALA  A  1  3  ?  28.082  32.412  25.646  1.00  0.00  ?  3  ALA  A  HB3   1  
HETATM  27  N  N     .  NME  A  1  4  ?  31.379  34.139  23.857  1.00  0.00  ?  4  NME  A  N     1  
HETATM  28  H  H     .  NME  A  1  4  ?  30.550  34.685  23.686  1.00  0.00  ?  4  NME  A  H     1  
HETATM  29  C  C     .  NME  A  1  4  ?  32.677  34.798  23.722  1.00  0.00  ?  4  NME  A  CH3   1  
HETATM  30  H  H1    .  NME  A  1  4  ?  33.332  34.199  23.086  1.00  0.00  ?  4  NME  A  HH31  1  
HETATM  31  H  H2    .  NME  A  1  4  ?  32.556  35.788  23.280  1.00  0.00  ?  4  NME  A  HH32  1  
HETATM  32  H  H3    .  NME  A  1  4  ?  33.144  34.897  24.704  1.00  0.00  ?  4  NME  A  HH33  1  
#
"""


@pytest.fixture
def ALA_ALA(tmpdir):
    p = str(tmpdir / 'ala.cif')

    with open(p, 'w') as fout:
        fout.write(pdbx_ALA_ALA)

    return p


def test_ALA_ALA(ALA_ALA):
    m = pdbinf.load_pdbx_file(ALA_ALA, templates=[pdbinf.STANDARD_AA_DOC])

    assert m.GetNumAtoms() == 32
    assert m.GetNumBonds() == 31
