# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
import gemmi


_STANDARD_AA = """\
data_ALA
#
_chem_comp.id ALA
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
ALA  N    N    N  N  N  N
ALA  CA   CA   C  N  N  S
ALA  C    C    C  N  N  N
ALA  O    O    O  N  N  N
ALA  CB   CB   C  N  N  N
ALA  OXT  OXT  O  N  Y  N
ALA  HA   HA   H  N  N  N
ALA  HB1  1HB  H  N  N  N
ALA  HB2  2HB  H  N  N  N
ALA  HB3  3HB  H  N  N  N
ALA  H    H    H  N  N  N
ALA  H2   H2   H  N  N  N
ALA  H3   H3   H  N  Y  N
ALA  HXT  HXT  H  N  Y  N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
N   H3  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  HB1 SING N N
CB  HB2 SING N N
CB  HB3 SING N N
OXT HXT SING N N
##
data_NME
#
_chem_comp.id NME
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
NME  N    N    N  N  N  N
NME  C    C    C  N  N  N
NME  H    H    H  N  N  N
NME  H1   H1   H  N  N  N
NME  H2   H2   H  N  N  N
NME  H3   H3   H  N  N  N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   C   SING N N
N   H   SING N N
C   H1  SING N N
C   H2  SING N N
C   H3  SING N N
##
data_ACE
_chem_comp.id ACE
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
ACE  CH3  CH3  C  N  N  N
ACE  C    C    C  N  N  N
ACE  O    O    O  N  N  N
ACE  H1   H1   H  N  N  N
ACE  H2   H2   H  N  N  N
ACE  H3   H3   H  N  N  N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
CH3 C   SING N N
C   O   DOUB N N
CH3 H1  SING N N
CH3 H2  SING N N
CH3 H3  SING N N
##
data_ARG
#
_chem_comp.id ARG
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
ARG N    N    N N N N
ARG CA   CA   C N N S
ARG C    C    C N N N
ARG O    O    O N N N
ARG CB   CB   C N N N
ARG CG   CG   C N N N
ARG CD   CD   C N N N
ARG NE   NE   N N N N
ARG CZ   CZ   C N N N
ARG NH1  NH1  N N N N
ARG NH2  NH2  N N N N
ARG OXT  OXT  O N Y N
ARG H    H    H N N N
ARG H2   HN2  H N Y N
ARG HA   HA   H N N N
ARG HB2  1HB  H N N N
ARG HB3  2HB  H N N N
ARG HG2  1HG  H N N N
ARG HG3  2HG  H N N N
ARG HD2  1HD  H N N N
ARG HD3  2HD  H N N N
ARG HE   HE   H N N N
ARG HH11 1HH1 H N N N
ARG HH12 2HH1 H N N N
ARG HH21 1HH2 H N N N
ARG HH22 2HH2 H N N N
ARG HXT  HXT  H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA   SING N N
N   H    SING N N
N   H2   SING N N
N   H3   SING N N
CA  C    SING N N
CA  CB   SING N N
CA  HA   SING N N
C   O    DOUB N N
C   OXT  SING N N
CB  CG   SING N N
CB  HB2  SING N N
CB  HB3  SING N N
CG  CD   SING N N
CG  HG2  SING N N
CG  HG3  SING N N
CD  NE   SING N N
CD  HD2  SING N N
CD  HD3  SING N N
NE  CZ   SING N N
NE  HE   SING N N
CZ  NH1  SING N N
CZ  NH2  DOUB N N
NH1 HH11 SING N N
NH1 HH12 SING N N
NH2 HH21 SING N N
NH2 HH22 SING N N
OXT HXT  SING N N
##
data_TRP
#
_comp_chem.id TRP
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C Y N N
CD1 CD1 C Y N N
CD2 CD2 C Y N N
NE1 NE1 N Y N N
CE2 CE2 C Y N N
CE3 CE3 C Y N N
CZ2 CZ2 C Y N N
CZ3 CZ3 C Y N N
CH2 CH2 C Y N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HD1 HD1 H N N N
HE1 HE1 H N N N
HE3 HE3 H N N N
HZ2 HZ2 H N N N
HZ3 HZ3 H N N N
HH2 HH2 H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  CD1 DOUB Y N
CG  CD2 SING Y N
CD1 NE1 SING Y N
CD1 HD1 SING N N
CD2 CE2 DOUB Y N
CD2 CE3 SING Y N
NE1 CE2 SING Y N
NE1 HE1 SING N N
CE2 CZ2 SING Y N
CE3 CZ3 DOUB Y N
CE3 HE3 SING N N
CZ2 CH2 DOUB Y N
CZ2 HZ2 SING N N
CZ3 CH2 SING Y N
CZ3 HZ3 SING N N
CH2 HH2 SING N N
OXT HXT SING N N
##
data_PRO
#
_chem_comp.id PRO
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C N N N
CD  CD  C N N N
OXT OXT O N Y N
H   HT1 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HG2 1HG H N N N
HG3 2HG H N N N
HD2 1HD H N N N
HD3 2HD H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   CD  SING N N
N   H   SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  CD  SING N N
CG  HG2 SING N N
CG  HG3 SING N N
CD  HD2 SING N N
CD  HD3 SING N N
OXT HXT SING N N
##
data_CYS
#
_chem_comp.id CYS
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N R
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
SG  SG  S N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HG  HG  H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  SG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
SG  HG  SING N N
OXT HXT SING N N
##
data_LYS
#
_chem_comp.id LYS
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C N N N
CD  CD  C N N N
CE  CE  C N N N
NZ  NZ  N N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HG2 1HG H N N N
HG3 2HG H N N N
HD2 1HD H N N N
HD3 2HD H N N N
HE2 1HE H N N N
HE3 2HE H N N N
HZ1 1HZ H N N N
HZ2 2HZ H N N N
HZ3 3HZ H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  CD  SING N N
CG  HG2 SING N N
CG  HG3 SING N N
CD  CE  SING N N
CD  HD2 SING N N
CD  HD3 SING N N
CE  NZ  SING N N
CE  HE2 SING N N
CE  HE3 SING N N
NZ  HZ1 SING N N
NZ  HZ2 SING N N
NZ  HZ3 SING N N
OXT HXT SING N N
##
data_ASP
#
_chem_comp.id ASP
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C N N N
OD1 OD1 O N N N
OD2 OD2 O N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 HB1 H N N N
HB3 HB2 H N N N
HD2 HD2 H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N  
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  OD1 DOUB N N
CG  OD2 SING N N
OD2 HD2 SING N N
OXT HXT SING N N
##
data_ASN
#
_comp_chem.id ASN
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N    N    N N N N
CA   CA   C N N S
C    C    C N N N
O    O    O N N N
CB   CB   C N N N
CG   CG   C N N N
OD1  OD1  O N N N
ND2  ND2  N N N N
OXT  OXT  O N Y N
H    H    H N N N
H2   HN2  H N Y N
HA   HA   H N N N
HB2  1HB  H N N N
HB3  2HB  H N N N
HD21 1HD2 H N N N
HD22 2HD2 H N N N
HXT  HXT  H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA   SING N N
N   H    SING N N
N   H2   SING N N
CA  C    SING N N
CA  CB   SING N N
CA  HA   SING N N
C   O    DOUB N N
C   OXT  SING N N
CB  CG   SING N N
CB  HB2  SING N N
CB  HB3  SING N N
CG  OD1  DOUB N N
CG  ND2  SING N N
ND2 HD21 SING N N
ND2 HD22 SING N N
OXT HXT  SING N N
##
data_GLU
#
_chem_comp.id GLU
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C N N N
CD  CD  C N N N
OE1 OE1 O N N N
OE2 OE2 O N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 HB1 H N N N
HB3 HB2 H N N N
HG2 HG1 H N N N
HG3 HG2 H N N N
HE2 HE2 H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  CD  SING N N
CG  HG2 SING N N
CG  HG3 SING N N
CD  OE1 DOUB N N
CD  OE2 SING N N
OE2 HE2 SING N N
OXT HXT SING N N
##
data_GLN
#
_comp_chem.id GLN
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N    N    N N N N
CA   CA   C N N S
C    C    C N N N
O    O    O N N N
CB   CB   C N N N
CG   CG   C N N N
CD   CD   C N N N
OE1  OE1  O N N N
NE2  NE2  N N N N
OXT  OXT  O N Y N
H    H    H N N N
H2   HN2  H N Y N
HA   HA   H N N N
HB2  1HB  H N N N
HB3  2HB  H N N N
HG2  1HG  H N N N
HG3  2HG  H N N N
HE21 1HE2 H N N N
HE22 2HE2 H N N N
HXT  HXT  H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA   SING N N
N   H    SING N N
N   H2   SING N N
CA  C    SING N N
CA  CB   SING N N
CA  HA   SING N N
C   O    DOUB N N
C   OXT  SING N N
CB  CG   SING N N
CB  HB2  SING N N
CB  HB3  SING N N
CG  CD   SING N N
CG  HG2  SING N N
CG  HG3  SING N N
CD  OE1  DOUB N N
CD  NE2  SING N N
NE2 HE21 SING N N
NE2 HE22 SING N N
OXT HXT  SING N N
##
data_GLY
#
_comp_chem.id GLY
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N N
C   C   C N N N
O   O   O N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA2 HA1 H N N N
HA3 HA2 H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  HA2 SING N N
CA  HA3 SING N N
C   O   DOUB N N
C   OXT SING N N
OXT HXT SING N N
##
data_ILE
#
_comp_chem.id ILE
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N    N    N N N N
CA   CA   C N N S
C    C    C N N N
O    O    O N N N
CB   CB   C N N S
CG1  CG1  C N N N
CG2  CG2  C N N N
CD1  CD1  C N N N
OXT  OXT  O N Y N
H    H    H N N N
H2   HN2  H N Y N
HA   HA   H N N N
HB   HB   H N N N
HG12 1HG1 H N N N
HG13 2HG1 H N N N
HG21 1HG2 H N N N
HG22 2HG2 H N N N
HG23 3HG2 H N N N
HD11 1HD1 H N N N
HD12 2HD1 H N N N
HD13 3HD1 H N N N
HXT  HXT  H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA   SING N N
N   H    SING N N
N   H2   SING N N
CA  C    SING N N
CA  CB   SING N N
CA  HA   SING N N
C   O    DOUB N N
C   OXT  SING N N
CB  CG1  SING N N
CB  CG2  SING N N
CB  HB   SING N N
CG1 CD1  SING N N
CG1 HG12 SING N N
CG1 HG13 SING N N
CG2 HG21 SING N N
CG2 HG22 SING N N
CG2 HG23 SING N N
CD1 HD11 SING N N
CD1 HD12 SING N N
CD1 HD13 SING N N
OXT HXT  SING N N
##
data_LEU
#
_chem_comp.id LEU
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N    N    N N N N
CA   CA   C N N S
C    C    C N N N
O    O    O N N N
CB   CB   C N N N
CG   CG   C N N N
CD1  CD1  C N N N
CD2  CD2  C N N N
OXT  OXT  O N Y N
H    H    H N N N
H2   HN2  H N Y N
HA   HA   H N N N
HB2  1HB  H N N N
HB3  2HB  H N N N
HG   HG   H N N N
HD11 1HD1 H N N N
HD12 2HD1 H N N N
HD13 3HD1 H N N N
HD21 1HD2 H N N N
HD22 2HD2 H N N N
HD23 3HD2 H N N N
HXT  HXT  H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA   SING N N
N   H    SING N N
N   H2   SING N N
CA  C    SING N N
CA  CB   SING N N
CA  HA   SING N N
C   O    DOUB N N
C   OXT  SING N N
CB  CG   SING N N
CB  HB2  SING N N
CB  HB3  SING N N
CG  CD1  SING N N
CG  CD2  SING N N
CG  HG   SING N N
CD1 HD11 SING N N
CD1 HD12 SING N N
CD1 HD13 SING N N
CD2 HD21 SING N N
CD2 HD22 SING N N
CD2 HD23 SING N N
OXT HXT  SING N N
##
data_MET
#
_chem_comp.id MET
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C N N N
SD  SD  S N N N
CE  CE  C N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HG2 1HG H N N N
HG3 2HG H N N N
HE1 1HE H N N N
HE2 2HE H N N N
HE3 3HE H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  SD  SING N N
CG  HG2 SING N N
CG  HG3 SING N N
SD  CE  SING N N
CE  HE1 SING N N
CE  HE2 SING N N
CE  HE3 SING N N
OXT HXT SING N N
##
data_PHE
#
_chem_comp.id PHE
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C Y N N
CD1 CD1 C Y N N
CD2 CD2 C Y N N
CE1 CE1 C Y N N
CE2 CE2 C Y N N
CZ  CZ  C Y N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HD1 HD1 H N N N
HD2 HD2 H N N N
HE1 HE1 H N N N
HE2 HE2 H N N N
HZ  HZ  H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  CD1 DOUB Y N
CG  CD2 SING Y N
CD1 CE1 SING Y N
CD1 HD1 SING N N
CD2 CE2 DOUB Y N
CD2 HD2 SING N N
CE1 CZ  DOUB Y N
CE1 HE1 SING N N
CE2 CZ  SING Y N
CE2 HE2 SING N N
CZ  HZ  SING N N
OXT HXT SING N N
##
data_SER
#
_chem_comp.id SER
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
OG  OG  O N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HG  HG  H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  OG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
OG  HG  SING N N
OXT HXT SING N N
##
data_THR
_chem_comp.id THR
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N    N    N N N N
CA   CA   C N N S
C    C    C N N N
O    O    O N N N
CB   CB   C N N R
OG1  OG1  O N N N
CG2  CG2  C N N N
OXT  OXT  O N Y N
H    H    H N N N
H2   HN2  H N Y N
HA   HA   H N N N
HB   HB   H N N N
HG1  HG1  H N N N
HG21 1HG2 H N N N
HG22 2HG2 H N N N
HG23 3HG2 H N N N
HXT  HXT  H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA   SING N N
N   H    SING N N
N   H2   SING N N
CA  C    SING N N
CA  CB   SING N N
CA  HA   SING N N
C   O    DOUB N N
C   OXT  SING N N
CB  OG1  SING N N
CB  CG2  SING N N
CB  HB   SING N N
OG1 HG1  SING N N
CG2 HG21 SING N N
CG2 HG22 SING N N
CG2 HG23 SING N N
OXT HXT  SING N N
##
data_TYR
#
_chem_comp.id TYR
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C Y N N
CD1 CD1 C Y N N
CD2 CD2 C Y N N
CE1 CE1 C Y N N
CE2 CE2 C Y N N
CZ  CZ  C Y N N
OH  OH  O N N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HD1 HD1 H N N N
HD2 HD2 H N N N
HE1 HE1 H N N N
HE2 HE2 H N N N
HH  HH  H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  CD1 DOUB Y N
CG  CD2 SING Y N
CD1 CE1 SING Y N
CD1 HD1 SING N N
CD2 CE2 DOUB Y N
CD2 HD2 SING N N
CE1 CZ  DOUB Y N
CE1 HE1 SING N N
CE2 CZ  SING Y N
CE2 HE2 SING N N
CZ  OH  SING N N
OH  HH  SING N N
OXT HXT SING N N
##
data_VAL
#
_chem_comp.id VAL
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N    N    N N N N
CA   CA   C N N S
C    C    C N N N
O    O    O N N N
CB   CB   C N N N
CG1  CG1  C N N N
CG2  CG2  C N N N
OXT  OXT  O N Y N
H    H    H N N N
H2   HN2  H N Y N
HA   HA   H N N N
HB   HB   H N N N
HG11 1HG1 H N N N
HG12 2HG1 H N N N
HG13 3HG1 H N N N
HG21 1HG2 H N N N
HG22 2HG2 H N N N
HG23 3HG2 H N N N
HXT  HXT  H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA   SING N N
N   H    SING N N
N   H2   SING N N
CA  C    SING N N
CA  CB   SING N N
CA  HA   SING N N
C   O    DOUB N N
C   OXT  SING N N
CB  CG1  SING N N
CB  CG2  SING N N
CB  HB   SING N N
CG1 HG11 SING N N
CG1 HG12 SING N N
CG1 HG13 SING N N
CG2 HG21 SING N N
CG2 HG22 SING N N
CG2 HG23 SING N N
OXT HXT  SING N N
##
data_HIS
#
_chem_comp.id HIS
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C Y N N
ND1 ND1 N Y N N
CD2 CD2 C Y N N
CE1 CE1 C Y N N
NE2 NE2 N Y N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HD1 HD1 H N N N
HD2 HD2 H N N N
HE1 HE1 H N N N
HE2 HE2 H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  ND1 SING Y N
CG  CD2 DOUB Y N
ND1 CE1 DOUB Y N
ND1 HD1 SING N N
CD2 NE2 SING Y N
CD2 HD2 SING N N
CE1 NE2 SING Y N
CE1 HE1 SING N N
NE2 HE2 SING N N
OXT HXT SING N N
##
data_HID
#
_chem_comp.id HID
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
N   N   N N N N
CA  CA  C N N S
C   C   C N N N
O   O   O N N N
CB  CB  C N N N
CG  CG  C Y N N
ND1 ND1 N Y N N
CD2 CD2 C Y N N
CE1 CE1 C Y N N
NE2 NE2 N Y N N
OXT OXT O N Y N
H   H   H N N N
H2  HN2 H N Y N
HA  HA  H N N N
HB2 1HB H N N N
HB3 2HB H N N N
HD1 HD1 H N N N
HD2 HD2 H N N N
HE1 HE1 H N N N
HE2 HE2 H N N N
HXT HXT H N Y N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
N   CA  SING N N
N   H   SING N N
N   H2  SING N N
CA  C   SING N N
CA  CB  SING N N
CA  HA  SING N N
C   O   DOUB N N
C   OXT SING N N
CB  CG  SING N N
CB  HB2 SING N N
CB  HB3 SING N N
CG  ND1 SING Y N
CG  CD2 DOUB Y N
ND1 CE1 SING Y N
ND1 HD1 SING N N
CD2 NE2 SING Y N
CD2 HD2 SING N N
CE1 NE2 DOUB Y N
CE1 HE1 SING N N
NE2 HE2 SING N N
OXT HXT SING N N
##
data_HOH
#
_chem_comp.id HOH
#
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
O  O  O N N N
H1 H1 H N N N
H2 H2 H N N N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
O  H1 SING N N
O  H2 SING N N
##
"""

STANDARD_AA_DOC = gemmi.cif.read_string(_STANDARD_AA)
