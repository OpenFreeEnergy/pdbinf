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
"""

STANDARD_AA_DOC = gemmi.cif.read_string(_STANDARD_AA)
