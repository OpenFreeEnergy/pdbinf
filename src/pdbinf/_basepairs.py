# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
import gemmi

_BASEPAIRS = """\
data_A
_chem_comp.id A
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.pdbx_aromatic_flag
OP3 N
P N
OP1 N
OP2 N
"O5'" N
"C5'" N
"C4'" N
"O4'" N
"C3'" N
"O3'" N
"C2'" N
"O2'" N
"C1'" N
N9 Y
C8 Y
N7 Y
C5 Y
C6 Y
N6 N
N1 Y
C2 Y
N3 Y
C4 Y
HOP3 N
HOP2 N
"H5'" N
"H5''" N
"H4'" N
"H3'" N
"HO3'" N
"H2'" N
"HO2'" N
"H1'" N
H8 N
H61 N
H62 N
H2 N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
OP3 P SING N
OP3 HOP3 SING N
P OP1 DOUB N
P OP2 SING N
P "O5'" SING N
OP2 HOP2 SING N
"O5'" "C5'" SING N
"C5'" "C4'" SING N
"C5'" "H5'" SING N
"C5'" "H5''" SING N
"C4'" "O4'" SING N
"C4'" "C3'" SING N
"C4'" "H4'" SING N
"O4'" "C1'" SING N
"C3'" "O3'" SING N
"C3'" "C2'" SING N
"C3'" "H3'" SING N
"O3'" "HO3'" SING N
"C2'" "O2'" SING N
"C2'" "C1'" SING N
"C2'" "H2'" SING N
"O2'" "HO2'" SING N
"C1'" N9 SING N
"C1'" "H1'" SING N
N9 C8 SING Y
N9 C4 SING Y
C8 N7 DOUB Y
C8 H8 SING N
N7 C5 SING Y
C5 C6 SING Y
C5 C4 DOUB Y
C6 N6 SING N
C6 N1 DOUB Y
N6 H61 SING N
N6 H62 SING N
N1 C2 SING Y
C2 N3 DOUB Y
C2 H2 SING N
N3 C4 SING Y
##
data_C
_chem_comp.id C
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.pdbx_aromatic_flag
OP3 N
P N
OP1 N
OP2 N
"O5'" N
"C5'" N
"C4'" N
"O4'" N
"C3'" N
"O3'" N
"C2'" N
"O2'" N
"C1'" N
N1 N
C2 N
O2 N
N3 N
C4 N
N4 N
C5 N
C6 N
HOP3 N
HOP2 N
"H5'" N
"H5''" N
"H4'" N
"H3'" N
"HO3'" N
"H2'" N
"HO2'" N
"H1'" N
H41 N
H42 N
H5 N
H6 N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
OP3 P SING N
OP3 HOP3 SING N
P OP1 DOUB N
P OP2 SING N
P "O5'" SING N
OP2 HOP2 SING N
"O5'" "C5'" SING N
"C5'" "C4'" SING N
"C5'" "H5'" SING N
"C5'" "H5''" SING N
"C4'" "O4'" SING N
"C4'" "C3'" SING N
"C4'" "H4'" SING N
"O4'" "C1'" SING N
"C3'" "O3'" SING N
"C3'" "C2'" SING N
"C3'" "H3'" SING N
"O3'" "HO3'" SING N
"C2'" "O2'" SING N
"C2'" "C1'" SING N
"C2'" "H2'" SING N
"O2'" "HO2'" SING N
"C1'" N1 SING N
"C1'" "H1'" SING N
N1 C2 SING N
N1 C6 SING N
C2 O2 DOUB N
C2 N3 SING N
N3 C4 DOUB N
C4 N4 SING N
C4 C5 SING N
N4 H41 SING N
N4 H42 SING N
C5 C6 DOUB N
C5 H5 SING N
C6 H6 SING N
##
data_G
_chem_comp.id G
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.pdbx_aromatic_flag
OP3 N
P N
OP1 N
OP2 N
"O5'" N
"C5'" N
"C4'" N
"O4'" N
"C3'" N
"O3'" N
"C2'" N
"O2'" N
"C1'" N
N9 Y
C8 Y
N7 Y
C5 Y
C6 N
O6 N
N1 N
C2 N
N2 N
N3 N
C4 Y
HOP3 N
HOP2 N
"H5'" N
"H5''" N
"H4'" N
"H3'" N
"HO3'" N
"H2'" N
"HO2'" N
"H1'" N
H8 N
H1 N
H21 N
H22 N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
OP3 P SING N
OP3 HOP3 SING N
P OP1 DOUB N
P OP2 SING N
P "O5'" SING N
OP2 HOP2 SING N
"O5'" "C5'" SING N
"C5'" "C4'" SING N
"C5'" "H5'" SING N
"C5'" "H5''" SING N
"C4'" "O4'" SING N
"C4'" "C3'" SING N
"C4'" "H4'" SING N
"O4'" "C1'" SING N
"C3'" "O3'" SING N
"C3'" "C2'" SING N
"C3'" "H3'" SING N
"O3'" "HO3'" SING N
"C2'" "O2'" SING N
"C2'" "C1'" SING N
"C2'" "H2'" SING N
"O2'" "HO2'" SING N
"C1'" N9 SING N
"C1'" "H1'" SING N
N9 C8 SING Y
N9 C4 SING Y
C8 N7 DOUB Y
C8 H8 SING N
N7 C5 SING Y
C5 C6 SING N
C5 C4 DOUB Y
C6 O6 DOUB N
C6 N1 SING N
N1 C2 SING N
N1 H1 SING N
C2 N2 SING N
C2 N3 DOUB N
N2 H21 SING N
N2 H22 SING N
N3 C4 SING N
##
data_T
_chem_comp.id T
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.pdbx_aromatic_flag
OP3 N
P N
OP1 N
OP2 N
"O5'" N
"C5'" N
"C4'" N
"O4'" N
"C3'" N
"O3'" N
"C2'" N
"C1'" N
N1 N
C2 N
O2 N
N3 N
C4 N
O4 N
C5 N
C7 N
C6 N
HOP3 N
HOP2 N
"H5'" N
"H5''" N
"H4'" N
"H3'" N
"HO3'" N
"H2'" N
"H2''" N
"H1'" N
H3 N
H71 N
H72 N
H73 N
H6 N
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
OP3 P SING N
OP3 HOP3 SING N
P OP1 DOUB N
P OP2 SING N
P "O5'" SING N
OP2 HOP2 SING N
"O5'" "C5'" SING N
"C5'" "C4'" SING N
"C5'" "H5'" SING N
"C5'" "H5''" SING N
"C4'" "O4'" SING N
"C4'" "C3'" SING N
"C4'" "H4'" SING N
"O4'" "C1'" SING N
"C3'" "O3'" SING N
"C3'" "C2'" SING N
"C3'" "H3'" SING N
"O3'" "HO3'" SING N
"C2'" "C1'" SING N
"C2'" "H2'" SING N
"C2'" "H2''" SING N
"C1'" N1 SING N
"C1'" "H1'" SING N
N1 C2 SING N
N1 C6 SING N
C2 O2 DOUB N
C2 N3 SING N
N3 C4 SING N
N3 H3 SING N
C4 O4 DOUB N
C4 C5 SING N
C5 C7 SING N
C5 C6 DOUB N
C7 H71 SING N
C7 H72 SING N
C7 H73 SING N
C6 H6 SING N
##
"""

BASEPAIRS_DOC = gemmi.cif.read_string(_BASEPAIRS)
