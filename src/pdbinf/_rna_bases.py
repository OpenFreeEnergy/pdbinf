# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf

import gemmi


_RNA_BASES = """\
data_A
_chem_comp.id A
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
OP3 O N
P P N
OP1 O N
OP2 O N
"O5'" O N
"C5'" C N
"C4'" C N
"O4'" O N
"C3'" C N
"O3'" O N
"C2'" C N
"O2'" O N
"C1'" C N
N9 N Y
C8 C Y
N7 N Y
C5 C Y
C6 C Y
N6 N N
N1 N Y
C2 C Y
N3 N Y
C4 C Y
HOP3 H N
HOP2 H N
"H5'" H N
"H5''" H N
"H4'" H N
"H3'" H N
"HO3'" H N
"H2'" H N
"HO2'" H N
"H1'" H N
H8 H N
H61 H N
H62 H N
H2 H N
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
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
OP3 O N
P P N
OP1 O N
OP2 O N
"O5'" O N
"C5'" C N
"C4'" C N
"O4'" O N
"C3'" C N
"O3'" O N
"C2'" C N
"O2'" O N
"C1'" C N
N1 N N
C2 C N
O2 O N
N3 N N
C4 C N
N4 N N
C5 C N
C6 C N
HOP3 H N
HOP2 H N
"H5'" H N
"H5''" H N
"H4'" H N
"H3'" H N
"HO3'" H N
"H2'" H N
"HO2'" H N
"H1'" H N
H41 H N
H42 H N
H5 H N
H6 H N
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
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
OP3 O N
P P N
OP1 O N
OP2 O N
"O5'" O N
"C5'" C N
"C4'" C N
"O4'" O N
"C3'" C N
"O3'" O N
"C2'" C N
"O2'" O N
"C1'" C N
N9 N Y
C8 C Y
N7 N Y
C5 C Y
C6 C N
O6 O N
N1 N N
C2 C N
N2 N N
N3 N N
C4 C Y
HOP3 H N
HOP2 H N
"H5'" H N
"H5''" H N
"H4'" H N
"H3'" H N
"HO3'" H N
"H2'" H N
"HO2'" H N
"H1'" H N
H8 H N
H1 H N
H21 H N
H22 H N
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
data_U
_chem_comp.id U
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
OP3 O N
P P N
OP1 O N
OP2 O N
"O5'" O N
"C5'" C N
"C4'" C N
"O4'" O N
"C3'" C N
"O3'" O N
"C2'" C N
"O2'" O N
"C1'" C N
N1 N N
C2 C N
O2 O N
N3 N N
C4 C N
O4 O N
C5 C N
C6 C N
HOP3 H N
HOP2 H N
"H5'" H N
"H5''" H N
"H4'" H N
"H3'" H N
"HO3'" H N
"H2'" H N
"HO2'" H N
"H1'" H N
H3 H N
H5 H N
H6 H N
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
N3 C4 SING N
N3 H3 SING N
C4 O4 DOUB N
C4 C5 SING N
C5 C6 DOUB N
C5 H5 SING N
C6 H6 SING N
##
"""

RNA_DOC = gemmi.cif.read_string(_RNA_BASES)