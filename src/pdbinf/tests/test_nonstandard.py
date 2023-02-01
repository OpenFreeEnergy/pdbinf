import gemmi
import pdbinf
from rdkit import Chem
import pytest

# copied straight from ccd
tpo_template = """\
data_TPO
#

_chem_comp.id                                   TPO
_chem_comp.name                                 PHOSPHOTHREONINE
_chem_comp.type                                 "L-PEPTIDE LINKING"
_chem_comp.pdbx_type                            ATOMP
_chem_comp.formula                              "C4 H10 N O6 P"
_chem_comp.mon_nstd_parent_comp_id              THR
_chem_comp.pdbx_synonyms                        PHOSPHONOTHREONINE
_chem_comp.pdbx_formal_charge                   0
_chem_comp.pdbx_initial_date                    1999-07-08
_chem_comp.pdbx_modified_date                   2020-06-17
_chem_comp.pdbx_ambiguous_flag                  N
_chem_comp.pdbx_release_status                  REL
_chem_comp.pdbx_replaced_by                     ?
_chem_comp.pdbx_replaces                        ?
_chem_comp.formula_weight                       199.099
_chem_comp.one_letter_code                      T
_chem_comp.three_letter_code                    TPO
_chem_comp.pdbx_model_coordinates_details       ?
_chem_comp.pdbx_model_coordinates_missing_flag  N
_chem_comp.pdbx_ideal_coordinates_details       ?
_chem_comp.pdbx_ideal_coordinates_missing_flag  N
_chem_comp.pdbx_model_coordinates_db_code       1FMO
_chem_comp.pdbx_subcomponent_list               ?
_chem_comp.pdbx_processing_site                 EBI
#   #
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_align
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
_chem_comp_atom.pdbx_component_atom_id
_chem_comp_atom.pdbx_component_comp_id
_chem_comp_atom.pdbx_ordinal
TPO  N     N     N  0  1  N  N  N  21.891  2.133  -14.748   1.153  -1.040   2.377  N     TPO   1  
TPO  CA    CA    C  0  1  N  N  S  22.318  2.994  -13.673   0.572   0.199   1.844  CA    TPO   2  
TPO  CB    CB    C  0  1  N  N  R  21.313  4.075  -13.361   1.111   0.449   0.434  CB    TPO   3  
TPO  CG2   CG2   C  0  1  N  N  N  21.837  5.045  -12.302   2.634   0.580   0.485  CG2   TPO   4  
TPO  OG1   OG1   O  0  1  N  N  N  20.898  4.716  -14.523   0.755  -0.645  -0.412  OG1   TPO   5  
TPO  P     P     P  0  1  N  N  N  19.424  4.424  -14.993  -0.142  -0.039  -1.603  P     TPO   6  
TPO  O1P   O1P   O  0  1  N  N  N  19.358  5.014  -16.321   0.644   0.968  -2.350  O1P   TPO   7  
TPO  O2P   O2P   O  0  1  N  N  N  19.243  2.986  -14.834  -0.580  -1.224  -2.601  O2P   TPO   8  
TPO  O3P   O3P   O  0  1  N  N  N  18.506  5.082  -14.021  -1.456   0.656  -0.985  O3P   TPO   9  
TPO  C     C     C  0  1  N  N  N  22.539  2.278  -12.384  -0.927   0.070   1.794  C     TPO  10  
TPO  O     O     O  0  1  N  N  N  21.778  1.390  -12.005  -1.435  -1.012   1.626  O     TPO  11  
TPO  OXT   OXT   O  0  1  N  Y  N  23.582  2.721  -11.720  -1.700   1.159   1.935  OXT   TPO  12  
TPO  H     H     H  0  1  N  N  N  22.570  1.402  -14.958   2.154  -0.949   2.296  H     TPO  13  
TPO  H2    2HN   H  0  1  N  Y  N  21.663  2.673  -15.582   0.877  -1.782   1.751  H2    TPO  14  
TPO  HA    HA    H  0  1  N  N  N  23.275  3.418  -14.056   0.844   1.034   2.490  HA    TPO  15  
TPO  HB    HB    H  0  1  N  N  N  20.410  3.593  -12.916   0.680   1.369   0.039  HB    TPO  16  
TPO  HG21  1HG2  H  0  0  N  N  N  21.094  5.844  -12.071   3.065  -0.339   0.881  HG21  TPO  17  
TPO  HG22  2HG2  H  0  0  N  N  N  22.154  4.506  -11.378   3.018   0.758  -0.518  HG22  TPO  18  
TPO  HG23  3HG2  H  0  0  N  N  N  22.821  5.477  -12.598   2.906   1.415   1.131  HG23  TPO  19  
TPO  HOP2  2HOP  H  0  0  N  N  N  18.353  2.809  -15.117  -1.114  -0.819  -3.298  HOP2  TPO  20  
TPO  HOP3  3HOP  H  0  0  N  N  N  17.616  4.905  -14.304  -1.938  -0.033  -0.509  HOP3  TPO  21  
TPO  HXT   HXT   H  0  1  N  Y  N  23.722  2.264  -10.898  -2.662   1.076   1.902  HXT   TPO  22  
#   #
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
TPO  N    CA    SING  N  N   1  
TPO  N    H     SING  N  N   2  
TPO  N    H2    SING  N  N   3  
TPO  CA   CB    SING  N  N   4  
TPO  CA   C     SING  N  N   5  
TPO  CA   HA    SING  N  N   6  
TPO  CB   CG2   SING  N  N   7  
TPO  CB   OG1   SING  N  N   8  
TPO  CB   HB    SING  N  N   9  
TPO  CG2  HG21  SING  N  N  10  
TPO  CG2  HG22  SING  N  N  11  
TPO  CG2  HG23  SING  N  N  12  
TPO  OG1  P     SING  N  N  13  
TPO  P    O1P   DOUB  N  N  14  
TPO  P    O2P   SING  N  N  15  
TPO  P    O3P   SING  N  N  16  
TPO  O2P  HOP2  SING  N  N  17  
TPO  O3P  HOP3  SING  N  N  18  
TPO  C    O     DOUB  N  N  19  
TPO  C    OXT   SING  N  N  20  
TPO  OXT  HXT   SING  N  N  21  
#   #
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
TPO  SMILES            ACDLabs               10.04  "O=P(O)(O)OC(C(N)C(=O)O)C"  
TPO  SMILES_CANONICAL  CACTVS                3.341  "C[C@@H](O[P](O)(O)=O)[C@H](N)C(O)=O"  
TPO  SMILES            CACTVS                3.341  "C[CH](O[P](O)(O)=O)[CH](N)C(O)=O"  
TPO  SMILES_CANONICAL  "OpenEye OEToolkits"  1.5.0  "C[C@H]([C@@H](C(=O)O)N)OP(=O)(O)O"  
TPO  SMILES            "OpenEye OEToolkits"  1.5.0  "CC(C(C(=O)O)N)OP(=O)(O)O"  
TPO  InChI             InChI                 1.03   "InChI=1S/C4H10NO6P/c1-2(3(5)4(6)7)11-12(8,9)10/h2-3H,5H2,1H3,(H,6,7)(H2,8,9,10)/t2-,3+/m1/s1"  
TPO  InChIKey          InChI                 1.03   USRGIUJOYOXOQJ-GBXIJSLDSA-N  
#   #
loop_
_pdbx_chem_comp_identifier.comp_id
_pdbx_chem_comp_identifier.type
_pdbx_chem_comp_identifier.program
_pdbx_chem_comp_identifier.program_version
_pdbx_chem_comp_identifier.identifier
TPO  "SYSTEMATIC NAME"  ACDLabs               10.04  O-phosphono-L-threonine  
TPO  "SYSTEMATIC NAME"  "OpenEye OEToolkits"  1.5.0  "(2S,3R)-2-amino-3-phosphonooxy-butanoic acid"  
#   #
loop_
_pdbx_chem_comp_audit.comp_id
_pdbx_chem_comp_audit.action_type
_pdbx_chem_comp_audit.date
_pdbx_chem_comp_audit.processing_site
TPO  "Create component"   1999-07-08  EBI   
TPO  "Modify descriptor"  2011-06-04  RCSB  
TPO  "Modify synonyms"    2020-06-05  PDBE  
#
_pdbx_chem_comp_synonyms.ordinal     1
_pdbx_chem_comp_synonyms.comp_id     TPO
_pdbx_chem_comp_synonyms.name        PHOSPHONOTHREONINE
_pdbx_chem_comp_synonyms.provenance  ?
_pdbx_chem_comp_synonyms.type        ?
##
"""

tpo = """\
REMARK   1 PDBFIXER FROM: tpo.pdb
REMARK   1 CREATED WITH OPENMM 7.7, 2023-01-19
HETATM    1  H1  ACE     1      -1.861   4.051   1.454  1.00  0.00           H  
HETATM    2  CH3 ACE     1      -1.252   3.158   1.245  1.00  0.00           C  
HETATM    3  H2  ACE     1      -1.461   2.423   2.061  1.00  0.00           H  
HETATM    4  H3  ACE     1      -0.170   3.423   1.346  1.00  0.00           H  
HETATM    5  C   ACE     1      -1.522   2.585  -0.102  1.00  0.00           C  
HETATM    6  O   ACE     1      -2.666   2.764  -0.584  1.00  0.00           O  
HETATM    7  N   TPO     2      -0.365   1.748  -0.928  1.00  0.00           N  
HETATM    8  CA  TPO     2      -0.159   0.351  -0.158  1.00  0.00           C  
HETATM    9  CB  TPO     2       1.315  -0.018   0.021  1.00  0.00           C  
HETATM   10  CG2 TPO     2       1.942   0.984   0.953  1.00  0.00           C  
HETATM   11  OG1 TPO     2       1.924  -0.131  -1.190  1.00  0.00           O  
HETATM   12  P   TPO     2       2.522  -1.636  -1.565  1.00  0.00           P  
HETATM   13  O1P TPO     2       2.285  -2.666  -0.492  1.00  0.00           O  
HETATM   14  O2P TPO     2       4.193  -1.597  -1.883  1.00  0.00           O  
HETATM   15  O3P TPO     2       1.769  -2.179  -2.981  1.00  0.00           O  
HETATM   16  C   TPO     2      -0.881  -0.709  -0.837  1.00  0.00           C  
HETATM   17  O   TPO     2      -1.555  -0.570  -1.856  1.00  0.00           O  
HETATM   18  H   TPO     2       0.320   2.242  -1.638  1.00  0.00           H  
HETATM   19  HA  TPO     2      -0.554   0.439   0.896  1.00  0.00           H  
HETATM   20  HB  TPO     2       1.344  -1.004   0.536  1.00  0.00           H  
HETATM   21 HG21 TPO     2       2.873   0.507   1.359  1.00  0.00           H  
HETATM   22 HG22 TPO     2       2.124   1.962   0.496  1.00  0.00           H  
HETATM   23 HG23 TPO     2       1.262   1.078   1.828  1.00  0.00           H  
HETATM   24  N   NME     3      -0.837  -2.343  -0.220  1.00  0.00           N  
HETATM   25  H   NME     3      -0.418  -3.247  -1.189  1.00  0.00           H  
HETATM   26  C   NME     3      -1.923  -2.794   0.887  1.00  0.00           C  
HETATM   27  H1  NME     3      -2.865  -1.683   1.480  1.00  0.00           H  
HETATM   28  H2  NME     3      -3.399  -3.145  -0.439  1.00  0.00           H  
HETATM   29  H3  NME     3      -1.984  -3.991   1.501  1.00  0.00           H  
TER      30      NME     3
END
"""


def test_tpo():
    template = gemmi.cif.read_string(tpo_template)

    m = Chem.MolFromPDBBlock(tpo, proximityBonding=False, removeHs=False)

    m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC, template])

    assert m.GetNumAtoms() == 29
    assert m.GetNumBonds() == 28
    for i, atom in enumerate(m.GetAtoms()):
        if i in [13, 14]:  # this is P-[O-]
            assert atom.GetFormalCharge() == -1
        else:
            assert atom.GetFormalCharge() == 0

        if i == 12:  # this is the P=O
            b = atom.GetBonds()

            assert len(b) == 1
            assert b[0].GetBondType() == Chem.BondType.DOUBLE


def test_tpo_fail():
    m = Chem.MolFromPDBBlock(tpo, proximityBonding=False, removeHs=False)

    with pytest.raises(ValueError, match="Failed to find template"):
        m = pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])

