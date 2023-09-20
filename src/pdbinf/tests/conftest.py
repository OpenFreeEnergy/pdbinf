# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
import gemmi
import pooch
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


@pytest.fixture(scope='session')
def tpo_doc():
    return gemmi.cif.read_string(tpo_template)


@pytest.fixture()
def openmm_nucleic():
    return pooch.retrieve(
        url="https://github.com/openmm/openmm/raw/8.0.0/wrappers/python/tests/systems/nucleic.pdb",
        known_hash="97eae94bffff1b8e524957477cd841f114ed58a99bbcdf8f518de31380892907",
    )


PDBS = pooch.create(
    path=pooch.os_cache('pdbinf'),
    base_url='https://github.com/openforcefield/protein-ligand-benchmark/raw/d3387602bbeb0167abf00dfb81753d8936775dd2/data/',
    version=None,
    registry={
        'p38/01_protein/crd/protein.pdb': '3f0bf718644e7c29f5200cd3def4240ac25ef5fb1948b2e64deb5015d8a45aa4',
        'mcl1/01_protein/crd/protein.pdb': 'f80ff9dd93a5d9dd6e90091e9631a8ce7fe0dc931e16543e22c1f92009660306',
        'cdk2/01_protein/crd/protein.pdb': '15d1e509d7951ca45ea266d51a627d5f452dcf0bb5bd48751ae57eb29e28ab69',
        'shp2/01_protein/crd/protein.pdb': 'd6759cbd135aaddaa658446064df4095d978d3681c014a0528b542d60b2c8770',
        'pde2/01_protein/crd/protein.pdb': '3b7967c1717789215452cdf919520625602d5438a9d2a18620726b8b1b3a8ef0',
        'cmet/01_protein/crd/protein.pdb': '155ec32941a9082dbdbbfde460ff97c88d4fe7e100e9a9577edb5a9e7b6467ae',
        'ptp1b/01_protein/crd/protein.pdb': 'bfa0f9204e96aa463b80946b788c4153cd24701291007eb77638a16fd156634e',
        'thrombin/01_protein/crd/protein.pdb': 'eb4ea18bef9c4c71dcdc922616d6719ee918112be87a0bd6b274c856eff1dd59',
        'cdk8/01_protein/crd/protein.pdb': 'b058774526a19775d8f438b14e9d6da331b6de74e0ef9e96db575f6c0bb067b2',
        'pfkfb3/01_protein/crd/protein.pdb': '4367710db0dbf284cc715ae9a8dd82d06bd77dcc3fb0885678e16632a2732dcc',
        'tyk2/01_protein/crd/protein.pdb': '9090684f4bdae90afbe5f2698a14c778396c024c19ceb6333de4808d9e29fae6',
        'syk/01_protein/crd/protein.pdb': 'f6199d0c1818eb5bb24e164426789cf39cae7aa32c8ca2e98f5f44d299a6f82f',
        'tnks2/01_protein/crd/protein.pdb': 'fc7681a05dbf07590aa8de133f981b6d8ae9cebcc23d54addc2c4fe80be80299',
        'eg5/01_protein/crd/protein.pdb': 'f2964a785c922502dc86fb4e2e5295d32d41d5b68b8c3246e989de5234c3fd0f',
        'hif2a/01_protein/crd/protein.pdb': '5bbf520e7c102a65cc7ba0253fd66f43562f77284c82b3b9613e997b7ac76c93',

    },
)


@pytest.fixture(params=['p38', 'mcl1', 'cdk2', 'shp2', 'pde2', 'cmet', 'ptp1b',
                        'thrombin', 'cdk8', 'pfkfb3', 'tyk2', 'syk', 'tnks2',
                        'eg5', 'hif2a'])
def plb_proteins(request):
    return PDBS.fetch('{}/01_protein/crd/protein.pdb'.format(request.param))


@pytest.fixture
def cdk2_pdb():
    return PDBS.fetch('cdk2/01_protein/crd/protein.pdb')
