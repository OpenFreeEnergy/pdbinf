# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
import pdbinf
import pooch
import pytest

app = pytest.importorskip('openmm.app')

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


def test_openmm_regression(plb_proteins, tpo_doc, request):
    ref = app.PDBFile(plb_proteins)

    if 'cdk2' in plb_proteins:
        t = [pdbinf.STANDARD_AA_DOC, tpo_doc]
    else:
        t = [pdbinf.STANDARD_AA_DOC]

    m = pdbinf.load_pdb_file(plb_proteins, templates=t)

    for b in ref.topology.bonds():
        ix1, ix2 = b.atom1.index, b.atom2.index

        x = m.GetBondBetweenAtoms(ix1, ix2)

        assert x is not None
