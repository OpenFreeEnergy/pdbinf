[![Logo](https://img.shields.io/badge/OSMF-OpenFreeEnergy-%23002f4a)](https://openfree.energy/)
[![build](https://github.com/OpenFreeEnergy/pdbinf/actions/workflows/ci.yaml/badge.svg)](https://github.com/OpenFreeEnergy/pdbinf/actions/workflows/ci.yaml)
[![coverage](https://codecov.io/gh/OpenFreeEnergy/pdbinf/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenFreeEnergy/pdbinf)

PDBInf allows you to read PDB files using RDKit and assign bond orders from .cif (PDBx) template files.

```python
# provide a PDBx defining the TPO nonstandard residue
nonstandard_residues = gemmi.cif.read('./tpo.cif')

# load a PDB file using the built-in standard amino acid doc and providing information on the nonstandard residues
pdbinf.load_pdb_file('./tpo.pdb', templates=[pdbinf.STANDARD_AA_DOC, nonstandard_residues])

```
