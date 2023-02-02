# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
from importlib import resources
import pdbinf
import pooch
import pytest


@pytest.fixture
def pdb_181l():
    file_path = pooch.retrieve(
        url="https://github.com/OpenFreeEnergy/gufe/raw/6f59135b72dc9620b7dbc2884556309dea8498bb/gufe/tests/data/181l_openmmClean.pdb",
        known_hash="4e7d857a1c55891286c8350ce9ffa18e1428dab971a70a24a745cd20de1c5eb3",
    )
    yield file_path


def test_read_181l(pdb_181l):
    m = pdbinf.load_pdb_file(pdb_181l, templates=[pdbinf.STANDARD_AA_DOC])

    assert m
