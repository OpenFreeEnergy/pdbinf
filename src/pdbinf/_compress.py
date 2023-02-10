# This code is part of pdbinf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/pdbinf
"""
Utilities for slimming down cif files to only fields required for templates

Useful for hacking down the chemical components dictionary

e.g.::

import pdbinf

with open('components.cif', 'r') as instream:
    with open('slim_components.cif', 'w') as outstream:
        for line in slimdown_cif(instream):
            outstream.write(line)

"""
from typing import Generator, Iterable


REQUIRED_ATOM_FIELDS = [
    'atom_id', 'pdbx_aromatic_flag',
]
REQUIRED_BOND_FIELDS = [
    'atom_id_1', 'atom_id_2', 'value_order', 'pdbx_aromatic_flag',
]


def strip_comments(thing: Iterable[str]) -> Generator[str, None, None]:
    for line in thing:
        yield line.partition('#')[0].strip()


def slimdown_cif(inp: Iterable[str]) -> Generator[str, None, None]:
    """Take cif contents and return only what is required for templates

    Parameters
    ----------
    inp : iterable of strings
       a stream of the input cif

    Returns
    -------
    out : generator of strings
       line by line, the file to be written
    """
    sup = strip_comments(inp)

    for line in sup:
        if line.startswith('data'):
            yield line + '\n'
            # yield '#\n'
        elif line.startswith('_chem_comp.id'):
            a, b = line.split()
            yield a + ' ' + b + '\n'
            # yield '#\n'
        elif line.startswith('loop_'):
            line = next(sup)
            if line.startswith('_chem_comp_atom'):
                yield 'loop_\n'
                for t in REQUIRED_ATOM_FIELDS:
                    yield f'_chem_comp_atom.{t}\n'
                # grab atom columns to keep
                i = 0
                reqd_atom = {t: None for t in REQUIRED_ATOM_FIELDS}
                while line.startswith('_chem_comp_atom'):
                    cat = line.split('.')[1].strip()
                    if cat in REQUIRED_ATOM_FIELDS:
                        reqd_atom[cat] = i
                    i += 1
                    line = next(sup)
                # spit out only required fields for the following loop
                while line:
                    tokens = line.split()
                    pieces = [tokens[reqd_atom[cat]]
                              for cat in REQUIRED_ATOM_FIELDS]
                    yield ' '.join(pieces) + '\n'
                    line = next(sup)
                yield '#\n'
            elif line.startswith('_chem_comp_bond'):
                yield 'loop_\n'
                for t in REQUIRED_BOND_FIELDS:
                    yield f'_chem_comp_bond.{t}\n'
                i = 0

                reqd_bond = {t: None for t in REQUIRED_BOND_FIELDS}
                while line.startswith('_chem_comp_bond'):
                    cat = line.split('.')[1].strip()
                    if cat in REQUIRED_BOND_FIELDS:
                        reqd_bond[cat] = i
                    i += 1
                    line = next(sup)

                while line:
                    tokens = line.split()
                    pieces = [tokens[reqd_bond[cat]] for cat in
                              REQUIRED_BOND_FIELDS]
                    yield ' '.join(pieces) + '\n'
                    line = next(sup)
                yield '##\n'
