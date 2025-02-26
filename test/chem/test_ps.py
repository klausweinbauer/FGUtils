import pytest

from fgutils.chem.ps import (
    atomic_sym2num,
    get_electronegativity,
    get_num_valence_electrons,
    find_atom_numbers,
)


def test_atom2symbol_dict():
    tmpk = []
    tmpv = []
    for k, v in atomic_sym2num.items():
        assert k not in tmpk, "Found symbol '{}' twice.".format(k)
        assert v not in tmpv, "Found number '{}' twice.".format(v)
        tmpk.append(k)
        tmpv.append(v)


@pytest.mark.parametrize("atom,exp_value", [(6, 2.55), ("O", 3.44)])
def test_electonegativity(atom, exp_value):
    value = get_electronegativity(atom)
    assert exp_value == value


@pytest.mark.parametrize("atom,exp_value", [(6, 4), ("O", 6)])
def test_valence_electrons(atom, exp_value):
    value = get_num_valence_electrons(atom)
    assert exp_value == value


@pytest.mark.parametrize(
    "exp_value, params",
    [([8], (6, 3.44)), ([21, 93], (None, 1.36)), ([9, 17, 35, 53, 85, 117], (7, None))],
)
def test_find_atom_numbers(exp_value, params):
    value = find_atom_numbers(*params)
    assert all(v in value for v in exp_value)
