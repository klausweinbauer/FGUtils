import pytest

from fgutils.const import BOND_KEY
from fgutils.parse import parse
from fgutils.chem.valence import _check_mol_valence, _check_its_valence, check_valence
from fgutils.its import ITS


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C", True),
        ("HC(H)(H)H", True),
        ("HC(H)(H)(H)H", False),
        ("N", True),
        ("HN(H)H", True),
        ("HN(H)(H)H", False),
        ("c1cccc2c1cccc2", True)
    ],
)
def test_mol_valence_check(smiles, expected):
    g = parse(smiles)
    result = _check_mol_valence(g)
    assert expected == result


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C", False),
        ("HC(H)(H)H", True),
        ("HC(H)(H)(H)H", False),
        ("N", False),
        ("HN(H)H", True),
        ("HN(H)(H)H", False),
    ],
)
def test_mol_valence_check_exact(smiles, expected):
    g = parse(smiles)
    result = _check_mol_valence(g, exact=True)
    assert expected == result


def test_mol_valence_check_with_charge():
    g = parse("HN(H)(H)H")
    print(g.nodes(data=True))
    g.add_edge(1, 1, **{BOND_KEY: -0.5})
    result = _check_mol_valence(g)
    assert result is True


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C<2,1>C", True),
        ("C<1,2>C(<0,1>H)<1,2>C", False),
        ("C<2,1>C(<1,0>H)<2,1>C", False),
    ],
)
def test_its_valence_check(smiles, expected):
    its = parse(smiles)
    result = _check_its_valence(its)
    assert expected == result


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C<2,1>C", False),
        ("HC(H)<2,1>C(H)H", False),
        ("HC(H)<1,2>C(H)H", False),
        ("HC(H)<2,2>C(H)H", True),
    ],
)
def test_its_valence_check_exact(smiles, expected):
    its = parse(smiles)
    result = _check_its_valence(its, exact=True)
    assert expected == result


@pytest.mark.parametrize(
    "smiles,charge,expected",
    [
        ("H<1,1>N(H)H", [0, 0], True),
        ("H<1,1>N(H)(H)H", [0, 0], False),
        ("H<1,1>N(H)(H)H", [0, -0.5], False),
        ("H<1,1>N(H)(H)H", [-0.5, 0], False),
        ("H<1,1>N(H)(H)H", [-0.5, -0.5], True),
    ],
)
def test_its_valence_check_with_charge(smiles, charge, expected):
    its = parse(smiles)
    its.add_edge(1, 1, **{BOND_KEY: charge})
    result = _check_its_valence(its)
    assert expected == result


@pytest.mark.parametrize("smiles", [("C<1,1>C"), ("C")])
def test_check_valence(smiles):
    its = parse(smiles)
    result = check_valence(its)
    assert result is True


def test_check_valence_with_its():
    its = ITS(parse("C<2,1>C"))
    result = check_valence(its)
    assert result is True
