import pytest

from fgutils.parse import parse
from fgutils.rdkit import graph_to_smiles, smiles_to_graph


def test_simple_graph():
    exp_smiles = "[CH3][CH2][OH]"
    g = parse("CCO")
    smiles = graph_to_smiles(g, implicit_h=True)
    assert exp_smiles == smiles


def test_with_Si():
    g = parse("CSi(C)(C)C")
    smiles = graph_to_smiles(g, implicit_h=True)
    assert "[CH3][Si]([CH3])([CH3])[CH3]" == smiles


def test_aromaticity():
    g = parse("c1ccccc1")
    smiles = graph_to_smiles(g, implicit_h=True)
    assert "[cH]1[cH][cH][cH][cH][cH]1" == smiles


def test_aromaticity2():
    input_smiles = "c1cc[nH]c1"
    g = smiles_to_graph(input_smiles)
    out_smiles = graph_to_smiles(g, implicit_h=True)
    assert "[cH]1[cH][cH][nH][cH]1" == out_smiles


def test_parse_invalid():
    with pytest.raises(ValueError):
        smiles_to_graph("CP(=O)(=O)C")


def test_remove_hydrogen():
    input_smiles = "[CH3][CH2][OH]"
    g = smiles_to_graph(input_smiles, implicit_h=False)
    out_smiles = graph_to_smiles(g, implicit_h=False)
    assert "CCO" == out_smiles


@pytest.mark.parametrize(
    "smiles", [("CCO"), ("c1cc[nH]c1"), ("c1ccccc1"), ("c1ccncc1"), ("c1c[nH]cn1")]
)
def test_ignore_implicit_hydrogen(smiles):
    g = smiles_to_graph(smiles, implicit_h=False)
    result_smiles = graph_to_smiles(g, implicit_h=False)
    assert smiles == result_smiles
