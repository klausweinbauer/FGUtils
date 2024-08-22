from fgutils.parse import parse
from fgutils.rdkit import graph_to_smiles


def test_simple_graph():
    exp_smiles = "CCO"
    g = parse(exp_smiles)
    smiles = graph_to_smiles(g)
    assert exp_smiles == smiles


def test_with_Si():
    g = parse("CSi(C)(C)C")
    smiles = graph_to_smiles(g)
    assert "C[Si](C)(C)C" == smiles


def test_aromaticity():
    g = parse("c1ccccc1")
    smiles = graph_to_smiles(g)
    assert "c1ccccc1" == smiles
