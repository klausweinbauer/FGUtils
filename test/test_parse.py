import pytest

from fgutils.parse import *


def _assert_graph(g, exp_nodes, exp_edges):
    assert len(exp_nodes) == g.number_of_nodes()
    assert len(exp_edges) == g.number_of_edges()
    for i, sym in exp_nodes.items():
        assert sym == g.nodes[i]["symbol"]
    for i1, i2, order in exp_edges:
        assert order == g.edges[i1, i2]["bond"]

def test_tokenize():
    def _ct(token, exp_type, exp_value, exp_col):
        return (
            token[0] == exp_type and token[1] == exp_value and token[2] == exp_col
        )

    it = tokenize("RC(=O)OR")
    assert True == _ct(next(it), "WILDCARD", "R", 0)
    assert True == _ct(next(it), "ATOM", "C", 1)
    assert True == _ct(next(it), "BRANCH_START", "(", 2)
    assert True == _ct(next(it), "BOND", "=", 3)
    assert True == _ct(next(it), "ATOM", "O", 4)
    assert True == _ct(next(it), "BRANCH_END", ")", 5)
    assert True == _ct(next(it), "ATOM", "O", 6)
    assert True == _ct(next(it), "WILDCARD", "R", 7)

def test_branch():
    exp_nodes = {0: "R", 1: "C", 2: "O", 3: "O", 4: "R"}
    exp_edges = [(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)]
    g = parse("RC(=O)OR")
    _assert_graph(g, exp_nodes, exp_edges)

def test_multi_branch():
    exp_nodes = {0: "C", 1: "C", 2: "C", 3: "O", 4: "O", 5: "C"}
    exp_edges = [(0, 1, 1), (1, 2, 2), (2, 3, 1), (2, 4, 1), (1, 5, 1)]
    g = parse("CC(=C(O)O)C")
    _assert_graph(g, exp_nodes, exp_edges)

def test_ring_3():
    exp_nodes = {0: "C", 1: "C", 2: "C"}
    exp_edges = [(0, 1, 1), (1, 2, 1), (0, 2, 1)]
    g = parse("C1CC1")
    _assert_graph(g, exp_nodes, exp_edges)

def test_ring_4():
    exp_nodes = {0: "C", 1: "C", 2: "C", 3: "C"}
    exp_edges = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (0, 3, 1)]
    g = parse("C1CCC1")
    _assert_graph(g, exp_nodes, exp_edges)

def test_multi_ring():
    exp_nodes = {0: "C", 1: "C", 2: "C", 3: "C"}
    exp_edges = [(0, 1, 1), (1, 2, 1), (0, 2, 1), (1, 3, 1), (2, 3, 1)]
    g = parse("C1C2C1C2")
    _assert_graph(g, exp_nodes, exp_edges)

def test_aromatic_ring():
    exp_nodes = {i: "c" for i in range(6)}
    exp_edges = [(0, 5, 1.5), *[(i, i + 1, 1.5) for i in range(5)]]
    g = parse("c1ccccc1")
    _assert_graph(g, exp_nodes, exp_edges)

def test_aromatic_ring_syntax_error():
    with pytest.raises(SyntaxError):
        parse("c1ccccC1")

def test_complex_aromatic_ring():
    exp_nodes = {i: "c" for i in range(9)}
    exp_nodes[0] = "C"
    exp_nodes[3] = "C"
    exp_nodes[5] = "C"
    exp_edges = [
        (0, 1, 1),
        (1, 2, 1.5),
        (2, 3, 1),
        (2, 4, 1.5),
        (4, 5, 2),
        (4, 6, 1.5),
        (6, 7, 1.5),
        (7, 8, 1.5),
        (8, 1, 1.5),
    ]
    g = parse("Cc1c(C)c(=C)ccc1")
    _assert_graph(g, exp_nodes, exp_edges)

def test_symtax_errir():
    with pytest.raises(SyntaxError):
        parse("X")

def test_parse_explicit_hydrogen():
    exp_nodes = {0: "H", 1: "O"}
    exp_edges = [(0, 1, 1)]
    g = parse("HO")
    _assert_graph(g, exp_nodes, exp_edges)
