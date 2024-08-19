import pytest

from fgutils.parse import parse, tokenize
from fgutils.utils import print_graph


def _assert_graph(g, exp_nodes, exp_edges):
    assert len(exp_nodes) == g.number_of_nodes()
    assert len(exp_edges) == g.number_of_edges()
    for i, sym in exp_nodes.items():
        assert sym == g.nodes[i]["symbol"]
    for i1, i2, order in exp_edges:
        assert order == g.edges[i1, i2]["bond"]


def _ct(token, exp_type, exp_value, exp_col):
    return token[0] == exp_type and token[1] == exp_value and token[2] == exp_col


def test_tokenize():
    it = tokenize("RC(=O)OR")
    assert True is _ct(next(it), "WILDCARD", "R", 0)
    assert True is _ct(next(it), "ATOM", "C", 1)
    assert True is _ct(next(it), "BRANCH_START", "(", 2)
    assert True is _ct(next(it), "BOND", "=", 3)
    assert True is _ct(next(it), "ATOM", "O", 4)
    assert True is _ct(next(it), "BRANCH_END", ")", 5)
    assert True is _ct(next(it), "ATOM", "O", 6)
    assert True is _ct(next(it), "WILDCARD", "R", 7)


def test_tokenize_multichar():
    it = tokenize("RClR")
    assert True is _ct(next(it), "WILDCARD", "R", 0)
    assert True is _ct(next(it), "ATOM", "Cl", 1)
    assert True is _ct(next(it), "WILDCARD", "R", 3)


@pytest.mark.parametrize(
    "pattern,exp_value", (("C<1,2>C", "<1,2>"), ("C<,1>C", "<,1>"), ("C<1,>C", "<1,>"))
)
def test_tokenize_rcbonds(pattern, exp_value):
    it = tokenize(pattern)
    assert True is _ct(next(it), "ATOM", "C", 0)
    assert True is _ct(next(it), "RC_BOND", exp_value, 1)
    assert True is _ct(next(it), "ATOM", "C", len(pattern) - 1)


@pytest.mark.parametrize(
    "pattern,exp_value,exp_col",
    (("C{group}C", "{group}", 1), ("CR{pattern_1}C", "{pattern_1}", 2)),
)
def test_tokenize_node_labels(pattern, exp_value, exp_col):
    it = tokenize(pattern)
    for _ in range(exp_col):
        next(it)
    assert True is _ct(next(it), "NODE_LABEL", exp_value, exp_col)


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


def test_parse_disconnected_graphs():
    exp_nodes = {0: "C", 1: "O"}
    exp_edges = []
    g = parse("C.O")
    _assert_graph(g, exp_nodes, exp_edges)


def test_parse_disconnected_in_ring():
    exp_nodes = {0: "C", 1: "C", 2: "C"}
    exp_edges = [(0, 1, 1), (1, 2, 1)]
    g = parse("C1CC.1")
    _assert_graph(g, exp_nodes, exp_edges)


def test_syntax_error():
    with pytest.raises(SyntaxError):
        parse("X")


def test_syntax_error_invalid_ring_start():
    with pytest.raises(SyntaxError):
        parse("1CCC1")


def test_parse_explicit_hydrogen():
    exp_nodes = {0: "H", 1: "O"}
    exp_edges = [(0, 1, 1)]
    g = parse("HO")
    _assert_graph(g, exp_nodes, exp_edges)


def test_parse_its():
    exp_nodes = {i: "C" for i in range(6)}
    exp_edges = [
        (0, 1, (2, 1)),
        (0, 5, (0, 1)),
        (1, 2, (1, 2)),
        (2, 3, (2, 1)),
        (3, 4, (0, 1)),
        (4, 5, (2, 1)),
    ]
    g = parse("C1<2,>C<,2>C<2,>C<0,1>C<2,>C<0,1>1")
    _assert_graph(g, exp_nodes, exp_edges)


def test_parse_labled_graph():
    exp_nodes = {0: "C", 1: "#"}
    exp_edges = [(0, 1, 1)]
    g = parse("C{group}")
    print_graph(g)
    _assert_graph(g, exp_nodes, exp_edges)
    assert not g.nodes(data=True)[0]["is_labeled"]
    assert g.nodes(data=True)[1]["is_labeled"]


def test_parse_multi_labled_graph():
    exp_nodes = {0: "#"}
    exp_edges = []
    g = parse("{group1,group2}")
    _assert_graph(g, exp_nodes, exp_edges)
    assert g.nodes(data=True)[0]["is_labeled"]
    assert g.nodes(data=True)[0]["labels"] == ["group1", "group2"]


def test_bond_tpye_in_double_closing():
    exp_nodes = {i: "C" for i in range(5)}
    exp_edges = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 2), (3, 4, 1), (4, 2, 1)]
    g = parse("C1CC2C=1C2")
    _assert_graph(g, exp_nodes, exp_edges)
