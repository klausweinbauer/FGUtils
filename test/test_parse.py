import pytest
import networkx as nx

from fgutils.parse import parse, tokenize, Parser
from fgutils.utils import print_graph
from fgutils.const import SYMBOL_KEY, BOND_KEY, AAM_KEY, IS_LABELED_KEY, LABELS_KEY


def _assert_graph(g, exp_nodes, exp_edges):
    assert len(exp_nodes) == g.number_of_nodes()
    assert len(exp_edges) == g.number_of_edges()
    for i, sym in exp_nodes.items():
        assert sym == g.nodes[i][SYMBOL_KEY]
    for i1, i2, order in exp_edges:
        assert order == g.edges[i1, i2][BOND_KEY]


def _ct(token, exp_type, exp_value, exp_col):
    return token[0] == exp_type and token[1] == exp_value and token[2] == exp_col


def assert_multi_graph_eq(exp_graph, act_graph, ignore_keys=[AAM_KEY]):
    def _nm(n1, n2):
        for k, v in n1.items():
            if k in ignore_keys:
                continue
            if k not in n2.keys() or n2[k] != v:
                return False
        for k, v in n2.items():
            if k in ignore_keys:
                continue
            if k not in n1.keys() or n1[k] != v:
                return False
        return True

    def _em(e1, e2):
        e1_bonds = sorted([d[BOND_KEY] for j, d in e1.items()])
        e2_bonds = sorted([d[BOND_KEY] for j, d in e2.items()])
        return e1_bonds == e2_bonds

    is_isomorphic = nx.is_isomorphic(
        exp_graph, act_graph, node_match=_nm, edge_match=_em
    )
    assert is_isomorphic, "Graphs are not isomorphic."


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
    assert not g.nodes(data=True)[0][IS_LABELED_KEY]
    assert g.nodes(data=True)[1][IS_LABELED_KEY]


def test_parse_multi_labled_graph():
    exp_nodes = {0: "#"}
    exp_edges = []
    g = parse("{group1,group2}")
    _assert_graph(g, exp_nodes, exp_edges)
    assert g.nodes(data=True)[0][IS_LABELED_KEY]
    assert g.nodes(data=True)[0][LABELS_KEY] == ["group1", "group2"]


def test_bond_tpye_in_double_closing():
    exp_nodes = {i: "C" for i in range(5)}
    exp_edges = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 2), (3, 4, 1), (4, 2, 1)]
    g = parse("C1CC2C=1C2")
    _assert_graph(g, exp_nodes, exp_edges)


def test_parse_multigraph():
    exp_g = nx.MultiGraph()
    exp_g.add_nodes_from([0, 1], **{SYMBOL_KEY: "C"})
    exp_g.add_edge(0, 1, **{BOND_KEY: 1})
    exp_g.add_edge(0, 1, **{BOND_KEY: 2})
    parser = Parser(use_multigraph=True)
    g = parser.parse("C1=C1")
    assert_multi_graph_eq(exp_g, g, ignore_keys=[LABELS_KEY, IS_LABELED_KEY, AAM_KEY])


def test_doc_example1():
    # example for fgutils.parse.Parser()
    parser = Parser()
    g = parser("CC(O)=O")
    assert "Graph with 4 nodes and 3 edges" == str(g)


def test_bond_type():
    parser = Parser()
    g = parser("CCC")
    for u, v, d in g.edges(data=True):
        assert isinstance(d[BOND_KEY], int), "Bond {}-{}:{} is of type {}.".format(
            u, v, d[BOND_KEY], type(d[BOND_KEY])
        )


def test_bond_type_in_its():
    parser = Parser()
    g = parser("CC<2,1>C")
    for u, v, d in g.edges(data=True):
        assert isinstance(d[BOND_KEY], tuple), "Bond {}-{}:{} is of type {}.".format(
            u, v, d[BOND_KEY], type(d[BOND_KEY])
        )
