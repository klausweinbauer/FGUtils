from fgutils.parse import parse
from fgutils.utils import add_implicit_hydrogens


def _assert_Hs(graph, idx, h_cnt):
    atom_sym = graph.nodes[idx]["symbol"]
    h_neighbors = [
        n_id for n_id in graph.neighbors(idx) if graph.nodes[n_id]["symbol"] == "H"
    ]
    assert h_cnt == len(
        h_neighbors
    ), "Expected atom {} to have {} hydrogens but found {} instead.".format(
        atom_sym, h_cnt, len(h_neighbors)
    )


def test_add_implicit_hydrogens_1():
    graph = parse("C=O")
    graph = add_implicit_hydrogens(graph)
    assert 4 == len(graph)
    _assert_Hs(graph, 0, 2)
    _assert_Hs(graph, 1, 0)


def test_add_implicit_hydrogens_2():
    graph = parse("CO")
    graph = add_implicit_hydrogens(graph)
    assert 6 == len(graph)
    _assert_Hs(graph, 0, 3)
    _assert_Hs(graph, 1, 1)


def test_add_implicit_hydrogens_3():
    graph = parse("HC(H)(H)OH")
    graph = add_implicit_hydrogens(graph)
    assert 6 == len(graph)
    _assert_Hs(graph, 1, 3)
    _assert_Hs(graph, 4, 1)
