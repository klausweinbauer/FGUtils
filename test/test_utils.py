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


def test_sulfur_ring():
    graph = parse("C:1N:C:S:C:1")
    graph = add_implicit_hydrogens(graph)
    assert 8 == len(graph)
    _assert_Hs(graph, 0, 1)
    _assert_Hs(graph, 1, 0)
    _assert_Hs(graph, 2, 1)
    _assert_Hs(graph, 3, 0)
    _assert_Hs(graph, 4, 1)


def test_nitrogen_5ring():
    graph = parse("C:1C:N(H):C:C:1")
    graph = add_implicit_hydrogens(graph)
    assert 10 == len(graph)
    _assert_Hs(graph, 0, 1)
    _assert_Hs(graph, 1, 1)
    _assert_Hs(graph, 2, 1)
    _assert_Hs(graph, 3, 0)
    _assert_Hs(graph, 4, 1)
    _assert_Hs(graph, 5, 1)


def test_nitrogen_6ring():
    graph = parse("C:1C:C:N:C:C:1")
    graph = add_implicit_hydrogens(graph)
    assert 11 == len(graph)
    _assert_Hs(graph, 0, 1)
    _assert_Hs(graph, 1, 1)
    _assert_Hs(graph, 2, 1)
    _assert_Hs(graph, 3, 0)
    _assert_Hs(graph, 4, 1)
    _assert_Hs(graph, 5, 1)


def test_boric_acid():
    graph = parse("OB(O)O")
    graph = add_implicit_hydrogens(graph)
    assert 7 == len(graph)
    _assert_Hs(graph, 0, 1)
    _assert_Hs(graph, 1, 0)
    _assert_Hs(graph, 2, 1)
    _assert_Hs(graph, 3, 1)


def test_selenium_dioxide():
    graph = parse("O=Se=O")
    graph = add_implicit_hydrogens(graph)
    assert 3 == len(graph)
    _assert_Hs(graph, 0, 0)
    _assert_Hs(graph, 1, 0)
    _assert_Hs(graph, 2, 0)


def test_tin_tetrachloride():
    graph = parse("ClSn(Cl)(Cl)Cl")
    graph = add_implicit_hydrogens(graph)
    assert 5 == len(graph)
    _assert_Hs(graph, 0, 0)
    _assert_Hs(graph, 1, 0)
    _assert_Hs(graph, 2, 0)
    _assert_Hs(graph, 3, 0)
    _assert_Hs(graph, 4, 0)


# SeO2 SnCL4
