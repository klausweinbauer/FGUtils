from fgutils.parse import parse
from fgutils.utils import add_implicit_hydrogens, split_its
from fgutils.const import SYMBOL_KEY

from .test_parse import _assert_graph


def _assert_Hs(graph, idx, h_cnt):
    atom_sym = graph.nodes[idx][SYMBOL_KEY]
    h_neighbors = [
        n_id for n_id in graph.neighbors(idx) if graph.nodes[n_id][SYMBOL_KEY] == "H"
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


def test_add_implicit_hydrogens_4():
    graph = parse("C")
    graph = add_implicit_hydrogens(graph)
    assert 5 == len(graph)
    _assert_Hs(graph, 0, 4)


# def test_add_implicit_hydrogens_to_its_1():
#     exp_its = parse("HC1(=O)<1,0>O(<0,1>H<1,0>O<0,1>1)C(H)(H)H")
#     its = parse("C(=O)(<0,1>O)<1,0>OC", init_aam=True)
#     g, h = split_its(its)
#     print(g.nodes(data=True))
#     print("{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h)))
#     its_h = add_implicit_hydrogens(its)
#     g, h = split_its(its_h)
#     assert_graph_eq(exp_its, its_h)
#     assert False


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


def test_split_its():
    its = parse("C1<2,>C<,2>C<2,>C(C)<0,1>C<2,>C(C(=O)O)<0,1>1")
    exp_nodes = {i: "C" for i in range(10)}
    exp_nodes[8] = "O"
    exp_nodes[9] = "O"
    exp_edges_g = [
        (0, 1, 2),
        (1, 2, 1),
        (2, 3, 2),
        (3, 4, 1),
        (5, 6, 2),
        (6, 7, 1),
        (7, 8, 2),
        (7, 9, 1),
    ]
    exp_edges_h = [
        (0, 1, 1),
        (0, 6, 1),
        (1, 2, 2),
        (2, 3, 1),
        (3, 4, 1),
        (3, 5, 1),
        (5, 6, 1),
        (6, 7, 1),
        (7, 8, 2),
        (7, 9, 1),
    ]
    g, h = split_its(its)
    _assert_graph(g, exp_nodes, exp_edges_g)
    _assert_graph(h, exp_nodes, exp_edges_h)
