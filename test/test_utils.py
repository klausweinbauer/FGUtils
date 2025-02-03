import numpy as np
from numpy.testing import assert_array_equal

from fgutils.rdkit import mol_smiles_to_graph, graph_to_smiles
from fgutils.parse import parse
from fgutils.utils import add_implicit_hydrogens, complete_aam, mol_compare
from fgutils.const import SYMBOL_KEY


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


def test_aam_complete():
    in_smiles = "[C:1][C:3]CO"
    exp_smiles = "[CH3:1][CH2:3][CH2:2][OH:4]"
    g = mol_smiles_to_graph(in_smiles)
    complete_aam(g)
    out_smiles = graph_to_smiles(g, canonical=False)
    assert out_smiles == exp_smiles


def test_aam_complete_with_offset():
    in_smiles = "[C:1][C:3]CO"
    exp_smiles = "[CH3:1][CH2:3][CH2:10][OH:11]"
    g = mol_smiles_to_graph(in_smiles)
    complete_aam(g, offset=10)
    out_smiles = graph_to_smiles(g, canonical=False)
    assert out_smiles == exp_smiles


def test_aam_complete_with_offset_min():
    in_smiles = "[C:5][C:9]CO"
    exp_smiles = "[CH3:5][CH2:9][CH2:6][OH:7]"
    g = mol_smiles_to_graph(in_smiles)
    complete_aam(g, offset="min")
    out_smiles = graph_to_smiles(g, canonical=False)
    assert out_smiles == exp_smiles


def test_aam_complete_empty_mapping_with_offset_min():
    in_smiles = "CO"
    exp_smiles = "[CH3:1][OH:2]"
    g = mol_smiles_to_graph(in_smiles)
    complete_aam(g, offset="min")
    out_smiles = graph_to_smiles(g, canonical=False)
    assert out_smiles == exp_smiles


def test_mol_compare():
    target_smiles = "CC(=O)O"
    c1 = "O=C(C)O"
    c2 = "CC(=O)OC"
    c3 = "CC(=O)O"
    c4 = "[CH3:1][C:2](=[O:3])[O:4]"
    exp_result = np.array([1, 0, 1, 1])
    target = mol_smiles_to_graph(target_smiles)
    candidates = [mol_smiles_to_graph(c) for c in [c1, c2, c3, c4]]
    output = mol_compare(candidates, target)
    assert_array_equal(output, exp_result)


def test_mol_compare_disconnected():
    target_smiles = "CC.O"
    c1 = "O.CC"
    c2 = "[OH2].[CH3][CH3]"
    exp_result = np.array([1, 1])
    target = mol_smiles_to_graph(target_smiles)
    candidates = [mol_smiles_to_graph(c) for c in [c1, c2]]
    output = mol_compare(candidates, target)
    assert_array_equal(output, exp_result)
