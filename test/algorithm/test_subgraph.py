from fgutils.permutation import PermutationMapper
from fgutils.parse import parse
from fgutils.algorithm.subgraph import map_anchored_subgraph, map_subgraph
from fgutils.rdkit import smiles_to_graph
from fgutils.chem.its import get_its

default_mapper = PermutationMapper(wildcard="R", ignore_case=True)


def _assert_anchored_mapping(mapping, valid, exp_mapping=[]):
    assert mapping[0] == valid
    for emap in exp_mapping:
        assert emap in mapping[1]


def _assert_mapping(mapping, valid, exp_mapping=[], index=0):
    assert mapping[index][0] == valid
    for emap in exp_mapping:
        assert emap in mapping[index][1]


def test_simple_match():
    exp_mapping = [(1, 0), (2, 1)]
    g = parse("CCO")
    p = parse("RO")
    m = map_anchored_subgraph(g, 2, p, 1, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_branched_match():
    exp_mapping = [(0, 0), (1, 1), (2, 2), (3, 3)]
    g = parse("CC(=O)O")
    p = parse("RC(=O)O")
    m = map_anchored_subgraph(g, 2, p, 2, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_ring_match():
    exp_mapping = [(0, 2), (1, 1), (2, 0)]
    g = parse("C1CO1")
    p = parse("R1CC1")
    m = map_anchored_subgraph(g, 1, p, 1, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_not_match():
    g = parse("CC=O")
    p = parse("RC(=O)NR")
    m = map_anchored_subgraph(g, 2, p, 2, mapper=default_mapper)
    _assert_anchored_mapping(m, False)


def test_1():
    g = parse("CC=O")
    p = parse("RC(=O)R")
    m = map_anchored_subgraph(g, 0, p, 3, mapper=default_mapper)
    _assert_anchored_mapping(m, False)


def test_2():
    exp_mapping = [(0, 0), (1, 1), (2, 2)]
    g = parse("CC=O")
    p = parse("RC=O")
    m = map_anchored_subgraph(g, 2, p, 2, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_ignore_aromaticity():
    exp_mapping = [(1, 0), (2, 1)]
    g = parse("c1c(=O)cccc1")
    p = parse("C=O")
    m = map_anchored_subgraph(g, 2, p, 1, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_3():
    exp_mapping = [(0, 4), (1, 3), (2, 1), (4, 2), (3, 0)]
    g = parse("COC(C)=O")
    p = parse("RC(=O)OR")
    m = map_anchored_subgraph(g, 4, p, 2, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_explore_wrong_branch():
    exp_mapping = [(0, 2), (1, 1), (2, 0), (3, 3)]
    g = parse("COCO")
    p = parse("C(OR)O")
    m = map_anchored_subgraph(g, 1, p, 1, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_match_subgraph_to_mol():
    exp_mapping = [(0, 2), (1, 0), (2, 1)]
    g = parse("NC(=O)C")
    p = parse("C(=O)N")
    m = map_anchored_subgraph(g, 2, p, 1, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_match_hydrogen():
    # H must be explicit
    g = parse("C=O")
    p = parse("C(H)=O")
    m = map_anchored_subgraph(g, 1, p, 2, mapper=default_mapper)
    _assert_anchored_mapping(m, False)


def test_match_implicit_hydrogen():
    exp_mapping = [(0, 0), (1, 2)]
    g = parse("C=O")
    p = parse("C(H)=O")
    mapper = PermutationMapper(can_map_to_nothing=["H"])
    m = map_anchored_subgraph(g, 1, p, 2, mapper=mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_invalid_bond_match():
    g = parse("C=O")
    p = parse("CO")
    m = map_anchored_subgraph(g, 0, p, 0, mapper=default_mapper)
    _assert_anchored_mapping(m, False)


def test_match_not_entire_subgraph():
    g = parse("C=O")
    p = parse("C(=O)C")
    m = map_anchored_subgraph(g, 0, p, 0, mapper=default_mapper)
    _assert_anchored_mapping(m, False)


def test_start_with_match_to_nothing():
    g = parse("CCO")
    p = parse("HO")
    m = map_anchored_subgraph(g, 2, p, 0, mapper=default_mapper)
    _assert_anchored_mapping(m, False)


def test_match_explicit_hydrogen():
    exp_mapping = [(2, 1), (3, 0)]
    g = parse("CCOH")
    p = parse("HO")
    m = map_anchored_subgraph(g, 2, p, 1, mapper=default_mapper)
    _assert_anchored_mapping(m, True, exp_mapping)


def test_map_subgraph_with_anchor():
    exp_mapping = [(2, 1), (1, 0)]
    g = parse("CCO")
    p = parse("CO")
    m = map_subgraph(g, 2, p, subgraph_anchor=1, mapper=default_mapper)
    _assert_mapping(m, True, exp_mapping)


def test_map_subgraph_without_anchor():
    exp_mapping = [(2, 1), (1, 0)]
    g = parse("CCO")
    p = parse("CO")
    m = map_subgraph(g, 2, p, mapper=default_mapper)
    _assert_mapping(m, True, exp_mapping, index=1)


def test_map_empty_subgraph():
    exp_mapping = []
    g = parse("CCO")
    p = parse("")
    m = map_subgraph(g, 2, p, mapper=default_mapper)
    _assert_mapping(m, True, exp_mapping)


def test_map_invalid_subgraph():
    g = parse("CCO")
    p = parse("Cl")
    m = map_subgraph(g, 2, p, mapper=default_mapper)
    _assert_mapping(m, False)


def test_map_specific_subgraph_to_general_graph():
    g = parse("R")
    p = parse("C")
    m = map_subgraph(g, 0, p, mapper=default_mapper)
    _assert_mapping(m, False)


def test_map_with_bond_tuple():
    smiles = "[C:1](=[O:2])=[O:3]>>[C:1](=[O:2])[O:3]"
    pattern = "R(=O)<2,1>O"
    exp_mapping = [(1, 0), (2, 1), (3, 2)]
    g, h = smiles_to_graph(smiles)
    its = get_its(g, h)
    p = parse(pattern)
    m = map_subgraph(its, 1, p, mapper=default_mapper)
    _assert_mapping(m, True, exp_mapping)


def test_map_with_its_with_wildcard():
    smiles = "[C:1][O:2].[O:3]>>[C:1][O:3].[O:2]"
    pattern = "R(<0,1>R)<1,0>R"
    exp_mapping = [(1, 0), (2, 2), (3, 1)]
    g, h = smiles_to_graph(smiles)
    its = get_its(g, h)
    p = parse(pattern)
    m = map_subgraph(its, 1, p, mapper=default_mapper)
    _assert_mapping(m, True, exp_mapping)


def test_doc_example_2():
    # example for functional_groups.rst:Get changing groups in reaction
    smiles = "[C:1][C:2](=[O:3])[O:4][C:5].[O:6]>>[C:1][C:2](=[O:3])[O:6].[O:4][C:5]"
    pattern = "C(=O)(<0,1>R)<1,0>R"
    exp_mapping = [(2, 0), (3, 1), (4, 3), (6, 2)]
    g, h = smiles_to_graph(smiles)
    its = get_its(g, h)
    p = parse(pattern)
    m = map_subgraph(its, 2, p, mapper=default_mapper)
    _assert_mapping(m, True, exp_mapping)
