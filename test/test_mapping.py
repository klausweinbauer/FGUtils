import pytest

from fgutils.parse import parse
from fgutils.mapping import *


def _assert_mapping(mapping, valid, exp_mapping=[]):
    assert mapping[0] == valid
    for emap in exp_mapping:
        assert emap in mapping[1]


def test_simple_match():
    exp_mapping = [(1, 0), (2, 1)]
    g = parse("CCO")
    p = parse("RO")
    m = map_anchored_pattern(g, 2, p, 1)
    _assert_mapping(m, True, exp_mapping)


def test_branched_match():
    exp_mapping = [(0, 0), (1, 1), (2, 2), (3, 3)]
    g = parse("CC(=O)O")
    p = parse("RC(=O)O")
    m = map_anchored_pattern(g, 2, p, 2)
    _assert_mapping(m, True, exp_mapping)


def test_ring_match():
    exp_mapping = [(0, 2), (1, 1), (2, 0)]
    g = parse("C1CO1")
    p = parse("R1CC1")
    m = map_anchored_pattern(g, 1, p, 1)
    _assert_mapping(m, True, exp_mapping)


def test_not_match():
    g = parse("CC=O")
    p = parse("RC(=O)NR")
    m = map_anchored_pattern(g, 2, p, 2)
    _assert_mapping(m, False)


def test_1():
    g = parse("CC=O")
    p = parse("RC(=O)R")
    m = map_anchored_pattern(g, 0, p, 3)
    _assert_mapping(m, False)


def test_2():
    exp_mapping = [(0, 0), (1, 1), (2, 2)]
    g = parse("CC=O")
    p = parse("RC=O")
    m = map_anchored_pattern(g, 2, p, 2)
    _assert_mapping(m, True, exp_mapping)


def test_ignore_aromaticity():
    exp_mapping = [(1, 0), (2, 1)]
    g = parse("c1c(=O)cccc1")
    p = parse("C=O")
    m = map_anchored_pattern(g, 2, p, 1)
    _assert_mapping(m, True, exp_mapping)


def test_3():
    exp_mapping = [(0, 4), (1, 3), (2, 1), (4, 2), (3, 0)]
    g = parse("COC(C)=O")
    p = parse("RC(=O)OR")
    m = map_anchored_pattern(g, 4, p, 2)
    _assert_mapping(m, True, exp_mapping)


def test_explore_wrong_branch():
    exp_mapping = [(0, 2), (1, 1), (2, 0), (3, 3)]
    g = parse("COCO")
    p = parse("C(OR)O")
    m = map_anchored_pattern(g, 1, p, 1)
    _assert_mapping(m, True, exp_mapping)


def test_match_pattern_to_mol():
    exp_mapping = [(0, 2), (1, 0), (2, 1)]
    g = parse("NC(=O)C")
    p = parse("C(=O)N")
    m = map_anchored_pattern(g, 2, p, 1)
    _assert_mapping(m, True, exp_mapping)


def test_match_hydrogen():
    # H must be explicit
    g = parse("C=O")
    p = parse("C(H)=O")
    m = map_anchored_pattern(g, 1, p, 2)
    _assert_mapping(m, False)


def test_invalid_bond_match():
    g = parse("C=O")
    p = parse("CO")
    m = map_anchored_pattern(g, 0, p, 0)
    _assert_mapping(m, False)


def test_match_not_entire_pattern():
    g = parse("C=O")
    p = parse("C(=O)C")
    m = map_anchored_pattern(g, 0, p, 0)
    _assert_mapping(m, False)


def test_start_with_match_to_nothing():
    g = parse("CCO")
    p = parse("HO")
    m = map_anchored_pattern(g, 2, p, 0)
    _assert_mapping(m, False)


def test_match_explicit_hydrogen():
    exp_mapping = [(2, 1), (3, 0)]
    g = parse("CCOH")
    p = parse("HO")
    m = map_anchored_pattern(g, 2, p, 1)
    _assert_mapping(m, True, exp_mapping)


def test_map_pattern_with_anchor():
    exp_mapping = [(2, 1), (1, 0)]
    g = parse("CCO")
    p = parse("CO")
    m = map_pattern(g, 2, p, pattern_anchor=1)
    _assert_mapping(m, True, exp_mapping)


def test_map_pattern_without_anchor():
    exp_mapping = [(2, 1), (1, 0)]
    g = parse("CCO")
    p = parse("CO")
    m = map_pattern(g, 2, p)
    _assert_mapping(m, True, exp_mapping)


def test_map_empty_pattern():
    exp_mapping = []
    g = parse("CCO")
    p = parse("")
    m = map_pattern(g, 2, p)
    _assert_mapping(m, True, exp_mapping)


def test_map_invalid_pattern():
    g = parse("CCO")
    p = parse("Cl")
    m = map_pattern(g, 2, p)
    _assert_mapping(m, False)


def test_map_specific_pattern_to_general_graph():
    g = parse("R")
    p = parse("C")
    m = map_pattern(g, 0, p)
    _assert_mapping(m, False)
