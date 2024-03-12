import pytest
import networkx as nx

from fgutils.fgconfig import *
from fgutils.mapping import map_pattern


def test_init():
    config = {
        "name": "carbonyl",
        "pattern": "C(=O)",
        "group_atoms": [1, 2],
        "anti_pattern": "RC(=O)R",
    }

    fgc = FGConfig(**config)
    assert "carbonyl" == fgc.name
    assert True == isinstance(fgc.group_atoms, list)
    assert [1, 2] == fgc.group_atoms
    assert isinstance(fgc.pattern, nx.Graph)
    assert 1 == len(fgc.anti_pattern)
    assert 4 == len(fgc.anti_pattern[0])


def map_full(graph, pattern):
    for i in range(len(graph)):
        r, _ = map_pattern(graph, i, pattern)
        if r is True:
            return True
    return False


def _assert_structure(
    config: FGConfig,
    exp_parent: None | FGConfig | list[FGConfig],
    exp_subgroups: None | FGConfig | list[FGConfig] = None,
):
    exp_parent = (
        [exp_parent]
        if isinstance(exp_parent, FGConfig)
        else []
        if exp_parent is None
        else exp_parent
    )
    exp_subgroups = (
        [exp_subgroups]
        if isinstance(exp_subgroups, FGConfig)
        else []
        if exp_subgroups is None
        else exp_subgroups
    )
    assert (
        exp_parent == config.parents
    ), "Expected parent of {} to be {} but got {}.".format(
        config.name,
        [g.name for g in exp_parent] if isinstance(exp_parent, list) else None,
        [g.name for g in config.parents] if isinstance(config.parents, list) else None,
    )
    exp_subgroups = sort_by_pattern_len(exp_subgroups, reverse=True)
    assert (
        exp_subgroups == config.subgroups
    ), "Expected sugroups of {} to be {} but got {}.".format(
        config.name, [g.name for g in exp_subgroups], [g.name for g in config.subgroups]
    )


def test_search_parent():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    parents = search_parents([fg1], fg2)
    assert parents == [fg1]


def test_get_no_parent():
    fg1 = FGConfig(name="1", pattern="RO")
    fg2 = FGConfig(name="2", pattern="RC")
    parents = search_parents([fg1], fg2)
    assert parents == None


def test_get_correct_parent():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RO")
    fg3 = FGConfig(name="3", pattern="ROR")
    parents = search_parents([fg1, fg2], fg3)
    assert parents == [fg2]


def test_get_multiple_parents():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RO")
    fg3 = FGConfig(name="3", pattern="RCO")
    parents = search_parents([fg1, fg2], fg3)
    assert parents == [fg1, fg2]


def test_get_multiple_unique_parents():
    fg11 = FGConfig(name="11", pattern="RC")
    fg12 = FGConfig(name="12", pattern="RO")
    fg2 = FGConfig(name="2", pattern="RCOR")
    fg3 = FGConfig(name="3", pattern="RCOCR")
    fg11.subgroups = [fg2]
    fg12.subgroups = [fg2]
    fg2.parents = [fg11, fg12]
    parents = search_parents([fg11, fg12], fg3)
    assert parents == [fg2]


def test_get_parent_recursive():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    fg3 = FGConfig(name="3", pattern="RCO")
    fg1.subgroups = [fg2]
    fg2.parents = [fg1]
    parents = search_parents([fg1], fg3)
    assert parents == [fg2]


def test_insert_child_between():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg1, fg3, fg2])
    assert tree == [fg1]
    _assert_structure(fg1, None, fg2)
    _assert_structure(fg2, fg1, fg3)
    _assert_structure(fg3, fg2)


def test_insert_child_after():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg1, fg2, fg3])
    assert tree == [fg1]
    _assert_structure(fg1, None, fg2)
    _assert_structure(fg2, fg1, fg3)
    _assert_structure(fg3, fg2)


def test_insert_new_root():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg2, fg3, fg1])
    assert tree == [fg1]
    _assert_structure(fg1, None, fg2)
    _assert_structure(fg2, fg1, fg3)
    _assert_structure(fg3, fg2)


def test_insert_child_in_between_multiple():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCOR")
    fg31 = FGConfig(name="31", pattern="RCOH")
    fg32 = FGConfig(name="32", pattern="RC=O")
    tree = build_config_tree_from_list([fg1, fg31, fg32, fg2])
    assert tree == [fg1]
    _assert_structure(fg1, None, [fg2, fg32])
    _assert_structure(fg2, fg1, fg31)
    _assert_structure(fg31, fg2)
    _assert_structure(fg32, fg1)


def test_multiple_parents():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RO")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg1, fg2, fg3])
    assert tree == [fg1, fg2]
    _assert_structure(fg1, None, fg3)
    _assert_structure(fg2, None, fg3)
    _assert_structure(fg3, [fg1, fg2])


def test_tree_structure():
    def _check_fg(fg):
        for sfg in fg.subgroups:
            print("Test {} -> {}.".format(fg.name, sfg.name))
            assert fg.pattern_len <= sfg.pattern_len
            assert True == map_full(sfg.pattern, fg.pattern)
            assert False == map_full(
                fg.pattern, sfg.pattern
            ), "Parent pattern {} contains child pattern {}.".format(
                fg.pattern_str, sfg.pattern_str
            )
            _check_fg(sfg)

    for root_fg in build_FG_tree():
        _check_fg(root_fg)


@pytest.mark.parametrize(
    "fg_group,exp_pattern_len",
    [("carbonyl", 2), ("aldehyde", 3), ("ketone", 2), ("carboxylic_acid", 4)],
)
def test_pattern_len(fg_group, exp_pattern_len):
    fg = get_FG_by_name(fg_group)
    assert exp_pattern_len == fg.pattern_len


def test_config_name_uniqueness():
    name_list = []
    for fg in get_FG_list():
        assert fg.name not in name_list, "Config name '{}' already exists.".format(
            fg.name
        )
        name_list.append(fg.name)
    assert len(functional_group_config) == len(name_list)


def test_config_pattern_uniqueness():
    pattern_list = []
    for fg in get_FG_list():
        assert (
            fg.pattern_str not in pattern_list
        ), "Config pattern '{}' already exists with name '{}'.".format(
            fg.pattern_str, fg.name
        )
        pattern_list.append(fg.pattern_str)
    assert len(functional_group_config) == len(pattern_list)


def test_build_tree():
    build_FG_tree()
    assert False
