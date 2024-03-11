import pytest
import networkx as nx

from fgutils.fgconfig import *
from fgutils.mapping import map_pattern


def test_init():
    config = {
        "name": "carbonyl",
        "pattern": "C(=O)",
        "subgroups": [
            {
                "name": "aldehyde",
                "pattern": "RC(=O)",
                "group_atoms": [1, 2],
                "anti_pattern": "RC(=O)R",
            },
            {"name": "ketone", "pattern": "RC(=O)R", "group_atoms": [1, 2]},
        ],
    }

    fgc = FGConfig(**config)
    assert 2 == len(fgc.subgroups)
    sfgc1 = fgc.subgroups[0]
    sfgc2 = fgc.subgroups[1]
    assert "carbonyl" == fgc.name
    assert "ketone" == sfgc1.name
    assert "aldehyde" == sfgc2.name
    assert "carbonyl" == sfgc1.parent.name
    assert "carbonyl" == sfgc2.parent.name
    assert True == isinstance(fgc.group_atoms, list)
    assert True == isinstance(sfgc1.group_atoms, list)
    assert [[1 == 2]], sfgc1.group_atoms
    assert [[1 == 2]], sfgc2.group_atoms
    assert isinstance(fgc.pattern, nx.Graph)
    assert isinstance(sfgc1.pattern, nx.Graph)
    assert isinstance(sfgc2.pattern, nx.Graph)


def map_full(graph, pattern):
    for i in range(len(graph)):
        r, _ = map_pattern(graph, i, pattern)
        if r is True:
            return True
    return False


def test_fg_tree_init():
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


def test_fg_config_uniqueness():
    name_list = []
    pattern_list = []
    for fg in get_FG_list():
        assert fg.name not in name_list
        name_list.append(fg.name)
        assert fg.pattern_str not in pattern_list
        pattern_list.append(fg.pattern_str)
    assert len(functional_group_config) == len(name_list)
    assert len(functional_group_config) == len(pattern_list)


def test_build_FG_tree():
    def _print(fg, indent=0):
        print(
            "{}{:<20}{:<30} {}".format(
                indent * " ",
                fg.name,
                "[Parent: {}]".format(fg.parent.name if fg.parent is not None else "ROOT"),
                fg.pattern_str,
            )
        )
        for sfg in fg.subgroups:
            _print(sfg, indent + 2)

    for fg in build_FG_tree():
        _print(fg)
    assert False
