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


def test_carbonyl_init():
    fg_parent = get_FG_by_name("carbonyl")
    for fg_child in fg_parent.subgroups:
        assert len(fg_parent.pattern.nodes) <= len(fg_child.pattern.nodes)
        assert True == map_full(fg_child.pattern, fg_parent.pattern)


@pytest.mark.parametrize(
    "fg_group,exp_pattern_len",
    [("carbonyl", 2), ("aldehyde", 3), ("ketone", 2), ("carboxylic_acid", 4)],
)
def test_pattern_len(fg_group, exp_pattern_len):
    fg = get_FG_by_name(fg_group)
    assert exp_pattern_len == fg.pattern_len


def test_subgroup_order():
    fg_parent = get_FG_by_name("carbonyl")
    exp_names = ["carboxylic_acid", "amide", "carboxylate_ester", "aldehyde", "ketone"]
    fg_subgroup_names = [fg.name for fg in fg_parent.subgroups]
    for exp_name, name in zip(exp_names, fg_subgroup_names):
        assert (
            exp_name == name
        ), "Subgroup order missmatch. Expected {} but got {}.".format(exp_name, name)
