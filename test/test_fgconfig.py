import pytest
import networkx as nx

from fgutils.algorithm.subgraph import map_subgraph_to_graph
from fgutils.permutation import PermutationMapper
from fgutils.fgconfig import (
    FGConfigProvider,
    FGConfig,
    FGTreeNode,
    search_parents,
    is_subgroup,
    build_config_tree_from_list,
    _default_fg_config,
)
from fgutils.const import SYMBOL_KEY


def test_init():
    config = {
        "name": "carbonyl",
        "pattern": "C(=O)",
        "group_atoms": [1, 2],
        "anti_pattern": "RC(=O)R",
    }

    fgc = FGConfig(**config)
    assert "carbonyl" == fgc.name
    assert True is isinstance(fgc.group_atoms, list)
    assert [1, 2] == fgc.group_atoms
    assert isinstance(fgc.pattern, nx.Graph)
    assert 1 == len(fgc.anti_pattern)
    assert 4 == len(fgc.anti_pattern[0])


def _init_fgnode(name, pattern) -> FGTreeNode:
    return FGTreeNode(FGConfig(name=name, pattern=pattern))


default_mapper = PermutationMapper(wildcard="R", ignore_case=True)


def test_search_parent():
    fg1 = _init_fgnode("1", "RC")
    fg2 = _init_fgnode("2", "RCR")
    parents = search_parents([fg1], fg2, mapper=default_mapper)
    assert parents == [fg1]


def test_get_no_parent():
    fg1 = _init_fgnode("1", "RO")
    fg2 = _init_fgnode("2", "RC")
    parents = search_parents([fg1], fg2, mapper=default_mapper)
    assert parents is None


def test_get_correct_parent():
    fg1 = _init_fgnode("1", "RC")
    fg2 = _init_fgnode("2", "RO")
    fg3 = _init_fgnode("3", "ROR")
    parents = search_parents([fg1, fg2], fg3, mapper=default_mapper)
    assert parents == [fg2]


def test_get_multiple_parents():
    fg1 = _init_fgnode("1", "RC")
    fg2 = _init_fgnode("2", "RO")
    fg3 = _init_fgnode("3", "RCO")
    parents = search_parents([fg1, fg2], fg3, mapper=default_mapper)
    assert parents is not None
    assert 2 == len(parents)
    assert all([fg in parents for fg in [fg1, fg2]])


def test_get_multiple_unique_parents():
    fg11 = _init_fgnode("11", "RC")
    fg12 = _init_fgnode("12", "RO")
    fg2 = _init_fgnode("2", "RCOR")
    fg3 = _init_fgnode("3", "RCOCR")
    fg11.children = [fg2]
    fg12.children = [fg2]
    fg2.parents = [fg11, fg12]
    parents = search_parents([fg11, fg12], fg3, mapper=default_mapper)
    assert parents == [fg2]


def test_get_parent_recursive():
    fg1 = _init_fgnode("1", "RC")
    fg2 = _init_fgnode("2", "RCR")
    fg3 = _init_fgnode("3", "RCO")
    fg1.children = [fg2]
    fg2.parents = [fg1]
    parents = search_parents([fg1], fg3, mapper=default_mapper)
    assert parents == [fg2]


def _assert_structure(
    roots: list[FGTreeNode],
    config: FGConfig,
    exp_parents: FGConfig | list[FGConfig],
    exp_children: FGConfig | list[FGConfig] = [],
):
    def _get_node(
        nodes: list[FGTreeNode], config: FGConfig, node: FGTreeNode | None
    ) -> FGTreeNode | None:
        for n in nodes:
            if n.fgconfig == config:
                assert (
                    node is None or node == n
                ), "Found config multipel times in the tree."
                node = n
            node = _get_node(n.children, config, node)
        return node

    def _assert_get_node(config: FGConfig) -> FGTreeNode:
        node = _get_node(roots, config, None)
        assert node is not None, "Could not find config {} in tree.".format(config.name)
        return node

    node = _assert_get_node(config)

    pnodes = []
    for p in exp_parents if isinstance(exp_parents, list) else [exp_parents]:
        pnodes.append(_assert_get_node(p))
    exp_parent_nodes = sorted(pnodes, key=lambda x: x.order_id())

    cnodes = []
    for c in exp_children if isinstance(exp_children, list) else [exp_children]:
        cnodes.append(_assert_get_node(c))
    exp_children_nodes = sorted(cnodes, key=lambda x: x.order_id(), reverse=True)

    assert len(exp_parent_nodes) == len(node.parents), "Unequal number of parents."
    assert all(
        n in node.parents for n in exp_parent_nodes
    ), "Expected parent of {} to be {} but got {}.".format(
        node.fgconfig.name,
        [p.name for p in exp_parents] if isinstance(exp_parents, list) else None,
        [p.fgconfig.name for p in node.parents]
        if isinstance(node.parents, list)
        else None,
    )

    assert len(exp_children_nodes) == len(node.children), "Unequal number of children."
    assert all(
        n in node.children for n in exp_children_nodes
    ), "Expected sugroups of {} to be {} but got {}.".format(
        node.fgconfig.name,
        [c.name for c in exp_children] if isinstance(exp_children, list) else None,
        [c.fgconfig.name for c in node.children],
    )


def test_insert_child_between():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg1, fg3, fg2], mapper=default_mapper)
    _assert_structure(tree, fg1, [], fg2)
    _assert_structure(tree, fg2, fg1, fg3)
    _assert_structure(tree, fg3, fg2)


def test_insert_child_after():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg1, fg2, fg3], mapper=default_mapper)
    _assert_structure(tree, fg1, [], fg2)
    _assert_structure(tree, fg2, fg1, fg3)
    _assert_structure(tree, fg3, fg2)


def test_insert_new_root():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCR")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg2, fg3, fg1], mapper=default_mapper)
    _assert_structure(tree, fg1, [], fg2)
    _assert_structure(tree, fg2, fg1, fg3)
    _assert_structure(tree, fg3, fg2)


def test_insert_child_in_between_multiple():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RCOR")
    fg31 = FGConfig(name="31", pattern="RCOH")
    fg32 = FGConfig(name="32", pattern="RC=O")
    tree = build_config_tree_from_list([fg1, fg31, fg32, fg2], mapper=default_mapper)
    _assert_structure(tree, fg1, [], [fg2, fg32])
    _assert_structure(tree, fg2, fg1, fg31)
    _assert_structure(tree, fg31, fg2)
    _assert_structure(tree, fg32, fg1)


def test_insert_child_in_between_multiple_2():
    fg1 = FGConfig(name="1", pattern="RCR")
    fg2 = FGConfig(name="2", pattern="RCRCR")
    fg31 = FGConfig(name="31", pattern="RCCCR")
    fg32 = FGConfig(name="32", pattern="RCOCR")
    fg4 = FGConfig(name="4", pattern="RCCCCR")
    tree = build_config_tree_from_list(
        [fg1, fg31, fg32, fg4, fg2], mapper=default_mapper
    )
    _assert_structure(tree, fg1, [], [fg2])
    _assert_structure(tree, fg2, fg1, [fg31, fg32])
    _assert_structure(tree, fg31, fg2, fg4)
    _assert_structure(tree, fg32, fg2)


def test_multiple_parents():
    fg1 = FGConfig(name="1", pattern="RC")
    fg2 = FGConfig(name="2", pattern="RO")
    fg3 = FGConfig(name="3", pattern="RCOH")
    tree = build_config_tree_from_list([fg1, fg2, fg3], mapper=default_mapper)
    _assert_structure(tree, fg1, [], fg3)
    _assert_structure(tree, fg2, [], fg3)
    _assert_structure(tree, fg3, [fg1, fg2])


def test_tree_structure():
    def _check_fg(node: FGTreeNode):
        for c in node.children:
            print("Test {} -> {}.".format(node.fgconfig.name, c.fgconfig.name))
            assert node.fgconfig.pattern_len <= c.fgconfig.pattern_len
            assert True is map_subgraph_to_graph(
                c.fgconfig.pattern, node.fgconfig.pattern, mapper=default_mapper
            )
            assert False is map_subgraph_to_graph(
                node.fgconfig.pattern, c.fgconfig.pattern, mapper=default_mapper
            ), "Parent pattern {} contains child pattern {}.".format(
                node.fgconfig.pattern_str, c.fgconfig.pattern_str
            )
            _check_fg(c)

    provider = FGConfigProvider()
    for root_fg in provider.get_tree():
        _check_fg(root_fg)


@pytest.mark.parametrize(
    "fg_group,exp_pattern_len",
    [("carbonyl", 2), ("aldehyde", 3), ("ketone", 2), ("carboxylic_acid", 4)],
)
def test_pattern_len(fg_group, exp_pattern_len):
    provider = FGConfigProvider()
    fg = provider.get_by_name(fg_group)
    assert exp_pattern_len == fg.pattern_len


def test_config_name_uniqueness():
    name_list = []
    provider = FGConfigProvider()
    for fg in provider.config_list:
        assert fg.name not in name_list, "Config name '{}' already exists.".format(
            fg.name
        )
        name_list.append(fg.name)
    assert len(_default_fg_config) == len(name_list)


def test_config_pattern_validity():
    provider = FGConfigProvider()
    for c in provider.config_list:
        valid = False
        for _, sym in c.pattern.nodes(data=SYMBOL_KEY):  # type: ignore
            if sym != "C":
                valid = True
                break
        assert valid, "Carbon only groups are not allowed (Config: {})".format(c.name)


def test_config_pattern_uniqueness():
    pattern_list = []
    provider = FGConfigProvider()
    for fg in provider.config_list:
        assert (
            fg.pattern_str not in pattern_list
        ), "Config pattern '{}' already exists with name '{}'.".format(
            fg.pattern_str, fg.name
        )
        pattern_list.append(fg.pattern_str)
    assert len(_default_fg_config) == len(pattern_list)


def test_is_not_subgroup_if_matched_with_anti_pattern():
    fp = FGConfig(name="p", pattern="ROR", anti_pattern=["ROH"])
    fc = FGConfig(name="c", pattern="COH")
    result = is_subgroup(fp, fc, mapper=default_mapper)
    assert result is False


# def test_build_tree():
#     from fgutils.fgconfig import print_tree
#     provider = FGConfigProvider()
#     tree = provider.get_tree()
#     print_tree(tree)
#     assert False
