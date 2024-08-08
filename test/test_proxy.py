import pytest
import networkx as nx

from fgutils.parse import parse as pattern_to_graph
from fgutils.proxy import ReactionProxy, parse_group_dict, parse


def assert_graph_eq(exp_graph, act_graph):
    def _nm(n1, n2):
        for k, v in n1.items():
            if k not in n2.keys() or n2[k] != v:
                return False
        for k, v in n2.items():
            if k not in n1.keys() or n1[k] != v:
                return False
        return True

    def _em(e1, e2):
        for k, v in e1.items():
            if k not in e2.keys() or e2[k] != v:
                return False
        for k, v in e2.items():
            if k not in e1.keys() or e1[k] != v:
                return False
        return True

    return nx.is_isomorphic(exp_graph, act_graph, node_match=_nm, edge_match=_em)


@pytest.mark.parametrize(
    "conf",
    (
        {"core": "A", "groups": {"test": "B"}},
        {"core": "A", "groups": {"test": ["B"]}},
        {"core": "A", "groups": {"test": {"pattern": "B"}}},
    ),
)
def test_init(conf):
    proxy = ReactionProxy(conf)
    assert "A" == proxy.core
    assert 1 == len(proxy.groups)
    assert "test" in proxy.groups.keys()
    assert ["B"] == proxy.groups["test"].pattern


@pytest.mark.parametrize(
    "core,group_conf,exp_result",
    (
        ("C#{group}", {"group": "N"}, "C#N"),
        ("C{group2}={group1}", {"group1": "O", "group2": "C"}, "CC=O"),
        ("C{group1}", {"group1": "{group2}C", "group2": "C=O"}, "CC(=O)C"),
    ),
)
def test_insert_groups(core, group_conf, exp_result):
    exp_graph = pattern_to_graph(exp_result)
    groups = parse_group_dict(group_conf)
    result = parse(core, groups)
    assert_graph_eq(exp_graph, result)


def test_reaction_generation():
    conf = {
        "core": "CC(<2,1>O)<0,1>{nucleophile}",
        "groups": {"nucleophile": "C#N"},
    }
    exp_g = pattern_to_graph("CC(=O).C#N")
    exp_h = pattern_to_graph("CC(O)C#N")
    proxy = ReactionProxy(conf)
    g, h = next(proxy)
    assert_graph_eq(g, exp_g)
    assert_graph_eq(h, exp_h)
