import pytest
import networkx as nx

from fgutils.utils import print_graph
from fgutils.parse import parse as pattern_to_graph
from fgutils.proxy import Proxy, ProxyGraph, ReactionProxy, ProxyGroup, build_graph


def assert_graph_eq(exp_graph, act_graph, ignore_keys=["aam"]):
    def _nm(n1, n2):
        for k, v in n1.items():
            if k in ignore_keys:
                continue
            if k not in n2.keys() or n2[k] != v:
                return False
        for k, v in n2.items():
            if k in ignore_keys:
                continue
            if k not in n1.keys() or n1[k] != v:
                return False
        return True

    def _em(e1, e2):
        for k, v in e1.items():
            if k in ignore_keys:
                continue
            if k not in e2.keys() or e2[k] != v:
                return False
        for k, v in e2.items():
            if k in ignore_keys:
                continue
            if k not in e1.keys() or e1[k] != v:
                return False
        return True

    is_isomorphic = nx.is_isomorphic(
        exp_graph, act_graph, node_match=_nm, edge_match=_em
    )
    assert is_isomorphic, "Graphs are not isomorphic."


@pytest.mark.parametrize(
    "conf",
    (
        {
            "core": "A",
            "groups": {"test": {"graphs": [{"pattern": "BB", "anchor": [0]}]}},
        },
        {
            "core": "A",
            "groups": {"test": {"graphs": {"pattern": "BB", "anchor": [0]}}},
        },
        {
            "core": "A",
            "groups": {"test": {"graphs": ["BB"]}},
        },
        {
            "core": "A",
            "groups": {"test": {"graphs": "BB"}},
        },
        {
            "core": "A",
            "groups": {"test": ["BB"]},
        },
        {
            "core": "A",
            "groups": {"test": "BB"},
        },
    ),
)
def test_init(conf):
    proxy = ReactionProxy.from_json(conf)
    assert "A" == proxy.core
    assert 1 == len(proxy.groups)
    assert "test" in proxy.groups.keys()
    assert 1 == len(proxy.groups["test"].graphs)
    assert "BB" == proxy.groups["test"].graphs[0].pattern
    assert [0] == proxy.groups["test"].graphs[0].anchor


@pytest.mark.parametrize(
    "core,config,exp_graph",
    (
        ("{g}1CC1", {"g": "C"}, "C1CC1"),
        ("C{g}C", {"g": {"graphs": {"pattern": "OC", "anchor": [1]}}}, "CC(O)C"),
        ("{g}1CC1", {"g": "CC"}, "C1(C)CC1"),
        (
            "C1{g}C1",
            {"g": {"graphs": {"pattern": "CCC", "anchor": [0, 2]}}},
            "C1CCCC1",
        ),
    ),
)
def test_build_graph(core, config, exp_graph):
    exp_graph = pattern_to_graph(exp_graph)
    proxy = Proxy.from_json({"core": core, "groups": config})
    result = next(proxy)
    assert_graph_eq(exp_graph, result)


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
    proxy = Proxy.from_json({"core": core, "groups": group_conf})
    result = next(proxy)
    assert_graph_eq(exp_graph, result)


def test_reaction_generation():
    group = ProxyGroup("nucleophile", pattern="C#N")
    proxy = ReactionProxy("CC(<2,1>O)<0,1>{nucleophile}", groups=group)
    exp_g = pattern_to_graph("CC(=O).C#N")
    exp_h = pattern_to_graph("CC(O)C#N")
    g, h = next(proxy)
    assert_graph_eq(g, exp_g)
    assert_graph_eq(h, exp_h)
