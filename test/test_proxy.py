import pytest
import networkx as nx

from fgutils.parse import parse as pattern_to_graph, Parser
from fgutils.proxy import Proxy, ProxyGraph, ReactionProxy, ProxyGroup


def _node_match(n1, n2, ignore_keys):
    for k, v in n1.items():
        if k in ignore_keys:
            continue
        if k not in n2.keys() or n2[k] != v:
            print("unequal 1", k)
            return False
    for k, v in n2.items():
        if k in ignore_keys:
            continue
        if k not in n1.keys() or n1[k] != v:
            print("unequal 2", k)
            return False
    return True


def _edge_match(e1, e2, ignore_keys):
    for k, v in e1.items():
        if k in ignore_keys:
            continue
        if k not in e2.keys() or e2[k] != v:
            print("unequal 3", k)
            return False
    for k, v in e2.items():
        if k in ignore_keys:
            continue
        if k not in e1.keys() or e1[k] != v:
            print("unequal 4", k)
            return False
    return True


def assert_graph_eq(exp_graph, act_graph, ignore_keys=["aam"]):
    def _nm(n1, n2):
        return _node_match(n1, n2, ignore_keys)

    def _em(e1, e2):
        return _edge_match(e1, e2, ignore_keys)

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
    assert "A" == proxy.core.graphs[0].pattern
    assert 1 == len(proxy.groups)
    assert "test" == proxy.groups[0].name
    assert 1 == len(proxy.groups[0].graphs)
    assert "BB" == proxy.groups[0].graphs[0].pattern
    assert [0] == proxy.groups[0].graphs[0].anchor


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


def test_multigraph_reaction_generation():
    group_diene = ProxyGroup("diene", ProxyGraph("C<2,1>C<1,2>C<2,1>C", anchor=[0, 3]))
    group_dienophile = ProxyGroup("dienophile", ProxyGraph("C<2,1>C", anchor=[0, 1]))
    proxy = ReactionProxy(
        "{diene}1<0,1>{dienophile}<0,1>1",
        groups=[group_diene, group_dienophile],
        parser=Parser(use_multigraph=True),
    )
    exp_g = pattern_to_graph("C=CC=C.C=C")
    exp_h = pattern_to_graph("C1C=CCCC1")
    g, h = next(proxy)
    assert_graph_eq(g, exp_g)
    assert_graph_eq(h, exp_h)


def test_graph_sampling():
    n = 5
    proxy = ProxyGroup("group", pattern=["A", "B", "C"], sampler=lambda x: x[0])
    patterns = [next(proxy).pattern for _ in range(n)]
    assert ["A"] * n == patterns


def test_group_sampling():
    pattern = ["C", "O", "N"]
    proxy = Proxy("{g}C", ProxyGroup("g", pattern=pattern))
    graphs = [graph for graph in proxy]
    assert 3 == len(graphs)
    for g, p in zip(graphs, pattern):
        assert 2 == len(g.nodes)
        assert 1 == len(g.edges)
        assert "C" == g.nodes(data=True)[0]["symbol"]  # type: ignore
        assert p == g.nodes(data=True)[1]["symbol"]  # type: ignore


def test_multiple_cores():
    cores = ["C{g0}", "O{g1}", "N{g2}"]
    core_group = ProxyGroup("core", [ProxyGraph(p) for p in cores])
    proxy = Proxy(
        core_group,
        [ProxyGroup("g{}".format(i), pattern="C") for i in range(3)],
    )
    graphs = [graph for graph in proxy]
    assert 3 == len(graphs)
    for g, p in zip(graphs, cores):
        assert 2 == len(g.nodes)
        assert 1 == len(g.edges)
        print(g.nodes(data=True))
        assert p[0] == g.nodes(data=True)[0]["symbol"]  # type: ignore
        assert "C" == g.nodes(data=True)[1]["symbol"]  # type: ignore
