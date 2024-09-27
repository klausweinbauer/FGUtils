import pytest
import numpy as np

from fgutils.parse import parse as pattern_to_graph, Parser
from fgutils.proxy import (
    Proxy,
    ProxyGraph,
    ReactionProxy,
    ProxyGroup,
    GraphSampler,
    build_group_tree,
    build_graphs,
)
from fgutils.const import SYMBOL_KEY

from test.my_asserts import assert_graph_eq


@pytest.mark.parametrize(
    "conf",
    (
        {
            "core": "A",
            "groups": {
                "test": {"graphs": [{"pattern": "BB", "anchor": [0], "order": 7}]}
            },
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
            "core": ["A"],
            "groups": {"test": "BB"},
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
    proxy = ReactionProxy.from_dict(conf)
    assert isinstance(proxy.core.graphs[0], ProxyGraph)
    assert "A" == proxy.core.graphs[0].pattern
    assert 1 == len(proxy.groups)
    group = proxy.groups[0]
    assert isinstance(group, ProxyGroup)
    assert "test" == group.name
    assert 1 == len(group.graphs)
    graph = group.graphs[0]
    assert isinstance(graph, ProxyGraph)
    assert "BB" == graph.pattern
    assert [0] == graph.anchor
    if "order" in graph.properties.keys():
        assert 7 == graph["order"]


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
    proxy = Proxy.from_dict({"core": core, "groups": config})
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
    proxy = Proxy.from_dict({"core": core, "groups": group_conf})
    result = next(proxy)
    assert_graph_eq(exp_graph, result)


def test_reaction_generation():
    group = ProxyGroup("nucleophile", "C#N")
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
    group = ProxyGroup("group", ["A", "B", "C"], sampler=lambda x: [x[0]])
    patterns = [group.sample_graphs()[0].pattern for _ in range(n)]
    assert ["A"] * n == patterns


def test_doc_example1():
    # example for fgutils.proxy.Proxy()
    pattern = ["C", "O", "N"]
    proxy = Proxy("C{g}", ProxyGroup("g", pattern))
    graphs = [graph for graph in proxy]
    assert 3 == len(graphs)
    for g, p in zip(graphs, pattern):
        assert 2 == len(g.nodes)
        assert 1 == len(g.edges)
        assert "C" == g.nodes(data=True)[0][SYMBOL_KEY]  # type: ignore
        assert p == g.nodes(data=True)[1][SYMBOL_KEY]  # type: ignore


def test_multiple_cores():
    cores = ["C{g0}", "O{g1}", "N{g2}"]
    core_group = ProxyGroup("core", [ProxyGraph(p) for p in cores])
    proxy = Proxy(
        core_group,
        [
            ProxyGroup("g{}".format(i), "C", sampler=GraphSampler(unique=True))
            for i in range(3)
        ],
    )
    graphs = [graph for graph in proxy]
    assert 3 == len(graphs)
    for g, p in zip(graphs, cores):
        assert 2 == len(g.nodes)
        assert 1 == len(g.edges)
        print(g.nodes(data=True))
        assert p[0] == g.nodes(data=True)[0][SYMBOL_KEY]  # type: ignore
        assert "C" == g.nodes(data=True)[1][SYMBOL_KEY]  # type: ignore


def test_create_proxy_tree():
    core_group = ProxyGroup("core", "{g1,g2,g3}")
    groups = [
        ProxyGroup("g1", "O"),
        ProxyGroup("g2", "C{g1}"),
        ProxyGroup("g3", "N"),
    ]
    tree = build_group_tree(core_group, groups)
    assert 5 == len(tree.nodes)
    assert 4 == len(tree.edges)
    assert 3 == np.sum([1 for n, d in tree.degree() if d == 1])


def test_proxy_tree_with_two_groups():
    core_group = ProxyGroup("core", "{g1,g2}{g1,g3}")
    groups = [
        ProxyGroup("g1", "O"),
        ProxyGroup("g2", "C{g1}"),
        ProxyGroup("g3", "N"),
    ]
    tree = build_group_tree(core_group, groups)
    print(tree.nodes)
    assert 6 == len(tree.nodes)
    assert 5 == len(tree.edges)
    assert 4 == len(list(tree.neighbors("core_#0")))


def test_sample_graphs_from_group():
    pattern = ["A", "B", "C"]
    group = ProxyGroup("group", pattern)
    graphs = group.sample_graphs()
    assert graphs is not None
    assert len(graphs) == len(pattern)
    for graph, p in zip(graphs, pattern):
        assert graph.pattern == p


def test_build_with_empyt_pattern_group():
    core = ProxyGraph("C{H}")
    group = {"H": ProxyGroup("H", "")}
    graphs = build_graphs(core, group, Parser())
    assert 1 == len(graphs)
    assert 1 == len(graphs[0].nodes)
    assert 0 == len(graphs[0].edges)


def test_insert_multiple_groups():
    parser = Parser()
    core = ProxyGraph("{g12}C{g3}")
    groups = {
        "g12": ProxyGroup("g12", ["C", "O"]),
        "g3": ProxyGroup("g3", "N"),
    }
    graphs = build_graphs(core, groups, parser)
    assert 2 == len(graphs)
    assert_graph_eq(parser("CCN"), graphs[0])
    assert_graph_eq(parser("OCN"), graphs[1])


def test_single_result_from_sampler():
    parser = Parser()
    core = ProxyGraph("{g12}C{g3}")
    groups = {
        "g12": ProxyGroup("g12", ["C", "O"], sampler=lambda g: g[0]),
        "g3": ProxyGroup("g3", "N"),
    }
    graphs = build_graphs(core, groups, parser)
    assert 1 == len(graphs)
    assert_graph_eq(parser("CCN"), graphs[0])


def test_build_with_multiple_graphs():
    parser = Parser()
    core = ProxyGraph("{g1}")
    groups = {
        "g1": ProxyGroup("g1", ["C", "O", "N"]),
    }
    graphs = build_graphs(core, groups, parser)
    assert 3 == len(graphs)
    assert_graph_eq(parser("C"), graphs[0])
    assert_graph_eq(parser("O"), graphs[1])
    assert_graph_eq(parser("N"), graphs[2])


def test_insert_aromatic_ring():
    parser = Parser()
    core = ProxyGraph("C1CCC{g}1")
    groups = {
        "g": ProxyGroup(
            "g",
            [ProxyGraph("CC", anchor=[0, 1]), ProxyGraph("c1ccccc1", anchor=[0, 5])],
        )
    }
    graphs = build_graphs(core, groups, parser)
    assert 2 == len(graphs)
    assert_graph_eq(parser("C1CCCCC1"), graphs[0])
    assert_graph_eq(parser("c1ccc2c(c1)CCCC2"), graphs[1])


def test_graph_dependency_with_custom_sampler():
    class CustomSampler:
        def __init__(self, group1, group2):
            if group1 == group2:
                raise ValueError()
            self.group1 = group1
            self.group2 = group2
            self.group1_order = 0
            self.group2_order = np.inf

        def __call__(
            self, graphs: list[ProxyGraph], group_name=None
        ) -> ProxyGraph | None:
            if group_name not in [self.group1, self.group2]:
                raise RuntimeError()
            result_group = None
            if group_name == self.group1:
                result_group = [
                    g for g in graphs if g["order"] < self.group2_order - 1
                ][0]
                self.group1_order = result_group["order"]
            if group_name == self.group2:
                result_group = [
                    g for g in graphs if g["order"] > self.group1_order + 1
                ][0]
                self.group2_order = result_group["order"]
            return result_group

    parser = Parser()
    sampler = CustomSampler("g1", "g2")
    graphs = [
        ProxyGraph("C", order=1),
        ProxyGraph("O", order=2),
        ProxyGraph("C=O", order=3),
    ]
    groups = [
        ProxyGroup("g1", graphs, sampler=sampler),
        ProxyGroup("g2", graphs, sampler=sampler),
    ]
    proxy_1 = Proxy("{g1}C{g2}", groups)
    proxy_2 = Proxy("{g2}C{g1}", groups)
    result_1 = next(proxy_1)
    result_2 = next(proxy_2)
    assert_graph_eq(parser("CCC=O"), result_1)
    assert_graph_eq(parser("CCC=O"), result_2)
