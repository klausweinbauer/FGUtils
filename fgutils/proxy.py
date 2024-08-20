from __future__ import annotations

import random
import networkx as nx

from fgutils.utils import split_its
from fgutils.parse import parse as pattern_to_graph


def relabel_graph(g):
    mapping = {}
    for i, u in enumerate(sorted(g.nodes)):
        mapping[u] = i
    return nx.relabel_nodes(g, mapping)


class ProxyGraph:
    """
    ProxyGraph is a essentially a subgraph used to expand molecules. If the
    node that is replaced by the pattern has more edges than the pattern has
    anchors, the last anchor will be used multiple times. In the default case
    the first node in the pattern will connect to all the neighboring nodes of
    the replaced node.

    :param pattern: String representation of the graph.
    :param anchor: A list of indices in the pattern that are used to
        connect to the parent graph. (Default = [0])
    """

    def __init__(self, pattern: str, anchor: list[int] = [0]):
        self.pattern = pattern
        if self.pattern is None:
            raise ValueError("Missing config 'pattern'.")
        self.anchor = anchor

    def __str__(self):
        return "{} Anchors: {}".format(self.pattern, self.anchor)


class ProxyGroup:
    """

    ProxyGroup is a collection of patterns that can be replaced for a labeled
    node in a graph. The node label is the respective group name where one of
    the patterns will be replaced.

    """

    def __init__(
        self,
        name,
        graphs: ProxyGraph | list[ProxyGraph] | None = None,
        pattern: str | list[str] | None = None,
    ):
        self.name = name
        if graphs is None:
            self.graphs = []
        else:
            self.graphs = graphs if isinstance(graphs, list) else [graphs]
        if pattern is not None:
            pattern = [pattern] if isinstance(pattern, str) else pattern
            for p in pattern:
                self.graphs.append(ProxyGraph(p))
        if len(self.graphs) == 0:
            raise ValueError("Group '{}' has no graphs.".format(name))

    def __str__(self):
        s = "ProxyGroup {}\n".format(self.name)
        for g in self.graphs:
            s += "  {}\n".format(g)
        return s

    @staticmethod
    def from_json_single(name, config: dict) -> ProxyGroup:
        graph_configs = config.get("graphs", [])
        graphs = []
        if isinstance(graph_configs, str):
            graph_configs = [graph_configs]
        if isinstance(graph_configs, dict):
            graph_configs = [graph_configs]
        print("HERE", graph_configs)
        for graph_config in graph_configs:
            if isinstance(graph_config, str):
                graph_config = {"pattern": graph_config}
            graphs.append(ProxyGraph(**graph_config))
        return ProxyGroup(name, graphs)

    @staticmethod
    def from_json(config: dict) -> dict[str, ProxyGroup]:
        groups = {}
        for name, group_config in config.items():
            if isinstance(group_config, str):
                group_config = [group_config]
            if isinstance(group_config, list):
                group_config = {"graphs": group_config}
            group = ProxyGroup.from_json_single(name, group_config)
            groups[name] = group
        return groups

    def get_graph(self):
        return random.sample(self.graphs, 1)[0]

    def __iter__(self):
        return self

    def __next__(self):
        return self.get_graph()


def _is_group_node(g: nx.Graph, idx: int, groups: dict[str, ProxyGroup]) -> bool:
    d = g.nodes[idx]
    return d["is_labeled"] and any([lbl in groups.keys() for lbl in d["labels"]])


def _has_group_nodes(g: nx.Graph, groups: dict[str, ProxyGroup]) -> bool:
    for u in g.nodes:
        if _is_group_node(g, u, groups):
            return True
    return False


def insert_groups(core: nx.Graph, groups: dict[str, ProxyGroup]) -> nx.Graph:
    _core = core.copy()
    idx_offset = len(core.nodes)
    for anchor, d in _core.nodes(data=True):
        if _is_group_node(_core, anchor, groups):
            sym = random.sample(d["labels"], 1)[0]
            group = next(groups[sym])
            h = pattern_to_graph(group.pattern, idx_offset=idx_offset)
            core = nx.compose(core, h)
            for i, (_, v, d) in enumerate(core.edges(anchor, data=True)):
                anchor_idx = i if len(group.anchor) > i else len(group.anchor) - 1
                core.add_edge(idx_offset + group.anchor[anchor_idx], v, **d)
            core.remove_node(anchor)
            idx_offset += len(h.nodes)
    core = relabel_graph(core)
    return core


def build_graph(pattern: str, groups: dict[str, ProxyGroup] = {}) -> nx.Graph:
    core = pattern_to_graph(pattern)
    while _has_group_nodes(core, groups):
        core = insert_groups(core, groups)
    return core


class Proxy:
    """
    Proxy is a generator class. It randomly extends a specific core graph by a
    set of subgraphs (groups).

    :param core: String representation of core graph. For example a
        specific functional group or a reaction center.
    :param groups: A list of groups to expand the core graph with.
    :param enable_aam: Flag to specify if the 'aam' label is set in the graph.
    """

    def __init__(
        self,
        core: str,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup] = [],
        enable_aam: bool = True,
    ):
        self.enable_aam = enable_aam
        self.core = core
        self.groups = {}
        if isinstance(groups, ProxyGroup):
            self.groups[groups.name] = groups
        elif isinstance(groups, list):
            for group in groups:
                self.groups[group.name] = group
        else:
            self.groups = groups

    def __str__(self):
        s = "ReactionProxy | Core: {} Enable AAM: {}\n".format(
            self.core, self.enable_aam
        )
        for group in self.groups.values():
            group_s = str(group)
            group_s = "\n  ".join(group_s.split("\n"))
            s += "  {}\n".format(group_s)
        return s

    @staticmethod
    def from_json(config: dict) -> Proxy:
        config["groups"] = ProxyGroup.from_json(config["groups"])
        return Proxy(**config)

    def generate(self):
        """Generate a new sample.

        :returns: A new graph.
        """
        graph = build_graph(self.core, self.groups)
        if self.enable_aam:
            for n in graph.nodes:
                graph.nodes[n]["aam"] = n + 1
        return graph

    def __iter__(self):
        return self

    def __next__(self):
        return self.generate()


class ReactionProxy(Proxy):
    def __init__(
        self,
        core: str,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup] = [],
        enable_aam: bool = True,
    ):
        super().__init__(core, groups, enable_aam)

    def generate(self):
        return split_its(super().generate())


class MolProxy(Proxy):
    def __init__(
        self,
        core: str,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup] = [],
    ):
        super().__init__(core, groups, False)

    def generate(self):
        return super().generate()
