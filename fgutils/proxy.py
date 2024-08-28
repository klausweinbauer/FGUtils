from __future__ import annotations

import random
import numpy as np
import networkx as nx

from fgutils.utils import split_its
from fgutils.parse import Parser
from fgutils.const import IS_LABELED_KEY, LABELS_KEY


def relabel_graph(g):
    mapping = {}
    for i, u in enumerate(sorted(g.nodes)):
        mapping[u] = i
    return nx.relabel_nodes(g, mapping)


class GraphSampler:
    """
    Base class for sampling ProxyGraphs.
    """

    def __init__(self, unique=False):
        self.unique = unique
        self.__hist = []

    def sample(self, graphs: list[ProxyGraph]) -> list[ProxyGraph] | None:
        """
        Method to retrive a new sample from a list of graphs.

        :param graphs: A list of graphs to sample from.

        :returns: Returns one or more graphs from the list or None if sampling
            should stop.
        """
        if self.unique:
            for g in graphs:
                if g not in self.__hist:
                    self.__hist.append(g)
                    return [g]
            return None
        else:
            return graphs

    def __call__(self, graphs: list[ProxyGraph]) -> list[ProxyGraph] | None:
        return self.sample(graphs)


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
    the patterns will be replaced

    :param name: The name of the group.
    :param graphs: (optional) A list of subgraphs in this group.
    :param pattern: (optional) A list of graph descriptions. The patterns are
        converted to ProxyGraphs with one anchor at index 0. If you need more
        control over how the subgraphs are inserted use the ``graphs``
        argument.
    :param sampler: (optional) An object or a function to retrive individual
        graphs from the list. The expected function interface is:
        ``func(list[ProxyGraph]) -> list[ProxyGraph]``. Implement the ``__call__``
        method if you use a class.
    :param unique: Argument to specify if graphs can be returned multiple
        times. This only takes effect if sampler is not set. (Default = False)
    """

    def __init__(
        self,
        name,
        graphs: str | list[str] | ProxyGraph | list[ProxyGraph],
        sampler=None,
        unique=False,
    ):
        self.name = name
        self.sampler = GraphSampler(unique=unique) if sampler is None else sampler
        self.graphs = graphs

    def __str__(self):
        s = "ProxyGroup {}\n".format(self.name)
        for g in self.graphs:
            s += "  {}\n".format(g)
        return s

    @property
    def graphs(self) -> list[ProxyGraph]:
        return self.__graphs

    @graphs.setter
    def graphs(self, value):
        self.__graphs: list[ProxyGraph] = []
        if isinstance(value, str) or isinstance(value, ProxyGraph):
            value = [value]
        if isinstance(value, list):
            for graph in value:
                if isinstance(graph, str):
                    self.graphs.append(ProxyGraph(graph))
                elif isinstance(graph, ProxyGraph):
                    self.graphs.append(graph)
                else:
                    raise TypeError(
                        "Invalid type '{}' for graphs argument.".format(
                            type(self.__graphs)
                        )
                    )
        else:
            raise TypeError(
                "Invalid type '{}' for graphs argument.".format(type(self.__graphs))
            )
        if 0 == len(self.__graphs):
            raise ValueError("Group '{}' has no graphs.".format(self.name))

    @staticmethod
    def from_json_single(name, config: dict) -> ProxyGroup:
        graph_configs = config.get("graphs", [])
        graphs = []
        if isinstance(graph_configs, str):
            graph_configs = [graph_configs]
        if isinstance(graph_configs, dict):
            graph_configs = [graph_configs]
        for graph_config in graph_configs:
            if isinstance(graph_config, str):
                graph_config = {"pattern": graph_config}
            graphs.append(ProxyGraph(**graph_config))  # type: ignore
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

    def sample_graphs(self) -> list[ProxyGraph] | None:
        """Sample graphs. This method uses the sampler to select a list of
        graphs.

        :returns: A list of ProxyGraph objects or None if there is nothing more
        to sample.
        """
        result = self.sampler(self.graphs)
        if result is not None and not isinstance(result, list):
            result = [result]
        return result


def _is_group_node(g: nx.Graph, idx: int, groups: dict[str, ProxyGroup]) -> bool:
    d = g.nodes[idx]
    return d["is_labeled"] and any([lbl in groups.keys() for lbl in d["labels"]])


def _get_next_group_node(g: nx.Graph, groups: dict[str, ProxyGroup]) -> int | None:
    for anchor, d in g.nodes(data=True):
        if d is None:
            raise ValueError("Expected a labeled graph.")
        if _is_group_node(g, anchor, groups):
            return anchor
    return None


def replace_node(graph, node, replacement_graph: ProxyGraph, parser: Parser):
    """Replace ``node`` in ``graph`` with ``replacement_graph``.

    :param graph: The graph where a node should be replace by a subgraph.
    :param node: The node to replace in graph.
    :param replacement_graph: The subgraph that is inserted instead of the
        node.
    :param parser: The parser to use to convert the pattern into structure.

    :returns: Returns a new graph with ``node`` replace by
        ``replacement_graph``.
    """
    idx_offset = len(graph.nodes)
    h = parser.parse(replacement_graph.pattern, idx_offset=idx_offset)
    graph = nx.compose(graph, h)
    if len(h.nodes) > 0:
        for i, (_, v, d) in enumerate(graph.edges(node, data=True)):
            anchor_idx = i
            if len(replacement_graph.anchor) <= i:
                anchor_idx = len(replacement_graph.anchor) - 1
            graph.add_edge(idx_offset + replacement_graph.anchor[anchor_idx], v, **d)
    graph.remove_node(node)
    graph = relabel_graph(graph)
    return graph


def replace_next_node(graph, groups, parser: Parser):
    """Replace the next labeled node in ``graph`` with the respective group.

    :param graph: The graph where a node should be replace by a subgraph.
    :param groups: A list of groups to replace the labeled nodes in the parent
        with. The dictionary keys must be the group names.
    :param parser: The parser to use to convert the pattern into structure.

    :returns: Returns a list of new graphs with ``node`` replace with the
        groups selected by the label sampler. Or None if nothing is left to
        replace.
    """
    result_graphs = []
    anchor = _get_next_group_node(graph, groups)
    if anchor is None:
        return None
    anchor_labels = graph.nodes[anchor][LABELS_KEY]
    group_labels = []
    for anchor_label in anchor_labels:
        if anchor_label in groups.keys():
            group_labels.append(anchor_label)
    if len(group_labels) != 1:
        raise RuntimeError(
            "Multiple group labels found on node ({}).".format(group_labels)
        )
    group_name = group_labels[0]
    if groups[group_name].name != group_name:
        raise ValueError(
            "Dictionary key '{}' does not match group name '{}'.".format(
                group_name, groups[group_name].name
            )
        )
    sub_graphs = groups[group_name].sample_graphs()
    if sub_graphs is None:
        raise StopIteration()
    for sub_graph in sub_graphs:
        _graph = graph.copy()
        _graph = replace_node(_graph, anchor, sub_graph, parser)
        result_graphs.append(_graph)
    return result_graphs


def build_graphs(
    core: ProxyGraph,
    groups: dict[str, ProxyGroup],
    parser: Parser,
):
    """
    Replace labeled nodes in the core graph with groups. For each labeled node
    the label sampler will select a set of labels (group names) used to replace
    the node by the specified subgraphs. Selecting multiple labels per labeled
    node results in multiple result graphs.

    :param core: The parent graph with labeled nodes.
    :param groups: A list of groups to replace the labeled nodes in the parent
        with. The dictionary keys must be the group names.
    :param parser: The parser that is used to convert subgraph patterns into
        graphs.

    :returns: Returns a list of graphs with replaced nodes.
    """
    result_set = []
    working_set = [parser(core.pattern)]
    while len(working_set) > 0:
        _working_set = []
        for ws_graph in working_set:
            result_graphs = replace_next_node(ws_graph, groups, parser)
            if result_graphs is None:
                result_set.append(ws_graph)
            else:
                _working_set.extend(result_graphs)
        working_set = _working_set
    return result_set


class Proxy:
    """
    Proxy is a generator class. It randomly extends a specific core graph by a
    set of subgraphs (groups). This class implements the iterator interface so
    it can be used in a for loop to generate samples::

        >>> proxy = Proxy("C{g}", ProxyGroup("g", ["C", "O", "N"], unique=True))
        >>> for graph in proxy:
        >>>    print([d["symbol"] for n, d in graph.nodes(data=True)])
        ['C', 'C']
        ['C', 'O']
        ['C', 'N']

    :param core: A pattern string or ProxyGroup representing the core graph.
        For example a specific functional group or a reaction center.
    :param groups: A list of groups to expand the core graph with.
    :param enable_aam: Flag to specify if the 'aam' label is set in the graph.
    :param parser: The parser to use to convert patterns into structures.
    """

    def __init__(
        self,
        core: str | list[str] | ProxyGroup,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup],
        enable_aam: bool = True,
        parser: Parser | None = None,
    ):
        self.enable_aam = enable_aam
        self.core = (
            ProxyGroup("__core__", core) if not isinstance(core, ProxyGroup) else core
        )
        self.groups = groups
        if parser is None:
            self.parser = Parser(use_multigraph=True)
        else:
            self.parser = parser

        self.__active_generator = self.__generate()

    @property
    def groups(self) -> list[ProxyGroup]:
        return list(self.__groups.values())

    @groups.setter
    def groups(self, value):
        self.__groups = {}
        if isinstance(value, ProxyGroup):
            self.__groups[value.name] = value
        elif isinstance(value, list):
            for group in value:
                self.__groups[group.name] = group
        else:
            self.__groups = value

    def __str__(self):
        s = "ReactionProxy | Core: {} Enable AAM: {}\n".format(
            self.core, self.enable_aam
        )
        for group in self.__groups.values():
            group_s = str(group)
            group_s = "\n  ".join(group_s.split("\n"))
            s += "  {}\n".format(group_s)
        return s

    @staticmethod
    def from_json(config: dict) -> Proxy:
        config["groups"] = ProxyGroup.from_json(config["groups"])
        return Proxy(**config)

    def __generate(self):
        core_graphs = self.core.sample_graphs()
        while core_graphs is not None:
            try:
                for core_graph in core_graphs:
                    graphs = build_graphs(core_graph, self.__groups, self.parser)
                    for graph in graphs:
                        if self.enable_aam:
                            for n in graph.nodes:
                                graph.nodes[n]["aam"] = n + 1
                        if isinstance(graph, nx.MultiGraph):
                            graph = nx.Graph(graph)
                        yield graph
                core_graphs = self.core.sample_graphs()
            except StopIteration:
                return

    def get_next(self):
        return next(self.__active_generator)

    def __iter__(self):
        return self

    def __next__(self):
        return self.get_next()


class ReactionProxy(Proxy):
    """
    Proxy to generate reactions.
    """

    def __init__(
        self,
        core: str | list[str] | ProxyGroup,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup] = [],
        enable_aam: bool = True,
        parser: Parser | None = None,
    ):
        super().__init__(core, groups, enable_aam, parser)

    def get_next(self):
        """
        Generate a new reaction sample. The reaction proxy returns two graphs G
        and H. G is the reactant graph and H is the product graph.

        :returns: A tuple of two graphs (G, H) representing the reaction G
            \u2192 H.
        """
        return split_its(super().get_next())


class MolProxy(Proxy):
    """
    Proxy to generate molecules.
    """

    def __init__(
        self,
        core: str | list[str] | ProxyGroup,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup] = [],
        parser: Parser | None = None,
    ):
        super().__init__(core, groups, False, parser)


def build_group_tree(
    core: ProxyGroup, groups: ProxyGroup | list[ProxyGroup], parser=None
) -> nx.Graph:
    """
    Constructs a tree of all possible graph instantiations. The number of leave
    nodes in this tree is the number of possible samples.

    :param core: The ProxyGroup that serves as core group.
    :param groups: A list of groups to replace labeled nodes.
    :param parser: (optional) A parser to use for conversion from pattern to
        structures.

    :returns: Returns the construction tree as nx.Graph object.
    """

    def _add_node(
        tree: nx.Graph, group: ProxyGroup, groups: dict[str, ProxyGroup], parser: Parser
    ) -> str:
        node_name = "{}_#{}".format(group.name, len(tree.nodes))
        tree.add_node(node_name)
        for g in group.graphs:
            graph = parser(g.pattern)
            for _, d in graph.nodes(data=True):
                if d is None:
                    raise ValueError("Expected labeled graph.")
                if d[IS_LABELED_KEY]:
                    keys = d[LABELS_KEY]
                    for k in keys:
                        _group = groups[k]
                        _node_name = _add_node(tree, _group, groups, parser)
                        tree.add_edge(node_name, _node_name)
        return node_name

    if parser is None:
        parser = Parser()
    if not isinstance(groups, list):
        groups = [groups]
    group_dict = {}
    for g in groups:
        group_dict[g.name] = g
    tree = nx.Graph()
    _add_node(tree, core, group_dict, parser)
    return tree
