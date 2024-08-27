from __future__ import annotations

import random
import networkx as nx

from fgutils.utils import split_its
from fgutils.parse import Parser


def relabel_graph(g):
    mapping = {}
    for i, u in enumerate(sorted(g.nodes)):
        mapping[u] = i
    return nx.relabel_nodes(g, mapping)


class GraphSampler:
    """
    Base class for sampling ProxyGraphs.
    """

    def __init__(self, unique=True):
        self.__used_graphs = []
        self.unique = unique
        self.__index = 0

    def sample(self, graphs: list[ProxyGraph]) -> ProxyGraph | None:
        """
        Method to retrive a new sample from a list of graphs. The default
        behaviour is that each graph is sample at most once in the order they
        are provided. If ``unique`` is set to False the sample method will
        iterate over the list of graphs infinately.

        :param graphs: A list of graphs to sample from.

        :returns: Returns one of the graphs or None if sampling should stop.
        """
        if self.unique:
            for graph in graphs:
                if graph not in self.__used_graphs:
                    self.__used_graphs.append(graph)
                    return graph
            return None
        else:
            if self.__index >= len(graphs):
                self.__index = 0
            result = graphs[self.__index]
            self.__index += 1
            return result

    def __call__(self, graphs: list[ProxyGraph]) -> ProxyGraph | None:
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
    the patterns will be replaced. This class implements the iterator interface
    so it can be directly used in a for loop to retrive graphs::

        >>> proxy = ProxyGroup("group", pattern=["A", "B", "C"])
        >>> for graph in proxy:
        >>>     print(graph.pattern)
        A
        B
        C

    :param name: The name of the group.
    :param graphs: (optional) A list of subgraphs in this group.
    :param pattern:

        (optional) A list of graph descriptions. The patterns are converted to
        ProxyGraphs with one anchor at index 0. If you need more control over
        how the subgraphs are inserted use the ``graphs`` argument.

    :param sampler:

        (optional) An object or a function to retrive individual graphs from
        the list. The expected function interface is: ``func(list[ProxyGraph])
        -> ProxyGraph``. Implement the ``__call__`` method if you use a class.
    """

    def __init__(
        self,
        name,
        graphs: ProxyGraph | list[ProxyGraph] | None = None,
        pattern: str | list[str] | None = None,
        sampler=None,
    ):
        self.name = name
        self.sampler = GraphSampler() if sampler is None else sampler
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

    def get_graph(self) -> ProxyGraph | None:
        """Get the next graph. This method uses the sampler to select one of
        the specified graphs.

        :returns: A ProxyGraph or None if there is nothing more to sample.
        """
        return self.sampler(self.graphs)

    def __iter__(self):
        return self

    def __next__(self):
        graph = self.get_graph()
        if graph is None:
            raise StopIteration()
        return graph


def _is_group_node(g: nx.Graph, idx: int, groups: dict[str, ProxyGroup]) -> bool:
    d = g.nodes[idx]
    return d["is_labeled"] and any([lbl in groups.keys() for lbl in d["labels"]])


def _has_group_nodes(g: nx.Graph, groups: dict[str, ProxyGroup]) -> bool:
    for u in g.nodes:
        if _is_group_node(g, u, groups):
            return True
    return False


def insert_groups(
    core: nx.Graph, groups: dict[str, ProxyGroup], parser: Parser
) -> nx.Graph:
    """

    Replace labeled nodes in the core graph with groups. For each labeled node
    one label is chosen at random and replaced by the identically named group.
    This function does not resolve recursive labeled nodes. If a group has
    again labeled nodes they will be part of the result graph.

    :param core: The parent graph with labeled nodes.
    :param groups: A list of groups to replace the labeled nodes in the parent
        with.
    :param parser: The parser that is used to convert subgraph patterns into
        graphs.

    :returns: Returns the core graph with replaced nodes.
    """
    _core = core.copy()
    idx_offset = len(core.nodes)
    for anchor, d in _core.nodes(data=True):
        if d is None:
            raise ValueError("Expected a labeled graph.")
        if _is_group_node(_core, anchor, groups):
            sym = random.sample(d["labels"], 1)[0]
            graph = next(groups[sym])
            h = parser.parse(graph.pattern, idx_offset=idx_offset)
            core = nx.compose(core, h)
            for i, (_, v, d) in enumerate(core.edges(anchor, data=True)):  # type: ignore
                anchor_idx = i if len(graph.anchor) > i else len(graph.anchor) - 1
                core.add_edge(idx_offset + graph.anchor[anchor_idx], v, **d)  # type: ignore
            core.remove_node(anchor)
            idx_offset += len(h.nodes)
    core = relabel_graph(core)
    return core


def build_graph(pattern: str, parser: Parser, groups: dict[str, ProxyGroup] = {}):
    """
    Construct a graph from a pattern and replace all labeled nodes by the
    structure defined in the list of groups. This function resolves recursive
    labeling. The result graph has no labeled noes as long as a group is given
    for each label.

    :param pattern: The graph description for the parent graph.
    :param parser: The parser to use to convert patterns into structures.
    :param groups: A list of groups to replace the labeled nodes in the parent
        with.

    :returns: Returns the resulting graph with replaced nodes.
    """
    core = parser.parse(pattern)
    while _has_group_nodes(core, groups):
        core = insert_groups(core, groups, parser)
    return core


class Proxy:
    """
    Proxy is a generator class. It randomly extends a specific core graph by a
    set of subgraphs (groups). This class implements the iterator interface so
    it can be used in a for loop to generate samples::

        >>> proxy = Proxy("C{g}", ProxyGroup("g", pattern=["C", "O", "N"]))
        >>> for graph in proxy:
        >>>    print([d["symbol"] for n, d in graph.nodes(data=True)])
        ['C', 'C']
        ['C', 'O']
        ['C', 'N']

    :param core: A pattern string or ProxyGroup representing the core graph.
        For example a specific functional group or a reaction center.
    :param groups: A list of groups to expand the core graph with.
    :param enable_aam: Flag to specify if the 'aam' label is set in the graph.
    """

    def __init__(
        self,
        core: str | ProxyGroup,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup],
        enable_aam: bool = True,
        parser: Parser | None = None,
    ):
        self.enable_aam = enable_aam
        self.core = (
            ProxyGroup("__core__", pattern=core, sampler=GraphSampler(unique=False))
            if isinstance(core, str)
            else core
        )
        self.groups = groups
        if parser is None:
            self.parser = Parser(use_multigraph=True)
        else:
            self.parser = parser

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

    def generate(self):
        """Generates a new sample.

        :returns: A new graph.
        """
        core_graph = next(self.core)
        graph = build_graph(core_graph.pattern, self.parser, self.__groups)
        if self.enable_aam:
            for n in graph.nodes:
                graph.nodes[n]["aam"] = n + 1
        if isinstance(graph, nx.MultiGraph):
            graph = nx.Graph(graph)
        return graph

    def __iter__(self):
        return self

    def __next__(self):
        return self.generate()


class ReactionProxy(Proxy):
    """
    Proxy to generate reactions.
    """

    def __init__(
        self,
        core: str,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup] = [],
        enable_aam: bool = True,
        parser: Parser | None = None,
    ):
        super().__init__(core, groups, enable_aam, parser)

    def generate(self):
        return split_its(super().generate())


class MolProxy(Proxy):
    """
    Proxy to generate molecules.
    """

    def __init__(
        self,
        core: str,
        groups: ProxyGroup | list[ProxyGroup] | dict[str, ProxyGroup] = [],
        parser: Parser | None = None,
    ):
        super().__init__(core, groups, False, parser)

    def generate(self):
        return super().generate()
