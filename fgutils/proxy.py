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
    def __init__(self, **kwargs):
        self.pattern = kwargs.get("pattern", None)
        if self.pattern is None:
            raise ValueError("Missing config 'pattern'.")
        self.anchor = kwargs.get("anchor", [0])


class ProxyGroup:
    def __init__(self, name, **kwargs):
        self.name = name
        self.graphs = []
        graph_configs = kwargs.get("graphs", [])
        if isinstance(graph_configs, str):
            graph_configs = [graph_configs]
        if isinstance(graph_configs, dict):
            graph_configs = [graph_configs]
        for graph_config in graph_configs:
            if isinstance(graph_config, str):
                graph_config = {"pattern": graph_config}
            self.graphs.append(ProxyGraph(**graph_config))
        if len(self.graphs) == 0:
            raise ValueError("Group '{}' has no graphs.".format(name))

    def get_graph(self):
        return random.sample(self.graphs, 1)[0]

    def __iter__(self):
        return self

    def __next__(self):
        return self.get_graph()

    @staticmethod
    def parse(config: dict):
        _groups = {}
        for name, group_conf in config.items():
            if isinstance(group_conf, str):
                group_conf = [group_conf]
            if isinstance(group_conf, list):
                group_conf = {"graphs": group_conf}
            group = ProxyGroup(name, **group_conf)
            _groups[name] = group
        return _groups


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


class ReactionProxy:
    def __init__(self, config: dict, enable_aam=True):
        self.enable_aam = enable_aam
        self.core = config.get("core", "")
        self.groups = ProxyGroup.parse(config.get("groups", {}))

    def generate(self):
        its = build_graph(self.core, self.groups)
        if self.enable_aam:
            for n in its.nodes:
                its.nodes[n]["aam"] = n + 1
        return split_its(its)

    def __iter__(self):
        return self

    def __next__(self):
        return self.generate()
