import random
import networkx as nx

from fgutils.utils import split_its
from fgutils.parse import parse as pattern_to_graph


def relabel_graph(g):
    mapping = {}
    for i, u in enumerate(sorted(g.nodes)):
        mapping[u] = i
    return nx.relabel_nodes(g, mapping)


class ProxyGroup:
    def __init__(self, name, **kwargs):
        self.name = name
        self.pattern = kwargs.get("pattern", [])
        if not isinstance(self.pattern, list):
            self.pattern = [self.pattern]

    def get_pattern(self):
        return random.sample(self.pattern, 1)[0]

    def __iter__(self):
        return self

    def __next__(self):
        return self.get_pattern()


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
            group_pattern = next(groups[sym])
            h = pattern_to_graph(group_pattern, idx_offset=idx_offset)
            core = nx.compose(core, h)
            for _, v, d in core.edges(anchor, data=True):
                core.add_edge(idx_offset, v, **d)
            core.remove_node(anchor)
            idx_offset += len(h.nodes)
    core = relabel_graph(core)
    return core


def parse(pattern: str, groups: dict[str, ProxyGroup] = {}) -> nx.Graph:
    core = pattern_to_graph(pattern)
    while _has_group_nodes(core, groups):
        core = insert_groups(core, groups)
    return core


def parse_group_dict(config: dict) -> dict[str, ProxyGroup]:
    groups = {}
    for name, group_conf in config.items():
        if name in groups.keys():
            raise ValueError("Group name '{}' is already in use.".format(name))
        if isinstance(group_conf, str) or isinstance(group_conf, list):
            group_conf = {"pattern": group_conf}
        group = ProxyGroup(name, **group_conf)
        groups[name] = group
    return groups


class ReactionProxy:
    def __init__(self, config: dict, enable_aam=True):
        self.enable_aam = enable_aam
        self.core = config.get("core", "")
        self.groups = parse_group_dict(config.get("groups", {}))

    def generate(self):
        its = parse(self.core, self.groups)
        if self.enable_aam:
            for n in its.nodes:
                its.nodes[n]["aam"] = n
        return split_its(its)

    def __iter__(self):
        return self

    def __next__(self):
        return self.generate()
