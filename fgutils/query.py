import copy
import collections

from fgutils.utils import add_implicit_hydrogens
from fgutils.mapping import map_pattern
from fgutils.fgconfig import FGConfig, get_FG_tree


def get_functional_group(graph, index: int, config: FGConfig):
    is_fg = False
    max_id = len(graph)
    graph = add_implicit_hydrogens(copy.deepcopy(graph))

    is_match, mapping = map_pattern(graph, index, config.pattern)
    fg_indices = []
    if is_match:
        fg_indices = [
            m_id
            for m_id, fg_id in mapping
            if fg_id in config.group_atoms and m_id < max_id
        ]
        is_fg = is_fg or index in fg_indices

    last_len = config.max_pattern_size
    for apattern, apattern_size in sorted(
        [(m, m.number_of_nodes()) for m in config.anti_pattern],
        key=lambda x: x[1],
        reverse=True,
    ):
        if not is_fg:
            break
        if last_len > apattern_size:
            last_len = apattern_size
        is_match, _ = map_pattern(graph, index, apattern)
        is_fg = is_fg and not is_match
    return is_fg, sorted(fg_indices)


def check_functional_group(graph, idx, fg):
    def _get_recursive(fg, groups, indices):
        is_fg, fg_indices = get_functional_group(graph, idx, fg)
        if is_fg:
            groups.append(fg.name)
            indices.append(fg_indices)
            for sfg in fg.subgroups:
                if _get_recursive(sfg, groups, indices):
                    break
            return True
        return False

    groups = []
    indices = []
    _get_recursive(fg, groups, indices)
    return groups, indices


def get_functional_groups_raw(graph) -> tuple[dict, list[str]]:
    fg_nodes = [
        n_id for n_id, n_sym in graph.nodes(data="symbol") if n_sym not in ["H", "C"]
    ]
    idx_map = collections.defaultdict(lambda: [])
    groups = []
    for fg_node in fg_nodes:
        for fg in get_FG_tree():
            fg_groups, fg_indices = check_functional_group(graph, fg_node, fg)
            if len(fg_groups) > 0:
                for _group, _indices in zip(fg_groups, fg_indices):
                    _i = len(groups)
                    groups.append(_group)
                    for _idx in _indices:
                        idx_map[_idx].append(_i)
    return dict(idx_map), groups
