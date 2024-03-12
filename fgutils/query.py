import copy
import collections

from fgutils.utils import add_implicit_hydrogens
from fgutils.mapping import map_pattern
from fgutils.fgconfig import FGConfig, build_FG_tree, FGTreeNode


def is_functional_group(graph, index: int, config: FGConfig):
    max_id = len(graph)
    graph = add_implicit_hydrogens(copy.deepcopy(graph))

    is_fg, mapping = map_pattern(graph, index, config.pattern)
    fg_indices = []
    if is_fg:
        fg_indices = [
            m_id
            for m_id, fg_id in mapping
            if fg_id in config.group_atoms and m_id < max_id
        ]
        is_fg = index in fg_indices

    if is_fg:
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


def get_functional_groups_raw(graph) -> tuple[dict, list[str]]:
    def _query(nodes: list[FGTreeNode], graph, idx, checked_groups = []):
        fg_groups = []
        fg_indices = []
        for node in nodes:
            if node.fgconfig.name in checked_groups:
                continue
            is_fg, indices = is_functional_group(graph, idx, node.fgconfig)
            if is_fg:
                checked_groups.append(node.fgconfig.name)
                fg_groups.append(node.fgconfig.name)
                fg_indices.append(indices)
                _fg_groups, _fg_indices = _query(node.children, graph, idx, checked_groups)
                fg_groups.extend(_fg_groups)
                fg_indices.extend(_fg_indices)
        return fg_groups, fg_indices

    fg_candidate_ids = [
        n_id for n_id, n_sym in graph.nodes(data="symbol") if n_sym not in ["H", "C"]
    ]
    roots = build_FG_tree()
    idx_map = collections.defaultdict(lambda: [])
    groups = []
    for atom_id in fg_candidate_ids:
        fg_groups, fg_indices = _query(roots, graph, atom_id)
        print("Check atom {} | Groups: {} Indices: {}".format(atom_id, fg_groups, fg_indices))
        if len(fg_groups) > 0:
            for _group, _indices in zip(fg_groups, fg_indices):
                assert atom_id in _indices
                _i = len(groups)
                groups.append(_group)
                for _idx in _indices:
                    idx_map[_idx].append(_i)
    return dict(idx_map), groups
