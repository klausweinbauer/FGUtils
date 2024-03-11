import copy
import networkx as nx

from fgutils.permutation import Mapper


def _get_neighbors(graph, idx, excluded_nodes=set()):
    return [
        (nidx, graph.nodes[nidx]["symbol"])
        for nidx in graph.neighbors(idx)
        if nidx not in excluded_nodes
    ]


def _get_symbol(graph, idx):
    return graph.nodes[idx]["symbol"]


def map_anchored_pattern(
    graph: nx.Graph, anchor: int, pattern: nx.Graph, pattern_anchor: int
):
    mapper = Mapper(wildcard="R", ignore_case=True)

    def _fit(idx, pidx, visited_nodes=set(), visited_pnodes=set(), indent=0):
        visited_nodes = copy.deepcopy(visited_nodes)
        visited_nodes.add(idx)
        visited_pnodes = copy.deepcopy(visited_pnodes)
        visited_pnodes.add(pidx)

        node_neighbors = _get_neighbors(graph, idx, visited_nodes)
        pnode_neighbors = _get_neighbors(pattern, pidx, visited_pnodes)

        nn_syms = [n[1] for n in node_neighbors]
        pnn_syms = [n[1] for n in pnode_neighbors]

        is_valid = False
        mappings = [(idx, pidx)]
        if len(pnn_syms) == 0:
            is_valid = True
        else:
            for n_mapping in mapper.permute(pnn_syms, nn_syms):
                _is_valid = True
                _mapping = set()
                _vnodes = set()
                _vpnodes = set()
                for pnn_i, nn_i in n_mapping:
                    pnn_idx = pnode_neighbors[pnn_i][0]
                    if nn_i == -1:
                        _vpnodes.add(pnn_idx)
                        continue
                    pnn_bond = pattern.edges[pidx, pnn_idx]["bond"]
                    nn_idx = node_neighbors[nn_i][0]
                    nn_bond = graph.edges[idx, nn_idx]["bond"]
                    if nn_bond == pnn_bond:
                        r_fit, r_mapping, r_vnodes = _fit(
                            nn_idx,
                            pnn_idx,
                            visited_nodes,
                            visited_pnodes,
                            indent=indent + 2,
                        )
                        if r_fit:
                            _vnodes.update(r_vnodes[0])
                            _vpnodes.update(r_vnodes[1])
                            _mapping.update(r_mapping)
                        else:
                            _is_valid = False
                    else:
                        _is_valid = False
                    if not _is_valid:
                        break
                if _is_valid:
                    is_valid = True
                    visited_nodes.update(_vnodes)
                    visited_pnodes.update(_vpnodes)
                    mappings.extend(_mapping)
                    break

        return is_valid, mappings, (visited_nodes, visited_pnodes)

    fit = False
    mapping = []
    sym = _get_symbol(graph, anchor)
    psym = _get_symbol(pattern, pattern_anchor)
    init_mapping = mapper.permute([psym], [sym])
    if init_mapping == [[(0, 0)]]:
        fit, mapping, _ = _fit(anchor, pattern_anchor)

    return fit, mapping


def map_pattern(
    graph: nx.Graph, anchor: int, pattern: nx.Graph, pattern_anchor: None | int = None
):
    if pattern_anchor is None:
        if len(pattern) == 0:
            return True, []
        for pidx in pattern.nodes:
            result = map_anchored_pattern(graph, anchor, pattern, pidx)
            if result[0]:
                return result
        return False, []
    else:
        return map_anchored_pattern(graph, anchor, pattern, pattern_anchor)
