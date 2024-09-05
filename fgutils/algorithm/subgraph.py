import copy
import networkx as nx

from fgutils.permutation import PermutationMapper

from fgutils.const import SYMBOL_KEY, BOND_KEY


def _get_neighbors(graph, idx, excluded_nodes=set()):
    return [
        (nidx, graph.nodes[nidx][SYMBOL_KEY])
        for nidx in graph.neighbors(idx)
        if nidx not in excluded_nodes
    ]


def _get_symbol(graph, idx):
    return graph.nodes[idx][SYMBOL_KEY]


def map_anchored_subgraph(
    graph: nx.Graph,
    anchor: int,
    subgraph: nx.Graph,
    subgraph_anchor: int,
    mapper: PermutationMapper,
):
    def _fit(idx, pidx, visited_nodes=set(), visited_pnodes=set(), indent=0):
        visited_nodes = copy.deepcopy(visited_nodes)
        visited_nodes.add(idx)
        visited_pnodes = copy.deepcopy(visited_pnodes)
        visited_pnodes.add(pidx)

        node_neighbors = _get_neighbors(graph, idx, visited_nodes)
        pnode_neighbors = _get_neighbors(subgraph, pidx, visited_pnodes)

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
                    pnn_bond = subgraph.edges[pidx, pnn_idx][BOND_KEY]
                    nn_idx = node_neighbors[nn_i][0]
                    nn_bond = graph.edges[idx, nn_idx][BOND_KEY]
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
    psym = _get_symbol(subgraph, subgraph_anchor)
    init_mapping = mapper.permute([psym], [sym])
    if init_mapping == [[(0, 0)]]:
        fit, mapping, _ = _fit(anchor, subgraph_anchor)

    return fit, mapping


def map_subgraph(
    graph: nx.Graph,
    anchor: int,
    subgraph: nx.Graph,
    mapper: PermutationMapper,
    subgraph_anchor: None | int = None,
):
    if subgraph_anchor is None:
        if len(subgraph) == 0:
            return [(True, [])]
        results = []
        for pidx in subgraph.nodes:
            result = map_anchored_subgraph(graph, anchor, subgraph, pidx, mapper)
            results.append(result)
        if len(results) > 0:
            return results
        else:
            return [(False, [])]
    else:
        return [map_anchored_subgraph(graph, anchor, subgraph, subgraph_anchor, mapper)]


def map_subgraph_to_graph(
    graph: nx.Graph, subgraph: nx.Graph, mapper: PermutationMapper
):
    for i in range(len(graph)):
        mappings = map_subgraph(graph, i, subgraph, mapper)
        for r, _ in mappings:
            if r is True:
                return True
    return False
