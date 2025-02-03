import copy
import networkx as nx

from fgutils.permutation import PermutationMapper, MappingMatrix

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
) -> tuple[bool, list[tuple[int, int]], tuple[set[int], set[int]]]:
    """Map an anchored subgraph into an anchored parent graph. The anchor in
    the subgraph and the anchor in the parent graph is the first fixed mapping.
    The subgraph is aligned to the parent graph based on this anchor. This
    function will only find one alignment and not all possible alignments. It
    will return once a solution is found.

    :param graph: The paraent graph in which to align the subgraph.
    :param anchor: The node index for the anchor in the parent graph.
    :param subgraph: The subgraph to align in the parent.
    :param subgraph_anchor: The node index for the anchor in the subgraph.
    :param mapper: A :py:class:`~fgutils.permutation.PermutationMapper`
        instance. The permutation mapper specifies which nodes in the subgraph
        can align with nodes in the parent.

    :returns: The function returns a 3-tuple of the form ``(is_valid, mapping,
        visited_nodes)``. The first entry ``is_valid`` is a boolean indicating
        if a mapping was found. The second entry ``mapping`` is a list of
        integer tuples. Each integer tuple encodes one node mapping where the
        first element is the node index in the parent graph and the second
        entry is the node index in the subgraph. The third element is a tuple
        of sets with node indices for the visited nodes during graph alignment
        in the subgraph and parent graph, respectively.
    """

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
    visited_nodes = set([anchor]), set([subgraph_anchor])
    if init_mapping == [[(0, 0)]]:
        fit, mapping, visited_nodes = _fit(anchor, subgraph_anchor)

    return fit, mapping, visited_nodes


def map_subgraph(
    graph: nx.Graph,
    anchor: int,
    subgraph: nx.Graph,
    mapper: PermutationMapper,
    subgraph_anchor: None | int = None,
) -> list[tuple[bool, list[tuple[int, int]]]]:
    """Finds an anchored map of a subgraph in a parent graph.

    :param anchor: The node index in the parent graph where the subgraph is
        anchored.
    :param subgraph: The subgraph to align at the anchor node.
    :param mapper: A :py:class:`~fgutils.permutation.PermutationMapper`
        instance. The permutation mapper specifies which nodes in the subgraph
        can align with nodes in the parent.
    :param subgraph_anchor: (Optional) A node index in the subgraph for a fixed
        anchor. If provided this is the first fixed mapping between subgraph
        and parent.

    :returns: Returns a list of tuples: [(is_valid, mapping)]. Each list entry
        corresponds to one subgraph mapping. The first element in the tuple is
        a boolean indicating if a mapping was found. The second element is the
        mapping. The mapping is a list of integer 2-tuples where the first
        element is the node index in the subgraph and the second element is the
        node index in the parent graph.
    """
    if subgraph_anchor is None:
        if len(subgraph) == 0:
            return [(True, [])]
        results = []
        for pidx in subgraph.nodes:
            is_valid, mapping, _ = map_anchored_subgraph(
                graph, anchor, subgraph, pidx, mapper
            )
            results.append((is_valid, mapping))
        if len(results) > 0:
            return results
        else:
            return [(False, [])]
    else:
        is_valid, mapping, _ = map_anchored_subgraph(
            graph, anchor, subgraph, subgraph_anchor, mapper
        )
        return [(is_valid, mapping)]


def map_subgraph_to_graph(
    graph: nx.Graph, subgraph: nx.Graph, mapper: PermutationMapper
):
    for i in range(len(graph)):
        mappings = map_subgraph(graph, i, subgraph, mapper)
        for r, _ in mappings:
            if r is True:
                return True
    return False


def map_subgraph2(
    graph: nx.Graph,
    subgraph: nx.Graph,
    mapper: PermutationMapper,
    matrix: MappingMatrix | None = None,
):
    if len(list(nx.connected_components(subgraph))) > 1:
        raise ValueError("Do not use map_subgraph2 for disconnected subgraphs.")
    g_labels = [d[SYMBOL_KEY] for _, d in graph.nodes(data=True)]
    sub_labels = [d[SYMBOL_KEY] for _, d in subgraph.nodes(data=True)]
    if matrix is None:
        matrix = MappingMatrix(sub_labels, g_labels, mapper)
    mapping_sym = matrix.min_mapping_symbol(sub_labels, g_labels)
    assert mapping_sym is not None
    possible_mappings = []
    for u, ud in subgraph.nodes(data=True):
        if ud[SYMBOL_KEY] != mapping_sym[0]:
            continue
        for v, vd in graph.nodes(data=True):
            if vd[SYMBOL_KEY] != mapping_sym[1]:
                continue
            print("Map {} {} {}".format(u, v, mapping_sym))
            is_valid, mapping, (sub_vnodes, g_vnodes) = map_anchored_subgraph(
                graph, v, subgraph, u, mapper
            )
            if is_valid:
                possible_mappings.append((is_valid, mapping))
    return possible_mappings
