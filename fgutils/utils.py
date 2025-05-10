import numpy as np
import networkx as nx

from fgutils.const import SYMBOL_KEY, BOND_KEY, AAM_KEY


def print_graph(graph):
    print(
        "Graph Nodes: {}".format(
            " ".join(
                [
                    "{}[{}]".format(n[1][SYMBOL_KEY], n[0])
                    for n in graph.nodes(data=True)
                ]
            )
        )
    )
    print(
        "Graph Edges: {}".format(
            " ".join(
                [
                    "[{}]-[{}]:{}".format(n[0], n[1], n[2][BOND_KEY])
                    for n in graph.edges(data=True)
                ]
            )
        )
    )


def add_implicit_hydrogens(graph: nx.Graph, inplace=True) -> nx.Graph:
    """Add implicit hydrogens as dedicated nodes to graph.

    :param graph: The mol graph to add hydrogens to.
    :param inplace: (optional) Inplace will change the input instance. When set
        to False the returned instance is a copy.

    :returns: Returns the molecular graph with implicit hydrogens or a copy if
        inplace is False.
    """
    if not inplace:
        graph = graph.copy()
    valence_dict = {
        2: ["Be", "Mg", "Ca", "Sr", "Ba"],
        3: ["B", "Al", "Ga", "In", "Tl"],
        4: ["C", "Si", "Sn", "Pb", "Pb"],
        5: ["N", "P", "As", "Sb", "Bi"],
        6: ["O", "S", "Se", "Te", "Po"],
        7: ["F", "Cl", "Br", "I", "At"],
    }
    valence_table = {}
    for v, elmts in valence_dict.items():
        for elmt in elmts:
            valence_table[elmt] = v
    nodes = [
        (n_id, n_sym)
        for n_id, n_sym in graph.nodes(data=SYMBOL_KEY)  # type: ignore
        if n_sym not in ["R", "H"]
    ]
    for n_id, n_sym in nodes:
        if n_sym not in valence_table.keys():
            # No hydrogens are added if element is not in dict. These atoms
            # are most likley not part of a functional group anyway so skipping
            # hydrogens is fine
            continue
        bond_cnt = sum([b for _, _, b in graph.edges(n_id, data=BOND_KEY)])  # type: ignore
        # h_cnt can be negative; aromaticity is complicated, we just ignore that
        valence = valence_table[n_sym]
        h_cnt = int(np.min([8, 2 * valence]) - valence - bond_cnt)
        for h_id in range(len(graph), len(graph) + h_cnt):
            node_attributes = {SYMBOL_KEY: "H"}
            edge_attributes = {BOND_KEY: 1}
            graph.add_node(h_id, **node_attributes)
            graph.add_edge(n_id, h_id, **edge_attributes)
    return graph


def remove_implicit_hydrogens(graph: nx.Graph, inplace=True) -> nx.Graph:
    """Remove implicit hydrogens as dedicated nodes to graph.

    :param graph: The mol graph to remove hydrogens to.
    :param inplace: (optional) Inplace will change the input instance. When set
        to False the returned instance is a copy.

    :returns: Returns the molecular graph with implicit hydrogens removed or a
        copy if inplace is False.
    """
    if not inplace:
        graph = graph.copy()
    nodes_to_remove = [n for n, d in graph.nodes(data=True) if d[SYMBOL_KEY] == "H"]
    graph.remove_nodes_from(
        nodes_to_remove,
    )
    return graph


# TODO Remove in future version. Deprecated since v0.1.6
def split_its(graph: nx.Graph) -> tuple[nx.Graph, nx.Graph]:
    """
    .. warning::
        Deprecated since v0.1.6. This function was moved to module fgutils.its.
        It will be removed from fgutils.utils in the future.

    Split an ITS graph into reactant graph G and product graph H. Required
    labels on the ITS graph are BOND_KEY.

    :param graph: ITS graph to split up.

    :returns: Tuple of two graphs (G, H).
    """
    print(
        "[WARNING] Function fgutils.split_its() is deprecated and will "
        + "be removed in a future version. Use function "
        + "fgutils.its.split_its() instead."
    )

    def _set_rc_edge(g, u, v, b):
        if b == 0:
            g.remove_edge(u, v)
        else:
            g[u][v][BOND_KEY] = b

    g = graph.copy()
    h = graph.copy()
    for u, v, d in graph.edges(data=True):  # type: ignore
        if d is None:
            raise ValueError("No edge labels found.")
        bond = d[BOND_KEY]
        if isinstance(bond, tuple) or isinstance(bond, list):
            _set_rc_edge(g, u, v, bond[0])
            _set_rc_edge(h, u, v, bond[1])
    return g, h


def initialize_aam(graph: nx.Graph, offset=1):
    """Initialize atom-atom map on a graph based on node indices.

    :param graph: The graph where to initialize the atom-atom map.

    :param offset: (optional) The mapping offset. Offset is the first value
        used for numbering.
    """
    for n, d in graph.nodes(data=True):
        if AAM_KEY in d:
            raise RuntimeError(
                "Graph has already an atom-atom map. "
                + "The original atom-atom map would be overwritten. "
                + "You might want to use fgutils.utils.complete_aam()."
            )
        d[AAM_KEY] = n + offset


def complete_aam(graph: nx.Graph, offset: None | int | str = None):
    """Complete the atom-atom map on a graph based on node indices. This
    function does not override an existing atom-atom map. It extends the
    existing atom-atom map to all nodes. The numbering of the new nodes starts
    at 1 or ``offset`` and skipps all existing mapping numbers.

    :param graph: The graph where to complete the atom-atom map.
    :param offset: (optional) The mapping offset. Offset is the first value
        used for numbering. If set to ``"min"`` the offset will be set to the
        lowest existing number.
    """
    mappings = [d[AAM_KEY] for _, d in graph.nodes(data=True) if AAM_KEY in d]
    next_mapping = 1
    if offset is not None:
        if isinstance(offset, int):
            next_mapping = offset
        elif offset == "min":
            if len(mappings) > 0:
                next_mapping = int(np.min(mappings))
        else:
            raise ValueError(
                (
                    "Unknown value '{}' for offset. " + 'Use integer or "min" instead.'
                ).format(offset)
            )
    for n, d in graph.nodes(data=True):
        if AAM_KEY in d or d[SYMBOL_KEY] == "H":
            continue
        while next_mapping in mappings:
            next_mapping += 1
        graph.nodes[n][AAM_KEY] = next_mapping
        mappings.append(next_mapping)


def mol_equal(
    candidate: nx.Graph,
    target: nx.Graph,
    compare_mode="both",
    iterations=3,
    ignore_hydrogens=False,
    min_atoms=0,
) -> bool:
    """Compare a candidate molecule to a specific target molecule. The
    result is true if the candidate matches the target. The
    comparison is done with 3 iterations WL.

    :param candidate: A candidate molecule.
    :param target: A target molecule to compare to.
    :param compare_mode: (optional) The compare mode: "target", "candidate",
        "largest_target", "largest_candidate", or "both". For "target" its
        sufficient that the candidate matches all components in the target, for
        "candidate" its sufficient that the target matches all candidate
        components and for "both" a complete match must be found. For the two
        largest_* compare modes only the respective largest component needs to
        find a match. This is helpful to ignore additional small compounds.
        (Default: both)
    :param iterations: (optional) The number of Weisfeiler-Leman iterations.
        (Default: 3)
    :param ignore_hydrogens: (optional) Flag to ignore hydrogens in the
        comparison. Molecules with unequal hydrogen configurations will count
        as identical. (Default: False)
    :param min_atoms: (optional) The minimum number of atoms a component needs
        to have to be considered in the comparison. This can be useful to
        exlude small molecules like water. (Default: 0)

    :returns: True if the candidate matches the target and else otherwise.
    """

    def _get_hash_list(g, iterations, largest_only=False, min_atoms=0):
        connected_node_sets = sorted(nx.connected_components(g), key=len, reverse=True)
        if largest_only:
            connected_node_sets = [connected_node_sets[0]]
        connected_components = [
            g.subgraph(c).copy() for c in connected_node_sets if len(c) >= min_atoms
        ]
        hash_list = [
            nx.weisfeiler_lehman_graph_hash(
                component,
                edge_attr=BOND_KEY,
                node_attr=SYMBOL_KEY,
                iterations=iterations,
            )
            for component in connected_components
        ]
        return hash_list

    if ignore_hydrogens:
        target = target.subgraph(
            [n for n, d in target.nodes(data=True) if d[SYMBOL_KEY] != "H"]
        ).copy()
        candidate = candidate.subgraph(
            [n for n, d in candidate.nodes(data=True) if d[SYMBOL_KEY] != "H"]
        ).copy()
    target_hash_list = _get_hash_list(
        target, iterations, compare_mode == "largest_target", min_atoms
    )
    candidate_hash_list = _get_hash_list(
        candidate, iterations, compare_mode == "largest_candidate", min_atoms
    )
    target_match = [h in candidate_hash_list for h in target_hash_list]
    candidate_match = [h in target_hash_list for h in candidate_hash_list]
    if "target" in compare_mode:
        return all(target_match)
    elif "candidate" in compare_mode:
        return all(candidate_match)
    elif compare_mode == "both":
        return all(target_match) and all(candidate_match)
    else:
        raise ValueError(
            "For 'compare_mode' use 'target', 'candidate' or 'both'. Value '{}' is unknown.".format(
                compare_mode
            )
        )


def get_unreachable_nodes(g, start_nodes, radius=1):
    """Get the list of nodes that can not be reached from start nodes within a
    given distance.

    :param g: The graph for which to get unreachable nodes.
    :param start_nodes: A list of nodes to start from. Start nodes count as
        radius 0. For a reachable node there must exist a path from a start
        node with at most radius number of steps.
    :param radius: The maximum number of hops from start_nodes. (Default: 1)

    :returns: Returns the list of unreachable nodes.
    """
    A = nx.adjacency_matrix(g, nodelist=range(len(g.nodes))).toarray()
    if radius == 0:
        D_sum = np.identity(A.shape[0])
    else:
        D = A.copy()
        D_sum = A.copy()
        for _ in range(radius - 1):
            D = np.matmul(D, A)
            D_sum += D
    center_paths = D_sum[start_nodes].sum(axis=0)
    return np.where(center_paths == 0)[0]


def relabel_graph(g, offset=0):
    mapping = {}
    for i, u in enumerate(sorted(g.nodes)):
        mapping[u] = i + offset
    return nx.relabel_nodes(g, mapping)


def to_non_aromatic_symbol(symbol):
    sym_map = {"c": "C", "n": "N", "b": "B", "o": "O", "p": "P", "s": "S"}
    return sym_map.get(symbol, symbol)
