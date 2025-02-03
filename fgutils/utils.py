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


def add_implicit_hydrogens(graph: nx.Graph) -> nx.Graph:
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
        if AAM_KEY in d:
            continue
        while next_mapping in mappings:
            next_mapping += 1
        graph.nodes[n][AAM_KEY] = next_mapping
        mappings.append(next_mapping)


def mol_compare(candidates: list[nx.Graph], target: nx.Graph) -> np.ndarray:
    """Compare a set of candidate molecules to a specific target molecule. The
    result is a binary vector of the same length as the candidates with the
    i-th entry beeing 1 if the candidate structurally matches the target. The
    comparison is done with 3 iterations WL.

    :param candidates: A list of candidate molecules.
    :param target: A target molecule to compare to.

    :returns: A numpy array of the same length as the candidates list with the
        i-th entry set to 1 if the candidate matches the target.
    """
    target_hash = nx.weisfeiler_lehman_graph_hash(
        target, edge_attr=BOND_KEY, node_attr=SYMBOL_KEY, iterations=3
    )
    result = np.zeros(len(candidates))
    for i, candidate in enumerate(candidates):
        candidate_hash = nx.weisfeiler_lehman_graph_hash(
            candidate, edge_attr=BOND_KEY, node_attr=SYMBOL_KEY, iterations=3
        )
        if candidate_hash == target_hash:
            result[i] = 1
    return result
