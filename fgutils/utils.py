import numpy as np
import networkx as nx

from fgutils.const import SYMBOL_KEY, BOND_KEY


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


def split_its(graph: nx.Graph) -> tuple[nx.Graph, nx.Graph]:
    """Split an ITS graph into reactant graph G and product graph H.

    :param graph: ITS graph to split up.

    :returns: Tuple of two graphs (G, H).
    """

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
        if isinstance(bond, tuple):
            _set_rc_edge(g, u, v, bond[0])
            _set_rc_edge(h, u, v, bond[1])
    return g, h
