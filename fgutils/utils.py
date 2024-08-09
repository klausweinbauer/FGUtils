import numpy as np
import networkx as nx


def print_graph(graph):
    print(
        "Graph Nodes: {}".format(
            " ".join(
                ["{}[{}]".format(n[1]["symbol"], n[0]) for n in graph.nodes(data=True)]
            )
        )
    )
    print(
        "Graph Edges: {}".format(
            " ".join(
                [
                    "[{}]-[{}]:{}".format(n[0], n[1], n[2]["bond"])
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
        for n_id, n_sym in graph.nodes(data="symbol")  # type: ignore
        if n_sym not in ["R", "H"]
    ]
    for n_id, n_sym in nodes:
        if n_sym not in valence_table.keys():
            # No hydrogens are added if element is not in dict. These atoms
            # are most likley not part of a functional group anyway so skipping
            # hydrogens is fine
            continue
        bond_cnt = sum([b for _, _, b in graph.edges(n_id, data="bond")])  # type: ignore
        # h_cnt can be negative; aromaticity is complicated, we just ignore that
        valence = valence_table[n_sym]
        h_cnt = int(np.min([8, 2 * valence]) - valence - bond_cnt)
        for h_id in range(len(graph), len(graph) + h_cnt):
            graph.add_node(h_id, symbol="H")
            graph.add_edge(n_id, h_id, bond=1)
    return graph
