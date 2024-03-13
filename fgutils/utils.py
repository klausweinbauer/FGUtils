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
        4: ["C", "Si"],
        5: ["N", "P"],
        6: ["O", "S"],
        7: ["F", "Cl", "Br", "I"],
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
        assert (
            n_sym in valence_table.keys()
        ), "Element {} not found in valence table.".format(n_sym)
        bond_cnt = sum([b for _, _, b in graph.edges(n_id, data="bond")])  # type: ignore
        h_cnt = int(8 - valence_table[n_sym] - bond_cnt)
        assert h_cnt >= 0, "Negative hydrogen count."
        for h_id in range(len(graph), len(graph) + h_cnt):
            graph.add_node(h_id, symbol="H")
            graph.add_edge(n_id, h_id, bond=1)
    return graph
