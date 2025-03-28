import torch

from torch_geometric.data import Data


def node_induced_subgraph(graph, nodes):
    """Get a node induced subgraph form the torch graph representation. Keep
    in mind that the returned node indices (e.g. in the edge_index tensor) are
    not the same as in the input graph. This is because the node indices are
    implicit in the row number of the node features and hence a relabeling is
    necessary.

    :param graph: The graph from which to get a node induced subgraph.
    :param nodes: A list of node indices for the induced subgraph.

    :returns: Node induced subgraph.
    """
    if not isinstance(nodes, torch.Tensor):
        nodes = torch.tensor(nodes)
    node_map = {n: i for i, n in enumerate(nodes.tolist())}
    new_edge_index = []
    selected_edge_indices = []
    for i, (u, v) in enumerate(graph.edge_index.T.tolist()):
        if u in nodes and v in nodes:
            selected_edge_indices.append(i)
            new_edge_index.append([node_map[u], node_map[v]])
    selected_edge_indices = torch.tensor(selected_edge_indices)
    new_x = graph.x[nodes].detach().clone()
    new_edge_index = torch.tensor(new_edge_index).T
    subgraph = Data(x=new_x, edge_index=new_edge_index)
    if graph.edge_attr is not None:
        subgraph.edge_attr = graph.edge_attr[selected_edge_indices].detach().clone()
    return subgraph


def edge_induced_subgraph(graph, edges):
    """Get an edge induced subgraph form the torch graph representation. Keep
    in mind that the returned node indices (e.g. in the edge_index tensor) are
    not the same as in the input graph. This is because the node indices are
    implicit in the row number of the node features and hence a relabeling is
    necessary.

    :param graph: The graph from which to get a edge induced subgraph.
    :param edges: A list of edge indices for the induced subgraph.

    :returns: Node induced subgraph.
    """
    if not isinstance(edges, torch.Tensor):
        edges = torch.tensor(edges)
    new_edge_index = graph.edge_index.T[edges].detach().clone().T
    selected_node_indices = new_edge_index.unique()
    node_map = {n: i for i, n in enumerate(selected_node_indices.tolist())}
    new_edge_index.apply_(lambda x: node_map[x])
    new_x = graph.x[selected_node_indices].detach().clone()
    subgraph = Data(x=new_x, edge_index=new_edge_index)
    if graph.edge_attr is not None:
        subgraph.edge_attr = graph.edge_attr[edges].detach().clone()
    return subgraph
