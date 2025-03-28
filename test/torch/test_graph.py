import torch

from test.torch.test_utils import get_torch_sample

from fgutils.torch.graph import node_induced_subgraph, edge_induced_subgraph


def test_node_induced_subgraph():
    x = [1, 2, 3, 4]
    edge_index = [[0, 0, 0, 1, 1, 2], [1, 2, 3, 2, 3, 3]]
    graph = get_torch_sample(x, edge_index)
    subgraph = node_induced_subgraph(graph, [0, 2, 3])
    assert torch.equal(subgraph.x, torch.tensor([1, 3, 4]).unsqueeze(1))
    assert torch.equal(
        subgraph.edge_index, torch.tensor([[0, 0, 1, 1, 2, 2], [1, 2, 2, 0, 0, 1]])
    )


def test_node_induced_subgraph_with_edge_attr():
    x = [1, 2, 3, 4]
    edge_index = [[0, 0, 0, 1, 1, 2], [1, 2, 3, 2, 3, 3]]
    edge_attr = [[1], [2], [3], [4], [5], [6]]
    graph = get_torch_sample(x, edge_index, edge_attr)
    subgraph = node_induced_subgraph(graph, [0, 2, 3])
    assert torch.equal(subgraph.x, torch.tensor([1, 3, 4]).unsqueeze(1))
    assert torch.equal(
        subgraph.edge_index, torch.tensor([[0, 0, 1, 1, 2, 2], [1, 2, 2, 0, 0, 1]])
    )
    assert torch.equal(subgraph.edge_attr, torch.tensor([[2], [3], [6], [2], [3], [6]]))


def test_edge_induced_subgraph():
    x = [1, 2, 3, 4]
    edge_index = [[0, 0, 0, 1, 1, 2], [1, 2, 3, 2, 3, 3]]
    graph = get_torch_sample(x, edge_index)
    subgraph = edge_induced_subgraph(graph, [0, 6, 3])
    assert torch.equal(subgraph.x, torch.tensor([1, 2, 3]).unsqueeze(1))
    assert torch.equal(subgraph.edge_index, torch.tensor([[0, 1, 1], [1, 0, 2]]))


def test_edge_induced_subgraph_node_relabeling():
    x = [1, 2, 3, 4]
    edge_index = [[0, 0, 0, 1, 1, 2], [1, 2, 3, 2, 3, 3]]
    graph = get_torch_sample(x, edge_index)
    subgraph = edge_induced_subgraph(graph, [2, 8])
    assert torch.equal(subgraph.x, torch.tensor([1, 4]).unsqueeze(1))
    assert torch.equal(subgraph.edge_index, torch.tensor([[0, 1], [1, 0]]))


def test_edge_induced_subgraph_unsorted_edges():
    x = [1, 2, 3, 4]
    edge_index = [[2, 0, 1, 0, 1, 0], [3, 2, 2, 3, 3, 1]]
    graph = get_torch_sample(x, edge_index)
    subgraph = edge_induced_subgraph(graph, [5, 2, 11])
    assert torch.equal(subgraph.x, torch.tensor([1, 2, 3]).unsqueeze(1))
    assert torch.equal(subgraph.edge_index, torch.tensor([[0, 1, 1], [1, 2, 0]]))


def test_node_induced_subgraph_with_edge_attr2():
    x = [1, 2, 3, 4]
    edge_index = [[0, 0, 0, 1, 1, 2], [1, 2, 3, 2, 3, 3]]
    edge_attr = [[1], [2], [3], [4], [5], [6]]
    graph = get_torch_sample(x, edge_index, edge_attr)
    subgraph = edge_induced_subgraph(graph, [0, 6, 3])
    assert torch.equal(subgraph.x, torch.tensor([1, 2, 3]).unsqueeze(1))
    assert torch.equal(subgraph.edge_index, torch.tensor([[0, 1, 1], [1, 0, 2]]))
    assert torch.equal(subgraph.edge_attr, torch.tensor([[1], [1], [4]]))
