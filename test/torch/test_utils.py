import torch
import torch_geometric.data

from test.my_asserts import assert_graph_eq
from fgutils.const import LABELS_KEY, IS_LABELED_KEY
from fgutils.parse import parse as pattern_to_graph

from fgutils.torch.utils import its_from_torch


def get_torch_sample(x, edge_index, edge_attr=None, y=0):
    x = torch.tensor(x)
    if len(x.size()) == 1:
        x = x.unsqueeze(1)
    assert x.size(1) == 1, "x.size(1) must be 1 but got {}".format(x.size())
    edge_index = torch.tensor(edge_index)
    assert edge_index.size(0) == 2, "edge_index.size(0) must be 2 but got {}".format(
        edge_index.size()
    )
    if edge_attr is not None:
        edge_attr = torch.tensor(edge_attr)
        assert edge_attr.size(0) == edge_index.size(
            1
        ), "edge_attr.size(0) must be equal to edge_index.size(1) but got {} != {}".format(
            edge_attr.size(), edge_index.size()
        )
    return torch_geometric.data.Data(
        x=x, edge_index=edge_index, edge_attr=edge_attr, y=y
    )


def test_its_from_torch_data():
    x = [6, 6, 6, 6]
    edge_index = [[0, 1, 2, 3], [1, 2, 3, 0]]
    edge_attr = [[0, 1], [1, 0], [0, 1], [1, 0]]
    sample = get_torch_sample(x, edge_index, edge_attr)
    G = its_from_torch(sample)
    exp_G = pattern_to_graph("C1<0,1>C<1,0>C<0,1>C<1,0>1")
    assert_graph_eq(G, exp_G, ignore_keys=[LABELS_KEY, IS_LABELED_KEY])


def test_its_from_torch_databatch():
    x1 = [6, 6, 6, 6]
    edge_index1 = [[0, 1, 2, 3], [1, 2, 3, 0]]
    edge_attr1 = [[0, 1], [1, 0], [0, 1], [1, 0]]
    sample1 = get_torch_sample(x1, edge_index1, edge_attr1)
    x2 = [6, 6, 8]
    edge_index2 = [[0, 1], [1, 2]]
    edge_attr2 = [[1, 1], [2, 1]]
    sample2 = get_torch_sample(x2, edge_index2, edge_attr2)
    batch = torch_geometric.data.Batch.from_data_list([sample1, sample2])
    graphs = its_from_torch(batch)
    exp_G1 = pattern_to_graph("C1<0,1>C<1,0>C<0,1>C<1,0>1")
    exp_G2 = pattern_to_graph("CC<2,1>O")
    assert len(graphs) == 2
    assert_graph_eq(graphs[0], exp_G1, ignore_keys=[LABELS_KEY, IS_LABELED_KEY])
    assert_graph_eq(graphs[1], exp_G2, ignore_keys=[LABELS_KEY, IS_LABELED_KEY])
