import torch

from test.my_asserts import assert_edge_index, assert_edge_attr
from fgutils.parse import parse
from fgutils.torch.ITSDataset import ITSDataset


def test_ITSDataset():
    its = parse("C1<1,0>O<0,1>H<1,0>N<0,1>1")
    exp_x = torch.Tensor([6, 8, 1, 7]).unsqueeze(1)
    exp_edge_index = torch.Tensor([[0, 1, 2, 3, 1, 2, 3, 0], [1, 2, 3, 0, 0, 1, 2, 3]])
    exp_edge_attr = torch.Tensor([[1, 0], [0, 1], [1, 0], [0, 1]] * 2)

    ds = ITSDataset([its], [0])
    sample = ds.__getitem__(0)

    assert torch.equal(exp_x, sample.x)
    assert_edge_index(exp_edge_index, sample.edge_index)
    assert_edge_attr(exp_edge_attr, exp_edge_index, sample)
    assert 0 == sample.y


def test_ITSDataset_without_targets():
    its = parse("C1<1,0>O<0,1>H<1,0>N<0,1>1")
    exp_x = torch.Tensor([6, 8, 1, 7]).unsqueeze(1)
    exp_edge_index = torch.Tensor([[0, 1, 2, 3, 1, 2, 3, 0], [1, 2, 3, 0, 0, 1, 2, 3]])
    exp_edge_attr = torch.Tensor([[1, 0], [0, 1], [1, 0], [0, 1]] * 2)

    ds = ITSDataset([its])
    sample = ds.__getitem__(0)

    assert torch.equal(exp_x, sample.x)
    assert_edge_index(exp_edge_index, sample.edge_index)
    assert_edge_attr(exp_edge_attr, exp_edge_index, sample)
    assert sample.y is None
