import torch
import torch_geometric.data

from test.my_asserts import assert_graph_eq, assert_edge_index, assert_edge_attr
from fgutils.const import (
    LABELS_KEY,
    IS_LABELED_KEY,
    AAM_KEY,
    IDX_MAP_KEY,
)
from fgutils.parse import parse as pattern_to_graph
from fgutils.its import ITS

from fgutils.torch.utils import (
    _build_its,
    its_from_torch,
    its_to_torch,
    get_adjacency_matrix,
    prune,
)


def get_torch_sample(x, edge_index, edge_attr=None, y=0):
    x = torch.tensor(x)
    if len(x.size()) == 1:
        x = x.unsqueeze(1)
    assert x.size(1) == 1, "x.size(1) must be 1 but got {}".format(x.size())
    bidirectional_edge_index = [
        [*edge_index[0], *edge_index[1]],
        [*edge_index[1], *edge_index[0]],
    ]
    edge_index = torch.tensor(bidirectional_edge_index)
    assert edge_index.size(0) == 2, "edge_index.size(0) must be 2 but got {}".format(
        edge_index.size()
    )
    if edge_attr is not None:
        bidirectional_edge_attr = [*edge_attr, *edge_attr]
        edge_attr = torch.tensor(bidirectional_edge_attr)
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


def test_its_to_torch_data():
    exp_x = torch.tensor([6, 6, 6, 6]).unsqueeze(1)
    exp_edge_index = [[0, 1, 2, 3], [1, 2, 3, 0]]
    exp_edge_index = torch.tensor(
        [
            [*exp_edge_index[0], *exp_edge_index[1]],
            [*exp_edge_index[1], *exp_edge_index[0]],
        ]
    )
    exp_edge_attr = torch.tensor([[0, 1], [1, 0], [0, 1], [1, 0]] * 2)

    its = pattern_to_graph("C1<0,1>C<1,0>C<0,1>C<1,0>1")
    sample = its_to_torch(its)

    assert torch.equal(exp_x, sample.x)
    assert_edge_index(exp_edge_index, sample.edge_index)
    assert_edge_attr(exp_edge_attr, exp_edge_index, sample)


def test_its2torch2its():
    G1 = pattern_to_graph("C1<0,1>C<1,0>C<0,1>C<1,0>1")
    G2 = pattern_to_graph("CC<2,1>O")
    batch = its_to_torch([G1, G2])
    its_list = its_from_torch(batch)
    assert_graph_eq(its_list[0], G1, ignore_keys=[LABELS_KEY, IS_LABELED_KEY])
    assert_graph_eq(its_list[1], G2, ignore_keys=[LABELS_KEY, IS_LABELED_KEY])


def test_its2torch2its_2():
    smiles = (
        "[OH:1][c:2]1[cH:3][cH:4][cH:5][c:6]([OH:7])[cH:8]1."
        + "[c:9]1([NH2:19])[cH:10][cH:11][cH:12][c:13]2[cH:14][cH:15][cH:16][cH:17][c:18]12>>"
        + "[OH2:1].[c:2]1([NH:19][c:9]2[cH:10][cH:11][cH:12][c:13]3[cH:14][cH:15][cH:16]"
        + "[cH:17][c:18]23)[cH:3][cH:4][cH:5][c:6]([OH:7])[cH:8]1"
    )
    in_its = ITS.from_smiles(smiles).graph
    torch_its = its_to_torch(in_its)
    out_its = its_from_torch(torch_its)
    assert_graph_eq(in_its, out_its, ignore_keys=[AAM_KEY, IDX_MAP_KEY])


def test_get_adjacency_matrix():
    x = [6, 6, 6, 6, 6]
    edge_index = [[0, 1, 2, 3, 2], [1, 2, 3, 4, 0]]
    edge_attr = [[0, 1], [1, 1], [2, 2], [1, 1], [1, 1]]
    sample = get_torch_sample(x, edge_index, edge_attr)
    exp_A = torch.zeros((len(x), len(x)))
    exp_A[edge_index] = 1
    exp_A[torch.tensor(edge_index)[torch.LongTensor([1, 0])].tolist()] = 1
    A = get_adjacency_matrix(sample)
    assert torch.equal(exp_A, A)


def test_prune():
    x = [6, 6, 6, 6, 6]
    edge_index = [[0, 1, 2, 3, 2], [1, 2, 3, 4, 0]]
    edge_attr = [[0, 1], [1, 1], [2, 2], [1, 1], [1, 1]]
    sample = get_torch_sample(x, edge_index, edge_attr)

    exp_x = torch.tensor([6, 6, 6, 6]).unsqueeze(1)
    exp_edge_index = [[0, 1, 2, 2], [1, 2, 3, 0]]
    exp_edge_index = torch.tensor(
        [
            [*exp_edge_index[0], *exp_edge_index[1]],
            [*exp_edge_index[1], *exp_edge_index[0]],
        ]
    )
    exp_edge_attr = torch.tensor(
        [[0, 1], [1, 1], [2, 2], [1, 1], [0, 1], [1, 1], [2, 2], [1, 1]]
    )

    pruned_sample = prune(sample, torch.tensor([1]), radius=2)
    assert torch.equal(exp_x, pruned_sample.x)
    assert torch.equal(exp_edge_index, pruned_sample.edge_index)
    assert torch.equal(exp_edge_attr, pruned_sample.edge_attr)


def test_prune_without_edge_attr():
    x = [6, 6, 6, 6, 6]
    edge_index = [[0, 1, 2, 3, 2], [1, 2, 3, 4, 0]]
    sample = get_torch_sample(x, edge_index)

    exp_x = torch.tensor([6, 6, 6, 6]).unsqueeze(1)
    exp_edge_index = [[0, 1, 2, 2], [1, 2, 3, 0]]
    exp_edge_index = torch.tensor(
        [
            [*exp_edge_index[0], *exp_edge_index[1]],
            [*exp_edge_index[1], *exp_edge_index[0]],
        ]
    )

    pruned_sample = prune(sample, torch.tensor([1]), radius=2)
    assert torch.equal(exp_x, pruned_sample.x)
    assert torch.equal(exp_edge_index, pruned_sample.edge_index)
    assert pruned_sample.edge_attr is None


def test_build_its_with_unordered_edge_index():
    x = [1, 2, 3, 4]
    edge_index = [[0, 0, 1], [1, 3, 2]]
    edge_attr = [[1, 1], [1, 1], [1, 1]]
    sample = get_torch_sample(x, edge_index, edge_attr)
    its = _build_its(sample.x, sample.edge_index, sample.edge_attr)
    new_x = [d["symbol"].item() for _, d in its.nodes(data=True)]
    assert all([x[i] == new_x[i] for i in range(len(x))])
