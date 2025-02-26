import torch
import networkx as nx

from typing import Callable
from fgutils.chem.ps import atomic_num2sym


def _build_its(node_attrs, edge_index, edge_attrs):
    assert edge_index.size(0) == 2
    assert edge_attrs.size(0) == edge_index.size(1)
    its = nx.Graph()
    for edge, eattr in zip(edge_index.T, edge_attrs):
        u = edge[0].item()
        v = edge[1].item()
        its.add_edge(u, v, bond=tuple(*eattr.tolist()))
        its.nodes[u]["symbol"] = node_attrs[u]
        its.nodes[v]["symbol"] = node_attrs[v]
    return its


def _default_node_feature_trans(x):
    return atomic_num2sym[int(x)]


def _its_from_torch_data(sample, node_feature_transform=None):
    if node_feature_transform is None:
        node_feature_transform = _default_node_feature_trans
    edge_attrs = torch.tensor(
        [[sample.edge_attr[i].tolist()] for i in range(sample.edge_attr.size(0))]
    )
    node_attrs = [node_feature_transform(sample.x[i]) for i in range(sample.x.size(0))]
    assert sample.edge_index.size(1) == len(edge_attrs)
    its = _build_its(node_attrs, sample.edge_index, edge_attrs)
    return its


def _its_from_torch_databatch(sample, node_feature_transform=None):
    if node_feature_transform is None:
        node_feature_transform = _default_node_feature_trans
    batch_size = sample.batch.max() + 1
    batch_indices = sample.batch.unique()
    assert (
        len(batch_indices) == batch_size
    ), "Batch size is not euqal to the number of indices? ({} != {})".format(
        batch_size, len(batch_indices)
    )
    graphs = []
    node_idx_offset = 0
    for batch_idx in batch_indices:
        node_indices = (sample.batch == batch_idx).nonzero().squeeze()
        edge_attrs = [
            [sample.edge_attr[i].tolist()] for i in range(sample.edge_attr.size(0))
        ]
        node_attrs = [node_feature_transform(sample.x[i]) for i in node_indices]
        assert sample.edge_index.size(1) == len(edge_attrs)
        smpl_edge_indices = []
        smpl_edge_attrs = []
        for edge, eattr in zip(sample.edge_index.T, edge_attrs):
            u = edge[0].item()
            if u not in node_indices:
                continue
            smpl_edge_indices.append(edge.tolist())
            smpl_edge_attrs.append(eattr)
        smpl_edge_indices = torch.tensor(smpl_edge_indices).T - node_idx_offset
        smpl_edge_attrs = torch.tensor(smpl_edge_attrs)
        its = _build_its(node_attrs, smpl_edge_indices, smpl_edge_attrs)
        graphs.append(its)
        node_idx_offset = node_indices.max() + 1
    return graphs


def its_from_torch(
    data, node_feature_transform: Callable[[torch.Tensor], str] | None = None
):
    """Convert a ITS from the PyTorch data format back into the NetworkX ITS
    format.

    :param data: The PyTorch data to convert. This can be a single sample or a
        batch of samples.
    :param node_feature_transform: (optional) A transform function to convert
        the node features back into the atom symbol. The function gets the node
        feature tensor as argument and is expected to return the atom symbol
        for that node.

    :returns: The sample ITS graph or a list of ITS graph if data is a batch.
    """
    if hasattr(data, "batch") and data.batch is not None:
        return _its_from_torch_databatch(
            data, node_feature_transform=node_feature_transform
        )
    else:
        return _its_from_torch_data(data, node_feature_transform=node_feature_transform)
