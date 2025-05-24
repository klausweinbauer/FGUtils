import torch
import networkx as nx

from typing import Callable
from torch_geometric.data import Batch

from fgutils.const import BOND_KEY, SYMBOL_KEY
from fgutils.chem.ps import atomic_num2sym, atomic_sym2num
from fgutils.its import ITS
from fgutils.utils import relabel_graph

from torch_geometric.data import Data


def _build_its(node_attrs, edge_index, edge_attrs):
    assert edge_index.size(0) == 2
    assert edge_attrs.size(0) == edge_index.size(1)
    its = nx.Graph()
    for i in range(len(node_attrs)):
        its.add_node(i, **{SYMBOL_KEY: node_attrs[i]})
    for edge, eattr in zip(edge_index.T, edge_attrs):
        u = edge[0].item()
        v = edge[1].item()
        assert u in its.nodes
        assert v in its.nodes
        its.add_edge(u, v, bond=tuple(eattr.tolist()))
    return its


def _default_node_feature_trans_torch2its(x):
    return atomic_num2sym[int(x[0])]


def _its_from_torch_data(sample, node_feature_transform=None):
    if node_feature_transform is None:
        node_feature_transform = _default_node_feature_trans_torch2its
    node_attrs = [node_feature_transform(sample.x[i]) for i in range(sample.x.size(0))]
    its = _build_its(node_attrs, sample.edge_index, sample.edge_attr)
    return its


def _its_from_torch_databatch(sample, node_feature_transform=None):
    if node_feature_transform is None:
        node_feature_transform = _default_node_feature_trans_torch2its
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
        node_attrs = [node_feature_transform(sample.x[i]) for i in node_indices]
        smpl_edge_indices = []
        smpl_edge_attrs = []
        for edge, eattr in zip(sample.edge_index.T, sample.edge_attr):
            u = edge[0].item()
            if u not in node_indices:
                continue
            smpl_edge_indices.append(edge.tolist())
            smpl_edge_attrs.append(eattr.tolist())
        smpl_edge_indices = torch.tensor(smpl_edge_indices).T - node_idx_offset
        smpl_edge_attrs = torch.tensor(smpl_edge_attrs)
        its = _build_its(node_attrs, smpl_edge_indices, smpl_edge_attrs)
        graphs.append(its)
        node_idx_offset = node_indices.max() + 1
    return graphs


def its_from_torch(
    data, node_feature_transform: Callable[[torch.Tensor], str] | None = None
):
    """Convert an ITS in PyTorch data format back into the NetworkX ITS format.

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


def _default_node_feature_trans_its2torch(d):
    return [atomic_sym2num[d[SYMBOL_KEY]]]


def _default_edge_feature_trans_its2torch(d):
    g_b, h_b = d[BOND_KEY]
    g_b = 0 if g_b is None else g_b
    h_b = 0 if h_b is None else h_b
    return [g_b, h_b]


def _its_to_torch(
    its: nx.Graph | ITS, node_feature_transform=None, edge_feature_transform=None
) -> Data:
    if isinstance(its, ITS):
        its = its.graph
    its = relabel_graph(its)

    if node_feature_transform is None:
        node_feature_transform = _default_node_feature_trans_its2torch
    if edge_feature_transform is None:
        edge_feature_transform = _default_edge_feature_trans_its2torch

    node_attrs = [None] * len(its.nodes)
    for u, d in its.nodes(data=True):
        node_attr = node_feature_transform(d)
        node_attrs[u] = node_attr
    assert None not in node_attrs
    x = torch.tensor(node_attrs)
    n_attr_cnt = len(node_attrs)

    edge_attrs = []
    edge_indices = []
    for u, v, d in its.edges(data=True):
        assert u < n_attr_cnt
        assert v < n_attr_cnt
        edge_indices.extend([[u, v], [v, u]])
        edge_attr = edge_feature_transform(d)
        edge_attrs.extend([edge_attr, edge_attr])
    edge_indices = torch.tensor(edge_indices).T
    edge_attrs = torch.tensor(edge_attrs)
    return Data(x=x, edge_index=edge_indices, edge_attr=edge_attrs)


def its_to_torch(
    its: nx.Graph | ITS | list[nx.Graph | ITS],
    node_feature_transform=None,
    edge_feature_transform=None,
) -> Data | Batch:
    """Convert an ITS graph or a list of ITS graphs to the PyTorch data format.

    :param its: The ITS graph to convert. If it is a single ITS graph the
        return type will be Data. For a list of ITS graphs a Batch is returned.
    :param node_feature_transform: (optional) A transform function to convert
        the node annotations into a node feature vector. The function gets the
        node dict as argument and is expected to return a torch.Tensor node
        feature vector.
    :param edge_feature_transform: (optional) A transform function to convert
        the edge annotations into an edge feature vector. The function gets the
        edge dict as argument and is expected to return a torch.Tensor edge
        feature vector.

    :returns: The sample ITS graph or a list of ITS graph if data is a batch.
    """
    if isinstance(its, list):
        torch_samples = []
        for _its in its:
            torch_samples.append(_its_to_torch(_its))
        return Batch.from_data_list(torch_samples)
    else:
        return _its_to_torch(
            its,
            node_feature_transform=node_feature_transform,
            edge_feature_transform=edge_feature_transform,
        )


def get_adjacency_matrix(sample: Data) -> torch.Tensor:
    """Get the adjacency matrix from a torch data sample. The sample needs
    property ``x`` and ``edge_index``.

    :param sample: The torch data sample.

    :returns: Returns the adjacency matrix as tensor.
    """
    if sample.x is None:
        raise ValueError("Sample is missing argument x (node features).")
    if sample.edge_index is None:
        raise ValueError("Sample is missing argument edge_index.")

    edges = sample.edge_index
    n = sample.x.size(0)
    A = torch.zeros((n, n))
    A[edges[0, :], edges[1, :]] = 1
    return A


def prune(sample: Data, start_nodes: torch.Tensor, radius=1) -> Data:
    """Prune the graph represented as torch data sample. This function removes
    all nodes that are ``radius`` farther away from a set of ``start_nodes``.

    :param sample: The torch data sample.
    :param start_nodex: A list of start node indices.
    :param radius: (optional) The maximum path length of nodes to keep around
        the start nodes. Every node with a longer shortes path to any start
        node is removed from the graph. Default: 1

    :returns: Returns the pruned graph as torch data sample.
    """
    if sample.x is None:
        raise ValueError("Sample is missing argument x (node features).")
    if sample.edge_index is None:
        raise ValueError("Sample is missing argument edge_index.")

    A = get_adjacency_matrix(sample)
    if radius == 0:
        D_sum = torch.eye(A.size(0))
    else:
        D = A.detach().clone()
        D_sum = A.detach().clone()
        for _ in range(radius - 1):
            D = torch.matmul(D, A)
            D_sum += D
    center_paths = D_sum[start_nodes].sum(axis=0)  # type: ignore
    reachable_nodes = torch.where(center_paths > 0)[0]
    node_map = {u: i for i, u in enumerate(reachable_nodes.tolist())}
    new_edge_index = []
    new_edge_attr = []
    edge_attrs = (
        sample.edge_attr.tolist()
        if sample.edge_attr is not None
        else [None] * sample.edge_index.size(1)
    )
    for (u, v), attr in zip(sample.edge_index.T.tolist(), edge_attrs):
        if u in reachable_nodes and v in reachable_nodes:
            new_edge_index.append([node_map[u], node_map[v]])
            if attr is not None:
                new_edge_attr.append(attr)
    new_edge_index = torch.tensor(new_edge_index).T
    new_edge_attr = torch.tensor(new_edge_attr)
    pruned_sample = Data(
        edge_index=new_edge_index,
        x=sample.x[reachable_nodes],
        y=sample.y,
    )
    if sample.edge_attr is not None:
        pruned_sample.edge_attr = new_edge_attr
    if "id" in sample:
        pruned_sample.id = sample.id
    return pruned_sample


def prune_rc(sample, radius=1):
    """Prune the torch data ITS graph represented. This function removes
    all nodes that are ``radius`` farther away from the reaction center.

    :param sample: The torch data sample.
    :param radius: (optional) The maximum path length of nodes to keep around
        the reaction center. Every node with a longer shortes path to any node
        in the reaction center is removed from the graph. Default: 1

    :returns: Returns the pruned graph as torch data sample.
    """
    if sample.edge_attr.size(-1) != 2:
        raise NotImplementedError(
            "Edge attributes should be explicit, i.e., "
            + "only the bond change attributes. "
            + "Found tensor of size {} instead.".format(sample.edge_attr.size())
        )
    rc_edge_idx = (sample.edge_attr[:, 0] != sample.edge_attr[:, 1]).nonzero().squeeze()
    rc_node_idx = sample.edge_index[0, rc_edge_idx].unique()
    return prune(sample, rc_node_idx, radius=radius)
