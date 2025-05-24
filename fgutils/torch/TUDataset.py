# Copied and modified from
# https://github.com/pyg-team/pytorch_geometric/blob/master/torch_geometric/datasets/tu_dataset.py

import os
import sys
import torch
import os.path as osp
import logging

from typing import Dict, Tuple
from torch import Tensor
from torch_geometric.data import Data, InMemoryDataset
from torch_geometric.data.dataset import _repr, files_exist
from torch_geometric.io import fs
from torch_geometric.io.tu import read_file, cat
from torch_geometric.utils import coalesce, cumsum, one_hot, remove_self_loops

names = [
    "A",
    "graph_indicator",
    "node_labels",
    "node_attributes" "edge_labels",
    "edge_attributes",
    "graph_labels",
    "graph_attributes",
    "graph_ids",
]


def split(data: Data, batch: Tensor) -> Tuple[Data, Dict[str, Tensor]]:
    node_slice = cumsum(torch.bincount(batch))

    assert data.edge_index is not None
    row, _ = data.edge_index
    edge_slice = cumsum(torch.bincount(batch[row]))

    # Edge indices should start at zero for every graph.
    data.edge_index -= node_slice[batch[row]].unsqueeze(0)

    slices = {"edge_index": edge_slice}
    if data.x is not None:
        slices["x"] = node_slice
    else:
        # Imitate `collate` functionality:
        data._num_nodes = torch.bincount(batch).tolist()
        data.num_nodes = batch.numel()
    if data.edge_attr is not None:
        slices["edge_attr"] = edge_slice
    if data.y is not None:
        assert isinstance(data.y, Tensor)
        if data.y.size(0) == batch.size(0):
            slices["y"] = node_slice
        else:
            slices["y"] = torch.arange(0, int(batch[-1]) + 2, dtype=torch.long)
    if data.id is not None:
        assert isinstance(data.id, Tensor)
        assert data.id.size(0) != batch.size(0)
        slices["id"] = torch.arange(0, int(batch[-1]) + 2, dtype=torch.long)

    return data, slices


def read_tu_data_with_ids(
    folder: str,
    prefix: str,
) -> Tuple[Data, Dict[str, Tensor], Dict[str, int]]:
    files = fs.glob(osp.join(folder, f"{prefix}_*.txt"))
    prefix_i = len(prefix) + 1
    names = [osp.basename(f)[prefix_i:-4] for f in files]

    edge_index = read_file(folder, prefix, "A", torch.long).t() - 1
    batch = read_file(folder, prefix, "graph_indicator", torch.long) - 1

    node_attribute = torch.empty((batch.size(0), 0))
    if "node_attributes" in names:
        node_attribute = read_file(folder, prefix, "node_attributes")
        if node_attribute.dim() == 1:
            node_attribute = node_attribute.unsqueeze(-1)

    node_label = torch.empty((batch.size(0), 0))
    if "node_labels" in names:
        node_label = read_file(folder, prefix, "node_labels", torch.long)
        if node_label.dim() == 1:
            node_label = node_label.unsqueeze(-1)
        node_label = node_label - node_label.min(dim=0)[0]
        node_labels = list(node_label.unbind(dim=-1))
        node_labels = [one_hot(x) for x in node_labels]
        if len(node_labels) == 1:
            node_label = node_labels[0]
        else:
            node_label = torch.cat(node_labels, dim=-1)

    edge_attribute = torch.empty((edge_index.size(1), 0))
    if "edge_attributes" in names:
        edge_attribute = read_file(folder, prefix, "edge_attributes")
        if edge_attribute.dim() == 1:
            edge_attribute = edge_attribute.unsqueeze(-1)

    edge_label = torch.empty((edge_index.size(1), 0))
    if "edge_labels" in names:
        edge_label = read_file(folder, prefix, "edge_labels", torch.long)
        if edge_label.dim() == 1:
            edge_label = edge_label.unsqueeze(-1)
        edge_label = edge_label - edge_label.min(dim=0)[0]
        edge_labels = list(edge_label.unbind(dim=-1))
        edge_labels = [one_hot(e) for e in edge_labels]
        if len(edge_labels) == 1:
            edge_label = edge_labels[0]
        else:
            edge_label = torch.cat(edge_labels, dim=-1)

    x = cat([node_attribute, node_label])
    edge_attr = cat([edge_attribute, edge_label])

    y = None
    if "graph_attributes" in names:  # Regression problem.
        y = read_file(folder, prefix, "graph_attributes")
    elif "graph_labels" in names:  # Classification problem.
        y = read_file(folder, prefix, "graph_labels", torch.long)
        _, y = y.unique(sorted=True, return_inverse=True)

    id = None
    if "ids" in names:
        id = read_file(folder, prefix, "ids", torch.int64)

    num_nodes = int(edge_index.max()) + 1 if x is None else x.size(0)
    edge_index, edge_attr = remove_self_loops(edge_index, edge_attr)
    edge_index, edge_attr = coalesce(edge_index, edge_attr, num_nodes)

    if id is not None:
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y, id=id)
    else:
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y)
    data, slices = split(data, batch)

    sizes = {
        "num_node_attributes": node_attribute.size(-1),
        "num_node_labels": node_label.size(-1),
        "num_edge_attributes": edge_attribute.size(-1),
        "num_edge_labels": edge_label.size(-1),
    }

    return data, slices, sizes


class TUDataset(InMemoryDataset):
    """Class for reading datasets in the
    [TUDataset](https://chrsmrrs.github.io/datasets/docs/format/) file format.
    This implementation is copied and modified from
    [torch_geometric.datasets.TUDataset](https://pytorch-geometric.readthedocs.io/en/latest/generated/torch_geometric.datasets.TUDataset.html)
    to also support instance IDs for easier identification of samples. IDs are
    loaded from an additional ids.txt file where the i-th line is the id of the
    i-th graph.

    :param root: Root directory where the dataset is saved.
    :param name: The name of the dataset.
    :param transform: A function/transform that takes in a data object and
        returns a transformed version. The data object will be transformed
        before every access. (default: None)
    :param pre_transform: A function/transform that takes in
        a data object and returns a transformed version. The data object
        will be transformed before being saved to disk. (default: None)
    :param pre_filter: A function that takes in an data object and returns
        a boolean value, indicating whether the data object should be
        included in the final dataset. (default: None)
    :param force_reload: Whether to re-process the dataset. (default:
        False)
    :param use_node_attr: If True, the dataset will contain additional
        continuous node attributes (if present). (default: False)
    :param use_edge_attr: If True, the dataset will contain additional
        continuous edge attributes (if present). (default: False)
    """

    def __init__(
        self,
        root,
        name,
        transform=None,
        pre_transform=None,
        pre_filter=None,
        use_node_attr: bool = True,
        use_edge_attr: bool = True,
    ):
        self.ds_name = name
        super().__init__(root, transform, pre_transform, pre_filter)

        out = fs.torch_load(self.processed_paths[0])
        data, self.slices, self.sizes, data_cls = out
        self.data = data_cls.from_dict(data)

        assert isinstance(self._data, Data)
        if self._data.x is not None and not use_node_attr:
            num_node_attributes = self.num_node_attributes
            self._data.x = self._data.x[:, num_node_attributes:]
        if self._data.edge_attr is not None and not use_edge_attr:
            num_edge_attrs = self.num_edge_attributes
            self._data.edge_attr = self._data.edge_attr[:, num_edge_attrs:]

    @property
    def raw_dir(self) -> str:
        return osp.join(self.root, self.ds_name)

    @property
    def processed_dir(self) -> str:
        return osp.join(self.root, self.ds_name)

    @property
    def raw_file_names(self):
        return [
            "{}_{}".format(self.ds_name, n)
            for n in [
                "A",
                "edge_attributes",
                "node_attributes",
                "graph_indicator",
                "graph_labels",
                "ids",
            ]
        ]

    @property
    def processed_file_names(self):
        return ["data.pt", "pre_transform.pt", "pre_filter.pt"]

    def _process(self):
        f = osp.join(self.processed_dir, "pre_transform.pt")
        if osp.exists(f) and torch.load(f, weights_only=False) != _repr(
            self.pre_transform
        ):
            logging.info("Dataset pre_transform changed.")
            os.remove(f)

        f = osp.join(self.processed_dir, "pre_filter.pt")
        if osp.exists(f) and torch.load(f, weights_only=False) != _repr(
            self.pre_filter
        ):
            logging.info("Dataset pre_filter changed.")
            os.remove(f)

        if files_exist(self.processed_paths):
            return

        if self.log and "pytest" not in sys.modules:
            logging.info("Preprocessing dataset.")

        fs.makedirs(self.processed_dir, exist_ok=True)
        self.process()

        path = osp.join(self.processed_dir, "pre_transform.pt")
        fs.torch_save(_repr(self.pre_transform), path)
        path = osp.join(self.processed_dir, "pre_filter.pt")
        fs.torch_save(_repr(self.pre_filter), path)

    def download(self):
        pass
        # Download to `self.raw_dir`.
        # raise NotImplementedError()
        # download_url(url, self.raw_dir)
        # ...

    def process(self) -> None:
        self.data, self.slices, sizes = read_tu_data_with_ids(
            self.raw_dir, self.ds_name
        )

        if self.pre_filter is not None or self.pre_transform is not None:
            data_list = [self.get(idx) for idx in range(len(self))]

            if self.pre_filter is not None:
                data_list = [d for d in data_list if self.pre_filter(d)]

            if self.pre_transform is not None:
                data_list = [self.pre_transform(d) for d in data_list]

            self.data, self.slices = self.collate(data_list)
            self._data_list = None  # Reset cache.

        assert isinstance(self._data, Data)
        fs.torch_save(
            (self._data.to_dict(), self.slices, sizes, self._data.__class__),
            self.processed_paths[0],
        )
