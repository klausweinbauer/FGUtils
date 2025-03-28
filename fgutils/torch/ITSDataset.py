import torch
import networkx as nx
from fgutils.its import ITS

from fgutils.torch.utils import its_to_torch


class ITSDataset(torch.utils.data.Dataset):
    """Instantiates a torch.Dataset from a list of ITS graphs.

    :param data: The list of ITS graphs to create the dataset from.
    :param y: (optional) The list of ITS labels. These are the prediction
        targets. This list must be of the same length as the ITS graphs.
    :param ids: (optional) A list of IDs. The number of IDs must be equal to
        the number of graphs. Default: None
    :param node_feature_transform: (optional) A transform function to convert
        node attributes to a node feature vector. Default: None
    :param edge_feature_transform: (optional) A transform function to convert
        edge attributes to an edge feature vector. Default: None
    :param pre_transform: (optional) A callback function to transform the
        Tensor representation of the ITS graph before it is stored. This
        function is executed only once on initialization.
    :param transform: (optional) A callback function to transform the
        Tensor representation of the ITS graph when accessed. This function is
        executed on every access.
    """

    def __init__(
        self,
        data: list[nx.Graph | ITS],
        y: list[int] | None = None,
        ids: list[int] | None = None,
        node_feature_transform=None,
        edge_feature_transform=None,
        pre_transform=None,
        transform=None,
    ):
        super(ITSDataset, self).__init__()
        self.has_ids = True
        if y is None:
            self.has_ids = False
            y = [0] * len(data)
        if len(data) != len(y):
            raise ValueError(
                "Number of ITS graphs must be equal to the number of labels. ({} != {})".format(
                    len(data), len(y)
                )
            )
        if ids is not None and len(data) != len(ids):
            raise ValueError(
                "Number of IDs must be equal to the number of ITS graphs. ({} != {})".format(
                    len(ids), len(data)
                )
            )

        torch_data = []
        for its, y, id in zip(data, y, ids if ids is not None else [None] * len(y)):
            its_torch = its_to_torch(
                its, node_feature_transform, edge_feature_transform
            )
            if self.has_ids:
                its_torch.y = torch.tensor(y)
            if id is not None:
                its_torch.id = torch.tensor(id)
            if pre_transform is not None:
                its_torch = pre_transform(its_torch)
            torch_data.append(its_torch)

        self.data = torch_data
        self.pre_transform = pre_transform
        self.transform = transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        sample = self.data[index]
        if self.transform is not None:
            sample = self.transform(sample)
        return sample
