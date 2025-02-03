import os
import numpy as np

from fgutils.const import SYMBOL_KEY, BOND_KEY
from fgutils.chem.ps import get_atomic_number
from fgutils.its import ITS


class RawTUDataset:
    """Helper class for TUDataset construction.

    :param name: The name of the new dataset.
    :param A: The adjacency matrix list.
    :param graph_indicator: Graph indicator list.
    :param node_attributes: Node attribute list.
    :param edge_attributes: Edge attribute list.
    :param ids: Id list.
    """

    def __init__(
        self,
        name,
        A: list,
        graph_indicator: list,
        graph_labels: None | list = None,
        node_attributes: None | list = None,
        edge_attributes: None | list = None,
        ids: None | list = None,
    ):
        self.name = name
        self.A = A
        self.ids = ids
        self.graph_indicator = graph_indicator
        self.graph_labels = graph_labels
        self.node_attributes = node_attributes
        self.edge_attributes = edge_attributes

    @staticmethod
    def from_reaction_dict(
        data: list[dict], id_col: str, its_col: str, label_col: str, dataset_name: str
    ):
        """Create a RawTUDataset instance from a reaction dict. Each entry in
        the data list contains one reaction. Each entry (dict) needs to provide
        a reaction id, an ITS graph, and a label. The label is the reaction
        class for example. In the ML context the label is the
        classification/regression target.

        :param data: List of dictionaries. Each dict requires the keys:
            ``id_col``, ``its_col``, ``label_col``
        :param id_col: The dict key where to find the reaction id.
        :param its_col: The dict key where to find the ITS graph.
        :param label_col: The dict key where to find the graph label.

        :returns: A instance of RawTUDataset
        """
        DS_A = []
        DS_ids = []
        DS_graph_indicator = []
        DS_graph_labels = []
        DS_node_attributes = []
        DS_edge_attributes = []

        node_idx_offset = 1

        node_cnt = 0
        edge_cnt = 0
        graph_cnt = len(data)

        for i, entry in enumerate(data):
            its = entry[its_col]
            if isinstance(its, ITS):
                its = its.graph
            node_cnt += len(its.nodes)
            edge_cnt += 2 * len(its.edges)
            node_idx_offset += len(its.nodes)
            DS_graph_labels.append(entry[label_col])
            DS_ids.append(entry[id_col])
            for u, d in its.nodes(data=True):
                DS_graph_indicator.append(i + 1)
                DS_node_attributes.append(get_atomic_number(d[SYMBOL_KEY]))
            for u, v, d in its.edges(data=True):
                DS_A.extend([(u, v), (v, u)])
                g_b, h_b = d[BOND_KEY]
                g_b = 0 if g_b is None else g_b
                h_b = 0 if h_b is None else h_b
                DS_edge_attributes.extend([(g_b, h_b), (g_b, h_b)])
        assert edge_cnt == len(DS_A)
        assert node_cnt == len(DS_graph_indicator)
        assert graph_cnt == len(DS_graph_labels)
        assert node_cnt == len(DS_node_attributes)
        assert edge_cnt == len(DS_edge_attributes)
        assert graph_cnt == len(DS_ids)

        return RawTUDataset(
            dataset_name,
            DS_A,
            DS_graph_indicator,
            graph_labels=DS_graph_labels,
            node_attributes=DS_node_attributes,
            edge_attributes=DS_edge_attributes,
            ids=DS_ids,
        )

    @property
    def graph_cnt(self) -> int:
        """The number of graphs in the dataset."""
        return np.max(self.graph_indicator)

    def print_stats(self):
        """Prints a summary of the dataset."""
        graph_cnt = self.graph_cnt
        node_cnt = len(self.graph_indicator)
        edge_cnt = len(self.A)
        print("Dataset stats:")
        print("  {:<15}: {}".format("Graphs", graph_cnt))
        print("  {:<15}: {}".format("Nodes", node_cnt))
        print("  {:<15}: {}".format("Edges", edge_cnt))
        print("  {:<15}: {:.2f}".format("Avg. Nodes", node_cnt / graph_cnt))
        print("  {:<15}: {:.2f}".format("Avg. Edges", edge_cnt / graph_cnt))
        if self.node_attributes is not None:
            node_attr = self.node_attributes[0]
            node_attr_cnt = (
                len(node_attr)
                if isinstance(node_attr, list) or isinstance(node_attr, tuple)
                else 1
            )
            print("  {:<15}: {}".format("Node Features", node_attr_cnt))
        if self.edge_attributes is not None:
            edge_attr = self.edge_attributes[0]
            edge_attr_cnt = (
                len(edge_attr)
                if isinstance(edge_attr, list) or isinstance(edge_attr, tuple)
                else 1
            )
            edge_attr_cnt = len(self.edge_attributes[0])
            print("  {:<15}: {}".format("Edge Features", edge_attr_cnt))

    def save(self, directory):
        """Write the TUDataset to a directory. This step converts the helper
        instance RawTUDataset into a valid TUDataset. You can load the
        TUDataset from the directory where you saved the RawTUDataset.

        :param directory: The directory where the TUDataset should be stored.
        """

        def _write_to_file(data: list | None, file_name: str):
            if data is None:
                return
            if not os.path.exists(directory):
                os.mkdir(directory)
            with open("{}/{}_{}.txt".format(directory, self.name, file_name), "w") as f:
                for entry in data:
                    entry_str = str(entry).replace(" ", "")
                    if isinstance(entry, tuple):
                        entry_str = entry_str.strip("(").strip(")")
                    elif isinstance(entry, list):
                        entry_str = entry_str.strip("[").strip("]")
                    f.write("{}\n".format(entry_str))

        _write_to_file(self.A, "A")
        _write_to_file(self.graph_indicator, "graph_indicator")
        _write_to_file(self.graph_labels, "graph_labels")
        _write_to_file(self.node_attributes, "node_attributes")
        _write_to_file(self.edge_attributes, "edge_attributes")
        _write_to_file(self.ids, "ids")
