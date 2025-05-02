import os
import numpy as np

from fgutils.proxy import relabel_graph
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
        self.failed_conversions = []

    @staticmethod
    def from_reaction_dict(
        data: list[dict],
        id_col: str,
        its_col: str,
        label_col: str,
        dataset_name: str,
        skip_failed_conversion=False,
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
        :param skip_failed_conversion: (optional) If set to true failed
            conversions will not raise an exception. Failed conversions are
            then stored in the self.failed_conversions property. (Default:
            False)

        :returns: A instance of RawTUDataset
        """
        DS_A = []
        DS_ids = []
        DS_graph_indicator = []
        DS_graph_labels = []
        DS_node_attributes = []
        DS_edge_attributes = []
        failed_conversions = []

        node_idx_offset = 1

        node_cnt = 0
        edge_cnt = 0
        graph_cnt = len(data)

        for entry in data:
            try:
                _DS_A = []
                _DS_ids = []
                _DS_graph_indicator = []
                _DS_graph_labels = []
                _DS_node_attributes = []
                _DS_edge_attributes = []

                its = entry[its_col]
                if isinstance(its, ITS):
                    its = its.graph
                its = relabel_graph(its, offset=node_idx_offset)
                _DS_graph_labels.append(entry[label_col])
                _DS_ids.append(entry[id_col])
                for u, d in its.nodes(data=True):
                    _DS_graph_indicator.append(len(DS_ids) + 1)
                    _DS_node_attributes.append(get_atomic_number(d[SYMBOL_KEY]))
                for u, v, d in its.edges(data=True):
                    _DS_A.extend([(u, v), (v, u)])
                    g_b, h_b = d[BOND_KEY]
                    g_b = 0 if g_b is None else g_b
                    h_b = 0 if h_b is None else h_b
                    _DS_edge_attributes.extend([(g_b, h_b), (g_b, h_b)])

                node_cnt += len(its.nodes)
                edge_cnt += 2 * len(its.edges)
                node_idx_offset += len(its.nodes)
                DS_A.extend(_DS_A)
                DS_ids.extend(_DS_ids)
                DS_graph_indicator.extend(_DS_graph_indicator)
                DS_graph_labels.extend(_DS_graph_labels)
                DS_node_attributes.extend(_DS_node_attributes)
                DS_edge_attributes.extend(_DS_edge_attributes)
            except Exception as e:
                if skip_failed_conversion:
                    failed_conversions.append((entry, e))
                else:
                    raise e
        assert edge_cnt == len(DS_A)
        assert node_cnt == len(DS_graph_indicator)
        assert graph_cnt == len(DS_graph_labels) + len(failed_conversions)
        assert node_cnt == len(DS_node_attributes)
        assert edge_cnt == len(DS_edge_attributes)
        assert graph_cnt == len(DS_ids) + len(failed_conversions)

        rawds = RawTUDataset(
            dataset_name,
            DS_A,
            DS_graph_indicator,
            graph_labels=DS_graph_labels,
            node_attributes=DS_node_attributes,
            edge_attributes=DS_edge_attributes,
            ids=DS_ids,
        )
        rawds.failed_conversions = failed_conversions
        return rawds

    @staticmethod
    def write_to_disk(
        data: list[dict],
        id_col: str,
        its_col: str,
        label_col: str,
        dataset_name: str,
        directory: str,
        skip_failed_conversion=False,
    ):
        """Generate a TUDataset instance from a reaction dict and write it to
        disk. Each entry in the data list contains one reaction. Each entry
        (dict) needs to provide a reaction id, an ITS graph, and a label. The
        label is the reaction class for example. In the ML context the label is
        the classification/regression target.

        :param data: List of dictionaries. Each dict requires the keys:
            ``id_col``, ``its_col``, ``label_col``
        :param id_col: The dict key where to find the reaction id.
        :param its_col: The dict key where to find the ITS graph.
        :param label_col: The dict key where to find the graph label.
        :param directory: The directory where the dataset files are stored.
        :param skip_failed_conversion: (optional) If set to true failed
            conversions will not raise an exception. Failed conversions are
            then stored in the self.failed_conversions property. (Default:
            False)

        :returns: A list of failed conversions.
        """

        def _write_to_file(data: list | None, file):
            if data is None:
                return
            if not os.path.exists(directory):
                os.mkdir(directory)
            for entry in data:
                entry_str = str(entry).replace(" ", "")
                if isinstance(entry, tuple):
                    entry_str = entry_str.strip("(").strip(")")
                elif isinstance(entry, list):
                    entry_str = entry_str.strip("[").strip("]")
                file.write("{}\n".format(entry_str))

        f_A = open("{}/{}_{}.txt".format(directory, dataset_name, "A"), "w")
        f_ids = open("{}/{}_{}.txt".format(directory, dataset_name, "ids"), "w")
        f_graph_indicator = open(
            "{}/{}_{}.txt".format(directory, dataset_name, "graph_indicator"), "w"
        )
        f_graph_labels = open(
            "{}/{}_{}.txt".format(directory, dataset_name, "graph_labels"), "w"
        )
        f_node_attributes = open(
            "{}/{}_{}.txt".format(directory, dataset_name, "node_attributes"), "w"
        )
        f_edge_attributes = open(
            "{}/{}_{}.txt".format(directory, dataset_name, "edge_attributes"), "w"
        )
        failed_conversions = []

        node_idx_offset = 1
        graph_i = 0

        for entry in data:
            try:
                _DS_A = []
                _DS_ids = []
                _DS_graph_indicator = []
                _DS_graph_labels = []
                _DS_node_attributes = []
                _DS_edge_attributes = []

                its = entry[its_col]
                if isinstance(its, ITS):
                    its = its.graph
                its = relabel_graph(its, offset=node_idx_offset)
                _DS_graph_labels.append(entry[label_col])
                _DS_ids.append(entry[id_col])
                for u, d in its.nodes(data=True):
                    _DS_graph_indicator.append(graph_i + 1)
                    _DS_node_attributes.append(get_atomic_number(d[SYMBOL_KEY]))
                for u, v, d in its.edges(data=True):
                    _DS_A.extend([(u, v), (v, u)])
                    g_b, h_b = d[BOND_KEY]
                    g_b = 0 if g_b is None else g_b
                    h_b = 0 if h_b is None else h_b
                    _DS_edge_attributes.extend([(g_b, h_b), (g_b, h_b)])

                node_idx_offset += len(its.nodes)
                graph_i += 1
                _write_to_file(_DS_A, f_A)
                _write_to_file(_DS_ids, f_ids)
                _write_to_file(_DS_graph_indicator, f_graph_indicator)
                _write_to_file(_DS_graph_labels, f_graph_labels)
                _write_to_file(_DS_node_attributes, f_node_attributes)
                _write_to_file(_DS_edge_attributes, f_edge_attributes)
            except Exception as e:
                if skip_failed_conversion:
                    failed_conversions.append((entry, e))
                else:
                    raise e

        f_A.close()
        f_ids.close()
        f_graph_indicator.close()
        f_graph_labels.close()
        f_node_attributes.close()
        f_edge_attributes.close()

        return failed_conversions

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
