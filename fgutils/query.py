import copy
import networkx as nx
import numpy as np


from fgutils.utils import add_implicit_hydrogens
from fgutils.permutation import PermutationMapper
from fgutils.algorithm.subgraph import map_subgraph
from fgutils.fgconfig import FGConfig, FGConfigProvider, FGTreeNode
from fgutils.rdkit import smiles_to_graph
from fgutils.const import SYMBOL_KEY
from fgutils.chem.its import get_its


def is_functional_group(
    graph, index: int, config: FGConfig, mapper: PermutationMapper, max_id=None
):
    if max_id is None:
        max_id = np.max(list(graph.nodes))

    is_fg = False
    fg_indices = []
    mappings = map_subgraph(graph, index, config.pattern, mapper)
    for _is_fg, _mapping in mappings:
        if _is_fg:
            fg_indices = [
                m_id
                for m_id, fg_id in _mapping
                if fg_id in config.group_atoms and m_id <= max_id
            ]
            is_fg = index in fg_indices
            if is_fg:
                break

    if is_fg:
        last_len = config.max_pattern_size
        for apattern, apattern_size in sorted(
            [(m, m.number_of_nodes()) for m in config.anti_pattern],
            key=lambda x: x[1],
            reverse=True,
        ):
            if last_len > apattern_size:
                last_len = apattern_size
            mappings = map_subgraph(graph, index, apattern, mapper)
            for _is_fg, _ in mappings:
                is_fg = is_fg and not _is_fg
                if not is_fg:
                    break
    return is_fg, sorted(fg_indices)


class FGQuery:
    """
    Class to get functional groups from a molecule.

    :param mapper: (optional) The permutation mapper to use.
    :param config_provider: (optional) A functional group config provider. If
        not specified the default functional group collection will be used.
    :param require_implicit_hydrogen: Flag to specify if implicit hydrogens are
        required for the query. This usually depends on the provided FGConfigs.
        If the configured functional group patterns do not require hydrogens
        this can be set to false. (Default = True)
    """

    def __init__(
        self,
        mapper: PermutationMapper | None = None,
        config: FGConfig | list[FGConfig] | FGConfigProvider | None = None,
        require_implicit_hydrogen: bool = True,
    ):
        self.mapper = (
            mapper
            if mapper is not None
            else PermutationMapper(wildcard="R", ignore_case=True)
        )

        if isinstance(config, FGConfig) or isinstance(config, list):
            self.config_provider = FGConfigProvider(config, mapper=self.mapper)
        elif isinstance(config, FGConfigProvider):
            self.config_provider = config
        elif config is None:
            self.config_provider = FGConfigProvider(mapper=self.mapper)
        else:
            raise ValueError(
                "Invalid type '{}' for argument config.".format(type(config))
            )

        self.require_implicit_hydrogen = require_implicit_hydrogen

    def __find_best_node_rec(self, nodes: list[FGTreeNode], graph, idx, max_id):
        best_node = None
        node_indices = []
        for node in nodes:
            is_fg, fg_indices = is_functional_group(
                graph, idx, node.fgconfig, mapper=self.mapper, max_id=max_id
            )
            if is_fg:
                r_node, r_indices = self.__find_best_node_rec(
                    node.children, graph, idx, max_id
                )
                if r_node is None:
                    best_node = node
                    node_indices = fg_indices
                else:
                    best_node = r_node
                    node_indices = r_indices
        return best_node, node_indices

    def __get_functional_groups(self, graph: nx.Graph) -> list[tuple[str, list[int]]]:
        fg_candidate_ids = [
            n_id
            for n_id, n_sym in graph.nodes(data=SYMBOL_KEY)  # type: ignore
            if n_sym not in ["H", "C"]
        ]
        roots = self.config_provider.get_tree()
        groups = []
        unidentified_ids = []
        max_id = None
        if self.require_implicit_hydrogen:
            max_id = np.max(list(graph.nodes))
            graph = add_implicit_hydrogens(copy.deepcopy(graph))
        while len(fg_candidate_ids) > 0:
            atom_id = fg_candidate_ids.pop(0)
            node, indices = self.__find_best_node_rec(roots, graph, atom_id, max_id)
            if node is None:
                unidentified_ids.append(atom_id)
            else:
                assert atom_id in indices
                for i in indices:
                    if i in fg_candidate_ids:
                        fg_candidate_ids.remove(i)
                    elif i in unidentified_ids:
                        unidentified_ids.remove(i)
                groups.append((node.fgconfig.name, indices))
        return groups

    def get(self, value) -> list[tuple[str, list[int]]]:
        """
        Get all functional groups from a molecule. The query returns two
        functional groups for acetylsalicylic acid::

            >>> smiles = "O=C(C)Oc1ccccc1C(=O)O" # acetylsalicylic acid
            >>> query = FGQuery()
            >>> query.get(smiles)
            [('ester', [0, 1, 3]), ('carboxylic_acid', [10, 11, 12])]

        :param value: This is either a graph or SMILES as string.

        :returns: Returns a list of tuples. The first element in a tuple is the
            functional group name and the second element is a list of node
            indices that belong to this functional group
            ``(functional_group_name, [idx_1, idx_2, ...])``.
        """
        mol_graph = None
        if isinstance(value, nx.Graph):
            mol_graph = value
        elif isinstance(value, str):
            mol_graph = smiles_to_graph(value)
            if isinstance(mol_graph, tuple):
                mol_graph = get_its(*mol_graph)
        else:
            raise ValueError(
                "Can not interpret '{}' (type: {}) as mol graph.".format(
                    value, type(value)
                )
            )
        return self.__get_functional_groups(mol_graph)  # type: ignore
