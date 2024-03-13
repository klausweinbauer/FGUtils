import copy
import networkx as nx

from fgutils.utils import add_implicit_hydrogens
from fgutils.permutation import PermutationMapper
from fgutils.mapping import map_pattern
from fgutils.fgconfig import FGConfig, FGConfigProvider, FGTreeNode


def is_functional_group(graph, index: int, config: FGConfig, mapper: PermutationMapper):
    max_id = len(graph)
    graph = add_implicit_hydrogens(copy.deepcopy(graph))

    is_fg, mapping = map_pattern(graph, index, config.pattern, mapper)
    fg_indices = []
    if is_fg:
        fg_indices = [
            m_id
            for m_id, fg_id in mapping
            if fg_id in config.group_atoms and m_id < max_id
        ]
        is_fg = index in fg_indices

    if is_fg:
        last_len = config.max_pattern_size
        for apattern, apattern_size in sorted(
            [(m, m.number_of_nodes()) for m in config.anti_pattern],
            key=lambda x: x[1],
            reverse=True,
        ):
            if last_len > apattern_size:
                last_len = apattern_size
            is_match, _ = map_pattern(graph, index, apattern, mapper)
            is_fg = is_fg and not is_match
            if not is_fg:
                break
    return is_fg, sorted(fg_indices)


class FGQuery:
    def __init__(
        self,
        use_smiles=False,
        mapper: PermutationMapper | None = None,
        config_provider: FGConfigProvider | None = None,
    ):
        self.use_smiles = use_smiles
        self.mapper = (
            mapper
            if mapper is not None
            else PermutationMapper(wildcard="R", ignore_case=True)
        )
        self.config_provider = (
            config_provider
            if config_provider is not None
            else FGConfigProvider(mapper=self.mapper)
        )

    def __find_best_node_rec(self, nodes: list[FGTreeNode], graph, idx):
        best_node = None
        node_indices = []
        for node in nodes:
            is_fg, fg_indices = is_functional_group(
                graph, idx, node.fgconfig, mapper=self.mapper
            )
            if is_fg:
                r_node, r_indices = self.__find_best_node_rec(
                    node.children,
                    graph,
                    idx,
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
            for n_id, n_sym in graph.nodes(data="symbol")  # type: ignore
            if n_sym not in ["H", "C"]
        ]
        roots = self.config_provider.get_tree()
        groups = []
        unidentified_ids = []
        while len(fg_candidate_ids) > 0:
            atom_id = fg_candidate_ids.pop(0)
            node, indices = self.__find_best_node_rec(roots, graph, atom_id)
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
        if len(unidentified_ids) > 0:
            raise RuntimeError(
                "Could not find a functional group for atom(s) {}.".format(
                    [
                        "{}@{}".format(graph.nodes[i]["symbol"], i)
                        for i in unidentified_ids
                    ]
                )
            )
        return groups

    def get(self, value) -> list[tuple[str, list[int]]]:
        mol_graph = None
        if isinstance(value, nx.Graph):
            mol_graph = value
        elif self.use_smiles:
            import rdkit.Chem.rdmolfiles as rdmolfiles
            from fgutils.rdkit import mol_to_graph

            mol = rdmolfiles.MolFromSmiles(value)
            mol_graph = mol_to_graph(mol)
        else:
            raise ValueError(
                "Can not interpret '{}' (type: {}) as mol graph.".format(
                    value, type(value)
                )
            )
        return self.__get_functional_groups(mol_graph)  # type: ignore
