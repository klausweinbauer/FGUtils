import re
import copy
import itertools
import numpy as np
import networkx as nx

from fgutils.parse import parse
from fgutils.permutation import Mapper

def print_graph(graph):
    print(
        "Graph Nodes: {}".format(
            " ".join(
                ["{}[{}]".format(n[1]["symbol"], n[0]) for n in graph.nodes(data=True)]
            )
        )
    )
    print(
        "Graph Edges: {}".format(
            " ".join(
                [
                    "[{}]-[{}]:{}".format(n[0], n[1], n[2]["bond"])
                    for n in graph.edges(data=True)
                ]
            )
        )
    )


class FGConfig:
    def __init__(self, **kwargs):
        self.parent = None

        pattern = kwargs.get("pattern", [])
        pattern = pattern if isinstance(pattern, list) else [pattern]
        self.pattern = [parse(p) for p in pattern]

        self.name = kwargs.get("name", None)
        if self.name is None:
            raise ValueError(
                "Functional group config requires a name. Add 'name' property to config."
            )

        group_atoms = kwargs.get("group_atoms", None)
        if group_atoms is None:
            group_atoms = []
            for p in self.pattern:
                group_atoms.append(list(p.nodes))
        else:
            group_atoms = (
                group_atoms if isinstance(group_atoms[0], list) else [group_atoms]
            )
        self.group_atoms = group_atoms

        anti_pattern = kwargs.get("anti_pattern", [])
        anti_pattern = (
            anti_pattern if isinstance(anti_pattern, list) else [anti_pattern]
        )
        self.anti_pattern = sorted(
            [parse(p) for p in anti_pattern],
            key=lambda x: x.number_of_nodes(),
            reverse=True,
        )

        depth = kwargs.get("depth", None)
        self.max_pattern_size = (
            depth
            if depth is not None
            else np.max([p.number_of_nodes() for p in self.pattern + self.anti_pattern])
        )

        subgroups = []
        for sgc in kwargs.get("subgroups", []):
            fgc = FGConfig(**sgc)
            fgc.parent = self
            subgroups.append(fgc)
        self.subgroups = subgroups


def pattern_match(graph, anchor, pattern, pattern_anchor=None, verbose=False):
    if pattern_anchor is None:
        for pidx in pattern.nodes:
            result = anchored_pattern_match(
                graph, anchor, pattern, pidx, verbose=verbose
            )
            if result[0]:
                return result
        return False, [[]]
    else:
        return anchored_pattern_match(
            graph, anchor, pattern, pattern_anchor, verbose=verbose
        )


functional_group_config = [
    {
        "name": "carbonyl",
        "pattern": "C=O",
        "group_atoms": [0, 1],
        "subgroups": [
            {
                "name": "aldehyde",
                "pattern": "RC=O",
                "anti_pattern": ["NC=O"],
                "group_atoms": [1, 2],
            },
            {
                "name": "ketone",
                "pattern": "RC(=O)R",
                "group_atoms": [1, 2],
                "anti_pattern": ["RC(=O)O"],
            },
            {
                "name": "carboxylic_acid",
                "pattern": "RC(=O)O",
                "group_atoms": [1, 2, 3],
                "anti_pattern": ["RC(=O)OR"],
            },
            {
                "name": "carboxylate_ester",
                "pattern": "RC(=O)OR",
                "group_atoms": [1, 2, 3],
            },
            {"name": "amide", "pattern": "C(=O)N"},
        ],
    },
]

fg_configs = None


def get_FGs() -> list[FGConfig]:
    global fg_configs
    if fg_configs is None:
        c = []
        for fgc in functional_group_config:
            c.append(FGConfig(**fgc))
        fg_configs = c
    return fg_configs


def get_FGs_flat() -> list[FGConfig]:
    def _get(conf: FGConfig) -> list[FGConfig]:
        configs = [conf]
        for sg in conf.subgroups:
            configs.extend(_get(sg))
        return configs

    configs = []
    for fg in get_FGs():
        configs.extend(_get(fg))
    return configs


def get_FG_by_name(name: str) -> FGConfig:
    for fg in get_FGs_flat():
        if fg.name == name:
            return fg
    raise KeyError("No functional group config with name '{}' found.".format(name))


def get_FG_names() -> list[str]:
    return [c.name for c in get_FGs_flat()]


def get_FG_root_chain(name: str) -> list[FGConfig]:
    fg = get_FG_by_name(name)
    chain = [fg]
    while fg.parent is not None:
        chain.insert(0, fg.parent)
        fg = fg.parent
    return chain


def check_functional_group(graph, config: FGConfig, index: int, verbose=False) -> bool:
    is_func_group = False

    for pattern, group_indices in zip(config.pattern, config.group_atoms):
        is_match, match = pattern_match(graph, index, pattern, verbose=verbose)
        if is_match:
            m_idx = [m[1] for m in match if m[0] == index][0]
            is_func_group = is_func_group or m_idx in group_indices

    last_len = config.max_pattern_size
    for apattern, apattern_size in sorted(
        [(m, m.number_of_nodes()) for m in config.anti_pattern],
        key=lambda x: x[1],
        reverse=True,
    ):
        if not is_func_group:
            break
        if last_len > apattern_size:
            last_len = apattern_size
        is_match, _ = pattern_match(graph, index, apattern, verbose=verbose)
        is_func_group = is_func_group and not is_match
    return is_func_group


def mol_to_graph(mol) -> nx.Graph:
    bond_order_map = {
        "SINGLE": 1,
        "DOUBLE": 2,
        "TRIPLE": 3,
        "QUADRUPLE": 4,
        "AROMATIC": 1.5,
    }
    g = nx.Graph()
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx(), symbol=atom.GetSymbol())
    for bond in mol.GetBonds():
        bond_type = str(bond.GetBondType()).split(".")[-1]
        g.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond=bond_order_map[bond_type],
        )
    return g
