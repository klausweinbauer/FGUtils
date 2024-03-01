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
