import copy
import itertools
import numpy as np
import networkx as nx

from fgutils.algorithm.subgraph_enumeration import (
    node_induced_connected_subgraphs,
    is_valid_extension,
)

expected_extensions = [
    [0],
    [0, 1],
    [0, 1, 2],
    [0, 1, 2, 3],
    [0, 1, 2, 3, 4],
    [0, 1, 2, 3, 4, 5],
    [0, 1, 2, 3, 5],
    [0, 1, 2, 5],
    [0, 1, 2, 4],
    [0, 1, 2, 4, 5],
    [0, 1, 3],
    [0, 1, 3, 5],
    [0, 1, 3, 5, 4],
    [0, 2],
    [0, 2, 4],
    [0, 2, 4, 5],
    [0, 2, 4, 5, 3],
    [0, 2, 5],
    [0, 2, 5, 3],
]


def found_all_extensions(subgraphs, ext_extensions=None):
    if ext_extensions is None:
        ext_extensions = expected_extensions
    found = [0] * len(ext_extensions)
    for subgraph in subgraphs:
        for i, ext in enumerate(ext_extensions):
            if subgraph == ext:
                found[i] = 1
    return np.sum(found) == len(ext_extensions)


def assert_extensions(subgraphs, ext_extensions=None):
    assert found_all_extensions(subgraphs, ext_extensions)


def test_enumeration():
    g = nx.Graph()
    g.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (3, 5), (4, 5)])
    subgraphs = list(node_induced_connected_subgraphs(g, 0))
    assert_extensions(subgraphs)


def test_enumeration2():
    g = nx.Graph()
    g.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (3, 5), (4, 5)])
    g = nx.relabel_nodes(g, {0: 3, 3: 0}, copy=True)
    subgraphs = list(node_induced_connected_subgraphs(g, 3))
    for sg in subgraphs:
        print(sg)
    relabeled_ext = copy.deepcopy(expected_extensions)
    for i, a in enumerate(relabeled_ext):
        for j, b in enumerate(a):
            if b == 0:
                b = 3
            elif b == 3:
                b = 0
            relabeled_ext[i][j] = b
    assert_extensions(subgraphs, relabeled_ext)


def test_indipendence_to_node_numberings():
    g = nx.Graph()
    g.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (3, 5), (4, 5)])
    failed = 0
    for permute in itertools.permutations(range(6)):
        nmap = {i: n for i, n in enumerate(permute)}
        _g = nx.relabel_nodes(g, nmap, copy=True)
        subgraphs = list(node_induced_connected_subgraphs(_g, nmap[0]))
        if len(subgraphs) != 19:
            failed += 1
    assert failed == 0


def test_is_valid_extension():
    D = [1, 0, 1, 2, 1, 2, 3, 4, 3, 3, np.inf, np.inf, np.inf, 4]
    U = [2, 3, 5, 4, 7, 9]
    v = 10
    assert is_valid_extension(U, v, D)
