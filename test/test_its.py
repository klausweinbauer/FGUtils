import networkx as nx

from fgutils.torch import its_from_torch, its_to_torch
from fgutils.its import get_its, split_its, ITS
from fgutils.parse import parse
from fgutils.rdkit import smiles_to_graph
from fgutils.const import (
    LABELS_KEY,
    IS_LABELED_KEY,
    IDX_MAP_KEY,
    AAM_KEY,
    BOND_KEY,
    SYMBOL_KEY,
)

from test.my_asserts import assert_graph_eq
from .test_parse import _assert_graph


def test_split_its():
    its = parse("C1<2,>C<,2>C<2,>C(C)<0,1>C<2,>C(C(=O)O)<0,1>1")
    exp_nodes = {i: "C" for i in range(10)}
    exp_nodes[8] = "O"
    exp_nodes[9] = "O"
    exp_edges_g = [
        (0, 1, 2),
        (1, 2, 1),
        (2, 3, 2),
        (3, 4, 1),
        (5, 6, 2),
        (6, 7, 1),
        (7, 8, 2),
        (7, 9, 1),
    ]
    exp_edges_h = [
        (0, 1, 1),
        (0, 6, 1),
        (1, 2, 2),
        (2, 3, 1),
        (3, 4, 1),
        (3, 5, 1),
        (5, 6, 1),
        (6, 7, 1),
        (7, 8, 2),
        (7, 9, 1),
    ]
    g, h = split_its(its)
    _assert_graph(g, exp_nodes, exp_edges_g)
    _assert_graph(h, exp_nodes, exp_edges_h)


def test_get_its():
    g, h = smiles_to_graph("[C:1][O:2].[C:3]>>[C:1].[O:2][C:3]")
    exp_its = parse("C<1,0>O<0,1>C")
    its = get_its(g, h)
    assert_graph_eq(
        exp_its, its, ignore_keys=[LABELS_KEY, IS_LABELED_KEY, IDX_MAP_KEY, AAM_KEY]
    )


def test_get_its_unbalanced_but_mapped():
    smiles = "[C:1][O:2].[N:3]>>[C:1][N:3]"
    exp_its = parse("C<0,1>N")
    g, h = smiles_to_graph(smiles)
    its = get_its(g, h)
    assert_graph_eq(
        exp_its, its, ignore_keys=[LABELS_KEY, IS_LABELED_KEY, IDX_MAP_KEY, AAM_KEY]
    )


def test_its_isomorphism():
    g1 = nx.Graph()
    g1.add_node(0, **{SYMBOL_KEY: "C"})
    g1.add_node(1, **{SYMBOL_KEY: "O"})
    g1.add_edge(0, 1, **{BOND_KEY: (1, 2)})

    g2 = nx.Graph()
    g2.add_node(0, **{SYMBOL_KEY: "C"})
    g2.add_node(1, **{SYMBOL_KEY: "O"})
    g2.add_edge(0, 1, **{BOND_KEY: (1.0, 2.0)})

    g3 = nx.Graph()
    g3.add_node(0, **{SYMBOL_KEY: "C"})
    g3.add_node(1, **{SYMBOL_KEY: "O"})
    g3.add_edge(0, 1, **{BOND_KEY: [1, 2]})

    its1 = ITS(g1)
    its2 = ITS(g2)
    its3 = ITS(g3)

    its1.standardize()
    its2.standardize()
    its3.standardize()

    h1 = its1.wl_hash
    h2 = its2.wl_hash
    h3 = its3.wl_hash

    assert h1 == h2
    assert h2 == h3


def test_its_isomorphism2():
    smiles = (
        "[c:3]1[c:4][c:5][c:6][c:7][c:8]1[C:1][O:2]"
        + ">>[c:3]1[c:4][c:5][c:6][c:7][c:8]1[C:1]=[O:2]"
    )
    its1 = ITS.from_smiles(smiles)
    torch_its = its_to_torch(its1.graph)
    its2 = ITS(its_from_torch(torch_its))

    its1.standardize()
    its2.standardize()
    assert its1.wl_hash == its2.wl_hash
