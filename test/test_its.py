from fgutils.its import get_its, split_its
from fgutils.parse import parse
from fgutils.rdkit import smiles_to_graph
from fgutils.const import LABELS_KEY, IS_LABELED_KEY, IDX_MAP_KEY, AAM_KEY

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
