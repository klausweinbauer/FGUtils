from fgutils.chem.its import get_its
from fgutils.parse import parse
from fgutils.rdkit import smiles_to_graph
from fgutils.const import LABELS_KEY, IS_LABELED_KEY, IDX_MAP_KEY

from test.my_asserts import assert_graph_eq


def test_get_its():
    g, h = smiles_to_graph("[C:1][O:2].[C:3]>>[C:1].[O:2][C:3]")
    exp_its = parse("C<1,0>O<0,1>C")
    its = get_its(g, h)
    assert_graph_eq(exp_its, its, ignore_keys=[LABELS_KEY, IS_LABELED_KEY, IDX_MAP_KEY])
