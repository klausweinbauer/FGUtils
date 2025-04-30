import pytest
import matplotlib.pyplot as plt

from fgutils.rdkit import smiles_to_graph
from fgutils.vis import plot_as_mol


@pytest.mark.parametrize(
    "smiles",
    [("c1cc[nH]c1"), ("[CH3][CH2][OH]"), ("CC(=O)O")],
)
def test_plot_as_mol(smiles):
    _, ax = plt.subplots()
    g = smiles_to_graph(smiles)
    plot_as_mol(g, ax=ax)
