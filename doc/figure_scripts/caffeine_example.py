import matplotlib.pyplot as plt
from fgutils import Parser
from fgutils.vis import plot_as_mol

parser = Parser()
mol = parser("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

fig, ax = plt.subplots(1, 1, dpi=200)
plot_as_mol(mol, ax)
plt.savefig("doc/figures/caffeine_example.png", bbox_inches="tight", transparent=True)
