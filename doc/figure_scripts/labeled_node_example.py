import matplotlib.pyplot as plt
from fgutils import Parser
from fgutils.proxy import MolProxy, ProxyGroup
from fgutils.vis import plot_graph

pattern = "CC(=O)O{propyl}"
propyl_group = ProxyGroup("propyl", pattern="CCC")
parser = Parser()
proxy = MolProxy(pattern, propyl_group, parser=parser)

g = parser(pattern)
mol = next(proxy)

fig, ax = plt.subplots(1, 2, dpi=100, figsize=(12, 4))
plot_graph(g, ax[0], show_labels=True)
plot_graph(mol, ax[1])
plt.savefig(
    "doc/figures/labeled_node_example.png", bbox_inches="tight", transparent=True
)
plt.show()
