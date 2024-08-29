import matplotlib.pyplot as plt
from fgutils import Parser
from fgutils.proxy import MolProxy, ProxyGroup
from fgutils.vis import GraphVisualizer

pattern = "CC(=O)O{propyl}"
propyl_group = ProxyGroup("propyl", "CCC")
parser = Parser()
proxy = MolProxy(pattern, propyl_group, parser=parser)

g = parser(pattern)
mol = next(proxy)

vis = GraphVisualizer()
fig, ax = plt.subplots(1, 2, dpi=100, figsize=(12, 3))
vis.plot(g, ax[0], title="Core Pattern")
vis.plot(mol, ax[1], title="Generated Molecule")
plt.savefig(
    "doc/figures/labeled_node_example.png", bbox_inches="tight", transparent=True
)
plt.show()
