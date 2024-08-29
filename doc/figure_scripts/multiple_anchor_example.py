import matplotlib.pyplot as plt
from fgutils.proxy import MolProxy, ProxyGroup, ProxyGraph, Parser
from fgutils.vis import GraphVisualizer

pattern = "N{g}C{g}(O)C"
g_1 = ProxyGroup("g", ProxyGraph("C1CCCCC1", anchor=[1, 3]))
g_2 = ProxyGroup("g", ProxyGraph("C1CCCCC1", anchor=[1, 3, 4]))

parser = Parser()
proxy1 = MolProxy(pattern, g_1)
proxy2 = MolProxy(pattern, g_2)

parent_graph = parser(pattern)
mol1 = next(proxy1)
mol2 = next(proxy2)

vis = GraphVisualizer(show_edge_labels=False)
fig, ax = plt.subplots(1, 3, dpi=200, figsize=(20, 3))
vis.plot(parent_graph, ax[0], title="parent")
vis.plot(mol1, ax[1], title="2 anchor nodes")
vis.plot(mol2, ax[2], title="3 anchor nodes")
plt.savefig(
    "doc/figures/multiple_anchor_example.png", bbox_inches="tight", transparent=True
)
plt.show()
