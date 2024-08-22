import matplotlib.pyplot as plt
from fgutils.proxy import MolProxy, ProxyGroup, ProxyGraph, Parser
from fgutils.vis import plot_graph

pattern = "C{g}(C)C{g}(C)(C)C"
g_2 = ProxyGroup("g", ProxyGraph("c1ccccc1", anchor=[1, 3]))
g_3 = ProxyGroup("g", ProxyGraph("c1ccccc1", anchor=[1, 3, 4]))
g_4 = ProxyGroup("g", ProxyGraph("c1ccccc1", anchor=[1, 3, 4, 5]))

parser = Parser()
proxy1 = MolProxy(pattern, g_2)
proxy2 = MolProxy(pattern, g_3)
proxy3 = MolProxy(pattern, g_4)

parent_graph = parser(pattern)
mol1 = next(proxy1)
mol2 = next(proxy2)
mol3 = next(proxy3)

fig, ax = plt.subplots(2, 2, dpi=200, figsize=(16, 9))
plot_graph(
    parent_graph, ax[0, 0], show_labels=True, show_edge_labels=False, title="parent"
)
plot_graph(mol1, ax[0, 1], show_edge_labels=False, title="2 anchor nodes")
plot_graph(mol2, ax[1, 0], show_edge_labels=False, title="3 anchor nodes")
plot_graph(mol3, ax[1, 1], show_edge_labels=False, title="4 anchor nodes")
plt.savefig(
    "doc/figures/multiple_anchor_example.png", bbox_inches="tight", transparent=True
)
plt.show()
