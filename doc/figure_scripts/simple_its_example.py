import matplotlib.pyplot as plt
from fgutils import Parser
from fgutils.vis import plot_graph


pattern = "C1<2,1>C<1,2>C<2,1>C(C)<0,1>C<2,1>C(O)<0,1>1"
parser = Parser()
g = parser(pattern)

fig, ax = plt.subplots(1, 1)
plot_graph(g, ax)
plt.savefig("doc/figures/simple_its_example.png", bbox_inches="tight", transparent=True)
plt.show()
