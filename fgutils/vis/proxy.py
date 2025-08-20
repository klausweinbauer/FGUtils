from fgutils.parse import Parser
from fgutils.proxy import ProxyGraph

from .common import GraphVisualizer, AnchorNodeLabelFormatter


class ProxyVisualizer:
    def __init__(self, parser=None, use_mol_coords=True):
        self.parser = Parser(use_multigraph=True) if parser is None else parser
        self.graph_visualizer = GraphVisualizer(
            node_label_formatter=AnchorNodeLabelFormatter(),
            use_mol_coords=use_mol_coords,
        )

    def plot_graph(self, proxy_graph: ProxyGraph, ax, title=None):
        graph = self.parser(proxy_graph.pattern)
        self.graph_visualizer.node_label_formatter.anchor = proxy_graph.anchor  # type: ignore
        self.graph_visualizer.plot(graph, ax, title=title)
