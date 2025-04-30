import io
import tqdm
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import rdkit.Chem.rdmolfiles as rdmolfiles
import rdkit.Chem.Draw.rdMolDraw2D as rdMolDraw2D
import rdkit.Chem.rdChemReactions as rdChemReactions

from PIL import Image
from matplotlib.backends.backend_pdf import PdfPages

from fgutils.rdkit import graph_to_smiles, get_mol_coords
from fgutils.const import (
    SYMBOL_KEY,
    AAM_KEY,
    BOND_KEY,
    IS_LABELED_KEY,
    LABELS_KEY,
    BOND_CHAR_MAP,
)
from fgutils.parse import Parser
from fgutils.proxy import ProxyGraph


def plot_its(its, ax=None, use_mol_coords=True, title=None):
    if not isinstance(its, nx.Graph):
        its = its.graph

    if ax is None:
        _, _ax = plt.subplots()
    else:
        _ax = ax

    its = its.copy()
    positions = get_node_positions(its, use_mol_coords=use_mol_coords)
    its.remove_nodes_from([n for n in its.nodes if n not in positions.keys()])

    _ax.axis("equal")
    _ax.axis("off")
    if title is not None:
        _ax.set_title(title)

    nx.draw_networkx_edges(its, positions, edge_color="#000000", ax=_ax)
    nx.draw_networkx_nodes(its, positions, node_color="#FFFFFF", node_size=500, ax=_ax)

    labels = {n: "{}:{}".format(d[SYMBOL_KEY], n) for n, d in its.nodes(data=True)}
    edge_labels = {}
    for u, v, d in its.edges(data=True):
        bc1 = d[BOND_KEY][0]
        bc2 = d[BOND_KEY][1]
        if bc1 == bc2:
            continue
        if bc1 in BOND_CHAR_MAP.keys():
            bc1 = BOND_CHAR_MAP[bc1]
        if bc2 in BOND_CHAR_MAP.keys():
            bc2 = BOND_CHAR_MAP[bc2]
        edge_labels[(u, v)] = "({},{})".format(bc1, bc2)

    nx.draw_networkx_labels(its, positions, labels=labels, ax=_ax)
    nx.draw_networkx_edge_labels(its, positions, edge_labels=edge_labels, ax=_ax)

    if ax is None:
        plt.show()


def plot_as_mol(g: nx.Graph, ax=None, use_mol_coords=True):
    if ax is None:
        _, _ax = plt.subplots()
    else:
        _ax = ax

    g = g.copy()
    positions = get_node_positions(g, use_mol_coords=use_mol_coords)
    g.remove_nodes_from([n for n in g.nodes if n not in positions.keys()])

    _ax.axis("equal")
    _ax.axis("off")

    nx.draw_networkx_edges(g, positions, edge_color="#909090", ax=_ax)
    nx.draw_networkx_nodes(g, positions, node_color="#FFFFFF", node_size=500, ax=_ax)

    labels = {n: "{}".format(d[SYMBOL_KEY]) for n, d in g.nodes(data=True)}
    edge_labels = {}
    for u, v, d in g.edges(data=True):
        bc = d[BOND_KEY]
        if bc in BOND_CHAR_MAP.keys():
            bc = BOND_CHAR_MAP[bc]
        edge_labels[(u, v)] = "{}".format(bc)

    nx.draw_networkx_labels(g, positions, labels=labels, ax=_ax)
    nx.draw_networkx_edge_labels(g, positions, edge_labels=edge_labels, ax=_ax)

    if ax is None:
        plt.show()


def get_rxn_img(smiles, background_colour=None):
    drawer = rdMolDraw2D.MolDraw2DCairo(1600, 900)
    opts = drawer.drawOptions()
    if background_colour is None:
        background_colour = (0.0, 0.0, 0.0, 0.0)
    opts.setBackgroundColour(background_colour)
    if ">>" in smiles:
        rxn = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
        drawer.DrawReaction(rxn)
    else:
        mol = rdmolfiles.MolFromSmiles(smiles)
        if mol is None:
            mol = rdmolfiles.MolFromSmarts(smiles)
        drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img = Image.open(io.BytesIO(drawer.GetDrawingText()))
    nonwhite_positions = [
        (x, y)
        for x in range(img.size[0])
        for y in range(img.size[1])
        if img.getdata()[x + y * img.size[0]] != background_colour  # type: ignore
    ]
    rect = (
        min([x - 10 for x, _ in nonwhite_positions]),
        min([y - 10 for _, y in nonwhite_positions]),
        max([x + 10 for x, _ in nonwhite_positions]),
        max([y + 10 for _, y in nonwhite_positions]),
    )
    return img.crop(rect)


def plot_reaction(g: nx.Graph, h: nx.Graph, ax=None, title=None):
    if ax is None:
        _, _ax = plt.subplots()
    else:
        _ax = ax

    _ax.axis("off")
    if title is not None:
        _ax.set_title(title)
    rxn_smiles = "{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h))
    _ax.imshow(get_rxn_img(rxn_smiles))

    if ax is None:
        plt.show()


class EdgeLabelFormatter:
    def __init__(self, show_single_bonds=False, rc_only=False):
        self.rc_only = rc_only
        self.show_single_bonds = show_single_bonds
        self.bond_chars = BOND_CHAR_MAP

    def __call__(self, e, d):
        bc = d[BOND_KEY]
        if isinstance(bc, tuple):
            bc1 = bc[0]
            bc2 = bc[1]
            if bc1 in self.bond_chars.keys():
                bc1 = self.bond_chars[bc1]
            if bc2 in self.bond_chars.keys():
                bc2 = self.bond_chars[bc2]
            assert bc1 != bc2
            return "({},{})".format(bc1, bc2)
        elif not self.rc_only:
            if bc == 1 and not self.show_single_bonds:
                return ""
            if bc in self.bond_chars.keys():
                bc = self.bond_chars[bc]
            return "{}".format(bc)
        else:
            return ""


class NodeLabelFormatter:
    def __init__(self, show_aam=False):
        self.show_aam = show_aam

    def __call__(self, n, d):
        if self.show_aam:
            return "{}:{}".format(d[SYMBOL_KEY], d.get(AAM_KEY, n))
        else:
            return "{}".format(d[SYMBOL_KEY])


class AnchorNodeLabelFormatter(NodeLabelFormatter):
    def __init__(self, show_aam=False, anchor=[], format_str="[{}]"):
        self.anchor = anchor
        self.format_str = format_str
        super().__init__(show_aam=show_aam)

    def __call__(self, n, d):
        lbl = super().__call__(n, d)
        if n in self.anchor:
            lbl = self.format_str.format(lbl)
        return lbl


class LabelLegendFormatter:
    def __init__(self):
        pass

    def get_key(self, i, n, d):
        return "{{g{}}}".format(i + 1)

    def __call__(self, lbl_dict):
        lines = []
        for k, v in lbl_dict.items():
            list_str = "[{}]".format(",".join([e for e in v]))
            line = "{}: {}".format(k, list_str)
            lines.append(line)
        return "\n".join(lines)


def get_node_positions(g: nx.Graph, use_mol_coords: bool = True):
    if use_mol_coords:
        positions = get_mol_coords(g)
    else:
        positions = nx.spring_layout(g)
    return positions


class GraphVisualizer:
    def __init__(
        self,
        use_mol_coords=True,
        edge_color="#909090",
        node_color="#FFFFFF",
        node_size=500,
        node_label_formatter=None,
        edge_label_formatter=None,
        label_legend_formatter=None,
        show_node_labels=True,
        show_edge_labels=True,
        show_label_legend=True,
        label_legend_position=(0.5, 0.05),
    ):
        self.use_mol_coords = use_mol_coords
        self.edge_color = edge_color
        self.node_color = node_color
        self.node_size = node_size
        self.node_label_formatter = (
            NodeLabelFormatter()
            if node_label_formatter is None
            else node_label_formatter
        )
        self.edge_label_formatter = (
            EdgeLabelFormatter()
            if edge_label_formatter is None
            else edge_label_formatter
        )
        self.label_legend_formatter = (
            LabelLegendFormatter()
            if label_legend_formatter is None
            else label_legend_formatter
        )
        self.show_node_labels = show_node_labels
        self.show_edge_labels = show_edge_labels
        self.show_label_legend = show_label_legend
        self.label_legend_position = label_legend_position
        self.connectionstyle = ["arc3,rad={}".format(0.3 * i) for i in range(4)]

    @staticmethod
    def build_label_dict(graph: nx.Graph, formatter: LabelLegendFormatter):
        g = graph.copy()
        lbl_node_idx = 0
        lbl_dict = {}
        for n, d in g.nodes(data=True):
            assert d is not None
            if d[IS_LABELED_KEY]:
                _sym = formatter.get_key(lbl_node_idx, n, d)
                d[SYMBOL_KEY] = _sym
                lbl_node_idx += 1
                lbl_dict[_sym] = d[LABELS_KEY]
        return lbl_dict, g

    @staticmethod
    def get_legend_position(graph: nx.Graph, offset=(0, 0), use_mol_coords=True):
        node_positions = []
        for _, v in get_node_positions(graph, use_mol_coords).items():
            node_positions.append(v)
        node_positions = np.array(node_positions)
        return np.min(node_positions, axis=0) + np.array(offset)

    def plot(self, graph, ax, title=None, **kwargs):
        label_legend_position = kwargs.get(
            "label_legend_position", self.label_legend_position
        )
        positions = get_node_positions(graph, self.use_mol_coords)

        if self.show_label_legend:
            lbl_dict, graph = self.build_label_dict(
                graph, formatter=self.label_legend_formatter
            )
            ax.text(
                label_legend_position[0],
                label_legend_position[1],
                self.label_legend_formatter(lbl_dict),
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
            )

        ax.axis("equal")
        ax.axis("off")

        if title is not None:
            ax.set_title(title)

        nx.draw_networkx_edges(
            graph,
            positions,
            edge_color=self.edge_color,
            ax=ax,
            connectionstyle=self.connectionstyle,
        )
        nx.draw_networkx_nodes(
            graph,
            positions,
            node_color=self.node_color,
            node_size=self.node_size,
            ax=ax,
        )

        labels = {}
        for n, d in graph.nodes(data=True):
            labels[n] = self.node_label_formatter(n, d)

        edge_labels = {}
        if isinstance(graph, nx.MultiGraph):
            for u, v, i, d in graph.edges(data=True, keys=True):
                edge_labels[(u, v, i)] = self.edge_label_formatter((u, v), d)
        else:
            for u, v, d in graph.edges(data=True):
                edge_labels[(u, v)] = self.edge_label_formatter((u, v), d)

        if self.show_node_labels:
            nx.draw_networkx_labels(graph, positions, labels=labels, ax=ax)
        if self.show_edge_labels:
            nx.draw_networkx_edge_labels(
                graph,
                positions,
                edge_labels=edge_labels,
                ax=ax,
                connectionstyle=self.connectionstyle,
            )


class ProxyVisualizer:
    def __init__(self, parser=None, use_mol_coords=True):
        self.parser = Parser(use_multigraph=True) if parser is None else parser
        self.graph_visualizer = GraphVisualizer(
            node_label_formatter=AnchorNodeLabelFormatter(),
            use_mol_coords=use_mol_coords,
        )

    def plot_graph(self, proxy_graph: ProxyGraph, ax, title=None):
        graph = self.parser(proxy_graph.pattern)
        self.graph_visualizer.node_label_formatter.anchor = proxy_graph.anchor
        self.graph_visualizer.plot(graph, ax, title=title)


class PdfWriter:
    """Class to create PDF report.

    :param file: The file name of the pdf.
    :param plot_fn: Function to create the plot for a single data entry. The
        expected interface is func(data_entry, axis, **kwargs).
    :param plot_per_row: If True the plot function will be called for an entire
        row an not for individual subplots.
    :param max_pages: The maximal number of pages to create.
    :param rows: Plot rows per page.
    :param cols: Plot columns per page.
    :param pagesize: Size of a single page.
    :param width_ratios: Column width ratios.
    """

    def __init__(
        self,
        file,
        plot_fn,
        plot_per_row=False,
        max_pages=999,
        rows=7,
        cols=2,
        pagesize=(21, 29.7),
        width_ratios=None,
    ):
        self.pdf_pages = PdfPages(file)
        self.plot_fn = plot_fn
        self.plot_per_row = plot_per_row
        self.max_pages = max_pages
        self.rows = rows
        self.cols = cols
        self.pagesize = pagesize
        self.width_ratios = width_ratios

    def plot(self, data, **kwargs):
        if not isinstance(data, list):
            data = [data]
        if self.plot_per_row:
            pp_cnt = self.rows
        else:
            pp_cnt = self.rows * self.cols
        n_print = self.max_pages * pp_cnt
        step = np.max([len(data) / n_print, 1])
        pages = int(np.ceil(len(data) / step / pp_cnt))
        done = False
        for p in tqdm.tqdm(range(pages)):
            # print("Pring page {} of {}".format(p + 1, pages))
            fig, ax = plt.subplots(
                self.rows,
                self.cols,
                figsize=self.pagesize,
                width_ratios=self.width_ratios,
            )
            for r in range(self.rows):
                if self.plot_per_row:
                    _idx = int((p * self.rows + r) * step)
                    if _idx >= len(data):
                        done = True
                        break
                    entry = data[_idx]
                    self.plot_fn(entry, ax[r, :], index=_idx, **kwargs)
                else:
                    for c in range(self.cols):
                        _idx = int(
                            (p * self.rows * self.cols + r * self.cols + c) * step
                        )
                        if _idx >= len(data):
                            done = True
                            break
                        entry = data[_idx]
                        self.plot_fn(entry, ax[r, c], index=_idx, **kwargs)
            plt.tight_layout()
            self.pdf_pages.savefig(fig, bbox_inches="tight", pad_inches=1)
            plt.close()
            if done:
                break

    def close(self):
        self.pdf_pages.close()
