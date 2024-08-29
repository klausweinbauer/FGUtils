import io
import numpy as np
import networkx as nx

import rdkit.Chem as Chem
import rdkit.Chem.rdmolfiles as rdmolfiles
import rdkit.Chem.rdDepictor as rdDepictor
import rdkit.Chem.Draw.rdMolDraw2D as rdMolDraw2D
import rdkit.Chem.rdChemReactions as rdChemReactions

from PIL import Image

from fgutils.rdkit import graph_to_mol, graph_to_smiles
from fgutils.const import SYMBOL_KEY, AAM_KEY, BOND_KEY, IS_LABELED_KEY, LABELS_KEY
from fgutils.parse import Parser
from fgutils.proxy import ProxyGraph


def _get_its_as_mol(its: nx.Graph) -> Chem.rdchem.Mol:
    _its = its.copy()
    for n in _its.nodes:
        _its.nodes[n][AAM_KEY] = n
    for u, v in _its.edges():
        _its[u][v][BOND_KEY] = 1
    return graph_to_mol(_its)


def _get_graph_as_mol(g: nx.Graph) -> Chem.rdchem.Mol:
    _g = g.copy()
    for n, d in _g.nodes(data=True):
        if d[IS_LABELED_KEY]:
            _g.nodes[n][SYMBOL_KEY] = "C"
            _g.nodes[n][IS_LABELED_KEY] = False
        _g.nodes[n][AAM_KEY] = n
    for u, v in _g.edges():
        _g[u][v][BOND_KEY] = 1
    return graph_to_mol(_g)


def plot_its(its, ax, use_mol_coords=True, title=None):
    bond_char = {None: "∅", 0: "∅", 1: "—", 2: "=", 3: "≡"}

    if use_mol_coords:
        mol = _get_its_as_mol(its)
        positions = {}
        conformer = rdDepictor.Compute2DCoords(mol)
        for i, atom in enumerate(mol.GetAtoms()):
            aam = atom.GetAtomMapNum()
            apos = mol.GetConformer(conformer).GetAtomPosition(i)
            positions[aam] = [apos.x, apos.y]
    else:
        positions = nx.spring_layout(its)

    ax.axis("equal")
    ax.axis("off")
    if title is not None:
        ax.set_title(title)

    nx.draw_networkx_edges(its, positions, edge_color="#000000", ax=ax)
    nx.draw_networkx_nodes(its, positions, node_color="#FFFFFF", node_size=500, ax=ax)

    labels = {n: "{}:{}".format(d[SYMBOL_KEY], n) for n, d in its.nodes(data=True)}
    edge_labels = {}
    for u, v, d in its.edges(data=True):
        bc1 = d[BOND_KEY][0]
        bc2 = d[BOND_KEY][1]
        if bc1 == bc2:
            continue
        if bc1 in bond_char.keys():
            bc1 = bond_char[bc1]
        if bc2 in bond_char.keys():
            bc2 = bond_char[bc2]
        edge_labels[(u, v)] = "({},{})".format(bc1, bc2)

    nx.draw_networkx_labels(its, positions, labels=labels, ax=ax)
    nx.draw_networkx_edge_labels(its, positions, edge_labels=edge_labels, ax=ax)


def plot_as_mol(g: nx.Graph, ax, use_mol_coords=True):
    bond_char = {None: "∅", 1: "—", 2: "=", 3: "≡"}

    if use_mol_coords:
        mol = graph_to_mol(g)
        positions = {}
        conformer = rdDepictor.Compute2DCoords(mol)
        for i, atom in enumerate(mol.GetAtoms()):
            aidx = atom.GetIdx()
            apos = mol.GetConformer(conformer).GetAtomPosition(i)
            positions[aidx] = [apos.x, apos.y]
    else:
        positions = nx.spring_layout(g)

    ax.axis("equal")
    ax.axis("off")

    nx.draw_networkx_edges(g, positions, edge_color="#909090", ax=ax)
    nx.draw_networkx_nodes(g, positions, node_color="#FFFFFF", node_size=500, ax=ax)

    labels = {n: "{}".format(d[SYMBOL_KEY]) for n, d in g.nodes(data=True)}
    edge_labels = {}
    for u, v, d in g.edges(data=True):
        bc = d[BOND_KEY]
        if bc in bond_char.keys():
            bc = bond_char[bc]
        edge_labels[(u, v)] = "{}".format(bc)

    nx.draw_networkx_labels(g, positions, labels=labels, ax=ax)
    nx.draw_networkx_edge_labels(g, positions, edge_labels=edge_labels, ax=ax)


def get_rxn_img(smiles):
    drawer = rdMolDraw2D.MolDraw2DCairo(1600, 900)
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
        if img.getdata()[x + y * img.size[0]] != (255, 255, 255)  # type: ignore
    ]
    rect = (
        min([x - 10 for x, _ in nonwhite_positions]),
        min([y - 10 for _, y in nonwhite_positions]),
        max([x + 10 for x, _ in nonwhite_positions]),
        max([y + 10 for _, y in nonwhite_positions]),
    )
    return img.crop(rect)


def plot_reaction(g: nx.Graph, h: nx.Graph, ax, title=None):
    ax.axis("off")
    if title is not None:
        ax.set_title(title)
    rxn_smiles = "{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h))
    ax.imshow(get_rxn_img(rxn_smiles))


class EdgeLabelFormatter:
    def __init__(self, show_single_bonds=False, rc_only=False):
        self.rc_only = rc_only
        self.show_single_bonds = show_single_bonds
        self.bond_chars = {None: "∅", 0: "∅", 1: "—", 2: "=", 3: "≡"}

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
    def __init__(self):
        pass

    def __call__(self, n, d):
        lbl = "{}".format(d[SYMBOL_KEY])
        return lbl


class AnchorNodeLabelFormatter(NodeLabelFormatter):
    def __init__(self, anchor=[], format_str="[{}]"):
        self.anchor = anchor
        self.format_str = format_str

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
        mol = _get_graph_as_mol(nx.Graph(g))
        positions = {}
        conformer = rdDepictor.Compute2DCoords(mol)
        for i, atom in enumerate(mol.GetAtoms()):
            aidx = atom.GetIdx()
            apos = mol.GetConformer(conformer).GetAtomPosition(i)
            positions[aidx] = [apos.x, apos.y]
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
        label_legend_offset=(0, -0.5),
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
        self.label_legend_offset = label_legend_offset
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

    def plot(self, graph, ax, title=None):
        positions = get_node_positions(graph, self.use_mol_coords)

        if self.show_label_legend:
            lbl_dict, graph = self.build_label_dict(
                graph, formatter=self.label_legend_formatter
            )
            leg_pos = self.get_legend_position(
                graph,
                offset=self.label_legend_offset,
                use_mol_coords=self.use_mol_coords,
            )
            ax.text(leg_pos[0], leg_pos[1], self.label_legend_formatter(lbl_dict))

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
