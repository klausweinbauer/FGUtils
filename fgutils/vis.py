import io
import numpy as np
import networkx as nx
from numpy import who

import rdkit.Chem as Chem
import rdkit.Chem.rdmolfiles as rdmolfiles
import rdkit.Chem.rdDepictor as rdDepictor
import rdkit.Chem.Draw.rdMolDraw2D as rdMolDraw2D
import rdkit.Chem.rdChemReactions as rdChemReactions
import matplotlib.pyplot as plt

from PIL import Image
from matplotlib.backends.backend_pdf import PdfPages

from fgutils.rdkit import graph_to_mol, graph_to_smiles
from fgutils.const import SYMBOL_KEY, AAM_KEY, BOND_KEY, IS_LABELED_KEY, LABELS_KEY
from fgutils.parse import Parser
from fgutils.proxy import ProxyGroup


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


def plot_its(its, ax, use_mol_coords=True):
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


class AutoEdgeLabelFormatter:
    def __init__(self, rc_only=False):
        self.rc_only = rc_only
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
            if bc1 != bc2:
                return "({},{})".format(bc1, bc2)
            elif not self.rc_only:
                return "{}".format(bc1)
            else:
                return ""
        else:
            if bc in self.bond_chars.keys():
                bc = self.bond_chars[bc]
            return "{}".format(bc)


def plot_graph(
    g: nx.Graph,
    ax,
    use_mol_coords=True,
    show_labels=False,
    show_edge_labels=True,
    title=None,
    fmt_node_label=None,
    fmt_edge_label=None,
):
    if fmt_node_label is None:
        fmt_node_label = lambda n, d: "{}".format(d[SYMBOL_KEY])
    if fmt_edge_label is None:
        fmt_edge_label = AutoEdgeLabelFormatter()

    if use_mol_coords:
        mol = _get_graph_as_mol(g)
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

    if title is not None:
        ax.set_title(title)

    nx.draw_networkx_edges(g, positions, edge_color="#909090", ax=ax)
    nx.draw_networkx_nodes(g, positions, node_color="#FFFFFF", node_size=500, ax=ax)

    labels = {}  # {n: "{}".format(d[SYMBOL_KEY]) for n, d in g.nodes(data=True)}
    for n, d in g.nodes(data=True):
        lbl = fmt_node_label(n, d)
        if d[IS_LABELED_KEY] and show_labels:
            lbl = "{}".format(d[LABELS_KEY])
        labels[n] = lbl

    edge_labels = {}
    for u, v, d in g.edges(data=True):
        edge_labels[(u, v)] = fmt_edge_label((u, v), d)

    nx.draw_networkx_labels(g, positions, labels=labels, ax=ax)
    if show_edge_labels:
        nx.draw_networkx_edge_labels(g, positions, edge_labels=edge_labels, ax=ax)


def plot_reaction(g: nx.Graph, h: nx.Graph, ax):
    rxn_smiles = "{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h))
    ax.imshow(get_rxn_img(rxn_smiles))


def create_pdf_group_report(file: str, groups: list[ProxyGroup], parser: Parser = None):
    def _fmtnode(n, d):
        lbl = ""
        if n in graph.anchor:
            lbl = "[{}]".format(d[SYMBOL_KEY])
        else:
            lbl = "{}".format(d[SYMBOL_KEY])
        # if d[IS_LABELED_KEY]:
        #     lbl = "{}".format(d[LABELS_KEY])
        return lbl

    if parser is None:
        parser = Parser()
    pp = PdfPages(file)
    rows = 6
    cols = 4
    p_idx = 0
    figsize = (21, 29.7)
    dpi = 100
    fig, axs = plt.subplots(rows, cols, figsize=figsize, dpi=dpi)
    for group in groups:
        for i, graph in enumerate(group.graphs):
            ax = axs[int(p_idx / cols), p_idx % cols]
            g = parser(graph.pattern)
            plot_graph(
                g,
                ax,
                fmt_node_label=_fmtnode,
                fmt_edge_label=AutoEdgeLabelFormatter(rc_only=True),
            )
            ax.axis("off")
            ax.set_title("{} Graph {}".format(group.name, i + 1))
            p_idx += 1
            if p_idx == rows * cols:
                p_idx = 0
                plt.tight_layout()
                pp.savefig(fig)
                fig, axs = plt.subplots(rows, cols, figsize=figsize, dpi=dpi)
    pp.close()
