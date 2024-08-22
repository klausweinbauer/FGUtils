from fgutils.proxy_collection.diels_alder_proxy import OutOfSampleError
from fgutils.proxy_collection import DielsAlderProxy
from fgutils.rdkit import graph_to_smiles, graph_to_mol

import io
import rdkit.Chem.rdmolfiles as rdmolfiles
import rdkit.Chem.Draw.rdMolDraw2D as rdMolDraw2D
import rdkit.Chem.rdChemReactions as rdChemReactions
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

from PIL import Image
from synrbl.SynUtils.chem_utils import remove_atom_mapping, normalize_smiles

import collections
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import rdkit.Chem as Chem
import rdkit.Chem.rdmolfiles as rdmolfiles
import rdkit.Chem.rdDepictor as rdDepictor


def _add_its_nodes(ITS, G, H, eta, symbol_key):
    eta_G, eta_G_inv, eta_H, eta_H_inv = eta[0], eta[1], eta[2], eta[3]
    for n, d in G.nodes(data=True):
        n_ITS = eta_G[n]
        n_H = eta_H_inv[n_ITS]
        if n_ITS is not None and n_H is not None:
            ITS.add_node(n_ITS, symbol=d[symbol_key], idx_map=(n, n_H))
    for n, d in H.nodes(data=True):
        n_ITS = eta_H[n]
        n_G = eta_G_inv[n_ITS]
        if n_ITS is not None and n_G is not None and n_ITS not in ITS.nodes:
            ITS.add_node(n_ITS, symbol=d[symbol_key], idx_map=(n_G, n))


def _add_its_edges(ITS, G, H, eta, bond_key):
    eta_G, eta_G_inv, eta_H, eta_H_inv = eta[0], eta[1], eta[2], eta[3]
    for n1, n2, d in G.edges(data=True):
        if n1 > n2:
            continue
        e_G = d[bond_key]
        n_ITS1 = eta_G[n1]
        n_ITS2 = eta_G[n2]
        n_H1 = eta_H_inv[n_ITS1]
        n_H2 = eta_H_inv[n_ITS2]
        e_H = None
        if H.has_edge(n_H1, n_H2):
            e_H = H[n_H1][n_H2][bond_key]
        if not ITS.has_edge(n_ITS1, n_ITS2) and n_ITS1 > 0 and n_ITS2 > 0:
            ITS.add_edge(n_ITS1, n_ITS2, bond=(e_G, e_H))

    for n1, n2, d in H.edges(data=True):
        if n1 > n2:
            continue
        e_H = d[bond_key]
        n_ITS1 = eta_H[n1]
        n_ITS2 = eta_H[n2]
        n_G1 = eta_G_inv[n_ITS1]
        n_G2 = eta_G_inv[n_ITS2]
        if n_G1 is None or n_G2 is None:
            continue
        if not G.has_edge(n_G1, n_G2) and n_ITS1 > 0 and n_ITS2 > 0:
            ITS.add_edge(n_ITS1, n_ITS2, bond=(None, e_H))


def get_its(
    G: nx.Graph, H: nx.Graph, aam_key="aam", symbol_key="symbol", bond_key="bond"
) -> nx.Graph:
    """Get the ITS graph of reaction G \u2192 H.

    :param G: Educt molecular graph.
    :param H: Product molecular graph.
    :param aam_key: (optional) The node label that encodes atom indices on the
        molecular graph.
    :param symbol_key: (optional) The node label that encodes atom symbols on
        the molecular graph.
    :param bond_key: (optional) The edge label that encodes bond order on the
        molecular graph.

    :returns: Returns the ITS graph.
    """
    eta_G = collections.defaultdict(lambda: None)
    eta_G_inv = collections.defaultdict(lambda: None)
    eta_H = collections.defaultdict(lambda: None)
    eta_H_inv = collections.defaultdict(lambda: None)
    eta = (eta_G, eta_G_inv, eta_H, eta_H_inv)

    for n, d in G.nodes(data=True):
        if aam_key in d.keys() and d[aam_key] >= 0:
            eta_G[n] = d[aam_key]
            eta_G_inv[d[aam_key]] = n
    for n, d in H.nodes(data=True):
        if aam_key in d.keys() and d[aam_key] >= 0:
            eta_H[n] = d[aam_key]
            eta_H_inv[d[aam_key]] = n

    ITS = nx.Graph()
    _add_its_nodes(ITS, G, H, eta, symbol_key)
    _add_its_edges(ITS, G, H, eta, bond_key)

    return ITS


def its2mol(its: nx.Graph, aam_key="aam", bond_key="bond") -> Chem.rdchem.Mol:
    _its = its.copy()
    for n in _its.nodes:
        _its.nodes[n][aam_key] = n
    for u, v in _its.edges():
        _its[u][v][bond_key] = 1
    return graph_to_mol(_its)


def plot_its(
    its, ax, bond_key="bond", aam_key="aam", symbol_key="symbol", use_mol_coords=True
):
    bond_char = {None: "∅", 1: "—", 2: "=", 3: "≡"}
    mol = its2mol(its, aam_key=aam_key, bond_key=bond_key)

    if use_mol_coords:
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

    labels = {n: "{}:{}".format(d[symbol_key], n) for n, d in its.nodes(data=True)}
    edge_labels = {}
    for u, v, d in its.edges(data=True):
        bc1 = d[bond_key][0]
        bc2 = d[bond_key][1]
        if bc1 == bc2:
            continue
        if bc1 in bond_char.keys():
            bc1 = bond_char[bc1]
        if bc2 in bond_char.keys():
            bc2 = bond_char[bc2]
        edge_labels[(u, v)] = "({},{})".format(bc1, bc2)

    nx.draw_networkx_labels(its, positions, labels=labels, ax=ax)
    nx.draw_networkx_edge_labels(its, positions, edge_labels=edge_labels, ax=ax)

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
        min([x - 10 for x, y in nonwhite_positions]),
        min([y - 10 for x, y in nonwhite_positions]),
        max([x + 10 for x, y in nonwhite_positions]),
        max([y + 10 for x, y in nonwhite_positions]),
    )
    return img.crop(rect)

proxy = DielsAlderProxy(enable_aam=True)

n = 6
fig, ax = plt.subplots(n, 2, dpi=50, figsize=(16, 9))
for i in range(n):
    try:
        g, h = next(proxy)
        rxn_smiles = "{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h))
        print(rxn_smiles)
        its = get_its(g, h)
        ax[i, 0].imshow(get_rxn_img(rxn_smiles))
        ax[i, 0].axis("off")
        plot_its(its, ax[i, 1])
    except OutOfSampleError:
        break
plt.tight_layout()
plt.show()
