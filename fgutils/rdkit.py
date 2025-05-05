import networkx as nx
import rdkit.Chem as Chem
import rdkit.Chem.rdmolfiles as rdmolfiles
import rdkit.Chem.rdDepictor as rdDepictor
import rdkit.Chem.rdmolops as rdmolops

from fgutils.utils import to_non_aromatic_symbol
from fgutils.const import (
    IS_LABELED_KEY,
    SYMBOL_KEY,
    AAM_KEY,
    LABELS_KEY,
    BOND_KEY,
    BOND_ORDER_MAP,
)


RDKIT_BOND_ORDER_MAP = {
    1: Chem.rdchem.BondType.SINGLE,
    2: Chem.rdchem.BondType.DOUBLE,
    3: Chem.rdchem.BondType.TRIPLE,
    4: Chem.rdchem.BondType.QUADRUPLE,
    1.5: Chem.rdchem.BondType.AROMATIC,
}


def mol_to_graph(mol: Chem.rdchem.Mol, implicit_h=False) -> nx.Graph:
    """Convert an RDKit molecule to a graph.

    :param mol: An RDKit molecule.
    :param implicit_h: Flag to add all Hydrogen atoms. (Default: False)

    :returns: The molecule as node and edge labeled graph.
    """
    g = nx.Graph()
    h_idx = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        aam = atom.GetAtomMapNum()
        node_attributes = {SYMBOL_KEY: atom.GetSymbol()}
        if aam > 0:
            node_attributes[AAM_KEY] = aam
        atom_idx = atom.GetIdx()
        g.add_node(atom_idx, **node_attributes)

        # Add hydrogens
        h_cnt = atom.GetNumExplicitHs()
        if implicit_h:
            h_cnt += atom.GetNumImplicitHs()
        for _ in range(h_cnt):
            h_attributes = {SYMBOL_KEY: "H"}
            h_edge_attributes = {BOND_KEY: 1}
            g.add_node(h_idx, **h_attributes)
            g.add_edge(atom_idx, h_idx, **h_edge_attributes)
            h_idx += 1

        # Add Charge
        charge = atom.GetFormalCharge()
        if charge != 0:
            edge_attributes = {BOND_KEY: 0.5 * -charge}
            g.add_edge(atom_idx, atom_idx, **edge_attributes)

    for bond in mol.GetBonds():
        bond_type = str(bond.GetBondType()).split(".")[-1]
        edge_attributes = {BOND_KEY: 1}
        if bond_type in BOND_ORDER_MAP.keys():
            edge_attributes[BOND_KEY] = BOND_ORDER_MAP[bond_type]
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), **edge_attributes)

    return g


def get_mol_coords(g: nx.Graph, scale=1) -> dict[int, tuple[float, float]]:
    """Try to get a molecule like coordinate representation of the graph.

    :param g: The graph to get the coordinates for.
    :param scale: (optional) A scale for the coordinates. (Default: 1)

    :returns: Returns a dict of coordinates. The keys are the node indices and
        the values are the 2 coordinates x and y.
    """
    _g = g.copy()
    for n, d in _g.nodes(data=True):
        if IS_LABELED_KEY in d and d[IS_LABELED_KEY]:
            _g.nodes[n][SYMBOL_KEY] = "C"
            _g.nodes[n][IS_LABELED_KEY] = False
        _g.nodes[n][AAM_KEY] = n
    for u, v in _g.edges():
        _g[u][v][BOND_KEY] = 1
    positions = {}
    mol = graph_to_mol(_g)
    conformer = rdDepictor.Compute2DCoords(mol)
    for i, atom in enumerate(mol.GetAtoms()):
        aam = atom.GetAtomMapNum()
        apos = mol.GetConformer(conformer).GetAtomPosition(i)
        positions[aam] = [scale * apos.x, scale * apos.y]
    return positions


def _get_node_H_count(g, v):
    H_cnt = 0
    for n in g.neighbors(v):
        d = g.nodes[n]
        if d[SYMBOL_KEY] == "H":
            H_cnt += 1
    return H_cnt


def _graph_to_smiles_node_check(n, d):
    if d is None:
        raise ValueError("Graph node {} has no data.".format(n))

    if IS_LABELED_KEY in d.keys() and d[IS_LABELED_KEY]:
        raise ValueError(
            "Graph contains labeled nodes. Node {} with label [{}].".format(
                n, ",".join(d[LABELS_KEY])
            )
        )


def graph_to_mol(g: nx.Graph, ignore_aam=False) -> Chem.rdchem.Mol:
    """Convert a graph to an RDKit molecule.

    :param g: The molecule as node and edge labeled graph. The graph requires
        ``SYMBOL`` node labels and ``BOND_KEY`` edge labels. The node label
        ``AAM_KEY`` is optional to annotate the molecule with an atom-atom map.

    :param ignore_aam: If set to true the atom-atom map will not be
        initialized.

    :returns: Returns the graph as RDKit molecule.
    """
    rw_mol = Chem.rdchem.RWMol()
    idx_map = {}
    H_nodes = set()
    for n, d in g.nodes(data=True):
        _graph_to_smiles_node_check(n, d)
        atom_symbol = to_non_aromatic_symbol(d[SYMBOL_KEY])
        if atom_symbol == "H":
            H_nodes.add(n)
            continue
        atom = Chem.rdchem.Atom(atom_symbol)
        atom.SetNumExplicitHs(_get_node_H_count(g, n))
        idx = rw_mol.AddAtom(atom)
        idx_map[n] = idx
        if not ignore_aam and AAM_KEY in d.keys() and d[AAM_KEY] >= 0:
            rw_mol.GetAtomWithIdx(idx).SetAtomMapNum(d[AAM_KEY])

    for n1, n2, d in g.edges(data=True):
        if n1 in H_nodes or n2 in H_nodes:
            continue
        if d is None:
            raise ValueError("Graph edge {} has no data.".format((n1, n2)))
        idx1 = idx_map[n1]
        idx2 = idx_map[n2]
        if n1 == n2:
            rw_mol.GetAtomWithIdx(idx_map[n1]).SetFormalCharge(int(-2 * d[BOND_KEY]))
        else:
            rw_mol.AddBond(idx1, idx2, RDKIT_BOND_ORDER_MAP[d[BOND_KEY]])

    mol = rw_mol.GetMol()
    rdmolops.SanitizeMol(mol)
    return mol


def graph_to_smiles(g: nx.Graph, implicit_h=False, ignore_aam=False) -> str:
    """Convert a molecular graph into a SMILES string. This function uses
    RDKit for SMILES generation.

    :param g: Graph to convert to SMILES representation.
    :param implicit_h: Flag to add all Hydrogen atoms. (Default: False)
    :param ignore_aam: If set to True the returned SMILES has no atom-atom map.

    :returns: Returns the SMILES.
    """
    mol = graph_to_mol(g, ignore_aam=ignore_aam)
    return rdmolfiles.MolToSmiles(mol, allHsExplicit=implicit_h)


def reaction_smiles_to_graph(
    smiles: str, implicit_h=False
) -> tuple[nx.Graph, nx.Graph]:
    """Converts a reaction SMILES to the graph representation G \u2192 H,
    where G is the reactant graph and H is the product graph.

    :param smiles: Reaction SMILES to convert to graph tuple.
    :param implicit_h: Flag to add all Hydrogen atoms. (Default: False)

    :returns: Returns the graphs G and H as tuple.
    """
    rxn_tokens = smiles.split(">>")
    if len(rxn_tokens) != 2:
        raise ValueError("Expected reaction SMILES but found '{}'.".format(smiles))
    r_smiles, p_smiles = rxn_tokens
    g = smiles_to_graph(r_smiles, implicit_h=implicit_h)
    h = smiles_to_graph(p_smiles, implicit_h=implicit_h)
    assert isinstance(g, nx.Graph)
    assert isinstance(h, nx.Graph)
    return g, h


def mol_smiles_to_graph(smiles: str, implicit_h=False) -> nx.Graph:
    """Converts a SMILES to a graph.

    :param smiles: SMILES to convert to graph(s).
    :param implicit_h: Flag to add all Hydrogen atoms. (Default: False)

    :returns: A node and edge labeled molecular graph.
    """
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = rdmolfiles.MolFromSmiles(smiles, params)
    if mol is None:
        raise ValueError("RDKit was unable to parse SMILES '{}'.".format(smiles))
    return mol_to_graph(mol, implicit_h=implicit_h)


def smiles_to_graph(
    smiles: str, implicit_h=False
) -> nx.Graph | tuple[nx.Graph, nx.Graph]:
    """Converts a SMILES to a graph. If the SMILES encodes a reaction a graph
    tuple is returned.

    :param smiles: SMILES to convert to graph(s).
    :param implicit_h: Flag to add all Hydrogen atoms. (Default: False)

    :returns: A molecular graph or graph tuple if SMILES is a reaction SMILES.
    """
    if ">>" in smiles:
        return reaction_smiles_to_graph(smiles, implicit_h=implicit_h)
    else:
        return mol_smiles_to_graph(smiles, implicit_h=implicit_h)
