import networkx as nx
import rdkit.Chem as Chem
import rdkit.Chem.rdmolfiles as rdmolfiles

from fgutils.const import IS_LABELED_KEY, SYMBOL_KEY, AAM_KEY, LABELS_KEY, BOND_KEY


def mol_to_graph(mol: Chem.rdchem.Mol) -> nx.Graph:
    """Convert an RDKit molecule to a graph.

    :param mol: An RDKit molecule.

    :returns: The molecule as node and edge labeled graph.
    """
    bond_order_map = {
        "SINGLE": 1,
        "DOUBLE": 2,
        "TRIPLE": 3,
        "QUADRUPLE": 4,
        "AROMATIC": 1.5,
    }
    g = nx.Graph()
    for atom in mol.GetAtoms():
        aam = atom.GetAtomMapNum()
        node_attributes = {SYMBOL_KEY: atom.GetSymbol()}
        if aam > 0:
            node_attributes[AAM_KEY] = aam
        g.add_node(atom.GetIdx(), **node_attributes)
    for bond in mol.GetBonds():
        bond_type = str(bond.GetBondType()).split(".")[-1]
        edge_attributes = {BOND_KEY: 1}
        if bond_type in bond_order_map.keys():
            edge_attributes[BOND_KEY] = bond_order_map[bond_type]
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), **edge_attributes)
    return g


def _get_rdkit_atom_sym(symbol):
    sym_map = {"c": "C", "n": "N", "b": "B", "o": "O", "p": "P", "s": "S"}
    return sym_map.get(symbol, symbol)


def graph_to_mol(g: nx.Graph, ignore_aam=False) -> Chem.rdchem.Mol:
    """Convert a graph to an RDKit molecule.

    :param g: The molecule as node and edge labeled graph. The graph requires
        ``SYMBOL`` node labels and ``BOND_KEY`` edge labels. The node label
        ``AAM_KEY`` is optional to annotate the molecule with an atom-atom map.

    :param ignore_aam: If set to true the atom-atom map will not be
        initialized.

    :returns: Returns the graph as RDKit molecule.
    """
    bond_order_map = {
        1: Chem.rdchem.BondType.SINGLE,
        2: Chem.rdchem.BondType.DOUBLE,
        3: Chem.rdchem.BondType.TRIPLE,
        4: Chem.rdchem.BondType.QUADRUPLE,
        1.5: Chem.rdchem.BondType.AROMATIC,
    }
    rw_mol = Chem.rdchem.RWMol()
    idx_map = {}
    for n, d in g.nodes(data=True):
        if d is None:
            raise ValueError("Graph node {} has no data.".format(n))
        atom_symbol = _get_rdkit_atom_sym(d[SYMBOL_KEY])
        if IS_LABELED_KEY in d.keys() and d[IS_LABELED_KEY]:
            raise ValueError(
                "Graph contains labeled nodes. Node {} with label [{}].".format(
                    n, ",".join(d[LABELS_KEY])
                )
            )
        idx = rw_mol.AddAtom(Chem.rdchem.Atom(atom_symbol))
        idx_map[n] = idx
        if not ignore_aam and AAM_KEY in d.keys() and d[AAM_KEY] >= 0:
            rw_mol.GetAtomWithIdx(idx).SetAtomMapNum(d[AAM_KEY])
    for n1, n2, d in g.edges(data=True):
        if d is None:
            raise ValueError("Graph edge {} has no data.".format((n1, n2)))
        idx1 = idx_map[n1]
        idx2 = idx_map[n2]
        rw_mol.AddBond(idx1, idx2, bond_order_map[d[BOND_KEY]])
    return rw_mol.GetMol()


def graph_to_smiles(g: nx.Graph, ignore_aam=False, canonical=True) -> str:
    """Convert a molecular graph into a SMILES string. This function uses
    RDKit for SMILES generation.

    :param g: Graph to convert to SMILES representation.
    :param ignore_aam: If set to True the returned SMILES has no atom-atom map.
    :param canonical: If set to False no attempt is made to canonicalize the
        SMILES.

    :returns: Returns the SMILES.
    """
    mol = graph_to_mol(g, ignore_aam=ignore_aam)
    return rdmolfiles.MolToSmiles(mol, canonical=canonical)


def reaction_smiles_to_graph(smiles: str) -> tuple[nx.Graph, nx.Graph]:
    """Converts a reaction SMILES to the graph representation G \u2192 H,
    where G is the reactant graph and H is the product graph.

    :param smiles: Reaction SMILES to convert to graph tuple.

    :returns: Returns the graphs G and H as tuple.
    """
    rxn_tokens = smiles.split(">>")
    if len(rxn_tokens) != 2:
        raise ValueError("Expected reaction SMILES but found '{}'.".format(smiles))
    r_smiles, p_smiles = rxn_tokens
    g = smiles_to_graph(r_smiles)
    h = smiles_to_graph(p_smiles)
    assert isinstance(g, nx.Graph)
    assert isinstance(h, nx.Graph)
    return g, h


def mol_smiles_to_graph(smiles: str) -> nx.Graph:
    """Converts a SMILES to a graph.

    :param smiles: SMILES to convert to graph(s).

    :returns: A node and edge labeled molecular graph.
    """
    mol = rdmolfiles.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit was unable to parse SMILES '{}'.".format(smiles))
    return mol_to_graph(mol)


def smiles_to_graph(smiles: str) -> nx.Graph | tuple[nx.Graph, nx.Graph]:
    """Converts a SMILES to a graph. If the SMILES encodes a reaction a graph
    tuple is returned.

    :param smiles: SMILES to convert to graph(s).

    :returns: A molecular graph or graph tuple if SMILES is a reaction SMILES.
    """
    if ">>" in smiles:
        return reaction_smiles_to_graph(smiles)
    else:
        return mol_smiles_to_graph(smiles)
