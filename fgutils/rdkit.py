import networkx as nx
import rdkit.Chem as Chem
import rdkit.Chem.rdmolfiles as rdmolfiles

from fgutils.const import IS_LABELED_KEY, SYMBOL_KEY, AAM_KEY, LABELS_KEY, BOND_KEY


def mol_to_graph(mol: Chem.rdchem.Mol) -> nx.Graph:
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


def graph_to_mol(g: nx.Graph) -> Chem.rdchem.Mol:
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
        if AAM_KEY in d.keys() and d[AAM_KEY] >= 0:
            rw_mol.GetAtomWithIdx(idx).SetAtomMapNum(d[AAM_KEY])
    for n1, n2, d in g.edges(data=True):
        if d is None:
            raise ValueError("Graph edge {} has no data.".format((n1, n2)))
        idx1 = idx_map[n1]
        idx2 = idx_map[n2]
        rw_mol.AddBond(idx1, idx2, bond_order_map[d[BOND_KEY]])
    return rw_mol.GetMol()


def graph_to_smiles(g: nx.Graph) -> str:
    mol = graph_to_mol(g)
    return rdmolfiles.MolToSmiles(mol)


def smiles_to_graph(smiles: str) -> nx.Graph | tuple[nx.Graph, nx.Graph]:
    if ">>" in smiles:
        r_smiles, p_smiles = smiles.split(">>")
        g = smiles_to_graph(r_smiles)
        h = smiles_to_graph(p_smiles)
        assert isinstance(g, nx.Graph)
        assert isinstance(h, nx.Graph)
        return g, h
    else:
        mol = rdmolfiles.MolFromSmiles(smiles)
        return mol_to_graph(mol)
