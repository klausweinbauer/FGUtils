import networkx as nx
import rdkit.Chem as Chem
import rdkit.Chem.rdmolfiles as rdmolfiles


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
        g.add_node(atom.GetIdx(), symbol=atom.GetSymbol())
    for bond in mol.GetBonds():
        bond_type = str(bond.GetBondType()).split(".")[-1]
        bond_order = 1
        if bond_type in bond_order_map.keys():
            bond_order = bond_order_map[bond_type]
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond=bond_order)
    return g


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
        idx = rw_mol.AddAtom(Chem.rdchem.Atom(d["symbol"]))
        idx_map[n] = idx
        if "aam" in d.keys() and d["aam"] >= 0:
            rw_mol.GetAtomWithIdx(idx).SetAtomMapNum(d["aam"])
    for n1, n2, d in g.edges(data=True):
        idx1 = idx_map[n1]
        idx2 = idx_map[n2]
        rw_mol.AddBond(idx1, idx2, bond_order_map[d["bond"]])
    return rw_mol.GetMol()


def graph_to_smiles(g: nx.Graph) -> str:
    mol = graph_to_mol(g)
    return rdmolfiles.MolToSmiles(mol)


def smiles_to_graph(smiles: str) -> nx.Graph:
    mol = rdmolfiles.MolFromSmiles(smiles)
    return mol_to_graph(mol)
