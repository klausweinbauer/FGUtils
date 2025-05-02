# Modified from https://github.com/klausweinbauer/AAMUtils/blob/main/aamutils/algorithm/aaming.py

import collections
import networkx as nx

from fgutils.const import SYMBOL_KEY, AAM_KEY, BOND_KEY, IDX_MAP_KEY
from fgutils.rdkit import smiles_to_graph, graph_to_smiles
from fgutils.utils import complete_aam, get_unreachable_nodes, relabel_graph


def _add_its_nodes(ITS, G, H, eta):
    eta_G, eta_G_inv, eta_H, eta_H_inv = eta[0], eta[1], eta[2], eta[3]
    for n, d in G.nodes(data=True):
        n_ITS = eta_G[n]
        n_H = eta_H_inv[n_ITS]
        if n_ITS is not None and n_H is not None:
            node_attributes = {
                SYMBOL_KEY: d[SYMBOL_KEY],
                IDX_MAP_KEY: (n, n_H),
                AAM_KEY: n_ITS,
            }
            ITS.add_node(n_ITS, **node_attributes)
    for n, d in H.nodes(data=True):
        n_ITS = eta_H[n]
        n_G = eta_G_inv[n_ITS]
        if n_ITS is not None and n_G is not None and n_ITS not in ITS.nodes:
            node_attributes = {
                SYMBOL_KEY: d[SYMBOL_KEY],
                IDX_MAP_KEY: (n_G, n),
                AAM_KEY: n_ITS,
            }
            ITS.add_node(n_ITS, **node_attributes)


def _add_its_edges(ITS, G, H, eta):
    eta_G, eta_G_inv, eta_H, eta_H_inv = eta[0], eta[1], eta[2], eta[3]
    for n1, n2, d in G.edges(data=True):
        if n1 > n2:
            continue
        e_G = d[BOND_KEY]
        n_ITS1 = eta_G[n1]
        n_ITS2 = eta_G[n2]
        n_H1 = eta_H_inv[n_ITS1]
        n_H2 = eta_H_inv[n_ITS2]
        e_H = 0
        if H.has_edge(n_H1, n_H2):
            e_H = H[n_H1][n_H2][BOND_KEY]
        if (
            not ITS.has_edge(n_ITS1, n_ITS2)
            and n_ITS1 is not None
            and n_ITS1 > 0
            and n_ITS2 is not None
            and n_ITS2 > 0
        ):
            edge_attributes = {BOND_KEY: (e_G, e_H)}
            ITS.add_edge(n_ITS1, n_ITS2, **edge_attributes)

    for n1, n2, d in H.edges(data=True):
        if n1 > n2:
            continue
        e_H = d[BOND_KEY]
        n_ITS1 = eta_H[n1]
        n_ITS2 = eta_H[n2]
        n_G1 = eta_G_inv[n_ITS1]
        n_G2 = eta_G_inv[n_ITS2]
        if n_G1 is None or n_G2 is None:
            continue
        if not G.has_edge(n_G1, n_G2) and n_ITS1 > 0 and n_ITS2 > 0:
            edge_attributes = {BOND_KEY: (0, e_H)}
            ITS.add_edge(n_ITS1, n_ITS2, **edge_attributes)


def get_its(G: nx.Graph, H: nx.Graph) -> nx.Graph:
    """Get the ITS graph of reaction G \u2192 H. G and H must be molecular
    graphs with node labels 'aam' and 'symbol' and bond label 'bond'.

    :param G: Reactant molecular graph.
    :param H: Product molecular graph.

    :returns: Returns the ITS graph.
    """
    eta_G = collections.defaultdict(lambda: None)
    eta_G_inv = collections.defaultdict(lambda: None)
    eta_H = collections.defaultdict(lambda: None)
    eta_H_inv = collections.defaultdict(lambda: None)
    eta = (eta_G, eta_G_inv, eta_H, eta_H_inv)

    for n, d in G.nodes(data=True):
        if d is None:
            raise ValueError("Graph node {} has no data.".format(n))
        if AAM_KEY in d.keys() and d[AAM_KEY] >= 0:
            eta_G[n] = d[AAM_KEY]
            eta_G_inv[d[AAM_KEY]] = n
    for n, d in H.nodes(data=True):
        if d is None:
            raise ValueError("Graph node {} has no data.".format(n))
        if AAM_KEY in d.keys() and d[AAM_KEY] >= 0:
            eta_H[n] = d[AAM_KEY]
            eta_H_inv[d[AAM_KEY]] = n

    ITS = nx.Graph()
    _add_its_nodes(ITS, G, H, eta)
    _add_its_edges(ITS, G, H, eta)

    return ITS


def get_rc(ITS: nx.Graph) -> nx.Graph:
    """Get the reaction center (RC) graph from an ITS graph.

    :param ITS: The ITS graph to get the RC from.
    :param symbol_key: (optional) The node label that encodes atom symbols on
        the molecular graph.
    :param bond_key: (optional) The edge label that encodes bond order on the
        molecular graph.

    :returns: Returns the reaction center (RC) graph.
    """
    rc = nx.Graph()
    for n1, n2, d in ITS.edges(data=True):
        edge_label = d[BOND_KEY]
        if edge_label[0] != edge_label[1]:
            rc.add_node(n1, **{SYMBOL_KEY: ITS.nodes[n1][SYMBOL_KEY]})
            rc.add_node(n2, **{SYMBOL_KEY: ITS.nodes[n2][SYMBOL_KEY]})
            rc.add_edge(n1, n2, **{BOND_KEY: edge_label})
    return rc


def split_its(graph: nx.Graph) -> tuple[nx.Graph, nx.Graph]:
    """Split an ITS graph into reactant graph G and product graph H. Required
    labels on the ITS graph are BOND_KEY.

    :param graph: ITS graph to split up.

    :returns: Tuple of two graphs (G, H).
    """

    def _set_rc_edge(g, u, v, b):
        if b == 0:
            g.remove_edge(u, v)
        else:
            g[u][v][BOND_KEY] = b

    g = graph.copy()
    h = graph.copy()
    for u, v, d in graph.edges(data=True):  # type: ignore
        if d is None:
            raise ValueError("No edge labels found.")
        bond = d[BOND_KEY]
        if isinstance(bond, tuple) or isinstance(bond, list):
            _set_rc_edge(g, u, v, bond[0])
            _set_rc_edge(h, u, v, bond[1])
    return g, h


def prune_its_to_rc(its, radius=0, insert_hydrogens=True):
    """Prune an ITS graph to its reaction center and some context. The context
    is defined by a radius. Radius 0 gives only the reaction center.

    :param its: The ITS graph to prune.
    :param radius: The size of the context around the reaction center.
        (Default: 0)
    :param insert_hydrogens: If true, removed nodes will be replaced by
        hydrogen atoms if were adjacent to kept nodes. This ensures that the
        result is still a valid molecular graph. (Default: True)

    :returns: Returns the pruned ITS graph.
    """
    rc = get_rc(its)
    unreachable_nodes = get_unreachable_nodes(its, rc.nodes, radius=radius)
    its_pruned = its.copy()
    new_node_id = len(its.nodes)
    for u in unreachable_nodes:
        if insert_hydrogens:
            for v in its.neighbors(u):
                if v not in unreachable_nodes:
                    its_pruned.add_node(new_node_id, **{SYMBOL_KEY: "H"})
                    its_pruned.add_edge(new_node_id, v, **{BOND_KEY: (1, 1)})
                    new_node_id += 1
        its_pruned.remove_node(u)
    return its_pruned


def remove_reagents(its):
    """Remove all reagents from the ITS graph. These are all compounds that
    are not connected to the reaction center.

    :returns: Returns the new ITS graph.
    """
    for component_nodes in sorted(nx.connected_components(its), key=len, reverse=True):
        component = its.subgraph(component_nodes).copy()
        for _, _, d in component.edges(data=True):
            if d["bond"][0] != d["bond"][1]:
                return relabel_graph(component)
    raise RuntimeError("No reaction center found.")


class ITS:
    """Imaginary Transition State graph class. Superposition graph of
    reactants and products in a chemical reaction.

    :param graph: The raw ITS graph.
    """

    def __init__(self, graph: nx.Graph):
        self.graph = graph
        complete_aam(graph, offset="min")

    @classmethod
    def from_smiles(cls, smiles: str):
        """Construct an ITS graph from an atom-atom mapped reaction smiles.

        :param smiles: An atom-atom mapped reaction smiles.

        :returns: Returns the ITS graph of the reaction.
        """
        g, h = smiles_to_graph(smiles)
        its = get_its(g, h)
        return cls(its)

    def to_smiles(self, ignore_aam=False, implicit_h=False) -> str:
        """Convert the ITS graph into a reaction smiles.

        :param ignore_aam: If set to True the returned SMILES has no atom-atom
            map.
        :param implicit_h: Flag to add all Hydrogen atoms. (Default: False)

        :returns: Returns the reaction smiles.
        """
        g, h = split_its(self.graph)
        smiles = "{}>>{}".format(
            graph_to_smiles(g, ignore_aam=ignore_aam, implicit_h=implicit_h),
            graph_to_smiles(h, ignore_aam=ignore_aam, implicit_h=implicit_h),
        )
        return smiles

    def split(self) -> tuple[nx.Graph, nx.Graph]:
        """Split the ITS graph into reactant graph G and product graph H.

        :returns: Returns G and H as tuple.
        """
        g, h = split_its(self.graph)
        return g, h

    def prune(self, radius=1, insert_hydrogens=True):
        """Prune the ITS graph to its reaction center and some context. The
        context is defined by a radius. Radius 0 gives only the reaction
        center.

        :param radius: The size of the context around the reaction center.
            (Default: 0)
        :param insert_hydrogens: If true, removed nodes will be replaced by
            hydrogen atoms if were adjacent to kept nodes. This ensures that the
            result is still a valid molecular graph. (Default: True)
        """
        self.graph = prune_its_to_rc(
            self.graph, radius=radius, insert_hydrogens=insert_hydrogens
        )

    def remove_reagents(self):
        """Remove all reagents from the ITS graph. These are all compounds that
        are not connected to the reaction center.

        :returns: Returns the new ITS graph.
        """
        self.graph = remove_reagents(self.graph)
