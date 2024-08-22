import collections
import networkx as nx

from fgutils.const import SYMBOL_KEY, AAM_KEY, BOND_KEY


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


def get_its(G: nx.Graph, H: nx.Graph) -> nx.Graph:
    """

    Get the ITS graph of reaction G \u2192 H. G and H must be molecular graphs
    with node labels 'aam' and 'symbol' and bond label 'bond'.

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
    _add_its_nodes(ITS, G, H, eta, SYMBOL_KEY)
    _add_its_edges(ITS, G, H, eta, BOND_KEY)

    return ITS
