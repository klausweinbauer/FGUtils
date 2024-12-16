# Modified from https://github.com/klausweinbauer/AAMUtils/blob/main/aamutils/algorithm/aaming.py

import collections
import networkx as nx

from fgutils.const import SYMBOL_KEY, AAM_KEY, BOND_KEY, IDX_MAP_KEY


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
        if not ITS.has_edge(n_ITS1, n_ITS2) and n_ITS1 > 0 and n_ITS2 > 0:
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
