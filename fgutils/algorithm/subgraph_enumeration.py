# This algorithm implementation is based on the paper 'A linear delay algorithm
# for enumerating all connected induced subgraphs' by Mohammed Alokshiya, Saeed
# Salem and Fidaa Abed
# https://doi.org/10.1186/s12859-019-2837-y

import copy
import networkx as nx
import numpy as np


def is_valid_extension(U, v, D):
    s = U[0]  # anchor(U)
    x = U[-1]  # last_added(U)
    if v < s:
        return False
    if D[v] > D[x]:  # distance(s, v) > distance(s, x):
        return True
    return D[v] == D[x] and v > x  # distance(s, v) == distance(s, x) and v > x


def enumerateCIS(G, U, C, D):
    yield U
    for v in C:
        if is_valid_extension(U, v, D):
            new_C = [u for u in G.neighbors(v) if u not in C and u not in U]
            _D = copy.deepcopy(D)
            for u in new_C:
                assert (
                    D[v] + 1 <= _D[u]
                ), "Wanted to increase dist @ node {} from {} to {}.".format(
                    u, _D[u], D[v] + 1
                )
                _D[u] = D[v] + 1
            yield from enumerateCIS(G, U + [v], C + new_C, _D)


def _node_induced_connected_subgraphs(G, anchor):
    U = [anchor]
    C = list(G.neighbors(anchor))
    D = [np.inf] * len(G.nodes)
    P = [-np.inf] * len(G.nodes)
    D[anchor] = 0
    for c in C:
        D[c] = 1
    P[anchor] = -1
    return enumerateCIS(G, U, C, D)


def node_induced_connected_subgraphs(G, anchor):
    """Iterates over all node induced anchored connected subgraphs. The
    generation of new subgraphs has linear delay.

    :param G: The networkx.Graph to get connected subgraphs from.
    :param anchor: The node index to use as anchor.

    :returns: Returns a generator to get all connected induced subgraphs. The
        return type of a single iteration is a list of node indices that form a
        connected induced subgraph.
    """
    # anchor vertext needs to be lexicographically the smallest
    nmap = {anchor: 0}
    for n in G.nodes:
        if n == anchor:
            continue
        else:
            nmap[n] = len(nmap.keys())
    nmap_inv = {v: k for k, v in nmap.items()}
    G = nx.relabel_nodes(G, nmap, copy=True)

    for subgraph in _node_induced_connected_subgraphs(G, 0):
        relabeled_subgraph = [nmap_inv[u] for u in subgraph]
        yield relabeled_subgraph
