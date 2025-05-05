import networkx as nx

from fgutils.its import ITS
from fgutils.const import SYMBOL_KEY, BOND_KEY
from fgutils.utils import to_non_aromatic_symbol
from .ps import get_num_valence_electrons


def _check_mol_valence(g: nx.Graph, exact=False) -> bool:
    for n, d in g.nodes(data=True):
        exp_valence = 8
        if d[SYMBOL_KEY] == "H":
            exp_valence = 2
        sym = to_non_aromatic_symbol(d[SYMBOL_KEY])
        valence_electrons = get_num_valence_electrons(sym)
        covalent_bonds = 0
        for neighbor in g.neighbors(n):
            bond = g.edges[n, neighbor][BOND_KEY]
            if neighbor == n:
                bond = int(2 * bond)
            covalent_bonds += bond
        valence = valence_electrons + covalent_bonds
        if exact:
            valid = valence == exp_valence
        else:
            # Need 0.5 tolerance for aromatic systems
            valid = valence <= (exp_valence + 0.5)
        if not valid:
            return False
    return True


def _check_its_valence(its: nx.Graph | ITS, exact=False) -> bool:
    if isinstance(its, ITS):
        its = its.graph
    for n, d in its.nodes(data=True):
        exp_valence = 8
        if d[SYMBOL_KEY] == "H":
            exp_valence = 2
        sym = to_non_aromatic_symbol(d[SYMBOL_KEY])
        valence_electrons = get_num_valence_electrons(sym)
        g_covalent_bonds = 0
        h_covalent_bonds = 0
        for neighbor in its.neighbors(n):
            g_bond = its.edges[n, neighbor][BOND_KEY][0]
            h_bond = its.edges[n, neighbor][BOND_KEY][1]
            if neighbor == n:
                g_bond = int(2 * g_bond)
                h_bond = int(2 * h_bond)
            g_covalent_bonds += g_bond
            h_covalent_bonds += h_bond
        g_valence = valence_electrons + g_covalent_bonds
        h_valence = valence_electrons + h_covalent_bonds
        if exact:
            valid = g_valence == exp_valence and h_valence == exp_valence
        else:
            # Need 0.5 tolerance for aromatic systems
            valid = g_valence <= (exp_valence + 0.5) and h_valence <= (
                exp_valence + 0.5
            )
        if not valid:
            return False
    return True


def check_valence(graph: nx.Graph | ITS, exact=False) -> bool:
    """Check if the valence in a molecular graph or an ITS graph is valid. Use
    exact if you also want to check if all hydrogens are correct.

    :param graph: The molecular graph or an ITS graph to check the valence for.
    :param exact: (optional) If set to true the valence must be exact for all
        atoms. (Default: False)

    :returns: Returns true if the valence is valid for all atoms and false
        otherwise.
    """
    if isinstance(graph, ITS):
        return _check_its_valence(graph, exact=exact)
    else:
        edges = list(graph.edges(data=True))
        if len(edges) == 0:
            return _check_mol_valence(graph, exact=exact)
        bond = edges[0][2][BOND_KEY]
        if isinstance(bond, tuple) or isinstance(bond, list):
            return _check_its_valence(graph, exact=exact)
        elif isinstance(bond, int) or isinstance(bond, float):
            return _check_mol_valence(graph, exact=exact)
        else:
            raise NotImplementedError(
                "Bond type '{}' not implemented in valence check.".format(type(bond))
            )
