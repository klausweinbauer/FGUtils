import re
import numpy as np
import networkx as nx

from fgutils.const import SYMBOL_KEY


token_specification = [
    ("ATOM", r"H|Br|Cl|Se|Sn|Si|Mg|Li|C|N|O|P|S|F|B|I|b|c|n|o|p|s"),
    ("BOND", r"\.|-|=|#|$|:|/|\\"),
    ("BRANCH_START", r"\("),
    ("BRANCH_END", r"\)"),
    ("RING_NUM", r"\d+"),
    ("WILDCARD", r"R"),
    ("RC_BOND", r"<\d*,\d*>"),
    ("NODE_LABEL", r"\{[a-zA-Z0-9_,-]+\}"),
    ("MISMATCH", r"."),
]


def tokenize(pattern):
    token_re = "|".join("(?P<%s>%s)" % pair for pair in token_specification)
    for m in re.finditer(token_re, pattern):
        ttype = m.lastgroup
        value = m.group()
        if value == "":
            break
        column = m.start()
        yield ttype, value, column


class Parser:
    """

    Class to convert a SMILES like graph description into a NetworkX graph.

    Example for parsing acetic acid::

      >>> parser = Parser()
      >>> g = parser("CC(O)=O")
      Graph with 4 nodes and 3 edges

    :param use_multigraph:

        Flag to specify if the resulting graph object should be of type
        networkx.MultiGraph or networkx.Graph. The difference is that a
        MultiGraph can have more than one edge between two nodes. For parsing
        molecule like graphs this is not necessary because bond types are
        encoded as edge labels. (Default = False)

    :param verbose:

        Flag to print information during parsing. (Default = False)

    """

    def __init__(self, use_multigraph=False, verbose=False):
        self.bond_to_order_map = {"-": 1, "=": 2, "#": 3, "$": 4, ":": 1.5, ".": None}
        self.verbose = verbose
        self.use_multigraph = use_multigraph
        self.__clear()

    def __clear(self):
        if self.use_multigraph:
            self.graph = nx.MultiGraph()
        else:
            self.graph = nx.Graph()
        self.anchor = None
        self.branches = []
        self.rings = {}
        self.bond_order = 1

    def __print_process_token(self, ttype, value):
        if self.verbose:
            print(
                "Process Token: {:>15}={} | Anchor: {}@{} Bond: {}".format(
                    ttype,
                    value,
                    self.graph.nodes[self.anchor][SYMBOL_KEY]
                    if self.anchor is not None
                    else "None",
                    self.anchor,
                    self.bond_order,
                )
            )

    def __process_token_add_node(self, ttype, value, idx):
        is_labeled = False
        labels = []
        if ttype == "NODE_LABEL":
            is_labeled = True
            labels = value.lstrip("{").rstrip("}").split(",")
            value = "#"
        self.graph.add_node(idx, symbol=value, labels=labels, is_labeled=is_labeled)
        if self.anchor is not None:
            anchor_sym = self.graph.nodes[self.anchor][SYMBOL_KEY]
            if self.bond_order == 1 and anchor_sym.islower() and value.islower():
                self.bond_order = 1.5
            if self.bond_order is not None:
                self.graph.add_edge(self.anchor, idx, bond=self.bond_order)
            self.bond_order = 1
        self.anchor = idx

    def __process_token_rc_bond(self, value):
        rc_bonds = value.replace("<", "").replace(">", "").split(",")
        if len(rc_bonds) != 2:
            raise SyntaxError(
                "Reaction center bond should be of form <g_bond, h_bond>."
            )
        g_bond, h_bond = rc_bonds[0], rc_bonds[1]
        g_bond = 1 if g_bond == "" else int(g_bond)
        h_bond = 1 if h_bond == "" else int(h_bond)
        self.bond_order = (g_bond, h_bond)

    def __process_token_ring(self, value):
        if value in self.rings.keys():
            anchor_sym = self.graph.nodes[self.anchor][SYMBOL_KEY]
            ring_anchor = self.rings[value]
            if anchor_sym.islower():
                self.bond_order = 1.5
            if self.bond_order is not None:
                self.graph.add_edge(self.anchor, ring_anchor, bond=self.bond_order)
            self.bond_order = 1
            del self.rings[value]
        else:
            if self.anchor is None:
                raise SyntaxError("Invalid ring anchor. Did you start with '1...'?")
            self.rings[value] = self.anchor

    def __process_token(self, ttype, value, idx) -> bool:
        if ttype == "ATOM" or ttype == "WILDCARD" or ttype == "NODE_LABEL":
            self.__process_token_add_node(ttype, value, idx)
        elif ttype == "BOND":
            self.bond_order = self.bond_to_order_map[value]
        elif ttype == "RC_BOND":
            self.__process_token_rc_bond(value)
        elif ttype == "BRANCH_START":
            self.branches.append(self.anchor)
        elif ttype == "BRANCH_END":
            self.anchor = self.branches.pop()
        elif ttype == "RING_NUM":
            self.__process_token_ring(value)
        else:
            return False
        return True

    def parse(self, pattern: str, idx_offset: int = 0):
        """

        Method to parse a SMILES like graph pattern.

        :param pattern:

            The pattern to convert into a graph. The pattern is a tree-like
            description of the graph. It is strongly oriented at the SMILES
            notation.

        :param idx_offset:

            The index offset argument provides the starting value for the
            consecutive node numbering. (Default = 0)

        :returns:

            Returns the converted graph object.

        """
        self.__clear()
        for ttype, value, col in tokenize(pattern):
            self.__print_process_token(ttype, value)
            idx = self.graph.number_of_nodes() + idx_offset
            if not self.__process_token(ttype, value, idx):
                selection = pattern[
                    col - np.min([col, 4]) : col + np.min([len(pattern) - col + 1, 5])
                ]
                raise SyntaxError(
                    "Invalid character '{}' found in column {} near '{}'.".format(
                        pattern[col], col, selection
                    )
                )
        return self.graph

    def __call__(self, pattern, idx_offset=0):
        return self.parse(pattern, idx_offset)


def parse(pattern, verbose=False, idx_offset=0):
    parser = Parser(verbose=verbose)
    return parser(pattern, idx_offset)
