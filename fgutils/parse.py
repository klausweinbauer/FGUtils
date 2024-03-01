import re
import numpy as np
import networkx as nx

def tokenize(pattern):
    token_specification = [
        ("ATOM", r"H|Br|Cl|C|N|O|P|S|F|B|I|b|c|n|o|p|s"),
        ("BOND", r"\.|-|=|#|$|:|/|\\"),
        ("BRANCH_START", r"\("),
        ("BRANCH_END", r"\)"),
        ("RING_NUM", r"\d+"),
        ("WILDCARD", r"R"),
        ("MISMATCH", r"."),
    ]
    token_re = "|".join("(?P<%s>%s)" % pair for pair in token_specification)
    for m in re.finditer(token_re, pattern):
        ttype = m.lastgroup
        value = m.group()
        if value == "":
            break
        column = m.start()
        yield ttype, value, column


def parse(pattern, verbose=False):
    bond_to_order = {"-": 1, "=": 2, "#": 3, "$": 4, ":": 1.5, ".": 0}
    g = nx.Graph()
    anchor = None
    branches = []
    rings = {}
    bond_order = 1
    for ttype, value, col in tokenize(pattern):
        if verbose:
            print(
                "Process Token: {:>15}={} | Anchor: {}@{} Bond: {}".format(
                    ttype,
                    value,
                    g.nodes[anchor]["symbol"] if anchor is not None else "None",
                    anchor,
                    bond_order,
                )
            )
        idx = g.number_of_nodes()
        if ttype == "ATOM" or ttype == "WILDCARD":
            g.add_node(idx, symbol=value)
            if anchor is not None:
                anchor_sym = g.nodes[anchor]["symbol"]
                if bond_order == 1 and anchor_sym.islower() and value.islower():
                    bond_order = 1.5
                g.add_edge(anchor, idx, bond=bond_order)
                bond_order = 1
            anchor = idx
        elif ttype == "BOND":
            bond_order = bond_to_order[value]
        elif ttype == "BRANCH_START":
            branches.append(anchor)
        elif ttype == "BRANCH_END":
            anchor = branches.pop()
        elif ttype == "RING_NUM":
            if value in rings.keys():
                anchor_sym = g.nodes[anchor]["symbol"]
                ring_anchor = rings[value]
                ring_anchor_sym = g.nodes[ring_anchor]["symbol"]
                if anchor_sym.islower() != ring_anchor_sym.islower():
                    raise SyntaxError(
                        (
                            "Ring {} must be of same aromaticity type. "
                            + "Started with {} and ended with {}."
                        ).format(value, ring_anchor_sym, anchor_sym)
                    )
                if anchor_sym.islower():
                    bond_order = 1.5
                g.add_edge(anchor, ring_anchor, bond=bond_order)
                del rings[value]
            else:
                rings[value] = anchor
        else:
            selection = pattern[
                col - np.min([col, 4]) : col + np.min([len(pattern) - col + 1, 5])
            ]
            raise SyntaxError(
                "Invalid character found in column {} near '{}'.".format(col, selection)
            )

    return g
