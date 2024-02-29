import re
import copy
import itertools
import numpy as np
import networkx as nx


def print_graph(graph):
    print(
        "Graph Nodes: {}".format(
            " ".join(
                ["{}[{}]".format(n[1]["symbol"], n[0]) for n in graph.nodes(data=True)]
            )
        )
    )
    print(
        "Graph Edges: {}".format(
            " ".join(
                [
                    "[{}]-[{}]:{}".format(n[0], n[1], n[2]["bond"])
                    for n in graph.edges(data=True)
                ]
            )
        )
    )


def tokenize(pattern):
    token_specification = [
        ("ATOM", r"Br|Cl|C|N|O|P|S|F|B|I|b|c|n|o|p|s"),
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


class FGConfig:
    def __init__(self, **kwargs):
        self.parent = None

        pattern = kwargs.get("pattern", [])
        pattern = pattern if isinstance(pattern, list) else [pattern]
        self.pattern = [parse(p) for p in pattern]

        self.name = kwargs.get("name", None)
        if self.name is None:
            raise ValueError(
                "Functional group config requires a name. Add 'name' property to config."
            )

        group_atoms = kwargs.get("group_atoms", None)
        if group_atoms is None:
            group_atoms = []
            for p in self.pattern:
                group_atoms.append(list(p.nodes))
        else:
            group_atoms = (
                group_atoms if isinstance(group_atoms[0], list) else [group_atoms]
            )
        self.group_atoms = group_atoms

        anti_pattern = kwargs.get("anti_pattern", [])
        anti_pattern = (
            anti_pattern if isinstance(anti_pattern, list) else [anti_pattern]
        )
        self.anti_pattern = sorted(
            [parse(p) for p in anti_pattern],
            key=lambda x: x.number_of_nodes(),
            reverse=True,
        )

        depth = kwargs.get("depth", None)
        self.max_pattern_size = (
            depth
            if depth is not None
            else np.max([p.number_of_nodes() for p in self.pattern + self.anti_pattern])
        )

        subgroups = []
        for sgc in kwargs.get("subgroups", []):
            fgc = FGConfig(**sgc)
            fgc.parent = self
            subgroups.append(fgc)
        self.subgroups = subgroups


def get_mapping_permutations(match_symbols, sym_dict, wildcard=None, ignore_case=True):
    def _compare(match_sym, sym, wildcard, ignore_case) -> bool:
        if ignore_case:
            match_sym = match_sym.lower()
            wildcard = wildcard.lower() if wildcard is not None else None
            sym = sym.lower() if sym is not None else None
        return match_sym == wildcard or match_sym == sym

    mappings = []
    mapping_sets = []
    added_indices = []
    sym_dict = copy.deepcopy(sym_dict)
    while len(sym_dict) < len(match_symbols):
        added_indices.append(len(sym_dict))
        sym_dict.append(None)
    sym_map = [(i, s) for i, s in enumerate(sym_dict)]
    for sym_permut in itertools.permutations(sym_map):
        mapping = []
        is_match = True
        for i, s1 in enumerate(match_symbols):
            if _compare(s1, sym_permut[i][1], wildcard, ignore_case):
                sym_idx = sym_permut[i][0]
                mapping.append((i, -1 if sym_idx in added_indices else sym_idx))
            else:
                is_match = False
                break
        if is_match:
            mapping_set = set(mapping)
            if not mapping_set in mapping_sets:
                mappings.append(mapping)
                mapping_sets.append(set(mapping))

    return mappings


def anchored_pattern_match(graph, anchor, pattern, pattern_anchor, verbose=False):
    def _fits(
        idx, pidx, visited_nodes=set(), visited_pnodes=set(), indent=0, verbose=False
    ):
        assert idx not in visited_nodes
        assert pidx not in visited_pnodes

        visited_nodes = copy.deepcopy(visited_nodes)
        visited_nodes.add(idx)
        visited_pnodes = copy.deepcopy(visited_pnodes)
        visited_pnodes.add(pidx)

        node_neighbors = [
            (nidx, graph.nodes[nidx]["symbol"])
            for nidx in graph.neighbors(idx)
            if nidx not in visited_nodes
        ]
        pnode_neighbors = [
            (nidx, pattern.nodes[nidx]["symbol"])
            for nidx in pattern.neighbors(pidx)
            if nidx not in visited_pnodes
        ]

        fits = False
        match = []
        sym = graph.nodes[idx]["symbol"]
        psym = pattern.nodes[pidx]["symbol"]
        if verbose:
            print(
                "{}_fit({}@{}, {}@{}, visited_nodes=[{}], visited_pnodes=[{}]".format(
                    indent * " ", sym, idx, psym, pidx, visited_nodes, visited_pnodes
                )
            )
        if (sym == "R" or psym == "R" or sym.lower() == psym.lower()) and (
            len(node_neighbors) >= len([n for n in pnode_neighbors if n[1] != "R"])
        ):
            if verbose:
                print("{}Atom Match: {} -> {}".format(indent * " ", sym, psym))
            match.append((idx, pidx))
            if len(pnode_neighbors) > 0:
                an_syms = [a[1] for a in node_neighbors]
                pn_syms = [a[1] for a in pnode_neighbors]
                mappings = get_mapping_permutations(pn_syms, an_syms, wildcard="R")
                if verbose:
                    print(
                        "{}Check Mappings: (Pattern) {} -> (Mol) {} ==> {}".format(
                            indent * " ", pn_syms, an_syms, mappings
                        )
                    )
                for mapping in mappings:
                    valid_mapping = True
                    n_matches = set()
                    n_vnodes = set()
                    n_vpnodes = set()
                    for pn_i, an_i in mapping:
                        pnn = pnode_neighbors[pn_i]
                        if an_i == -1:
                            n_vpnodes.add(pnn[0])
                            continue
                        pnn_bond = pattern.edges[pidx, pnn[0]]["bond"]
                        nn = node_neighbors[an_i]
                        nn_bond = graph.edges[idx, nn[0]]["bond"]
                        if nn_bond != pnn_bond:
                            valid_mapping = False
                            if verbose:
                                print(
                                    "{}Invalid mapping: {} | Mol: {}-{} P: {}-{} Bond {} != {}".format(
                                        indent * " ",
                                        mapping,
                                        idx,
                                        nn[0],
                                        pidx,
                                        pnn[0],
                                        nn_bond,
                                        pnn_bond,
                                    )
                                )
                            break
                        n_fit, n_match, vnodes = _fits(
                            nn[0],
                            pnn[0],
                            visited_nodes,
                            visited_pnodes,
                            indent=indent + 2,
                            verbose=verbose,
                        )
                        if verbose:
                            print(
                                "{}Fits: {}  Match: {}".format(
                                    indent * " ", n_fit, n_match
                                )
                            )
                        if not n_fit:
                            valid_mapping = False
                            break
                        else:
                            n_vnodes.update(vnodes[0])
                            n_vpnodes.update(vnodes[1])
                            n_matches.update(n_match)
                    if valid_mapping:
                        fits = True
                        visited_nodes.update(n_vnodes)
                        visited_pnodes.update(n_vpnodes)
                        match.extend(n_matches)
                        if verbose:
                            print(
                                    "{}Valid mapping: {} ==> new matches: {}".format(
                                    indent * " ", mapping, n_matches
                                )
                            )
                        break
                    else:
                        if verbose:
                            print(
                                "{}{} is not a valid mapping.".format(
                                    indent * " ", mapping
                                )
                            )
            else:
                fits = True
        else:
            if verbose:
                print(
                    "{}No Match: {} -> {} neighbors: {} {}".format(
                        indent * " ", sym, psym, node_neighbors, pnode_neighbors
                    )
                )
        if verbose:
            print(
                "{}Visited Nodes: graph={} pattern={}".format(
                    indent * " ", visited_nodes, visited_pnodes
                )
            )
        return fits, match, (visited_nodes, visited_pnodes)

    if verbose:
        print(
            "anchored_pattern_match(gnodes: {}, {}, pnodes: {}, {})".format(
                " ".join(
                    [
                        "{}[{}]".format(n[1]["symbol"], n[0])
                        for n in graph.nodes(data=True)
                    ]
                ),
                anchor,
                " ".join(
                    [
                        "{}[{}]".format(n[1]["symbol"], n[0])
                        for n in pattern.nodes(data=True)
                    ]
                ),
                pattern_anchor,
            )
        )
    fits, match, vnodes = _fits(anchor, pattern_anchor, indent=2, verbose=verbose)
    all_matched = set(pattern.nodes).issubset(vnodes[1])
    if not all_matched:
        fits = False
        if verbose:
            print("NOT ALL PATTERN NODES MATCHED!")

    if verbose:
        print(
            "RETURN: {}, {} (visited nodes: {}, visited pattern nodes: {})".format(
                fits, match, vnodes[0], vnodes[1]
            )
        )
    if not fits:
        match = []
    return fits, match


def pattern_match(graph, anchor, pattern, pattern_anchor=None, verbose=False):
    if pattern_anchor is None:
        for pidx in pattern.nodes:
            result = anchored_pattern_match(
                graph, anchor, pattern, pidx, verbose=verbose
            )
            if result[0]:
                return result
        return False, [[]]
    else:
        return anchored_pattern_match(
            graph, anchor, pattern, pattern_anchor, verbose=verbose
        )


functional_group_config = [
    {
        "name": "carbonyl",
        "pattern": "C=O",
        "group_atoms": [0, 1],
        "subgroups": [
            {
                "name": "aldehyde",
                "pattern": "RC=O",
                "anti_pattern": ["NC=O"],
                "group_atoms": [1, 2],
            },
            {
                "name": "ketone",
                "pattern": "RC(=O)R",
                "group_atoms": [1, 2],
                "anti_pattern": ["RC(=O)O"],
            },
            {
                "name": "carboxylic_acid",
                "pattern": "RC(=O)O",
                "group_atoms": [1, 2, 3],
                "anti_pattern": ["RC(=O)OR"],
            },
            {
                "name": "carboxylate_ester",
                "pattern": "RC(=O)OR",
                "group_atoms": [1, 2, 3],
            },
            {"name": "amide", "pattern": "C(=O)N"},
        ],
    },
]

fg_configs = None


def get_FGs() -> list[FGConfig]:
    global fg_configs
    if fg_configs is None:
        c = []
        for fgc in functional_group_config:
            c.append(FGConfig(**fgc))
        fg_configs = c
    return fg_configs


def get_FGs_flat() -> list[FGConfig]:
    def _get(conf: FGConfig) -> list[FGConfig]:
        configs = [conf]
        for sg in conf.subgroups:
            configs.extend(_get(sg))
        return configs

    configs = []
    for fg in get_FGs():
        configs.extend(_get(fg))
    return configs


def get_FG_by_name(name: str) -> FGConfig:
    for fg in get_FGs_flat():
        if fg.name == name:
            return fg
    raise KeyError("No functional group config with name '{}' found.".format(name))


def get_FG_names() -> list[str]:
    return [c.name for c in get_FGs_flat()]


def get_FG_root_chain(name: str) -> list[FGConfig]:
    fg = get_FG_by_name(name)
    chain = [fg]
    while fg.parent is not None:
        chain.insert(0, fg.parent)
        fg = fg.parent
    return chain


def check_functional_group(graph, config: FGConfig, index: int, verbose=False) -> bool:
    is_func_group = False

    for pattern, group_indices in zip(config.pattern, config.group_atoms):
        is_match, match = pattern_match(graph, index, pattern, verbose=verbose)
        if is_match:
            m_idx = [m[1] for m in match if m[0] == index][0]
            is_func_group = is_func_group or m_idx in group_indices

    last_len = config.max_pattern_size
    for apattern, apattern_size in sorted(
        [(m, m.number_of_nodes()) for m in config.anti_pattern],
        key=lambda x: x[1],
        reverse=True,
    ):
        if not is_func_group:
            break
        if last_len > apattern_size:
            last_len = apattern_size
        is_match, _ = pattern_match(graph, index, apattern, verbose=verbose)
        is_func_group = is_func_group and not is_match
    return is_func_group


def mol_to_graph(mol) -> nx.Graph:
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
        g.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond=bond_order_map[bond_type],
        )
    return g
