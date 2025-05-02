import re
import networkx as nx

from fgutils.chem.valence import _check_its_valence
from fgutils.its import ITS, split_its
from fgutils.const import SYMBOL_KEY, BOND_KEY

_BOND_MAP = {"-": 1, "=": 2, ":": 1.5}


class DPORule:
    """Double Pushout Rule class.

    :param rule_id: The id of the rule. This can also be seen as the rule name.

    :param left: The left graph L in the DPO rule ``L \u2190 C \u2192 R``.
    :param context: The context graph C in the DPO rule ``L \u2190 C \u2192 R``.
    :param right: The right graph R in the DPO rule ``L \u2190 C \u2192 R``.
    """

    def __init__(
        self, rule_id: str, left: nx.Graph, context: nx.Graph, right: nx.Graph
    ):
        self.rule_id = rule_id
        self.left = left
        self.context = context
        self.right = right

    def to_rc_graph(self) -> nx.Graph:
        """Method to convert the DPO rule representation to a reaction center
        graph, i.e., the superposition of L and R on C. The reaction center
        graph is a node and edge labeled graph. The ``BOND_KEY`` edge label
        encodes the bond change in the reaction.

        :returns: Returns the reaction center graph.
        """
        rc = nx.Graph()
        for n, d in self.context.nodes(data=True):
            rc.add_node(n, **d)
        assert len(self.context.edges) == 0
        for n, d in self.left.nodes(data=True):
            if not rc.has_node(n):
                raise ValueError(
                    "Node {}:{} is not in context.".format(n, d[SYMBOL_KEY])
                )
        for u, v, d in self.left.edges(data=True):
            rc.add_edge(u, v, **{BOND_KEY: [d[BOND_KEY], 0]})
        for u, v, d in self.right.edges(data=True):
            if rc.has_edge(u, v):
                rc.edges[u, v][BOND_KEY][1] = d[BOND_KEY]
            else:
                rc.add_edge(u, v, **{BOND_KEY: [0, d[BOND_KEY]]})
        return rc


def _is_start(line):
    p = re.compile(r"^rule \[$")
    return p.match(line) is not None


def _match_id(line):
    p = re.compile(r"^\s*ruleID \"(?P<ruleID>[a-zA-Z\d]+)\"$")
    m = p.match(line)
    if m:
        return m.group("ruleID")
    return None


def _match_graph(line):
    p = re.compile(r"^\s*(?P<name>[a-z]+) \[$")
    m = p.match(line)
    if m:
        return m.group("name")
    return None


def _match_edge(line):
    p = re.compile(
        r"^\s+edge \[ source (?P<source>\d+) target (?P<target>\d+) label \"(?P<label>[-=:])\" ]$"
    )
    m = p.match(line)
    if m:
        return int(m.group("source")), int(m.group("target")), m.group("label")
    return None


def _match_node(line):
    p = re.compile(r"^\s+node \[ id (?P<id>\d+) label \"(?P<label>[a-zA-Z\d+-]+)\" ]$")
    m = p.match(line)
    if m:
        return int(m.group("id")), m.group("label")
    return None


def _is_end(line):
    p = re.compile(r"^\s*]$")
    return p.match(line) is not None


def _parse_graph(lines: list[str]):
    graph_name = _match_graph(lines.pop(0))
    edges = []
    nodes = []
    while not _is_end(lines[0]):
        line = lines.pop(0)
        edge = _match_edge(line)
        if edge is None:
            node = _match_node(line)
            if node is None:
                raise ValueError(
                    "Expected node or edge in graph not '{}'.".format(line)
                )
            nodes.append(node)
        else:
            edges.append(edge)
    lines.pop(0)
    g = nx.Graph()
    for n, sym in nodes:
        g.add_node(n, **{SYMBOL_KEY: sym})
    for u, v, bond in edges:
        g.add_edge(u, v, **{BOND_KEY: _BOND_MAP[bond]})
    return graph_name, g


def parse_gml_dpo_rule(lines: list[str]) -> DPORule:
    """Parse a GML string into a DPO rule object.

    :param lines: The lines of the GML string.

    :returns: Returns the parsed DPO rule.
    """
    start_line = lines.pop(0)
    if not _is_start(start_line):
        raise ValueError("Expected GML start not '{}'.".format(start_line))
    id_line = lines.pop(0)
    rule_id = _match_id(id_line)
    if rule_id is None:
        raise ValueError("Expected ruleID not '{}'.".format(id_line))
    g_name, left = _parse_graph(lines)
    if g_name != "left":
        raise ValueError("Expected graph 'left' but got '{}'.".format(g_name))
    g_name, context = _parse_graph(lines)
    if g_name != "context":
        raise ValueError("Expected graph 'context' but got '{}'.".format(g_name))
    g_name, right = _parse_graph(lines)
    if g_name != "right":
        raise ValueError("Expected graph 'right' but got '{}'.".format(g_name))
    end_line = lines.pop(0)
    if not _is_end(end_line):
        raise ValueError("Expected end line but got '{}'.".format(end_line))
    if len(lines) > 0:
        raise ValueError(
            "Expected no more lines but found {} more lines.".format(len(lines))
        )
    rule = DPORule(rule_id, left, context, right)
    return rule


class ReactionRule:
    """Class describing a reaction rule.

    :param rc_graph: A reaction center graph.
    :param name: (optional) The name of the reaction rule.
    """

    def __init__(self, rc_graph: nx.Graph, name="unnamed rule"):
        self.rc = rc_graph
        self.l, self.r = split_its(rc_graph)
        self.name = name

    @staticmethod
    def from_gml(src: str):
        """Load reaction rule from a GML file. Supported format is the M\u00D8D GML
        rule format.

        :param data: This can be either a file path or a GML string.

        :returns: Returns an instance of the ReactionRule class.
        """
        if src.endswith(".gml"):
            with open(src, "r") as f:
                lines = f.readlines()
        else:
            lines = src.split("\n")
        dpo_rule = parse_gml_dpo_rule(lines)
        return ReactionRule(dpo_rule.to_rc_graph(), name=dpo_rule.rule_id)


def _has_bonding_violation(its, rule, ctx2its_mapping):
    # check non-bonding condition of new edges in the reactant
    new_edges = [(u, v) for u, v, d in rule.rc.edges(data=True) if d[BOND_KEY][0] == 0]
    for ur, vr in new_edges:
        u = ctx2its_mapping[ur]
        v = ctx2its_mapping[vr]
        if its.has_edge(u, v):
            return True
    return False


def apply_rule(
    g: nx.Graph,
    rule: ReactionRule,
    n: int | None = None,
    unique=True,
    connected_only=False,
) -> list[ITS]:
    """Apply a reaction rule to a reactant graph G. The function returns the
    generated ITS graphs.

    :param g: The reactant graph. This can be a disconnected graph for multiple
        reactant molecules.

    :param rule: The reaction rule to apply to the reactant graph.

    :param n: (optional) Limits the maximum number of solutions. If n is set to
        an integer it will return once n solutions are found. (Default: None)

    :param unique: (optional) Flag to specify if isomorphic solutions should be
        returned as one solution. Isomorphism is checked with 3 iterations WL.

    :param connected_only: (optional) Flag to specify if the ITS graph must be
        connected. Setting this to true means that all reactants in g must take
        part in the reaction.

    :returns: Returns a list of ITS graphs.
    """
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        g,
        rule.l,
        node_match=lambda d1, d2: d1[SYMBOL_KEY] == d2[SYMBOL_KEY],
        edge_match=lambda d1, d2: d1[BOND_KEY] == d2[BOND_KEY],
    )
    its_graphs = {}
    for its2ctx_mapping in matcher.subgraph_monomorphisms_iter():
        ctx2its_mapping = {v: k for k, v in its2ctx_mapping.items()}
        its = g.copy()
        if _has_bonding_violation(its, rule, ctx2its_mapping):
            continue

        its_edge_attrs = {
            (u, v): [d[BOND_KEY], d[BOND_KEY]] for u, v, d in its.edges(data=True)
        }
        for ur, vr, dr in rule.rc.edges(data=True):
            u = ctx2its_mapping[ur]
            v = ctx2its_mapping[vr]
            if not its.has_edge(u, v):
                its.add_edge(u, v)
            its_edge_attrs[u, v] = dr[BOND_KEY]

        nx.set_edge_attributes(its, its_edge_attrs, BOND_KEY)
        if connected_only and not nx.is_connected(its):
            continue
        if not _check_its_valence(its):
            continue
        if unique is True:
            wl_hash = nx.weisfeiler_lehman_graph_hash(
                its, edge_attr=BOND_KEY, node_attr=SYMBOL_KEY, iterations=3
            )
            if wl_hash not in its_graphs:
                its_graphs[wl_hash] = ITS(its)
        else:
            its_graphs[len(its_graphs)] = ITS(its)
        if n is not None and len(its_graphs) == n:
            break

    return list(its_graphs.values())
