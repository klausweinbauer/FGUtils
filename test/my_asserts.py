import networkx as nx

from fgutils.const import AAM_KEY


def _node_match(n1, n2, ignore_keys):
    for k, v in n1.items():
        if k in ignore_keys:
            continue
        if k not in n2.keys() or n2[k] != v:
            print("unequal 1", k)
            return False
    for k, v in n2.items():
        if k in ignore_keys:
            continue
        if k not in n1.keys() or n1[k] != v:
            print("unequal 2", k)
            return False
    return True


def _edge_match(e1, e2, ignore_keys):
    for k, v in e1.items():
        if k in ignore_keys:
            continue
        if k not in e2.keys() or e2[k] != v:
            print("unequal 3", k)
            return False
    for k, v in e2.items():
        if k in ignore_keys:
            continue
        if k not in e1.keys() or e1[k] != v:
            print("unequal 4", k)
            return False
    return True


def assert_graph_eq(exp_graph, act_graph, ignore_keys=[AAM_KEY]):
    def _nm(n1, n2):
        return _node_match(n1, n2, ignore_keys)

    def _em(e1, e2):
        return _edge_match(e1, e2, ignore_keys)

    is_isomorphic = nx.is_isomorphic(
        exp_graph, act_graph, node_match=_nm, edge_match=_em
    )
    assert is_isomorphic, "Graphs are not isomorphic."
