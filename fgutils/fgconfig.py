from __future__ import annotations
import numpy as np

from fgutils.parse import parse
from fgutils.mapping import map_pattern

functional_group_config = [
    {
        "name": "carbonyl",
        "pattern": "C(=O)",
    },
    {
        "name": "aldehyde",
        "pattern": "RC(=O)H",
        "group_atoms": [1, 2],
    },
    {
        "name": "ketone",
        "pattern": "RC(=O)R",
        "group_atoms": [1, 2],
    },
    {
        "name": "carboxylic_acid",
        "pattern": "RC(=O)OH",
        "group_atoms": [1, 2, 3],
    },
    {"name": "amide", "pattern": "RC(=O)N(R)R", "group_atoms": [1, 2, 3]},
    {"name": "alcohol", "pattern": "ROH", "group_atoms": [1, 2]},
    {"name": "enol", "pattern": "C=COH"},
    {"name": "acetal", "pattern": "RC(OC)(OC)H", "group_atoms": [1, 2, 4, 6]},
    {"name": "ketal", "pattern": "RC(OR)(OR)R", "group_atoms": [1, 2, 4]},
    {"name": "hemiacetal", "pattern": "RC(OC)(OH)H", "group_atoms": [1, 2, 4, 5, 6]},
    {"name": "ether", "pattern": "ROR", "group_atoms": [1]},
    {"name": "thioether", "pattern": "RSR", "group_atoms": [1]},
    {"name": "ester", "pattern": "RC(=O)OR", "group_atoms": [1, 2, 3]},
    {"name": "thioester", "pattern": "RC(=O)SR", "group_atoms": [1, 2, 3]},
    {"name": "anhydride", "pattern": "RC(=O)OC(=O)R", "group_atoms": [1, 2, 3, 4, 5]},
    {"name": "amine", "pattern": "RN(R)R", "group_atoms": [1]},
    {"name": "nitrile", "pattern": "RC#N", "group_atoms": [1, 2]},
    {"name": "nitrose", "pattern": "RN=O", "group_atoms": [1, 2]},
    {"name": "nitro", "pattern": "RN(=O)O", "group_atoms": [1, 2, 3]},
    {"name": "peroxy_acid", "pattern": "RC(=O)OOH", "group_atoms": [1, 2, 3, 4, 5]},
    {"name": "hemiketal", "pattern": "RC(OH)(OR)R", "group_atoms": [1, 2, 3, 4]},
    {"name": "phenol", "pattern": "C:COH", "group_atoms": [2, 3]},
    {"name": "anilin", "pattern": "C:CN(R)R", "group_atoms": [2]},
    {"name": "ketene", "pattern": "RC(R)=C=O", "group_atoms": [1, 3, 4]},
    {"name": "carbamate", "pattern": "ROC(=O)N(R)R", "group_atoms": [1, 2, 3, 4]},
]


class FGConfig:
    len_exclude_nodes = ["R"]

    def __init__(self, **kwargs):
        self.pattern_str = kwargs.get("pattern", None)
        if self.pattern_str is None:
            raise ValueError("Expected value for argument pattern.")
        self.pattern = parse(self.pattern_str)

        self.name = kwargs.get("name", None)
        if self.name is None:
            raise ValueError(
                "Functional group config requires a name. Add 'name' property to config."
            )

        group_atoms = kwargs.get("group_atoms", None)
        if group_atoms is None:
            group_atoms = list(self.pattern.nodes)
        if not isinstance(group_atoms, list):
            raise ValueError("Argument group_atoms must be a list.")
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
            else np.max(
                [p.number_of_nodes() for p in [self.pattern] + self.anti_pattern]
            )
        )

    @property
    def pattern_len(self) -> int:
        return len(
            [
                _
                for _, n_sym in self.pattern.nodes(data="symbol")  # type: ignore
                if n_sym not in self.len_exclude_nodes
            ]
        )


def is_subgroup(parent: FGConfig, child: FGConfig) -> bool:
    p2c = map_full(child.pattern, parent.pattern)
    c2p = map_full(parent.pattern, child.pattern)
    if p2c:
        assert c2p == False, "{} ({}) -> {} ({}) matches in both directions.".format(
            parent.name, parent.pattern_str, child.name, child.pattern_str
        )
        return True
    return False


class TreeNode:
    def __init__(self, is_child_callback):
        self.parents: list[TreeNode] = []
        self.children: list[TreeNode] = []
        self.is_child_callback = is_child_callback

    def is_child(self, parent: TreeNode) -> bool:
        return self.is_child_callback(parent, self)

    def add_child(self, child: TreeNode):
        child.parents.append(self)
        self.children.append(child)


class FGTreeNode(TreeNode):
    def __init__(self, fgconfig: FGConfig):
        self.fgconfig = fgconfig
        self.parents: list[FGTreeNode]
        self.children: list[FGTreeNode]
        super().__init__(lambda a, b: is_subgroup(a.fgconfig, b.fgconfig))

    def order_id(self):
        return (
            self.fgconfig.pattern_len,
            len(self.fgconfig.pattern),
            hash(self.fgconfig.pattern_str),
        )

    def add_child(self, child: FGTreeNode):
        super().add_child(child)
        self.parents = sorted(self.parents, key=lambda x: x.order_id(), reverse=True)
        self.children = sorted(self.children, key=lambda x: x.order_id(), reverse=True)


fg_configs = None


def get_FG_list() -> list[FGConfig]:
    global fg_configs
    if fg_configs is None:
        c = []
        for fgc in functional_group_config:
            c.append(FGConfig(**fgc))
        fg_configs = c
    return fg_configs


def get_FG_by_name(name: str) -> FGConfig:
    for fg in get_FG_list():
        if fg.name == name:
            return fg
    raise KeyError("No functional group config with name '{}' found.".format(name))


def get_FG_names() -> list[str]:
    return [c.name for c in get_FG_list()]


def sort_by_pattern_len(configs: list[FGConfig], reverse=False) -> list[FGConfig]:
    return list(
        sorted(
            configs,
            key=lambda x: (x.pattern_len, len(x.pattern), hash(x.pattern_str)),
            reverse=reverse,
        )
    )


def map_full(graph, pattern):
    for i in range(len(graph)):
        r, _ = map_pattern(graph, i, pattern)
        if r is True:
            return True
    return False


def search_parents(roots: list[TreeNode], child: TreeNode) -> None | list[TreeNode]:
    parents = set()
    for root in roots:
        if child.is_child(root):
            _parents = search_parents(root.children, child)
            if _parents is None:
                parents.add(root)
            else:
                parents.update(_parents)
    return None if len(parents) == 0 else list(parents)


def print_tree(roots: list[FGTreeNode]):
    def _print(node: FGTreeNode, indent=0):
        print(
            "{}{:<{width}}{:<40} {}".format(
                indent * " ",
                node.fgconfig.name,
                "[Parents: {}]".format(
                    ", ".join([p.fgconfig.name for p in node.parents])
                    if len(node.parents) > 0
                    else "ROOT"
                ),
                node.fgconfig.pattern_str,
                width=30 - indent,
            )
        )
        for child in node.children:
            _print(child, indent + 2)

    for root in roots:
        _print(root)


def build_config_tree_from_list(config_list: list[FGConfig]) -> list[FGTreeNode]:
    roots = []
    for config in sort_by_pattern_len(config_list):
        node = FGTreeNode(config)
        parents = search_parents(roots, node)
        if parents is None:
            roots.append(node)
        else:
            for parent in parents:
                parent.add_child(node)
    return roots


_fg_tree_roots = None


def build_FG_tree() -> list[FGTreeNode]:
    global _fg_tree_roots
    if _fg_tree_roots is None:
        _fg_tree_roots = build_config_tree_from_list(get_FG_list())
    return _fg_tree_roots
