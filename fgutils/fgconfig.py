from __future__ import annotations
import numpy as np

from fgutils.permutation import PermutationMapper
from fgutils.parse import parse
from fgutils.mapping import map_to_entire_graph

_default_fg_config = [
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
    {"name": "alcohol", "pattern": "COH", "group_atoms": [1, 2]},
    {
        "name": "primary_alcohol",
        "pattern": "CCOH",
        "group_atoms": [2, 3],
        "anti_pattern": ["CC(O)O"],
    },
    {
        "name": "secondary_alcohol",
        "pattern": "C(C)(C)OH",
        "group_atoms": [3, 4],
        "anti_pattern": ["CC(O)O"],
    },
    {
        "name": "tertiary_alcohol",
        "pattern": "C(C)(C)(C)OH",
        "group_atoms": [4, 5],
        "anti_pattern": ["CC(O)O"],
    },
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
    {"name": "peroxide", "pattern": "ROOR", "group_atoms": [1, 2]},
    {"name": "peroxy_acid", "pattern": "RC(=O)OOH", "group_atoms": [1, 2, 3, 4, 5]},
    {"name": "hemiketal", "pattern": "RC(OH)(OR)R", "group_atoms": [1, 2, 3, 4]},
    {"name": "phenol", "pattern": "C:COH", "group_atoms": [2, 3]},
    {"name": "anilin", "pattern": "C:CN(R)R", "group_atoms": [2]},
    {"name": "ketene", "pattern": "RC(R)=C=O", "group_atoms": [1, 3, 4]},
    {"name": "carbamate", "pattern": "ROC(=O)N(R)R", "group_atoms": [1, 2, 3, 4]},
    {"name": "acyl_chloride", "pattern": "RC(=O)Cl", "group_atoms": [1, 2, 3]},
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


def is_subgroup(parent: FGConfig, child: FGConfig, mapper: PermutationMapper) -> bool:
    p2c = map_to_entire_graph(child.pattern, parent.pattern, mapper)
    c2p = map_to_entire_graph(parent.pattern, child.pattern, mapper)
    if p2c:
        assert c2p is False, "{} ({}) -> {} ({}) matches in both directions.".format(
            parent.name, parent.pattern_str, child.name, child.pattern_str
        )
        for anti_pattern in parent.anti_pattern:
            p2c_anti = map_to_entire_graph(child.pattern, anti_pattern, mapper)
            if p2c_anti:
                return False
        return True
    return False


class FGTreeNode:
    def __init__(self, fgconfig: FGConfig):
        self.fgconfig = fgconfig
        self.parents: list[FGTreeNode] = []
        self.children: list[FGTreeNode] = []

    def order_id(self):
        return (
            self.fgconfig.pattern_len,
            len(self.fgconfig.pattern),
            hash(self.fgconfig.pattern_str),
        )

    def add_child(self, child: FGTreeNode):
        child.parents.append(self)
        self.children.append(child)
        self.parents = sorted(self.parents, key=lambda x: x.order_id(), reverse=True)
        self.children = sorted(self.children, key=lambda x: x.order_id(), reverse=True)


def sort_by_pattern_len(configs: list[FGConfig], reverse=False) -> list[FGConfig]:
    return list(
        sorted(
            configs,
            key=lambda x: (x.pattern_len, len(x.pattern), hash(x.pattern_str)),
            reverse=reverse,
        )
    )


def search_parents(
    roots: list[FGTreeNode], child: FGTreeNode, mapper: PermutationMapper
) -> None | list[FGTreeNode]:
    parents = set()
    for root in roots:
        if is_subgroup(root.fgconfig, child.fgconfig, mapper):
            _parents = search_parents(root.children, child, mapper)
            if _parents is None:
                parents.add(root)
            else:
                parents.update(_parents)
    return None if len(parents) == 0 else list(parents)


def tree2str(roots: list[FGTreeNode]):
    sym = {"branch": "├── ", "skip": "│   ", "end": "└── ", "empty": "    "}

    def _print(node: FGTreeNode, indent, is_last=False, width=(35, 30)):
        result = "{}{}{:<{width1}}{:<{width2}}{}\n".format(
            indent,
            sym["end"] if is_last else sym["branch"],
            node.fgconfig.name,
            "[{}]".format(
                ", ".join([p.fgconfig.name for p in node.parents])
                if len(node.parents) > 0
                else "ROOT"
            ),
            node.fgconfig.pattern_str,
            width1=width[0] - len(indent) - len(sym["branch"]),
            width2=width[1],
        )

        for i, child in enumerate(node.children):
            _is_last = i == len(node.children) - 1
            if is_last:
                _indent = indent + sym["empty"]
            else:
                _indent = indent + "{}".format(sym["skip"])
            result += _print(child, _indent, _is_last, width=width)
        return result

    in_const = len(sym["skip"])
    width = (40, 25)
    tree_str = "{}{:<{width1}}{:<{width2}}{:<}\n".format(
        " " * in_const,
        "Functional Group",
        "Parents",
        "Pattern",
        width1=width[0] - in_const,
        width2=width[1],
    )
    for root in roots:
        tree_str += _print(root, "", width=width)

    result = []
    max_l = 0
    for line in tree_str.split("\n"):
        _line = line[4:]
        if len(_line) > max_l:
            max_l = len(_line)
        result.append(_line)

    result.insert(1, "{}".format("-" * max_l))

    return "\n".join(result)


def print_tree(roots: list[FGTreeNode]):
    print(tree2str(roots))


def build_config_tree_from_list(
    config_list: list[FGConfig], mapper: PermutationMapper
) -> list[FGTreeNode]:
    roots = []
    for config in sort_by_pattern_len(config_list):
        node = FGTreeNode(config)
        parents = search_parents(roots, node, mapper)
        if parents is None:
            roots.append(node)
        else:
            for parent in parents:
                parent.add_child(node)
    return roots


class FGConfigProvider:
    def __init__(
        self,
        config: list[dict] | list[FGConfig] | None = None,
        mapper: PermutationMapper | None = None,
    ):
        self.config_list: list[FGConfig] = []
        if config is None:
            config = _default_fg_config
        if isinstance(config, list) and len(config) > 0:
            if isinstance(config[0], dict):
                for fgc in config:
                    self.config_list.append(FGConfig(**fgc))  # type: ignore
            elif isinstance(config[0], FGConfig):
                self.config_list = config  # type: ignore
            else:
                raise ValueError("Invalid config value.")
        else:
            raise ValueError("Invalid config value.")

        self.mapper = (
            mapper
            if mapper is not None
            else PermutationMapper(wildcard="R", ignore_case=True)
        )

        self.__tree_roots = None

    def get_tree(self) -> list[FGTreeNode]:
        if self.__tree_roots is None:
            self.__tree_roots = build_config_tree_from_list(
                self.config_list, self.mapper
            )
        return self.__tree_roots

    def get_by_name(self, name: str) -> FGConfig:
        for fg in self.config_list:
            if fg.name == name:
                return fg
        raise KeyError("No functional group config with name '{}' found.".format(name))
