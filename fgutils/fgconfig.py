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
    {"name": "amine", "pattern": "RN(R)R", "group_atoms": [1, 2, 3]},
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
        self.parents: list[FGConfig] = []

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

        self.subgroups: list[FGConfig] = []

    @property
    def pattern_len(self) -> int:
        return len(
            [
                _
                for _, n_sym in self.pattern.nodes(data="symbol")
                if n_sym not in self.len_exclude_nodes
            ]
        )


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


def get_FG_root_chain(name: str) -> list[FGConfig]:
    # TODO fix
    fg = get_FG_by_name(name)
    chain = [fg]
    while fg.parent is not None:
        chain.insert(0, fg.parent)
        fg = fg.parent
    return chain


def sort_by_pattern_len(configs: list[FGConfig], reverse=False) -> list[FGConfig]:
    return list(
        sorted(
            configs,
            key=lambda x: (x.pattern_len, len(x.pattern)),
            reverse=reverse,
        )
    )


def map_full(graph, pattern):
    for i in range(len(graph)):
        r, _ = map_pattern(graph, i, pattern)
        if r is True:
            return True
    return False


def is_subgroup(parent, child):
    p2c = map_full(child.pattern, parent.pattern)
    c2p = map_full(parent.pattern, child.pattern)
    if p2c:
        assert c2p == False, "{} ({}) -> {} ({}) matches in both directions.".format(
            parent.name, parent.pattern_str, child.name, child.pattern_str
        )
        return True
    return False


def search_parents(tree: list[FGConfig], child) -> None | list[FGConfig]:
    parents = set()
    for group in tree:
        is_child = is_subgroup(group, child)
        if is_child:
            sub_parents = search_parents(group.subgroups, child)
            if sub_parents is None:
                parents.add(group)
            else:
                parents.update(sub_parents)
    return None if len(parents) == 0 else list(parents)


def insert_child(parent: FGConfig, child: FGConfig):
    assert len(child.subgroups) == 0
    new_subgroups = []
    for sg in parent.subgroups:
        assert parent in sg.parents, "Inconsistency in parent child structure."
        if is_subgroup(child, sg):
            sg.parents.remove(parent)
            sg.parents.append(child)
            child.subgroups.append(sg)
        else:
            new_subgroups.append(sg)
    child.parents.append(parent)
    new_subgroups.append(child)
    parent.subgroups = new_subgroups

    parent.subgroups = sort_by_pattern_len(parent.subgroups, reverse=True)
    child.subgroups = sort_by_pattern_len(child.subgroups, reverse=True)


def print_tree(tree: list[FGConfig]):
    def _print(fg, indent=0):
        print(
            "{}{:<{width}}{:<40} {}".format(
                indent * " ",
                fg.name,
                "[Parents: {}]".format(
                    ", ".join([g.name for g in fg.parents])
                    if len(fg.parents) > 0
                    else "ROOT"
                ),
                fg.pattern_str,
                width=30 - indent,
            )
        )
        for sfg in fg.subgroups:
            _print(sfg, indent + 2)

    for fg in tree:
        _print(fg)


def build_config_tree_from_list(config_list: list[FGConfig]) -> list[FGConfig]:
    roots = []
    for fg in sort_by_pattern_len(config_list):
        print("Add {} to tree.".format(fg.name))
        fg.subgroups = []
        parents = search_parents(roots, fg)
        print(
            "Found parents: {}".format(
                [g.name for g in parents] if parents is not None else []
            )
        )
        if parents is None:
            roots.append(fg)
        else:
            for parent in parents:
                print("Insert child {} to parent {}.".format(fg.name, parent.name))
                insert_child(parent, fg)
        print_tree(roots)
        print("--- DONE ---")
    return roots


def build_FG_tree() -> list[FGConfig]:
    return build_config_tree_from_list(get_FG_list())
