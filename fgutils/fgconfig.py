import numpy as np

from fgutils.parse import parse

functional_group_config = [
    {
        "name": "carbonyl",
        "pattern": "RC(=O)R",
        "group_atoms": [1, 2],
        "subgroups": [
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
            {
                "name": "carboxylate_ester",
                "pattern": "RC(=O)OR",
                "group_atoms": [1, 2, 3],
            },
            {"name": "amide", "pattern": "RC(=O)N(R)R", "group_atoms": [1, 2, 3]},
        ],
    },
]


class FGConfig:
    len_exclude_nodes = ["R"]

    def __init__(self, **kwargs):
        self.parent = None

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

        subgroups = []
        for sgc in kwargs.get("subgroups", []):
            fgc = FGConfig(**sgc)
            fgc.parent = self
            subgroups.append(fgc)
        self.subgroups = sorted(
            subgroups, key=lambda x: (x.pattern_len, len(x.pattern.nodes)), reverse=True
        )

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


def get_FG_tree() -> list[FGConfig]:
    global fg_configs
    if fg_configs is None:
        c = []
        for fgc in functional_group_config:
            c.append(FGConfig(**fgc))
        fg_configs = c
    return fg_configs


def get_FG_list() -> list[FGConfig]:
    def _get(conf: FGConfig) -> list[FGConfig]:
        configs = [conf]
        for sg in conf.subgroups:
            configs.extend(_get(sg))
        return configs

    configs = []
    for fg in get_FG_tree():
        configs.extend(_get(fg))
    return configs


def get_FG_by_name(name: str) -> FGConfig:
    for fg in get_FG_list():
        if fg.name == name:
            return fg
    raise KeyError("No functional group config with name '{}' found.".format(name))


def get_FG_names() -> list[str]:
    return [c.name for c in get_FG_list()]


def get_FG_root_chain(name: str) -> list[FGConfig]:
    fg = get_FG_by_name(name)
    chain = [fg]
    while fg.parent is not None:
        chain.insert(0, fg.parent)
        fg = fg.parent
    return chain
