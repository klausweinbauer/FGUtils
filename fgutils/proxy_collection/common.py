from fgutils.proxy import ProxyGroup, ProxyGraph


simple_groups = []
complex_groups = []

hydrogen_group = ProxyGroup("H", "")
simple_groups += [hydrogen_group]

alkyl_group_names = [
    "methyl",
    "ethyl",
    "propyl",
    # "butyl",
    # "pentyl",
    # "hexyl",
    # "heptyl",
    # "octyl",
    # "nonyl",
    # "decyl",
    # "undecyl",
    # "dodecyl",
    # "tridecyl",
    # "tetradecyl",
]

alkyl_groups = [
    ProxyGroup(name, "C" * (i + 1)) for i, name in enumerate(alkyl_group_names)
]
simple_groups += alkyl_groups

carbon_chain_max_length = 4 #10
carbon_chains = [ProxyGroup("CC1", "C")] + [
    ProxyGroup("CC{}".format(i), ProxyGraph("C" * i, anchor=[0, i - 1]))
    for i in range(2, carbon_chain_max_length + 1)
]
simple_groups += carbon_chains

aryl_groups = [
    ProxyGroup("phenyl", "c1ccccc1"),
    ProxyGroup("o-tolyl", "c1c(C)cccc1"),
    ProxyGroup("m-tolyl", "c1cc(C)ccc1"),
    ProxyGroup("p-tolyl", "c1ccc(C)cc1"),
]
simple_groups += aryl_groups

amine_groups = [
    ProxyGroup("1-amine", ProxyGraph("N")),
    ProxyGroup("2-amine", ProxyGraph("N{alkyl}")),
    ProxyGroup("3-amine", ProxyGraph("N({alkyl}){alkyl}")),
]
simple_groups += amine_groups


def add_multi_group(complex_groups, groups, name):
    graphs = ["{{{}}}".format(g.name) for g in groups]
    complex_groups += [ProxyGroup(name, graphs)]


add_multi_group(complex_groups, alkyl_groups, "alkyl")
add_multi_group(complex_groups, carbon_chains, "carbon_chain")
add_multi_group(complex_groups, aryl_groups, "aryl")
add_multi_group(complex_groups, amine_groups, "amine")

simple_groups += [ProxyGroup("alkohol", "CO")]
simple_groups += [ProxyGroup("trimethylsilanol", "OSi(C)(C)C")]
simple_groups += [ProxyGroup("aldehyde", "C=O")]
simple_groups += [ProxyGroup("acid", ProxyGraph("C(=O)O"))]
simple_groups += [
    ProxyGroup(
        "ester",
        [
            ProxyGraph("C(=O)O{aryl}", anchor=[0, 3]),
            ProxyGraph("C(=O)O{alkyl}", anchor=[0, 3]),
        ],
    )
]
simple_groups += [
    ProxyGroup(
        "ether",
        [
            ProxyGraph("O{aryl}", anchor=[0, 1]),
            ProxyGraph("O{alkyl}", anchor=[0, 1]),
        ],
    )
]
simple_groups += [ProxyGroup("nitrogen_dioxide", "N(O)O")]
simple_groups += [ProxyGroup("halogen", ["F", "Cl", "Br", "I"])]
simple_groups += [ProxyGroup("nitrile", "C#N")]
simple_groups += [ProxyGroup("alkenyl", ["C=CC", "C=C(C)C"])]
simple_groups += [ProxyGroup("hydrogen_sulfite", "S(=O)(=O)O")]
simple_groups += [ProxyGroup("carbonyl", ["C(=O){aryl}", "C(=O){alkyl}"])]

common_groups = simple_groups + complex_groups

any_group = ProxyGroup("any", ["{{{}}}".format(g.name) for g in common_groups])
common_groups += [any_group]
