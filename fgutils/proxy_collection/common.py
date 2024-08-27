from fgutils.proxy import ProxyGroup, ProxyGraph


simple_groups = []
complex_groups = []

alkyl_group_names = [
    "methyl",
    "ethyl",
    "propyl",
    "butyl",
    "pentyl",
    "hexyl",
    "heptyl",
    "octyl",
    "nonyl",
    "decyl",
    "undecyl",
    "dodecyl",
    "tridecyl",
    "tetradecyl",
]

alkyl_groups = [
    ProxyGroup(name, pattern="C" * (i + 1)) for i, name in enumerate(alkyl_group_names)
]
simple_groups += alkyl_groups

carbon_chain_max_length = 10
carbon_chains = [ProxyGroup("CC1", pattern="C")] + [
    ProxyGroup("CC{}".format(i), ProxyGraph("C" * i, anchor=[0, i - 1]))
    for i in range(2, carbon_chain_max_length + 1)
]
simple_groups += carbon_chains

aryl_groups = [
    ProxyGroup("phenyl", pattern="c1ccccc1"),
    ProxyGroup("o-tolyl", pattern="c1c(C)cccc1"),
    ProxyGroup("m-tolyl", pattern="c1cc(C)ccc1"),
    ProxyGroup("p-tolyl", pattern="c1ccc(C)cc1"),
]
simple_groups += aryl_groups

amine_groups = [
    ProxyGroup("1-amine", ProxyGraph("N")),
    ProxyGroup("2-amine", ProxyGraph("N{alkyl}")),
    ProxyGroup("3-amine", ProxyGraph("N({alkyl}){alkyl}")),
]
simple_groups += amine_groups


def add_multi_group(complex_groups, groups, name):
    pattern = "{" + ",".join([g.name for g in groups]) + "}"
    complex_groups += [ProxyGroup(name, pattern=pattern)]


add_multi_group(complex_groups, alkyl_groups, "alkyl")
add_multi_group(complex_groups, carbon_chains, "carbon_chain")
add_multi_group(complex_groups, aryl_groups, "aryl")
add_multi_group(complex_groups, amine_groups, "amine")

simple_groups += [ProxyGroup("alkohol", pattern="CO")]
simple_groups += [ProxyGroup("trimethylsilanol", pattern="OSi(C)(C)C")]
simple_groups += [ProxyGroup("aldehyde", pattern="C=O")]
simple_groups += [ProxyGroup("acid", ProxyGraph("C(=O)O"))]
simple_groups += [ProxyGroup("ester", ProxyGraph("C(=O)O{aryl,alkyl}", anchor=[0, 3]))]
simple_groups += [ProxyGroup("ether", ProxyGraph("O{aryl,alkyl}", anchor=[0, 1]))]
simple_groups += [ProxyGroup("nitrogen_dioxide", pattern="N(O)O")]
simple_groups += [ProxyGroup("halogen", pattern=["F", "Cl", "Br", "I"])]
simple_groups += [ProxyGroup("nitrile", pattern="C#N")]
simple_groups += [ProxyGroup("alkenyl", pattern=["C=CC", "C=C(C)C", "C=C(C)(C)C"])]
simple_groups += [ProxyGroup("hydrogen_sulfite", pattern="S(=O)(=O)O")]
simple_groups += [ProxyGroup("carbonyl", pattern="C(=O){aryl,alkyl}")]

common_groups = simple_groups + complex_groups
