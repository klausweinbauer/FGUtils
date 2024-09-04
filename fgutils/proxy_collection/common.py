from fgutils.proxy import ProxyGroup, ProxyGraph


simple_groups = []
complex_groups = []

hydrogen_group = ProxyGroup("H", "")
simple_groups += [hydrogen_group]


alkyl_carbon_name_map = {
    1: ["methyl"],
    2: ["ethyl"],
    3: ["propyl", "isopropyl"],
    4: ["butyl", "isobutyl", "sec-butyl", "tert-butyl"],
    5: [
        "pentyl",
        "sec-pentyl",
        "3-pentyl",
        "isopentyl",
        "sec-isopentyl",
        "tert-pentyl",
        "neopentyl",
        "active_pentyl",
    ],
}

alkyl_groups = [
    ProxyGroup("methyl", "C"),
    ProxyGroup("ethyl", "CC"),
    ProxyGroup("propyl", "CCC"),
    ProxyGroup("isopropyl", "C(C)C"),
    ProxyGroup("butyl", "CCCC"),
    ProxyGroup("isobutyl", "CC(C)C"),
    ProxyGroup("sec-butyl", "C(C)CC"),
    ProxyGroup("tert-butyl", "C(C)(C)C"),
    ProxyGroup("pentyl", "CCCCC"),
    ProxyGroup("sec-pentyl", "C(C)CCC"),
    ProxyGroup("3-pentyl", "C(CC)CC"),
    ProxyGroup("isopentyl", "CCC(C)C"),
    ProxyGroup("sec-isopentyl", "C(C)C(C)C"),
    ProxyGroup("tert-pentyl", "C(C)(C)CC"),
    ProxyGroup("neopentyl", "CC(C)(C)C"),
    ProxyGroup("active_pentyl", "CC(C)CC"),
]


def get_alkyl_group(name, carbon_counts=None, max_carbons=5):
    graphs = []
    if carbon_counts is None:
        carbon_counts = range(1, max_carbons + 1)
    for cc in carbon_counts:
        for alkyl_name in alkyl_carbon_name_map[cc]:
            graphs.append(ProxyGraph("{{{}}}".format(alkyl_name)))
    group = ProxyGroup(name, graphs)
    return group


simple_groups += alkyl_groups

carbon_chain_max_length = 4  # 10
carbon_chains = [ProxyGroup("CC1", "C")] + [
    ProxyGroup("CC{}".format(i), ProxyGraph("C" * i, anchor=[0, i - 1]))
    for i in range(2, carbon_chain_max_length + 1)
]
simple_groups += carbon_chains

# alkyl_2_groups = [
#     ProxyGroup("methyl_2", ProxyGraph("C")),
#     ProxyGroup("ethyl_2", ProxyGraph("CC", anchor=[0, 1])),
#     ProxyGroup("propyl_2", ProxyGraph("CCC", anchor=[0, 2])),
#     ProxyGroup("butyl_2", ProxyGraph("CCCC", anchor=[0, 3])),
#     ProxyGroup("isobutyl_2", ProxyGraph("CC(C)C", anchor=[0, 3])),
#     ProxyGroup("pentyl_2", ProxyGraph("CCCCC", anchor=[0, 4])),
#     ProxyGroup("3-pentyl_2", ProxyGraph("CC(CC)C", anchor=[0, 4])),
#     ProxyGroup("active_pentyl_2", ProxyGraph("CC(C)CC", anchor=[0, 4])),
#     ProxyGroup("isopentyl_2", ProxyGraph("CCC(C)C", anchor=[0, 4])),
#     ProxyGroup("neopentyl_2", ProxyGraph("CC(C)(C)C", anchor=[0, 4])),
# ]
# simple_groups += alkyl_2_groups
#
#
# def get_alkyl_2_group(name, carbon_counts=None, max_carbons=5):
#     group_names = [g.name for g in simple_groups]
#     graphs = []
#     if carbon_counts is None:
#         carbon_counts = range(1, max_carbons + 1)
#     for cc in carbon_counts:
#         for alkyl_name in alkyl_carbon_name_map[cc]:
#             alkyl_2_name = "{}_2".format(alkyl_name)
#             if alkyl_2_name not in group_names:
#                 continue
#             graphs.append(ProxyGraph("{{{}}}".format(alkyl_2_name)))
#     group = ProxyGroup(name, graphs)
#     return group


aryl_groups = [
    ProxyGroup("phenyl", "c1ccccc1"),
    ProxyGroup("o-tolyl", "c1c(C)cccc1"),
    ProxyGroup("m-tolyl", "c1cc(C)ccc1"),
    ProxyGroup("p-tolyl", "c1ccc(C)cc1"),
    ProxyGroup("2-3-xylene", "c1c(C)c(C)ccc1"),
    ProxyGroup("2-4-xylene", "c1c(C)cc(C)cc1"),
    ProxyGroup("2-5-xylene", "c1c(C)ccc(C)c1"),
    ProxyGroup("2-6-xylene", "c1c(C)cccc(C)1"),
    ProxyGroup("3-4-xylene", "c1cc(C)c(C)cc1"),
    ProxyGroup("3-5-xylene", "c1cc(C)cc(C)c1"),
    ProxyGroup("naphthalene", "c1cccc2c1cccc2"),
    ProxyGroup("beta-naphthalene", ProxyGraph("c1cccc2c1cccc2", anchor=[1])),
]
simple_groups += aryl_groups

simple_groups += [
    get_alkyl_group("__amine_end1", carbon_counts=[1, 3]),
    get_alkyl_group("__amine_end2", max_carbons=2),
]

amine_groups = [
    ProxyGroup("1-amine", ProxyGraph("N")),
    ProxyGroup("2-amine", ProxyGraph("N{__amine_end2}")),
    ProxyGroup("3-amine", ProxyGraph("N({__amine_end1}){__amine_end2}")),
]
simple_groups += amine_groups


def add_multi_group(complex_groups, groups, name):
    graphs = ["{{{}}}".format(g.name) for g in groups]
    complex_groups += [ProxyGroup(name, graphs)]


add_multi_group(complex_groups, alkyl_groups, "alkyl")
# add_multi_group(complex_groups, alkyl_2_groups, "alkyl_2")
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
            ProxyGraph("C(=O)O", anchor=[0, 2]),
        ],
    ),
]
simple_groups += [
    ProxyGroup(
        "ether",
        [
            ProxyGraph("O"),
        ],
    ),
]

simple_groups += [
    get_alkyl_group("__alkene_end", max_carbons=4),
    get_alkyl_group("__alkene_intermediate", max_carbons=2),
]
simple_groups += [ProxyGroup("alkene", "{__alkene_intermediate}C=C{__alkene_end}")]

simple_groups += [
    get_alkyl_group("__allyl_end", carbon_counts=[1, 2, 3, 4]),
]
simple_groups += [ProxyGroup("allyl", ["CC={__allyl_end}"])]

simple_groups += [ProxyGroup("nitrogen_dioxide", "N(O)O")]
simple_groups += [ProxyGroup("halogen", ["F", "Cl", "Br", "I"])]
simple_groups += [ProxyGroup("nitrile", "C#N")]
simple_groups += [ProxyGroup("hydrogen_sulfite", "S(=O)(=O)O")]
simple_groups += [
    ProxyGroup("carbonyl", "C(=O)"),
]

common_groups = simple_groups + complex_groups

any_group = ProxyGroup("any", ["{{{}}}".format(g.name) for g in common_groups])
common_groups += [any_group]
