from fgutils import ReactionProxy
from fgutils.proxy import ProxyGroup, ProxyGraph

from fgutils.proxy_collection.common import (
    common_groups,
    get_alkyl_group,
)


group_collection = [
    # override alkene group to prevent formation of diene in dienophiles
    get_alkyl_group("__alkene_end", carbon_counts=[1, 2, 3]),
    get_alkyl_group("__alkene_intermediate", carbon_counts=[1, 2]),
    ProxyGroup("alkene", "C{__alkene_intermediate}C=C{__alkene_end}"),
    # =================================
    ProxyGroup(
        "s-trans_diene_bridge",
        [
            ProxyGraph("{CC1}"),
            ProxyGraph("{CC2}"),
            ProxyGraph("CO", anchor=[0, 1]),
            ProxyGraph("CN", anchor=[0, 1]),
        ],
    ),
    ProxyGroup(
        "ring_bridge",
        [
            ProxyGraph("C"),
            ProxyGraph("O"),
            ProxyGraph("C({fg_col3})"),
            ProxyGraph("CC", anchor=[0, 1]),
            ProxyGraph("C1CCCCC1", anchor=[0, 5]),
        ],
    ),
    ProxyGroup(
        "s-cis_diene",
        [
            ProxyGraph(
                "C<2,1>C1<1,2>C(<2,1>C)C2OC1C(=C)C(=C)2",
                anchor=[0, 3],
                name="tetramethyldiene",
            ),
            ProxyGraph(
                "C1(Cl)<2,1>C(Cl)<1,2>C(Cl)<2,1>C(Cl)C(Cl)(Cl)1",
                anchor=[0, 6],
                name="perchlorocyclopentadiene",
            ),
            ProxyGraph(
                "C<2,1>C<1,2>C1<2,1>CCCC2=C1C=CC=C2",
                anchor=[0, 3],
                name="x-naphthalene",
            ),
            ProxyGraph("C<2,1>C<1,2>C<2,1>C", anchor=[0, 3], name="butadiene"),
            ProxyGraph(
                "C1CCC2<2,1>C<1,2>C<2,1>C3CCCCC3C2C1",
                anchor=[3, 6],
                name="dodecahydrophenanthrene",
            ),
            ProxyGraph(
                "C1<2,1>C<1,2>C<2,1>C{ring_bridge}1",
                anchor=[0, 3],
                name="hydronaphthalene_1",
            ),
            ProxyGraph(
                "C<2,1>C1<1,2>C(C{ring_bridge}C1)<2,1>C",
                anchor=[0, 6],
                name="dimethyl_ring_diene",
            ),
            ProxyGraph(
                "C1<2,1>C2<1,2>C(C{ring_bridge}C2)<2,1>CCC1",
                anchor=[0, 6],
                name="hydronaphthalene_2",
            ),
            ProxyGraph("C<2,1>C<1,2>C<2,1>C{electron_donating_group}", anchor=[0, 3]),
            ProxyGraph("C<2,1>C<1,2>C({electron_donating_group})<2,1>C", anchor=[0, 4]),
            ProxyGraph("C<2,1>C<1,2>C<2,1>C{s-trans_diene_mol}", anchor=[0, 3]),
            ProxyGraph("C<2,1>C<1,2>C({s-trans_diene_mol})<2,1>C", anchor=[0, 4]),
        ],
    ),
    ProxyGroup(
        "s-trans_diene_mol",
        [
            ProxyGraph(
                "C1=C2C=CCCC2CCC1", anchor=[5], name="hexahydronaphthalene_asym"
            ),
            ProxyGraph(
                "C1=C2C=CCCC2CCC1", anchor=[7], name="hexahydronaphthalene_asym"
            ),
            ProxyGraph(
                "C1=C2C=CCCC2CCC1", anchor=[8], name="hexahydronaphthalene_asym"
            ),
            ProxyGraph("C1CCC=C2C1=CCCC2", anchor=[1], name="hexahydronaphthalene_sym"),
            ProxyGraph("CC=C1C=CCCC1", anchor=[0]),
            ProxyGraph("C=C1C=CCCC1", anchor=[5], name="methylenecyclohexene"),
        ],
    ),
    ProxyGroup(
        "s-trans_diene",
        [
            ProxyGraph(
                "C1<2,1>C2<1,2>C<2,1>CCCC2CCC1",
                anchor=[0, 3],
                name="hexahydronaphthalene_asym",
            ),
            ProxyGraph(
                "C1CCC<2,1>C2<1,2>C1<2,1>CCCC2",
                anchor=[3, 6],
                name="hexahydronaphthalene_sym",
            ),
            ProxyGraph(
                "C<2,1>C1<1,2>C<2,1>CCCC1", anchor=[0, 3], name="methylenecyclohexene"
            ),
            ProxyGraph(
                "C1C<2,1>C2{s-trans_diene_bridge}C1C<2,1>C<1,2>2", anchor=[1, 5]
            ),
            ProxyGraph(
                "C1({electron_donating_group})<2,1>C2<1,2>C<2,1>CCCC2CCC1",
                anchor=[0, 4],
            ),
            ProxyGraph(
                "C1<2,1>C2<1,2>C({electron_donating_group})<2,1>CCCC2CCC1",
                anchor=[0, 4],
            ),
            ProxyGraph(
                "C1<2,1>C2<1,2>C<2,1>C({electron_donating_group})CCC2CCC1",
                anchor=[0, 3],
            ),
            ProxyGraph(
                "C1CCC({electron_donating_group})<2,1>C2<1,2>C1<2,1>CCCC2",
                anchor=[3, 7],
            ),
            ProxyGraph(
                "C<2,1>C1<1,2>C<2,1>C({electron_donating_group})CCC1", anchor=[0, 3]
            ),
            ProxyGraph(
                "{electron_donating_group}C<2,1>C1<1,2>C<2,1>CCCC1", anchor=[1, 4]
            ),
        ],
    ),
    ProxyGroup(
        "electron_donating_group",
        ["{aryl}", "{alkyl}", "{trimethylsilanol}", "{amine}"],
    ),
    ProxyGroup(
        "electron_withdrawing_group",
        [
            "{ether}",
            "{halogen}",
            "{aryl}",
            "{alkene}",
            "{acid}",
            "{ester}",
            "{carbonyl}",
            "{aldehyde}",
            "{hydrogen_sulfite}",
            "{nitrile}",
            "{nitrogen_dioxide}",
        ],
    ),
    ProxyGroup(
        "dienophile_bridge",
        [
            ProxyGraph("C(=O)OC(=O)", anchor=[0, 3], name="maleic_anhydride"),
            ProxyGraph("S(=O)(=O)CC", anchor=[0, 4], name="sulfolene"),
            ProxyGraph("CCc3ccccc3", anchor=[0, 7], name="dihydronaphthalene"),
            ProxyGraph("C{ring_bridge}C", anchor=[0, 2]),
        ],
    ),
    get_alkyl_group("fg_col1_alkyl", max_carbons=3),
    get_alkyl_group("fg_col2_alkyl", max_carbons=2),
    get_alkyl_group("fg_col3_alkyl", max_carbons=2),
    ProxyGroup("ester_fg", ["{ester}C", "{ester}CC"]),
    ProxyGroup("fg_col1", ["{ester_fg}", "{acid}", "{fg_col1_alkyl}"]),
    ProxyGroup("fg_col2", ["{ester_fg}", "{acid}", "{fg_col2_alkyl}", "{aryl}"]),
    ProxyGroup("fg_col3", ["{ester_fg}", "{acid}", "{fg_col3_alkyl}"]),
    ProxyGroup(
        "CC12",
        [
            ProxyGraph("C"),
            ProxyGraph("CC", anchor=[0, 1]),
            ProxyGraph("C({fg_col1})C", anchor=[0, 2]),
        ],
    ),
    ProxyGroup(
        "CC34",
        [
            ProxyGraph("CCC", anchor=[0, 2]),
            ProxyGraph("CCCC", anchor=[0, 3]),
            ProxyGraph("C({fg_col1})CC", anchor=[0, 3]),
        ],
    ),
    ProxyGroup(
        "intra_mol_bridge",
        [
            ProxyGraph("{CC34}"),
            ProxyGraph("{CC34}(=O)", anchor=[0]),
            ProxyGraph("{CC34}(=N)", anchor=[0]),
            ProxyGraph("{CC12}{ether}C", anchor=[0, 2]),
            ProxyGraph("{CC12}NC", anchor=[0, 2]),
            ProxyGraph("{CC12}{ester}", anchor=[0, 1]),
        ],
    ),
    ProxyGroup(
        "intra_mol_bridge_invalid",
        [
            ProxyGraph("{CC12}"),
            ProxyGraph("C({fg_col1})"),
        ],
    ),
    ProxyGroup(
        "dienophile",
        [
            ProxyGraph("C1<2,1>CC2CC1C=C2", anchor=[0, 1], name="norbomadiene"),
            ProxyGraph("C<2,1>C", anchor=[0, 1], name="ethylene"),
            ProxyGraph("C<2,1>C(Cl)OC(=O)C", anchor=[0, 1], name="chlorovinylacetate"),
            ProxyGraph("CC<2,1>CC#N", anchor=[1, 2], name="but-2-enenitrile"),
            ProxyGraph("C<3,2>C", anchor=[0, 1], name="acetylene"),
            ProxyGraph("C1<2,1>C{dienophile_bridge}1", anchor=[0, 1]),
            ProxyGraph("C<2,1>C{electron_withdrawing_group}", anchor=[0, 1]),
        ],
    ),
    ProxyGroup(
        "diene",
        [
            ProxyGraph("{s-cis_diene}"),
        ],
    ),
]


class DielsAlderProxy(ReactionProxy):
    """A proxy for the generation of Diels-Alder reaction samples and
    counter-samples. The proxy returns two graphs G and H as tuple. G is the
    reactant graph and H is the product graph. For a comprehensive description
    of the proxy configuration read section :ref:`diels-alder_reaction_proxy`.

    :param enable_aam: Flag to specify if the ``aam`` label is set in the
        result graphs. (Default = True)
    :param neg_sample: If set to true the proxy will exclusively generate
        negative samples, i.e., reactions where a Diels-Alder graph
        transformation rule is theoretically applicable but the reaction will
        never happen in reality. (Default = False)
    """

    core_graphs = [
        ProxyGraph(
            "{fg_col2}C1<2,1>C<1,2>C<2,1>C2{intra_mol_bridge}C<0,1>2<2,1>C<0,1>1{fg_col3}",
            name="Intra Molecular Center",
        ),
        ProxyGraph("{diene}1<0,1>{dienophile}<0,1>1", name="Inter Molecular Center"),
    ]

    def __init__(self, enable_aam=True, neg_sample=False):
        _groups = group_collection
        if neg_sample:
            _groups = group_collection.copy()
            _groups += [ProxyGroup("diene", ProxyGraph("{s-trans_diene}"))]
            _groups += [
                ProxyGroup("intra_mol_bridge", ProxyGraph("{intra_mol_bridge_invalid}"))
            ]
        core_group = ProxyGroup(
            "__DA_core__",
            self.core_graphs,
            unique=True,
        )
        super().__init__(
            core_group,
            common_groups + _groups,
            enable_aam=enable_aam,
        )
