from fgutils import ReactionProxy
from fgutils.proxy import ProxyGroup, ProxyGraph

from fgutils.proxy_collection.common import common_groups


group_collection = [
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
        "s-cis_diene",
        [
            ProxyGraph(
                "C1(Cl)<2,1>C(Cl)<1,2>C(Cl)<2,1>C(Cl)C(Cl)(Cl)1",
                anchor=[0, 6],
                name="perchlorocyclopentadiene",
            ),
            ProxyGraph("C<2,1>C<1,2>C<2,1>C", anchor=[0, 3], name="butadiene"),
            ProxyGraph("C1<2,1>C<1,2>C<2,1>CC1", anchor=[0, 3], name="cyclopentadiene"),
            ProxyGraph("C1<2,1>C<1,2>C<2,1>CO1", anchor=[0, 3], name="furan"),
            ProxyGraph(
                "C<2,1>C1<1,2>C(CCCC1)<2,1>C",
                anchor=[0, 7],
                name="dimethylenecyclohexane",
            ),
            ProxyGraph("C1<2,1>C<1,2>C<2,1>CCC1", anchor=[0, 3], name="cyclohexadiene"),
            ProxyGraph(
                "C1<2,1>C2<1,2>C(CCCC2)<2,1>CCC1",
                anchor=[0, 7],
                name="hexahydronaphthalene",
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
            "C{alkenyl}",
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
            ProxyGraph("CCC", anchor=[0, 2], name="cyclopentene"),
            ProxyGraph("CCCC", anchor=[0, 2], name="cyclohexene"),
            ProxyGraph("COC", anchor=[0, 2], name="dihydrofuran"),
            ProxyGraph("S(=O)(=O)CC", anchor=[0, 4], name="sulfolene"),
            ProxyGraph("CCc3ccccc3", anchor=[0, 7], name="dihydronaphthalene"),
        ],
    ),
    ProxyGroup(
        "intra_mol_bridge",
        [
            ProxyGraph("{CC2}"),
            ProxyGraph("{CC3}"),
            ProxyGraph("CC(=O)", anchor=[0, 1]),
            ProxyGraph("CCC(=O)", anchor=[0, 2]),
            ProxyGraph("CC(=N)", anchor=[0, 1]),
            ProxyGraph("CCC(=N)", anchor=[0, 2]),
        ],
    ),
    ProxyGroup(
        "intra_mol_bridge_invalid",
        [
            ProxyGraph("{CC1}"),
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
    core_patterns = [
        "C1<2,1>C<1,2>C(C)<2,1>C2C{intra_mol_bridge}C<0,1>2<2,1>C<0,1>1",
        "{diene}1<0,1>{dienophile}<0,1>1",
    ]

    def __init__(self, enable_aam=True, unique=True, neg_sample=False):
        _groups = group_collection
        if neg_sample:
            _groups = group_collection.copy()
            _groups += [ProxyGroup("diene", ProxyGraph("{s-trans_diene}"))]
            _groups += [
                ProxyGroup("intra_mol_bridge", ProxyGraph("{intra_mol_bridge_invalid}"))
            ]
        core_graphs = []
        for core_pattern in self.core_patterns:
            core_graphs.append(ProxyGraph(core_pattern))
        core_group = ProxyGroup(
            "__DA_core__",
            core_graphs,
            unique=unique,
        )
        super().__init__(
            core_group,
            common_groups + _groups,
            enable_aam=enable_aam,
        )
