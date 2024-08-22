import random
import collections

from fgutils import ReactionProxy
from fgutils.proxy import ProxyGroup, ProxyGraph
from fgutils.rdkit import graph_to_smiles

import fgutils.proxy_collection.common


class OutOfSampleError(Exception):
    def __init__(self, n):
        super().__init__()
        self.n = n

    def __str__(self):
        return "Could not find a new sample after {} tries.".format(self.n)


_groups = [
    ProxyGroup(
        "s-cis_diene_sample",
        [
            ProxyGraph("C1<2,1>C<1,2>C<2,1>CC1", anchor=[0, 3]),
            ProxyGraph("C1<2,1>C<1,2>C<2,1>CO1", anchor=[0, 3]),
            ProxyGraph("C<2,1>C1<1,2>C(CCCC1)<2,1>C", anchor=[0, 7]),
            ProxyGraph("C1<2,1>C<1,2>C<2,1>CCC1", anchor=[0, 3]),
            ProxyGraph("C1<2,1>C2<1,2>C(CCCC2)<2,1>CCC1", anchor=[0, 7]),
        ],
    ),
    ProxyGroup(
        "s-trans_diene_sample",
        [
            ProxyGraph("C1<2,1>C2<1,2>C<2,1>CCCC2CCC1", anchor=[0, 3]),
            ProxyGraph("C1CCC<2,1>C2<1,2>C1<2,1>CCCC2", anchor=[3, 6]),
            ProxyGraph("C<2,1>C1<1,2>C<2,1>CCCC1", anchor=[0, 3]),
        ],
    ),
    ProxyGroup("ethylene", ProxyGraph("C<2,1>C", anchor=[0, 1])),
    ProxyGroup(
        "electron_donating_group", pattern="{alkyl,aryl,trimethylsilanol,amine}"
    ),
    ProxyGroup(
        "electron_withdrawing_group",
        pattern="{alkohol,ether,aldehyde,ester,nitrile,nitrogen_dioxide,halogen,alkenyl,aryl}",
    ),
    ProxyGroup(
        "penta_diene_bridge",
        [
            ProxyGraph("C(=O)OC(=O)", anchor=[0, 3]),
            ProxyGraph("CCC", anchor=[0, 2]),
            ProxyGraph("COC", anchor=[0, 2]),
            ProxyGraph("S(=O)(=O)CC", anchor=[0, 4]),
            ProxyGraph("CCc3ccccc3", anchor=[0, 7]),
        ],
    ),
    ProxyGroup(
        "intra_mol_bridge",
        [
            ProxyGraph("{CC2,CC3}"),
            ProxyGraph("CC(=O)", anchor=[0, 1]),
            ProxyGraph("CCC(=O)", anchor=[0, 2]),
        ],
    ),
    ProxyGroup(
        "dienophile_mol",
        [
            ProxyGraph("{ethylene}"),
            ProxyGraph("C<2,1>C(Cl)OC(=O)C", anchor=[0, 1]),
            ProxyGraph("CC<2,1>CC#N", anchor=[1, 2]),
            ProxyGraph("C1<2,1>CCC=C1", anchor=[0, 1]),
            ProxyGraph("C1<2,1>CC=CC1", anchor=[0, 1]),
            ProxyGraph("C<3,2>C", anchor=[0, 1]),
        ],
    ),
    ProxyGroup(
        "diene_mol",
        [
            ProxyGraph("{s-cis_diene_sample}"),
            # ProxyGraph("C<2,1>C<1,2>C<2,1>C", anchor=[0, 3]),
        ],
    ),
    ProxyGroup(
        "dienophile",
        [
            ProxyGraph("{dienophile_mol}"),
            ProxyGraph("C1<2,1>C{penta_diene_bridge}1", anchor=[0, 1]),
            ProxyGraph("C<2,1>C{electron_withdrawing_group}", anchor=[0, 1]),
        ],
    ),
    ProxyGroup(
        "diene",
        [
            # ProxyGraph("{s-trans_diene_sample}"),
            ProxyGraph("{diene_mol}"),
            ProxyGraph("C<2,1>C<1,2>C<2,1>C{electron_donating_group}", anchor=[0, 3]),
            ProxyGraph("C<2,1>C<1,2>C({electron_donating_group})<2,1>C", anchor=[0, 4]),
        ],
    ),
]


class DielsAlderProxy:
    cores = [
        "C1<2,1>C<1,2>C(C)<2,1>C2C{intra_mol_bridge}C<0,1>2<2,1>C<0,1>1",
        "{diene}1<0,1>{dienophile}<0,1>1",
    ]

    def __init__(self, ignore_duplicates=True, enable_aam=True):
        self.proxies = []
        for core in self.cores:
            proxy = ReactionProxy(
                core, groups=fgutils.proxy_collection.common.common_groups + _groups
            )
            proxy.enable_aam = enable_aam
            self.proxies.append(proxy)
        self.ignore_duplicates = ignore_duplicates
        self.__dup_cache = collections.defaultdict(lambda: None)
        self.duplicate_retries = 100

    def generate(self):
        for _ in range(self.duplicate_retries):
            proxy = random.sample(self.proxies, 1)[0]
            g, h = next(proxy)
            if not self.ignore_duplicates:
                return g, h
            g_smiles = graph_to_smiles(g)
            h_smiles = graph_to_smiles(h)
            smpl_key = hash("{}>>{}".format(g_smiles, h_smiles))
            if self.__dup_cache[smpl_key] is None:
                self.__dup_cache[smpl_key] = 1
                return g, h
        raise OutOfSampleError(self.duplicate_retries)

    def __iter__(self):
        return self

    def __next__(self):
        return self.generate()
