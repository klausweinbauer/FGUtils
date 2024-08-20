import copy
import random
import collections
import rdkit.Chem as Chem

from fgutils import ReactionProxy
from fgutils.proxy import ProxyGroup, ProxyGraph
from fgutils.utils import print_graph
from fgutils.rdkit import graph_to_smiles, graph_to_mol

import fgutils.proxy_collection.common


class OutOfSampleError(Exception):
    def __init__(self, n):
        super().__init__()
        self.n = n

    def __str__(self):
        return "Could not find a new sample after {} tries.".format(self.n)


_groups = [
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
        "simple_dienophile",
        [
            ProxyGraph("C<2,1>C(Cl)OC(=O)C", anchor=[0, 1]),
            ProxyGraph("CC<2,1>CC#N", anchor=[1, 2]),
            ProxyGraph("C1<2,1>CCC=C1", anchor=[0, 1]),
        ],
    ),
    ProxyGroup(
        "dienophile",
        [
            ProxyGraph("{simple_dienophile}"),
            ProxyGraph("C1<2,1>C{penta_diene_bridge}1", anchor=[0, 1]),
            ProxyGraph("C1<2,1>C{electron_withdrawing_group}", anchor=[0, 1]),
        ],
    ),
    ProxyGroup("diene", [ProxyGraph("C<2,1>C<1,2>C<2,1>C", anchor=[0, 3])]),
    ProxyGroup(
        "intra_mol_diels_alder_bridge",
        [
            ProxyGraph("{CC2,CC3}"),
            ProxyGraph("CC(=O)", anchor=[0, 1]),
            ProxyGraph("CCC(=O)", anchor=[0, 2]),
        ],
    ),
]


class DielsAlderProxy:
    cores = [
        "{electron_donating_group}C1<2,1>C<1,2>C<2,1>C<0,1>{dienophile}<0,1>1",
        "C1<2,1>C<1,2>C({electron_donating_group})<2,1>C<0,1>{dienophile}<0,1>1",
        "C1<2,1>C<1,2>C(C)<2,1>C2C{intra_mol_diels_alder_bridge}C<0,1>2<2,1>C<0,1>1",
        # Needs MultiGraph "{dienes}1<0,1>{dienophiles}<0,1>1",
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
