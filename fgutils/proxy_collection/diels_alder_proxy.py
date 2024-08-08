import copy
import random
import collections
import rdkit.Chem as Chem

from fgutils import ReactionProxy
from fgutils.utils import print_graph
from fgutils.rdkit import graph_to_smiles, graph_to_mol


class DielsAlderProxy:
    standard_groups = {
        "alkyl": ["C" for _ in range(4)],
        "aryl": ["c1ccccc1", "c1(C)ccccc1", "c1c(C)cccc1", "c1cc(C)ccc1"],
        "alkohol": ["CO"],
        "trimethylsilanol": "OSi(C)(C)C",
        "amine": ["NC", "N(C)C"],
        "alkenyl": ["C=CC", "C=C(C)C", "C=C(C)(C)C"],
        "ether": ["COC{}".format("C" * i) for i in range(2)] + ["C(Cl)OC(=O)C"],
        "aldehyde": "C(=O)O",
        "ester": ["C(=O)OC{}".format("C" * i) for i in range(2)],
        "nitrile": "C#N",
        "nitrogen_dioxide": "N(O)O",
        "halogen": ["F", "Cl", "Br", "I"],
    }

    super_groups = {
        "electron_donating_group": "{alkyl,aryl,trimethylsilanol,amine}",
        "electron_withdrawing_group": "{alkohol,ether,aldehyde,ester,"
        + "nitrile,nitrogen_dioxide,halogen,alkenyl,aryl}",
        "penta_diene_bridge": [
            "C(=O)OC(=O)",
            "CCC",
            "COC",
            "S(=O)(=O)CC",
            "Cc3cccccc3",
        ],
    }

    cores = [
        "{electron_donating_group}C1<2,1>C<1,2>C<2,1>C<0,1>"
        + "C<2,1>C({electron_withdrawing_group})<0,1>1",
        "C1<2,1>C<1,2>C({electron_donating_group})<2,1>C<0,1>"
        + "C<2,1>C({electron_withdrawing_group})<0,1>1",
        # Pentadiene as dienophile
        # "{electron_donating_group}C1<2,1>C<1,2>C<2,1>C<0,1>"
        # + "C2<2,1>C<0,1>1{penta_diene_bridge}2",
        #"C1<2,1>C<1,2>C<2,1>C<0,1>"
        #+ "C2<2,1>C<0,1>1{penta_diene_bridge}2",
    ]

    def __init__(self, ignore_duplicates=True, enable_aam=True):
        self.proxies = []
        for core in self.cores:
            config = {
                "core": core,
                "groups": {**self.standard_groups, **self.super_groups},
            }
            self.proxies.append(ReactionProxy(config, enable_aam=enable_aam))
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
        raise RuntimeError(
            "Could not find a new sample after {} tries.".format(self.duplicate_retries)
        )

    def __iter__(self):
        return self

    def __next__(self):
        return self.generate()


proxy = DielsAlderProxy(enable_aam=False)

for _ in range(100):
    try:
        g, h = next(proxy)
        print("{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h)))
    except RuntimeError:
        break
