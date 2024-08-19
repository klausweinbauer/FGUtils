import copy
import random
import collections
import rdkit.Chem as Chem

from fgutils import ReactionProxy
from fgutils.utils import print_graph
from fgutils.rdkit import graph_to_smiles, graph_to_mol

alkyl_group_names = [
    "methyl",
    "ethyl",
    "propyl",
    "butyl",
    "pentyl",
    "hexyl",
    "heptyl",
    "octyl",
    # "nonyl",
    # "decyl",
    # "undecyl",
    # "dodecyl",
    # "tridecyl",
    # "tetradecyl",
]

alkyl_groups = {name: "C" * (i + 1) for i, name in enumerate(alkyl_group_names)}
n_carbs = {
    "1C": "C",
    **{"{}C".format(i): [{"pattern": "C" * i, "anchor": [0, i - 1]}] for i in range(2, 7)},
}

base_groups = {**alkyl_groups, **n_carbs}


class DielsAlderProxy:
    standard_groups = {
        "alkyl": "{" + ",".join(alkyl_group_names) + "}",
        "aryl": ["c1ccccc1", "c1c(C)cccc1", "c1cc(C)ccc1"],
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
            {"pattern": "C(=O)OC(=O)", "anchor": [0, 3]},
            {"pattern": "CCC", "anchor": [0, 2]},
            {"pattern": "COC", "anchor": [0, 2]},
            {"pattern": "S(=O)(=O)CC", "anchor": [0, 4]},
            {"pattern": "CCc3ccccc3", "anchor": [0, 7]},
        ],
        "dienophile": [
            # {"pattern": "C<2,1>C(Cl)OC(=O)C", "anchor": [0, 1]},
            # {"pattern": "CC<2,1>CC#N", "anchor": [1, 2]},
            # {"pattern": "C1<2,1>CCC=C1", "anchor": [0, 1]},
            {"pattern": "C1<2,1>C{penta_diene_bridge}1", "anchor": [0, 1]},
            {"pattern": "C1<2,1>C{electron_withdrawing_group}", "anchor": [0, 1]},
        ],
        "diene": [{"pattern": "C<2,1>C<1,2>C<2,1>C", "anchor": [0, 3]}],
    }

    cores = [
        "{electron_donating_group}C1<2,1>C<1,2>C<2,1>C<0,1>{dienophile}<0,1>1",
        "C1<2,1>C<1,2>C({electron_donating_group})<2,1>C<0,1>{dienophile}<0,1>1",
        "C1<2,1>C<1,2>C(C)<2,1>C2C{2C,3C}C<0,1>2<2,1>C<0,1>1",
        # Needs MultiGraph "{dienes}1<0,1>{dienophiles}<0,1>1",
    ]

    def __init__(self, ignore_duplicates=True, enable_aam=True):
        self.proxies = []
        for core in self.cores:
            config = {
                "core": core,
                "groups": {**base_groups, **self.standard_groups, **self.super_groups},
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
