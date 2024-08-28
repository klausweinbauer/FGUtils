from fgutils.proxy import Proxy, ProxyGroup, LabelSampler
from fgutils.proxy_collection.common import common_groups
from fgutils.vis import plot_its, plot_reaction
from fgutils.chem.its import get_its
from fgutils.rdkit import graph_to_smiles

import matplotlib.pyplot as plt

core = ProxyGroup("core", pattern="{alkyl}(=O)O", sample_unique=True)
proxy = Proxy(core, groups=common_groups)
proxy.label_sampler = LabelSampler()

for g in proxy:
    print(g)
