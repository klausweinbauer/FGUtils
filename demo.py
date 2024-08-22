from fgutils.proxy_collection.diels_alder_proxy import OutOfSampleError
from fgutils.proxy_collection import DielsAlderProxy
from fgutils.vis import plot_its, plot_reaction
from fgutils.chem.its import get_its
from fgutils.rdkit import graph_to_smiles

import matplotlib.pyplot as plt

proxy = DielsAlderProxy(enable_aam=True)

n = 6
fig, ax = plt.subplots(n, 2, dpi=50, figsize=(16, 9))
for i in range(n):
    try:
        g, h = next(proxy)
        rxn_smiles = "{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h))
        print(rxn_smiles)
        its = get_its(g, h)
        plot_reaction(g, h, ax[i, 0])
        ax[i, 0].axis("off")
        plot_its(its, ax[i, 1])
    except OutOfSampleError:
        break
plt.tight_layout()
plt.show()
