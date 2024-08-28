import numpy as np
from fgutils.proxy import Proxy, ProxyGroup
from fgutils.proxy_collection.common import common_groups
from fgutils.proxy_collection.diels_alder_proxy import DielsAlderProxy
from fgutils.vis import GraphVisualizer, plot_reaction, plot_its
from fgutils.chem.its import get_its
from fgutils.rdkit import graph_to_smiles

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

neg_sample = True

if neg_sample:
    file_name = "DA_reactions_neg.pdf"
else:
    file_name = "DA_reactions.pdf"
proxy = DielsAlderProxy(neg_sample=neg_sample)
vis = GraphVisualizer()

data = [r for r in proxy]
print("Created {} reactions.".format(len(data)))

pages = 10
rows, cols = 7, 2
n_print = pages * rows
step = np.max([len(data) / n_print, 1])
pp = PdfPages(file_name)
done = False
for p in range(pages):
    print("Pring page {} of {}".format(p + 1, pages))
    fig, ax = plt.subplots(rows, cols, figsize=(21, 29.7), width_ratios=[2, 1])
    for r in range(rows):
        _idx = int((p * rows + r) * step)
        if _idx >= len(data):
            done = True
            break
        g, h = data[_idx]
        its = get_its(g, h)
        plot_reaction(g, h, ax[r, 0])
        ax[r, 0].set_title("Reaction | Index: {}".format(_idx))
        plot_its(its, ax[r, 1])
        ax[r, 1].set_title("ITS | Index: {}".format(_idx))
    plt.tight_layout()
    pp.savefig(fig, bbox_inches="tight", pad_inches=1)
    plt.close()
    if done:
        break
pp.close()
