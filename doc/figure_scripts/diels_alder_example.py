import random
import matplotlib.pyplot as plt
from fgutils.proxy import ProxyGroup, ProxyGraph, ReactionProxy
from fgutils.proxy_collection.common import common_groups
from fgutils.vis import plot_reaction, plot_its
from fgutils.chem.its import get_its


electron_donating_group = ProxyGroup(
    "electron_donating_group",
    ["{methyl}", "{ethyl}", "{propyl}", "{aryl}", "{amine}"],
)
electron_withdrawing_group = ProxyGroup(
    "electron_withdrawing_group",
    ["{alkohol}", "{ether}", "{aldehyde}", "{ester}", "{nitrile}"],
)
diene_group = ProxyGroup(
    "diene",
    ProxyGraph("C<2,1>C<1,2>C<2,1>C{electron_donating_group}", anchor=[0, 3]),
)
dienophile_group = ProxyGroup(
    "dienophile",
    ProxyGraph("C<2,1>C{electron_withdrawing_group}", anchor=[0, 1]),
)
groups = common_groups + [
    electron_donating_group,
    electron_withdrawing_group,
    diene_group,
    dienophile_group,
]

proxy = ReactionProxy("{diene}1<0,1>{dienophile}<0,1>1", groups)

n = 4
fig, ax = plt.subplots(n, 2, width_ratios=[2, 1], figsize=(20, n * 4))
for i, (g, h) in enumerate(random.sample(list(proxy), n)):
    plot_reaction(g, h, ax[i, 0], title="Reaction")
    plot_its(get_its(g, h), ax[i, 1], title="ITS Graph")
plt.tight_layout()
plt.savefig(
    "doc/figures/diels_alder_example.png", bbox_inches="tight", transparent=False
)
plt.show()
