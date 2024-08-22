import matplotlib.pyplot as plt
from fgutils.proxy import ProxyGroup, ProxyGraph, ReactionProxy
from fgutils.proxy_collection.common import common_groups
from fgutils.vis import plot_reaction


electron_donating_group = ProxyGroup(
    "electron_donating_group", pattern="{alkyl,aryl,amine}"
)
electron_withdrawing_group = ProxyGroup(
    "electron_withdrawing_group",
    pattern="{alkohol,ether,aldehyde,ester,nitrile}",
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

r, c = 3, 2
fig, ax = plt.subplots(r, c, dpi=400)
for ri in range(r):
    for ci in range(c):
        g, h = next(proxy)
        ax[ri, ci].axis("off")
        plot_reaction(g, h, ax[ri, ci])

plt.tight_layout()
plt.savefig("doc/figures/diels_alder_example.png", bbox_inches="tight", transparent=True)
plt.show()
