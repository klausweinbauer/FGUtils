from fgutils.proxy_collection.common import common_groups
from fgutils.proxy import MolProxy

from fgutils.vis import PdfWriter, GraphVisualizer

top_level_groups = [
    "alkyl",
    "aryl",
    "allyl",
    "alkene",
    "amine",
    "ester",
    "ether",
    "carbonyl",
]

vis = GraphVisualizer()


def plot(value, ax, length=0, group="", index=0):
    vis.plot(value, ax, title="{} {}/{}".format(group, index + 1, length))


pdf_writer = PdfWriter("common_groups.pdf", plot_fn=plot, cols=3)

data = common_groups
for tlg in top_level_groups:
    proxy = MolProxy("{{{}}}".format(tlg), common_groups)
    data = list(proxy)
    print("Group '{}' has {} molecules.".format(tlg, len(data)))
    pdf_writer.plot(data, group=tlg, length=len(data))

pdf_writer.close()
