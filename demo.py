from fgutils.proxy_collection.diels_alder_proxy import DielsAlderProxy
from fgutils.vis import GraphVisualizer, plot_reaction, plot_its, PdfWriter
from fgutils.chem.its import get_its


def plot(data, ax, index=0, **kwargs):
    g, h = data
    its = get_its(g, h)
    plot_reaction(g, h, ax[0])
    ax[0].set_title("Reaction | Index: {}".format(index))
    plot_its(its, ax[1])
    ax[1].set_title("ITS | Index: {}".format(index))


neg_sample = True

if neg_sample:
    file_name = "DA_reactions_neg.pdf"
else:
    file_name = "DA_reactions.pdf"

proxy = DielsAlderProxy(neg_sample=neg_sample)
vis = GraphVisualizer()
pdf_writer = PdfWriter(
    file_name, plot, max_pages=50, plot_per_row=True, width_ratios=[2, 1]
)

data = list(proxy)
print("Generated {} reactions.".format(len(data)))

pdf_writer.plot(data)
pdf_writer.close()
