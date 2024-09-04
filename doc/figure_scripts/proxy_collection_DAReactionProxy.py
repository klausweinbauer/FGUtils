import numpy as np
import matplotlib.pyplot as plt

from fgutils.chem.its import get_its
from fgutils.vis import (
    GraphVisualizer,
    AnchorNodeLabelFormatter,
    plot_reaction,
    plot_its,
)
from fgutils.const import IS_LABELED_KEY, LABELS_KEY
from fgutils.parse import Parser
from fgutils.proxy_collection.diels_alder_proxy import DielsAlderProxy

plot = False


def get_group(groups, name):
    for g in groups:
        if g.name == name:
            return g
    raise ValueError("Group '{}' not found.".format(name))


def plot_core_graphs():
    parser = Parser(use_multigraph=True)
    vis = GraphVisualizer()
    vis.node_label_formatter.show_aam = True

    n = len(DielsAlderProxy.core_graphs)
    fig, ax = plt.subplots(1, n, figsize=(15, 6))
    for i, g in enumerate(DielsAlderProxy.core_graphs):
        vis.plot(parser(g.pattern), ax[i], title=g.name)
    plt.tight_layout()
    plt.savefig("doc/figures/proxy_collection_DAReactionProxy.png", bbox_inches="tight")
    if plot:
        plt.show()


def plot_graphs(graphs, name, cols=4, show_aam=False):
    n = len(graphs)
    rows = int(np.ceil(n / cols))
    parser = Parser()
    nlformatter = AnchorNodeLabelFormatter(show_aam=show_aam)
    vis = GraphVisualizer(node_label_formatter=nlformatter)
    fig, axs = plt.subplots(rows, cols, figsize=(20, 5 * rows))
    for r in range(rows):
        for c in range(cols):
            i = r * cols + c
            ax = axs[r, c]
            ax.axis("off")
            if i < len(graphs):
                graph = graphs[i]
                nlformatter.anchor = graph.anchor
                vis.plot(parser(graph.pattern), ax, title=graph.name)
                ax.axis("on")
                ax.set_xticks([])
                ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(
        "doc/figures/proxy_collection_DAReactionProxy_{}.png".format(name),
        bbox_inches="tight",
    )
    if plot:
        plt.show()


def plot_group(group, name, cols=4):
    plot_graphs(group.graphs, name, cols)


def plot_diene_and_dienophile_graphs():
    proxy = DielsAlderProxy()
    diene_group = get_group(proxy.groups, "s-cis_diene")
    dienophile_group = get_group(proxy.groups, "dienophile")
    plot_group(diene_group, "diene")
    plot_group(dienophile_group, "dienophile")


def resolve_group(group, groups):
    def _resolve(group, groups, graphs, parser):
        for g in group.graphs:
            _g = parser(g.pattern)
            if len(_g.nodes) == 1 and _g.nodes(data=True)[0][IS_LABELED_KEY]:
                _resolve(
                    get_group(groups, _g.nodes(data=True)[0][LABELS_KEY][0]),
                    groups,
                    graphs,
                    parser,
                )
            else:
                graphs.append(g)
        return

    graphs = []
    parser = Parser()
    _resolve(group, groups, graphs, parser)
    return graphs


def plot_edg_ewg():
    proxy = DielsAlderProxy()
    ewg = get_group(proxy.groups, "electron_withdrawing_group")
    edg = get_group(proxy.groups, "electron_donating_group")
    ewg_graphs = resolve_group(ewg, proxy.groups)
    edg_graphs = resolve_group(edg, proxy.groups)
    plot_graphs(ewg_graphs, "EWG")
    plot_graphs(edg_graphs, "EDG")


def plot_samples(neg_sample=False):
    if neg_sample:
        file_name = "DA_reactions_neg"
    else:
        file_name = "DA_reactions"
    proxy = DielsAlderProxy(neg_sample=neg_sample)

    data = [r for r in proxy]
    print("Created {} reactions.".format(len(data)))

    rows, cols = 10, 2
    step = np.max([len(data) / rows, 1])
    fig, ax = plt.subplots(rows, cols, figsize=(21, 4 * rows), width_ratios=[2, 1])
    for r in range(rows):
        _idx = int(r * step)
        if _idx >= len(data):
            break
        g, h = data[_idx]
        its = get_its(g, h)
        plot_reaction(g, h, ax[r, 0])
        ax[r, 0].set_title("Reaction")
        plot_its(its, ax[r, 1])
        ax[r, 1].set_title("ITS")
    plt.tight_layout()
    plt.savefig(
        "doc/figures/proxy_collection_DAReactionProxy_{}.png".format(file_name),
        bbox_inches="tight",
    )
    if plot:
        plt.show()


def plot_s_trans_diene():
    proxy = DielsAlderProxy()
    diene_group = get_group(proxy.groups, "s-trans_diene")
    plot_group(diene_group, "s-trans_diene")


plot_core_graphs()
plot_diene_and_dienophile_graphs()
plot_edg_ewg()
plot_samples()
plot_samples(True)
plot_s_trans_diene()
