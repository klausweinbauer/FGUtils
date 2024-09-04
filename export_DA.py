import json

from fgutils.proxy_collection import DielsAlderProxy
from fgutils.rdkit import graph_to_smiles


def export(file_name, neg_sample):
    idx_fmt = "DA_{}"
    if neg_sample:
        idx_fmt = "DA_neg_{}"

    dataset = []
    proxy = DielsAlderProxy(neg_sample=neg_sample)
    for i, r in enumerate(proxy):
        g, h = r
        smiles = "{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h))
        dataset.append({"index": idx_fmt.format(i), "reaction": smiles})

    with open(file_name, "w") as f:
        json.dump(dataset, f, indent=4)


export("Diels-Alder_synthetic_data.json", neg_sample=False)
export("Diels-Alder_synthetic_negative_data.json", neg_sample=True)
