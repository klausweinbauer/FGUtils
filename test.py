from fgutils.rdkit import smiles_to_graph, graph_to_smiles
from fgutils.chem.its import get_its
from fgutils.utils import split_its

smiles = "[C:1][C:2](=[O:3])[Cl:4].[O:5][C:6]>>[C:1][C:2](=[O:3])[O:5][C:6].[Cl:4]"
g, h = smiles_to_graph(smiles)
its = get_its(g, h)
print(next(iter(its.nodes(data=True))))
print(next(iter(its.edges(data=True))))

g, h = split_its(its)
smiles = "{}>>{}".format(graph_to_smiles(g), graph_to_smiles(h))
print(smiles)
