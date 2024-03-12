FGUtils is a collection of utility functions for querying functional groups in molecules from their molecular graph representation.

## Dependencies
- Python (>= 3.11)
- numpy (>= 1.26.4)
- networkx (>= 3.2.1)
- rdkit (>= 2023.09.4 optional)

## Installation
You can install FGUtils using pip.
```
pip install fgutils
```

## Getting Started
A simple example querying the functional groups for acetylsalicylic acid.
```
import fgutils
import rdkit.Chem.rdmolfiles as rdmolfiles
from fgutils.utils import mol_to_graph

smiles = "O=C(C)Oc1ccccc1C(=O)O" # acetylsalicylic acid
mol_graph = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
index_map, groups = fgutils.get_functional_groups(mol_graph)
print(index_map, groups)
```

The output is an index map and a list of functional group names. 

```
{0: [0, 1, 2], 1: [0, 1, 2], 3: [2, 3]}
['carbonyl', 'ketone', 'ester', 'ether']
```

The index map maps atom numbers to functional groups. The key of the map is the atom index and the value is the list of indices from the functional group list. E.g. atom 0 (oxygen) is in functional groups carbonyl, ketone, and ester.
