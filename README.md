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

smiles = "O=C(C)Oc1ccccc1C(=O)O"  # acetylsalicylic acid
query = fgutils.FGQuery(use_smiles=True) # requires rdkit to be installed
groups = query.get(smiles)
print(groups)
```

The output is a list of tuples containing the functional group name and the corresponding atom indices.

```
[('ester', [0, 1, 3]), ('carboxylic_acid', [10, 11, 12])]
```
