FGUtils is a collection of utility functions for querying functional groups in
molecules from their molecular graph representation.

## Dependencies
- Python (>= 3.11)
- numpy (>= 1.26.4)
- networkx (>= 3.2.1)
- rdkit (>= 2023.09.4 optional)

## Installation
You can install [FGUtils](https://pypi.org/project/fgutils/) using pip.
```
pip install fgutils
```

## Getting Started
For a comprehensive description of FGUtils features read through the
[documentation](https://klausweinbauer.github.io/FGUtils/). However, querying
the functional groups for a molecule like acetylsalicylic acid is as simple as
running the following:
```
>>> from fgutils import FGQuery
>>> 
>>> smiles = "O=C(C)Oc1ccccc1C(=O)O" # acetylsalicylic acid
>>> query = FGQuery()
>>> query.get(smiles)
[('ester', [0, 1, 3]), ('carboxylic_acid', [10, 11, 12])]
```

The output is a list of tuples containing the functional group name and the
corresponding atom indices.
