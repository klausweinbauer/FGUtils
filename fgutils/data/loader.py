import csv

from fgutils.its import ITS
from fgutils.rdkit import reaction_smiles_to_graph

from rdkit import RDLogger


class CSVLoader:
    """The base class for loading CSV files. This class is a generator and can
    directly be used in a for loop. The returned column values are in the order
    in which they were specified in the columns argument.

    :param file: File to the CSV file.
    :param columns: Single column or list of columns to read from CSV file.
    :param delimiter: (optional) The CSV delimiter. Default: ","
    """

    def __init__(self, file, columns: str | list[str], delimiter=","):
        self.file = file
        self.columns = columns if isinstance(columns, list) else [columns]
        self.delimiter = delimiter

    def __get_col_indices(self, headers: list[str]) -> list[int]:
        indices = []
        for col in self.columns:
            found = False
            for i, header in enumerate(headers):
                if header == col:
                    indices.append(i)
                    found = True
            if not found:
                raise ValueError(
                    "Column '{}' not found. Available columns are {}.".format(
                        col, headers
                    )
                )
        return indices

    def __iter__(self):
        with open(self.file, "r") as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            header = next(reader)
            col_indices = self.__get_col_indices(header)
            single_col = len(self.columns) == 1
            for line in reader:
                sel = [line[i] for i in col_indices]
                yield sel[0] if single_col else sel


class CSVReactionLoader(CSVLoader):
    """Load SMILES reactions from CSV file.

    :param file: The path to the CSV file.
    :param reaction_col: The column where to find the SMILES reaction.
    :param id_col: (optional) An id column. Default: None

    :param mode: (optional) Defines how reactions and molecules are loaded.
        Available modes are:

        * ``"GH"``: Loads the reactant and product molecules as graphs. The
          reaction is returned in the form G \u2192 H.
        * ``"SMILES"``: Loads the reaction as SMILES.
        * ``"ITS"``: Loads the reaction as :py:class:`~fgutils.its.ITS` graph.
          The SMILES must have annotated atom-atom maps.
    """

    def __init__(self, file, reaction_col, id_col=None, mode="GH", delimiter=","):
        RDLogger.DisableLog("rdApp.*")
        columns = []
        if id_col is not None:
            columns.append(id_col)
        columns.append(reaction_col)
        self.mode = mode
        super(CSVReactionLoader, self).__init__(file, columns, delimiter=delimiter)

    def __iter__(self):  # type: ignore
        for row in super().__iter__():
            if self.mode == "GH":
                if isinstance(row, list):
                    try:
                        g, h = reaction_smiles_to_graph(row.pop())
                        row.extend([g, h])  # type: ignore
                    except ValueError:
                        row.extend([None, None])
                    yield row
                else:
                    try:
                        g = reaction_smiles_to_graph(row)
                    except ValueError:
                        g = None
                    yield g
            elif self.mode == "SMILES":
                yield row
            elif self.mode == "ITS":
                if isinstance(row, list):
                    try:
                        its = ITS.from_smiles(row.pop())
                        row.append(its)  # type: ignore
                    except ValueError:
                        row.append(None)
                    yield row
                else:
                    try:
                        g = ITS.from_smiles(row)
                    except ValueError:
                        g = None
                    yield g
            else:
                raise ValueError(
                    "Reaction reading mode '{}' is not implemented.".format(self.mode)
                )
