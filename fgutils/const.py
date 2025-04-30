import collections

AAM_KEY = "aam"
SYMBOL_KEY = "symbol"
BOND_KEY = "bond"
IS_LABELED_KEY = "is_labeled"
LABELS_KEY = "labels"
IDX_MAP_KEY = "idx_map"


BOND_ORDER_MAP = {
    "SINGLE": 1,
    "DOUBLE": 2,
    "TRIPLE": 3,
    "QUADRUPLE": 4,
    "AROMATIC": 1.5,
}


BOND_CHAR_MAP = {None: "∅", 0: "∅", 1: "—", 1.5: ":", 2: "=", 3: "≡"}


ATOM_COLORS = collections.defaultdict(
    lambda: "#000000",
    {
        "N": "#333399",
        "O": "#e61919",
        "H": "#555555",
        "S": "#666600",
        "F": "#996600",
        "I": "#660099",
        "P": "#996600",
        "Cl": "#008901",
        "Br": "#663333",
        "Na": "#0000ff",
        "K": "#008383",
        "Zn": "#663333",
        "Cu": "#663333",
        "Sn": "#336699",
        "Mg": "#006600",
        "B": "#008901",
    },
)
