import collections

AAM_KEY = "aam"
SYMBOL_KEY = "symbol"
BOND_KEY = "bond"
IS_LABELED_KEY = "is_labeled"
LABELS_KEY = "labels"
IDX_MAP_KEY = "idx_map"


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
