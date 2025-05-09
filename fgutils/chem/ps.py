import re

from fgutils.utils import to_non_aromatic_symbol

# Electronegativity data from:
# https://sciencenotes.org/list-of-electronegativity-values-of-the-elements/
# Valence Electron data from:
# https://www.schoolmykids.com/learn/periodic-table/valence-electrons-of-all-the-elements

atomic_data = {
    "H": {
        "num": 1,
        "name": "Hydrogen",
        "electronegativity": 2.2,
        "valence_electrons": ["1s1"],
    },
    "He": {
        "num": 2,
        "name": "Helium",
        "electronegativity": 0,
        "valence_electrons": ["1s2"],
    },
    "Li": {
        "num": 3,
        "name": "Lithium",
        "electronegativity": 0.98,
        "valence_electrons": ["2s1"],
    },
    "Be": {
        "num": 4,
        "name": "Beryllium",
        "electronegativity": 1.57,
        "valence_electrons": ["2s2"],
    },
    "B": {
        "num": 5,
        "name": "Boron",
        "electronegativity": 2.04,
        "valence_electrons": ["2s2", "2p1"],
    },
    "C": {
        "num": 6,
        "name": "Carbon",
        "electronegativity": 2.55,
        "valence_electrons": ["2s2", "2p2"],
    },
    "N": {
        "num": 7,
        "name": "Nitrogen",
        "electronegativity": 3.04,
        "valence_electrons": ["2s2", "2p3"],
    },
    "O": {
        "num": 8,
        "name": "Oxygen",
        "electronegativity": 3.44,
        "valence_electrons": ["2s2", "2p4"],
    },
    "F": {
        "num": 9,
        "name": "Fluorine",
        "electronegativity": 3.98,
        "valence_electrons": ["2s2", "2p5"],
    },
    "Ne": {
        "num": 10,
        "name": "Neon",
        "electronegativity": 0,
        "valence_electrons": ["2s2", "2p6"],
    },
    "Na": {
        "num": 11,
        "name": "Sodium",
        "electronegativity": 0.93,
        "valence_electrons": ["3s1"],
    },
    "Mg": {
        "num": 12,
        "name": "Magnesium",
        "electronegativity": 1.31,
        "valence_electrons": ["3s2"],
    },
    "Al": {
        "num": 13,
        "name": "Aluminum",
        "electronegativity": 1.61,
        "valence_electrons": ["3s2", "3p1"],
    },
    "Si": {
        "num": 14,
        "name": "Silicon",
        "electronegativity": 1.9,
        "valence_electrons": ["3s2", "3p2"],
    },
    "P": {
        "num": 15,
        "name": "Phosphorus",
        "electronegativity": 2.19,
        "valence_electrons": ["3s2", "3p3"],
    },
    "S": {
        "num": 16,
        "name": "Sulfur",
        "electronegativity": 2.58,
        "valence_electrons": ["3s2", "3p4"],
    },
    "Cl": {
        "num": 17,
        "name": "Chlorine",
        "electronegativity": 3.16,
        "valence_electrons": ["3s2", "3p5"],
    },
    "Ar": {
        "num": 18,
        "name": "Argon",
        "electronegativity": 0,
        "valence_electrons": ["3s2", "3p6"],
    },
    "K": {
        "num": 19,
        "name": "Potassium",
        "electronegativity": 0.82,
        "valence_electrons": ["4s1"],
    },
    "Ca": {
        "num": 20,
        "name": "Calcium",
        "electronegativity": 1.0,
        "valence_electrons": ["4s2"],
    },
    "Sc": {
        "num": 21,
        "name": "Scandium",
        "electronegativity": 1.36,
        "valence_electrons": ["3d1", "4s2"],
    },
    "Ti": {
        "num": 22,
        "name": "Titanium",
        "electronegativity": 1.54,
        "valence_electrons": ["3d2", "4s2"],
    },
    "V": {
        "num": 23,
        "name": "Vanadium",
        "electronegativity": 1.63,
        "valence_electrons": ["3d3", "4s2"],
    },
    "Cr": {
        "num": 24,
        "name": "Chromium",
        "electronegativity": 1.66,
        "valence_electrons": ["3d5", "4s1"],
    },
    "Mn": {
        "num": 25,
        "name": "Manganese",
        "electronegativity": 1.55,
        "valence_electrons": ["3d5", "4s2"],
    },
    "Fe": {
        "num": 26,
        "name": "Iron",
        "electronegativity": 1.83,
        "valence_electrons": ["3d6", "4s2"],
    },
    "Co": {
        "num": 27,
        "name": "Cobalt",
        "electronegativity": 1.88,
        "valence_electrons": ["3d7", "4s2"],
    },
    "Ni": {
        "num": 28,
        "name": "Nickel",
        "electronegativity": 1.91,
        "valence_electrons": ["3d8", "4s2"],
    },
    "Cu": {
        "num": 29,
        "name": "Copper",
        "electronegativity": 1.9,
        "valence_electrons": ["3d10", "4s1"],
    },
    "Zn": {
        "num": 30,
        "name": "Zinc",
        "electronegativity": 1.65,
        "valence_electrons": ["3d10", "4s2"],
    },
    "Ga": {
        "num": 31,
        "name": "Gallium",
        "electronegativity": 1.81,
        "valence_electrons": ["4s2", "4p1"],
    },
    "Ge": {
        "num": 32,
        "name": "Germanium",
        "electronegativity": 2.01,
        "valence_electrons": ["4s2", "4p2"],
    },
    "As": {
        "num": 33,
        "name": "Arsenic",
        "electronegativity": 2.18,
        "valence_electrons": ["4s2", "4p3"],
    },
    "Se": {
        "num": 34,
        "name": "Selenium",
        "electronegativity": 2.55,
        "valence_electrons": ["4s2", "4p4"],
    },
    "Br": {
        "num": 35,
        "name": "Bromine",
        "electronegativity": 2.96,
        "valence_electrons": ["4s2", "4p5"],
    },
    "Kr": {
        "num": 36,
        "name": "Krypton",
        "electronegativity": 3.0,
        "valence_electrons": ["4s2", "4p6"],
    },
    "Rb": {
        "num": 37,
        "name": "Rubidium",
        "electronegativity": 0.82,
        "valence_electrons": ["5s1"],
    },
    "Sr": {
        "num": 38,
        "name": "Strontium",
        "electronegativity": 0.95,
        "valence_electrons": ["5s2"],
    },
    "Y": {
        "num": 39,
        "name": "Yttrium",
        "electronegativity": 1.22,
        "valence_electrons": ["4d1", "5s2"],
    },
    "Zr": {
        "num": 40,
        "name": "Zirconium",
        "electronegativity": 1.33,
        "valence_electrons": ["4d2", "5s2"],
    },
    "Nb": {
        "num": 41,
        "name": "Niobium",
        "electronegativity": 1.6,
        "valence_electrons": ["4d4", "5s1"],
    },
    "Mo": {
        "num": 42,
        "name": "Molybdenum",
        "electronegativity": 2.16,
        "valence_electrons": ["4d5", "5s1"],
    },
    "Tc": {
        "num": 43,
        "name": "Technetium",
        "electronegativity": 1.9,
        "valence_electrons": ["4d5", "5s2"],
    },
    "Ru": {
        "num": 44,
        "name": "Ruthenium",
        "electronegativity": 2.2,
        "valence_electrons": ["4d7", "5s1"],
    },
    "Rh": {
        "num": 45,
        "name": "Rhodium",
        "electronegativity": 2.28,
        "valence_electrons": ["4d8", "5s1"],
    },
    "Pd": {
        "num": 46,
        "name": "Palladium",
        "electronegativity": 2.2,
        "valence_electrons": ["4d10"],
    },
    "Ag": {
        "num": 47,
        "name": "Silver",
        "electronegativity": 1.93,
        "valence_electrons": ["4d10", "5s1"],
    },
    "Cd": {
        "num": 48,
        "name": "Cadmium",
        "electronegativity": 1.69,
        "valence_electrons": ["4d10", "5s2"],
    },
    "In": {
        "num": 49,
        "name": "Indium",
        "electronegativity": 1.78,
        "valence_electrons": ["5s2", "5p1"],
    },
    "Sn": {
        "num": 50,
        "name": "Tin",
        "electronegativity": 1.96,
        "valence_electrons": ["5s2", "5p2"],
    },
    "Sb": {
        "num": 51,
        "name": "Antimony",
        "electronegativity": 2.05,
        "valence_electrons": ["5s2", "5p3"],
    },
    "Te": {
        "num": 52,
        "name": "Tellurium",
        "electronegativity": 2.1,
        "valence_electrons": ["5s2", "5p4"],
    },
    "I": {
        "num": 53,
        "name": "Iodine",
        "electronegativity": 2.66,
        "valence_electrons": ["5s2", "5p5"],
    },
    "Xe": {
        "num": 54,
        "name": "Xenon",
        "electronegativity": 2.6,
        "valence_electrons": ["5s2", "5p6"],
    },
    "Cs": {
        "num": 55,
        "name": "Cesium",
        "electronegativity": 0.79,
        "valence_electrons": ["6s1"],
    },
    "Ba": {
        "num": 56,
        "name": "Barium",
        "electronegativity": 0.89,
        "valence_electrons": ["6s2"],
    },
    "La": {
        "num": 57,
        "name": "Lanthanum",
        "electronegativity": 1.1,
        "valence_electrons": ["5d1", "6s2"],
    },
    "Ce": {
        "num": 58,
        "name": "Cerium",
        "electronegativity": 1.12,
        "valence_electrons": ["4f1", "5d1", "6s2"],
    },
    "Pr": {
        "num": 59,
        "name": "Praseodymium",
        "electronegativity": 1.13,
        "valence_electrons": ["4f3", "6s2"],
    },
    "Nd": {
        "num": 60,
        "name": "Neodymium",
        "electronegativity": 1.14,
        "valence_electrons": ["4f4", "6s2"],
    },
    "Pm": {
        "num": 61,
        "name": "Promethium",
        "electronegativity": 1.13,
        "valence_electrons": ["4f5", "6s2"],
    },
    "Sm": {
        "num": 62,
        "name": "Samarium",
        "electronegativity": 1.17,
        "valence_electrons": ["4f6", "6s2"],
    },
    "Eu": {
        "num": 63,
        "name": "Europium",
        "electronegativity": 1.2,
        "valence_electrons": ["4f7", "6s2"],
    },
    "Gd": {
        "num": 64,
        "name": "Gadolinium",
        "electronegativity": 1.2,
        "valence_electrons": ["4f7", "5d1", "6s2"],
    },
    "Tb": {
        "num": 65,
        "name": "Terbium",
        "electronegativity": 1.22,
        "valence_electrons": ["4f9", "6s2"],
    },
    "Dy": {
        "num": 66,
        "name": "Dysprosium",
        "electronegativity": 1.23,
        "valence_electrons": ["4f10", "6s2"],
    },
    "Ho": {
        "num": 67,
        "name": "Holmium",
        "electronegativity": 1.24,
        "valence_electrons": ["4f11", "6s2"],
    },
    "Er": {
        "num": 68,
        "name": "Erbium",
        "electronegativity": 1.24,
        "valence_electrons": ["4f12", "6s2"],
    },
    "Tm": {
        "num": 69,
        "name": "Thulium",
        "electronegativity": 1.25,
        "valence_electrons": ["4f13", "6s2"],
    },
    "Yb": {
        "num": 70,
        "name": "Ytterbium",
        "electronegativity": 1.1,
        "valence_electrons": ["4f14", "6s2"],
    },
    "Lu": {
        "num": 71,
        "name": "Lutetium",
        "electronegativity": 1.27,
        "valence_electrons": ["5d1", "6s2"],
    },
    "Hf": {
        "num": 72,
        "name": "Hafnium",
        "electronegativity": 1.3,
        "valence_electrons": ["5d2", "6s2"],
    },
    "Ta": {
        "num": 73,
        "name": "Tantalum",
        "electronegativity": 1.5,
        "valence_electrons": ["5d3", "6s2"],
    },
    "W": {
        "num": 74,
        "name": "Tungsten",
        "electronegativity": 2.36,
        "valence_electrons": ["5d4", "6s2"],
    },
    "Re": {
        "num": 75,
        "name": "Rhenium",
        "electronegativity": 1.9,
        "valence_electrons": ["5d5", "6s2"],
    },
    "Os": {
        "num": 76,
        "name": "Osmium",
        "electronegativity": 2.2,
        "valence_electrons": ["5d6", "6s2"],
    },
    "Ir": {
        "num": 77,
        "name": "Iridium",
        "electronegativity": 2.2,
        "valence_electrons": ["5d7", "6s2"],
    },
    "Pt": {
        "num": 78,
        "name": "Platinum",
        "electronegativity": 2.28,
        "valence_electrons": ["5d9", "6s1"],
    },
    "Au": {
        "num": 79,
        "name": "Gold",
        "electronegativity": 2.54,
        "valence_electrons": ["5d10", "6s1"],
    },
    "Hg": {
        "num": 80,
        "name": "Mercury",
        "electronegativity": 2.0,
        "valence_electrons": ["5d10", "6s2"],
    },
    "Tl": {
        "num": 81,
        "name": "Thallium",
        "electronegativity": 1.62,
        "valence_electrons": ["6s2", "6p1"],
    },
    "Pb": {
        "num": 82,
        "name": "Lead",
        "electronegativity": 2.33,
        "valence_electrons": ["6s2", "6p2"],
    },
    "Bi": {
        "num": 83,
        "name": "Bismuth",
        "electronegativity": 2.02,
        "valence_electrons": ["6s2", "6p3"],
    },
    "Po": {
        "num": 84,
        "name": "Polonium",
        "electronegativity": 2.0,
        "valence_electrons": ["6s2", "6p4"],
    },
    "At": {
        "num": 85,
        "name": "Astatine",
        "electronegativity": 2.2,
        "valence_electrons": ["6s2", "6p5"],
    },
    "Rn": {
        "num": 86,
        "name": "Radon",
        "electronegativity": 0,
        "valence_electrons": ["6s2", "6p6"],
    },
    "Fr": {
        "num": 87,
        "name": "Francium",
        "electronegativity": 0.7,
        "valence_electrons": ["7s1"],
    },
    "Ra": {
        "num": 88,
        "name": "Radium",
        "electronegativity": 0.89,
        "valence_electrons": ["7s2"],
    },
    "Ac": {
        "num": 89,
        "name": "Actinium",
        "electronegativity": 1.1,
        "valence_electrons": ["6d1", "7s2"],
    },
    "Th": {
        "num": 90,
        "name": "Thorium",
        "electronegativity": 1.3,
        "valence_electrons": ["6d2", "7s2"],
    },
    "Pa": {
        "num": 91,
        "name": "Protactinium",
        "electronegativity": 1.5,
        "valence_electrons": ["5f2", "6d1", "7s2"],
    },
    "U": {
        "num": 92,
        "name": "Uranium",
        "electronegativity": 1.38,
        "valence_electrons": ["5f3", "6d1", "7s2"],
    },
    "Np": {
        "num": 93,
        "name": "Neptunium",
        "electronegativity": 1.36,
        "valence_electrons": ["5f4", "6d1", "7s2"],
    },
    "Pu": {
        "num": 94,
        "name": "Plutonium",
        "electronegativity": 1.28,
        "valence_electrons": ["5f6", "7s2"],
    },
    "Am": {
        "num": 95,
        "name": "Americium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f7", "7s2"],
    },
    "Cm": {
        "num": 96,
        "name": "Curium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f7", "6d1", "7s2"],
    },
    "Bk": {
        "num": 97,
        "name": "Berkelium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f9", "7s2"],
    },
    "Cf": {
        "num": 98,
        "name": "Californium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f10", "7s2"],
    },
    "Es": {
        "num": 99,
        "name": "Einsteinium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f11", "7s2"],
    },
    "Fm": {
        "num": 100,
        "name": "Fermium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f12", "7s2"],
    },
    "Md": {
        "num": 101,
        "name": "Mendelevium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f13", "7s2"],
    },
    "No": {
        "num": 102,
        "name": "Nobelium",
        "electronegativity": 1.3,
        "valence_electrons": ["5f14", "7s2"],
    },
    "Lr": {
        "num": 103,
        "name": "Lawrencium",
        "electronegativity": 0,
        "valence_electrons": ["7s2", "7p1"],
    },
    "Rf": {
        "num": 104,
        "name": "Rutherfordium",
        "electronegativity": 0,
        "valence_electrons": ["6d2", "7s2"],
    },
    "Db": {
        "num": 105,
        "name": "Dubnium",
        "electronegativity": 0,
        "valence_electrons": ["6d3", "7s2"],
    },
    "Sg": {
        "num": 106,
        "name": "Seaborgium",
        "electronegativity": 0,
        "valence_electrons": ["6d4", "7s2"],
    },
    "Bh": {
        "num": 107,
        "name": "Bohrium",
        "electronegativity": 0,
        "valence_electrons": ["6d5", "7s2"],
    },
    "Hs": {
        "num": 108,
        "name": "Hassium",
        "electronegativity": 0,
        "valence_electrons": ["6d6", "7s2"],
    },
    "Mt": {
        "num": 109,
        "name": "Meitnerium",
        "electronegativity": 0,
        "valence_electrons": ["6d7", "7s2"],
    },
    "Ds": {
        "num": 110,
        "name": "Darmstadtium",
        "electronegativity": 0,
        "valence_electrons": ["6d9", "7s1"],
    },
    "Rg": {
        "num": 111,
        "name": "Roentgenium",
        "electronegativity": 0,
        "valence_electrons": ["6d10", "7s1"],
    },
    "Cn": {
        "num": 112,
        "name": "Copernicium",
        "electronegativity": 0,
        "valence_electrons": ["6d10", "7s2"],
    },
    "Nh": {
        "num": 113,
        "name": "Nihonium",
        "electronegativity": 0,
        "valence_electrons": ["7s2", "7p1"],
    },
    "Fl": {
        "num": 114,
        "name": "Flerovium",
        "electronegativity": 0,
        "valence_electrons": ["7s2", "7p2"],
    },
    "Mc": {
        "num": 115,
        "name": "Moscovium",
        "electronegativity": 0,
        "valence_electrons": ["7s2", "7p3"],
    },
    "Lv": {
        "num": 116,
        "name": "Livermorium",
        "electronegativity": 0,
        "valence_electrons": ["7s2", "7p4"],
    },
    "Ts": {
        "num": 117,
        "name": "Tennessine",
        "electronegativity": 0,
        "valence_electrons": ["7s2", "7p5"],
    },
    "Og": {
        "num": 118,
        "name": "Oganesson",
        "electronegativity": 0,
        "valence_electrons": ["7s2", "7p6"],
    },
}

atomic_sym2num = {sym: d["num"] for sym, d in atomic_data.items()}
atomic_num2sym = {num: sym for sym, num in atomic_sym2num.items()}


def get_atomic_number(atomic_symbol) -> int:
    if atomic_symbol not in atomic_sym2num.keys():
        raise ValueError("Unknown atomic symbol '{}'.".format(atomic_symbol))
    return atomic_sym2num[atomic_symbol]


def get_electronegativity(atom: str | int) -> float:
    """Get the electronegativity of an atom.

    :param atom: The atom number (int) or the atom symbol (str) to get the
        electronegativity for.
    :returns: The electronegativity.
    """
    if isinstance(atom, int):
        atom = atomic_num2sym[atom]
    return atomic_data[atom]["electronegativity"]


def get_num_valence_electrons(atom: str | int) -> int:
    """Get the number of valence electrons of an atom.

    :param atom: The atom number (int) or the atom symbol (str) to get the
        valence electrons for.
    :returns: The number of valence electrons.
    """
    if isinstance(atom, int):
        atom = atomic_num2sym[atom]
    atom = to_non_aromatic_symbol(atom)
    valence_electron_conf = atomic_data[atom]["valence_electrons"]
    num = 0
    p = re.compile(r"(?P<n>[0-9])(?P<t>[spfd])(?P<e>[0-9]{1,2})")
    for e_conf in valence_electron_conf:
        m = p.match(e_conf)
        assert m is not None
        num += int(m["e"])
    return num


def find_atom_numbers(num_valence_electrons=None, electronegativity=None) -> list[int]:
    """Get all atom numbers that match certain criteria.

    :param num_valence_electrons: (optional) The number of valence electrons.
    :param electronegativity: (optional) The electronegativity.

    :returns: Returns a list of atom numbers that match the search criteria.
    """
    subset = atomic_data.values()
    if electronegativity is not None:
        _subset = []
        for entry in subset:
            if entry["electronegativity"] == electronegativity:
                _subset.append(entry)
        subset = _subset
    if num_valence_electrons is not None:
        _subset = []
        for entry in subset:
            if get_num_valence_electrons(entry["num"]) == num_valence_electrons:
                _subset.append(entry)
        subset = _subset
    return [entry["num"] for entry in subset]
