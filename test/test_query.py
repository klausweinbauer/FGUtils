import pytest
import rdkit.Chem.rdmolfiles as rdmolfiles

from fgutils.query import *
from fgutils.parse import parse
from fgutils.fgconfig import *
from fgutils.utils import mol_to_graph
from fgutils.mapping import map_pattern


def test_get_functional_groups_raw():
    mol = parse("C=O")
    idx_map, groups = get_functional_groups_raw(mol)
    assert 2 == len(groups)
    assert "carbonyl" == groups[0]
    assert "aldehyde" == groups[1]
    assert [0, 1] == idx_map[0]
    assert [0, 1] == idx_map[1]


def test_check_fg():
    fg = get_FG_by_name("carbonyl")
    mol = parse("CC(=O)O")
    fgs, indices = check_functional_group(mol, 2, fg)
    assert 2 == len(fgs)
    assert 2 == len(indices)
    assert "carbonyl" == fgs[0]
    assert "carboxylic_acid" == fgs[1]
    assert [1, 2] == indices[0]
    assert [1, 2, 3] == indices[1]


@pytest.mark.parametrize(
    "name,mol,anchor,exp_indices",
    [
        ("carbonyl", "CC(=O)O", 2, [1, 2]),
        ("carboxylic_acid", "CC(=O)O", 2, [1, 2, 3]),
        ("amide", "C(=O)N", 2, [0, 1, 2]),
    ],
)
def test_get_functional_group(name, mol, anchor, exp_indices):
    fg = get_FG_by_name(name)
    mol = parse(mol)
    is_fg, indices = get_functional_group(mol, anchor, fg)
    print(is_fg, indices)
    assert True == is_fg
    assert len(exp_indices) == len(indices)
    assert exp_indices == indices


def _test_fg(smiles, group_name, fg_anchor, fg_indices=None):
    group_names = get_FG_names()
    assert group_name in group_names, "Functional group {} is not implemented.".format(
        group_name
    )

    mol = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
    fg_indices = [i for i in mol.nodes()] if fg_indices is None else fg_indices
    idx_map, groups = get_functional_groups_raw(mol)
    print("Index map: {}  Groups: {}".format(idx_map, groups))
    fg_root_chain = [fg.name for fg in get_FG_root_chain(group_name)]
    _groups = {}
    for _id in idx_map[fg_anchor]:
        _group = groups[_id]
        _indices = []
        for k, v in idx_map.items():
            if _id in v:
                _indices.append(k)
        _groups[_group] = sorted(_indices)
    print(_groups)
    assert group_name in _groups.keys()
    assert len(_groups[group_name]) == len(fg_indices)
    assert fg_indices == _groups[group_name]
    for _group, _indices in _groups.items():
        if _group == group_name:
            continue
        assert (
            _group in fg_root_chain
        ), "Functional group {} is not in root chane {}.".format(_group, fg_root_chain)
        for _i in _indices:
            assert (
                _i in _groups[group_name]
            ), "Group {} must be a subset of group {}.".format(_group, group_name)


def test_amid():
    # Formamide
    _test_fg("C(=O)N", "amide", 2, [0, 1, 2])
    # Asparagine
    _test_fg("NC(=O)CC(N)C(=O)O", "amide", 2, [0, 1, 2])

    #    # def test_acyl(self):
    #    #    # Acetyl cloride
    #    #    self.__test_fg("CC(=O)[Cl]", "acyl", [1, 2])
    #
    #    def test_diol(self):
    #        self.__test_fg("OCO", "diol")
    #
    #    def test_hemiacetal(self):
    #        self.__test_fg("COCO", "hemiacetal")
    #
    #    def test_acetal(self):
    #        self.__test_fg("COCOC", "acetal")
    #
    #    def test_urea(self):
    #        self.__test_fg("O=C(N)O", "urea")
    #
    #    def test_carbonat(self):
    #        self.__test_fg("O=C(O)O", "carbonat")
    #        self.__test_fg("COC(=O)O", "carbonat", [1, 2, 3, 4])


def test_ester():
    # Methyl acetate
    _test_fg("COC(C)=O", "carboxylate_ester", 1, [1, 2, 4])

    #    def test_anhydrid(self):
    #        self.__test_fg("O=C(C)OC=O", "anhydrid")


def test_acid():
    # Acetic acid
    _test_fg("CC(=O)O", "carboxylic_acid", 2, [1, 2, 3])

    #    def test_anilin(self):
    #        self.__test_fg("Nc1ccccc1", "anilin")
    #
    #    def test_amin(self):
    #        # Glycin
    #        self.__test_fg("NCC(=O)O", "amin", [0, 1])
    #        # Methcatione
    #        self.__test_fg("CNC(C)C(=O)c1ccccc1", "amin", [0, 1, 2])
    #
    #    def test_nitril(self):
    #        self.__test_fg("C#N", "nitril")
    #
    #    def test_hydroxylamin(self):
    #        self.__test_fg("NO", "hydroxylamin")
    #
    #    def test_nitrose(self):
    #        self.__test_fg("N=O", "nitrose")
    #
    #    def test_nitro(self):
    #        self.__test_fg("ON=O", "nitro")
    #
    #    def test_thioether(self):
    #        # Diethylsulfid
    #        self.__test_fg("CCSCC", "thioether", [1, 2, 3])
    #
    #    def test_thioester(self):
    #        # Methyl thionobenzonat
    #        self.__test_fg("CSC(=O)c1ccccc1", "thioester", [1, 2, 3])
    #        self.__test_fg("COC(=S)c1ccccc1", "thioester", [1, 2, 3])


def test_keton():
    # Methcatione
    _test_fg("CNC(C)C(=O)c1ccccc1", "ketone", 5, [4, 5])


def test_aldehyde():
    # Formaldehyde
    _test_fg("C=O", "aldehyde", 1)
    _test_fg("CC=O", "aldehyde", 2, [1, 2])
