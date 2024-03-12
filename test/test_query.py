import pytest
import collections
import rdkit.Chem.rdmolfiles as rdmolfiles

from fgutils.query import *
from fgutils.parse import parse
from fgutils.fgconfig import *
from fgutils.utils import mol_to_graph


def _get_group_map(index_map):
    group_map = collections.defaultdict(lambda: [])
    for idx, group_indices in index_map.items():
        for group_idx in group_indices:
            assert idx not in group_map[group_idx]
            group_map[group_idx].append(idx)
    return dict(group_map)


def _assert_fg(raw_result, fg_name, indices):
    r_idxmap, r_groups = raw_result
    assert (
        indices[0] in r_idxmap.keys()
    ), "No functional group found for atom {}.".format(indices[0])
    groups = [r_groups[i] for i in r_idxmap[indices[0]]]
    assert (
        fg_name in groups
    ), "Could not find functional group '{}' for for index {}.".format(
        fg_name, indices[0]
    )
    group_i = -1
    for i, g in enumerate(groups):
        if g == fg_name:
            assert group_i == -1, "Found group multipel times ({}).".format(groups)
            group_i = i
    group_map = _get_group_map(r_idxmap)
    assert len(group_map[group_i]) == len(
        indices
    ), "Expected group '{}' to have {} atoms but found {}.".format(
        fg_name, len(indices), len(group_map[group_i])
    )
    for i in group_map[group_i]:
        assert (
            i in indices
        ), "Could not find index {} in functional group {} (Indices: {}).".format(
            i, fg_name, group_map[group_i]
        )


def _assert_not_fg(raw_result, fg_name, indices):
    r_idxmap, r_groups = raw_result
    for idx in indices:
        groups = [r_groups[i] for i in r_idxmap[idx]]
        assert (
            fg_name not in groups
        ), "Wrongly identified functional group '{}' for for index {}.".format(
            fg_name, indices[0]
        )


def test_get_functional_groups_raw():
    mol = parse("C=O")
    r = get_functional_groups_raw(mol)
    _assert_fg(r, "carbonyl", [0, 1])
    _assert_fg(r, "ketone", [0, 1])
    _assert_fg(r, "aldehyde", [0, 1])


@pytest.mark.parametrize(
    "name,smiles,anchor,exp_indices",
    [
        ("carbonyl", "CC(=O)O", 2, [1, 2]),
        ("carboxylic_acid", "CC(=O)O", 2, [1, 2, 3]),
        ("amide", "C(=O)N", 2, [0, 1, 2]),
    ],
)
def test_get_functional_group(name, smiles, anchor, exp_indices):
    fg = get_FG_by_name(name)
    mol = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
    is_fg, indices = is_functional_group(mol, anchor, fg)
    assert True == is_fg
    assert len(exp_indices) == len(indices)
    assert exp_indices == indices


@pytest.mark.parametrize(
    "smiles,functional_groups,exp_indices",
    [
        pytest.param("C=O", ["aldehyde"], [[0, 1]], id="Formaldehyde"),
        pytest.param("C(=O)N", ["amide"], [[0, 1, 2]], id="Formamide"),
        pytest.param("NC(=O)CC(N)C(=O)O", ["amide"], [[0, 1, 2]], id="Asparagine"),
        pytest.param("CC(=O)[Cl]", ["carbonyl"], [[1, 2]], id="Acetyl cloride"),
        pytest.param("COC(C)=O", ["ester"], [[1, 2, 4]], id="Methyl acetate"),
        pytest.param("CC(=O)O", ["carboxylic_acid"], [[1, 2, 3]], id="Acetic acid"),
        pytest.param("NCC(=O)O", ["amine"], [[0]], id="Glycin"),
        pytest.param(
            "CNC(C)C(=O)c1ccccc1",
            ["amine", "ketone"],
            [[1], [4, 5]],
            id="Methcatione",
        ),
        pytest.param("CCSCC", ["thioether"], [[2]], id="Diethylsulfid"),
        pytest.param(
            "CSC(=O)c1ccccc1", ["thioester"], [[1, 2, 3]], id="Methyl thionobenzonat"
        ),
        #pytest.param("", [""], [[]], id=""),
    ],
)
def test_functional_group_on_compound(smiles, functional_groups, exp_indices):
    assert len(functional_groups) == len(exp_indices)
    mol = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
    r = get_functional_groups_raw(mol)
    for fg, indices in zip(functional_groups, exp_indices):
        _assert_fg(r, fg, indices)
