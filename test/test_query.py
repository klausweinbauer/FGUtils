import pytest
import rdkit.Chem.rdmolfiles as rdmolfiles

from fgutils.permutation import PermutationMapper
from fgutils.fgconfig import FGConfigProvider
from fgutils.query import FGQuery, is_functional_group
from fgutils.parse import parse
from fgutils.rdkit import mol_to_graph

default_mapper = PermutationMapper(wildcard="R", ignore_case=True)
default_config_provider = FGConfigProvider(mapper=default_mapper)
default_query = FGQuery(mapper=default_mapper, config_provider=default_config_provider)


@pytest.mark.parametrize(
    "name,smiles,anchor,exp_indices",
    [
        ("carbonyl", "CC(=O)O", 2, [1, 2]),
        ("carboxylic_acid", "CC(=O)O", 2, [1, 2, 3]),
        ("amide", "C(=O)N", 2, [0, 1, 2]),
        ("acyl_chloride", "CC(=O)[Cl]", 3, [1, 2, 3]),
    ],
)
def test_get_functional_group(name, smiles, anchor, exp_indices):
    fg = default_config_provider.get_by_name(name)
    mol = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
    is_fg, indices = is_functional_group(mol, anchor, fg, mapper=default_mapper)
    assert is_fg
    assert len(exp_indices) == len(indices)
    assert exp_indices == indices


def test_get_functional_groups():
    mol = parse("C=O")
    groups = default_query.get(mol)
    assert ("aldehyde", [0, 1]) in groups


def test_get_functional_group_once():
    mol = parse("CC(=O)OC")
    groups = default_query.get(mol)
    assert 1 == len(groups)
    assert ("ester", [1, 2, 3]) in groups


@pytest.mark.parametrize(
    "smiles,functional_groups,exp_indices",
    [
        pytest.param("C=O", ["aldehyde"], [[0, 1]], id="Formaldehyde"),
        pytest.param("C(=O)N", ["amide"], [[0, 1, 2]], id="Formamide"),
        pytest.param("NC(=O)CC(N)C(=O)O", ["amide"], [[0, 1, 2]], id="Asparagine"),
        pytest.param("[Cl]C(=O)C", ["acyl_chloride"], [[0, 1, 2]], id="Acetyl cloride"),
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
        pytest.param(
            "O=C(C)Oc1ccccc1C(=O)O",
            ["ester", "carboxylic_acid"],
            [[0, 1, 3], [10, 11, 12]],
            id="Acetylsalicylic acid",
        ),
        # pytest.param("", [""], [[]], id=""),
    ],
)
def test_functional_group_on_compound(smiles, functional_groups, exp_indices):
    assert len(functional_groups) == len(exp_indices)
    mol = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
    groups = default_query.get(mol)
    for fg, indices in zip(functional_groups, exp_indices):
        assert (fg, indices) in groups
