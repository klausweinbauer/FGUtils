from fgutils.parse import parse
from fgutils.synthesis import ReactionRule, apply_rule
from fgutils.utils import add_implicit_hydrogens
from fgutils.rdkit import mol_smiles_to_graph


def test_apply_rule_form_break():
    reactant_smiles = "[C:1][C:2](=[O:3])[O:4].[N:5]"
    reactant = mol_smiles_to_graph(reactant_smiles)
    exp_reaction = (
        "[CH3:1][C:2](=[O:3])[OH:4].[NH3:5]>>" + "[CH3:1][C:2](=[O:3])[NH2:5].[OH2:4]"
    )
    rule = "C(<0,1>N)<1,0>O"
    its_graphs = apply_rule(reactant, ReactionRule(parse(rule)), unique=False)
    assert len(its_graphs) == 1
    assert its_graphs[0].to_smiles() == exp_reaction


def test_apply_rule_form():
    reactant_smiles = "C=C.C"
    reactant = mol_smiles_to_graph(reactant_smiles, implicit_h=True)
    rule = "C1<2,1>C<0,1>C<1,0>H<0,1>1"
    its_graphs = apply_rule(reactant, ReactionRule(parse(rule)), unique=True)
    assert len(its_graphs) == 1
    assert (
        its_graphs[0].to_smiles(ignore_aam=True, implicit_h=True)
        == "[CH2]=[CH2].[CH4]>>[CH3][CH2][CH3]"
    )


def test_apply_rule_break():
    reactant = "CCC"
    rule = "C<1,2>C<1,0>C"
    its_graphs = apply_rule(parse(reactant), ReactionRule(parse(rule)), unique=False)
    assert len(its_graphs) == 2
    assert its_graphs[0].to_smiles(ignore_aam=True) == "CCC>>C.C=C"
    assert its_graphs[1].to_smiles(ignore_aam=True) == "CCC>>C.C=C"


def test_apply_rule_unique_argument():
    reactant = "CC(=O)O.N"
    reactant_g = add_implicit_hydrogens(parse(reactant))
    rule = "C1<0,1>N<1,0>H<0,1>O<1,0>1"
    its_graphs = apply_rule(reactant_g, ReactionRule(parse(rule)), unique=True)
    assert len(its_graphs) == 1


def test_apply_rule_without_disconnected():
    reactant = "NC(=O)O.N"
    reactant_g = add_implicit_hydrogens(parse(reactant))
    rule = "C1<0,1>N<1,0>H<0,1>O<1,0>1"
    its_graphs = apply_rule(
        reactant_g, ReactionRule(parse(rule)), unique=True, connected_only=True
    )
    assert len(its_graphs) == 1
