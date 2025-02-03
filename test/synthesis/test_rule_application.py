from fgutils.parse import parse
from fgutils.synthesis import ReactionRule, apply_rule
from fgutils.utils import add_implicit_hydrogens


def test_apply_rule_form_break():
    reactant = "CC(=O)O.N"
    rule = "C(<0,1>N)<1,0>O"
    its_graphs = apply_rule(parse(reactant), ReactionRule(parse(rule)), unique=False)
    assert len(its_graphs) == 1
    assert its_graphs[0].to_smiles(ignore_aam=True) == "CC(=O)O.N>>CC(N)=O.O"


def test_apply_rule_reduce_form():
    reactant = "C=C.C"
    rule = "C<2,1>C<0,1>C"
    its_graphs = apply_rule(parse(reactant), ReactionRule(parse(rule)), unique=False)
    assert len(its_graphs) == 2
    assert its_graphs[0].to_smiles(ignore_aam=True) == "C.C=C>>CCC"
    assert its_graphs[1].to_smiles(ignore_aam=True) == "C.C=C>>CCC"


def test_apply_rule_increase_brea():
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
