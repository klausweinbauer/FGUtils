import rdkit.Chem.rdmolfiles as rdmolfiles

from fgutils.query import *
from fgutils.fgconfig import get_FG_list, get_FG_names, get_FG_root_chain
from fgutils.utils import mol_to_graph


def _test_fg(smiles, group_name, indices=None):
    groups = get_FG_list()
    group_names = get_FG_names()
    assert group_name in group_names, "Functional group {} is not implemented.".format(
        group_name
    )

    mgraph = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
    indices = [i for i in mgraph.nodes()] if indices is None else indices
    for idx in mgraph.nodes():
        pos_groups = []
        for g in groups:
            result = is_functional_group(mgraph, g, idx)  # type: ignore
            if result:
                pos_groups.append(g.name)
        if idx in indices:
            valid_result = len(pos_groups) > 0 and group_name in pos_groups
            root_chain = get_FG_root_chain(group_name)
            root_chain_names = [g.name for g in root_chain]
            for pg in pos_groups:
                if pg not in root_chain_names:
                    valid_result = False
            assert True == valid_result, (
                "Atom {} ({}) in {} was expected to be in functional "
                + "group {} but was identified as {}. "
                + "Acceptable groups according to root chain are {}."
            ).format(
                mgraph.nodes[idx]["symbol"],
                idx,
                smiles,
                group_name,
                pos_groups,
                root_chain_names,
            )
        else:
            assert True == (group_name not in pos_groups), (
                "Atom {} ({}) in {} was expected to be not part of "
                + "functional group {} but was identified as {}."
            ).format(mgraph.nodes[idx]["symbol"], idx, smiles, group_name, pos_groups)


def test_amid():
    # Formamide
    _test_fg("C(=O)N", "amide", [0, 1, 2])
    # Asparagine
    _test_fg("NC(=O)CC(N)C(=O)O", "amide", [0, 1, 2])

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
    _test_fg("COC(C)=O", "carboxylate_ester", [1, 2, 4])

    #    def test_anhydrid(self):
    #        self.__test_fg("O=C(C)OC=O", "anhydrid")


def test_acid():
    # Acetic acid
    _test_fg("CC(=O)O", "carboxylic_acid", [1, 2, 3])

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
    _test_fg("CNC(C)C(=O)c1ccccc1", "ketone", [4, 5])


def test_aldehyde():
    # Formaldehyde
    _test_fg("C=O", "aldehyde")
    # Methcatione
    # self.__test_fg("CC=O", "aldehyde", [1, 2])


def test_xxx():
    g = parse("C=O")
    c = get_FG_by_name("aldehyde")
    r = check_functional_group(g, c, 0, verbose=True)
    print(r)
