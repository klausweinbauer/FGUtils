import unittest
import rdkit.Chem.rdmolfiles as rdmolfiles
from fgutils.functional_group import *


class TestPatternParse(unittest.TestCase):
    def _assert_graph(self, g, exp_nodes, exp_edges):
        self.assertEqual(len(exp_nodes), g.number_of_nodes())
        self.assertEqual(len(exp_edges), g.number_of_edges())
        for i, sym in exp_nodes.items():
            self.assertEqual(sym, g.nodes[i]["symbol"])
        for i1, i2, order in exp_edges:
            self.assertEqual(order, g.edges[i1, i2]["bond"])

    def test_tokenize(self):
        def _ct(token, exp_type, exp_value, exp_col):
            return (
                token[0] == exp_type and token[1] == exp_value and token[2] == exp_col
            )

        it = tokenize("RC(=O)OR")
        self.assertTrue(_ct(next(it), "WILDCARD", "R", 0))
        self.assertTrue(_ct(next(it), "ATOM", "C", 1))
        self.assertTrue(_ct(next(it), "BRANCH_START", "(", 2))
        self.assertTrue(_ct(next(it), "BOND", "=", 3))
        self.assertTrue(_ct(next(it), "ATOM", "O", 4))
        self.assertTrue(_ct(next(it), "BRANCH_END", ")", 5))
        self.assertTrue(_ct(next(it), "ATOM", "O", 6))
        self.assertTrue(_ct(next(it), "WILDCARD", "R", 7))

    def test_branch(self):
        exp_nodes = {0: "R", 1: "C", 2: "O", 3: "O", 4: "R"}
        exp_edges = [(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)]
        g = parse("RC(=O)OR")
        self._assert_graph(g, exp_nodes, exp_edges)

    def test_multi_branch(self):
        exp_nodes = {0: "C", 1: "C", 2: "C", 3: "O", 4: "O", 5: "C"}
        exp_edges = [(0, 1, 1), (1, 2, 2), (2, 3, 1), (2, 4, 1), (1, 5, 1)]
        g = parse("CC(=C(O)O)C")
        self._assert_graph(g, exp_nodes, exp_edges)

    def test_ring_3(self):
        exp_nodes = {0: "C", 1: "C", 2: "C"}
        exp_edges = [(0, 1, 1), (1, 2, 1), (0, 2, 1)]
        g = parse("C1CC1")
        self._assert_graph(g, exp_nodes, exp_edges)

    def test_ring_4(self):
        exp_nodes = {0: "C", 1: "C", 2: "C", 3: "C"}
        exp_edges = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (0, 3, 1)]
        g = parse("C1CCC1")
        self._assert_graph(g, exp_nodes, exp_edges)

    def test_multi_ring(self):
        exp_nodes = {0: "C", 1: "C", 2: "C", 3: "C"}
        exp_edges = [(0, 1, 1), (1, 2, 1), (0, 2, 1), (1, 3, 1), (2, 3, 1)]
        g = parse("C1C2C1C2")
        self._assert_graph(g, exp_nodes, exp_edges)

    def test_aromatic_ring(self):
        exp_nodes = {i: "c" for i in range(6)}
        exp_edges = [(0, 5, 1.5), *[(i, i + 1, 1.5) for i in range(5)]]
        g = parse("c1ccccc1")
        self._assert_graph(g, exp_nodes, exp_edges)

    def test_aromatic_ring_syntax_error(self):
        with self.assertRaises(SyntaxError):
            parse("c1ccccC1")

    def test_complex_aromatic_ring(self):
        exp_nodes = {i: "c" for i in range(9)}
        exp_nodes[0] = "C"
        exp_nodes[3] = "C"
        exp_nodes[5] = "C"
        exp_edges = [
            (0, 1, 1),
            (1, 2, 1.5),
            (2, 3, 1),
            (2, 4, 1.5),
            (4, 5, 2),
            (4, 6, 1.5),
            (6, 7, 1.5),
            (7, 8, 1.5),
            (8, 1, 1.5),
        ]
        g = parse("Cc1c(C)c(=C)ccc1")
        self._assert_graph(g, exp_nodes, exp_edges)


class TestMappingPermutations(unittest.TestCase):
    def test_1(self):
        self.assertEqual([[]], get_mapping_permutations([], ["C"]))

    def test_2(self):
        self.assertEqual([[(0, 1)]], get_mapping_permutations(["C"], ["S", "C"]))

    def test_3(self):
        self.assertEqual(
            [[(0, 0)], [(0, 1)]], get_mapping_permutations(["C"], ["C", "C"])
        )

    def test_4(self):
        self.assertEqual(
            [[(0, 0), (1, 1)], [(0, 1), (1, 0)]],
            get_mapping_permutations(["C", "C"], ["C", "C"]),
        )

    def test_5(self):
        self.assertEqual([], get_mapping_permutations(["C", "O"], ["C"]))

    def test_wildcard(self):
        self.assertEqual(
            [[(0, 0)]], get_mapping_permutations(["R"], ["C"], wildcard="R")
        )

    def test_wildcard_to_nothing(self):
        self.assertEqual([[(0, -1)]], get_mapping_permutations(["R"], [], wildcard="R"))

    def test_wildcard_to_nothing_double(self):
        self.assertEqual(
            [[(0, -1), (1, -1)]], get_mapping_permutations(["R", "R"], [], wildcard="R")
        )

    def test_wildcard_match_and_match_to_nothing(self):
        self.assertEqual(
            [[(0, 0), (1, -1)], [(0, -1), (1, 0)]],
            get_mapping_permutations(["R", "R"], ["C"], wildcard="R"),
        )

    def test_ignore_case(self):
        self.assertEqual([[(0, 0)]], get_mapping_permutations(["c"], ["C"]))


class TestAnchoredPatternMatch(unittest.TestCase):
    def _assert_mapping(self, mapping, valid, exp_mapping=[]):
        self.assertEqual(mapping[0], valid)
        for emap in exp_mapping:
            self.assertIn(emap, mapping[1])

    def test_simple_match(self):
        exp_mapping = [(1, 0), (2, 1)]
        g = parse("CCO")
        p = parse("RO")
        m = anchored_pattern_match(g, 2, p, 1)
        self._assert_mapping(m, True, exp_mapping)

    def test_branched_match(self):
        exp_mapping = [(0, 0), (1, 1), (2, 2), (3, 3)]
        g = parse("CC(=O)O")
        p = parse("RC(=O)O")
        m = anchored_pattern_match(g, 2, p, 2)
        self._assert_mapping(m, True, exp_mapping)

    def test_ring_match(self):
        exp_mapping = [(0, 2), (1, 1), (2, 0)]
        g = parse("C1CO1")
        p = parse("R1CC1")
        m = anchored_pattern_match(g, 1, p, 1)
        self._assert_mapping(m, True, exp_mapping)

    def test_not_match(self):
        g = parse("CC=O")
        p = parse("RC(=O)NR")
        m = anchored_pattern_match(g, 2, p, 2)
        self._assert_mapping(m, False)

    def test_1(self):
        g = parse("CC=O")
        p = parse("RC(=O)R")
        m = anchored_pattern_match(g, 0, p, 3)
        self._assert_mapping(m, False)

    def test_2(self):
        exp_mapping = [(0, 0), (1, 1), (2, 2)]
        g = parse("CC=O")
        p = parse("RC=O")
        m = anchored_pattern_match(g, 2, p, 2)
        self._assert_mapping(m, True, exp_mapping)

    def test_ignore_aromaticity(self):
        exp_mapping = [(1, 0), (2, 1)]
        g = parse("c1c(=O)cccc1")
        p = parse("C=O")
        m = anchored_pattern_match(g, 2, p, 1)
        self._assert_mapping(m, True, exp_mapping)

    def test_3(self):
        exp_mapping = [(0, 4), (1, 3), (2, 1), (4, 2), (3, 0)]
        g = parse("COC(C)=O")
        p = parse("RC(=O)OR")
        m = anchored_pattern_match(g, 4, p, 2)
        self._assert_mapping(m, True, exp_mapping)

    def test_explore_wrong_branch(self):
        exp_mapping = [(0, 2), (1, 1), (2, 0), (3, 3)]
        g = parse("COCO")
        p = parse("C(OR)O")
        m = anchored_pattern_match(g, 1, p, 1)
        self._assert_mapping(m, True, exp_mapping)

    def test_match_pattern_to_mol(self):
        exp_mapping = [(0, 2), (1, 0), (2, 1)]
        g = parse("NC(=O)C")
        p = parse("C(=O)N")
        m = anchored_pattern_match(g, 2, p, 1)
        self._assert_mapping(m, True, exp_mapping)

    def test_match_hydrogen(self):
        exp_mapping = [(0, 1), (1, 2)]
        g = parse("C=O")
        p = parse("RC=O")
        m = anchored_pattern_match(g, 1, p, 2)
        self._assert_mapping(m, True, exp_mapping)

    def test_match_hydrogen2(self):
        exp_mapping = [(0, 0), (1, 1)]
        g = parse("C=O")
        p = parse("C(=O)R")
        m = anchored_pattern_match(g, 1, p, 1)
        self._assert_mapping(m, True, exp_mapping)

    def test_match_hydrogen3(self):
        exp_mapping = [(0, 1), (1, 2)]
        g = parse("C=O")
        p = parse("RC=O")
        m = anchored_pattern_match(g, 0, p, 1, verbose=True)
        self._assert_mapping(m, True, exp_mapping)

class TestFGTreeGeneration(unittest.TestCase):
    def test1(self):
        config = {
            "name": "carbonyl",
            "pattern": "C(=O)",
            "subgroups": [
                {
                    "name": "aldehyde",
                    "pattern": "RC(=O)",
                    "group_atoms": [1, 2],
                    "anti_pattern": "RC(=O)R",
                },
                {"name": "ketone", "pattern": "RC(=O)R", "group_atoms": [1, 2]},
            ],
        }

        fgc = FGConfig(**config)
        self.assertEqual(2, len(fgc.subgroups))
        sfgc1 = fgc.subgroups[0]
        sfgc2 = fgc.subgroups[1]
        self.assertEqual("carbonyl", fgc.name)
        self.assertEqual("aldehyde", sfgc1.name)
        self.assertEqual("ketone", sfgc2.name)
        self.assertEqual("carbonyl", sfgc1.parent.name)
        self.assertEqual("carbonyl", sfgc2.parent.name)
        self.assertTrue(isinstance(fgc.group_atoms, list))
        self.assertTrue(isinstance(fgc.group_atoms[0], list))
        self.assertTrue(isinstance(sfgc1.group_atoms, list))
        self.assertTrue(isinstance(sfgc1.group_atoms[0], list))
        self.assertEqual([[1, 2]], sfgc1.group_atoms)
        self.assertEqual([[1, 2]], sfgc2.group_atoms)


class TestFunctionalGroupCheck(unittest.TestCase):
    def __test_fg(self, smiles, group_name, indices=None):
        groups = get_FGs_flat()
        group_names = get_FG_names()
        assert (
            group_name in group_names
        ), "Functional group {} is not implemented.".format(group_name)

        mgraph = mol_to_graph(rdmolfiles.MolFromSmiles(smiles))
        # print(
        #    "Graph: {}".format(
        #        " ".join(
        #            [
        #                "{}[{}]".format(n[1]["symbol"], n[0])
        #                for n in mgraph.nodes(data=True)
        #            ]
        #        )
        #    )
        # )
        indices = [i for i in mgraph.nodes()] if indices is None else indices
        for idx in mgraph.nodes():
            pos_groups = []
            for g in groups:
                result = check_functional_group(mgraph, g, idx, verbose=False)  # type: ignore
                if result:
                    pos_groups.append(g.name)
            if idx in indices:
                valid_result = len(pos_groups) > 0 and group_name in pos_groups
                root_chain = get_FG_root_chain(group_name)
                root_chain_names = [g.name for g in root_chain]
                for pg in pos_groups:
                    if pg not in root_chain_names:
                        valid_result = False
                self.assertTrue(
                    valid_result,
                    (
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
                    ),
                )
            else:
                self.assertTrue(
                    group_name not in pos_groups,
                    (
                        "Atom {} ({}) in {} was expected to be not part of "
                        + "functional group {} but was identified as {}."
                    ).format(
                        mgraph.nodes[idx]["symbol"], idx, smiles, group_name, pos_groups
                    ),
                )

            #    def test_phenol(self):
            #        # 2,3 Xylenol
            #        self.__test_fg("Cc1cccc(O)c1C", "phenol", [1, 2, 3, 4, 5, 6, 7])
            #        self.__test_fg("Oc1ccc[nH]1", "phenol")
            #
            #    def test_alcohol(self):
            #        # Ethanol
            #        self.__test_fg("CCO", "alcohol", [1, 2])
            #
            #    def test_ether(self):
            #        # Methylether
            #        self.__test_fg("COC", "ether", [0, 1, 2])
            #        self.__test_fg("COc1cccc(OC)c1OC", "ether", [0, 1, 2, 9, 10, 11, 6, 7, 8])
            #        self.__test_fg("COc1ccc[nH]1", "ether", [0, 1, 2])
            #        self.__test_fg("COc1ccccc1", "ether", [0, 1, 2])
            #
            #    def test_enol(self):
            #        # 3-pentanone enol
            #        self.__test_fg("CCC(O)=CC", "enol", [2, 3, 4])

    def test_amid(self):
        # Formamide
        self.__test_fg("C(=O)N", "amide", [0, 1, 2])
        # Asparagine
        self.__test_fg("NC(=O)CC(N)C(=O)O", "amide", [0, 1, 2])

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

    def test_ester(self):
        # Methyl acetate
        self.__test_fg("COC(C)=O", "carboxylate_ester", [1, 2, 4])


        #    def test_anhydrid(self):
        #        self.__test_fg("O=C(C)OC=O", "anhydrid")

    def test_acid(self):
        # Acetic acid
        self.__test_fg("CC(=O)O", "carboxylic_acid", [1, 2, 3])

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

    def test_keton(self):
        # Methcatione
        self.__test_fg("CNC(C)C(=O)c1ccccc1", "ketone", [4, 5])

    def test_aldehyde(self):
        # Formaldehyde
        self.__test_fg("C=O", "aldehyde")
        # Methcatione
        # self.__test_fg("CC=O", "aldehyde", [1, 2])

    def test_xxx(self):
        g = parse("C=O")
        c = get_FG_by_name("aldehyde")
        r = check_functional_group(g, c, 0, verbose=True)
        print(r)
