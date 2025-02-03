import copy
import itertools
import collections
import numpy as np


def generate_mapping_permutations(pattern, structure, wildcard=None):
    mappings = []
    struct_map = [(i, s) for i, s in enumerate(structure)]
    for struct_permut in itertools.permutations(struct_map):
        mapping = []
        is_match = len(pattern) > 0
        for i, pattern_sym in enumerate(pattern):
            if len(struct_permut) > i and (
                pattern_sym == wildcard or pattern_sym == struct_permut[i][1]
            ):
                mapping.append((i, struct_permut[i][0]))
            else:
                is_match = False
                break
        if is_match:
            mappings.append(mapping)
    return mappings


class PermutationMapper:
    """The Permutation Mapper class specifies how nodes can match. The wildcard
    is specified for the pattern characters and not the structure characters.
    This means that the wildcard will only match in the direction Pattern ->
    Structure. The ``R`` as wildcard and ``C`` as normal label for example.
    ``R -> C`` will match but ``C -> R`` will not.

    :param wildcard: (optional) The wildcard is a specific pattern node label
        that can match to any other node. In the chemical notation this would
        be the R group, e.g., RC(=O)OH for an acid.
    :param ignore_case: (optional) Flag if the matching of node labels is case
        insensitive. Default: False
    :param can_map_to_nothing: (optional) A list of characters that can map to
        nothing. A node with a label in this list does not need a mapping at
        all. For example optional nodes that can or can not have a mapping.
        This is different to the wildcard in the sense that a wildcard must map
        to node (with an arbitrary label though).
    """

    def __init__(self, wildcard=None, ignore_case=False, can_map_to_nothing=[]):
        self.wildcard = wildcard
        self.ignore_case = ignore_case
        self.can_map_to_nothing = sorted(
            (
                can_map_to_nothing
                if isinstance(can_map_to_nothing, list)
                else [can_map_to_nothing]
            ),
            key=lambda x: 1 if wildcard is not None and x in wildcard else 0,
        )

    def permute(self, pattern, structure):
        wildcard = self.wildcard
        can_map_to_nothing = self.can_map_to_nothing
        if self.ignore_case:
            wildcard = None if wildcard is None else wildcard.lower()
            pattern = [p.lower() for p in pattern]
            structure = [s.lower() for s in structure]
            can_map_to_nothing = [cmtn.lower() for cmtn in can_map_to_nothing]

        struct_additions = []
        if len(can_map_to_nothing) > 0:
            structure = copy.deepcopy(structure)
            for cmtn in can_map_to_nothing:
                if cmtn == wildcard:
                    num_to_add = len(pattern) - len(structure)
                else:
                    pattern_ref = [(i, p) for i, p in enumerate(pattern) if p == cmtn]
                    struct_ref = [(i, s) for i, s in enumerate(structure) if s == cmtn]
                    num_to_add = len(pattern_ref) - len(struct_ref)
                for i in range(len(structure), len(structure) + num_to_add):
                    structure.append(cmtn)
                    struct_additions.append(i)

        mappings = generate_mapping_permutations(pattern, structure, wildcard=wildcard)

        if len(struct_additions) > 0:
            for mapping in mappings:
                for i, (pi, si) in enumerate(mapping):
                    if si in struct_additions:
                        mapping[i] = (pi, -1)

        unique_mappings = []
        mapping_sets = []
        for mapping in mappings:
            mapping_set = set(mapping)
            if mapping_set not in mapping_sets:
                unique_mappings.append(mapping)
                mapping_sets.append(mapping_set)

        return unique_mappings


class MappingMatrix:
    """A class for fast label mapping evaluation according to a
    :py:class:`~fgutils.permutation.PermutationMapper`.

    :param pattern_symbols: A complete list of all possible pattern labels.
    :param structure_symbols: A complete list of all possible symbol labels.
    :param mapper: A PermutationMapper instance. The matrix is constructed
        based on these mappings.
    """

    def __init__(
        self,
        pattern_symbols: list[str],
        structure_symbols: list[str],
        mapper: PermutationMapper,
    ):
        self.__s2i = {
            sym: idx
            for idx, sym in enumerate(set([*pattern_symbols, *structure_symbols]))
        }
        self.__i2s = {idx: sym for sym, idx in self.__s2i.items()}
        self.__valid_mappings = np.zeros((len(self.__s2i), len(self.__s2i)))
        for ps in pattern_symbols:
            for ss in structure_symbols:
                mappings = mapper.permute(ps, ss)
                assert len(mappings) <= 1
                if len(mappings) > 0:
                    assert len(mappings[0]) == 1
                    self.__valid_mappings[
                        self.__s2i[ps],
                        self.__s2i[ss],
                    ] = 1

    def is_mapping(self, pattern_symbol: str, structure_symbol: str) -> bool:
        """Check if two symbols can match, i.e., two nodes with respective
        label can be mapped.

        :param pattern_symbol: The pattern symbol to check if it matches the
            structure symbol.
        :param structure_symbol: The structure symbol to check if it matches
            the pattern symbol.

        :returns: True if the symbols can be mapped and False otherwise.
        """
        i = self.__s2i[pattern_symbol]
        j = self.__s2i[structure_symbol]
        return self.__valid_mappings[i, j] == 1

    def min_mapping_symbol(
        self,
        pattern_symbols: list[str],
        structure_symbols: list[str],
    ) -> tuple[str, str] | None:
        """This function finds the symbol with the lowest permutation count.
        Take ``[A, A, B]`` and ``[A, A, A, B, B]`` as pattern and symbol lists
        for example. The function will return ``B``, because there are only 2
        possible mappings whereas for ``A`` there are 6. This is useful to find
        suitable anchor nodes in pattern matching. It is 3x as fast to anchor
        ``B`` and check the alighment then anchoring ``A``. The min mapping
        excludes zero mappings. If there is no mapping possible at all None
        is returned.

        :param pattern_symbols: The list of pattern symbols.
        :param structure_symbols: The list of structure symbols.

        :returns: The tuple of symbols with the lowest permutation count. The
            first entry is the pattern symbol and the second entry is the
            structure symbol. If no mapping is found at all None is returned.
        """
        if len(pattern_symbols) > len(structure_symbols):
            raise ValueError("Pattern has more symbols than structure.")

        def _setup_vec(symbols):
            counter = collections.Counter(symbols)
            vec = np.zeros((len(self.__s2i), 1))
            for s, n in counter.items():
                vec[self.__s2i[s]] = n
            return vec

        ps_vec = _setup_vec(pattern_symbols)
        ss_vec = _setup_vec(structure_symbols)
        m_cnt = np.matmul(
            np.matmul(self.__valid_mappings, ss_vec),
            np.matmul(ps_vec.T, self.__valid_mappings),
        ) * self.__valid_mappings
        non_zero = m_cnt[np.nonzero(m_cnt)]
        if len(non_zero) == 0:
            return None
        index = np.where(m_cnt == np.min(non_zero))
        return (self.__i2s[index[0][0]], self.__i2s[index[1][0]])
