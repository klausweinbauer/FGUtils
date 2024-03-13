import copy
import itertools


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
    def __init__(self, wildcard=None, ignore_case=False, can_map_to_nothing=[]):
        self.wildcard = wildcard
        self.ignore_case = ignore_case
        self.can_map_to_nothing = sorted(
            can_map_to_nothing
            if isinstance(can_map_to_nothing, list)
            else [can_map_to_nothing],
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
