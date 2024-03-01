import re
import copy
import itertools
import numpy as np
import networkx as nx

from fgutils.parse import parse
from fgutils.permutation import Mapper


def check_functional_group(graph, config: FGConfig, index: int, verbose=False) -> bool:
    is_func_group = False

    for pattern, group_indices in zip(config.pattern, config.group_atoms):
        is_match, match = pattern_match(graph, index, pattern, verbose=verbose)
        if is_match:
            m_idx = [m[1] for m in match if m[0] == index][0]
            is_func_group = is_func_group or m_idx in group_indices

    last_len = config.max_pattern_size
    for apattern, apattern_size in sorted(
        [(m, m.number_of_nodes()) for m in config.anti_pattern],
        key=lambda x: x[1],
        reverse=True,
    ):
        if not is_func_group:
            break
        if last_len > apattern_size:
            last_len = apattern_size
        is_match, _ = pattern_match(graph, index, apattern, verbose=verbose)
        is_func_group = is_func_group and not is_match
    return is_func_group


