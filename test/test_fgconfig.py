import pytest

from fgutils.fgconfig import *

def test_init():
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
    assert 2 == len(fgc.subgroups)
    sfgc1 = fgc.subgroups[0]
    sfgc2 = fgc.subgroups[1]
    assert "carbonyl" == fgc.name
    assert "aldehyde" == sfgc1.name
    assert "ketone" == sfgc2.name
    assert "carbonyl" == sfgc1.parent.name
    assert "carbonyl" == sfgc2.parent.name
    assert True == isinstance(fgc.group_atoms, list)
    assert True == isinstance(fgc.group_atoms[0], list)
    assert True == isinstance(sfgc1.group_atoms, list)
    assert True == isinstance(sfgc1.group_atoms[0], list)
    assert [[1 == 2]], sfgc1.group_atoms
    assert [[1 == 2]], sfgc2.group_atoms
