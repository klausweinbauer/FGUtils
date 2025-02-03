from fgutils.chem.ps import atomic_sym2num


def test_atom2symbol_dict():
    tmpk = []
    tmpv = []
    for k, v in atomic_sym2num.items():
        assert k not in tmpk, "Found symbol '{}' twice.".format(k)
        assert v not in tmpv, "Found number '{}' twice.".format(v)
        tmpk.append(k)
        tmpv.append(v)
