import numpy as np

from fgutils.const import BOND_KEY, SYMBOL_KEY, ATOM_COLORS
from fgutils.rdkit import get_mol_coords


def delay_import_cairo():
    try:
        import cairosvg
    except OSError as e:
        raise ImportError("Failed to load CairoSVG.") from e
    except ImportError as e:
        raise ImportError("Install CairoSVG to rasterize images.") from e
    return cairosvg


def delay_import_drawsvg():
    try:
        import drawsvg
    except OSError as e:
        raise ImportError("Failed to load Drawsvg.") from e
    except ImportError as e:
        raise ImportError("Install Drawsvg to create SVG images.") from e
    return drawsvg


def _rotmat(angle):
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


class BondLine:

    break_color = "#c26a77ff"
    form_color = "#5da899ff"
    default_color = "#000000"
    aromatic_dash = "8,4"
    bond_dist = 5
    stroke_width = 5

    def __init__(
        self, bond_type="default", length=1, loffset=0, noffset=0, is_aromatic=False
    ):
        if bond_type == "default":
            self.color = BondLine.default_color
        elif bond_type == "break":
            self.color = BondLine.break_color
        elif bond_type == "form":
            self.color = BondLine.form_color
        else:
            raise ValueError(
                "Argument bond_type can only be 'default', 'break' or 'form'."
            )
        self.length = length
        self.noffset = noffset
        self.loffset = loffset
        self.is_aromatic = is_aromatic

    def draw(self, p1, p2):
        draw = delay_import_drawsvg()
        dp_len = np.linalg.norm(p2 - p1)
        dp = (p2 - p1) / dp_len
        rmat = _rotmat(0.5 * np.pi)
        dpn = np.matmul(dp.T, rmat)
        start_p = p1 + BondLine.bond_dist * self.noffset * dpn + self.loffset * dp_len
        end_p = start_p + dp * dp_len * self.length
        pattern = None
        if self.is_aromatic:
            pattern = BondLine.aromatic_dash
        return draw.Line(
            *start_p,
            *end_p,
            stroke_width=BondLine.stroke_width,
            stroke=self.color,
            stroke_dasharray=pattern,
        )


bond_config = {
    "0,1": [BondLine("form")],
    "0,1.5": [BondLine("form", noffset=1), BondLine("form", noffset=-1, is_aromatic=True)],
    "0,2": [BondLine("form", noffset=1), BondLine("form", noffset=-1)],
    "0,3": [
        BondLine("form", noffset=2),
        BondLine("form"),
        BondLine("form", noffset=-2),
    ],
    "1,0": [BondLine("break")],
    "1,1": [BondLine()],
    "1,1.5": [BondLine(noffset=1), BondLine("form", noffset=-1, is_aromatic=True)],
    "1,2": [BondLine(noffset=1), BondLine("form", noffset=-1)],
    "1,3": [BondLine(noffset=2), BondLine("form"), BondLine("form", noffset=-2)],
    "1.5,1": [BondLine(noffset=1), BondLine("break", noffset=-1, is_aromatic=True)],
    "1.5,1.5": [BondLine(noffset=1), BondLine(noffset=-1, is_aromatic=True)],
    "1.5,2": [BondLine(noffset=1), BondLine("form", noffset=-1, is_aromatic=True)],
    "2,0": [BondLine("break", noffset=1), BondLine("break", noffset=-1)],
    "2,1": [BondLine(noffset=1), BondLine("break", noffset=-1)],
    "2,1.5": [BondLine(noffset=1), BondLine("break", noffset=-1, is_aromatic=True)],
    "2,2": [BondLine(noffset=1), BondLine(noffset=-1)],
    "2,3": [BondLine(noffset=2), BondLine(), BondLine("form", noffset=-2)],
    "3,0": [
        BondLine("break", noffset=2),
        BondLine("break"),
        BondLine("break", noffset=-2),
    ],
    "3,1": [BondLine(noffset=2), BondLine("break"), BondLine("break", noffset=-2)],
    "3,2": [BondLine(noffset=2), BondLine(), BondLine("break", noffset=-2)],
    "3,3": [BondLine(noffset=2), BondLine(), BondLine(noffset=-2)],
}


def graph_to_svg_drawing(g):
    draw = delay_import_drawsvg()

    border = 30
    positions = get_mol_coords(g, scale=66)

    p_array = np.array([[x, y] for x, y in positions.values()])
    p_min = np.min(p_array, axis=0)
    p_max = np.max(p_array, axis=0)

    size_x, size_y = p_max[0] - p_min[0] + 2 * border, p_max[1] - p_min[1] + 2 * border
    d = draw.Drawing(size_x, size_y, origin=(0, 0))
    # d.append(draw.Rectangle(-1000, -1000, 2000, 2000, fill="lightgray"))

    offset = -p_min + border
    for u, v, edata in g.edges(data=True):
        p1 = np.array(positions[u]) + offset
        p2 = np.array(positions[v]) + offset

        if isinstance(edata[BOND_KEY], float) or isinstance(edata[BOND_KEY], int):
            b_key = "{},{}".format(edata[BOND_KEY], edata[BOND_KEY])
        else:
            b_key = "{},{}".format(edata[BOND_KEY][0], edata[BOND_KEY][1])
        bonds = bond_config[b_key]
        for bond in bonds:
            d.append(bond.draw(p1, p2))
    for k, p in positions.items():
        p = np.array(p) + offset
        atom_symbol = g.nodes(data=True)[k][SYMBOL_KEY]
        d.append(draw.Circle(p[0], p[1], 24, fill="white"))
        d.append(
            draw.Text(
                atom_symbol,
                font_size=32,
                x=p[0],
                y=p[1],
                text_anchor="middle",
                dominant_baseline="middle",
                fill=ATOM_COLORS[atom_symbol],
            )
        )

    return d
