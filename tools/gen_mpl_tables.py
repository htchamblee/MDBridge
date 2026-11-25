"""Write Matplotlib's colour data as C++ tables: named colours, the default colour cycle, colormaps.

    python3 -I -B tools/gen_mpl_tables.py

Imports Matplotlib 3.11.2 (the version of tests/parity/requirements-regenerate.txt; another version
is refused) and writes src/render/core/mpl_colors.inc, included by src/render/core/color.cpp and
colormap.cpp (namespace mdbridge::render::detail):

* kNamedColors: every name matplotlib.colors.to_rgba knows (its full map after the "grey" / "gray"
  aliases and the priority order of matplotlib/colors.py: base colours over Tableau over CSS4 over
  xkcd), sorted by UTF-8 bytes, with the RGBA that to_rgba gives for the name (alpha 1);
* kDefaultCycle: the colours of rcParamsDefault["axes.prop_cycle"] ("C0" ... "C9");
* kListedMaps / kListedMapColors: the colormaps MDBridge bundles (owner decision D6 of
  docs/specs/m4-plan.md): the qualitative ColorBrewer / Tableau maps and the perceptually uniform
  viridis family, each a ListedColormap of RGB rows, in this order.

Numbers are written as C++ hexadecimal floating literals (float.hex()), so the build reads exactly
the doubles Matplotlib computes. The output is deterministic; the build needs no Python.

The data is Matplotlib's (Matplotlib licence, PSF-style; ColorBrewer schemes under the Apache 2.0
licence; xkcd colour survey CC0; CSS4 names from the W3C; Tableau 10 / 20 as distributed with
Matplotlib); the notices are in the generated file, THIRD_PARTY_NOTICES.md and
packaging/licenses/Matplotlib-LICENSE.txt, ColorBrewer-LICENSE.txt.
"""

from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INC = os.path.join(ROOT, "src", "render", "core", "mpl_colors.inc")
MATPLOTLIB_VERSION = "3.11.2"
WIDTH = 100

# owner decision D6: Matplotlib's qualitative and perceptually uniform maps (and their "_r" reverses,
# made at run time)
QUALITATIVE = ["tab10", "tab20", "tab20b", "tab20c", "Set1", "Set2", "Set3", "Pastel1", "Pastel2", "Dark2",
               "Accent", "Paired"]
PERCEPTUAL = ["viridis", "plasma", "inferno", "magma", "cividis"]
BUNDLED_MAPS = QUALITATIVE + PERCEPTUAL

NOTICE = [
    "// Colour data of Matplotlib 3.11.2 (lib/matplotlib/_color_data.py, _cm.py, _cm_listed.py,",
    "// rcsetup.py). Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved.",
    "// Copyright (c) 2002-2011 John D. Hunter; All Rights Reserved. Used under the Matplotlib licence",
    "// (packaging/licenses/Matplotlib-LICENSE.txt). Change from Matplotlib: the data is converted to",
    "// RGBA doubles (the values matplotlib.colors.to_rgba gives) and laid out as C++ tables.",
    "// Included data: the xkcd colour survey names (Randall Munroe, CC0 1.0); the CSS4 colour names",
    "// (W3C); the Tableau 10 / 20 colours as distributed with Matplotlib; the ColorBrewer schemes",
    "// Set1-3, Pastel1-2, Dark2, Accent, Paired: Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and",
    "// The Pennsylvania State University, Apache License 2.0",
    "// (packaging/licenses/ColorBrewer-LICENSE.txt); viridis, plasma, inferno, magma (Nathaniel J.",
    "// Smith, Stefan van der Walt, CC0 1.0) and cividis (Jamie R. Nunez, Christopher R. Anderton, Ryan",
    "// S. Renslow) as distributed with Matplotlib.",
]


def cpp_string(text: str) -> str:
    out = []
    for ch in text:
        if ch in "\\\"":
            out.append("\\" + ch)
        elif " " <= ch <= "~":
            out.append(ch)
        else:
            out.append("".join(f"\\{b:03o}" for b in ch.encode("utf-8")))
    return '"' + "".join(out) + '"'


def hexf(value: float) -> str:
    value = float(value)
    assert value == value and abs(value) != float("inf"), value
    text = value.hex()  # '0x1.8000000000000p-1', '-0x0.0p+0', '0x0.0p+0'
    mantissa, exponent = text.split("p")
    whole, _, fraction = mantissa.partition(".")
    fraction = fraction.rstrip("0")
    return f"{whole}.{fraction}p{exponent}" if fraction else f"{whole}p{exponent}"


def main(argv: list[str]) -> int:
    if argv:
        print(__doc__, file=sys.stderr)
        return 2
    import matplotlib
    import matplotlib.colors as mcolors

    if matplotlib.__version__ != MATPLOTLIB_VERSION:
        print(f"gen_mpl_tables.py needs Matplotlib {MATPLOTLIB_VERSION} (found {matplotlib.__version__})",
              file=sys.stderr)
        return 1

    names = dict(mcolors.get_named_colors_mapping())
    rows = []
    for name in sorted(names, key=lambda n: n.encode("utf-8")):
        value = names[name]
        assert isinstance(value, tuple) or (isinstance(value, str) and len(value) == 7), (name, value)
        rgba = mcolors.to_rgba(name)
        assert rgba == mcolors.to_rgba(value) and rgba[3] == 1.0, name
        assert name == name.lower(), name  # the case-folded lookup of to_rgba relies on it
        row = f"    {{{cpp_string(name)}, {{{', '.join(hexf(v) for v in rgba[:3])}}}}},"
        if len(row) > WIDTH:
            row = f"    {{{cpp_string(name)},\n     {{{', '.join(hexf(v) for v in rgba[:3])}}}}},"
        rows.append(row)

    cycle = matplotlib.rcParamsDefault["axes.prop_cycle"].by_key()["color"]
    cycle_rows = []
    for color in cycle:
        rgba = mcolors.to_rgba(color)
        assert rgba[3] == 1.0
        cycle_rows.append(f"    {{{', '.join(hexf(v) for v in rgba[:3])}}},  // {mcolors.to_hex(rgba)}")

    map_rows, color_rows, first = [], [], 0
    for name in BUNDLED_MAPS:
        cmap = matplotlib.colormaps[name]
        assert type(cmap) is mcolors.ListedColormap, name
        lut = mcolors.to_rgba_array(cmap.colors)
        assert lut.shape == (cmap.N, 4) and (lut[:, 3] == 1.0).all(), name
        assert (cmap.N <= 20) == (name in QUALITATIVE), name  # MDATK's palette() tells them apart by N
        map_rows.append(f"    {{{cpp_string(name)}, {first}, {cmap.N}}},")
        color_rows.append(f"    // {name}")
        for r, g, b, _a in lut:
            color_rows.append(f"    {{{hexf(r)}, {hexf(g)}, {hexf(b)}}},")
        first += cmap.N

    lines = [
        "// Matplotlib's named colours, default colour cycle and the bundled colormaps.",
        f"// Generated by tools/gen_mpl_tables.py from Matplotlib {matplotlib.__version__}; do not edit.",
        "// Included by src/render/core/color.cpp and colormap.cpp inside namespace mdbridge::render::detail.",
        "//",
        *NOTICE,
        "",
        f"// {len(rows)} names, sorted by bytes; RGB of to_rgba(name) (alpha 1)",
        "inline constexpr NamedColor kNamedColors[] = {",
        *rows,
        "};",
        "",
        "// rcParamsDefault[\"axes.prop_cycle\"]: the colours of \"C0\", \"C1\", ...",
        "inline constexpr Rgb kDefaultCycle[] = {",
        *cycle_rows,
        "};",
        "",
        "// name, first row in kListedMapColors, colours (N)",
        "inline constexpr ListedMapRow kListedMaps[] = {",
        *map_rows,
        "};",
        "",
        "inline constexpr Rgb kListedMapColors[] = {",
        *color_rows,
        "};",
        "",
    ]
    for line in "\n".join(lines).split("\n"):
        assert len(line) <= WIDTH, line
    os.makedirs(os.path.dirname(INC), exist_ok=True)
    with open(INC, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines))
    print(f"wrote {os.path.relpath(INC, ROOT)}: {len(rows)} named colours, {len(cycle)} cycle colours, "
          f"{len(BUNDLED_MAPS)} colormaps")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
