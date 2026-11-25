"""What Matplotlib's colour functions and MDATK's themes give, for src/render/core
(tests/render/test_colors.cpp).

    python3 -I -B make_color_cases.py MDATK_CHECKOUT OUT_DIR

Imports Matplotlib 3.11.2 and MDATK (a checkout at v0.0.4.1; nothing is written into it:
mdatk/render/theme.py and coloring.py) and writes OUT_DIR/color_cases.json:

  to_rgba        [color, alpha, rgba | {"error": message}]: every named colour, case variants,
                 "none", "C<n>", hex of every length and case (and invalid ones), grey levels in
                 Python's float() syntax, Unicode spellings, alpha overrides and invalid alphas
  to_rgba_tuple  [components, alpha, rgba | {"error"}] (NaN components: improvements_colors)
  to_hex         [rgba, keep_alpha, hex]
  names          the number of named colours
  colormaps      maps as specs (bundled name, from_list of colours or (value, colour) pairs,
                 segment data; reversed; with_extremes) with N, under / over / bad, the whole
                 lookup table for N <= 16, and samples: [x, alpha, rgba] (float input),
                 [i, alpha, rgba] (integer input), [x, alpha, bytes] (bytes=True)
  colormap_errors  [spec, message] of invalid maps (Matplotlib's message, raised at first use)
  palettes       [name, n, colours] of MDATK's coloring.palette for every bundled map (and "_r")
  colormap_names [name, problem]: MDBridge's message for names that are not bundled
                 (improvements_colors.colormap_problem), None for bundled ones
  themes         per theme: the tokens, series, series_color(k), the sequential map (whole lookup
                 table, and with_extremes(bad=grid) as the heatmap panel uses it), ordinal(n)
                 (n > 5: improvements_colors.ordinal_colors)
  get_theme      [name, theme name]
  style_axes     [theme, grid, font_size, properties read back from a Matplotlib Axes]
  style_legend   [theme, font_size, properties read back from a Legend]

Floats are written with repr (exact); NaN and infinities as json.dumps writes them. The output is
deterministic (fixed random seed) and does not depend on the CPU (no BLAS; numpy's elementwise
float64 multiply / add / divide only).
"""

from __future__ import annotations

import json
import math
import os
import random
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_color_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))
sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

import improvements_colors as improved  # noqa: E402
from mdatk.render import coloring, theme as mtheme  # noqa: E402

assert matplotlib.__version__ == "3.11.2", matplotlib.__version__
RNG = random.Random(20261125)
NAMES = sorted(mcolors.get_named_colors_mapping(), key=lambda n: n.encode("utf-8"))


def plain(rgba):
    return [float(v) for v in rgba]


def rgba_or_error(fn):
    try:
        return plain(fn())
    except (ValueError, TypeError) as exc:
        return {"error": str(exc)}


# ---------------------------------------------------------------------------------------------
# to_rgba

def to_rgba_cases():
    inputs = []  # (color, alpha)
    for name in NAMES:
        inputs.append((name, None))
    for name in RNG.sample(NAMES, 200):
        inputs.append((name.upper(), None))
    for name in RNG.sample(NAMES, 80):
        inputs.append((name.title(), None))
        inputs.append((name.swapcase(), RNG.choice([0.0, 0.25, 1.0, None])))
    for name in RNG.sample(NAMES, 120):
        inputs.append((name, RNG.choice([0.0, 0.5, 1.0, RNG.random()])))
    inputs += [(c, None) for c in "RGBCMYKWxXzZ"]
    inputs += [("blacK", None), ("K", None), ("BLACK", None), ("Black", None), ("tab:Blue", None),
               ("TAB:BLUE", None), ("tab:grey", None), ("tab:gray", None), ("xkcd:Sky Blue", None),
               ("XKCD:SKY BLUE", None), ("xkcd:grey", None), ("xkcd:gray", None), ("redd", None), ("re d", None),
               ("red ", None), (" red", None), ("İ", None), ("gréy", None), ("", None), (" ", None),
               ("grey", None), ("GREY", None), ("Gray", None), ("darkgrey", None), ("darkgray", None)]
    # "none"
    for c in ["none", "None", "NONE", "nOnE", " none", "none "]:
        for alpha in [None, 0.3, 1.0, 1.5]:
            inputs.append((c, alpha))
    # colour cycle
    inputs += [(f"C{k}", None) for k in range(25)]
    inputs += [(c, None) for c in ["C00", "C01", "C09", "C010", "C123456789012345678901234567890", "c0", "c1", "C",
                                   "C-1", "C+1", "C1.5", " C1", "C1 ", "C١", "C1\n", "CC1", "C 1", "C1a"]]
    inputs += [("C3", 0.5), ("C3", 0.0), ("C3", 2.0), ("C12", 0.75)]
    # hex
    for digits in (3, 4, 6, 8):
        for _ in range(60):
            text = "".join(RNG.choice("0123456789abcdefABCDEF") for _ in range(digits))
            inputs.append(("#" + text, RNG.choice([None, None, None, 0.5, 0.0, 1.0])))
    for v in range(16):
        inputs.append((f"#{v:x}{v:x}{v:x}", None))
        inputs.append((f"#{v:X}{v:x}{15 - v:x}{v:x}", None))
    for v in range(0, 256, 5):
        inputs.append((f"#{v:02x}{255 - v:02x}{v:02X}", None))
        inputs.append((f"#{v:02x}{v:02x}{v:02x}{255 - v:02x}", None))
    for digits in (1, 2, 5, 7, 9, 10, 12):
        inputs.append(("#" + "".join(RNG.choice("0123456789abcdef") for _ in range(digits)), None))
    inputs += [(c, None) for c in ["#", "#ggg", "#12345g", "##123", "#fff ", " #fff", "#fff\n", "fff", "ffffff",
                                   "#١٢٣", "0x123", "#-12", "#+123"]]
    inputs += [("#11223344", 0.5), ("#1234", 0.0), ("#1234", None), ("#11223344", None), ("#abc", 1.0)]
    # grey levels: Python float() syntax
    greys = ["0", "1", "0.5", "0.25", ".5", "5.", " 0.5 ", "\t.5\n", "+0.5", "-0", "-0.0", "+0", "1e-1", "5E-1",
             "0.2_5", "1_0", "1__0", "_1", "1_", "0._5", "nan", "NaN", "-nan", "inf", "-inf", "Infinity",
             "1.0000001", "1.0000000000000002", "0.9999999999999999", "-0.1", "-1e-300", "1e-320", "2",
             "0.5.5", ".", "e1", "1e", "0x1", "１", "０.５", "٠.٥", " 0.5",
             " 0.5 ", "　0.5", "0.5\x00", "1.", "00.5", "0e0", "1E0", "0.75 ", "0,5"]
    inputs += [(g, None) for g in greys]
    inputs += [(g, 0.5) for g in ["0", "0.5", "1", "2", "nan"]]
    for _ in range(80):
        inputs.append((repr(RNG.random()), None))
    for _ in range(20):
        inputs.append((repr(RNG.uniform(-2.0, 3.0)), None))
    # invalid alphas: checked before the colour
    for c in ["red", "#123", "0.5", "none", "nonsense", "C1", "#12"]:
        for alpha in [-0.1, 1.1, float("nan"), float("inf"), -0.0]:
            inputs.append((c, alpha))
    out = []
    for color, alpha in inputs:
        mcolors.get_named_colors_mapping().cache.clear()
        out.append([color, alpha, rgba_or_error(lambda: mcolors.to_rgba(color, alpha))])
    return out


def tuple_cases():
    inputs = [((0, 0, 0), None), ((1, 1, 1), None), ((0.2, 0.4, 0.6), None), ((0.2, 0.4, 0.6, 0.8), None),
              ((0.2, 0.4, 0.6, 0.8), 0.1), ((0.2, 0.4, 0.6), 0.1), ((1.1, 0, 0), None), ((0, -0.1, 0), None),
              ((0, 0, 0, 1.5), None), ((0, 0, 0, 1.5), 0.5), ((0, 0), None), ((0, 0, 0, 0, 0), None), ((), None),
              ((float("nan"), 0, 0), None), ((0, 0, 0, float("nan")), None), ((0, 0, 0, float("nan")), 0.5),
              ((-0.0, 0, 0), None), ((float("inf"), 0, 0), None), ((0.5, 0.5, 0.5), 1.5), ((0.5, 0.5), 1.5)]
    for _ in range(20):
        inputs.append((tuple(RNG.random() for _ in range(RNG.choice([3, 4]))), RNG.choice([None, RNG.random()])))
    out = []
    for comps, alpha in inputs:
        # a list, not a tuple: to_rgba reads a 2-tuple as (colour, alpha); and no cached result
        # (the cache cannot tell -0.0 from 0.0)
        mcolors.get_named_colors_mapping().cache.clear()
        result = rgba_or_error(lambda: mcolors.to_rgba([float(v) for v in comps], alpha))
        out.append([[float(v) for v in comps], alpha, improved.tuple_result(comps, alpha, result)])
    return out


def hex_cases():
    values = [k / 255 for k in range(256)] + [(k + 0.5) / 255 for k in range(255)] + [0.5, 0.25, 1 / 3, 0.002]
    out = []
    for v in values:
        rgba = (v, 1 - v, v / 2, 1 - v / 3)
        for keep in (False, True):
            out.append([plain(rgba), keep, mcolors.to_hex(rgba, keep_alpha=keep)])
    for _ in range(60):
        rgba = tuple(RNG.random() for _ in range(4))
        out.append([plain(rgba), RNG.choice([False, True]), None])
        out[-1][2] = mcolors.to_hex(rgba, keep_alpha=out[-1][1])
    return out


# ---------------------------------------------------------------------------------------------
# colormaps

def float_inputs(n):
    xs = [0.0, -0.0, 1.0, 0.5, 0.25, 0.75, 1 / 3, 2 / 3, -1e-300, -0.001, -1.0, 1.0000000000000002, 1.5,
          0.9999999999999999, float("nan"), float("inf"), float("-inf"), 5e-324, 0.05, 0.95, 1e9, -1e9]
    for k in sorted({1, n // 2, n - 1, max(n - 2, 0)}):
        b = k / n
        xs += [b, math.nextafter(b, -math.inf), math.nextafter(b, math.inf)]
    xs += [RNG.random() for _ in range(12)]
    return xs


def int_inputs(n):
    return sorted({-5, -1, 0, 1, n // 2, n - 1, n, n + 1, 10 ** 12})


def build(spec):
    kind = spec["kind"]
    if kind == "bundled":
        cmap = matplotlib.colormaps[spec["name"]]
    elif kind == "from_list":
        colors = spec["colors"]
        if spec.get("values") is not None:
            colors = list(zip(spec["values"], colors))
        cmap = mcolors.LinearSegmentedColormap.from_list(spec["name"], colors, N=spec["N"], gamma=spec["gamma"])
    elif kind == "segments":
        data = {k: [tuple(r) for r in v] for k, v in spec["data"].items()}
        cmap = mcolors.LinearSegmentedColormap(spec["name"], data, N=spec["N"], gamma=spec["gamma"])
    elif kind == "listed":
        cmap = mcolors.ListedColormap(spec["colors"], name=spec["name"])
    else:
        raise AssertionError(kind)
    if spec.get("reversed"):
        cmap = cmap.reversed()
    ext = spec.get("extremes")
    if ext:
        cmap = cmap.with_extremes(**ext)
    cmap._init()
    return cmap


def map_record(spec, alphas=(None,)):
    cmap = build(spec)
    n = cmap.N
    rec = {"spec": spec, "name": cmap.name, "N": n, "under": plain(cmap.get_under()), "over": plain(cmap.get_over()),
           "bad": plain(cmap.get_bad())}
    if n <= 16:
        rec["lut"] = [plain(row) for row in cmap._lut]
    samples, ints, raw = [], [], []
    for x in float_inputs(n):
        for alpha in alphas:
            samples.append([x, alpha, plain(cmap(x, alpha=alpha))])
            if not (alpha is not None and math.isnan(alpha)):
                raw.append([x, alpha, [int(v) for v in cmap(x, alpha=alpha, bytes=True)]])
    for i in int_inputs(n):
        for alpha in alphas:
            ints.append([i, alpha, plain(cmap(i, alpha=alpha))])
            if not (alpha is not None and math.isnan(alpha)):
                raw.append([i, alpha, [int(v) for v in cmap(np.int64(i), alpha=alpha, bytes=True)], "int"])
    rec["samples"], rec["int_samples"], rec["bytes"] = samples, ints, raw
    return rec


def random_colors(k):
    pool = ["red", "#123456", "#abc", "0.5", "tab:orange", "xkcd:sky blue", "C3", "#12345678", "none", "k",
            "#fcfcfb", "navy"]
    return [RNG.choice(pool) for _ in range(k)]


def colormap_cases():
    specs = []
    for name in improved.BUNDLED_MAPS:
        specs.append(({"kind": "bundled", "name": name}, (None,)))
        specs.append(({"kind": "bundled", "name": name + "_r"}, (None,)))
    for name in ("viridis", "tab10", "Set2"):
        specs.append(({"kind": "bundled", "name": name}, (0.0, 0.5, 1.0, 1.5, -0.5, float("nan"))))
        specs.append(({"kind": "bundled", "name": name,
                       "extremes": {"bad": "0.5", "under": "r", "over": "#00ff0080"}}, (None, 0.5)))
        specs.append(({"kind": "bundled", "name": name, "extremes": {"bad": "none"}}, (None, 0.5)))
        specs.append(({"kind": "bundled", "name": name, "extremes": {"under": "C1"}}, (None,)))
    for n in (1, 2, 3, 5, 16, 256, 1000):
        for k in (2, 3, 6):
            specs.append(({"kind": "from_list", "name": f"fl_{n}_{k}", "colors": random_colors(k), "N": n,
                           "gamma": 1.0}, (None,)))
    specs.append(({"kind": "from_list", "name": "pairs", "colors": ["red", "green", "blue", "white"],
                   "values": [0.0, 0.3, 0.3, 1.0], "N": 11, "gamma": 1.0}, (None,)))
    specs.append(({"kind": "from_list", "name": "pairs2", "colors": ["#000000", "#ff000080", "#00ff00", "w"],
                   "values": [0.0, 0.2, 0.9, 1.0], "N": 256, "gamma": 1.0}, (None, 0.25)))
    specs.append(({"kind": "from_list", "name": "pairs3", "colors": ["k", "w"], "values": [0.0, 1.0], "N": 7,
                   "gamma": 1.0, "reversed": True}, (None,)))
    specs.append(({"kind": "from_list", "name": "rev", "colors": random_colors(5), "N": 256, "gamma": 1.0,
                   "reversed": True, "extremes": {"under": "k", "over": "w", "bad": "C2"}}, (None,)))
    saw = {"red": [[0.0, 1.0, 0.0], [0.2, 1.0, 0.0], [0.4, 1.0, 0.0], [0.6, 1.0, 0.0], [0.8, 1.0, 0.0],
                   [1.0, 1.0, 0.0]],
           "green": [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.75, 1.0, 1.0], [1.0, 1.0, 1.0]],
           "blue": [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 1.0], [1.0, 1.0, 1.0]]}
    for n in (1, 2, 11, 256):
        specs.append(({"kind": "segments", "name": f"saw{n}", "data": saw, "N": n, "gamma": 1.0}, (None,)))
    specs.append(({"kind": "segments", "name": "saw_r", "data": saw, "N": 11, "gamma": 1.0, "reversed": True},
                  (None,)))
    alpha_saw = dict(saw, alpha=[[0.0, 0.0, 0.0], [1.0, 2.0, 2.0]])  # clipped to 1
    specs.append(({"kind": "segments", "name": "alpha_saw", "data": alpha_saw, "N": 9, "gamma": 1.0}, (None,)))
    specs.append(({"kind": "listed", "name": "lst", "colors": ["r", "#00ff0080", "0.25"]}, (None, 0.5)))
    specs.append(({"kind": "listed", "name": "lst1", "colors": ["navy"], "reversed": True}, (None,)))
    return [map_record(spec, alphas) for spec, alphas in specs]


def colormap_errors():
    specs = [
        {"kind": "from_list", "name": "one", "colors": ["red"], "N": 256, "gamma": 1.0},
        {"kind": "from_list", "name": "dec", "colors": ["red", "blue", "green"], "values": [0.0, 0.6, 0.4],
         "N": 256, "gamma": 1.0},
        {"kind": "from_list", "name": "neg", "colors": ["red", "blue"], "values": [-0.1, 1.0], "N": 256, "gamma": 1.0},
        {"kind": "from_list", "name": "big", "colors": ["red", "blue"], "values": [0.0, 1.1], "N": 256, "gamma": 1.0},
        {"kind": "from_list", "name": "short", "colors": ["red", "blue"], "values": [0.0, 0.9], "N": 256,
         "gamma": 1.0},
        {"kind": "from_list", "name": "late", "colors": ["red", "blue"], "values": [0.1, 1.0], "N": 256, "gamma": 1.0},
        {"kind": "segments", "name": "segdec", "N": 16, "gamma": 1.0,
         "data": {"red": [[0.0, 0.0, 0.0], [0.7, 1.0, 1.0], [0.5, 1.0, 1.0], [1.0, 1.0, 1.0]],
                  "green": [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], "blue": [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]}},
        {"kind": "segments", "name": "segend", "N": 16, "gamma": 1.0,
         "data": {"red": [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], "green": [[0.0, 0.0, 0.0], [0.9, 1.0, 1.0]],
                  "blue": [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]}},
    ]
    out = []
    for spec in specs:
        try:
            build(spec)(0.5)
        except ValueError as exc:
            out.append([spec, str(exc)])
        else:
            raise AssertionError(spec)
    return out


def palette_cases():
    out = []
    for base in improved.BUNDLED_MAPS:
        for name in (base, base + "_r"):
            for n in (0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 20, 21, 25, 40):
                out.append([name, n, [plain(c) for c in coloring.palette(name, n)]])
    return out


def colormap_name_cases():
    names = list(improved.BUNDLED_MAPS) + [n + "_r" for n in improved.BUNDLED_MAPS]
    names += ["jet", "Blues", "set2", "SET2", "Viridis", "TAB10_R", "viridis_r_r", "_r", "", "turbo", "tab10 ",
              "x_r", "okabe_ito", "grey", "Greys", "coolwarm", "Spectral", "gist_rainbow"]
    out = []
    for name in names:
        known = name in matplotlib.colormaps
        out.append([name, known, improved.colormap_problem(name)])
    return out


# ---------------------------------------------------------------------------------------------
# MDATK's themes

def theme_cases():
    out = {}
    for t in (mtheme.LIGHT, mtheme.DARK):
        rec = {k: getattr(t, k) for k in ("name", "surface", "page", "ink", "ink2", "muted", "grid", "axis")}
        rec["series"] = list(t.series)
        rec["series_color"] = [[k, t.series_color(k)] for k in range(-10, 21)]
        seq = t.sequential()
        seq._init()
        rec["sequential"] = {"name": seq.name, "N": seq.N, "lut": [plain(row) for row in seq._lut]}
        bad = seq.with_extremes(bad=t.grid)
        bad._init()
        rec["sequential_bad"] = {"bad": plain(bad.get_bad()),
                                 "samples": [[x, plain(bad(x))] for x in (float("nan"), 0.0, 0.5, 1.0, -1.0, 2.0)]}
        ordinals = []
        for n in range(0, 14):
            cmap = t.ordinal(n)
            cmap._init()
            mdatk_colors = [plain(row) for row in cmap._lut[: cmap.N]]
            colors = improved.ordinal_colors(t.name, n, mdatk_colors)
            ordinals.append({"n": n, "name": cmap.name, "colors": colors, "improved": colors != mdatk_colors})
        rec["ordinal"] = ordinals
        out[t.name] = rec
    return out


def get_theme_cases():
    return [[name, mtheme.get_theme(name).name] for name in ("light", "dark", "Dark", "DARK", "", "x", " dark")]


def tick_props(axis):
    tick = axis.get_major_ticks()[0]
    return {"tick_color": plain(mcolors.to_rgba(tick.tick1line.get_color())),
            "tick_width": float(tick.tick1line.get_markeredgewidth()),
            "tick_length": float(tick.tick1line.get_markersize()),
            "label_color": plain(mcolors.to_rgba(tick.label1.get_color())),
            "label_size": float(tick.label1.get_fontsize()),
            "grid_visible": bool(tick.gridline.get_visible()),
            "grid_color": plain(mcolors.to_rgba(tick.gridline.get_color())),
            "grid_width": float(tick.gridline.get_linewidth()),
            "grid_linestyle": tick.gridline.get_linestyle(),
            "axis_label_color": plain(mcolors.to_rgba(axis.label.get_color())),
            "axis_label_size": float(axis.label.get_fontsize())}


def style_cases():
    axes_out, legend_out = [], []
    for t in (mtheme.LIGHT, mtheme.DARK):
        for grid in (True, False):
            for size in (10.0, 11.0, 7.5):
                fig = Figure()
                ax = fig.add_subplot(1, 1, 1)
                mtheme.style_axes(ax, t, grid=grid, font_size=size)
                props = {"face": plain(ax.get_facecolor()),
                         "spines": {side: [bool(s.get_visible()), plain(s.get_edgecolor()), float(s.get_linewidth())]
                                    for side, s in ax.spines.items()},
                         "x": tick_props(ax.xaxis), "y": tick_props(ax.yaxis),
                         "title_color": plain(mcolors.to_rgba(ax.title.get_color())),
                         "title_size": float(ax.title.get_fontsize()),
                         "axisbelow": str(ax.get_axisbelow())}
                axes_out.append([t.name, grid, size, props])
        for size in (10.0, 11.0):
            fig = Figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot([0, 1], [0, 1], label="a")
            leg = ax.legend()
            mtheme.style_legend(leg, t, size)
            frame = leg.get_frame()
            text = leg.get_texts()[0]
            legend_out.append([t.name, size, {"face": plain(frame.get_facecolor()),
                                              "edge": plain(frame.get_edgecolor()),
                                              "edge_width": float(frame.get_linewidth()),
                                              "text_color": plain(mcolors.to_rgba(text.get_color())),
                                              "text_size": float(text.get_fontsize())}])
    return axes_out, legend_out


def main(out_dir: str) -> None:
    axes, legends = style_cases()
    doc = {
        "matplotlib": matplotlib.__version__,
        "to_rgba": to_rgba_cases(),
        "to_rgba_tuple": tuple_cases(),
        "to_hex": hex_cases(),
        "names": len(NAMES),
        "colormaps": colormap_cases(),
        "colormap_errors": colormap_errors(),
        "palettes": palette_cases(),
        "colormap_names": colormap_name_cases(),
        "themes": theme_cases(),
        "get_theme": get_theme_cases(),
        "style_axes": axes,
        "style_legend": legends,
    }
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "color_cases.json")
    # one line per list entry (readable diffs), compact inside
    dumps = lambda value: json.dumps(value, ensure_ascii=False, separators=(",", ":"))  # noqa: E731
    parts = []
    for key, value in doc.items():
        if isinstance(value, list):
            body = ",\n".join(dumps(item) for item in value)
            parts.append(f"{dumps(key)}:[\n{body}\n]")
        elif isinstance(value, dict):
            body = ",\n".join(f"{dumps(k)}:{dumps(v)}" for k, v in value.items())
            parts.append(f"{dumps(key)}:{{\n{body}\n}}")
        else:
            parts.append(f"{dumps(key)}:{dumps(value)}")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("{\n" + ",\n".join(parts) + "\n}\n")
    count = (len(doc["to_rgba"]) + len(doc["to_rgba_tuple"]) + len(doc["to_hex"]) + len(doc["palettes"])
             + sum(len(m["samples"]) + len(m["int_samples"]) + len(m["bytes"]) for m in doc["colormaps"]))
    print(f"wrote {path}: {count} cases")


if __name__ == "__main__":
    main(sys.argv[2])
