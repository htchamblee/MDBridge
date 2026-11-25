"""What Matplotlib 3.11.2 (as MDATK v0.0.4.1 uses it) and MDATK compute for axis limits, tick
locations and tick labels, for src/render/core/ticks.cpp (tests/render/test_ticks.cpp).

    /opt/mdbridge-venv/bin/python -I -B make_tick_cases.py MDATK_CHECKOUT OUT_DIR

Writes OUT_DIR/tick_cases.json, one case per line in each section:

  nonsingular   transforms._nonsingular(vmin, vmax, expander, tiny, increasing)
  locator       MaxNLocator / AutoLocator tick_values and view_limits for many intervals (signs,
                magnitudes 1e-300..1e300, widths down to a few ulps, large offsets), nbins (fixed or
                'auto' with an axis' tick space 0..12), steps, integer, symmetric, prune,
                min_n_ticks; Matplotlib's ValueError messages for invalid steps
  fixed         FixedLocator(locs, nbins).tick_values
  formatter     ScalarFormatter.format_ticks + get_offset of located ticks for a view interval:
                automatic / disabled / fixed offsets, offset thresholds, power limits (also fixed
                exponents), scientific off, mathtext, ASCII minus; empty and single locations
  format_data   ScalarFormatter.format_data of numpy floats (computed offsets) and Python floats
  axes          real Axes of a Figure (sizes, dpi, axes rectangles, tick label sizes, x and y):
                limits set or autoscaled from plotted data (margins 5 %, sticky edges of a line,
                tight), MDATK's integer_axis; the axis length in points, get_tick_space, the view
                interval, get_majorticklocs, the tick labels and offset text as Axis._update_ticks
                makes them, and which ticks it draws
  set_lim       set_xlim(v0, v1) of an Axes (orientation kept), imshow's extents
  colorbar      fig.colorbar(imshow(vmin, vmax), fraction 0.04, pad 0.02) as MDATK's heatmap makes
                it (label size 0.75 x font, integer_axis for counts): the long axis' limits, tick
                space, ticks (get_ticks), labels and offset text; vmin <= vmax only (MDBridge
                clamps vmax first, RN-04; Matplotlib's colour bar of a reversed norm depends on
                the order its norm callbacks run in)
  thin_ticks    mdatk.render.panels.thin_ticks(n, max_ticks); max_ticks < 1 is MDBridge's
                improvement (taken as 1; MDATK raises ZeroDivisionError for 0 and labels every
                row for negative values), docs/DIFFERENCES.md

Floats are JSON numbers (repr; NaN / Infinity as json writes them). Matplotlib's rcParams are its
defaults (MDATK changes none that matter here). Results depend on libm's pow / log10 and numpy's
scalar loops only, so numpy's dispatch is pinned to its x86-64 baseline as in
make_analysis_cases.py (np.log10 may otherwise use SIMD kernels); the output is the same on every
machine and run.
"""

import os
import platform
import sys
import tempfile

if platform.machine().lower() in ("x86_64", "amd64"):  # before numpy is imported
    os.environ["NPY_ENABLE_CPU_FEATURES"] = "X86_V2"
os.environ["MPLBACKEND"] = "Agg"
os.environ["MPLCONFIGDIR"] = tempfile.mkdtemp(prefix="mdbridge-mpl-")  # no user matplotlibrc

if len(sys.argv) != 3:
    sys.exit("usage: make_tick_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import json  # noqa: E402
import logging  # noqa: E402
import math  # noqa: E402
import warnings  # noqa: E402

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
mpl.rcdefaults()
logging.getLogger("matplotlib").setLevel(logging.ERROR)  # MAXTICKS notes
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from matplotlib import ticker  # noqa: E402
from matplotlib import transforms as mtransforms  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from mdatk.render.panels import integer_axis, thin_ticks  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
rng = np.random.default_rng(20261125)


def f(x):
    return float(x)


def fl(xs):
    return [float(x) for x in xs]


# ---------------------------------------------------------------------------------------------
# intervals

def intervals():
    out = [
        (0.0, 1.0), (0.0, 0.0), (1.0, 1.0), (-1.0, -1.0), (-1.0, 1.0), (1.0, 0.0), (0.0, 10.0),
        (0.0, 100.0), (0.0, 7.0), (0.0, 13.0), (-0.5, 9.5), (0.5, 119.5), (-0.5, 15.5), (0.0, 1e-9),
        (1e6, 1e6 + 10), (1_000_000.0, 1_000_010.0), (-5e-7, 3e-7), (1e-300, 2e-300),
        (12345.678, 12345.679), (99.9999, 100.0001), (-100.0001, -99.9999), (0.1, 0.3),
        (1e15, 1e15 + 1), (2.0, 2.0 + 4.4e-16), (5.0, 5.0 + 1e-12), (-3.2, 47.1), (0.0, 0.03),
        (0.0, 3e5), (0.0, 2.4e6), (1.5e-5, 3.5e-5), (-2e7, 3e7), (10.0, 10.5), (1e300, 1e301),
        (-1e300, 1e300), (0.0, 120.0), (0.0, 24.0), (0.0, 1.2), (0.002, 0.012), (-7.5, -2.5),
        (0.95, 1.05), (19999.5, 20000.5), (0.0, 1e100), (123456789.0, 123456789.5),
        (float("inf"), 1.0), (float("nan"), 0.0), (-float("inf"), float("inf")),
    ]
    for _ in range(140):
        mag = 10 ** rng.uniform(-12, 12)
        center = rng.choice([0.0, 1.0, -1.0]) * mag * rng.uniform(0, 3)
        rel = 10 ** rng.uniform(-13, 2)
        width = max(abs(center), mag) * rel
        a = center - width * rng.uniform(0, 1)
        b = a + width
        if rng.uniform() < 0.1:
            a, b = b, a
        out.append((float(a), float(b)))
    for _ in range(40):  # "nice" data ranges of time axes and counts
        a = float(rng.integers(-50, 50)) * float(rng.choice([1, 0.1, 0.5, 2, 10, 0.001]))
        b = a + float(rng.integers(1, 400)) * float(rng.choice([1, 0.25, 0.02, 3, 100]))
        out.append((a, b))
    return out


class StubAxis:
    """An axis object whose only use by MaxNLocator is get_tick_space (no 'axes' attribute)."""

    def __init__(self, space):
        self.space = space

    def get_tick_space(self):
        return self.space


# ---------------------------------------------------------------------------------------------
# sections

def nonsingular_cases():
    pairs = [(0.0, 0.0), (1.0, 1.0), (-1.0, -1.0), (1e-300, 1e-300), (1e-290, 2e-290), (2.0, 1.0),
             (float("inf"), 1.0), (float("nan"), 0.0), (-float("inf"), float("inf")), (1.0, 1.0 + 1e-16),
             (1e10, 1e10 + 1e-6), (-5.0, -5.0), (0.0, 1e-310), (3.0, 3.0 + 1e-14), (-1e-14, 1e-14),
             (1e308, 1e308), (-2.5, 7.25), (100.0, -100.0), (0.0, -0.0), (5e-324, 0.0)]
    params = [(0.001, 1e-15), (0.05, 1e-15), (1e-13, 1e-14), (1e-12, 1e-13), (0.1, 1e-15)]
    out = []
    for vmin, vmax in pairs:
        for expander, tiny in params:
            for increasing in (True, False):
                r = mtransforms._nonsingular(vmin, vmax, expander=expander, tiny=tiny, increasing=increasing)
                out.append({"vmin": vmin, "vmax": vmax, "expander": expander, "tiny": tiny,
                            "increasing": increasing, "result": fl(r)})
    return out


LOCATOR_PARAMS = [
    {"kind": "auto"},
    {"kind": "integer"},
    {"nbins": 10},
    {"nbins": 5},
    {"nbins": 1},
    {"nbins": 3, "steps": [1, 2, 5, 10]},
    {"nbins": 7, "steps": [1, 2, 4, 5, 10]},
    {"nbins": 4, "steps": [2, 5]},
    {"nbins": 6, "steps": [1, 3, 10]},
    {"nbins": 12, "integer": True},
    {"nbins": 4, "integer": True, "min_n_ticks": 3},
    {"nbins": 8, "symmetric": True},
    {"nbins": 9, "prune": "lower"},
    {"nbins": 9, "prune": "upper"},
    {"nbins": 9, "prune": "both"},
    {"nbins": 2, "min_n_ticks": 4},
    {"nbins": 10, "min_n_ticks": 0},
    {"nbins": "auto", "steps": [1, 2, 2.5, 5, 10], "integer": True},
]


def make_locator(p, space):
    if p.get("kind") == "auto":
        loc = ticker.AutoLocator()
    elif p.get("kind") == "integer":
        loc = ticker.MaxNLocator(integer=True)
    else:
        kw = {k: v for k, v in p.items() if k != "kind"}
        loc = ticker.MaxNLocator(**kw)
    if space is not None:
        loc.set_axis(StubAxis(space))
    return loc


def locator_cases():
    out = []
    ivals = intervals()
    for i, (vmin, vmax) in enumerate(ivals):
        for j, p in enumerate(LOCATOR_PARAMS):
            auto = p.get("kind") == "auto" or p.get("nbins") == "auto"
            spaces = [None]
            if auto:
                spaces = [None, int(rng.integers(0, 13))]
            elif (i + j) % 3:  # every interval with every parameter set would be too many
                continue
            for space in spaces:
                loc = make_locator(p, space)
                case = {"params": p, "tick_space": space, "vmin": vmin, "vmax": vmax}
                try:
                    case["ticks"] = fl(loc.tick_values(vmin, vmax))
                    case["view_limits"] = fl(loc.view_limits(vmin, vmax))
                except Exception as exc:  # noqa: BLE001
                    case["error"] = f"{type(exc).__name__}: {exc}"
                out.append(case)
    for steps in ([], [0.5, 1], [1, 1], [2, 1], [1, 11], [5, 2, 10]):
        try:
            ticker.MaxNLocator(steps=steps)
            msg = None
        except ValueError as exc:
            msg = str(exc)
        except IndexError as exc:
            msg = f"IndexError: {exc}"
        out.append({"params": {"nbins": 10, "steps": steps}, "invalid_steps": msg})
    return out


def fixed_cases():
    out = []
    for n in (0, 1, 2, 3, 5, 9, 10, 11, 17, 40):
        base = list(range(n))
        for locs in (base, [x - n // 2 for x in base], [0.5 * x - 3.25 for x in base]):
            for nbins in (None, 1, 3, 4, 10, n if n else None):
                loc = ticker.FixedLocator(locs, nbins=nbins)
                out.append({"locs": fl(locs), "nbins": nbins, "ticks": fl(loc.tick_values(None, None))})
    return out


FORMATTER_OPTIONS = [
    {},
    {"useOffset": False},
    {"useOffset": 1000.0},
    {"useOffset": -2.5},
    {"useMathText": True},
    {"unicode_minus": False},
    {"powerlimits": [-3, 4]},
    {"powerlimits": [0, 0]},
    {"powerlimits": [3, 3]},
    {"powerlimits": [-2, -2]},
    {"scientific": False},
    {"offset_threshold": 2},
    {"offset_threshold": 6},
    {"useMathText": True, "useOffset": 1000.0},
]


def make_formatter(o):
    rc = {"axes.unicode_minus": o.get("unicode_minus", True),
          "axes.formatter.offset_threshold": o.get("offset_threshold", 4)}
    with mpl.rc_context(rc):
        fmt = ticker.ScalarFormatter(useOffset=o.get("useOffset"), useMathText=o.get("useMathText"))
    if "powerlimits" in o:
        fmt.set_powerlimits(tuple(o["powerlimits"]))
    if "scientific" in o:
        fmt.set_scientific(o["scientific"])
    fmt.create_dummy_axis()
    return fmt, rc


def format_case(o, locs, view):
    fmt, rc = make_formatter(o)
    fmt.axis.set_view_interval(*view)
    try:
        with mpl.rc_context(rc):  # fix_minus reads axes.unicode_minus when called
            labels = fmt.format_ticks(locs)
            offset_text = fmt.get_offset()
    except ValueError as exc:  # a fixed offset and an empty view interval: log10(0)
        return {"options": o, "locs": fl(locs), "view": fl(view), "error": str(exc)}
    sigfigs = None
    if fmt._format:
        sigfigs = int(fmt._format.split("%1.")[1].split("f")[0])
    return {"options": o, "locs": fl(locs), "view": fl(view), "labels": labels, "offset_text": offset_text,
            "offset": f(fmt.offset), "order_of_magnitude": int(fmt._orderOfMagnitude), "decimals": sigfigs}


def formatter_cases():
    out = []
    auto = ticker.AutoLocator()
    ivals = [iv for iv in intervals() if all(math.isfinite(v) for v in iv)]
    for i, (a, b) in enumerate(ivals):
        view = mtransforms._nonsingular(a, b, expander=0.05)
        if a > b:
            view = (view[1], view[0])
        locs = auto.tick_values(*view)
        for j, o in enumerate(FORMATTER_OPTIONS):
            if j and (i + j) % 4:
                continue
            out.append(format_case(o, locs, view))
    # ticks set by hand: one, two equal, outside the view, zero, negative
    special = [([], (0.0, 1.0)), ([0.5], (0.0, 1.0)), ([1e6], (999990.0, 1000010.0)), ([0.0], (-1.0, 1.0)),
               ([0.0], (0.0, 0.0)), ([3.0, 3.0], (2.0, 4.0)), ([10.0, 20.0, 30.0], (100.0, 200.0)),
               ([-1e-7, 0.0, 1e-7], (-1e-7, 1e-7)), ([1.0000001, 1.0000002, 1.0000003], (1.0, 1.0000004)),
               ([-1000002.0, -1000001.0, -1000000.0], (-1000002.5, -999999.5)), ([0.1, 0.2, 0.30000000000000004], (0.1, 0.3)),
               ([1e-30, 2e-30, 3e-30], (0.0, 4e-30)), ([2e20, 4e20], (0.0, 5e20)), ([0, 1, 2, 3], (0.0, 3.0)),
               ([-2.0, -1.0], (-2.0, -1.0)), ([5.0], (5.0, 5.0))]
    for locs, view in special:
        for o in FORMATTER_OPTIONS:
            out.append(format_case(o, np.asarray(locs, dtype=np.float64), view))
    return out


def format_data_cases():
    out = []
    values = [1000.0, 1500.0, -2.5e-4, 7.0, 1e6, 123456.0, -3.0e12, 0.1, 9.99999999999, 1.00048828125,
              2.0 ** -20, 6.02e23, 1e23, 1e-5, 4.5, 12.0, -0.001]
    values += [float(10 ** rng.uniform(-20, 20) * rng.choice([1, -1])) for _ in range(40)]
    for math_text in (False, True):
        for minus in (True, False):
            with mpl.rc_context({"axes.unicode_minus": minus}):
                fmt = ticker.ScalarFormatter(useMathText=math_text)
                for v in values:
                    out.append({"value": v, "numpy": True, "mathtext": math_text, "unicode_minus": minus,
                                "text": fmt.format_data(np.float64(v))})
                    out.append({"value": v, "numpy": False, "mathtext": math_text, "unicode_minus": minus,
                                "text": fmt.format_data(v)})
    return out


def axis_record(ax, axis, fig):
    ends = mtransforms.Bbox.unit().transformed(ax.transAxes - fig.dpi_scale_trans)
    name = axis.axis_name
    length = (ends.width if name == "x" else ends.height) * 72
    view = axis.get_view_interval()
    locs = axis.get_majorticklocs()
    fmt = axis.get_major_formatter()
    labels = fmt.format_ticks(locs)
    drawn_locs = [t.get_loc() for t in axis._update_ticks()]
    drawn = [i for i, x in enumerate(locs) if x in drawn_locs]
    loc = axis.get_major_locator()
    return {"length_points": f(length), "label_size": f(axis._get_tick_label_size(name)),
            "tick_space": int(axis.get_tick_space()), "integer": isinstance(loc, ticker.MaxNLocator)
            and not isinstance(loc, ticker.AutoLocator), "view": fl(view), "locs": fl(locs), "labels": labels,
            "offset_text": fmt.get_offset(), "drawn": drawn}


def axes_cases():
    out = []
    figures = [((16.0, 9.0), 120), ((8.0, 4.5), 120), ((6.4, 4.8), 100), ((2.6666666666666665, 1.5), 120),
               ((10.0, 3.0), 72), ((4.0, 8.0), 150)]
    rects = [(0.1, 0.1, 0.8, 0.8), (0.05, 0.55, 0.4, 0.4), (0.06, 0.08, 0.27, 0.18), (0.5, 0.1, 0.45, 0.05)]
    sizes = [8.5, 10.2, 7.0, 0.85 * 14, 0.85 * 6, 20.0, 0.5, 3.4]
    ivals = [iv for iv in intervals() if all(math.isfinite(v) for v in iv)]
    k = 0
    for fi, (figsize, dpi) in enumerate(figures):
        for ri, rect in enumerate(rects):
            for si, size in enumerate(sizes):
                if (fi + ri + si) % 2:
                    continue
                for mode in ("limits", "data", "sticky", "tight", "integer"):
                    fig = Figure(figsize=figsize, dpi=dpi)
                    ax = fig.add_axes(rect)
                    ax.tick_params(labelsize=size)
                    xa, xb = ivals[k % len(ivals)]
                    ya, yb = ivals[(7 * k + 3) % len(ivals)]
                    k += 1
                    case = {"figsize": fl(figsize), "dpi": dpi, "rect": fl(rect), "label_size": size, "mode": mode}
                    if mode in ("limits", "integer"):
                        ax.set_xlim(xa, xb)
                        ax.set_ylim(ya, yb)
                        case["xlim"] = [xa, xb]
                        case["ylim"] = [ya, yb]
                        if mode == "integer":
                            integer_axis(ax.yaxis)
                    else:
                        n = int(rng.integers(1, 30))
                        xs = np.sort(rng.uniform(min(xa, xb), max(xa, xb), n))
                        ys = rng.uniform(min(ya, yb), max(ya, yb), n)
                        (line,) = ax.plot(xs, ys)
                        if mode == "sticky":
                            sx = [float(xs[0]), float(xs[-1])]
                            line.sticky_edges.x[:] = sx
                            line.sticky_edges.y[:] = [float(min(ys))]
                            case["sticky_x"] = sx
                            case["sticky_y"] = [float(min(ys))]
                        if mode == "tight":
                            ax.autoscale_view(tight=True)
                        case["data_x"] = fl(ax.dataLim.intervalx)
                        case["data_y"] = fl(ax.dataLim.intervaly)
                    case["x"] = axis_record(ax, ax.xaxis, fig)
                    case["y"] = axis_record(ax, ax.yaxis, fig)
                    out.append(case)
    return out


def set_lim_cases():
    out = []
    pairs = [(0.0, 1.0), (1.0, 0.0), (5.0, 5.0), (0.0, 0.0), (-3.0, -3.0), (119.5, -0.5), (1e-300, 1e-300),
             (2.0, 2.0 + 1e-15), (1e9, 1e9), (-0.5, 0.5), (15.5, -0.5)]
    for v0, v1 in pairs:
        fig = Figure()
        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
        ax.set_xlim(v0, v1)
        out.append({"v0": v0, "v1": v1, "result": fl(ax.get_xlim())})
    for extent in ([-0.5, 99.5, 15.5, -0.5], [0.0, 1.0, 0.0, 1.0], [-1.25, 3.75, 2.5, -0.5], [3.0, 3.0, 1.0, 1.0]):
        fig = Figure()
        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
        ax.imshow(np.zeros((3, 4)), extent=extent, aspect="auto", interpolation="nearest")
        out.append({"v0": extent[0], "v1": extent[1], "result": fl(ax.get_xlim()), "imshow": True})
        out.append({"v0": extent[2], "v1": extent[3], "result": fl(ax.get_ylim()), "imshow": True})
    return out


def colorbar_cases():
    out = []
    ranges = [(0.0, 1.0), (0.0, 1e-9), (0.0, 37.0), (0.0, 12.0), (0.0, 3.0), (0.0, 0.0), (5.0, 5.0),
              (0.0, 0.73), (0.0, 250000.0), (1e-6, 4e-6), (-2.0, 1e-9), (-0.5, 4.5), (1000.0, 1000.5),
              (0.0, 2.0), (0.0, 1.5), (0.0, 113.0)]
    for figsize, dpi in (((8.0, 4.5), 120), ((16.0, 9.0), 120), ((4.0, 2.0), 72)):
        for font in (10.0, 14.0, 7.0):
            for vmin, vmax in ranges:
                for integer in (False, True):
                    fig = Figure(figsize=figsize, dpi=dpi)
                    ax = fig.add_axes((0.08, 0.1, 0.6, 0.8))
                    im = ax.imshow(np.zeros((3, 4)), aspect="auto", vmin=vmin, vmax=vmax)
                    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
                    if integer:
                        integer_axis(cb.ax.yaxis)
                    cb.ax.tick_params(labelsize=font * 0.75)
                    axis = cb.ax.yaxis
                    rec = axis_record(cb.ax, axis, fig)
                    rec["ticks"] = fl(cb.get_ticks())
                    out.append({"figsize": fl(figsize), "dpi": dpi, "font_size": font, "vmin": vmin, "vmax": vmax,
                                "integer": integer, "limits": fl(cb.ax.get_ylim()), "axis": rec})
    return out


def thin_tick_cases():
    out = []
    for n in list(range(0, 50)) + [64, 99, 100, 101, 120, 160, 399, 400, 401, 1000, 2047, 5000]:
        for max_ticks in (1, 2, 3, 7, 14, 16, 20, 40, 45, 200):
            out.append({"n": n, "max_ticks": max_ticks, "ticks": [int(i) for i in thin_ticks(n, max_ticks)]})
        for max_ticks in (0, -1, -16):  # improvement: as max_ticks 1
            out.append({"n": n, "max_ticks": max_ticks, "improved": True,
                        "ticks": [int(i) for i in thin_ticks(n, 1)]})
    return out


def write(path, sections):
    parts = []
    for name, cases in sections:
        lines = ",\n".join("  " + json.dumps(c, ensure_ascii=False, sort_keys=True) for c in cases)
        parts.append(f'"{name}": [\n{lines}\n]')
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("{\n" + ",\n".join(parts) + "\n}\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    sections = [
        ("nonsingular", nonsingular_cases()),
        ("locator", locator_cases()),
        ("fixed", fixed_cases()),
        ("formatter", formatter_cases()),
        ("format_data", format_data_cases()),
        ("axes", axes_cases()),
        ("set_lim", set_lim_cases()),
        ("colorbar", colorbar_cases()),
        ("thin_ticks", thin_tick_cases()),
    ]
    write(os.path.join(OUT, "tick_cases.json"), sections)
    print(", ".join(f"{name} {len(cases)}" for name, cases in sections))


if __name__ == "__main__":
    main()
