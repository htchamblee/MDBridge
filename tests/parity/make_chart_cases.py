"""What MDATK's standalone charts (mdatk/render/charts.py render_charts) draw and write, for
MDBridge's charts (src/render charts; docs/specs/m4-plan.md sections 1.6 and 3).

    /opt/mdbridge-venv/bin/python -I -B make_chart_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies; nothing is written into the checkout. MDATK runs in-process, patched by
render_common.setup() as in make_figure_cases.py. Reads the input datasets make_figure_cases.py
wrote (OUT_DIR/inputs: demo_dry, demo_dry_1, demo_dry_nan, demo_dry_times; run that generator first)
and its own charts/inputs/tall (a 700-row matrix, for RN-05's height cap). Writes
OUT_DIR/charts/<case>.json and charts/<case>.npz (the images' arrays).

A case is the demo project (make_figure_cases.demo_base) with `--set`-style overrides of
output.charts / visual / inputs, then render_charts(dataset, project, out_dir, reporter, stamp) into
an empty temporary folder. Every chart figure is recorded where MDATK saves it (Figure.savefig is
wrapped: the figure is drawn, walked by mpl_walker.Walker, then saved as MDATK saves it). Recorded:

  case, description, input, overrides, project (tagged JSON), stamp, improvements
  paths     what render_charts returned, the temporary folder as {out}
  log       the Reporter's log lines ("Wrote N chart file(s) to {out}/out")
  progress  [fraction (float.hex), message] of the Reporter's progress calls ("Chart <key>")
  messages  MDATK's log records, "LEVEL message" (the missing-series and unknown-style warnings)
  error     {type, message} when render_charts raises (e.g. a format Matplotlib does not know);
            the charts saved before stay recorded
  charts    one per figure in the order saved: key, kind ("series" / "matrix", from the
            function that drew it), files [{path, format, size_px or size_pt (PNG: IHDR; JPEG:
            the frame header; TIFF: the first IFD; SVG: width / height; PDF: MediaBox; EPS / PS:
            the BoundingBox's width and height), magic (first bytes, hex)}], figure (the walk, see
            mpl_walker.py; array members "c<k>/<artist id>.<what>" in the .npz)
  files     every file left in the folder (names, sorted)

Image parity is contents only (owner decision D10): file bytes are not recorded.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import render_common as rc  # noqa: E402  (before numpy: pins its kernels)

CHECKOUT, OUT = rc.args("make_chart_cases.py MDATK_CHECKOUT OUT_DIR")
rc.setup(CHECKOUT)
rc.add_mask(OUT, "{data}/render")

import json  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import struct  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

import improvements_render  # noqa: E402
import mpl_walker  # noqa: E402
from mdatk import cli  # noqa: E402
from mdatk.data import AnalysisResult, Dataset  # noqa: E402
from mdatk.demo import demo_project  # noqa: E402
from mdatk.progress import Reporter  # noqa: E402
from mdatk.render import charts as charts_module  # noqa: E402

INPUTS = os.path.join(OUT, "inputs")
CHARTS = os.path.join(OUT, "charts")
CHART_INPUTS = os.path.join(CHARTS, "inputs")  # this generator's own inputs


def demo_base():
    p = demo_project(rc.DEMO_TOP, rc.DEMO_NC, "{work}/out", 12, cleaning=False, mid_time_ps=120.0)
    p.output.out_dir = "{work}/out"
    return p


def s(key, value):
    return f"{key}={json.dumps(value, ensure_ascii=False)}"


LOW = [s("output.charts.dpi", 80)]  # most cases: 640 x 360 px charts
THREE = [s("output.charts.series", ["rmsd/Complex", "hbonds/count", "energy/total"])]

CASES = [
    # (name, input, overrides, stamp, description)
    ("default", "demo_dry", [], "", "every series of the demo (no series option) at 300 dpi, 8 x 4.5 in"),
    ("three", "demo_dry", THREE + LOW, "", "three series: Å, count (integer ticks), energy with its target"),
    ("smoothing", "demo_dry", THREE + LOW + [s("output.charts.smoothing", 5)], "",
     "raw line and running mean (legend)"),
    ("no_target", "demo_dry", THREE + LOW + [s("output.charts.show_target", False)], "", "target line off"),
    ("dark", "demo_dry", THREE + LOW + [s("visual.theme", "dark"), s("output.charts.smoothing", 3)], "",
     "dark theme"),
    ("style_ggplot", "demo_dry", THREE + LOW + [s("output.charts.style", "ggplot")], "",
     "a Matplotlib style: light theme, Matplotlib's chrome"),
    ("style_default", "demo_dry", THREE + LOW + [s("output.charts.style", "default")], "",
     "Matplotlib's default style"),
    ("style_unknown", "demo_dry", THREE + LOW + [s("output.charts.style", "no-such-style")], "",
     "unknown style: warning, MDATK chrome"),
    ("labels", "demo_dry", [s("output.charts.series", ["rmsd/Backbone"])] + LOW + [
        s("output.charts.color", "#e34948"), s("output.charts.line_width", 3.0),
        s("output.charts.title", "{label} of {key}"), s("output.charts.xlabel", "Simulation time"),
        s("output.charts.ylabel", "Deviation")], "", "colour, width, title template, axis labels"),
    ("heatmaps", "demo_dry", [s("output.charts.series", ["rmsd/Complex"]), s("output.charts.include_heatmaps", True)]
     + LOW, "", "matrix charts: counts and Mg occupancy classes (proximity/ excluded)"),
    ("heatmaps_dark_style", "demo_dry", [s("output.charts.series", ["contacts/atom_pairs"]),
                                         s("output.charts.include_heatmaps", True), s("visual.theme", "dark"),
                                         s("output.charts.style", "ggplot")] + LOW, "",
     "matrix charts with a Matplotlib style"),
    ("formats", "demo_dry", [s("output.charts.series", ["rmsd/Complex"]),
                             s("output.charts.formats", ["png", "svg", "PDF", ".png"])] + LOW, "",
     "several formats (case and dot ignored)"),
    ("format_jpg_eps", "demo_dry", [s("output.charts.series", ["rmsd/Complex"]),
                                    s("output.charts.formats", ["jpg", "eps"])] + LOW, "",
     "formats MDBridge does not offer (RN-11)"),
    ("format_unknown", "demo_dry", [s("output.charts.series", ["rmsd/Complex", "hbonds/count"]),
                                    s("output.charts.formats", ["png", "bogus"])] + LOW, "",
     "an unknown format: Matplotlib's error after the first file (RN-11)"),
    ("missing_and_matrix", "demo_dry", [s("output.charts.series", ["nokey/x", "rmsd/Complex", "mg/occupancy",
                                                                   "contacts/map", "rmsd@Nope/x"])] + LOW, "",
     "missing keys warn; matrix and pair-map keys are dropped silently (RN-10)"),
    ("stamped", "demo_dry", THREE[:1] + LOW + [s("output.charts.prefix", "run1"),
                                               s("output.charts.include_heatmaps", True)], "2026-10-15_143205",
     "time-stamped names, own prefix"),
    ("time_frame", "demo_dry", THREE + LOW + [s("inputs.time_unit", "frame")], "", "time axis in frames"),
    ("time_us", "demo_dry", THREE[:1] + LOW + [s("inputs.time_unit", "us")], "", "time axis in µs"),
    ("one_frame", "demo_dry_1", THREE + LOW + [s("output.charts.include_heatmaps", True)], "",
     "one frame: no x limits set; dt = 1 for the matrix"),
    ("small", "demo_dry", THREE[:1] + [s("output.charts.dpi", 50), s("output.charts.width", 4.0),
                                       s("output.charts.height", 3.0)], "", "small chart"),
    ("nan_smoothing", "demo_dry_nan", [s("output.charts.series", ["rmsd/Complex"]), s("output.charts.smoothing", 3)]
     + LOW, "", "a NaN value with smoothing (RN-14)"),
    ("format_tif_ps", "demo_dry", [s("output.charts.series", ["rmsd/Complex"]),
                                   s("output.charts.formats", ["tif", "TIFF", "ps"])] + LOW, "",
     "TIFF and PostScript (pillow / Matplotlib's PS backend in MDATK)"),
    ("format_unknown_two", "demo_dry", [s("output.charts.series", ["rmsd/Complex"]),
                                        s("output.charts.formats", ["gif", "png", ".WebP"])] + LOW, "",
     "two formats MDBridge does not write (RN-11)"),
    ("matrix_series", "demo_dry", [s("output.charts.series", ["mg/occupancy", "rmsd/Complex", "contacts/per_residue_a",
                                                              "mg/occupancy"])] + LOW, "",
     "matrix keys in the series list: matrix charts after the series, each once (RN-10)"),
    ("matrix_series_heatmaps", "demo_dry", [s("output.charts.series", ["hbonds/per_residue_b", "energy/total"]),
                                            s("output.charts.include_heatmaps", True)] + LOW, "",
     "a requested matrix first, then the other heatmaps (RN-10)"),
    ("uneven_times", "demo_dry_times", [s("output.charts.series", ["rmsd/Complex"]),
                                        s("output.charts.include_heatmaps", True)] + LOW, "",
     "matrix cells of two runs with different time steps (RN-01)"),
    ("style_default_smoothing", "demo_dry", THREE + LOW + [s("output.charts.style", "default"),
                                                          s("output.charts.smoothing", 4),
                                                          s("output.charts.include_heatmaps", True)], "",
     "Matplotlib's default look: colour cycle, legend, target, matrix charts"),
    ("style_default_color", "demo_dry", [s("output.charts.series", ["energy/total"]), s("output.charts.style", "default"),
                                         s("output.charts.color", "tab:green"), s("output.charts.smoothing", 3)]
     + LOW, "", "Matplotlib's default look with a colour"),
    ("dark_heatmaps", "demo_dry", [s("output.charts.series", ["rmsd/Complex"]), s("output.charts.include_heatmaps", True),
                                   s("visual.theme", "dark")] + LOW, "",
     "dark theme matrix charts: title and colour bar in the theme's colours (RN-50)"),
    ("matrix_tall", "tall", [s("output.charts.series", ["tall/occupancy"]), s("output.charts.width", 1.0)], "",
     "a 700-row matrix at 300 dpi: 38 250 px tall, capped at 32 768 px (RN-05)"),
]


def tall_input():
    """demo_dry's frames with one synthetic 700 x 24 count matrix (labels R1 ... R700)."""
    demo = rc.load_dataset(os.path.join(INPUTS, "demo_dry"))
    ds = Dataset(demo.frames, demo.time_ps, time_unit=demo.time_unit)
    rows, frames = 700, len(demo.frames)
    r = np.arange(rows)[:, None]
    f = np.arange(frames)[None, :]
    res = AnalysisResult("tall", "synthetic", {}, {"occupancy": ((r * 7 + f * 3) % 5).astype(np.float64)},
                         {"rows": [f"R{k + 1}" for k in range(rows)]})
    res.add_matrix("occupancy", "occupancy", "rows", "Tall occupancy", kind="count")
    ds.add(res)
    folder = os.path.join(CHART_INPUTS, "tall")
    rc.save_dataset(ds, folder)
    return rc.load_dataset(folder)


# --- recording ---------------------------------------------------------------------------------------

def file_info(path, folder):
    with open(path, "rb") as fh:
        data = fh.read()
    info = {"path": rc.mask_text(path.replace(folder, "{out}")), "format": os.path.splitext(path)[1][1:],
            "magic": data[:8].hex(), "bytes_nonzero": len(data) > 0}
    if data[:8] == b"\x89PNG\r\n\x1a\n":
        info["size_px"] = list(struct.unpack(">II", data[16:24]))
    elif b"<svg" in data[:2000]:
        m = re.search(rb'<svg[^>]*width="([0-9.]+)pt"[^>]*height="([0-9.]+)pt"', data)
        info["size_pt"] = [float(m.group(1)), float(m.group(2))] if m else None
    elif data[:5] == b"%PDF-":
        m = re.search(rb"/MediaBox \[ *([0-9. ]+)\]", data)
        info["size_pt"] = [float(v) for v in m.group(1).split()] if m else None
    elif data[:2] == b"%!":
        m = re.search(rb"%%BoundingBox: *(-?[0-9]+) +(-?[0-9]+) +(-?[0-9]+) +(-?[0-9]+)", data)
        info["size_pt"] = [int(m.group(3)) - int(m.group(1)), int(m.group(4)) - int(m.group(2))] if m else None
    elif data[:2] == b"\xff\xd8":
        info["size_px"] = jpeg_size(data)
    elif data[:4] in (b"II*\x00", b"MM\x00*"):
        info["size_px"] = tiff_size(data)
    return info


def jpeg_size(data):
    """[width, height] from a JPEG's first SOF marker."""
    at = 2
    while at + 9 < len(data) and data[at] == 0xFF:
        marker, length = data[at + 1], struct.unpack(">H", data[at + 2:at + 4])[0]
        if 0xC0 <= marker <= 0xCF and marker not in (0xC4, 0xC8, 0xCC):
            h, w = struct.unpack(">HH", data[at + 5:at + 9])
            return [w, h]
        at += 2 + length
    return None


def tiff_size(data):
    """[width, height] from a TIFF's first IFD (tags 256 / 257)."""
    e = "<" if data[:2] == b"II" else ">"
    ifd = struct.unpack(e + "I", data[4:8])[0]
    n = struct.unpack(e + "H", data[ifd:ifd + 2])[0]
    size = {}
    for k in range(n):
        tag, typ, _count = struct.unpack(e + "HHI", data[ifd + 2 + 12 * k:ifd + 10 + 12 * k])
        raw = data[ifd + 10 + 12 * k:ifd + 14 + 12 * k]
        value = struct.unpack(e + "H", raw[:2])[0] if typ == 3 else struct.unpack(e + "I", raw)[0]
        if tag in (256, 257):
            size[tag] = value
    return [size.get(256), size.get(257)]


class Recorder:
    def __init__(self, folder):
        self.folder = folder
        self.charts = []
        self.arrays = {}
        self.current = {}

    def savefig(self, fig, original, path, *a, **kw):
        k = len(self.charts)
        if self.current.get("fig") is not fig:
            fig.canvas.draw()
            walker = mpl_walker.Walker(fig, prefix=f"c{k}/")
            self.current = {"fig": fig, "record": {"key": None, "kind": None, "files": [],
                                                   "figure": walker.figure()}}
            self.arrays.update(walker.arrays)
            self.charts.append(self.current["record"])
        try:
            return original(fig, path, *a, **kw)
        finally:
            if os.path.exists(path):
                self.current["record"]["files"].append(file_info(path, self.folder))


def run_case(name, inp, overrides, stamp, description, inputs):
    p = demo_base()
    for assignment in overrides:
        cli.apply_override(p, assignment)
    ds = inputs[inp] if inp in inputs else rc.load_dataset(os.path.join(INPUTS, inp))
    folder = tempfile.mkdtemp(prefix="mdbridge-charts-")
    out_dir = os.path.join(folder, "out")
    rec = {"case": name, "description": description, "input": inp, "overrides": overrides,
           "project": rc.project_dict(p), "stamp": stamp, "improvements": improvements_render.applied()}
    recorder = Recorder(folder)
    original = Figure.savefig

    def savefig(fig, path, *a, **kw):
        return recorder.savefig(fig, original, path, *a, **kw)

    series_fn, matrix_fn = charts_module.render_series_chart, charts_module.render_matrix_chart

    def tagged(fn, kind):
        def run(dataset, project, key, *a, **kw):
            start = len(recorder.charts)
            try:
                return fn(dataset, project, key, *a, **kw)
            finally:
                for c in recorder.charts[start:]:
                    c["key"], c["kind"] = key, kind
        return run

    logs, progress = [], []
    reporter = Reporter(log_cb=logs.append, progress_cb=lambda f, m: progress.append([float(f).hex(), m]))
    rc.reset_log()
    Figure.savefig = savefig
    charts_module.render_series_chart = tagged(series_fn, "series")
    charts_module.render_matrix_chart = tagged(matrix_fn, "matrix")
    try:
        paths = charts_module.render_charts(ds, p, out_dir, reporter, stamp)
        rec["paths"] = [path.replace(folder, "{out}") for path in paths]
    except Exception as exc:  # noqa: BLE001 - MDATK's error is the expected behaviour
        rec["error"] = {"type": type(exc).__name__, "message": rc.mask_text(str(exc).replace(folder, "{out}"))}
    finally:
        Figure.savefig = original
        charts_module.render_series_chart, charts_module.render_matrix_chart = series_fn, matrix_fn
    rec["log"] = [rc.mask_text(line.replace(folder, "{out}")) for line in logs]
    rec["progress"] = progress
    rec["messages"] = [rc.mask_text(line.replace(folder, "{out}")) for line in rc.LOG_LINES]
    rec["charts"] = recorder.charts
    rec["files"] = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
    shutil.rmtree(folder, ignore_errors=True)
    return improvements_render.adjust("chart", name, rec), recorder.arrays


def main():
    if os.path.isdir(CHARTS):
        shutil.rmtree(CHARTS)
    os.makedirs(CHARTS)
    inputs = {"tall": tall_input()}
    for case in CASES:
        rec, arrays = run_case(*case, inputs)
        rc.write_json(os.path.join(CHARTS, f"{case[0]}.json"), rec)
        if arrays:
            rc.write_npz(os.path.join(CHARTS, f"{case[0]}.npz"), arrays)
        print(case[0], len(rec["charts"]), rec.get("error", ""), flush=True)


if __name__ == "__main__":
    main()
