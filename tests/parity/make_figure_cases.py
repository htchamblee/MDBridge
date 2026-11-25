"""What MDATK's multi-panel figure (mdatk/render/video.py build_figure, panels.py, seqpanels.py)
shows, for MDBridge's figure and panels (src/render; docs/specs/m4-plan.md section 3).

    /opt/mdbridge-venv/bin/python -I -B make_figure_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies (numpy, scipy, matplotlib 3.11.2, MDAnalysis 2.11.0); nothing is written
into the checkout. MDATK runs in-process, patched by render_common.setup() (MDBridge's independent
tables and analysis improvements, as make_pipeline_cases.py; MDBridge's rendering improvements,
improvements_render.py), with numpy's and OpenBLAS' kernels pinned. Reads
tests/data/topology/demo_dry.parm7, tests/data/io/trajectories/demo_dry.nc and
tests/data/analysis/systems/aromatic.pdb. Writes into OUT_DIR (tests/data/render):

  inputs/<name>/    the datasets the figures (and make_chart_cases.py's charts) are built from
                    (render_common.save_dataset: dataset.json + <result>.npz), read back before
                    use, so the C++ tests start from exactly these bytes:
                      demo_dry      the demo project's analyses on demo_dry: every analysis plus
                                    trace / atoms / proximity / basepairs (24 frames, 3 Mg ions)
                      demo_dry_1    the same on the first frame only
                      demo_dry_nan  demo_dry's trace, proximity and rmsd with a non-finite trace
                                    point (frame 3) and RMSD value (frame 2): MDATK's crash and NaN
                                    windows (RN-03, RN-14)
                      protein_rna   trace, atoms and proximity of aromatic.pdb (RNA chain A,
                                    protein chain B; 3 models)
                      three_chains  a synthetic trace of three chains (A and B 30 nucleotides, LIG
                                    9; 40 frames) with base pairs within A and between A and B
                                    that change every 10 frames, and proximity
                      example_<case> MDATK's example dataset (mdatk/render/example.py) of the case's
                                    project: the GUI's layout preview (120 frames; 24 for
                                    example_rep, whose atoms are repeated per frame)
  inputs/demo_dry_map.mrc   a density map simulated from demo_dry's first frame (atoms, unit
                    weights, resolution 6 Å; mdatk.density.simulate_map / mdatk.io.mrc.write_map),
                    placed with the identity transform
  figures/<case>.json  one figure per case (CASES below), figures/<case>.npz its arrays

A case is a project (the input's base project with `--set`-style overrides,
mdatk.cli.apply_override), a size, dpi and watermark: build_figure(dataset, project, w, h, dpi,
watermark), then, for each position (0, the middle and the last analysed frame), every panel's
update(pos) and canvas.draw(). Recorded (render_common.dumps; floats as repr):

  case, description, input, overrides, project (Project.to_dict(), tagged JSON, paths {data}),
  size [w, h] (as passed; the canvas is figure_size(w, h)), dpi, watermark, positions,
  improvements (improvements_render.applied())
  error     {stage ("build" or "pos <n>"), type, message} when MDATK raises; what was recorded
            before stays
  log       "LEVEL message" of MDATK's log records while the case ran (build_figure's warnings,
            the panels' warn_once, the 3D layer's notes)
  panels    per panel in slot order: type, ok, auto_title, axes (ids of its main / twin /
            colour-bar axes), animated (ids of animated_artists(): what the video loop redraws per
            frame) and the panel attributes it has of: radius, data_lo, data_hi (the structure
            window before the aspect fit), view_note, stage_note, density_note, stage {renderer,
            quality, failed}; per_row, n_rows, cell_size (sequence); source, span (secondary)
  frames    [{pos, figure}]: mpl_walker.Walker(fig).figure() (see mpl_walker.py) at the first
            position; at the later ones {pos, artists}: the animated artists only (the rest of the
            figure does not change between frames; the video loop draws it once), each as its id
            and the fields whose value differs from the first position's record of that artist

The .npz holds the walker's arrays, named "p<pos>/<artist id>.<what>": images' data (float64, NaN
where masked), the RGBA they show (uint8) and colour-bar lookup tables (float64 RGBA).

Image parity is contents only (owner decision D10): there are no reference images. The 3D layer's
pixels are make_molecule_cases.py's subject (fixed viewports); here only its size and coverage.
The output is reproducible (numpy / OpenBLAS kernels pinned by render_common, matplotlib's bundled
DejaVu fonts, inputs rounded as render_common.save_dataset describes).
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import render_common as rc  # noqa: E402  (before numpy: pins its kernels)

CHECKOUT, OUT = rc.args("make_figure_cases.py MDATK_CHECKOUT OUT_DIR")
rc.setup(CHECKOUT)
rc.add_mask(OUT, "{data}/render")

import copy  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402

import improvements_render  # noqa: E402
import mpl_walker  # noqa: E402
from mdatk import cli  # noqa: E402
from mdatk.config import Project  # noqa: E402
from mdatk.data import AnalysisResult, CSRBuilder, Dataset  # noqa: E402
from mdatk.demo import demo_project  # noqa: E402
from mdatk.pipeline import Pipeline  # noqa: E402
from mdatk.progress import Reporter  # noqa: E402
from mdatk.render import residues as residues_module  # noqa: E402
from mdatk.render.example import example_dataset  # noqa: E402
from mdatk.render.video import build_figure  # noqa: E402

INPUTS = os.path.join(OUT, "inputs")
FIGURES = os.path.join(OUT, "figures")
MAP = os.path.join(INPUTS, "demo_dry_map.mrc")
AROMATIC = os.path.join(rc.SYSTEMS, "aromatic.pdb")
ALL_ANALYSES = ["trace", "atoms", "proximity", "basepairs", "rmsd", "contacts", "hbonds", "energy", "mg",
                "stacking"]
WATERMARK = "Example data – layout preview"


# --- base projects ----------------------------------------------------------------------------------

def demo_base():
    """The demo project as `mdatk demo` writes it for demo_dry (its layout: large_plus_4)."""
    p = demo_project(rc.DEMO_TOP, rc.DEMO_NC, "{work}/out", 12, cleaning=False, mid_time_ps=120.0)
    p.output.out_dir = "{work}/out"
    return p


def protein_base():
    p = Project()
    p.inputs.topology = AROMATIC
    p.inputs.trajectories = []
    p.inputs.time_unit = "ps"
    p.inputs.dt_ps = 2.0
    p.groups[0].name, p.groups[0].residues = "RNA", "1-4"
    p.groups[1].name, p.groups[1].residues = "Protein", "5-10"
    p.visual.proximity_set = "Protein"
    p.output.out_dir = "{work}/out"
    p.layout.preset = "single"
    p.layout.panels = []
    return p


def synthetic_base():
    p = Project()
    p.inputs.time_unit = "ns"
    p.groups[0].name, p.groups[0].residues = "Strand A", "A:1-30"
    p.groups[1].name, p.groups[1].residues = "Strand B", "B:31-60"
    p.output.out_dir = "{work}/out"
    p.visual.proximity_set = "Strand B"
    p.layout.preset = "single"
    p.layout.panels = []
    return p


BASES = {"demo_dry": demo_base, "demo_dry_1": demo_base, "demo_dry_nan": demo_base, "protein_rna": protein_base,
         "three_chains": synthetic_base, "example": demo_base}


# --- input datasets ---------------------------------------------------------------------------------

def analysed(project, names):
    work = tempfile.mkdtemp(prefix="mdbridge-render-inputs-")
    try:
        p = copy.deepcopy(project)
        p.output.out_dir = os.path.join(work, "out")
        return Pipeline(p, Reporter()).analyze(names)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def three_chains():
    """Chains A (resid 1-30) and B (31-60) as a duplex and a 9-residue chain LIG (101-109), 40 frames,
    10 ps apart; base pairs A-B (fraying in frames 10-19 and 30-39) and a hairpin within A's first
    ten residues from frame 20; every residue drifts in and out of proximity."""
    n_frames, dt = 40, 10.0
    coords, chain, names, gidx, resid, letters, links = [], [], [], [], [], [], []
    seq = "GCAUGGCAUCGAGCAUGGCAUCGAGCAUGG"
    comp = {"G": "C", "C": "G", "A": "U", "U": "A"}
    strands = [("A", seq, list(range(1, 31)), 0),
               ("B", "".join(comp[b] for b in reversed(seq)), list(range(31, 61)), 1),
               ("LIG", "AAGGCCUUA", list(range(101, 110)), -1)]
    for s, (name, letters_s, ids, group) in enumerate(strands):
        for k, r in enumerate(ids):
            level = k if s != 1 else len(ids) - 1 - k
            ang = np.radians(level * 32.7 + (150 if s == 1 else 0))
            if s < 2:
                xyz = [9.0 * np.cos(ang), 9.0 * np.sin(ang), 2.8 * level]
            else:
                xyz = [25.0 + 3.0 * k, 4.0 * np.sin(k / 2.0), 10.0]
            coords.append(xyz)
            chain.append(s)
            names.append(name)
            gidx.append(group)
            resid.append(r)
            letters.append(letters_s[k])
            if k:
                links.append([len(coords) - 2, len(coords) - 1])
    xyz = np.asarray(coords, dtype=np.float64)
    wobble = 0.6 * np.sin(np.arange(n_frames)[:, None, None] / 5.0 + np.arange(len(xyz))[None, :, None] / 7.0)
    frames_xyz = (xyz[None] + wobble).astype(np.float32)
    n_res = len(xyz)
    member = np.zeros((2, n_res), bool)
    member[0, :30], member[1, 30:60] = True, True
    trace = AnalysisResult("trace", "synthetic", {}, {}, {})
    trace.arrays.update(coords=frames_xyz, ion_coords=np.zeros((n_frames, 0, 3), np.float32),
                        chain=np.asarray(chain, np.int64), group_index=np.asarray(gidx, np.int64),
                        group_member=member, is_point=np.ones(n_res, bool),
                        links=np.asarray(links, np.int64).reshape(-1, 2), resindex=np.arange(n_res, dtype=np.int64),
                        resid=np.asarray(resid, np.int64), resname=np.asarray(letters, dtype="U8"),
                        chain_name=np.asarray(names, dtype="U3"))
    b = CSRBuilder()
    for f in range(n_frames):
        fray = 3 if (f // 10) % 2 else 0
        i = np.arange(fray, 30)
        j = 59 - i
        if f >= 20:  # a hairpin in A replaces A's first pairs
            i = np.concatenate([np.arange(0, 4), i[i >= 10]])
            j = np.concatenate([np.arange(9, 5, -1), 59 - i[4:]])
        b.add(i, j, np.ones(len(i)))
    bp = AnalysisResult("basepairs", "synthetic", {}, {}, {})
    bp.arrays.update(b.arrays("map"))
    labels = [str(k) for k in range(n_res)]
    bp.meta["map_r"], bp.meta["map_c"] = labels, labels
    bp.add_pairmap("map", "map", "map_r", "map_c", "Base pairs")
    k = np.arange(n_res)[:, None]
    dist = 5.0 + 4.0 * np.sin(k / 3.0 + np.arange(n_frames)[None, :] / 6.0)
    prox = AnalysisResult("proximity", "synthetic", {}, {}, {})
    prox.arrays.update(min_distance=dist.astype(np.float32), member=member[1].copy(),
                       resindex=np.arange(n_res, dtype=np.int64))
    ds = Dataset(np.arange(n_frames), np.arange(n_frames) * dt, time_unit="ns")
    for r in (trace, bp, prox):
        ds.add(r)
    return ds


def with_nan(ds):
    """demo_dry's trace, proximity and rmsd results with trace point 5 of frame 3 and the first RMSD
    series' frame 2 not finite."""
    out = Dataset(ds.frames, ds.time_ps, time_unit=ds.time_unit)
    for name in ("trace", "proximity", "rmsd"):
        out.add(copy.deepcopy(ds.results[name]))
    out.results["trace"].arrays["coords"][3, 5] = np.nan
    out.results["rmsd"].arrays["rmsd_0"][2] = np.nan
    return out


def write_inputs():
    """Every input dataset (and the map), written and read back."""
    demo = demo_base()
    datasets = {}
    datasets["demo_dry"] = analysed(demo, ALL_ANALYSES)
    one = copy.deepcopy(demo)
    one.inputs.stop = 1
    datasets["demo_dry_1"] = analysed(one, ALL_ANALYSES)
    datasets["protein_rna"] = analysed(protein_base(), ["trace", "atoms", "proximity"])
    datasets["three_chains"] = three_chains()
    for name, ds in datasets.items():
        rc.save_dataset(ds, os.path.join(INPUTS, name))
    loaded = {name: rc.load_dataset(os.path.join(INPUTS, name)) for name in datasets}
    rc.save_dataset(with_nan(loaded["demo_dry"]), os.path.join(INPUTS, "demo_dry_nan"))
    loaded["demo_dry_nan"] = rc.load_dataset(os.path.join(INPUTS, "demo_dry_nan"))

    from mdatk.density import simulate_map
    from mdatk.io.mrc import write_map

    xyz = np.asarray(loaded["demo_dry"].results["atoms"].arrays["coords"][0], dtype=np.float64)
    write_map(MAP, simulate_map(xyz, np.ones(len(xyz)), 6.0))
    return loaded


# --- cases ------------------------------------------------------------------------------------------

def P(type_, series=(), title="", **options):
    return {"type": type_, "series": list(series), "title": title, "options": options}


def layout(preset, *panels):
    return [f"layout.preset={json.dumps(preset)}", f"layout.panels={json.dumps(list(panels), ensure_ascii=False)}"]


def s(key, value):
    return f"{key}={json.dumps(value, ensure_ascii=False)}"


S = P("structure")
RMSD = P("timeseries", ["rmsd/Complex"])
HB = P("timeseries", ["hbonds/count"])
EN = P("timeseries", ["energy/total"])
CM = P("contact_map", ["contacts/map"])
HM = P("heatmap", ["contacts/per_residue_a"])
MG = P("heatmap", ["mg/occupancy"])
DA = P("dual_axis", ["rmsd/Complex", "energy/total"])
SEQ = P("sequence")
SEC = P("secondary")
PRESETS = ["single", "large_plus_2", "large_plus_3", "large_plus_4", "large_top_2", "large_top_3", "equal_2",
           "equal_3", "equal_4", "equal_5", "equal_6", "single_strip", "large_plus_2_strip", "large_plus_3_strip",
           "large_top_3_strip", "equal_4_strip"]
SLOTS = {"single": 1, "large_plus_2": 3, "large_plus_3": 4, "large_plus_4": 5, "large_top_2": 3, "large_top_3": 4,
         "equal_2": 2, "equal_3": 3, "equal_4": 4, "equal_5": 5, "equal_6": 6, "single_strip": 2,
         "large_plus_2_strip": 4, "large_plus_3_strip": 5, "large_top_3_strip": 5, "equal_4_strip": 5}
FILL = [S, RMSD, CM, HM, DA, EN]
NINE = ["rmsd/Complex", "rmsd/Backbone", "energy/total", "energy/elec", "energy/vdw", "hbonds/count",
        "contacts/atom_pairs", "contacts/residue_pairs", "basepairs/count"]
DOTBRACKET = "(" * 12 + ")" * 12
SMALL = {"size": (960, 540), "dpi": 60}


def preset_cases():
    out = []
    for k, name in enumerate(PRESETS):
        n = SLOTS[name]
        panels = [FILL[(k + i) % len(FILL)] if i else S for i in range(n)]
        if name.endswith("_strip"):
            panels[-1] = SEQ if k % 2 else SEC
        out.append((f"preset_{name}", "demo_dry", layout(name, *panels), {}, f"layout preset {name}"))
    return out


def rep(*specs):
    return s("visual.representations", list(specs))


CASES = [
    # (name, input, overrides, options (size, dpi, watermark), description)
    ("demo_video", "demo_dry", [], {"size": (1920, 1080), "dpi": 120},
     "the demo project at the video's default size (1920 x 1080, dpi 120)"),
    ("demo", "demo_dry", [], {}, "the demo project (large_plus_4) at 960 x 540, dpi 60"),
    ("demo_dark", "demo_dry", [s("visual.theme", "dark")], {}, "dark theme"),
    ("demo_title_watermark", "demo_dry", [s("visual.title", "Aptamer and target, 24 frames")],
     {"watermark": WATERMARK}, "figure title (1.3 x font size) and the GUI's watermark"),
    ("demo_font_16", "demo_dry", [s("visual.font_size", 16.0)], {}, "larger base font"),
    ("demo_odd_size", "demo_dry", [], {"size": (1001, 563), "dpi": 77}, "odd size: rounded down to even pixels"),
    ("demo_small_dpi", "demo_dry", [], {"size": (640, 360), "dpi": 40}, "small figure"),
    ("extra_panels", "demo_dry", layout("single", S, RMSD, HB), {}, "more panels than slots: warning, extras ignored"),
    *preset_cases(),
    # time series
    ("ts_two_same_unit", "demo_dry", layout("single", P("timeseries", ["rmsd/Complex", "rmsd/Backbone"])), {},
     "two series in Å: legend, y label = unit"),
    ("ts_mixed_units", "demo_dry", layout("single", P("timeseries", ["rmsd/Complex", "hbonds/count", "energy/total"])),
     {}, "three units: no y label; energy's target line"),
    ("ts_nine_series", "demo_dry", layout("single", P("timeseries", NINE)), {}, "nine series: the palette repeats"),
    ("ts_progressive", "demo_dry", layout("single", P("timeseries", ["rmsd/Complex"], progressive=True)), {},
     "progressive lines (animated)"),
    ("ts_smooth", "demo_dry", layout("single", P("timeseries", ["rmsd/Complex"], smooth=5)), {}, "running mean of 5"),
    ("ts_target_ylim", "demo_dry", layout("single", P("timeseries", ["rmsd/Complex"], "RMSD (Å)", target=1.5,
                                                      ylim=[0.0, 4.0])), {}, "target option, ylim, own title"),
    ("ts_energy_target", "demo_dry", layout("single", EN), {}, "the series' own target (energy.target_value)"),
    ("ts_counts", "demo_dry", layout("single", HB), {}, "counts: integer ticks"),
    ("ts_no_series", "demo_dry", layout("single", P("timeseries")), {}, "No data: no series selected"),
    ("ts_missing_key", "demo_dry", layout("single", P("timeseries", ["nokey/x"])), {}, "No data: unknown key"),
    ("ts_time_ps", "demo_dry", layout("single", RMSD) + [s("inputs.time_unit", "ps")], {}, "time in ps"),
    ("ts_time_frame", "demo_dry", layout("single", RMSD) + [s("inputs.time_unit", "frame")], {}, "time in frames"),
    ("ts_nan_smooth", "demo_dry_nan", layout("single", P("timeseries", ["rmsd/Complex"], smooth=3)), {},
     "a NaN value with smoothing blanks whole windows (RN-14)"),
    # dual axis
    ("dual_axis", "demo_dry", layout("single", DA), {}, "RMSD (left) and energy (right)"),
    ("dual_axis_counts", "demo_dry", layout("single", P("dual_axis", ["hbonds/count", "basepairs/count"])), {},
     "two count axes"),
    ("dual_axis_one", "demo_dry", layout("single", P("dual_axis", ["rmsd/Complex"])), {},
     "No data: needs two series"),
    # heatmaps
    ("heatmap_counts", "demo_dry", layout("single", HM), {}, "per-residue contacts, colour bar 'count'"),
    ("heatmap_classes", "demo_dry", layout("single", MG), {}, "Mg occupancy classes (ordinal map)"),
    ("heatmap_values", "demo_dry", layout("single", P("heatmap", ["proximity/min_distance"])), {},
     "distances (value matrix)"),
    ("heatmap_options", "demo_dry", layout("single", P("heatmap", ["hbonds/per_residue_b"], vmin=1.0, vmax=3.0,
                                                       max_ticks=4)), {}, "vmin, vmax, max_ticks"),
    ("heatmap_vmin_above_vmax", "demo_dry", layout("single", P("heatmap", ["contacts/per_residue_a"], vmin=5.0,
                                                               vmax=1.0)), {}, "vmin > vmax (RN-04)"),
    ("heatmap_dark", "demo_dry", layout("equal_2", HM, MG) + [s("visual.theme", "dark")], {}, "dark ramps"),
    ("heatmap_one_frame", "demo_dry_1", layout("equal_2", HM, RMSD), {}, "one frame: dt = 1, xlim t0..t0+1"),
    # contact maps
    ("contact_instant", "demo_dry", layout("single", CM), {}, "instant map (animated)"),
    ("contact_frame", "demo_dry", layout("single", P("contact_map", ["contacts/map"], mode="frame", frame=5)), {},
     "fixed frame 5"),
    ("contact_frame_clamped", "demo_dry", layout("single", P("contact_map", ["hbonds/map"], mode="frame",
                                                             frame=999)), {}, "frame beyond the end: clamped"),
    ("contact_frequency", "demo_dry", layout("single", P("contact_map", ["contacts/map"], mode="frequency")), {},
     "fraction of frames"),
    ("contact_transpose", "demo_dry", layout("single", P("contact_map", ["stacking/map"], transpose=True,
                                                         max_ticks=5)), {}, "transposed, 5 ticks"),
    # structure view
    ("structure_groups", "demo_dry", layout("single", S) + [s("visual.glow_mode", "none")], {},
     "group colours, no glow"),
    ("structure_chain", "demo_dry", layout("single", S) + [s("visual.color_mode", "chain")], {},
     "chain colours (tab10), legend Chain 1 / 2 (RN-09)"),
    ("structure_chain_palette", "demo_dry", layout("single", S) + [s("visual.color_mode", "chain"),
                                                                    s("visual.chain_palette", "viridis")], {},
     "chain colours from a continuous map"),
    ("structure_proximity", "demo_dry", layout("single", S) + [s("visual.color_mode", "proximity"),
                                                                s("visual.glow_mode", "current")], {},
     "proximity colour and current glow"),
    ("structure_glow_ever", "demo_dry", layout("single", S) + [s("visual.glow_mode", "ever")], {}, "glow ever"),
    ("structure_glow_fade", "demo_dry", layout("single", S) + [s("visual.glow_mode", "fade"),
                                                                s("visual.fade_frames", 3),
                                                                s("visual.glow_strength", 0.5)], {},
     "fade over 3 frames"),
    ("structure_zoom", "demo_dry", layout("single", S) + [s("visual.zoom", 2.5)], {}, "zoom 2.5"),
    ("structure_zoom_tiny", "demo_dry", layout("single", S) + [s("visual.zoom", 0.01)], {}, "zoom below 0.05"),
    ("structure_centre_residue", "demo_dry", layout("single", S) + [s("visual.center_on", "B:20")], {},
     "centred on B:20"),
    ("structure_centre_origin", "demo_dry", layout("single", S) + [s("visual.center_on", "origin")], {},
     "centred on the origin"),
    ("structure_centre_unknown", "demo_dry", layout("single", S) + [s("visual.center_on", "Z:99")], {},
     "a centre that matches nothing: warning, model centre"),
    ("structure_centre_follow", "demo_dry", layout("single", S) + [s("visual.center_on", "B:20"),
                                                                    s("visual.center_follow", True)], {},
     "followed centre"),
    ("structure_tags", "demo_dry", layout("single", S) + [s("visual.label_every", 4),
                                                           s("visual.label_residues", "1, 12-13")], {},
     "residue tags"),
    ("structure_plain", "demo_dry", layout("single", S) + [s("visual.show_time_label", False),
                                                            s("visual.depth_cue", False),
                                                            s("visual.rotation", [30.0, 45.0, 60.0]),
                                                            s("visual.point_size", 40.0), s("visual.line_width", 4.0)],
     {}, "no time label, no depth cue, rotated, larger points"),
    ("structure_time_us", "demo_dry", layout("single", S) + [s("inputs.time_unit", "us")], {}, "time label in µs"),
    ("structure_nan", "demo_dry_nan", layout("single", S), {}, "a non-finite trace point (RN-03)"),
    ("structure_one_frame", "demo_dry_1", layout("equal_2", S, CM), {}, "one frame"),
    # 3D representations (software renderer)
    ("rep_cartoon", "demo_dry", layout("equal_2", S, RMSD) + [rep({"style": "cartoon"})], {}, "cartoon"),
    ("rep_cartoon_surf", "demo_dry", layout("equal_2", S, RMSD) + [
        rep({"style": "cartoon", "groups": ["Aptamer"], "sugars": True},
            {"style": "surf", "groups": ["Target"], "opacity": 0.5})], {},
     "cartoon with sugars and a translucent surface"),
    ("rep_licorice_hide", "demo_dry", layout("equal_2", S, RMSD) + [
        rep({"style": "licorice", "groups": ["Target"], "color": "element"}), s("visual.rep_rest", "hide")], {},
     "licorice, the rest hidden"),
    ("rep_vdw_tube", "demo_dry", layout("equal_2", S, RMSD) + [
        rep({"style": "vdw", "groups": ["Aptamer"], "size": 0.6}, {"style": "tube", "groups": ["Target"],
                                                                   "color": "nucleotide"})], {}, "vdw and tube"),
    ("rep_ribbon_glow", "demo_dry", layout("equal_2", S, RMSD) + [
        rep({"style": "ribbon"}), s("visual.glow_mode", "current")], {}, "ribbon with glow"),
    ("rep_trace_vtk", "demo_dry", layout("equal_2", S, RMSD) + [
        rep({"style": "trace"}), s("visual.renderer", "vtk")], {}, "renderer vtk without VTK: note, built-in renderer"),
    # density map shadow
    ("density_behind", "demo_dry", layout("equal_2", S, RMSD) + [
        s("density.map_path", MAP), s("density.transform", [1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0, 0, 0, 0]),
        s("density.level", 1.0)], {}, "density shadow behind the trace"),
    ("density_front_rep", "demo_dry", layout("equal_2", S, RMSD) + [
        s("density.map_path", MAP), s("density.transform", [1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0, 0, 0, 0]),
        s("density.front", True), s("density.opacity", 0.8), s("density.color", "#e34948"),
        rep({"style": "cartoon"})], {}, "density shadow in front of a cartoon (animated)"),
    ("density_follow", "demo_dry", layout("single", S) + [
        s("density.map_path", MAP), s("density.transform", [1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0, 0, 0, 0]),
        s("visual.center_on", "A:3"), s("visual.center_follow", True)], {}, "the shadow follows a followed centre"),
    ("density_not_fitted", "demo_dry", layout("single", S) + [s("density.map_path", MAP),
                                                               s("density.auto_fit", False)], {},
     "a map without placement: note"),
    # sequence and secondary structure
    ("sequence", "demo_dry", layout("single_strip", S, SEQ), {}, "sequence strip"),
    ("sequence_wrap", "demo_dry", layout("single_strip", S, P("sequence", per_row=5)), {},
     "5 per row: 5 rows, note and warning"),
    ("sequence_range_tags", "demo_dry", layout("single_strip", S, P("sequence", residues="B:13-18")) + [
        s("visual.label_every", 2), s("visual.glow_mode", "current")], {}, "a range, tags, glow"),
    ("sequence_bad_range", "demo_dry", layout("single_strip", S, P("sequence", residues="Z:1")), {},
     "No data: range matches nothing"),
    ("secondary_dynamic", "demo_dry", layout("single_strip", S, SEC), {}, "pairs per frame"),
    ("secondary_first", "demo_dry", layout("single_strip", S, P("secondary", source="first",
                                                                 residues="A:1-6, B:19-24")), {},
     "first frame's pairs of a range"),
    ("secondary_dotbracket", "demo_dry", layout("single_strip", S, P("secondary", source="dotbracket",
                                                                      dotbracket=DOTBRACKET)) +
     [s("visual.label_every", 6)], {}, "pairs from a dot-bracket"),
    ("secondary_dotbracket_length", "demo_dry", layout("single_strip", S, P("secondary", source="dotbracket",
                                                                             dotbracket="((..))")), {},
     "No data: dot-bracket of the wrong length"),
    ("secondary_dotbracket_empty", "demo_dry", layout("single_strip", S, P("secondary", source="dotbracket")), {},
     "No data: no dot-bracket"),
    # three chains (synthetic)
    ("chains_sequence", "three_chains", layout("single_strip", S, SEQ) + [s("visual.color_mode", "chain")], {},
     "three chains: chain labels in the gaps, wrapping, warning"),
    ("chains_secondary", "three_chains", layout("single_strip", S, SEC) + [s("visual.glow_mode", "ever")], {},
     "intra-chain arcs above, inter-chain arcs below; pairs change every 10 frames"),
    ("chains_sequence_range", "three_chains", layout("equal_4_strip", S, P("secondary", residues="A:1-12"),
                                                     P("sequence", residues="A:1-10, LIG:101-109"),
                                                     P("secondary", residues="LIG:101-109"),
                                                     P("sequence", per_row=12)), {}, "ranges across chains"),
    # protein and RNA
    ("protein_trace", "protein_rna", layout("single_strip", S, SEQ) + [s("visual.color_mode", "chain")], {},
     "protein CA trace next to RNA; amino-acid letters"),
    ("protein_tube", "protein_rna", layout("single", S) + [rep({"style": "tube"}), s("visual.glow_mode", "none")],
     {}, "protein tube and RNA tube"),
    ("protein_cartoon_licorice", "protein_rna", layout("single", S) + [
        rep({"style": "cartoon", "groups": ["RNA"]}, {"style": "licorice", "groups": ["Protein"]})], {},
     "RNA cartoon, protein licorice"),
    # the GUI's layout preview: MDATK's example dataset
    ("example_demo", "example", [], {"watermark": WATERMARK}, "example data of the demo layout"),
    ("example_strip", "example", layout("large_plus_2_strip", S, RMSD, HM, SEC), {"watermark": WATERMARK},
     "example data with a secondary-structure strip"),
    ("example_rep", "example", layout("equal_2", S, DA) + [rep({"style": "cartoon"})],
     {"watermark": WATERMARK, "example_frames": 24}, "schematic duplex with a cartoon (24 example frames)"),
]


# --- recording --------------------------------------------------------------------------------------

def error(exc):
    if isinstance(exc, KeyError) and exc.args:
        return {"type": "KeyError", "message": str(exc.args[0])}
    return {"type": type(exc).__name__, "message": rc.mask_text(str(exc))}


def project_for(name, inp, overrides):
    p = BASES[inp]()
    for assignment in overrides:
        cli.apply_override(p, assignment)
    return p


def positions(n):
    return sorted({0, n // 2, n - 1}) if n else []


def panel_record(k, panel, walker):
    rec = {"type": panel.type_name, "ok": bool(panel.ok), "auto_title": panel.auto_title,
           "axes": {"main": walker.id_of(panel.ax)}}
    if getattr(panel, "ax2", None) is not None:
        rec["axes"]["twin"] = walker.id_of(panel.ax2)
    for ax in panel.fig.axes:
        cb = getattr(ax, "_colorbar", None)
        if cb is not None and cb.mappable.axes is panel.ax:
            rec["axes"]["colorbar"] = walker.id_of(ax)
    try:
        rec["animated"] = [walker.id_of(a) for a in panel.animated_artists()] if panel.ok else []
    except AttributeError as exc:  # a panel whose setup failed half-way
        rec["animated_error"] = str(exc)
    for attr in ("radius", "data_lo", "data_hi", "view_note", "stage_note", "density_note", "per_row", "n_rows",
                 "cell_size", "source", "span"):
        if hasattr(panel, attr):
            rec[attr] = getattr(panel, attr)
    stage = getattr(panel, "stage", None)
    if stage is not None:
        rec["stage"] = {"renderer": type(stage.renderer).__name__, "quality": stage.quality,
                        "failed": bool(stage.failed)}
    return rec


def roles(panels):
    out = {}
    for k, panel in enumerate(panels):
        out[id(panel.ax)] = {"panel": k, "role": "main"}
        if getattr(panel, "ax2", None) is not None:
            out[id(panel.ax2)] = {"panel": k, "role": "twin"}
        for ax in panel.fig.axes:
            cb = getattr(ax, "_colorbar", None)
            if cb is not None and cb.mappable.axes is panel.ax:
                out[id(ax)] = {"panel": k, "role": "colorbar"}
    return out


def changed(record, first):
    """The id and the fields that differ from the first position's record of the same artist."""
    old = first.get(record["id"], {})
    return {k: v for k, v in record.items() if k == "id" or rc.plain(old.get(k)) != rc.plain(v)}


def run_case(name, inp, overrides, opts, description, inputs):
    p = project_for(name, inp, overrides)
    if inp == "example":
        folder = os.path.join(INPUTS, f"example_{name}")
        rc.save_dataset(example_dataset(p, n_frames=opts.get("example_frames", 120)), folder)
        ds, input_name = rc.load_dataset(folder), f"example_{name}"
    else:
        ds, input_name = inputs[inp], inp
    w, h = opts.get("size", SMALL["size"])
    dpi = opts.get("dpi", SMALL["dpi"])
    watermark = opts.get("watermark", "")
    rec = {"case": name, "description": description, "input": input_name,
           "overrides": [rc.mask_text(o) for o in overrides],
           "project": rc.project_dict(p), "size": [w, h], "dpi": dpi, "watermark": watermark,
           "positions": positions(ds.n_frames), "improvements": improvements_render.applied()}
    arrays = {}
    rc.reset_log()
    residues_module._WARNED.clear()  # warn_once: once per case, as in a fresh process
    try:
        fig, canvas, panels = build_figure(ds, p, w, h, dpi, watermark)
    except Exception as exc:  # noqa: BLE001 - MDATK's error is the expected behaviour
        rec["error"] = {"stage": "build", **error(exc)}
        rec["log"] = list(rc.LOG_LINES)
        return rec, arrays
    frames = []
    rec["frames"] = frames
    animated = [a for pn in panels if pn.ok for a in pn.animated_artists()]
    try:
        for n, pos in enumerate(rec["positions"]):
            for pn in panels:
                pn.update(pos)
            canvas.draw()
            walker = mpl_walker.Walker(fig, prefix=f"p{pos}/")
            if n == 0:
                figure = walker.figure(roles(panels))
                frames.append({"pos": pos, "figure": figure})
                first = {a["id"]: a for ax in figure["axes"] for a in ax["artists"]}
                rec["panels"] = [panel_record(k, pn, walker) for k, pn in enumerate(panels)]
            else:
                walker.register()
                frames.append({"pos": pos, "artists": [changed(walker.artist(a), first) for a in animated]})
            arrays.update(walker.arrays)
    except Exception as exc:  # noqa: BLE001
        rec["error"] = {"stage": f"pos {pos}", **error(exc)}
    rec["log"] = [rc.mask_text(line) for line in rc.LOG_LINES]
    return improvements_render.adjust("figure", name, rec), arrays


def main():
    for sub in (INPUTS, FIGURES):
        if os.path.isdir(sub):
            shutil.rmtree(sub)
        os.makedirs(sub)
    inputs = write_inputs()
    names = [c[0] for c in CASES]
    assert len(names) == len(set(names)), "duplicate case names"
    for case in CASES:
        rec, arrays = run_case(*case, inputs)
        rc.write_json(os.path.join(FIGURES, f"{case[0]}.json"), rec)
        if arrays:
            rc.write_npz(os.path.join(FIGURES, f"{case[0]}.npz"), arrays)
        print(case[0], rec.get("error", ""), flush=True)


if __name__ == "__main__":
    main()
