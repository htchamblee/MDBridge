#!/usr/bin/env python3
"""MDATK v0.0.4.1's analysis window driven by action scripts (tests/data/README.md "gui").

    python -I -B tests/parity/make_gui_cases.py MDATK_CHECKOUT DATA_DIR [--only NAME ...]

Writes DATA_DIR/gui/<topic>/<case>.json (topics: roundtrip, pages, ap, canvas) with the case's
input files (<case>/in/), the files MDBridge's window must write (<case>/out/) and MDATK's own,
where they differ (<case>/mdatk/). DATA_DIR is a tests/data tree: the cases read its committed
files (projects/, project-file/mdatk-0.1.1/, the demo files topology/demo_dry.parm7,
io/trajectories/demo_dry.nc, render/inputs/demo_dry_map.mrc) and point at them.

Every case runs twice, each time in a child process of this script (--worker): as MDATK v0.0.4.1
("mdatk") and with MDBridge's window improvements (improvements_gui.py, "improved"). Both import
MDATK's GTK 3 window (mdatk/gui/app.py) under the display of DISPLAY, or of an Xvfb this script
starts (-displayfd), with numpy's and OpenBLAS' kernels pinned (chem_common.pin_kernels, set in
this script's environment before numpy is imported and inherited by the workers; the density fits'
files are compared across machines as regenerate_all.py's FITTED_TEXT / FITTED_GUI) and with
render_common.setup(): Matplotlib with MDBridge's bundled fonts only (improvements_fonts.py),
MDBridge's earlier improvements patched in. A case gets a new work folder (the current folder) and XDG_CONFIG_HOME / HOME folders of its
own; the window is a MainWindow of a Gtk.Application, driven through its own methods and widgets
(the binder's widgets, the slot editors, the canvas' press / motion / release / scroll handlers with
synthetic events, save_project with the file dialog answered by the script), its canvas widget fixed
at the case's size. Pending canvas redraws and settle timers run before each action (MDATK's 120 /
250 ms timers are waited out, as a user would). MDATK's clock (projectfile._now) is fixed: inputs
are written at 2026-11-25T11:00:00, the window's saves count from the case's `clock` one second per
reading. Paths of the work folder and the config folder become {case_dir} / {config_dir}; cpptraj
and conda are never found.

Needs GTK 3 introspection data, PyGObject and pycairo for this Python
(tests/parity/requirements-gui.txt) and Xvfb when DISPLAY is unset. The output is reproducible.
"""

from __future__ import annotations

import base64
import datetime as real_dt
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import types

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import chem_common  # noqa: E402

# numpy's and OpenBLAS' kernels pinned before numpy is imported, as the chem generators pin them (the
# density fits' Powell path runs through BLAS dot products); the MDATK workers inherit the environment
chem_common.pin_kernels()

CLOCK = "2026-11-25T12:00:00"           # the window's clock of every case
INPUT_CLOCK = real_dt.datetime(2026, 11, 25, 11, 0, 0)  # input projects written by MDATK
CANVAS = {"width": 960, "height": 540}
TOPICS = ("roundtrip", "pages", "ap", "canvas", "tools")
DEMO = {"demo_dry.parm7": "topology/demo_dry.parm7", "demo_dry.nc": "io/trajectories/demo_dry.nc"}
MAP = {"map.mrc": "render/inputs/demo_dry_map.mrc"}
PROJECT = "p.mdatk.toml"

# ================================================================================================
# the cases

CASES = []


def case(topic, name, description, actions, *, files=None, build=None, start=PROJECT, items=(),
         mdatk_test="", settings=None, canvas=None, fitted=None):
    """files: name -> path under tests/data; build: name -> builder function name (an input project
    written by MDATK in the worker); start: the project the window starts with (None: none);
    fitted: the case runs a density fit (owner decision D10): expected file -> absolute tolerance of
    its decimal numbers (at least one unit of the last printed digit; a project file's [density]
    table within D10's angle, shift and CC tolerances; the rest of the text exact), and its log,
    reports and snapshot texts with their numbers within the last printed digit (the case needs
    "fitted")."""
    files = dict(files or {})
    if start and start not in files and start not in (build or {}):
        files[start] = "projects/default.mdatk.toml"  # what `mdatk init` writes
    CASES.append({"topic": topic, "name": name, "description": description, "items": list(items),
                  "mdatk_test": mdatk_test, "files": files, "build": dict(build or {}),
                  "start": start, "settings": settings, "canvas": canvas or CANVAS, "actions": actions,
                  "fitted": dict(fitted or {})})


def a(action, **kw):
    return {"action": action, **kw}


def snap(*keys, **kw):
    return a("snapshot", keys=list(keys), **kw)


DRAG_ROTATE = [[4, 0], [10, 2], [20, 5], [40, 8], [60, 10]]

# ---- round trips: open, save unchanged (lossless: AP-01 ... AP-06 for what the files hold) -------
for fixture in ("default", "full", "handwritten", "datetimes", "missing_log"):
    case("roundtrip", f"open_save_{fixture}", f"Start with projects/{fixture}.mdatk.toml and save it unchanged",
         [a("wait"), a("save")], files={PROJECT: f"projects/{fixture}.mdatk.toml"},
         items=["AP-06"] if fixture in ("full", "handwritten") else [],
         mdatk_test="tests/test_gui.py::test_widgets_roundtrip_project")
case("roundtrip", "open_save_split", "A project whose step log is in a separate file: open, save unchanged",
     [a("open", path="split.mdatk.toml"), a("save")], start=None,
     files={"split.mdatk.toml": "projects/split.mdatk.toml", "split.log.toml": "projects/split.log.toml"})
case("roundtrip", "open_params_log", "A run's parameter log opens as an unsaved project; Save asks for a name",
     [a("open", path="MDATK_params_log.toml"), a("save", answers={"name": "replay"})], start=None,
     files={"MDATK_params_log.toml": "projects/MDATK_params_log.toml"}, items=["AP-06", "AP-11"])
case("roundtrip", "open_save_mdatk_0_1_1", "A project file of MDATK 0.1.1 (planner keys 0.0.4.1 does not know)",
     [a("save")], files={PROJECT: "project-file/mdatk-0.1.1/plan.mdatk.toml"}, items=["D9"])
case("roundtrip", "lossy_values", "Values the widgets of MDATK cannot show survive open and save",
     [a("save")], build={PROJECT: "lossy_values"},
     items=["AP-01", "AP-02", "AP-03", "AP-04", "AP-05", "AP-19"])
case("roundtrip", "relative_paths", "Relative paths stay relative on Save; Save as into a sub-folder rebases them",
     [a("save"), a("save_as", name="sub/copy"), snap("inputs.topology", "inputs.trajectories", "output.out_dir")],
     build={PROJECT: "relative_paths"}, items=["AP-06", "D8"])

# ---- per page edits ------------------------------------------------------------------------------
case("pages", "input_page", "Input & cleaning: notes, files, frames, time, external data, cleaning", [
    a("type", key="notes", text="Edited in the window.\nSecond line"),
    a("type", key="inputs.topology", text="other.parm7"),
    a("table", key="inputs.trajectories", do="add"),
    a("set", key="inputs.trajectories.2", value="md3.nc"),
    a("table", key="inputs.trajectories", do="move", row=2, by=-1),
    a("table", key="inputs.trajectories", do="remove", row=0),
    a("set", key="inputs.start", value=3),
    a("type", key="inputs.stop", text="250"),
    a("set", key="inputs.step", value=5),
    a("type", key="inputs.dt_ps", text="0.002"),
    a("set", key="inputs.time_unit", value="ns"),
    a("table", key="inputs.external", do="add"),
    a("set", key="inputs.external.1.path", value="energy.dat"),
    a("set", key="inputs.external.1.label", value="Energy"),
    a("set", key="inputs.external.1.x_unit", value="ns"),
    a("set", key="inputs.external.1.columns", value="2, 3"),
    a("set", key="cleaning.enabled", value=True),
    a("type", key="cleaning.solvent_mask", text=":WAT,HOH"),
    a("set", key="cleaning.fit", value=True),
    a("set", key="cleaning.fit_reference", value="frame"),
    a("set", key="cleaning.reference_frame", value=10),
    a("type", key="cleaning.stop", text="500"),
    a("type", key="cleaning.extra_commands", text="center :1-24\nimage"),
    snap("inputs", "cleaning.stop"),
    a("save")], files={PROJECT: "projects/full.mdatk.toml"})
case("pages", "groups_page", "Residue groups: add, edit, recolour, move, remove; the group lists follow", [
    a("table", key="groups", do="add"),
    a("set", key="groups.3.name", value="Loop"),
    a("set", key="groups.3.residues", value="A:20-25"),
    a("set", key="groups.3.color", value="#aa22cc"),
    a("table", key="groups", do="move", row=3, by=-1),
    a("set", key="groups.0.selection", value="segid A and not name H*"),
    a("table", key="groups", do="remove", row=1),
    a("set", key="analyses.contacts.group_b", value="Loop"),
    snap("groups", slots=True),
    a("save")], files={PROJECT: "projects/full.mdatk.toml"},
    mdatk_test="tests/test_gui.py::test_swatch_survives_selection")
case("pages", "analyses_page", "Analyses: every section's fields", [
    a("set", key="analyses.use_pbc", value=True),
    a("set", key="analyses.proximity_search_radius", value=12.5),
    a("table", key="analyses.rmsd.targets", do="add"),
    a("set", key="analyses.rmsd.targets.1.label", value="Loop RMSD"),
    a("set", key="analyses.rmsd.targets.1.fit_selection", value="group:Target and name P"),
    a("set", key="analyses.rmsd.reference_frame", value=4),
    a("set", key="analyses.rmsd.mass_weighted", value=True),
    a("set", key="analyses.contacts.cutoff", value=4.25),
    a("set", key="analyses.contacts.group_a", value="Target"),
    a("set", key="analyses.hbonds.distance_cutoff", value=3.3),
    a("set", key="analyses.hbonds.use_angle", value=False),
    a("type", key="analyses.hbonds.donor_elements", text="N,O,S"),
    a("set", key="analyses.mg.inner_cutoff", value=2.4),
    a("set", key="analyses.mg.bridge_shell", value="outer"),
    a("set", key="analyses.stacking.scope", value="inter"),
    a("set", key="analyses.stacking.angle_cutoff", value=25),
    a("set", key="analyses.basepairs.stride", value=3),
    a("set", key="analyses.energy.dielectric", value=4),
    a("set", key="analyses.energy.exclude_bonded", value=2),
    a("type", key="analyses.energy.target_value", text="-55.5"),
    snap("analyses"),
    a("save")])
case("pages", "structure_page", "Structure view: orientation, representations, density, tags, appearance", [
    a("set", key="visual.rotation", value=[15.5, -30.0, 120.25]),  # slider values are floats
    a("set", key="visual.zoom", value=1.75),
    a("type", key="visual.center_on", text="B:12"),
    a("table", key="visual.representations", do="add"),
    a("set", key="visual.representations.0.style", value="licorice"),
    a("set", key="visual.representations.0.groups", value="Target"),
    a("set", key="visual.representations.0.color", value="element"),
    a("set", key="visual.representations.0.size", value=0.8),
    a("table", key="visual.representations", do="add"),
    a("set", key="visual.representations.1.opacity", value=0.5),
    a("set", key="visual.representations.1.sugars", value=True),
    a("set", key="visual.rep_rest", value="hide"),
    a("set", key="visual.backbone_atom", value="C4'"),
    a("set", key="density.level", value=2.5),
    a("set", key="density.color", value="#336699"),
    a("set", key="density.opacity", value=0.4),
    a("set", key="visual.label_every", value=5),
    a("type", key="visual.label_residues", text="1, 12, 24"),
    a("set", key="visual.structure_groups", value=["Aptamer"]),
    a("set", key="visual.base_color", value="#ffaa00"),
    a("set", key="visual.glow_mode", value="fade"),
    a("set", key="visual.theme", value="dark"),
    a("type", key="visual.title", text="Edited figure"),
    snap("visual", "density"),
    a("save")], mdatk_test="tests/test_gui.py::test_representation_editor_and_3d_preview")
case("pages", "layout_page", "Layout & panels: preset, slot types, data, groups, titles, more keys, options", [
    a("set", key="layout.preset", value="equal_4"),
    a("set", key="layout.panels.1.type", value="timeseries"),
    a("set", key="layout.panels.1.data", value="contacts/atom_pairs"),
    a("set", key="layout.panels.1.group_a", value="Target"),
    a("set", key="layout.panels.1.group_b", value="Aptamer"),
    snap(slots=True),
    a("type", key="layout.panels.1.title", text="My title"),
    a("set", key="layout.panels.1.group_b", value=""),
    a("set", key="layout.panels.2.type", value="secondary"),
    a("set", key="layout.panels.3.type", value="heatmap"),
    a("type", key="layout.panels.3.options", text='{"max_ticks": 5}'),
    a("type", key="layout.panels.3.more", text="hbonds/count"),
    snap("layout", slots=True),
    a("save")], mdatk_test="tests/test_gui.py::test_panel_group_dropdowns_and_auto_title",
    items=["AP-07"])
case("pages", "outputs_page", "Outputs: locations, video, charts, frame PDBs", [
    a("type", key="output.out_dir", text="results"),
    a("set", key="output.write_params_log", value=False),
    a("set", key="output.video.width", value=1280),
    a("set", key="output.video.height", value=720),
    a("set", key="output.video.fps", value=24),
    a("type", key="output.video.codec", text="libx265"),
    a("set", key="output.charts.formats", value=["png", "pdf"]),
    a("set", key="output.charts.width", value=6.5),
    a("type", key="output.charts.style", text="ggplot"),
    a("type", key="output.charts.color", text="tab:green"),
    a("type", key="output.charts.series", text="rmsd/Complex, hbonds/count"),
    a("set", key="output.pdb.enabled", value=True),
    a("type", key="output.pdb.frames", text="0, 10, last"),
    snap("output"),
    a("save")])

# ---- AP items --------------------------------------------------------------------------------------
case("ap", "ap01_set_out_of_range", "Values typed beyond the spin ranges are kept (validation flags them)", [
    a("set", key="analyses.contacts.cutoff", value=25.0),
    a("set", key="analyses.hbonds.angle_cutoff", value=60.0),
    a("set", key="output.video.width", value=7682),
    snap("analyses.contacts.cutoff", "analyses.hbonds.angle_cutoff", "output.video.width"),
    a("save")], items=["AP-01"])
case("ap", "ap02_typed_numbers", "Optional numbers keep every digit typed", [
    a("type", key="inputs.dt_ps", text="0.0123456789"),
    a("type", key="analyses.energy.target_value", text="-1234567.5"),
    a("save")], items=["AP-02"])
case("ap", "ap05_unknown_choice", "A stored renderer the drop-down does not list stays after another project",
     [a("open", path="other.mdatk.toml"), snap("visual.renderer"), a("save")],
     build={PROJECT: "renderer_software", "other.mdatk.toml": "renderer_vulkan"}, items=["AP-05", "AP-19"])
case("ap", "ap07_empty_slots", "A 4-slot preset with 2 panels is saved with 2", [
    a("set", key="visual.zoom", value=1.2), a("save")], build={PROJECT: "two_panels_equal_4"}, items=["AP-07"])
case("ap", "ap08_preset_round_trip", "equal_4 -> single -> equal_4 keeps every panel", [
    a("set", key="layout.preset", value="single"), snap("layout.panels"),
    a("set", key="layout.preset", value="equal_4"), snap("layout.panels", slots=True), a("save")],
    files={PROJECT: "projects/full.mdatk.toml"}, items=["AP-08"])
case("ap", "ap09_swap_bad_options", "A swap or preset change while slot options are no JSON object", [
    a("type", key="layout.panels.0.options", text="{bad"),
    a("canvas.drag", **{"from": {"slot": 1, "at": [0.5, 0.5]}, "to": {"slot": 2, "at": [0.5, 0.5]}}),
    a("set", key="layout.preset", value="equal_4"),
    a("type", key="layout.panels.0.options", text="[1, 2]"),
    a("save"),
    a("type", key="layout.panels.0.options", text=""),
    a("save")], files={PROJECT: "projects/full.mdatk.toml"}, items=["AP-09", "AP-10"])
case("ap", "ap10_form_errors", "Settings that do not parse: the save is refused with located messages", [
    a("type", key="inputs.stop", text="x"),
    a("type", key="analyses.energy.target_value", text="low"),
    a("save"), snap("inputs.stop"),
    a("type", key="inputs.stop", text=""), a("type", key="analyses.energy.target_value", text=""),
    a("save")], items=["AP-10", "AP-18"])
case("ap", "ap11_save_as_suffix", "Save as with a name typed without .mdatk.toml over an existing file", [
    a("set", key="notes", value="new text"),
    a("save_as", name="existing", answers={"overwrite": "cancel"}),
    a("save_as", name="existing", answers={"overwrite": "replace"})],
    build={PROJECT: "default_project", "existing.mdatk.toml": "notes_keep_me"}, items=["AP-11"])
case("ap", "ap12_start_broken", "The window started with a project that does not load",
     [a("save", answers={"name": "fresh"})], files={"broken.mdatk.toml": "projects/broken.mdatk.toml"},
     start="broken.mdatk.toml", items=["AP-12"])
case("ap", "ap13_unsaved_changes", "New, Open and closing with unsaved changes ask first", [
    a("set", key="analyses.contacts.cutoff", value=6.0), snap(modified=True),
    a("new", answers={"unsaved": "cancel"}), snap("analyses.contacts.cutoff"),
    a("open", path="other.mdatk.toml", answers={"unsaved": "save"}), snap("notes", modified=True),
    a("set", key="notes", value="edited after opening"),
    a("close", answers={"unsaved": "discard"})],
    build={PROJECT: "default_project", "other.mdatk.toml": "notes_keep_me"}, items=["AP-13", "D6"])
case("ap", "ap17_recent_keeps_settings", "Opening a project keeps the other keys of settings.toml", [
    a("open", path="other.mdatk.toml")],
    build={PROJECT: "default_project", "other.mdatk.toml": "notes_keep_me"}, items=["AP-17", "D16"],
    settings='# MDATK user settings.\n\n# Recently used project files, newest first.\nrecent = ["{case_dir}/p.mdatk.toml"]\n'
             "# RNAdo's Cat-proof mode: the keyboard is ignored (untick it with the mouse).\ncat_proof = true\n")

# ---- canvas -------------------------------------------------------------------------------------
case("canvas", "drag_rotate_swap", "Drag in the structure panel rotates; dragging panel 2 onto 3 swaps them", [
    a("wait"),
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.5]}, "path": [[40, 0]]}),
    snap("visual.rotation"),
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.6]}, "path": DRAG_ROTATE}),
    a("canvas.drag", **{"from": {"slot": 1, "at": [0.5, 0.5]}, "to": {"slot": 2, "at": [0.5, 0.5]}}),
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.05]}, "to": {"slot": 3, "at": [0.4, 0.6]}}),
    a("canvas.drag", **{"from": {"slot": 1, "at": [0.5, 0.5]}, "path": [[1, 1], [2, 0]]}),
    snap("visual.rotation", "layout.panels", slots=True), a("save")],
    files=DEMO, build={PROJECT: "demo_project"}, items=[],
    mdatk_test="tests/test_gui.py::test_drag_rotates_structure_and_swaps_panels")
case("canvas", "shift_drag_moves_structure", "Shift + drag moves even the structure panel", [
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.5]}, "to": {"slot": 1, "at": [0.5, 0.5]}, "shift": True}),
    snap("layout.panels"), a("save")], files=DEMO, build={PROJECT: "demo_project"})
case("canvas", "wheel_zoom", "Mouse wheel over the structure panel zooms; elsewhere it does nothing", [
    a("canvas.wheel", at={"slot": 0, "at": [0.5, 0.5]}, notches=[1, 1, 1, -1]),
    snap("visual.zoom"),
    a("canvas.wheel", at={"slot": 1, "at": [0.5, 0.5]}, notches=[1]),
    snap("visual.zoom"), a("save")], files=DEMO, build={PROJECT: "demo_project"},
    mdatk_test="tests/test_gui.py::test_slot_defaults_for_residue_panels_and_wheel_zoom")
case("canvas", "wheel_zoom_from_tiny", "The wheel zooms the stored zoom, not the clamped one shown", [
    a("canvas.wheel", at={"slot": 0, "at": [0.5, 0.5]}, notches=[1]), snap("visual.zoom"), a("save")],
    files=DEMO, build={PROJECT: "demo_tiny_zoom"}, items=["AP-01"])
case("canvas", "double_click_opens_slot", "Double-click a panel: the layout page shows its editor", [
    a("canvas.double_click", at={"slot": 2, "at": [0.5, 0.5]}), snap(page=True)],
    files=DEMO, build={PROJECT: "demo_project"})
case("canvas", "swap_into_empty_slot", "A panel dragged onto an empty slot", [
    a("canvas.drag", **{"from": {"slot": 1, "at": [0.5, 0.5]}, "to": {"slot": 3, "at": [0.5, 0.5]}}),
    snap("layout.panels"), a("save")], files=DEMO, build={PROJECT: "demo_two_panels"}, items=["AP-07"])
case("canvas", "map_move_rotate", "Map mode: drag moves the density map, Shift + drag rotates it", [
    a("wait"),
    a("click", key="density.centre_on_model"),
    snap("density", density_status=True),
    a("toggle", key="density.move_with_mouse", value=True),
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.5]}, "path": [[10, 0], [20, 5], [30, 10]]}),
    snap("density.fit_info", density_status=True),
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.5]}, "path": [[10, 0], [20, 0]], "shift": True}),
    a("toggle", key="density.move_with_mouse", value=False),
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.5]}, "path": [[0, 20]]}),
    snap("visual.rotation", "density"), a("save"),
    a("click", key="density.clear"), snap("density", density_status=True)],
    files={**DEMO, **MAP}, build={PROJECT: "demo_map"},
    mdatk_test="tests/test_gui.py::test_density_map_fit_shadow_and_drag")
# the fit part of test_density_map_fit_shadow_and_drag (milestone 6 W3-A): Global search, a drag in
# map mode, Refine from there; the jobs are waited for (MDATK: pump until the worker is done)
case("canvas", "density_fit", "Global search, a drag in map mode, Refine from the placement", [
    a("wait"),
    snap("density", density_status=True, modified=True),
    a("click", key="density.global_search"),
    snap("density", density_status=True, modified=True),
    a("toggle", key="density.move_with_mouse", value=True),
    a("canvas.drag", **{"from": {"slot": 0, "at": [0.5, 0.5]}, "path": [[10, 0], [20, 5], [30, 10]]}),
    a("toggle", key="density.move_with_mouse", value=False),
    snap("density.fit_info", density_status=True),
    a("click", key="density.refine"),
    snap("density", density_status=True),
    a("save"), snap(modified=True)],
    files={**DEMO, **MAP}, build={PROJECT: "demo_map"}, fitted={PROJECT: 0.0},
    mdatk_test="tests/test_gui.py::test_density_map_fit_shadow_and_drag")


# ---- the workflow tools (milestone 6 W3-B, owner decision D1): the PDB clean-up and structure
# builder windows of the Workflow menu (mdatk/gui/pdbdialog.py, builderdialog.py), driven through
# their widgets (tool.set: an entry's text, a drop-down's id, a check box, a spin button's value)
# and their buttons (tool.click: an insensitive or hidden button does nothing, as for a user); the
# files they write (and the project's step log) and their reports compared.

def tool(action, name, **kw):
    return a(f"tool.{action}", tool=name, **kw)


def tsnap(name, report=True, fields=(), sensitive=(), **kw):
    """A snapshot of a tool window: its report text, widget values, which widgets are sensitive."""
    return a("snapshot", keys=[], tool={"name": name, "report": report, "fields": list(fields),
                                        "sensitive": list(sensitive)}, **kw)


PDB_IN = {"in.pdb": "pdb_four_residues"}
case("tools", "pdb_tools_check_and_clean", "PDB clean-up: the output suggested, Check input, Amber names and "
     "renumbering, then REMARK 950 annotation (builtin engine, no stacking); the window closed", [
    tool("open", "pdb_tools"),
    tsnap("pdb_tools", report=False, fields=["output"], sensitive=["ann_engine", "reference", "chain_ids", "start"]),
    tool("set", "pdb_tools", key="input", value="in.pdb"),
    tsnap("pdb_tools", report=False, fields=["output"]),
    tool("click", "pdb_tools", key="check"),
    tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="scheme", value="amber"),
    tsnap("pdb_tools", report=False, fields=["strip_5p"]),
    tool("set", "pdb_tools", key="renumber", value=True),
    tsnap("pdb_tools", report=False, sensitive=["start", "per_chain", "keep_gaps"]),
    tool("click", "pdb_tools", key="write"),
    tsnap("pdb_tools", sensitive=["ann_engine"]),
    tool("set", "pdb_tools", key="annotate", value=True),
    tool("set", "pdb_tools", key="ann_engine", value="builtin"),
    tool("set", "pdb_tools", key="ann_stacking", value=False),
    tsnap("pdb_tools", report=False, sensitive=["ann_engine"]),
    tool("click", "pdb_tools", key="write"),
    tsnap("pdb_tools"),
    tool("click", "pdb_tools", key="close"),
    a("snapshot", keys=[], tool={"name": "pdb_tools", "report": False}),
    a("save")], build=PDB_IN, items=["TW-01"],
    mdatk_test="tests/test_gui.py::test_pdb_tools_window_checks_and_cleans")
case("tools", "pdb_tools_chains_and_sequence", "PDB clean-up: chain ids assigned, CCD names with forced L, SEQRES "
     "and a FASTA file, renumbering per chain; ids per segment, split, combine and sort; a reference PDB", [
    tool("open", "pdb_tools"),
    tool("set", "pdb_tools", key="input", value="blank.pdb"),
    tool("set", "pdb_tools", key="chain_mode", value="assign"),
    tool("set", "pdb_tools", key="scheme", value="ccd"),
    tool("set", "pdb_tools", key="l_ranges", value="B:1"),
    tool("set", "pdb_tools", key="seqres", value=True),
    tool("set", "pdb_tools", key="fasta", value="blank.fasta"),
    tool("set", "pdb_tools", key="renumber", value=True),
    tool("set", "pdb_tools", key="per_chain", value=True),
    tool("set", "pdb_tools", key="start", value=10),
    tool("click", "pdb_tools", key="check"),
    tsnap("pdb_tools"),
    tool("click", "pdb_tools", key="write"),
    tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="output", value="lists.pdb"),
    tool("set", "pdb_tools", key="fasta", value=""),
    tool("set", "pdb_tools", key="chain_mode", value="list"),
    tool("set", "pdb_tools", key="chain_ids", value="X, Y"),
    tool("set", "pdb_tools", key="combine", value="X,Y"),
    tool("set", "pdb_tools", key="combine_into", value="Q"),
    tool("set", "pdb_tools", key="records", value="pdb"),
    tool("set", "pdb_tools", key="terminals", value=False),
    tool("set", "pdb_tools", key="keep_gaps", value=True),
    tool("click", "pdb_tools", key="write"),
    tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="input", value="in.pdb"),
    tsnap("pdb_tools", report=False, fields=["output"]),
    tool("set", "pdb_tools", key="output", value="ref.pdb"),
    tool("set", "pdb_tools", key="combine", value=""),
    tool("set", "pdb_tools", key="combine_into", value=""),
    tool("set", "pdb_tools", key="chain_mode", value="reference"),
    tool("set", "pdb_tools", key="reference", value="in.pdb"),
    tool("set", "pdb_tools", key="copy_numbers", value=True),
    tool("set", "pdb_tools", key="split_at", value="B:4"),
    tool("set", "pdb_tools", key="sort_chains", value=True),
    tool("set", "pdb_tools", key="renumber", value=False),
    tool("click", "pdb_tools", key="write"),
    tsnap("pdb_tools")], build={**PDB_IN, "blank.pdb": "pdb_blank_chains"}, items=["TW-01"])
case("tools", "pdb_tools_errors", "PDB clean-up: no input, a missing input, no output, the output over the input, "
     "a reversed L range, a file that is not UTF-8 text, a missing reference PDB", [
    tool("open", "pdb_tools"),
    tool("click", "pdb_tools", key="check"), tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="input", value="missing.pdb"),
    tool("click", "pdb_tools", key="write"), tsnap("pdb_tools", fields=["output"]),
    tool("set", "pdb_tools", key="input", value="in.pdb"),
    tool("set", "pdb_tools", key="output", value=""),
    tool("click", "pdb_tools", key="write"), tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="output", value="in.pdb"),
    tool("click", "pdb_tools", key="write"), tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="output", value="out.pdb"),
    tool("set", "pdb_tools", key="l_ranges", value="A:4-1"),
    tool("click", "pdb_tools", key="check"), tsnap("pdb_tools"),
    tool("click", "pdb_tools", key="write"), tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="l_ranges", value=""),
    tool("set", "pdb_tools", key="input", value="latin1.pdb"),
    tool("set", "pdb_tools", key="output", value="out.pdb"),
    tool("click", "pdb_tools", key="check"), tsnap("pdb_tools"),
    tool("click", "pdb_tools", key="write"), tsnap("pdb_tools"),
    tool("set", "pdb_tools", key="input", value="in.pdb"),
    tool("set", "pdb_tools", key="chain_mode", value="reference"),
    tool("set", "pdb_tools", key="reference", value="nowhere.pdb"),
    tool("click", "pdb_tools", key="write"), tsnap("pdb_tools")],
    build={**PDB_IN, "latin1.pdb": "pdb_latin1"}, items=["TW-02", "CH-13", "CH-15"])
case("tools", "builder_helix", "Structure builder: Check and Write of Watson-Crick helices (D, L with Amber names and "
     "hydrogens, an explicit partner with a mismatch, a single strand); the checks of the other pairings and a "
     "triplex; the sensitivity rules", [
    tool("open", "builder"),
    tsnap("builder", report=False, sensitive=["partner", "chain2", "pairing", "triplex", "third", "chain3", "check"]),
    tool("set", "builder", key="seq", value="GGCAUGCC"),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="helix_out", value="helix_d.pdb"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="chirality", value="L"),
    tool("set", "builder", key="scheme", value="amber"),
    tool("set", "builder", key="hydrogens", value=True),
    tool("set", "builder", key="phosphate5", value=True),
    tool("set", "builder", key="chain1", value="P"),
    tool("set", "builder", key="chain2", value="Q"),
    tool("set", "builder", key="start1", value=101),
    tool("set", "builder", key="start2", value=-3),
    tool("set", "builder", key="helix_out", value="helix_l.pdb"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="partner_mode", value="explicit"),
    tsnap("builder", report=False, sensitive=["partner"]),
    tool("set", "builder", key="partner", value="GGCAUGUC"),
    tool("set", "builder", key="scheme", value=""),
    tool("set", "builder", key="hydrogens", value=False),
    tool("set", "builder", key="phosphate5", value=False),
    tool("set", "builder", key="chirality", value="D"),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="helix_out", value="mismatch.pdb"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="partner_mode", value="none"),
    tsnap("builder", report=False, sensitive=["partner", "chain2", "start2", "pairing", "triplex"]),
    tool("set", "builder", key="seq", value="5'-acgt-3'"),
    tool("set", "builder", key="helix_out", value="single.pdb"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="partner_mode", value="auto"),
    tool("set", "builder", key="seq", value="AAGG"),
    tool("set", "builder", key="pairing", value="hh"),
    tsnap("builder", report=False, sensitive=["triplex"]),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="pairing", value="reverse-hoogsteen"),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="pairing", value="wc"),
    tool("set", "builder", key="triplex", value="hoogsteen"),
    tsnap("builder", report=False, sensitive=["third", "chain3", "start3"]),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="third", value="UUCX"),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="helix_out", value=""),
    tool("click", "builder", key="write"), tsnap("builder")], items=["TW-01"],
    mdatk_test="tests/test_gui.py::test_builder_non_wc_helix_builds_in_the_background")
case("tools", "builder_mutate", "Structure builder: a duplex written, then mutated with and without compensatory "
     "partners (Check lists the residues, NOT FOUND, a spec that cannot be read); the output suggested", [
    tool("open", "builder"),
    tool("set", "builder", key="seq", value="GGCAUGCC"),
    tool("set", "builder", key="helix_out", value="duplex.pdb"),
    tool("click", "builder", key="write"),
    tool("page", "builder", page="mutate"),
    tool("set", "builder", key="mut_input", value="duplex.pdb"),
    tsnap("builder", report=False, fields=["mut_out"], sensitive=["check", "wobble"]),
    tool("set", "builder", key="mutations", value="A:3G, B:20C"),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="mutations", value="A:3G"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="mutations", value="A:3G, Z9"),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="mutations", value="A:4C"),
    tool("set", "builder", key="compensatory", value=False),
    tsnap("builder", report=False, sensitive=["wobble"]),
    tool("set", "builder", key="mut_out", value="nocomp.pdb"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="mut_input", value="in.pdb"),
    tsnap("builder", report=False, fields=["mut_out"]),
    tool("set", "builder", key="compensatory", value=True),
    tool("set", "builder", key="wobble", value=True),
    tool("set", "builder", key="mut_out", value="in_mut.pdb"),
    tool("set", "builder", key="mutations", value="A:1A"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="mut_input", value="latin1.pdb"),
    tool("click", "builder", key="check"), tsnap("builder"),
    tool("set", "builder", key="mut_input", value="duplex.pdb"),
    tool("set", "builder", key="mut_out", value="duplex.pdb"),
    tool("click", "builder", key="write"), tsnap("builder")],
    build={**PDB_IN, "latin1.pdb": "pdb_latin1"}, items=["TW-02", "CH-08", "CH-09"])
case("tools", "builder_fit", "Structure builder: Fit to map, a local refinement of a model moved off the map it was "
     "simulated from (in the background); the inputs checked first", [
    tool("open", "builder"),
    tool("page", "builder", page="fit"),
    tsnap("builder", report=False, sensitive=["check", "step"]),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="fit_input", value="in.pdb"),
    tsnap("builder", report=False, fields=["fit_out"]),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="fit_map", value="model.mrc"),
    tool("set", "builder", key="fit_out", value="in.pdb"),
    tool("click", "builder", key="write"), tsnap("builder"),
    tool("set", "builder", key="fit_out", value="in_fit.pdb"),
    tool("set", "builder", key="resolution", value=6.0),
    tool("set", "builder", key="search", value="local"),
    tsnap("builder", report=False, sensitive=["step"]),
    tool("click", "builder", key="write"),
    tool("wait", "builder"),
    tsnap("builder")], build={**PDB_IN, "model.mrc": "map_of_shifted_model"}, items=["TW-01"],
    fitted={"in_fit.pdb": 0.01, "p.mdatk.toml": 0.0})


# ================================================================================================
# main process: run the workers, merge, write

def start_xvfb():
    r, w = os.pipe()
    proc = subprocess.Popen(["Xvfb", "-displayfd", str(w), "-screen", "0", "1600x1200x24", "-nolisten", "tcp",
                             "-noreset"], pass_fds=(w,), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.close(w)
    data = b""
    while not data.endswith(b"\n"):
        chunk = os.read(r, 64)
        if not chunk:
            raise RuntimeError("Xvfb did not start")
        data += chunk
    os.close(r)
    return proc, ":" + data.decode().strip()


def run_worker(kind, mdatk, data, cases_file, out_file, env):
    cmd = [sys.executable, "-I", "-B", os.path.abspath(__file__), "--worker", kind, mdatk, data, cases_file, out_file]
    subprocess.run(cmd, check=True, env=env)
    with open(out_file, encoding="utf-8") as fh:
        return json.load(fh)


def dump(obj):
    return json.dumps(obj, indent=1, ensure_ascii=False) + "\n"


def write_bytes(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)


def content(entry):
    return base64.b64decode(entry["b64"]) if "b64" in entry else entry["text"].encode("utf-8")


OUTCOME_KEYS = ("log", "stderr", "dialogs", "questions", "snapshots", "canvas", "uncaught")


def needs_of(c, outcome, lossless):
    base = {"open", "wait", "save", "set", "click", "screenshot"}
    out = set()
    for act in c["actions"]:
        name = act["action"]
        if name.startswith("canvas."):
            out.add("canvas")
        elif name.startswith("tool."):
            out.add("tools")
        elif name not in base:
            out.add(name)
        if "answers" in act:
            out.add("answers")
        if name == "open" and outcome.get("dialogs"):
            out.add("open_failure")
    if c["start"]:
        out.add("start")
    if lossless:  # the expected files are not what MDATK wrote: MDBridge's edit model (W0-B, AP items)
        out.add("lossless")
    if c["settings"]:
        out.add("settings")
    if c["fitted"]:
        out.add("fitted")
    return sorted(out)


def merge(c, inputs, mdatk, improved, data):
    rel = f"gui/{c['topic']}/{c['name']}"
    root = os.path.join(data, rel)
    files = {}
    for name, entry in inputs.items():
        if name in c["files"]:
            files[name] = c["files"][name]
        else:
            write_bytes(os.path.join(root, "in", name), content(entry))
            files[name] = f"{rel}/in/{name}"
    expected = {}
    for name, entry in improved["files"].items():
        write_bytes(os.path.join(root, "out", name), content(entry))
        expected[name] = f"{rel}/out/{name}"
    outcome = {k: improved[k] for k in OUTCOME_KEYS if improved[k]}
    differs = {}
    mfiles = {}
    for name in sorted(set(mdatk["files"]) | set(improved["files"])):
        m, i = mdatk["files"].get(name), improved["files"].get(name)
        if m == i:
            continue
        if m is None:
            mfiles[name] = None
        else:
            write_bytes(os.path.join(root, "mdatk", name), content(m))
            mfiles[name] = f"{rel}/mdatk/{name}"
    if mfiles:
        differs["expected"] = mfiles
    for key in OUTCOME_KEYS:
        if mdatk[key] != improved[key]:
            differs[key] = mdatk[key]
    if improved["settings"] is not None:
        write_bytes(os.path.join(root, "settings.after.toml"), improved["settings"].encode("utf-8"))
        outcome["settings"] = f"{rel}/settings.after.toml"
    if mdatk["settings"] != improved["settings"]:
        differs["settings"] = None
        if mdatk["settings"] is not None:
            write_bytes(os.path.join(root, "mdatk", "settings.after.toml"), mdatk["settings"].encode("utf-8"))
            differs["settings"] = f"{rel}/mdatk/settings.after.toml"
    out = {"format": "mdbridge-gui-case", "version": 1, "description": c["description"],
           "source": "MDATK v0.0.4.1's GTK window driven by tests/parity/make_gui_cases.py, with MDBridge's "
                     "window improvements (tests/parity/improvements_gui.py)",
           "clock": CLOCK, "items": c["items"], "mdatk_test": c["mdatk_test"],
           "needs": needs_of(c, outcome, bool(mfiles) or "settings" in differs),
           "canvas": c["canvas"], "files": files}
    if c["fitted"]:  # expected files with fitted numbers: compared within their last printed digit
        out["fitted"] = c["fitted"]
    if c["settings"] is not None:
        write_bytes(os.path.join(root, "settings.before.toml"), c["settings"].encode("utf-8"))
        out["settings"] = f"{rel}/settings.before.toml"
    if c["start"]:
        out["start"] = {"project": c["start"]}
    out["actions"] = c["actions"]
    out["expected"] = expected
    out["log"] = outcome.pop("log", [])
    out["outcome"] = outcome
    out["mdatk"] = differs
    write_bytes(os.path.join(data, rel + ".json"), dump(out).encode("utf-8"))


def main(argv):
    only = []
    if "--only" in argv:
        k = argv.index("--only")
        argv, only = argv[:k], argv[k + 1:]
    if len(argv) != 2:
        sys.exit(__doc__.splitlines()[2].strip())
    mdatk, data = (os.path.abspath(p) for p in argv)
    cases = [c for c in CASES if not only or c["name"] in only]
    if only and len(cases) != len(only):
        sys.exit(f"unknown case in {only}")
    env = dict(os.environ)
    xvfb = None
    if not env.get("DISPLAY"):
        xvfb, env["DISPLAY"] = start_xvfb()
    tmp = tempfile.mkdtemp(prefix="gui-cases-")
    try:
        cases_file = os.path.join(tmp, "cases.json")
        with open(cases_file, "w", encoding="utf-8") as fh:
            json.dump(cases, fh)
        mres = run_worker("mdatk", mdatk, data, cases_file, os.path.join(tmp, "mdatk.json"), env)
        ires = run_worker("improved", mdatk, data, cases_file, os.path.join(tmp, "improved.json"), env)
    finally:
        if xvfb is not None:
            xvfb.terminate()
            xvfb.wait()
    if not only:
        for topic in TOPICS:
            shutil.rmtree(os.path.join(data, "gui", topic), ignore_errors=True)
    for c in cases:
        m, i = mres[c["name"]], ires[c["name"]]
        if m["inputs"] != i["inputs"]:
            raise RuntimeError(f"{c['name']}: the two passes wrote different inputs")
        if only:
            shutil.rmtree(os.path.join(data, "gui", c["topic"], c["name"]), ignore_errors=True)
        merge(c, m["inputs"], m, i, data)
        print(f"gui/{c['topic']}/{c['name']}.json")
    shutil.rmtree(tmp, ignore_errors=True)
    return 0


# ================================================================================================
# worker: one pass over the cases, in this process

class Ev:
    def __init__(self, x, y, state=0, button=1, typ=None, direction=None):
        self.x, self.y, self.state, self.button, self.type, self.direction = x, y, state, button, typ, direction


class Harness:
    """What the window asks and shows during one case."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.log, self.stderr, self.dialogs, self.questions = [], [], [], []
        self.snapshots, self.canvas, self.uncaught = [], [], []
        self.action = -1
        self.current = {}
        self.capturing = False

    def typed(self):
        act = self.current
        if act.get("action") == "save_as":
            return act["name"]
        if act.get("action") == "open":
            return act["path"]
        return act.get("answers", {}).get("name")

    def ask(self, kind, path):
        answer = self.current.get("answers", {}).get(kind, {"overwrite": "replace", "unsaved": "save"}[kind])
        self.questions.append({"action": self.action, "kind": kind, "path": path, "answer": answer})
        return answer

    def record_question(self, q):
        self.questions.append({"action": self.action, **q})


def worker(kind, mdatk, data, cases_file, out_file):
    display = os.environ.get("DISPLAY")
    import render_common  # pins the kernels before numpy; removes DISPLAY (Agg only)

    os.environ["DISPLAY"] = display
    os.environ["NO_AT_BRIDGE"] = "1"
    os.environ["GSETTINGS_BACKEND"] = "memory"
    os.environ["PATH"] = "/usr/bin:/bin"  # no cpptraj
    for var in ("CONDA_EXE", "CONDA_PREFIX", "MAMBA_ROOT_PREFIX", "CONDA_ROOT"):
        os.environ.pop(var, None)
    render_common.setup(mdatk)
    import gi

    gi.require_version("Gtk", "3.0")
    gi.require_version("Gdk", "3.0")
    from gi.repository import Gdk, GLib, Gtk

    if not Gtk.init_check([])[0]:
        raise RuntimeError("cannot initialise GTK")
    import mdatk as mdatk_pkg
    from mdatk import config, projectfile as PF
    from mdatk.cleaning import conda
    from mdatk.config import PanelSpec, Project, RepSpec, _to_plain
    from mdatk.gui import app as A
    from mdatk.gui import widgets as W
    from mdatk.render.layout import n_slots

    improved = kind == "improved"
    H = Harness()

    # -- fixed clock
    class Clock:
        t = INPUT_CLOCK

    class FixedDateTime(real_dt.datetime):
        @classmethod
        def now(cls, tz=None):
            v = Clock.t
            Clock.t = Clock.t + real_dt.timedelta(seconds=1)
            return cls(v.year, v.month, v.day, v.hour, v.minute, v.second, 123456)

    PF._dt = types.SimpleNamespace(datetime=FixedDateTime)
    conda.detect_ambertools_env = lambda *a_, **k_: None

    # -- the file dialogs answered by the script
    class FakeChooser:
        def __init__(self, *args, **kw):
            self.confirm = False
            self.path = None

        def __getattr__(self, name):
            return lambda *a_, **k_: None

        def set_do_overwrite_confirmation(self, value):
            self.confirm = bool(value)

        def run(self):
            typed = H.typed()
            if not typed:
                return Gtk.ResponseType.CANCEL
            self.path = os.path.abspath(typed)
            if self.confirm and os.path.exists(self.path) and H.ask("overwrite", self.path) != "replace":
                return Gtk.ResponseType.CANCEL
            return Gtk.ResponseType.OK

        def get_filename(self):
            return self.path

        def get_filenames(self):
            return [self.path] if self.path else []

    class GtkProxy:
        FileChooserDialog = FakeChooser

        def __getattr__(self, name):
            return getattr(Gtk, name)

    A.Gtk = GtkProxy()

    # -- what the window shows
    ready = f"{mdatk_pkg.PROJECT_TITLE} {mdatk_pkg.__version__} ready."
    append_log = A.MainWindow.append_log

    def append_log_(self, text):
        if H.capturing and not text.startswith(ready):
            H.log.extend(str(text).split("\n"))
        append_log(self, text)

    A.MainWindow.append_log = append_log_

    def dialog(kind_):
        def show(self, text, detail=""):
            self.append_log(f"{'ERROR: ' if kind_ == 'error' else ''}{text}{(': ' + detail) if detail else ''}")
            H.dialogs.append({"action": H.action, "kind": kind_, "text": text,
                              "detail": detail if improved else detail[:3000]})  # AP-18
        return show

    A.MainWindow.error = dialog("error")
    A.MainWindow.info = dialog("info")

    def excepthook(typ, value, tb):
        if not H.capturing:
            sys.__excepthook__(typ, value, tb)
            return
        H.uncaught.append({"action": H.action, "error": f"{typ.__name__}: {value}"})
        if os.environ.get("GUI_CASES_DEBUG"):
            sys.__excepthook__(typ, value, tb)

    sys.excepthook = excepthook

    if improved:
        import improvements_gui as IG

        IG.apply(A, W, config, PF, n_slots, H)

    gapp = Gtk.Application(application_id=None)
    gapp.register(None)

    # -- input projects written by MDATK
    def groups_ab():
        return [config.ResidueGroup("Aptamer", "1-12", "", "#2a78d6"), config.ResidueGroup("Target", "13-24", "", "#eb6834")]

    def demo_layout(p):
        p.layout.preset = "large_plus_3"
        p.layout.panels = [PanelSpec("structure", [], "Structure"), PanelSpec("timeseries", ["rmsd/Complex"], "RMSD"),
                           PanelSpec("timeseries", ["hbonds/count"], ""),
                           PanelSpec("heatmap", ["contacts/per_residue_a"], "", {"max_ticks": 6})]

    def demo_base():
        p = Project()
        p.inputs.topology = "demo_dry.parm7"
        p.inputs.trajectories = ["demo_dry.nc"]
        p.groups = groups_ab()
        p.output.video.width, p.output.video.height, p.output.video.dpi = 320, 180, 40
        demo_layout(p)
        return p

    def b_demo_project():
        return demo_base()

    def b_demo_tiny_zoom():
        p = demo_base()
        p.visual.zoom = 0.01
        return p

    def b_demo_two_panels():
        p = demo_base()
        p.layout.preset = "equal_4"
        p.layout.panels = p.layout.panels[:2]
        return p

    def b_demo_map():
        p = demo_base()
        p.density.map_path = "map.mrc"
        p.density.resolution = 6.0
        p.layout.preset = "single"
        p.layout.panels = p.layout.panels[:1]
        return p

    def b_default_project():
        return Project()

    def b_notes_keep_me():
        p = Project()
        p.notes = "keep me"
        return p

    def b_renderer_software():
        p = Project()
        p.visual.renderer = "software"
        return p

    def b_renderer_vulkan():
        p = Project()
        p.visual.renderer = "vulkan"
        return p

    def b_two_panels_equal_4():
        p = Project()
        p.layout.preset = "equal_4"
        p.layout.panels = [PanelSpec("structure", [], "Structure"), PanelSpec("timeseries", ["rmsd/Target"], "RMSD")]
        return p

    def b_lossy_values():
        p = Project()
        p.analyses.contacts.cutoff = 25.0          # spin 0.5 .. 20
        p.analyses.hbonds.angle_cutoff = 60.0      # spin 90 .. 180 (validation: 0 .. 180)
        p.analyses.energy.cutoff = 60.0            # spin 4 .. 50
        p.inputs.dt_ps = 0.0123456789              # "%g"
        p.analyses.energy.target_value = -1234567.5
        p.visual.glow_color = "tab:orange"         # Matplotlib colours GDK does not parse
        p.visual.proximity_color = "#ff000080"
        p.visual.base_color = "gold"
        p.visual.renderer = "vulkan"               # not in the drop-down (RN-20)
        p.analyses.basepairs.engine = "barnaba"
        p.output.charts.formats = ["jpg", "svg"]
        p.output.video.width = 7682
        p.visual.zoom = 0.01
        p.visual.rotation = [12.34567, -200.0, 0.0]
        p.density.opacity = 0.12345
        p.visual.representations = [RepSpec("cartoon", [], "auto", 1.23456, 0.5)]
        return p

    def b_relative_paths():
        p = Project()
        p.inputs.topology = "data/system.parm7"
        p.inputs.trajectories = ["data/md1.nc", "/abs/md2.nc"]
        p.output.out_dir = "results"
        p.density.map_path = "maps/emd.mrc"
        return p

    # -- input files of the tool cases (file_<name>(path) writes them; the same bytes in both passes)
    def test_pdbtools():
        import importlib.util

        if "mdatk_test_pdbtools" not in sys.modules:
            spec = importlib.util.spec_from_file_location("mdatk_test_pdbtools",
                                                          os.path.join(mdatk, "tests", "test_pdbtools.py"))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            sys.modules["mdatk_test_pdbtools"] = module
        return sys.modules["mdatk_test_pdbtools"]

    FOUR = [("A", [("G", 1, "G", False), ("C", 2, "C", False)]), ("B", [("0G", 3, "G", True), ("0U", 4, "U", True)])]

    def write_text(path, text):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)

    def b_file_pdb_four_residues(path):
        """MDATK's test_pdb_tools_window_checks_and_cleans input: D G C (chain A), L 0G 0U (chain B)."""
        write_text(path, test_pdbtools().make_pdb(FOUR))

    def b_file_pdb_blank_chains(path):
        write_text(path, test_pdbtools().make_pdb([(" ", [("G", 1, "G", False), ("C", 2, "C", False)]),
                                                   (" ", [("A", 1, "A", False), ("U", 2, "U", False)])]))

    def b_file_pdb_latin1(path):
        text = test_pdbtools().make_pdb(FOUR, header=["REMARK   1 MOD\u00c8LE SANS ACCENT"])
        with open(path, "wb") as fh:
            fh.write(text.encode("latin-1"))

    def b_file_map_of_shifted_model(path):
        import numpy as np

        from mdatk import builder as B, density as D, pdbtools as PT
        from mdatk.io.mrc import write_map

        xyz, elements = B.model_atoms(PT.read_pdb(test_pdbtools().make_pdb(FOUR)))
        write_map(path, D.simulate_map(np.asarray(xyz) + np.array([0.4, -0.3, 0.2]), D.atom_weights(elements), 6.0))

    builders = {k[2:]: v for k, v in locals().items() if k.startswith("b_")}

    # -- helpers on the window
    def flush(win):
        for _ in range(1000):
            while Gtk.events_pending():
                Gtk.main_iteration_do(False)
            c = win.canvas
            if c._pending is not None:
                GLib.source_remove(c._pending)
                c._pending = None
                c.render()
                continue
            if getattr(c, "_settle", None) is not None:
                GLib.source_remove(c._settle)
                c._settle_now()
                continue
            return

    def wait_jobs(win, timeout=300):
        t0 = time.time()
        while win.worker is not None and time.time() - t0 < timeout:
            Gtk.main_iteration_do(False)
            time.sleep(0.002)
        if win.worker is not None:
            raise RuntimeError("job did not finish")
        flush(win)

    def binding(win, key):
        found = [b for b in win.binder.bindings if b.path == key]
        if len(found) != 1:
            raise KeyError(f"no field {key}")
        return found[0]

    TABLES = {"groups": ("group_editor", ["name", "residues", "selection", "color"]),
              "inputs.trajectories": ("traj_editor", [None]),
              "inputs.external": ("ext_editor", ["path", "label", "x_column", "x_unit", "columns"]),
              "analyses.rmsd.targets": ("rmsd_editor", ["label", "selection", "fit_selection"])}
    SLOT_FIELDS = {"type": "type", "data": "data", "group_a": "ga", "group_b": "gb", "more": "more",
                   "title": "title", "options": "opts"}

    def table_cell(win, key):
        for prefix, (attr, cols) in TABLES.items():
            if key.startswith(prefix + "."):
                rest = key[len(prefix) + 1:].split(".")
                row = int(rest[0])
                col = cols.index(rest[1] if len(rest) > 1 else None)
                return getattr(win, attr), row, col, prefix
        return None

    def commit_entry(w, text):
        w.set_text(text)
        w.emit("activate")

    def type_into(widget, text):
        if isinstance(widget, Gtk.SpinButton):
            widget.set_text(text)
            widget.update()
        elif isinstance(widget, Gtk.TextView):
            widget.get_buffer().set_text(text)
        elif isinstance(widget, Gtk.ComboBoxText) and widget.get_has_entry():
            commit_entry(widget.get_child(), text)
        elif isinstance(widget, Gtk.Entry):
            commit_entry(widget, text)
        elif hasattr(widget, "entry"):
            commit_entry(widget.entry, text)
        else:
            raise TypeError(f"cannot type into {type(widget).__name__}")

    def set_field(win, key, value=None, text=None):
        cell = table_cell(win, key)
        if cell:
            editor, row, col, prefix = cell
            editor._edited(None, str(row), str(value if text is None else text), col)
            return
        parts = key.split(".")
        if key.startswith("visual.representations."):
            row = win.rep_rows[int(parts[2])]
            w = row[parts[3]]
            if isinstance(w, Gtk.ComboBoxText):
                w.set_active_id(str(value))
            elif isinstance(w, Gtk.SpinButton):
                w.set_value(float(value)) if text is None else type_into(w, text)
            elif isinstance(w, Gtk.CheckButton):
                w.set_active(bool(value))
            else:
                w.set_text(str(value if text is None else text))
            return
        if key.startswith("layout.panels."):
            k, field = int(parts[2]), parts[3]
            if improved:
                IG.slot_edited(win, k)
            w = win.slot_rows[k][SLOT_FIELDS[field]]
            if isinstance(w, Gtk.ComboBoxText):
                w.set_active_id(str(value))
            else:
                commit_entry(w, str(value if text is None else text))
            return
        if key == "layout.preset":
            win.preset_combo.set_active_id(str(value))
            return
        b = binding(win, key)
        if text is not None:
            type_into(b.widget, text)
            return
        if key == "visual.rotation":
            win._set_rotation(value)
        else:
            b.set(value)
        if improved:
            IG.on_set(b)(value)

    def plain(value):
        return json.loads(json.dumps(_to_plain(value), ensure_ascii=False, default=str))

    def slots_state(win):
        out = []
        for row in win.slot_rows:
            model = row["data"].get_model()
            out.append({"frame": row["frame"].get_label(), "type": row["type"].get_active_id(),
                        "data": row["data"].get_active_id(), "choices": [[r[1], r[0]] for r in model],
                        "group_a": row["ga"].get_active_id(), "group_b": row["gb"].get_active_id(),
                        "group_a_enabled": row["ga"].get_sensitive(), "group_b_enabled": row["gb"].get_sensitive(),
                        "more": row["more"].get_text(), "title": row["title"].get_text(),
                        "options": row["opts"].get_text()})
        return out

    def snapshot(win, act):
        rec = {"action": H.action}
        if act.get("keys"):
            try:
                project = win.collect(validate=False)
                rec["values"] = {k: plain(W.get_path(project, k)) for k in act["keys"]}
            except ValueError as exc:
                rec["error"] = str(exc)
        if act.get("slots"):
            rec["slots"] = slots_state(win)
        if act.get("page"):
            rec["page"] = win.side_stack.get_visible_child_name()
        if act.get("density_status"):
            rec["density_status"] = win.density_status.get_text()
        if act.get("modified") and improved:
            rec["modified"] = IG.modified(win)
        if act.get("tool"):
            rec["tool"] = tool_snapshot(win, act["tool"])
        H.snapshots.append(rec)

    def point(win, spec):
        x, y, w, h = win.canvas.slots[spec["slot"]]
        fx, fy = spec["at"]
        return x + fx * w, y + fy * h

    def canvas_record(win):
        rec = {"action": H.action, "slots": [list(s) for s in win.canvas.slots],
               "structure_slot": win.canvas.structure_slot, "events": []}
        H.canvas.append(rec)
        return rec["events"]

    def canvas_action(win, act):
        c = win.canvas
        flush(win)
        events = canvas_record(win)
        name = act["action"]
        if name == "canvas.wheel":
            x, y = point(win, act["at"])
            for n in act["notches"]:
                events.append(["scroll", x, y, n])
                c._on_scroll(None, Ev(x, y, direction=Gdk.ScrollDirection.UP if n > 0 else Gdk.ScrollDirection.DOWN))
                flush(win)
            return
        if name == "canvas.double_click":
            x, y = point(win, act["at"])
            for kind_, typ in (("press", Gdk.EventType.BUTTON_PRESS), ("release", None),
                               ("press", Gdk.EventType.BUTTON_PRESS), ("double", Gdk.EventType._2BUTTON_PRESS),
                               ("release", None)):
                events.append([kind_, x, y])
                if kind_ == "release":
                    c._on_release(None, Ev(x, y))
                else:
                    c._on_press(None, Ev(x, y, typ=typ))
            flush(win)
            return
        state = Gdk.ModifierType.SHIFT_MASK if act.get("shift") else 0
        x0, y0 = point(win, act["from"])
        moves = [point(win, act["to"])] if "to" in act else [(x0 + dx, y0 + dy) for dx, dy in act["path"]]
        events.append(["press", x0, y0])
        c._on_press(None, Ev(x0, y0, state, typ=Gdk.EventType.BUTTON_PRESS))
        for x, y in moves:
            events.append(["motion", x, y])
            c._on_motion(None, Ev(x, y, state))
        x, y = moves[-1]
        events.append(["release", x, y])
        try:
            c._on_release(None, Ev(x, y, state))
        except Exception as exc:  # noqa: BLE001 - in GTK the signal handler's exception is printed
            H.uncaught.append({"action": H.action, "error": f"{type(exc).__name__}: {exc}"})
        flush(win)

    def topology_exists(win):
        try:
            top = win.collect(validate=False).inputs.topology
        except ValueError:
            return False
        return bool(top) and os.path.isfile(top)

    def new_window(project, path):
        win = A.MainWindow(gapp, project, path)
        alloc = types.SimpleNamespace(x=0, y=0, width=CANVAS_SIZE[0], height=CANVAS_SIZE[1])
        win.canvas.get_allocation = lambda: alloc
        win.resize(1320, 920)
        win.show_all()
        if improved:
            win._mdb_base = os.path.dirname(os.path.abspath(path)) if path else os.getcwd()
        flush(win)
        return win

    # -- the tool windows (W3-B)
    TOOL_ACTIONS = {"pdb_tools": {"check": "run_check", "write": "run_clean", "close": "destroy"},
                    "builder": {"check": "check_current", "write": "run_current", "cancel": "_cancel",
                                "close": "destroy"}}
    tools_ready = []

    def tools_setup():
        """Once, before the first tool case: the step log's durations fixed (0.0: logged() times the
        step with time.monotonic), and in the improved pass MDBridge's chemistry and tool window
        improvements (improvements_chem.py, improvements_gui.apply_tools)."""
        if tools_ready:
            return
        tools_ready.append(True)
        PF.time = types.SimpleNamespace(monotonic=lambda: 0.0)
        if improved:
            import improvements_chem
            from mdatk.gui import builderdialog, pdbdialog

            improvements_chem.apply_mdatk()
            IG.apply_tools(pdbdialog, builderdialog)  # after: it binds the patched annotate

    def tool_window(win, name):
        dlg = {"pdb_tools": win.pdb_tools, "builder": win.builder}[name]
        if dlg is None:
            raise RuntimeError(f"the {name} window is not open")
        return dlg

    def tool_value(w):
        if isinstance(w, Gtk.Box) and "entry" in vars(w):
            return w.entry.get_text()
        if isinstance(w, Gtk.ComboBoxText):
            return w.get_active_id()
        if isinstance(w, Gtk.CheckButton):
            return w.get_active()
        if isinstance(w, Gtk.SpinButton):
            return w.get_value()
        return w.get_text()

    def tool_set(w, value):
        if isinstance(w, Gtk.Box) and "entry" in vars(w):
            w.entry.set_text(str(value))
        elif isinstance(w, Gtk.ComboBoxText):
            if not w.set_active_id(str(value)):
                raise ValueError(f"no item {value!r}")
        elif isinstance(w, Gtk.CheckButton):
            w.set_active(bool(value))
        elif isinstance(w, Gtk.SpinButton):
            w.set_value(float(value))
        else:
            w.set_text(str(value))

    def tool_action(win, act):
        tools_setup()
        name, what = act["tool"], act["action"][len("tool."):]
        if what == "open":
            {"pdb_tools": win.open_pdb_tools, "builder": win.open_builder}[name]()
            return
        dlg = tool_window(win, name)
        if what == "set":
            tool_set(getattr(dlg, act["key"]), act["value"])
        elif what == "page":
            dlg.notebook.set_current_page(("helix", "mutate", "fit").index(act["page"]))
        elif what == "click":
            button = getattr(dlg, f"{act['key']}_button", None)
            if button is not None and not (button.get_sensitive() and button.get_visible()):
                return  # a user cannot press it
            try:
                getattr(dlg, TOOL_ACTIONS[name][act["key"]])()
            except Exception as exc:  # noqa: BLE001 - in GTK the button handler's exception is printed
                H.uncaught.append({"action": H.action, "error": f"{type(exc).__name__}: {exc}"})
        elif what == "wait":
            t0 = time.time()
            while dlg.worker is not None and time.time() - t0 < 300:
                Gtk.main_iteration_do(False)
                time.sleep(0.002)
            if dlg.worker is not None:
                raise RuntimeError("the tool's job did not finish")
        else:
            raise ValueError(f"unknown action {act['action']}")

    def tool_snapshot(win, spec):
        dlg = {"pdb_tools": win.pdb_tools, "builder": win.builder}[spec["name"]]
        out = {"name": spec["name"], "open": dlg is not None}
        if dlg is None:
            return out
        if spec.get("report"):
            buf = dlg.report.get_buffer()
            out["report"] = buf.get_text(buf.get_start_iter(), buf.get_end_iter(), False)
        if spec.get("fields"):
            out["values"] = {k: tool_value(getattr(dlg, k)) for k in spec["fields"]}
        if spec.get("sensitive"):
            out["sensitive"] = {k: (getattr(dlg, k) if hasattr(dlg, k) else getattr(dlg, f"{k}_button")).get_sensitive()
                                for k in spec["sensitive"]}
        return out

    def run_action(state, act):
        win = state["win"]
        name = act["action"]
        H.current = act
        answers = act.get("answers", {})
        if name == "wait":
            if topology_exists(win):
                win.load_preview()
                wait_jobs(win)
            flush(win)
        elif name == "set":
            set_field(win, act["key"], value=act["value"])
        elif name == "type":
            set_field(win, act["key"], text=act["text"])
        elif name == "table":
            if act["key"] == "visual.representations":
                if act["do"] == "add":
                    win._add_rep_row(RepSpec())
                elif act["do"] == "remove":
                    win._remove_rep_row(win.rep_rows[act["row"]])
                else:
                    raise ValueError("representations do not move")
            else:
                editor = getattr(win, TABLES[act["key"]][0])
                if act["do"] == "add":
                    editor.add_row()
                else:
                    editor.view.get_selection().select_path(Gtk.TreePath(act["row"]))
                    if act["do"] == "remove":
                        editor.remove_selected()
                    else:
                        editor.move(act["by"])
        elif name == "click":
            def fit(search):
                win.fit_density(search)
                wait_jobs(win)
            {"density.centre_on_model": win.centre_density, "density.clear": win.clear_density_fit,
             "density.global_search": lambda: fit("global"), "density.refine": lambda: fit("manual"),
             "structure.reset_orientation": lambda: win._set_rotation([0, 0, 0]),
             "structure.reset_view": win._reset_view, "toolbar.refresh": win.refresh_canvas}[act["key"]]()
        elif name == "toggle":
            {"density.move_with_mouse": win.map_toggle}[act["key"]].set_active(bool(act["value"]))
        elif name.startswith("canvas."):
            canvas_action(win, act)
        elif name.startswith("tool."):
            tool_action(win, act)
        elif name == "snapshot":
            snapshot(win, act)
        elif name == "save":
            win.save_project()
        elif name == "save_as":
            win.save_project(save_as=True)
        elif name == "new":
            if not improved or IG.confirm_discard(win, answers, H.record_question):
                win.new_project()
                if improved:
                    win._mdb_base = os.getcwd()
        elif name == "open":
            if not improved or IG.confirm_discard(win, answers, H.record_question):
                win.open_project()
                if improved and win.project_path:
                    win._mdb_base = os.path.dirname(os.path.abspath(win.project_path))
        elif name == "close":
            if not improved or IG.confirm_discard(win, answers, H.record_question):
                win._on_close()
                win.destroy()
                state["win"] = None
        else:
            raise ValueError(f"unknown action {name}")
        if state["win"] is not None:
            if improved:
                IG.after_action(state["win"])
            flush(state["win"])

    def snapshot_files(folder):
        out = {}
        for here, dirs, names in os.walk(folder):
            dirs.sort()
            for n in sorted(names):
                if n.startswith(".") and n.endswith(".lock"):
                    continue
                full = os.path.join(here, n)
                with open(full, "rb") as fh:
                    out[os.path.relpath(full, folder).replace(os.sep, "/")] = fh.read()
        return out

    def entry(data_, subst):
        try:
            text = data_.decode("utf-8")
        except UnicodeDecodeError:
            return {"b64": base64.b64encode(data_).decode("ascii")}
        return {"text": subst(text)}

    with open(cases_file, encoding="utf-8") as fh:
        cases = json.load(fh)
    results = {}
    start_dir = os.getcwd()
    global CANVAS_SIZE
    for c in cases:
        H.reset()
        CANVAS_SIZE = (c["canvas"]["width"], c["canvas"]["height"])
        work = os.path.realpath(tempfile.mkdtemp(prefix="gui-case-"))
        cfg = os.path.realpath(tempfile.mkdtemp(prefix="gui-config-"))
        home = os.path.realpath(tempfile.mkdtemp(prefix="gui-home-"))
        os.environ["XDG_CONFIG_HOME"], os.environ["HOME"] = cfg, home

        def subst(text, work=work, cfg=cfg):
            return text.replace(work, "{case_dir}").replace(cfg, "{config_dir}")

        os.chdir(work)
        try:
            for name, rel in c["files"].items():
                os.makedirs(os.path.dirname(os.path.join(work, name)) or work, exist_ok=True)
                shutil.copyfile(os.path.join(data, rel), os.path.join(work, name))
            Clock.t = INPUT_CLOCK
            for name, builder in c["build"].items():
                if "file_" + builder in builders:
                    builders["file_" + builder](os.path.join(work, name))
                else:
                    builders[builder]().save(os.path.join(work, name))
            for n in os.listdir(work):
                if n.startswith(".") and n.endswith(".lock"):
                    os.remove(os.path.join(work, n))
            if c["settings"] is not None:
                os.makedirs(os.path.join(cfg, "mdatk"))
                with open(os.path.join(cfg, "mdatk", "settings.toml"), "w", encoding="utf-8") as fh:
                    fh.write(c["settings"].replace("{case_dir}", work))
            before = snapshot_files(work)
            inputs = {n: entry(v, subst) for n, v in before.items()}
            Clock.t = real_dt.datetime.fromisoformat(CLOCK)
            H.capturing = True
            project, path = None, None
            if c["start"]:
                start = os.path.join(work, c["start"])
                try:
                    project, path = Project.load(start), start
                    if PF.is_params_log(path):
                        path = None
                except (OSError, ValueError) as exc:
                    problem = f"Could not open {start}: {exc}"
                    H.stderr.append(f"warning: could not load {start}: {exc}")
            win = new_window(project, path)
            handler = A.GuiLogHandler(win)
            import logging
            logging.getLogger("mdatk").addHandler(handler)
            if c["start"] and project is None and improved:  # AP-12 (MDATK 0.1.0)
                win.append_log(f"ERROR: {problem}\n{IG.START_PROBLEM}")
            flush(win)
            state = {"win": win}
            for k, act in enumerate(c["actions"]):
                if state["win"] is None:
                    raise RuntimeError(f"{c['name']}: action after the window closed")
                H.action = k
                run_action(state, act)
            if state["win"] is not None:
                flush(state["win"])
                t0 = time.time()  # timers of the window (open_section's 60 ms scroll) fire before it goes
                while time.time() - t0 < 0.4:
                    flush(state["win"])
                    time.sleep(0.01)
                state["win"].destroy()
            logging.getLogger("mdatk").removeHandler(handler)
            while Gtk.events_pending():
                Gtk.main_iteration_do(False)
            H.capturing = False
            after = snapshot_files(work)
            files = {n: entry(v, subst) for n, v in after.items() if before.get(n) != v}
            spath = os.path.join(cfg, "mdatk", "settings.toml")
            settings = None
            if os.path.isfile(spath):
                with open(spath, encoding="utf-8") as fh:
                    settings = subst(fh.read())
            fix = lambda o: json.loads(subst(json.dumps(o, ensure_ascii=False)))  # noqa: E731
            results[c["name"]] = {"inputs": inputs, "files": files, "settings": settings,
                                  **{k: fix(getattr(H, k)) for k in OUTCOME_KEYS}}
        finally:
            os.chdir(start_dir)
            for folder in (work, cfg, home):
                shutil.rmtree(folder, ignore_errors=True)
    with open(out_file, "w", encoding="utf-8") as fh:
        json.dump(results, fh, ensure_ascii=False)
    return 0


CANVAS_SIZE = (CANVAS["width"], CANVAS["height"])

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        try:
            sys.exit(worker(*sys.argv[2:7]))
        except BaseException:
            sys.excepthook = sys.__excepthook__  # the window's hook records handler errors only
            raise
    sys.exit(main(sys.argv[1:]))
