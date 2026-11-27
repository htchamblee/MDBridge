"""MDBridge's improvements of MDATK v0.0.4.1's analysis window, applied on the test side
(docs/specs/m5-plan.md section 4, docs/DIFFERENCES.md "Interface", AP-xx).

make_gui_cases.py drives MDATK's GTK window; with apply() the window behaves as MDBridge's does
where MDBridge deliberately differs, so the recorded files and messages are MDBridge's expected
behaviour, never edited by hand. Each rule mirrors what MDBridge's window does (W0-B edit model,
W1 / W2 pages and canvas); where a later MDATK fixed the item, its design is followed and named.

  AP-01  no clamping on load: a field keeps the value it was loaded with until the user changes it
         (spin boxes and sliders clamp only what they show). Generally: a field the user did not
         change saves the loaded value verbatim (lossless round trip; also AP-02 ... AP-05 and the
         rotation sliders' 3-decimal rounding, which applies to values the user or the canvas set)
  AP-02  optional numbers are shown with the shortest round-trip text (repr), never "%g"
  AP-03  colour fields keep their text (Matplotlib colours MDATK accepts: "tab:orange", "#ff000080",
         "gold"); the picker writes "#rrggbb"
  AP-04  chart formats: the list as set, unknown entries kept
  AP-05  drop-downs keep a stored value that is not in their list (also AP-19)
  AP-06  paths are written as the file had them (relative stays relative; D8); Save as into another
         folder rebases relative paths onto the new folder
  AP-07  empty slots stay empty: layout.panels keeps its length; a slot becomes a panel when the user
         edits its editor or swaps a panel into it
  AP-08  a preset with fewer slots keeps the panels beyond them (render's "ignoring N extra
         panel(s)" warning) and a larger preset shows them again
  AP-09  a swap or preset change while a slot's options are no JSON object is refused with the
         located message (MDATK 0.1.0, c6d8d15: "Panels not swapped; fix this first: ...")
  AP-10  form errors are located: "Input & cleaning › Frames › Stop: 'x' is not a whole number"
         (the locations of MDATK 0.0.10's problems.py, ported below)
  AP-11  Save as: the overwrite question sees the final name (suffix appended first; MDATK 0.1.0,
         172df88)
  AP-12  a project that fails to load at the start is named in the window's log (MDATK 0.1.0,
         13d4fea: "ERROR: Could not open ...", "These are the default settings instead (unsaved).")
  AP-13  unsaved changes: New, Open and closing ask Save / Discard / Cancel (D6)
  PRJ-09 a later MDATK's file keeps its mdatk_version and its comments when the window saves it
         (W0-B's edit model writes it back byte for byte); PRJ-10 a step-log
         marker inside the notes does not start the log (MDATK 0.1.1 _log_start)
  AP-17  the recent list keeps the settings file's other keys (MDATK 0.0.14 add_recent; MDBridge
         also writes atomically, which a file comparison cannot show)
  AP-18  dialog details are not cut at 3 000 characters
  W2-A   the slot editors' data labels follow the settings they name after every change (MDATK
         relabelled them only when it rebuilt the editors: a group chosen for an analysis left
         "Aptamer–Target …" in the Data drop-down); the reference frame loads by itself (its summary
         goes to the status bar), so the window's log has no "--- Loading preview ---" and
         "Reference frame loaded: …" lines for it

The PDB clean-up and structure builder windows (milestone 6 W3-B, apply_tools; the chemistry itself
is improvements_chem.py's):

  TW-02  *Check input* (PDB clean-up) and *Check* of the builder's *Mutate* page report a file that
         is not UTF-8 text ("ERROR: 'utf-8' codec can't decode ...") instead of failing in the
         button's handler (MDATK 0.1.0: `except (OSError, ValueError)`)
  CH-13  the PDB clean-up names a reference PDB that is not a file ("ERROR: [Errno 2] No such file
         or directory: 'x.pdb'"), as the command line does; MDATK read the name as PDB text ("no
         ATOM/HETATM records found ...")

Not observable in the recorded files and messages, tested by the Qt tests of their packages:
AP-14 (close during a job), AP-15 (table edits commit on focus loss), AP-16 (About text with the
run-time citations: `mdbridge citations`, cli_parity.py), AP-19's offered choices, AP-20 (drawing
off the GUI thread).
"""

from __future__ import annotations

import copy
import json
import os
import types

PROJECT_SUFFIX = ".mdatk.toml"
_ERR = object()  # a getter that raised

# ---------------------------------------------------------------------------------------------
# AP-10: locations of MDATK 0.0.10+ (mdatk/problems.py of 0.1.1: TABS, SECTION_TABS, _LABELS,
# _LEAF_LABELS, PANEL_NAMES, describe; MDATK's own code, ported)

TABS = {"input": "Input & cleaning", "groups": "Residue groups", "analyses": "Analyses",
        "structure": "Structure view", "layout": "Layout & panels", "outputs": "Outputs"}
SECTION_TABS = {"notes": "input", "inputs": "input", "cleaning": "input", "groups": "groups",
                "analyses": "analyses", "visual": "structure", "density": "structure", "layout": "layout",
                "output": "outputs"}
_LABELS = {
    "notes": "Project notes",
    "inputs.topology": "Input files › Topology",
    "inputs.trajectories": "Input files › Trajectories",
    "inputs.start": "Frames › Start", "inputs.stop": "Frames › Stop", "inputs.step": "Frames › Step",
    "inputs.dt_ps": "Frames › Time step override",
    "inputs.external": "External data",
    "inputs.external.#": "External file {n}",
    "cleaning": "Trajectory cleaning",
    "groups.#": "Group {n}{name}",
    "groups.#.selection": "Selection (overrides)",
    "groups.#.color": "Colour",
    "analyses.rmsd": "RMSD",
    "analyses.rmsd.targets.#": "Target {n}{name}",
    "analyses.contacts": "Pairwise contacts",
    "analyses.hbonds": "Hydrogen bonds",
    "analyses.mg": "Mg²⁺ coordination and occupancy",
    "analyses.mg.occupancy_selection": "Occupancy residues",
    "analyses.stacking": "π–π stacking",
    "analyses.basepairs": "Base pairs",
    "analyses.energy": "Intermolecular energy",
    "visual.representations.#": "Representation {n}",
    "visual.structure_groups": "Appearance › Groups shown",
    "visual.proximity_set": "Appearance › Proximity set",
    "visual.color_mode": "Appearance › Colouring",
    "visual.glow_mode": "Appearance › Glow",
    "visual.center_on": "Orientation › Centre on",
    "visual.zoom": "Orientation › Zoom",
    "visual.label_residues": "Residue number tags › Label residues",
    "visual.rep_rest": "Representations › Other residues",
    "visual.backbone_atom": "Representations › Backbone atom",
    "visual.renderer": "Representations › Renderer",
    "density": "Density map",
    "density.map_path": "Map file",
    "layout.preset": "Grid preset",
    "layout.panels.#": "Panel {n}{name}",
    "layout.panels.#.series": "Data",
    "output.video": "Video",
    "output.charts": "Time charts",
    "output.pdb": "Frame PDBs",
    "output.pdb.selection": "Atom subset",
    "output.video.fps": "FPS",
}
_LEAF_LABELS = {
    "group_a": "Group A", "group_b": "Group B", "groups": "Groups", "selection": "Selection",
    "fit_selection": "Fit selection", "mg_selection": "Mg selection", "ligand_selection": "Ligand selection",
    "residues": "Residues", "name": "Name", "options": "Options", "style": "Style", "color": "Colouring",
}
PANEL_NAMES = {
    "structure": "Structure", "timeseries": "Time series", "dual_axis": "Dual-axis time series",
    "heatmap": "Residue heatmap", "contact_map": "Contact map", "sequence": "Sequence",
    "secondary": "Secondary structure",
}
SEP = " › "


def tab_of(field):
    return SECTION_TABS.get(field.split(".", 1)[0], "") if field else ""


def _item_name(parts, k, project):
    if project is None:
        return ""
    try:
        obj = project
        for part in parts[:k + 1]:
            obj = obj[int(part)] if part.isdigit() else getattr(obj, part)
    except (AttributeError, IndexError, TypeError, ValueError):
        return ""
    if parts[0] == "layout":
        kind = getattr(obj, "type", "")
        name = PANEL_NAMES.get(kind, kind)
    else:
        name = getattr(obj, "name", "") or getattr(obj, "label", "")
    return f" ({name})" if name else ""


def describe(field, project=None):
    """"Layout & panels › Panel 3 (Contact map) › Options" for "layout.panels.2.options"."""
    if not field:
        return ""
    parts = field.split(".")
    out = [TABS.get(tab_of(field), parts[0].capitalize())]
    for k in range(len(parts)):
        pattern = ".".join("#" if p.isdigit() else p for p in parts[:k + 1])
        label = _LABELS.get(pattern)
        if label is None and k == len(parts) - 1 and k > 0 and not parts[k].isdigit():
            label = _LEAF_LABELS.get(parts[k], parts[k].replace("_", " ").capitalize())
        if label:
            n = int(parts[k]) + 1 if parts[k].isdigit() else 0
            out.append(label.format(n=n, name=_item_name(parts, k, project) if "{name}" in label else ""))
    return SEP.join(out)


def field_message(kind, text, exc):
    """AP-10: the reason of a field that does not parse (MDATK: Python's int() / float() text)."""
    if kind == "int":
        return f"'{text}' is not a whole number"
    if kind == "float":
        return f"'{text}' is not a number"
    return str(exc)


def options_problem(k, panel_type, text, exc):
    """AP-10: located message of slot k's options (MDATK: "panel options '{bad': <json error>")."""
    project = types.SimpleNamespace(layout=types.SimpleNamespace(
        panels=[types.SimpleNamespace(type=panel_type) if j == k else None for j in range(k + 1)]))
    where = describe(f"layout.panels.{k}.options", project)
    if isinstance(exc, json.JSONDecodeError):
        return f"{where}: '{text}' is not valid JSON: {exc}"
    return f"{where}: '{text}' is not a JSON object"


# ---------------------------------------------------------------------------------------------
# AP-06: paths relative to the project folder

def load_as_written(project_cls, PF, path):
    """Project.load without resolve_paths (AP-06, D8): the project with its paths as the file has
    them, what MDBridge's edit model opens. Also make_gui_harness_cases.py's open."""
    return project_cls.from_dict(PF.read(path))


def rebase(project, old_dir, new_dir):
    """Save as into another folder: every relative path (the fields MDATK's resolve_paths makes
    absolute, the planner steps' paths included) is rewritten relative to the new folder; absolute
    paths, "~" paths and empty ones stay."""
    from mdatk.workflows import resolve_node_paths

    def fix(p):
        if not p or os.path.isabs(p) or p.startswith("~"):
            return p
        return os.path.relpath(os.path.normpath(os.path.join(old_dir, p)), new_dir)

    out = project.copy()
    out.workflow = project.workflow
    out.inputs.topology = fix(out.inputs.topology)
    out.inputs.trajectories = [fix(p) for p in out.inputs.trajectories]
    for ext in out.inputs.external:
        ext.path = fix(ext.path)
    out.cleaning.reference_file = fix(out.cleaning.reference_file)
    out.cleaning.output_dir = fix(out.cleaning.output_dir)
    if out.cleaning.cpptraj_path and os.sep in out.cleaning.cpptraj_path:
        out.cleaning.cpptraj_path = fix(out.cleaning.cpptraj_path)
    out.output.out_dir = fix(out.output.out_dir)
    out.output.cache_dir = fix(out.output.cache_dir)
    out.density.map_path = fix(out.density.map_path)
    out.workflow = copy.deepcopy(project.workflow)
    out.workflow.out_dir = fix(out.workflow.out_dir)
    for node in list(out.workflow.nodes) + list(out.workflow.scratch):
        resolve_node_paths(node, fix)
    return out


# ---------------------------------------------------------------------------------------------
# PRJ-09 / PRJ-10 in the window's saves

def later_version(path):
    """The mdatk_version of a project file a later MDATK wrote (MDBridge keeps it on save, PRJ-09),
    else None."""
    if not path or not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh.read(8192).splitlines()[:40]:
            if line.startswith("mdatk_version = "):
                version = line.split("=", 1)[1].strip().strip('"')
                try:
                    later = tuple(int(x) for x in version.split(".")) > (0, 0, 4, 1)
                except ValueError:
                    return None
                return version if later else None
    return None


def _open_after(line, in_triple, depth):
    """(in_triple, bracket depth) after `line` of a TOML value."""
    k = 0
    while k < len(line):
        if line.startswith('"""', k):
            in_triple = not in_triple
            k += 3
            continue
        ch = line[k]
        if in_triple:
            k += 2 if ch == "\\" else 1
            continue
        if ch == "#":
            break
        if ch in "\"'":
            end = k + 1
            while end < len(line) and line[end] != ch:
                end += 2 if (ch == '"' and line[end] == "\\") else 1
            k = end + 1
            continue
        if ch in "[{":
            depth += 1
        elif ch in "]}":
            depth -= 1
        k += 1
    return in_triple, depth


def _entries(text):
    """The entries of a project file's text: [signature, comment block, body lines], and the lines
    after the last one."""
    items, pending, table = [], [], ""
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip() or line.startswith("#"):
            pending.append(line)
            i += 1
            continue
        if line.startswith("["):
            table = line.strip()
            sig = ("table", table)
        else:
            sig = ("key", table, line.split("=", 1)[0].strip())
        body = [line]
        in_triple, depth = _open_after(line.split("=", 1)[1] if "=" in line and not line.startswith("[") else "",
                                       False, 0)
        i += 1
        while (in_triple or depth > 0) and i < len(lines):
            body.append(lines[i])
            in_triple, depth = _open_after(lines[i], in_triple, depth)
            i += 1
        items.append([sig, pending, body])
        pending = []
    return items, pending


def keep_file_comments(new_text, old_text):
    """MDBridge's save of a file a later MDATK wrote (W0-B, ProjectDocument::comments): each entry the
    file had keeps the comment the file had above it. Entries are matched in order by their table and
    key (so the steps of a plan, whose tables repeat, match one by one)."""
    new, tail = _entries(new_text)
    old, _ = _entries(old_text)
    j = 0
    for item in new:
        for k in range(j, min(len(old), j + 200)):
            if old[k][0] == item[0]:
                item[1] = old[k][1]
                j = k + 1
                break
    out = []
    for _sig, comments, body in new:
        out += comments + body
    return "\n".join(out + tail)


def log_start(PF, text):
    """MDATK 0.1.1's _log_start (projectfile.py; MDATK's own code, ported): the marker line that ends
    the settings; one inside a string (the notes) does not count."""
    first = at = text.find("\n" + PF.LOG_MARKER)
    while at >= 0:
        try:
            PF.tomlio.loads(text[:at + 1])
            return at
        except PF.tomlio.TOMLDecodeError:
            at = text.find("\n" + PF.LOG_MARKER, at + 1)
    if first < 0:
        return -1
    try:
        PF.tomlio.loads(text)
        return -1
    except PF.tomlio.TOMLDecodeError:
        return first


# ---------------------------------------------------------------------------------------------
# AP-13: unsaved changes

def modified(win):
    try:
        now = win.collect(validate=False).to_dict()
    except ValueError:
        return True
    return now != getattr(win, "_mdb_saved", None)


def mark_saved(win):
    try:
        win._mdb_saved = win.collect(validate=False).to_dict()
    except ValueError:
        win._mdb_saved = None


def confirm_discard(win, answers, record):
    """Before New, Open or closing: True when the action may go on."""
    if not modified(win):
        return True
    answer = answers.get("unsaved", "save")
    record({"kind": "unsaved", "path": win.project_path or "", "answer": answer})
    if answer == "cancel":
        return False
    if answer == "save":
        return bool(win.save_project())
    return True


# ---------------------------------------------------------------------------------------------
# apply

def apply(app, W, config, PF, n_slots, harness):
    """Patch MDATK's window (module `app`), its widgets (`W`), Project (`config`), projectfile
    (`PF`). `harness` is make_gui_cases.py's object: it answers the save dialog (`harness.typed`,
    `harness.ask(kind, path, answers)`) and records."""
    Binder = W.Binder

    # --- AP-01 ... AP-05: lossless fields ------------------------------------------------------
    def load(self, project):
        for b in self.bindings:
            value = W.get_path(project, b.path)
            b.set(value)
            b.loaded = copy.deepcopy(value)
            b.baseline = _ERR

    def take_baseline(self):
        for b in self.bindings:
            try:
                b.baseline = b.get()
            except (ValueError, TypeError):
                b.baseline = _ERR

    def store(self, project, skip=()):
        errors = []
        for b in self.bindings:
            if b.path in skip:
                continue
            try:
                current = b.get()
            except (ValueError, TypeError) as exc:
                text = getattr(b, "text", lambda: "")()
                message = field_message(getattr(b, "kind", ""), text, exc)
                if not any(message.startswith(t + SEP) for t in TABS.values()):  # not located yet
                    message = f"{describe(b.path, project)}: {message}"
                errors.append(message)
                continue
            unchanged = hasattr(b, "loaded") and b.baseline is not _ERR and current == b.baseline
            W.set_path(project, b.path, copy.deepcopy(b.loaded) if unchanged else current)
        return errors

    Binder.load = load
    Binder.take_baseline = take_baseline
    Binder.store = store

    def optional(kind, parse):
        def bind(self, path, w):
            def get():
                t = w.get_text().strip()
                return parse(t) if t else None
            # AP-02: shortest round-trip text, not "%g"
            show = (lambda v: "" if v is None else str(v)) if kind == "int" else \
                (lambda v: "" if v is None else repr(float(v)))
            self._add(path, w, get, lambda v: w.set_text(show(v)))
            binding = self.bindings[-1]
            binding.kind = kind
            binding.text = lambda: w.get_text().strip()
            return w
        return bind

    Binder.optional_int = optional("int", int)
    Binder.optional_float = optional("float", float)

    load_project = app.MainWindow.load_project

    def load_project_(self, project, path=None):
        self._mdb_version = later_version(path)
        self._mdb_source = None
        if self._mdb_version:  # the settings text of the file, for its comments
            with open(path, "r", encoding="utf-8") as fh:
                text = fh.read()
            at = log_start(PF, text)
            self._mdb_source = text[:at + 1] if at >= 0 else text
        self._mdb_count = len(project.layout.panels)
        self._mdb_extra = []
        self._mdb_preset = project.layout.preset
        load_project(self, project, path)
        self.binder.take_baseline()
        mark_saved(self)

    app.MainWindow.load_project = load_project_

    # --- AP-06: no absolute paths on load; Save as rebases ---------------------------------------
    def project_load(cls, path):
        return load_as_written(cls, PF, path)

    config.Project.load = classmethod(project_load)

    def save_project(self, save_as=False):
        try:
            project = self.collect(validate=False)
        except ValueError as exc:
            self.error("Some settings are invalid", str(exc))
            return False
        path = self.project_path
        if save_as or not path:
            typed = harness.typed()
            if not typed:
                return False
            path = os.path.abspath(typed)
            if not path.endswith(".toml"):  # AP-11: before the question, not after it
                path += PROJECT_SUFFIX
            if os.path.exists(path) and harness.ask("overwrite", path) != "replace":
                return False
        base = getattr(self, "_mdb_base", None)
        new_base = os.path.dirname(os.path.abspath(path))
        rebased = bool(base) and os.path.normpath(base) != new_base
        if rebased:
            project = rebase(project, base, new_base)
        version, project_text = PF.__version__, PF.project_text
        if getattr(self, "_mdb_version", None):  # PRJ-09: a later MDATK's file keeps its mdatk_version
            PF.__version__ = self._mdb_version   # and the comments it had
            source = self._mdb_source
            PF.project_text = lambda p, log_file="": keep_file_comments(project_text(p, log_file), source)
        try:
            project.save(path)
        except OSError as exc:
            self.error("Could not save the project", str(exc))
            return False
        finally:
            PF.__version__, PF.project_text = version, project_text
        self._mdb_base = new_base
        if rebased:  # the fields show the rebased paths
            self.load_project(project, path)
        self.project = project
        self.project_path = path
        self.header.set_subtitle(path)
        PF.add_recent(path)
        self.append_log(f"Saved project {path}")
        self.binder.take_baseline()
        for b in self.binder.bindings:  # what was saved is the new loaded state
            b.loaded = copy.deepcopy(W.get_path(project, b.path))
        mark_saved(self)
        return True

    app.MainWindow.save_project = save_project

    # --- PRJ-10: the step log starts at the marker that ends the settings (MDATK 0.1.1 _log_start:
    # a marker pasted into the notes is part of a string and starts nothing) ---------------------
    def existing_log(path):
        if not os.path.isfile(path):
            return "", ""
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read()
        log_file = ""
        for line in text.splitlines()[:60]:
            if line.startswith("log_file = "):
                try:
                    log_file = PF.tomlio.loads(line)["log_file"]
                except PF.tomlio.TOMLDecodeError:
                    log_file = ""
                break
        at = log_start(PF, text)
        return (text[at + 1:] if at >= 0 else ""), log_file

    PF._existing_log = existing_log

    # --- AP-07 / AP-08 / AP-09: slots ------------------------------------------------------------
    def read_slots(self):
        specs = []
        project = None
        for k, row in enumerate(self.slot_rows):
            text = row["opts"].get_text().strip()
            ptype = row["type"].get_active_id() or "timeseries"
            try:
                options = json.loads(text) if text else {}
            except json.JSONDecodeError as exc:
                raise ValueError(options_problem(k, ptype, text, exc)) from exc
            if not isinstance(options, dict):
                raise ValueError(options_problem(k, ptype, text, None))
            keys = []
            if app.PANEL_TYPES.get(ptype) is None or app.PANEL_TYPES[ptype].max_series > 0:
                project = project or self._project_for_choices()
                primary = self._primary_key(row, project)
                keys = ([primary] if primary else []) + [
                    s.strip() for s in row["more"].get_text().split(",") if s.strip()]
            specs.append(app.PanelSpec(ptype, keys, row["title"].get_text().strip(), options))
        return specs

    app.MainWindow._read_slots = read_slots

    def real_panels(self):
        return self._read_slots()[:self._mdb_count] + list(self._mdb_extra)

    def rebuild_slots(self, specs=None):
        if specs is None:
            specs = real_panels(self)
        specs = list(specs)
        preset = self.preset_combo.get_active_id() or "single"
        n = n_slots(preset)
        self._mdb_count = min(len(specs), n)
        self._mdb_extra = specs[n:]
        for child in self.slot_box.get_children():
            self.slot_box.remove(child)
        self.slot_rows = []
        project = self._project_for_choices()
        loading, self._loading = self._loading, True
        try:
            for k in range(n):
                spec = specs[k] if k < len(specs) else app.PanelSpec("timeseries", [], "")
                row = self._slot_editor(k, spec, preset, project)
                self.slot_rows.append(row)
                self.slot_box.pack_start(row["frame"], False, False, 0)
        finally:
            self._loading = loading
        self.slot_box.show_all()

    app.MainWindow._rebuild_slots = rebuild_slots

    def get_layout(self):
        from mdatk.config import LayoutParams
        return LayoutParams(preset=self.preset_combo.get_active_id() or "single", panels=real_panels(self))

    app.MainWindow._get_layout = get_layout

    set_layout = app.MainWindow._set_layout

    def set_layout_(self, layout):
        self._mdb_preset = layout.preset
        set_layout(self, layout)

    app.MainWindow._set_layout = set_layout_

    def on_preset_changed(self):
        if self._loading:
            return
        try:
            specs = real_panels(self)
        except ValueError as exc:
            self.append_log(f"Grid preset not changed; fix this first: {exc}")
            loading, self._loading = self._loading, True
            try:
                self.preset_combo.set_active_id(self._mdb_preset)
            finally:
                self._loading = loading
            return
        self._mdb_preset = self.preset_combo.get_active_id() or "single"
        self._rebuild_slots(specs)
        self.refresh_canvas()

    app.MainWindow._on_preset_changed = on_preset_changed

    def on_canvas_swap(self, i, j):
        try:
            specs = self._read_slots()
        except ValueError as exc:
            self.append_log(f"Panels not swapped; fix this first: {exc}")
            return
        if max(i, j) >= len(specs):
            return
        specs[i], specs[j] = specs[j], specs[i]
        count = self._mdb_count
        if min(i, j) < count:  # a panel moved into an empty slot: that slot holds a panel now
            count = max(count, i + 1, j + 1)
        self._rebuild_slots(specs[:count] + list(self._mdb_extra))
        self.append_log(f"Swapped panels {i + 1} and {j + 1}")
        self.refresh_canvas()

    app.MainWindow._on_canvas_swap = on_canvas_swap

    # the canvas' wheel zoom works on the project's value (AP-01: the loaded value, not the clamped
    # one the spin box shows), within MDATK's wheel limits
    def on_canvas_zoom(self, factor):
        b = next(b for b in self.binder.bindings if b.path == "visual.zoom")
        try:
            current = b.get()
        except (ValueError, TypeError):
            current = b.baseline
        zoom = b.loaded if current == b.baseline else current
        value = min(50.0, max(0.05, zoom * factor))
        b.loaded = value
        self.zoom_spin.set_value(value)
        b.baseline = b.get()

    app.MainWindow._on_canvas_zoom = on_canvas_zoom

    # --- W2-A: the reference frame loads by itself (status bar, not the log) ---------------------
    logged = app.MainWindow.append_log

    def append_log(self, text):
        if text == "--- Loading preview ---" or str(text).startswith("Reference frame loaded: "):
            return
        logged(self, text)

    app.MainWindow.append_log = append_log

    # --- AP-17: the settings file's other keys stay (MDATK 0.0.14 add_recent) --------------------
    docs = {"recent": "Recently used project files, newest first.",
            "cat_proof": "RNAdo's Cat-proof mode: the keyboard is ignored (untick it with the mouse)."}

    def add_recent(path):
        path = os.path.abspath(path)
        recent = [path] + [p for p in PF.recent_projects() if p != path]
        try:
            data = PF.tomlio.load(PF._settings_path())
        except (OSError, PF.tomlio.TOMLDecodeError):
            data = {}
        if not isinstance(data, dict):
            data = {}
        data["recent"] = recent[:PF.RECENT_MAX]
        target = PF._settings_path()
        try:
            os.makedirs(os.path.dirname(target), exist_ok=True)
            text = PF.tomlio.dumps(data, lambda p, s: docs.get(p[-1] if p else "", None),
                                   [f"{PF.PROJECT_TITLE} user settings."])
            with open(target, "w", encoding="utf-8") as fh:
                fh.write(text)
        except OSError as exc:
            PF.log.warning("could not update the recent projects list: %s", exc)

    PF.add_recent = add_recent


def on_set(binding):
    """A value the user chose in a field: it is what the field holds now (AP-01: no clamping)."""
    def after(value):
        try:
            shown = binding.get()
        except (ValueError, TypeError):
            shown = None
        if isinstance(shown, float) and isinstance(value, int) and not isinstance(value, bool):
            value = float(value)  # a number field of a float setting stores a float (spin_float), as MDBridge's
        binding.loaded = copy.deepcopy(value)
        try:
            binding.baseline = binding.get()
        except (ValueError, TypeError):
            binding.baseline = _ERR
    return after


def after_action(win):
    """W2-A: the slot editors' data labels name the settings as they are now (MDBridge's editors
    follow every change; MDATK relabelled them only when it rebuilt them). Rows' texts only: the
    choices, the selection and the titles stay."""
    rows = getattr(win, "slot_rows", None)
    if not rows:
        return
    from mdatk.analysis.base import key_label
    project = win._project_for_choices()
    for row in rows:
        for r in row["data"].get_model():
            if r[1]:
                r[0] = f"{key_label(r[1], project)}   [{r[1]}]"


def slot_edited(win, k):
    """AP-07: a slot whose editor the user changes holds a panel (and every slot before it)."""
    win._mdb_count = max(win._mdb_count, k + 1)


START_PROBLEM = "These are the default settings instead (unsaved)."  # AP-12 (MDATK 0.1.0)


def apply_tools(pdbdialog, builderdialog):
    """MDBridge's tool window improvements (TW-02, CH-13 above) on MDATK's two windows' modules."""
    import errno

    real = pdbdialog.pdbtools

    def read_pdb(source, *args, **kwargs):
        try:
            return real.read_pdb(source, *args, **kwargs)
        except UnicodeDecodeError as exc:  # TW-02 (MDATK 0.1.0 catches ValueError)
            raise real.PDBError(str(exc)) from None

    def clean(source, *args, chain_reference=None, **kwargs):
        if isinstance(chain_reference, str) and not os.path.isfile(chain_reference):  # CH-13
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), chain_reference)
        return real.clean(source, *args, chain_reference=chain_reference, **kwargs)

    proxy = types.SimpleNamespace(**{k: getattr(real, k) for k in dir(real) if not k.startswith("__")})
    proxy.read_pdb, proxy.clean = read_pdb, clean
    pdbdialog.pdbtools = proxy
    builderdialog.pdbtools = proxy
    from mdatk import annotate

    pdbdialog.annotate = annotate.annotate  # improvements_chem's (CH-14, D13): the window imported the original
