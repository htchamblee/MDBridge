# MDATK bug report: issues still present in 0.1.1

- **Date:** 2026-11-25
- **Source:** found while building MDBridge, the C++ successor of MDATK; checked against MDATK 0.1.1 at
  commit `1efa0ef44e9cb0d04b9c04a409919f7913f59fdb` (tag-less release commit "0.1.1", 2026-11-25).
  Every item was first seen in MDATK 0.0.4.1 and is still present at that commit.
- **Test environment:** Python 3.13.18, MDAnalysis 2.11.0, numpy 2.5.5, scipy 1.18.1, matplotlib 3.11.2;
  the GUI items with PyGObject 3.48.2 / GTK 3 under Xvfb (`xvfb-run -a`). MDATK was run from the checkout
  (not installed), with `python -I -B` and `sys.path.insert(0, "<checkout>")`.

## How to use this report

- Each item has a stable id (`MDATK-BUG-001` …) and the MDBridge register id it came from (`AP-xx` GUI,
  `AN-xx` analyses, `SEL-xx` selections, `RD-xx` / `LD-xx` readers and loader, `PRJ-xx` project files,
  `CLI-xx` command line, `RN-xx` rendering, `MRC-x` maps). Use the MDBridge id only for cross-reference.
- `file:line` is MDATK 0.1.1 (the commit above), relative to the repository root.
- **Status:** *Verified* = reproduced against 0.1.1 with the command shown (output quoted as observed);
  *Confirmed by reading* = the code at the cited lines still does what is described (not run).
- **Inherited from MDAnalysis** marks items whose root cause is MDAnalysis behaviour that reaches users
  through MDATK. They describe behaviour only; the suggested fix is a workaround in MDATK's own code.
- "Suggested fix" says what MDBridge did, in words, and why. Prefer whatever fits MDATK's design.
- Reproductions assume a shell in the repository root, `mdatk` meaning `python -m mdatk`, and Python
  snippets run with the checkout on `sys.path`. Several use the demo system:
  `mdatk demo demo --n-frames 6 --quiet` (writes `demo/demo.parm7`, `demo/demo.nc`,
  `demo/demo.mdatk.toml`).
- Severity scale: **critical data loss** (user settings or results silently lost or rewritten) >
  **wrong results** > **crash** (traceback or abort) > **misleading message** > **minor**.

## Summary

85 items: 9 critical data loss, 26 wrong results, 17 crash, 15 misleading message, 18 minor.
54 Verified, 31 Confirmed by reading, 0 Unverified (the appendix lists further reader-message items that
were not re-run).

| Id | Area | Severity | Title | Status |
|---|---|---|---|---|
| MDATK-BUG-001 | GUI save | critical data loss | Spin buttons clamp loaded values, saving rewrites valid settings | Verified |
| MDATK-BUG-002 | GUI save | critical data loss | Optional numbers rounded to 6 significant digits on save | Verified |
| MDATK-BUG-003 | GUI save | critical data loss | Colour buttons turn valid Matplotlib colours into black or rewrite them | Verified |
| MDATK-BUG-004 | GUI save | critical data loss | Chart formats other than png/svg/pdf are dropped on save | Verified |
| MDATK-BUG-005 | GUI save | critical data loss | Combo boxes replace an unknown value with the previous project's choice | Verified |
| MDATK-BUG-006 | GUI save | critical data loss | GUI save makes every path absolute | Verified |
| MDATK-BUG-007 | GUI save | critical data loss | Empty layout slots are saved as real panels | Verified |
| MDATK-BUG-008 | GUI save | critical data loss | Switching to a smaller preset deletes panels for good | Verified |
| MDATK-BUG-009 | GUI | critical data loss | New / Open / Close / quit discard unsaved edits without asking | Confirmed by reading |
| MDATK-BUG-011 | Analyses | wrong results | Element of an atom depends on which group it is evaluated in | Verified |
| MDATK-BUG-012 | Project validation | wrong results | NaN / negative / zero cutoffs, radii and dielectric are accepted | Verified |
| MDATK-BUG-013 | RMSD | wrong results | Two RMSD targets with the same data key: one series silently lost | Verified |
| MDATK-BUG-014 | Loader | wrong results | `inputs.dt_ps` negative or NaN accepted: negative / NaN time axis | Verified |
| MDATK-BUG-015 | RMSD | wrong results | RMSD of 1–2 atoms is NaN or off by up to 8 Å (inherited from MDAnalysis) | Verified |
| MDATK-BUG-016 | RMSD | wrong results | Fit selection of < 3 or collinear atoms: silent no rotation / arbitrary rotation | Verified |
| MDATK-BUG-017 | Periodic geometry | wrong results | Periodic cutoff decisions depend on box size; far points get wrong distances (inherited from MDAnalysis) | Verified |
| MDATK-BUG-018 | Analyses | wrong results | Selections evaluated on whatever frame is current at setup, always periodic | Confirmed by reading |
| MDATK-BUG-019 | Energy | wrong results | No bonded exclusions (silently) when the topology has no bonds | Confirmed by reading |
| MDATK-BUG-020 | Energy | wrong results | Bonded exclusions miss pairs linked only through atoms outside both groups | Verified |
| MDATK-BUG-021 | Energy | wrong results | Generic LJ gives virtual sites (TIP4P EP/M) carbon parameters | Verified |
| MDATK-BUG-022 | H-bonds | wrong results | Guessed hydrogens pick their donor by non-periodic distance | Confirmed by reading |
| MDATK-BUG-023 | Base pairs | wrong results | Base pairs ignore the periodic cell (0 pairs on wrapped duplexes) | Confirmed by reading |
| MDATK-BUG-024 | Stacking | wrong results | Ring geometry from wrapped atoms: split rings never stack | Confirmed by reading |
| MDATK-BUG-025 | Proximity | wrong results | Proximity members stay `inf` when no structure atom is left | Confirmed by reading |
| MDATK-BUG-045 | GRO reader | wrong results | Wide GRO coordinates silently read from single characters (inherited from MDAnalysis) | Verified |
| MDATK-BUG-046 | GRO reader | wrong results | Multi-frame GRO: only the first frame is read, silently (inherited from MDAnalysis) | Verified |
| MDATK-BUG-047 | mdcrd reader | wrong results | One-atom mdcrd with a box: box lines read as frames (inherited from MDAnalysis) | Verified |
| MDATK-BUG-048 | PDB reader | wrong results | Partial element column: atoms get mass 0 silently (inherited from MDAnalysis) | Verified |
| MDATK-BUG-052 | Selections | wrong results | `prop 6 < x` means `x >= 6` (reversed operands negated) (inherited from MDAnalysis) | Verified |
| MDATK-BUG-053 | Selections | wrong results | Malformed `prop` operators accepted (`><`) (inherited from MDAnalysis) | Verified |
| MDATK-BUG-054 | Selections | wrong results | `atom A 1 and` silently selects nothing (inherited from MDAnalysis) | Verified |
| MDATK-BUG-057 | Rendering | wrong results | Heatmap columns all `t[1] - t[0]` wide: misplaced with uneven times | Confirmed by reading |
| MDATK-BUG-058 | Rendering | wrong results | Running mean: one NaN blanks the whole window | Verified |
| MDATK-BUG-059 | Rendering | wrong results | Class matrices with > 5 classes reuse colours | Verified |
| MDATK-BUG-060 | Density | wrong results | One NaN voxel removes the density shadow | Confirmed by reading |
| MDATK-BUG-028 | CLI | crash | Many errors escape as Python tracebacks (`info` with parm7 alone, reader errors) | Verified |
| MDATK-BUG-029 | Analyses | crash | Non-finite coordinates abort open-boundary analyses (silently dropped with PBC) | Verified |
| MDATK-BUG-030 | Density maps | crash | `density info`: MemoryError / EOFError tracebacks, no file name | Verified |
| MDATK-BUG-031 | Rendering | crash | No analysed frames: IndexError in panels, preview frame -1 | Confirmed by reading |
| MDATK-BUG-032 | Rendering | crash | A NaN coordinate makes the structure view crash (`Axis limits cannot be NaN`) | Confirmed by reading |
| MDATK-BUG-033 | Rendering | crash | Surface grid of a blown-up frame allocates hundreds of GiB | Verified |
| MDATK-BUG-034 | Charts | crash | Matrix chart of thousands of residues exceeds Agg's size limit | Confirmed by reading |
| MDATK-BUG-035 | Rendering | crash | Panel option `max_ticks = 0` → ZeroDivisionError | Verified |
| MDATK-BUG-036 | Rendering | crash | Panel options of the wrong type (`ylim = 5`, `per_row = "x"`) abort the render | Confirmed by reading |
| MDATK-BUG-037 | GUI preview | crash | Example data of a single-strand structure is corrupt (IndexError) | Verified |
| MDATK-BUG-038 | Software renderer | crash | Equal view limits → ZeroDivisionError | Confirmed by reading |
| MDATK-BUG-039 | PDB reader | crash | Hybrid-36 serials in CONECT records rejected (inherited from MDAnalysis) | Verified |
| MDATK-BUG-040 | mdcrd reader | crash | One-frame mdcrd without box: empty StopIteration (inherited from MDAnalysis) | Verified |
| MDATK-BUG-041 | NetCDF reader | crash | Cell lengths (0, 0, 0): TypeError (inherited from MDAnalysis) | Verified |
| MDATK-BUG-042 | parm7 reader | crash | Unknown %FLAG right after POINTERS: UnboundLocalError (inherited from MDAnalysis) | Verified |
| MDATK-BUG-043 | DCD reader | crash | NSAVC = 0: ZeroDivisionError (inherited from MDAnalysis) | Verified |
| MDATK-BUG-044 | Loader | crash | Trajectory chain with mixed extensions: wrong reader chosen | Verified |
| MDATK-BUG-010 | GUI | misleading message | Form errors are Python texts (`invalid literal for int()`) | Verified |
| MDATK-BUG-026 | Energy | misleading message | `use_parm7` with a non-parm7 topology silently uses generic LJ | Confirmed by reading |
| MDATK-BUG-027 | Periodic geometry | misleading message | `use_pbc` without a usable cell silently computes open distances | Confirmed by reading |
| MDATK-BUG-049 | Loader | misleading message | All reader warnings are discarded | Verified |
| MDATK-BUG-050 | Loader | misleading message | A GROMACS `.top` is read as parm7: `%VE Missing in header` | Verified |
| MDATK-BUG-051 | Readers | misleading message | Reader errors are Python internals naming MDAnalysis classes (inherited from MDAnalysis) | Verified |
| MDATK-BUG-055 | Selections | misleading message | Negative `around` radius: error or empty result depending on atom counts (inherited from MDAnalysis) | Verified |
| MDATK-BUG-056 | Selections | misleading message | Selection errors at end of input and in float ranges are Python texts (inherited from MDAnalysis) | Verified |
| MDATK-BUG-061 | Rendering | misleading message | Ion legend always says `Mg²⁺` | Confirmed by reading |
| MDATK-BUG-062 | Rendering | misleading message | Chain legend numbers chains although letters are known | Confirmed by reading |
| MDATK-BUG-063 | Charts | misleading message | A matrix key in `output.charts.series` is dropped silently | Confirmed by reading |
| MDATK-BUG-064 | Charts / render | misleading message | Bad chart formats and colours fail after the analyses, without naming the setting | Verified |
| MDATK-BUG-065 | CLI | misleading message | KeyErrors printed with Python's quotes | Verified |
| MDATK-BUG-066 | CLI | misleading message | `--frames x`, `--set groups.--1.name=x`: `invalid literal for int()` | Verified |
| MDATK-BUG-067 | RMSD | misleading message | RMSD reference frame beyond the run clamped silently | Verified |
| MDATK-BUG-068 | GUI | minor | Quitting during a job kills the worker mid-write | Confirmed by reading |
| MDATK-BUG-069 | User settings | minor | `settings.toml` written non-atomically; an unreadable file is replaced by an empty one | Confirmed by reading |
| MDATK-BUG-070 | Energy | minor | Net force direction inconsistent with the 0.5 Å distance clamp | Confirmed by reading |
| MDATK-BUG-071 | Atoms analysis | minor | Bond guess misses S–S (second half withdrawn 2026-11-27) | Confirmed by reading |
| MDATK-BUG-072 | GRO reader | minor | GRO box fallback loses minus signs (inherited from MDAnalysis) | Verified |
| MDATK-BUG-073 | Selections | minor | Duplicate group names: `group:` uses the last, `Project.group()` the first | Verified |
| MDATK-BUG-074 | Selections | minor | `group:` token stops at invalid characters, misleading error | Verified |
| MDATK-BUG-075 | Rendering | minor | `wrap_degrees` folds every angle within 0.0018° of −180 to +180 | Verified |
| MDATK-BUG-076 | Rendering | minor | Secondary panel: unknown `source` silently means `dynamic` | Confirmed by reading |
| MDATK-BUG-077 | Video | minor | Codec fallback ignores `crf` (fixed mpeg4 quality, no rate for openh264) | Confirmed by reading |
| MDATK-BUG-078 | Video | minor | Videos encoded with BT.601 and left untagged: shifted hues | Confirmed by reading |
| MDATK-BUG-079 | GUI preview | minor | Example data expands whole residue ranges (8 s for `1-300000000`) | Verified |
| MDATK-BUG-080 | Software renderer | minor | `fast` quality at odd sizes: blank last row/column, view scaled | Confirmed by reading |
| MDATK-BUG-081 | Software renderer | minor | Glow halo seam at the molecule box; halo strength depends on earlier renders | Confirmed by reading |
| MDATK-BUG-082 | Charts | minor | Chart title and colour bar stay black on the dark theme | Confirmed by reading |
| MDATK-BUG-083 | CLI | minor | `stop=None` in the empty frame-selection message | Verified |
| MDATK-BUG-084 | Pipeline | minor | External data files are read only after all analyses ran | Confirmed by reading |
| MDATK-BUG-085 | PDB export | minor | A frame beyond the PDB limits aborts the export, leaving earlier files | Confirmed by reading |

---------------------------------------------------------------------------------------------------

## 1. Data loss on save and in the GUI

Reproduction for 001–008 (one probe; run headless with `xvfb-run -a python probe.py`):

```python
import os, tempfile, gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from mdatk.config import Project, PanelSpec, LayoutParams
from mdatk.gui.app import MainWindow
app = Gtk.Application(application_id=None); app.register(None)
d = tempfile.mkdtemp(); src = os.path.join(d, "p.mdatk.toml")
p = Project()
p.analyses.contacts.cutoff = 25.0; p.analyses.hbonds.angle_cutoff = 60.0; p.analyses.energy.cutoff = 60.0
p.inputs.dt_ps = 0.0123456789; p.analyses.energy.target_value = -1234567.5
p.visual.glow_color = "tab:orange"; p.visual.proximity_color = "#ff000080"; p.visual.base_color = "gold"
p.output.charts.formats = ["jpg", "svg"]; p.output.video.width = 7682
p.save(src); q = Project.load(src)
win = MainWindow(app, q, src); r = win.collect(validate=False)
for name, a, b in [("contacts.cutoff", q.analyses.contacts.cutoff, r.analyses.contacts.cutoff),
                   ("hbonds.angle_cutoff", q.analyses.hbonds.angle_cutoff, r.analyses.hbonds.angle_cutoff),
                   ("energy.cutoff", q.analyses.energy.cutoff, r.analyses.energy.cutoff),
                   ("inputs.dt_ps", q.inputs.dt_ps, r.inputs.dt_ps),
                   ("energy.target_value", q.analyses.energy.target_value, r.analyses.energy.target_value),
                   ("glow_color", q.visual.glow_color, r.visual.glow_color),
                   ("proximity_color", q.visual.proximity_color, r.visual.proximity_color),
                   ("base_color", q.visual.base_color, r.visual.base_color),
                   ("charts.formats", q.output.charts.formats, r.output.charts.formats),
                   ("video.width", q.output.video.width, r.output.video.width)]:
    print(("SAME " if a == b else "DIFF ") + f"{name}: {a!r} -> {b!r}")
q2 = Project.load(src); q2.visual.renderer = "vulkan"          # any value not in the combo
win.load_project(q2, src); print("renderer:", win.collect(validate=False).visual.renderer)
p = Project(); p.layout = LayoutParams(preset="equal_4", panels=[PanelSpec("structure", [], "Structure"),
                                                                 PanelSpec("timeseries", ["rmsd/Target"], "RMSD")])
p.save(src); win.load_project(Project.load(src), src)
print("collect:", [(s.type, s.series, s.title) for s in win.collect(validate=False).layout.panels])
win.preset_combo.set_active_id("single"); win.preset_combo.set_active_id("equal_4")
print("after single->equal_4:", [(s.type, s.series) for s in win.collect(validate=False).layout.panels])
```

Observed at 0.1.1:

```
DIFF contacts.cutoff: 25.0 -> 20.0
DIFF hbonds.angle_cutoff: 60.0 -> 90.0
DIFF energy.cutoff: 60.0 -> 50.0
DIFF inputs.dt_ps: 0.0123456789 -> 0.0123457
DIFF energy.target_value: -1234567.5 -> -1234570.0
DIFF glow_color: 'tab:orange' -> '#000000'
DIFF proximity_color: '#ff000080' -> '#000000'
DIFF base_color: 'gold' -> '#ffd700'
DIFF charts.formats: ['jpg', 'svg'] -> ['svg']
DIFF video.width: 7682 -> 7680
renderer: software
collect: [('structure', [], 'Structure'), ('timeseries', ['rmsd/Target'], 'RMSD'),
          ('timeseries', ['contacts/atom_pairs'], ''), ('timeseries', ['contacts/atom_pairs'], '')]
after single->equal_4: [('structure', []), ('timeseries', ['contacts/atom_pairs']),
                        ('timeseries', ['contacts/atom_pairs']), ('timeseries', ['contacts/atom_pairs'])]
```

### MDATK-BUG-001 — Spin buttons clamp loaded values; saving rewrites valid settings

- **MDBridge ref:** AP-01 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/gui/widgets.py:81-87` (`spin`: the value is clamped to the adjustment), `:373-378`
  (`spin_float` / `spin_int`); ranges in `mdatk/gui/app.py:1092` (contacts 0.5–20), `:1103` (H-bond angle
  90–180), `:1168` (energy cutoff 4–50), `:1793` (video width 64–7680).
- **Description:** loading a project puts each value into a `Gtk.SpinButton` whose range is narrower than
  what `validate()` allows; GTK clamps, and *Save* writes the clamped value. Opening and saving a project
  silently changes its science settings.
- **Reproduce:** probe above.
- **Expected:** a value that `validate()` accepts survives open + save unchanged.
- **Actual:** `contacts.cutoff 25 → 20`, `hbonds.angle_cutoff 60 → 90` (validate allows 0–180),
  `energy.cutoff 60 → 50`, `video.width 7682 → 7680`.
- **Suggested fix:** MDBridge's form fields hold any value of the setting's type; the ranges only come from
  validation, which flags an out-of-range value in the form but never rewrites it. In GTK: set the
  adjustment's bounds wide (or to the validation limits) and flag values outside the suggested range.
- **Test idea:** for every spin-bound setting, load a project with a valid value outside the widget's
  current range, `collect()`, and assert the value is unchanged.

### MDATK-BUG-002 — Optional numbers are rounded to 6 significant digits on save

- **MDBridge ref:** AP-02 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Where:** `mdatk/gui/widgets.py:403-407` (`optional_float` shows `f"{v:g}"` and reads the text back);
  the same `%g` in `mdatk/gui/wizard.py` (MDBridge M7).
- **Description:** `%g` keeps 6 significant digits; the shown text is what `collect()` parses.
- **Actual:** `inputs.dt_ps 0.0123456789 → 0.0123457`, `energy.target_value −1234567.5 → −1234570`.
- **Expected:** unchanged.
- **Suggested fix:** show the shortest round-trip text (`repr(float)`), as the project writer does.
- **Test idea:** round-trip random doubles through every `optional_float` binding.

### MDATK-BUG-003 — Colour buttons turn valid colours black or rewrite them

- **MDBridge ref:** AP-03 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Where:** `mdatk/gui/widgets.py:38-49` (`hex_to_rgba` falls back to `#000000` when GDK cannot parse),
  `:393-395` (`color` binding writes `rgba_to_hex`).
- **Description:** colours are parsed with GDK, not Matplotlib. Matplotlib colours MDATK itself accepts
  (`tab:orange`, `#rrggbbaa`, `C0`, `0.5` grey) become black; names GDK knows are rewritten.
- **Actual:** `tab:orange → #000000`, `#ff000080 → #000000`, `gold → #ffd700`.
- **Expected:** the stored text is kept unless the user picks a colour.
- **Suggested fix:** keep the setting's text in the field; only the picker writes `#rrggbb`; validate the
  text with `matplotlib.colors.is_color_like` and flag invalid entries.
- **Test idea:** load each of `tab:orange`, `C3`, `#ff000080`, `0.5`, `gold`; collect; assert unchanged.

### MDATK-BUG-004 — Chart formats other than png / svg / pdf are dropped on save

- **MDBridge ref:** AP-04 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Where:** `mdatk/gui/app.py:1811-1816` (three check boxes; nothing ticked → `["png"]`).
- **Actual:** `formats = ["jpg", "svg"]` is saved as `["svg"]`.
- **Expected:** formats the charts support but the GUI does not show are kept.
- **Suggested fix:** MDBridge offers every supported format and keeps (and flags) unknown entries. At
  minimum, carry over entries that have no check box.
- **Test idea:** load `["jpg", "svg", "tiff"]`, collect, assert equal.

### MDATK-BUG-005 — Combo boxes replace an unknown value with the previous project's choice

- **MDBridge ref:** AP-05 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Where:** `mdatk/gui/widgets.py:384-386` (`combo`: `set_active_id` of an unknown id leaves the old
  selection; the getter returns it, or `""` the first time).
- **Description:** a value not in the list (a newer option, a typo, a hand edit) is silently replaced with
  whatever the previous project had selected, or with `""` in a fresh window.
- **Actual:** project with `visual.renderer = "vulkan"` opened after one with `software` → saved as `software`.
- **Expected:** the unknown value stays, shown as e.g. `vulkan (unknown)`, and the check flags it.
- **Suggested fix:** when `set_active_id` returns False, append a temporary entry for the value (MDBridge
  shows `<value> (unknown)`), never substitute.
- **Test idea:** load A (`software`), then B (`unknownvalue`); collect B; assert `unknownvalue`.

### MDATK-BUG-006 — GUI save makes every path absolute

- **MDBridge ref:** AP-06 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Where:** the GUI loads with path resolution (`mdatk/config.py:557-578`, `resolve`) and saves the
  resolved project (`mdatk/gui/app.py:1970`); the CLI's `init` keeps relative paths
  (`mdatk/cli.py:228`, `resolve=False`).
- **Reproduce:** open `demo/demo.mdatk.toml` in the GUI and save.
- **Actual (diff of the file):** `topology = "demo.parm7"` → `"/tmp/…/demo/demo.parm7"`,
  `trajectories = ["demo.nc"]` → `["/tmp/…/demo/demo.nc"]`, `out_dir = "MDATK_output"` →
  `"/tmp/…/demo/MDATK_output"`.
- **Expected:** paths written as the file had them. A project folder that is moved or shared breaks after
  one GUI save.
- **Suggested fix:** MDBridge writes paths as they were given (relative to the project folder); *Save as*
  to another folder rebases relative paths. In MDATK: keep the unresolved strings next to the resolved
  ones, or `os.path.relpath` paths under the project folder on save.
- **Test idea:** open + save a project with relative paths; assert the text of those keys is unchanged.

### MDATK-BUG-007 — Empty layout slots are saved as real panels

- **MDBridge ref:** AP-07 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Where:** `mdatk/gui/app.py:1721-1722` (`_rebuild_slots`: a missing slot becomes
  `PanelSpec("timeseries", [], "")`, whose editor then picks the first data key).
- **Actual:** a 4-slot preset with 2 panels is saved with 4 (two untitled `contacts/atom_pairs` series).
- **Expected:** `layout.panels` keeps its length; empty slots stay empty.
- **Suggested fix:** keep empty slots as placeholders (MDBridge draws them as such) and save only panels
  the user configured.
- **Test idea:** probe above, assert two panels after `collect()`.

### MDATK-BUG-008 — Switching to a smaller preset deletes panels for good

- **MDBridge ref:** AP-08 · **Area:** GUI save · **Severity:** critical data loss · **Status:** Verified
- **Where:** `mdatk/gui/app.py:1562-1566` → `:1710-1725` (`_rebuild_slots` keeps only `n_slots(preset)`
  editors; `_read_slots` reads only those).
- **Actual:** `equal_4 → single → equal_4` loses the RMSD panel (see probe output).
- **Expected:** panels beyond the slot count are kept until the user removes them (the renderer already
  ignores extra panels with a warning).
- **Suggested fix:** keep the full spec list in the window and show only the first N; MDBridge does this.
- **Test idea:** probe above, assert panel 2 is still `rmsd/Target`.

### MDATK-BUG-009 — New / Open / Close / quit discard unsaved edits without asking

- **MDBridge ref:** AP-13 · **Area:** GUI · **Severity:** critical data loss · **Status:** Confirmed by reading
- **Where:** `mdatk/gui/app.py:1936-1946` (`new_project`, `close_project`), `:1953-1968` (`open_path`),
  `:2551-2553` (`_on_close` only stops jobs and returns False).
- **Description:** the window knows about unsaved edits (it uses that for co-editing reloads, `:2078-2110`),
  but none of these actions asks; nothing marks a modified project in the title.
- **Expected:** a modified marker and a Save / Discard / Cancel question.
- **Suggested fix:** compare `collect(validate=False)` with the last loaded/saved project (MDBridge keeps a
  modified flag); ask before replacing it.
- **Test idea:** change a field, call `new_project()` with a stubbed dialog; assert the dialog was shown.

### MDATK-BUG-010 — Form errors are Python texts

- **MDBridge ref:** AP-10 · **Area:** GUI · **Severity:** misleading message · **Status:** Verified
- **Where:** `mdatk/gui/widgets.py:397-401` (`optional_int`: `int(t)`); panel options JSON in
  `mdatk/gui/app.py` `_read_slots`.
- **Reproduce:** in the probe, set `inputs.stop` to `x`, or a slot's options to `{bad`, then `collect()`.
- **Actual:** `inputs.stop: invalid literal for int() with base 10: 'x'`;
  `layout.panels.0.options: panel options '{bad': Expecting property name enclosed in double quotes: line 1 column 2 (char 1)`.
- **Expected:** a located message in words, e.g. `Input & cleaning › Frames › Stop: 'x' is not a whole
  number` (0.0.10's `problems.py` already locates selection problems; extend it to every field).
- **Test idea:** type a non-number into each int/float field; assert the message names the field and
  contains no `invalid literal`.

---------------------------------------------------------------------------------------------------

## 2. Wrong numerical results

### MDATK-BUG-011 — The element of an atom depends on which group it is evaluated in

- **MDBridge ref:** SEL-13 · **Area:** analyses (H-bonds, Mg, contacts, energy, atoms) · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/selections.py:360-378` (`atom_elements`: if *any* atom of the group has an empty
  element, every atom's element is guessed from its name), `:345-357` (`guess_element`: first letter).
- **Reproduce:**
  ```python
  import MDAnalysis as mda
  from mdatk.selections import atom_elements
  u = mda.Universe.empty(5, n_residues=2, atom_resindex=[0, 0, 0, 0, 1], trajectory=True)
  u.add_TopologyAttr("names", ["C1", "CL1", "BR1", "MG", "EPW"]); u.add_TopologyAttr("resnames", ["LIG", "WAT"])
  u.add_TopologyAttr("elements", ["C", "Cl", "Br", "Mg", ""])
  print(list(atom_elements(u.select_atoms("resname LIG"))), list(atom_elements(u.atoms)))
  ```
- **Actual:** `['C', 'CL', 'BR', 'MG']` for the ligand alone, `['C', 'C', 'B', 'M', 'E']` when the group also
  holds one atom without an element (e.g. an Amber extra point with Z = 0).
- **Expected:** the same element for an atom wherever it is evaluated.
- **Suggested fix:** fall back per atom: use the topology element when it is non-empty, otherwise the
  guess for that atom alone. Bump the versions of analyses that use elements so caches recompute.
- **Test idea:** the snippet above; assert both lists agree on the shared atoms.

### MDATK-BUG-012 — NaN / negative / zero cutoffs, radii and dielectric are accepted

- **MDBridge ref:** PRJ-01 · **Area:** project validation · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/config.py:741-744` (`validate()` checks only `contacts.cutoff <= 0` and Mg inner >
  outer; `nan <= 0` is False).
- **Reproduce:**
  ```
  mdatk analyze -p demo/demo.mdatk.toml -o o1 --only hbonds --set analyses.hbonds.distance_cutoff=-3.5 -q
  sed 's/^distance_cutoff = 3.5/distance_cutoff = nan/' demo/demo.mdatk.toml > demo/nan.mdatk.toml
  mdatk analyze -p demo/nan.mdatk.toml -o o2 --only hbonds -q
  ```
  then read `count` from `o*/MDATK_cache/hbonds-*.npz`.
- **Actual:** cutoff −3.5 → `[30, 31, 21, 21, 31, 31]` (identical to +3.5: the open path uses |cutoff|,
  the periodic path gives no pairs or an error); `nan` → `[0, 0, 0, 0, 0, 0]`. No message in either case.
  The same holds for energy / Mg / stacking cutoffs, `energy.dielectric = 0` (infinite energies), angle
  cutoffs outside 0–180, `proximity_search_radius <= 0`.
- **Expected:** validation errors.
- **Suggested fix:** in `validate()`: every cutoff, radius and tolerance finite and > 0; angles in
  [0, 180]; dielectric > 0; Mg inner ≤ outer with both finite. Messages in the existing style
  (`H-bond distance cutoff must be positive.`, `… must be finite.`).
- **Test idea:** parametrise over each numeric field with `nan`, `inf`, `-1`, `0`.

### MDATK-BUG-013 — Two RMSD targets with the same data key: one series silently lost

- **MDBridge ref:** PRJ-03 · **Area:** RMSD · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/analysis/rmsd.py:14-16` (`_out_key` maps runs of other characters to `_`), `:116-119`
  (`add_series` under the same key overwrites).
- **Reproduce:**
  ```python
  from mdatk.config import Project, RMSDTarget
  from mdatk.analysis.rmsd import RMSDAnalysis
  p = Project.load("demo/demo.mdatk.toml")
  p.analyses.rmsd.targets = [RMSDTarget("A B", "group:Aptamer and name P"), RMSDTarget("A_B", "group:Target and name P")]
  print(p.validate(), RMSDAnalysis.expected_outputs(p)["series"])
  ```
- **Actual:** `[] ['A_B', 'A_B']`.
- **Expected:** `validate()` reports `RMSD targets 1 and 2 give the same data key 'A_B'; rename one.`
- **Suggested fix:** that validation (MDBridge's choice); alternatively disambiguate keys as 0.1.0 did for
  external data (`-2` suffix), but a silent rename changes panel keys users wrote.
- **Test idea:** the snippet; assert one problem.

### MDATK-BUG-014 — `inputs.dt_ps` negative or NaN gives a negative / NaN time axis

- **MDBridge ref:** LD-04 · **Area:** loader / runner · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/analysis/base.py:283-286` and `mdatk/pipeline.py:235-237` (`if dt_ps:` is true for NaN
  and negatives), `mdatk/io/loader.py:61-62, 70-71`; nothing in `validate()`.
- **Reproduce:** `mdatk analyze -p demo/demo.mdatk.toml -o o5 --only hbonds --set inputs.dt_ps=-2 -q`, then
  read `_time_ps` from the cache `.npz`.
- **Actual:** `_time_ps: [-0.0, -2.0, -4.0, -6.0, -8.0, -10.0]`; with NaN every time is NaN.
- **Expected:** `inputs.dt_ps must be positive` (0 / unset keeps meaning "from the file").
- **Suggested fix:** the validation above.
- **Test idea:** `validate()` with `-2`, `nan`, `inf`.

### MDATK-BUG-015 — RMSD of one or two atoms is NaN or off by up to 8 Å (inherited from MDAnalysis)

- **MDBridge ref:** AN-16, AN-17 · **Area:** RMSD · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/analysis/rmsd.py:99-101` (`rms.rmsd(..., superposition=True)` when no fit selection).
  Reached by the default *Backbone* target (`name C4'`, `rmsd.py:33`) of a 1–2 residue group.
- **Description (behaviour of MDAnalysis' QCP superposition):** for 1 atom the result is NaN; for 2 atoms
  (rigid copies, even identical sets, of float32-rounded coordinates) the iteration does not converge in
  many cases and returns NaN or values of several Å that change with the last input bit. Identical
  structures of ≥ 3 atoms can give up to ~5e-6 Å instead of 0.
- **Reproduce:**
  ```python
  import numpy as np
  from MDAnalysis.analysis import rms
  from scipy.spatial.transform import Rotation
  rng = np.random.default_rng(3); f32 = lambda a: np.asarray(a, np.float32).astype(np.float64)
  print(rms.rmsd(np.array([[1., 2, 3]]), np.array([[4., 5, 6]]), center=True, superposition=True))
  for fam in ("rotated copy", "identical"):
      vals = []
      for t in range(100):
          a = f32(rng.normal(size=(2, 3)) * 5)
          b = a.copy() if fam == "identical" else f32(Rotation.random(random_state=t).apply(a) + 2)
          vals.append(rms.rmsd(a, b, center=True, superposition=True))
      v = np.array(vals); print(fam, np.isnan(v).sum(), (v > 1e-3).sum(), np.nanmax(v))
  ```
- **Actual:** `nan` for one atom; 2 atoms, exact answer 0: *rotated copy* 24 NaN and 11 values > 1e-3 Å
  (worst 6.01 Å) out of 100; *identical* 34 NaN and 14 > 1e-3 Å (worst 8.12 Å).
- **Expected:** 0 for one atom; the exact superposition RMSD otherwise.
- **Suggested fix (in MDATK's code):** detect degenerate input (n < 3, or nearly rank-deficient
  inner-product matrix) and compute the RMSD yourself: n = 1 → 0; otherwise superpose with
  `mdatk.spatial.kabsch` (SVD) and compute the residual after rotating. MDBridge uses an exact eigen
  solution there and keeps the QCP value elsewhere; it also returns exactly 0 when E0 − λ is at the
  rounding level (recomputing the residual with the optimal rotation).
- **Test idea:** the families above (plus collinear sets) against an SVD reference within 1e-6 Å.

### MDATK-BUG-016 — Fit selections of fewer than 3 or collinear atoms: no rotation or an arbitrary one

- **MDBridge ref:** AN-18 · **Area:** RMSD (and the trace/atoms alignment) · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/spatial.py:89-90` (`len(mobile) < 3` → identity rotation, silently); `:91-97` for
  collinear sets the rotation about the line is whatever LAPACK returns.
- **Reproduce:** `kabsch(np.array([[0,0,0],[1,0,0.]]), np.array([[0,0,0],[0,1,0.]]))[0]`
- **Actual:** the identity matrix, although the two sets differ by a 90° rotation; no warning.
- **Expected:** a warning naming the target, and a deterministic result.
- **Suggested fix:** MDBridge warns once (`RMSD target 'X': fit selection has N atoms (collinear);
  superposition is translation only`) and uses translation only for both cases, so results do not
  depend on the platform's LAPACK.
- **Test idea:** 2-atom and collinear fits, assert the warning and identical results across runs.

### MDATK-BUG-017 — Periodic cutoff decisions depend on the box size; far points get wrong distances (inherited from MDAnalysis)

- **MDBridge ref:** AN-01, AN-02, AN-03 · **Area:** periodic geometry (contacts, H-bonds, Mg, stacking, energy, proximity) · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/spatial.py:36-43` (`capped_distance` on float32 input whenever a box is given).
- **Description:** with a box, each distance is computed with float32 scale factors, so a pair exactly at
  the cutoff is kept or dropped depending on the box length, and points far outside the primary cell get
  errors of 1e-4 Å and more. Which internal method runs (and whether `Cutoff … too large for box` or
  `Cutoff must be positive` is raised) depends on the atom counts. The pair order differs between the
  open and periodic paths.
- **Reproduce:**
  ```python
  import numpy as np; from mdatk.spatial import pairs_within
  a, b = np.array([[1., 1, 1]]), np.array([[4.5, 1, 1]])
  for L in (25, 30, 31.7, 40, 60): print(L, len(pairs_within(a, b, 3.5, box=[L, L, L, 90, 90, 90])[0]))
  print(pairs_within(np.array([[10000., 0, 0]]), np.array([[0., 0, 0]]), 5.0, box=[20, 20, 20, 90, 90, 90])[2])
  ```
- **Actual:** pair at exactly 3.5 Å kept for L = 25 and 31.7, dropped for L = 30, 40, 60; the point at
  10000 Å (500 box lengths, distance 0) gives `1.49e-04`.
- **Expected:** one decision independent of the box (`d² <= c²` on exact minimum-image distances, as the
  open path does).
- **Suggested fix (MDATK side):** compute periodic pairs yourself in float64: wrap each difference with
  `δ − L·round(δ/L)` (orthorhombic) or the shortest of the 27 neighbour images after reduction
  (triclinic), use a cell list for speed, sort pairs by (i, j). This also fixes MDATK-BUG-029's periodic
  half and removes the size-dependent errors.
- **Test idea:** the boundary pair over many box lengths; random far-outside points against a float64
  reference.

### MDATK-BUG-018 — Selections are evaluated on whatever frame is current at setup, and always periodically

- **MDBridge ref:** SEL-17, AN-13 · **Area:** analyses · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/base.py:273-276` (all `setup()` calls run before the frame loop, with no seek);
  `mdatk/analysis/rmsd.py:85` (RMSD's setup seeks to its reference frame and leaves it current);
  `mdatk/analysis/structure.py:344` (`atom_bonds` guesses bonds from the current frame);
  `mdatk/selections.py:326-330` (`select_atoms` uses the universe's box regardless of `use_pbc`).
- **Description:** distance keywords (`around`, `sphzone`, …) in group selections and the atoms analysis'
  bond guess use frame 0 of the trajectory, or the RMSD reference frame for analyses set up after RMSD
  (trace, proximity, atoms, base pairs). They also apply the cell even with `analyses.use_pbc = false`.
- **Expected:** evaluated on the first analysed frame (`frame_indices[0]`), whatever the setup order, and
  honouring `use_pbc`.
- **Suggested fix:** seek to `frame_indices[0]` before each `setup()` (MDBridge passes that frame to
  every analysis); for `use_pbc = false`, evaluate distance keywords without the box (e.g. clear
  `ts.dimensions` temporarily, or `select_atoms(..., periodic=False)`).
- **Test idea:** a group `around 5 resname MG` on a trajectory where the answer changes between frame 0
  and the start frame; with and without RMSD enabled.

### MDATK-BUG-019 — Energy: no bonded exclusions, silently, when the topology has no bonds

- **MDBridge ref:** AN-07 · **Area:** energy · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/energy.py:36-41` (`_bonded_exclusions` returns nothing when `u.bonds` is
  missing or empty).
- **Description:** a PDB without CONECT or a GRO file: covalently linked atoms of two groups (e.g. two
  halves of one strand) are not excluded and give huge LJ terms (with the 0.5 Å clamp), without a
  message.
- **Expected:** at least `energy: the topology has no bonds; bonded pairs are not excluded`.
- **Suggested fix:** MDBridge warns and excludes along distance-guessed bonds (covalent radii, see
  MDATK-BUG-071) when `exclude_bonded > 0` and the topology has none.
- **Test idea:** energy between residues 1–6 and 7–12 of one strand from a PDB without CONECT; assert the
  warning and a finite vdW.

### MDATK-BUG-020 — Energy: exclusions miss pairs linked only through atoms outside both groups

- **MDBridge ref:** AN-08 · **Area:** energy · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/analysis/energy.py:46-48, 53` (search starts only from A atoms in a direct A–B bond).
- **Reproduce:**
  ```python
  import MDAnalysis as mda; from mdatk.analysis.energy import _bonded_exclusions
  u = mda.Universe.empty(3, n_residues=3, atom_resindex=[0, 1, 2], trajectory=True); u.add_bonds([(0, 1), (1, 2)])
  print(_bonded_exclusions(u, u.atoms[[0]], u.atoms[[2]], 3).tolist())
  ```
- **Actual:** `[]` — A–X–B (two bonds apart) is not excluded with `exclude_bonded = 3`. With a direct A–B
  bond it is (`[0]`).
- **Expected:** every A–B pair within `exclude_bonded` bonds along any path is excluded.
- **Suggested fix:** BFS from every A atom over the whole bond graph up to `depth`, collecting B atoms
  (MDBridge does this).
- **Test idea:** the snippet; also a linker residue between the groups.

### MDATK-BUG-021 — Generic Lennard-Jones gives virtual sites (TIP4P EP / M) carbon parameters

- **MDBridge ref:** AN-26 (part) · **Area:** energy (non-parm7 topologies) · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/analysis/generic_lj.py:97-108, 166-185` (`lookup` of an unknown element → `DEFAULT`,
  carbon); `mdatk/selections.py:345-357` (`EPW` → element `E`, see MDATK-BUG-011); also
  `mdatk/analysis/energy.py:167-168` passes `ts.dimensions` to the hydrogen search regardless of
  `use_pbc`.
- **Description:** a TIP4P water's extra point (names `EP`, `EPW`, `M`, `LP`; Z = 0) is a massless charge
  site with no LJ sphere, but gets R* 1.908 Å / ε 0.086 with the warning `no generic Lennard-Jones
  parameters for element(s) E: carbon's are used`.
- **Expected:** extra points get (0, 0) and no warning.
- **Suggested fix:** recognise extra points (empty element with Z ≤ 0 in the topology, or the usual names
  in a water residue) and give them zero LJ, as MDBridge does (type `EP`). Use `self.ctx.box(ts)` for the
  hydrogen search so it follows `use_pbc`.
- **Test idea:** a TIP4P PDB (O, H1, H2, EPW): `parameters()` gives 0 for EPW and lists no missing element.

### MDATK-BUG-022 — H-bonds: guessed hydrogens pick their donor by non-periodic distance

- **MDBridge ref:** AN-09 · **Area:** H-bonds (topologies without bonds, `use_pbc`) · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/hbonds.py:141-150` (candidates within 1.25 Å come from the periodic pair
  search, then the nearest is chosen with `np.linalg.norm` of the raw difference; the search's own
  distances are discarded).
- **Description:** a hydrogen across a face from its donor loses to any non-periodic candidate (or keeps
  a far-away "nearest").
- **Expected:** the choice uses the same minimum-image distance as the search.
- **Suggested fix:** use the `d` returned by `pairs_within` (ties: lowest donor index).
- **Test idea:** an N–H split across the box face with another heavy atom 1.2 Å away in plain distance.

### MDATK-BUG-023 — Base pairs ignore the periodic cell

- **MDBridge ref:** AN-11 · **Area:** base pairs · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/basepairs.py:102-114` (the detector gets raw positions; no box is used).
- **Description:** a duplex whose strands sit on opposite sides of a face (each strand imaged as a
  molecule, which `autoimage` does) has 0 base pairs per frame, with `use_pbc` on or off.
- **Expected:** with `use_pbc`, pairs found across the face.
- **Suggested fix:** with `use_pbc` and a usable cell, compute the detector's distances with minimum
  images (MDBridge), or make each residue pair whole before detection.
- **Test idea:** the demo duplex with strand B shifted by one box vector; same pairs as unshifted.

### MDATK-BUG-024 — Stacking: rings split by wrapping never stack

- **MDBridge ref:** AN-12 · **Area:** stacking · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/stacking.py:146-158` (`ring_geometry(pos[idx])` on raw positions; only the
  centroid offsets later use minimum images).
- **Description:** on per-atom-wrapped trajectories a ring crossing a face has a centroid inside the box
  and is "non-planar", so it is missed.
- **Expected:** with `use_pbc`, ring atoms made whole first.
- **Suggested fix:** unwrap each ring to the minimum image of its first atom before the geometry.
- **Test idea:** a stacked pair with one ring split across a face.

### MDATK-BUG-025 — Proximity: members stay `inf` when no structure atom is left

- **MDBridge ref:** AN-15 · **Area:** proximity · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/structure.py:288-297` (`process` returns before `self.dist[self.member, k] = 0.0`
  when `self.struct` is empty).
- **Description:** when every structure-view residue is in the proximity set, the members' distance stays
  `inf` instead of 0, and the proximity colouring treats them as far away.
- **Suggested fix:** set the members to 0 before the early return.
- **Test idea:** proximity set = all structure residues; assert all zeros.

### MDATK-BUG-026 — `use_parm7` with a non-parm7 topology silently uses generic LJ

- **MDBridge ref:** AN-06 · **Area:** energy · **Severity:** misleading message · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/energy.py:143-156` (the warning is given only when reading a parm7 fails).
- **Expected:** `energy: use_parm7 is set but the topology <file> is not a parm7 file; using generic
  element parameters`.
- **Test idea:** energy on a PDB with `use_parm7 = true`; assert the warning.

### MDATK-BUG-027 — `use_pbc` without a usable cell silently computes open distances

- **MDBridge ref:** AN-19 · **Area:** periodic geometry · **Severity:** misleading message · **Status:** Confirmed by reading
- **Where:** `mdatk/spatial.py:21-27` (`_usable_box` returns None for a missing cell or a zero length);
  `mdatk/analysis/base.py:137-140`.
- **Expected:** one warning per run: `use_pbc is on but <file> frame k has no usable periodic cell;
  distances are not periodic`.

---------------------------------------------------------------------------------------------------

## 3. Crashes

### MDATK-BUG-028 — Many errors escape as Python tracebacks

- **MDBridge ref:** CLI-05, LD-06 · **Area:** CLI · **Severity:** crash · **Status:** Verified
- **Where:** `mdatk/cli.py:1641` catches only `OSError, ValueError, RuntimeError, KeyError`; reader and
  loader errors of other classes (`AttributeError`, `StopIteration`, `TypeError`, `IndexError`,
  `UnboundLocalError`, `EOFError`, `ZeroDivisionError`, `MemoryError`) print a traceback.
- **Reproduce:**
  ```
  mdatk info --top demo/demo.parm7
  ```
  and the reader cases of section 4 (e.g. `mdatk info --top no_box.gro`).
- **Actual:** `AttributeError: No trajectory loaded into Universe. Did you mean: '_trajectory'?` (traceback);
  `UnboundLocalError: cannot access local variable 'unitcell' …` (traceback), and so on.
- **Expected:** `error: demo/demo.parm7 has no coordinates; give a trajectory with --traj`, and
  `error: <file>: <message>` for every reader failure, exit status 1.
- **Suggested fix:** in `main()`, catch `Exception` after the specific handlers (keep `MDATK_DEBUG`), and
  wrap `load_universe` so its errors name the file. For `info` without a trajectory, check
  `hasattr(u, "trajectory")` first.
- **Test idea:** run each failing reader case through `main()`; assert exit 1 and no `Traceback`.

### MDATK-BUG-029 — Non-finite coordinates abort open-boundary analyses (and are silently dropped with PBC)

- **MDBridge ref:** AN-20 · **Area:** analyses · **Severity:** crash · **Status:** Verified
- **Where:** `mdatk/spatial.py:44-46` (scipy `cKDTree` raises), `:36-43` (periodic path drops the point).
- **Reproduce:**
  ```python
  import numpy as np; from mdatk.spatial import pairs_within
  a = np.array([[0, 0, 0], [np.nan, 0, 0]], float); b = np.array([[1.0, 0, 0]])
  pairs_within(a, b, 3.0)                                   # open
  print(pairs_within(a, b, 3.0, box=[20, 20, 20, 90, 90, 90]))
  ```
- **Actual:** open: `ValueError: data must be finite, check for nan or inf values` (the whole analysis
  aborts); periodic: the NaN atom is silently in no pair.
- **Expected:** the same behaviour on both paths: non-finite points take part in no pair, with one
  warning `frame k: N atoms have non-finite coordinates`.
- **Suggested fix:** mask non-finite rows before building the trees and map indices back.
- **Test idea:** the snippet; a blown-up frame in a short trajectory.

### MDATK-BUG-030 — `density info`: MemoryError / EOFError tracebacks without the file name

- **MDBridge ref:** MRC-1, MRC-2, CLI-04 · **Area:** density maps · **Severity:** crash · **Status:** Verified
- **Where:** `mdatk/io/mrc.py` (header sizes are trusted; gzip errors are not wrapped).
- **Reproduce:**
  ```python
  import struct, gzip
  h = bytearray(1024)
  struct.pack_into("<3i", h, 0, 100000, 100000, 100000); struct.pack_into("<i", h, 12, 2)
  struct.pack_into("<3i", h, 28, 100000, 100000, 100000); struct.pack_into("<6f", h, 40, 10, 10, 10, 90, 90, 90)
  struct.pack_into("<3i", h, 64, 1, 2, 3); h[208:212] = b"MAP "; h[212:216] = bytes([0x44, 0x44, 0, 0])
  open("huge.mrc", "wb").write(bytes(h) + b"\0" * 16)
  struct.pack_into("<3i", h, 0, 4, 4, 4); struct.pack_into("<3i", h, 28, 4, 4, 4)
  d = gzip.compress(bytes(h) + b"\0" * 256); open("trunc.mrc.gz", "wb").write(d[:-30])
  ```
  `mdatk density info huge.mrc`; `mdatk density info trunc.mrc.gz`
- **Actual:** traceback ending `MemoryError`; traceback ending
  `EOFError: Compressed file ended before the end-of-stream marker was reached`.
- **Expected:** `error: huge.mrc: file is shorter than its header says` (checked before allocating);
  `error: trunc.mrc.gz: …` naming the file for every gzip / zlib error. When MODE or the cell angles only
  make sense byte-swapped, say `cannot determine the byte order of the map header (mode …)` instead of
  `unsupported MRC mode 33554432`.
- **Test idea:** the two files; assert exit 1, the path in the message, no traceback.

### MDATK-BUG-031 — No analysed frames: IndexError in the panels, preview frame −1

- **MDBridge ref:** RN-02 · **Area:** rendering · **Severity:** crash · **Status:** Confirmed by reading
- **Where:** `mdatk/render/panels.py:604, 614-615, 666-679, 720-722` (`t[0]` on an empty time axis);
  `mdatk/cli.py:372` (`pos = min(max(0, frame), n_frames - 1)` = −1); `mdatk/render/video.py:180-203`
  (starts FFmpeg on zero frames).
- **Description:** reached with `--from-cache` of a run whose frame selection was empty, or an empty cached
  dataset.
- **Expected:** `preview` / `render` / `charts` stop with `no analysed frames to render` (exit 1); panels
  built for the GUI show `No data: no frames`.
- **Test idea:** a `Dataset` with zero frames through `build_figure` and `render_video`.

### MDATK-BUG-032 — A NaN coordinate makes the structure view crash

- **MDBridge ref:** RN-03 · **Area:** rendering · **Severity:** crash · **Status:** Confirmed by reading
- **Where:** `mdatk/render/panels.py:176-180` (radius from `norm(...).max()` is NaN), then axis limits.
- **Actual (per MDBridge's parity run on 0.0.4.1; code unchanged):** `Axis limits cannot be NaN or Inf`,
  render aborted.
- **Expected:** frame the view by finite points and warn once naming the frames; leave non-finite residues
  out of those frames.
- **Test idea:** a trace with one NaN row in one frame.

### MDATK-BUG-033 — Surface grid of a blown-up frame allocates hundreds of GiB

- **MDBridge ref:** RN-19 · **Area:** rendering (surface representation) · **Severity:** crash · **Status:** Verified
- **Where:** `mdatk/render/molecule/surface.py:22-42` (`density_grid` sizes the grid from the extent).
- **Reproduce:** `density_grid(np.array([[0, 0, 0], [5000, 5000, 5000.]]))`
- **Actual:** `MemoryError Unable to allocate 939. GiB for an array with shape (125977537197,)` (or, below
  the allocator's limit, gigabytes and minutes per frame).
- **Expected:** grids above a limit (MDBridge: 2^26 points) refused up front; the surface left out of that
  frame with a warning.
- **Test idea:** the call above returns quickly with the warning.

### MDATK-BUG-034 — Matrix chart of thousands of residues exceeds Agg's size limit

- **MDBridge ref:** RN-05 · **Area:** charts · **Severity:** crash · **Status:** Confirmed by reading
- **Where:** `mdatk/render/charts.py:123` (`height = max(cp.height, 0.18 * rows + 1.5)` inches).
- **Description:** 4000 rows at 150 dpi → ~108 000 px, beyond Agg's 2^16 px limit; the chart step fails
  after the analyses ran. Below the limit, images of hundreds of MB.
- **Expected:** height capped (MDBridge: 2^15 px; labels are thinned anyway) with a warning naming the cap.
- **Test idea:** a 5000-row matrix chart.

### MDATK-BUG-035 — Panel option `max_ticks = 0` raises ZeroDivisionError

- **MDBridge ref:** (rendering ticks, no register id; "thin_ticks") · **Area:** rendering · **Severity:** crash · **Status:** Verified
- **Where:** `mdatk/render/panels.py:60-62` (`n / max_ticks`), used with the panel option at `:744`.
- **Reproduce:** `from mdatk.render.panels import thin_ticks; thin_ticks(10, 0)`
- **Actual:** `ZeroDivisionError: division by zero`; `thin_ticks(10, -3)` labels every row.
- **Expected:** any value below 1 treated as 1.
- **Test idea:** `thin_ticks(10, 0)` and `(10, -3)` give `[0]`.

### MDATK-BUG-036 — Panel options of the wrong type abort the render

- **MDBridge ref:** RN-15, RN-40 · **Area:** rendering · **Severity:** crash · **Status:** Confirmed by reading
- **Where:** `mdatk/render/panels.py:617` (`set_ylim(*options["ylim"])` with `ylim = 5` → TypeError);
  `mdatk/render/seqpanels.py:186` (`int(options["per_row"])` with `"x"` → ValueError; a negative value
  inverts the axis and draws giant boxes).
- **Expected:** `check()`/`validate()` reports `layout.panels[k].options.ylim: expected two numbers`; a
  panel with an unusable option shows `No data: option 'per_row' must be a whole number …`.
- **Test idea:** each panel type with options of the wrong type through `check()` and `build_figure`.

### MDATK-BUG-037 — Example data of a single-strand structure is corrupt

- **MDBridge ref:** RN-21 · **Area:** GUI layout preview · **Severity:** crash · **Status:** Verified
- **Where:** `mdatk/render/example.py:230-238` (`m = 0` for one strand; `a[:m - 2]` is `a[:-2]`, not empty).
- **Reproduce:**
  ```python
  import numpy as np; from mdatk.render import example
  r = example._example_basepairs(np.zeros(10, dtype=int), 12)
  print({k: v.shape for k, v in r.arrays.items()})
  ```
- **Actual:** `map_i (48,)`, `map_v (48,)` but `map_j (0,)`: a corrupt CSR map; drawing it raises
  `IndexError`, so the GUI's layout preview breaks for a single-strand structure.
- **Expected:** no example pairs for one strand.
- **Suggested fix:** `if m == 0: return` an empty map; clamp `m - open_end` at 0.
- **Test idea:** the snippet; assert equal lengths.

### MDATK-BUG-038 — Software renderer: equal view limits raise ZeroDivisionError

- **MDBridge ref:** RN-28 · **Area:** software 3D renderer · **Severity:** crash · **Status:** Confirmed by reading
- **Where:** `mdatk/render/molecule/software.py:105-106` (`W / (xlim[1] - xlim[0])`). Reversed or NaN
  limits / rotation, and coordinates beyond float32, give garbage; a sphere larger than the image can
  raise `MemoryError`.
- **Expected:** `ValueError("3D view: …")` for unusable views; unplaceable primitives skipped; sprites
  clipped to the image.
- **Test idea:** render one sphere with `xlim = (1, 1)`.

---------------------------------------------------------------------------------------------------

## 4. File readers and the loader

All files below are tiny and written inline; `load` means:

```python
from mdatk.io.loader import load_universe
def load(top, *trajs):
    u = load_universe(top, trajs); print(u.atoms.n_atoms, u.trajectory.n_frames, u.dimensions, u.atoms.positions[0], u.atoms.masses)
```

### MDATK-BUG-039 — Hybrid-36 serials in CONECT records are rejected (inherited from MDAnalysis)

- **MDBridge ref:** RD-B01 · **Area:** PDB reader · **Severity:** crash · **Status:** Verified
- **Reproduce:**
  ```
  ATOM  A0000  O1  HOH A   1       0.000   0.000   0.000  1.00  0.00           O
  ATOM  A0001  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H
  CONECTA0000A0001
  END
  ```
  `load("b01.pdb")`
- **Actual:** `ValueError: Failed to construct topology from file b01.pdb with parser <class 'MDAnalysis.topology.PDBParser.PDBParser'>. Error: invalid literal for int() with base 10: 'A0000'`.
  ATOM serials accept hybrid-36; CONECT does not, so PDB files of > 99 999 atoms with CONECT cannot be loaded.
- **Expected:** loads, with the bond.
- **Suggested fix (MDATK side):** pre-scan the PDB; if CONECT records hold hybrid-36 serials, decode them
  and pass the bonds yourself (load without CONECT, then `u.add_bonds`), or rewrite the CONECT records to
  decimal indices in a temporary copy.
- **Test idea:** the file above loads with one bond.

### MDATK-BUG-040 — One-frame mdcrd without a box: empty StopIteration (inherited from MDAnalysis)

- **MDBridge ref:** RD-A01 · **Area:** Amber mdcrd reader · **Severity:** crash · **Status:** Verified
- **Reproduce:** `printf 'title\n%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n' 0 0 0 0.96 0 0 > one_frame.mdcrd`
  with a 2-atom PDB `two.pdb`; `load("two.pdb", "one_frame.mdcrd")`.
- **Actual:** `StopIteration` with an empty message (CLI: traceback).
- **Expected:** a non-periodic trajectory of 1 frame.
- **Suggested fix:** detect the case before handing the file to MDAnalysis (one frame and EOF) and read
  the coordinates yourself, or at least turn the error into `<file>: …` (MDATK-BUG-028).

### MDATK-BUG-041 — NetCDF frames with cell lengths (0, 0, 0): TypeError (inherited from MDAnalysis)

- **MDBridge ref:** RD-N01 · **Area:** NetCDF reader · **Severity:** crash · **Status:** Verified
- **Reproduce:** write an Amber NetCDF (scipy `netcdf_file`, version 2) with `cell_lengths = [0, 0, 0]` and
  `cell_angles = [90, 90, 90]` per frame; `load("two.pdb", "zero_box.nc")`.
- **Actual:** `TypeError: 'NoneType' object is not subscriptable`. Some tools write such zero cells for
  non-periodic runs.
- **Expected:** frames with zero lengths have no cell.
- **Suggested fix:** detect all-zero `cell_lengths` and open the file with a reader option or a
  transformation that drops the cell; or read such files without the cell variables.

### MDATK-BUG-042 — parm7 with an unknown %FLAG right after POINTERS (inherited from MDAnalysis)

- **MDBridge ref:** RD-P01 · **Area:** parm7 topology reader · **Severity:** crash · **Status:** Verified
- **Reproduce:**
  ```python
  t = open("demo/demo.parm7").read(); i = t.index("%FLAG POINTERS"); j = t.index("%FLAG", i + 5)
  open("unknown_first.parm7", "w").write(t[:j] + "%FLAG MY_EXTRA_SECTION\n%FORMAT(10I8)\n       1\n" + t[j:])
  ```
  `load("unknown_first.parm7")`
- **Actual:** `UnboundLocalError: cannot access local variable 'numlines' where it is not associated with a value`.
  The same section later in the file is fine.
- **Expected:** unknown sections skipped wherever they are (newer tleap / ParmEd write extra sections).
- **Suggested fix:** MDATK already parses parm7 sections itself (`mdatk/io/parm7.py`); move unknown
  sections after POINTERS to the end in a temporary copy, or build the topology from MDATK's own parser.

### MDATK-BUG-043 — DCD with NSAVC = 0: ZeroDivisionError (inherited from MDAnalysis)

- **MDBridge ref:** RD-D01 · **Area:** DCD reader · **Severity:** crash · **Status:** Verified
- **Reproduce:** write a 2-frame DCD with MDAnalysis' writer, then set the int32 at byte offset 16 (NSAVC)
  to 0; `load("two.pdb", "nsavc0.dcd")`.
- **Actual:** `ZeroDivisionError: division by zero`.
- **Expected:** dt from the header's DELTA with a warning, or `<file>: DCD header has frame interval 0`.
- **Suggested fix:** check the header before opening and pass `dt=` explicitly.

### MDATK-BUG-044 — Trajectory chain with mixed extensions picks the wrong reader

- **MDBridge ref:** LD-01 · **Area:** loader · **Severity:** crash · **Status:** Verified
- **Where:** `mdatk/io/loader.py:56-60` (the extension table is applied only when all trajectories share
  one extension).
- **Reproduce:** two copies of a valid 2-frame mdcrd named `two.mdcrd` and `two.crd`;
  `load("two.pdb", "two.mdcrd", "two.crd")`.
- **Actual:** `TypeError: Unable to read two.crd with <class 'MDAnalysis.coordinates.CRD.CRDReader'>`
  (each file alone loads). `.rst7` in a chain gives `Unknown coordinate trajectory format 'RST7'`.
- **Expected:** the format chosen per file.
- **Suggested fix:** pass a list of `(filename, format)` tuples to `mda.Universe` (MDAnalysis' chain reader
  accepts per-file formats), from `_TRAJECTORY_FORMATS`.
- **Test idea:** the mixed chain loads with 4 frames.

### MDATK-BUG-045 — Wide GRO coordinates are read from single characters, silently (inherited from MDAnalysis)

- **MDBridge ref:** RD-G01 · **Area:** GRO reader · **Severity:** wrong results · **Status:** Verified
- **Reproduce:**
  ```
  wide
      1
      1SOL     OW    1 10000.000   2.000   3.000
     3.00000   3.00000   3.00000
  ```
  `load("wide.gro")`
- **Actual:** position `[10.0, 0.0, 0.0]` Å (expected `[100000, 20, 30]`); no error. The field width is
  taken from the first `.` at or after column 26 of the first atom line.
- **Expected:** width = distance between the decimal points of the y and z fields, checked against x;
  otherwise `<file> line 3: cannot determine coordinate field width`.
- **Suggested fix (MDATK side):** pre-check the first atom line of GRO inputs and refuse inconsistent
  widths with that message.

### MDATK-BUG-046 — Multi-frame GRO: only the first frame is read, silently (inherited from MDAnalysis)

- **MDBridge ref:** RD-G07 · **Area:** GRO reader · **Severity:** wrong results · **Status:** Verified
- **Reproduce:** a GRO file holding two frames (two title/count/atoms/box blocks) as topology alone:
  `load("two_frames.gro")`.
- **Actual:** `n_frames 1`; the second frame is ignored without a message.
- **Expected:** every frame read (GRO is a valid GROMACS trajectory), or at least
  `<file>: only the first of N frames is read`.
- **Suggested fix:** count blocks; if > 1, warn, or read the frames yourself into a `MemoryReader`.

### MDATK-BUG-047 — One-atom mdcrd with a box: box lines become extra frames (inherited from MDAnalysis)

- **MDBridge ref:** RD-A03 · **Area:** mdcrd reader · **Severity:** wrong results · **Status:** Verified
- **Reproduce:** a one-atom PDB and an mdcrd of 3 frames, each `x y z` line followed by a box line
  `30.000 30.000 30.000`; `load("one_atom.pdb", "n1.mdcrd")`.
- **Actual:** `n_frames 6`, no box: every box line is a frame.
- **Expected:** 3 frames with a 30 Å cell.
- **Suggested fix:** for N = 1, detect a box line (exactly 3 fields followed by the next frame or EOF) and
  open with `format="MDCRD"` and the box option, or read the file yourself.

### MDATK-BUG-048 — Partial PDB element column: atoms get mass 0 silently (inherited from MDAnalysis)

- **MDBridge ref:** RD-B10 · **Area:** PDB reader · **Severity:** wrong results · **Status:** Verified
- **Reproduce:** two ATOM lines, the first with element `C`, the second with a blank element column;
  `load("b10.pdb")`.
- **Actual:** masses `[12.011, 0.0]`; the reader's warning is discarded (MDATK-BUG-049).
  Mass-weighted RMSD and heavy-atom logic then use mass 0.
- **Expected:** element and mass of blank-element atoms guessed from the name, as for files without the
  column; one warning with the count of atoms whose mass stays unknown.
- **Suggested fix:** after loading, guess masses for atoms whose mass is 0 and whose element is blank.

### MDATK-BUG-049 — All reader warnings are discarded

- **MDBridge ref:** LD-03 (also RD-B13) · **Area:** loader · **Severity:** misleading message · **Status:** Verified
- **Where:** `mdatk/io/loader.py:62-63` (`warnings.simplefilter("ignore")` around `mda.Universe`);
  `mdatk/cli.py:1632`.
- **Reproduce:** `mdatk info --top b10.pdb --traj two.mdcrd` prints no warning, while loading the same
  file with warnings on shows `Unknown element  found for some atoms` and `Unknown masses …`. The same
  `info` reports `"dt_ps": 1.0` for an mdcrd that has no time information.
- **Expected:** reader warnings collected and printed once (`WARNING: <file>: …`; GUI: the log): placeholder
  CRYST1, missing elements, unknown masses, skipped sections, missing times (with `dt_ps: null` or a note).
- **Suggested fix:** `warnings.catch_warnings(record=True)` and forward each unique message through the
  reporter.

### MDATK-BUG-050 — A GROMACS `.top` is read as parm7

- **MDBridge ref:** LD-02 · **Area:** loader · **Severity:** misleading message · **Status:** Verified
- **Where:** `mdatk/io/loader.py:18` (`.top` → `PARM7`).
- **Reproduce:** `printf '[ defaults ]\n1 2\n' > gmx.top; mdatk info --top gmx.top --traj two.mdcrd`
- **Actual:** `error: Failed to construct topology from file …/gmx.top with parser <class 'MDAnalysis.topology.TOPParser.TOPParser'>. Error: …/gmx.top is not a valid TOP file. %VE Missing in header`.
- **Expected:** keep the mapping, but when line 1 is not `%VERSION`: `<file>: not an Amber parm7 topology
  (GROMACS .top files are not supported)`.

### MDATK-BUG-051 — Reader errors are Python internals naming MDAnalysis classes (inherited from MDAnalysis)

- **MDBridge ref:** RD-B05, RD-B06, RD-B07, RD-G02, RD-P11 (more in the appendix) · **Area:** readers · **Severity:** misleading message · **Status:** Verified
- **Reproduce and actual:**
  - empty `empty.pdb`: `EOFError: Compressed file ended before the end-of-stream marker was reached`
    (the file is not compressed; same for GRO and mdcrd);
  - PDB with only `REMARK` and `END`: `IndexError: index 0 is out of bounds for axis 0 with size 0`;
  - PDB with a Latin-1 byte in a REMARK: `Failed to construct topology … with parser <class
    'MDAnalysis.topology.PDBParser.PDBParser'>. Error: 'utf-8' codec can't decode byte 0xe9 in position 14`;
  - GRO without the box line: `UnboundLocalError: cannot access local variable 'unitcell' …`.
- **Expected:** `<file>: file is empty`, `<file>: no ATOM or HETATM records`, `<file>: missing box line
  after atom N`; non-UTF-8 bytes read as Latin-1 with one warning naming the line; the format named in
  words, not by MDAnalysis class.
- **Suggested fix:** pre-checks for the cheap cases (empty file, no atom records, GRO box line) before
  calling MDAnalysis, and one wrapper that rewrites `Failed to construct topology … with parser <class …>`
  to `<file>: cannot read <format> topology: <message>`. Latin-1: decode and re-encode into a temporary
  UTF-8 copy.

---------------------------------------------------------------------------------------------------

## 5. Selections

All through `mdatk.selections.select(u, text, allow_empty=True)` on the demo universe
(`u = mda.Universe("demo/demo.parm7", "demo/demo.nc")`; masses of atoms 0–3 are 14.01, 12.01, 14.01, 12.01).

### MDATK-BUG-052 — `prop 6 < x` means `x >= 6` (inherited from MDAnalysis)

- **MDBridge ref:** SEL-01 · **Severity:** wrong results · **Status:** Verified
- **Actual:** `prop mass > 14.01` → 4263 atoms, but `prop 14.01 < mass` → 4354 (= `mass >= 14.01`);
  `prop 14.01 <= mass` → 4263 (= `mass > 14.01`). Reversed operands get the *negated* operator.
- **Expected:** `V < P` means `P > V` (mirror, not negate).
- **Suggested fix (MDATK side):** before handing the text to MDAnalysis, rewrite `prop [abs] <number> <op>
  <name>` to `prop [abs] <name> <mirrored op> <number>` (a regex in `resolve_selection`).
- **Test idea:** the four forms select identical sets pairwise.

### MDATK-BUG-053 — Malformed `prop` operators are accepted (inherited from MDAnalysis)

- **MDBridge ref:** SEL-02 · **Severity:** wrong results · **Status:** Verified
- **Actual:** `prop mass ><14` → 8587 atoms (text before the operator is dropped); `prop mass <> 14` →
  `could not convert string to float: '>'`.
- **Expected:** `invalid operator '><' in prop (use <, >, <=, >=, ==, !=)`.
- **Suggested fix:** validate the operator token in a pre-pass and raise `SelectionError`.

### MDATK-BUG-054 — `atom` swallows keywords (inherited from MDAnalysis)

- **MDBridge ref:** SEL-05 · **Severity:** wrong results · **Status:** Verified
- **Actual:** `atom A 1 and` → 0 atoms, no error (`and` read as the atom name); likewise `atom and 1 P`.
- **Expected:** `atom: expected segid resid name`.
- **Suggested fix:** pre-pass: an `atom` argument that is a selection keyword or parenthesis is an error
  (a backslash escape still allows the literal).

### MDATK-BUG-055 — Negative `around` radius: error or empty result depending on atom counts (inherited from MDAnalysis)

- **MDBridge ref:** SEL-06 · **Severity:** misleading message · **Status:** Verified
- **Actual:** `around -2 resid 1` → `invalid selection …: MDAnalysis says: Cutoff must be positive`, but
  `around -2 (resid 1 and name P)` → 0 atoms, no error. Which happens depends on the internal search method
  (atom counts); `Cutoff … too large for box` errors are size-dependent in the same way.
- **Expected:** any negative radius: `around: distance must not be negative`; any non-negative radius
  answered without size-dependent errors.
- **Suggested fix:** check the radius in a pre-pass; for the box errors see MDATK-BUG-017.

### MDATK-BUG-056 — Python texts at end of input and in float ranges (inherited from MDAnalysis)

- **MDBridge ref:** SEL-03, SEL-09 (SEL-04 is partly addressed by 0.0.10's hint) · **Severity:** misleading message · **Status:** Verified
- **Actual:** `around` / `point 1 2 3` → `float() argument must be a string or a real number, not 'NoneType'`;
  `mass 12:+13`, `mass -inf:inf`, `charge :1` → `could not convert string to float: ''`;
  `resid 1 AND name P` → `Failed to parse value: AND`.
- **Expected:** one message per keyword with its signature (`around: expected a distance and a selection,
  e.g. 'around 3.5 resname MG'`, `point: expected x y z radius`); range bounds accept the same float
  syntax as single values (sign, `E`, `inf`); upper-case keywords get `keywords are lower case (did you
  mean 'and'?)`.
- **Suggested fix:** extend `_syntax_hint` with these patterns (the `'None'` hint of 0.0.10 is the model);
  pre-normalise float range bounds.

---------------------------------------------------------------------------------------------------

## 6. Rendering and video

### MDATK-BUG-057 — Heatmap columns are all `t[1] - t[0]` wide

- **MDBridge ref:** RN-01 · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/render/panels.py:720-722`; `mdatk/render/charts.py:131`.
- **Description:** with non-uniform times (chained runs with different time steps, time resets) the cells
  sit at the wrong times and the moving cursor leaves its frame's cell.
- **Expected:** column edges at the midpoints between frame times, half a step outward at both ends
  (identical to today for uniform times).
- **Test idea:** 12 frames 10 ps apart then 12 frames 2 ps apart; check cell edges.

### MDATK-BUG-058 — Running mean: one NaN blanks the whole window

- **MDBridge ref:** RN-14 · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/render/panels.py:569-577` (`np.convolve` of padded data).
- **Reproduce:** `_running_mean(np.array([1,2,3,np.nan,5,6,7,8,9,10.]), 5)`
- **Actual:** `[1.6 nan nan nan nan nan 7. 8. 8.8 9.4]` (an empty RMSD target or NaN frames erase whole
  stretches of a smoothed curve).
- **Expected:** NaN-aware mean (sum of finite values / count of finite values); NaN only where the whole
  window is NaN.

### MDATK-BUG-059 — Class matrices with more than 5 classes reuse colours

- **MDBridge ref:** RN-18 · **Severity:** wrong results · **Status:** Verified
- **Where:** `mdatk/render/theme.py:51-63` (`ordinal`: for n − 1 > 4 the 4 ramp colours are kept).
- **Reproduce:** `LIGHT.ordinal(8)` → 5 colours for 8 classes, while the colour bar labels 8 classes.
- **Expected:** n − 1 evenly spaced ramp steps (`LinearSegmentedColormap.from_list(ramp, N=n-1)`),
  today's colours for n ≤ 5.

### MDATK-BUG-060 — One NaN voxel removes the density shadow

- **MDBridge ref:** RN-30, RN-33 · **Severity:** wrong results · **Status:** Confirmed by reading
- **Where:** `mdatk/io/mrc.py:44-47` (`stats` with `mean` / `std` → NaN), `mdatk/density.py:257-258`
  (`map_level`), `:704` (`np.maximum(map - level, 0)` keeps NaN), `:688-691` (peak NaN). An infinite voxel
  gives NaN or a blank shadow too.
- **Expected:** non-finite voxels carry no density; the level from the finite voxels.
- **Suggested fix:** `np.nanmean` / `np.nanstd` on finite data; `np.where(np.isfinite(d), ...)`.

### MDATK-BUG-061 — The ion legend always says `Mg²⁺`

- **MDBridge ref:** RN-08 · **Severity:** misleading message · **Status:** Confirmed by reading
- **Where:** `mdatk/render/panels.py:468`.
- **Description:** the ions drawn come from `analyses.mg.mg_selection`; with Na⁺, K⁺ or Ca²⁺ the legend
  still says Mg²⁺.
- **Expected:** a label from the ions' element(s) (`K⁺`, `Na⁺`, …; `Ions` for mixtures).

### MDATK-BUG-062 — The chain legend numbers chains although letters are known

- **MDBridge ref:** RN-09 · **Severity:** misleading message · **Status:** Confirmed by reading
- **Where:** `mdatk/render/panels.py:455` (`f"Chain {int(c) + 1}"`).
- **Description:** next to residue tags like `B:12` and a sequence panel showing `A, B`, the legend says
  `Chain 1, 2`.
- **Expected:** `Chain A` (the trace's chain names; the index only for chains without a name).

### MDATK-BUG-063 — A matrix key in `output.charts.series` is dropped silently

- **MDBridge ref:** RN-10 · **Severity:** misleading message · **Status:** Confirmed by reading
- **Where:** `mdatk/render/charts.py:31-37` (keys of another kind are filtered out without a warning;
  only missing keys warn).
- **Expected:** a requested matrix is charted as a matrix chart; a pair-map key warns
  `chart series K is a pair map; …`.

### MDATK-BUG-064 — Bad chart formats and colours fail after the analyses, without naming the setting

- **MDBridge ref:** RN-11, RN-12 · **Severity:** misleading message · **Status:** Verified
- **Where:** `mdatk/render/charts.py:108-110, 154-156`; `validate()` (`mdatk/config.py:694-`) checks no
  format and no colour.
- **Reproduce:**
  ```
  mdatk charts -p demo/demo.mdatk.toml -o o6 --set 'output.charts.formats=["png","pnng"]'
  mdatk preview -p demo/demo.mdatk.toml -o o0 --set 'visual.glow_color="notacolour"'
  ```
- **Actual:** the analyses run, one png is written, then `error: Format 'pnng' is not supported (supported
  formats: avif, eps, gif, …)`; the preview fails with `error: Invalid RGBA argument: 'notacolour'` (no
  setting named; in a video the error comes mid-render).
- **Expected:** both reported by `validate()` before anything runs, naming the field
  (`output.charts.formats: …`, `visual.glow_color: 'notacolour' is not a colour`).
- **Suggested fix:** check formats against `FigureCanvasBase.get_supported_filetypes()` and every colour
  field with `matplotlib.colors.is_color_like` in `validate()` (as warnings when the step does not run,
  per 0.0.10's rule).

---------------------------------------------------------------------------------------------------

## 7. Command line and messages

### MDATK-BUG-065 — KeyErrors printed with Python's quotes

- **MDBridge ref:** CLI-01 · **Severity:** misleading message · **Status:** Verified
- **Where:** `mdatk/cli.py:1641-1645` (`print(f"error: {exc}")` of a KeyError prints its repr).
- **Reproduce:** `mdatk analyze -p demo/demo.mdatk.toml -o o3 --only nope -q`
- **Actual:** `error: "unknown analysis 'nope' (known: contacts, energy, …)"` (also in
  `WARNING: skipping data source 'bogus@Nope': "unknown analysis …"`).
- **Expected:** no surrounding double quotes (use `exc.args[0]` for KeyError).

### MDATK-BUG-066 — `--frames x`, `--set groups.--1.name=x`: `invalid literal for int()`

- **MDBridge ref:** PRJ-04, PRJ-05 · **Severity:** misleading message · **Status:** Verified
- **Where:** `mdatk/cli.py:121-123` (`--frames`), `:69-71` (`'--1'.lstrip('-').isdigit()` passes, then `int('--1')`).
- **Reproduce:** `mdatk init x.toml --frames x`; `mdatk init x.toml -p demo/demo.mdatk.toml --set 'groups.--1.name=x'`
- **Actual:** `error: invalid literal for int() with base 10: 'x'` / `… '--1'`.
- **Expected:** `--frames: 'x' is not a frame number or range`; `groups.--1.name: no list element '--1'`.

### MDATK-BUG-067 — RMSD reference frame beyond the run is clamped silently

- **MDBridge ref:** PRJ-02 · **Severity:** misleading message · **Status:** Verified
- **Where:** `mdatk/analysis/rmsd.py:71`.
- **Reproduce:** `mdatk analyze -p demo/demo.mdatk.toml -o o4 --only rmsd --set analyses.rmsd.reference_frame=999 -q`
- **Actual:** runs with the last frame as reference, no message.
- **Expected:** keep the clamp, warn `RMSD: reference frame 999 is beyond the 6 analysed frames; using the
  last` (`… -3 is negative; using the first`).

---------------------------------------------------------------------------------------------------

## 8. Minor

### MDATK-BUG-068 — Quitting during a job kills the worker mid-write

- **MDBridge ref:** AP-14 · **Area:** GUI · **Status:** Confirmed by reading
- **Where:** `mdatk/gui/app.py:2551-2553` (`jobs.stop()` with the default `wait=0`), `mdatk/gui/jobs.py:95-101`,
  `mdatk/gui/worker.py:36` (`daemon=True`).
- **Description:** the analysis window cancels the job but does not wait; the daemon thread dies with the
  process, possibly mid-way through a cache archive or a temporary video file. (RNAdo asks since 0.0.11;
  this window does not.)
- **Suggested fix:** ask; on confirmation cancel and join the worker (it already cleans up its temporary
  files on cancellation) before quitting.

### MDATK-BUG-069 — `settings.toml` written non-atomically; an unreadable file is replaced

- **MDBridge ref:** AP-17 (remaining part; 0.0.14 fixed the loss of other keys) · **Area:** user settings · **Status:** Confirmed by reading
- **Where:** `mdatk/projectfile.py:623-629` (`open(target, "w")`), `:611-617` (an unparseable file reads as
  `{}`), `:643-650` (`add_recent` then writes only what it read).
- **Description:** a crash or a concurrent write (MDATK GUI and RNAdo both update the file) can leave a
  truncated file, and the next `add_recent` replaces it with `recent` alone, losing `cat_proof` and any
  other setting.
- **Suggested fix:** write to a temporary file in the same folder and `os.replace` (as the project writer
  already does at `projectfile.py:358-362`); do not write back over a file that failed to parse.

### MDATK-BUG-070 — Energy: net-force direction inconsistent with the 0.5 Å clamp

- **MDBridge ref:** AN-05 · **Area:** energy · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/energy.py:190, 211-213` (`d` is clamped to 0.5, `rvec` keeps the real length r,
  so the force is scaled by r / 0.5 for r < 0.5 Å).
- **Suggested fix:** unit vector `rvec / |rvec|` times the scalar force at the clamped distance.

### MDATK-BUG-071 — `atoms` analysis bond guess misses S–S

- **MDBridge ref:** AN-14 · **Area:** atoms analysis (3D view) · **Status:** Confirmed by reading
- **Where:** `mdatk/analysis/structure.py:376-403` (fixed 1.95 Å heavy-atom limit: disulfides at 2.05 Å are
  missed).
- **Suggested fix:** covalent radii (sum + 0.4 Å).
- **Withdrawn (2026-11-27):** the second half of the original item ("links chains end to end": the
  residue test `abs(resi[i] - resi[j]) <= 1` ignores the chain) was a false report. A bond also needs the
  two atoms within 1.95 Å, so chain ends that are not already in bonding distance are never joined, and
  chains covalently joined across a chain break are drawn joined. MDBridge's own one-chain rule (AN-14)
  misses those links, and disulfides between residues that are not consecutive; MDBridge adopts 0.1.2's
  rule.

### MDATK-BUG-072 — GRO box fallback loses minus signs (inherited from MDAnalysis)

- **MDBridge ref:** RD-G05 · **Area:** GRO reader · **Status:** Verified
- **Reproduce:** a one-atom GRO whose box line is `-1.00000-2.00000-3.00000`; `load("neg_box.gro")`.
- **Actual:** box `[10, 20, 30, 90, 90, 90]` Å (signs dropped; a malformed box becomes a valid-looking one).
- **Expected:** fixed-width 10-character fields with signs kept, or `GRO box line has neither 3 nor 9 values`.

### MDATK-BUG-073 — Duplicate group names: `group:` uses the last, `Project.group()` the first

- **MDBridge ref:** SEL-15 · **Area:** selections · **Status:** Verified
- **Where:** `mdatk/selections.py:160` (`by_name` dict, last wins) vs `mdatk/config.py:583-588` (first wins).
- **Reproduce:** groups `[ResidueGroup("Dup", "1-2"), ResidueGroup("Dup", "3-12")]`:
  `Project(groups=…).group("Dup").residues` → `1-2`; `select(u, "group:Dup", groups)` → residues 3–12.
- **Expected:** one rule everywhere (first wins, or refuse duplicates on load). `validate()` flags the
  duplicate, but `info` and selections still resolve it.

### MDATK-BUG-074 — `group:` token stops at invalid characters with a misleading error

- **MDBridge ref:** SEL-14 · **Area:** selections · **Status:** Verified
- **Where:** `mdatk/selections.py:20` (`_GROUP_TOKEN` stops at a space; it also allows `+` and `'`, which
  group names may not contain).
- **Actual:** `group:Bad name` → `unknown residue group 'Bad' (defined: Aptamer, Target)`;
  `group:Ünï` → `Unknown selection token: 'group:Ünï' (did you mean 'group'?)`.
- **Expected:** `group name 'Bad name' may only contain letters, digits, '_', '.' and '-'`.

### MDATK-BUG-075 — `wrap_degrees` folds angles near −180 to +180

- **MDBridge ref:** RN-22 · **Area:** rendering (orientation) · **Status:** Verified
- **Where:** `mdatk/render/orientation.py:39-41` (`np.isclose(a, -180.0)` has a relative tolerance).
- **Reproduce:** `wrap_degrees(-179.999)` → `180.0`.
- **Expected:** only −180 itself (within ~1e-9) becomes +180.

### MDATK-BUG-076 — Secondary-structure panel: an unknown `source` silently means `dynamic`

- **MDBridge ref:** RN-41 · **Area:** rendering · **Status:** Confirmed by reading
- **Where:** `mdatk/render/seqpanels.py:264, 277-284`.
- **Expected:** `No data: unknown source 'x' (option 'source': dynamic, first or dotbracket)`.

### MDATK-BUG-077 — Codec fallback ignores `crf`

- **MDBridge ref:** RN-07 · **Area:** video · **Status:** Confirmed by reading
- **Where:** `mdatk/render/video.py:57-64, 80-83`: without `libx264` the encoder falls back to `libopenh264`
  (no rate option: FFmpeg's low default bitrate) or `mpeg4` (fixed `-q:v 3`); `crf` is ignored and the
  only notice is a `log.warning`.
- **Suggested fix:** map `crf` onto each fallback encoder's quality scale (e.g. `-b:v` for openh264,
  `-q:v` from crf for mpeg4) and print the fallback as a `WARNING:` line.

### MDATK-BUG-078 — Videos are encoded with BT.601 and left untagged

- **MDBridge ref:** RVE-1 (DIFFERENCES "Colours tagged BT.709") · **Area:** video · **Status:** Confirmed by reading
- **Where:** `mdatk/render/video.py:77-84` (RGBA → `yuv420p` with FFmpeg's default BT.601 matrix, no colour
  tags). Players read untagged HD video as BT.709, so greens and reds look shifted.
- **Suggested fix:** add `-vf scale=out_color_matrix=bt709:out_range=tv` and
  `-colorspace bt709 -color_primaries bt709 -color_trc bt709 -color_range tv`.

### MDATK-BUG-079 — Example data expands whole residue ranges before keeping 400

- **MDBridge ref:** RN-25 · **Area:** GUI layout preview · **Status:** Verified
- **Where:** `mdatk/render/example.py:29-38`.
- **Reproduce:** a group `Big = "1-300000000"`; `example.group_resids(p, "Big")` returns 400 ids after
  8.3 s (gigabytes for larger ranges). The schematic sequence has 120 letters, so larger groups can hit
  an IndexError.
- **Suggested fix:** stop expanding at 400; repeat the sequence.

### MDATK-BUG-080 — Software renderer `fast` quality at odd sizes

- **MDBridge ref:** RN-27 · **Area:** software 3D renderer · **Status:** Confirmed by reading
- **Where:** `mdatk/render/molecule/software.py:104` (`round(width * 0.5)` rounds half to even),
  `:425-431` (upsampled by 2 and copied into the full size).
- **Description:** at odd widths/heights the last column/row stays transparent and the view is scaled by
  up to 1 px (GUI drags).
- **Suggested fix:** `ceil(size / 2)` with the view limits widened accordingly, then crop.

### MDATK-BUG-081 — Glow halo: seam at the molecule box; strength depends on earlier renders

- **MDBridge ref:** RN-26, RN-29 · **Area:** software 3D renderer · **Status:** Confirmed by reading
- **Where:** `mdatk/render/molecule/software.py:667-671` (`uniform_filter(..., mode="constant")` on the
  molecule's bounding-box crop treats pixels outside as 0: a darker rectangle in the halo);
  `:569-577` (`_disc_norm` cached under sizes rounded to 0.01 px, computed from whichever size first
  rounded there).
- **Suggested fix:** filter with the real neighbours (crop one pixel larger, or filter the full alpha);
  compute the normalisation for the exact size (or key the cache by the exact floats).

### MDATK-BUG-082 — Chart title and colour bar stay black on the dark theme

- **MDBridge ref:** RN-50 · **Area:** charts · **Status:** Confirmed by reading
- **Where:** `mdatk/render/charts.py:104, 152` (left title set after `style_axes`, which styles the centre
  title), `:138-144` (colour bar not styled).
- **Suggested fix:** title in the theme's ink; colour-bar outline, ticks and labels in the theme's colours.

### MDATK-BUG-083 — `stop=None` in the empty frame-selection message

- **MDBridge ref:** CLI-03 · **Area:** CLI · **Status:** Verified
- **Where:** `mdatk/io/loader.py:84-86`.
- **Reproduce:** `mdatk analyze -p demo/demo.mdatk.toml -o o3 --frames 30: -q`
- **Actual:** `error: frame selection start=30 stop=None step=1 selects no frames (trajectory has 6)`.
- **Expected:** `stop=end`.

### MDATK-BUG-084 — External data files are read only after all analyses ran

- **MDBridge ref:** CLI-02 · **Area:** pipeline · **Status:** Confirmed by reading
- **Where:** `mdatk/pipeline.py:241-245` (after `run_analyses` at `:225`); `validate()` does not check
  `inputs.external` paths.
- **Description:** a missing or malformed external file fails the command after the whole frame pass.
- **Suggested fix:** load and validate external data before the frame loop (same messages).

### MDATK-BUG-085 — PDB export: a frame beyond the PDB limits aborts the export after writing earlier files

- **MDBridge ref:** CLI-08 · **Area:** PDB export · **Status:** Confirmed by reading
- **Where:** `mdatk/io/pdb_export.py:68-78` (one file per frame; the writer raises at the first coordinate
  outside −999.999…9999.999 Å, leaving the earlier files and possibly a partial one).
- **Suggested fix:** check every requested frame's coordinate range first, then write; message
  `… frame 2 is outside this range; no files were written.`

---------------------------------------------------------------------------------------------------

## MDATK's response (2026-11-27)

MDATK 0.1.2 (commit `2d93338`, 2026-11-27) fixed the other items. It kept 023, 026 and 061 by design and
found the second half of 071 not to be a bug; MDBridge keeps its own behaviour for the three design
choices (its DIFFERENCES entries stand). 099 was first answered "by design"; the owner decided otherwise
(D11) and 0.1.2 fixed it.

| Item | MDATK's verdict | MDATK's reason | MDBridge's note |
|---|---|---|---|
| 023 | By design | Base pairs ignore periodic images; the base-pair analysis code documents this (`basepairs.py`: paired bases belong to one molecule, which cleaned trajectories keep whole). | Holds when the cleaning keeps a duplex's two strands together, as cpptraj's `autoimage` normally does around its anchor; it does not hold when each strand is imaged on its own (for example `gmx trjconv -pbc mol`). Our "which `autoimage` does" overstated how often that happens. MDBridge keeps minimum-image distances (AN-11). |
| 026 | By design | `use_parm7` means "use the parm7 when one is available", as the GUI label and README say. | Accepted. MDBridge keeps its warning when the topology is not a parm7 (AN-06); the option's meaning is unchanged. |
| 061 | By design | The ion selection is documented as the Mg²⁺ ions, so the legend says Mg²⁺; labelling from the selected ions would need the trace analysis to store their elements. | Accepted for MDATK. MDBridge stores the ions' labels and names them (RN-08). |
| 099 (addendum) | Fixed in 0.1.2 | First answered by design (the README documented the dropped hydrogens); changed after the owner's decision D11. | When the starting model has hydrogens, the mutated result has them too, in MDBridge (CH-09) and in 0.1.2. |
| 071, second half | False report | Two chain ends bond only when their atoms are within bonding distance, so separate chains are not joined. | Agreed; withdrawn in the item above. |

## Appendix A: further reader-message items (not itemised, not re-run against 0.1.1)

MDATK 0.1.1 still reads every format through MDAnalysis with the same loader (`mdatk/io/loader.py`
changed only for the time-step override), so these are expected to remain; all are Python-text or
misleading errors (severity: misleading message), MDBridge refs in brackets:

- parm7: file ends early → `StopIteration` [RD-P02]; POINTERS with too few values → `list index out of
  range` [RD-P03]; missing required section → `KeyError 'bonda'` [RD-P04]; a blank field mid-line shifts
  later values, surfacing as `Length of charges does not match …` [RD-P05]; Fortran `D` exponents
  rejected [RD-P07]; RESIDUE_POINTER not validated (silently wrong residues) [RD-P08]; ATOMIC_NUMBER ≥ 119
  → `KeyError 119` [RD-P09]; bond index out of range / not a multiple of 3 [RD-P10].
- PDB: bad hybrid-36 serial → `KeyError 'b'` [RD-B02]; models without CRYST1 take the cell of the
  previously *read* frame, so random access gives order-dependent cells [RD-B03]; indented ATOM lines
  counted in the topology but not in frames [RD-B04]; bad coordinate field without a line number [RD-B08].
- GRO: header count larger than the atom lines; non-integer count [RD-G03, RD-G04]; nine-zero triclinic
  box silently "no cell" [RD-G06].
- NetCDF: atom-count mismatch message gives MDAnalysis API advice [RD-N02]; missing `units` →
  `AttributeError` [RD-N03]; CDF-5 / truncated files → numpy texts [RD-N04]; streaming record count
  rejected [RD-N05]; dt ignores the `time` scale factor [RD-N06]; 1-D `cell_*` variables misread [RD-N07];
  CDF-1 files rejected [RD-N08]; "no dt" warning on every access [RD-N10].
- mdcrd: atom count never checked (misframed data or `cannot reshape`) [RD-A02]; trailing blank lines
  count as a frame [RD-A04]; `********` fields without explanation [RD-A05].
- DCD: 4-D DCD → numpy broadcast error [RD-D02]; XDR errors without file / frame [RD-D03].

## Appendix B: items checked and left out

- **Fixed in MDATK after 0.0.4.1 (dropped):** AP-09 (swap with invalid options JSON, 0.1.0), AP-11 (*Save
  as* suffix after the overwrite question, 0.1.0), AP-12 (startup load errors in the window, 0.1.0),
  AP-15 (table cells commit on focus loss, 0.1.0), AP-18 (selectable dialog texts, 0.1.0), AP-17's loss
  of other keys (0.0.14; the non-atomic write remains as 069), AN-04 (triclinic minimum-image vectors,
  0.1.0), AN-10 (N with a guessed H is no acceptor, 0.1.0), AN-21 / RN-23 (unknown structure groups and
  unreadable `label_residues` now warned at check time, 0.0.10), AN-23 (contacts of a group with
  itself, 0.1.0), RN-24 (chain-qualified tags, 0.1.0), RD-P06 (parm7 LJ-table format lines, 0.1.0),
  LD-05 (time-step override, 0.1.0: the override now replaces stored times by design), SEL-04 (the
  `'None'` token now has a hint, 0.0.10), AN-25 (generic LJ table replaced, 0.1.1).
- **Not reproduced (dropped):** RD-B09 (one corrupt CONECT record dropping all bonds: two malformed
  records tried, bonds kept); RN-06 (canvas size differing from the request: 723 size/DPI combinations
  tried, none differs); RN-04 (reversed heatmap limits: MDATK does not raise, as MDBridge's own notes
  record); RN-13 (unknown `chain_palette`: the demo preview did not reach it).
- **Not bugs in MDATK (excluded):** MDBridge-only design choices and conventions (AP-16, AP-19, RN-16,
  RN-17, RN-20, Barnaba / VTK / x264 replacements, AV1, citations CIT-*, data tables DATA-*, KEEP items
  SEL-11/12/16/18, RD-B12, RD-G08, RD-N09, RD-A07, RD-D04, AN-22/24, CLI-06/09, PRJ-06/08, MRC-3/4);
  AP-20 (canvas drawn on the GUI thread) is performance, not correctness; AN-26's other refinements
  (CT / OH / OS / OW types) are accuracy improvements over a documented table, not bugs.
