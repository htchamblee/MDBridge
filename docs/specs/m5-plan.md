# Milestone 5 plan: the `app` module

Goal (ARCHITECTURE "Migration order" 5 + "Definition of done"): MDATK v0.0.4.1's application layer
in C++ and Qt 6 Widgets: the analysis window (`mdatk/gui/app.py`, `widgets.py`, `canvas.py`,
`worker.py`, `desktop.py`; 3 130 lines) with its six settings pages, the WYSIWYG canvas layout
editor, the run / progress / log views, the About and citations dialogs, and the command-line
commands of MDATK's analysis side that MDBridge does not have yet. The planner, script wizard and
runner (`planner.py`, `wizard.py`; milestone 7) and the PDB clean-up and structure builder windows
(`pdbdialog.py`, `builderdialog.py`; milestone 6) are inventoried here only to fix the boundary.

Licensing: MDATK code (owner's) is ported freely. The app layer involves no MDAnalysis behaviour, so
it has no clean-room specs. Qt Widgets, Qt Test and the platform plugins are Qt (LGPL v3, already
bundled); no new third-party dependency is planned (section 7, D-items name the only candidates).

References: MDATK v0.0.4.1 = `htchamblee/MDATK` tag `v0.0.4.1` (9875bb6); "0.1.1" = its later
`main` (1efa0ef, version 0.1.1, tags up to v0.1.0), read for section 5 only. Line numbers are
v0.0.4.1's unless a version is named.

---------------------------------------------------------------------------------------------------

## 1. Inventory of MDATK's app layer

### 1.1 Main window (`gui/app.py`)

| Element | MDATK (file:line) | Milestone |
|---|---|---|
| Window 1320 x 920, title `MDATK`, subtitle = project path or `unsaved project` | app.py:83-133, 1431 | 5 |
| Header bar: *Open project* + drop-down (recent projects, *Demo project…*, *Browse…*), *Save*, *Save as…*, *Workflow* menu (Project Planner…, Script wizard…, PDB clean-up…, Structure builder…), main menu (New project (defaults), Create demo system…, Project log…, Export citations (BibTeX)…, About and citations) | app.py:135-224 | 5; Workflow items: shells (D2) |
| Toolbar: six page toggle buttons (icons, accessible names), hint text, *Load real frame*, *Refresh* | app.py:226-254 | 5 (*Load real frame* is automatic in MDBridge since M4) |
| Side bar: fixed 520 px, page title, stacked pages, `open_section(page, frame)` scrolls to a section | app.py:53, 256-321 | 5 |
| Canvas (right): WYSIWYG figure, section 1.3 | app.py:323-367 | 5 |
| Log expander `Log (pipeline and cpptraj output)` (monospace, Clear, keeps 5 000-6 000 lines, follows the end), logging handler for WARNING and above | app.py:72-81, 369-382, 2003-2010 | 5 |
| Action bar: *Render* (suggested action), *Cancel*, progress bar (`Idle`, `<title> …`, `NN%  message`, `Failed`, `Cancelled`, `<title>: done`), *Open output folder*, *Open video* | app.py:384-405, 1911-1995 | 5 |
| Shortcuts: Ctrl+O open, Ctrl+S save, Ctrl+Shift+S save as, Ctrl+R render, Ctrl+P planner, Ctrl+Shift+W wizard (no others; no undo / redo anywhere in MDATK) | app.py:407-420 | 5 (Ctrl+P / Ctrl+Shift+W reserved for M7) |
| Start: `mdatk gui [PROJECT]`; a project that fails to load prints `warning: could not load …` on stderr and opens defaults; a params log opens as an unsaved project | app.py:2046-2075, cli.py:182 | 5 |

### 1.2 Settings pages (one per part of `.mdatk.toml`; two-way `Binder`, widgets.py:314-406)

| Page (id) | Sections and fields (project key) | MDATK | Milestone |
|---|---|---|---|
| Input & cleaning (`input`) | Project notes (`notes`, text); Input files (`inputs.topology` file chooser with filter, `inputs.trajectories` list editor + Browse multi-select); Frames, time and external data (`start`/`stop`/`step`, `dt_ps`, `time_unit`, `inputs.external` table: path, label, x col, x unit combo, columns); Trajectory cleaning (`cleaning.*`: enabled, cpptraj path, conda env + *Detect*, conda exe, masks, autoimage + anchor, fit + mask, fit reference + frame, reference file, remove box, trajin start/stop/offset, extra commands); *Run cleaning now*, *Inspect system* | app.py:422-559, 1440-1460, 1786-1831 | 5 (cleaning per D12) |
| Residue groups (`groups`) | table name / residues / selection / colour + swatch column, *Pick colour…*, add / remove / move, hint on `group:` | app.py:561-637 | 5 |
| Analyses (`analyses`) | General (`use_pbc`, `proximity_search_radius`); RMSD (targets table, reference frame, superpose, mass weighted); contacts, H-bonds, Mg²⁺, stacking, base pairs (engine `builtin` / `barnaba`), energy (target value) — each section's fields and group combos | app.py:639-756 | 5 |
| Structure view (`structure`) | Orientation (3 sliders −180..180, angle label, zoom, centre on, follow, reset buttons); Representations (rows: style, groups, colouring, size, opacity, sugars, remove; rest / backbone atom / hydrogens; renderer `auto` / `software` / `vtk`); Density map (file, show / colour / opacity / front, level / resolution, fit frame / search / step, auto fit, *Global search* / *Refine* / *Centre on model* / *Clear*, *Move the map with the mouse* toggle, status line); Residue number tags; Appearance (groups shown, trace atoms, align / ions / depth cue, colour mode, base / ion colours, chain palette, proximity set / cutoff / colour, glow mode / colour, fade / strength, line width / point size, theme / font size, time label / title) | app.py:758-1107 | 5; density fit buttons M6 (m4 D9) |
| Layout & panels (`layout`) | grid preset; one slot editor per slot (panel type, data key with label, group A / B combos with `(analysis default)`, more keys, title (automatic until typed), options JSON); *Render still* of frame N | app.py:1109-1334 | 5 |
| Outputs (`outputs`) | Locations (out dir, cache dir, params log); Video (enabled, name, size, dpi, fps / stride, codec / CRF, ffmpeg path); Time charts (enabled, series, heatmaps, formats PNG / SVG / PDF check boxes, size, dpi, style (mdatk + Matplotlib styles), colour, line width, title, axis labels, smoothing, target, prefix); Frame PDBs (enabled, frames, times, subset, prefix, *Export PDBs now*, *Re-render video + charts from cache*) | app.py:1336-1419 | 5 |

`collect()` copies the project, stores every binding and (optionally) validates (app.py:1462-1474);
unbound keys keep their values. Panel-type slot defaults: app.py:66-69.

### 1.3 Canvas layout editor (`gui/canvas.py`)

`set_content` / debounced `schedule_render` 120 ms (canvas.py:95-112); `render` builds the figure at
the widget size with the video's aspect and `dpi * scale`, updates panels at 35 % of the example
time axis, blits the moving artists (114-149); slot hit areas from each panel's tight bounding box,
empty slots from the grid (151-178). Mouse: press in the structure panel below the title strip
(12 %) = rotate (0.5°/px trackball, `on_rotate` writes the sliders), elsewhere / title strip / Shift =
move panel (drag threshold 4 px, target slot highlighted with dashed outline and type label, release
= swap), map mode: drag moves the density map, Shift-drag rotates it (`nudge_density`, `on_map_move`
stores transform + anchor + `fit_info`), double-click = open the slot's editor, wheel over the
structure = zoom × 1.15 per notch; cursors `grab` / `all-scroll` / `move`; fast quality while dragging,
full quality 250 ms after (208-345). Messages `Cannot draw this layout: …` (130).

### 1.4 Dialogs, jobs and desktop

| Item | MDATK | Milestone |
|---|---|---|
| File dialogs and filters: topology, trajectory, reference, data, map, project (`*.mdatk.toml`), parameter logs, BibTeX | app.py:56-63, widgets.py:246-312 | 5 |
| Save / Save as: overwrite confirmation, default name `<workflow.project_name or MDATK>.mdatk.toml`, `.mdatk.toml` appended when missing | app.py:1509-1543 | 5 |
| Project log dialog (`time  step  [status]  title` + summary) | app.py:1545-1572 | 5 |
| Demo system (folder chooser, background job, opens the project) | app.py:1574-1593 | 5 (D13) |
| About and citations (text of `project_resources` grouped, *Export BibTeX…*) | app.py:1648-1720 | 5 |
| Open output folder (FileManager1 D-Bus, `inode/directory`, xdg-open; Explorer), open latest video (`latest_video` by mtime) | desktop.py, app.py:1722-1752 | 5 |
| Jobs: one at a time (`Another job is still running.`), worker thread, throttled progress (0.5 %), log lines, error dialog + traceback in the log, `Cancelled` | worker.py, app.py:1911-1942 | 5 |
| Render / Run cleaning / Inspect system / Load preview / Render still (temp PNG shown in a dialog) / Export PDBs / Re-render from cache; result dialogs (`Render finished`, file counts) | app.py:1758-1909, 2030-2044 | 5 |
| Error / info: modal message dialogs, detail cut at 3 000 characters | app.py:2012-2028 | 5 |
| Recent projects: `$XDG_CONFIG_HOME/mdatk/settings.toml` (else `~/.config`), 10 entries, existing files only | projectfile.py:528-554 | 5 (core part exists: `core/settings.hpp`) |
| Run generated scripts (confirmation, runs logged in the project) | app.py:1944-1985 | 7 |
| Project Planner window (canvas of steps, ports, templates, Delete key, zoom, details pane, write / run) | planner.py (1 056) | 7 |
| Script wizard (schema-built forms, records editor, live preview, write / run) | wizard.py (910) | 7 |
| PDB clean-up window / structure builder window | pdbdialog.py, builderdialog.py | 6 |
| `apply_inputs` (planner's files as analysis inputs) | app.py:1611-1630 | 7 |

No preferences dialog, no settings beyond the recent list, no undo, no auto-save, no translations
in MDATK.

### 1.5 Command line: MDATK commands vs MDBridge

MDBridge has `init` (`init-params`), `info`, `list`, `analyze`, `export-pdb`, `render`, `charts`,
`preview`, `density info` (MDATK's names and options; `cli_parity.py` holds about 120 cases over
them) and its own `gui`,
`project show`, `graphics`, `fonts`, `probe-video` (src/app/cli.cpp).

| MDATK command (cli.py) | Milestone |
|---|---|
| `gui [params]` (794) | 5: MDBridge's `gui [project]` already takes it; MDATK's positional name kept |
| `demo OUT_DIR [--n-frames 120] [--dt 10] [--no-solvent] [--run] [-q]` (798) | 5 (D13) |
| `clean [common] [--force]` (820) and cleaning inside `analyze` / `render` | 5 or later (D12) |
| `density fit`, `density simulate` (869-887) | 6 (m4 D9) |
| `pdb clean/check/annotate`, `build …`, `dict …` (892-932) | 6 |
| `workflow tools/check/write/run/generate/template/restraints` (934-992) | 7 |

Citations are not a command in MDATK (GUI only: About, *Export citations (BibTeX)…*); D14.

### 1.6 MDATK tests that specify the app layer

`test_gui.py` (13: widgets round trip, canvas shows layout, drag rotates and swaps, panel group
drop-downs and automatic titles, swatch, real-frame preview and render, density fit / shadow / drag,
representation editor and 3D preview, PDB tools window (M6), citation text, slot defaults and wheel
zoom, value widgets ignore the wheel, builder helix (M6)); `test_workflow_gui.py` (5, M7);
`test_cli.py` (`test_apply_override`, `test_cli_demo_render_and_tools`, `test_cli_errors`). Each
M5 test gets a Qt Test counterpart (section 3.3).

---------------------------------------------------------------------------------------------------

## 2. What MDBridge already has

- `src/app/gui` (M1 / M4): `MainWindow` (splitter: `ProjectPanel` sidebar with a read-only summary,
  settings tree and step log; right: the 3D view above `PreviewWidget`), File (Open, Quit), View
  (software / Vulkan view), Help (Graphics, About); drag-and-drop of project files; status bar with
  progress; the reference frame loads on a worker thread (`Background`, `background.hpp`) with
  cancellation by generation; `PreviewWidget` draws the figure on a worker thread, debounced.
- `src/app/commands/project_view.*`: `loadReferenceFrame`, `structurePreview`, `structureView`,
  `previewFigure` (MDATK's canvas content, no Qt widgets) — the canvas editor builds on it.
- `src/render`: figure model with `tightBox(ax)` (slot hit areas), `StructurePanel::setRotation`,
  `setFast`, `nudgeDensity`, `orientation.hpp` `drag()`, example dataset, Vulkan and software views
  (`render::View`), figure → QImage (`figure_qt.cpp`).
- `src/core`: the project model and file writer (byte parity), `applyOverride`, `validate`,
  `citations` (`keysFor`, `withRunKeys`, BibTeX texts), `settings.hpp` (recent list in MDATK's
  file), step log, `Reporter` / `Cancelled`.
- CLI framework (`cli.cpp`, `commands.cpp`) with MDATK's exit codes and messages; `cli_parity.py`.
- CI: offscreen-platform Qt tests, `xvfb-run` for Vulkan window tests, `packaging/smoke_test.py --gui`
  screenshots on Ubuntu, bare Rocky and Windows.

Missing for M5: every editable page, the editor canvas, jobs from the window, dialogs, recent
projects in the GUI, demo / clean / citations commands, GUI tests that compare project files.

---------------------------------------------------------------------------------------------------

## 3. Parity strategy

### 3.1 What is compared

| Artifact | Same means |
|---|---|
| Project file written after a GUI action script (open, edit fields, add rows, swap / rotate / zoom / map drag, preset change, save) | byte for byte with MDATK's for the same actions (time stamps fixed as `make_project_file_golden.py` does), after the AP patches (`improvements_gui.py`) |
| `collect()` result without edits | the opened project (lossless round trip; MDATK's `test_widgets_roundtrip_project` only checks the demo) |
| Layout editor results | `layout` table of the saved file (panels, series, titles, options), `visual.rotation` / `zoom`, `density.transform` / `anchor` / `fit_info` |
| Automatic titles, data-key choices and labels of slot editors, group choices | equal strings per (panel type, project) |
| Messages shown (log lines, dialog titles / texts, status) | MDATK's text where MDATK has one, MDBridge's documented text where an AP changes it |
| New CLI commands (`demo`, `clean`, …) | as `cli_parity.py`: stdout, stderr, exit code, files byte for byte (demo trajectory: float32 bytes) |
| Pixels of the window | not compared (D10 of M4: no golden images); screenshots are CI artifacts and the smoke test checks content (preview drawn, structure shown) |

### 3.2 Generators (`tests/parity/`, `regenerate_all.py`)

- `make_gui_cases.py`: imports MDATK v0.0.4.1's `mdatk.gui.app` under Xvfb (GTK 3, PyGObject,
  pycairo), opens each fixture project (`tests/data/projects`, demo, MDATK's examples), replays action
  scripts through the window's own methods and widgets (`binder` setters, `_on_canvas_swap`,
  `canvas.set_rotation`, `_on_map_moved`, `_rebuild_slots`, `save_project`), writes the files and a
  JSON of `(action script, file text, log lines)` to `tests/data/gui/`. Checked during planning: the
  v0.0.4.1 GUI runs headless this way (probes in section 4).
- `improvements_gui.py`: test-side patches of MDATK for the AP items (as `improvements_*.py`).
- `cli_parity.py` gains `demo` (and `clean` with D12) cases.
- Regenerate job: needs GTK 3 introspection data and PyGObject / pycairo for Python 3.13
  (pip-pinned in `requirements-regenerate.txt`; Ubuntu's `python3-gi` is built for 3.12) and `xvfb`.

### 3.3 GUI tests in CI

- Qt Test (`Qt6::Test`, LGPL, test-only) in `mdbridge_gui_tests`, on the `offscreen` platform on Linux
  and Windows (no display); Vulkan-view cases keep `xvfb-run` + lavapipe as now.
- Every widget is reachable by `objectName` = project key (`analyses.contacts.cutoff`), so tests set
  values as a user would (`QTest::keyClicks`, mouse events on the canvas at slot centres) and compare
  the saved file with `tests/data/gui/*.json`.
- Counterparts of MDATK's tests (section 1.6) plus one per AP item; `smoke_test.py --gui` extended to
  open, edit, save and render the demo from the AppImage / installer; screenshots of every page and
  of a canvas drag uploaded as artifacts.
- The core edit model (W0-B) is tested without Qt (GoogleTest), so most logic runs on every
  compiler without a display.

---------------------------------------------------------------------------------------------------

## 4. MDATK app-layer bugs and fragile behaviour to fix (ids AP-xx)

Found by reading v0.0.4.1 and running its GUI headless (Xvfb, probe scripts that open a project,
call `collect()` / `save_project()` and compare). "0.1.1" says whether MDATK's later `main` still has
it (same probes run against it) and how it was fixed; where it was, MDBridge follows that design.

| Id | MDATK v0.0.4.1 behaviour (verified) | MDBridge | 0.1.1 | Status (DIFFERENCES ids) |
|---|---|---|---|---|
| AP-01 | Spin buttons clamp loaded values to their widget ranges, so opening and saving rewrites valid settings: `contacts.cutoff` 25 → 20, `hbonds.angle_cutoff` 60 → 90 (validation allows 0-180), `energy.cutoff` 60 → 50, `video.width` 7682 → 7680, `visual.zoom` 0.01 → 0.05 (widgets.py:81, 342, 345; ranges app.py:639-756, 1336-1419) | no clamping on load: fields hold any value of the setting's type; ranges come from `validate()` (PRJ-01) and out-of-range values are flagged, never rewritten | not fixed (same probe) | done: ED-01, FW-01, PG-01 |
| AP-02 | Optional numbers shown with `%g`: `inputs.dt_ps` 0.0123456789 → 0.0123457, `energy.target_value` −1234567.5 → −1234570 on save (widgets.py:375; wizard.py:126 likewise, M7) | shortest round-trip text (the project writer's float format) | not fixed | done: ED-01, FW-02 |
| AP-03 | Colour buttons parse with GDK and write `#rrggbb`: Matplotlib colours MDATK itself accepts become black (`tab:orange` → `#000000`, `#ff000080` → `#000000`) and names are rewritten (`gold` → `#ffd700`) (widgets.py:38-49, 361) | colour field keeps the text; the picker writes `#rrggbb`; validity by MDBridge's Matplotlib colour grammar (RN-12) | not fixed | done: ED-01, FW-02 |
| AP-04 | Chart formats are three check boxes: `formats = ["jpg", "svg"]` saves as `["svg"]`, none ticked saves `["png"]` (app.py:1371-1376) | every format the charts support (RN-11) selectable; unknown entries kept and flagged | not fixed | done: FW-03, PG-02 |
| AP-05 | Fixed combos keep the previous project's choice when a value is not in the list (`visual.renderer = "vulkan"` saved as `software` after a project with `software`); first time: `""` (widgets.py:353) | an unknown value stays, shown as `<value> (unknown)` and flagged by the check; never replaced | not fixed | done: FW-03 |
| AP-06 | Saving from the GUI writes every path absolute (`topology = "demo.parm7"` → `/…/demo/demo.parm7`, `out_dir`), so a project folder that is moved or shared breaks after one save; `init` keeps them relative (cli.py:214) | paths written as the file had them (relative to the project folder); Save as to another folder rebases relative paths — D8 | not fixed | done: ED-02, FW-08 |
| AP-07 | Empty slots become panels on save: a 4-slot preset with 2 panels is saved with 4 (two `contacts/atom_pairs` series without titles) (app.py:1282-1320) | empty slots stay empty (drawn as placeholders); `layout.panels` keeps its length | not fixed | done: LE-05 |
| AP-08 | Changing the preset to fewer slots drops the extra panels at once and for good (equal_4 → single → equal_4 loses panel 2) (app.py:1137-1141, 1282) | panels beyond the slots are kept (render's existing `ignoring M extra panel(s)` warning) until the user removes them | not fixed | done: LE-05 |
| AP-09 | A panel swap or preset change while a slot's options JSON is invalid raises out of the signal handler (swap silently not done) (app.py:356-363) | the action is refused with the located message | fixed in 0.1.0 (c6d8d15, gui/app.py): adopt | done: LE-05 |
| AP-10 | Form errors are Python texts: `inputs.stop: invalid literal for int() with base 10: 'x'`, `layout: panel options '{bad': Expecting property name …` | located problems (`Input & cleaning › Frames › Stop: 'x' is not a whole number`), the field flagged | partly (0.0.10 `problems.py` locates selection / group problems; the int text remains): adopt the 0.0.10 design for every field | done: ED-03, ED-04, FW-04, LE-05 |
| AP-11 | *Save as* with a name typed without `.mdatk.toml` appends the suffix after the overwrite confirmation: an existing file is replaced unasked (app.py:1531) | the confirmation sees the final name | fixed in 0.1.0 (172df88, gui/app.py): adopt | done: WN-03 |
| AP-12 | A project that fails to load at start is reported on stderr only; the window opens defaults (app.py:2062) | named in the window's log and status bar (MDBridge M4 already prints it) | fixed in 0.1.0 (13d4fea): adopt | done: WN-09 |
| AP-13 | Closing the window, *New project* or opening another project discards unsaved edits without asking; nothing marks a modified project (app.py:1476-1507, 1997) | modified marker in the title; ask (Save / Discard / Cancel) — D6 | not fixed (`_on_close` only stops jobs) | done: WN-03 |
| AP-14 | Closing during a job cancels the worker but does not wait for it (daemon thread, worker.py:39); the process ends mid-write | ask; on confirmation cancel and wait for the job's own cleanup (temporary video files) before quitting | RNAdo asks (0.0.11); the analysis window does not | done: JB-01 |
| AP-15 | Table cells (groups, trajectories, external data, RMSD targets) drop text typed but not confirmed with Enter when the cell is left by a click | edits commit on focus loss (Qt delegates; tested) | fixed in 0.1.0 (608f322, gui/widgets.py): adopt | done: FW-05 |
| AP-16 | About and *Export citations* cite the settings' resources only | settings' citations followed by the run-time ones of the cached results (`core::citations::withRunKeys`, CIT-06) | MDBridge-specific | done: WN-07, WA-09 |
| AP-17 | Updating the recent list rewrites `settings.toml` with `recent` only, non-atomically (projectfile.py:542-554); MDBridge's `core::addRecent` (src/core/settings.cpp:78) does the same, so it deletes settings newer MDATK keeps there (0.0.14 cat-proof mode) | keep every other key; write through a temporary file and rename | fixed in 0.0.14 (`projectfile.py` `_user_settings`, other keys kept; still a plain write): adopt and make atomic | done: PRJ-12 |
| AP-18 | Dialog details cut at 3 000 characters, not selectable; the still preview scaled to a fixed 1 200 px (app.py:2017, 2026, 2036) | full text in the log, selectable scrollable details; the still fits the screen | selectable texts in 0.1.0 (`dialogtext.py`): adopt | done: WN-07, WA-07 |
| AP-19 | Choices MDBridge does not have are offered: base-pair engine `barnaba`, renderer `vtk` as VTK, all Matplotlib chart styles (app.py:728-731, 1086-1090, 1382-1388) | MDBridge's choices (renderer adds `vulkan`, RN-20; D5 of M4 styles; Barnaba never bundled) with stored values kept (AP-05) | n/a | done: FW-03, PG-03, GP-05 |
| AP-20 | The canvas builds and draws the figure on the GUI thread; each wheel notch rebuilds the whole figure through the zoom spin (canvas.py:114, app.py:338) | drawing on a worker (done in M4's preview); zoom and rotation update the structure panel only, full rebuild debounced | not fixed | done: LE-01, LE-02 |

Not bugs, kept: page order and names, field order, hint texts, default file names, the 35 %
preview position, the 120 ms / 250 ms timings, the `(analysis default)` group entry, automatic
titles that follow the data until typed.

---------------------------------------------------------------------------------------------------

## 5. MDATK after 0.0.4.1 (to 0.1.1): what the app layer may want

Reported, not ported. Parity tests stay on v0.0.4.1 (D18); items adopted in M5 are AP fixes or
features whose expected data the test side patches in, naming the MDATK version.

### 5.1 App-layer features added after 0.0.4.1

| Feature | Version (CHANGELOG / file) | Class |
|---|---|---|
| Value fields ignore the mouse wheel | 0.0.4 (already in v0.0.4.1: widgets.py:55-79) | M5 (parity) |
| Settings checked as you type: page badges, flagged fields, problems list that jumps to the field; located messages | 0.0.10 (`problems.py`, `gui/problemview.py`) | adopt in M5 (the form layer is built once) |
| "Only what is about to run can block it": errors vs warnings by use | 0.0.10 | later (post-parity: changes `validate()` output the CLI parity compares) — D7 |
| Residue-group renames cascade through the project; delete says where a group is used | 0.0.10 (`Project.rename_group`, `group_references`) | adopt in M5 (group editor; otherwise reworked) |
| Groups from chains (on choosing a topology, *Groups from chains*) | 0.0.10 (`chaingroups.py`) | adopt in M5 if cheap (needs the topology's chains, which `io` has) — D15 |
| *Close project* | 0.0.13 | adopt in M5 |
| Open drop-down fix, selectable dialog texts, About dialog shared | 0.0.13, 0.1.0 (`about.py`, `dialogtext.py`) | adopt in M5 |
| Output log panel: per-step marks, stderr red, ERROR / WAITING styles, bounded batched inserts, Ctrl+L fold | 0.0.11 (`logpanel.py`) | M5 for the bounded / batched view and styles; per-step filters M7 |
| RNAdo as a program of its own (planner, wizard, PDB clean-up, builder); `mdatk gui --rnado/--open`, `rnado` | 0.0.11 (`gui/rnado.py`, `launch.py`) | M7 — D1 |
| Co-editing: each program notices the other's saves (file monitor), merges under a file lock, change bar | 0.0.11 (`projectui.py`) | M5: notice external changes and offer reload (QFileSystemWatcher); the merge M7 |
| *Use the plan's files* | 0.0.11 | M7 |
| Desktop entries: `mdatk desktop install/uninstall/status`, menu item | 0.0.12 (`launchers.py`) | later (packaging; AppImage users have no menu entry until then) |
| Cat-proof mode (keyboard ignored), stored in `settings.toml` | 0.0.14 (`gui/catproof.py`) | later (RNAdo, M7) |
| Cluster profiles window, `mdatk cluster list/learn`, `<config>/mdatk/clusters.toml` | 0.1.0 (`clusterdialog.py`) | separate post-parity phase |
| CLI `workflow status/jobs/manifest`, `--parallel`, `--include-separate` | 0.0.7-0.1.0 | M7 / post-parity |
| CLI `pdb qc`, `decoys`, `rosetta`, `frame`, `ensemble`, `map`; `build chimera`, B-DNA | 0.0.12-0.0.18 | separate post-parity phase (M6 / M7 domains) |
| Windows cursor workaround | 0.0.10 (`cursors.py`) | n/a (Qt) |

The analysis-side commands (`init`, `info`, `analyze`, `render`, `charts`, `preview`, `export-pdb`,
`clean`, `demo`, `density`) keep v0.0.4.1's options in 0.1.1 (cli.py:1184-1330 there).

### 5.2 Project files after 0.0.4.1

`format_version` stays 1; only keys were added, all in the workflow part: `workflow.supersede_runs`,
`parallel` (missing = 1), `cpu_cores`, `gpus`, `last_node`, per node `adopted` and `folder_tag`,
plus new options of step types (cluster job, NAMD MDFF, …). The analysis settings have no new keys
(types only annotated). `mdatk_version` says `0.1.1`. MDBridge reads such a file today with
`WARNING: project.workflow: ignoring unknown key 'parallel'` (one per key) and drops the keys when it
writes (checked with `init -p`): a GUI save would destroy a 0.1.1 user's plan settings. M5 keeps
unknown keys and writes them back unchanged (D9); full support of the new keys comes with M7.

### 5.3 Behaviour changes outside the app layer the owner may want to map

- 0.1.1 generic Lennard-Jones table (`docs/GENERIC_LJ.md`, `analysis/generic_lj.py`): one parm10
  type per element, hydrogens typed by the bonded atom, lone atoms as Joung-Cheatham / Li-Merz 12-6
  ions. MDBridge's AN-25 chose OpenMM amber14 residue templates instead; the two now differ for every
  non-parm7 energy (demo duplex, strand against strand, vdW in kcal/mol: MDATK 0.1.1 −34.1, MDBridge
  −34.62 from the templates (DIFFERENCES, energy), the parm7's own parameters −34.67; v0.0.4.1's
  table +145.4) — D17.
- 0.1.0 analysis fixes overlapping MDBridge's register: triclinic minimum-image vectors (= AN-04),
  H-bonds by distance counted twice, `inputs.dt_ps` ignored for time-stamped trajectories, contacts of
  a group with itself, CONECT-only-ligand chains, chain-qualified tags (= RN-24), CHAMBER / CMAP parm7,
  external files of one name. Each needs a check against MDBridge (not M5 work).

### 5.4 MDATK's `docs/ARCHITECTURE_CPP.md` (0.1.1) vs MDBridge's `docs/ARCHITECTURE.md`

It is the plan as of 0.0.5 and conflicts with MDBridge's current plan in: "no C++ code yet";
"MDATK frozen, bug fixes only" (MDATK has had feature releases 0.0.5-0.1.1); chemfiles for reading
(dropped in M2); "e.g. QCustomPlot" (GPL, excluded in M4); no "parity is the floor" rule, no Windows
path rules; one GUI with "Widgets for the sidebar, wizard and planner" (MDATK has since split RNAdo
into its own program, D1). No conflict in modules, delivery, compilers or migration order.

---------------------------------------------------------------------------------------------------

## 6. Work breakdown

Code layout: `src/core` gains the edit model (no Qt); `src/app/gui/` gains `pages/` (one file per
page), `forms/` (bound field widgets), `canvas/` (editor), `jobs/` (job runner, log view), `dialogs/`;
`src/app/commands/` gains `demo.cpp`, `clean.cpp` (D12), `citations.cpp` (D14).

### Wave 0 (independent)

| Chunk | Content | Acceptance |
|---|---|---|
| **W0-A Generators** | `make_gui_cases.py` (action scripts over the fixture projects: per page, per AP, canvas actions), `improvements_gui.py`, `cli_parity.py` demo / clean cases; regenerate job with GTK + PyGObject for 3.13 + xvfb | `regenerate_all.py --only gui` byte-identical twice |
| **W0-B Edit model** (core) | `ProjectDocument`: the project value, its file path, field paths (`analyses.contacts.cutoff`, `layout.panels.2.options`), set / get by path with type checks, modified flag vs last saved, undo / redo of snapshots (D5), located problems (path → page, section, label; 0.0.10 texts), lossless numbers / colours / unknown enum values (AP-01..05), unknown keys kept (D9), paths kept relative (AP-06, D8), group rename cascade (0.0.10), settings file keeps keys + atomic write (AP-17) | GoogleTest: round trip of every fixture project and a 0.1.1 file byte for byte; AP-01..06, AP-17 cases |
| **W0-C Field widgets** (Qt) | bound widgets: number (no clamping, range hint), optional number, text, path + file dialog with MDATK's filters, combo (+ unknown value), editable combo, colour (text + picker + swatch), check, multi-line text, CSV, list table (QTableView, commit on focus loss AP-15, add / remove / move), problem flag + tooltip; wheel ignored (0.0.4); `objectName` = key | Qt Test per widget: load → save is identity; AP-15 |
| **W0-D Jobs and log** | job runner on a worker thread (one at a time, `Another job is still running.`), `core::Reporter` progress throttled as MDATK, cancel, result / error delivery, close-while-running (AP-14); log view (bounded 5 000-6 000 lines, batched, styles of 0.0.11, Clear, Copy, Save, fold Ctrl+L) | tests: 10^6 lines keep the event loop responsive (< 0.2 s stalls); cancel paths |
| **W0-E Commands** | `demo` (port of `demo.py`: numpy PCG64 normal exists, `shuffle` to port; parm7 / NetCDF writers exist), `citations` (D14), `clean` (D12: cpptraj script writer, `conda run`, cleaning key and reuse, MDATK's messages) | `cli_parity.py` demo / clean cases; MDATK's `test_cli_demo_render_and_tools` counterpart |
| **W0-F GUI test harness** | `mdbridge_gui_tests` (Qt Test, offscreen, Linux + Windows), helpers to run an action script against the window and compare with `tests/data/gui`, screenshot artifacts | a trivial open → save case green on both platforms |

### Wave 1

| Chunk | Needs | Content | Acceptance |
|---|---|---|---|
| **W1-A Window shell** | B, C, D | MDATK's layout (D4): header / menus (Open + recent + demo, Save, Save as (AP-11), New, Close project, Project log, citations, About), page buttons with problem badges, side bar with scroll-to-section, log, action bar, status; shortcuts (D11); modified marker and close / open / new questions (AP-13, D6); start-up load errors (AP-12); external-change notice (5.1); drag-and-drop | Qt Test: each menu action; AP-11..14 |
| **W1-B Pages 1** | B, C | Input & cleaning, Residue groups (swatches, colour picker, rename cascade), Analyses | per-page action scripts equal MDATK's files; `test_panel_group_dropdowns…` parts |
| **W1-C Pages 2** | B, C | Structure view (orientation sliders + angle label, representations rows, density section with fit buttons as shells (M6), tags, appearance), Outputs (formats AP-04, styles AP-19) | action scripts; `test_representation_editor…`, `test_value_widgets_ignore_the_mouse_wheel` |
| **W1-D Problems** | B, A (badges) | check as you type (debounced), problems list jumping to fields (0.0.10 presentation, MDATK's v0.0.4.1 messages where they exist) | every `validate_cases.json` problem lands on its field |

### Wave 2

| Chunk | Needs | Content | Acceptance |
|---|---|---|---|
| **W2-A Layout page + canvas editor** | W1-A, W0-B; M4 render APIs | slot editors (type, data key + label, group combos, more keys, automatic titles, options JSON with AP-09), presets (AP-07, AP-08); editor canvas on `PreviewWidget`'s worker drawing: slot hit areas from `Figure::tightBox`, drag-to-swap with highlight and label, rotate (`orientation::drag`, `StructurePanel::setRotation`, fast quality, settle 250 ms), wheel zoom (AP-20), map move / rotate (`nudgeDensity`, transform + anchor + `fit_info` texts), double-click → slot editor, cursors; 3D view placement (D4) | `test_canvas_shows_layout…`, `test_drag_rotates_structure_and_swaps_panels`, `test_density_map_fit_shadow_and_drag` (drag part), `test_slot_defaults…`; saved `layout` / `rotation` / `zoom` / `density` equal MDATK's for the same mouse paths |
| **W2-B Actions** | W0-D, W1-A..C | Render, Cancel, Render still (preview dialog), Export PDBs, Re-render from cache, Inspect system, Run cleaning (D12), Create demo, Open output folder / video (`QDesktopServices`, `latest_video`), result dialogs, About + BibTeX export through `withRunKeys` (AP-16, AP-18) | `test_real_frame_preview_and_render`, `test_citation_text…` counterparts; smoke test renders the demo from the window |

### Wave 3

Docs (DIFFERENCES "Interface" with every AP item and test, README, BUILDING (Qt Test), LICENSING if
anything is added), accessibility pass (D10: keyboard-only walk through every page, accessible names,
high-DPI at 100 / 150 / 200 %), Windows checks (paths in dialogs, `%APPDATA%` per D16, installer
start-menu entry opens a project), `smoke_test.py --gui` editing case on bare Rocky and Windows,
screenshots of every page as artifacts.

Status of W3-A: done. Paths in interface text natively on Windows (`gui::displayPath`, WN-12);
accessibility (D10: Tab reaches and names every control, Tab order per page, canvas keys, Esc,
150 % / 200 % runs; WN-13); settings locations checked (D16, WN-14); DIFFERENCES "Interface" with
an index of the ids and every AP item, README, BUILDING, ARCHITECTURE; the AP table's status column.

### Critical path and parallelism

`W0-B → W0-C → W1-B/C → W2-A` is the critical path (the edit model decides every widget's round
trip). Second chain `W0-D → W2-B`. W0-A (generators, GTK in the regenerate job) has long lead time;
start it first. Six agents in wave 0, four in wave 1, two in wave 2, one or two in wave 3.

Chunks needing owner decisions first: W0-B (D5, D6, D8, D9), W1-A (D1, D2, D3, D4, D11, D16), W0-E
(D12, D13, D14), W1-B (D15), W1-D (D7), W3 (D10).

---------------------------------------------------------------------------------------------------

## 7. Owner decisions (one line each, with recommendation)

All decided by the owner on 2026-11-25: every recommendation below is accepted; D4 as amended.
Also accepted: the pre-M5 fix wave from the MDATK 0.1.1 comparison (section 5).

- **D1 RNAdo split**: recommend M5 builds the analysis window only; in M7 the planner, wizard, PDB
  clean-up and builder become a second program `rnado` (AppImage action, Start-menu entry) as MDATK
  0.0.11 did.
- **D2 Shells for later milestones**: recommend menu entries for planner, wizard, PDB clean-up,
  builder and the density fit buttons present but disabled, tooltip "comes with milestone 6 / 7"; no
  empty windows.
- **D3 Look**: recommend Qt's Fusion style on both platforms with the system light / dark scheme
  (Qt 6.5+) and native file dialogs (XDG portal on Linux).
- **D4 Window layout**: MDATK's (page buttons, 520 px side bar, canvas right, log fold, action bar).
  **Decided:** the separate live 3D view of M1/M4 is removed; the canvas's structure panel is the
  3D view (rotate, zoom, map move in place, as MDATK), drawn by the Vulkan renderer while dragging
  when a GPU is available, else the software renderer. A dedicated 3D window may return later if it
  earns a purpose (builder, M6; very large systems).
- **D5 Undo**: recommend undo / redo (Ctrl+Z, Ctrl+Shift+Z) of whole-project snapshots for every edit
  and canvas action (a drag = one step), not persisted (MDATK has none).
- **D6 Unsaved changes and auto-save**: recommend a modified marker, the Save / Discard / Cancel
  question (AP-13), and a recovery copy every 2 minutes in the user config folder offered at the next
  start; the project file itself is never written without Save.
- **D7 Validation semantics**: recommend the 0.0.10 presentation (badges, located messages) in M5 and
  keeping v0.0.4.1's blocking rules until the post-parity phase (they change CLI output).
- **D8 Paths on save (AP-06)**: recommend writing paths as the file had them (relative stays
  relative; Save as rebases), differing from MDATK's GUI.
- **D9 Newer project files**: recommend keeping unknown keys and writing them back unchanged, with
  one warning on open; no refusal to save.
- **D10 Accessibility and high-DPI**: recommend keyboard access to every control, accessible names
  on icon buttons (as MDATK), Qt's per-monitor DPI scaling with the preview drawn at the device pixel
  ratio; no screen-reader testing in CI.
- **D11 Shortcuts**: recommend MDATK's (Ctrl+O, Ctrl+S, Ctrl+Shift+S, Ctrl+R; Ctrl+P / Ctrl+Shift+W
  reserved for M7) plus Ctrl+N, Ctrl+W (close project), Ctrl+Q, Ctrl+Z / Ctrl+Shift+Z, Ctrl+L (log,
  as RNAdo), F5 (refresh preview), Esc (cancel a drag), Ctrl+1..6 (pages).
- **D12 Cleaning (cpptraj) in M5**: recommend yes: `clean`, cleaning in `analyze` / `render` and the
  GUI's *Run cleaning now* / *Detect*, with a process runner M7 reuses (today MDBridge only reuses
  MDATK's cleaned files).
- **D13 `demo` in M5**: recommend yes (GUI's *Create demo system* and *Demo project…*; writers and
  the numpy random port exist).
- **D14 Citations on the command line**: recommend a new `mdbridge citations [-p PROJECT] [--bibtex
  FILE] [--json]` citing settings and run-time keys (`withRunKeys`, CIT-06); name not used by MDATK.
  **Format** (W0-A; the expected output of `cli_parity.py`'s `citations_*` cases is written from
  this specification, `tests/parity/citations_command.py`): options of every analysis command
  (`-p`, `--top`, `--traj`, `-o`, `--cache-dir`, `--frames`, `-s` ...) plus `--bibtex FILE` and
  `--json`. Keys = the settings' keys, then the `citation_keys` of the analyses the cache manifest's
  last run lists (none without a manifest), through `withRunKeys`. Plain: MDATK's About text
  (`Resources this project uses (N references; each generated script lists only its own):`, per
  resource its name and ` - note`, references indented by four spaces, empty lines between), then
  `Recorded by the cached results of the last run:` with the references only the cache adds; no key:
  `This project uses no external resource with a reference yet.` `--bibtex FILE`: MDATK's BibTeX of
  the keys, `Wrote N references to FILE`; no key: `Nothing to cite yet: The project uses no external
  resource with a reference.` and no file. `--json`: `{"keys", "settings_keys", "run_keys",
  "resources": [{"id", "name", "note", "keys"}]}` (indent 2, non-ASCII as is) instead of the text or
  the `Wrote` line. Errors: `error: <message>`, exit 1; usage errors exit 2.
- **D15 Groups from chains (0.0.10) in M5**: recommend yes, with the rename cascade; both are
  additions beyond v0.0.4.1, listed in DIFFERENCES.
- **D16 Settings location**: recommend Linux: keep the shared `$XDG_CONFIG_HOME/mdatk/settings.toml`
  (recent list shared with MDATK) plus MDBridge-only settings in `$XDG_CONFIG_HOME/mdbridge/`;
  Windows: `%APPDATA%\MDBridge\` (MDATK's Windows beta uses `%USERPROFILE%\.config\mdatk`, not shared).
- **D17 Generic Lennard-Jones (not M5)**: MDATK 0.1.1 and MDBridge AN-25 now use different documented
  tables; recommend keeping AN-25 and deciding before the post-parity phase.
- **D18 Parity reference**: recommend v0.0.4.1 stays the reference for all tests; adopted later fixes
  (AP-09, -11, -12, -15, -17, -18 and section 5.1 items) are applied test-side with the version named.
- **D19 Translations**: recommend English only, all strings through `tr()` so they can be added later.
- **D20 Crash reports**: recommend local only: on a crash write a report (stack trace; minidump on
  Windows) to the user config folder and offer it at the next start for sending by e-mail; no
  network reporting.

---------------------------------------------------------------------------------------------------

## 8. Risks

1. **GTK in the regenerate job**: MDATK's GUI needs GTK 3 typelibs and PyGObject / pycairo built for
   the job's Python 3.13 (Ubuntu ships them for 3.12); building PyGObject from pip needs the GTK
   development packages. Mitigation: start W0-A first; fall back to a 3.12 venv with system packages
   for the GUI generator only (the GUI's outputs do not depend on numpy versions beyond the pinned
   ones).
   **Solved (W0-A):** the regenerate job builds PyGObject 3.52.3 and pycairo 1.28.0 from source for
   its Python 3.13 venv against Ubuntu 24.04's GTK (`libgirepository-2.0-dev`, `libcairo2-dev`,
   `gir1.2-gtk-3.0`; pinned build tools, `tests/parity/requirements-gui-build.txt` /
   `requirements-gui.txt`); `make_gui_cases.py` starts Xvfb itself. Its output is the same with the
   system's PyGObject 3.48 on Python 3.12 (checked byte for byte).
2. **Round-trip completeness**: one field bound lossily breaks byte parity of every saved file.
   Mitigation: W0-B tests every fixture and a 0.1.1 file byte for byte before widgets exist.
3. **Two programs, one settings file** (D16): MDATK and MDBridge both write `settings.toml`; without
   AP-17 they erase each other's settings, and without locking concurrent writes can still race.
4. **Canvas responsiveness**: interactive rotation at 1080p with the 3D layer and density shadow on
   the software renderer; keep fast quality while dragging and draw off the GUI thread (AP-20).
5. **Scope creep from 0.1.1**: RNAdo, co-editing, cluster profiles are large; held by D1, D7 and the
   class column of section 5.1.
6. **Windows GUI tests**: offscreen platform differences (fonts, DPI) and file dialogs; tests drive
   widgets by object name and never open native dialogs.
7. **Qt Test on GCC 11 / MSVC `/W4 /WX`**: Qt's test macros warn on some compilers; isolate in the test
   target's own flags if needed.
