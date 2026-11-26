# Milestone 5: status (2026-11-26)

Milestone 5 is complete pending CI on the wave-3 merge: every wave is done and checked locally;
CI runs on the pushed branch decide.

Pushed and green on CI: everything up to the pre-M5 fix wave (`9abefbf`) and the M5 plan with the
owner's decisions (`docs/specs/m5-plan.md`, all D1-D20 accepted, D4 amended: no separate 3D view).

Reports kept for reference (`docs/reports/`):

- `2026-11-25-mdatk-0.1.1-fix-comparison.md`: every MDBridge fix compared with MDATK 0.1.1 (the
  adopted items are in the fix wave).
- `2026-11-25-mdatk-bug-report.md`: 85 bugs still present in MDATK 0.1.1, handed to the MDATK side.

## Wave 0: done (CI green at 183f698)

Edit model (`core::ProjectDocument`, problems), field widgets (`src/app/gui/forms/`), job runner
and log view (`src/app/gui/jobs/`), `demo` / `clean` / `citations` (with the process runner in
`src/core/process`, cleaning in `src/analysis/cleaning`), GUI parity generators
(`make_gui_cases.py`, 34 cases, skipped until their pages exist) and the GUI test harness
(`mdbridge_gui_tests`, `gui.forms`, `gui.jobs`). Owner decision: cpptraj scripts name MDBridge
(CMD-08).

Known issue for W3 (packaging), resolved in wave 3 (PKG-01): when the host lacked one of the
AppImage's fallback libraries, AppRun put the whole `usr/lib/fallback` folder first on
LD_LIBRARY_PATH, so the bundled fontconfig, freetype, X11 and GL replaced the host's.

## Wave 1: done (CI green at 70b8a42)

Window shell in MDATK's layout (no separate 3D view, D4), the Input & cleaning, Residue groups,
Analyses, Structure view and Outputs pages, problems panel and page badges; every page fits the
side bar (WN-11). GUI tests use the bundled font and report through stderr on Windows. 31 of the
42 recorded MDATK GUI cases run and pass; the other 11 need the canvas or the slot editors (W2-A).

Notes for wave 3: on Windows, show paths in interface text with native '\\' separators
consistently (title, recent list, change notice, log); storage unchanged.

Owner decisions since: fonts switch to IBM Plex Serif everywhere (interface, figures, charts,
videos) with IBM Plex Math as the only fallback for missing symbols (branch `plex-serif`).

## Wave 2: done (CI green at d1a4310)

Layout & panels page and canvas editor (the canvas's structure panel is the 3D view; Vulkan while
dragging with a real GPU; one undo step per drag and per wheel burst), the window's actions
(Render, still, Export PDBs, Re-render, Inspect, cleaning, Centre on model, result dialogs,
About / BibTeX with run citations). All 42 recorded MDATK GUI cases pass. Fonts: IBM Plex Serif
everywhere, Plex Math fallback, Plex Mono for fixed-width interface text and sequence text
(FONT-01..03); DejaVu removed. Windows: project saves retry a rename the system refuses for up
to 2 s (PRJ-13).

## Wave 3: done (local checks at the wave-3 merge; CI pending)

W3-A: on Windows every path the interface shows uses native '\\' separators (`gui::displayPath` /
`displayPathsIn`: title, recent list, change notice, log lines including the pipeline's, dialogs,
status bar, problems, field flags; storage and the command line unchanged, WN-12); keyboard and
screen-reader access (D10, WN-13: Tab reaches every control in visual order, every focusable widget
named, the canvas takes the focus with a focus frame, arrow keys rotate and + / - zoom the
structure, one undo step per burst); every page fits at 150 % and 200 % (`gui.window.scale150` /
`scale200`); settings locations checked (`%APPDATA%\\MDBridge`, D16, WN-14); the docs (README,
BUILDING, DIFFERENCES interface index, ARCHITECTURE status).

W3-B: the AppImage supplies only the fallback libraries the host lacks, the GL dispatch libraries
(glvnd) as one unit (PKG-01: the known issue above); `graphics --json` and the script report list
the libraries loaded; `mdbridge gui PROJECT --script FILE` drives the window from a JSON action
list (PKG-02; `gui.script`), used by `packaging/smoke_test.py --gui` to edit a field, save (exactly
one line changes), render a still and photograph every page on the AppImage (Ubuntu, bare Rocky 9
under Xvfb) and on the installed Windows program; CI checks the Windows Start-menu entry and starts
its target on a project (PKG-03; no file association, as MDATK).

Integration: the script runner's messages (its report's "problem", stderr) show paths as the
interface does (`gui::displayPath`), its report's path fields keep the stored form
(`ScriptTests::problemShowsPathsAsTheInterface`).
