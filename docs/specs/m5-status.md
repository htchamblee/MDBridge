# Milestone 5: status (2026-11-25)

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

Known issue for W3 (packaging): when the host lacks one of the AppImage's fallback libraries,
AppRun puts the whole `usr/lib/fallback` folder first on LD_LIBRARY_PATH, so Rocky's fontconfig,
freetype, harfbuzz and GL replace the host's; only the missing libraries should be supplied.

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

## Wave 3: in progress

Docs, accessibility (D10), Windows checks (native '\\' in every path shown in the interface;
%APPDATA%, D16; Start-menu entry), GUI smoke tests on bare Rocky and Windows, AppRun fallback fix.
