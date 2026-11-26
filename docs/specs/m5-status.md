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

## Wave 1: in progress

W1-A window shell, W1-B pages 1 (Input & cleaning, Residue groups, Analyses), W1-C pages 2
(Structure view, Outputs), W1-D problems. Then wave 2 (layout page + canvas editor, actions) and
wave 3 (docs, accessibility, Windows checks, packaging fix above).
