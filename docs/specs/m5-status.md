# Milestone 5: status (2026-11-25)

Pushed and green on CI: everything up to the pre-M5 fix wave (`9abefbf`) and the M5 plan with the
owner's decisions (`docs/specs/m5-plan.md`, all D1-D20 accepted, D4 amended: no separate 3D view).

Reports kept for reference (`docs/reports/`):

- `2026-11-25-mdatk-0.1.1-fix-comparison.md`: every MDBridge fix compared with MDATK 0.1.1 (the
  adopted items are in the fix wave).
- `2026-11-25-mdatk-bug-report.md`: 85 bugs still present in MDATK 0.1.1, handed to the MDATK side.

## Wave 0 (started, then stopped for a usage limit before any code was written)

All six packages restart from the plan (section 6); nothing was lost. Notes the packages gathered:

- **W0-B edit model:** commit the documented header first. Lossless numbers need the original text
  of each loaded field kept until the user edits it (the typed structs hold doubles). The rename
  cascade ports MDATK 0.1.1's `Project.rename_group`.
- **W0-C field widgets:** start after W0-B's header.
- **W0-D jobs and log:** follow MDATK 0.1.1 (`gui/jobs.py`, `logpanel.py`, close handling in
  `app.py` / `rnado.py`): progress throttled to 0.5 % steps; the log keeps at most 6000 lines and
  trims to 5000, appends in batches every 40 ms of at most 2000 lines with a note counting skipped
  lines; on close ask, stop the job, wait up to 15 s.
- **W0-E commands:** byte-identical `demo` output depends on numpy's random numbers, scipy's
  nearest-neighbour search, numpy rounding and the NetCDF writer layout; check the existing
  `numpy_random` and parm7 / NetCDF writers first.
- **W0-A / W0-F:** W0-A defines the action-script JSON first (actions open, set, save, canvas
  paths) so W0-F can consume it. `MainWindow` has no save yet: the harness needs a save hook calling
  `core::saveProject`, with `core::setClock` for fixed time stamps.

Next: run the packages, verify each against its acceptance column, merge, regenerate, full ctest
with app + FFmpeg + Vulkan, cli_parity, push; then waves 1-3.
