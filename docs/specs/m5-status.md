# Milestone 5: status (2026-11-25)

Pushed and green on CI: everything up to the pre-M5 fix wave (`9abefbf`) and the M5 plan with the
owner's decisions (`docs/specs/m5-plan.md`, all D1-D20 accepted, D4 amended: no separate 3D view).

Reports kept for reference (`docs/reports/`):

- `2026-11-25-mdatk-0.1.1-fix-comparison.md`: every MDBridge fix compared with MDATK 0.1.1 (the
  adopted items are in the fix wave).
- `2026-11-25-mdatk-bug-report.md`: 85 bugs still present in MDATK 0.1.1, handed to the MDATK side.

## Wave 0 (in progress, not pushed)

Started from `9abefbf`, each package on its own local branch and worktree under the session
scratchpad (`m5-w0a`, `m5-w0b`, `m5-w0d`, `m5-w0e`, `m5-w0f`). Work was stopped early for a usage
limit; each branch ends in a `WIP:` commit that says what is done and untested. The worktrees live
in the session container only: if it is reclaimed they are lost and the packages restart from the
plan (section 6).

| Package | Content | State |
|---|---|---|
| W0-A | GUI parity generators (`make_gui_cases.py`, `improvements_gui.py`, demo / clean / citations CLI cases, regenerate job with GTK) | WIP |
| W0-B | Core edit model (`ProjectDocument`, paths, undo, located problems, rename cascade, recovery) | WIP |
| W0-C | Field widgets | not started (waits for W0-B's API) |
| W0-D | Job runner and log view | WIP |
| W0-E | `demo`, `clean`, `citations` commands | WIP |
| W0-F | GUI test harness (`mdbridge_gui_tests`, CI) | WIP |

Next: finish and verify each package (the acceptance column of the plan), merge, regenerate,
full ctest with app + FFmpeg + Vulkan, cli_parity, push; then W0-C, then waves 1-3.
