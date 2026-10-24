# Test data

`projects/` holds project files; `expected/` holds what MDATK v0.0.4.1 (the reference
implementation) reads from each of them. `tests/core/test_project_file.cpp` checks that MDBridge
reads every file in `projects/` the same way, byte for byte in the JSON form of the result.

| File | Made by |
|---|---|
| `default.mdatk.toml`, `full.mdatk.toml`, `split.mdatk.toml` + `split.log.toml`, `MDATK_params_log.toml` | MDATK v0.0.4.1, `tests/parity/make_projects.py` |
| `handwritten.mdatk.toml` (CRLF line endings, alias keys, no format marker), `datetimes.mdatk.toml`, `empty_params_log.toml`, `missing_log.mdatk.toml`, `other_format.toml`, `broken.mdatk.toml`, `preset.json` | by hand |

To regenerate (for example after adding a fixture), with a checkout of MDATK at `v0.0.4.1`:

```bash
git clone --branch v0.0.4.1 --depth 1 https://github.com/htchamblee/MDATK /tmp/mdatk
export PYTHONPATH=/tmp/mdatk          # MDATK's dependencies must be installed (pip install /tmp/mdatk)
python tests/parity/make_projects.py tests/data/projects
cd tests/data/projects && python ../../parity/reference.py --out ../expected *.toml *.json
```

`make_projects.py` time-stamps what it writes, so regenerated files differ from the committed ones
in their `saved` / `time` entries; commit the project files and their references together. CI's
parity job (`.github/workflows/ci.yml`) does the same with fresh files on every push.

`text/` holds the expected results of Python's text functions that project files depend on:
`textwrap_cases.toml` (textwrap.wrap and str.splitlines over string literals of MDATK's sources and
edge cases; `python tests/parity/textwrap_cases.py <MDATK checkout> tests/data/text/textwrap_cases.toml`)
and `tomlio_golden.toml` (a document written by MDATK's `tomlio.dumps`; the data and comments are
repeated in `tests/core/test_text.cpp`).

`project-file/` holds what MDATK writes when it saves projects (`tests/core/test_project_io.cpp`,
`test_step_log.cpp`, `test_toml_writer.cpp` compare MDBridge's output with it byte for byte):

| Folder | Content | Made by |
|---|---|---|
| `golden/` | a default and a rich project, its step log (every failure kind), a truncated output, the log moving to a side file (`split*`), a parameter log, the recent-projects list (`settings.toml`, work folder shown as `/home/user/projects`) and MDATK's messages | `python3 -I -B tests/parity/make_project_file_golden.py <MDATK checkout> tests/data/project-file/golden` (then delete `manifest.txt`, which lists sizes that depend on the temporary folder) |
| `vectors/` | `tomlio` / `textwrap` / `_entry_text` input-output pairs; tagged values explained at the top of the generator | `python3 -I -B tests/parity/make_writer_vectors.py <MDATK checkout> tests/data/project-file/vectors` (keep the six `*_vectors.json` files the tests use) |
| `workflow/` | projects with planner steps of every step type, edge cases and the projects of two MDATK tests, saved with `saved = 2026-10-15T14:32:00` | `golden_projects()` of the step-type generator in `tests/parity/` |

The generators fix the clocks (`saved`, entry times, durations), so their output is reproducible
and equal to the committed files.

Never put research structures or trajectories of real projects here.
