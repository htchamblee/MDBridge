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

`python/` holds what Python and MDATK compute for the functions of `mdbridge/core/python.hpp`,
`json.hpp`, `hash.hpp` and `naming.hpp` (`tests/core/test_python.cpp`, `test_hash.cpp`,
`test_naming.cpp`). All of it is written by `tests/parity/python_cases.py`; the first six files are
the vectors of the milestone-2 porting specs, which the script rebuilds byte for byte:

| File | Content |
|---|---|
| `python_repr_cases.json` | `repr()`, `str()` and `type().__name__` of every TOML value type (strings, ints, floats, datetimes ...) |
| `stable_hash_cases.json` | `config.stable_hash` of 19 inputs: the exact JSON blob, its SHA-256, hashes of length 0 / 8 / 16 / 64 / 100; `file_fingerprint` of a file with a set mtime and of a missing file |
| `resolve_paths_cases.json` | `Project.resolve_paths` with `HOME=/home/tester` and base `/base/proj` (expanduser, join, normpath; root's home was `/root`) |
| `dirs_cases.json` | `out_dir` / `cache_dir` / `workflow_dir` (working directory `{cwd}`), `sanitize_project`, `group()` |
| `round2_vectors.json` | `repr(round(x, 2))` of 732 doubles given as `float.hex()` (step-log durations) |
| `naming_cases.json` | every function of `mdatk/workflows/naming.py` |
| `stable_hash_inputs.json` | the inputs of `stable_hash_cases.json` as data (tagged JSON: `{"$float": "nan"}`, `{"$datetime": iso}`, `{"$date": iso}`, `{"$int": "..."}` beyond 64 bits) |
| `posixpath_cases.json` | `normpath`, `isabs`, `dirname`, `basename`, `splitext`, `join`, `abspath` of awkward paths; `expanduser` with several `HOME` values |
| `string_repr_cases.json` | `repr()` and `json.dumps()` of the characters on both sides of every printable / non-printable boundary (Unicode 15.1) and of undecodable bytes (surrogateescape) |
| `json_loads_cases.json` | `json.loads` of valid, invalid and randomly mutated texts: `repr()` of the result or the error message (`lossy` marks results MDBridge cannot hold: integers beyond 64 bits, lone surrogates) |
| `round_cases.json` | `repr(round(x, n))` for n from -309 to 324, with the OverflowError cases |
| `sha256_cases.json` | `hashlib.sha256` of `"a" * n` for n = 0..200 (every padding case) and of longer texts |
| `naming_extra_cases.json` | `parse_name`, `sanitize_project` / `sanitize_tag` and `new_stamp` of further edge cases (Unicode digits, final newlines, greedy matching, years below 1000) |

To regenerate, with Python 3.13 (json's messages and the Unicode tables depend on the version) on
Linux and an MDATK checkout at `v0.0.4.1` (numpy importable; nothing is written into the checkout):

```bash
python3 -I -B tests/parity/python_cases.py /tmp/mdatk tests/data/python
```

Never put research structures or trajectories of real projects here.
