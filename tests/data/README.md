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

`project-file/` holds what MDATK writes when it saves projects (`tests/core/test_project_io.cpp`,
`test_step_log.cpp`, `test_toml_writer.cpp` compare MDBridge's output with it byte for byte):

| Folder | Content | Made by |
|---|---|---|
| `golden/` | a default and a rich project, its step log (every failure kind), a truncated output, the log moving to a side file (`split*`), a parameter log, the recent-projects list (`settings.toml`, work folder shown as `/home/user/projects`) and MDATK's messages | `python3 -I -B tests/parity/make_project_file_golden.py <MDATK checkout> tests/data/project-file/golden` (then delete `manifest.txt`, which lists sizes that depend on the temporary folder) |
| `vectors/` | `tomlio` / `textwrap` / `_entry_text` input-output pairs; tagged values explained at the top of the generator | `python3 -I -B tests/parity/make_writer_vectors.py <MDATK checkout> tests/data/project-file/vectors` (keep the six `*_vectors.json` files the tests use) |
| `workflow/` | projects with planner steps of every step type, edge cases and the projects of two MDATK tests, saved with `saved = 2026-10-15T14:32:00` | `python3 -I -B tests/parity/make_workflow_goldens.py <MDATK checkout> tests/data/project-file/workflow` (MDATK's dependencies importable) |

The generators fix the clocks (`saved`, entry times, durations), so their output is reproducible
and equal to the committed files.

`project-model/` holds what MDATK's parameter model (`mdatk/config.py`) computes, for
`mdbridge/core/project.hpp` (`tests/core/test_project_model.cpp`, which also reuses
`python/resolve_paths_cases.json` and `python/dirs_cases.json`). Values are tagged as in
`project-file/vectors/`; results are compared as Python `repr()` texts.

| File | Content |
|---|---|
| `from_dict_cases.json` | `Project.from_dict(input, base_dir)`: hand-made inputs (aliases, unknown keys, error order, lists of objects replacing the defaults, free-form tables, path resolution), every field given a value of every TOML type (with None, NaN, infinity, 64-bit limits), integral floats in integer fields; `repr()` of the result or the error message, and the warnings logged |
| `validate_cases.json` | `validate()`, `density_problems()`, `density_active()`, `layout_warnings()`, `group_names()` of 103 projects (every check, NaN values, stored placements, fit frames), run in a folder holding the empty files listed in it |
| `config_doc.json` | `projectfile.config_doc` of every path of the model |
| `default_project.json` | `repr(Project().to_dict())` |

To regenerate, with Python 3.13 on Linux and an MDATK checkout at `v0.0.4.1` (numpy and scipy
importable; nothing is written into the checkout):

```bash
python3 -I -B tests/parity/make_model_cases.py /tmp/mdatk tests/data/project-model

`step-types/` holds the step types of MDATK's Project Planner and what MDATK computes with them
(`tests/core/test_step_types.cpp`). Values are tagged JSON (`{"t": "float", "hex": float.hex()}`,
`{"t": "str", "v": text}`, ...; see the top of the generator):

| File | Content |
|---|---|
| `step_types.json` | every registered step type (label, category, tag, description, ports, sections, options with their defaults' Python types, templates), `CATEGORIES` and `FILE_KINDS`; the same data as `src/core/step_types.inc` |
| `coerce_cases.json` | `Option.coerce` of about 220 inputs (numbers, `float()` spellings with whitespace, underscores, inf / nan, subnormals and overflow, Unicode digits, other types) for an option of every shape the step types use, and synthetic ones |
| `digit_cases.json` | `float()` of every Unicode decimal digit and of `5` between every whitespace character (an int option) |
| `normalise_cases.json` | `Tool.normalise` of defaults, valid and invalid values, unknown keys and malformed records per step type: result, issues, and the log lines with and without an issue list |
| `requirements_cases.json` | `schema.requirements_met` (truthiness, `==` across bool / int / float, tuple membership) |
| `setting_comments.json` | `projectfile._settings_doc` for every key of every step type, unknown keys and step types, and two synthetic step types that pin down its search order |
| `option_text_cases.json` | `projectfile._option_text` of synthetic labels, tooltips and choices (the `strip(": ")` edge cases) |
| `registry_cases.json` | `get_tool`'s error text, `register()`'s errors, `tools_by_category`, `Section.template`'s error, `guess_kind`, `kind_label` |

To regenerate, with Python 3.13 and an MDATK checkout at `v0.0.4.1` (nothing is written into it;
the build itself needs no Python):

```bash
python3 -I -B tools/gen_step_types.py /tmp/mdatk          # src/core/step_types.inc, step_types.json
python3 -I -B tests/parity/make_step_type_cases.py /tmp/mdatk tests/data/step-types
```

`topology/` holds the topology fixtures and what MDAnalysis 2.11.0 builds from them:
`fixture_chains.pdb`, `fixture_noelements.pdb`, `fixture_segids.pdb`, `fixture_chains_ba.pdb` and
`demo_dry.parm7` (MDATK's dry demo system) from
`python3 -I -B tests/parity/make_selection_fixtures.py <MDATK checkout> tests/data/topology`, and their
`*.topology.json` dumps from
`python3 -I tests/parity/dump_topologies.py tests/data/topology tests/data/topology/*.pdb tests/data/topology/demo_dry.parm7`
(both with MDAnalysis installed). `golden_model.json` holds MDATK's chain numbering and Reporter events
for them.

`selections/` holds what MDAnalysis' `select_atoms` and MDATK's selection layer compute, for
`src/core/atom_selection.cpp` (`tests/core/test_atom_selection.cpp`) and `src/core/selections.cpp`
(`tests/core/test_selections.cpp`). Selected atoms are `[first, last]` runs of indices; errors are
`str(exception)` (the text MDATK puts into its messages).

| File | Content |
|---|---|
| `atom_selection_cases.json` | per fixture of `topology/` (rebuilt from the dumps, first frame included): about 6000 selection strings each with the selected atoms or the error: every keyword with matching, glob and malformed values, resid ranges with insertion codes, numbers beyond 64 bits, Unicode digits and whitespace, `prop` / `same` on every property, macros, syntax errors, random boolean combinations and random token sequences. `unsupported` marks the cases MDBridge rejects on purpose (`smarts`, `cyzone` / `cylayer` with coordinates, `type_index`, bond / angle objects of topologies with bonds) |
| `atom_selection_distance_cases.json` | frames with random positions for `fixture_chains.pdb` and `demo_dry.parm7`: no box, orthorhombic, triclinic, skewed, zero and invalid boxes, atoms inside and far outside the cell; `around`, `sphzone`, `sphlayer`, `isolayer`, `point`, `prop`, `same x`, `cyzone` on each |
| `selections_cases.json` | every function of `mdatk/selections.py` (residue ranges, groups nested / cyclic / duplicated / 8 and 9 levels deep, `select` on every fixture, element guessing, labels, frame and time lists) and `resolve_export_selection` / `frames_from_specs` of `mdatk/io/pdb_export.py` with the warnings they log; `python_numbers` holds `int()` / `float()` of awkward strings |

To regenerate, with MDAnalysis 2.11.0 (and, for the second, MDATK's dependencies) installed and an
MDATK checkout at `v0.0.4.1` (nothing is written into it):

```bash
python3 -I -B tests/parity/make_atom_selection_cases.py /tmp/mdatk tests/data/selections
python3 -I -B tests/parity/make_selections_cases.py /tmp/mdatk tests/data/selections
python3 -I tools/gen_selection_tables.py     # src/core/selection_tables.inc (macro name lists)
```

`io/parm7/` holds what MDATK's parm7 module (`mdatk/io/parm7.py`) and MDAnalysis 2.11.0's
TOPParser make of parm7 files, for `src/io/parm7.cpp` (`tests/io/test_parm7.cpp`). Floats are
`float.hex()` strings; errors are `{"type", "message"}` with the file's path written `{path}`.

| File | Content |
|---|---|
| `read_cases.json` | `read_sections` of about 75 crafted files (bytes as hex, with the `wanted` set): `%FORMAT` spellings (spaces, lower case, Unicode digits, two parentheses, width 0, bad lines), D exponents, blank fields, long and short lines, repeated flags, `%VERSION` / `%COMMENT` lines, CRLF / CR / mixed line ends, Unicode and invalid UTF-8, bad and overflowing numbers, random sections; and of `topology/demo_dry.parm7` |
| `writer_cases.json` | the bytes `write_parm7` writes for sections of every kind (or its error and what it wrote before it: non-ASCII text, NaN / infinity as integers, values of another type), and `_format_values` of 300 random doubles and edge values |
| `lj_cases.json` | `LJTable.from_parm7(demo_dry.parm7)`, `pair_coefficients` of about 1800 atom pairs (negative indices included), the types and `mdatk.demo.LJ_TYPES` MDATK's own LJ test uses, and the "missing parm7 sections" errors |
| `topology/*.parm7`, `*.topology.json` | variants of a small 40-atom system written by `write_parm7` (line ends, `RESIDUE_CHAINID` sections that are used or skipped, atomic numbers 0 / -1 / missing, other formats, `%COMMENT` lines, unknown sections, empty bond sections, residue pointers) with the dump of `tests/parity/dump_topologies.py` |
| `topology_errors.json` | the variants TOPParser (or the Universe's length check) rejects, with `str(exception)` |
| `is_parm7_cases.json` | the extension test MDATK uses for parm7 topologies |

`io/external/` holds what MDATK's external-series reader (`mdatk/io/external.py`) and numpy
compute, for `src/io/external.cpp` (`tests/io/test_external.cpp`). Floats are `repr()` strings;
errors are `str(exception)` with the folder written `{dir}`.

| File | Content |
|---|---|
| `table_cases.json` | `read_table` of about 60 files (bytes as hex): CSV, tabs, whitespace and per-line detection, explicit delimiters, comment lines and the header taken from them, xvg files, ragged rows, Unicode whitespace and digits, CR / CRLF, invalid UTF-8, errors |
| `column_cases.json` | `_column_index` (names before indices, Unicode digits, errors) and `_sanitize` |
| `interp_cases.json` | `numpy.interp(x, xp, fp, left, right)` of 2000 random small inputs (duplicates, unsorted points, NaN and infinities everywhere) and 16 long ones; its errors |
| `load_cases.json` | `load_external` of files and specs on datasets: array names and values and `repr(meta)`, or the error |

`dataset/dataset_cases.json` holds what MDATK's data model (`mdatk/data.py`) computes, for
`src/core/dataset.cpp` (`tests/core/test_dataset.cpp`): results built with `add_series` /
`add_matrix` / `add_pairmap` (their `repr(meta)`), datasets of them and about 65 queries each
(`time`, `time_label`, `keys`, `has`, `kind_of`, `series`, `matrix`, `pairmap` with `dense` and
`frequency`, `describe`) with their values or KeyError / ValueError messages, `CSRBuilder` runs
and the TypeErrors of `add_*`.

To regenerate, with Python 3.13, numpy and (for parm7) MDAnalysis 2.11.0 installed and an MDATK
checkout at `v0.0.4.1` (nothing is written into it):

```bash
python3 -I -B tests/parity/make_parm7_cases.py /tmp/mdatk tests/data/io/parm7
python3 -I -B tests/parity/make_external_cases.py /tmp/mdatk tests/data/io/external
python3 -I -B tests/parity/make_dataset_cases.py /tmp/mdatk tests/data/dataset
```

`io/trajectories/` holds coordinate files and what MDAnalysis 2.11.0 reads from them, for
`src/io/netcdf.cpp`, `amber_ascii.cpp`, `dcd.cpp`, `gromacs.cpp` and `trajectory.cpp`
(`tests/io/test_trajectory.cpp`, `test_netcdf.cpp`). Times and dt are `float.hex()` strings, boxes
and positions big-endian float32 bits; errors are `{"type", "message"}` with the files' paths
written `{path0}`, `{path1}` ...

| File | Content |
|---|---|
| `files/*` | small synthetic trajectories (7 to 90 atoms): NetCDF from MDAnalysis' NCDFWriter and crafted with scipy (CDF-1 / CDF-2, every netCDF type, scale factors in float32 / float64, no time, no box, fixed frame dimension, broken Amber conventions and units, corrupt files), Amber mdcrd (box lines, 1 / 4 / 10 atoms, CRLF, long titles, bad fields) and rst7 / inpcrd (with and without time), DCD (MDAnalysis' writer and crafted: both byte orders, X-PLOR, unit cells as cosines / degrees / box vectors, fixed atoms, 4D, truncation, broken headers), XTC (compressed and uncompressed, precisions 10^2 to 10^5, integer sizes beyond 2^24 and beyond 64 bits) and TRR (single and double precision), truncated and broken XDR files; `demo_dry.xtc` / `.trr` / `.dcd`: frames of the demo trajectory |
| `demo_dry.nc` | `mdatk.demo.build_demo(n_frames=24, solvent=False)`'s trajectory (topology: `topology/demo_dry.parm7`) |
| `cases.json` | `cases`: per case (one per line) the files (several: MDAnalysis' ChainReader), `format`, `dt` and atom count given to the Universe, and `n_frames`, `dt`, each frame's `time`, `box`, positions (`pos`, or `pos_sha256` of the little-endian float32 bytes for the demo), `dt` where it changes along a chain, the TRJ `header`, frame errors; or the error. `writer`: frames of a source file and the file MDAnalysis' NCDFWriter writes from them |
| `writer/*.nc` | those NCDFWriter files (`demo_dry.nc` rewritten is `demo_dry.nc` itself) |
| `formats.json` | MDATK's extension table (`mdatk/io/loader.py`), MDAnalysis' per-file format guess and its "Unknown coordinate trajectory format" errors |

To regenerate (MDAnalysis 2.11.0 with scipy and without netCDF4; the output is reproducible):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_trajectory_cases.py /tmp/mdatk tests/data/io/trajectories
```

Never put research structures or trajectories of real projects here.
