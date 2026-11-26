# Test data

## MDAnalysis with MDBridge's independent tables

Wherever expected data comes from MDAnalysis 2.11.0, it is MDAnalysis 2.11.0's behaviour with
MDBridge's independent residue-name, atom-name, element and mass tables patched in, not
MDAnalysis' own tables. The owner decided that MDBridge carries no data tables extracted from
MDAnalysis: its tables are rebuilt from independent, attributable public sources (the wwPDB
Chemical Component Dictionary, the PDB format v3.3 and IUPAC-IUB residue names, the CIAAW 2024
abridged standard atomic weights and NUBASE2020 most stable isotopes, the Amber 2024 and GROMACS
manuals, OpenMM 8.1.2's force-field files; see `data/sources/ATTRIBUTION.md`).
`tools/gen_data_tables.py` writes them to
`src/core/include/mdbridge/core/detail/data_tables_independent.inc` from the extracts committed
under `data/sources/`. MDBridge's results therefore differ from stock MDAnalysis
wherever these tables matter, and the expected data describes that behaviour: MDAnalysis'
algorithms, unchanged, reading MDBridge's tables.

`tests/parity/independent_tables.py` loads the tables with `tools/gen_data_tables.py`'s own
functions (from `data/sources/`, not from the `.inc`) and replaces MDAnalysis' data in place
before a generator uses it: the residue sets of the `protein`, `nucleic` and `water` macros and
the atom-name sets of `backbone`, `nucleicbackbone`, `nucleicbase` and `nucleicsugar`; the
guessers' element symbols (`Z2SYMB` / `SYMB2Z`), masses (every element, under its symbol and its
upper-case spelling), special atom-name -> element table (the force-field name -> element table,
which is by name alone like MDAnalysis' lookup; virtual sites map to no element, mass 0) and
guessable one- and two-letter elements; PDBWriter's ion names and Hg / Cl residue
lists, and its include / exclude pairs (empty: no independent counterpart). Only data is
replaced, no MDAnalysis function. The generators that depend on these tables apply the patch:
`dump_topologies.py` (every dump, also when imported), `make_atom_selection_cases.py`,
`make_selections_cases.py`, `make_parm7_cases.py`, `make_loader_cases.py`, `make_export_cases.py`,
`make_rmsd_cases.py`, `make_universe_cases.py`, `make_analysis_cases.py`,
`make_analysis_framework_cases.py`, `make_pipeline_cases.py`, `make_secstruct_cases.py`, the
render generators (`render_common.py`: `make_figure_cases.py`, `make_chart_cases.py`,
`make_molecule_cases.py`, `make_render_cli_cases.py`) and `cli_parity.py` (its MDATK side); the
others do not touch them. `tools/gen_selection_tables.py` and
`tools/gen_pdb_tables.py` no longer write any of these tables, only the property names of the
selection language and Python's Unicode data.

MDAnalysis' algorithms read atom names alone, so the tables are filtered by a collision rule in
`tools/gen_data_tables.py` (documented there): a force-field name -> element entry, a guessable
two-letter element or a PDB-writer ion name is dropped when a CCD standard-residue atom name (or
its PDB v2.3 spelling, or an atom name of the CCD peptide-linking components HIP, HSE, HYP, LYN of
the protein list) of another element would reach it (CA, CD, CE, HG, ND names; Ca and Cd as
guessable elements; CA, CD, CE, HG, HO, ND ion names). The generator fails if a collision is
left, and `/opt/mdbridge-venv/bin/python -I -B tests/parity/independent_tables.py` guesses the
element of every such reference name, and the column the PDB writer gives it, with the patched
MDAnalysis and lists any that differ from the CCD (none at present). Ions whose atom name is one of
the dropped names (a calcium ion named `CA`, cadmium `CD`, mercury `HG`) are therefore guessed as
carbon / hydrogen and written from column 14 by the name-only rules; MDBridge's residue rule
below handles them.

**MDBridge's residue-scoped ion rule (a deliberate improvement over MDAnalysis).** The generated
table `kIndResidueAtomElements` lists (residue name, atom name, element) for atoms whose element
the residue decides: every CCD one-atom component (`CA`/`CA` -> Ca, `CD`/`CD` -> Cd, `HG`/`HG` ->
Hg, `ZN`/`ZN` -> Zn, `NA`/`NA` -> Na ...), the `HG` / `CL` atoms of the CCD components in which
they are mercury / chlorine (`EMC`/`HG` -> Hg ...), and the one-atom ion residues of the OpenMM
Amber and CHARMM36 ion libraries (`SOD`/`SOD` -> Na, `CAL`/`CAL` -> Ca, `CLA`/`CLA` -> Cl,
`POT`/`POT` -> K, `CD2`/`CD` -> Cd, `ZN2`/`ZN` -> Zn ...); upper case, without the residue names of
the protein, nucleic and water lists. When an atom's type is guessed from its name and its
(residue name, atom name), compared in upper case, is in the table, the type is that element
symbol in upper case and, where the mass follows the type (no elements column, no masses in the
file), the mass is the element's; otherwise the name-only procedure applies unchanged. The PDB
writer starts the name of such an atom in column 13 when the element has two letters; otherwise
its rules apply unchanged. MDAnalysis cannot do this, so `tests/parity/independent_tables.py`
post-processes patched MDAnalysis' results with exactly these semantics:
`residue_rule_universe` (guessed types and masses, applied in every dump of `dump_topologies.py`)
and `residue_rule_pdb_text` (the name column of the files the export generator records).
`io/loader/files/ions.pdb` and `ions.gro` exercise it: a protein fragment (CA, CD1, HG1, CE stay C
/ H), lone ions named the CCD way and the force-field way, `CD2`/`CD`, `EMC`/`HG`, and `K+`/`K+`,
which no source names (name rule: K).

Where MDBridge deliberately improves on MDATK (docs/DIFFERENCES.md), the generators apply the
improvement to MDATK's results: `tests/parity/improvements_cli.py` (project validation, parameter
and command-line messages, external data, PDB export, map errors; items PRJ-* and CLI-*) is used by
cli_parity.py and the model, loader, export, maps, analysis and pipeline generators.

Reader improvements (docs/DIFFERENCES.md, "Amber, DCD and GROMACS readers"): the expected data of
the parm7, NetCDF, mdcrd / rst7, DCD and XTC / TRR readers is MDBridge's improved behaviour.
`tests/parity/improvements_readers_amber.py` (documented per item id) adjusts what the generators
record from MDAnalysis: `apply()` (called by `independent_tables.apply()`) gives every parm7
Universe float64 charges and MDATK's analyses MDBridge's versions
(`tests/parity/analysis_versions.py`, checked against src/analysis/*.cpp); make_parm7_cases.py,
make_trajectory_cases.py, make_loader_cases.py and make_export_cases.py call its adjust functions.

**Selection improvements (docs/DIFFERENCES.md "Selections").** `tests/parity/improvements_selections.py`
turns MDAnalysis' and MDATK's selection answers into MDBridge's: `make_atom_selection_cases.py`
passes every case through its model of MDBridge's grammar (checked against MDAnalysis wherever
MDBridge keeps its answer), and `make_selections_cases.py`, `make_analysis_framework_cases.py` and
`cli_parity.py` replace MDATK's group expansion, `select()`, `atom_elements` and the cache keys'
group lookup with MDBridge's behaviour.

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
parity job (`.github/workflows/ci.yml`) does the same with fresh files on every push. With
`--fixed-clock` it uses the clock of the run that wrote the committed files and writes them byte
for byte.

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
| `mdatk-0.1.1/` | project files a later MDATK writes (D9; DIFFERENCES PRJ-09, PRJ-10): a plan with every key 0.1.x added (`workflow.parallel` ..., per step `adopted`, `folder_tag`), step types and options 0.0.4.1 lacks, the log marker in the notes and a script, a step log with `waiting`, `stopped` and `workflow supersede`; the same project with its log in a side file; `read.json`, what MDATK 0.1.1's `read()` / `read_log()` return for them | `python3 -I -B tests/parity/make_mdatk_latest_projects.py <MDATK 0.1.1 checkout> tests/data/project-file/mdatk-0.1.1` (MDATK at the commit `regenerate_all.py` names, `--mdatk-latest`; standard library only) |

The generators fix the clocks (`saved`, entry times, durations), so their output is reproducible
and equal to the committed files.

`edit-model/edit_model_cases.json` holds what MDATK 0.1.1 does where the GUI's edit model
(`mdbridge/core/project_document.hpp`, `problems.hpp`; `tests/core/test_project_document.cpp`) follows
it: `Project.rename_group` and `group_references` on projects that refer to their groups in every way
(changed locations, errors, `to_dict()` afterwards), `problems()` of every input of
`project-model/validate_cases.json` (field, message, location; 0.1.1's selection checks left out, D7) and
`problems.describe()` of fields with and without a project. Made by
`python3 -I -B tests/parity/make_edit_model_cases.py <MDATK 0.1.1 checkout> tests/data/project-model/validate_cases.json tests/data/edit-model`
(`--mdatk-latest`; numpy and scipy importable).

`project-model/` holds what MDATK's parameter model (`mdatk/config.py`) computes, for
`mdbridge/core/project.hpp` (`tests/core/test_project_model.cpp`, which also reuses
`python/resolve_paths_cases.json` and `python/dirs_cases.json`). Values are tagged as in
`project-file/vectors/`; results are compared as Python `repr()` texts.

| File | Content |
|---|---|
| `from_dict_cases.json` | `Project.from_dict(input, base_dir)` with MDBridge's unknown keys kept (PRJ-09, `tests/parity/improvements_core_io.py`): hand-made inputs (aliases, unknown keys, error order, lists of objects replacing the defaults, free-form tables, path resolution), every field given a value of every TOML type (with None, NaN, infinity, 64-bit limits), integral floats in integer fields; `repr()` of the result or the error message, and the warnings logged |
| `validate_cases.json` | `validate()`, `density_problems()`, `density_active()`, `layout_warnings()`, `group_names()` of 103 projects (every check, NaN values, stored placements, fit frames), run in a folder holding the empty files listed in it |
| `config_doc.json` | `projectfile.config_doc` of every path of the model |
| `default_project.json` | `repr(Project().to_dict())` |

To regenerate, with Python 3.13 on Linux and an MDATK checkout at `v0.0.4.1` (numpy and scipy
importable; nothing is written into the checkout):

```bash
python3 -I -B tests/parity/make_model_cases.py /tmp/mdatk tests/data/project-model
```

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
`python3 -I -B tests/parity/dump_topologies.py tests/data/topology tests/data/topology/*.pdb tests/data/topology/demo_dry.parm7`
(both with MDAnalysis installed; the dumps use MDBridge's independent tables, see above). `golden_model.json` holds MDATK's chain numbering and Reporter events
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
python3 -I tools/gen_selection_tables.py     # src/core/selection_tables.inc (prop / same property names)
python3 -I tools/gen_data_tables.py          # src/core/include/mdbridge/core/detail/data_tables_independent.inc (macro names, elements, masses)
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

`io/pdb/pdb_cases.json` holds what MDATK's PDB tools compute, for `src/io/pdb.cpp`
(`tests/io/test_pdb.cpp`): `read` — `pdbtools.read_pdb` of about 270 texts (MDATK's test fixtures,
ChimeraX damage, models, ANISOU / SIGATM / SIGUIJ / orphan ANISOU, MASTER, records after END, CONECT
variants, serial and residue-number wraps, lower-case records, Unicode, every line break of
`str.splitlines`, broken fields and random mutations of a normal and a damaged file): the structure
(residues as `[name, chain, resseq, icode, label, coords(), atoms]`, atoms as `[record, name, alt,
tail, element, charge, serial, anisou]`) or the error, and `write_pdb` of it with CONECT stripped,
kept (`write_conect`, when it differs) and with SEQRES records of `seqres_chains`; `files` —
`read_pdb` of files (bytes as hex: line ends, BOM, invalid UTF-8); `residue_names` —
`parse_residue_name` of every name of up to three characters over an alphabet; `insert` —
`annotate._insert_records`; `transform` — `builder.transform_pdb_text` (rotations, translations,
ANISOU, line ends, limits); `tail` — `builder._tail`. Floats are `float.hex()` strings. MDATK runs
with MDBridge's chem improvements (`improvements_chem.py`): CH-01 makes a residue number below -999
a write error (`write_error`) and names it when a line MDATK wrote that way is read; CH-07 gives
`seqres_chains` the three-character SEQRES names of `chem::seqresNames`.

`io/pdb/export/` holds the PDB files MDATK's frame export writes (`mdatk.io.pdb_export.export_frames`
through MDAnalysis 2.11.0's PDBWriter), for `src/io/pdb_export.cpp` (`tests/io/test_pdb_export.cpp`):
`export_cases.json` lists per case the system (topology, trajectories, dt), the selection, frames,
prefix and stamp, and the files written (text, or SHA-256 and first lines when longer than 16 KiB),
the progress calls, or the error and the files left behind. The systems are the dry demo
(`topology/demo_dry.parm7` with `io/trajectories/demo_dry.nc` and its DCD / XTC copies), the
7-atom files of `io/loader/files` with the mdcrd / rst7 / XTC / TRR / DCD / NetCDF files of
`io/trajectories/files` (HEADER from mdcrd titles, cells and no cells), PDB and GRO files exported
from their own coordinates (HEADER / COMPND / REMARK, models, formal charges, occupancies, insertion
codes, segids, CONECT), the parm7 chain-ID variants of `io/parm7/topology` with `files/forty.nc`, and
the inputs of `files/`: chain IDs MDAnalysis keeps or replaces by `X`, atom names of every alignment
rule, coordinates at and beyond the PDB limits (`limit*`, `edge*.gro`; a failing frame removes its
file). `big` holds the SHA-256 of three exports of a 100001-atom system that the generator and the
C++ test build from the same formulas (serials beyond 99999, resids beyond 9999 and below -999, no
CONECT records once a bonded serial reaches 100000).

`io/loader/` holds what MDATK's loader (`mdatk/io/loader.py`) and MDAnalysis 2.11.0 make of
topology files and file combinations, for `src/io/loader.cpp` and `src/io/pdb_topology.cpp`
(`tests/io/test_loader.cpp`, which also reads the PDB fixtures of `topology/` and `demo_dry.parm7`):

| File | Content |
|---|---|
| `files/*.pdb`, `*.gro` | crafted PDB / GRO files: elements (valid, lower case, invalid, blank, none), types and masses MDAnalysis guesses, segids and chain IDs, formal charges (valid, unknown entries), resid wraps, hybrid-36 and `*****` serials, CONECT variants, models with one, every or some CRYST1 records, placeholder / broken / zero cells, HEADER / COMPND / REMARK, CRLF, Unicode and invalid UTF-8, GRO boxes (3 / 9 values, zero, packed, broken), GRO resid wraps, precision and broken lines; `ions.pdb` / `ions.gro` for the residue-scoped ion rule (see the first section) |
| `files/*.topology.json` | the topology MDAnalysis builds from each readable file (`tests/parity/dump_topologies.py`) |
| `topology_cases.json` | per file: its dump or the Universe's error, and the PDB / GRO reader (frame count, dt, header / compound / remarks, every frame's float32 positions, cell and time, in order and shuffled, with and without dt) |
| `loader_cases.json` | `load_universe` + `system_summary` over a matrix of the demo and crafted files (formats by extension, chains, dt 0 / NaN / given, missing files, reader errors) with the residue-name counts; `frame_indices` of 1620 slices; the JSON `mdatk info --top ... --traj ...` prints |

`system_summary` lists the 12 most frequent residue names in the order of numpy's `argsort`, which
for equal counts depends on the CPU's sort kernel; MDBridge orders equal counts by name, and the
tests compare with that order (and with MDATK's count sequence). Self-bonds (`CONECT 1 1`), which
MDAnalysis keeps, cannot be held by `core::Topology` and are dropped.

To regenerate, with Python 3.13, MDAnalysis 2.11.0 (with scipy) in `/opt/mdbridge-venv` and an MDATK
checkout at `v0.0.4.1` (nothing is written into it; the output is reproducible; the export cases
use the loader's files, so run the loader generator first):

```bash
python3 -I -B tests/parity/make_pdb_cases.py /tmp/mdatk tests/data/io/pdb         # numpy and pytest
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_loader_cases.py /tmp/mdatk tests/data/io/loader
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_export_cases.py /tmp/mdatk tests/data/io/pdb/export
python3 -I tools/gen_pdb_tables.py                       # src/io/pdb_tables.inc (Python 3.13)
```

MDBridge's improvements of the PDB / GRO readers and the loader (docs/DIFFERENCES.md) are applied
to the loader, export, pipeline and command-line expectations by `tests/parity/improvements_readers_pdb.py`
(MDAnalysis runs on equivalent rewrites of the files; MDBridge's own error and warning texts).

`src/io/pdb_tables.inc` holds data the PDB code needs: Python's `str.upper()` mapping and the
characters `str.isalpha()` excludes from `\w` (Unicode 15.1). The element, mass and atom-name
tables of the guessers, PDBParser and the PDB writer are MDBridge's independent tables
(`src/core/include/mdbridge/core/detail/data_tables_independent.inc`, see the first section).

`io/maps/` holds density maps and what MDATK's map code (`mdatk/io/mrc.py`, the map functions of
`mdatk/density.py`, `mdatk density info`) gives for them, for `src/io/mrc.cpp`
(`tests/io/test_mrc.cpp`). Floats are JSON numbers, float32 statistics their bit patterns, map data
the SHA-256 of its little-endian float32 bytes ([z][y][x]); file paths in messages are `{path}`.

| File | Content |
|---|---|
| `files/*` | small crafted MRC / CCP4 files (mostly 5 x 4 x 3 voxels): every axis order in both byte orders, modes 0 (signed, unsigned without the `MAP ` tag, IMOD stamps and flags), 1, 2, 6 and 12 with extreme values, machine stamps and byte-order guessing, NSTART / ORIGIN / MXYZ / cell / angle rules, extended headers, every MapError; gzip files (members, empty members, FEXTRA / FNAME / FCOMMENT / FHCRC, zero padding, stored blocks, trailing data, truncation, corrupt blocks, bad CRC / size / magic / method) and a 40 x 32 x 24 map in one and in many members |
| `read_cases.json` | `read_map` of each file: grid, origin, voxel, data hash, `stats()`, first values; or the error. Where Python raises EOFError, zlib.error, BadGzipFile, MemoryError or OverflowError, `error` is MDBridge's documented MapError and `python` what Python raised |
| `info_cases.json` | what `mdatk density info` prints for some of the files |
| `stats_cases.json` | `stats()` and `map_level` of generated float32 arrays of 1 to 2^24 + 1 values (numpy's pairwise float32 summation, around its block sizes) |
| `write_cases.json` | `write_map` of generated maps (origins on and off whole voxels, beyond the int32 range, infinite): SHA-256 of the file, its first 224 header bytes, what `read_map` reads back; the error for an empty map |
| `simulate_cases.json` | `atom_weights` of element spellings; `simulate_map` of atom sets (MDATK's test geometry, random models, spacings, sigma factors, padding): grid, origin, voxel, data hash, statistics; its FitErrors |

The arrays and maps the C++ test builds itself come from the formulas in the generator's docstring.
numpy on AVX-512 CPUs computes `exp` with Intel's SVML, which MDBridge does not reproduce: the
generator runs `simulate_map` with the C library's `exp` (numpy's scalar path, what C++ `std::exp`
gives with glibc) and records the hashes numpy's own `exp` gives next to them (`numpy_sha256`,
identical for every committed case). Off glibc the test compares only the grid of simulated maps.

`rigid/rigid_cases.json` holds what MDATK's `RigidTransform` and rotation helpers compute
(`mdatk/density.py`), for `src/core/rigid_transform.cpp` (`tests/core/test_rigid_transform.cpp`):
`from_list` of rotations rounded to 1 to 12 digits, perturbed, scaled, reflected or singular, with
huge, tiny, non-finite entries and wrong lengths (the re-orthonormalised rotation or the message),
`to_list` at several digits, `apply` / `inverse` / `compose` / `about`, `rotvec_matrix`,
`quaternion_matrix`, `rotation_angle` and `anchor_points`. numpy multiplies with BLAS and takes the
SVD from LAPACK, whose polar factor is up to 3e-15 off the exact one for coarse inputs (MDBridge's is
within an ulp), so the test compares those results within a few 1e-15; `to_list`, `anchor_points`
and the messages exactly.

To regenerate, with Python 3.13, numpy 2.5.5 and scipy in `/opt/mdbridge-venv` and an MDATK checkout
at `v0.0.4.1` (nothing is written into it; the output is reproducible):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_map_cases.py /tmp/mdatk tests/data/io/maps
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_rigid_cases.py /tmp/mdatk tests/data/rigid
```

`citations/citation_cases.json` holds what MDATK's citation module (`mdatk/citations.py`)
computes, for `src/core/citations.cpp` (`tests/core/test_citations.cpp`): `REFERENCES` and
`RESOURCES` with every field (the same data as `src/core/citation_data.inc`), `unique`,
`keys_for`, `references`, `bibtex`, `text_lines` and `grouped` of about 160 id and key lists
(single entries, all, duplicates, unknown and empty entries, random subsets), `comment_lines` of
about 1200 resource lists, prefixes (non-ASCII included), widths (down to textwrap's "invalid
width" errors) and headings (`None`, empty, long), `pipeline_resources` of about 110 projects made
with `Project.from_dict` (cleaning, video, renderers, base-pair engines, secondary panels with
every kind of `source` option, RMSD targets with and without fit selections) and
`project_resources` with stand-in plans (`generating` is what MDBridge's caller passes: the
resources of the steps whose tool generates a script, in plan order). MDATK runs with MDBridge's
citation changes (`tests/parity/improvements_citations.py`, CIT-01 / CIT-02 in
docs/DIFFERENCES.md), which `tools/gen_citation_data.py` applies to the table too; `cli_parity.py`
and the render generators apply them before MDATK writes its params and step logs.

To regenerate, with Python 3.13 and an MDATK checkout at `v0.0.4.1` (nothing is written into it;
the output is reproducible and the build needs no Python):

```bash
python3 -I -B tools/gen_citation_data.py /tmp/mdatk        # src/core/citation_data.inc
python3 -I -B tests/parity/make_citation_cases.py /tmp/mdatk tests/data/citations
```

`io/cif/cif_cases.json` holds what MDATK's mmCIF code (the mmCIF part of `mdatk/dictionaries.py`)
gives, for `src/io/cif.cpp` (`tests/io/test_cif.cpp`). Paths in results and messages are written
`{root}`.

| Key | Content |
|---|---|
| `parse` | about 140 texts: MDATK's parser fixture, the synthetic acedrg / eLBOW dictionaries of its `tests/test_dictionaries.py` and what `finish_dictionary` makes of them, crafted texts (quoting rules, text fields, loops, pairs, several blocks, comments, keywords in any case, `global_` sections, line breaks of `str.splitlines`, Unicode whitespace and tags, invalid UTF-8 as `hex`, every error) and random token soups. For each: `list(_cif_tokens)` or its error, `parse_cif` or its error, `write_cif` of the result and `parse_cif` of that, and `block` / `get` / `col` / `column` lookups in several cases |
| `write` | `write_cif` of documents built directly (every quoting branch, multi-row and ragged pair items, ragged loops, line breaks in cells, empty items, non-ASCII widths) and `parse_cif` of the text |
| `cif_value` | `cif_value` of about 50 values |
| `lower` | `str.lower()` of strings around final sigma, U+0130, the Kelvin sign and other cased / case-ignorable characters, and SHA-256 hashes of `lower()` over every code point alone and in sigma contexts (`core::pythonLower`, which the case-insensitive lookups use) |
| `dictionary` | `is_dictionary_file` of files given as byte parts (hex, repeat count) around the 2,000,000-character limit (CRLF, CR, 4-byte characters, invalid UTF-8), of a folder and of a missing file |
| `find` | `find_dictionary` scenarios (MDATK's test plus code / `tmp` / age / path ordering, hidden names, glob patterns and sets, folders named `*.cif`, names with glob characters, relative and absolute sources and excludes); files carry their age in seconds, which the test recreates |

`src/core/case_tables.inc` (the Unicode data of `core::pythonLower`) is written by
`python3 -I tools/gen_case_tables.py` from the running Python's `str.lower()`.

To regenerate, with Python 3.13 and pytest installed and an MDATK checkout at `v0.0.4.1` (nothing
is written into it; the output is reproducible):

```bash
python3 -I tools/gen_case_tables.py                       # src/core/case_tables.inc
python3 -I -B tests/parity/make_cif_cases.py /tmp/mdatk tests/data/io/cif
```

`cli/cli_cases.json` holds what MDATK's command line (`mdatk/cli.py`) prints and writes for the
commands MDBridge has, for `src/app/commands/` (`tests/cli/test_commands.cpp`, which runs the same
command lines in-process): `init` / `init-params` (defaults, every common option, `--set` with list
indices, JSON values, free-form tables, Project properties and every error, projects and a parameter
log as `--project`, existing files, JSON presets), `info` (the demo system with NetCDF / XTC / DCD
trajectories, PDB and GRO topologies, residue groups that fail in every way, missing files), `list`
(without conda, and with a fake conda installation: environments with and without an executable
cpptraj, roots from `CONDA_EXE` / `CONDA_PREFIX` / `MAMBA_ROOT_PREFIX`, `~/.conda/environments.txt`),
`analyze` (the demo system with the default project, `--only` with instances, duplicates, spaces,
internal analyses and an empty list, cache reuse after an earlier run, partial reuse, a changed
parameter, `--force`, frame selections, output and cache folders, `use_pbc`, layout warnings and the
panels' and charts' analyses, no analyses at all, external data, a project file, every error),
`export-pdb` (frame lists, times, groups and selections, time stamps, `--quiet` / `--verbose`,
layout warnings, a project's `[output.pdb]`, every error), `density info` of every map of
`io/maps/files`, and `render`, `charts` and `preview` (`render: true`; the demo project of
`make_render_cli_cases.py`: full runs, `--from-cache`, strides, codecs, odd sizes, `ffmpeg_path`, PDBs,
a parameter log as project, every format, unknown formats and invalid settings, frames out of range),
`demo` (dry, solvated and default systems byte for byte, the project's cleaning with and without
cpptraj on PATH, impossible arguments, `--run` as a render case) and `clean` and cleaning in
`analyze` / `export-pdb` with a fake cpptraj (a shell script that copies the inputs to the outputs)
run directly, from `--cpptraj`, `$AMBERHOME` or a fake `conda run` (`--conda-env`, detected), reuse,
`--force`, changed settings, `--verbose`, every error; the script's time stamp, the marker's time and
the cleaning key are normalised. The cases run with `PATH=/usr/bin:/bin` and must find no conda,
mamba or micromamba there (cli_parity.py refuses to run otherwise; CI removes the link GitHub's
Ubuntu runners put in /usr/bin); per-user library caches under the cases' `home/.cache` (fontconfig)
are not compared. The render cases' figures and videos are recorded by metadata only (PNG size, JPEG / TIFF kind, SVG
and PDF page size in points, PostScript bounding box, ffprobe's video metadata without the codec
name and faststart; the parameter log as its run's `outputs` and `timestamp`; the step log in the
project file with `time` / `duration_s` masked), with MDBridge's rendering improvements patched into
MDATK (`improvements_render.py`); `tests/cli/test_commands.cpp` runs those that draw nothing, the
others need the Qt program (CI's `parity` job; FFmpeg's `ffprobe` on PATH). Each case names the command line, the inputs (files of `tests/data` copied into an
empty work folder, files written there: `setup`, with `{work}` for the folder; `executable` ones
made executable), command lines run first in the same folder (`before`, with their exit statuses
`before_exits`), environment variables (`env`) and MDATK's exit status, stdout, stderr and the files
it left (text, or SHA-256 and size beyond 16 KiB). The analysis cache is MDBridge's own
(`src/analysis/include/mdbridge/analysis/cache.hpp`), so its files are compared by meaning: names
with `{manifest}` for `MDATK_cache.json` / `MDBridge_cache.json` and `<name>-{key}.npz` / `.json`
for the entries (`.mdbridge` dropped), the manifest and the entries' descriptions as `json` (project,
version, time stamps, keys and `mtime_ns` masked; `all_entries` in any order), the archives as `npz`
(per array name, dtype, shape and the flat values as one line of JSON text; floats compared within
1e-9 relative, pair maps per frame as sorted entries) and several entries of one analysis as
`several` (matched in any order). MDATK runs in-process with the independent tables applied (first
section). Normalised: the work folder is `{work}`, `saved` time stamps `{time}`, file-name time
stamps `{stamp}`, cache keys `{key}` (they hash the inputs' mtimes); `posix_only` marks cases with
absolute paths, OS error texts, project paths or environment variables (skipped on Windows; the
cache manifest's path printed by `analyze` and the paths in cache files are compared there in
generic form); `exit_only` cases compare the exit status alone (usage errors, and the six gzip
errors `src/io/include/mdbridge/io/mrc.hpp` words differently); where MDATK crashes with a
traceback (`stderr` null, `mdatk_crash` the exception), MDBridge must print `error: ...` and exit
with 1. The chem commands of milestone 6 (`pdb clean / check / to-lrna / annotate`, `build helix /
mutate / fit`, `dict list / smiles / finish / pubchem`, `density fit / simulate`, and `render` with a
map and `auto_fit`) have their cases too, marked `pending` with the commands they `needs` (cli_parity.py
`PENDING`; skipped when comparing and by the C++ test until the commands land): their inputs come from
`build helix` run first, `chem/maps` and `chem/dictionaries`; MDATK runs with `improvements_chem.py`
and the command line in `sys.argv` (the step log records it), the step log's `time` / `duration_s`
are masked, and `dict pubchem` asks a stand-in PubChem on 127.0.0.1 (`MDBRIDGE_PUBCHEM_URL`, owner
decision D7; `{pubchem}` in `env`). Cases with `fitted` (rebuilt helices: 0.05; density fits: 0.6)
compare their texts with their decimal numbers within that absolute tolerance (D10) and keep files of
any size as text. numpy's and OpenBLAS' kernels are pinned (`chem_common.pin_kernels`, as in the chem
generators): with other kernels the least-squares closures of the non-Watson-Crick and triplex builds
land on other solutions, up to 0.8 Å apart (CI run 45882807857; `OPENBLAS_CORETYPE=Haswell` on the
reference machine), which no tolerance covers; pinning changed only the fitted cases. The same script compares a build with MDATK directly (CI's `parity` job, against the
AppImage):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/cli_parity.py --mdatk /tmp/mdatk --write tests/data/cli/cli_cases.json
/opt/mdbridge-venv/bin/python -I -B tests/parity/cli_parity.py --mdatk /tmp/mdatk --mdbridge build/tests/mdbridge_commands_cli
```

`list` prints the panel types and layout presets of MDATK's figure from `src/core/panel_tables.inc`
(also the pipeline's panel types, later the renderer's grid), which `tools/gen_panel_tables.py`
writes from `mdatk.render.panels.PANEL_TYPES` and `mdatk.render.layout` (`PRESETS`, `PRESET_LABELS`,
`HEIGHT_RATIOS`); it imports MDATK, so it needs MDATK's dependencies (matplotlib):

```bash
/opt/mdbridge-venv/bin/python -I -B tools/gen_panel_tables.py /tmp/mdatk   # src/core/panel_tables.inc
```

`chem/secstruct_cases.json` holds what MDATK's secondary-structure code (`mdatk/secstruct.py`)
computes, for `src/chem/secstruct.cpp` (`tests/chem/test_secstruct.cpp`), one case per line:

| Key | Content |
|---|---|
| `constants` | `OPENERS`, `CLOSERS`, `UNPAIRED`, `CHAIN_SEPARATOR`, `PAIR_RULES`, `C1_WINDOW`, `C1_IDEAL`, `SCORE_TIE`, `BASE_PAIR_ATOMS`, `KIND_LABELS` |
| `residue_names` | `base_letter` / `canonical_name` of every name of up to three characters over the letters of the naming schemes (and others), lower case, aliases, padded and non-nucleotide names |
| `lw_kind` | `lw_kind` of Leontis-Westhof codes |
| `parse` | `dotbracket_to_pairs` of MDATK's test texts, every bracket level, `&` and Unicode whitespace separators, invalid characters (repr() in the message), random balanced, mutated and random texts: pairs, n, chain breaks or the message |
| `write` | `pairs_to_dotbracket` of MDATK's test pairs, random nested and crossing sets, 26 to 32 mutually crossing pairs (30 levels fit), reversed, duplicated and invalid pairs, chain breaks in and out of range, n <= 0: the text or the message |
| `stems`, `crosses`, `unique_pairs` | `_stems`, `_crosses`, `unique_pairs` of pair lists |
| `detector` | `PairDetector(resnames, atom_index, cutoff, wobble)` on the geometries of MDATK's `tests/test_basepairs.py` (its own helpers build them), C1'-C1' distances and H-bond distances exactly at the limits and one ulp beyond, missing atoms, `C1*` names, an infinite cutoff, non-finite coordinates, near-tie sweeps, L-RNA and every residue-name scheme, perturbed helices, dense random clouds: per frame the coordinates, `candidates` (i, j, kind, score, deviation) and `detect` |
| `select` | `detect` with `candidates` replaced (as MDATK's tie test does): scores on, around and half-way between points of the 0.05 grid, equal deviations |
| `lattice` | two large detector cases the C++ test rebuilds from a formula (multiples of 0.25 A): counts, the first rows and SHA-256 of the candidate array (float64) and of `detect` (int64) |

Floats are `float.hex()` strings (`nan`, `inf`, `-inf`); the test compares everything exactly, the
scores and deviations bit for bit (MDBridge computes the distances in numpy's order,
`sqrt((dx² + dy²) + dz²)`, without fused multiply-adds). The geometries built by MDATK's test helpers
go through numpy's LAPACK / BLAS, which differ in the last bits between CPUs, so the generator rounds
them to float32 (as MDBridge's frames are) before use; perturbations are added in float64 afterwards.
MDATK's tests import pytest, which the generator replaces with a minimal stand-in when it is not
installed (only decorators run at import). To regenerate, with numpy and MDAnalysis in
`/opt/mdbridge-venv` and an MDATK checkout at `v0.0.4.1` (nothing is written into it; the output is
reproducible):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_secstruct_cases.py /tmp/mdatk tests/data/chem
```

The chem module of milestone 6 (PDB hygiene, annotation, the structure builder, restraint
dictionaries, the density fit; docs/specs/m6-plan.md 4.2) is checked against the other files of
`chem/`, written by five generators that import MDATK's chem code with MDBridge's improvements patched
in (`tests/parity/improvements_chem.py`: CH-01 ... CH-15, D7, D11 ... D13, NA-01's names, and the reference-data differences CHD-01, CHD-02 from `tools/gen_chem_tables.py`; each item
names the MDATK version whose design it follows). Their contents are described in each generator's
docstring:

| File | Generator | Content |
|---|---|---|
| `pdbtools_cases.json` | `make_pdbtools_cases.py` (pytest: MDATK's `tests/test_pdbtools.py` fixtures) | inputs; `scheme_name` of every base x hand x terminal x scheme; `chiral_volume` / `sugar_chirality`; `parse_ranges`; each PDB-tools function alone on the inputs (notes, residues, `write_pdb` where atoms or records change); `clean` over MDATK's test calls and an option grid. Read by `tests/chem/test_pdbtools.cpp` (every case byte for byte) |
| `annotate_cases.json` | `make_annotate_cases.py` | built helices and crafted files; `annotate` with every option and engine (records, notes, text) |
| `dict_cases.json`, `dictionaries/*.cif` | `make_dict_cases.py` (pytest: `tests/test_dictionaries.py`'s synthetic dictionaries) | `COMPONENTS`; `parse_smiles` (test strings, the components, a fuzzed set); `match_graph`; `check_smiles` / `identify_smiles`; `finish_dictionary` of synthetic acedrg / eLBOW dictionaries with every option; two dictionaries for `cli_parity.py`'s `dict finish` cases |
| `builder_cases.json` | `make_builder_cases.py` | `check_helix` and the sequence functions; the A-form templates; Watson-Crick `build_helix` over its options; `mutate` and `parse_mutations` |
| `builder_fitted_cases.json`, `maps/{helix_moved,strand_moved,strand_nan}.mrc` | `make_builder_cases.py` | the reference helix fits of the other pairings (start, the refined parameters `least_squares` reaches, cost, screw, torsions, contacts), `_pair_fit` of every canonical pair, the rebuilt helices (other pairings, triplexes: text, notes, heavy-atom coordinates), `fit_to_map` on the maps |
| `density_fit_cases.json`, `maps/fit_*.mrc` | `make_density_fit_cases.py` (MDAnalysis: the demo frame) | `rotation_samples`, `_fft_size`, per model x map the `FitProblem`, CCs, flips, the coarse search's grid and peaks, `refine`, `fit_model` global and manual; errors (read by `tests/chem/test_density_fit.cpp`: fits within D10's tolerance) |

numpy's and OpenBLAS' kernels are pinned (`chem_common.py`, as `make_analysis_cases.py`). The
least-squares closures, Powell's path in the density fits, Kabsch's SVD and the helices' 3 x 3
products still depend on the CPU in the last bits, so `regenerate_all.py` compares these JSON files
(not `pdbtools_cases.json`, pure Python) within the chem tolerance where their bytes differ
(`FITTED`: numbers within a relative 1e-6, decimals inside strings within one unit of their last
digit); the maps are made without BLAS products and compared byte for byte. The C++ tests compare
fitted results with D10's tolerances (non-Watson-Crick helices: heavy-atom RMSD <= 0.01 Å, every
atom <= 0.05 Å; fits: 0.5°, 0.2 Å, CC 1e-3; plan 4.1). To regenerate, with an MDATK checkout at
`v0.0.4.1` (nothing is written into it):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_pdbtools_cases.py /tmp/mdatk tests/data/chem   # needs pytest
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_annotate_cases.py /tmp/mdatk tests/data/chem
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_dict_cases.py /tmp/mdatk tests/data/chem       # needs pytest
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_builder_cases.py /tmp/mdatk tests/data/chem
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_density_fit_cases.py /tmp/mdatk tests/data/chem
```

The chemistry reference tables of milestone 6 are not test data but are regenerated (and compared)
by the same job: `tools/gen_chem_tables.py` writes `src/chem/helix_tables.inc` (the Arnott A-RNA
monomer from `data/sources/helix`, the 4JRD nucleotide and screw from `data/sources/pdb`) and
`src/chem/component_tables.inc` (the CCD's L-ribonucleotides 0A 0C 0G 0U from
`data/sources/ccd/ribonucleotide_*.tsv`), pure Python and deterministic, and compares them with
MDATK 0.0.4.1's embedded tables (`builder.py`, `dictionaries.py`): a difference that is not one of
its documented ones (docs/DIFFERENCES.md CHD-01, CHD-02) fails the step. `tests/chem/test_reference_tables.cpp`
checks the tables' chemistry (bond lengths, links, compositions).

```bash
python3 -I -B tools/gen_chem_tables.py /tmp/mdatk   # src/chem/helix_tables.inc, component_tables.inc
```

The analysis module of milestone 3 (`src/analysis/`, tested by `tests/analysis/`) is checked against
the files of `analysis/`, written by six generators (in the order `regenerate_all.py` runs them):

| Generator | Files | Tests |
|---|---|---|
| `make_analysis_cases.py` | `micro_cases.json`, `demo_dry_cases.json`, `demo_dry_variants.json`, `systems/`, `systems_cases.json`, `spatial_cases.json` | `test_contacts_hbonds.cpp`, `test_mg_stacking.cpp`, `test_energy.cpp`, `test_rmsd_analysis.cpp`, `test_structure.cpp`, `test_basepairs.cpp`, `test_parallel.cpp` |
| `make_distance_cases.py` | `distance_cases.json` | `test_periodic.cpp` |
| `make_rmsd_cases.py` | `rmsd_cases.json` | `test_rmsd_qcp.cpp` |
| `make_universe_cases.py` | `universe_cases.json`, `universe_files/` | none: the record behind `docs/specs/m3-universe-notes.md` |
| `make_analysis_framework_cases.py` | `framework_cases.json` | `test_framework.cpp` |
| `make_spatial_cases.py` | `spatial_extra_cases.json` | `test_spatial.cpp` |
| `make_pipeline_cases.py` | `pipeline_cases.json` | `test_pipeline.cpp` |

MDBridge's improvements of the analyses (docs/DIFFERENCES.md, "Analyses") are applied on the test
side by `tests/parity/improvements_analyses.py`: `make_analysis_cases.py`,
`make_analysis_framework_cases.py`, `make_pipeline_cases.py` and `cli_parity.py` (its MDATK child)
patch MDATK with it, and `make_rmsd_cases.py` records MDBridge's values for the improved RMSD cases
("improved", "near", "mdanalysis"). The energy step of `tests/parity/improvements_periodic.py` takes
integer powers through libm's `pow` (`math.pow`, as `energy.cpp` does) and the net force's length as
`sqrt((x*x + y*y) + z*z)`: numpy's `**` uses SVML's `pow` on AVX-512 CPUs and `np.linalg.norm` of a
vector BLAS' `ddot`, both of which made `cli/cli_cases.json` depend on the CPU. The covalent radii of the bond guess come from
`data/sources/covalent_radii` through `tools/gen_covalent_radii.py` (`src/analysis/covalent_radii.inc`);
the energy analysis' parameters without parm7 (AN-25 - AN-27) from `data/sources/forcefields` through
`tools/gen_energy_lj.py` (`src/analysis/energy_lj_tables.inc`), whose `load()` `tests/parity/energy_parameters.py`
uses to apply the same rules on MDATK's side.

The tests share one helper header, `tests/analysis/analysis_test_support.hpp` (case files, micro
systems, case runs, tagged JSON, array comparison, margins). Every test that goes through the
analysis runner runs twice under ctest, with the frame-parallel runner (`MDBRIDGE_THREADS=4`) and
sequentially (`<test>/threads1`), and must pass both ways; the spatial, periodic-search and QCP
tests run once.

`make_analysis_cases.py` records what MDATK's ten analyses (`mdatk/analysis/*.py`: contacts, energy,
hbonds, mg, rmsd, stacking, trace, proximity, atoms, basepairs) compute (`docs/specs/m3-plan.md`,
4.1 / 4.2). MDATK runs in-process with the independent tables applied and the residue-scoped ion
rule on every topology file (first section).

| File | Content |
|---|---|
| `micro_cases.json` | the micro systems of MDATK's `tests/test_analyses.py` and `test_basepairs.py` as data (`systems`: atoms with name, resid, resname, element, mass, charge; segid `SYS`; bonds; dt; float32 coordinates per frame, as `conftest.make_universe` builds them; the C++ test builds a `Topology` and an in-memory trajectory from them) and the analyses run on them: contacts and the native fraction, H-bonds linear / bent / far with and without bonds and with the angle off, Mg bridges / octahedral geometry / outer shell / occupancy / no ions, stacking stacked / offset / tilted / inter / no rings, analytic Coulomb + LJ (generic parameters, eps 4r), bonded exclusions 0-4 (and `_bonded_exclusions` itself, depths 0-5), RMSD superposition / fit selection / none / mass-weighted / `rmsd@G` instances, the base-pair duplex (strides 1 and 2, wobble off, cutoff 3.0, every analysis), structure residues only (`visual.structure_groups`), an unknown proximity set (error), L-RNA names, wobble and cutoff |
| `demo_dry_cases.json` | the demo project (`mdatk.demo.demo_project`, as `mdatk demo` writes it for `build_demo(n_frames=24, solvent=False)`) on `topology/demo_dry.parm7` + `io/trajectories/demo_dry.nc`: every analysis, all 24 frames, in one pass |
| `demo_dry_variants.json` | the demo project with settings changed: the DCD / XTC / TRR copies of `io/trajectories/files` (every analysis), `use_pbc` (the demo's orthorhombic cell), frame slices (start / stop / step, negative start, dt given), RMSD mass-weighted / fit selections / no superposition / reference frames 5, 999, -3 (clamped) / odd target labels and an empty selection / no targets, contacts A = A / all atoms / cutoffs 3 and 8, H-bonds without angle / strict / other element lists, Mg geometry / tolerance / outer shell / occupancy selection / cutoffs / no ions, stacking inter / strict / no rings / one group, energy generic LJ / exclusions 0, 1, 5 / plain Coulomb / cutoff 8, base pairs stride 1 and 7 / no wobble / cutoff 3.0, trace without alignment / ions / with `trace_atoms` P or none, structure groups (Target only, unknown), proximity sets `sel:resid 3` / a group with radius 4 / unknown, atoms with hydrogens, instances `X@Target+Aptamer`, `rmsd@G`, and errors (instance groups unknown / wrong count / not taken, unknown analysis, group B inside A for contacts / H-bonds / energy, unknown group, empty group) |
| `systems/` | new small systems (generated, no research data): `periodic_ortho` / `periodic_tric` (`.parm7` written by `mdatk.demo.write_topology`, `.nc` by MDAnalysis' NCDFWriter: `build_system("GCAU", n_frames=6, seed=11)`'s 4-bp duplex and its 3 Mg ions moved to straddle a cell corner and wrapped into a 25 Å cube / a 25 x 27 x 29 Å cell with angles 60 / 70 / 80), `pdb_nobonds.pdb` (demo frame 0: no CONECT, no element column, chain IDs A / B / C, CRYST1), `.gro` (MDAnalysis' GROWriter) and `.dcd` (8 frames), `noh.pdb` (the heavy atoms, element column) and `.dcd` (frames 0, 6, 12, 18), `noh_conect.pdb` (`noh.pdb` without chain IDs, CONECT records for one bond of residue 24 only: chains from bonds, TOP-01), `aromatic.pdb` (3 models: G / C / A / U and HSD with PHE / TYR / TRP / HIE / HIS stacked on them, offset, tilted, too far, one ring atom lifted out of plane) |
| `systems_cases.json` | the analyses on `systems/`: periodic systems with `use_pbc` off and on, cutoffs up to 12 Å against the 25 Å cell, Mg geometry and H-bond angles across faces (minimum image only for orthorhombic cells, by MDATK's design); PDB / GRO topologies without bonds or charges, alone and with the DCD (elements and masses guessed, hydrogens assigned from distances, bonds guessed, generic LJ, no-charge warning, chain IDs), mass-weighted RMSD, `use_pbc`; heavy atoms only (the H-bond "no hydrogens" fallback); the aromatic rings (ring sizes 5, 6 and 9, all / inter / strict / proteins only) |
| `spatial_cases.json` | (not read by the C++ tests, which use `spatial_extra_cases.json` and `distance_cases.json` below) `mdatk.spatial`: `_usable_box` of valid and unusable boxes; `pairs_within` (i, j, d in MDATK's order, and their dtypes) of random points without a box, in orthorhombic and triclinic cells and with unusable boxes, empty inputs, a == b, an exact cutoff distance, pairs across a face, points far outside the cell, duplicated points, cutoffs at and beyond half the cell; `minimum_image_vectors` with half-length ties (rounding half to even), angles within and beyond `np.allclose`, triclinic and unusable boxes; `kabsch` of 1, 2, 3, 5, 20 points with and without weights, identical, reflected, collinear and planar sets, zero weights |

Each case is one run of `run_analyses` built as `Pipeline.analyze` builds it: the base project with
`--set`-style overrides (`mdatk.cli.apply_override`; listed per case), the instances with
`instance_project`'s projects and their own `AnalysisContext`s sharing one Reporter, the frames of
`frame_indices(inputs.start, stop, step)` and a fresh Universe. Recorded per case: the project
sections the analyses read (tagged JSON), the analysed frames and their times, the Reporter's log
lines (`Preparing ...`, `Analysing ...`, `WARNING: ...`) and progress calls, Python warnings, the
trajectory frame current at each setup, or the error (`{type, message}`); per result: params and
meta (tagged JSON; the meta holds the outputs and labels), the order of the outputs, the cache key
for a fixed input fingerprint (recorded: the fingerprint's shape with paths `{data}/...` /
`{micro}/name`, the real sizes and `mtime_ns` 0), and every array as `{dtype, shape, sha256,
values}` (numpy `dtype.str`; integers and booleans as JSON, floats as `float.hex()`, strings as
strings). `margins` holds per result and analysed frame the smallest distance of a compared
quantity from its threshold (cutoff distances, cosines of angle cutoffs, planarity, lateral offset,
H-guess distances and ties, the base-pair C1' window, score grid and score ties, trace links, guessed
bonds), and the energy's sums of |term| for relative tolerances: a count that differs where the
margin is tiny can be a last-bit flip, anywhere else it is a bug (a margin of exactly 0 is an exact
tie built into the geometry, e.g. MDATK's contacts test at 3.5 Å, and must match exactly). The
cosines MDATK compares with are recorded as numpy computes them and as libm's
`cos(x * pi / 180)`. The full list is in the generator's docstring.

numpy's SIMD kernels (`power`, `cos`) and OpenBLAS / LAPACK's kernels (SVD of `kabsch` and of the
ring planes, BLAS products) are chosen per CPU and differ in the last bits (relative 1e-16 to 1e-13
on the reference machine: energies, aligned coordinates, RMSD with a fit selection, stacking
geometry). The generator therefore pins them before numpy loads (`OPENBLAS_CORETYPE=Nehalem`, one
thread, numpy's dispatch limited to its x86-64 baseline with `NPY_ENABLE_CPU_FEATURES=X86_V2`), so
the data is the same on almost every x86-64 machine (one CI runner still moved RMSDs with a fit
selection by a few ulps, so `regenerate_all.py` compares these files within a relative 1e-12 where
their bytes differ; below); `--native` skips the pinning to show what a machine's own kernels
change (not the committed data). The C++ tests compare such values within the
tolerances of the plan (4.3).

To regenerate (after the topology fixtures and the trajectory cases, whose `demo_dry` files it
reads; MDATK's test helpers are imported, with a stand-in for pytest when it is missing; the output
is reproducible):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_analysis_cases.py /tmp/mdatk tests/data/analysis
```

`analysis/framework_cases.json` holds what MDATK's analysis framework (`mdatk/analysis/base.py`,
`__init__.py`, `common.py`) computes, for `src/analysis/` (`tests/analysis/test_framework.cpp`):
the registry (names, labels, versions, `internal`, group slots, output-label templates); for five
projects given as the plain data `Project.from_dict` reads (the defaults, nested and self-referencing
groups with unusual names, disabled analyses, repeated group names, no groups) every analysis's
`enabled`, `repr(params_for)`, `selection_texts`, `expected_outputs`, `groups`, `output_label` of its
outputs and of unknown names, `cache_key` for two fixed input fingerprints (no file is read),
`referenced_groups` and `new_result`; `expected_keys` with and without disabled analyses;
`referenced_groups` of extra texts; per instance id `parse_instance`, `instance_id` and what
`instance_project` gives (the same description, or the KeyError / ValueError); `key_label` of data
keys; `split_key`; `residue_rows` and `without` on index lists of a small fake topology;
`per_row_counts` of random CSR triplets (float32 values as `float.hex()`); and `run_analyses` of
recording analyses over a MemoryReader trajectory (log lines, progress calls, setup / process /
finalize calls, times). MDATK runs in-process with the independent tables applied (first section).
The data keys must equal MDATK's exactly; the analysis cache itself is MDBridge's own.

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_analysis_framework_cases.py /tmp/mdatk tests/data/analysis
```

`analysis/pipeline_cases.json` holds what MDATK's analyze pipeline (`mdatk/pipeline.py`
`Pipeline.required_analyses` / `analyze` / `load_cached_dataset`, `mdatk/io/cache.py`) does, for
`src/analysis/pipeline.cpp` and `cache.cpp` (`tests/analysis/test_pipeline.cpp`): `required` lists
`required_analyses(video, charts)` of projects covering every panel type (structure with and without
representations and proximity, sequence, secondary with each `source`, time series, dual axis,
heatmap, contact map), charted series, analyses switched on and off and invalid instances (the
"skipping data source" lines, sorted: MDATK logs them in the order of a Python set), and an unknown
panel type; `sequences` are runs of `Pipeline(project).analyze(names, force)` one after the other on
one cache folder (cache reuse, invalidation by a parameter and by the frame selection, a display-only
change, `force`, instance lists, a repeated instance, no analyses, frame selections with external
data, a reused cleaned trajectory, and errors: `check()`, unknown analyses and groups, an empty frame
selection, bad external data). Per run: the project, the log and progress lines, the cache keys, the
manifest (`MDATK_cache.json` without `project`, `mdatk_version` and `updated`; `all_entries` sorted by
the masked keys, since MDATK sorts by the real ones), each entry's `.json`
description, the entries in the cache folder, the dataset (frames, times, keys, labels, units,
energy target) and what `load_cached_dataset` rebuilds. Paths are `{data}` / `{work}` placeholders;
cache keys are those of the fingerprints with placeholder paths and `mtime_ns` 0 (the C++ test masks
its own fingerprints the same way), and every `mtime_ns` is written as 0. Only the demo files'
contents are read, through MDATK's loader with the independent tables and the residue-scoped ion
rule; the analyses run are contacts, energy, H-bonds, RMSD, trace and proximity. The cache files
themselves are MDBridge's own (`src/analysis/include/mdbridge/analysis/cache.hpp`) and not compared.

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_pipeline_cases.py /tmp/mdatk tests/data/analysis
```

`analysis/spatial_extra_cases.json` holds what MDATK's spatial helpers (`mdatk/spatial.py`) compute,
for `src/analysis/spatial.cpp` and the numerics under it (`core/linalg.hpp`,
`core/detail/numpy_sum.hpp`, `core/detail/neighbors.hpp`; `tests/analysis/test_spatial.cpp`),
written by `tests/parity/make_spatial_cases.py`, which runs MDATK's functions in-process. Every float
is `float.hex()`. `pairs`: `pairs_within` without a box (scipy's cKDTree) on random float64 and
float32-widened points, clusters, duplicates, a lattice at distances equal to the cutoff, crafted
pairs within an ulp of the cutoff (including pairs where `sqrt(d2) <= r` and `d2 <= r * r` disagree:
scipy compares squares), a spread that defeats a cell grid, special cutoffs (0, negative, infinite,
NaN) and non-finite points (scipy's error); the pairs in MDATK's (tree) order, which the test sorts:
MDBridge returns them sorted by (i, j) and must match the set and every distance exactly.
`min_image` / `usable_box`: `minimum_image_vectors` and `_usable_box` for boxes near and away from 90
degrees (numpy's `allclose` margins), unusable boxes and half-box ties. `kabsch`: random, weighted,
mirrored, coplanar and tiny sets, fewer than 3 points (compared within a relative 1e-14).
`rings`: `ring_geometry` of the stacking analysis (batched `np.linalg.svd`; planarity within 1e-14,
normals up to sign within 1e-13). `sums`: `np.sum` / `np.mean` of float64 arrays of 0 to 300 and
some larger lengths, whose values the test computes itself (`value(k)` in the generator), bit for
bit. `tests`: the inputs of MDATK's `test_pairs_within_matches_brute_force_and_pbc` (open and
periodic halves) and `test_kabsch_recovers_rotation`, which the test ports.

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_spatial_cases.py /tmp/mdatk tests/data/analysis
```

`analysis/distance_cases.json` holds what MDAnalysis 2.11.0's `capped_distance` returns for the calls
MDATK's `pairs_within` makes with a periodic box (`mdatk/spatial.py`; every case is also run through
`pairs_within` and must agree), for MDBridge's periodic pair search
(`docs/specs/m3-capped-distance.md`; `tests/analysis/test_periodic.cpp`, lattice points
x fastest): the automatic method rule (sizes 9 / 10, cutoff around
0.3 x the cell, n x m >= 1e8 on generated lattices), points at and one float32 step around the
cutoff along axes and across faces, random sets inside, near and far outside orthogonal and
triclinic cells, cutoffs beyond half a cell, the same array on both sides, duplicates, dense sets for
the grid method's order, and degenerate cutoffs (0, negative, NaN, infinite), angles, lengths and
coordinates. Each case gives the method chosen, the pairs in MDAnalysis' order and the distances as
`float.hex` (or the error and warnings); lattice cases give SHA-256 hashes; `band` lists the pairs
within the spec's tolerance of the cutoff. That is the reference behaviour; MDBridge's own periodic
search (docs/DIFFERENCES.md AN-01 to AN-03) is tested against each case's `exact` column, computed
independently by `tests/parity/improvements_periodic.py` (numpy). That module is also patched into
MDATK by `make_analysis_cases.py`, `make_analysis_framework_cases.py`, `make_pipeline_cases.py` and
`cli_parity.py`, so their `use_pbc` results, energy forces and the bumped analysis versions are
MDBridge's (AN-04, AN-05, AN-11, AN-12), and adds the `improved` column of the `min_image` cases of
`spatial_extra_cases.json`. Format in the generator's docstring:

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_distance_cases.py /tmp/mdatk tests/data/analysis
```

`make_rmsd_cases.py` and `make_universe_cases.py` write the probe data of the milestone-3 behaviour
specifications `docs/specs/m3-rmsd.md` and `m3-universe-notes.md`:

| File | Content |
|---|---|
| `rmsd_cases.json` | `docs/specs/m3-rmsd.md`: MDAnalysis 2.11.0's `rms.rmsd(a, b, weights, center=True, superposition=True)` of about 240 coordinate pairs (general structures with noise from 0 to 5 Å, identical / translated / rotated copies, planar sets, one and two atoms, collinear and coincident atoms, weights leaving one or two atoms, NaN / infinite inputs and weights, negative, zero-sum, tiny and huge weights, empty inputs, shape and weight errors) with the result, error or numpy warnings, `scale` (E0 / N, the scale of the tolerance) and a `degenerate` flag; and the series MDATK's rmsd analysis computes on the demo trajectory (targets as atom-index runs, masses or no weights, reference frames 0 / 5 / 23) |
| `universe_cases.json` | `docs/specs/m3-universe-notes.md`: iteration over `trajectory[int64 array]` on the demo NetCDF / XTC / DCD / TRR files, a multi-model PDB and a GRO file (order, repeats, negative and out-of-range indices, stopping early; the current frame afterwards); cells read back from files with zero lengths / angles; per topology (parm7, PDB with and without CONECT, GRO) whether charges and chain IDs exist, the errors when not, bonds in order, fragment indices and, per atom group, its bonds, fragment indices, residues, residue atoms, boolean-mask subsets and unions |
| `universe_files/*` | the crafted inputs of `universe_cases.json`: `frag_conect.pdb` (CONECT records out of order, four fragments), `conect_lonely.pdb` (a CONECT record without partners), `cell_*.nc` / `.dcd` / `.xtc` / `.trr` (cells with a zero length, zero lengths, all zeros, a zero angle, orthorhombic) |

Floats are `float.hex()` strings. Both generators run MDAnalysis with the independent tables (first
section); neither imports MDATK (its checkout is only the common first argument). Nothing recorded
goes through LAPACK or BLAS, whose kernels are chosen per CPU: the random rotations come from
quaternions and are applied elementwise, and the exact superposition RMSD of the improved values is
computed in plain floats (`improvements_analyses.exact_rmsd`, Horn's quaternion by Jacobi rotations
as in `rmsd_qcp.cpp`), so the files do not depend on the CPU (`make_rmsd_cases.py --report` prints the
spread between MDAnalysis and that exact superposition that the spec's tolerances rest on). To regenerate, with
Python 3.13 and MDAnalysis 2.11.0 in `/opt/mdbridge-venv` (the output is reproducible; the universe
cases read `io/loader/files` and `io/trajectories`, so run those generators first):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_rmsd_cases.py /tmp/mdatk tests/data/analysis
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_universe_cases.py /tmp/mdatk tests/data/analysis
```

## Rendering (milestone 4)

`render/` holds what MDATK's presentation layer (`mdatk/render/`) shows and writes, for MDBridge's
`src/render` and the `render`, `charts` and `preview` commands (docs/specs/m4-plan.md section 3).
Image parity is **contents only** (owner decision D10): figures are compared with MDATK's by what
they contain - the data shown, colours, ticks and their labels, legends, texts and the layout
boxes - never by pixels, and there are no reference images. The 3D layer of the structure view is
the exception: MDATK's numpy software renderer is compared pixel by pixel in fixed small
viewports, within the plan's tolerance (99.5 % of pixels equal, all within 2/255).

The four generators run MDATK in-process through `tests/parity/render_common.py`, which pins
numpy's and OpenBLAS' kernels (as `make_analysis_cases.py` does), gives Matplotlib MDBridge's figure
fonts (`tests/parity/improvements_fonts.py`: the bundled IBM Plex Serif and Plex Mono faces of
`data/fonts` only, IBM Plex Math as their fallback, residue tick labels in the sequence font; never a machine's fonts or the wheel's DejaVu; docs/DIFFERENCES.md
FONT-01 ... FONT-03), and patches MDATK like the pipeline
generator (independent tables, analysis / selection / CLI improvements) and with MDBridge's
rendering improvements, `tests/parity/improvements_render.py` (items RN-01 to RN-17 of the plan;
patched so far: RN-01, RN-04 and RN-14, the plot panels, and RN-40 ... RN-42, the sequence /
secondary panels; the other items are MDATK's own
behaviour, including the cases that show its bugs, until their chunk fills in its patch with its
C++ code and docs/DIFFERENCES.md entry and regenerates). Every record lists the active items under
`improvements`. The Matplotlib artist walker `tests/parity/mpl_walker.py` turns a drawn figure into
data (its docstring describes the records). JSON is written by `render_common.dumps`: floats as
Python's repr, non-finite ones as `{"$float": "nan"}`, short objects on one line. Arrays are in
`.npz` archives with stored (not deflated) members and fixed time stamps.

| Folder | Generator | Content |
|---|---|---|
| `inputs/` | `make_figure_cases.py` | the datasets every figure and chart is built from, so the C++ tests read exactly what MDATK drew: `<name>/dataset.json` (frames, time_ps, time_unit, per result name, params, meta, array names; cache keys as `{key}`) and `<name>/<result>.npz`; float64 arrays rounded to 10 significant digits. `demo_dry` (the demo project's analyses on the demo files: every analysis plus trace, atoms, proximity, base pairs; 24 frames, 3 Mg ions), `demo_dry_1` (first frame only), `demo_dry_nan` (non-finite trace point and RMSD value), `demo_dry_times` (demo_dry's rmsd and contacts at the times of two concatenated runs: 10 ps, then 2 ps steps; RN-01), `protein_rna` (`analysis/systems/aromatic.pdb`: RNA and protein chains), `three_chains` (synthetic: chains A, B and LIG, changing base pairs), `example_<case>` (MDATK's example dataset of a case's project, the GUI preview); `demo_dry_map.mrc`, a map simulated from the demo's first frame |
| `figures/` | `make_figure_cases.py` | 100 figures (`<case>.json` + `<case>.npz`): every panel type and option, every layout preset, themes, title and watermark, sizes and DPIs, missing data ("No data: ..."), the structure view's colour, glow, centre, zoom, tags, 3D representations and density shadow, sequence wrapping, arcs within and between chains, protein and RNA, MDATK's errors (RN-03); `build_figure` and `update(pos)` at the first, middle and last frame, the figure walked after drawing (later positions: the animated artists' changes only) |
| `charts/` | `make_chart_cases.py` | `render_charts` cases (`<case>.json` + `.npz`): series and matrix charts, smoothing, targets, themes, styles (D5), labels and title templates, formats (with the files' first bytes and sizes), missing / matrix / pair-map series, stamps, time units, one frame, uneven times, a 700-row matrix (`inputs/tall`, the generator's own input); per chart the walked figure, the paths returned, log, progress and warnings; MDBridge's improvements RN-01, RN-05, RN-10, RN-11, RN-14, RN-50, D5 applied (improvements_render.py) |
| `molecule/` | `make_molecule_cases.py` | `molecule_cases.json`: topologies (backbone control atoms, rings, letters, segments, element radii and colours) of the demo, protein / RNA and schematic systems; scenes (`build_scene` + `frame()` per style and colour scheme: specs, residue masks, per item kind / count) with their geometry in `scenes.npz` (float32 coordinates, int32 indices); software-renderer images (`renders.npz`, uint8 RGBA) of scenes and of the synthetic items of MDATK's `tests/test_molecule_render.py`, with viewport, shading and the items rendered; MDATK's renderer patched with MDBridge's RN-26 (anti-aliasing with the real neighbours) and RN-27 (fast quality at odd sizes), `improvements_render.SOFTWARE_RENDERER`, listed in the file's `improvements` |
| `cli/` | `make_render_cli_cases.py` | `render_cli_cases.json`: `render` (full run, `--from-cache`, `--no-video`, stride, codecs, odd sizes, missing ffmpeg, PDBs), `charts` and `preview` (formats, frames out of range) on the demo files: exit status, stdout, stderr, the files written (PNG / SVG / PDF sizes, ffprobe's metadata and faststart of the videos, the params log's outputs), the step log; normalised like `cli_parity.py` |

`tests/render/test_figure.cpp` rebuilds the first frame of every `figures/` case in MDBridge's figure
model (grid, roles, chrome, texts, legends and the recorded artists) and compares its layout,
ticks, texts and legends with the record (the figure framework, W1-A). `tests/render/test_plot_panels.cpp`
builds every case whose panels are implemented from its input dataset and project, as
build_figure does, and compares the whole contents (`figure_compare.hpp`: data artists, images and
colour-bar tables of the `.npz`, legends, ticks, the panels' records, the animated artists at the
later positions, the log); `tests/render/test_seqpanels.cpp` does the same for the sequence and
secondary panels in MDATK's frozen axes boxes (the structure view, until it lands, as a stand-in).
`tests/render/test_structure_panel.cpp` builds the structure view (W2-B) of every case whose first
panel is one from the case's input and project and compares its axes with the record at every
recorded position (the other panels are placeholders there; the structure axes is put at
Matplotlib's box first). The structure view's improvements RN-03 (non-finite trace points: a
warning instead of MDATK's error, case `structure_nan`), RN-08 (ion legend), RN-09 (`Chain A`),
RN-16 (`vtk` means the GPU renderer) and RN-33 (density level from the finite voxels) are patched
in `improvements_render.py`.

`make_figure_cases.py` runs first (it writes `inputs/`); the cases are listed, with what each
shows, in each generator. The command cases need `ffmpeg` and `ffprobe` (CI's regenerate job
installs Ubuntu 24.04's package); MDATK's video bytes are not kept, only their metadata (the codec
is MDATK's, libx264: MDBridge encodes with its own, RN-07). To regenerate (MDATK's checkout as the
first argument; the output is reproducible):

```bash
for g in figure chart molecule render_cli; do
    /opt/mdbridge-venv/bin/python -I -B tests/parity/make_${g}_cases.py /tmp/mdatk tests/data/render
done
```

The JSON records of `figures/`, `charts/` and `molecule/` hold view rotations and projected
coordinates that go through BLAS products of 3 x 3 matrices; with numpy's AVX-512 kernels
(`--native`) they move by an ulp, so `regenerate_all.py` compares them, where their bytes differ,
like the machine-dependent analysis files below (numbers within a relative 1e-12). Their `.npz`
files did not change.

`render/colors/color_cases.json` holds what Matplotlib 3.11.2's colour functions and MDATK's
themes (`mdatk/render/theme.py`, `coloring.py`) give, for `src/render/core`
(`tests/render/test_colors.cpp`): about 2300 `to_rgba` strings (every named colour, case variants,
`none`, `C<n>`, hex of every length and case, grey levels in Python's `float()` syntax, Unicode
spellings, alpha overrides and invalid alphas, with Matplotlib's messages), RGB(A) tuples, `to_hex`,
79 colormaps (the 17 bundled maps of owner decision D6 and their reverses, `with_extremes`,
`from_list` of colours and of (value, colour) pairs, segment data with jumps; lookup tables for
N <= 16 and float, integer and `bytes=True` samples around every boundary), the errors of invalid
maps, MDATK's chain palettes `palette(name, n)`, the themes' tokens, series colours, sequential and
ordinal ramps, and the properties `style_axes` / `style_legend` leave on a Matplotlib axes and
legend. MDBridge's improvements (ordinal ramp for more than five classes, NaN tuples, the message
for unbundled map names) come from `tests/parity/improvements_colors.py` (docs/DIFFERENCES.md,
"Rendering: colours and themes"). Matplotlib's colour cache is cleared before every `to_rgba` (it cannot tell
`-0.0` from `0.0`). `tools/gen_mpl_tables.py` writes the named colours, default colour cycle and
bundled colormaps that the C++ code uses to `src/render/core/mpl_colors.inc` (Matplotlib's data
with its notices). To regenerate, with Python 3.13 and Matplotlib 3.11.2 in `/opt/mdbridge-venv`
(the output is reproducible and does not depend on the CPU):

```bash
/opt/mdbridge-venv/bin/python -I -B tools/gen_mpl_tables.py                    # src/render/core/mpl_colors.inc
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_color_cases.py /tmp/mdatk tests/data/render/colors

## render/text: the figure fonts and text metrics (milestone 4, W0-D)

`render/text/text_cases.json` holds what Matplotlib 3.11.2 does with the IBM Plex fonts MDBridge
bundles (`data/fonts`, IBM's release files, `data/sources/fonts/README.md`), configured as MDBridge
draws (`tests/parity/improvements_fonts.py`: Plex Serif first in the sans-serif / serif lists, Plex
Mono first in the monospace list, IBM Plex Math as the fallback font), for `src/render/core/font.cpp` and `text.cpp`
(`tests/render/test_text.cpp`), written by `tests/parity/make_text_cases.py` (no MDATK needed): the
seven files' roles, checksums, names and metrics (the committed `data/fonts` must have the release's
checksums); the face Matplotlib's `findfont` picks among the six Plex Serif and Plex Mono faces for
1080 font-property sets (family lists, weights, styles); and `get_text_width_height_descent` of 73
strings (MDATK's labels, tick labels, superscripts, ion charges from the fallback font, Greek,
ligatures, kerning pairs, sequences and dot-bracket, blank text) on every figure face at four of ten (points, dpi) sizes, with
`text.hinting: no_hinting` (width / height / descent and the glyph x positions in 26.6
fixed point, which MDBridge reproduces exactly) and with Matplotlib's default hinting (MDATK's
figures; MDBridge stays within a bound, docs/DIFFERENCES.md "Rendering: fonts and text"). Matplotlib's
rcParams are its defaults with MDBridge's font lists and the machine's fonts are left out, so the
output depends only on the Matplotlib wheel (its FreeType 2.14.3 and HarfBuzz) and the bundled
fonts, and is the same everywhere:

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_text_cases.py tests/data/render/text
```

## render/molecule_geometry: the 3D representations' geometry (milestone 4, W0-F)

`tests/parity/make_molecule_geometry_cases.py` runs MDATK's `mdatk/render/molecule/` (topology.py,
geometry.py, surface.py, scene.py) in-process, for `src/render/molecule` and
`tests/render/test_molecule_geometry.cpp` (vertices within 1e-5 Å, indices exact; docs/specs/m4-plan.md,
section 3):

| File | Content |
|---|---|
| `systems.json` | the inputs of each `Topology`: per-atom names, elements, residue names and rows, bonds, per-residue chain and resindex, groups, and the frames used. `demo_dry` comes from MDATK's `trace` and `atoms` analyses on `topology/demo_dry.parm7` + `io/trajectories/demo_dry.nc` (independent tables, residue rule, analysis improvements, `visual.align_frames = false`, so the coordinates are the trajectory's float32 values); frame 100 is frame 0 with three atoms set to NaN. `mixed` (a protein chain with a numbering gap and a far jump, four demo nucleotides, a ligand and an Mg ion), `charmm` (bases told by their atoms) and `odd_names` (old atom names, blank / lower-case elements, pseudouridine, a residue without atoms) are built in the generator |
| `topology.json` | `normalise_name` / `base_letter` of name lists; per system and backbone atom (P, C4', and an unknown one) every role of `Topology` |
| `styles.json` | one case per line: a style (cartoon, ribbon, tube, licorice, vdw, surf, trace) on a system, frame, residue mask, size, sugars, spline samples and backbone atom: `drawn()` and the `DrawList` (meshes, spheres, sticks) |
| `surface.json` | `density_grid` (shape, origin, the float32 values for small grids), `surface_nets` and `gaussian_surface` (with each vertex's nearest atom and its distance, for ties) of small point sets and demo residues |
| `scenes.json` | `MoleculeScene` with `residue_masks`, `fill_rest` ("trace" / "hide") and `frame(xyz, residue_rgb)`: the specs and masks, `drawn()`, per item its kind, owners, atoms / roles, radii, colours and opacity, and `depth_range` |
| `ndimage.json` | scipy's `gaussian_filter1d` / `gaussian_filter` (mode "constant") and `map_coordinates` (order 1, mode "nearest") on float32 arrays, for the port in `src/numerics/ndimage.cpp` |

Floats are `repr()` (exact). numpy's and OpenBLAS' kernels are pinned as in
`make_analysis_cases.py` (the spline weights use a small matmul, the Gaussian kernel `np.exp`);
`regenerate_all.py` compares `styles.json`, `scenes.json` and `ndimage.json` within a relative
1e-12 where their bytes differ. The output is reproducible; it reads the demo files written by the
topology and trajectory generators:

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_molecule_geometry_cases.py /tmp/mdatk tests/data/render/molecule_geometry
```

## render/density_shadow: the density shadow (milestone 4, W0-I)

`tests/parity/make_density_shadow_cases.py` runs MDATK's `mdatk/density.py` (`_above_level`,
`project_density`, `shadow_rgba`, `anchored`) and `scipy.ndimage.affine_transform` in-process, for
`src/render/molecule/density_shadow.cpp` and the affine resampling of `ndimage.cpp`
(`tests/render/test_density_shadow.cpp`: images within 1e-6, placements within 1e-9). The maps are
committed files read with MDATK's `read_map`: `render/inputs/demo_dry_map.mrc` (figure step) and
`io/maps/files/*` (map step: permuted axes, NSTART origins, constant, NaN, inf, float32 extremes,
NaN cell). MDBridge's improvements RN-30 ... RN-32 (docs/DIFFERENCES.md) are applied test-side;
records keep MDATK's own outcome under `mdatk` where they change it.

| File | Content |
|---|---|
| `affine.json` | `affine_transform(order=1, mode="constant", prefilter=False)` of float32 arrays: shape, values, 3 x 3 matrix, offset, output shape, cval and the float32 output (last grid lines exactly, outside, NaN offset, empty output) |
| `shadows.json` | `kind` `above_level`: `_above_level(map, level)` of small maps; `kind` `project`: `project_density` inputs (map, placement, view, centre, level, window, max_pixels, max_depth) and the image (rows bottom -> top), extent, `null` or the error, and `shadow_rgba` of some images |
| `anchored.json` | `anchored(values, anchor, reference).to_list(12)`, the error, or the RN-32 note |

Floats are exact (float32 values widened). numpy's and OpenBLAS' kernels are pinned as in
`make_analysis_cases.py`; `regenerate_all.py` compares these files as JSON within a relative 1e-12
where their bytes differ (3 x 3 products and kabsch's SVD):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_density_shadow_cases.py /tmp/mdatk tests/data/render/density_shadow

## render/shading: residue shading, orientation and the example dataset (milestone 4, W0-H)

`tests/parity/make_shading_cases.py` runs MDATK's `mdatk/render/coloring.py`, `residues.py`,
`orientation.py` and `example.py` in-process (patched by `render_common.setup()`, including
`improvements_shading.py`: RN-21 ... RN-25) and numpy's `default_rng`, for `src/render/core`
(`tests/render/test_shading.cpp`, bit for bit; BLAS rotations within 1e-15). No trajectory is read.

| File | Content |
|---|---|
| `shading_cases.json` | crc32 of keys; per seed the SeedSequence state, PCG64 raw outputs and a scripted series of `uniform` / `integers` / `normal` calls; `glow_intensity`, `static_colors` and `frame_colors` cases; `one_letter`; orientation (`rx` / `ry` / `rz`, `euler_to_matrix`, `wrap_degrees`, `matrix_to_euler` incl. gimbal lock, `drag`, `project`, `format_angles`); `group_resids`; per example trace `chain_names`, `resnames`, `residue_rows`, `tag_rows` (with their log lines) and `residue_shading` settings; a proximity result with other rows than the trace; the schematic trace and duplex |
| `examples.json` | per case (default layout, MDATK's test_v002 layout, every data kind, representations, a loaded preview with and without atoms, single strand, large / odd / no groups, one frame): the project, frame count, preview, and the dataset `example_dataset` makes (frames, times, unit, per result its name, key, params, meta and array names) |
| `arrays.npz` | the arrays (stored members): `<case>.<result>.<array>`, the long random sequences `rng.<seed>.<op>`, the shading arrays `shading.<case>.<visual>.<need>.<what>` and `remap.*` |

Floats are `repr()` (exact); numpy's kernels are pinned (`render_common`), so its sin / cos / exp are
the C library's. The output is reproducible:

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_shading_cases.py /tmp/mdatk tests/data/render/shading
```

## gui: action scripts for the window (milestone 5, W0-F)

`gui/` holds GUI cases: an action script run against MDBridge's main window and the files the
window must have written when it ends (`tests/gui`, program `mdbridge_gui_tests`, Qt Test on the
`offscreen` platform; docs/specs/m5-plan.md section 3). The test runs every `gui/**/*.json` file.
A case is one JSON object (format proposed with the harness; the generators of the per-page cases,
`tests/parity/make_gui_cases.py`, are being written in parallel and either use it or the runner is
adapted, `tests/gui/gui_harness.cpp`):

| Key | Content |
|---|---|
| `format`, `version` | `"mdbridge-gui-case"`, `1` (other values: the case fails) |
| `description`, `source` | what the case does and what wrote its expected files (text only) |
| `clock` | optional local time `YYYY-MM-DDTHH:MM:SS`: the project file clock (`core::setClock`) starts there and goes one second further on every reading, as MDATK's generators fix `_now()` |
| `files` | `{name in the case folder: path under tests/data}`: copied into a new, empty case folder before the first action |
| `actions` | the steps, in order (below) |
| `expected` | `{name in the case folder: path under tests/data}`: after the last action each file must equal it byte for byte, after two normalisations: line ends (CR LF on Windows) become LF, and the case folder's absolute path (generic form, `/` separators) in the written file becomes `{case_dir}` |
| `log` | optional, not checked yet: the log lines MDATK showed |

Actions (`"action"` names the step; an unknown one fails the case, naming it):

| Action | Keys | What the harness does |
|---|---|---|
| `open` | `path` (case folder) | the window's *Open project* with `path` chosen in its file dialog; it must succeed (unless the case needs `open_failure`) |
| `wait` | `timeout_s` (default 60) | until the window has loaded the project's structure and preview data (`MainWindow::loading()` false) and the preview has settled |
| `set` | `key`, `value` | finds the one widget whose `objectName` is `key` (field widgets use their project key, `analyses.contacts.cutoff`) and enters `value` as a user would: keys typed into line edits, spin boxes and editable combo boxes (select all, type, Enter), check boxes toggled with Space when their state differs, the item whose project value (item data) is `value` chosen in other combo boxes, plain-text edits replaced by typing; a key no widget has is a list editor's cell (`groups.3.name`): its editor opened, the value entered, Enter (check cells: Space) |
| `click` | `key` | clicks the button whose `objectName` is `key` |
| `save` | `path` (optional, case folder) | the window's *Save* (Ctrl+S); with `path`, *Save as…* onto that file |
| `screenshot` | `name`, `widget` (optional objectName) | saves the window (or that widget) as `<case>-<name>.png` in the screenshot folder (below); never compared (D10 of M4) |

Screenshots go to `$MDBRIDGE_GUI_SCREENSHOTS` (default: `gui-screenshots/` in the build tree);
CI uploads them as the `gui-test-screenshots-linux` / `-windows` artifacts. A failing case keeps
its case folder and names it.

| Folder | Content | Made by |
|---|---|---|
| `harness/` | `open_save_default`: the default project (`projects/default.mdatk.toml`) opened and saved unchanged, and the file MDBridge's window writes for that: MDATK v0.0.4.1's (its window's Open and Save without edits are `Project.load` + `Project.save`; checked once against the GTK window itself under Xvfb: same bytes) with the paths kept as written (AP-06, D8: `improvements_gui.load_as_written`; MDATK writes `out_dir` absolute) | `python3 -I -B tests/parity/make_gui_harness_cases.py <MDATK checkout> tests/data` |
| `roundtrip/`, `pages/`, `ap/`, `canvas/` | MDATK's GTK window driven by action scripts with MDBridge's window fixes (next section): open / save of the projects, edits of every page, the AP items, canvas drags, swaps, wheel zoom and map moves | `python -I -B tests/parity/make_gui_cases.py <MDATK checkout> tests/data` (step `gui`) |

### gui: the per-page, AP and canvas cases (milestone 5, W0-A)

`tests/parity/make_gui_cases.py` (step `gui` of `regenerate_all.py`) drives MDATK v0.0.4.1's GTK 3
window (`mdatk/gui/app.py`) headless under Xvfb through the same kind of action script and records
what it wrote and showed. The cases are in `gui/<topic>/` (`roundtrip`, `pages`, `ap`, `canvas`):
`<case>.json`, its input files in `<case>/in/`, the expected files in `<case>/out/` and MDATK's own
files, where they differ, in `<case>/mdatk/`. MDBridge's deliberate fixes of the window (AP-xx of
`docs/specs/m5-plan.md` section 4, `docs/DIFFERENCES.md` "Interface") are applied on the test side
by `tests/parity/improvements_gui.py`, so `expected` and `outcome` are MDBridge's behaviour (as in
the other generators, MDATK also runs with MDBridge's earlier improvements: `render_common.setup`,
unknown project keys kept by PRJ-09 / D9). The format is the one above, extended (a runner that
lacks something a case `needs` skips the case and names what is missing):

| Key | Content |
|---|---|
| `clock` | always set: every case starts at `2026-11-25T12:00:00` (MDATK reads the clock once per save) |
| `items` | ids the case covers (`AP-01`, `D8` ...; `[]`: plain parity with MDATK) |
| `mdatk_test` | the MDATK test (`tests/test_gui.py::...`) the case mirrors, or `""` |
| `needs` | the actions and keys of this section the case uses (`canvas`, `save_as`, `new`, `close`, `table`, `type`, `toggle`, `snapshot`, `start`, `settings`, `answers`; `open_failure` when an `open` is expected to fail), so a runner can skip what it does not have yet (`mdbridge_gui_tests` runs them all, below) |
| `canvas` | `{"width": W, "height": H}`: the figure canvas widget's size in pixels (device pixel ratio 1); the figure is drawn at `scale = min(W / video.width, H / video.height)`, `int(video.width * scale)` x `int(video.height * scale)` pixels at `video.dpi * scale`, centred |
| `settings` | path under tests/data of `$XDG_CONFIG_HOME/mdatk/settings.toml` before the start (`{case_dir}` in it: the case folder), or absent: no such file (D16, AP-17) |
| `start` | `{"project": name}`: the window is started with `{case_dir}/name` (as `gui PATH`), before the first action; absent: started without a project |
| `outcome` | MDBridge's expected outcome besides files (below) |
| `mdatk` | the keys of `expected` / `outcome` whose MDATK value differs, with MDATK's value (`expected` paths point into `<case>/mdatk/`); `{}`: MDATK did the same. Never compared |

`mdbridge_gui_tests` runs every case against the window (`app::MainWindow`, W1-A, W2-A); the figure
canvas gets the case's `canvas` size, and the canvas actions send their mouse events at the positions
of MDATK's record of the same action (`outcome.canvas`), resolving points from this canvas's hit areas
only for an action without a record (hit areas that differ from MDATK's by more than 3 px are
printed, not failed). Each case gets an
environment of its own (`XDG_CONFIG_HOME` and `APPDATA` = `{config_dir}`, `HOME` an empty folder,
the case folder current), the `settings` file is written first, the window is started with
`start.project` (`MainWindow::startWith`, as `gui PATH`), questions are answered from `answers`
(else the first choice) and recorded, file dialogs get `name` (Save as), `path` (Open) or
`answers.name` (Save of a project without a file), and message dialogs are recorded. After every
action the harness waits until the preview has redrawn (the generator's `flush`). Compared: the
files of `expected`; `log` (the window's log lines; the TOML parser's own words in `not valid TOML:
<words> (at line ...` are masked, tests/core/test_project_file.cpp); `outcome.questions`,
`dialogs`, `snapshots` (`values`, `error`, `page`, `density_status`, `modified`; numbers within
1e-9 relative), `stderr` (an absent key: none expected); and `outcome.settings` (comments line for
line, values as TOML, since arrays wrap by the length of the paths in them). A case without
`outcome` (`harness/`) is compared by its files only.

Placeholders, in expected files and in `outcome` texts: `{case_dir}` (the case folder, absolute)
and `{config_dir}` (`$XDG_CONFIG_HOME`, a folder outside it). Paths are written with `/`: on
Windows, where the case folder is Qt's `/` path joined with `\` and the window shows a path as
given or normalised (MDATK's `os.path` would write `\`), the folder matches with
either separator and the separators of the path after it become `/` (up to a space, quote, colon,
comma, semicolon or closing bracket); the settings file is compared by its parsed values, so a
TOML string's escaped `\\` is a path's `\` there. The environment of a run:
`XDG_CONFIG_HOME` = `{config_dir}`, `HOME` an empty folder, no cpptraj or conda installation found,
the current folder = the case folder.

Fields are addressed by project key (the widget's `objectName`): `analyses.contacts.cutoff`. List
items take their index: `groups.1.color`, `visual.representations.0.size`,
`inputs.trajectories.2`, `inputs.external.0.columns`, `analyses.rmsd.targets.0.selection`; the
slot editors of the layout page are `layout.panels.N.<field>`, fields `type`, `data`, `group_a`,
`group_b`, `more`, `title`, `options` (slot N, 0-based, whether or not it holds a panel yet).

Actions added or extended:

| Action | Keys | What the user does |
|---|---|---|
| `set` | `key`, `value` | as above; `value` is the project value: drop-downs choose the item whose project value (item data) is `value` (ids outside the list are allowed: AP-05), spin boxes and sliders take numbers (integers for integer fields), colour buttons a colour picked (`"#rrggbb"`, or any text: AP-03), check boxes `true` / `false`, list-valued fields a list of strings (`output.charts.formats`, `visual.structure_groups`), `visual.rotation` `[x, y, z]` (the three sliders), table cells their text, slot editors their choice (`.group_a` `""` = "(analysis default)") |
| `type` | `key`, `text` | replaces a text field's text with `text` and commits it (Enter / leaving the field): line edits, optional numbers (`inputs.stop`), comma lists, notes, table cells, the slot editors' `more`, `title`, `options` |
| `table` | `key`, `do`, `row`, `by` | list editors (`groups`, `inputs.trajectories`, `inputs.external`, `analyses.rmsd.targets`, `visual.representations`): `"add"` (a row with the editor's defaults at the end), `"remove"` (row `row`), `"move"` (row `row` by `by` = -1 / +1) |
| `click` | `key` | also `density.centre_on_model`, `density.clear`, `structure.reset_orientation`, `structure.reset_view` |
| `toggle` | `key`, `value` | a toggle button: `density.move_with_mouse` (map mode of the canvas) |
| `wait` | | as above; MDATK's side presses *Load real frame* when the project's topology exists (MDBridge loads it by itself) |
| `canvas.drag` | `from`, `path` or `to`, `shift` | press the left button at `from`, move, release. Points are `{"slot": N, "at": [fx, fy]}`: fractions of slot N's hit rectangle (the panel's drawn extent: axes, labels, colour bars; an empty slot's grid cell) from its top-left corner, y down. `path`: the motion events, each `[dx, dy]` pixels from the press point, released at the last; `to`: one motion to that point, released there. `shift`: Shift held throughout |
| `canvas.wheel` | `at`, `notches` | wheel turns over `at`, one event per entry: `1` up (zoom x 1.15), `-1` down |
| `canvas.double_click` | `at` | opens the slot's editor on the layout page |
| `new`, `close` | `answers` | *New project (defaults)*; closing the window |
| `open`, `save` | `answers` | as above (`open` may now fail: the error is in `outcome.dialogs`) |
| `save_as` | `name`, `answers` | *Save as…*, `name` typed in the file dialog (relative to the case folder, with or without `.mdatk.toml`: AP-11) |
| `snapshot` | `keys`, `slots`, `page`, `density_status`, `modified` | records the state into `outcome.snapshots` (below) |

`answers`: `"overwrite": "replace" | "cancel"` (a chosen file exists), `"unsaved": "save" |
"discard" | "cancel"` (unsaved changes before New, Open or closing: AP-13), `"name"` (typed when
*Save* asks for a file). A question without an answer gets the first choice. Before a canvas action
the canvas shows the figure of all earlier actions (pending redraws done; MDATK's 120 / 250 ms
timers do not matter). A press moved at least 4 pixels (|dx| + |dy|) is a drag; below the title
strip (top 12 %) of the structure panel it rotates (0.5° per pixel), in map mode it moves the map
(Shift: rotates it), anywhere else or with Shift it moves the panel onto the slot under the release.

`log` (the key above) holds MDBridge's expected log: the lines the window's log got from the start
on (MDATK's `... ready.` line left out; `{case_dir}` for the case folder); MDATK's own under
`mdatk.log` where they differ. `outcome`:

| Key | Content |
|---|---|
| `stderr` | lines printed on stderr (MDATK's start-up load warning, AP-12) |
| `dialogs` | `{"action": k, "kind": "error" or "info", "text", "detail"}` per message dialog (`action`: index into `actions`, -1 the start) |
| `questions` | `{"action": k, "kind": "overwrite" or "unsaved", "path", "answer"}` per question asked |
| `snapshots` | `{"action": k, "values": {key: value}}` (project values as the window would save them, JSON with the project file's types; `"error"` instead when the settings are invalid), with `slots`: per slot editor `{"frame", "type", "data", "choices": [[id, label], ...], "group_a", "group_b", "group_a_enabled", "group_b_enabled", "more", "title", "options"}`; `page`: the page shown (`input`, `groups`, `analyses`, `structure`, `layout`, `outputs`); `density_status`; `modified` |
| `canvas` | per canvas action `{"action": k, "slots": [[x, y, w, h], ...], "structure_slot": n, "events": [[kind, x, y], ...]}` (kinds `press`, `motion`, `release`, `scroll`, `double`): MDATK's hit rectangles and the pointer positions the action resolved to, for checking a runner's geometry |
| `uncaught` | `{"action": k, "error"}` per exception a handler let escape (MDATK: AP-09) |
| `settings` | path under tests/data of the expected settings file afterwards (absent: none) |

Files are compared as text; snapshot values that come from the canvas (rotation, zoom, density
placement) within 1e-9 relative, everything else exactly. The generator needs GTK 3 with its
introspection data, PyGObject and pycairo for its Python (pinned in
`tests/parity/requirements-gui.txt`, built against the system's GTK: `libgirepository-2.0-dev`,
`libcairo2-dev`, `gir1.2-gtk-3.0`) and Xvfb (started by the generator when `DISPLAY` is unset). It
runs every case twice in child processes, as MDATK and with `improvements_gui.py`:

```bash
python -I -B tests/parity/make_gui_cases.py /tmp/mdatk tests/data
```

## numerics: SciPy's numerics (milestone 6, W0-C)

`tests/parity/make_numerics_cases.py` runs SciPy 1.18.1 and numpy 2.5.5 alone (no MDATK), for
`src/numerics` (`tests/numerics/test_fft.cpp`, `test_ndimage.cpp`, `test_optimize.cpp`,
`test_pairs.cpp`; docs/DIFFERENCES.md "Numerics"):

| File | Content |
|---|---|
| `numerics_cases.json` | `fft`: SHA-256 of `rfftn` / `irfftn` of formula arrays (and of MDATK's correlation product and a raw half spectrum), every value for small shapes (C++: bit for bit, any thread count); `next_fast_len`; `gaussian` (per-axis sigma, modes "reflect" / "constant") and `map_coords` (order 1, mode "constant", cval 0 / 1.5 / NaN) of float32 arrays; `powell` and `brent` runs of named objectives with every iterate (C++: within 1e-9; they agree bit for bit); `least_squares` runs (dense and sparse, x_scale, diff_step, max_nfev, LSMR options) with every iterate and `ulp_spread`, SciPy's own movement when its residuals move by an ulp (C++: within 1e-9 or four times that); `chain_sparsity`; `permutations` (`RandomState(0).permutation`) and `groups` (`group_columns`); `lsmr`; `roots` (np.roots' real roots); `pairs` (`cKDTree.query_pairs`, sorted) |

The arrays come from `value(i) = ((i * 2654435761) mod 2^32) / 2^32 - 0.5` (exact), the objectives
from `math` and fixed-order arithmetic, so the C++ test rebuilds them. `least_squares`, `lsmr` and
`roots` go through BLAS and LAPACK: the generator pins their kernels as `make_analysis_cases.py` does,
and `regenerate_all.py` compares the file as JSON within a relative 1e-12 where its bytes differ:

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_numerics_cases.py tests/data/numerics
```

## Checking that the data is reproducible

`tests/parity/regenerate_all.py` runs every generator above (and `tools/gen_*.py`) in this order
into a fresh tree and compares it with the committed files byte for byte; CI's `regenerate` job
runs it on every push with Python 3.13, MDATK `v0.0.4.1` and the environment pinned in
`tests/parity/requirements-regenerate.txt`:

```bash
python3 tests/parity/regenerate_all.py --mdatk /tmp/mdatk --out /tmp/fresh --python /tmp/regen/bin/python
```

The step `gui` (make_gui_cases.py) needs PyGObject and pycairo too (`tests/parity/requirements-gui.txt`,
built against the system's GTK 3; `--gui-python` names another interpreter for it) and Xvfb.

`tools/gen_data_tables.py` runs first, on the committed `data/sources/` extracts copied into the
fresh tree (they are inputs, cut from downloads by `tools/extract_ccd.py` and
`tools/extract_openmm.py`, and not regenerated), so `src/core/include/mdbridge/core/detail/data_tables_independent.inc` is
checked like the other tables, and the patched generators read the same tables.

It handles the special cases noted above (`--fixed-clock`, `manifest.txt` and the two extra vector
files deleted, the loader before the export, MDATK's lock files and MDAnalysis' hidden
`.*_offsets.npz` / `.lock` files removed) and masks the patch number of Python's version
(`3.13.N`) where a header records it. Not compared, because no generator writes them: the
hand-made project files, `text/tomlio_golden.toml` and `topology/golden_model.json`. Some files may
differ by CPU: `rigid/rigid_cases.json` and `analysis/spatial_extra_cases.json` (its `kabsch` and `rings`
entries) come from OpenBLAS / LAPACK, whose kernels are chosen per CPU (the committed files were
written with the SkylakeX ones; Haswell / Zen kernels differ in the last bits), so where their bytes
differ they are compared as JSON with numbers (and `float.hex()` strings) within a relative 1e-12.
The same holds for `analysis/demo_dry_cases.json`, `demo_dry_variants.json`, `micro_cases.json` and
`systems_cases.json` despite the generator's pinning (CI run 56 moved RMSDs with a fit selection
by 12 to 33 ulps), where an array's `sha256` is ignored when its values agree within the tolerance;
not `analysis/spatial_cases.json`, whose degenerate `kabsch` inputs get a different rotation from
other kernels. The JSON files under `render/figures/`, `render/charts/` and `render/molecule/` are
compared the same way (section "Rendering").

other kernels. Likewise `render/molecule_geometry/styles.json`, `scenes.json` and `ndimage.json`
(pinned kernels; a small matmul and `np.exp`), and `numerics/numerics_cases.json` (least_squares'
BLAS sums and LAPACK's SVD, QR and eigenvalues).
The random rotations of `io/pdb/pdb_cases.json` are written out in `make_pdb_cases.py` for the same
reason. `--list` shows the steps; `--only STEP ...` reruns some of them.

`render/ticks/tick_cases.json` holds what Matplotlib 3.11.2 (MDATK v0.0.4.1's plotting library,
rcParams at their defaults as MDATK leaves them) and MDATK's `thin_ticks` compute, for
`src/render/core/ticks.cpp` (`tests/render/test_ticks.cpp`, compared bit for bit and string for
string): `_nonsingular`; `MaxNLocator` / `AutoLocator` `tick_values` and `view_limits` over about 230
intervals (magnitudes 1e-300 to 1e300, widths down to a few ulps, large offsets, non-finite ends)
with fixed and 'auto' bins, steps, integer, symmetric, prune and min_n_ticks, and the invalid-steps
messages; `FixedLocator` with nbins; `ScalarFormatter` labels, offset texts and formats (automatic,
disabled and fixed offsets, offset thresholds, power limits, scientific off, mathtext, ASCII minus)
and `format_data`; real Axes (figure sizes, dpi, rectangles, tick label sizes 0.5 to 20 pt; limits
set, autoscaled from plotted data with the 5 % margins, sticky edges or tight; MDATK's
`integer_axis`) with the axis length in points, `get_tick_space`, view interval, located and drawn
ticks, labels and offset text; `set_xlim` and imshow extents; MDATK's heatmap colour bar (fraction
0.04, pad 0.02, labels at 0.75 x font, integer ticks for counts) with its limits and ticks; and
`thin_ticks`, where `max_ticks` < 1 is MDBridge's improvement (docs/DIFFERENCES.md). The generator
pins numpy's dispatch to its x86-64 baseline; the output is reproducible:

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/make_tick_cases.py /tmp/mdatk tests/data/render/ticks
```

Never put research structures or trajectories of real projects here.
