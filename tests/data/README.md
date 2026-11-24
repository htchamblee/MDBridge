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
`make_selections_cases.py`, `make_parm7_cases.py`, `make_loader_cases.py`, `make_export_cases.py` and
`cli_parity.py` (its MDATK side); the others do not touch them. `tools/gen_selection_tables.py` and
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
ANISOU, line ends, limits); `tail` — `builder._tail`. Floats are `float.hex()` strings.

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
width" errors) and headings (`None`, empty, long), `pipeline_resources` of about 80 projects made
with `Project.from_dict` (cleaning, video, renderers, base-pair engines, secondary panels with
every kind of `source` option) and `project_resources` with stand-in plans (`generating` is what
MDBridge's caller passes: the resources of the steps whose tool generates a script, in plan order).

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
milestone-2 commands, for `src/app/commands/` (`tests/cli/test_commands.cpp`, which runs the same
command lines in-process): `init` / `init-params` (defaults, every common option, `--set` with list
indices, JSON values, free-form tables, Project properties and every error, projects and a parameter
log as `--project`, existing files, JSON presets), `info` (the demo system with NetCDF / XTC / DCD
trajectories, PDB and GRO topologies, residue groups that fail in every way, missing files),
`export-pdb` (frame lists, times, groups and selections, time stamps, `--quiet` / `--verbose`,
layout warnings, a project's `[output.pdb]`, every error) and `density info` of every map of
`io/maps/files`. Each case names the command line, the inputs (files of `tests/data` copied into an
empty work folder, project files written there) and MDATK's exit status, stdout, stderr and the files
it left (text, or SHA-256 and size beyond 16 KiB). MDATK runs in-process with the independent
tables applied (first section). Normalised: the work folder is `{work}`, `saved` time stamps
`{time}`, file-name time stamps `{stamp}`; `posix_only` marks cases with absolute paths, OS error
texts or project paths (skipped on Windows); `exit_only` cases compare the exit status alone (usage
errors, and the six gzip errors `src/io/include/mdbridge/io/mrc.hpp` words differently); where MDATK
crashes with a traceback (`stderr` null, `mdatk_crash` the exception), MDBridge must print
`error: ...` and exit with 1. The same script compares a build with MDATK directly (CI's `parity`
job, against the AppImage):

```bash
/opt/mdbridge-venv/bin/python -I -B tests/parity/cli_parity.py --mdatk /tmp/mdatk --write tests/data/cli/cli_cases.json
/opt/mdbridge-venv/bin/python -I -B tests/parity/cli_parity.py --mdatk /tmp/mdatk --mdbridge build/tests/mdbridge_commands_cli
```

## Checking that the data is reproducible

`tests/parity/regenerate_all.py` runs every generator above (and `tools/gen_*.py`) in this order
into a fresh tree and compares it with the committed files byte for byte; CI's `regenerate` job
runs it on every push with Python 3.13, MDATK `v0.0.4.1` and the environment pinned in
`tests/parity/requirements-regenerate.txt`:

```bash
python3 tests/parity/regenerate_all.py --mdatk /tmp/mdatk --out /tmp/fresh --python /tmp/regen/bin/python
```

`tools/gen_data_tables.py` runs first, on the committed `data/sources/` extracts copied into the
fresh tree (they are inputs, cut from downloads by `tools/extract_ccd.py` and
`tools/extract_openmm.py`, and not regenerated), so `src/core/include/mdbridge/core/detail/data_tables_independent.inc` is
checked like the other tables, and the patched generators read the same tables.

It handles the special cases noted above (`--fixed-clock`, `manifest.txt` and the two extra vector
files deleted, the loader before the export, MDATK's lock files and MDAnalysis' hidden
`.*_offsets.npz` / `.lock` files removed) and masks the patch number of Python's version
(`3.13.N`) where a header records it. Not compared, because no generator writes them: the
hand-made project files, `text/tomlio_golden.toml` and `topology/golden_model.json`. One file may
differ by CPU: `rigid/rigid_cases.json` comes from OpenBLAS / LAPACK, whose kernels are chosen per
CPU (the committed file was written with the SkylakeX ones; Haswell / Zen kernels differ in the
last bits), so where its bytes differ it is compared as JSON with numbers within a relative 1e-12.
The random rotations of `io/pdb/pdb_cases.json` are written out in `make_pdb_cases.py` for the same
reason. `--list` shows the steps; `--only STEP ...` reruns some of them.

Never put research structures or trajectories of real projects here.
