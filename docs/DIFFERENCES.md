# Deliberate differences from MDATK v0.0.4.1 and MDAnalysis 2.11.0

MDBridge reproduces MDATK v0.0.4.1, and through it MDAnalysis 2.11.0, wherever users rely on the
behaviour. Where the reference has a bug, a crash, a misleading message or a numerically fragile
result, MDBridge does better instead (docs/ARCHITECTURE.md, "Definition of done"). Each entry below
says what differs, why, and how it is tested. The parity tests expect MDBridge's behaviour. Their
expected data is adjusted on the test side; generated files are never edited by hand.

## Data and licensing

- **Independent data tables.** Residue names of the `protein`, `nucleic` and `water` macros, the
  atom-name sets of `backbone`, `nucleicbackbone`, `nucleicbase` and `nucleicsugar`, element symbols,
  masses, the force-field name → element table and the PDB writer's ion names come from public
  sources (wwPDB CCD, PDB v3.3, IUPAC-IUB, CIAAW 2024, NUBASE2020, the Amber and GROMACS manuals,
  OpenMM 8.1.2; data/sources/ATTRIBUTION.md), not from MDAnalysis. Names that would collide with a
  standard-residue atom name of another element are dropped (tools/gen_data_tables.py). Masses are
  found under the element symbol and its upper-case spelling, so `LI` gets lithium's mass where
  MDAnalysis gives 0. *Why:* licensing (LGPL tables are not carried). *Tested:*
  tests/parity/independent_tables.py patches MDAnalysis' data before every generator that depends on
  it, and checks every CCD reference name.
- **Residue-scoped ion rule.** An atom whose (residue, atom name) is a one-atom ion or a known
  mercury/chlorine site (`CA`/`CA` → Ca, `CD`/`CD` → Cd, `EMC`/`HG` → Hg, `SOD`/`SOD` → Na, ...;
  table `kIndResidueAtomElements`) gets that element's type and mass. The PDB writer starts its name
  in column 13 when the element has two letters. MDAnalysis decides by name alone, so a calcium ion
  named CA is carbon there. *Tested:* io/loader/files/ions.pdb and ions.gro; residue_rule_universe and
  residue_rule_pdb_text post-process MDAnalysis' results.
- **No Barnaba.** Barnaba is GPL. With `analyses.basepairs.engine = "barnaba"`, MDBridge logs
  `WARNING: basepairs: Barnaba is not available in MDBridge; using the builtin detector` and runs its
  builtin detector (meta `engine = "builtin"`, no `lw_first_frame`, `noncanonical` has no effect).
  (MDATK without Barnaba installed stops at setup with an installation hint.) The params, and so the
  cache key, keep the project's values; the parity cases leave the barnaba engine out. *Tested:*
  tests/analysis/test_basepairs.cpp.

## Analysis cache

- **Own cache files.** MDBridge writes its own archives and manifest in the same cache folder
  (`MDBridge_cache.json`, `<name>-<key>.mdbridge.npz` / `<name>-<key>.mdbridge.json`;
  src/analysis/include/mdbridge/analysis/cache.hpp), so they never collide with MDATK's
  (`MDATK_cache.json`, `<name>-<key>.npz` / `.json`). It does not read MDATK's cache: results are
  recomputed once. The archives are numpy `.npz` files; data keys, cache keys and the JSON
  descriptions keep MDATK's meaning (the manifest's version field is `mdbridge_version`). *Why:*
  the stored format is MDBridge's to evolve, and a project shared between MDATK and MDBridge keeps
  both caches intact. *Tested:* tests/analysis/test_pipeline.cpp; cli_parity.py compares cache files
  by meaning.

## Frame-parallel runner

MDATK processes the frames of a pass one after the other. MDBridge processes them on several
threads (`MDBRIDGE_THREADS`, default: the processor's threads, at most 16; 1 runs sequentially).
The trajectory is read in order on one thread, analyses that need the frames in order see them in
order, and results are retired in frame order, so arrays, metadata, logs, progress calls, frame
times, the frame left current and the first error are identical to a sequential run, bit for bit.
*Why:* speed; MDATK's single-threaded loop is the bottleneck on long trajectories. *Tested:*
tests/analysis/test_parallel.cpp (every parity case with 1, 2, 4 and 8 threads, errors, reader
failures, cancellation); ctest runs every analysis parity test with 4 threads and again with 1
(`<test>/threads1`).

## Order effects of the neighbour search

MDATK's pairs without a periodic box come in the traversal order of scipy's k-d tree, which depends on
how scipy was built. MDBridge returns them sorted by (i, j). The pair set and every distance are
identical. With a periodic box (`use_pbc`) MDBridge returns MDAnalysis' order, so nothing below
applies there. *Why:* a defined order that does not depend on the build of a library. Visible
consequences without a box:
- energy sums over pairs (elec, vdw, net force) can differ in the last bits;
- H-bond hydrogen guess: a hydrogen exactly as far from two donors goes to the donor with the lowest
  index (MDATK: the first in tree order);
- stacking map entries within a frame are in (i, j) order;
- `atoms` analysis: distance-guessed bond rows are sorted by (i, j).

*Tested:* test_contacts_hbonds.cpp (tie case), test_mg_stacking.cpp (sorted comparison),
test_structure.cpp, test_energy.cpp (tolerances).

## Readers and topologies

- **Self-bonds** (`CONECT 1 1`, kept by MDAnalysis) are dropped; the topology model has none.
- **`info` residue-name ranking:** equal counts are ordered by name (MDATK: numpy's argsort order,
  which depends on the CPU's sort kernel).
- **parm7 end of file:** the StopIteration carries "<file>: the parm7 file ends unexpectedly" instead
  of an empty message.
- **Periodic pair search where MDAnalysis is undefined:** a point more than 2^31 cells away is not
  moved along that vector (MDAnalysis overflows a 32-bit integer); a triclinic box vector with an
  infinite component is not periodic; a point with a NaN or infinite coordinate is in no pair for
  triclinic cells too (MDAnalysis reports it at distance 0 from every point, or crashes).
- **Distance selections** (`around`, `sphzone`, ...) compute distances directly. The grid-search
  artefact of MDAnalysis with cutoff 0 and no box (atoms on opposite faces reported at distance 0) does
  not occur. Pairs at exactly the cutoff may differ.
- **Unsupported selection features** are rejected with an error instead of being evaluated: `smarts`,
  `cyzone`/`cylayer` with coordinates, `type_index`, and bond / angle objects of topologies with bonds
  (marked `unsupported` in atom_selection_cases.json).
- **RMSD:** numpy's RuntimeWarnings for non-finite or empty input are not reproduced (MDATK ignores
  them).

## Density maps (src/io/include/mdbridge/io/mrc.hpp)

- Truncated or corrupt gzip data raises MapError naming the file ("file too short for an MRC header",
  "file is shorter than its header says", "Error -3 while decompressing data: ...", "CRC check failed
  ...", ...) where Python raises EOFError / zlib.error / BadGzipFile without the path. A member's
  trailer is checked only when the map needs data beyond it.
- Header sizes beyond the file raise "file is shorter than its header says" without allocating
  (Python: MemoryError or OverflowError).
- min/max of maps whose extreme is ±0 may differ in the sign of the zero; a NaN result keeps the first
  NaN's bits.
- Simulated maps use the C library's exp(), as numpy does without SIMD kernels.

## Project files

- An integer field given an integral number beyond 64 bits is an "expected an integer" error (Python
  keeps a big integer). NaN or infinity in an integer field gives "cannot convert float NaN to integer"
  / "cannot convert float infinity to integer" as a ProjectError.
- A planner node's `files` value that is not a string gives "<where>.files.<port>: expected a file
  path, got <repr>" (MDATK crashes).
- **Windows paths** (MDATK never ran on Windows; src/core/include/mdbridge/core/ntpath.hpp): `/` and `\`
  both separate; drive, UNC and root-relative paths are absolute (root-relative takes the project
  folder's drive); `C:x` resolves against the project folder when that is on drive C, otherwise against
  `C:/`; `~` expands to %USERPROFILE% (else %HOMEDRIVE%%HOMEPATH%); results are normalised and written
  with `/` (`C:\data\..\a.pdb` → `C:/a.pdb`); device paths (`\\?\`, `\\.\`) are kept as written;
  comparison is case-sensitive. On Linux and macOS paths follow posixpath exactly, as MDATK does.
  *Tested:* tests/core/test_ntpath.cpp on every platform.
- **PRJ-01: every numeric analysis parameter is validated.** MDATK checks only `contacts.cutoff <= 0`
  and Mg inner > outer, so NaN passes everything, and negative or NaN H-bond, Mg, stacking and
  energy cutoffs, `energy.dielectric = 0` (infinite energies), angle cutoffs outside 0-180 and
  `proximity_search_radius <= 0` give silent zero, NaN or inf series. validate() now requires every
  cutoff, radius and tolerance to be a finite number > 0 ("H-bond distance cutoff must be
  positive." / "... must be finite."), the H-bond and stacking angle cutoffs to lie in [0, 180]
  ("... must be between 0 and 180 degrees."), and keeps the Mg inner <= outer check for valid
  cutoffs; the checks stand where MDATK's two did (core::analysisParameterProblems). *Why:* silent
  nonsense results. *Tested:* project-model/validate_cases.json (new cases), cli_cases.json
  (info_invalid_parameters, analyze_invalid_parameters); tests/parity/improvements_cli.py applies
  the same checks to MDATK's validate().
- **PRJ-02: a clamped RMSD reference frame is reported.** `analyses.rmsd.reference_frame` beyond
  the analysed frames (or negative) is still clamped to the last (first) frame, now with
  `WARNING: RMSD: reference frame 999 is beyond the 24 analysed frames; using the last` (`... -3 is
  negative; using the first`). Results and cache keys are unchanged. *Tested:*
  analysis/demo_dry_variants.json (rmsd_reference_999 / -3), cli_cases.json
  (analyze_reference_beyond / _negative).
- **PRJ-03: RMSD targets with the same data key are a validation problem.** Labels `A B` and `A_B`
  both give the key `rmsd/A_B`; MDATK silently overwrites the first series. validate() reports
  `RMSD targets 1 and 2 give the same data key 'A_B'; rename one.` *Tested:* validate_cases.json
  (rmsd keys), cli_cases.json.
- **PRJ-04 / PRJ-05: index and frame texts.** `--set groups.--1.name=x` reports
  `groups.--1.name: no list element '--1'` and `--frames 1:x` reports
  `--frames: '1:x' is not a frame number or range` (MDATK: `invalid literal for int() with base
  10: ...`). *Tested:* cli_cases.json (init_set_index_dashes, init_frames_bad, init_frames_bad_part),
  test_commands.cpp.
- **PRJ-07: a Project property given as a key gets a hint.** `--set out_dir=x` (or a top-level
  `out_dir` / `cache_dir` / `workflow_dir` key in a project file) warns
  `project: ignoring unknown key 'out_dir' (did you mean output.out_dir?)`
  (`output.cache_dir`, `workflow.out_dir`). *Tested:* from_dict_cases.json (properties as keys),
  cli_cases.json (init_set_property, init_set_cache_dir).

## Command line

- **No tracebacks.** Where MDATK crashes with a Python traceback (`info` with a parm7 file alone,
  corrupt gzip maps, map dimensions too large to allocate), MDBridge prints `error: <message>` and exits
  with 1. *Tested:* cli_cases.json (`mdatk_crash` cases), tests/cli/test_commands.cpp.
- **CLI-01: KeyErrors without Python's quotes.** `error: unknown analysis 'nope' (known: ...)` and
  `WARNING: skipping data source 'bogus@Nope': unknown analysis ...` (MDATK prints str(KeyError),
  the message's repr in quotes). The step log keeps `KeyError: '...'` (a project-file convention).
  *Tested:* cli_cases.json (analyze_unknown, analyze_layout), pipeline_cases.json.
- **CLI-02: external data is read before the frame loop.** A missing or unreadable
  `inputs.external` file fails before the analyses run (same message, `external data: ...`), so no
  run time is wasted and nothing is cached for the failed run; MDATK finds it after the pass.
  *Tested:* cli_cases.json (analyze_external_bad / _missing), pipeline_cases.json (errors).
- **CLI-03: `stop=end`** in `frame selection start=30 stop=end step=1 selects no frames ...`
  (MDATK: Python's `None`). *Tested:* cli_cases.json, io/loader/loader_cases.json,
  pipeline_cases.json.
- **CLI-04: map errors always name the map.** `density info` prints MapError texts with the path,
  gzip errors included (MDATK: `Not a gzipped file (b'XY')`, `CRC check failed ...` without it), and
  `{path}: cannot determine the byte order of the map header`; a MODE or cell angles that only make
  sense in the other byte order give `... byte order of the map header (mode 33554432)` / `(cell
  angles [...])` instead of `unsupported MRC mode 33554432` / `non-orthogonal cell angles [6.47e-41,
  ...]`. *Tested:* io/maps/read_cases.json, cli_cases.json (density_info_* compared in full).
- **CLI-08: PDB export checks every frame first.** A frame with a coordinate beyond the PDB limits
  stops the export before any file is written: `PDB files must have coordinate values between
  -1000.000 and 9999.999 Angstroem: frame 2 is outside this range; no files were written.` (MDATK
  aborts at that frame and leaves the earlier files.) *Tested:* io/pdb/export/export_cases.json,
  cli_cases.json (export_coordinate_limits).
- **A single-value option given twice takes the last value**, as argparse does (`--top`, `-o`,
  `--cache-dir`, `--frames`, `--project`, `--cpptraj`, `--conda-env`, `--only`, `--frame-list`,
  `--times`, `--selection`; a repeated `--traj` replaces the earlier list). This is MDATK's
  behaviour, listed because CLI11 rejects repeats by default. *Tested:* cli_cases.json
  (init_options_twice, analyze_options_twice, export_options_twice).
- JSON presets: `init` with a `.json` preset reports that MDBridge uses TOML project files.

## Planned (from the improvement register; add each entry here when it lands)

Periodic distances and minimum images for every cell (AN-01–AN-04); degenerate RMSD (AN-16); base
pairs and stacking on wrapped trajectories (AN-11, AN-12); setup-order independence (AN-13, SEL-17);
`prop` operand reversal (SEL-01, SEL-02); per-atom element fallback (SEL-13); reader fixes (RD-*,
LD-*); readable messages in place of Python's internal texts (SEL-03/04/08, RD-*).
