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

### PDB, GRO and the loader (src/io/pdb_topology.cpp, loader.cpp, trajectory.cpp)

Expected data: tests/parity/improvements_readers_pdb.py, called by make_loader_cases.py,
make_export_cases.py, make_pipeline_cases.py and cli_parity.py. MDAnalysis runs on *equivalent* rewrites
of the files (a rewrite MDAnalysis reads the way MDBridge reads the original), and the errors MDBridge
words itself replace MDAnalysis'. Errors marked RM follow the reader-message rule
(docs/IMPROVEMENTS.md): they name the file and the line, say what was expected, are ValueErrors and are
never rewrapped in MDAnalysis' "Failed to construct topology ..." / ChainReader texts.

- **Reader warnings (LD-03).** MDATK loads with every Python warning ignored, so defaulted cells,
  elements, masses, serials and skipped records went unnoticed. MDBridge's readers report them through
  `io::readerWarning()` (src/io/include/mdbridge/io/warnings.hpp): one warning per kind and file, with
  a count and the first line. `loadSystem()` collects them in `core::System::warnings`
  (`readTopologyFile()`: `TopologyFile::warnings`); `info` prints each as `WARNING: <message>` on
  stderr, `analyze` (Pipeline::system()) and `export-pdb` log them as `WARNING: <message>` (hidden by
  `-q`, like their other log lines); a warning issued outside a load goes to `core::logWarning()`.
  *Tested:* test_loader.cpp (ReaderWarnings, SystemWarnings, WarningScopes); cli_cases.json (`info_pdb`,
  `export_pdb_trajectory`, the `inputs.dt_ps` cases).
- **Hybrid-36 CONECT serials (RD-B01).** CONECT serials are decoded like ATOM serials (`A0000` =
  100000), so files with more than 99,999 atoms keep their bonds (MDAnalysis: `invalid literal for
  int() ...`). *Tested:* io/loader/files/hy36.pdb.
- **Malformed CONECT records (RD-B09, RD-B11).** A record MDAnalysis calls corrupt (its length
  rule), or with a field that is no serial, is skipped alone with a warning (MDAnalysis drops every
  CONECT bond of the file, or fails); a blank partner field is no partner. The per-pair "CONECT record
  to TER entry" and the single-entry warnings become one counted warning each. *Tested:*
  conect_one_bad.pdb, conect_bad_int.pdb, conect_variants.pdb.
- **Bad hybrid-36 atom serial (RD-B02), RM:** `<file> line 1: invalid hybrid-36 atom serial 'Aa000'`
  (MDAnalysis: KeyError `'a'` or UnboundLocalError). *Tested:* hy36_bad.pdb.
- **One atom-line rule (RD-B04).** Both passes take a line as an atom when columns 1-6 are `ATOM  ` or
  `HETATM` (MDAnalysis' topology pass strips the line first and tests a prefix, its coordinate pass does
  not, so an indented line gave `Inconsistency in file ...`). Other lines that start with ATOM/HETATM
  (indented, `ATOMS`, `ATOM 123456`) are ignored with one warning. *Tested:* leading_space.pdb.
- **Model and record errors (RD-B05), RM:** `<file>: no ATOM or HETATM records`, `<file> line 5: ATOM
  or HETATM record after END (line 4)`, `<file>: model 2 has 6 atoms, model 1 has 7 (line 11)`
  (MDAnalysis: numpy IndexErrors, or the Inconsistency error when the model is read). Every model's
  atom count is checked when the file is opened, so a file with inconsistent models fails at load
  instead of in the middle of an analysis. *Tested:* remarks_only.pdb, endmdl_first.pdb,
  models_mismatch.pdb, models_more.pdb.
- **Empty files (RD-B06), RM:** `<file>: file is empty` for PDB and GRO (MDAnalysis: `Compressed file
  ended before the end-of-stream marker was reached`). *Tested:* empty.pdb.
- **Text that is not UTF-8 (RD-B07).** Bytes that are not part of valid UTF-8 are read as Latin-1, with
  one warning naming the first line (MDAnalysis rejects the file). HEADER / COMPND / REMARK texts are
  cut by characters of the decoded line (MDAnalysis cuts bytes, so a cut through a character failed).
  *Tested:* latin1_name.pdb, invalid_utf8.pdb, header_cut_utf8.pdb.
- **Bad coordinates (RD-B08), RM:** `<file> line 3: x coordinate '   1.0x0' is not a number` (with
  `(model k)` in multi-model files; GRO likewise with `(frame k)`). *Tested:* bad_coordinate.pdb.
- **Cells do not depend on the read order (RD-B03).** Every model's cell is fixed when the file is
  opened: its last CRYST1, else the cell of the nearest earlier model in the file (MDAnalysis: of the
  model read before, so random access gave other cells). A CRYST1 before the first MODEL counts as
  earlier, so the common layout "one CRYST1, then MODEL records" (cpptraj) keeps its box (MDAnalysis
  never reads it). CRYST1 placeholder / unreadable-record warnings come once, with the line.
  *Tested:* models_some_cryst.pdb and models_one_cryst.pdb (in and out of order), export cases.
- **Partial element column (RD-B10).** Atoms whose columns 77-78 are blank, in a file where others are
  filled, get type, element and mass guessed from the name (with the residue rule), as in a file without
  the column (MDAnalysis: type `""`, mass 0, silently). One warning counts the atoms whose mass stays
  unknown (types that are no element). *Tested:* elements.pdb, segids.pdb, conect_duplicate_serial.pdb
  (dumps and export cases).
- **Per-atom warnings (RD-B11)** (unknown elements, missing resids, made-up serials) are one counted
  warning each.
- **No time information (RD-B13).** PDB and multi-frame GRO files keep 1 ps between frames, and warn
  that the file has no time information (unless inputs.dt_ps is given).
- **GRO field width (RD-G01).** The width is the distance between the decimal points of y and z, and the
  x field's point must agree (GROMACS' own rule); otherwise `<file> line 3: cannot determine the
  coordinate field width ...` (MDAnalysis takes the first point from column 26 and read 11-character
  fields as 1-character ones, silently). *Tested:* wide.gro.
- **GRO structure errors (RD-G02, RD-G03, RD-G04), RM:** `<file>: missing box line after atom 7 (line
  9)`, `<file>: the header (line 2) says 7 atoms but the file has 4 atom lines`, `<file>: no atoms (the
  atom count on line 2 is 0)`,
  `<file> line 2: atom count ' seven' is not an integer` (MDAnalysis: UnboundLocalError, IndexError,
  `invalid literal for int()`). *Tested:* no_box.gro, too_few.gro, bad_count.gro.
- **GRO box line (RD-G05).** When splitting on blanks fails, the values are read as 10-character fields
  (GROMACS' %10.5f, which run together beyond 9999 nm), else as run-together numbers with five decimals,
  **signs kept** (MDAnalysis' fallback pattern dropped them); otherwise `<file> line 10: GRO box line
  '   3.0   4.0' has neither 3 nor 9 values`. *Tested:* neg_box.gro, packed_box.gro, bad_box.gro.
- **GRO boxes that are no cell (RD-G06)** warn: `<file>: 1 box is zero; no periodic box (line 10)`,
  `... box has vectors that do not form a cell ...`. *Tested:* zero_box.gro, bad_cell.gro.
- **Multi-frame GRO (RD-G07).** Every frame of a GRO file is read (GROMACS trajectories are such frames
  one after another; MDAnalysis reads the first only). A later frame with another atom count is an
  error; text after the last frame that is no complete frame is ignored with a warning. *Tested:*
  frames.gro.
- **Mixed chains (LD-01).** Each trajectory gets the format of MDATK's extension table for its own
  extension (MDATK applied the table only when all files shared one extension: `.rst7` + `.inpcrd`
  failed, `.crd` became CHARMM CRD). *Tested:* loader_cases.json (time.rst7 + time.inpcrd).
- **GROMACS .top (LD-02).** `.top` stays parm7 (MDATK's mapping), but a file whose first line is not
  `%VE` gives `<file>: not an Amber parm7 topology (GROMACS .top files are not supported)`. *Tested:*
  test_loader.cpp GromacsTopFile.
- **Time step validation (LD-04).** validate() reports `inputs.dt_ps must be positive (0: the
  trajectory's own time step).` for NaN, infinite and negative values, and loadSystem() rejects them
  (MDATK passed them on: NaN or decreasing times). *Tested:* test_project_model.cpp ValidateTimeStep,
  test_loader.cpp InvalidTimeStep, loader_cases.json.
- **Given vs stored time step (LD-05).** When inputs.dt_ps is given for a file that keeps its own times
  (NetCDF with a time variable) and the steps differ: `inputs.dt_ps (2.5 ps) differs from the time step
  in <file> (10.0 ps); the file's times are used`. *Tested:* SystemWarnings, pipeline_cases.json,
  cli_cases.json.
- **No coordinates (LD-06):** `<file> has no coordinates; give a trajectory with --traj` (MDATK:
  AttributeError `No trajectory loaded into Universe`, a traceback in `info`). *Tested:*
  loader_cases.json, cli_cases.json `info_parm7_alone`.

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

## Command line

- **No tracebacks.** Where MDATK crashes with a Python traceback (`info` with a parm7 file alone,
  corrupt gzip maps, map dimensions too large to allocate), MDBridge prints `error: <message>` and exits
  with 1. *Tested:* cli_cases.json (`mdatk_crash` cases), tests/cli/test_commands.cpp.
- Six gzip error texts of `density info` are worded as in mrc.hpp (compared by exit status only).
- JSON presets: `init` with a `.json` preset reports that MDBridge uses TOML project files.

## Planned (from the improvement register; add each entry here when it lands)

Periodic distances and minimum images for every cell (AN-01–AN-04); degenerate RMSD (AN-16); base
pairs and stacking on wrapped trajectories (AN-11, AN-12); setup-order independence (AN-13, SEL-17);
`prop` operand reversal (SEL-01, SEL-02); per-atom element fallback (SEL-13); reader fixes (RD-*,
LD-*); validation of every cutoff (PRJ-01); readable messages in place of Python's internal texts
(SEL-03/04/08, RD-*, CLI-01/03/04).
