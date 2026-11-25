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
- **Analysis versions.** An analysis' version (part of its cache key) is bumped whenever an
  improvement below can change its results for real inputs, so caches written before the change are
  not reused. MDBridge has not been released, so one bump (MDATK's version + 1) covers every
  improvement of this release: energy 2 (AN-07, AN-08, RD-P12's float64 charges, SEL-13), hbonds 2
  (AN-09, AN-10, SEL-13), rmsd 2, proximity 2 (AN-15, SEL-13), atoms 3 (AN-13, AN-14, SEL-13),
  contacts 2, mg 2, basepairs 2 and trace 5 (SEL-13). The test side gives MDATK's analyses these
  versions from one table, tests/parity/analysis_versions.py, which regenerate_all.py checks against
  the C++ descriptors.

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

## Amber, DCD and GROMACS readers (docs/IMPROVEMENTS.md 2.1, 2.4-2.6)

Expected data: the generators record MDAnalysis; tests/parity/improvements_readers_amber.py turns
it into MDBridge's behaviour (per item, below). Error texts follow the reader-message rule: they name
the file and, where known, the line, columns, section or frame. Warnings go through the reader-warning
channel (LD-03, `io::readerWarning`); tests/data/io/trajectories/cases.json lists each case's
warnings (`warnings`, computed by improvements_readers_amber.trajectory_warnings(); test_trajectory.cpp
expects exactly these, and none for the other cases).

parm7 topologies (src/io/parm7.cpp; tested by io/parm7/topology*, test_parm7.cpp
`ImprovedTopologyRules`):
- **RD-P01** An unknown section right after POINTERS is skipped like one anywhere else (MDAnalysis:
  `UnboundLocalError ... 'numlines'`). unknown_first.parm7 loads.
- **RD-P02 / RD-P03 / RD-P04** Python texts replaced: `<file>: unexpected end of file in %FLAG MASS
  (expected 8 data lines, found 1)`; `<file>: required section %FLAG POINTERS is missing` (was an
  empty StopIteration); `<file>: required section %FLAG ANGLES_WITHOUT_HYDROGEN is missing` (was
  KeyError `'anga'`); `<file> line L: expected a %FLAG line after the 8 data lines of %FLAG MASS,
  found '...'`; `<file> line L: expected %FORMAT or %COMMENT after %FLAG MASS, found '...'`;
  `<file> line L, %FLAG RESIDUE_POINTER field 2: '1x' is not an integer`; a second header line
  without a second word is the "'TITLE' missing in header" error (was IndexError). POINTERS is read
  up to the next %FLAG: any count of at least 12 values is accepted (`<file>: %FLAG POINTERS needs
  at least 12 integers (found N)`, `<file> line L: %FLAG POINTERS value 'x' is not an integer`), so
  pointer_lines_short.parm7 (30 values in 3 lines) loads.
- **RD-P05** A blank numeric field before the last non-blank field of a line is an error (`<file>
  line L: blank field F in %FLAG MASS`) instead of shifting every later value. A blank text field
  there is the empty string (chain_blank.parm7: residue 3 has chain ID "", MDAnalysis drops it and
  ignores RESIDUE_CHAINID). Trailing blank fields stay allowed.
- **RD-P06** Data lines are counted from each section's %FORMAT (fields per line), not from a fixed
  per-section constant: charge_three_per_line.parm7 (`3E24.16`) loads.
- **RD-P07** `D` / `d` exponents are accepted in real fields.
- **RD-P08** RESIDUE_POINTER must start at 1, increase strictly and stay within the atoms
  (`<file>: %FLAG RESIDUE_POINTER is not increasing at residue 3 (11, then 6)`, `... residue 1 starts
  at atom 3; the first residue must start at atom 1`); MDAnalysis accepted such files with wrong
  residue membership (resptr_zero, resptr_unsorted, residue_pointer_offset are now errors).
- **RD-P09** `<file>: ATOMIC_NUMBER 119 of atom 1 is not an element` (was KeyError `119`).
- **RD-P10** Bond, angle and dihedral entries are checked: `<file>: %FLAG BONDS_INC_HYDROGEN entry k:
  atom index 4 is not a multiple of 3` / `refers to atom 41 of 40` / `is negative`, `... entry 5 has 2
  of 3 values` (negative third and fourth dihedral indices are Amber flags and stay accepted).
- **RD-P11** Topology errors never name MDAnalysis' parser classes. parm7: every error names the
  file (MDBridge never wrapped them). PDB: an error that is not one of the reader's own RM messages
  reads `<file>: cannot read PDB topology: <message>`. GRO: `<file> line 6: residue number '    x'
  (columns 1-5) is not an integer` (or `atom number ... (columns 16-20)`), `<file> line 5: no atom
  name in columns 11-15` (MDAnalysis: `Failed to load from the topology file ... with parser <class
  'MDAnalysis.topology.GROParser.GROParser'>.\nError: Couldn't read the following line ...` / `Missing
  atom name on line: 5`). *Tested:* bad_line.gro, missing_name.gro (loader topology_cases.json).
- **RESIDUE_CHAINID with the wrong count (LD-03)** is ignored as before (one segment SYSTEM), with
  `<file>: %FLAG RESIDUE_CHAINID has 7 values for 8 residues; chain IDs are not read`. *Tested:*
  test_parm7.cpp TopologyVariants (chain_short, chain_long).
- **RD-P12** Charges are the field divided by 18.2223 in float64 (MDAnalysis: float32 field divided
  by float32(18.2223) in float32, ~1e-8 relative error). Every expected value built from a parm7
  topology (topology dumps, selections, energies, the command line) uses these charges:
  improvements_readers_amber.apply() (called by independent_tables.apply()) wraps TOPParser.parse.
- **RD-P14** Count mismatches: `<file>: %FLAG MASS has 39 values for 40 atoms` (was "Length of masses
  does not match number of atoms.\nExpect: 40 Have: 39").

Amber NetCDF (src/io/netcdf.cpp; io/trajectories/cases.json):
- **RD-N01** Cell lengths (0, 0, 0) mean no cell for that frame (MDAnalysis: TypeError
  `'NoneType' object is not subscriptable`). zero_box.nc loads.
- **RD-N02** Atom-count mismatch: the message of the other readers, `The topology and NCDF trajectory
  files don't have the same number of atoms!\nTopology number of atoms 8\nTrajectory: <file> Number
  of atoms 7` (was advice about MDAnalysis' n_atoms argument); in a chain, the chain's checks report
  it (chain_atoms_differ_nc).
- **RD-N03** `<file>: variable 'coordinates' has no 'units' attribute (expected 'angstrom')`, `...
  'units' attribute of variable 'x' is not text ...` (were AttributeErrors).
- **RD-N04** `<file>: file is empty`, `<file>: unsupported NetCDF variant (CDF-5)`, `<file>: truncated
  header`, `<file>: corrupt NetCDF header (...)`, `<file>: file ends inside frame 3 (expected N bytes,
  found M)` (were numpy / mmap texts).
- **RD-N05** A streaming record count (0xFFFFFFFF) takes the records the file holds (new case
  streaming.nc); no records: `<file>: trajectory has no frames`.
- **RD-N06** dt is computed from the scaled times, as the frames' times are (MDAnalysis ignored the
  time scale_factor for dt): mda_scaled, coords_int_f64scale, coords_float_f64scale,
  time_short_scaled.
- **RD-N07** A cell variable without the frame dimension, (cell_spatial) / (cell_angular), is the
  cell of every frame (new case cell_constant.nc); any other shape is an error (cell_per_frame_scalar).
- **RD-N08** CDF-1 (32-bit offsets) files are read, with `WARNING: <file>: NetCDF file uses 32-bit
  offsets (CDF-1); the Amber convention asks for 64-bit offsets (CDF-2)`. cdf1, cdf1_dt and
  chain_bad_first load (expected data: MDAnalysis on a CDF-2 copy).
- **RD-N10** A file of more than one frame without a time variable (and no inputs.dt_ps) warns once,
  when it is opened, in the PDB / GRO readers' words (RD-B13): `<file>: the file has no time
  information; frames are taken to be 1 ps apart (set inputs.dt_ps for the real time step)`
  (MDAnalysis: `Reader has no dt information, set to 1.0 ps` on every dt / time access).

Amber ASCII (src/io/amber_ascii.cpp):
- **RD-A01** One frame without a box line is one non-periodic frame (MDAnalysis: empty
  StopIteration).
- **RD-A02** The layout is checked against the atom count at open: every coordinate line has 10
  fields, the frame's last 3N mod 10 (or 10); `<file> line L: expected 10 coordinates for 7 atoms,
  found 5 (does the topology match?)`; the line after the first frame is a box line (3 values) or the
  next frame's first line, else `<file> line 5: expected a box line (3 values) or the first 10
  coordinates of the next frame, found 4 values ...` (MDAnalysis mis-framed the data or gave numpy's
  reshape text).
- **RD-A03** One atom: box lines are recognised (every second line, when there is an even number of
  data lines and all of them hold three positive values); MDAnalysis read 3 periodic frames as 6.
- **RD-A04** Trailing blank lines are ignored; an incomplete last frame is ignored with `<file>: the
  last frame is incomplete (lines 18-18) and is ignored` (partial.mdcrd). Frame count = complete frames.
- **RD-A05** `<file> line 2, columns 1-8: '********' is not a number (Fortran overflow: the value does
  not fit the field)`; rst7: `<file> line 2: time 'abc' is not a number`, `... atom count 'abc' is
  not a non-negative integer`, `<file> line 6: expected 6 coordinates, found 3 (the file has 8
  atoms)`, `<file>: the file ends at line 3; 7 atoms need 6 lines`.
- **RD-A06, RD-B06** `<file>: file is empty` (mdcrd and rst7: cases mdcrd_empty, rst7_empty) /
  `<file>: no frames after the title` (were EOFErrors).

DCD / XTC / TRR (src/io/dcd.cpp, gromacs.cpp):
- **RD-D01** NSAVC = 0 is read as 1: dt = DELTA, times (frame + ISTART) * DELTA, with `<file>: the
  DCD header's frame interval (NSAVC) is 0; it is taken as 1 step` (MDAnalysis: ZeroDivisionError). Expected data: MDAnalysis on a copy with NSAVC 1.
- **RD-D02** `<file>: 4-dimensional DCD files are not supported` (was numpy's broadcast text).
- **RD-D03** XDR errors name the file and, while reading a frame, the frame: `<file> frame 4: XTC read
  error = compression`, `<file>: XDR read error = magic`.

Error classes: the new reader errors (NetCDF, mdcrd, rst7, DCD 4-D) are OSErrors, so a chain of files
reports them as they are (a ValueError would be rewrapped as "Unable to read <file> with <class
...>"); parm7 errors are ValueErrors.

## Selections (src/core/atom_selection.cpp, src/core/selections.cpp)

The selection language stays MDAnalysis' (keywords, precedence, matching). Where MDAnalysis or MDATK
selects the wrong atoms, accepts a typo silently or reports an internal Python text, MDBridge does
the following instead. *Tested* (all items): tests/data/selections/atom_selection_cases.json,
atom_selection_distance_cases.json and selections_cases.json, whose MDAnalysis / MDATK answers
tests/parity/improvements_selections.py turns into MDBridge's (its model of the grammar is checked
against MDAnalysis on every case whose answer is kept); tests/core/test_atom_selection.cpp and
test_selections.cpp; cli_cases.json (`info_groups`); framework_cases.json (cache keys).

- **SEL-01 Reversed `prop` operands.** `prop 6 > x` means x < 6 (the operator is mirrored).
  MDAnalysis negates it (x <= 6), so the boundary was included or excluded the wrong way round.
- **SEL-02 Malformed `prop` operators** are errors: `invalid operator '><' in prop (use <, >, <=,
  >=, ==, !=)`. MDAnalysis dropped text before the operator (`prop x ><5` meant x < 5) or failed with
  Python texts (`too many values to unpack`, `could not convert string to float: '>'`). The
  operator is the run of `<>=!` characters after the property; anything else is reported.
- **SEL-03 Arguments that run out** name what the keyword expects: `around: expected a distance
  and a selection, e.g. 'around 3.5 resname MG'` (likewise isolayer, sphzone, sphlayer, cyzone,
  cylayer), `point: expected x y z radius`, `atom: expected segid resid name`, `same: expected
  'same <property> as <selection>'`, `prop: expected 'prop [abs] <property> <op> <value>'`.
  MDAnalysis: `float() argument must be ... not 'NoneType'`, `pop from an empty deque`, `deque
  index out of range`, `argument of type 'NoneType' is not iterable`, ...
- **SEL-04 End of the selection.** `Unexpected end of selection after 'and'` (any token: `name`,
  `not`, `(`, `group`, ...) and `Unexpected end of selection (expected ')')` replace MDAnalysis'
  `Unknown selection token: 'None'` / `Unexpected token 'None'` / `Unexpected token: 'None'
  Expected: ')'`. A selection of blanks only is `empty selection` (an empty string still selects
  nothing).
- **SEL-05 `atom` arguments** may not be keywords, `and`, `or` or parentheses (`atom A 1 and` is
  `atom: expected segid resid name`; MDAnalysis read `and` as the atom name and selected nothing). A
  leading backslash escapes, as for other values (`atom \and 1 P`: segid `and`).
- **SEL-06 Distance keywords** no longer raise the grid search's errors (`Cutoff ... too large for
  box (max ...)`, `Cutoff must be positive`, `Any of the box dimensions cannot be 0`), which
  depended on the atom counts and the box; every cutoff gets exact minimum-image distances. A
  negative distance or radius is always an error (`around: distance must not be negative`,
  `sphzone: radius must not be negative`, also isolayer, sphlayer, point); MDAnalysis selected
  nothing for a negative cutoff on the brute-force route. *Also tested:* AtomSelection.CutoffsBeyondTheBox.
- **SEL-08 Messages.** `Unknown 'same' property 'RESID' (choose one of: fragment, x, y, ...)`
  instead of `Unknown same property : RESIDChoose one of : dict_keys([...])`; an unknown token whose
  lower-case spelling is a keyword gets a hint: `Unknown selection token: 'AND' (keywords are lower
  case: did you mean 'and'?)`. New messages are not wrapped in `Selection failed: '...'`.
- **SEL-09 Float ranges** (`mass`, `charge`, ...) accept in each bound what a single value accepts
  (`mass 12:+13`, `mass -inf:inf`, `charge 1E-1:2`); a piece that is no number and no range of two
  numbers is `Failed to parse number: <piece>` (`mass :3`, `mass 1:2abc`). MDAnalysis failed with
  `could not convert string to float: ''` and ignored text after the upper bound.
- **SEL-10 Reversed `resid` ranges** (lower resid above the upper one) select nothing, with or
  without insertion codes. With insertion codes MDAnalysis selected both end residues (`resid 3:1`
  on a PDB file: residues 3 and 1).
- **SEL-13 Element fallback per atom** (`atomElements`, used by heavy-atom filters, H-bonds, Mg
  ligands, generic LJ in energy and the structure panel): an atom without an element gets the guess
  from its own name; atoms with one keep it. MDATK guessed every atom of the group from its name as
  soon as one atom had no element, so an atom's element depended on the group it was in (one Amber
  extra point switched a whole system to name guessing; `CL1` became C). No recorded case changes
  (their topologies have all elements or none), but results can change for topologies with partly
  missing elements, so the analyses that use elements are bumped (contacts, energy, hbonds, mg,
  basepairs, trace, proximity, atoms; "Analysis cache" above).
- **SEL-14 Invalid group names in `group:`.** The name after `group:` runs to the next blank or
  parenthesis; a name a group may not have (validate()'s `[A-Za-z0-9_.-]+`), or a defined group with
  such a name, is `group name 'bad name' may only contain letters, digits, '_', '.' and '-'`. MDATK
  resolved `group:bad name` as the group `bad` ("unknown residue group 'bad'") and left `group:Ünï`
  to MDAnalysis (`Unknown selection token: 'group:Ünï'`); `group:A+B` and `group:a'b` were names.
- **SEL-15 Repeated group names:** the first group of a name wins in `group:` expansion and in the
  groups recorded in cache keys, as in `core::group()` and the GUI (MDATK: the last one in the
  expansion). validate() still reports the repetition.
- **SEL-16 (kept)** A residue number in a group's ranges means `resid N`, as the `resid` keyword:
  `45` excludes `45A` on files with insertion codes; ranges `44-46` include them. **SEL-11 (kept)**
  `and` and `or` have the same precedence and combine left to right (`A or B and C` is
  `(A or B) and C`); parenthesise to group otherwise.

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

## Analyses (non-periodic improvements)

Expected data: `tests/parity/improvements_analyses.py` patches MDATK (and the RMSD probes) to the
behaviour below before `make_analysis_cases.py`, `make_analysis_framework_cases.py`,
`make_pipeline_cases.py`, `make_rmsd_cases.py` and `cli_parity.py --write` record it; hand-made
checks are in tests/analysis/test_improvements_analyses.cpp. Version bumps (new cache keys): energy 2,
hbonds 2, rmsd 2, proximity 2, atoms 3.

- **Setup order (SEL-17, AN-13).** Every analysis' setup starts on the first analysed frame
  (MDATK: on frame 0 of the trajectory, or on the rmsd reference frame for analyses set up after
  rmsd), so selections with distance keywords and the `atoms` analysis' guessed bonds no longer
  depend on the analysis order or the rmsd settings; `atoms` guesses its bonds in prepare() from the
  first analysed frame. *Tested:* setup_frame of every recorded case (all the first analysed frame),
  test_rmsd_analysis.cpp `EverySetupStartsOnTheFirstAnalysedFrame`.
- **Periodicity of selections (SEL-17).** `around`, `sphzone`, ... in the analyses' selections use the
  cell only with `analyses.use_pbc` (MDATK: whenever the frame has one). *Tested:*
  `SelectionsHonourUsePbc`.
- **use_parm7 with another topology (AN-06).** `WARNING: energy: use_parm7 is set but the topology
  <file> is not a parm7 file; using generic element parameters` (MDATK: silent). *Tested:* systems_cases
  logs (pdb / gro / noh / aromatic), test_energy.cpp.
- **Energy exclusions without topology bonds (AN-07).** With `exclude_bonded > 0` and no bonds:
  `WARNING: energy: the topology has no bonds; bonded pairs are excluded along bonds guessed from
  distances on the first analysed frame`, and the exclusions follow the bonds the `atoms` guesser
  (below) finds on the whole system (MDATK: nothing excluded, silently). *Tested:* micro and
  systems_cases energy logs and values.
- **Exclusions along any path (AN-08).** Pairs joined by at most `exclude_bonded` bonds through any
  atoms are excluded (MDATK: only paths that start at a direct A-B bond and run through A). *Tested:*
  `ExclusionsThroughOtherAtoms`; the recorded chains give the same keys.
- **Guessed hydrogens (AN-09, AN-10).** Without topology bonds a hydrogen goes to the donor nearest by
  the pair search's distance (minimum image with use_pbc; MDATK: the plain distance), and a nitrogen
  that gets a hydrogen is not an acceptor, as with topology bonds (MDATK: kept as acceptor, so PDB
  and parm7 files of one structure gave different counts). *Tested:* systems_cases pdb_nobonds /
  gro_nobonds H-bond counts and maps.
- **Bond guess (AN-14).** The `atoms` analysis (and AN-07) guesses bonds from covalent radii
  (Cordero et al. 2008, data/sources/covalent_radii; distance <= r_i + r_j + 0.4 Å, carbon's radius
  for unknown elements), never between two hydrogens, never to an atom of a one-atom residue (ions:
  covalent radii would bond Mg to its water), and between residues only for consecutive residues of
  one chain (same segment and chain ID). MDATK: 1.95 Å (a 2.05 Å disulfide missed), hydrogens within
  1.25 Å, adjacent residue indices across chains. *Tested:* systems_cases `atoms` bonds,
  `CovalentRadiusBondGuess`.
- **Proximity members (AN-15).** Residues of the proximity set are at distance 0 also when no other
  structure atoms remain (MDATK: left at inf). *Tested:* demo_dry_variants structure_target_only.
- **Degenerate RMSD (AN-16).** For fewer than 3 atoms, collinear or coincident sets (by the weighted
  inner-product matrix: q <= 1e-10 p^2) or no convergence in 50 Newton steps, the RMSD is the exact
  superposition RMSD: the optimal rotation from a Jacobi eigen-solution of the 4 x 4 key matrix and
  the residual summed after rotating; one atom gives 0 (MDAnalysis: NaN, or values off by up to 10 Å
  that change chaotically with the last input bit). Reached by the default "Backbone" target of a
  1-2 residue group. *Tested:* rmsd_cases.json (expected: numpy's SVD value, within 1e-10 of the
  coordinates' size; the nearly collinear float32 sets are that ill-conditioned) and its demo series.
- **RMSD near 0 (AN-17).** Where E0 - lambda < 4096 eps E0 the residual is recomputed with the
  optimal rotation; identical structures (the reference frame against itself) give exactly 0
  (MDAnalysis: up to 5e-6 Å). The register proposed 64 eps; the Newton iteration stops at a relative
  step of 1e-14, and identical sets were seen at 250 eps, so 4096 eps (about 1e-12 E0) is used. The
  QCP value is kept bit for bit everywhere else. *Tested:* rmsd_cases.json (identical, translated,
  rotated, noise <= 1e-6 families), the rmsd series of every recorded case, cli_cases.json.
- **Fits of fewer than 3 or collinear atoms (AN-18).** Translation only, with `WARNING: RMSD target
  'X': fit selection has N atoms; superposition is translation only` (or `N collinear atoms`) (MDATK:
  translation only silently below 3 atoms, an arbitrary LAPACK-dependent rotation about the line for
  collinear ones). *Tested:* `TranslationOnlyFit`.
- **use_pbc without a usable cell (AN-19).** One warning per run: `use_pbc is on but trajectory frame
  K has no usable periodic cell; distances are not periodic` (missing cell, a length or angle not
  finite or not positive, an angle of 180° or more; MDATK: silent). The trajectory file is not named:
  the runner sees one concatenated trajectory. *Tested:* `UsePbcWithoutUsableCell`.
- **Non-finite coordinates (AN-20).** Atoms with a NaN or infinite coordinate are in no pair on the
  open path too (MDATK: the analysis aborted with scipy's `data must be finite, check for nan or inf
  values`); one warning per run of consecutive analysed frames: `frame K: N atoms have non-finite
  coordinates; they are left out of distance searches` (`frames K1-K2: up to N atoms ...`).
  *Tested:* test_spatial.cpp (the recorded error cases), `NonFiniteCoordinates`.
- **Unknown structure groups (AN-21).** `WARNING: structural view: unknown residue group 'X' ignored`
  from the trace analysis' setup (once per run of the structural view; MDATK: silent). *Tested:*
  demo_dry_variants structure_unknown_group.

## Planned (from the improvement register; add each entry here when it lands)

Periodic distances and minimum images for every cell (AN-01–AN-04); base pairs and stacking on
wrapped trajectories (AN-11, AN-12).
