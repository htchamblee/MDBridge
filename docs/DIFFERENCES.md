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

## Amber, DCD and GROMACS readers (docs/IMPROVEMENTS.md 2.1, 2.4-2.6)

Expected data: the generators record MDAnalysis; tests/parity/improvements_readers_amber.py turns
it into MDBridge's behaviour (per item, below). Error texts follow the reader-message rule: they name
the file and, where known, the line, columns, section or frame. Warnings marked "(warning pending)"
are not printed yet: MDBridge has no reader-warning channel until LD-03; the code carries the text
as `TODO(LD-03)`.

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
- **RD-P11** (parm7 part) every parm7 topology error names the file; MDBridge never wrapped them in
  MDAnalysis' parser-class text. The PDB / GRO wrapper (pdb_topology.cpp) belongs to the PDB chunk.
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
- **RD-N08** CDF-1 (32-bit offsets) files are read (warning pending: "uses 32-bit offsets (CDF-1)").
  cdf1, cdf1_dt and chain_bad_first load (expected data: MDAnalysis on a CDF-2 copy).
- **RD-N10** The "no time information; dt set to 1.0 ps" warning is to be issued once per file
  (warning pending).

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
- **RD-A04** Trailing blank lines are ignored; an incomplete last frame is ignored (warning pending:
  "the last frame is incomplete"). Frame count = complete frames.
- **RD-A05** `<file> line 2, columns 1-8: '********' is not a number (Fortran overflow: the value does
  not fit the field)`; rst7: `<file> line 2: time 'abc' is not a number`, `... atom count 'abc' is
  not a non-negative integer`, `<file> line 6: expected 6 coordinates, found 3 (the file has 8
  atoms)`, `<file>: the file ends at line 3; 7 atoms need 6 lines`.
- **RD-A06** `<file>: file is empty` / `<file>: no frames after the title` (were EOFErrors).

DCD / XTC / TRR (src/io/dcd.cpp, gromacs.cpp):
- **RD-D01** NSAVC = 0 is read as 1: dt = DELTA, times (frame + ISTART) * DELTA (warning pending)
  (MDAnalysis: ZeroDivisionError). Expected data: MDAnalysis on a copy with NSAVC 1.
- **RD-D02** `<file>: 4-dimensional DCD files are not supported` (was numpy's broadcast text).
- **RD-D03** XDR errors name the file and, while reading a frame, the frame: `<file> frame 4: XTC read
  error = compression`, `<file>: XDR read error = magic`.

Error classes: the new reader errors (NetCDF, mdcrd, rst7, DCD 4-D) are OSErrors, so a chain of files
reports them as they are (a ValueError would be rewrapped as "Unable to read <file> with <class
...>"); parm7 errors are ValueErrors.

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
`prop` operand reversal (SEL-01, SEL-02); per-atom element fallback (SEL-13); reader fixes (RD-B*,
RD-G*, LD-*); validation of every cutoff (PRJ-01); readable messages in place of Python's internal texts
(SEL-03/04/08, RD-*, CLI-01/03/04).
