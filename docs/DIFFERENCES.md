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
  improvement of this release: energy 2 (AN-01-AN-05, AN-07, AN-08, RD-P12's float64 charges,
  SEL-13, AN-25-AN-27), hbonds 2 (AN-01-AN-04, AN-09, AN-10, SEL-13, distance-only pairs counted once),
  rmsd 2 (AN-16-AN-18), proximity 2 (AN-01-AN-04,
  AN-15, SEL-13), atoms 3 (AN-13, AN-14, SEL-13), contacts 2 (AN-01-AN-04, SEL-13, a group with itself, TOP-01)
  and mg 2 (AN-01-AN-04, SEL-13),
  basepairs 2 (AN-01-AN-04, AN-11, SEL-13, NA-01), stacking 2 (AN-01-AN-04, AN-12, NA-01) and trace 5 (SEL-13, TOP-01, NA-01; its key
  now covers every group). The fixes taken from MDATK 0.1.x ("Fixes adopted from MDATK 0.1.x" below) are
  in this release's bump; MDATK 0.1.1 numbers its own versions independently (energy 2, hbonds 2,
  contacts 2, trace 5 mean other results there), which matters only when the parity reference moves. The test side gives MDATK's analyses these
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
identical. With a periodic box (`use_pbc`) the pairs are sorted by (i, j) too (AN-03, below).
*Why:* a defined order that does not depend on the build of a library. Visible consequences:
- energy sums over pairs (elec, vdw, net force) can differ in the last bits;
- H-bond hydrogen guess: a hydrogen exactly as far from two donors goes to the donor with the lowest
  index (MDATK: the first in tree order);
- stacking map entries within a frame are in (i, j) order;
- `atoms` analysis: distance-guessed bond rows are sorted by (i, j).

*Tested:* test_contacts_hbonds.cpp (tie case), test_mg_stacking.cpp (sorted comparison),
test_structure.cpp, test_energy.cpp (tolerances).

## Periodic geometry (`use_pbc`; src/analysis/spatial.cpp)

MDATK computes periodic pairs with MDAnalysis' `capped_distance` and corrects difference vectors for
orthorhombic cells only. MDBridge has one float64 minimum image for every cell shape
(`PeriodicCell` in spatial.hpp) and uses it for the pairs and for every vector. The analyses whose
periodic results change have version 2 (contacts, energy, hbonds, mg, stacking, proximity,
basepairs; "Analysis cache" above), so caches written before recompute. *Tested:* the expected data of every case with
`use_pbc` (and of the distance cases) comes from tests/parity/improvements_periodic.py, an
independent numpy implementation of the same geometry that the generators patch into MDATK
(bit-identical results); test_periodic.cpp also checks the cases below directly, and the cell list
against comparing every pair.

- **AN-01 One exact periodic distance.** MDAnalysis' brute-force route scaled each axis difference
  by L * fl32(1 / L), so every distance was off by up to 2^-24 relative, even for pairs that need no
  image: a pair at exactly 3.5 Å was kept for L = 25 and dropped for L = 30, and points far outside
  the cell got errors of up to 0.5 Å (1e10 Å away: 1.48 instead of 1.0). The grid route wrapped in
  float32. MDBridge takes the difference of the float32 positions in float64 and reduces it by
  whole cells (v - L * nearbyint(v / L), repeated while |v| > L / 2; triclinic cells along c, b, a
  and then the shortest of the 27 neighbouring images). A pair needing no image has exactly the
  open path's distance, and a pair is kept iff d² <= c², the open path's test (MDAnalysis' brute
  force tested sqrt(d²) <= c). Periodic energies change in about the sixth digit (MDAnalysis' grid
  wrapped in float32), counts only at boundary pairs. *Why:* results must not depend on the box
  size or route.
  *Tested:* distance_cases.json `exact` column (495 cases); test_periodic.cpp
  ExactBoundaryForEveryBoxSize, FarOutsidePointsAreExact.
- **AN-02 No size-dependent errors.** `Cutoff must be positive`, `Cutoff <c> too large for box (max
  <m>)` and `Any of the box dimensions cannot be 0` came only from MDAnalysis' grid route (both sets
  of at least 10 points and a small cutoff, or n·m >= 1e8), so a project could fail on a large
  system and work on a small one. MDBridge never raises them: a cell list with bins at least one
  cutoff wide prunes candidates (all candidates along a direction narrower than three bins; points
  more than 10^4 cells away are compared with everything), and the result does not depend on the
  sizes. A cell without periodic geometry (angles that admit no third box vector, an angle outside
  (0, 180), a non-finite triclinic length) gives plain distances; an orthogonal axis with an
  infinite or zero length is not periodic. A point with a NaN or infinite coordinate is in no pair.
  *Tested:* NoSizeDependentErrorsAndSameResultForAnySize, CellListMatchesAllPairs, the former error
  cases of distance_cases.json.
- **AN-03 Pair order.** Periodic pairs come sorted by (i, j), as without a box (MDAnalysis: brute
  force row-major, grid cell by cell with j descending). The consequences listed under "Order effects
  of the neighbour search" now apply with `use_pbc` too. *Tested:* test_periodic.cpp (sorted order),
  the `use_pbc` analysis cases.
- **AN-04 Minimum-image vectors for triclinic cells.** MDATK's `minimum_image_vectors` corrected
  orthorhombic cells only (angles within `np.allclose` of 90), so H-bond angles, Mg geometry,
  stacking offsets and the energy net force were wrong across a face of a triclinic cell (net force
  195 vs 26 kcal/mol/Å for the same system), and for angles such as 90.0005° the pairs used
  triclinic images while the vectors were not corrected. MDBridge uses the pair search's image for
  every vector and one orthogonality test everywhere: the three angles exactly 90 after rounding to
  float32. *Tested:* spatial_extra_cases.json `improved`; systems_cases.json `periodic_tric_*`.
- **AN-05 Energy net force direction.** MDATK multiplied the real difference vector by the force
  divided by the clamped distance max(r, 0.5), so for r < 0.5 Å the force was scaled by r / 0.5.
  MDBridge uses the unit vector times the force at the clamped distance (zero for r = 0); for
  r >= 0.5 nothing changes, bit for bit. *Tested:* test_energy.cpp
  NetForceAtClampedDistanceHasFullMagnitude (open and across a face); the energy cases (expected
  values from the patched MDATK step in improvements_periodic.py).
- **AN-11 Base pairs on wrapped trajectories.** MDATK ignored the cell, so a duplex wrapped across a
  face had no base pairs (every `periodic_*` case: 0 per frame). With `use_pbc` and a periodic cell
  every donor-acceptor and C1'-C1' distance the detector reads is a minimum-image distance; without
  `use_pbc` nothing changes. The register proposed making each residue pair whole by its first
  atoms; per-distance images are used instead because the distances read are at most 12.5 Å while
  reference atoms across a pair (P-P about 18 Å) can be more than half a 25 Å cell apart, where whole
  residues would pick the wrong image. *Tested:* systems_cases.json `periodic_*_pbc` (basepairs).
- **AN-12 Stacking rings split by wrapping.** Centroid, normal and planarity came from raw positions,
  so a ring wrapped atom by atom was "non-planar" and missed. With `use_pbc` each ring atom is put at
  the minimum image of the ring's first atom first. *Tested:* systems_cases.json `periodic_*_pbc`
  (stacking).

## Readers and topologies

- **Self-bonds** (`CONECT 1 1`, kept by MDAnalysis) are dropped; the topology model has none.
- **TOP-01: chains from bonds only when the bonds hold the residues together** (MDATK 0.1.0,
  analysis/structure.py `chain_ids`, lines 75-84; core::chainIds). Without chain IDs or segments that
  tell chains apart, the bonded fragments do; but a PDB file whose only CONECT records are its
  ligand's leaves every other atom a fragment of its own, and MDATK 0.0.4.1 made every residue its own
  chain (no backbone links in the structure view, wrong chain colours and legend). When some residue's
  atoms lie in more than one fragment, every residue is in chain 0. Part of the trace analysis'
  version 5; the contacts of a group with itself take their strand ends from these chains ("Contacts
  of a group with itself"). *Tested:* test_topology.cpp (ChainsFromBondsOnlyWhenTheyHoldTheResidues);
  systems_cases.json `noh_conect_self` (the heavy-atom demo without chain IDs, CONECT for one bond of
  residue 24 only: one chain, so the contacts of residues 1-24 against themselves leave out every
  sequence neighbour, 12-13 included, and the trace has one chain); the generators' chain_ids is 0.1.0's
  (tests/parity/improvements_core_io.py).
- **`info` residue-name ranking:** equal counts are ordered by name (MDATK: numpy's argsort order,
  which depends on the CPU's sort kernel).
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
- **The time step override replaces the files' times (LD-05).** With `inputs.dt_ps` set, frame k of
  the trajectory (of the chain, counted from 0) is at k x dt_ps whatever the files' time stamps, for
  every consumer (analysis time axes, the frame axis without analyses, frames chosen by time); a file
  that stores its own times with another step says so: `inputs.dt_ps (2.5 ps) replaces the time step in
  <file> (10.0 ps)`. The cache keys' input fingerprint gets `"time_axis": "frame x dt_ps"` when the
  override is set, so results cached with the file's times are recomputed. MDATK 0.0.4.1 (MDAnalysis'
  `dt`) kept a NetCDF file's stored times, so its time axis disagreed with the PDB export by time (frame x
  dt) and, for a chain of NetCDF files each starting at its own first time, was not monotonic; MDBridge
  first kept the file's times with a warning. This follows MDATK 0.1.0 (`mdatk/io/loader.py` _FrameTimes,
  `analysis/base.py` run_analyses, `pipeline.py` input_fingerprint): MDBridge applies it in
  `io::openTrajectory` / `openTrajectories` (`io::withFrameTimes`) and to coordinates read from the
  topology file (MDATK 0.1.0's run_analyses uses frame x dt_ps there too). *Tested:* `Loader.SystemWarnings`,
  trajectory cases.json (every case with a dt), pipeline_cases.json, cli_cases.json, the analysis cases
  with `inputs.dt_ps`.
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
  case: did you mean 'and'?)`, and since SEL-19 a misspelt one too: `Unknown selection token:
  'resnam' (did you mean 'resname'?)`. New messages are not wrapped in `Selection failed: '...'`.
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
- **SEL-19 Selection errors: MDATK 0.1.x's frame, MDBridge's messages, hints** (follows MDATK
  0.1.1 `mdatk/selections.py` `syntax_error`, `_syntax_hint`, `suggest`; also covers SEL-04,
  SEL-08 and CLI-07). `select()` reports `invalid selection 'X' (expanded: 'Y'): <message>
  (<hint>)`: the text as written (stripped), its expansion only when it differs (also from the text
  after `sel:`), the parser's own message (SEL-03, SEL-04, SEL-08 ...; MDATK 0.1.x shows
  `MDAnalysis says: <MDAnalysis' message>`, which MDBridge has no use for) and a hint where one
  helps. MDATK 0.0.4.1 printed `invalid selection 'X' -> 'Y': ...`, the expansion even when it was
  the text itself (`invalid selection 'name (' -> 'name (': ...`). The hints, as 0.1.x's:
  - a token that is (by case or a typo) a group's name: `Unknown selection token: 'Target' (to use
    the residue group write group:Target)`, in place of a keyword hint;
  - a misspelt keyword: `(did you mean 'resname'?)` (`keywords are lower case: did you mean 'and'?`
    when only the case differs, as SEL-08); no hint for `and` / `or` themselves where a selection
    should start, or when nothing is close (0.1.x: `'foo' is not a selection keyword`, which repeats
    the message);
  - `group Target` (MDAnalysis' own group keyword): `Selection failed: 'Failed to find group:
    Target' (write group:Target, with a colon, to use the residue group)`, or `(residue groups are
    written group:<Name>)` when no group is close;
  - `Unexpected token: 'x' Expected: ')'`: `(a parenthesis is not closed)`. The end of the selection
    needs no hint: SEL-04's `Unexpected end of selection (expected ')')` says it (0.1.x hints at
    MDAnalysis' `'None'` token instead).

  "Close" is 0.1.x's `suggest`: a candidate equal ignoring case, else Python difflib's
  `get_close_matches(..., n=1, cutoff=0.6)` on the lower-case spellings, reimplemented in
  `core::suggestName` / `core::similarity` (src/core/suggest.cpp; no suggestion for a name of 200 or
  more characters, where difflib's junk heuristic would start). The errors carry their kind and
  offending token (`SelectionError::kind()`, `token()`), as 0.1.x's do, so that a later change can
  put the setting's location in front of the message. *Also tested:* `Selections.SuggestLikeDifflib`
  (ratios and suggestions against difflib, selections_cases.json `similarity` / `suggest`),
  `Selections.MessagesNameTheProblem`; cli_cases.json (`info_groups`, `export_bad_selection`).
- **SEL-20 Unknown groups suggest the group meant** (follows MDATK 0.1.1 `selections.py`
  `unknown_group_error`). `unknown residue group 'Tagret' (defined: Aptamer, Target); did you mean
  'Target'?`, and when a group's own selection names it, `unknown residue group 'Tagret' in the
  selection of group 'Uses' (defined: ...)` (MDATK 0.0.4.1 named neither). The defined groups are
  listed once each, in order. *Tested:* selections_cases.json (groups set `typos`).
- **SEL-21 Circular group definitions named** (follows MDATK 0.1.1 `selections.py`
  `expand_groups`). A group whose expansion reaches itself is `group definitions are nested in a
  circle: group:Circle1 → group:Circle2 → group:Circle3 → group:Circle1`, from the group where the
  expansion started; MDATK 0.0.4.1 expanded until the depth limit and reported `group definitions
  are nested too deeply (circular reference?)`, which now only nine nested distinct groups give.
  *Tested:* selections_cases.json (`Loop1` / `Loop2`, `Circle1` ...).
- **SEL-22 Why a selection matched nothing** (follows MDATK 0.1.1 `selections.py`
  `explain_empty`). `selection 'group:Far' matched no atoms: group 'Far' (residues '900-910')
  matches no atoms - the structure's residues are numbered 1-302`: each group the text uses is tried
  on its own, the empty ones are named with their definition (residue ranges add the structure's
  residue numbers); otherwise each top-level and / or part is tried: `selection 'name P and resid
  999' matched no atoms: the part 'resid 999' matches no atoms`, or `each part matches atoms, but
  not in this combination (check 'and' / 'or')`. Parts that do not stand alone (`around 3 ...`) are
  skipped. The tries use the same structure and frame as the selection (with `use_pbc` off, an open
  box, SEL-17). The text is shown stripped. MDATK 0.0.4.1: `selection 'X' matched no atoms`.
  *Tested:* selections_cases.json, cli_cases.json (`info_groups`), analysis micro / variant cases.
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
- **Windows paths** (src/core/include/mdbridge/core/ntpath.hpp). MDATK 0.0.4.1 has a Windows build
  (a beta: `packaging/windows`, MSYS2's Python bundled with PyInstaller); there `resolve_paths` uses
  Python's `ntpath` as it is: `expanduser`, the path itself when `isabs`, else
  `normpath(join(folder, path))`. MDBridge's rules: `/` and `\` both separate; drive, UNC and
  root-relative paths are absolute (root-relative takes the project folder's drive); `C:x` resolves
  against the project folder when that is on drive C, otherwise against `C:/`; `~` expands to
  %USERPROFILE% (else %HOMEDRIVE%%HOMEPATH%) and `~user` as `ntpath.expanduser` guesses it (the
  profile folder's sibling `C:\Users\user` when the profile folder is named %USERNAME%, else
  unchanged); results are normalised and written with `/` (`C:\data\..\a.pdb` → `C:/a.pdb`); device
  paths (`\\?\`, `\\.\`) are kept as written; comparison is case-sensitive (MDATK never calls
  `normcase`; the same-drive test of `C:x` ignores the letter's case, as `ntpath.join` does). On
  Linux and macOS paths follow posixpath exactly, as MDATK does. Deliberate differences from MDATK's
  Windows build (checked with Python 3.12's and 3.13's `ntpath`), each for one reason: a stored path
  means one file, whatever the process's per-drive current folders, the Python version bundled or
  how the path was typed.
  1. *Separators and normalisation.* MDATK writes `\` in resolved relative paths
     (`C:\projects\p\sub\a.pdb`) and keeps absolute ones exactly as given (`C:\data\..\a.pdb`,
     `C:\Users\me/a.pdb` from `~/a.pdb`, `\\server\share\x\..\a.pdb`); MDBridge normalises every
     result to one generic form with `/` (the same file; one text to compare, log and show).
  2. *A drive-relative path on another drive.* `D:a.pdb` with the project on C: MDATK keeps `D:a.pdb`,
     which Windows resolves against the process's current folder on drive D when the file is opened
     (a GUI program does not control it); MDBridge uses `D:/a.pdb` (`D:` alone: `D:/`).
  3. *Root-relative paths and the Python version.* `\data\a.pdb`: Python 3.13 (and MDBridge) put it on
     the project folder's drive; Python 3.12 and older count it as absolute and keep `\data\a.pdb`
     (the current drive when opened). Which Python MDATK's bundle has depends on MSYS2 at build time;
     MDBridge's rule does not.
  4. *`cleaning.cpptraj_path`* counts as a path when it holds `/` or `\` (MDATK: `os.sep`, `\` only),
     so `bin/cpptraj` is resolved against the project folder as `bin\cpptraj` is, instead of against
     the current folder when cpptraj is started.
  5. *Output folders against the current directory* (`output.out_dir` relative or empty): MDATK's
     `abspath` asks Windows (`GetFullPathNameW`: per-drive current folders; trailing dots and spaces
     of names removed); MDBridge applies the string rules above to the current directory.
  6. The recent-projects list is in `%APPDATA%\MDBridge\settings.toml` (owner decision D16, PRJ-12;
     MDATK's Windows build: `%USERPROFILE%\.config\mdatk\settings.toml`, not shared).
  *Tested:* tests/core/test_ntpath.cpp on every platform (`~user` against Python 3.13's
  `ntpath.expanduser`).
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
  `out_dir` / `cache_dir` / `workflow_dir` key in a project file) names the field meant in the
  unknown-key warning (PRJ-09): `project.out_dir (did you mean output.out_dir?)`
  (`output.cache_dir`, `workflow.out_dir`). *Tested:* from_dict_cases.json (properties as keys),
  cli_cases.json (init_set_property, init_set_cache_dir).
- **PRJ-09: unknown keys are kept (owner decision D9 of milestone 5).** MDATK warns about every key
  its model does not know (`project.workflow: ignoring unknown key 'parallel'`, one line each) and
  drops it, so saving a project of a later MDATK lost its settings (0.1.x's plan budgets, adopted
  runs, folder tags). MDBridge keeps every unknown key of every table, at any depth (lists of tables
  included), with its value and type as read, and writes it back after the table's own fields; one
  warning names them all: `project: 2 unknown keys kept and saved back unchanged (from a newer
  MDATK, or misspelt?): project.bogus, project.inputs.nope` (`core::unknownKeysWarning`; none when
  the file cannot be read). `--set` warns only about the keys it adds and, as MDATK, cannot set a
  key the model does not know (`unknown parameter`). Unknown keys inside an analysis section are part
  of that analysis' cache key, as everything in the section is. The keys MDATK 0.1.x added to the
  plan are modelled with 0.1.1's comments (config.py `WorkflowParams`, `PlanNode`):
  `workflow.supersede_runs`, `parallel`, `cpu_cores`, `gpus`, `last_node` and per step `adopted`,
  `folder_tag`; they are written only when the file had them, so files of MDATK 0.0.4.1 come out as
  before (unset means what 0.1.x reads for a file without them: `parallel` 1). Step logs of 0.1.x
  (statuses `waiting` and `stopped`, the step `workflow supersede`) are read and kept as they are.
  MDBridge still writes its own comments (0.0.4.1's, so options of step types 0.0.4.1 lacks have
  none). `mdatk_version` stays what a later MDATK wrote (`0.1.1` is saved as `0.1.1`); files without
  it, or of 0.0.4.1 and older, get `0.0.4.1` (core::writtenMdatkVersion). MDATK 0.1.1 only writes
  the key and drops it when reading (config.py `from_dict`, projectfile.py `read`); nothing branches
  on it. *Tested:* project-file/mdatk-0.1.1 (files MDATK 0.1.1 writes,
  tests/parity/make_mdatk_latest_projects.py run by regenerate_all.py with `--mdatk-latest`): what
  MDBridge reads equals what 0.1.1 reads, and a save keeps every key, value and order, the step log
  byte for byte and `mdatk_version` (test_project_io.cpp `ProjectCompat`); from_dict_cases.json and
  cli_cases.json with tests/parity/improvements_core_io.py applying the same rule to MDATK.
- **PRJ-10: the step-log marker inside the settings is no log** (MDATK 0.1.0, projectfile.py
  `_log_start`, lines 399-418, and `append_log`, 488-492). MDATK 0.0.4.1 took the first marker line
  anywhere in a project file for the start of its log, so the marker pasted into the notes or a
  custom script (with a piece of another project file) made the next entry go without a marker of its
  own and the next save cut the settings at the notes: the file was no TOML. A marker starts the log
  only where the text before it is complete TOML (else the first one, for a damaged file), for saving,
  appending and moving the log out. A separate log file (`log_file`) that is gone starts again with the
  next entry instead of failing it. *Tested:* test_project_io.cpp (LogMarkerInsideTheSettingsIsNoLog,
  MissingSeparateLogStartsAgain, MdatkLatestProjectsSurviveASave: 0.1.1's file holds the marker in
  its notes and a script).
- **PRJ-11: comments stay valid TOML** (MDATK 0.0.13, tomlio.py `_comment_text`, lines 59-65 and
  171-174). Control characters other than tab in a comment (a free-form key such as a panel option
  named with one, which its comment quotes) become spaces, and a header item with line breaks is one
  comment line per line; MDATK 0.0.4.1 wrote them as they were, and the file did not load. Text
  without them is unchanged. *Tested:* test_toml_writer.cpp (CommentsHoldNoControlCharacters,
  MDATK 0.1.1's texts).
- **PRJ-12: the user settings file keeps its other entries** (AP-17; MDATK 0.0.14, projectfile.py
  `_user_settings` / `add_recent`). Updating the recent list rewrote the settings file MDATK and
  MDBridge share with `recent` alone, deleting what MDATK keeps there (RNAdo's `cat_proof`). Every
  other entry now stays in its place, with 0.0.14's comments, and the file is written to `.tmp` and
  renamed (MDATK writes it in place). Locations per owner decision D16: Linux keeps the shared
  `$XDG_CONFIG_HOME/mdatk/settings.toml` (else `~/.config/mdatk`) for the list, MDBridge's own
  settings go to `$XDG_CONFIG_HOME/mdbridge/` (`core::mdbridgeConfigDir`); Windows keeps both in
  `%APPDATA%\MDBridge\`. *Tested:* test_project_io.cpp (RecentListKeepsOtherSettings: MDATK 0.1.1's
  text for the same file; ConfigFolders).
- **RN-20: `visual.renderer = "vulkan"` is accepted** (owner decision D7, docs/specs/m4-plan.md).
  MDATK knows `auto`, `software` and `vtk` (VTK, optional); MDBridge's GPU renderer is Vulkan, so
  `vulkan` is a valid value besides them (`vtk` keeps meaning "the GPU renderer", RN-16), and the
  message for an unknown value reads `Unknown renderer 'opengl' (auto, software, vtk or vulkan).`
  Values stay case-sensitive as in MDATK; project files and their comments are unchanged. *Tested:*
  project-model/validate_cases.json ("renderer", "renderer vulkan", "renderer case", "everything");
  tests/parity/improvements_cli.py applies the same rule to MDATK's validate().

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
- **CLI-07: selection errors** (`info`'s group entries, `analyze`, `export-pdb`) are SEL-19 ...
  SEL-22's: `invalid selection 'group:Broken' (expanded: '(name ()'): Selection failed: 'Unexpected
  token '(''`, `selection 'group:Target' matched no atoms: group 'Target' (residues '13-24') matches
  no atoms - the structure's residues are numbered 1-1`. *Tested:* cli_cases.json.
- **CLI-08: PDB export checks every frame first.** A frame with a coordinate beyond the PDB limits
  stops the export before any file is written: `PDB files must have coordinate values between
  -1000.000 and 9999.999 Angstroem: frame 2 is outside this range; no files were written.` (MDATK
  aborts at that frame and leaves the earlier files.) *Tested:* io/pdb/export/export_cases.json,
  cli_cases.json (export_coordinate_limits).
- **EXP-01: the PDB export's atom subset** (MDATK 0.0.10, io/pdb_export.py lines 19-40). A subset
  of commas only (`,`, ` , , `) is `',' names no residue group (leave the atom subset empty to export
  every atom)`, and a comma list of plain names that are not all groups names the first unknown one
  with the message `group:<name>` gets (`unknown residue group 'Nope' (defined: Aptamer, Target)`,
  with SEL-20's suggestion for a typo: `...; did you mean 'Target'?`);
  MDATK 0.0.4.1 parsed both as a selection (`Unknown selection token: 'Aptamer,'`).
  *Tested:* selections_cases.json (export cases, `Aptamer, Tagret` among them), with
  tests/parity/improvements_core_io.py.
- **EXT-01: two external series of one key keep both** (MDATK 0.1.0, io/external.py lines 111-121).
  Two files of one name (`run1/rmsd.dat`, `run2/rmsd.dat`) or two columns whose names read alike gave
  the same series key, and the later series replaced the earlier one. The later one now gets a key of
  its own, `rmsd.rmsd-2`, `-3` ..., labelled `rmsd: rmsd (2)` (the same column named twice still
  replaces). *Tested:* io/external/load_cases.json (duplicate_names, same_key_twice), with
  tests/parity/improvements_core_io.py.
- **A single-value option given twice takes the last value**, as argparse does (`--top`, `-o`,
  `--cache-dir`, `--frames`, `--project`, `--cpptraj`, `--conda-env`, `--only`, `--frame-list`,
  `--times`, `--selection`; a repeated `--traj` replaces the earlier list). This is MDATK's
  behaviour, listed because CLI11 rejects repeats by default. *Tested:* cli_cases.json
  (init_options_twice, analyze_options_twice, export_options_twice).
- JSON presets: `init` with a `.json` preset reports that MDBridge uses TOML project files.

## Demo, cleaning and citations commands (src/app/commands/demo.cpp, citations.cpp; src/analysis/cleaning.cpp; src/core/process.cpp)

`demo` is a port of MDATK's `mdatk/demo.py` and writes MDATK's files byte for byte for the same
arguments (demo.parm7, demo.nc; demo.mdatk.toml apart from `saved`): numpy's `default_rng(7)`
(normal, `shuffle`), numpy's evaluation order without fused multiply-adds, scipy's KD-tree queries
as exact searches. `clean` and the cleaning inside `analyze`, `render` and `export-pdb` are a port of
`mdatk/cleaning/cpptraj.py` and `conda.py`: the same cpptraj script (byte for byte apart from the
time stamp of its first line), folder (`<cache>/cleaned/<key>` or `cleaning.output_dir`), file
names, marker `MDATK_clean.json`, `conda run` command lines and messages, so MDATK and MDBridge reuse
each other's cleaned files. MDATK 0.1.1's demo.py, cpptraj.py and conda.py are unchanged from
0.0.4.1; its pipeline change is adopted (CMD-05). *Tested:* cli_cases.json (`demo_*`, `clean_*`,
`analyze_clean*`, `export_clean`, `demo_run` against the Qt program; tests/parity/cli_parity.py runs
a fake cpptraj and a fake `conda run`), tests/cli/test_cleaning.cpp.

- **CMD-01: `demo` refuses impossible arguments before writing anything.** `--n-frames 0` (MDATK:
  `error: [Errno 5] trying to go over trajectory limit` after writing demo.parm7), a negative count
  (`negative dimensions are not allowed`) and a time step that is not a positive finite number
  (MDATK writes a trajectory whose frames are all at 0 ps, go back in time or are at NaN ps, and a
  project with the PDB time `nanps`) give `the demo needs at least one frame (--n-frames 0)` and
  `the time between frames must be a positive number of ps (--dt 0.0)`; the library's sequence
  (not a command option) must have at least four of A, C, G, U (MDATK: KeyError or IndexError).
  Not fixed in MDATK 0.1.1. *Tested:* cli_cases.json (demo_no_frames, demo_negative_frames,
  demo_bad_dt; tests/parity/improvements_cli.py patches MDATK's build_demo the same way),
  `Demo.RefusesImpossibleArguments`.
- **CMD-02: the demo's hint names MDBridge:** `Next: mdbridge render -p d/demo.mdatk.toml   (or:
  mdbridge gui d/demo.mdatk.toml)`. The files keep MDATK's texts (the parm7 title `MDATK synthetic
  RNA aptamer/target demo`, the notes `... written by 'mdatk demo' ...`) so that they are MDATK's
  byte for byte. *Tested:* cli_cases.json (NAMING in cli_parity.py), `Demo.Command`.
- **CMD-03: cancelling cleaning stops cpptraj at once, with its children.** MDATK checks the cancel
  flag between reads of cpptraj's output, so a cancel waits until cpptraj next prints (it is silent
  for long stretches of a large trajectory), and `terminate()` ends `conda run` but not the cpptraj
  it started. MDBridge's process runner (src/core/process.hpp, also for the planner's scripts in
  milestone 7) polls every 50 ms, starts the program in its own process group (Windows: a job
  object) and ends the whole group: SIGTERM, SIGKILL after 10 s (Windows: at once, as Python's
  terminate()). *Tested:* `Process.CancelEndsTheProcessGroup`, `Cleaning.Cancel`.
- **CMD-04: cpptraj's output is decoded line by line.** MDATK decodes each 4 KiB chunk on its own,
  so a UTF-8 character split between two reads becomes two U+FFFD; MDBridge decodes whole lines
  (invalid bytes still become U+FFFD) and writes cpptraj.log with the bytes cpptraj wrote (MDATK:
  the decoded text; the same file for UTF-8 output). *Tested:* `Process.LineSplitter`.
- **CMD-05: cleaned files are matched by their key everywhere** (MDATK 0.1.0, pipeline.py
  `cleaned_if_available`): previews and the interface's structure view use cleaned files only when
  the marker names the key of the current inputs and settings, also in a folder given as
  `cleaning.output_dir` (MDATK 0.0.4.1 used any cleaned files found there). Previews never start
  cpptraj. *Tested:* `Cleaning.OutputDirNeedsTheKey`, `ProjectView.CleaningWithoutCleanedFilesShowsTheRawFiles`.
- **CMD-06: `mdbridge citations [common options] [--bibtex FILE] [--json]`** (new; owner decision
  D14, AP-16; specification on MDATK's functions: tests/parity/citations_command.py). The project:
  `-p` and the overrides of every analysis command. Cited: the settings' resources (as MDATK's About
  dialog: `project_resources`, `keys_for`) followed by the run-time references of the cached results
  (meta `citation_keys` of the analyses the cache manifest's last run lists; none when one of them
  is missing; CIT-06, `core::citations::withRunKeys`). Text (default): MDATK's About text,
  `Resources this project uses (N references; each generated script lists only its own):`, a blank
  line, per resource `<name>[ - <note>]` and its references indented by four spaces and a blank
  line; then, when the cache adds references, `Recorded by the cached results of the last run:` and
  those references; without any: `This project uses no external resource with a reference yet.`
  `--json`: `{"keys", "settings_keys", "run_keys" (the cache's, unfiltered), "resources": [{"id",
  "name", "note", "keys"}]}` (indent 2, non-ASCII as it is). `--bibtex FILE` writes the BibTeX
  entries and prints `Wrote N references to FILE` (with `--json`, the JSON instead); nothing to cite:
  `Nothing to cite yet: The project uses no external resource with a reference.` and no file.
  *Tested:* cli_cases.json (`citations_*`), `Citations.SettingsAndRunKeys`, `Citations.Command`.
- **CMD-07: programs started from the AppImage get the user's environment** (as MDATK's frozen
  build: `frozen.child_env`): `LD_LIBRARY_PATH` and `PATH` entries inside `$APPDIR` and a
  `FONTCONFIG_FILE` there are left out, so cpptraj and conda do not load the AppImage's libraries.
  *Tested:* `Process.ChildEnvironmentLeavesTheAppImageOut`.
- The parameter log's `[run.inputs.cleaning]` holds `key`, `script`, `log` and `reused`, as MDATK's
  (MDBridge wrote `key` and `reused = true` alone while it could only reuse MDATK's cleaned files).

## Citations (src/core/citations.cpp, src/core/citation_data.inc)

Owner decision (2026-11-25): a method is cited by the people who originated it, not by one
implementation of it; an implementation of common knowledge, or of a method published by others,
is not a reason to cite the implementation. Every method MDBridge took over from MDAnalysis, NumPy,
SciPy and Matplotlib through MDATK was checked against that rule (the inventory is in
`src/core/include/mdbridge/core/citations.hpp` and below).

- **CIT-01: the software MDATK ran on is not cited; the RMSD's superposition methods are.**
  MDATK's `pipeline_resources` cites MDAnalysis (Michaud-Agrawal 2011, Gowers 2016), NumPy (Harris
  2020), SciPy (Virtanen 2020) and Matplotlib (Hunter 2007) for every project, because MDATK ran on
  them. MDBridge does not use these packages, so they are not cited and their resources and
  references are no longer in the table. Of the methods MDBridge reproduces, two decide a reported
  result and have a published originator outside those projects; they are cited when the project
  uses them (`analyses.rmsd` enabled with targets and `superposition`): resource `qcp`
  (Theobald 2005, doi:10.1107/S0108767305015266; Liu, Agrafiotis and Theobald 2010,
  doi:10.1002/jcc.21439) for a target superposed without a fit selection (MDAnalysis' `rms.rmsd`,
  which implements QCP; `rmsd_qcp.cpp`), and resource `kabsch` (Kabsch 1976,
  doi:10.1107/S0567739476001873; Kabsch 1978, doi:10.1107/S0567739478001680) for a target
  superposed on a fit selection (`kabsch` of `spatial.cpp`). Not cited, as common knowledge or
  standard algorithms whose implementation does not change a result: distance and neighbour searches
  (cell grids, k-d trees), minimum-image distances in triclinic cells, the geometric criteria of
  the H-bond, contact, stacking and Mg analyses, plane fits, bond guessing by distance, pairwise
  summation, linear interpolation, Gaussian and uniform filters and trilinear resampling of the
  surface and density shadow, axis ticks and figure layout; the selection language (MDAnalysis'
  syntax, in the tradition of CHARMM's and VMD's) is an interface convention, not a method. Colour
  maps (viridis family, cividis, ColorBrewer, Tableau) are styling and not cited (owner decision);
  their licences stay in `packaging/licenses`. The random numbers (PCG64, SeedSequence, ziggurat)
  only make the example data of the layout preview and are not cited. *Why:* citations should
  credit the originators of what a result depends on; citing packages MDBridge does not run
  misstates how a result was produced. *Tested:* citation_cases.json (`references`, `resources`,
  `pipeline` with RMSD variants, `project`; tables and functions from MDATK with
  tests/parity/improvements_citations.py applied, which tools/gen_citation_data.py applies to
  src/core/citation_data.inc too), cli_cases.json and render/cli/render_cli_cases.json (params and
  step logs of `render` / `charts`), `Citations.MethodsByTheirOriginators`.
- **CIT-02: Barnaba is not cited.** MDATK cites Barnaba (Bottaro 2019) when
  `analyses.basepairs.engine = "barnaba"` and a secondary-structure panel reads dynamic pairs.
  MDBridge never runs Barnaba (see "No Barnaba" above: its builtin detector runs instead), so it is
  not cited; the table keeps its entry. *Tested:* citation_cases.json (`pipeline`, base-pair
  engines), `Citations.MethodsByTheirOriginators`.
- **CIT-03: the bond-guess warning cites the covalent radii.** When the energy analysis excludes
  bonded pairs along guessed bonds (AN-07), its warning names the source of the radii the guess
  uses (AN-14): `WARNING: energy: the topology has no bonds; bonded pairs are excluded along bonds
  guessed from distances on the first analysed frame (covalent radii: Cordero et al., Dalton Trans.
  2008, 2832-2838, doi:10.1039/b801115j)` (owner decision). *Why:* the radii decide which pairs
  are excluded and so the reported energies; they are Cordero et al.'s data, not common knowledge.
  The method (a distance cut-off over summed radii) stays uncited (CIT-01); since CIT-06 the radii
  are in the citation table and a run that guessed bonds cites them in the project's citations. *Tested:* micro_cases.json and systems_cases.json energy logs
  (tests/parity/improvements_analyses.py, AN-07, records the message).
- **CIT-04: the energy analysis' parameters without parm7 are cited by what was used.** MDATK's
  per-element table (`GENERIC_LJ`: Amber parm99 / parm10 types, Joung-Cheatham ions cut to four
  decimals, and Ca (1.71, 0.12), Mn (1.2, 0.05) and a default (1.9, 0.1) that match no published
  set) was an uncited ad-hoc mix. MDBridge replaces it by published parameters (AN-25 - AN-27,
  under "Analyses") and cites their originators when a run uses them, by what it used: the
  analysis' log line (`energy: Lennard-Jones parameters without parm7: N atoms from Amber
  force-field residue templates (OL3, TIP3P), M atoms from Amber atom types by element; partial
  charges ... Cite: Cornell et al. 1995, doi:...; ...`) and its metadata `lj_references` (citation
  keys in first-use order). Owner policy: cite who originated the numbers actually used. The
  analysis uses only charges and Lennard-Jones parameters, so: the protein and nucleic templates
  (ff14SB, OL15, OL3) cite Cornell et al. 1995 alone, since their charges are ff94's and every
  Lennard-Jones class in the three files has the R* / ε of a parm94.dat type (checked class by
  class against parm94.dat, SHA-256 `34021e968833329f...` of ParmEd's copy; parm99 and parm10 kept
  those values); the torsion refinements that name these force fields (ff14SB, bsc0, χOL3, OL15
  and their parts) do not enter the energies and are not cited for them. TIP3P water: Jorgensen et
  al. 1983; ions: Joung and Cheatham 2008 (monovalent), Li et al. 2013 (divalent, 12-6 CM set), Li
  et al. 2015 (3+ / 4+, IOD set). The typed layer: parm94's types Cornell et al. 1995, Cl and Br
  (parm99) Wang et al. 2000, the TIP3P water oxygen Jorgensen et al. 1983, the ions Joung and
  Cheatham 2008 / Li et al. 2013; the extra point (no Lennard-Jones term) nobody. Li et al. 2015
  joins the citation table (tests/parity/improvements_citations.py, `verified` from Crossref); the
  others were in MDATK's and were checked against Crossref too. Not traced further: the charges of
  the residues Amber added after ff94 (ASH, GLH, LYN, CYM, HYP), which were fitted later in ff94's
  scheme; they are cited as ff94's. No project resource is added: which force fields a run uses
  depends on the residues it meets, which the settings do not tell, so `pipelineResources` stays as
  it was and the run names them. *Why:* owner policy, cite the originators of what a result
  depends on, only when it is used. *Tested:* `EnergyParameters.*` (test_energy.cpp), the energy
  logs and metadata of micro_cases.json, systems_cases.json, demo_dry_variants.json,
  pipeline_cases.json and cli_cases.json (improvements_analyses.py), citation_cases.json.
- **CIT-05: the `atoms` analysis' bond-guess warning cites the covalent radii.** When the topology
  has no bonds among the structure atoms, the `atoms` analysis (structural view) guesses them
  (AN-13, AN-14) and now says so, once per run: `WARNING: structural view: the topology has no bonds
  among the structure atoms; bonds are guessed from distances on the first analysed frame (covalent
  radii: Cordero et al., Dalton Trans. 2008, 2832-2838, doi:10.1039/b801115j)` (MDATK: silent), as
  the energy analysis' warning does (CIT-03). *Why:* the radii decide the drawn bonds; they are
  Cordero et al.'s data. *Tested:* `AnalysisImprovements.AtomsBondGuessWarningCitesCovalentRadii`,
  the `atoms` logs of systems_cases.json (improvements_analyses.py).
- **CIT-06: run-time citations reach the project's citations.** What a run used that its settings
  do not tell is recorded in the result's metadata `citation_keys` (citation keys in first-use
  order): the energy analysis' `lj_references`, then Cordero et al. 2008 when it guessed bonds
  (exclusions, AN-07, or the parameters' chain ends and neighbours, AN-25 / AN-26); the `atoms`
  analysis' Cordero et al. 2008 when it guessed its bonds (CIT-05). The parameter log
  (`citation_keys`, `citations`) and the render step's entry in the project's step log cite the
  settings' references (as before) followed by those the cached results of the run record
  (`analysis::Pipeline::runCitationKeys`: the entries of the run's analysis keys, in that order,
  each entry's keys in its order; `core::citations::withRunKeys`: without repeats and unknown
  keys). Without cached results (no analysis ran or the entries are gone; results cached before
  this change have no `citation_keys`) only the settings' citations appear. The About dialog, the
  BibTeX export and generated script headers are not in MDBridge yet; they are to cite through
  `withRunKeys` too. Cordero et al. 2008 joins the citation table (Crossref-verified). *Why:*
  owner decision: cite what a result depends on wherever the project is cited, not only in a log
  line. *Tested:* `Pipeline.RunCitationKeysFromTheCache`, the energy and atoms metadata of the
  analysis cases, the parameter and step logs of cli_cases.json (improvements_citations.py
  `_patch_run_citations`), citation_cases.json.
- VTK is not cited either (RN-16, under "Rendering: the render, charts and preview commands").

## Analyses (non-periodic improvements)

Expected data: `tests/parity/improvements_analyses.py` patches MDATK (and the RMSD probes) to the
behaviour below before `make_analysis_cases.py`, `make_analysis_framework_cases.py`,
`make_pipeline_cases.py`, `make_rmsd_cases.py` and `cli_parity.py --write` record it; hand-made
checks are in tests/analysis/test_improvements_analyses.cpp. Version bumps (new cache keys): energy 2,
hbonds 2, rmsd 2, proximity 2, atoms 3 ("Analysis cache" above has every analysis' version).

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
  <file> is not a parm7 file; Lennard-Jones parameters come from force-field residue templates and
  Amber atom types` (MDATK: silent). The warning when a parm7 file's parameters cannot be read ends
  the same way (MDATK: `...; using generic element parameters`), naming the assignment of AN-25 /
  AN-26 instead of MDATK's per-element table. *Tested:* systems_cases
  logs (pdb / gro / noh / aromatic), test_energy.cpp.
- **Energy exclusions without topology bonds (AN-07).** With `exclude_bonded > 0` and no bonds:
  `WARNING: energy: the topology has no bonds; bonded pairs are excluded along bonds guessed from
  distances on the first analysed frame (covalent radii: Cordero et al., Dalton Trans. 2008,
  2832-2838, doi:10.1039/b801115j)` (the citation: CIT-03), and the exclusions follow the bonds the `atoms` guesser
  (below) finds on the whole system (MDATK: nothing excluded, silently). *Tested:* micro and
  systems_cases energy logs and values.
- **Exclusions along any path (AN-08).** Pairs joined by at most `exclude_bonded` bonds through any
  atoms are excluded (MDATK: only paths that start at a direct A-B bond and run through A). *Tested:*
  `ExclusionsThroughOtherAtoms`; the recorded chains give the same keys.
- **Energy parameters without parm7: residue templates (AN-25).** Without a parm7 topology (or with
  `use_parm7 = false`, or a parm7 file that cannot be read) MDATK gave every atom one Lennard-Jones
  (R*, ε) per element (`GENERIC_LJ`, CIT-04): every hydrogen Amber's `HP` (1.1 Å), which made
  hydrogen-bonded interfaces strongly repulsive (MDATK's demo duplex, strand against strand, first
  frame: vdW +145.4 kcal/mol with the table against -34.7 with the topology's parameters).
  Measured on the generated `noh` system (the demo's heavy atoms, chain IDs, no bonds): 468 of
  511 atoms take template parameters and 43 the typed layer (the two 5' residues, which keep their
  phosphate and so match no Amber template); with chain ends decided by bonds alone it was 221 and
  290. MDBridge
  matches each residue of the analysed atoms against the residue templates of Amber's force fields
  as OpenMM 8.1.2 converts them (amber14: protein.ff14SB, DNA.OL15, RNA.OL3, tip3p with the
  Joung-Cheatham and Li-Merz TIP3P ions; MIT; data/sources/forcefields, tools/gen_energy_lj.py ->
  src/analysis/energy_lj_tables.inc, never typed by hand). Candidates: the template named like the
  residue (upper case), with its terminal variants when it is the middle form (Amber's N/C-prefixed
  amino acids, 5/3/N-suffixed nucleotides); a name that is not a template stands for the templates
  pdbNames.xml makes it an alias of (HIS, HSD -> HID / HIE / HIP; WAT, SOL, TIP3 -> HOH; ADE -> A
  ...); a one-atom residue is an ion template named like the residue or the atom (as written, upper
  case, without charge signs and digits: `Na+` -> NA) of the same element. An atom maps to the
  template atom of its name, else by pdbNames.xml's aliases (residue-specific first, then the
  Protein / Nucleic / HOH ones: `O1P` -> OP1, `HO5'` <- `H5T`, ...), one template atom per atom, and
  its element (when known) must be the template atom's. A residue *matches* a template when every
  atom maps; *completely* when every template atom is mapped too (hydrogens included). Exactly one
  complete match: that template, whatever the residue's position. Otherwise the templates that
  contain the residue are narrowed by its position in its chain when it has terminal variants: it
  starts a chain when it does not continue the residue before it (by index), ends one when the
  next does not continue it. Residue r + 1 continues residue r unless the evidence agrees on a
  break: a new chain (another segment, or another chain ID where the topology has chain IDs), or a
  gap or reset in the numbering (resid(r + 1) is neither resid(r) + 1 nor resid(r) with another
  insertion code) *and* no bond between the two (the topology's bonds where either residue has
  any, else those AN-14's guess finds on the first analysed frame). So consecutive numbering in
  one chain is never split by geometry (MDATK's synthetic demo has O3'-P distances of 0.8-3.9 Å),
  and a bond bridges a numbering gap. Without chain IDs (GRO) and with numbering that runs on
  across strands, two strands are not told apart this way; their ends still match by their atoms
  when hydrogens are present. Start -> N / 5', end -> C / 3', both -> the nucleic N form (no amino-acid
  template; such a residue falls through), neither -> the middle form. If the remaining templates
  give every atom the same R* and ε (HID / HIE / HIP without hydrogens do), the residue takes them;
  else (none left, or they disagree, e.g. a 5'-terminal nucleotide with its phosphate) it falls to
  AN-26. R* = σ 2^(1/6) / 2 and ε / 4.184 from OpenMM's σ (nm) and ε (kJ/mol): Amber's values to
  about 1e-15 relative (P: 2.1000000000000005). Lorentz-Berthelot combination, as before. The water
  model is unknown without a topology that says it: TIP3P's are assumed throughout, for the water
  template (O 1.7683 / 0.152, charges -0.834 / 0.417), the ion templates and the ion and `OW`
  types of AN-26 (the TIP3P sets of Joung and Cheatham and of Li et al.; their SPC/E or TIP4P-Ew
  sets differ, e.g. Na+ 1.212 / 0.3526 for SPC/E). A 4-site water (TIP4P, OPC) matches the TIP3P
  template by its O and H (its extra point takes AN-26's `EP`) and takes no charges.
  *Why:* the force field's own parameters per atom type reproduce Amber's energies (demo duplex:
  -34.62 kcal/mol from the templates with the parm7 file's charges, -34.63 from a PDB of the
  same frame, against -34.67 with the parm7 file's own parameters, whose base hydrogens and
  hydroxyl hydrogens use simplified types). *Tested:* `EnergyParameters.DemoMatchesTheTemplates`,
  `IncompleteResiduesAndTermini`, `ChainEndsFromChainIdsNumberingAndBonds`, the energy cases of micro_cases.json, systems_cases.json,
  demo_dry_variants.json (energy_generic_lj), pipeline_cases.json and cli_cases.json
  (tests/parity/energy_parameters.py, the same rules on MDATK's side).
- **Energy parameters without parm7: Amber atom types (AN-26).** An atom of a residue no template
  matches (ligands, modified nucleotides, other residues, ions with other names) takes the values
  of one Amber atom type, the design of MDATK's later generic tables (same owner; the values read
  from parm10.dat, frcmod.ionsjc_tip3p and frcmod.ions234lm_126_tip3p, whose SHA-256 match those
  that design records; data/sources/forcefields/amber_lj_types.tsv, tools/extract_amber_lj.py): a
  hydrogen by the element of the atom it is bonded to (topology bond, else the nearest non-hydrogen
  atom of its residue within 1.45 Å on the first analysed frame; none: as on nitrogen): on N `H`
  (0.6, 0.0157), on O `HO` (0, 0), on S `HS` (0.6, 0.0157), otherwise `H1` (1.387, 0.0157); an atom
  alone in its residue whose element has a TIP3P ion set: that ion (Joung-Cheatham Li+ ... I-,
  Li-Merz 12-6 CM Mg2+, Ca2+, Sr2+, Ba2+, Mn2+, Fe2+, Co2+, Ni2+, Cu2+, Zn2+, Cd2+); C, N, O, P,
  S, F, Cl, Br, I: `C`, `N`, `O`, `P`, `S`, `F`, `Cl`, `Br`, `I` of parm10.dat; a metal of that ion
  list inside a larger residue: its ion's; any other element: carbon's, with `WARNING: energy: no
  Lennard-Jones parameters for element(s) 'X' (N atoms); carbon's (Amber type C) are used`
  (metadata `lj_default_elements`). Extra points (virtual sites, lone pairs) take parm10.dat's `EP`
  (0, 0: no Lennard-Jones term), are no atom's bonded neighbour, are left out of template matching
  (their residue then takes no template charges) and are not reported as unknown elements: an atom
  whose element is no chemical element (CIAAW symbols) and whose name is a water virtual site of
  OpenMM's force fields (M, M1, M2, LP1, LP2, OM), GROMACS' TIP4P `MW`, or starts with `EP` or `LP`
  (Amber's `EPW`, `EP1`, lone pairs), and the virtual-site names of OpenMM's water residues in them
  whatever their guessed element (`OM` of CHARMM's `TP4E`). Refinements of that design, where the bonded neighbours are
  known (the topology's bonds where the residue has any, else AN-14's guess on the first analysed
  frame), all parm10 types: a carbon with four neighbours `CT` (1.908, 0.1094; sp3, as Amber types
  them) instead of `C` (1.908, 0.086); an oxygen with two neighbours `OS` (1.6837, 0.17; ether,
  ester, phosphate ester), `OH` (1.721, 0.2104) when one of them is a hydrogen, `OW` (1.7683, 0.152;
  TIP3P's) when both are. Without hydrogens a carbon rarely has four neighbours, so it stays `C`.
  No UFF or other generic force field (owner decision). *Why:* MDATK's table gave every hydrogen
  one radius and two elements values of unknown origin. Measured on the demo duplex (strand against
  strand, first frame, vdW, kcal/mol): -34.74 from the types alone (every residue renamed),
  -34.67 parm7, +145.4 MDATK's table. *Tested:* `EnergyParameters.TypedLayer`, `ExtraPoints`, the systems cases.
- **Energy charges without topology charges (AN-27).** MDATK computed no electrostatics without
  topology charges (PDB, GRO). MDBridge gives the atoms of a residue that matches a template
  *completely* (AN-25, hydrogens included) the template's partial charges; a residue with a missing
  or extra atom gets none (never part of a residue's charge set), and the run says so: `WARNING:
  energy: the topology has no partial charges; K of N atoms take the charges of the force-field
  residue template they match completely (hydrogens included); the other M atoms (R residues) have
  no charge, so their electrostatic terms are zero` (all or none matched: no warning, or MDATK's
  `topology has no partial charges; electrostatics will be zero`). Topology charges (parm7, PSF,
  ...) take precedence for every atom and are never mixed with template charges: the formats store
  a charge for every atom or for none (a 0 is a charge), and two charge models in one sum would be
  inconsistent. Metadata: `charge_source` `topology`, `templates` (every atom), `partial` or
  `none`; `charge_atoms` {topology, template, none}, `lj_atoms` {parm7, template, typed} (atoms of
  both groups), `lj_source` `parm7`, `templates`, `typed` or `templates+typed` (MDATK: `generic`).
  *Why:* a PDB or GRO of a standard system then gives the force field's electrostatics (demo duplex,
  PDB of the first frame: -48.36 kcal/mol from OL3's charges; the demo's own parm7 has synthetic
  charges, -43.76). *Tested:* `EnergyParameters.AnalysisChargesMetadataAndLog`, the systems and
  pipeline cases.
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

## Fixes adopted from MDATK 0.1.x (analysis)

MDATK fixed these after 0.0.4.1; MDBridge takes MDATK's design (same owner; its code ported), cited per
entry by version and file. Expected data: MDATK 0.0.4.1 patched to the same behaviour on the test side
(tests/parity/improvements_analyses.py, improvements_readers_amber.py, improvements_readers_pdb.py,
improvements_citations.py); the parity tests keep 0.0.4.1 as the reference. No analysis version is bumped
for them: MDBridge has not been released, and hbonds, contacts and trace already carry this release's one
bump ("Analysis cache" above).

- **H-bonds by distance only count each atom pair once.** Without the angle (no hydrogens in the
  topology, or `use_angle = false`) MDATK counted a donor-acceptor pair once per direction when both atoms
  can donate and accept (every N / O of a heavy-atom PDB, two hydroxyls), inflating the count, maps and
  per-residue series by up to 2x. Now a pair already counted with the donor in group A is not counted
  again with the donor in group B (MDATK 0.1.0, `mdatk/analysis/hbonds.py` process(); MDBridge
  src/analysis/hbonds.cpp process()). Distance-only runs with an angle are unchanged. Measured: the
  `noh` system (the demo's heavy atoms, no hydrogens; systems_cases.json) 114, 114, 72, 108 -> 57, 57, 36,
  54 H-bonds over its four frames (every N / O there both donates and accepts, so every pair was counted
  twice). *Tested:* systems_cases.json `noh` and `noh_no_angle`, the only cases whose results changed.
- **Contacts of a group with itself.** (MDATK 0.1.0, `mdatk/analysis/contacts.py` version 2; MDBridge
  src/analysis/contacts.cpp; register AN-23, which kept these as conventions.) (1) The map holds each
  residue pair once (i < j), and `per_res_a` / `per_res_b` counted it on one residue only, so a residue's
  count depended on index order; now both arrays are each residue's row plus its column (equal). (2)
  Residues next to each other in sequence (`|resindex difference| <= 1`) are left out only within one
  chain (`core::chainIds`: chain ID, segment, then bonded fragments when they hold every residue
  together, TOP-01): the last residue of a strand and
  the first of the next (a duplex's closing base pair) are no neighbours and count. (3) The `atom_pairs`
  series is labelled `A–B atom contacts`, not `heavy-atom contacts`, when hydrogens count
  (`heavy_only = false`); the output label (`outputLabel`, chart and panel titles) says the same. *Why:*
  bugs: order-dependent counts, missing contacts of the closing pair, a wrong label. *Tested:*
  demo_dry_variants.json `contacts_self_duplex` (the demo duplex against itself: per_res_a equals
  per_res_b; the closing pair 12-13, 14 atom contacts on the first frame, was left out before), the contacts_same_group / contacts_all_atoms cases of
  micro_cases.json and demo_dry_variants.json; systems_cases.json `noh_conect_self` (chains from
  TOP-01); `AnalysisFramework.InstanceKeysAndLabels`.
- **Structure trace cache key covers every group.** The trace result records every residue group's
  residues (`group_member`, `group_index`), which colour the view and choose the residues of the 3D
  representations, but its cache key covered only the groups shown (`visual.structure_groups`) and the
  Mg selection, so editing a group that is not shown reused a stale trace. Its selection texts (and so
  the key's `groups`) are now every group plus the Mg selection (MDATK 0.1.0,
  `mdatk/analysis/structure.py` TraceAnalysis.selection_texts; MDBridge src/analysis/structure.cpp
  TraceDescriptor). *Tested:* framework_cases.json (selection texts and keys), pipeline_cases.json and
  cli_cases.json cache keys.
- **parm7 LJ reader: CHAMBER and CMAP format lines.** MDATK's parm7 section reader (which reads the
  Lennard-Jones table of the energy analysis) required a repeat count and no inner parentheses in
  `%FORMAT`, and raised `cannot parse format line` for CHAMBER's `(a80)` (CTITLE) and `(i2,a78)`
  (FORCE_FIELD_TYPE) and ff19SB's `(8(F9.5))` (CMAP parameters), before checking whether the section was
  wanted; the energy analysis then warned and fell back to other parameters. Now `(a80)` reads with a
  repeat count of 1, `(8(F9.5))` as `8F9.5`, and a section whose format still cannot be parsed is
  skipped when it is not wanted (MDATK 0.1.0, `mdatk/io/parm7.py` `_FMT` and read_sections; MDBridge
  src/io/parm7.cpp searchFormat / readParm7Sections). The topology reader already skipped such
  sections. *Tested:* read_cases.json `format_no_count`, `format_inner_parens`,
  `format_unparsed_unwanted`, `format_unparsed_wanted` (tests/parity/make_parm7_cases.py),
  `Energy.Parm7WithChamberAndCmapFormatLines` (the demo's parm7 with CHAMBER and CMAP sections gives the
  same energies, `lj_source` parm7).
- **Energy metadata `lj_table`; citation resource `generic_lj`.** A run without parm7 parameters records
  `lj_table`, MDATK 0.1.1's key (`mdatk/analysis/energy.py` finalize, `generic_lj.py` TABLE_NAME): the
  name of the table its parameters come from, here `Amber force-field residue templates (OpenMM amber14:
  ff14SB, OL15, OL3, TIP3P with Joung-Cheatham and Li-Merz ions), then Amber parm10 atom types,
  Joung-Cheatham and Li-Merz (CM) ions for TIP3P (docs/GENERIC_LJ.md)`, after `n_excluded_pairs` as in
  MDATK; MDBridge's own `lj_source`, `lj_atoms`, `lj_references` and `lj_default_elements` stay (AN-25 -
  AN-27, CIT-04). `docs/GENERIC_LJ.md` records both layers as MDATK 0.1.1's `docs/GENERIC_LJ.md` records
  its table: rules, every number with its file, SHA-256 and publication, and a test
  (`EnergyParameters.GenericLjDocument`) that keeps the document equal to the extracts in
  data/sources/forcefields. MDATK 0.1.1's citation resource `generic_lj` (`mdatk/citations.py`,
  `_citation_data.py`) is in the citation table, so its params logs and project citations resolve; its
  keys are those of MDBridge's typed layer (Cornell 1995, Wang 2000, Jorgensen 1983, Joung 2008, Li
  2013: MDATK's four plus TIP3P's water oxygen `OW`). It is never a project resource: MDBridge cites the
  energy parameters per run (CIT-04, CIT-06). The values differ from MDATK 0.1.1's for every run without
  parm7 (demo duplex, strand against strand: -34.62 kcal/mol from the templates, MDATK 0.1.1 -34.1,
  parm7 -34.67). *Tested:* `EnergyParameters.AnalysisChargesMetadataAndLog`, the energy metadata of the
  analysis, pipeline and cli cases, citation_cases.json, `Citations.MethodsByTheirOriginators`.

## Nucleotide names: DNA, L-DNA and CCD components (src/io/pdb.cpp, src/chem/secstruct.cpp)

- **NA-01 L-DNA names, the CCD's own atom names and the backbone guard** (follows MDATK 0.1.1
  `mdatk/pdbtools.py` `parse_residue_name`, `has_backbone_atoms`, `standard_atom_name`, and their
  readers `secstruct.PairDetector`, `analysis/basepairs.py`, `analysis/stacking.py` `ring_atoms`,
  `analysis/structure.py`'s trace flag `backbone`, `render/molecule/topology.py`,
  `render/residues.py` `one_letter` and the sequence panel; MDATK 0.0.13 - 0.0.17). MDATK 0.0.4.1
  read D-DNA (DA, DA5, DA3, DAN, THY ...) but no L-DNA, took any residue named like a nucleotide for
  one, and looked for the standard atom names only. MDBridge:
  - reads the L-DNA names: the CCD's L-deoxyribonucleotides 0DA 0DC 0DG 0DT, and MDATK's own LDA LDC
    LDG LDT (LDA5 / LDA3 / LDAN in older files) as L-DNA A C G T (`io::parseResidueName`; the
    sequence panel, base pairs and the molecule view's base colours). MDATK 0.0.4.1: `X` / no base;
  - takes a residue for a nucleotide by its name only when it has a backbone atom (a primed name,
    an old-style `*` or P; `io::hasBackboneAtoms`): in the CCD, LA and RA are the lanthanum and
    radium ions and LDA a detergent, which neither pair (the base-pair detector) nor show as `A` in
    the sequence panel (through the trace's new per-row array `backbone`; a trace cached without it
    reads the names alone, as 0.1.x does);
  - reads the atoms that the CCD names its own way in 0DA (adenine N9A, C8A, N7A, C5A, C6A, N6A, N1A,
    C2A, N3A, C4A; H8A, H2A, HN61, HN62; H2P on a protonated OP2) and 0DT (C5M, the methyl carbon) by
    their standard names (N9 ... C7; `io::standardAtomName`) in the base-pair detector, the stacking
    rings and the molecule topology, so their pairs, rings and bonds are found. An atom of the
    standard name wins over a renamed one in the detector.

  The tables come from the CCD (data/sources/ccd/l_deoxyribonucleotide_atoms.tsv, cut by
  `tools/extract_ccd.py`: the released L-DNA linking components whose atoms and bonds match DA, DC,
  DG or DT, each atom with the matching standard atom by a graph match that keeps names where the
  bonds allow), not from MDATK's or MDAnalysis' tables; the result is the same as 0.1.x's
  `CCD_ATOM_NAMES`. Not adopted here: 0.1.x's PDB reading and writing of the components' names
  (`read_pdb` / `write_pdb`), its naming schemes and the `pdb clean` changes (the PDB tools). The
  cache versions of base pairs (2), stacking (2) and the trace (5) already differ from MDATK's
  ("Analysis cache"). *Tested:* chem/secstruct_cases.json (`residue_names`; detector cases `L-DNA
  0DA with CCD atom names`, `LDA without a backbone`, `LA without a backbone`), io/pdb/pdb_cases.json
  (`residue_names`), render/shading one_letter, `Pdb.BackboneAtomsAndComponentAtomNames`,
  `Shading.OneLetterNeedsABackbone`; the trace's `backbone` array in every analysis case. Expected
  data: tests/parity/improvements_nucleotides.py.

## Rendering (milestone 4): overview

Milestone 4 ports MDATK's presentation layer (`mdatk/render/`: figures, panels, charts, the
structure view and its 3D representations, the video, the `render` / `charts` / `preview` commands)
and Matplotlib's colour, tick and layout rules. Figures are compared by their contents, not their
pixels (owner decision D10): the expected data of `tests/data/render` records what MDATK's figures
contain (data, colours, ticks, labels, legends, images, colour bars, layout boxes) and what its
software renderer, molecule geometry, density shadow and commands produce, and the C++ tests compare
MDBridge's with it; there are no golden reference images. The sections below list the deliberate
differences. The plan's improvement items (docs/specs/m4-plan.md section 4, with their status) are
here:

| Items | Section |
|---|---|
| RN-01, RN-04, RN-14 | Rendering: plot panels (charts: RN-01, RN-14 too) |
| RN-02 | Rendering: structure view and sequence panels (the panels' `No data: no frames`), the video frame loop, the render commands |
| RN-03, RN-08, RN-09, RN-33 | Rendering: structure view |
| RN-05, RN-10, RN-11, RN-50 | Rendering: charts (RN-11 also: the render commands) |
| RN-06 | Rendering: the figure framework, the video frame loop |
| RN-07 | Rendering: video encoding, the render commands |
| RN-12, RN-13, RN-18 | Rendering: colours and themes (RN-12, RN-13 also: the render commands) |
| RN-15 | Rendering: the render commands |
| RN-16 | Rendering: structure view, the Vulkan 3D renderer, the render commands |
| RN-17 | kept as MDATK has it (Rendering: the video frame loop) |
| RN-19 | Rendering: molecule geometry |
| RN-20 | Project files |
| RN-21 ... RN-25 | Rendering: residue shading, orientation and the example dataset |
| RN-26 ... RN-29 | Rendering: software 3D renderer |
| RN-30 ... RN-32 | Rendering: density shadow |
| RN-40 ... RN-42 | Rendering: sequence and secondary-structure panels |
| D5, D6, D8, D9, D11, D13 (owner decisions) | charts (D5), colours (D6), the render commands (D8, D9), structure view (D11), interface (D13) |
| D15, D16 (owner decisions) | the Vulkan 3D renderer (D15), video encoding (D16) |
| CIT-01 ... CIT-06 | Citations (CIT-03 also: Analyses, AN-07; CIT-04 also: Analyses, AN-25 ... AN-27) |

## Rendering: colours and themes (src/render/core)

Colour strings follow Matplotlib 3.11.2's `to_rgba` grammar exactly (names in its priority order and
case rules, `C0`-`C9`, hex of 3 / 4 / 6 / 8 digits, grey levels in Python's `float()` syntax, its
error messages); colormaps, `with_extremes`, the chain palettes and MDATK's themes give Matplotlib's
and MDATK's numbers bit for bit. The expected data is Matplotlib's and MDATK's
(`tests/parity/make_color_cases.py`), with `tests/parity/improvements_colors.py` for the entries
below. *Tested:* `tests/render/test_colors.cpp` (`render/colors/color_cases.json`).

- **Ordinal colours for more than five classes (RN-18).** MDATK's `theme.ordinal(n)` kept its four
  ramp colours when a class matrix had more than five classes, so with the panel's `vmin -0.5` /
  `vmax n - 0.5` several classes shared one colour while the colour bar labelled them apart.
  MDBridge gives classes 1 .. n-1 the n - 1 evenly spaced steps of the ramp
  (`LinearSegmentedColormap.from_list(ramp, N = n - 1)`); up to five classes (every class matrix
  MDATK 0.0.4.1 makes has three) the colours are MDATK's. *Tested:* `Theme.Ramps` (n up to 13, every
  class a distinct colour).
- **NaN colour components (RN-12, part).** An RGB(A) tuple with a NaN component is rejected
  (`RGBA values should be within 0-1 range`); Matplotlib's range test lets NaN through, which drew
  nothing or black later. `colorProblem()` / `isColorLike()` give the message for `check()` (which
  names the project field, W3-B). *Tested:* `Color.ToRgbaTuples`.
- **Bundled colormaps only (RN-13, owner decision D6).** `visual.chain_palette` and the other map
  names accept the qualitative maps tab10, tab20, tab20b, tab20c, Set1-3, Pastel1-2, Dark2, Accent,
  Paired and the perceptual viridis, plasma, inferno, magma, cividis, each also reversed (`_r`), as
  in Matplotlib's registry. Any other name gives `'jet' is not a bundled colormap; choose from tab10,
  ... (or one of them reversed, with "_r")`, with `; did you mean 'Set2'?` for a name that differs in
  case only (MDATK: Matplotlib's registry error in the middle of the render, or its other maps).
  *Tested:* `Colormap.Names`.
- **Invalid colormap data is reported when the map is made.** Matplotlib builds the lookup table
  lazily and raised `data mapping points must start with x=0 and end with x=1` (a single colour in
  `from_list`) at the first use; MDBridge raises the same message when the map is created, and also
  rejects a map without colours or with fewer than one entry. *Tested:* `Colormap.Errors`.

## Rendering: axis ticks and numbers (src/render/core/ticks.cpp)

Tick locations, tick labels, offset / scientific-notation texts, autoscaled and colour-bar limits are
Matplotlib 3.11.2's bit for bit (ported with its notice); MDATK's `thin_ticks` likewise. Differences:

- **`thin_ticks` with `max_ticks` < 1.** MDATK: `max_ticks = 0` (the `max_ticks` option of the
  heatmap and contact-map panels) crashed the render with `ZeroDivisionError`; negative values
  labelled every row. MDBridge takes any value below 1 as 1 (one label, the first row).
  *Tested:* test_ticks.cpp `ThinTicks` (cases marked `improved` in
  `tests/data/render/ticks/tick_cases.json`, written by `make_tick_cases.py` from MDATK's
  `thin_ticks(n, 1)`).
- **Inputs Matplotlib cannot handle, unreachable from MDATK's figures:** `MaxNLocator` with fewer
  than one bin throws `std::invalid_argument` (Matplotlib divides by zero and returns garbage); a
  formatter with a fixed offset on an empty view interval uses exponent 0 (Matplotlib: `math domain
  error`); the offset search stops at 10^-400 instead of looping forever on intervals one ulp wide.
  An empty `steps` list means Matplotlib's default steps (Matplotlib: `IndexError`). *Tested:*
  test_ticks.cpp `MaxNLocatorRejectsNoBins`, `ScalarFormatter` (the recorded error cases).

## Rendering: fonts and text (src/render/core/font.cpp, text.cpp; src/render/text_qt.cpp)

MDATK draws text with Matplotlib's Agg backend in the DejaVu fonts Matplotlib ships. MDBridge
bundles the same five files (DejaVu 2.35 from Matplotlib 3.11.2's wheel, `data/fonts`: Sans, Sans
Bold, Sans Oblique, Sans Mono, Sans Mono Bold) compiled into the programs, chooses among them with
Matplotlib's font-manager scoring, and lays text out with its own TrueType reader. Differences:

- **Unhinted text metrics and glyphs.** With its default `text.hinting`, Matplotlib 3.11 measures and
  draws glyphs hinted by FreeType's TrueType interpreter: advances and heights snap to whole pixels,
  so a label's width depends on the dpi grid and on the FreeType build (up to 0.66 px per glyph and
  2 px in height against the unhinted outline at MDATK's sizes; 17 % narrower for small monospace
  text). MDBridge measures and draws the unhinted outlines, exactly as Matplotlib does with
  `text.hinting: no_hinting` (width, height, descent and every glyph position equal in 26.6 fixed
  point, kerning and ligatures included). *Why:* the same layout and text pixels on Linux, Windows
  and any server, at every dpi, independent of a font engine's version or the system's font
  settings; the figure layout tolerance (docs/specs/m4-plan.md, section 3) covers the difference.
  *Tested:* `tests/render/test_text.cpp` (`TextMetrics.AsMatplotlibUnhinted`: exact against
  `no_hinting`, and within 0.7 px per glyph / 2 px / 1.5 px of the default hinting, for every case of
  `tests/data/render/text/text_cases.json`).
- **Figures never use the machine's fonts.** MDATK's text follows the Matplotlib configuration of the
  machine it runs on (`font.family` of a matplotlibrc, other DejaVu releases or fonts installed
  there), so the same project renders differently on two machines. MDBridge always uses the bundled
  faces; font properties choose among them as Matplotlib's `findfont` scores (bold oblique and mono
  oblique, which MDATK does not use, are not bundled: the closest bundled face by that score).
  Missing glyphs are the font's `.notdef` box (Matplotlib falls back to other fonts with a warning).
  *Tested:* `FontChoice.AsMatplotlibsFontManager` (Matplotlib's choice among the bundled faces for
  792 property sets), `BundledFonts.AreMatplotlibsFiles` (the files' checksums).
- **Text pixels come from the outlines.** Raster figures fill the glyph outlines as paths with Qt's
  antialiasing instead of FreeType's glyph bitmaps (Agg); shapes and positions are the unhinted
  outline's, the antialiasing differs slightly (no gamma). A text image is identical wherever the
  same build runs: CI draws one on the build machine and on a bare Rocky Linux 9 without fonts,
  fontconfig or X11 and compares the pixels (`packaging/smoke_test.py --compare-text-report`).
- **No display needed.** `render`, `charts`, `preview` (and `fonts`) start Qt on its `offscreen`
  platform on Linux whatever `DISPLAY` or `QT_QPA_PLATFORM` say, like MDATK's Agg backend, which
  never needed a display. *Tested:* `TextQt.BatchPlatform`, the smoke test's `fonts` check.

## Rendering: the figure framework (src/render/core figure*.cpp, display_list.cpp; src/render/figure_qt.cpp)

MDATK composes its figures with Matplotlib (GridSpec, the "constrained" layout engine, Agg, SVG and
PDF backends). MDBridge has its own figure model with Matplotlib's layout rules ported (constrained
layout of one GridSpec with height ratios and MDATK's pads, tick / label / title / colour-bar space
from the text metrics, legends with "best" placement, colour bars right of their axes, twin axes),
a display list, and QPainter backends (QImage, QSvgGenerator, QPdfWriter). Differences:

- **Exact canvas size (RN-06, the figure part).** Matplotlib sizes a figure in inches
  (`w / dpi`) and truncates the canvas back to pixels, so some sizes and dpis give a canvas a pixel
  smaller than requested (MDATK then stops with `canvas size ... not even` / `choose a DPI that
  divides the resolution`). MDBridge's figure is exactly `w x h` pixels at any dpi; inches and
  points are derived from it. *Tested:* `FigureFramework.FrozenLayoutAndExactPixels`,
  `FigureRaster.ExactSizeAndColours` (1000 x 562 at dpi 77, 320 x 180 at dpi 120).
- **Layout solved in closed form.** Matplotlib's constrained layout feeds its margins to the
  kiwisolver constraint solver; MDBridge computes the same solution directly (columns and rows fill
  the figure, inner sizes in the grid's ratios), with Matplotlib's two passes, margins, colour-bar
  pads, submerged-margin matching and collapse warning. Positions differ from Matplotlib's by the
  text metrics only (the unhinted text above): every axes, twin and colour-bar box, tight box, axis
  label, title, text and legend of the 97 drawn figure cases lies within 1.25 % of the figure size of
  Matplotlib's (tolerance 1.5 %, docs/specs/m4-plan.md section 3), and every tick and tick label is
  equal (the automatic bin counts agree on every axis). *Tested:* `tests/render/test_figure.cpp`
  (`FigureCase.LayoutTicksAndTextsAsMatplotlib` on `tests/data/render/figures`).
- **"Best" legend locations follow MDBridge's text widths.** Matplotlib places a `loc="best"`
  legend where it covers the fewest line vertices, line segments, points and texts; its box width
  comes from the hinted text, up to 7 % narrower than the unhinted text at small pixel sizes (8.8 pt
  at 60 dpi). With Matplotlib's box size MDBridge's search picks Matplotlib's location in every case;
  with its own (wider) box, 3 of the 23 "best" legends of the figure cases go to another corner
  where the wider box fits better. *Tested:* `FigureCase.LayoutTicksAndTextsAsMatplotlib` (the
  location for Matplotlib's box in Matplotlib's axes box, and the box size).
- **Vector output keeps text as text in the bundled fonts.** SVG files have `<text>` elements in the
  DejaVu families (Matplotlib's default SVG output converts text to paths); PDF pages embed the
  bundled faces. Raster images fill the glyph outlines (above). The SVG size is given in mm (the
  same physical size as Matplotlib's pt). *Tested:* `FigureSvg.TextAsTextAndPhysicalSize`,
  `FigurePdf.PageOfTheFigureSize`.

## Rendering: plot panels (src/render/panels plots.cpp)

The time-series, dual-axis, heatmap and contact-map panels are MDATK's (`mdatk/render/panels.py`)
on MDBridge's figure model; every figure case of these panels has the same contents as MDATK's
(lines, dots, cursors, targets, legends, images and the RGBA they show, colour bars and their
tables, ticks, labels, titles, animated artists at the first, middle and last frame; layout within
the tolerance). Differences:

- **Heatmap cells at the midpoints between frame times (RN-01).** MDATK made every column of a
  residue heatmap `t[1] - t[0]` wide, so with non-uniform times (concatenated runs with different
  time steps) the cells sat at the wrong times and the moving cursor left the cell of its frame.
  MDBridge's column edges are the midpoints between frame times, half a step outward at both
  ends; uniform times give MDATK's cells, one frame `t0 - 0.5 .. t0 + 0.5` as MDATK. *Tested:*
  `tests/render/test_plot_panels.cpp` (`PlotArithmetic.ColumnEdgesAtMidpoints`; the case
  `heatmap_uneven_times`, 12 frames 10 ps apart then 12 frames 2 ps apart, whose expected data
  `improvements_render.py` makes MDATK's panel record with these edges; `PlotFrameLoop.
  UnevenColumnsDrawnAtTheirWidths`: the drawn image shows each frame's column at its time).
- **Heatmap `vmin` above `vmax` (RN-04).** The plan expected MDATK to stop with Matplotlib's
  `minvalue must be less than or equal to maxvalue`. It does not: the heatmap's colour bar runs
  `nonsingular()` on the reversed norm through Matplotlib's norm callbacks, which set the new vmin
  while the old vmax is still in place, so `vmin = 5, vmax = 1` silently became `0.9 .. 5` (an
  artefact of the callback order; Matplotlib without a colour bar does raise). A negative-valued
  matrix with the default `vmin = 0` was never an error (vmax is at least 1e-9: everything shows the
  lowest colour). MDBridge swaps reversed limits (the given ones, or a `vmin` above the largest
  value) and warns `heatmap <key>: vmin 5 is above vmax 1; using vmin 1, vmax 5`. *Tested:*
  `PlotArithmetic.HeatmapNormSwapsReversedLimits`, the case `heatmap_vmin_above_vmax` (norm,
  colour bar, pixels and the warning). Validating the options in `check()` is part of RN-15.
- **NaN-aware running mean (RN-14).** The time-series panel's `smooth` (and the series charts'
  smoothing, which use the same function) averaged with a convolution, so one NaN (an empty RMSD
  target, a NaN frame) blanked every window it fell in. MDBridge averages each window's finite
  values (MDATK's edge padding and window); a value is NaN only where its whole window is.
  *Tested:* `PlotArithmetic.RunningMeanNanAware`, the cases `ts_nan_smooth` and `ts_smooth`
  (and the chart records, regenerated with the same patch).

## Rendering: charts (src/render/panels/charts.cpp; src/render/core/figure_formats.cpp)

The standalone charts are MDATK's (`mdatk/render/charts.py`): the same charts (every series, or
`output.charts.series`; the matrices with `include_heatmaps`), sizes, chrome, smoothing, targets,
labels, title templates, file names, log (`Wrote N chart file(s) to D`), progress (`Chart <key>`)
and warnings. Every chart case has MDATK's paths, files, messages and figure contents
(`tests/render/test_charts.cpp` against `tests/data/render/charts`; the files' first bytes and
sizes from the Qt backend in `tests/render/test_charts_qt.cpp`). Differences:

- **Matrix-chart cells at the midpoints between frame times (RN-01), NaN-aware smoothing (RN-14).**
  The heatmap panel's improvements (above) apply to the charts, which share the code. *Tested:*
  the cases `uneven_times` (two runs with different time steps) and `nan_smoothing`.
- **Matrix-chart height capped (RN-05).** A matrix chart is 0.18 in per residue row tall, so a
  system of thousands of residues gave images beyond Agg's 2^16 px limit (a crash after the
  analyses ran) and, below it, images of hundreds of megabytes. MDBridge caps the height at 2^15 px
  (the row labels are thinned to 40 anyway) and warns `chart K: N rows need H px at D dpi; height
  capped at 32768 px`. *Tested:* `ChartHeight.CappedAtTwoToTheFifteen`, the case `matrix_tall`
  (700 rows at 300 dpi).
- **Matrix keys in `output.charts.series` (RN-10).** MDATK dropped them silently (only series were
  charted from that list). MDBridge charts a requested matrix as a matrix chart: after the series,
  each matrix once, before the other `include_heatmaps` matrices. A pair-map key (no chart type)
  warns `chart series K is a pair map; charts draw series and matrices (a contact_map panel shows
  it)`. *Tested:* `ChartRun.MatrixSeriesAndPairMaps`, the cases `missing_and_matrix`,
  `matrix_series`, `matrix_series_heatmaps`.
- **Formats checked first (RN-11).** MDATK wrote every format Matplotlib and pillow know and failed
  with Matplotlib's `ValueError` traceback at the first unknown one, after the analyses ran and
  after the files before it were written. MDBridge writes png, svg, pdf, jpg / jpeg (Qt's JPEG
  plugin), tif / tiff, eps and ps (its own writers: an uncompressed RGBA TIFF; PostScript Level 2
  with text as glyph outlines and, as Matplotlib's PS backend, no transparency) and checks the list
  before anything is drawn: `chart format 'bogus' is not supported (choose from png, svg, pdf, jpg,
  jpeg, tif, tiff, eps, ps)` (`check()` reports the same, W3-B). Matplotlib's gif, webp, avif,
  svgz, pgf, raw and rgba are no longer offered. The plan assumed jpg and eps failed in MDATK; they
  did not, and they still work. *Tested:* `ChartFormats.CheckedBeforeAnythingIsWritten`, the cases
  `formats`, `format_jpg_eps`, `format_tif_ps`, `format_unknown`, `format_unknown_two`;
  `ChartFiles.*` and `ChartPostScript.VectorsTextAndImage` (Qt).
- **Title and colour bar in the theme's colours (RN-50).** `style_axes` styled the axes' centre
  title, but charts set the left one, which stayed Matplotlib's black at 12 pt; the matrix charts'
  colour bar (outline, ticks, labels) stayed black too: black on the dark theme's near-black
  surface, unreadable. MDBridge draws the title in the theme's ink at 11 pt (as `style_axes`
  meant) and the colour bar's outline in the grid colour, its ticks muted and its tick labels and
  label in ink2 (sizes unchanged). Matplotlib's `default` style is unchanged. *Tested:* every chart
  case with MDATK's chrome; `dark_heatmaps`.
- **Styles (owner decision D5).** MDATK applied any Matplotlib style sheet (`ggplot`,
  `seaborn-v0_8`, a style file or URL) and warned only for unknown names, with Matplotlib's
  message. MDBridge draws MDATK's chrome (`mdatk`, or empty) and Matplotlib's `default` look; any
  other name warns, for each chart as MDATK warned, `unknown matplotlib style 'ggplot' (MDBridge
  draws the styles 'mdatk' and 'default'); using MDATK chrome` and draws MDATK's chrome. *Tested:*
  `style_default`, `style_default_smoothing`, `style_default_color`, `style_ggplot`,
  `style_unknown`, `heatmaps_dark_style`.

## Rendering: molecule geometry (src/render/molecule)

The 3D representations (MDATK's `mdatk/render/molecule/topology.py`, `geometry.py`, `surface.py`,
`scene.py`) give MDATK's topology roles, DrawLists, surfaces and colours (vertices within 1e-5 Å,
indices exact; `tests/render/test_molecule_geometry.cpp` against
`tests/parity/make_molecule_geometry_cases.py`), with these differences:

- **Oversized surface grids (RN-19).** A `surf` grid of more than 2^26 points (a blown-up frame with
  atoms thousands of Å apart, or an implausibly large system: 2^26 points is a 400 Å cube at 1 Å,
  600 Å at the 1.5 Å used above 5000 atoms) is refused before anything is allocated, with
  `surface grid of NX x NY x NZ points is too large (at most 67108864 points)`. MDATK allocated it
  (numpy needs over 30 bytes per point): a `MemoryError` that left the surface out of that frame,
  or, below what the allocator refuses, gigabytes of memory and minutes per frame. MDBridge's scene
  leaves the surface out of that frame with MDATK's warning, once per representation
  (`3D representation K (surf) not drawn on some frames: surface grid of ...`), and draws the other
  representations. A representation whose geometry fails for any other reason is skipped the same
  way (MDATK caught only `ValueError`, `MemoryError`, `IndexError` and `FloatingPointError`; others
  aborted the video). *Tested:* `MoleculeGeometry.OversizedSurfaceGridIsRefused`.
- **Nearest atom of a surface vertex.** scipy's `cKDTree` is replaced by a cell grid; where two atoms
  are exactly as near, MDBridge takes the lower atom index (cKDTree's choice depends on its tree).
  *Tested:* `MoleculeGeometryParity.SurfacesMatchMdatk` (distances where the indices differ),
  `MoleculeGeometry.NearestPointsAgreeWithBruteForce`.

## Rendering: software 3D renderer (src/render/molecule/software.cpp)

The built-in renderer (MDATK's `mdatk/render/molecule/software.py`: z-buffer with sphere and stick
impostors and mesh rasterisation, tile culling of hidden spheres, the front-most translucent layer,
fog, glow tint and halo, ink outlines, alpha anti-aliasing, `fast` quality) gives MDATK's RGBA bit
for bit on every recorded render (37 images of scenes and of the synthetic items of MDATK's
`tests/test_molecule_render.py`; `tests/render/test_molecule_render.cpp` against
`tests/parity/make_molecule_cases.py`, plan tolerance >= 99.5 % of pixels equal, all within 2/255).
It renders on several threads (bands of rows; the image does not depend on their number) where
MDATK used one. Differences (the generator applies RN-26 and RN-27 to MDATK test-side,
`improvements_render.SOFTWARE_RENDERER`):

- **Halo seam at the molecule's box (RN-26).** The edge anti-aliasing compares each alpha with its
  3 x 3 mean over the molecule's bounding box (plus 2 px). MDATK filtered that box alone, with zeros
  outside it, so inside a glow halo (which reaches beyond the box) every pixel along the box's edge
  was averaged with zeros: a rectangle about 17 % darker drawn into the halo around every glowing
  molecule (e.g. alpha 61, 48, 60 across the edge where the halo rises smoothly through 61, 66, 72).
  MDBridge's mean sees the real neighbours around the box (zeros only beyond the image); where the
  image around the box is clear the result is MDATK's to the bit. *Tested:*
  `MoleculeRender.HaloHasNoSeamAtTheMoleculeBox`, `MoleculeRender.MatchesMdatkRenders` (the glow
  renders).
- **Fast quality at odd sizes (RN-27).** `fast` renders at half size: MDATK used
  `int(round(size * 0.5))` (round half to even), so a 101 x 77 view became 50 x 38, doubled to
  100 x 76: the last column and row stayed transparent and the image was scaled by 100 / 101 (a
  103-px view became 52 px, doubled to 104 and cropped: scaled by 104 / 103). MDBridge renders
  `ceil(size / 2)` pixels with the view's limits widened to exactly two pixels per half pixel, then
  crops: the whole image is covered at the view's scale. Even sizes are unchanged. *Tested:*
  `MoleculeRender.FastQualityAtOddSizesCoversTheWholeImage`, `MatchesMdatkRenders` (case `fast`).
- **Unusable views and primitives (RN-28).** Equal x or y limits made MDATK fail with
  `ZeroDivisionError: float division by zero` (the panel then logged `3D representation not drawn
  (2D trace shown instead): float division by zero`); reversed limits failed for spheres with
  numpy's `negative dimensions are not allowed` (and drew meshes mirrored); NaN limits, rotation or
  depth range gave garbage. MDBridge throws
  `std::invalid_argument("3D view limits are empty, reversed or not finite")` (or `... rotation or
  depth range is not finite`, `3D view size must be at least 1 x 1 pixels`). Primitives whose
  screen position is not finite in float32 or lies beyond 2^40 px are skipped like MDATK's
  non-finite atoms (MDATK cast such values to integers: undefined results); non-finite glow
  values count as no glow (MDATK: NaN colours); a sphere far larger than the image is rasterised
  only where it is on the image (MDATK built its whole sprite of (2r + 3)^2 entries: gigabytes, or a
  `MemoryError`, for a sphere 10^4 pixels across). Inconsistent items (colour or owner counts, triangle indices out of range) throw
  `std::invalid_argument` instead of numpy's `IndexError`. *Tested:*
  `MoleculeRender.UnusableViewsAndPrimitives`.
- **Halo strength independent of earlier renders (RN-29).** MDATK cached the halo's normalisation
  (`_disc_norm`) under the halo size rounded to 0.01 px, computed from the first size that rounded
  there, so a render's halo depended on which renders the process had made before (by up to the
  normalisation's change over 0.005 px). MDBridge computes it for the exact size; no recorded render
  collides. *Tested:* `MoleculeRender.HaloDoesNotDependOnEarlierRenders`.
- **Fragment order.** MDATK numbered every fragment (`2^32` fragments per buffer at most: `scene too
  large for the built-in renderer`); MDBridge keys fragments by their primitive's rank in MDATK's
  fragment order, which picks the same front-most fragment and breaks depth ties the same way, with
  `2^32` primitives per buffer at most.

## Rendering: density shadow (src/render/molecule/density_shadow.cpp)

The projected density map of the structure view (MDATK's `project_density`, `_above_level`,
`shadow_rgba` and `anchored` in `mdatk/density.py`, with scipy's order-1 `affine_transform`) gives
MDATK's images, extents and placements (images and RGBA within 1e-6, bit-identical on the parity
cases; placements within 1e-9; `tests/render/test_density_shadow.cpp` against
`tests/parity/make_density_shadow_cases.py`, which applies the items below test-side and records
MDATK's own outcome next to them), with these differences:

- **Non-finite voxels (RN-30).** MDATK's `np.maximum(map - level, 0)` keeps NaN, so one NaN voxel
  made the integrated image NaN where it was sampled and the peak NaN, i.e. the whole shadow NaN
  (an infinite voxel gave inf / inf = NaN or a blank image). MDBridge counts non-finite voxels (and
  any non-finite difference to the level) as no density, and a non-finite level (the map level of a
  map holding NaN is NaN, as in MDATK) gives no shadow. *Tested:*
  `DensityShadow.AboveLevelMatchesMdatk`, `DensityShadow.ProjectDensityMatchesMdatk` (cases
  `nan_voxel`, `inf_voxel`, `*_map_level`).
- **Non-finite inputs (RN-31).** A NaN view, centre or voxel size (a map with NaN cell lengths reads
  fine) or an infinite placement made MDATK fail in `int(math.ceil(...))` with `cannot convert float
  NaN to integer` (or an `OverflowError`), which the structure panel logged as `density map not
  drawn: cannot convert float NaN to integer`. MDBridge raises `ValueError("density map placement,
  view or grid is not finite")` at the same point. *Tested:* `DensityShadow.ProjectDensityMatchesMdatk`
  (cases `nan_centre`, `nan_view`, `inf_translation`, `nan_voxel_size`).
- **Unusable anchors (RN-32).** `anchored` with a ragged stored anchor (a hand-edited project file)
  raised numpy's `inhomogeneous shape` `ValueError`, and with NaN in the anchor or in the reference
  frame `LinAlgError: SVD did not converge`; the map was not drawn. MDBridge uses the stored
  placement as it is (as MDATK does for an anchor of another shape) and returns a note
  (`density anchor ignored: ...`) for the caller's warning. *Tested:* `DensityShadow.AnchoredMatchesMdatk`
  (cases `ragged`, `nan_anchor`, `nan_reference`).

## Rendering: structure view (src/render/panels/structure.cpp)

The structure panel (MDATK's `StructurePanel` in `mdatk/render/panels.py`, with `stage.py` and
`renderers.py`) draws what MDATK draws: the projected trace (points, links, the two glow halos,
ions), residue tags, the time label, the legend, the bare equal-aspect window, the density shadow
behind or in front (following a followed centre) and the 3D layer of the representations (the
molecule library's scene, the built-in software renderer), with the 2D trace as the fallback for a
frame whose 3D render fails. `tests/render/test_structure_panel.cpp` compares every figure case
whose first panel is a structure view with MDATK's record (`tests/parity/make_figure_cases.py`;
positions, sizes, colours, widths, texts, legend, limits, the shadow's extent and RGBA, the 3D
layer's size and coverage, the panel's notes and animated artists, its warnings) at the first,
middle and last frame. The generator applies the items below to MDATK test-side
(`tests/parity/improvements_render.py`); differences:

- **Non-finite coordinates (RN-03, the structure part).** One NaN trace point (a blown-up frame)
  made MDATK's bounding-sphere radius NaN and the figure failed with Matplotlib's `Axis limits
  cannot be NaN or Inf`. MDBridge frames the view by the finite points (the model centre too),
  warns once (`structure view: non-finite coordinates in frame(s) 3; left out of the view (those
  residues are not drawn there)`, at most ten frame numbers then `and N more`) and leaves those
  points, their links, halos and tags out of the frames concerned. *Tested:* case `structure_nan`
  (was MDATK's error), `StructurePanel.NonFiniteFrameIsLeftOut`.
- **Ion legend (RN-08).** MDATK's legend always said `Mg²⁺`, whatever `analyses.mg.mg_selection`
  and `show_ions` selected. MDBridge names the element of the trace's ions (`meta.ion_labels`:
  residue name, charge signs and digits dropped): `Mg²⁺`, `Na⁺` (NA, SOD), `K⁺` (K, POT),
  `Ca²⁺` (CA, CAL), `Zn²⁺`, `Mn²⁺`, `Cl⁻` (CL, CLA), `Li⁺`, `Rb⁺`, `Cs⁺`, and `Ions` for a
  mixture, unknown names or no labels. *Tested:* `StructurePanel.IonAndChainLegend`; the demo's
  `Mg²⁺` in every case.
- **Chain legend (RN-09, owner decision D11).** MDATK's chain-colour legend numbered the chains
  (`Chain 1`, `Chain 2`) while the sequence panel and residue ranges use the chain names. MDBridge
  writes `Chain A` (the trace's chain names, letters from the chain index when the trace has none,
  as the sequence panel does; the number only for a chain without a name). *Tested:* cases
  `structure_chain`, `structure_chain_palette`, `protein_trace`, `chains_sequence`,
  `StructurePanel.IonAndChainLegend`.
- **`renderer = "vtk"` (RN-16, the panel part).** MDATK probed for VTK and, without it, noted `VTK
  is not installed (pip install vtk); using the built-in renderer`. In MDBridge `vtk` (like
  `vulkan`) means the GPU renderer: the frame loop and the GUI install the Vulkan renderer
  (`setStageRendererFactory`); where there is none the panel uses the built-in renderer with the
  note `the GPU renderer (Vulkan) is not available; using the built-in renderer` (logged as a
  `3D representations:` warning). *Tested:* case `rep_trace_vtk`, `StructurePanel.GpuRendererNote`.
- **Density level of maps with non-finite voxels (RN-33).** The shadow's contour level is `mean +
  level x rms` of the map; with one NaN voxel MDATK's level was NaN, and with RN-30 such a map
  showed no shadow at all. MDBridge computes the level from the finite voxels (the same value as
  MDATK's when all voxels are finite), so the map's shadow shows. *Tested:*
  `StructurePanel.DensityLevelFromFiniteVoxels`.
- **No analysed frames (RN-02, the panel part).** A trace without frames made MDATK fail with an
  `IndexError`; the panel shows `No data: no frames`.
- **The centre warning names its setting (RN-51; follows MDATK 0.1.1 `render/panels.py`
  `_view_centre`).** A `visual.center_on` that matches no residue warns `Structure view ›
  Orientation › Centre on: no residue of the structure view matches 'Z:99'; using the model centre`
  (MDATK 0.0.4.1: `structure view centre: ...`), the readable location MDATK 0.1.x gives that
  setting. *Tested:* case `structure_centre_unknown`, `StructurePanel` zoom-and-centre test.

## Rendering: residue shading, orientation and the example dataset (src/render/core)

The structure view's colours and glow (MDATK's `mdatk/render/coloring.py`, `residues.py`), its
orientation (`orientation.py`) and the GUI's example dataset (`example.py`, with the demo's
nucleotide geometry) give MDATK's numbers bit for bit (`tests/render/test_shading.cpp` against
`tests/parity/make_shading_cases.py`; rotations through numpy's BLAS within 1e-15, sin / cos / exp
bit for bit with glibc). The example series are seeded as MDATK seeds them: zlib's crc32 of the data
key into numpy's `default_rng` (SeedSequence, PCG64, `uniform` / `normal` / `integers`, ported from
numpy 2.5.5 with its notice), so the preview shows the same placeholder curves. The expected data
is MDATK's with `tests/parity/improvements_shading.py` for the entries below.

- **Example data of a single-strand structure (RN-21).** With one structure group (or one chain in
  the loaded reference frame), MDATK's example base pairs sliced `a[:m - 2]` with `m = 0` and the
  example dataset failed with `IndexError: shape mismatch`, so the GUI's layout preview broke.
  MDBridge gives such a structure no example pairs (a frayed end never goes below zero pairs).
  *Tested:* `ExampleDataset.SingleStrandAndLargeGroupsRN21RN25`, cases `single_strand`,
  `single_strand_atoms`.
- **Angles next to -180 degrees (RN-22).** `wrap_degrees` used `np.isclose(a, -180)`, whose relative
  tolerance turned every angle between -180 and -179.9982 into +180: a drag or a matrix read back
  near -180 jumped by up to 0.0018 degrees. MDBridge folds only -180 itself (within 1e-9) to +180.
  *Tested:* `Orientation.AnglesNextToMinus180StayRN22`, the `wrap` cases.
- **Unreadable `visual.label_residues` (RN-23).** MDATK ignored a spec it could not parse without a
  word (no tags from it, `label_every` tags kept). MDBridge keeps that behaviour and warns once:
  `visual.label_residues: invalid residue range item 'x'; no residue-number tags from it`.
  *Tested:* `Residues.RowsTagsAndShading` (tags with their log lines).
- **Chain-qualified `visual.label_residues` (RN-24).** MDATK dropped the chain of `B:12` and tagged
  residue 12 of every chain; MDBridge tags it on chain B only (chain names as `residue_rows` reads
  them, which already honoured the chain). *Tested:* `Residues.ChainQualifiedTagsRN24`.
- **Large residue groups in the example data (RN-25).** `group_resids` listed every number of a
  group's ranges before keeping the first 400 (a group `1-1000000000` took gigabytes and minutes),
  and the schematic trace's sequence had 120 letters, so a shown group of more than 120 residues
  raised `IndexError`. MDBridge stops at 400 residues and repeats the sequence. *Tested:* case
  `big_groups`, `ExampleDataset.GroupResids`.
- **Fewer than one example frame.** `exampleDataset` needs at least one frame
  (`the example dataset needs at least one frame`; MDATK: `IndexError` at `t[-1]`). One frame gives
  MDATK's NaN placeholders (`t / t[-1]`). *Tested:* case `one_frame`.

## Rendering: sequence and secondary-structure panels (src/render/panels/seqpanels.cpp)

The `sequence` and `secondary` panels (MDATK's `mdatk/render/seqpanels.py`) show what MDATK's show:
with MDATK's frozen axes boxes their records (`dumpFigure`) equal MDATK's at every recorded position
of every figure case that has one (cells, letters, boxes, tags, chain labels, the re-flow after the
layout with its cell size and rows, the note and its warning, arcs within and between chains,
halos, colours, fonts and limits; numbers within 1e-9, text extents within the layout tolerance),
and MDBridge's own layout places their axes within 1.5 % of the figure size
(`tests/render/test_seqpanels.cpp` against `tests/parity/make_figure_cases.py`, with the items
below applied by `tests/parity/improvements_render.py`). A dataset without frames gives
`No data: no frames` in both panels, as RN-02 asks for every panel (MDATK built them and failed
later). The differences:

- **Unusable `per_row` (RN-40).** The sequence panel read `int(options.per_row)`: a value `int()`
  cannot read (`"x"`, `"7.5"`, a list) raised `ValueError` / `TypeError` out of `build_figure`, which
  stopped the render or the preview with a traceback; a negative value inverted the x axis and gave
  boxes of tens of thousands of points squared over the whole figure. MDBridge shows
  `No data: option 'per_row' must be a whole number of residues per row (0: automatic), not -3` and
  carries on; 0, empty and values `int()` truncates (`7.9`, `"7"`, `0.5` = automatic) are read as
  MDATK reads them. *Tested:* `SeqPanels.UnusablePerRowIsAMessageRN40`, case
  `sequence_per_row_invalid`.
- **Unknown `source` (RN-41).** The secondary panel treated every `source` other than `first` and
  `dotbracket` as `dynamic`, so a typo (`dynamc`, `dot-bracket`) silently showed the base-pair
  analysis instead of the intended pairs. MDBridge shows `No data: unknown source 'dynamc' (option
  'source': dynamic, first or dotbracket)`; the analyses it requires are unchanged (`basepairs` for
  any source but `dotbracket`). *Tested:* `SeqPanels.UnknownSourceIsAMessageRN41`, case
  `secondary_source_unknown`.
- **"No data" frame (RN-42).** A "No data" panel shows its message on the axes' frame (left and
  bottom spines in the grid colour), but the secondary panel hid its spines before reading the
  dot-bracket or the base-pair map, so its messages (`no dot-bracket given`, a dot-bracket of the
  wrong length, no `basepairs/map`) floated without a frame. MDBridge's sequence and secondary
  panels keep the frame every other "No data" panel has. *Tested:*
  `SeqPanels.NoDataPanelsKeepTheirFrameRN42`, cases `secondary_dotbracket_empty`,
  `secondary_dotbracket_length`. The structure view likewise keeps the frame when its proximity
  (`No data: proximity`) or atoms (`No data: atoms`) result is missing, which it notices after
  hiding its spines (W2-B; *tested:* `StructurePanel.MissingResultsShowNoData`).

## Rendering: the Vulkan 3D renderer (src/render/vulkan/scene_renderer.cpp, src/render/panels/vulkan_stage.cpp)

MDATK's GPU path was VTK (`renderer: vtk`): OpenGL impostor spheres and sticks with VTK's own
lighting, no outlines, the glow tint and halo computed in numpy, and a crash risk without EGL
(RN-16). MDBridge's GPU renderer is Vulkan 1.3 and draws the software renderer's picture: the
primitives are set up on the CPU exactly as the software renderer sets them up (pixel boxes, the
spheres' quarter-pixel phase, the widened triangle edges, the capsule bound of the sticks), the
shaders repeat its per-pixel arithmetic in float32 (impostor spheres and cylinders, mesh
triangles, the depth key with ties to the later primitive, the same lighting), and compute passes
composite as the software renderer does (depth cue, glow tint, the halo from a quarter-resolution
grid with scipy's Gaussian and the disc normalisation, ink outlines, edge anti-aliasing, the
front-most translucent layer, fast quality). Differences:

- **Same look, not same pixels.** On lavapipe the images equal the software renderer's on every
  molecule case; on other drivers float32 rounding (division, square roots) may move single
  pixels. The tests require coverage IoU >= 0.99 and SSIM >= 0.97 against the software renderer on
  every render case and every scene of `tests/data/render/molecule` (also larger, rotated,
  glowing and fast), and the resolved G-buffer against `SoftwareRenderer::resolve` (>= 99.5 %
  of pixels, owners and depths). *Tested:* `VulkanScene.*`.
- **`renderer = "auto"` uses the GPU (RN-16).** MDATK's `auto` always meant its built-in
  renderer. In MDBridge `auto` uses the Vulkan renderer when it runs on a GPU, silently (same
  picture), and the software renderer otherwise; a Vulkan driver that runs on the CPU (Mesa's
  lavapipe) is left to the software renderer, which is about twice as fast there (`VulkanScene.
  Benchmark`; kept so by owner decision D15, docs/specs/m4-plan.md). `vtk` and `vulkan` use Vulkan on any device, the CPU's included, and fall back to
  the software renderer with the RN-16 note. The program installs this choice at start
  (`installGpuStageRenderers` in `main`), so figures, previews and video frames all use it; the
  device opens on first use only. *Tested:* `VulkanStage.*` (`AutoLeavesACpuDriverToTheSoftwareRenderer` on a CPU driver).
- **A GPU failure never stops a video.** When the GPU cannot draw a frame (device lost, out of
  memory, a scene or image beyond the device's limits) that frame and every later one are drawn by
  the software renderer, with one warning `3D representations: the GPU renderer failed (...);
  using the built-in renderer` (MDATK: VTK errors or crashes inside `Render()`). Unusable views and
  items raise the software renderer's errors (RN-28). *Tested:*
  `VulkanScene.UnusableViewsAndItemsThrowAsTheSoftwareRenderer`, `VulkanStage.InstallsAndUninstallsTheFactory`.

## Rendering: video encoding (src/render/video encoder.cpp)

MDATK pipes RGBA frames into an `ffmpeg` executable; MDBridge links its own LGPL-only build of
FFmpeg's libraries with SVT-AV1 (`packaging/ffmpeg`; owner decisions D1, D2) and writes the same
kind of file: MP4, `yuv420p`, `moov` in front (faststart), no audio, frame k at k / fps with the
frame rate read as `ffmpeg -r` reads it, the same temporary name `.MDATK_tmp_<pid>_<name>` in the
target folder renamed at the end and removed on failure. Stream metadata (frame count, size,
`r_frame_rate`, time base 1/15360 at 30 fps, duration) equals what `ffmpeg -r 30 ... -movflags
+faststart` writes; verified against the ffmpeg command line while writing the encoder.

- **Codec mapping and quality (RN-07).** MDATK passes `codec` to ffmpeg; with `libx264` (the project
  files' default, kept by `init`: D12) a build without x264 fell back silently to the first of
  libx264 / libopenh264 / mpeg4, and `crf` was used only for libx264 / libx265 (`mpeg4` always got
  `-q:v 3`). x264 and x265 are GPL and cannot be bundled, so `libx264`, `libx265`, `h264` and `hevc`
  map to the AV1 encoder `libsvtav1` with one warning, `video: codec libx264 is not bundled (GPL);
  using libsvtav1 (AV1)` (`h264` / `hevc`: without "(GPL)"); `av1`, `libsvtav1` and `mpeg4` are
  used as named; any other name is used if FFmpeg has that encoder (as MDATK passed it through), else
  `video: encoder X is not available; using libsvtav1` (or `mpeg4`). `crf` is always on the x264 /
  x265 scale (0-51, clamped) and is mapped: SVT-AV1 `crf = round(crf * 63 / 51)` (20 -> 25,
  23 -> 28, 51 -> 63; preset 8, kept by owner decision D16 until a side-by-side quality
  comparison with preset 10), mpeg4 `qscale = round(1 + crf / 10)` up to 20 and
  `round(3 + (crf - 20) * 28 / 31)` above (20 -> 3, MDATK's fixed value; 23 -> 6; 51 -> 31). When
  SVT-AV1 cannot encode the size (below 4x4) the encoder falls back to mpeg4 with
  `video: libsvtav1 cannot encode WxH (<reason>); using mpeg4`. *Tested:*
  tests/render/test_video_encoder.cpp (`VideoCodecChoice.*`, `SizesSvtAv1CannotEncodeFallBackToMpeg4`).
- **Colours tagged BT.709.** ffmpeg converts RGBA to `yuv420p` with the BT.601 matrix and leaves
  the stream untagged; players read untagged HD video as BT.709, so MDATK's videos show slightly
  shifted hues (greens and reds). MDBridge converts with the BT.709 matrix (limited range) and tags
  the stream `bt709` / `bt709` / `bt709` / `tv`, so players show the colours as drawn. *Tested:*
  `Av1Default24Frames320x180` (colour tags); decoded frames were checked at 44 dB PSNR (AV1) and
  42 dB (mpeg4) against the RGBA input while writing the encoder.
- **LGPL check at run time.** The encoder refuses FFmpeg libraries whose `avcodec_configuration()`
  contains `--enable-gpl` or `--enable-nonfree` or whose licence is not LGPL (MDBridge may only be
  combined with LGPL builds), with a message naming the configuration. *Tested:*
  `VideoFfmpeg.LgplOnlyBuild`; the configure-time check is in cmake/FindFFmpegLGPL.cmake.
- **Messages.** MDATK reports `FFmpeg exited with code N: <stderr>` or `FFmpeg stopped
  unexpectedly`; MDBridge's `VideoError`s say which step failed (`video: cannot write PATH: ...`,
  `video: cannot open the encoder libsvtav1 for WxH at 30/1 fps: ...`, `video size WxH is not even
  (yuv420p needs even sizes of at least 2x2)`) with FFmpeg's own error text, without pointer
  addresses. Without FFmpeg in the build: `video encoding is not available in this build`. FFmpeg
  and SVT-AV1 print nothing on the console. *Tested:* `VideoFiles.ClearErrors`,
  `VideoFfmpeg.NotAvailableMessage`.

## Rendering: the video frame loop (src/render/video/frame_loop.cpp, src/render/video_qt.cpp)

MDATK's `render_video` builds the figure once, marks the panels' animated artists, draws the
static rest once (Matplotlib blitting) and per frame updates the panels, restores the background,
draws the animated artists and pipes the RGBA buffer into ffmpeg. MDBridge does the same steps on
its own figure (the static and animated display-list layers, QPainter composition) with the same
positions (`0, s, 2s, ...`; RN-17 kept), log lines (`Rendering N video frame(s) at WxH, F fps ->
PATH`, `Wrote PATH (N frames)`), progress (`Video frame k/n` every `n // 200` frames and after the
last), cancellation between frames and temporary file (removed on failure or cancellation, an
existing file at the target kept; `VideoRender.CancelLeavesNoFile`,
`VideoRender.ExistingFileKeptOnFailure`). Differences:

- **Frames composed on several threads.** The panels are updated in order on the calling thread
  (they and the 3D stage renderer they call keep per-frame state), each frame's animated layer is
  composed over the background on a pool of threads (at most 8; the 3D renderer has its own
  threads) and the frames are encoded in order on a thread of their own, while the next frames are
  drawn; the encoder's RGBA -> yuv420p conversion runs on several threads too (whole rows per
  thread). MDATK did all of it one frame at a time (with ffmpeg in a second process). Composition
  is a pure function of the background and the frame's layer, and the conversion computes every row
  from the whole frame, so the frames and the encoder's input are the same bytes for any number of
  threads. *Tested:* `FrameLoop.FramesInOrderAndIdenticalForAnyThreadCount`,
  `QtVideo.ComposedFramesIdenticalForAnyThreadCount`,
  `QtVideo.RepresentationsThroughTheStageRendererFactory`, `VideoFiles.ConversionThreadsGiveTheSameVideo`
  (the same mpeg4 file for 1, 3 and 8 conversion threads and for padded rows).
- **No analysed frames (RN-02, the video part).** MDATK started ffmpeg on zero frames (or failed
  in the panels); MDBridge stops with `VideoError: no analysed frames to render` before anything is
  written. *Tested:* `VideoRender.NoFramesIsAnError`.
- **No canvas size errors (RN-06, the video part).** The canvas is exactly `figure_size(w, h)`
  (even) pixels for any dpi, so MDATK's `canvas size WxH is not even` and `canvas is WxH px but ...
  choose a DPI that divides the resolution` cannot occur (`render_odd_size`: 333 x 187 at dpi 37
  gives 332 x 186, as MDATK's). *Tested:* `RenderCliCases/VideoRun.AsMdatksVideo`.
- **Quality check.** Decoded mpeg4 frames are within 42 dB (RGB PSNR) of the composed ones at the
  default 1920 x 1080; at the parity cases' 320 x 180 yuv420p's chroma subsampling alone limits RGB
  to about 35 dB (thin coloured lines), and the luma is within 42 dB. AV1 videos are checked by
  their metadata only (the bundled FFmpeg has no AV1 decoder). *Tested:*
  `QtVideo.DecodedFramesWithin38dB`; metadata (frame count, size, frame rate, pixel format,
  streams, faststart; not the codec name, RN-07) against `tests/data/render/cli` by
  `RenderCliCases/VideoRun.AsMdatksVideo` and `QtVideo.FullRun24Frames320x180` (MDATK's
  `test_full_run` video).

## Rendering: the render, charts and preview commands (src/app/commands/render_pipeline.cpp, src/render/panels/checks.cpp)

`render`, `charts` and `preview` take MDATK's options (`--video` / `--no-video`, `--charts` /
`--no-charts`, `--pdb` / `--no-pdb`, `--force`, `--from-cache`; `--series`, `--formats`, `--dpi`;
`--frame`, `--output`) and print MDATK's JSON / path. A full run keeps MDATK's order and progress
ranges, one time stamp for all its files, the log lines (`<program> <version> - output directory D`,
`Finished in N s`), the params log (`MDATK_params_log.toml`: inputs, cache keys, time stamp,
`outputs`, citations) and the step log entry `render` / `Analyses and outputs` in the project file
(files, `N frames analysed; M file(s) written to D`, status and summary of a failure). Every command
line of `tests/data/cli/cli_cases.json`'s render cases prints and writes what MDATK does, figures and
videos compared by their metadata (sizes, ffprobe's stream data and faststart; the contents of the
figures are compared by tests/render). Differences:

- **Program line and environment.** The first log line names MDBridge and its version (MDATK:
  `Molecular Dynamic Analysis Tool Kit (MDATK) 0.0.4.1`); the params log's `environment` table is
  MDBridge's (`mdbridge`, `platform`, `qt`, `ffmpeg`: the bundled libraries' versions) instead of
  Python's packages and the `ffmpeg -version` line; `run.project` is `MDBridge`. *Tested:*
  `RenderCommands.FullRunStampsEveryFileLogsTheStepAndWritesTheParamsLog`.
- **Settings checked before the analyses (RN-11, RN-12, RN-13, RN-15; RN-04's option check).**
  MDATK found unknown chart formats, invalid colours, unknown chain palettes and panel options of the
  wrong type in the middle of a video or chart run, after the analyses, as Matplotlib tracebacks
  (`ValueError: Invalid RGBA argument`, `TypeError`), sometimes after files were written. MDBridge's
  `check()` (every command that checks a project: `info`'s `problems`, `analyze`, `export-pdb`,
  `render`, `charts`, `preview`) reports them with the other problems, each naming its field:
  `output.charts.formats: chart format 'gif' is not supported (choose from png, svg, pdf, jpg, jpeg,
  tif, tiff, eps, ps)`; colours with Matplotlib's grammar and message (`groups[k].color`,
  `visual.base_color` / `proximity_color` / `glow_color` / `ion_color`, `density.color` with a map,
  `output.charts.color` when set): `visual.glow_color: Invalid RGBA argument: 'gold-ish'`;
  `visual.chain_palette: 'Tab10' is not a bundled colormap; choose from ... ; did you mean 'tab10'?`
  (owner decision D6); panel options the panels cannot read: `layout.panels[1].options.ylim:
  expected two numbers, got 5`, `smooth` / `frame` / `max_ticks` (`expected a whole number`;
  `max_ticks` at least 1, MDATK divided by it), `target` / `vmin` / `vmax` (`expected a number`),
  the contact map's `mode` (`expected instant, frame or frequency, got 'mean'`; MDATK silently drew
  `instant`). Values the panels accept as MDATK did (`"5"`, `true`, `2.7` truncated) pass. The runs
  from the cache, which skip `check()` as in MDATK, check these settings alone before reading the
  cache. *Tested:* cli_cases.json (`render_invalid`, `render_from_cache_invalid`,
  `charts_bad_format_check`, `charts_bad_format`), `RenderCommands.InvalidSettingsNameTheirFields`,
  `RenderProblems.*`.
- **Preview formats (RN-11).** The preview's file format (from `--output`'s extension) is checked
  before the analyses run: `preview format 'gif' is not supported (choose from ...)` (MDATK: a
  Matplotlib `ValueError` after the analyses). A name without an extension is written as
  `name.png`, as Matplotlib's `savefig` writes it, and MDBridge prints that name (MDATK printed the
  name without `.png`). JPEG needs Qt's `qjpeg` image plugin, shipped in the AppImage and the Windows
  installer (both packaging scripts check it). *Tested:* cli_cases.json (`preview_bad_format`,
  `preview_no_extension`, `preview_jpg`), `RenderCommands.FormatsAreCheckedBeforeTheAnalyses`.
- **No analysed frames (RN-02).** The video, the charts and the preview stop with `no analysed
  frames to render` before a figure is built (MDATK: an `IndexError` in the panels, the preview at
  position -1). *Tested:* `RenderPipeline.NoAnalysedFramesToRender`.
- **`output.video.ffmpeg_path` (owner decision D8).** Not used: the bundled FFmpeg libraries encode,
  and the setting stays in project files. A project naming another executable gets the log line `Note: output.video.ffmpeg_path (P) is not
  used; MDBridge encodes with its bundled FFmpeg libraries` and its video (MDATK: `error: ffmpeg not
  executable at P`). *Tested:* cli_cases.json (`render_ffmpeg_path`),
  `RenderCommands.FfmpegPathIsNotUsedAndTheDensityMapIsNotFitted`.
- **Density maps without a placement (owner decision D9).** MDBridge has no density fit yet: a
  shown map without a stored placement is not fitted before the video or the preview but reported,
  `WARNING: density map not drawn: no placement (run density fit)` (with `auto_fit` off MDATK's own
  warning). *Tested:* `RenderCommands.FfmpegPathIsNotUsedAndTheDensityMapIsNotFitted`.
- **Codec warnings (RN-07).** The video's codec choice warns in MDBridge's words (see "Video
  encoding"): `WARNING: video: codec libx264 is not bundled (GPL); using libsvtav1 (AV1)` for the
  default codec, `WARNING: video: encoder nosuchcodec is not available; using libsvtav1` (MDATK:
  `FFmpeg encoder nosuchcodec unavailable; using libx264`). *Tested:* cli_cases.json (every video
  case).
- **No VTK citation (RN-16).** MDATK cites VTK for `visual.renderer = "vtk"`; MDBridge has no VTK
  (`vtk` means its Vulkan renderer), so the params log, the step log and the BibTeX export never cite
  it. *Tested:* citation_cases.json (`pipeline`, `project`), `Citations.VtkIsNeverCited`.
- **`probe-video` (a check command of MDBridge's own).** `mdbridge probe-video FILE` prints what a
  video holds as JSON (streams, codec, size, frame rate, time base, frame and packet counts, pixel
  format, colour tags, `faststart`), read back with the bundled FFmpeg libraries, so a machine
  without `ffprobe` can check a `render` (the smoke test does so on a bare Rocky Linux 9). MDATK has
  no such command. *Tested:* packaging/smoke_test.py.
- **Locale.** The batch commands start Qt on the offscreen platform; on a host without a UTF-8
  locale Qt's notice about switching to `C.UTF-8` is not printed (the switch is made silently), and
  numbers in logs and files keep the decimal point whatever the user's locale.

## Interface: 3D view and figure preview (src/app/gui, src/app/commands/project_view.cpp, src/render/view.cpp, software_view.cpp)

Milestone 4's interface (owner decision D13) shows an opened project's structure in the 3D view and a
read-only preview of its figure; the interactive layout editor (MDATK's canvas: dragging panels,
rotating the model in the figure, moving the density map) is milestone 5. What it draws is MDATK's:
the reference frame is `Pipeline.structure_preview()` (the first analysed frame of the cleaned files
MDATK left, else the raw files; cpptraj never runs), the 3D scene is the structure panel's 3D layer
of that frame (the project's representations, colours, `visual.rotation`, depth cue toward the panel
surface, glow, ions), and the preview is the GUI canvas' figure (`example_dataset` with that frame,
`build_figure` at `output.video`'s size and dpi scaled to the space available, the watermark
`Example data – layout preview (structure: reference frame)`, panels at 35 % of the example time
axis). Deliberate differences from MDATK's GUI:

- **A 3D view of the reference frame.** MDATK's GUI showed the structure only inside the figure
  canvas. MDBridge's main window has a 3D view (Vulkan, or the software view when Vulkan 1.3 is not
  usable) above the preview; drag to rotate, wheel to zoom, double-click to reset. Its rotation and
  zoom are the view's own: they do not change `visual.rotation` (MDATK's canvas rotation edited the
  project; that comes with the editor in milestone 5). A project without representations is drawn
  as a cartoon in the 3D view (the figure keeps its 2D trace, which is no 3D picture). *Tested:*
  `ProjectView.*` (tests/cli/test_project_view.cpp), `ViewCamera.*`, `SoftwareView.*`,
  packaging/smoke_test.py `--gui` (screenshots: the view must show the structure, not the
  placeholder, on Vulkan / lavapipe and on the software view of a bare Rocky Linux).
- **The reference frame loads by itself.** MDATK drew the schematic duplex until "Load real frame"
  was pressed; MDBridge loads the reference frame on a worker thread as soon as a project is opened
  (progress in the status bar), so the preview's structure panel shows the real trace and the
  watermark says so. When the frame cannot be loaded (missing or unreadable input files, a bad
  selection) the preview keeps the schematic structure and the reason is shown in the status bar
  (and printed on stderr), never in a dialog; the window stays usable. With cleaning on and no
  cleaned files yet, the raw files are shown, as MDATK, with a note in the status bar.
  *Tested:* `ProjectView.CleaningWithoutCleanedFilesShowsTheRawFiles`, `ProjectView.MissingInputsFail`,
  smoke_test.py (`gui runs 5 s and quits (a project without input files)`).
- **Software 3D view.** Without Vulkan the 3D view draws the molecule with the built-in software
  renderer (the structure panel's) on a worker thread, at fast quality while dragging; it shows the
  last image until the next one is ready, so the interface never waits for it. *Tested:*
  `SoftwareView.DrawsTheMoleculeAsTheSoftwareRenderer` (same pixels as the renderer, mouse rotation,
  zoom and reset).
- **Screenshots for checks.** `mdbridge gui PROJECT --screenshot PNG [--preview-screenshot PNG]`
  waits until the structure and the preview are drawn and reports what was shown (renderer,
  `molecule`, pixel `coverage`, the preview's size, load problems) as JSON (MDATK had no such check).

## Interface: the analysis window (milestone 5, AP items; tests/data/gui)

The expected behaviour of MDBridge's analysis window is MDATK v0.0.4.1's GTK window driven by action
scripts (`tests/parity/make_gui_cases.py`, tests/data/README.md "gui") with these fixes applied on
the test side (`tests/parity/improvements_gui.py`); each case keeps MDATK's own result under `mdatk`.
The window itself lands in milestone 5's waves 1 and 2 (docs/specs/m5-plan.md section 6); the C++ GUI
tests (`mdbridge_gui_tests`) replay the cases. *Tested:* the case named after each item, in
`tests/data/gui/<topic>/` (`roundtrip`, `pages`, `ap`, `canvas`).

- **AP-01: no clamping on load.** MDATK's spin boxes and sliders clamp a loaded value to their range,
  so opening and saving rewrote valid settings (`contacts.cutoff` 25 -> 20, `hbonds.angle_cutoff`
  60 -> 90, `energy.cutoff` 60 -> 50, `video.width` 7682 -> 7680, `visual.zoom` 0.01 -> 0.05, a
  rotation of -200 -> -180). MDBridge's fields hold any value of the setting's type; a field the user
  did not change saves its loaded value verbatim (the whole window round-trips losslessly), a value
  typed beyond the range is kept and left to validation, and the wheel zoom multiplies the stored
  zoom (within MDATK's wheel limits 0.05 ... 50). Values the user or the canvas sets keep MDATK's
  storage rounding (rotation, representation size / opacity, map opacity: 3 decimals). *Why:* opening
  a file must not change it. Not fixed in 0.1.1. *Tested:* `roundtrip/lossy_values`,
  `ap/ap01_set_out_of_range`, `canvas/wheel_zoom_from_tiny`, every `roundtrip/open_save_*`.
- **AP-02: optional numbers keep every digit.** MDATK showed `inputs.dt_ps` and
  `energy.target_value` with `%g` and saved what it showed (0.0123456789 -> 0.0123457, -1234567.5 ->
  -1234570). MDBridge shows and stores the shortest round-trip text. Not fixed in 0.1.1. *Tested:*
  `roundtrip/lossy_values`, `ap/ap02_typed_numbers`.
- **AP-03: colour fields keep their text.** MDATK's colour buttons parse with GDK and write
  `#rrggbb`: Matplotlib colours MDATK itself accepts became black (`tab:orange`, `#ff000080`) or were
  rewritten (`gold` -> `#ffd700`). MDBridge keeps the text; its picker writes `#rrggbb`. Not fixed in
  0.1.1. *Tested:* `roundtrip/lossy_values`.
- **AP-04: chart formats as set.** Three check boxes saved `["jpg", "svg"]` as `["svg"]`. MDBridge
  keeps the list (validation flags unknown entries). Not fixed in 0.1.1. *Tested:*
  `roundtrip/lossy_values`.
- **AP-05 / AP-19: unknown drop-down values stay.** A stored value outside a drop-down's list was
  saved as the previous project's choice (`visual.renderer = "vulkan"` -> `"software"`) or as `""`.
  MDBridge keeps it (shown as unknown and flagged). *Tested:* `ap/ap05_unknown_choice`,
  `roundtrip/lossy_values`.
- **AP-06: paths written as the file had them (owner decision D8).** MDATK's window made every path
  absolute on load and saved them so (`topology = "/home/.../system.parm7"`), which broke a project
  folder that is moved or shared. MDBridge writes relative paths relative (to the project's folder);
  *Save as* into another folder rebases them (`data/system.parm7` -> `../data/system.parm7` for a copy
  in `sub/`); absolute and `~` paths stay. The paths it reads are resolved against the project's
  folder as before. Not fixed in 0.1.1. *Tested:* `roundtrip/relative_paths`, every round trip.
- **AP-07: empty slots stay empty.** A 4-slot preset with 2 panels was saved with 4 (two untitled
  `contacts/atom_pairs` series). MDBridge saves the panels there are; a slot becomes a panel when its
  editor is changed (and so do the slots before it, which are positional) or a panel is dragged onto
  it. *Tested:* `ap/ap07_empty_slots`, `canvas/swap_into_empty_slot`, `pages/layout_page`.
- **AP-08: a smaller preset keeps the panels.** equal_4 -> single dropped panels 2 ... 4 at once
  and for good. MDBridge keeps them (rendering warns `layout single has 1 slots; ignoring 5 extra
  panel(s)`, as M4's render does) and a larger preset shows them again. *Tested:*
  `ap/ap08_preset_round_trip`, `roundtrip/open_save_full` (6 panels on 4 slots: MDATK dropped two).
- **AP-09: swaps and preset changes with invalid options are refused.** MDATK raised out of the
  signal handler (the swap silently not done, a traceback on the terminal). MDBridge logs `Panels not
  swapped; fix this first: <located message>` (MDATK 0.1.0, c6d8d15) and, for a preset change, `Grid
  preset not changed; fix this first: ...` with the preset put back. *Tested:* `ap/ap09_swap_bad_options`.
- **AP-10: located form errors.** MDATK named Python's errors with the field path (`inputs.stop:
  invalid literal for int() with base 10: 'x'`, `layout: panel options '{bad': Expecting ...`).
  MDBridge names the place as MDATK 0.0.10's `problems.py` does (`Input & cleaning › Frames › Stop:
  'x' is not a whole number`, `... › Target value: 'low' is not a number`, `Layout & panels › Panel 1
  (Structure) › Options: '{bad' is not valid JSON: <json's message>` or `... is not a JSON object`);
  0.1.1 still has the Python texts for numbers, so MDBridge follows its locations and words the
  reasons. *Tested:* `ap/ap10_form_errors`, `ap/ap09_swap_bad_options`.
- **AP-11: Save as asks about the final name.** MDATK appended `.mdatk.toml` after the overwrite
  question, so an existing `name.mdatk.toml` was replaced unasked. MDBridge asks about the final name
  (MDATK 0.1.0, 172df88). *Tested:* `ap/ap11_save_as_suffix`.
- **AP-12: a start-up load error is shown in the window.** MDATK printed `warning: could not load
  ...` on stderr only. MDBridge also logs `ERROR: Could not open <path>: <reason>` and `These are the
  default settings instead (unsaved).` (MDATK 0.1.0, 13d4fea). *Tested:* `ap/ap12_start_broken`.
- **AP-13: unsaved changes are not lost unasked (owner decision D6).** MDATK's New, Open and closing
  discarded edits. MDBridge asks Save / Discard / Cancel when the project differs from what was last
  loaded or saved. *Tested:* `ap/ap13_unsaved_changes`.
- **AP-17: the recent list keeps the settings file's other keys.** MDATK 0.0.4.1 rewrote
  `settings.toml` with `recent` only, deleting newer MDATK's keys (0.0.14's `cat_proof`). MDBridge
  keeps every other key (MDATK 0.0.14) and writes atomically. *Tested:* `ap/ap17_recent_keeps_settings`.
- **AP-18: dialog details in full.** MDATK cut message details at 3 000 characters. *Tested:* the
  recorded dialogs (`detail` uncut).

AP-14 (closing during a job), AP-15 (table edits committed on focus loss), AP-16 (About and BibTeX
with the run-time citations; the `citations` command's cases in `cli_parity.py` cover the keys),
AP-19's offered choices and AP-20 (drawing off the GUI thread) leave no trace in the recorded files
and messages; their packages' Qt tests cover them.

## Interface: the edit model (src/core/project_document.cpp, problems.cpp; milestone 5 W0-B)

MDATK's GUI binds its widgets straight to the `Project` dataclasses (gui/widgets.py `Binder`) and
saves whatever the widgets hold. MDBridge's GUI binds to `core::ProjectDocument` by field path
(`analyses.contacts.cutoff`, `layout.panels.2.options`), a Qt-free model tested on every compiler.
Its deliberate differences from MDATK v0.0.4.1:

- **ED-01: lossless fields** (milestone 5 AP-01 .. AP-05; MDATK 0.1.1 has the same bugs). MDATK's
  spin buttons clamped loaded values to their ranges (`contacts.cutoff` 25 saved as 20, `visual.zoom`
  0.01 as 0.05), optional numbers went through `%g` (`inputs.dt_ps` 0.0123456789 saved as 0.0123457),
  colour buttons rewrote colours (`tab:orange` and `#ff000080` saved as `#000000`, `gold` as
  `#ffd700`), chart formats were three check boxes (`["jpg", "svg"]` saved as `["svg"]`) and combos
  replaced values not in their lists. The document holds exactly what was read: no value is clamped,
  rounded or replaced; numbers show as the project writer writes them (shortest round trip), colours
  and choices are text, lists stay lists; giving a field the text it shows changes nothing, so a field
  the user does not edit is saved back as read. Ranges are the checks' business (ED-03). *Tested:*
  test_project_document.cpp `EveryFieldGivenItsOwnTextChangesNothing` (every field of every fixture
  project and MDATK 0.1.1 file: no change, no undo step), `ValuesMdatksWidgetsRewroteAreKept` (the
  values above, saved byte for byte).
- **ED-02: paths as written; Save As rebases them** (AP-06, owner decision D8). MDATK's GUI wrote
  every path absolute on save (`topology = "demo.parm7"` became `/home/.../demo/demo.parm7`), so a
  moved or shared project folder broke after one save. The document keeps paths as the file has them
  (relative to the project's folder, `~`, absolute) and gives the resolved project (what
  `loadProject` reads) to runs; Save As to another folder rewrites the relative paths (the ones
  `resolvePaths` resolves, planner files and path settings included) so that they still name the same
  files (`inputs/a.parm7` saved from /p into /q is `../p/inputs/a.parm7`; another Windows drive: the
  absolute path); the undo steps are rebased too. *Tested:* `PathsStayAsWritten`,
  `SaveAsRebasesRelativePaths`, `RelativePath`.
- **ED-03: located problems** (AP-10; the presentation of MDATK 0.0.10, `mdatk/problems.py` and
  `config.py` `problems()` of 0.1.1). Every line of `validate()` (which stays MDATK v0.0.4.1's
  blocking rules with MDBridge's fixes, owner decision D7) is also a `Problem` with its field path,
  page and location ("Analyses › Hydrogen bonds › Group A"), in 0.1.1's words where 0.1.1 has the same
  check (`unknown residue group 'X' (defined: Aptamer, Target); did you mean ...?` for v0.0.4.1's
  `contacts: unknown residue group 'X'.`; "Residue group names must be unique: 'A' is also the name of
  group 1."); the command line keeps printing `validate()`'s lines. Field and location are 0.1.1's
  for every problem both have. *Tested:* `ProblemsOfEveryValidateCaseAreLocated` (every
  validate_cases.json input: one problem per validate() line, each on an existing field; against
  MDATK 0.1.1's `problems()` of the same inputs, data/edit-model, tests/parity/make_edit_model_cases.py),
  `DescribeAsMdatkLatest` (describe() of 64 fields), `DocumentProblemsUseResolvedPaths`.
- **ED-04: form errors in words** (AP-10). MDATK showed Python's texts (`inputs.stop: invalid literal
  for int() with base 10: 'x'`, `panel options '{bad': Expecting property name ...`). A field given text
  it cannot hold is refused with a located message and keeps its value: "Input & cleaning › Frames ›
  Stop: 'x' is not a whole number", "'4,5' is not a number", "'maybe' is not true or false", "Options:
  not valid JSON: ...", "options must be a JSON object"; values of the wrong type set by path give the
  project file's messages ("expected a number, got 'far'"). *Tested:* `TextOfFields`, `GetAndSetByPath`.
- **ED-05: undo, redo and a modified flag** (owner decisions D5, D6; AP-13). MDATK has neither. Every
  edit and canvas action is one undo step holding a snapshot of the whole project (a drag is one step:
  `beginGroup` / `endGroup`, `cancelGroup` restores the state before it), at most 200 steps, not saved.
  The project counts as modified while it differs from the last opened or saved state (undoing back to
  it clears the flag). *Tested:* `UndoAndRedo`, `GroupsAreOneStep`, `UndoKeepsAtMostTheLimit`,
  `EditAndListeners`, `ListItems`.
- **ED-06: renaming a residue group renames its uses** (MDATK 0.0.10; follows 0.1.1's
  `Project.rename_group`, `group_references`, `selections.rename_group_tokens`). In v0.0.4.1 a renamed
  group left every analysis, selection, panel and chart that used the old name pointing at nothing.
  Fields holding a group name (analysis groups, representation groups, groups shown), selections that
  are the name or a list of group names, analysis instances in data keys
  (`contacts@Aptamer+Target/atom_pairs`) and `group:<name>` tokens in any other text (planner settings
  and panel options included) follow the rename; `group:<name>er`, other groups and the notes stay; the
  changed fields are returned (as field paths; 0.1.1 writes `groups[1].selection`). Keys MDBridge keeps
  without knowing them (PRJ-09) are written back unchanged, also by a rename. *Tested:*
  `RenameGroupAsMdatkLatest` (14 rename and 5 reference cases: the changed fields, errors and resulting
  project are 0.1.1's), `RenameGroupIsOneStep`.
- **ED-07: a later MDATK's file keeps its comments** (owner decision D9, with PRJ-09). MDBridge writes
  v0.0.4.1's comments, so saving a file of MDATK 0.1.1 replaced 0.1.1's comments (`# RNAdo's plan and
  script wizard`, the options of step types 0.0.4.1 lacks) with older or no ones. For a file whose
  `mdatk_version` is later than 0.0.4.1 (which a save keeps), each entry the file had keeps the comment
  above it in the file; entries it did not have get MDBridge's. Such a file comes back byte for byte;
  the only exception is the blank line MDATK's log move leaves at the end of a project file, which
  MDATK's own next save drops too (golden/split.resaved.mdatk.toml). *Tested:*
  `EveryMdatkProjectRoundTripsByteForByte` (all 18 project files MDATK 0.0.4.1 and 0.1.1 wrote, open ->
  save without an edit: byte for byte, step logs and log files untouched, the resolved project equal to
  `loadProject`'s), `LaterMdatkCommentsAreKept` (an edit changes its line only; recovery copies keep
  them), `HandMadeProjectsKeepTheirValues`.
- **ED-08: recovery copies and external changes** (owner decision D6; MDATK 0.0.11 notices other
  programs' saves). The document can write a recovery copy (a project file with `recovery_of` and
  `recovery_folder` entries, written atomically into `<mdbridge config>/recovery/`; the project file is
  never written without Save), list the copies and reopen one (modified, at its original path, paths
  rebased when its folder differs). `changedOnDisk()` tells whether the project file's settings changed
  since it was opened or saved (a step log appended by a run is no change). *Tested:* `RecoveryCopies`,
  `ChangedOnDisk`.

## Interface: the field widgets (src/app/gui/forms; milestone 5 W0-C)

MDATK's pages are GTK widgets read and written by its `Binder` (gui/widgets.py) when the project is
loaded and collected. MDBridge's field widgets are bound to `core::ProjectDocument` by field path:
each shows the document's value, again after every change (undo, redo, open) without making an edit,
and stores what the user enters as one undo step (text when editing finishes: Enter, leaving the
field, or `FormContext::commitPendingEdits()` before a save; a pick at once). The editor is named
after its field path. Deliberate differences from MDATK v0.0.4.1:

- **FW-01: no clamping, ranges as hints** (AP-01; MDATK 0.1.1 still clamps). MDATK's spin buttons
  rewrote loaded values to their ranges (`contacts.cutoff` 25 saved as 20, `video.width` 7682 as
  7680). A number field holds any text of the field's type; MDATK's spin range is shown in the tooltip
  ("Usual range: 0.5 to 20 (other values are kept)") and a value outside it is flagged as a warning,
  never changed; the arrows step by MDATK's step, rounded to its digits (4.5 + 0.1 is 4.6), and stop at
  the range's ends but never move a value that is already outside into it. Range errors of
  `validate()` flag the field as errors. *Tested:* test_forms.cpp `numbersAreNeverClamped`,
  `numberArrowsStepWithinTheRangeHint`.
- **FW-02: every digit and every colour kept** (AP-02, AP-03; not fixed in 0.1.1). Optional numbers
  show the project writer's shortest round-trip text (MDATK's `%g` saved `dt_ps` 0.0123456789 as
  0.0123457). Colour fields are text with a swatch and a picker (MDATK: a GTK colour button that wrote
  `tab:orange` and `#ff000080` as `#000000`, `gold` as `#ffd700`): any text is kept as typed, the
  swatch shows what Matplotlib's grammar (RN-12) reads, text it cannot read is flagged with its
  message (`Invalid RGBA argument: 'x'`); only a colour picked in the dialog is written `#rrggbb`.
  *Tested:* `optionalNumbersKeepEveryDigit`, `coloursAreKeptAsTyped`.
- **FW-03: lists and unknown choices kept** (AP-04, AP-05). A list field keeps every entry
  (`formats = ["jpg", "svg"]`; a page's check boxes come with W1-C). A fixed combo shows a value it
  does not offer as `<value> (unknown)`, flags it (with `validate()`'s message where it has one, else
  as a warning) and never replaces it, also not by the previous project's choice; picking an item
  stores the item's project value (the item data). Free-text combos (group names, palettes, codecs)
  keep the text as typed: no completion and no case-insensitive swap for an item. Base-pair engine and
  renderer offer MDBridge's choices (AP-19, RN-20). *Tested:* `chartFormatsAreKeptAsAList`,
  `unknownChoicesAreKeptAndShown`, `checkBoxesAndEditableCombos`.
- **FW-04: form errors and problems on the field** (AP-10; flags as MDATK 0.1.1's `set_flag`, colours
  of its `FLAG_COLORS`). Text a field cannot hold is refused with its located message (`Input &
  cleaning › Frames › Stop: 'x' is not a whole number`), stays shown and flagged (tinted, an icon at the
  end of line edits, the message as tooltip and accessible description above the field's hint); the
  document keeps its value; Esc shows the stored value again. The document's problems flag their
  fields (recomputed once per change for all widgets) and tint list cells. *Tested:*
  `formErrorsAreLocatedAndKeepTheValue`, `documentProblemsFlagTheirField`.
- **FW-05: table cells commit on focus loss** (AP-15; adopts MDATK 0.1.0 608f322
  `keep_edit_on_leave`). MDATK 0.0.4.1 dropped text typed into a table cell (groups, trajectories,
  external data, RMSD targets) when the cell was left by a click. A cell is stored on Enter, when it
  loses the focus (a click elsewhere, on a button of the editor) and before a save. Rows are the
  document's items: added (with MDATK's defaults for the row), removed and moved as one undo step each,
  a Browse… of several files as one step; rows without a name or selection are kept (MDATK dropped
  them on collect) and flagged by the checks. *Tested:* `tableCellsCommitOnFocusLossAndEnter`,
  `listEditorRowsAndUndo`, `pathFieldsBrowseWithMdatksFilters`.
- **FW-06: load and save without an edit is the identity.** A form of every field of each of the 20
  fixture projects (MDATK 0.0.4.1 and 0.1.1 files, hand-written ones), every field given Enter without
  typing, saves the file back byte for byte with no undo step (MDATK: `test_widgets_roundtrip_project`
  checks the demo only, and its widgets rewrote values, AP-01 .. AP-05). Only text the user types is
  stored; a multi-line field shows a lone carriage return as a line break but writes it back unchanged
  unless the text is edited. *Tested:* `everyFieldOfEveryFixtureSavesBackByteForByte`.
- **FW-07: the mouse wheel and the keyboard** (MDATK 0.0.4 `test_value_widgets_ignore_the_mouse_wheel`,
  owner decision D10). Spin boxes, combos and sliders ignore the wheel unless they have the keyboard
  focus (MDATK ignored it always) and no longer take the focus from it; every editor and button is
  reachable with Tab and has an accessible name. File dialogs are the platform's (D3) with MDATK's
  filters, each pattern also in lower and upper case (MDATK's match ignores case; Qt's dialogs do not on
  Linux). *Tested:* `valueWidgetsIgnoreTheMouseWheel`, `everyEditorIsNamedAndReachable`.

## Planned (from the improvement register; add each entry here when it lands)

None: every FIX and ASK item of the improvement register (docs/IMPROVEMENTS.md) has landed
(SEL-16 and SEL-11 are kept and documented under "Selections").
