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
basepairs), so caches written before recompute. *Tested:* the expected data of every case with
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
- **`info` residue-name ranking:** equal counts are ordered by name (MDATK: numpy's argsort order,
  which depends on the CPU's sort kernel).
- **parm7 end of file:** the StopIteration carries "<file>: the parm7 file ends unexpectedly" instead
  of an empty message.
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

## Command line

- **No tracebacks.** Where MDATK crashes with a Python traceback (`info` with a parm7 file alone,
  corrupt gzip maps, map dimensions too large to allocate), MDBridge prints `error: <message>` and exits
  with 1. *Tested:* cli_cases.json (`mdatk_crash` cases), tests/cli/test_commands.cpp.
- Six gzip error texts of `density info` are worded as in mrc.hpp (compared by exit status only).
- JSON presets: `init` with a `.json` preset reports that MDBridge uses TOML project files.

## Planned (from the improvement register; add each entry here when it lands)

Degenerate RMSD (AN-16); setup-order independence (AN-13, SEL-17);
`prop` operand reversal (SEL-01, SEL-02); per-atom element fallback (SEL-13); reader fixes (RD-*,
LD-*); validation of every cutoff (PRJ-01); readable messages in place of Python's internal texts
(SEL-03/04/08, RD-*, CLI-01/03/04).
