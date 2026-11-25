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
- **No Barnaba.** Barnaba is GPL. With `analyses.basepairs.engine = "barnaba"`, MDBridge warns and runs
  its builtin detector (meta `engine = "builtin"`, no `lw_first_frame`, `noncanonical` has no effect).
  The params, and so the cache key, keep the project's values. *Tested:* tests/analysis/test_basepairs.cpp.

## Analysis cache

- **Own cache format.** MDBridge writes its own archives and manifest (`MDBridge_cache.json`,
  `<name>-<key>.npz/.json`; src/analysis/include/mdbridge/analysis/cache.hpp). It does not read
  MDATK's cache: results are recomputed once. Data keys, cache keys and the JSON descriptions keep
  MDATK's meaning. *Tested:* tests/analysis/test_pipeline.cpp; cli_parity.py compares cache files by
  meaning.

## Order effects of the neighbour search

MDATK's pairs without a periodic box come in the traversal order of scipy's k-d tree, which depends on
how scipy was built. MDBridge returns them sorted by (i, j). The pair set and every distance are
identical. Visible consequences:
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
