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
  extra point switched a whole system to name guessing; `CL1` became C).
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

## Command line

- **No tracebacks.** Where MDATK crashes with a Python traceback (`info` with a parm7 file alone,
  corrupt gzip maps, map dimensions too large to allocate), MDBridge prints `error: <message>` and exits
  with 1. *Tested:* cli_cases.json (`mdatk_crash` cases), tests/cli/test_commands.cpp.
- Six gzip error texts of `density info` are worded as in mrc.hpp (compared by exit status only).
- JSON presets: `init` with a `.json` preset reports that MDBridge uses TOML project files.

## Planned (from the improvement register; add each entry here when it lands)

Periodic distances and minimum images for every cell (AN-01–AN-04); degenerate RMSD (AN-16); base
pairs and stacking on wrapped trajectories (AN-11, AN-12); setup-order independence (AN-13, SEL-17);
reader fixes (RD-*, LD-*); validation of every cutoff (PRJ-01); readable messages in place of
Python's internal texts (RD-*, CLI-01/03/04).
