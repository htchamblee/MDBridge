# MDBridge: C++ architecture plan

Status: agreed plan; milestones 1 (skeleton) and 2 (`core` + `io`, with the commands `init`, `info`,
`export-pdb` and `density info`) done, milestone 3 (`analysis`) next. The C++ successor of MDATK is
named **MDBridge** and lives in its own repository (`htchamblee/MDBridge`). The Python package MDATK (v0.0.4.1) is frozen as
the reference implementation and receives bug fixes only.

## Why

MDBridge must never ask a user to manage a language runtime. A chemist who only uses email should
download one file, start it, and work. Python cannot promise that (virtual environments, conda,
minor releases that break compatibility), and the roadmap (real-time Vulkan rendering of systems
with millions of atoms) needs native performance.

## Delivery

| Platform | Binary | Build from source |
|---|---|---|
| Linux | AppImage (one file, no install, no admin rights) | CMake, dependencies fetched at configure time |
| Windows | installer (.exe) | CMake, same |

macOS is not supported. Graphics: Vulkan 1.3 is required for the full 3D view; a simpler
software rasterizer is the fallback (machines without Vulkan 1.3: older GPUs, virtual machines,
remote desktops), with fewer representations and effects but the same panels and videos.
Standard: C++20; minimum compilers GCC 11 (RHEL 9) and MSVC 2022 (Windows; no MSYS2). Build: CMake; dependencies fetched (FetchContent / vcpkg manifest) so
a source build does not pollute the system.

## Toolkit and libraries

- GUI: Qt 6 (Widgets for the sidebar, wizard and planner; a Vulkan window for the 3D view;
  charts drawn in Qt, e.g. QCustomPlot, replacing Matplotlib early)
- Trajectories and topologies: own readers and writers (parm7, PDB, GRO, DCD, XTC/TRR, AMBER NetCDF
  and ASCII, MRC/CCP4, mmCIF), matched to what MDAnalysis reads; chemfiles was dropped in
  milestone 2 because it reads these formats differently and parity needs MDAnalysis' behaviour
- Linear algebra: Eigen; neighbour search: own cell list / nanoflann
- Video: FFmpeg libraries, LGPL-only build (no x264; OpenH264 or AV1 encoders)
- Scripting: embedded Lua for script/wizard generators
- Interim only: an embedded, private Python interpreter inside the bundle for pieces not yet
  ported (never visible to the user, never installed on the system), with permissively
  licensed packages only

Dependency rule: every bundled dependency is permissive (MIT, BSD, Apache, zlib, Boost, ...)
or LGPL, so the licence chosen after legal review (see `LICENSING.md`) is not constrained.
No GPL code is bundled or copied: Barnaba (GPL) is replaced by MDBridge's own base-pair
annotation, and the FFmpeg build is LGPL-only. MDATK's own Python code (same copyright holder)
may be ported freely.

## Modules

Each module is a library with a C++ API; only `app` depends on Qt widgets.

| Module | Contents | Depends on |
|---|---|---|
| `core` | data model (topology, coordinates, selections, residue groups), `.mdatk.toml` project file, step log, citations, settings | - |
| `io` | trajectory/topology readers, NetCDF/parm7, PDB read/write, MRC/CCP4 maps, mmCIF dictionaries | core |
| `chem` | PDB hygiene (L-RNA detection, naming schemes, chains, renumbering), base pairs, dot-bracket, annotation, helix builder, mutations, dictionary renaming | core, io |
| `analysis` | RMSD, contacts, H-bonds, Mg/O-Mg-O, stacking, energy, trace/proximity; parallel over frames; result cache | core, io, chem |
| `render` | Vulkan scene (trace, cartoon, licorice, surfaces, density shadow), panels and charts, video encoding; usable headless for batch videos | core, analysis |
| `workflow` | planner model, step types and checks; generators written in Lua and loaded from `.lua` files (users can add steps) | core, chem |
| `app` | Qt GUI and the CLI: two thin front ends over the same libraries | all |

Analyses run in one pass over the frames, frame-parallel (`runAnalyses`, `analysis.hpp`): one
thread reads the trajectory in order, worker threads run each analysis' per-frame step, and
results go to per-frame slots, so they are byte-identical for any thread count; first-frame steps
(rmsd reference, hbonds hydrogen guess, trace / atoms alignment reference and links) run in a
sequential `prepare` on the first frame. Threads: `RunOptions::threads` (`Pipeline::setThreads`),
else the environment variable `MDBRIDGE_THREADS` (a positive integer; `1` is the sequential run),
else the hardware concurrency capped at 16.

PDB hygiene stays in C++ (it shares chirality and base logic with the builder and analyses);
Lua generators call it rather than reimplement it.

External programs (cpptraj, Amber, GROMACS, Rosetta, SimRNA, RNAMake, Phenix, Coot, acedrg,
eLBOW, VMD) stay external: they are validated, published tools that keep users' work
publishable. MDBridge writes and (on request) runs their scripts.

## Compatibility

- MDATK's `.mdatk.toml` project files open in MDBridge, and the CLI commands keep their names
  and options.
- The analysis cache is not compatible (recomputed).
- Generated scripts must match the Python generators' output for the same settings.
- Paths stored in project files: relative ones resolve against the project file's folder. On
  Linux and macOS this is exactly MDATK's Python (posixpath; CI compares). MDATK never ran on
  Windows, so MDBridge sets its own Windows rules there (`src/core/include/mdbridge/core/ntpath.hpp`,
  pure string functions tested on every platform):
  - `/` and `\` both separate.
  - Absolute: drive paths (`C:\x`, `C:/x`), UNC paths (`\\server\share\x`, `//server/share/x`) and
    root-relative paths (`\x`, `/x`), which take the drive of the project's folder.
  - Drive-relative `C:x`: against the project's folder when that is on drive C (the letter's case
    ignored), otherwise against `C:/`.
  - `~` and `~\x` expand to `%USERPROFILE%` (else `%HOMEDRIVE%%HOMEPATH%`), as Python's ntpath;
    `~user` is left as it is.
  - Resolved paths, absolute ones included, are normalised and written in generic form: `/`
    separators, no `.` or `..` (never above the drive or UNC root), drive letter as written
    (`C:\data\..\a.pdb` becomes `C:/a.pdb`). The project's folder is in the same form
    (`C:/projects/p`), and so are the default output folders against the current directory.
  - Device paths (`\\?\C:\x`, `\\.\pipe\x`) are kept exactly as written.
  - Stored text is compared case-sensitively; nothing is case-folded.

## Migration order (speed first)

1. Skeleton: CMake, Qt 6 + Vulkan window, AppImage and Windows installer in CI; opens an existing
   `.mdatk.toml`. Proves the delivery promise before any science moves.
2. `core` + `io`: project file, trajectory and topology reading, PDB/maps.
3. `analysis`: all analyses, compared against Python on the test trajectories.
4. `render`: Vulkan structure view and representations; Qt charts (retires Matplotlib);
   video encoding.
5. `app`: GUI pages, canvas layout editor, CLI.
6. `chem`: hygiene, base pairs, annotation, builder, dictionaries.
7. `workflow`: Lua generators for every step type, planner, wizard, runner.
8. Remove the embedded Python.

## Definition of done (per module)

A module is done when it reproduces MDATK v0.0.4.1 on the existing test trajectories, PDB files
and workflows: the same numbers within floating-point tolerance, the same generated files
byte for byte (except time stamps), and the same project files. The Python test suite is the
specification; each Python test gets a C++ (or CLI-level) counterpart.

Parity is the floor, not the goal: MDBridge is meant to surpass MDATK. Where MDATK (or MDAnalysis,
whose behaviour it inherits) has a bug, a crash, a misleading message or a numerically fragile
result, MDBridge fixes it instead of reproducing it. Each such fix is a deliberate difference: it
is listed in `docs/DIFFERENCES.md` with the reason, and the parity tests check MDBridge's improved
behaviour there (the expected data is adjusted on the test side, as for the independent data tables,
`tests/parity/independent_tables.py`). Behaviour that is merely a convention (output formats, file
names, orders users rely on, option names) stays as MDATK has it.

## Repository and parity checks

Own repository, `htchamblee/MDBridge`, with its own licence, releases, CI and issues (MDATK
stays GPLv3; MDBridge's licence follows the legal review, all rights reserved until then).
Parity is checked by comparing outputs: MDBridge's CI checks out `htchamblee/MDATK` at the
frozen tag `v0.0.4.1`, runs both on the same inputs and compares numbers (floating-point
tolerance), generated scripts and files (byte for byte except time stamps) and project files.
Shared test inputs are copied into MDBridge or read from that checkout; the user's private
research structures are never committed to either repository.
