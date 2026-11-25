# MDBridge

Analysis and visualisation of molecular dynamics trajectories of RNA systems: the C++ successor of
[MDATK](https://github.com/htchamblee/MDATK) (whose Python version, v0.0.4.1, is frozen as the
reference implementation). MDBridge ships as one file per platform and never asks its users to
manage a language runtime: a Linux AppImage and a Windows installer.

The plan (modules, libraries, migration order, definition of done) is in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Status: milestone 3 done (`analysis`); next: milestone 4, rendering

**Milestone 1, the skeleton**, proved the delivery path before any science moved:

- **Builds**: CMake, C++20; Qt 6 is the only dependency to install, the rest is fetched (pinned and
  hash-checked) at configure time. CI builds a Linux AppImage on Rocky Linux 9 (glibc 2.34) and a
  Windows installer with MSVC 2022, and tests both as users get them (`.github/workflows/ci.yml`).
- **3D view**: a Vulkan 1.3 window inside the Qt interface, with a software view (QPainter) as
  the fallback on machines without Vulkan 1.3 (older GPUs, virtual machines, remote displays), or
  when Vulkan fails to start. It showed a placeholder scene until milestone 4, which draws an
  opened project's reference-frame structure there (mouse rotation and zoom) above a read-only
  preview of its figure.

**Milestone 2, `core` and `io`**, reads and writes everything MDATK v0.0.4.1 reads and writes
before its analyses start, checked against MDATK itself in CI:

- **Project files**: MDATK's typed parameter model; project files (`.mdatk.toml`) and parameter logs
  (`MDATK_params_log.toml`) read and written byte for byte as MDATK does (comments, step logs, logs
  moved to `<name>.log.toml`, locking), citations and the step log of every run.
- **Systems**: topologies and trajectories read as MDAnalysis reads them (Amber parm7, PDB, GRO;
  DCD, XTC, TRR, Amber NetCDF, Amber ASCII restart / mdcrd), MDAnalysis' selection language and
  MDATK's residue groups, PDB export of frames, MRC / CCP4 density maps (plain or gzip, read,
  written, simulated from a model), mmCIF restraint dictionaries.
- **Command line**: MDATK's `init` (`init-params`), `info`, `list`, `analyze`, `export-pdb` and
  `density info`, with MDATK's options, output, files, messages and exit codes
  (`tests/parity/cli_parity.py` runs both on the same inputs in CI; `analyze` keeps its own cache).
- **Clean room**: MDAnalysis' behaviour (LGPL) is reimplemented from its documented behaviour and
  observed results, never from its source; MDBridge's residue, element and mass tables are rebuilt
  from independent public sources (`data/sources/ATTRIBUTION.md`), so guessed masses and elements,
  some selection macros and PDB atom-name columns can differ from stock MDAnalysis by design
  (`tests/data/README.md`).

**Milestone 3, `analysis`**, computes what MDATK's analyses compute, checked against MDATK itself
in CI:

- **The ten analyses**: RMSD (optimal superposition, fit selections, mass weighting), contacts
  (with the native fraction), H-bonds, Mg2+ coordination (bridges, geometry, outer shell,
  occupancy), stacking, interaction energy (Coulomb + Lennard-Jones, parm7 or generic parameters,
  bonded exclusions), base pairs (MDBridge's own detector), and the structural view's trace,
  proximity and atoms; instances (`rmsd@G`, `contacts@Target+Aptamer`), frame selections and
  periodic cells (`use_pbc`) as MDATK has them.
- **Command line**: `analyze` (MDATK's options, `--only`, `--force`, cache reuse) and `list`
  (the analyses, panel types and layout presets), with MDATK's output and exit codes.
- **Cache**: MDBridge writes its own files (`MDBridge_cache.json`, `<name>-<key>.mdbridge.npz` /
  `.json`, readable with numpy) next to MDATK's, which it never reads or overwrites; the cache
  keys are MDATK's, so a change of input or parameter invalidates the same entries.
- **Frame-parallel runner**: the frames of one pass are processed on several threads
  (`MDBRIDGE_THREADS`, default: the processor's threads, at most 16); results, logs and progress are
  identical to a sequential run.
- **Parity method**: generators in `tests/parity/` run MDATK (and MDAnalysis) in-process on MDATK's
  own micro systems, the demo project and its variants and small generated systems, and record
  every array bit for bit with the decision margins; the C++ tests compare exactly, except for
  documented last-bit tolerances (`tests/data/README.md`, `docs/specs/m3-plan.md` 4.3).
- **Deliberate differences** from MDATK (no Barnaba, pair order of the neighbour search, own cache
  files, ...) are listed in [docs/DIFFERENCES.md](docs/DIFFERENCES.md); the inherited bugs being
  fixed next (periodic distances, degenerate RMSD, reader messages, ...) in
  [docs/IMPROVEMENTS.md](docs/IMPROVEMENTS.md).

Next, milestone 4: rendering (the Vulkan structure view and its representations, charts, figures
and videos). The interface pages, PDB hygiene and the workflow planner follow (milestones 5 to 7);
until then MDBridge's command line offers only the commands it has.

## Using it

```bash
chmod +x MDBridge-0.1.0-x86_64.AppImage
./MDBridge-0.1.0-x86_64.AppImage                          # the interface
./MDBridge-0.1.0-x86_64.AppImage project.mdatk.toml       # ... with a project
./MDBridge-0.1.0-x86_64.AppImage project show project.mdatk.toml   # summary on the command line
./MDBridge-0.1.0-x86_64.AppImage graphics                 # Vulkan devices, or why the software view is used
./MDBridge-0.1.0-x86_64.AppImage init run.mdatk.toml --top sys.parm7 --traj md.nc   # MDATK's commands
./MDBridge-0.1.0-x86_64.AppImage info -p run.mdatk.toml
./MDBridge-0.1.0-x86_64.AppImage export-pdb -p run.mdatk.toml --frame-list "0, last"
./MDBridge-0.1.0-x86_64.AppImage density info map.mrc
./MDBridge-0.1.0-x86_64.AppImage --help
```

On Windows the installer adds MDBridge to the Start menu (`mdbridge-gui.exe`); `mdbridge.exe` is
the command line. Builds of every push are attached to its CI run; tagged versions (`v*`) are
published as GitHub releases.

## Building

See [docs/BUILDING.md](docs/BUILDING.md). In short, with Qt 6.5 or newer (built with Vulkan
support, as the official Qt binaries are):

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/Qt/6.11.3/gcc_64
cmake --build build
ctest --test-dir build
./build/src/app/mdbridge --help
```

## Layout

| Path | Contents |
|---|---|
| `src/core` | data model and project files (no Qt): plain TOML data, MDATK's file formats, JSON |
| `src/io` | topology, trajectory, PDB, density-map and dictionary readers and writers (no Qt) |
| `src/chem` | secondary structure: dot-bracket notation, builtin base-pair detection (no Qt) |
| `src/analysis` | the analyses, their registry, the frame-parallel runner, the cache and the `analyze` pipeline (no Qt) |
| `src/render` | the 3D view: Vulkan 1.3 window, software view, Vulkan probe (Qt GUI, no widgets) |
| `src/app` | the Qt interface and the command line, two thin front ends over the libraries |
| `tests` | unit tests (GoogleTest), fixtures and MDATK's reference reading of them, parity scripts |
| `packaging` | AppImage and Windows installer builds, smoke test, licence texts, icon |
| `tools` | developer tools (shader compilation) |

## Licence

All rights reserved for now: the licence follows a legal review (see [LICENSE](LICENSE) and
[LICENSING.md](LICENSING.md)). Third-party software in the builds:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
