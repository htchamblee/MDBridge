# MDBridge

Analysis and visualisation of molecular dynamics trajectories of RNA systems: the C++ successor of
[MDATK](https://github.com/htchamblee/MDATK) (whose Python version, v0.0.4.1, is frozen as the
reference implementation). MDBridge ships as one file per platform and never asks its users to
manage a language runtime: a Linux AppImage and a Windows installer.

The plan (modules, libraries, migration order, definition of done) is in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Status: milestone 2 done (`core` + `io`); next: milestone 3, the analyses

**Milestone 1, the skeleton**, proved the delivery path before any science moved:

- **Builds**: CMake, C++20; Qt 6 is the only dependency to install, the rest is fetched (pinned and
  hash-checked) at configure time. CI builds a Linux AppImage on Rocky Linux 9 (glibc 2.34) and a
  Windows installer with MSVC 2022, and tests both as users get them (`.github/workflows/ci.yml`).
- **3D view**: a Vulkan 1.3 window inside the Qt interface, with a software view (QPainter) as
  the fallback on machines without Vulkan 1.3 (older GPUs, virtual machines, remote displays), or
  when Vulkan fails to start. It shows a placeholder scene until the `render` milestone.

**Milestone 2, `core` and `io`**, reads and writes everything MDATK v0.0.4.1 reads and writes
before its analyses start, checked against MDATK itself in CI:

- **Project files**: MDATK's typed parameter model; project files (`.mdatk.toml`) and parameter logs
  (`MDATK_params_log.toml`) read and written byte for byte as MDATK does (comments, step logs, logs
  moved to `<name>.log.toml`, locking), citations and the step log of every run.
- **Systems**: topologies and trajectories read as MDAnalysis reads them (Amber parm7, PDB, GRO;
  DCD, XTC, TRR, Amber NetCDF, Amber ASCII restart / mdcrd), MDAnalysis' selection language and
  MDATK's residue groups, PDB export of frames, MRC / CCP4 density maps (plain or gzip, read,
  written, simulated from a model), mmCIF restraint dictionaries.
- **Command line**: MDATK's `init` (`init-params`), `info`, `export-pdb` and `density info`, with
  MDATK's options, output, files, messages and exit codes (`tests/parity/cli_parity.py` runs both
  on the same inputs in CI).
- **Clean room**: MDAnalysis' behaviour (LGPL) is reimplemented from its documented behaviour and
  observed results, never from its source; MDBridge's residue, element and mass tables are rebuilt
  from independent public sources (`data/sources/ATTRIBUTION.md`), so guessed masses and elements,
  some selection macros and PDB atom-name columns can differ from stock MDAnalysis by design
  (`tests/data/README.md`).

Next, milestone 3: the analyses (RMSD, contacts, H-bonds, Mg2+ coordination, stacking, energy, base
pairs) and MDATK's `analyze` command, compared with MDATK on the test trajectories. Rendering,
charts and videos, the interface pages, PDB hygiene and the workflow planner follow (milestones 4
to 7); until then MDBridge's command line offers only the commands it has.

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
| `src/render` | the 3D view: Vulkan 1.3 window, software view, Vulkan probe (Qt GUI, no widgets) |
| `src/app` | the Qt interface and the command line, two thin front ends over the libraries |
| `tests` | unit tests (GoogleTest), fixtures and MDATK's reference reading of them, parity scripts |
| `packaging` | AppImage and Windows installer builds, smoke test, licence texts, icon |
| `tools` | developer tools (shader compilation) |

## Licence

All rights reserved for now: the licence follows a legal review (see [LICENSE](LICENSE) and
[LICENSING.md](LICENSING.md)). Third-party software in the builds:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
