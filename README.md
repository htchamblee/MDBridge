# MDBridge

Analysis and visualisation of molecular dynamics trajectories of RNA systems: the C++ successor of
[MDATK](https://github.com/htchamblee/MDATK) (whose Python version, v0.0.4.1, is frozen as the
reference implementation). MDBridge ships as one file per platform and never asks its users to
manage a language runtime: a Linux AppImage and a Windows installer.

The plan (modules, libraries, migration order, definition of done) is in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Status: milestone 1, the skeleton

Milestone 1 proves the delivery path before any science moves:

- **Builds**: CMake, C++20; Qt 6 is the only dependency to install, the rest is fetched (pinned and
  hash-checked) at configure time. CI builds a Linux AppImage on Rocky Linux 9 (glibc 2.34) and a
  Windows installer with MSVC 2022, and tests both as users get them (`.github/workflows/ci.yml`).
- **3D view**: a Vulkan 1.3 window inside the Qt interface, with a software view (QPainter) as
  the fallback on machines without Vulkan 1.3 (older GPUs, virtual machines, remote displays), or
  when Vulkan fails to start. It shows a placeholder scene until the `render` milestone.
- **Project files**: opens existing MDATK project files (`.mdatk.toml`) and parameter logs
  (`MDATK_params_log.toml`), including step logs moved to `<name>.log.toml`. They read exactly as
  MDATK v0.0.4.1 reads them: CI compares against MDATK itself on every push.

Not yet: the typed parameter model and writing project files (milestone 2), trajectories,
analyses, rendering, charts, videos and the MDATK command set (milestones 2 to 7).

## Using it

```bash
chmod +x MDBridge-0.1.0-x86_64.AppImage
./MDBridge-0.1.0-x86_64.AppImage                          # the interface
./MDBridge-0.1.0-x86_64.AppImage project.mdatk.toml       # ... with a project
./MDBridge-0.1.0-x86_64.AppImage project show project.mdatk.toml   # summary on the command line
./MDBridge-0.1.0-x86_64.AppImage graphics                 # Vulkan devices, or why the software view is used
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
