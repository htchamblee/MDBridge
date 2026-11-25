# MDBridge

Analysis and visualisation of molecular dynamics trajectories of RNA systems: the C++ successor of
[MDATK](https://github.com/htchamblee/MDATK) (whose Python version, v0.0.4.1, is frozen as the
reference implementation). MDBridge ships as one file per platform and never asks its users to
manage a language runtime: a Linux AppImage and a Windows installer.

The plan (modules, libraries, migration order, definition of done) is in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Status: milestone 4 done (`render`); next: milestone 5, the interface (`app`)

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

**Milestone 4, `render`**, draws what MDATK draws (charts, figures, videos and the structure view),
without Matplotlib, Python or an `ffmpeg` executable, on a desktop or a headless cluster node:

- **Figures**: MDATK's multi-panel figure (the structure view, time series, dual axis, heatmap,
  contact map, sequence and secondary-structure panels; the 16 layout presets; light and dark
  themes) on MDBridge's own figure model and QPainter backends (PNG, SVG, PDF, JPEG, TIFF, EPS,
  PS). Matplotlib's colours, colormaps, tick locators and formatters and its constrained layout are
  ported with their notice; text is drawn in the bundled DejaVu fonts (Matplotlib's), so a figure
  looks the same on every machine and needs no fonts, fontconfig or display (Qt's `offscreen`
  platform).
- **Standalone charts**: MDATK's series and matrix charts (`output.charts`: sizes, styles,
  smoothing, targets, title templates, file names).
- **Structure view and 3D representations**: MDATK's cartoon, ribbon, tube, licorice, vdW, surface
  and trace styles with its colour schemes, glow, ions, residue tags, centring and the projected
  density-map shadow, drawn by two interchangeable renderers: a port of MDATK's software
  rasteriser (multi-threaded, the same pixels as MDATK's) and a Vulkan 1.3 renderer (offscreen,
  no window needed) that draws the same picture on the GPU; `visual.renderer` chooses (`auto`: the
  GPU when there is one).
- **Video**: MP4 encoded in process with an own LGPL-only build of FFmpeg's libraries and the
  royalty-free AV1 encoder SVT-AV1 (owner decision D1; `mpeg4` as the fallback), bundled with the
  programs together with the corresponding FFmpeg source on every release. The project files'
  `codec = "libx264"` (GPL, never bundled) maps to AV1 with a warning; frames are composed on
  several threads and encoded in order.
- **Command line**: MDATK's `render` (full run: analyses, video, charts, PDBs, params log, step log),
  `charts` and `preview`, with MDATK's options, output, files, messages and exit codes; project
  settings are checked before any analysis runs. `probe-video FILE` reads back what a video holds.
- **Interface**: the 3D view shows an opened project's reference-frame structure with its
  representations (mouse rotation and zoom; Vulkan, or the software renderer without Vulkan 1.3),
  above a read-only preview of its figure on MDATK's example data. The layout editor comes with
  milestone 5.
- **Parity method** (owner decision D10, contents only): generators in `tests/parity/` record what
  MDATK's figures contain (the data, colours, ticks, labels, legends, images and colour bars of
  every panel at the first, middle and last frame, the layout's boxes), what its software renderer
  and molecule geometry produce, and what its `render` / `charts` / `preview` commands print and
  write; the C++ tests compare contents exactly (numbers to 1e-9, renderer pixels bit for bit) and
  layout boxes within 1.5 % of the figure size; videos by their metadata. There are no golden
  reference images.
- **Deliberate differences** (unhinted text metrics, exact canvas sizes, AV1 instead of H.264,
  BT.709 colour tags, heatmap cells at non-uniform times, NaN-aware smoothing, settings checked up
  front, the ion and chain legends, the renderer's halo seam, ...) are listed with their tests in
  [docs/DIFFERENCES.md](docs/DIFFERENCES.md); the plan and its RN-xx items in
  [docs/specs/m4-plan.md](docs/specs/m4-plan.md).

Next, milestone 5 (`app`, docs/ARCHITECTURE.md "Migration order"): the interface pages and the
canvas layout editor (dragging panels, rotating the model in the figure, moving the density map).
PDB hygiene and base-pair annotation (`chem`) and the workflow planner follow (milestones 6 and 7);
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
./MDBridge-0.1.0-x86_64.AppImage render -p run.mdatk.toml                 # analyses, video, charts, PDBs
./MDBridge-0.1.0-x86_64.AppImage charts -p run.mdatk.toml --formats png,svg
./MDBridge-0.1.0-x86_64.AppImage preview -p run.mdatk.toml --frame 10 --output frame10.png
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
| `src/render` | figures, charts, panels, the molecule's 3D geometry, the software and Vulkan renderers, video encoding and the 3D view (`core`, `molecule`, `video`: no Qt; the rest Qt GUI, no widgets) |
| `src/app` | the Qt interface and the command line, two thin front ends over the libraries |
| `tests` | unit tests (GoogleTest), fixtures and MDATK's reference reading of them, parity scripts |
| `packaging` | AppImage and Windows installer builds, the LGPL FFmpeg build, smoke test, licence texts, icon |
| `tools` | developer tools (shader compilation, Qt and lavapipe download, Matplotlib table generation) |

## Licence

All rights reserved for now: the licence follows a legal review (see [LICENSE](LICENSE) and
[LICENSING.md](LICENSING.md)). Third-party software in the builds:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
