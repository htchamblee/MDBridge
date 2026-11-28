# MDBridge

Analysis and visualisation of molecular dynamics trajectories of RNA systems: the C++ successor of
[MDATK](https://github.com/htchamblee/MDATK) (whose Python version, v0.0.4.1, is frozen as the
reference implementation). MDBridge ships as one file per platform and never asks its users to
manage a language runtime: a Linux AppImage and a Windows installer.

The plan (modules, libraries, migration order, definition of done) is in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Status: milestone 6 done (`chem`); next: milestone 7 (`workflow`)

**Milestone 1, the skeleton**, proved the delivery path before any science moved:

- **Builds**: CMake, C++20; Qt 6 is the only dependency to install, the rest is fetched (pinned and
  hash-checked) at configure time. CI builds a Linux AppImage on Rocky Linux 9 (glibc 2.34) and a
  Windows installer with MSVC 2022, and tests both as users get them (`.github/workflows/ci.yml`).
- **3D view**: a Vulkan 1.3 window inside the Qt interface, with a software view (QPainter) as
  the fallback on machines without Vulkan 1.3 (older GPUs, virtual machines, remote displays), or
  when Vulkan fails to start. It showed a placeholder scene until milestone 4, which drew an
  opened project's reference-frame structure there above a read-only preview of its figure; since
  milestone 5 the figure canvas's structure panel is the 3D view.

**Milestone 2, `core` and `io`**, reads and writes everything MDATK v0.0.4.1 reads and writes
before its analyses start, checked against MDATK itself in CI:

- **Project files**: MDATK's typed parameter model; project files (`.mdatk.toml`) and parameter logs
  (`MDATK_params_log.toml`) read and written byte for byte as MDATK does (comments, step logs, logs
  moved to `<name>.log.toml`, locking), citations and the step log of every run.
- **Systems**: topologies and trajectories read as MDAnalysis reads them (Amber parm7, PDB, GRO;
  DCD, XTC, TRR, Amber NetCDF, Amber ASCII restart / mdcrd), MDAnalysis' selection language and
  MDATK's residue groups, PDB export of frames, MRC / CCP4 density maps (plain or gzip, read,
  written, simulated from a model), mmCIF restraint dictionaries.
- **Command line**: MDATK's `demo`, `init` (`init-params`), `info`, `list`, `clean`, `analyze`,
  `export-pdb` and `density info`, with MDATK's options, output, files, messages and exit codes
  (`tests/parity/cli_parity.py` runs both on the same inputs in CI; `analyze` keeps its own cache).
  `demo` writes MDATK's synthetic system byte for byte; cleaning runs cpptraj (directly or through
  `conda run`) with MDATK's scripts and reuses MDATK's cleaned files. `citations` lists what a
  project cites: its settings' references and the run-time ones of its cached results (`--json`,
  `--bibtex FILE`).
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
  ported with their notice; text is drawn in the bundled IBM Plex Serif (with IBM Plex Math for the
  few symbols it lacks; the interface uses the same font), so a figure looks the same on every
  machine and needs no fonts, fontconfig or display (Qt's `offscreen`
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
- **Interface**: the 3D view showed an opened project's reference-frame structure with its
  representations above a read-only preview of its figure on MDATK's example data (replaced by
  milestone 5's window).
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

**Milestone 5, `app`**, is MDATK's analysis window in Qt 6 Widgets, checked against MDATK's own GTK
window driven by the same action scripts (`tests/parity/make_gui_cases.py`; 42 recorded cases):

- **Pages**: MDATK's six settings pages in its order and words (*Input & cleaning*, *Residue groups*,
  *Analyses*, *Structure view*, *Layout & panels*, *Outputs*), every field bound to the project by
  its key; a field you do not change is saved back as it was read (no clamping, every digit, colours
  and unknown choices kept), paths stay as the file has them (relative stays relative; *Save as*
  rebases them), files picked inside the project folder are stored relative. Settings are checked
  as you type: page badges, flagged fields and a problems list that jumps to the field (MDATK
  0.0.10's presentation). Group renames follow through the project; *Groups from chains* (D15).
- **Canvas layout editor**: the figure on MDATK's example data, drawn on a worker thread; drag a
  panel onto another to swap them, drag the structure to rotate it, the wheel (or + / - with the
  canvas focused) zooms, the arrow keys turn it, map mode moves and rotates the density map,
  double-click opens a slot's editor. The figure's structure panel is the 3D view (Vulkan while
  dragging on a GPU, else the software renderer).
- **Actions**: *Render* (analyses, video, charts, PDBs, params log; what `mdbridge render` writes),
  *Render still*, *Export PDBs*, *Re-render from cache*, *Run cleaning now* (cpptraj), *Inspect
  system*, *Centre on model*, *Create demo system…*, *Open output folder* / *Open video*, project
  log, About and citations (with the run's references) and BibTeX export; one job at a time with
  progress and Cancel, closing during a job asks and waits for its clean-up.
- **Undo and recovery**: undo / redo of every edit and canvas action (a drag is one step; Ctrl+Z,
  Ctrl+Shift+Z), a modified marker and Save / Discard / Cancel before anything discards changes, a
  recovery copy every 2 minutes offered at the next start (the project file is written only by
  Save), a notice when another program changes the file.
- **Keyboard and screens**: MDATK's shortcuts plus Ctrl+N, Ctrl+W, Ctrl+Q, Ctrl+L (log), F5 and
  Ctrl+1 ... Ctrl+6 (pages); every control reachable with Tab and named for screen readers; every
  page fits the side bar at 100, 150 and 200 % scale; on Windows paths are shown with `\`.
- **Fonts**: IBM Plex Serif in the interface, figures, charts and videos (IBM Plex Math for missing
  symbols, Plex Mono for fixed-width and sequence text), bundled, so text looks the same everywhere.
- **Settings**: the recent list is shared with MDATK on Linux (`~/.config/mdatk/settings.toml`, its
  other keys kept); MDBridge's own settings and recovery copies are in `~/.config/mdbridge`
  (`%APPDATA%\MDBridge` on Windows).
- **Packaging**: the AppImage uses the host's desktop libraries (GL, fontconfig, freetype, X11,
  xcb) and supplies its own copy only of one the host lacks (the GL dispatch libraries as one set);
  CI's smoke test drives the shipped programs' window (`mdbridge gui PROJECT --script FILE`): a
  field edited and saved, a still rendered, every page photographed, on Ubuntu, a bare Rocky Linux
  9 and the installed Windows program, whose Start-menu entry is checked.
- **Deliberate differences** (the AP items: MDATK window bugs fixed, and the window's additions) are
  listed by id with their tests in [docs/DIFFERENCES.md](docs/DIFFERENCES.md) ("Interface: index of
  the ids"); the plan in [docs/specs/m5-plan.md](docs/specs/m5-plan.md).

**Milestone 6, `chem`**, is MDATK's chemistry layer: PDB clean-up, base-pair annotation, the
structure builder, restraint dictionaries and the density fit, checked against MDATK itself in CI:

- **PDB clean-up** (`pdb clean`, `pdb check`, `pdb to-lrna`): ChimeraX repair, L-RNA detection from
  the sugar geometry (`--l-ranges` forces it), residue names in the `amber`, `ccd`, `rosetta` or
  `pdb` scheme (`--scheme`), chain ids from a reference PDB, a list or the TER segments, split,
  combined and sorted chains, renumbering, SEQRES and FASTA, 5'-phosphates removed (`--strip-5p`).
  `pdb check` prints the chains, sequences, names and sugar chirality as JSON.
- **Annotation** (`pdb annotate`): REMARK 950 base-pair and dot-bracket records from MDBridge's
  builtin detector. Barnaba (GPL) is not used: `--engine barnaba` runs the builtin detector with a
  warning.
- **Structure builder** (`build helix`, `build mutate`, `build fit`): ideal A-form RNA helices,
  right- or left-handed (`--chirality D` / `L`), single strands or duplexes with mismatches and G·U
  pairs, in any naming scheme, with hydrogens (`--hydrogens`) and a 5'-phosphate (`--phosphate5`) on
  request; base mutations that keep the glycosidic frame and change the paired partner too, with
  hydrogens on the new base when the model has them; rigid-body fits of a model into a density map
  (`--search global` or `local`). The other pairings (`--pairing hh`, `hoogsteen`) and the
  pyrimidine-motif triplex (`--triplex hoogsteen`) are built by a fitted builder: hh as PDB 4JRD, the
  Hoogsteen duplex and the triplex with the base triples, stacking and screw of RNA triple helices
  (PDB 3P22, 4PLX, 6SVS, 7JNH), the C+ of G·C+ Hoogsteen pairs as the CCD's CH; the kinds RNA is not
  known to form (`reverse-wc`, `reverse-hoogsteen`, the purine-motif triplex) and sequences a kind
  cannot pair (U or C on strand 1 of the Hoogsteen kinds, other bases than A in hh) are refused
  (docs/DIFFERENCES.md CH-16 ... CH-21).
- **Restraint dictionaries** (`dict list`, `dict smiles`, `dict finish`, `dict pubchem`): the SMILES
  of the L-ribonucleotides for acedrg and eLBOW, and their dictionaries renamed to the atom names of
  the wwPDB Chemical Component Dictionary (CCD) components 0A, 0C, 0G, 0U and their 5' / 3'
  variants. `dict pubchem` asks PubChem for a compound, only when run.
- **Density fit** (`density fit`, `density simulate`): MDATK's rigid-body fit of the structure into
  the project's map (a global search over rotations with FFT translations, then Powell refinement;
  `--save` stores the placement and keeps the project's paths as they were), and maps simulated from
  a frame. `render` and `preview` fit a shown map that has no placement first (`density.auto_fit`).
  The search runs on several threads; the result does not depend on their number.
- **Windows**: *Global search* and *Refine* on the Structure view page run the fit as a job with
  progress and Cancel; the placement found is one undo step. *PDB clean-up…* and *Structure
  builder…* (pages *Ideal helix*, *Mutate*, *Fit to map*) in the Workflow menu write what the
  commands write for the same options.
- **Reference data**: the builder's and the dictionaries' tables come from primary sources in
  `data/sources/` (the CCD components A, C, G, U and 0A, 0C, 0G, 0U, PDB 4JRD, the Arnott A-RNA fibre
  coordinates), checked against MDATK's embedded tables. Where they disagree MDBridge follows the
  source: guanine's H21 / H22 are named as the CCD names them.
- **Numerics**: the SciPy pieces MDATK uses (the FFT, ndimage filters, Powell's method,
  `least_squares`) in `src/numerics`: SciPy's own FFT (ducc0) vendored, the rest ported, each tested
  against SciPy.
- **Deliberate differences** (the CH items: residue numbers that do not fit the PDB columns refused,
  mutation specs without a chain, hydrogens of mutated bases, non-finite map voxels, `density fit
  --save` keeping relative paths, missing files named, ...) are listed by id with their tests in
  [docs/DIFFERENCES.md](docs/DIFFERENCES.md) ("Chemistry"); the plan in
  [docs/specs/m6-plan.md](docs/specs/m6-plan.md).

Next, milestone 7 (`workflow`, docs/ARCHITECTURE.md "Migration order";
[docs/specs/m7-plan.md](docs/specs/m7-plan.md)): the script generators of every step type in Lua,
the Project Planner, the script wizard and the runner (RNAdo). The planner's and the wizard's menu
entries are there, disabled, until then.

## Using it

```bash
chmod +x MDBridge-0.1.0-x86_64.AppImage
./MDBridge-0.1.0-x86_64.AppImage                          # the analysis window
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
./MDBridge-0.1.0-x86_64.AppImage pdb clean model.pdb --scheme amber --renumber -o clean.pdb
./MDBridge-0.1.0-x86_64.AppImage build helix --seq GGCAUGCC --chirality L -o l_helix.pdb
./MDBridge-0.1.0-x86_64.AppImage density fit -p run.mdatk.toml --map map.mrc --resolution 6 --save
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
| `src/numerics` | the SciPy pieces MDATK uses: FFT, ndimage filters and interpolation, Powell, `least_squares`, pair search (no Qt) |
| `src/chem` | PDB hygiene, base pairs and dot-bracket notation, annotation, the structure builder, mutations, restraint dictionaries, the density fit (no Qt) |
| `src/analysis` | the analyses, their registry, the frame-parallel runner, the cache and the `analyze` pipeline (no Qt) |
| `src/render` | figures, charts, panels, the molecule's 3D geometry, the software and Vulkan renderers, video encoding and the 3D view (`core`, `molecule`, `video`: no Qt; the rest Qt GUI, no widgets) |
| `src/app` | the Qt interface (with the PDB clean-up and structure builder windows, `gui/tools`) and the command line, two thin front ends over the libraries |
| `tests` | unit tests (GoogleTest), fixtures and MDATK's reference reading of them, parity scripts |
| `packaging` | AppImage and Windows installer builds, the LGPL FFmpeg build, smoke test, licence texts, icon |
| `third_party` | vendored source: SciPy's FFT (ducc0), unmodified and checked by SHA-256 |
| `data` | the bundled fonts, the icon and desktop entry, and the public-source extracts the generated tables are made from (`data/sources`) |
| `tools` | developer tools (shader compilation, Qt and lavapipe download, the extracts of `data/sources` and the tables generated from them, ducc0 vendoring) |

## Licence

All rights reserved for now: the licence follows a legal review (see [LICENSE](LICENSE) and
[LICENSING.md](LICENSING.md)). Third-party software in the builds:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
