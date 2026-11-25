# Milestone 4 plan: the `render` module

Goal (ARCHITECTURE "Migration order" 4 + "Definition of done"): MDATK v0.0.4.1's presentation
layer (`mdatk/render/`, 4 600 lines) in C++: the multi-panel figure (7 panel types, 16 layout
presets), the structure view with its 3D representations (Vulkan 1.3, software fallback ported from
MDATK's numpy renderer), standalone charts (PNG / SVG / PDF, retiring Matplotlib), the MP4 video
(bundled LGPL FFmpeg libraries instead of an `ffmpeg` executable), and the commands `render`,
`charts` and `preview` with MDATK's options, output and files. Usable headless (cluster nodes, the
bare Rocky system). GUI pages and the canvas layout editor stay milestone 5; M4 gives the GUI the
real structure in its 3D view and a read-only figure preview.

Licensing: MDATK code (owner's) ported freely. Matplotlib (its own licence, PSF-style, permissive)
and scipy (BSD-3) algorithms may be ported with notices (`THIRD_PARTY_NOTICES.md`,
`packaging/licenses/`, LICENSING table). No MDAnalysis code is involved in rendering (trace / atoms
come from the analyses of M3). Qt Charts and Qt Graphs are GPL v3 / commercial in open-source Qt 6,
and QCustomPlot is GPL / commercial: none of them may be used (ARCHITECTURE's "e.g. QCustomPlot"
must be corrected in this milestone).

---------------------------------------------------------------------------------------------------

## 1. Inventory of MDATK's rendering

### 1.1 Figure composition (`render/video.py: build_figure`, `layout.py`)

- **Size**: `figure_size(w, h)` rounds both down to even (≥ 2) for yuv420p; figure = `w/dpi x h/dpi`
  inches at `dpi` (video defaults 1920 x 1080, dpi 120). DPI sets text and line sizes: points ->
  pixels = `pt * dpi / 72`. Facecolor `theme.page`.
- **Grid**: `layout.preset` -> `(rows, cols, slots)`, 16 presets (`single`, `large_plus_2/3/4`,
  `large_top_2/3`, `equal_2..6`, `*_strip` x 5); `HEIGHT_RATIOS` for the strip presets; panels fill
  slots in order; extra panels: `WARNING layout X has N slots; ignoring M extra panel(s)`. Already in
  MDBridge as data (`core/panels.hpp`, `panel_tables.inc`).
- **Layout engine**: Matplotlib "constrained" (`h_pad=0.08, w_pad=0.08, hspace=0.06, wspace=0.04`
  inches / fractions), which sizes every axes around its tick labels, axis labels, title, colorbar
  and twin axis; `fig.suptitle(visual.title, fontsize=1.3 * font_size)`; optional watermark
  bottom-right (GUI "Example data – layout preview"). After the first draw the layout is frozen
  (`set_layout_engine("none")`), then each panel's `after_layout()` runs (structure: equal aspect by
  widening the data window; sequence: re-flow letters to the real pixel size).
- **Frame loop** (`render_video`): blitting: everything not in `animated_artists()` is drawn once
  into a background; per frame `update(pos)` moves the artists, the background is restored and only
  the animated artists are drawn; the RGBA buffer goes to FFmpeg.
- **RenderContext**: theme, time axis in `inputs.time_unit` (`dataset.time(unit)`,
  `time_label(unit)`), `font_size`, `format_time(pos)` = `frame N` or
  `t = 1.234 ns  (frame N)` (`us` shown as `µs`).

### 1.2 Theme and chrome (`theme.py`)

Two themes (light / dark) of 7 tokens (surface, page, ink, ink2, muted, grid, axis); 8-colour
categorical palette per theme assigned in fixed order (`k % 8`, never cycled to new colours);
13-step sequential blue ramp from the surface (`LinearSegmentedColormap.from_list`), reversed for
dark; `ordinal(n)` discrete ramp for class matrices (class 0 recedes). `style_axes`: no top/right
spines, spines `axis` colour 0.8 pt, ticks `muted` 0.6 pt x 3 pt, tick labels `ink2` at
0.85 x font, axis labels 0.9 x, title = font size, hairline grid 0.6 pt under the data.
`style_legend`: surface face, grid edge 0.6 pt, text `ink2` at 0.8 x font.

### 1.3 Panel types (`panels.py`, `seqpanels.py`; registry order = `list` order)

Common: title = `spec.title` or the panel's automatic title, left-aligned, `ink`; a `KeyError` in
setup shows centred `No data: <key>` in `muted` on a bare frame (the run continues).

| Type | Data | Axes / marks | Legend, colour bar | Options | Animated per frame |
|---|---|---|---|---|---|
| `structure` | `trace` (+ `atoms` with representations, + `proximity` for proximity colour / glow); density map | bare square data window (bounding sphere of all frames / zoom, pad `0.06*view + 2/max(zoom,1)`), equal aspect; points (size `ps*(0.55+0.9*depth)`, surface edge 0.6), links as round-capped segments (width by depth, colour mean of ends), depth fog 0.45 toward surface, two glow halos (`9*ps`, `3.8*ps`, alpha `0.20 / 0.45 * glow`), ions (`2.2*ps`, `ion_color`), residue tags pushed outward, time label top-left | lower right: `Chain k+1` (chain mode) or groups present, `≤ 4 Å of Target` (proximity), `Glow (fade)`, `Mg²⁺` | - (visual.* settings) | everything moving, 3D image, front shadow, legend on top |
| `timeseries` | 1-8 series | lines 1.6 pt in palette order, dot 6 pt with surface edge, dashed target line `(0,(4,3))` (`Target (v)`), cursor `axvline` ink alpha 0.7; xlim `t0..t-1` (1 frame: `t0..t0+1`); integer y ticks for counts; y label = axis label (1 series) or common unit | `best` when >1 series or a target | `progressive`, `smooth` (running mean, edge padding), `target`, `ylim` | cursor, dots (+ lines when progressive) |
| `dual_axis` | 2 series | twin y axes, right spine visible, labels `<label> (left axis)` / `(right axis)` | `best` | - (ignores smooth / target / ylim) | cursor, 2 dots |
| `heatmap` | 1 matrix (rows x frames) | `imshow` nearest, extent `t0-dt/2 .. t-1+dt/2` with `dt = t1-t0`; sequential map (NaN -> grid colour) `vmin` (0) .. `vmax` (max) or ordinal map for class matrices (tick labels = classes); residue ticks thinned to `max_ticks` (16); y tick labels 0.7 x font | colour bar `fraction 0.04, pad 0.02`, label = unit or `count`, integer ticks for counts | `vmin`, `vmax`, `max_ticks` | cursor |
| `contact_map` | 1 pairmap (CSR per frame) | residue x residue `imshow`, `instant` (frame pos), `frame` (fixed, clamped), `frequency` (fraction of frames, vmax 1); x labels rotated 90°; 0.65 x font | colour bar `count` / `fraction of frames` | `mode`, `frame`, `transpose`, `max_ticks` (14) | image (instant) |
| `sequence` | `trace` (+ proximity) | one-letter cells (bold monospace) coloured like the structure (text black / white by luminance 0.55), chain labels italic in gaps, tags under letters, largest cell that fits (shrink 0.92, min 6 px), note line `N residues need R rows: choose a smaller range ...` and one warning above 2 rows | - | `residues`, `per_row` | boxes, letters, halos, note |
| `secondary` | `trace` + `basepairs/map` (or dot-bracket) | residues on a line (chain gap 1.5), arcs (25 points, height ∝ span) above for intra-chain (`ink2` 1.4 pt) and transparent yellow below for inter-chain, letters ≤ 90 residues, tags, chain names | - | `source` (`dynamic`, `first`, `dotbracket`), `dotbracket`, `residues` | arcs when pairs change, points, halos, labels |

### 1.4 Structure view details (`panels.py`, `coloring.py`, `residues.py`, `orientation.py`)

- Orientation: `R = Rz Ry Rx` (degrees) about the reference-frame centroid; same matrix for the
  movie; `drag()` screen trackball 0.5°/px, `matrix_to_euler` with gimbal fold, `wrap_degrees` to
  (-180, 180]; `format_angles` (GUI). Centre: `visual.center_on` (`""` model centre, `origin`,
  residue spec, warning once and model centre when it matches nothing), `center_follow` shifts every
  frame by the residue's motion. Zoom `>= 0.05`.
- Colour modes `chain` (matplotlib palette `visual.chain_palette`, qualitative maps by index, else
  `linspace(0.05, 0.95)`), `groups`, `proximity` (base -> `proximity_color` by intensity, members keep
  their colour). Glow modes `none / ever / fade / current` (`fade_frames`, linear fade), strength,
  members never glow; renderers that draw the glow themselves get untinted colours.
- Tags: `label_every` (resid % n == 0) and `label_residues` ranges (bad spec ignored).
- Density shadow (`density.py: project_density / shadow_rgba`): map above `level` σ integrated along
  the view axis (`ndimage.affine_transform` order 1), ≤ 200 px (100 while dragging), normalised,
  `alpha = img^0.7 * opacity`, bilinear, behind (`zorder 0.5`) or in front of the 3D layer
  (`density.front`, 3.5); follows a followed centre; notes `not fitted yet`, `load the real reference
  frame`. Placement: stored 12 numbers + `anchor` (`anchored()`), else `ensure_density_fit` runs the
  FFT global search + Powell refinement first (`auto_fit`), a failed fit is a warning.
- 3D layer (`molecule/stage.py`): `PixelImage` drawn 1:1 into the axes box (re-rendered at the
  backend's pixel size, also for SVG / PDF magnification); the scene hides the 2D trace for every
  residue; a failed 3D render logs once and falls back to the 2D trace for that frame.

### 1.5 Molecule representations (`molecule/`)

- **Topology** (`topology.py`): backbone control atom per residue (`backbone_atom` P or C4', protein
  CA), segments broken at 8.5 Å, purine / pyrimidine / sugar rings, glycosidic atoms, C1'/C4',
  nucleotide letter from residue name or atoms (CHARMM names), element colours (CPK-like table),
  vdW radii table, nucleotide colours (A / G / C / U / T / ?).
- **Styles** (`geometry.py`): `cartoon` (Catmull-Rom backbone ellipse 1.6 x 0.45, rungs to C1' /
  C4', base slabs, optional `sugars` slabs, balls at lone residues and rung joints), `ribbon` (band
  1.2 x 0.12), `tube` (ellipse 0.45), `licorice` (half-bond sticks r 0.3, joins neighbour styles),
  `vdw` (spheres, radius x size), `surf` (QuickSurf-like Gaussian surface, σ 1.4, iso 0.5, grid
  1 Å or 1.5 Å above 5 000 atoms, surface nets, gradient normals, nearest-atom colour), `trace`
  (balls 0.85 + half-sticks 0.38). `RepSpec{style, groups, color, size, opacity, sugars}`;
  `rep_rest` = `trace` fills undrawn residues (ligands: licorice) or `hide`; `rep_hydrogens`.
- **Colour schemes**: `auto` (carbons keep the residue colour, others by element; bases by
  nucleotide), `residue`, `element`, `nucleotide`; sugar slabs lighter (`0.55 rgb + 0.45`).
- **Renderers** (`visual.renderer`): `auto` / `software` = numpy z-buffer renderer (`software.py`:
  sphere / stick impostor sprites, mesh rasterisation, sphere culling, translucency = front-most
  translucent layer blended once, depth fog, glow tint + Gaussian halo at 1/4 resolution, ink outlines
  at depth jumps and silhouettes, alpha AA, `fast` = half resolution, thread safe); `vtk` = optional
  VTK offscreen with a subprocess probe, falling back to the software renderer with a note.

### 1.6 Charts (`charts.py`, `[output.charts]`)

- Series: `output.charts.series` (missing keys: `WARNING chart series K is not available`; matrix
  keys silently dropped) or every series key except `trace/`, `proximity/`; heatmaps of every matrix
  key (same exclusion) when `include_heatmaps`.
- Series chart: `width x height` in (8 x 4.5) at `dpi` 300, constrained layout, theme of
  `visual.theme` with `style = "mdatk"`, else a Matplotlib style sheet (unknown style: warning, MDATK
  chrome) on the LIGHT theme; line `color` (or palette 0) `line_width`; `smoothing > 1`: raw line at
  0.6 width alpha 0.35 `Per frame` + `Running mean (n frames)`; target line when `show_target`;
  legend only with smoothing / target; labels `xlabel` / `ylabel` or automatic; integer ticks for
  counts; title `{label}` / `{key}` template, left; xlim `t0..t-1`.
- Matrix chart: height `max(height, 0.18 * rows + 1.5)`, as the heatmap panel but `Residue` y label,
  ≤ 40 row ticks at 8 pt, colour bar label unit / `count`.
- Files: `<out_dir>/<prefix>_<safe key>[_<stamp>].<fmt>` for each of `formats` (`png` default; any
  Matplotlib format works: svg, pdf, eps, jpg ...), `facecolor` kept; log
  `Wrote N chart file(s) to D`; progress `Chart <key>` per job.

### 1.7 Video (`video.py`, `[output.video]`)

FFmpeg executable (`ffmpeg_path`, PATH lookup; `VideoError` messages), raw RGBA on stdin, `-r fps`,
codec `libx264` (default; `-crf N -preset medium`), `libx265` (crf), `mpeg4` (`-q:v 3`), any other
name passed through; missing encoder -> first of libx264 / libopenh264 / mpeg4 with a warning;
`-pix_fmt yuv420p -movflags +faststart`, no audio. Positions `0, s, 2s, ...` (`frame_stride`); log
`Rendering N video frame(s) at WxH, F fps -> PATH`, progress every `n/200` frames, cancel between
frames, temporary `.MDATK_tmp_<pid>_<name>` in the target folder renamed at the end (removed on
failure), `Wrote PATH (N frames)`. File name `MDATK_video.mp4` time-stamped
(`MDATK_video_2026-10-15_143205.mp4`) when `output.timestamp_files`.

### 1.8 Commands and the pipeline's render path (`cli.py`, `pipeline.py`)

- `render [-p] [--video|--no-video] [--charts|--no-charts] [--pdb|--no-pdb] [--force] [--from-cache]`
  + common options: full run (inputs, required analyses, video 0.55-0.92, charts 0.92-0.97, PDBs,
  params log, step log entry `render` with files) and prints `summary()` JSON (`video`, `charts`,
  `pdbs`, `cache_manifest`, `params_log`); `--from-cache`: dataset from the manifest, prints
  `{"video": ..., "charts": [...]}` only for the enabled ones.
- `charts [--series a,b] [--formats png,svg] [--dpi N] [--from-cache]`: prints `{"charts": [...]}`.
- `preview [--frame N] [--output PATH] [--from-cache]`: one composite frame (video size and dpi) as
  PNG / SVG / PDF; default `<out_dir>/MDATK_preview_<pos:05d>.png`; prints the path.
- Density: video and preview run `ensure_density_fit` first. The params log's `environment` lists
  matplotlib and the `ffmpeg -version` line.

### 1.9 Interactive (GTK GUI) vs. batch

| GUI only (`gui/canvas.py`, `app.py`) | Batch (CLI / pipeline) |
|---|---|
| WYSIWYG canvas: the real panels drawn on the example dataset (`example.py`: 120 frames, seeded smooth series, schematic duplex or the loaded reference frame), at window size with the video's aspect and type scale | video, charts, preview from real data |
| drag = rotate (writes `visual.rotation`), scroll = zoom, map mode drag / Shift-drag moves / rotates the density map, drag panel onto slot = swap, double-click = edit panel; `fast` quality while dragging, settle after 250 ms | |
| "Load real frame" (`structure_preview`), "Preview frame", "Open video" (system player) | |

M4 ports what the canvas needs from `render/` (example dataset, rotation, fast quality, density
nudge API) so M5 only builds the widgets.

### 1.10 MDATK tests that specify rendering (each gets a C++ or CLI counterpart)

`test_render_cache.py` (glow modes, static colours, Euler round trip / drag, slots tile without
overlap), `test_seqpanels.py` (12: one-letter codes, residue rows, strip slots, wrapping and warning,
pairs change with frames, dot-bracket, basepairs requirement, zoom / centre, followed centre, long
chain names, animated labels, chain label never ends a row), `test_molecule_geometry.py` (11:
names, Catmull-Rom, sweep faces, ring slabs, roles, outward meshes per style, banded spline, CHARMM,
fall-backs, protein tube, rung balls), `test_molecule_render.py` (26: sphere area / depth exact,
stick width, no pinholes, order independence, depth range, handedness, fog, halo, translucency, fast
size, threads = serial, one-atom surface, cartoon independent of selection, rest as trace, headless,
video with representations, culling invariance, many items, translucent surfaces, sliver hidden,
non-finite atoms, vector output placement, straight colour boundary, glow helpers, fallback to 2D,
front shadow animated), `test_density.py` (shadow under the model, fit + video, followed shadow,
degraded placements), `test_pipeline.py` / `test_cli.py` / `test_v002.py` (full run: 24 frames at
320 x 180, 4 charts, stamped names; rerender from cache; `render`, `charts --formats svg`, `preview`;
example dataset matches layout).

---------------------------------------------------------------------------------------------------

## 2. External dependencies: MDATK -> MDBridge (licence-checked against LICENSING.md)

| Need | MDATK uses | MDBridge uses | Licence / rule |
|---|---|---|---|
| Figure, panels, charts | Matplotlib 3.11 (Agg) | own figure model + QPainter (`QImage` raster, `QSvgGenerator`, `QPdfWriter`) | Qt LGPL v3 (dynamic); QtSvg is one more Qt module (`get_qt.py --archive qtsvg`, AppImage / windeployqt). Qt Charts / Qt Graphs (GPL v3 in open-source Qt 6), QCustomPlot (GPL) excluded; JKQtPlotter (LGPL 2.1) allowed but not needed (we need Matplotlib-identical ticks, layout and colour maps) |
| Tick locator / formatter, autoscale margins, colormap interpolation, colour names | Matplotlib (`MaxNLocator`, `ScalarFormatter`, `to_rgba`, `LinearSegmentedColormap`, `colormaps`) | ported algorithms + tables with notice (`src/render/mpl_*`) | Matplotlib licence (PSF-style, permissive): port with notice; CSS4 names W3C; `tab10/tab20` (Tableau, distributed in Matplotlib's BSD-compatible terms), ColorBrewer (Apache-style), viridis family (CC0) |
| Fonts | DejaVu Sans / Mono shipped with Matplotlib | the same DejaVu fonts bundled as Qt resources, registered with `QFontDatabase::addApplicationFont`; no system fonts used for figures | Bitstream Vera / DejaVu licence (permissive), notice; identical metrics on every OS = reproducible layout, works on the bare Rocky system |
| Headless text / painting | Agg | `QT_QPA_PLATFORM=offscreen` set by the CLI for `render`/`charts`/`preview` when no display is needed (plugin already deployed in the AppImage) | Qt |
| 3D: software path | own numpy renderer (`software.py`) | C++ port (owner's code), also the fallback renderer of the GUI | - |
| 3D: GPU path | VTK (optional, BSD) | Vulkan 1.3: impostor spheres / cylinders, meshes, fog, outlines, glow, OIT-lite (front-most translucent layer); offscreen `VkImage` + readback for batch (no surface; works without X on lavapipe / GPUs) | Vulkan-Headers Apache/MIT; system loader; shaders via glslang (build tool only, header committed as now) |
| scipy pieces | `ndimage.gaussian_filter` / `uniform_filter` / `map_coordinates` / `affine_transform`, `cKDTree.query`, `fft`, `optimize.minimize(Powell)` | ported 1-D separable filters, order-1 interpolation, affine resampling, nearest-atom via own cell grid (exists in `analysis/spatial`); FFT / Powell only if the density fit is in scope (D9: pocketfft BSD port or own radix-2/3/5) | scipy BSD-3: port with notice |
| Video | `ffmpeg` executable (often GPL builds with x264) | FFmpeg libraries (`libavcodec`, `libavformat`, `libavutil`, `libswscale`), LGPL-only build (`--disable-gpl --disable-nonfree`), linked dynamically, shipped in the AppImage / installer with licence + source pointer; encoders per D1 | LGPL v2.1+ (dynamic link, relinkable, source offered). OpenH264 BSD-2 code but H.264 patent licence only covers Cisco's own binaries downloaded at run time; libaom BSD-2 / SVT-AV1 BSD-3-Clause-Clear + AOM patent licence (royalty-free); libvpx BSD-3; FFmpeg's native `mpeg4` (MPEG-4 Part 2) LGPL. Qt Multimedia (`QMediaRecorder` + `QVideoFrameInput`, Qt 6.8+) is an alternative using Qt's own LGPL FFmpeg, but gives no control of encoder, CRF, faststart or frame-exact timing: not recommended |
| Molecular surfaces | own Gaussian surface (numpy / scipy) | port of MDATK's QuickSurf-like surface nets (CPU, threaded per frame; GPU compute later). No SES / SAS library needed for parity: MSMS (proprietary), EDTSurf / NanoShaper (GPL) are excluded anyway; an own SES (e.g. distance-field based) is post-M4 | - |
| Density fit | numpy FFT + scipy Powell | see D9 | - |

Build and fetch: FFmpeg + encoder are built by a pinned, hash-checked script
(`packaging/ffmpeg/build_ffmpeg.sh` on Rocky 9 → glibc ≤ 2.34; Windows per D2) into a prefix that
CMake finds (`cmake/FindFFmpegLGPL.cmake`, checks `avcodec_configuration()` has no `--enable-gpl`
at configure time and in a unit test), cached in CI by script hash. A source build without the
libraries still builds; video then fails with `video encoding is not available in this build`.

---------------------------------------------------------------------------------------------------

## 3. Parity strategy

Pixels cannot match Matplotlib's Agg byte for byte. "Same output" is defined per artifact:

| Artifact | Same means | Tolerance |
|---|---|---|
| CLI stdout JSON, file names and folders, logs, warnings, progress messages, exit codes, step log and params-log `outputs` | byte for byte (time stamps fixed) | exact; `environment` block is MDBridge's own (documented difference) |
| Figure semantics: per panel the data shown (line x/y, scatter offsets, sizes, RGBA, segment colours / widths, image arrays and colour-map normalisation, texts and their strings, legend entries in order, colour-bar ticks / labels, axis limits, tick locations and tick labels, titles, labels, visibility, z-order classes) after `build_figure` and after `update(pos)` for chosen frames | equal | numbers 1e-9 relative (float64 math) / exact for strings and counts; RGBA ±1/255 |
| Layout geometry: each axes rectangle (and colour bar, twin axis) in figure pixels | MDATK's constrained-layout rectangles | ≤ 1.5 % of the figure size per edge (fonts and paddings identical; the solver is MDBridge's). Ticks are compared exactly whenever MDATK's auto `nbins` equals MDBridge's for that axis (both recorded); else only "nice and covering" |
| 3D layer image (software renderer) | MDATK's numpy `SoftwareRenderer.render` RGBA for the same items and viewport | ≥ 99.5 % of pixels equal, all within 2/255, coverage mask IoU ≥ 0.999 (rounding of float32 vs numpy order) |
| 3D layer image (Vulkan) | the software renderer's image | coverage IoU ≥ 0.99, SSIM ≥ 0.97 on lavapipe; checked in CI, documented as "same look, not same pixels" |
| Scene geometry (meshes, sticks, spheres per style, surface vertices / triangles, nearest atoms) | MDATK's `MoleculeScene.frame` DrawLists | vertices 1e-5 Å, indices exact |
| Composite PNG / SVG / PDF (preview, charts, video frames) | MDBridge's own golden images, rendered once, reviewed by the owner against MDATK's renders side by side (`docs/specs/m4-gallery.md`), then committed | PNG: ≥ 99 % pixels within 8/255 and mean abs ≤ 1/255 (Linux / Windows, bundled fonts); SVG / PDF: structure checks (page size, text runs, embedded 3D raster position) |
| Video | frame count, width / height, `r_frame_rate`, time base / duration, codec family, `yuv420p`, `moov` before `mdat` (faststart), no audio | exact; decoded frames vs MDBridge's composed RGBA: PSNR ≥ 38 dB (encoders are not bit-reproducible across thread counts, so no byte comparison) |

Generators (`tests/parity/`, run by `regenerate_all.py`, MDATK in-process with
`independent_tables.py` wherever trace / atoms / proximity / basepairs analyses run, as in
`tests/data/README.md`):

- `make_figure_cases.py` → `tests/data/render/figures/<case>.json`: walks MDATK's figure after
  `build_figure` + `update(pos)` (positions 0, mid, last) and dumps the semantics above (axes
  `get_position()` x figure size, `get_xlim/ylim`, `get_xticks` + tick label texts, `Line2D` data,
  `PathCollection` offsets / sizes / facecolors, `LineCollection` segments, `AxesImage` arrays as
  npz, texts, legend texts and handle colours, colour-bar ticks), plus `animated_artists` classes.
  Cases (≈ 40): every panel type x options, every layout preset, light / dark, title / watermark,
  1 frame, odd sizes, missing data, the demo system (`demo_dry`), MDATK's micro systems, a protein +
  RNA system (protein tube), a 3-chain system (sequence wrapping, inter-chain arcs), a density map
  case (`tests/data/io/maps`), proximity / glow modes, followed centre.
- `make_chart_cases.py` → charts' semantic dumps (same walker) + the list of paths `render_charts`
  returns.
- `make_molecule_cases.py` → scene DrawLists (npz) and software-renderer RGBA (npz) for
  `test_molecule_render`'s synthetic items and for each style on `demo_dry` (fixed viewports);
  topology roles and letters (json).
- `make_render_cli_cases.py` → MDATK's `render` / `charts` / `preview` stdout, file lists, logs and
  ffprobe metadata on the demo (needs `ffmpeg`/`ffprobe` in the regenerate job: apt install, pinned
  Ubuntu package; MDATK's own video bytes are not committed).
- `improvements_render.py`: test-side patches of MDATK that turn the dumps into MDBridge's improved
  behaviour (RN-xx below; e.g. per-frame heatmap edges, ion legend text), in the style of
  `improvements_*.py`; every patch named in DIFFERENCES.md.
- `cli_parity.py` gains `render --from-cache --no-video`, `charts`, `preview`, and a video run
  compared by metadata.

Reproducibility: Matplotlib 3.11.2 is already pinned in `requirements-regenerate.txt`; its bundled
DejaVu fonts make MDATK's dumps independent of the machine; dumps store floats with `repr`. MDATK's
video needs no x264 for metadata (`mpeg4` fallback is fine; codec family is MDBridge's own anyway).

Headless CI: Linux job runs the render tests with `QT_QPA_PLATFORM=offscreen` and lavapipe
(`mesa-vulkan-drivers` already installed) for the Vulkan tests; `appimage-rocky` (no X, no Vulkan)
runs `charts`, `preview` and a 24-frame `render` from the AppImage (software path, bundled fonts,
bundled FFmpeg) — this is the batch-on-a-cluster proof; `appimage-ubuntu` runs `preview` with
`--renderer vulkan` on lavapipe. Windows: render tests on the software path; Vulkan tests with
Mesa's lavapipe for Windows (`vulkan_lvp.dll`, MIT, test-only via `VK_DRIVER_FILES`, never
shipped) — WARP is Direct3D, not Vulkan (Mesa's Dozen on WARP is an option if lavapipe fails).
Golden-image tests run on both platforms.

---------------------------------------------------------------------------------------------------

## 4. MDATK rendering bugs and fragile behaviour to fix (CLAUDE.md rule; ids RN-xx)

Each lands in DIFFERENCES.md with a test; `improvements_render.py` adjusts the expected data.

| Id | MDATK behaviour | MDBridge |
|---|---|---|
| RN-01 | Heatmap panel / matrix chart cells use `dt = t[1]-t[0]` for every column: with non-uniform times (concatenated trajectories with different `dt`, time resets) cells are misplaced and the cursor misses them | column edges at the midpoints between frame times (first / last half-step outward) |
| RN-02 | No analysed frames (empty selection after `--from-cache` of an empty run) → `IndexError` at `t[0]` in every panel, crash with a traceback; `preview` computes `pos = -1` | panels show `No data: no frames`; `preview` / `render` stop with `no analysed frames to render` (exit 1) |
| RN-03 | Non-finite coordinates (blown-up frame) make the structure view's radius NaN → `Axis limits cannot be NaN or Inf` crash; non-finite series values crash y autoscaling | finite-only extents (the 3D layer already skips them); `WARNING` once naming the frames |
| RN-04 | `vmin > vmax` (heatmap options or a negative-valued matrix with default `vmin 0`) → Matplotlib `ValueError` at draw time, video aborted | clamp `vmax = max(vmax, vmin + 1e-9)` and warn once; `check()` validates numeric panel options and reports `layout.panels[k].options.vmin` |
| RN-05 | Matrix chart height `0.18 in` per row: thousands of residues → images beyond Agg's 2^16 px limit, crash | height capped at 2^15 px (rows thinned, already ≤ 40 labels); warning with the cap |
| RN-06 | Video size: `figsize = w/dpi` rounding can make the canvas differ from the request → `canvas size not even` / "choose a DPI that divides the resolution" errors | the canvas is exactly `figure_size(w, h)` pixels for any DPI (no error; message kept for impossible sizes < 2 px) |
| RN-07 | `codec = "libx264"` (the default) falls back silently to whatever exists; CRF ignored for `libopenh264` and `mpeg4` | codec names map to the bundled encoders (`libx264`, `libx265`, `h264`, `hevc` → the H.264 / AV1 encoder of D1; `mpeg4`, `libaom-av1`, `libsvtav1`, `libopenh264` as named); one line `video: codec libx264 is not bundled (GPL); using <x>`; `crf` mapped to the encoder's quality scale (documented table) |
| RN-08 | Structure legend always says `Mg²⁺`, whatever `analyses.mg.mg_selection` / `show_ions` selects (Na+, K+, Ca2+) | label from the ions' element(s): `Mg²⁺`, `K⁺`, `Na⁺`, `Ions` for mixtures |
| RN-09 | Chain legend `Chain 1, 2, ...` from the chain index although chain letters are known (sequence panel shows `A, B`) — misleading next to `B:12` ranges | `Chain A` (the trace's chain names; index only when a chain has no name) — ASK (D11) |
| RN-10 | `charts.series` naming a matrix key is dropped silently | charted as a heatmap (or warning `chart series K is a matrix; enable include_heatmaps`) |
| RN-11 | Unknown chart / preview format (`jpg`, `eps`, typo) → Matplotlib `ValueError` traceback after the analyses ran | formats validated in `check()`: png, svg, pdf (+ jpg, tif, webp via QImageWriter); others listed in the error |
| RN-12 | Invalid colour strings (group colours, `glow_color`, chart `color`) raise in the middle of a video | `check()` validates every colour with Matplotlib's grammar and names the field |
| RN-13 | `chain_palette` unknown → `ValueError` from Matplotlib's registry mid-render | validated in `check()`; error lists the bundled palettes |
| RN-14 | `_running_mean` with NaN (an empty RMSD target, NaN frames) blanks whole windows | NaN-aware window mean (NaN only where the whole window is NaN) |
| RN-15 | Panel options of the wrong type (`ylim = 5`, `frame = "x"`) → `TypeError` traceback | `check()` reports them (`layout.panels[k].options.ylim: expected two numbers`) |
| RN-16 | The VTK probe / fallback note is per process; renderer `vtk` on machines without EGL can crash inside VTK | `vtk` means "GPU": Vulkan when usable, else the software renderer with the note (no VTK) |
| RN-17 | Video positions with a stride skip the last frame (`0, s, 2s, ...`) | kept (convention: the frame count users expect); documented, not changed |
| RN-18 | `theme.ordinal(n)` for more than 5 classes keeps the 4 ramp colours: classes share colours under the colour bar's distinct labels | n - 1 evenly spaced ramp steps (`from_list(ramp, N = n - 1)`); MDATK's colours up to 5 classes (W0-B) |
| RN-19 | A `surf` grid for a blown-up frame (atoms thousands of Å apart) is allocated: `MemoryError` or gigabytes and minutes per frame | grids above 2^26 points refused up front; the surface is left out of that frame with the existing warning (W0-F) |
| RN-20 | `visual.renderer = "vulkan"` is rejected by `validate()` (MDATK knows auto / software / vtk) | accepted (D7); the unknown-renderer message lists it (W0-G) |
| RN-21 | Example dataset of a structure with one strand (one structure group / chain): `_example_basepairs` slices `a[:m - 2]` with `m = 0` → `IndexError`, the GUI's layout preview breaks | no example pairs for a single strand (W0-H) |
| RN-22 | `wrap_degrees` uses `np.isclose(a, -180)`: every angle within 0.0018° of -180 becomes +180 | only -180 itself (within 1e-9) folds to +180 (W0-H) |
| RN-23 | Unreadable `visual.label_residues` silently gives no tags | same tags, plus one warning naming the field and the parse error (W0-H) |
| RN-24 | `visual.label_residues = "B:12"` ignores the chain (residue 12 of every chain tagged) | the chain is honoured, as `residue_rows` does (W0-H) |
| RN-25 | Example data: `group_resids` expands whole ranges before keeping 400 (gigabytes for `1-1000000000`); the schematic sequence has 120 letters (`IndexError` for larger groups) | stop at 400 residues; the sequence repeats (W0-H) |
| RN-26 | Software renderer: the edge anti-aliasing filters the molecule's box alone (zeros outside): a darker rectangle in every glow halo | the 3 x 3 mean sees the real neighbours (W1-B) |
| RN-27 | Software renderer `fast` at odd sizes: half size rounded down (half to even), a transparent last column / row, the view scaled by up to 1 px | `ceil(size / 2)` with the limits widened, cropped (W1-B) |
| RN-28 | Software renderer: equal x / y limits → `ZeroDivisionError`; reversed or NaN limits / rotation, coordinates beyond float32 → garbage; a sphere larger than the image → `MemoryError` | `invalid_argument("3D view ...")`; unplaceable primitives skipped; sprites clipped to the image (W1-B) |
| RN-29 | Software renderer: the glow halo's normalisation cached under the size rounded to 0.01 px (depends on earlier renders) | computed for the exact size (W1-B) |
| RN-30 | Density shadow: one NaN voxel makes the whole shadow NaN (`np.maximum` keeps NaN, the peak is NaN); an infinite voxel makes it NaN or blank | non-finite voxels carry no density; a non-finite level gives no shadow (W0-I) |
| RN-31 | Non-finite placement, view, centre or map grid (NaN voxel size) → `cannot convert float NaN to integer` / `OverflowError` from `int(ceil(...))` | `ValueError("density map placement, view or grid is not finite")` (W0-I) |
| RN-32 | `anchored` with a ragged anchor (rows not all x, y, z) or non-finite anchor / reference points → numpy's `inhomogeneous shape` `ValueError` / `LinAlgError: SVD did not converge`, the map not drawn | the stored placement used as it is, with a note saying why (for a warning) (W0-I) |
| RN-40 | Sequence panel `per_row` that `int()` cannot read → `ValueError` out of `build_figure` (render aborted); a negative one → inverted axis, giant boxes | `No data: option 'per_row' must be a whole number ...`; 0 / truncating values as MDATK (W2-C) |
| RN-41 | Secondary panel: any unknown `source` silently means `dynamic` | `No data: unknown source 'x' (option 'source': dynamic, first or dotbracket)` (W2-C) |
| RN-42 | Secondary panel's "No data" messages after `_bare_axes` have no frame (other "No data" panels do) | the frame kept in both residue panels (W2-C) |

Not bugs, kept as conventions: DPI semantics (text size relative to the frame), `MDATK_` file names,
the example dataset's seeds (crc32 of the key → numpy PCG64: ported exactly so the GUI preview looks
the same), 120 example frames, the glow fade counted in analysed frames.

---------------------------------------------------------------------------------------------------

## 5. Work breakdown

Layout of the code (all in `src/render`, namespace `mdbridge::render`; Qt Gui, no widgets):
`color.*` (parsing, palettes, colormaps), `theme.*`, `ticks.*` (locators, formatters),
`text.*` (fonts, metrics), `figure/` (figure model, layout solver, display list, QPainter backends),
`panels/` (7 panels + registry), `charts.*`, `molecule/` (topology, geometry, surface, scene —
no Qt), `software/` (CPU renderer — no Qt), `vulkan/` (GPU scene renderer, offscreen + window),
`density_shadow.*`, `example.*`, `video/` (frame loop, encoder), `render_pipeline.*` (video /
charts / preview entry points used by `analysis::Pipeline` callers). Commands in
`src/app/commands/render.cpp`.

### Wave 0 (independent; start together)

| Chunk | Content (MDATK source) | Files | Acceptance |
|---|---|---|---|
| **W0-A Generators** | `make_figure_cases.py`, `make_chart_cases.py`, `make_molecule_cases.py`, `make_render_cli_cases.py`, the artist walker, `improvements_render.py` skeleton; regenerate job installs ffmpeg | `tests/parity/*`, `tests/data/render/**`, `tests/data/README.md` | `regenerate_all.py --only` reproduces byte for byte twice |
| **W0-B Colours** | `to_rgba` grammar (hex 3/4/6/8, CSS4 / base / `tab:` / `xkcd:` names, `C0`-`C9`, grey floats), colormaps (`from_list`, listed, `with_extremes`, bad colour), chain palettes, `theme.py` | `color.*`, `theme.*`, generated `mpl_colors.inc` (tool `tools/gen_mpl_tables.py` from Matplotlib's data, notice) | equality with a generated `to_rgba` / colormap sample table (≥ 2 000 cases) |
| **W0-C Ticks & numbers** | `MaxNLocator` (auto nbins from axis length / font, steps, integer), `AutoLocator`, `ScalarFormatter` (offset, sci, unicode minus), colour-bar locator, `thin_ticks`, autoscale margins 5 % | `ticks.*` (Matplotlib notice) | generated locator / formatter cases exact |
| **W0-D Fonts & headless text** | bundle DejaVu Sans / Sans Bold / Mono Bold; registry; pt→px; QPA offscreen selection for batch commands; licence files | `text.*`, `data/fonts/`, `packaging/licenses/`, AppImage / installer, `get_qt.py --archive qtsvg` in CI | bare-Rocky smoke: a text image renders identical to the Linux job's |
| **W0-E FFmpeg** (owner decisions D1, D2 first) | build script, CMake find module, GPL check, encoder wrapper (RGBA→yuv420p via swscale, codec map RN-07, CRF map, mp4 faststart, temp + rename, error messages), demux-based self-check for tests | `packaging/ffmpeg/*`, `cmake/FindFFmpegLGPL.cmake`, `src/render/video/encoder.*`, CI cache steps | encodes 24 frames 320 x 180 on Linux and Windows; metadata test; `avcodec_configuration()` has no GPL flag; AppImage glibc check passes |
| **W0-F Molecule geometry** | `topology.py`, `geometry.py`, `surface.py`, `scene.py` (DrawLists, styles, colour schemes, fill_rest), Gaussian filter port | `molecule/*` (no Qt) | `test_molecule_geometry` counterparts; DrawLists vs `make_molecule_cases` within 1e-5 Å |
| **W0-G Vulkan offscreen base** | device / queue without a window, render target + readback, pipeline cache, `compile_shaders.py` with more shaders, test harness that skips cleanly without Vulkan | `vulkan/context.*`, `vulkan/offscreen.*`, `tools/compile_shaders.py`, `tests/render/` | clear-colour image read back on lavapipe (Linux and Windows CI) |
| **W0-H Shading logic** | `coloring.py`, `residues.py` (one-letter codes via the existing `parseResidueName` of src/io/pdb.hpp, rows, tags), `orientation.py`, `example.py` (crc32 seeds + numpy PCG64 `normal` / `uniform` / `integers` ported bit-exact, numpy BSD notice) | `shading.*`, `orientation.*`, `example.*` (no Qt) | `test_render_cache` / `test_v002` example counterparts; example arrays bit-equal |
| **W0-I Density shadow** | `project_density`, `_above_level`, `shadow_rgba`, `anchored` use (core has RigidTransform), `affine_transform` order 1 port | `density_shadow.*` | shadow arrays vs MDATK within 1e-6 |

### Wave 1

| Chunk | Needs | Content | Acceptance |
|---|---|---|---|
| **W1-A Figure framework** | B, C, D | figure model (figure → axes, twin axes, colour bars, legends, texts, images, collections), constrained-style layout solver (GridSpec + height ratios + MDATK's pads, tick / label / title / colour-bar space from font metrics), freeze, display list, QPainter backends (QImage RGBA, SVG via QSvgGenerator with text kept as text, PDF via QPdfWriter, exact pixel size), static / dynamic layer split for the frame loop, `No data` panel | geometry within § 3 tolerance on every preset; SVG / PDF page sizes; semantic dump writer used by all later tests |
| **W1-B Software renderer** | F | port of `software.py` (z-pass, sprites, mesh raster, culling, translucency, fog, glow halo, outlines, AA, fast mode), threads | `test_molecule_render` counterparts; RGBA vs MDATK per § 3 |

### Wave 2

| Chunk | Needs | Content | Acceptance |
|---|---|---|---|
| **W2-A Plot panels** | W1-A, B, C | `timeseries`, `dual_axis`, `heatmap`, `contact_map` (+ RN-01, RN-04, RN-14) | figure dumps equal for their cases |
| **W2-B Structure panel** | W1-A, W1-B, H, I | 2D trace, halos, tags, legend (RN-08, RN-09), time label, centre / zoom / aspect, density layer (behind / front, follow), 3D stage with fallback to 2D, `set_rotation`, `nudge_density`, fast quality | dumps + `test_seqpanels` structure cases + density tests |
| **W2-C Sequence & secondary** | W1-A, H | both panels incl. re-flow after layout and the warning | `test_seqpanels` counterparts; dumps |
| **W2-D Charts** | W1-A (shares W2-A's heatmap code; can start on its own matrix path) | series / matrix charts, styles (D5), formats, names, logs (RN-05, RN-10, RN-11) | chart dumps; file lists exact |
| **W2-E Vulkan scene renderer** | F, G, W1-B (reference) | impostor spheres / cylinders, meshes, fog, outline pass, glow halo, translucency, AA; offscreen for batch, `VulkanWindow` path for the GUI | § 3 Vulkan tolerance against W1-B on every molecule case; frames/s benchmark on lavapipe and a real GPU |

### Wave 3

| Chunk | Needs | Content | Acceptance |
|---|---|---|---|
| **W3-A Video loop** | W2-A..C, W0-E | static background + dynamic layers per frame, positions / stride, progress / cancel / logs as MDATK, frames composed in parallel and encoded in order, temp + rename, renderer choice (auto: Vulkan offscreen if usable, else software) | `test_full_run` counterpart (24 frames 320 x 180); metadata exact; PSNR check; output identical for 1 and N threads before encoding |
| **W3-B Pipeline + CLI** | W3-A, W2-D | `render` (incl. PDBs, params log `outputs`, step log), `charts`, `preview` with MDATK's options / JSON; `check()` validations (RN-11..15); `required_analyses` already in M3; density placement missing → warning per D9 | `cli_parity.py` cases green; `test_cli` / `test_pipeline` counterparts |
| **W3-C GUI** | W2-B, W2-E | the 3D view shows the project's reference-frame structure (rotate / zoom with the mouse, the project's representations, software fallback), a read-only figure preview (example dataset, watermark) in the main window | smoke test screenshots (Vulkan on lavapipe, software on bare Rocky) |
| **W3-D Density fit** (if D9 = in M4) | W0-I | FFT global search + Powell refinement port, `density fit`, `auto_fit` | `test_density` fit cases |

### Wave 4 (after parity is green)

Docs (DIFFERENCES entries RN-xx, LICENSING rows: FFmpeg + encoder "since milestone 4", DejaVu,
Matplotlib / scipy ports, QtSvg; THIRD_PARTY_NOTICES; ARCHITECTURE: replace "e.g. QCustomPlot" by
"own QPainter charts (Qt Charts / QCustomPlot are GPL)"; README status), owner gallery review
(`m4-gallery.md`), performance pass (target: ≥ 10 x MDATK's frames/s for 1080p with cartoon on the
demo), smoke tests for video on the bare Rocky system and Windows.

### Critical path and parallelism

`W0-B/C/D → W1-A → W2-B → W3-A → W3-B` (figure framework is the bottleneck; start W1-A as soon as
B / C stubs exist). Second chain: `W0-F → W1-B → W2-B`. Long lead in parallel: `W0-E` (FFmpeg on
two platforms; blocked on D1 / D2) and `W0-G → W2-E` (Vulkan look). 9 agents in wave 0, 2 in wave 1
(the others continue tests / reviews), 5 in wave 2, 3-4 in wave 3.

Chunks needing owner decisions before they start: W0-E (D1, D2), W2-D (D5), W0-B (D6), W2-B (D11),
W3-B (D8, D9, D12), W3-C (D13), the golden-image review (D10).

---------------------------------------------------------------------------------------------------

## 6. Owner decisions (one line each, with recommendation)

- **D1 Video encoder**: H.264 plays everywhere (PowerPoint) but self-built OpenH264 carries patent
  royalties → recommend AV1 (SVT-AV1, royalty-free) in MP4 as the default for `libx264` requests,
  `mpeg4` as the universal fallback, H.264 only via the OS (Windows Media Foundation `h264_mf`) or
  Cisco's run-time OpenH264 download, after legal review.
- **D2 FFmpeg source**: recommend own pinned minimal LGPL build by script (Linux on Rocky 9;
  Windows with MSVC via vcpkg's `ffmpeg` port, which brings its own build shell) over third-party
  prebuilt LGPL binaries (simpler, but provenance and source-offer duties).
- **D3 Chart engine**: recommend own QPainter figure code (Matplotlib-identical ticks / layout,
  vector output) — not JKQtPlotter (LGPL) and never Qt Charts / Qt Graphs / QCustomPlot (GPL).
- **D4 Fonts**: recommend bundling DejaVu Sans / Mono (Matplotlib's fonts) for all figures; user
  fonts later as an option.
- **D5 `output.charts.style`**: recommend `mdatk` plus Matplotlib's `default` look; other style
  names warn and use MDATK chrome (as MDATK does for unknown styles).
- **D6 `visual.chain_palette`**: recommend porting Matplotlib's qualitative (tab10/20, Set1-3,
  Pastel1/2, Dark2, Accent, Paired) and perceptual (viridis, plasma, inferno, magma, cividis) maps;
  others error in `check()`.
- **D7 `visual.renderer`**: recommend keeping `auto` / `software` / `vtk` (vtk = GPU = Vulkan) and
  accepting `vulkan` as a new value.
- **D8 `output.video.ffmpeg_path`**: recommend ignoring it (bundled libraries) and keeping it in
  files; no external ffmpeg fallback.
- **D9 Density fit in M4**: recommend the shadow in M4 and the fit (`density fit`, `auto_fit`) as
  W3-D if capacity allows, else milestone 6 with the warning `density map not drawn: no placement
  (run density fit)`.
- **D10 Image parity bar**: recommend semantic dumps vs MDATK + MDBridge golden images approved in
  a side-by-side gallery, not pixel comparisons with Matplotlib.
- **D11 Chain legend (RN-09)**: recommend chain letters (`Chain A`).
- **D12 Defaults written by `init`**: recommend keeping `codec = "libx264"` and `crf = 20` in new
  project files (byte parity) and mapping at encode time (RN-07).
- **D13 GUI scope**: recommend M4 = real structure in the 3D view + read-only figure preview; the
  WYSIWYG editor (swap, double-click, map drag) in M5.
- **D14 Windows Vulkan CI**: recommend Mesa lavapipe for Windows as a test-only driver.

**Decided (2026-11-25).** Owner: D1 AV1 (SVT-AV1) in MP4 as the default, `mpeg4` fallback; D2 own
pinned LGPL build; D10 **contents only**: figure contents (data, colours, ticks, labels, legends)
compared with MDATK and layout within tolerance, no golden reference images; D13 3D view + read-only
preview, the editor in M5. The other decisions follow the recommendations above (D3-D9, D11, D12,
D14).

---------------------------------------------------------------------------------------------------

## 7. Risks

1. **Video codec / patents** (D1): AV1 in MP4 needs Windows' free AV1 extension for the built-in
   player and PowerPoint; users may see "cannot play". Mitigation: `mpeg4` fallback, clear message,
   legal review of H.264 paths.
2. **FFmpeg on Windows with MSVC** (no MSYS2 rule) and in the AppImage (glibc 2.34, size, encoder
   speed: libaom is slow, SVT-AV1 needs AVX2 for speed). Mitigation: start W0-E first; cache builds.
3. **Layout parity with Matplotlib's constrained layout**: an own solver will drift by pixels; tick
   counts can flip at `nbins` boundaries. Mitigation: identical fonts, tolerance in § 3, ticks
   compared when `nbins` agrees, owner gallery for the rest.
4. **Software renderer port fidelity**: numpy idioms (`np.maximum.at`, stable sorts on packed
   depth/id keys, float32 rounding) decide pixels; allow last-bit tolerance, port tests first.
5. **Vulkan look vs software look** (outlines, AA, glow) and lavapipe speed in CI; Vulkan and
   software must stay visually interchangeable because the fallback is silent.
6. **Headless reality**: Qt text needs a QGuiApplication (offscreen QPA) even in the CLI; fontconfig
   absent on bare systems; verify in `appimage-rocky` early (W0-D).
7. **Scope creep from the GUI** (canvas editor, map dragging) pulling M5 into M4 — held by D13.
8. **Density fit** (FFT + Powell) is a large numeric port; deferring it (D9) leaves `auto_fit`
   projects without a shadow until it lands.
