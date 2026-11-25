# Licensing plan (not yet in effect)

MDBridge is all rights reserved until its licence is chosen (see `LICENSE`). MDATK, the Python
predecessor, stays GPL v3; its own code (same copyright holder) may be ported into MDBridge
freely, but no third-party GPL code may be copied or bundled.

## Intended licensing

Not adopted yet; it takes effect only after legal review and a formal change of `LICENSE`.

- **Non-commercial use:** free, under a standard source-available non-commercial licence
  (candidate: PolyForm Noncommercial 1.0.0) with a notice prohibiting biological-weapons use.
- **Commercial use:** a paid commercial licence (draft held by the author; to be reviewed by a
  lawyer). Licensor to be the author's LLC once formed, with the copyright assigned to it.

## Rules until then

1. **Copyright:** all code must be owned by the author. Outside contributions require a
   contributor licence agreement that permits relicensing.
2. **Resources:** no third-party or institutional resources (equipment, networks, funded time)
   go into MDBridge.
3. **Dependencies bundled in binaries must allow proprietary use** (permissive or LGPL):

   | Dependency | Licence | Rule | Used since |
   |---|---|---|---|
   | Qt 6 (Qt Svg and the image-format plugins since milestone 4: SVG and JPEG figures) | LGPL v3 | link dynamically; ship the licence and a source pointer | milestone 1 |
   | libjpeg-turbo 3.2.0 (inside Qt's JPEG plugin, as Qt ships it) | IJG licence + BSD-3-Clause | OK (licence in `packaging/licenses/libjpeg-turbo-LICENSE.txt`; the IJG acknowledgement in `THIRD_PARTY_NOTICES.md`) | milestone 4 (JPEG figures) |
   | toml++ | MIT | OK | milestone 1 |
   | CLI11 | BSD-3-Clause | OK | milestone 1 |
   | Vulkan-Headers | Apache 2.0 / MIT | OK (headers only; the Vulkan loader is the user's) | milestone 1 |
   | zlib | zlib licence (permissive) | OK | milestone 2 (density maps) |
   | chemfiles | BSD-3 | not used: dropped in milestone 2 for MDBridge's own readers (its XTC decompression is ported, below) | - |
   | Eigen | MPL 2.0 | OK | planned |
   | nanoflann | BSD-2 | OK | planned |
   | FFmpeg (libavcodec, libavformat, libavutil, libswscale) | LGPL v2.1+ | own LGPL-only build (`packaging/ffmpeg`: `--disable-gpl --disable-nonfree --disable-version3`, no x264 / x265 or other GPL parts), linked dynamically, shipped with its licence and the corresponding source on each release; the GPL / non-free check runs at configure time (`cmake/FindFFmpegLGPL.cmake`) and at run time | milestone 4 (video) |
   | SVT-AV1 | BSD 3-Clause Clear + AOM Patent License 1.0 | OK (the AV1 encoder; shared library) | milestone 4 (video) |
   | DejaVu fonts 2.35 (as Matplotlib ships them) | Bitstream Vera / Arev font licences (permissive; modified fonts must be renamed), DejaVu changes public domain | OK (unmodified, compiled in; licence in `packaging/licenses/DejaVu-LICENSE.txt`) | milestone 4 (figures) |
   | Lua | MIT | OK | planned |
   | Python, NumPy, SciPy (interim, embedded) | PSF / BSD-style | OK | planned |
   | MDAnalysis (interim) | check the bundled version's licence | verify before bundling | planned |
   | Barnaba | GPL v3 | never bundled; replaced by MDBridge's own base-pair annotation | - |

   Code ported from permissive projects keeps its notice in the source file, in
   `THIRD_PARTY_NOTICES.md` and in `packaging/licenses/`:

   | Source | Licence | Ported into | Since |
   |---|---|---|---|
   | NumPy | BSD-3 | `numpy.interp` (`src/io/external.cpp`), `numpy.linspace` and `searchsorted` (`src/render/core/colormap.cpp`, milestone 4), pairwise summation of float32 / float64 arrays (`src/core/include/mdbridge/core/detail/numpy_sum.hpp`, used by `src/io/mrc.cpp` and the analyses (`src/analysis/spatial.cpp`, `rmsd.cpp`, `rmsd_qcp.cpp`, `energy.cpp`) and `src/render/molecule/ndimage.cpp`) | milestone 2 (summation shared since milestone 3) |
   | NumPy 2.5.5 `numpy.random` | BSD-3 / NCSA; PCG (MIT); ziggurat from Julia (MIT) | SeedSequence, PCG64, `Generator.uniform` / `normal` / `integers` and the ziggurat tables (`src/render/core/numpy_random.cpp`, `ziggurat_constants.inc`; `packaging/licenses/numpy-random-LICENSE.txt`) | milestone 4 |
   | SciPy 1.18.1 | BSD-3 (and the ndimage package's BSD-style notice) | `scipy.ndimage` Gaussian filter (mode "constant"), `uniform_filter1d` (mode "constant"), order-1 interpolation and order-1 `affine_transform` (mode "constant") (`src/render/molecule/ndimage.cpp`, used by the molecular surface `src/render/molecule/surface.cpp`, the density shadow `src/render/molecule/density_shadow.cpp`, the software renderer `src/render/molecule/software.cpp` and the Vulkan renderer's halo weights `src/render/vulkan/scene_renderer.cpp`) | milestone 4 |
   | CPython | PSF 2 | textwrap / splitlines (`src/core/text.cpp`), json.loads scanner (`src/core/json.cpp`), fnmatch / glob (`src/core/atom_selection.cpp`, `src/io/cif.cpp`) | milestone 2 |
   | Matplotlib 3.11.2 | Matplotlib licence (PSF-style); its ColorBrewer schemes Apache 2.0 | `to_rgba` / `to_hex`, colormaps (`src/render/core/color.cpp`, `colormap.cpp`), tick locators and formatters, autoscaling and colour-bar limits (`src/render/core/ticks.cpp`), the font manager's font choice (`src/render/core/text.cpp`), the figure layout (constrained layout, text boxes, axis decorations, legends, colour-bar placement; `src/render/core/figure.cpp`, `figure_layout.cpp`, `display_list.cpp`) and the colour tables of `src/render/core/mpl_colors.inc` (`tools/gen_mpl_tables.py`; xkcd names and viridis family CC0, CSS4 names W3C) | milestone 4 |
   | chemfiles, based on xdrfile 1.1.4 | BSD-3, BSD-2 | XTC decompression (`src/io/gromacs.cpp`); xdrfile releases before 1.1.4 were LGPL and must not be used | milestone 2 |
   | OpenMM 8.1.2 force-field files | MIT | residue / atom names and elements in `src/core/include/mdbridge/core/detail/data_tables_independent.inc` (with the wwPDB CCD, CC0, and CIAAW atomic weights; `data/sources/ATTRIBUTION.md`); the amber14 residue templates' charges and Lennard-Jones parameters in `src/analysis/energy_lj_tables.inc` (energy analysis, since AN-25) | milestone 2 |

   Formats implemented from their published descriptions, without porting code (no notice
   needed): the TrueType / OpenType tables of the bundled fonts (`src/render/core/font.cpp`: cmap,
   glyf, hmtx, name, GDEF, GSUB ligatures, GPOS kerning) follow the OpenType specification, with the
   fixed-point rounding FreeType (FTL / GPL 2) and HarfBuzz (MIT) document for scaling, checked
   against Matplotlib's output (no FreeType or HarfBuzz code used); the `.npy` / `.npz` files of
   the analysis cache (`src/analysis/cache.cpp`) follow numpy's documented format
   (`numpy.lib.format`, NEP 1) and the ZIP application note; zlib does the compression.

   MDAnalysis (LGPL) defines the behaviour MDBridge reproduces, but none of its code may be
   copied or translated: reimplement from observed behaviour or from behaviour-only descriptions
   written by someone else (clean room). Its data tables are not used either: residue, atom-name,
   element and mass tables are built from independent public sources (`data/sources/`), so
   MDBridge deliberately differs from MDAnalysis where those tables differ.

   Test-only tools are not bundled: GoogleTest (BSD-3); Mesa's lavapipe Vulkan driver (MIT), on
   Linux CI the distribution's package and on Windows CI MSYS2's build with the Khronos Vulkan
   loader (Apache 2.0) and the MSYS2 runtime libraries it needs (LLVM: Apache 2.0 with the LLVM
   exception; libstdc++ / libgcc with the GCC runtime exception; libiconv and gettext: LGPL; zlib,
   zstd, libxml2, libffi: permissive), fetched with their licence folders by
   `tools/get_lavapipe.py` (pinned by SHA-256) for the Vulkan tests only and never shipped; the
   `ffmpeg` / `ffprobe` executables of the parity and regenerate jobs (Ubuntu's package, which may
   be GPL) run MDATK's video and read metadata, and are never part of a build. External programs MDBridge writes scripts
   for (cpptraj, Amber, GROMACS, Rosetta, Phenix, Coot, ...) are not distributed with it and keep
   their own licences.
