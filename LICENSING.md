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
   | Qt 6 | LGPL v3 | link dynamically; ship the licence and a source pointer | milestone 1 |
   | toml++ | MIT | OK | milestone 1 |
   | CLI11 | BSD-3-Clause | OK | milestone 1 |
   | Vulkan-Headers | Apache 2.0 / MIT | OK (headers only; the Vulkan loader is the user's) | milestone 1 |
   | zlib | zlib licence (permissive) | OK | milestone 2 (density maps) |
   | chemfiles | BSD-3 | not used: dropped in milestone 2 for MDBridge's own readers (its XTC decompression is ported, below) | - |
   | Eigen | MPL 2.0 | OK | planned |
   | nanoflann | BSD-2 | OK | planned |
   | FFmpeg | LGPL (only) | no GPL components (no x264); use OpenH264 / AV1 | planned |
   | Lua | MIT | OK | planned |
   | Python, NumPy, SciPy (interim, embedded) | PSF / BSD-style | OK | planned |
   | MDAnalysis (interim) | check the bundled version's licence | verify before bundling | planned |
   | Barnaba | GPL v3 | never bundled; replaced by MDBridge's own base-pair annotation | - |

   Code ported from permissive projects keeps its notice in the source file, in
   `THIRD_PARTY_NOTICES.md` and in `packaging/licenses/`:

   | Source | Licence | Ported into | Since |
   |---|---|---|---|
   | NumPy | BSD-3 | `numpy.interp` (`src/io/external.cpp`), `numpy.linspace` (`src/render/core/colormap.cpp`, milestone 4), pairwise summation of float32 / float64 arrays (`src/core/include/mdbridge/core/detail/numpy_sum.hpp`, used by `src/io/mrc.cpp` and the analyses (`src/analysis/spatial.cpp`, `rmsd.cpp`, `rmsd_qcp.cpp`, `energy.cpp`)) | milestone 2 (summation shared since milestone 3) |
   | CPython | PSF 2 | textwrap / splitlines (`src/core/text.cpp`), json.loads scanner (`src/core/json.cpp`), fnmatch / glob (`src/core/atom_selection.cpp`, `src/io/cif.cpp`) | milestone 2 |
   | Matplotlib 3.11.2 | Matplotlib licence (PSF-style); its ColorBrewer schemes Apache 2.0 | `to_rgba` / `to_hex`, colormaps (`src/render/core/color.cpp`, `colormap.cpp`) and the colour tables of `src/render/core/mpl_colors.inc` (`tools/gen_mpl_tables.py`; xkcd names and viridis family CC0, CSS4 names W3C) | milestone 4 |
   | chemfiles, based on xdrfile 1.1.4 | BSD-3, BSD-2 | XTC decompression (`src/io/gromacs.cpp`); xdrfile releases before 1.1.4 were LGPL and must not be used | milestone 2 |
   | OpenMM 8.1.2 force-field files | MIT | residue / atom names and elements in `src/core/include/mdbridge/core/detail/data_tables_independent.inc` (with the wwPDB CCD, CC0, and CIAAW atomic weights; `data/sources/ATTRIBUTION.md`) | milestone 2 |

   Formats implemented from their published descriptions, without porting code (no notice
   needed): the `.npy` / `.npz` files of the analysis cache (`src/analysis/cache.cpp`) follow
   numpy's documented format (`numpy.lib.format`, NEP 1) and the ZIP application note; zlib does
   the compression.

   MDAnalysis (LGPL) defines the behaviour MDBridge reproduces, but none of its code may be
   copied or translated: reimplement from observed behaviour or from behaviour-only descriptions
   written by someone else (clean room). Its data tables are not used either: residue, atom-name,
   element and mass tables are built from independent public sources (`data/sources/`), so
   MDBridge deliberately differs from MDAnalysis where those tables differ.

   Test-only tools (GoogleTest, BSD-3) are not bundled. External programs MDBridge writes scripts
   for (cpptraj, Amber, GROMACS, Rosetta, Phenix, Coot, ...) are not distributed with it and keep
   their own licences.
