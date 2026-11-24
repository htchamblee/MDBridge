# Third-party software in MDBridge builds

MDBridge itself is all rights reserved (see `LICENSE`). Its builds contain the following
third-party software, each under its own licence; the licence texts are in the `licenses` folder of
a build (`share/doc/mdbridge/licenses` in the AppImage and on Linux, `licenses` next to the
programs on Windows).

| Component | Version | Licence | Use | Source |
|---|---|---|---|---|
| Qt (Core, Gui, Widgets, platform plugins) | 6.11 | LGPL v3 (`LGPL-3.0.txt`, which builds on `GPL-3.0.txt`) | user interface; shared libraries you may replace | https://download.qt.io/official_releases/qt/ |
| toml++ | 3.4.0 | MIT (`tomlplusplus-LICENSE.txt`) | reading project files (compiled in) | https://github.com/marzer/tomlplusplus |
| CLI11 | 2.7.2 | BSD-3-Clause (`CLI11-LICENSE.txt`) | command line parsing (compiled in) | https://github.com/CLIUtils/CLI11 |
| Vulkan-Headers | 1.4.363 | Apache 2.0 or MIT (`Vulkan-Headers-*.txt`) | Vulkan API declarations (compiled in) | https://github.com/KhronosGroup/Vulkan-Headers |
| zlib | 1.3.1 | zlib licence (`zlib-LICENSE.txt`) | reading gzip-compressed density maps (compiled in) | https://github.com/madler/zlib |
| CPython (algorithms only) | 3.13 | Python Software Foundation License 2 (`Python-LICENSE.txt`) | textwrap and str.splitlines rules in `src/core/text.cpp`, the json.loads scanner in `src/core/json.cpp`, fnmatch / glob pattern handling in `src/core/atom_selection.cpp` and `src/io/cif.cpp`, ported to C++ (compiled in) | https://github.com/python/cpython |
| NumPy (interpolation and summation routines only) | 2.5.5 | BSD-3-Clause (`numpy-LICENSE.txt`) | `numpy.interp` ported to C++ in `src/io/external.cpp`, the pairwise summation of float arrays (`loops_utils.h.src`) in `src/io/mrc.cpp` (compiled in) | https://github.com/numpy/numpy |
| chemfiles (XTC decompression only) | commit 920a1347 | BSD-3-Clause (`chemfiles-LICENSE.txt`) | XTC coordinate decompression ported to C++ in `src/io/gromacs.cpp` (compiled in) | https://github.com/chemfiles/chemfiles |
| xdrfile (via chemfiles' port) | as in chemfiles | BSD-2-Clause (`xdrfile-LICENSE.txt`) | the XTC decompression routines chemfiles took from xdrfile, in `src/io/gromacs.cpp` (compiled in) | https://github.com/chemfiles/chemfiles (src/files/XDRFile.cpp) |
| OpenMM (force-field residue and atom names only) | 8.1.2 | MIT (`OpenMM-LICENSE.txt`) | residue and atom names with their elements in `src/core/include/mdbridge/core/detail/data_tables_independent.inc`, used by the selection macros, element and mass guessing and the PDB writer (compiled in); the same tables also use the wwPDB Chemical Component Dictionary (CC0) and CIAAW/IUPAC atomic weights, see `data/sources/ATTRIBUTION.md` | https://github.com/openmm/openmm |

**Qt (LGPL v3).** MDBridge links Qt dynamically. The Qt libraries are separate files: in the
AppImage under `usr/lib` (`./MDBridge-*.AppImage --appimage-extract` unpacks them), on Windows next
to `mdbridge.exe`. You may replace them with your own build of the same Qt version. Qt's own
third-party components (e.g. ICU, HarfBuzz, FreeType, libpng, PCRE2, zlib) are listed with their
licences at https://doc.qt.io/qt-6/licenses-used-in-qt.html. The corresponding source of the Qt
version in these builds (6.11.3) is
https://download.qt.io/official_releases/qt/6.11/6.11.3/single/qt-everywhere-src-6.11.3.tar.xz.

**System libraries.** The Vulkan loader and graphics drivers belong to the user's system and are
not bundled. The AppImage carries libraries of its build system (Rocky Linux 9) that Qt needs and
not every Linux system has (XCB utilities, xkbcommon, ...), and fallback copies of the OpenGL / EGL
dispatch (libglvnd), X11, fontconfig and FreeType libraries, used only on systems without them.
All are separate shared libraries under permissive licences (MIT / X11 style, FreeType's FTL,
D-Bus under AFL 2.1) or the LGPL 2.1 (libsystemd, libgcrypt, needed by Qt's D-Bus support); their
licence texts are in `licenses/system`.

Unit-test tools (GoogleTest, BSD-3-Clause) are used only to build and test MDBridge and are not
part of any build.
