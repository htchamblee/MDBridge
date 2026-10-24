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
