# Building MDBridge

## Requirements

| | Linux | Windows |
|---|---|---|
| Compiler | GCC 11 or newer (RHEL 9's) | Visual Studio 2022 (MSVC), C++ workload |
| CMake | 3.24 or newer | 3.24 or newer (Visual Studio's own is fine) |
| Qt | 6.5 or newer with Vulkan support (releases use 6.11.3) | same, the `msvc2022_64` build |
| Build tool | Ninja (or Make) | Visual Studio's MSBuild |
| Python 3 | parity scripts, smoke test, shader compiler wrapper | same |

Everything else is fetched by CMake at configure time into the build folder, pinned to a release
and checked against its SHA-512 (`cmake/Dependencies.cmake`): toml++, CLI11, the Vulkan headers and
GoogleTest. Nothing is installed on the system. The Vulkan *loader* is not needed to build: the
program loads the user's at run time.

### Getting Qt

Qt is not fetched (building it per checkout would take an hour). Any of:

- the official binaries with [aqtinstall](https://github.com/miurahr/aqtinstall), which is how CI gets
  them:
  ```bash
  python3 -m pip install aqtinstall
  python3 -m aqt install-qt linux desktop 6.11.3 linux_gcc_64 -O ~/Qt --archives qtbase icu
  python  -m aqt install-qt windows desktop 6.11.3 win64_msvc2022_64 -O C:\Qt --archives qtbase
  ```
- the Qt online installer (Qt 6.11 for GCC 64-bit / MSVC 2022 64-bit);
- a distribution's Qt 6.5+ development packages with Vulkan enabled (e.g. Fedora's
  `qt6-qtbase-devel`). Ubuntu 24.04's Qt (6.4) is too old and RHEL / Rocky 9 has no Qt 6
  packages: use aqtinstall there.

## Linux

```bash
cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/Qt/6.11.3/gcc_64
cmake --build build/release
ctest --test-dir build/release --output-on-failure
./build/release/src/app/mdbridge gui tests/data/projects/full.mdatk.toml
```

The Qt interface needs Qt's X11 libraries on the build machine (most desktops have them; on RHEL /
Rocky: `xcb-util-cursor xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil
libxkbcommon-x11`). Useful options: `-DMDBRIDGE_WARNINGS_AS_ERRORS=ON` (as in CI),
`-DMDBRIDGE_BUILD_APP=OFF` (only the Qt-free `core` library and its tests: no Qt needed).

### AppImage

```bash
packaging/linux/build_appimage.sh --build build/release    # -> dist/MDBridge-<version>-x86_64.AppImage
packaging/linux/check_glibc.sh build/appimage/AppDir 2.34  # newest glibc the bundle needs
xvfb-run -a python3 packaging/smoke_test.py dist/MDBridge-*-x86_64.AppImage --gui
```

The script installs the build into an AppDir and lets [linuxdeploy](https://github.com/linuxdeploy/linuxdeploy)
and its Qt plugin bundle Qt and the libraries not every Linux system has, then packs it with
appimagetool (tool releases pinned by SHA-256 in the script). The AppImage needs the glibc of the
machine it was built on or newer: releases are built on Rocky Linux 9 (glibc 2.34), so they run on
RHEL 9 and Ubuntu 22.04 and newer. The Vulkan loader is never bundled; it comes with the user's
graphics driver.

## Windows

From a PowerShell prompt in the repository:

```powershell
packaging\windows\build_windows.ps1 -QtDir C:\Qt\6.11.3\msvc2022_64
```

It configures with the Visual Studio 2022 generator, builds, runs the tests, installs into
`build\windows\install` (Qt's DLLs via windeployqt and the MSVC runtime next to the programs, so the
installer needs no administrator rights) and writes `dist\MDBridge-<version>-windows-x64.zip` and,
with [Inno Setup 6](https://jrsoftware.org/isinfo.php), the installer `...-setup.exe`. Tests run
from the build folder need Qt's `bin` on `PATH` (the script sets it).

To work in Visual Studio instead: *Open a local folder* on the repository and set
`CMAKE_PREFIX_PATH` to the Qt folder in the CMake settings.

## Offline and mirrored builds

Point FetchContent at checkouts of the same tags instead of downloading:

```bash
cmake ... -DFETCHCONTENT_SOURCE_DIR_TOMLPLUSPLUS=/src/tomlplusplus   # v3.4.0
          -DFETCHCONTENT_SOURCE_DIR_CLI11=/src/CLI11                 # v2.7.2
          -DFETCHCONTENT_SOURCE_DIR_VULKANHEADERS=/src/Vulkan-Headers  # vulkan-sdk-1.4.363.0
          -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/src/googletest       # v1.18.0
```

or use installed packages of them with `-DMDBRIDGE_USE_SYSTEM_DEPS=ON`.

## Tests

- `ctest`: unit tests of `core` (every fixture in `tests/data/projects` must read exactly as MDATK
  reads it: `tests/data/expected`), of the Vulkan probe, and a check that the committed SPIR-V
  matches the shader sources.
- `packaging/smoke_test.py EXE [--gui]`: a finished build as users get it (AppImage, install tree
  or `mdbridge.exe`): command line, project files against MDATK's reading, and the interface with
  the software and the automatic 3D view (screenshots that must have content).
- Parity with MDATK (CI job `parity`): MDATK v0.0.4.1 writes fresh project files and records what it
  reads from them (`tests/parity/make_projects.py`, `reference.py`); `tests/parity/compare.py`
  checks that MDBridge reads the same. `tests/data/README.md` explains how to regenerate the
  committed fixtures.

The interface has check options for automation: `mdbridge gui PROJECT --renderer vulkan
--screenshot view.png` renders a few frames, saves the 3D view, prints a JSON report and quits;
`--quit-after SECONDS` quits after a while. `QT_QPA_PLATFORM=offscreen` runs without a display
(software view). `MDBRIDGE_VULKAN_VALIDATION=1` enables the Vulkan validation layer if installed.

## Shaders

The 3D view's GLSL shaders (`src/render/shaders`) are compiled to SPIR-V ahead of time into
`src/render/shaders/spirv.hpp`, which is committed, so builds need no shader compiler. After
changing a shader run `python tools/compile_shaders.py` (needs `glslangValidator` from the Vulkan
SDK or the `glslang-tools` / `glslang` package); the test `render.shaders_up_to_date` fails until then.
