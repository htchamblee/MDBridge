# Building MDBridge

## Requirements

| | Linux | Windows |
|---|---|---|
| Compiler | GCC 11 or newer (RHEL 9's) | Visual Studio 2022 (MSVC), C++ workload |
| CMake | 3.25 or newer | 3.25 or newer (Visual Studio's own is fine) |
| Qt | 6.5 or newer with Vulkan support (releases use 6.11.3) | same, the `msvc2022_64` build |
| Build tool | Ninja (or Make) | Visual Studio's MSBuild |
| Python 3 | parity scripts, smoke test, shader compiler wrapper | same |

GCC and MSVC are the only supported compilers. Clang is not supported: builds and CI do not use it,
and code need not compile with it.

Everything else is fetched by CMake at configure time into the build folder, pinned to a release
and checked against its SHA-512 (`cmake/Dependencies.cmake`): toml++, CLI11, the Vulkan headers,
zlib and GoogleTest. Nothing is installed on the system. The Vulkan *loader* is not needed to build: the
program loads the user's at run time.

### Getting Qt

Qt is not fetched (building it per checkout would take an hour). Any of:

- the official binaries with `tools/get_qt.py`, which is how CI gets them (it downloads the online
  installer's archives without an account and checks their SHA-1):
  ```bash
  python3 -m pip install py7zr
  python3 tools/get_qt.py --version 6.11.3 --arch linux_gcc_64 --out ~/Qt/6.11.3/gcc_64 --archive qtbase --archive icu --archive qtsvg
  python  tools\get_qt.py --version 6.11.3 --arch win64_msvc2022_64 --out C:\Qt\6.11.3\msvc2022_64 --archive qtbase --archive qtsvg
  ```
  Qt modules used: Qt Base (Core, Gui, Widgets; Gui also writes the PDF figures, and its JPEG, GIF
  and ICO image-format plugins come with it) and Qt SVG (`qtsvg`: the figures' SVG output and the
  SVG image plugin). Qt Image Formats (`qtimageformats`) is not needed: TIFF figures are written by
  MDBridge itself. The JPEG plugin (`plugins/imageformats/qjpeg`) must be deployed with the program
  for `.jpg` charts and previews; the AppImage and Windows scripts check it.
  ([aqtinstall](https://github.com/miurahr/aqtinstall) does the same, but its release 3.3.0 does not
  read the repository layout of Qt 6.11 for Windows);
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

Qt's CMake files need the OpenGL headers (RHEL / Rocky: `libglvnd-devel mesa-libGL-devel`; Debian /
Ubuntu: `libgl-dev`), and the interface needs Qt's X11 libraries (most desktops have them; on RHEL /
Rocky: `xcb-util-cursor xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil
libxkbcommon-x11`). Useful options: `-DMDBRIDGE_WARNINGS_AS_ERRORS=ON` (as in CI),
`-DMDBRIDGE_BUILD_APP=OFF` (only the Qt-free libraries - `core`, `io`, `analysis`, the render core
and panels, molecule geometry, software and offscreen Vulkan renderers, video encoder - the command
line's MDATK commands and their tests: no Qt needed; `build/tests/mdbridge_commands_cli` then offers
`init`, `info`, `list`, `analyze`, `export-pdb` and `density info` for `tests/parity/cli_parity.py
--mdbridge`; `render`, `charts` and `preview` need the Qt build).

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

## Video encoding (FFmpeg)

Videos are encoded in process with FFmpeg's libraries (`libavcodec`, `libavformat`, `libavutil`,
`libswscale`) and the AV1 encoder SVT-AV1, built by MDBridge's own scripts as an LGPL-only build
(no GPL or non-free parts; `LICENSING.md`). Versions are pinned in `packaging/ffmpeg/versions.env`
(FFmpeg 9.0.2, SVT-AV1 4.1.0; sources checked against their commit hashes). Build them once, then
point CMake at the result:

```bash
packaging/ffmpeg/build_ffmpeg.sh --prefix build/ffmpeg/prefix        # Linux: ~5 min; gcc, make, cmake, nasm, pkg-config
cmake ... -DMDBRIDGE_FFMPEG=ON -DMDBRIDGE_FFMPEG_ROOT=$PWD/build/ffmpeg/prefix
```

```powershell
packaging\ffmpeg\windows\build_ffmpeg.ps1 -Prefix build\ffmpeg\prefix -Work C:\mdbff   # Windows: vcpkg, ~30 min
packaging\windows\build_windows.ps1 -QtDir C:\Qt\6.11.3\msvc2022_64 -FFmpegRoot build\ffmpeg\prefix
```

Tests of the encoder and the video loop (`mdbridge_video_tests`, `QtVideo.*`) skip themselves
without the libraries; build with `-DMDBRIDGE_FFMPEG=ON` to run them (as CI does).
`MDBRIDGE_FFMPEG` is `AUTO` by default (used when `cmake/FindFFmpegLGPL.cmake` finds the libraries),
`ON` (required) or `OFF`. Without the libraries everything else builds and video encoding reports
`video encoding is not available in this build`. The find module refuses libraries configured with
`--enable-gpl` or `--enable-nonfree` (it reads the configuration string in libavcodec), and the
program checks the same at run time. Installs bundle the shared libraries (Linux: `lib/`, so the
AppImage's `usr/lib`; Windows: next to `mdbridge.exe`) with their licences in `licenses/`; each
script also writes the corresponding FFmpeg source archive next to the prefix
(`ffmpeg-<version>-source.tar.xz`, `ffmpeg-<version>-vcpkg-source.zip`), which CI attaches to
releases. CI caches the prefix by the hash of `packaging/ffmpeg/`. The Linux script runs on Rocky
Linux 9 in CI so that the libraries need glibc 2.34 at most; built elsewhere they need that
system's glibc.

## Offline and mirrored builds

Point FetchContent at checkouts of the same tags instead of downloading:

```bash
cmake ... -DFETCHCONTENT_SOURCE_DIR_TOMLPLUSPLUS=/src/tomlplusplus   # v3.4.0
          -DFETCHCONTENT_SOURCE_DIR_CLI11=/src/CLI11                 # v2.7.2
          -DFETCHCONTENT_SOURCE_DIR_VULKANHEADERS=/src/Vulkan-Headers  # vulkan-sdk-1.4.363.0
          -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/src/googletest       # v1.18.0
          -DFETCHCONTENT_SOURCE_DIR_ZLIB=/src/zlib                   # v1.3.1
```

A checkout of zlib's tag works as well as its release archive:
`git clone --branch v1.3.1 https://github.com/madler/zlib /src/zlib` (commit
`51b7f2abdade71cd9bb0e7a373ef2610ec6f9daf`). MDBridge compiles zlib's C sources itself (target
`mdbridge_zlib`) and never runs zlib's own CMake project, so the checkout stays unmodified.

Installed packages can be used instead with `-DMDBRIDGE_USE_SYSTEM_DEPS=ON` (zlib through CMake's
`find_package(ZLIB)`, at least 1.2.11).

## Tests

- `ctest`: unit tests of `core` (every fixture in `tests/data/projects` must read exactly as MDATK
  reads it: `tests/data/expected`), of the Vulkan probe, and a check that the committed SPIR-V
  matches the shader sources.
- `packaging/smoke_test.py EXE [--gui]`: a finished build as users get it (AppImage, install tree
  or `mdbridge.exe`): command line, project files against MDATK's reading, the figure fonts, a
  `render` of the demo project (an AV1 video read back with `mdbridge probe-video` and, when it is
  on `PATH`, `ffprobe`; charts), `charts` and `preview` without a display, and the interface with
  the software and the automatic 3D view (screenshots that must have content). CI runs it on the
  build machine, on Ubuntu 22.04, on a bare Rocky Linux 9 (no X11, fonts or Vulkan) and on Windows.
- Parity with MDATK (CI job `parity`): MDATK v0.0.4.1 writes fresh project files and records what it
  reads from them (`tests/parity/make_projects.py`, `reference.py`); `tests/parity/compare.py`
  checks that MDBridge reads the same. `tests/data/README.md` explains how to regenerate the
  committed fixtures.
- Reproducible test data (CI job `regenerate`): `tests/parity/regenerate_all.py` runs every
  generator of `tests/data/README.md` and `tools/gen_*.py` into a fresh folder and compares the
  result with the committed files byte for byte (unified-diff excerpts, exit status 1 on any
  difference). CI runs it on Ubuntu 24.04 with Python 3.13, MDATK v0.0.4.1 and the pinned
  environment of `tests/parity/requirements-regenerate.txt` (MDAnalysis 2.11.0, numpy 2.5.5,
  scipy 1.18.1, pytest). Locally, about 25 seconds:
  ```bash
  python3 -m venv /tmp/regen && /tmp/regen/bin/pip install --no-deps -r tests/parity/requirements-regenerate.txt
  python3 tests/parity/regenerate_all.py --mdatk /tmp/mdatk --out /tmp/fresh --python /tmp/regen/bin/python
  ```

The interface has check options for automation: `mdbridge gui PROJECT --renderer vulkan
--screenshot view.png` renders a few frames, saves the 3D view, prints a JSON report and quits;
`--quit-after SECONDS` quits after a while. `QT_QPA_PLATFORM=offscreen` runs without a display
(software view). `MDBRIDGE_VULKAN_VALIDATION=1` enables the Vulkan validation layer if installed.

Vulkan tests (`ctest -L vulkan`: the offscreen renderer and the Vulkan scene renderer against the
software renderer, `tests/render/test_vulkan_*.cpp`): on whatever Vulkan 1.3 driver the machine has;
without one they are skipped with the reason, and with `MDBRIDGE_REQUIRE_VULKAN=1` (set in CI) they
fail instead. Linux: Mesa's lavapipe (the `mesa-vulkan-drivers` package; RHEL / Rocky also
`vulkan-loader`) runs them on the CPU; the window test `mdbridge_vulkan_window_tests` needs a display
(`xvfb-run -a`, package `xorg-x11-server-Xvfb` / `xvfb`) and `MDBRIDGE_REQUIRE_VULKAN_WINDOW=1` to fail
rather than skip. Lavapipe is a test-only driver: it is never bundled. Windows:
`python tools/get_lavapipe.py --out build\lavapipe` fetches a test-only lavapipe (pinned, SHA-256
checked); then set `PATH=build\lavapipe;%PATH%` and `VK_DRIVER_FILES=build\lavapipe\lvp_icd.x86_64.json`.
In an elevated (administrator) shell the Vulkan loader ignores `VK_DRIVER_FILES`; register the
manifest instead (its full path as a DWORD value 0 under `HKLM\SOFTWARE\Khronos\Vulkan\Drivers`, as CI does).
`MDBRIDGE_VULKAN_DEVICE=<part of the name>` picks a device when there are several.

## Shaders

The Vulkan renderer's GLSL shaders (`src/render/shaders`: every `*.vert`, `*.frag`, `*.comp`, plus
shared `*.glsl` includes) are compiled to SPIR-V ahead of time into `src/render/shaders/spirv.hpp`,
which is committed, so builds need no shader compiler. After adding or changing a shader run
`python tools/compile_shaders.py` (needs `glslangValidator` from the Vulkan SDK or the
`glslang-tools` / `glslang` package; `--check` only compares); the test `render.shaders_up_to_date`
fails until then.
