# Building MDBridge

## Requirements

| | Linux | Windows |
|---|---|---|
| Compiler | GCC 11 or newer (RHEL 9's) | Visual Studio 2022 (MSVC), C++ workload |
| CMake | 3.25 or newer | 3.25 or newer (Visual Studio's own is fine) |
| Qt | 6.5 or newer with Vulkan support (releases use 6.11.3) | same, the `msvc2022_64` build |
| Build tool | Ninja (or Make) | Visual Studio's MSBuild |
| Python 3 | parity scripts, smoke test, shader compiler wrapper, table generators (none needed to build) | same |

GCC and MSVC are the only supported compilers. Clang is not supported: builds and CI do not use it,
and code need not compile with it.

Everything else is fetched by CMake at configure time into the build folder, pinned to a release
and checked against its SHA-512 (`cmake/Dependencies.cmake`): toml++, CLI11, the Vulkan headers,
zlib and GoogleTest. Nothing is installed on the system. The Vulkan *loader* is not needed to build: the
program loads the user's at run time.

One dependency is vendored instead: the FFT that SciPy 1.18.1 runs (ducc0 0.40.0, BSD-3-Clause),
copied unmodified from SciPy's source release into `third_party/ducc0` by `tools/vendor_ducc0.py`
and compiled with the numerics library. CMake checks every file against
`third_party/ducc0/SHA256SUMS` when it configures, so an edited copy stops the build;
`python -I -B tools/vendor_ducc0.py scipy-1.18.1.tar.gz --check` compares the copy with SciPy's
source release (PyPI), and without `--check` writes it again.

Lua 5.4.9 (MIT), which runs the workflow step generators, is vendored the same way: PUC-Rio's
release, unmodified, in `third_party/lua` (`tools/vendor_lua.py lua-5.4.9.tar.gz [--check]`, the
tarball from https://www.lua.org/ftp/, checked against its SHA-256; CMake checks
`third_party/lua/SHA256SUMS`). It is compiled as C++ (so Lua errors unwind through the host's C++
frames), only its core and the base, string, table, math and utf8 libraries, by
`src/workflow/lua/CMakeLists.txt`, with GCC's `-Wall -Wextra` (GCC 11 gives no warning) and
MSVC's default `/W3`; MDBridge's stricter warnings and warnings-as-errors apply to the host code
only. To update Lua, change the version and SHA-256 in `tools/vendor_lua.py` and the
version in `src/workflow/lua/CMakeLists.txt`, run the script, and check that
`src/workflow/lua/lua_seed.hpp` still matches `luai_makeseed` in `lstate.c`.

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
`demo`, `init`, `info`, `list`, `clean`, `analyze`, `export-pdb`, `citations`, `density info` / `fit`
/ `simulate`, `pdb`, `build` and `dict` for `tests/parity/cli_parity.py --mdbridge`; `render`,
`charts`, `preview` and `demo --run` need the Qt build, `dict pubchem` the `mdbridge` program's
network access).

### Sanitizers (AddressSanitizer, UndefinedBehaviorSanitizer)

`-DMDBRIDGE_SANITIZE=address,undefined` (any list GCC's `-fsanitize=` takes) compiles every target
with `-fsanitize=... -fno-omit-frame-pointer -fno-sanitize-recover=all`: MDBridge, its tests, the
fetched libraries and the vendored ducc0 alike (Qt, FFmpeg and the Vulkan driver are prebuilt and
stay uninstrumented, but their allocations are still checked for leaks). Configure the C compiler
too, so zlib is built by the same GCC as the rest (one sanitizer runtime); RHEL / Rocky need the
`libasan libubsan` packages. CI's `sanitizers` job does this on Rocky Linux 9 on every push, with the
tests in three shards of their own builds (`ctest ... -I 1,,3`, `2,,3`, `3,,3`: every third test):

```bash
cmake -S . -B build/sanitize -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=gcc-11 \
  -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_PREFIX_PATH=$HOME/Qt/6.11.3/gcc_64 -DMDBRIDGE_SANITIZE=address,undefined
cmake --build build/sanitize
export MDBRIDGE_REQUIRE_VULKAN=1 MDBRIDGE_TIME_SCALE=5 \
  VK_LOADER_DISABLE_DYNAMIC_LIBRARY_UNLOADING=1 \
  ASAN_OPTIONS=detect_leaks=1:abort_on_error=1:detect_stack_use_after_return=1:check_initialization_order=1 \
  UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 LSAN_OPTIONS=suppressions=$PWD/tests/sanitizers/lsan.supp
ctest --test-dir build/sanitize --output-on-failure --parallel 4
MDBRIDGE_REQUIRE_VULKAN_WINDOW=1 xvfb-run -a build/sanitize/tests/mdbridge_vulkan_window_tests
```

- Any finding (memory error, leak at exit, undefined behaviour) ends the program with an error, so
  the test fails. A finding in MDBridge's code is a bug to fix, never suppressed;
  `tests/sanitizers/lsan.supp` suppresses only leaks of prebuilt system libraries, each with its
  reason (today one: libdbus, reached through Qt's xcb platform in the Vulkan window test).
- `MDBRIDGE_TIME_SCALE` (a number of at least 1, default 1) multiplies the tests' time limits
  (`tests/core/time_scale.hpp`): the instrumented code runs several times slower, so a test that
  checks something is fast (`AtomSelection.MillionAtoms`: one second) or happens at once (a
  cancel, the log view never stalling the event loop) checks it against a proportionally larger
  budget instead of being skipped.
- `VK_LOADER_DISABLE_DYNAMIC_LIBRARY_UNLOADING=1` keeps the Vulkan driver (lavapipe) loaded until
  the program ends. Otherwise the loader unloads it at `vkDestroyInstance`, and memory its threads
  keep for the life of the process is reported as a leak from an `<unknown module>`.
- Some AddressSanitizer runtimes stop at start-up ("Shadow memory range interleaves with an
  existing memory mapping") when `/proc/sys/vm/mmap_rnd_bits` is above 28, as newer kernels set;
  Rocky 9's GCC 11 runs with 32 (CI). Where it does not: `sudo sysctl vm.mmap_rnd_bits=28`.
- About 40 minutes to build with two jobs and 13 minutes for `ctest --parallel 4` on four cores
  (the build is RelWithDebInfo: optimised, with line numbers in the reports).

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

Desktop libraries Qt needs at start (OpenGL / EGL dispatch, X11, xcb, fontconfig, freetype and what
they need) are copied from the build machine into `usr/lib/fallback/<name>/<name>`, one folder per
library. `AppRun` puts a copy's folder on the library path only when the host lacks that library
(not in `ldconfig -p`, the standard library folders or `LD_LIBRARY_PATH`), so on a desktop the
host's copies are used and a minimal or server install gets just the missing ones. The GL dispatch
libraries (glvnd: libGLdispatch, libGLX, libGL, libEGL, libOpenGL, libGLES*) count as one unit:
if the host lacks any member the AppImage carries, the whole bundled set is used, else the host's.
`MDBRIDGE_APPRUN_DEBUG=1 ./MDBridge-*.AppImage graphics --json` names the glvnd set and the copies
supplied on stderr and lists every library loaded (`libraries`). The smoke test checks this on every host it runs on
(`--expect-fallback LIBS|none` asserts which copies must be used there).

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
  matches the shader sources. The chemistry and numerics suites (`mdbridge_chem_tests`,
  `mdbridge_numerics_tests`) compare with the cases recorded from MDATK and SciPy (`tests/data/chem`,
  `tests/data/numerics`): byte for byte, or for fitted results within the tolerances of
  `docs/specs/m6-plan.md` (owner decision D10). The FFT test is bit for bit with glibc's libm; with
  another C library (MSVC's) it and the Powell tests allow the tolerances `docs/DIFFERENCES.md`
  "Numerics" gives.
- GUI tests (`ctest -L gui`; Qt Test, part of Qt's `qtbase`, test-only), one program per part of the
  interface, all on Qt's `offscreen` platform (no display needed; `MDBRIDGE_GUI_TEST_PLATFORM=xcb`
  shows the windows) with `QT_FORCE_STDERR_LOGGING=1` (so Windows prints their report). The window
  and page programs use Qt's Fusion style and the bundled IBM Plex Serif at 9 pt, as the
  application does (the fonts are compiled into the programs: no system fonts are needed, and the
  layout tests measure the application's own font):

  | ctest | program | what |
  |---|---|---|
  | `gui.cases` | `mdbridge_gui_tests` | the window driven by the action scripts of `tests/data/gui` (format: `tests/data/README.md`, "gui"), the saved files and messages compared with MDATK's; one case: `mdbridge_gui_tests cases:harness/open_save_default` |
  | `gui.window` | `mdbridge_gui_window_tests` | the window shell: layout, shortcuts, menus, questions, recovery, external changes, paths shown natively, keyboard access and accessible names, settings locations, pages fitting the side bar |
  | `gui.window.scale150`, `gui.window.scale200` | same | the side-bar fit, accessible names, Tab order and screenshots at `QT_SCALE_FACTOR` 1.5 and 2 |
  | `gui.canvas` | `mdbridge_gui_canvas_tests` | the canvas layout editor (drags, wheel, keyboard, undo) |
  | `gui.actions` | `mdbridge_gui_actions_tests` | Render, still, exports, cleaning, desktop actions, result dialogs (the video parts need FFmpeg) |
  | `gui.pages`, `gui.group_pages` | `mdbridge_gui_pages_tests`, `mdbridge_gui_group_pages_tests` | the settings pages |
  | `gui.forms`, `gui.problems`, `gui.jobs` | `mdbridge_gui_forms_tests`, `..._problems_tests`, `..._jobs_tests` | field widgets, problems as you type, job runner and log view |
  | `gui.tools` | `mdbridge_gui_tools_tests` | the PDB clean-up and structure builder windows: the files they write against the commands', *Check*, the fit on a worker thread and its Cancel, closing during a build, a host other than the analysis window, names and Tab order |
  | `gui.tools.scale150`, `gui.tools.scale200` | same | the tool windows' names, column fit and screenshots at `QT_SCALE_FACTOR` 1.5 and 2 |
  | `gui.script` | `mdbridge_gui_script_tests` | `mdbridge gui PROJECT --script FILE`, the action scripts of the smoke test |

  Run one function: `QT_QPA_PLATFORM=offscreen build/release/tests/mdbridge_gui_window_tests
  tabOrderVisitsEveryField`. Screenshots (every page, canvas drags, result dialogs; never compared)
  go to `$MDBRIDGE_GUI_SCREENSHOTS` (default `<build>/gui-screenshots`); CI uploads them as
  `gui-test-screenshots-linux` / `-windows`. The tests keep the recent list and recovery copies in
  a temporary folder (`XDG_CONFIG_HOME` / `APPDATA`).
- `packaging/smoke_test.py EXE [--gui]`: a finished build as users get it (AppImage, install tree
  or `mdbridge.exe`): command line, project files against MDATK's reading, the figure fonts, a
  `render` of the demo project (an AV1 video read back with `mdbridge probe-video` and, when it is
  on `PATH`, `ffprobe`; charts), `charts` and `preview` without a display, and the interface with
  the software and the automatic 3D view (screenshots that must have content), and an editing case
  run through the window with `mdbridge gui PROJECT --script FILE` (format in
  `src/app/gui/script_runner.hpp`): a field of the Outputs page typed in, saved (the file must
  change in exactly that line), *Render still*, a screenshot of every page. For an AppImage it
  also checks which libraries came from the AppImage's fallback copies (see "AppImage").
  `--screenshot-dir DIR` keeps the screenshots (CI: artifacts `gui-screenshots-<platform>`);
  `--gui-exe GUI_EXE` starts the program a desktop shortcut starts on a project file (Windows CI:
  the Start-menu entry's target, after checking the entry exists and starts MDBridge). CI runs it on
  the build machine, on Ubuntu 22.04, on a bare Rocky Linux 9 (no X11, fonts or Vulkan) and on the
  installed Windows program. The Windows installer registers no file type (as MDATK's; Windows
  could only claim every `.toml` file): projects open with *Open with*, by drag and drop, or from
  the window.
- Parity with MDATK (CI job `parity`): MDATK v0.0.4.1 writes fresh project files and records what it
  reads from them (`tests/parity/make_projects.py`, `reference.py`); `tests/parity/compare.py`
  checks that MDBridge reads the same. `tests/data/README.md` explains how to regenerate the
  committed fixtures.
- Reproducible test data (CI job `regenerate`): `tests/parity/regenerate_all.py` runs every
  generator of `tests/data/README.md` and `tools/gen_*.py` into a fresh folder and compares the
  result with the committed files byte for byte (unified-diff excerpts, exit status 1 on any
  difference). CI runs it on Ubuntu 24.04 with Python 3.13, MDATK v0.0.4.1, MDATK 0.1.1 (for the
  steps that record a later MDATK's project files) and the pinned environment of
  `tests/parity/requirements-regenerate.txt` (MDAnalysis 2.11.0, numpy 2.5.5, scipy 1.18.1, pytest);
  MDATK's render commands need the `ffmpeg` executable, and MDATK's GTK 3 window (the `gui` step)
  needs GTK 3, Xvfb and PyGObject / pycairo built against it (`tests/parity/requirements-gui.txt`;
  `.github/workflows/ci.yml` has the commands):
  ```bash
  python3 -m venv /tmp/regen && /tmp/regen/bin/pip install --no-deps -r tests/parity/requirements-regenerate.txt
  python3 tests/parity/regenerate_all.py --mdatk /tmp/mdatk --mdatk-latest /tmp/mdatk-latest --out /tmp/fresh \
    --python /tmp/regen/bin/python --pytest-python /tmp/regen/bin/python --gui-python /tmp/regen/bin/python
  ```
  `--list` names the steps, `--only STEP ...` runs some of them. The chemistry steps (milestone 6)
  are `chem_tables` (`tools/gen_chem_tables.py`: the builder's and the dictionaries' tables from
  `data/sources/`, compared with MDATK's embedded tables), `pdbtools`, `annotate`, `dict`, `builder`
  and `density_fit` (MDATK's chem layer with `tests/parity/improvements_chem.py`'s patches), the
  `pdb_*`, `build_*`, `dict_*` and `density_*` cases of `cli`, `numerics` (SciPy alone) and the
  `tools` cases of `gui`. All but `chem_tables` and `pdbtools` (pure Python) pin numpy's and
  OpenBLAS' kernels before numpy loads (`tests/parity/chem_common.py` `pin_kernels`: OpenBLAS's Nehalem kernels, one thread, numpy's
  x86-64 dispatch limited to X86_V2), so their numbers are the same on most machines. Where the
  bytes of a fitted result (density fits, Kabsch, least-squares fits) still differ, the job compares
  its numbers within the chem tolerance; on a runner whose pinned kernels give other bits than the
  reference machine's (its probe: the density fits' flips), the results whose path depends on them
  get the wider comparison the script describes (`KERNEL_PROBE`, `PATH_DEPENDENT`), and the job says
  so. The extracts in `data/sources/` are inputs, not regenerated: `tools/extract_ccd.py` cuts them
  from the wwPDB CCD (`components.cif.gz`), `tools/extract_pdb_entry.py` from a PDB entry's mmCIF and
  `tools/extract_arnott.py` from AmberTools' `nucgen.dat` (the last two check the download's
  SHA-256); each folder's README names the files and dates.

The interface has check options for automation: `mdbridge gui PROJECT --renderer vulkan
--screenshot view.png` renders a few frames, saves the 3D view, prints a JSON report and quits;
`--quit-after SECONDS` quits after a while; `--script FILE [--script-output DIR]` runs an action
script through the window (type values into fields, save, render a still, screenshots of every
page; format in `src/app/gui/script_runner.hpp`), prints a JSON report and quits with 0 when every
action succeeded. `QT_QPA_PLATFORM=offscreen` runs without a display
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
