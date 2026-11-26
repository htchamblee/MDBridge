"""Smoke test of a finished MDBridge build: AppImage, install tree or Windows mdbridge.exe.

    python packaging/smoke_test.py EXE [--gui] [--expect-renderer auto|vulkan|software] [--version X.Y.Z]
                                       [--text-report OUT.json] [--compare-text-report REF.json]
                                       [--screenshot-dir DIR] [--expect-fallback LIBS|none]
                                       [--gui-exe GUI_EXE]

Runs, in order:

* ``EXE --version``
* ``EXE project show`` on the fixtures in tests/data/projects: summaries, rejected files, and the
  JSON of every file against MDATK's reference reading (tests/data/expected)
* ``EXE graphics --json``; for an AppImage (or an AppDir's AppRun) also which libraries it loaded:
  a copy from the AppImage's fallback libraries (usr/lib/fallback) only where the host lacks the
  library, so fontconfig, freetype, OpenGL / EGL, X11 and xcb come from the host whenever it has
  them (packaging/linux/AppRun); ``--expect-fallback libGL.so.1,...`` names copies that must be
  used on this host (``none``: no copy at all)
* ``EXE fonts --json --render PNG``: the bundled figure fonts load from the program, and a text
  image is drawn without a display (Linux: the "offscreen" platform, also under X);
  ``--text-report`` writes the report, ``--compare-text-report`` checks that the image's pixels and
  the text metrics equal those of another machine's report (CI: the bare Rocky Linux system against
  the build machine)
* ``EXE info`` (demo system), ``EXE density info`` (a map) and ``EXE init`` (a new project file)
* rendering without a display, on a project of the demo system with a cartoon representation:
  ``EXE render`` (analyses, a 24-frame MP4 at 640 x 360 and the charts), the video read back with
  ``EXE probe-video`` (the bundled FFmpeg libraries) and, when it is on PATH, ``ffprobe``; ``EXE
  charts --from-cache`` (PNG, SVG, PDF) and ``EXE preview`` (PNG, PDF); ``--screenshot-dir`` keeps
  the preview image
* ``--gui``: ``EXE gui PROJECT --screenshot --preview-screenshot`` on a project of the demo system
  (tests/data) with the software view and the automatic choice (``--expect-renderer``: what the
  automatic choice must give; ``vulkan`` also forces a Vulkan run), checking that the 3D view shows
  the project's structure, not the placeholder (the program says so, and a good part of the pixels
  differ from the background), and that the figure preview was drawn; then
  ``EXE gui PROJECT --quit-after 5`` on a project whose input files are missing (the interface
  reports it and keeps running). On a headless Linux machine run it under ``xvfb-run -a``.
  Then the editing case: ``EXE gui PROJECT --renderer software --script FILE`` (an action script,
  src/app/gui/script_runner.hpp) saves the demo project, enters a new frame rate into the Outputs
  page's field (typed, as a user would), saves again, renders a still (*Render still*) and takes a
  screenshot of every page; the saved file must differ from the first save in exactly that line
  (and the time stamp of the save), and the still and the screenshots must have content. For an
  AppImage the interface's libraries (xcb, X11) are checked as above.
  ``--screenshot-dir``: keep the screenshots there (CI uploads them for review).
* ``--gui-exe GUI_EXE``: the program a desktop shortcut starts (Windows: the Start-menu entry's
  target, mdbridge-gui.exe), started as a double-click on a project or "Open with" does:
  ``GUI_EXE PROJECT`` (with a short script that checks the project was opened and takes a
  screenshot).

Uses the standard library only; exits with 1 when a check fails.
"""

from __future__ import annotations

import argparse
import difflib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tests", "parity"))
import compare  # noqa: E402  (tests/parity/compare.py)

PROJECTS = os.path.join(ROOT, "tests", "data", "projects")
EXPECTED = os.path.join(ROOT, "tests", "data", "expected")
FAILED = []
SCREENSHOT_DIR = None  # --screenshot-dir


def env(**extra):
    e = dict(os.environ, APPIMAGE_EXTRACT_AND_RUN="1")  # AppImages without FUSE (containers, CI)
    if os.name != "nt":
        e.setdefault("LC_ALL", "C.UTF-8")  # Qt warns on other locales
    e.update(extra)
    return e


def check(label: str, ok: bool, detail: str = "") -> bool:
    print(f"{'ok' if ok else 'FAIL':4s} {label}" + (f"\n     {detail}" if detail and not ok else ""))
    if not ok:
        FAILED.append(label)
    return ok


def run(exe, *args, timeout=180, extra_env=None):
    return subprocess.run([exe, *args], capture_output=True, env=env(**(extra_env or {})), encoding="utf-8",
                          errors="replace", timeout=timeout)


def last_json_line(text: str):
    lines = [line for line in text.splitlines() if line.startswith("{")]
    return read_json(lines[-1]) if lines else {}


# ---- the AppImage's fallback libraries (packaging/linux/AppRun) ------------------------------------

# the desktop libraries whose host copies must win (graphics driver, font set-up; libGL and libEGL
# with the rest of the glvnd unit)
HOST_FIRST = ("libGL.so.1", "libEGL.so.1", "libfontconfig.so.1", "libfreetype.so.6", "libX11.so.6", "libxcb.so.1")
EXPECT_FALLBACK = None  # --expect-fallback: None (not checked), set() for none, or the names


def is_appimage(exe: str) -> bool:
    return sys.platform.startswith("linux") and (exe.endswith(".AppImage") or os.path.basename(exe) == "AppRun")


def host_libraries() -> set:
    """The x86-64 libraries this host's dynamic loader finds: its cache, the standard folders and
    LD_LIBRARY_PATH (as AppRun decides, but found independently)."""
    names = set()
    for ldconfig in ("ldconfig", "/sbin/ldconfig"):
        try:
            out = subprocess.run([ldconfig, "-p"], capture_output=True, encoding="utf-8", errors="replace").stdout
        except OSError:
            continue
        for line in out.splitlines():
            parts = line.split()
            if parts and "x86-64" in line:
                names.add(parts[0])
        break
    folders = ["/usr/lib64", "/lib64", "/usr/lib/x86_64-linux-gnu", "/lib/x86_64-linux-gnu"]
    folders += [d for d in os.environ.get("LD_LIBRARY_PATH", "").split(":") if d]
    for folder in folders:
        try:
            names.update(os.listdir(folder))
        except OSError:
            pass
    return names


GLVND = re.compile(r"^lib(GLdispatch|GLX|GL|EGL|OpenGL|GLESv1_CM|GLESv2)\.so\.[0-9]+$")
APPRUN_DEBUG = {"MDBRIDGE_APPRUN_DEBUG": "1"}  # AppRun names the glvnd set it carries on stderr


def bundled_glvnd(stderr: str):
    """The glvnd set AppRun reported (MDBRIDGE_APPRUN_DEBUG), or None when it did not."""
    for line in stderr.splitlines():
        if line.startswith("AppRun: glvnd set:"):
            names = line.split(":", 2)[2].split(";")[0].split()
            return [n for n in names if n != "none"]
    return None


def library_check(label: str, libraries, stderr: str) -> None:
    """The fallback copies an AppImage loaded: only libraries the host lacks, the GL dispatch
    libraries (glvnd) as one unit (all the AppImage's when the host lacks a member of its set,
    else all the host's); fontconfig, freetype, X11 and xcb from the host whenever it has them;
    --expect-fallback: these copies are used."""
    if not libraries:
        check(f"{label}: the program lists the libraries it loaded", False, "no 'libraries' in the report")
        return
    host = host_libraries()
    loaded = {p.rsplit("/", 1)[-1]: p for p in libraries}
    fallback = {n: p for n, p in loaded.items() if "/usr/lib/fallback/" in p}
    print(f"     {label}: {len(libraries)} libraries loaded; copies from the AppImage's fallback: "
          f"{' '.join(sorted(fallback)) or 'none'}")
    wrong = sorted(n for n in fallback if n in host and not GLVND.match(n))
    check(f"{label}: fallback copies only of libraries the host lacks (glvnd as one unit)", not wrong,
          "the host has " + ", ".join(f"{n} (loaded {fallback[n]})" for n in wrong))
    glvnd = bundled_glvnd(stderr)
    if glvnd is None:
        check(f"{label}: AppRun names its glvnd set (MDBRIDGE_APPRUN_DEBUG)", False, stderr[-600:])
    else:
        missing = [n for n in glvnd if n not in host]
        loaded_glvnd = sorted(n for n in loaded if GLVND.match(n))
        mixed = [f"{n}: {loaded[n]}" for n in loaded_glvnd if (n in fallback) != bool(missing)]
        check(f"{label}: glvnd ({', '.join(loaded_glvnd) or 'none loaded'}) all from "
              + (f"the AppImage (the host lacks {', '.join(missing)})" if missing else "the host"),
              not mixed, "; ".join(mixed))
    present = [n for n in HOST_FIRST if n in host and not GLVND.match(n)]
    not_host = [f"{n}: {loaded[n]}" for n in present if n in fallback]
    check(f"{label}: the host's {', '.join(present) or 'desktop libraries (it has none)'} used",
          not not_host, "; ".join(not_host))
    if EXPECT_FALLBACK:
        absent = sorted(n for n in EXPECT_FALLBACK if n not in fallback)
        check(f"{label}: the AppImage's copies of {', '.join(sorted(EXPECT_FALLBACK))} used", not absent,
              "loaded from elsewhere or not at all: " + ", ".join(f"{n} ({loaded.get(n, 'not loaded')})" for n in absent))
    elif EXPECT_FALLBACK is not None:
        check(f"{label}: no fallback copy used", not fallback, ", ".join(sorted(fallback)))


# The demo system (tests/data) as a project: the 3D view shows its reference frame, the preview its
# figure with the example data
DEMO_PROJECT = """format = "mdatk-project"

[inputs]
topology = "demo_dry.parm7"
trajectories = ["demo_dry.nc"]

[[groups]]
name = "Aptamer"
residues = "1-12"

[[groups]]
name = "Target"
residues = "13-24"
color = "#eb6834"
"""
MIN_COVERAGE = 0.1  # the demo structure covers ~17 % of the view in the 1280 x 800 window


def demo_project(work: str) -> str:
    data = os.path.join(ROOT, "tests", "data")
    shutil.copy(os.path.join(data, "topology", "demo_dry.parm7"), work)
    shutil.copy(os.path.join(data, "io", "trajectories", "demo_dry.nc"), work)
    path = os.path.join(work, "demo.mdatk.toml")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(DEMO_PROJECT)
    return path


# The demo project, rendered from the command line: the structure as a cartoon (the 3D renderer),
# a small video, low-resolution charts
RENDER_SETTINGS = """
[visual]
representations = [{style = "cartoon"}]

[output]
out_dir = "out"
timestamp_files = false

[output.video]
width = 640
height = 360
dpi = 60

[output.charts]
dpi = 100
"""


def read_json(text: str):
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return {}


def file_starts_with(path: str, *signatures: bytes) -> bool:
    if not os.path.isfile(path) or os.path.getsize(path) < 200:
        return False
    with open(path, "rb") as fh:
        head = fh.read(64)
    return any(sig in head for sig in signatures)


def render_check(exe):
    """render (video + charts), probe-video / ffprobe, charts and preview on the demo project."""
    with tempfile.TemporaryDirectory() as work:
        project = demo_project(work)
        with open(project, "a", encoding="utf-8") as fh:
            fh.write(RENDER_SETTINGS)
        r = run(exe, "render", "-p", project, "--no-pdb", timeout=600)
        out = read_json(r.stdout)
        video = out.get("video") or ""
        charts = out.get("charts") or []
        detail = f"exit {r.returncode}; {r.stdout[-600:]} {r.stderr.strip()[-900:]}"
        check(f"render: a video and {len(charts)} charts of the demo project",
              r.returncode == 0 and os.path.isfile(video) and len(charts) > 0
              and all(file_starts_with(c, b"\x89PNG") for c in charts), detail)
        check("render: nothing left but the outputs (no temporary video file)",
              not [n for n in os.listdir(os.path.join(work, "out")) if n.startswith(".MDATK_tmp")],
              str(os.listdir(os.path.join(work, "out"))) if os.path.isdir(os.path.join(work, "out")) else "no out")

        r = run(exe, "probe-video", video)
        probe = read_json(r.stdout)
        check(f"probe-video: {probe.get('codec')} {probe.get('width')} x {probe.get('height')}, "
              f"{probe.get('packets')} frames at {probe.get('r_frame_rate')} fps, {probe.get('pix_fmt')}",
              r.returncode == 0 and probe.get("codec") == "av1" and probe.get("width") == 640
              and probe.get("height") == 360 and probe.get("packets") == 24 and probe.get("nb_frames") == 24
              and probe.get("pix_fmt") == "yuv420p" and probe.get("faststart") is True
              and probe.get("audio_streams") == 0 and probe.get("video_streams") == 1,
              f"exit {r.returncode}; {r.stdout[-800:]} {r.stderr[-400:]}")
        ffprobe = shutil.which("ffprobe")
        if ffprobe and os.path.isfile(video):
            f = subprocess.run([ffprobe, "-v", "error", "-select_streams", "v:0", "-count_packets", "-show_entries",
                                "stream=codec_name,width,height,pix_fmt,nb_read_packets", "-of", "json", video],
                               capture_output=True, encoding="utf-8", errors="replace", timeout=120)
            stream = (read_json(f.stdout).get("streams") or [{}])[0]
            check(f"ffprobe agrees: {stream}",
                  f.returncode == 0 and stream.get("codec_name") == "av1" and stream.get("width") == 640
                  and stream.get("height") == 360 and stream.get("nb_read_packets") == "24",
                  f.stderr[-400:])

        r = run(exe, "charts", "-p", project, "--from-cache", "--series", "rmsd/Complex", "--formats", "png,svg,pdf")
        files = read_json(r.stdout).get("charts") or []
        signatures = {".png": b"\x89PNG", ".svg": b"<?xml", ".pdf": b"%PDF"}
        check(f"charts --from-cache: {', '.join(os.path.basename(c) for c in files)}",
              r.returncode == 0 and len(files) == 3
              and all(file_starts_with(c, signatures.get(os.path.splitext(c)[1], b"?")) for c in files),
              f"exit {r.returncode}; {r.stdout[-600:]} {r.stderr.strip()[-600:]}")

        for name, signature in (("preview.png", b"\x89PNG"), ("preview.pdf", b"%PDF")):
            target = os.path.join(work, name)
            r = run(exe, "preview", "-p", project, "--frame", "12", "--output", target)
            check(f"preview --frame 12 --output {name}",
                  r.returncode == 0 and r.stdout.strip() == target and file_starts_with(target, signature)
                  and os.path.getsize(target) > 5000,
                  f"exit {r.returncode}; {r.stdout[-300:]} {r.stderr.strip()[-600:]}")
        if SCREENSHOT_DIR and os.path.isfile(os.path.join(work, "preview.png")):
            os.makedirs(SCREENSHOT_DIR, exist_ok=True)
            shutil.copy(os.path.join(work, "preview.png"), os.path.join(SCREENSHOT_DIR, "cli-preview.png"))


def screenshot(exe, renderer: str, work: str):
    png = os.path.join(work, f"view-{renderer}.png")
    preview_png = os.path.join(work, f"preview-{renderer}.png")
    r = run(exe, "gui", demo_project(work), "--renderer", renderer, "--screenshot", png,
            "--preview-screenshot", preview_png)
    lines = [line for line in r.stdout.splitlines() if line.startswith("{")]
    report = json.loads(lines[-1]) if lines else {}
    detail = f"exit {r.returncode}; {report}; {r.stderr.strip()[-600:]}"
    ok = r.returncode == 0 and report.get("colours", 0) > 100 and os.path.isfile(png) and os.path.getsize(png) > 1000
    check(f"gui --renderer {renderer}: {report.get('description', '?')}", ok, detail)
    check(f"gui --renderer {renderer}: the 3D view shows the project's structure ({report.get('structure')}, "
          f"{100 * report.get('coverage', 0):.0f} % of the pixels)",
          report.get("molecule") is True and report.get("coverage", 0) > MIN_COVERAGE, detail)
    preview = report.get("preview") or {}
    check(f"gui --renderer {renderer}: figure preview {preview.get('width')} x {preview.get('height')}",
          r.returncode == 0 and not preview.get("problem") and preview.get("width", 0) >= 100
          and preview.get("colours", 0) > 50 and os.path.isfile(preview_png) and os.path.getsize(preview_png) > 1000,
          detail)
    if SCREENSHOT_DIR:
        os.makedirs(SCREENSHOT_DIR, exist_ok=True)
        for f in (png, preview_png):
            if os.path.isfile(f):
                shutil.copy(f, SCREENSHOT_DIR)
        with open(os.path.join(SCREENSHOT_DIR, f"report-{renderer}.json"), "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=1, ensure_ascii=False)
    return report


# The editing case of --gui: a field of the Outputs page changed through the window, saved, a still
# rendered and every page photographed (gui --script; src/app/gui/script_runner.hpp)
EDIT_KEY, EDIT_OLD, EDIT_NEW = "output.video.fps", "fps = 24", "fps = 30"


def gui_config_env(work: str) -> dict:
    """A settings folder of the run's own (recent projects, recovery copies), not the user's."""
    folder = os.path.join(work, "config")
    os.makedirs(folder, exist_ok=True)
    return {"APPDATA": folder} if os.name == "nt" else {"XDG_CONFIG_HOME": folder}


def write_script(path: str, actions) -> str:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"format": "mdbridge-gui-script", "version": 1, "actions": actions}, fh, indent=1)
    return path


def saved_changes(before: str, after: str):
    """The lines removed and added between two saves, without the time stamp of the save."""
    diff = [line for line in difflib.ndiff(before.splitlines(), after.splitlines()) if line[:2] in ("- ", "+ ")]
    diff = [line for line in diff if not line[2:].startswith("saved = ")]
    return [line[2:] for line in diff if line.startswith("- ")], [line[2:] for line in diff if line.startswith("+ ")]


def gui_edit_check(exe, appimage: bool):
    with tempfile.TemporaryDirectory() as work:
        project = demo_project(work)
        shots = SCREENSHOT_DIR or os.path.join(work, "shots")
        script = write_script(os.path.join(work, "edit.json"), [
            {"action": "wait"},
            {"action": "save"},
            {"action": "copy_project", "to": "before.mdatk.toml"},
            {"action": "set", "key": EDIT_KEY, "value": 30},
            {"action": "save"},
            {"action": "render_still", "frame": 3, "to": "still.png", "screenshot": "edit-still-dialog"},
            {"action": "screenshot", "name": "edit-page", "pages": True},
        ])
        r = run(exe, "gui", project, "--renderer", "software", "--script", script, "--script-output", shots,
                timeout=900, extra_env={**gui_config_env(work), **APPRUN_DEBUG})
        report = last_json_line(r.stdout)
        detail = f"exit {r.returncode}; {report.get('problem')}; {r.stderr.strip()[-900:]}"
        check(f"gui --script: {len(report.get('actions', []))} actions through the window "
              f"(edit, save, render still, screenshots)", r.returncode == 0 and report.get("ok") is True, detail)
        actions = {a.get("action"): a for a in report.get("actions", [])}
        before = os.path.join(work, "before.mdatk.toml")
        if os.path.isfile(before):
            with open(before, encoding="utf-8") as fh:
                old = fh.read()
            with open(project, encoding="utf-8") as fh:
                new = fh.read()
            removed, added = saved_changes(old, new)
            check(f"gui --script: the saved file changed in exactly one line ({EDIT_OLD} -> {EDIT_NEW})",
                  removed == [EDIT_OLD] and added == [EDIT_NEW] and actions.get("set", {}).get("stored") == 30
                  and report.get("modified") is False,
                  f"removed {removed}, added {added}, stored {actions.get('set', {}).get('stored')}, "
                  f"modified {report.get('modified')}")
        else:
            check("gui --script: the project was saved", False, detail)
        still = os.path.join(work, "still.png")
        check("gui --script: Render still drew the frame",
              file_starts_with(still, b"\x89PNG") and os.path.getsize(still) > 5000, detail)
        if SCREENSHOT_DIR and os.path.isfile(still):
            shutil.copy(still, os.path.join(SCREENSHOT_DIR, "edit-still.png"))
        pages = report.get("pages", [])
        files = actions.get("screenshot", {}).get("files", [])
        check(f"gui --script: a screenshot of each of the {len(pages)} pages",
              len(pages) == 6 and len(files) == len(pages)
              and all(file_starts_with(f, b"\x89PNG") and os.path.getsize(f) > 5000 for f in files),
              f"pages {pages}, files {files}")
        if appimage:
            library_check("interface", report.get("libraries"), r.stderr)


def gui_exe_check(gui_exe):
    """The shortcut's program started on a project file, as a double-click or "Open with" does."""
    with tempfile.TemporaryDirectory() as work:
        project = demo_project(work)
        shots = SCREENSHOT_DIR or os.path.join(work, "shots")
        script = write_script(os.path.join(work, "open.json"), [
            {"action": "wait"},
            {"action": "screenshot", "name": "shortcut-window"},
        ])
        r = run(gui_exe, project, "--script", script, "--script-output", shots, timeout=600,
                extra_env=gui_config_env(work))
        report = last_json_line(r.stdout)
        opened = report.get("project") or ""
        same = bool(opened) and os.path.normcase(os.path.abspath(opened)) == os.path.normcase(os.path.abspath(project))
        check(f"{os.path.basename(gui_exe)} PROJECT: the window opens the project ({opened or '?'})",
              r.returncode == 0 and report.get("ok") is True and same and not report.get("messages"),
              f"exit {r.returncode}; {report}; {r.stderr.strip()[-600:]}")


BUNDLED_FONTS = {"IBMPlexSerif-Regular.ttf", "IBMPlexSerif-Bold.ttf", "IBMPlexSerif-Italic.ttf",
                 "IBMPlexSerif-BoldItalic.ttf", "IBMPlexMath-Regular.ttf", "IBMPlexMono-Regular.ttf",
                 "IBMPlexMono-Bold.ttf"}


def fonts_check(exe, report_path, reference_path):
    with tempfile.TemporaryDirectory() as work:
        png = os.path.join(work, "text.png")
        r = run(exe, "fonts", "--json", "--render", png)
        try:
            report = json.loads(r.stdout)
        except json.JSONDecodeError:
            report = {}
        image = report.get("image", {})
        files = {f.get("file") for f in report.get("faces", [])}
        ok = (r.returncode == 0 and files == BUNDLED_FONTS and os.path.isfile(png) and os.path.getsize(png) > 1000
              and len(image.get("pixels_sha256", "")) == 64)
        check(f"fonts: {len(files)} bundled faces, text image on the {report.get('platform')} platform "
              f"({image.get('pixels_sha256', '?')[:16]})", ok, f"exit {r.returncode}; {r.stdout[-600:]} {r.stderr[-600:]}")
        if os.name != "nt" and sys.platform != "darwin":
            check("fonts: drawn without a display (offscreen platform)", report.get("platform") == "offscreen",
                  str(report.get("platform")))
        if r.stderr.strip():
            check("fonts: nothing on stderr", False, r.stderr.strip()[-600:])
    if report_path:
        os.makedirs(os.path.dirname(os.path.abspath(report_path)), exist_ok=True)
        with open(report_path, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=1, ensure_ascii=False)
    if reference_path:
        with open(reference_path, encoding="utf-8") as fh:
            reference = json.load(fh)
        theirs = reference.get("image", {}).get("pixels_sha256")
        check(f"fonts: the text image has the same pixels as the reference machine's ({str(theirs)[:16]})",
              image.get("pixels_sha256") == theirs, f"here {image.get('pixels_sha256')}, reference {theirs}")
        check("fonts: the text metrics equal the reference machine's",
              report.get("metrics") == reference.get("metrics") and report.get("metrics"),
              "metrics differ")
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("exe")
    parser.add_argument("--gui", action="store_true", help="also start the interface")
    parser.add_argument("--expect-renderer", choices=["auto", "vulkan", "software"], default="auto",
                        help="renderer the automatic choice must give (auto: either)")
    parser.add_argument("--version", help="version the build must report")
    parser.add_argument("--text-report", help="write the report of `EXE fonts --json` to this file")
    parser.add_argument("--compare-text-report", help="a report of --text-report another machine wrote: the text "
                        "image's pixels and the metrics must be the same")
    parser.add_argument("--screenshot-dir", help="with --gui: keep the screenshots and reports in this folder")
    parser.add_argument("--expect-fallback", help="AppImage: the fallback copies that must be used on this host "
                        "(comma-separated library names), or 'none'")
    parser.add_argument("--gui-exe", help="the program a desktop shortcut starts (Windows: mdbridge-gui.exe): "
                        "check that it opens a project given as its argument")
    args = parser.parse_args(argv)
    exe = os.path.abspath(args.exe)
    global SCREENSHOT_DIR, EXPECT_FALLBACK
    SCREENSHOT_DIR = os.path.abspath(args.screenshot_dir) if args.screenshot_dir else None
    if args.expect_fallback is not None:
        EXPECT_FALLBACK = set() if args.expect_fallback == "none" else \
            {n.strip() for n in args.expect_fallback.split(",") if n.strip()}

    r = run(exe, "--version")
    check(f"--version: {r.stdout.strip()}", r.returncode == 0 and r.stdout.startswith("mdbridge ")
          and (not args.version or f" {args.version} " in r.stdout), r.stderr)

    r = run(exe, "project", "show", os.path.join(PROJECTS, "full.mdatk.toml"))
    check("project show (summary)", r.returncode == 0 and "3 steps, 2 connections" in r.stdout, r.stdout + r.stderr)
    r = run(exe, "project", "show", os.path.join(PROJECTS, "other_format.toml"))
    check("project show rejects another format", r.returncode == 1 and "not an MDATK project" in r.stderr, r.stderr)
    names = sorted(os.listdir(PROJECTS))
    problems = [p for name in names for p in compare.compare([exe], EXPECTED, os.path.join(PROJECTS, name), env())]
    check(f"project show --json reads {len(names)} files as MDATK 0.0.4.1 does", not problems, "\n".join(problems))

    r = run(exe, "graphics", "--json", extra_env=APPRUN_DEBUG)
    try:
        graphics = json.loads(r.stdout)
    except json.JSONDecodeError:
        graphics = {}
    check(f"graphics: {graphics.get('renderer')} ({graphics.get('device') or graphics.get('problem')})",
          r.returncode == 0 and graphics.get("renderer") in ("vulkan", "software"), r.stdout + r.stderr)
    if is_appimage(exe):
        library_check("command line", graphics.get("libraries"), r.stderr)

    fonts_check(exe, args.text_report, args.compare_text_report)

    # MDATK's commands (their output is compared with MDATK's by tests/parity/cli_parity.py)
    data = os.path.join(ROOT, "tests", "data")
    r = run(exe, "info", "--top", os.path.join(data, "topology", "demo_dry.parm7"),
            "--traj", os.path.join(data, "io", "trajectories", "demo_dry.nc"))
    try:
        info = json.loads(r.stdout)
    except json.JSONDecodeError:
        info = {}
    check(f"info: {info.get('n_atoms')} atoms, {info.get('n_frames')} frames of the demo system",
          r.returncode == 0 and info.get("n_atoms") == 773 and info.get("n_frames") == 24, r.stdout + r.stderr)
    r = run(exe, "density", "info", os.path.join(data, "io", "maps", "files", "gzip.mrc.gz"))
    check("density info", r.returncode == 0 and '"grid_xyz"' in r.stdout, r.stdout + r.stderr)
    with tempfile.TemporaryDirectory() as work:
        project = os.path.join(work, "new.mdatk.toml")
        r = run(exe, "init", project, "-s", "notes=smoke test")
        check("init writes a project file", r.returncode == 0 and os.path.exists(project), r.stdout + r.stderr)

    render_check(exe)

    if args.gui:
        with tempfile.TemporaryDirectory() as work:
            screenshot(exe, "software", work)
            auto = screenshot(exe, "auto", work)
            if args.expect_renderer != "auto":
                check(f"automatic choice is {args.expect_renderer}", auto.get("renderer") == args.expect_renderer,
                      f"got {auto.get('renderer')}: {auto.get('fallback_reason')}")
            if args.expect_renderer == "vulkan":
                vulkan = screenshot(exe, "vulkan", work)
                check("the Vulkan view ran", vulkan.get("renderer") == "vulkan", str(vulkan.get("fallback_reason")))
        r = run(exe, "gui", os.path.join(PROJECTS, "split.mdatk.toml"), "--quit-after", "5")
        check("gui runs 5 s and quits (a project without input files)", r.returncode == 0, r.stderr.strip()[-600:])
        gui_edit_check(exe, is_appimage(exe))
    if args.gui_exe:
        gui_exe_check(os.path.abspath(args.gui_exe))

    print(f"\n{'all checks passed' if not FAILED else str(len(FAILED)) + ' check(s) failed'}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
