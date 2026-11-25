"""Smoke test of a finished MDBridge build: AppImage, install tree or Windows mdbridge.exe.

    python packaging/smoke_test.py EXE [--gui] [--expect-renderer auto|vulkan|software] [--version X.Y.Z]
                                       [--text-report OUT.json] [--compare-text-report REF.json]
                                       [--screenshot-dir DIR]

Runs, in order:

* ``EXE --version``
* ``EXE project show`` on the fixtures in tests/data/projects: summaries, rejected files, and the
  JSON of every file against MDATK's reference reading (tests/data/expected)
* ``EXE graphics --json``
* ``EXE fonts --json --render PNG``: the bundled figure fonts load from the program, and a text
  image is drawn without a display (Linux: the "offscreen" platform, also under X);
  ``--text-report`` writes the report, ``--compare-text-report`` checks that the image's pixels and
  the text metrics equal those of another machine's report (CI: the bare Rocky Linux system against
  the build machine)
* ``EXE info`` (demo system), ``EXE density info`` (a map) and ``EXE init`` (a new project file)
* ``--gui``: ``EXE gui PROJECT --screenshot --preview-screenshot`` on a project of the demo system
  (tests/data) with the software view and the automatic choice (``--expect-renderer``: what the
  automatic choice must give; ``vulkan`` also forces a Vulkan run), checking that the 3D view shows
  the project's structure, not the placeholder (the program says so, and a good part of the pixels
  differ from the background), and that the figure preview was drawn; then
  ``EXE gui PROJECT --quit-after 5`` on a project whose input files are missing (the interface
  reports it and keeps running). On a headless Linux machine run it under ``xvfb-run -a``.
  ``--screenshot-dir``: keep the screenshots there (CI uploads them for review).

Uses the standard library only; exits with 1 when a check fails.
"""

from __future__ import annotations

import argparse
import json
import os
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


def env():
    e = dict(os.environ, APPIMAGE_EXTRACT_AND_RUN="1")  # AppImages without FUSE (containers, CI)
    if os.name != "nt":
        e.setdefault("LC_ALL", "C.UTF-8")  # Qt warns on other locales
    return e


def check(label: str, ok: bool, detail: str = "") -> bool:
    print(f"{'ok' if ok else 'FAIL':4s} {label}" + (f"\n     {detail}" if detail and not ok else ""))
    if not ok:
        FAILED.append(label)
    return ok


def run(exe, *args, timeout=180):
    return subprocess.run([exe, *args], capture_output=True, env=env(), encoding="utf-8", errors="replace",
                          timeout=timeout)


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


BUNDLED_FONTS = {"DejaVuSans.ttf", "DejaVuSans-Bold.ttf", "DejaVuSans-Oblique.ttf", "DejaVuSansMono.ttf",
                 "DejaVuSansMono-Bold.ttf"}


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
    args = parser.parse_args(argv)
    exe = os.path.abspath(args.exe)
    global SCREENSHOT_DIR
    SCREENSHOT_DIR = os.path.abspath(args.screenshot_dir) if args.screenshot_dir else None

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

    r = run(exe, "graphics", "--json")
    try:
        graphics = json.loads(r.stdout)
    except json.JSONDecodeError:
        graphics = {}
    check(f"graphics: {graphics.get('renderer')} ({graphics.get('device') or graphics.get('problem')})",
          r.returncode == 0 and graphics.get("renderer") in ("vulkan", "software"), r.stdout + r.stderr)

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

    print(f"\n{'all checks passed' if not FAILED else str(len(FAILED)) + ' check(s) failed'}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
