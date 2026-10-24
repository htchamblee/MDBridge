"""Smoke test of a finished MDBridge build: AppImage, install tree or Windows mdbridge.exe.

    python packaging/smoke_test.py EXE [--gui] [--expect-renderer auto|vulkan|software] [--version X.Y.Z]

Runs, in order:

* ``EXE --version``
* ``EXE project show`` on the fixtures in tests/data/projects: summaries, rejected files, and the
  JSON of every file against MDATK's reference reading (tests/data/expected)
* ``EXE graphics --json``
* ``--gui``: ``EXE gui PROJECT --screenshot`` with the software view and the automatic choice
  (``--expect-renderer``: what the automatic choice must give; ``vulkan`` also forces a Vulkan
  run), checking that each frame really has content; then ``EXE gui PROJECT --quit-after 5``.
  On a headless Linux machine run it under ``xvfb-run -a``.

Uses the standard library only; exits with 1 when a check fails.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tests", "parity"))
import compare  # noqa: E402  (tests/parity/compare.py)

PROJECTS = os.path.join(ROOT, "tests", "data", "projects")
EXPECTED = os.path.join(ROOT, "tests", "data", "expected")
FAILED = []


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


def screenshot(exe, renderer: str, work: str):
    png = os.path.join(work, f"view-{renderer}.png")
    r = run(exe, "gui", os.path.join(PROJECTS, "full.mdatk.toml"), "--renderer", renderer, "--screenshot", png)
    lines = [line for line in r.stdout.splitlines() if line.startswith("{")]
    report = json.loads(lines[-1]) if lines else {}
    ok = r.returncode == 0 and report.get("colours", 0) > 100 and os.path.getsize(png) > 1000
    check(f"gui --renderer {renderer}: {report.get('description', '?')}", ok,
          f"exit {r.returncode}; {report}; {r.stderr.strip()[-600:]}")
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("exe")
    parser.add_argument("--gui", action="store_true", help="also start the interface")
    parser.add_argument("--expect-renderer", choices=["auto", "vulkan", "software"], default="auto",
                        help="renderer the automatic choice must give (auto: either)")
    parser.add_argument("--version", help="version the build must report")
    args = parser.parse_args(argv)
    exe = os.path.abspath(args.exe)

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
        check("gui runs 5 s and quits", r.returncode == 0, r.stderr.strip()[-600:])

    print(f"\n{'all checks passed' if not FAILED else str(len(FAILED)) + ' check(s) failed'}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
