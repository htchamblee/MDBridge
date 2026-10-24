"""Compare what MDBridge reads from project files with what MDATK reads (reference.py's output).

    python tests/parity/compare.py --mdbridge EXE --expected DIR FILE...

Runs ``EXE project show --json FILE`` for every FILE and compares its "settings" and "log" with
``DIR/<name>.settings.json`` / ``.log.json`` byte for byte (after writing them the way
reference.py does), or checks that MDBridge rejects the file when MDATK did (``<name>.error.txt``).
Standard library only; exits with 1 on any difference.
"""

from __future__ import annotations

import argparse
import difflib
import json
import os
import subprocess
import sys
from typing import List, Sequence


def dumps(data) -> str:
    return json.dumps(data, indent=2, ensure_ascii=False) + "\n"


def read(path: str) -> str:
    with open(path, encoding="utf-8", newline="") as fh:
        return fh.read()


def diff(expected: str, actual: str, label: str) -> str:
    lines = difflib.unified_diff(expected.splitlines(), actual.splitlines(), f"MDATK {label}", f"MDBridge {label}",
                                 lineterm="", n=2)
    return "\n".join(list(lines)[:40])


def compare(exe: Sequence[str], expected_dir: str, path: str, env=None) -> List[str]:
    """Problems found for one file (empty: MDBridge reads it as MDATK does)."""
    name = os.path.basename(path)
    run = subprocess.run([*exe, "project", "show", "--json", path], capture_output=True, env=env,
                         encoding="utf-8", errors="replace", timeout=120)
    error_file = os.path.join(expected_dir, name + ".error.txt")
    if os.path.exists(error_file):
        if run.returncode == 0:
            return [f"{name}: MDATK rejects it ({read(error_file).strip()}), MDBridge read it"]
        return []
    if run.returncode != 0:
        return [f"{name}: MDBridge failed (exit {run.returncode}): {run.stderr.strip()}"]
    data = json.loads(run.stdout)
    problems = []
    for part in ("settings", "log"):
        reference = read(os.path.join(expected_dir, f"{name}.{part}.json"))
        actual = dumps(data[part])
        if actual != reference:
            problems.append(f"{name}: {part} differs\n{diff(reference, actual, part)}")
    return problems


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mdbridge", required=True, help="the mdbridge program (or AppImage)")
    parser.add_argument("--expected", required=True, help="folder written by reference.py")
    parser.add_argument("files", nargs="+")
    args = parser.parse_args(argv)
    env = dict(os.environ, APPIMAGE_EXTRACT_AND_RUN="1")
    failed = 0
    for path in args.files:
        problems = compare([args.mdbridge], args.expected, path, env)
        print(f"{'ok' if not problems else 'FAIL':4s} {path}")
        for problem in problems:
            print("     " + problem.replace("\n", "\n     "))
        failed += bool(problems)
    print(f"{len(args.files) - failed} of {len(args.files)} files read as MDATK reads them")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
