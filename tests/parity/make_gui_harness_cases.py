#!/usr/bin/env python3
"""The GUI test harness' own case: open a fixture project, save it, compare (tests/gui).

    python3 -I -B make_gui_harness_cases.py MDATK_CHECKOUT DATA_DIR

Reads DATA_DIR/projects/default.mdatk.toml (make_projects.py) and writes to DATA_DIR/gui/harness:

* ``open_save_default.json``: the action script (format: tests/data/README.md, "gui/"),
* ``open_save_default.mdatk.toml``: the file MDATK v0.0.4.1 writes when that project is opened and
  saved unchanged.

What MDATK's window does on *Open* and *Save* without edits (``open_path``, ``save_project``,
mdatk/gui/app.py) is ``Project.load(path)`` and ``project.save(path)`` (``collect()`` stores the
widgets' values, which for this project are the loaded ones); this generator calls those two, so it
needs no GTK, with MDBridge's window improvements of improvements_gui.py that apply to them: the
open keeps the paths as the file has them (AP-06, D8: ``improvements_gui.load_as_written``, where
MDATK's ``Project.load`` makes them absolute), so the expected file is what MDBridge's window
writes. The action-script cases that drive MDATK's widgets come from make_gui_cases.py.

The clock is fixed as in make_project_file_golden.py (2026-10-15T14:32:05, one second later on
every call). The project's folder is a temporary one; its absolute path, should the saved file
contain it, is replaced by ``{case_dir}``, the placeholder the harness puts in for its own folder.
MDATK's recent-projects list goes to a temporary XDG_CONFIG_HOME.
"""

from __future__ import annotations

import datetime as real_dt
import json
import os
import shutil
import sys
import tempfile
import types

sys.dont_write_bytecode = True  # never write __pycache__ into the MDATK checkout
if len(sys.argv) != 3:
    sys.exit("usage: make_gui_harness_cases.py MDATK_CHECKOUT DATA_DIR")
MDATK = sys.argv[1]
DATA = os.path.abspath(sys.argv[2])
sys.path.insert(0, MDATK)

from mdatk import projectfile as PF  # noqa: E402
from mdatk.config import Project  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import improvements_gui as IG  # noqa: E402

CLOCK = real_dt.datetime(2026, 10, 15, 14, 32, 5)
PLACEHOLDER = "{case_dir}"


class _Clock:
    t = CLOCK


class FixedDateTime(real_dt.datetime):
    @classmethod
    def now(cls, tz=None):  # noqa: D401 - mimic datetime.now
        value = _Clock.t
        _Clock.t = _Clock.t + real_dt.timedelta(seconds=1)
        return cls(value.year, value.month, value.day, value.hour, value.minute, value.second, 123456)


PF._dt = types.SimpleNamespace(datetime=FixedDateTime)  # _now() looks _dt up at call time


def write_text(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def open_save(out: str, name: str, fixture: str, description: str) -> None:
    _Clock.t = CLOCK
    base = os.path.basename(fixture)
    with tempfile.TemporaryDirectory() as tmp:
        work = os.path.join(tmp, "case")
        os.makedirs(work)
        path = os.path.join(work, base)
        shutil.copyfile(os.path.join(DATA, fixture), path)
        IG.load_as_written(Project, PF, path).save(path)  # the window's open_path, then save_project
        with open(path, encoding="utf-8") as f:
            text = f.read()
    expected = f"{name}.mdatk.toml"
    write_text(os.path.join(out, expected), text.replace(work, PLACEHOLDER))
    case = {
        "format": "mdbridge-gui-case",
        "version": 1,
        "description": description,
        "source": "MDATK v0.0.4.1: Project.load + Project.save, what the window's Open and Save do "
                  "without edits, with MDBridge's window improvements (paths as written, AP-06; "
                  "tests/parity/make_gui_harness_cases.py, improvements_gui.py)",
        "clock": CLOCK.isoformat(),
        "files": {base: fixture},
        "actions": [
            {"action": "open", "path": base},
            {"action": "wait"},
            {"action": "screenshot", "name": "opened"},
            {"action": "save"},
        ],
        "expected": {base: f"gui/harness/{expected}"},
    }
    write_text(os.path.join(out, f"{name}.json"), json.dumps(case, indent=2, ensure_ascii=False) + "\n")


def main() -> None:
    out = os.path.join(DATA, "gui", "harness")
    os.makedirs(out, exist_ok=True)
    config = tempfile.mkdtemp()
    os.environ["XDG_CONFIG_HOME"] = config  # add_recent of a real window would write there
    try:
        open_save(out, "open_save_default", "projects/default.mdatk.toml",
                  "Open the default project (what mdatk init writes) and save it unchanged")
    finally:
        shutil.rmtree(config, ignore_errors=True)


if __name__ == "__main__":
    main()
