"""Write MDATK project files with MDATK itself (the reference implementation).

    python tests/parity/make_projects.py OUT_DIR

MDATK (v0.0.4.1) must be importable: ``pip install`` a checkout of the tag, or set PYTHONPATH to it.
Writes, into OUT_DIR:

* ``default.mdatk.toml``    a project with default settings (what ``mdatk init`` writes)
* ``full.mdatk.toml``       every kind of entry: notes, relative input paths, unset values, special
                            floats, escapes and non-ASCII text, planner steps and connections, and a
                            step log with ok / failed / cancelled entries
* ``split.mdatk.toml``      a project whose step log moved to ``split.log.toml``
* ``MDATK_params_log.toml`` the parameter log of a run (opens like a project)

``reference.py`` then records what MDATK reads back from them.

The files carry the time they were written (``saved``, ``created``, log entry times and durations).
With ``--fixed-clock`` MDATK's clock is replaced by the one of the run that wrote the committed
files in tests/data/projects (2026-10-24T06:39:22 for the first time stamp, 06:39:23 for the
others, durations 0.0), so the output equals those files byte for byte
(tests/parity/regenerate_all.py uses it).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys
import types

import mdatk
from mdatk import projectfile as PF
from mdatk.config import ExternalSeries, PanelSpec, Project, RepSpec, ResidueGroup
from mdatk.progress import Cancelled
from mdatk.workflows.plan import add_node, connect


def default_project(out: str) -> None:
    Project().save(os.path.join(out, "default.mdatk.toml"))


def full_project(out: str) -> None:
    p = Project()
    p.notes = ('Aim: compare three "quoted" runs of the Mg²⁺-bound aptamer (5 Å cutoff).\n'
               "Tabs\there, a backslash \\ and a control character \x01.\n"
               "Unicode: ψ-angle, 🧬, naïve café.")
    p.inputs.topology = "system.parm7"
    p.inputs.trajectories = ["md1.nc", "md2.nc"]
    p.inputs.stop = None
    p.inputs.step = 2
    p.inputs.dt_ps = 2.5
    p.inputs.external.append(ExternalSeries(path="rmsd.dat", label="cpptraj RMSD", columns=["2", "rmsd"]))
    p.groups = [ResidueGroup("Aptamer", "A:1-40, 45", "", "#2a78d6"),
                ResidueGroup("Target", "B:1-30", "", "#eb6834"),
                ResidueGroup("Mg", "", "resname MG", "#4caf50")]
    p.analyses.energy.target_value = None
    p.visual.representations = [RepSpec("cartoon", ["Aptamer"], "nucleotide", 1.2, 0.85, True),
                                RepSpec("vdw", ["Mg"], "element", 0.5, 1.0)]
    p.layout.panels.append(PanelSpec("timeseries", ["rmsd/Complex", "external/cpptraj RMSD/rmsd"], "RMSD",
                                     {"smooth": 5, "progressive": True, "target": 1e-07, "vmax": 1e16,
                                      "vmin": 100000.0, "low": -0.0, "nan": float("nan"),
                                      "inf": float("inf"), "custom": [1, 2.5, "x", True]}))
    p.layout.panels.append(PanelSpec("contact_map", ["contacts/map"], "Contacts", {"transpose": True}))
    wf = p.workflow
    wf.project_name = "apt1"
    wf.stamp = "20261015-1432"
    build = add_node(wf, "amber-build", x=60, y=80, settings={"mg_count": 4})
    build.files["structure"] = "model.pdb"
    md = add_node(wf, "amber-md", x=320, y=80, title="Production MD")
    add_node(wf, "note", x=60, y=300, title="Why: compare with the cryo-EM map")
    connect(wf, build.id, "topology", md.id, "topology")
    connect(wf, build.id, "coordinates", md.id, "coordinates")
    path = os.path.join(out, "full.mdatk.toml")
    p.save(path)

    with PF.logged(path, "02", "tleap build", ["amber", "ol3"]) as entry:
        entry.command = "bash 02_amber-build.sh"
        entry.folder = "apt1_workflow/02_amber-build"
        entry.files = ["apt1_02_amber-build_20261015-1432.parm7", "café.rst7"]
        entry.summary = "4 Mg²⁺, 21 Na⁺\nbox 82.1 Å"
        entry.output = "x" * (PF.LOG_OUTPUT_CHARS + 50)
    try:
        with PF.logged(path, "03", "md"):
            raise RuntimeError("pmemd crashed: \"NaN\" in step 1200")
    except RuntimeError:
        pass
    try:
        with PF.logged(path, "render", "video"):
            raise Cancelled()
    except Cancelled:
        pass


def split_project(out: str) -> None:
    path = os.path.join(out, "split.mdatk.toml")
    p = Project()
    p.notes = "Long project: the step log lives in split.log.toml."
    p.save(path)
    with PF.logged(path, "pdb clean", "L-RNA naming"):
        pass
    limit = PF.LOG_SPLIT_BYTES
    PF.LOG_SPLIT_BYTES = 1000  # the project file is already larger: the next entry moves the log out
    try:
        with PF.logged(path, "04", "analysis") as entry:
            entry.summary = "RMSD, contacts, H-bonds"
    finally:
        PF.LOG_SPLIT_BYTES = limit
    p.notes = "Long project: the step log lives in split.log.toml (saved again after the split)."
    p.save(path)  # keeps log_file
    with PF.logged(path, "05", "video"):
        pass


def params_log(out: str) -> None:
    p = Project()
    p.inputs.topology = "/data/run1/system.parm7"
    p.inputs.trajectories = ["/data/run1/md1.nc"]
    p.output.out_dir = "/data/run1/MDATK_output"
    run = {"timestamp": "2026-10-15T14:32:05", "mdatk_version": mdatk.__version__,
           "inputs": {"topology": p.inputs.topology, "frames": 120},
           "outputs": ["MDATK_video.mp4"], "citations": ["mdanalysis2016"]}
    PF.write_params_log(os.path.join(out, "MDATK_params_log.toml"), p, run)


def fix_clock() -> None:
    """MDATK's clock as it ran when the committed files were written (see the module docstring)."""
    stamps = iter([_dt.datetime(2026, 10, 24, 6, 39, 22)])
    later = _dt.datetime(2026, 10, 24, 6, 39, 23)

    class FixedDateTime(_dt.datetime):
        @classmethod
        def now(cls, tz=None):
            value = next(stamps, later)
            return cls(value.year, value.month, value.day, value.hour, value.minute, value.second, 123456)

    PF._dt = types.SimpleNamespace(datetime=FixedDateTime)  # _now() looks _dt up at call time
    PF.time = types.SimpleNamespace(monotonic=lambda: 100.0)  # logged(): every duration 0.0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("out_dir")
    parser.add_argument("--fixed-clock", action="store_true",
                        help="the clock of the committed files (reproducible output)")
    args = parser.parse_args(argv)
    if args.fixed_clock:
        fix_clock()
    os.makedirs(args.out_dir, exist_ok=True)
    for make in (default_project, full_project, split_project, params_log):
        make(args.out_dir)
    for name in sorted(os.listdir(args.out_dir)):
        if name.endswith(".toml"):
            print(os.path.join(args.out_dir, name))
    print(f"(written by MDATK {mdatk.__version__})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
