"""Project files written by a later MDATK (0.1.1), for MDBridge's round-trip tests (owner decision D9 of
milestone 5; docs/DIFFERENCES.md PRJ-09, PRJ-10).

    python3 -I -B tests/parity/make_mdatk_latest_projects.py MDATK_LATEST_CHECKOUT OUT_DIR

MDATK_LATEST_CHECKOUT is MDATK at version 0.1.1 (the commit regenerate_all.py names); nothing is written
into it. Standard library only (MDATK's project files and planner import no numpy). The clock is fixed,
so the files are the same on every run. Writes into OUT_DIR (tests/data/project-file/mdatk-0.1.1):

* ``plan.mdatk.toml``   a project with every key 0.1.x added to the plan (workflow.supersede_runs,
                        parallel, cpu_cores, gpus, last_node; per step adopted and folder_tag), steps of
                        tool types and options 0.0.4.1 does not have (namd-mdff, user-script, the cluster
                        job options of amber-md), notes that hold the step-log marker line (pasted from
                        another project file), and a step log with the statuses waiting and stopped and the
                        step "workflow supersede", appended by append_log
* ``sidecar.mdatk.toml`` + ``sidecar.log.toml``  a project whose log moved to a separate file
* ``read.json``         per project: MDATK 0.1.1's read() (the settings, tagged as in make_model_cases.py)
                        and read_log() (the step log entries), and the byte offset where its _log_start
                        puts the log (-1: none)
"""

from __future__ import annotations

import datetime as dt
import json
import os
import sys

if len(sys.argv) != 3:
    sys.exit("usage: python3 -I -B tests/parity/make_mdatk_latest_projects.py MDATK_LATEST_CHECKOUT OUT_DIR")
MDATK, OUT = (os.path.abspath(a) for a in sys.argv[1:3])
sys.dont_write_bytecode = True
sys.path.insert(0, MDATK)

import mdatk  # noqa: E402
from mdatk import config as C  # noqa: E402
from mdatk import projectfile as PF  # noqa: E402
from mdatk.workflows import registry  # noqa: E402

if mdatk.__version__ != "0.1.1":
    sys.exit(f"{MDATK} is MDATK {mdatk.__version__}, not 0.1.1")

SAVED = dt.datetime(2026, 11, 25, 9, 30, 0)
PF._now = lambda: SAVED  # project_text's "saved"

os.makedirs(OUT, exist_ok=True)


def enc(v):
    """The tagged values of make_model_cases.py (tests/core/test_support.hpp)."""
    if v is None:
        return {"t": "none"}
    if isinstance(v, bool):
        return {"t": "bool", "v": v}
    if isinstance(v, int):
        return {"t": "int", "v": str(v)}
    if isinstance(v, float):
        return {"t": "float", "hex": v.hex(), "repr": repr(v)}
    if isinstance(v, str):
        return {"t": "str", "utf8": v.encode("utf-8").hex(), "v": v}
    if isinstance(v, dt.datetime):
        return {"t": "datetime", "iso": v.isoformat()}
    if isinstance(v, dt.date):
        return {"t": "date", "iso": v.isoformat()}
    if isinstance(v, dt.time):
        return {"t": "time", "iso": v.isoformat()}
    if isinstance(v, (list, tuple)):
        return {"t": "array", "v": [enc(x) for x in v]}
    if isinstance(v, dict):
        return {"t": "table", "v": [[k, enc(x)] for k, x in v.items()]}
    raise TypeError(type(v))


def node(nid, tool, step, x, **extra):
    n = C.PlanNode(id=nid, tool=tool, step=step, x=x, y=60.0, settings=registry.get_tool(tool).defaults())
    for key, value in extra.items():
        if key == "settings":
            n.settings.update(value)
        else:
            setattr(n, key, value)
    return n


def plan_project():
    p = C.Project()
    p.notes = ("Copied from another project file:\n"
               f"{PF.LOG_MARKER}\n"
               "[[log]]\nstep = \"render\"\n")
    p.inputs.topology = "inputs/system.parm7"
    p.inputs.trajectories = ["inputs/run1.nc"]
    w = p.workflow
    w.project_name = "apt1"
    w.stamp = "20261120-1015"
    w.supersede_runs = False
    w.parallel = 3
    w.cpu_cores = 8
    w.gpus = 1
    w.last_node = 7
    w.nodes = [
        node("n1", "input-file", 1, 40.0, settings={"path": "inputs/model.pdb", "kind": "pdb"}),
        node("n2", "amber-md", 2, 260.0, adopted="20261101-0900", folder_tag="model12",
             settings={"cluster": "", "job_header": "#SBATCH --time=24:00:00\n#SBATCH -p gpu",
                       "job_submit": "sbatch"}),
        node("n7", "namd-mdff", 3, 480.0),
        node("n5", "user-script", 4, 700.0, settings={"script": "echo done\n" + PF.LOG_MARKER + "\n"}),
    ]
    w.edges = [C.PlanEdge("n1", "file", "n2", "structure")]
    w.scratch = [node("s1", "user-script", 1, 40.0)]
    return p


def entry(step, status, minute, **kw):
    return PF.LogEntry(step=step, status=status, time=dt.datetime(2026, 11, 25, 10, minute, 0), **kw)


ENTRIES = [
    entry("workflow write", "ok", 0, title="Plan files written", duration_s=0.42,
          files=["apt1_workflow/02_amber-build_model12/apt1_02_amber-md_20261120-1015.sh"]),
    entry("02", "waiting", 5, title="Amber MD (n2)", duration_s=1.5, command="bash apt1_02_amber-md.sh",
          folder="apt1_workflow/02_amber-build_model12",
          summary="waiting for the separately run script (exit status 75)"),
    entry("03", "stopped", 9, title="NAMD MDFF (n7)", duration_s=12.0, summary="stopped by the user",
          output="line 1\nline 2\n"),
    entry("workflow supersede", "ok", 12, title="Earlier run of step 02 moved aside", duration_s=0.05,
          files=["apt1_workflow/02_amber-build_model12/superseded/20261120-1015_20261125-101200"]),
    entry("render", "failed", 20, duration_s=3.25, summary="ValueError: no frames", citations=["mdanalysis2016"]),
]


def write(name, project, entries, split=False):
    path = os.path.join(OUT, name)
    for stale in (path, path[:-len(PF.PROJECT_SUFFIX)] + ".log.toml"):
        if os.path.exists(stale):
            os.remove(stale)
    PF.save(path, project)
    for k, e in enumerate(entries):
        PF.append_log(path, e)
        if split and k == 1:
            PF._split_log(path)
    lock = os.path.join(OUT, f".{name}.lock")
    if os.path.exists(lock):
        os.remove(lock)
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    return {"settings": enc(PF.read(path)), "log": enc(PF.read_log(path)), "log_start": PF._log_start(text),
            "log_file": PF._existing_log(path)[1]}


def main() -> int:
    plan = plan_project()
    out = {"how": "MDATK 0.1.1: projectfile.read() and read_log() of the files it wrote here (save, "
                  "append_log, _split_log; clock fixed), and _log_start of the project file's text",
           "mdatk_version": mdatk.__version__,
           "files": {"plan.mdatk.toml": write("plan.mdatk.toml", plan, ENTRIES),
                     "sidecar.mdatk.toml": write("sidecar.mdatk.toml", plan, ENTRIES[:3], split=True)}}
    with open(os.path.join(OUT, "read.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print("wrote", ", ".join(sorted(os.listdir(OUT))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
