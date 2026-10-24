#!/usr/bin/env python3
"""Golden project files, step logs, parameter log and user settings written by MDATK v0.0.4.1.

Run (from anywhere; MDATK is imported from its checkout, nothing is written there):

    python3 -I -B make_project_file_golden.py MDATK_CHECKOUT OUT_DIR

Every output lands in OUT_DIR.  Time stamps and durations are made deterministic by
replacing ``mdatk.projectfile._dt`` (used by ``_now()``) and ``mdatk.projectfile.time`` (used by
``logged()``) with fakes, so the files can be compared byte for byte:

* ``_now()`` returns 2026-10-15T14:32:05, then one second later on every call.
* ``time.monotonic()`` returns 100.0, 101.23456, 102.0, 102.004, 103.0, 105.555, ... (pairs: start / end of
  one ``logged()`` block; see MONO below), giving the durations listed there.

Absolute paths of the temporary work folder are replaced by ``/home/user/projects`` in the saved copies
(``settings.toml``, ``manifest.txt``).
"""

from __future__ import annotations

import datetime as real_dt
import os
import shutil
import sys
import types

sys.dont_write_bytecode = True  # never write __pycache__ into the MDATK checkout
if len(sys.argv) != 3:
    sys.exit("usage: make_project_file_golden.py MDATK_CHECKOUT OUT_DIR")
MDATK = sys.argv[1]
sys.path.insert(0, MDATK)

from mdatk import projectfile as PF  # noqa: E402
from mdatk import tomlio  # noqa: E402
from mdatk.config import ExternalSeries, PanelSpec, Project, RepSpec, ResidueGroup  # noqa: E402
from mdatk.progress import Cancelled  # noqa: E402
from mdatk.workflows.plan import add_node, connect  # noqa: E402

HERE = os.path.abspath(sys.argv[2])
os.makedirs(HERE, exist_ok=True)
WORK = os.path.join(HERE, "work")
FAKE_ROOT = "/home/user/projects"


# ---------------------------------------------------------------- deterministic clock
class _Clock:
    t = real_dt.datetime(2026, 10, 15, 14, 32, 5)


class FixedDateTime(real_dt.datetime):
    @classmethod
    def now(cls, tz=None):  # noqa: D401 - mimic datetime.now
        value = _Clock.t
        _Clock.t = _Clock.t + real_dt.timedelta(seconds=1)
        return cls(value.year, value.month, value.day, value.hour, value.minute, value.second, 123456)


PF._dt = types.SimpleNamespace(datetime=FixedDateTime)  # _now() looks _dt up at call time

# start / end pairs of each logged() block (durations: 1.23456, 0.004, 2.555, 0.125, 0.135, 7.0, ...)
MONO = [100.0, 101.23456, 102.0, 102.004, 103.0, 105.555, 106.0, 106.125, 107.0, 107.135, 108.0, 115.0,
        116.0, 116.5, 117.0, 117.995, 118.0, 118.005, 119.0, 119.015, 120.0, 120.25, 121.0, 121.5,
        122.0, 122.75, 123.0, 123.1]
_mono = iter(MONO)
PF.time = types.SimpleNamespace(monotonic=lambda: next(_mono))


def fresh(name):
    path = os.path.join(WORK, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def keep(src, dst_name):
    shutil.copyfile(src, os.path.join(HERE, dst_name))


# ---------------------------------------------------------------- 1. default project
def default_project():
    path = fresh("default/default.mdatk.toml")
    Project().save(path)
    keep(path, "default.mdatk.toml")


# ---------------------------------------------------------------- 2. rich project
def rich_project() -> Project:
    p = Project()
    p.notes = ('Aim: compare three "quoted" runs\n'
               "Path C:\\runs\\apt1\twith a tab\n"
               "\n"
               "Unicode: 4.5 \u00c5 cut-off, Mg\u00b2\u207a ions, \u2264 2.6 \u00c5\n"
               "triple quotes \"\"\" inside, carriage\rreturn, bell\x07 and DEL\x7f\n")
    p.inputs.topology = "model.parm7"
    p.inputs.trajectories = ["md/prod1.nc", "md/prod2.nc", "md/prod3.nc"]
    p.inputs.stop = None  # -> "# stop =   (not set)"
    p.inputs.dt_ps = 2.5
    p.inputs.external = [ExternalSeries(path="rmsd.dat", label="cpptraj RMSD", x_column="0", columns=["2", "3"]),
                         ExternalSeries(path="energy.csv")]
    p.cleaning.enabled = True
    p.cleaning.stop = 5000
    p.groups = [ResidueGroup("Aptamer", "A:1-40", "", "#2a78d6"),
                ResidueGroup("Target", "B:1-12", "", "#eb6834"),
                ResidueGroup("Loop", "", "resid 15-20 and not name H*", "#7f7f7f")]
    p.analyses.energy.target_value = -150.5
    p.analyses.rmsd.reference_frame = 10
    p.analyses.contacts.cutoff = 4.0
    p.visual.representations = [RepSpec("cartoon", ["Aptamer"]),
                                RepSpec("licorice", ["Target", "Loop"], "element", 0.5, 0.8)]
    p.visual.rotation = [-90.0, 0.0, 12.5]
    p.visual.center_on = "B:12"
    p.visual.title = "A long figure title that is longer than ninety-six characters, to show that strings are never wrapped"
    p.density.map_path = "maps/emd_1234.mrc"
    p.density.transform = [0.999123, -0.012, 0.0405, 0.0119, 0.999927, 0.001, -0.0405, -0.000516, 0.99918,
                           12.345678, -3.5, 101.0]
    p.density.anchor = [[round(1.5 * k + 0.123, 3), round(-2.25 * k, 3), round(0.1 * k * k, 3)] for k in range(8)]
    p.density.fit_cc = 0.8123
    p.density.fit_info = "Global search: cc 0.8123 at 30 deg"
    p.layout.preset = "grid_2x3"
    p.layout.panels[1].options = {"smooth": 5, "progressive": True, "custom": 1, "target": None,
                                  "nested": {"a": 1, "b": None}}
    p.layout.panels.append(PanelSpec("heatmap", ["contacts/per_residue_a"], "Per-residue contacts",
                                     {"vmin": 0.0, "vmax": 0, "max_ticks": 20}))
    p.layout.panels.append(PanelSpec("secstruct", [], "", {"residues": "A:1-20, B:5-30", "source": "dotbracket",
                                                           "dotbracket": "((((....))))&((((....))))"}))
    p.output.out_dir = "results"
    p.output.charts.formats = ["png", "svg", "pdf"]
    p.output.charts.series = [f"rmsd/Target{k}" for k in range(12)]  # long list -> one element per line
    p.output.pdb.frames = "0, 10, 20-50:10, last"
    wf = p.workflow
    wf.project_name = "apt1"
    wf.stamp = "20261015-1432"
    build = add_node(wf, "amber-build", 40.0, 40.0, settings={"mg_count": 4})
    build.files["structure"] = "model.pdb"
    build.files["frcmod"] = "a.frcmod;b.frcmod"
    md = add_node(wf, "amber-md", 260.0, 40.0, title="Equilibration and production")
    md.settings["stages"] = md.settings["stages"][:2]  # two records keep the golden short
    add_node(wf, "note", 40.0, 220.0, settings={"text": "Remember: check the Mg2+ count.\nSecond line."})
    connect(wf, build.id, "topology", md.id, "topology")
    connect(wf, build.id, "coordinates", md.id, "coordinates")
    scratch = add_node(wf, "vmd", 0.0, 0.0)  # build the record with the planner's defaults ...
    wf.nodes.remove(scratch)  # ... and keep it as a stand-alone wizard setting instead
    wf.scratch.append(scratch)
    return p


def rich_and_logs():
    path = fresh("rich/rich.mdatk.toml")
    p = rich_project()
    p.save(path)
    keep(path, "rich.mdatk.toml")

    # step log: the first entry carries the comments, the others do not
    with PF.logged(path, "02", "tleap build", ["case-2023-amber", "key1"]) as entry:  # 1.23456 s
        entry.command = "bash apt1_02_tleap.sh"
        entry.folder = "/home/user/projects/apt1_workflow/02_amber-build"
        entry.files = ["apt1_02_tleap_20261015-1432.parm7", "apt1_02_tleap_20261015-1432.rst7"]
        entry.summary = "exit status 0"
        entry.output = "Welcome to LEaP!\nLoading parameters\n\tdone\n"
    try:
        with PF.logged(path, "03", "md"):  # 0.004 s
            raise RuntimeError("pmemd crashed")
    except RuntimeError:
        pass
    try:
        with PF.logged(path, "04", "slow run"):  # 2.555 s
            raise Cancelled("stopped by the user")
    except Cancelled:
        pass
    try:
        with PF.logged(path, "render", "Analyses and outputs"):  # 0.125 s
            raise KeyboardInterrupt
    except KeyboardInterrupt:
        pass
    try:
        with PF.logged(path, "pdb clean"):  # 0.135 s, no title
            raise ValueError()
    except ValueError:
        pass
    try:
        with PF.logged(path, "05", "summary set before the failure") as entry:  # 7.0 s
            entry.summary = "exit status 3"
            raise RuntimeError("ignored: the summary is already set")
    except RuntimeError:
        pass
    with PF.logged(path, "pdb clean", "PDB clean-up of caf\udce9.pdb") as entry:  # 0.5 s
        entry.files = ["caf\udce9.pdb", "clean/caf\udce9_clean.pdb"]
        entry.summary = "renamed 3 residues\nrenumbered chain B\n"
    keep(path, "rich.with-log.mdatk.toml")

    q = Project.from_dict(PF.read(path))  # no path resolution (Project.load would make paths absolute)
    q.notes = "edited notes"
    q.save(path)  # saving keeps the log
    keep(path, "rich.resaved.mdatk.toml")


def truncation_log():
    path = fresh("truncation/t.mdatk.toml")
    Project().save(path)
    with PF.logged(path, "06", "long output") as entry:  # 0.995 s
        entry.output = "".join(f"line {k:05d}\n" for k in range(2000))  # 22000 characters
    keep(path, "truncation.mdatk.toml")


def split_log():
    path = fresh("split/p.mdatk.toml")
    Project().save(path)
    with PF.logged(path, "02", "tleap build"):  # 0.005 s
        pass
    with PF.logged(path, "03", "md"):  # 0.015 s
        pass
    keep(path, "split.before.mdatk.toml")
    old = PF.LOG_SPLIT_BYTES
    PF.LOG_SPLIT_BYTES = 1000
    try:
        with PF.logged(path, "04", "analysis"):  # 0.25 s
            pass
    finally:
        PF.LOG_SPLIT_BYTES = old
    keep(path, "split.after.mdatk.toml")
    keep(os.path.join(os.path.dirname(path), "p.log.toml"), "split.after.log.toml")
    Project.from_dict(PF.read(path)).save(path)  # tests use Project.load (absolute paths); same layout
    with PF.logged(path, "05", "video"):  # 0.5 s
        pass
    keep(path, "split.resaved.mdatk.toml")
    keep(os.path.join(os.path.dirname(path), "p.log.toml"), "split.resaved.log.toml")

    # a project file over the limit that has no log yet: the side file gets marker + entry without comments
    path2 = fresh("split2/q.mdatk.toml")
    Project().save(path2)
    PF.LOG_SPLIT_BYTES = 1000
    try:
        with PF.logged(path2, "01", "first entry of a big project"):  # 0.75 s
            pass
    finally:
        PF.LOG_SPLIT_BYTES = old
    keep(path2, "split-nolog.after.mdatk.toml")
    keep(os.path.join(os.path.dirname(path2), "q.log.toml"), "split-nolog.after.log.toml")


def params_log():
    path = fresh("params/MDATK_params_log.toml")
    run = {
        "project": "Molecular Dynamic Analysis Tool Kit (MDATK)",
        "environment": {"mdatk": "0.0.4.1", "python": "3.13.17", "MDAnalysis": "2.9.0", "numpy": "2.3.1",
                        "platform": "Linux-6.18-x86_64", "cpptraj": None},
        "inputs": {"topology": None, "trajectories": None, "cleaning": None},
        "cache_dir": "/home/user/projects/results/MDATK_cache",
        "analysis_keys": {"rmsd": "0123456789abcdef", "contacts": "fedcba9876543210"},
        "timestamp": "2026-10-15_143205",
        "outputs": {"video": "", "charts": [], "pdbs": [], "params_log": None},
        "citation_keys": ["mdanalysis-2016", "mdanalysis-2011", "numpy-2020", "matplotlib-2007", "mdatk"],
        "citations": ["Gowers, R. J. et al. MDAnalysis: A Python Package for the Rapid Analysis of Molecular "
                      "Dynamics Simulations. Proc. 15th Python in Science Conf. (2016).",
                      "Harris, C. R. et al. Array programming with NumPy. Nature 585, 357-362 (2020)."],
    }
    PF.write_params_log(path, rich_project(), run)
    keep(path, "params_log.toml")
    # a params log is never logged into and cannot be saved over
    before = open(path, encoding="utf-8").read()
    with PF.logged(path, "render"):  # 0.75 s: no entry is written
        pass
    assert open(path, encoding="utf-8").read() == before
    try:
        Project().save(path)
    except ValueError as exc:
        with open(os.path.join(HERE, "messages.txt"), "a", encoding="utf-8") as fh:
            fh.write("save over a params log: " + str(exc).replace(WORK, "<work>") + "\n")


def recent():
    cfg = os.path.join(WORK, "config")
    os.environ["XDG_CONFIG_HOME"] = cfg
    paths = []
    for k in range(PF.RECENT_MAX + 2):
        path = fresh(f"recent/p{k}.mdatk.toml")
        Project().save(path)
        PF.add_recent(path)
        paths.append(path)
    text = open(os.path.join(cfg, "mdatk", "settings.toml"), encoding="utf-8").read()
    with open(os.path.join(HERE, "settings.toml"), "w", encoding="utf-8") as fh:
        fh.write(text.replace(WORK, FAKE_ROOT))
    assert PF.recent_projects()[0] == paths[-1]


def messages():
    out = []
    try:
        Project().save(os.path.join(WORK, "preset.json"))
    except ValueError as exc:
        out.append("save .json: " + str(exc).replace(WORK, "<work>"))
    try:
        PF.read(os.path.join(WORK, "preset.JSON"))
    except ValueError as exc:
        out.append("read .JSON: " + str(exc).replace(WORK, "<work>"))
    try:
        PF.append_log(os.path.join(WORK, "missing.mdatk.toml"), PF.LogEntry(step="x"))
    except FileNotFoundError as exc:
        out.append("append_log missing: " + repr(exc).replace(WORK, "<work>"))
    with open(os.path.join(HERE, "messages.txt"), "a", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")


def manifest():
    lines = []
    for root, dirs, files in os.walk(WORK):
        dirs.sort()
        for name in sorted(files):
            full = os.path.join(root, name)
            lines.append(f"{os.path.relpath(full, WORK)}\t{os.path.getsize(full)}")
    with open(os.path.join(HERE, "manifest.txt"), "w", encoding="utf-8") as fh:
        fh.write("# files left in the work folder (relative path, size in bytes): lock files stay, .tmp files do not\n")
        fh.write("\n".join(lines) + "\n")


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    if os.path.exists(os.path.join(HERE, "messages.txt")):
        os.remove(os.path.join(HERE, "messages.txt"))
    default_project()
    rich_and_logs()
    truncation_log()
    split_log()
    params_log()
    recent()
    messages()
    manifest()
    shutil.rmtree(WORK)
    print("golden files written to", HERE)


if __name__ == "__main__":
    main()
