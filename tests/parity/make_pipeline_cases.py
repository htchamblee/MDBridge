"""What MDATK's analyze pipeline (mdatk/pipeline.py Pipeline.required_analyses / analyze /
load_cached_dataset, mdatk/io/cache.py CacheStore) does, for src/analysis/pipeline.cpp and
cache.cpp (tests/analysis/test_pipeline.cpp).

    /opt/mdbridge-venv/bin/python -I -B make_pipeline_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies (numpy, scipy, MDAnalysis 2.11.0); nothing is written into the checkout.
MDATK runs in-process with MDBridge's independent tables patched into MDAnalysis
(independent_tables.apply(), and the residue-scoped ion rule on every Universe the pipeline loads;
tests/data/README.md). Reads tests/data/topology/demo_dry.parm7 and
tests/data/io/trajectories/demo_dry.nc. Writes OUT_DIR/pipeline_cases.json:

  required    Pipeline.required_analyses(video, charts) of projects with every panel type (structure
              with and without representations / proximity, sequence, secondary with each source,
              time series, dual axis, heatmap, contact map), chart series (instances, invalid ones,
              external data, keys without "/"), analyses switched on and off, an unknown panel
              type: the names or the error, and the log lines (skipped data sources;
              sorted, since MDATK logs them in the order of a Python set)
  sequences   runs of Pipeline(project).analyze(names, force) one after the other on one cache
              folder (a new Pipeline per run, as separate `mdatk analyze` calls make): the cache is
              reused, invalidated by parameter and frame changes, forced; frame selections,
              external data, explicit instance lists, no analyses at all, a reused cleaned
              trajectory, and errors (Pipeline.check(), unknown analyses and groups, empty frame
              selections, bad external data). Per run:
                project     Project.to_dict() of the run (tagged JSON, paths as placeholders)
                names, force, check
                log, progress   the Reporter's lines and [fraction (float.hex), message] calls
                error       {type, message} when the run raises (then nothing else)
                keys        [instance, cache key] in the order analysed
                manifest    MDATK_cache.json without project, mdatk_version and updated (its
                            all_entries sorted by the masked keys), and manifest_keys (the order
                            of its top-level keys)
                entries     per manifest entry, the result's .json file without project,
                            mdatk_version and created
                cached      the "<name>-<key>" entries in the cache folder after the run (sorted)
                dataset     frames, time_ps (float.hex), time unit, results (names in order),
                            describe(), and per data key the label / unit / target (series),
                            label / kind / unit / rows (matrices), label / shape (pair maps)
                reloaded    describe() and frames of Pipeline.load_cached_dataset() afterwards
              and `files`: small input files the case writes into its work folder first
              (external data; for the cleaning case, copies of the demo files as cleaned output
              plus a marker, written by the C++ test with its own cleaning key)

Placeholders: {data} for tests/data, {work} for the case's temporary folder. Cache keys hash the
input fingerprints (absolute paths, mtime_ns), so every key is replaced by the key MDATK computes for
the masked fingerprint (paths as placeholders, mtime_ns 0) - the same masking the C++ test applies -
and mtime_ns values are written as 0. The output is reproducible.
"""

import json
import logging
import os
import shutil
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

if len(sys.argv) != 3:
    sys.exit("usage: make_pipeline_cases.py MDATK_CHECKOUT OUT_DIR")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(sys.argv[1]))
sys.path.insert(0, HERE)

import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)
import improvements_analyses  # noqa: E402

improvements_analyses.apply_mdatk()  # MDATK with MDBridge's analysis improvements (docs/DIFFERENCES.md)

import improvements_cli  # noqa: E402

improvements_cli.apply()  # MDBridge's improvements (validation, external data before the frame loop, ...)

import numpy as np  # noqa: E402

import mdatk.cleaning.cpptraj as cpptraj  # noqa: E402
import mdatk.pipeline as pipeline_module  # noqa: E402
from mdatk import cli  # noqa: E402
from mdatk.analysis.base import analysis_class, instance_project  # noqa: E402
from mdatk.config import PanelSpec, Project  # noqa: E402
from mdatk.demo import demo_project  # noqa: E402
from mdatk.io import loader  # noqa: E402
from mdatk.progress import Reporter  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
DEMO_TOP = os.path.join(DATA, "topology", "demo_dry.parm7")
DEMO_NC = os.path.join(DATA, "io", "trajectories", "demo_dry.nc")

_load_universe = loader.load_universe


def load_universe(topology, trajectories=(), dt_ps=None):
    """load_universe with MDBridge's residue-scoped ion rule (as make_analysis_cases.py loads)."""
    u = _load_universe(topology, trajectories, dt_ps)
    fmt = loader._TOPOLOGY_FORMATS.get(os.path.splitext(topology)[1].lower())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        independent_tables.residue_rule_universe(u, topology, **({"topology_format": fmt} if fmt else {}))
    return u


pipeline_module.load_universe = load_universe


# --- masking ----------------------------------------------------------------------------------------

class Masker:
    def __init__(self, work):
        self.work = work
        self.keys = {}  # real cache key -> key of the masked fingerprint

    def text(self, s):
        s = s.replace(self.work, "{work}").replace(DATA, "{data}")
        for real, fixed in self.keys.items():
            s = s.replace(real, fixed)
        return s

    def __call__(self, o):
        if isinstance(o, dict):
            return {self.text(str(k)): (0 if k == "mtime_ns" else self(v)) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [self(v) for v in o]
        if isinstance(o, str):
            return self.text(o)
        if isinstance(o, np.ndarray):
            return self(o.tolist())
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, float) and not np.isfinite(o):
            return {"$float": "nan" if np.isnan(o) else ("inf" if o > 0 else "-inf")}
        return o

    def unmask(self, o):
        if isinstance(o, dict):
            return {k: self.unmask(v) for k, v in o.items()}
        if isinstance(o, list):
            return [self.unmask(v) for v in o]
        if isinstance(o, str):
            return o.replace("{work}", self.work).replace("{data}", DATA)
        return o


def error(exc):
    if isinstance(exc, KeyError) and exc.args:
        return {"type": "KeyError", "message": str(exc.args[0])}
    return {"type": type(exc).__name__, "message": str(exc)}


# --- projects ---------------------------------------------------------------------------------------

def base_project():
    """The demo project (as `mdatk demo` writes it for demo_dry) with the analyses MDBridge runs."""
    p = demo_project(DEMO_TOP, DEMO_NC, "{work}/out", 12, cleaning=False, mid_time_ps=120.0)
    p.analyses.mg.enabled = False
    p.analyses.stacking.enabled = False
    p.layout.preset = "large_plus_3"
    p.layout.panels = [PanelSpec("structure", [], "Structure"), PanelSpec("timeseries", ["rmsd/Complex"], "RMSD"),
                       PanelSpec("timeseries", ["energy/total"], "Energy"),
                       PanelSpec("contact_map", ["contacts/map"], "", {"mode": "instant"})]
    p.output.charts.series = ["hbonds/count"]
    return p


def with_overrides(make, overrides):
    p = make()
    for assignment in overrides:
        cli.apply_override(p, assignment)
    return p


# --- required_analyses ------------------------------------------------------------------------------

REQUIRED = [
    ("base", [], None, None),
    ("video_off", [], False, None),
    ("charts_off", [], None, False),
    ("nothing", ["output.video.enabled=false", "output.charts.enabled=false"], None, None),
    ("explicit_flags", ["output.video.enabled=false", "output.charts.enabled=false"], True, True),
    ("all_enabled", ["analyses.mg.enabled=true", "analyses.stacking.enabled=true"], None, None),
    ("none_enabled", ["analyses.contacts.enabled=false", "analyses.hbonds.enabled=false",
                      "analyses.energy.enabled=false", "analyses.rmsd.enabled=false"], None, False),
    ("rmsd_no_targets", ["analyses.rmsd.targets=[]", "output.video.enabled=false",
                         "output.charts.enabled=false"], None, None),
    ("structure_representations", ['visual.representations=[{"style": "cartoon"}]', "visual.glow_mode=\"none\""],
     None, False),
    ("structure_plain", ["visual.glow_mode=\"none\"", "visual.color_mode=\"groups\""], None, False),
    ("structure_proximity_colour", ["visual.glow_mode=\"none\"", "visual.color_mode=\"proximity\""], None, False),
    ("sequence", ['layout.panels=[{"type": "sequence"}]', "visual.glow_mode=\"none\""], None, False),
    ("sequence_glow", ['layout.panels=[{"type": "sequence"}]', "visual.glow_mode=\"ever\""], None, False),
    ("secondary_dynamic", ['layout.panels=[{"type": "secondary"}]', "visual.glow_mode=\"none\""], None, False),
    ("secondary_first", ['layout.panels=[{"type": "secondary", "options": {"source": "first"}}]'], None, False),
    ("secondary_dotbracket", ['layout.panels=[{"type": "secondary", "options": {"source": "dotbracket"}}]',
                              "visual.glow_mode=\"none\""], None, False),
    ("panels_series", ['layout.panels=[{"type": "heatmap", "series": ["mg/occupancy"]}, '
                       '{"type": "dual_axis", "series": ["stacking@Target/count", "external/x"]}, '
                       '{"type": "timeseries", "series": ["nokey", "trace/x", "contacts@Target+Aptamer/count"]}, '
                       '{"type": "contact_map", "series": ["hbonds@Aptamer+Target/map"]}]'], None, False),
    ("charts_series", ['output.charts.series=["external/col", "rmsd@Target/Complex", "contacts/count", '
                       '"contacts/count", "plain", "basepairs/x", "proximity/y"]', "output.video.enabled=false"],
     None, None),
    ("invalid_instances", ['output.charts.series=["bogus/x", "contacts@Nope/count", "rmsd@Aptamer+Target/a", '
                           '"trace@Aptamer/x", "contacts@Aptamer/count", "contacts@/count", "hbonds@Target+Aptamer/x"]',
                           "output.video.enabled=false"], None, None),
    ("unknown_panel", ['layout.panels=[{"type": "structure"}, {"type": "pie", "series": ["rmsd/x"]}]'], None, None),
]


def required_cases():
    work = tempfile.mkdtemp(prefix="mdbridge-pipeline-")
    mask = Masker(work)
    cases = []
    try:
        for name, overrides, video, charts in REQUIRED:
            p = with_overrides(base_project, overrides)
            p.output.out_dir = os.path.join(work, "out")
            logs = []
            case = {"name": name, "overrides": overrides, "video": video, "charts": charts,
                    "project": mask(p.to_dict())}
            pipe = pipeline_module.Pipeline(p, Reporter(log_cb=logs.append))
            try:
                case["names"] = pipe.required_analyses(video, charts)
            except Exception as exc:  # noqa: BLE001 - MDATK's error is the expected behaviour
                case["error"] = error(exc)
            case["log"] = sorted(logs)  # MDATK's order is that of a set of str (hash randomised)
            cases.append(case)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return cases


# --- analyze sequences ------------------------------------------------------------------------------

EXTERNAL = "time_ps,potential,temp K\n0,1.5,300\n100,2.5,301\n200,4.0,299.5\n300,3.0,300.25\n"


def mask_keys(mask, pipe, p):
    """The key of every instance analyze() computed keys for, for the masked fingerprint."""
    if pipe._inputs is None:
        return
    fp = mask(pipe.input_fingerprint(pipe._inputs))
    for inst, key in pipe.analysis_keys.items():
        mask.keys[key] = analysis_class(inst).cache_key(instance_project(p, inst), fp)


def run_record(mask, work, make, overrides, names=None, force=False, check=False):
    p = with_overrides(make, overrides)
    p = Project.from_dict(mask.unmask(p.to_dict()), None)
    logs, progress = [], []
    reporter = Reporter(log_cb=logs.append, progress_cb=lambda f, m: progress.append([float(f).hex(), m]))
    rec = {"overrides": overrides, "names": names, "force": force, "check": check}
    pipe = pipeline_module.Pipeline(p, reporter)
    try:
        if check:
            pipe.check()
        dataset = pipe.analyze(names, force=force)
    except Exception as exc:  # noqa: BLE001
        mask_keys(mask, pipe, p)
        rec["project"] = mask(p.to_dict())
        rec["error"] = mask(error(exc))
        rec["log"] = mask(logs)
        rec["progress"] = mask(progress)
        return rec
    mask_keys(mask, pipe, p)
    rec["project"] = mask(p.to_dict())
    rec["log"] = mask(logs)
    rec["progress"] = mask(progress)
    rec["keys"] = mask([[inst, key] for inst, key in pipe.analysis_keys.items()])
    with open(pipe.cache.manifest_path, encoding="utf-8") as fh:
        manifest = json.load(fh)
    rec["manifest_keys"] = list(manifest)
    entries = {}
    for entry in manifest["entries"]:
        with open(os.path.join(pipe.cache.cache_dir, f"{entry['name']}-{entry['key']}.json"), encoding="utf-8") as fh:
            info = json.load(fh)
        for k in ("project", "mdatk_version", "created"):
            info.pop(k)
        entries[entry["name"]] = mask(info)
    for k in ("project", "mdatk_version", "updated"):
        manifest.pop(k)
    rec["manifest"] = mask(manifest)
    # MDATK sorts all_entries by the real keys, which depend on paths and mtimes: by the masked ones
    rec["manifest"]["all_entries"].sort(key=lambda e: (e["name"], e["key"]))
    rec["entries"] = entries
    rec["cached"] = sorted(mask([f[:-4] for f in os.listdir(pipe.cache.cache_dir) if f.endswith(".npz")]))
    labels = {}
    for key in dataset.series_keys():
        s = dataset.series(key)
        labels[key] = {"label": s.label, "unit": s.unit, "target": s.target}
    for key in dataset.matrix_keys():
        m = dataset.matrix(key)
        labels[key] = {"label": m.label, "kind": m.kind, "unit": m.unit, "rows": len(m.row_labels)}
    for key in dataset.pairmap_keys():
        pm = dataset.pairmap(key)
        labels[key] = {"label": pm.label, "shape": list(pm.shape)}
    rec["dataset"] = mask({"frames": [int(f) for f in dataset.frames],
                           "time_ps": [float(t).hex() for t in dataset.time_ps],
                           "time_unit": dataset.time_unit, "results": list(dataset.results),
                           "describe": dataset.describe(), "labels": labels})
    try:
        reloaded = pipeline_module.Pipeline(p, Reporter()).load_cached_dataset()
        rec["reloaded"] = mask({"frames": [int(f) for f in reloaded.frames], "describe": reloaded.describe(),
                                "results": list(reloaded.results)})
    except Exception as exc:  # noqa: BLE001
        rec["reloaded"] = {"error": mask(error(exc))}
    return rec


def sequence(name, runs, files=None, prepare=None):
    work = tempfile.mkdtemp(prefix="mdbridge-pipeline-")
    mask = Masker(work)
    try:
        for rel, text in (files or {}).items():
            path = os.path.join(work, rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
        extra = prepare(work, mask) if prepare else {}
        records = [run_record(mask, work, base_project, *run) for run in runs]
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return dict({"name": name, "files": files or {}, "runs": records}, **extra)


EXT = ['inputs.external=[{"path": "{work}/ext.csv", "label": "MD"}]']


def cleaned_prepare(work, mask):
    """A cleaned trajectory as MDATK's cleaning leaves it (copies of the demo files + the marker)."""
    out = os.path.join(work, "cleaned")
    os.makedirs(out)
    shutil.copyfile(DEMO_TOP, os.path.join(out, cpptraj.CLEAN_PARM))
    shutil.copyfile(DEMO_NC, os.path.join(out, cpptraj.CLEAN_TRAJ))
    p = with_overrides(base_project, CLEANED)
    p = Project.from_dict(mask.unmask(p.to_dict()), None)
    key = cpptraj.cleaning_key(p.cleaning, p.inputs.topology, p.inputs.trajectories)
    with open(os.path.join(out, cpptraj.MARKER), "w", encoding="utf-8") as fh:
        json.dump({"key": key}, fh)
    # the key for masked fingerprints (what the C++ test computes from its own files)
    real_fp = cpptraj.file_fingerprint
    cpptraj.file_fingerprint = lambda path: mask(real_fp(path))
    try:
        fixed = cpptraj.cleaning_key(p.cleaning, p.inputs.topology, p.inputs.trajectories)
    finally:
        cpptraj.file_fingerprint = real_fp
    mask.keys[key] = fixed
    files = {"{work}/cleaned/" + cpptraj.CLEAN_PARM: "{data}/topology/demo_dry.parm7",
             "{work}/cleaned/" + cpptraj.CLEAN_TRAJ: "{data}/io/trajectories/demo_dry.nc"}
    return {"cleaning_key": fixed, "cleaned_files": files, "marker": "{work}/cleaned/" + cpptraj.MARKER}


CLEANED = ["cleaning.enabled=true", 'cleaning.output_dir="{work}/cleaned"']


def sequences():
    return [
        sequence("reuse", [
            ([],),                                                   # everything computed
            ([],),                                                   # everything from the cache
            (["analyses.contacts.cutoff=5.0"],),                     # contacts again, the rest cached
            (["analyses.contacts.cutoff=5.0", "analyses.energy.target_value=-75.5"],),  # display only
            (["analyses.contacts.cutoff=5.0"], None, True),          # forced
            (["analyses.contacts.cutoff=5.0", "inputs.step=2"],),    # other frames: all again
            (["inputs.step=2", "analyses.energy.target_value=null"], ["contacts@Target+Aptamer", "rmsd",
                                                                      "proximity"]),
            (["inputs.step=2"], ["rmsd", "rmsd", "trace"]),          # a repeated instance
            ([], []),                                                # no analyses: frames only
        ]),
        sequence("frames_external", [
            (["inputs.start=2", "inputs.stop=20", "inputs.step=3", 'inputs.time_unit="ps"', *EXT],
             ["hbonds", "rmsd@Target"]),
            (["inputs.start=-5", "inputs.dt_ps=2.5", *EXT, 'inputs.external.0.x_unit="frame1"'], ["contacts"]),
            ([*EXT], []),                                            # external data on the frame axis only
            (["inputs.start=2", "inputs.stop=20", "inputs.step=3", *EXT], ["hbonds", "rmsd@Target"]),
        ], files={"ext.csv": EXTERNAL}),
        sequence("errors", [
            ([], ["nope"]),
            ([], ["contacts@Nope"]),
            ([], ["contacts@Nope+Aptamer+Mg"]),
            ([], ["contacts@Nope+Aptamer"]),
            ([], ["trace@Aptamer"]),
            ([], ["rmsd@Aptamer+Target"]),
            (["inputs.step=0"], None, False, True),
            (["groups.0.residues=\"\"", "analyses.contacts.cutoff=-1"], None, False, True),
            (["inputs.start=100"], ["rmsd"]),
            ([*EXT, 'inputs.external.0.columns=["zz"]'], ["rmsd"]),
            (['inputs.external=[{"path": "{work}/missing.csv"}]'], ["rmsd"]),
            (['inputs.external=[{"path": "{work}/bad.csv"}]'], ["rmsd"]),
            ([], ["rmsd"], False, True),                             # valid in the end
        ], files={"ext.csv": EXTERNAL, "bad.csv": "a,b\n1,x\n"}),
        sequence("cleaned", [
            (CLEANED, ["rmsd", "contacts"]),
            (CLEANED, ["rmsd", "contacts"]),
        ], prepare=cleaned_prepare),
    ]


def main():
    # MDBridge's wording of KeyErrors and empty frame selections (improvements_cli.py, CLI-01 / CLI-03)
    data = improvements_cli.improve_strings({"required": required_cases(), "sequences": sequences()})
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "pipeline_cases.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False, allow_nan=False)
        fh.write("\n")


if __name__ == "__main__":
    main()
