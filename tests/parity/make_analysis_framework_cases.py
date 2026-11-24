"""What MDATK's analysis framework computes (mdatk/analysis/base.py, __init__.py, common.py), for
src/analysis (tests/analysis/test_framework.cpp).

    /opt/mdbridge-venv/bin/python -I -B make_analysis_framework_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies (numpy, MDAnalysis 2.11.0); nothing is written into the checkout. MDATK
runs in-process with MDBridge's independent tables patched into MDAnalysis
(independent_tables.apply(), tests/data/README.md). Writes OUT_DIR/framework_cases.json:

  registry       ANALYSES in order: name, label, version, internal, group_slots, OUTPUT_LABELS
  inputs         the fixed input fingerprints the cache keys are computed for (no file is read)
  projects       per project (given as the plain data Project.from_dict reads): for every analysis
                 enabled, repr(params_for), selection_texts, expected_outputs, groups,
                 output_label of its outputs and of unknown ones, cache_key per input fingerprint,
                 referenced_groups(selection_texts) and the repr of new_result(key, **meta);
                 expected_keys with and without disabled analyses; referenced_groups of extra texts;
                 per instance id: parse_instance, instance_id, and instance_project's outcome
                 (the instance's groups, params, expected outputs, labels and cache keys, or the
                 error); key_label of data keys
  split_key      split_key of keys
  residue_rows, without, per_row_counts
                 common.py on index lists of a small fake topology and on random CSR triplets
                 (float32 values as float.hex() of the float32 value)
  run            run_analyses of recording analyses over a MemoryReader trajectory: frames,
                 log lines, progress calls (fraction as float.hex()), the setup / process /
                 finalize calls (their number only for long runs), times (float.hex())

Error entries are {"type": Python class, "message": str(exc), or the KeyError's argument}.
The output is reproducible.
"""

import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_analysis_framework_cases.py MDATK_CHECKOUT OUT_DIR")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(sys.argv[1]))
sys.path.insert(0, HERE)

import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)

import numpy as np  # noqa: E402

from mdatk.analysis import ANALYSES, AnalysisContext, expected_keys, run_analyses  # noqa: E402
from mdatk.analysis import common  # noqa: E402
from mdatk.analysis.base import (Analysis, instance_id, instance_project, key_label,  # noqa: E402
                                 parse_instance, referenced_groups, split_key)
from mdatk.config import Project  # noqa: E402
from mdatk.progress import Reporter  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
rng = np.random.default_rng(20261124)

INPUTS = [
    {"topology": {"path": "/data/demo dry/sys\u00fc.parm7", "size": 123456, "mtime_ns": 1763900000123456789},
     "trajectories": [{"path": "/data/demo dry/run1.nc", "size": 98765, "mtime_ns": 1763900001000000000}],
     "frames": {"start": 0, "stop": None, "step": 1}, "dt_ps": None},
    {"topology": {"path": "C:/x/top.pdb", "size": 1, "mtime_ns": 0},
     "trajectories": [{"path": "a.xtc", "size": 2, "mtime_ns": 3}, {"path": "b.xtc", "size": 4, "mtime_ns": 5}],
     "frames": {"start": 5, "stop": 100, "step": 3}, "dt_ps": 2.5},
]

GROUPS_RICH = [
    {"name": "Aptamer", "residues": "1-12", "selection": "", "color": "#2a78d6"},
    {"name": "Target", "residues": "13-24", "selection": "", "color": "#eb6834"},
    {"name": "Mg", "residues": "", "selection": "resname MG", "color": "#00aa00"},
    {"name": "Pocket", "residues": "", "selection": "group:Aptamer and around 5 group:Mg", "color": "#111111"},
    {"name": "Weird.name-1'", "residues": "3", "selection": "", "color": "#222222"},
    {"name": "Both", "residues": "", "selection": " group:Pocket or group:Target or group:Nope ", "color": "#333333"},
    {"name": "My group", "residues": "4-5", "selection": "", "color": "#444444"},
    {"name": "Gruppe-\u00e4", "residues": "6", "selection": "group:Self", "color": "#555555"},
    {"name": "Self", "residues": "", "selection": "group:Self or group:Gruppe-\u00e4", "color": "#666666"},
]

PROJECTS = {
    "default": {},
    "rich": {
        "groups": GROUPS_RICH,
        "analyses": {
            "use_pbc": True,
            "proximity_search_radius": 6,
            "contacts": {"group_a": "Pocket", "group_b": "Both", "cutoff": 5},
            "hbonds": {"group_a": "Weird.name-1'", "group_b": "My group", "use_angle": False},
            "mg": {"ligand_selection": "group:Weird.name-1' or group:Both", "occupancy_selection": "Gruppe-\u00e4",
                   "mg_selection": "group:Mg", "geometry": True},
            "stacking": {"scope": "inter", "selection": "My group", "group_a": "Self", "group_b": "Nope"},
            "energy": {"target_value": -60.0, "cutoff": 10, "distance_dependent": True, "exclude_bonded": 0},
            "rmsd": {"targets": [{"label": "H\u00e9avy atoms!", "selection": "group:Pocket and not name H*"},
                                 {"label": "___", "selection": "group:Target", "fit_selection": "group:Aptamer"},
                                 {"label": "", "selection": "name P"},
                                 {"label": "a b/c..d", "selection": "group:My group"},
                                 {"label": "Backbone", "selection": "group:Self"}],
                     "reference_frame": 5, "mass_weighted": True, "superposition": False},
            "basepairs": {"engine": "barnaba", "stride": 1, "cutoff": 3, "wobble": False, "noncanonical": True},
        },
        "visual": {"structure_groups": ["Target", "Nope", "Mg", "Target"], "show_ions": False,
                   "proximity_set": "sel:resid 3", "trace_atoms": "P, C4'", "align_frames": False,
                   "rep_hydrogens": True},
    },
    "disabled": {
        "analyses": {"contacts": {"enabled": False}, "hbonds": {"enabled": False},
                     "rmsd": {"enabled": True, "targets": []}, "mg": {"enabled": False},
                     "stacking": {"enabled": True, "scope": "all"}},
        "visual": {"show_ions": True},
    },
    "duplicates": {
        "groups": [{"name": "Aptamer", "residues": "1-6"}, {"name": "Target", "residues": "13-24"},
                   {"name": "Aptamer", "residues": "7-12", "selection": "group:Target"}],
        "analyses": {"rmsd": {"enabled": False}, "energy": {"group_a": "Target", "group_b": "Target"}},
    },
    "no_groups": {"groups": [], "analyses": {"contacts": {"group_a": "", "group_b": ""}}},
}

INSTANCES = [
    "contacts", "contacts@Target+Aptamer", "contacts@Aptamer", "contacts@Nope+Target", "contacts@X+X",
    "contacts@", "contacts@Target++Aptamer", "contacts@+", "contacts@Target+Aptamer+Mg", "rmsd@Target",
    "rmsd@Target+Aptamer", "rmsd@Pocket", "rmsd@Weird.name-1'", "trace@Target", "trace", "trace@",
    "basepairs@Aptamer", "proximity@Nope", "bogus", "bogus@A+B", "", "@Target", "energy@Target+Target",
    "stacking@Target+Aptamer", "stacking@Both+Self", "mg@Aptamer+Target", "hbonds@Weird.name-1'+Both",
    "hbonds@My group+Gruppe-\u00e4", "atoms@Aptamer+Target", "Contacts@Target+Aptamer", "contacts@Target@Aptamer",
    "contacts@Aptamer+Target/x",
]

EXTRA_TEXTS = [
    ["  Aptamer  "], ["group:Pocket and not group:"], ["notgroup:Target"], ["group:Nope"],
    ["group:Target)or(group:Mg"], [""], ["   "], ["group:group:Aptamer"], ["Both", "Self"],
    ["group:Weird.name-1'"], ["group:Weird.name-1'x"], ["Gruppe-\u00e4"], ["group:Gruppe-\u00e4"],
    ["My group"], ["group:My group"], ["(group:Target)", "resname MG"], [], ["\u00a0Target\u3000"],
]


def error(exc):
    if isinstance(exc, KeyError):
        return {"type": "KeyError", "message": exc.args[0]}
    return {"type": type(exc).__name__, "message": str(exc)}


def output_names(cls, project):
    names = []
    for kind, outs in cls.expected_outputs(project).items():
        names += [o for o in outs if o not in names]
    names += [k for k in cls.OUTPUT_LABELS if k not in names]
    return names + ["bogus", ""]


def describe(cls, project):
    out = {
        "enabled": bool(cls.enabled(project)),
        "params": repr(cls.params_for(project)),
        "selection_texts": list(cls.selection_texts(project)),
        "expected_outputs": cls.expected_outputs(project),
        "groups": list(cls.groups(project)),
        "labels": [[o, cls.output_label(project, o)] for o in output_names(cls, project)],
        "cache_keys": [cls.cache_key(project, fp) for fp in INPUTS],
        "referenced_groups": referenced_groups(cls.selection_texts(project), project),
    }
    return out


def new_result_repr(cls, project, instance, key, meta):
    obj = cls.__new__(cls)
    obj.project = project
    obj.instance = instance or cls.name
    res = obj.new_result(key, **meta)
    return {"name": res.name, "key": res.key, "params": repr(res.params), "meta": repr(res.meta)}


def project_case(name, plain):
    p = Project.from_dict(plain)
    case = {"name": name, "plain": plain, "analyses": {}}
    for n, cls in ANALYSES.items():
        d = describe(cls, p)
        d["new_result"] = [new_result_repr(cls, p, "", "", {}),
                           new_result_repr(cls, p, f"{n}@X", "0123abcd", {"x": 1, "label": "old", "z": [1.5]})]
        case["analyses"][n] = d
    case["expected_keys"] = expected_keys(p)
    case["expected_keys_enabled"] = expected_keys(p, include_disabled=False)
    case["referenced_groups"] = [[texts, referenced_groups(texts, p)] for texts in EXTRA_TEXTS]
    instances = []
    for inst in INSTANCES:
        entry = {"instance": inst}
        base, groups = parse_instance(inst)
        entry["parse"] = [base, list(groups)]
        entry["id"] = instance_id(base, groups)
        try:
            q = instance_project(p, inst)
            cls = ANALYSES[base]
            entry["project"] = describe(cls, q)
            entry["unchanged"] = q is p
        except (KeyError, ValueError) as exc:
            entry["error"] = error(exc)
        instances.append(entry)
    case["instances"] = instances
    keys = []
    for inst in INSTANCES:
        base, _ = parse_instance(inst)
        outs = output_names(ANALYSES[base], p) if base in ANALYSES else ["x"]
        keys += [f"{inst}/{o}" for o in outs] + [inst]
    keys += ["x/y/z", "rmsd/Complex", "rmsd@Target/Heavy", "rmsd@Target/zzz", "contacts@Target+Aptamer/map/extra",
             "/", "//"]
    case["key_labels"] = [[k, key_label(k, p)] for k in keys]
    return case


def instance_id_cases():
    cases = []
    for base, groups in [("contacts", ["A", "B"]), ("contacts", ["", "B", ""]), ("rmsd", []), ("rmsd", [""]),
                         ("x", ["A+B"]), ("", ["A"]), ("", [])]:
        cases.append([base, groups, instance_id(base, groups)])
    return cases


class FakeResidues:
    def __getitem__(self, uniq):
        return uniq


class FakeUniverse:
    residues = FakeResidues()


class FakeAtoms:
    """The part of an AtomGroup that common.py uses."""

    universe = FakeUniverse()

    def __init__(self, ix, atom_residue):
        self.ix = np.asarray(ix, dtype=np.int64)
        self._atom_residue = np.asarray(atom_residue, dtype=np.int64)
        self.resindices = self._atom_residue[self.ix]
        self.n_atoms = len(self.ix)

    def __getitem__(self, mask):
        return FakeAtoms(self.ix[mask], self._atom_residue)


ATOM_RESIDUE = [0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 5, 6, 7, 7, 7, 7]  # 20 atoms, 8 residues


def atom_lists():
    lists = [[], [0], list(range(20)), [19, 3, 3, 0, 12], [5, 6, 7, 8], [15], [16, 9, 2, 17, 10]]
    for _ in range(8):
        k = int(rng.integers(1, 15))
        lists.append([int(x) for x in rng.integers(0, 20, size=k)])
    return lists


def common_cases():
    rows = []
    lists = atom_lists()
    for ix in lists:
        residues, r = common.residue_rows(FakeAtoms(ix, ATOM_RESIDUE))
        rows.append({"atoms": ix, "residues": [int(x) for x in residues], "rows": [int(x) for x in r]})
    without = []
    for a in lists[:8]:
        for b in lists[:8]:
            got = common.without(FakeAtoms(a, ATOM_RESIDUE), FakeAtoms(b, ATOM_RESIDUE))
            without.append({"atoms": a, "other": b, "result": [int(x) for x in got.ix]})
    counts = []
    specs = [(0, 0, 3), (1, 0, 2), (3, 5, 4), (5, 40, 6), (12, 300, 9), (2, 0, 0)]
    for n_frames, n_entries, n_rows in specs:
        cuts = np.sort(rng.integers(0, n_entries + 1, size=max(n_frames - 1, 0)))
        offsets = np.concatenate([[0], cuts, [n_entries]]).astype(np.int64) if n_frames else np.zeros(1, np.int64)
        if n_rows:
            r = rng.integers(0, n_rows, size=n_entries).astype(np.int32)
        else:
            r = np.zeros(0, np.int32)
        v = (rng.random(n_entries) * rng.choice([1.0, 3.7, 1e-3, 1e4], size=n_entries)).astype(np.float32)
        v[::7] = 1.0  # whole counts too
        out = common.per_row_counts(offsets, r, v, n_rows)
        assert out.dtype == np.float32
        counts.append({"offsets": [int(x) for x in offsets], "rows": [int(x) for x in r],
                       "values": [float(x).hex() for x in v], "n_rows": n_rows,
                       "shape": list(out.shape), "result": [float(x).hex() for x in out.ravel()]})
    return rows, without, counts


class Recorder(Analysis):
    """A streaming analysis that records what the runner calls."""

    name = "recorder"

    def __init__(self, ctx, label, calls):
        super().__init__(ctx)
        self.label = label
        self.calls = calls

    def setup(self, frame_indices):
        super().setup(frame_indices)
        self.calls.append(["setup", self.label, [int(f) for f in frame_indices]])

    def process(self, k, ts):
        self.calls.append(["process", self.label, int(k), int(ts.frame)])

    def finalize(self):
        self.calls.append(["finalize", self.label, int(self.n_frames)])
        return self.label


def run_cases():
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    cases = []
    n_total = 450
    for frames, labels, dt in [([0], ["A"], 1.0), ([2, 0, 1, 1], ["First", "Second"], 0.5),
                               (list(range(199)), ["Only"], 2.0), (list(range(200)), ["A", "", "C"], 1.0),
                               (list(range(401)), ["X"], 10.0), (list(range(449, -1, -1)), ["Rev"], 0.1),
                               ([], ["None"], 1.0), ([-1, -450], ["Neg"], 1.0)]:
        u = mda.Universe.empty(1, trajectory=False)
        coords = np.zeros((n_total, 1, 3), dtype=np.float32)
        u.load_new(coords, format=MemoryReader, dt=dt)
        logs, progress, calls = [], [], []
        rep = Reporter(log_cb=logs.append, progress_cb=lambda f, m: progress.append([float(f).hex(), m]))
        ctx = AnalysisContext(u, Project(), "", rep)
        objs = [Recorder(ctx, label, calls) for label in labels]
        results, times = run_analyses(ctx, objs, np.asarray(frames, dtype=np.int64))
        case = {"n_total": n_total, "dt": float(dt).hex(), "frames": frames, "labels": labels,
                "logs": logs, "progress": progress, "results": results, "times": [float(t).hex() for t in times]}
        if len(frames) <= 10:
            case["calls"] = calls
        else:  # long runs: the call sequence is checked by the C++ test from the frames
            case["n_calls"] = len(calls)
        cases.append(case)
    return cases


def main():
    registry = []
    for n, cls in ANALYSES.items():
        registry.append({"name": n, "label": cls.label, "version": cls.version, "internal": bool(cls.internal),
                         "group_slots": cls.group_slots, "output_labels": [[k, v] for k, v in cls.OUTPUT_LABELS.items()]})
    rows, without, counts = common_cases()
    data = {
        "registry": registry,
        "inputs": INPUTS,
        "projects": [project_case(n, plain) for n, plain in PROJECTS.items()],
        "instance_id": instance_id_cases(),
        "split_key": [[k, list(split_key(k))] for k in ["contacts@A+B/atom_pairs", "rmsd", "a/b/c", "/x", "x/", ""]],
        "residue_rows": rows,
        "without": without,
        "per_row_counts": counts,
        "run": run_cases(),
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "framework_cases.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
        fh.write("\n")


if __name__ == "__main__":
    main()
