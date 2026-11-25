"""What MDATK's residue shading, orientation and example dataset (mdatk/render/coloring.py,
residues.py, orientation.py, example.py) and numpy's default_rng give, for src/render/core
(shading.cpp, orientation.cpp, example.cpp, numpy_random.cpp; tests/render/test_shading.cpp;
milestone 4, chunk W0-H of docs/specs/m4-plan.md).

    /opt/mdbridge-venv/bin/python -I -B make_shading_cases.py MDATK_CHECKOUT OUT_DIR

MDATK runs in-process, patched by render_common.setup() (MDBridge's independent tables and
analysis improvements, improvements_render.py and improvements_shading.py: RN-21 ... RN-25). No
trajectory is read. Writes into OUT_DIR:

  shading_cases.json  crc32 of keys; numpy's PCG64 (raw outputs, SeedSequence states, scripted
                      uniform / integers / normal calls); glow_intensity, static_colors and
                      frame_colors cases; one_letter; orientation (euler_to_matrix, wrap_degrees,
                      matrix_to_euler, drag, project, format_angles); chain_names, resnames,
                      residue_rows, tag_rows and residue_shading of the example traces
  examples.json       per case: the project (Project.to_dict), n_frames, the preview passed (if
                      any) and the dataset example_dataset makes: frames, time, unit and per result
                      its name, key, params, meta and array names
  arrays.npz          the arrays: "<case>.<result>.<array>" of the examples, "rng.<seed>.<op>" of
                      the long random sequences, "shading.<case>.<what>" of residue_shading

Numbers are written as repr() (exact); the npz members are stored (render_common.write_npz). numpy's
kernels are pinned (render_common): sin / cos / exp are then the C library's, as in MDBridge.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import render_common as rc  # noqa: E402  (pins numpy's kernels before numpy is imported)

CHECKOUT, OUT = rc.args("make_shading_cases.py MDATK_CHECKOUT OUT_DIR")
rc.setup(CHECKOUT)

import copy  # noqa: E402
import types  # noqa: E402
import zlib  # noqa: E402

import numpy as np  # noqa: E402
from numpy.random.bit_generator import SeedSequence  # noqa: E402

from mdatk.config import PanelSpec, Project, RepSpec  # noqa: E402
from mdatk.render import coloring, example, orientation, residues  # noqa: E402
from mdatk.selections import SelectionError  # noqa: E402

ARRAYS = {}


def error(fn):
    try:
        fn()
    except Exception as exc:  # noqa: BLE001
        return {"type": type(exc).__name__, "message": str(exc)}
    return None


# --- crc32 and numpy's random generator -----------------------------------------------------------

KEYS = ["rmsd/Complex", "energy/total", "hbonds/count", "stacking@Aptamer+Target/count", "", "a",
        "external/cpp.rmsd", "q/Åµ", "x" * 300]

SEEDS = [0, 1, 42, 12345, 2 ** 32 - 1, 2 ** 32, 2 ** 40 + 7, 2 ** 64 - 1]
LONG_SEEDS = (0, 42, 2 ** 40 + 7)  # long sequences in arrays.npz

# scripted calls on one generator: (op, args); long ones go to arrays.npz
SCRIPT = [("uniform", (0.0, 2 * np.pi)), ("uniform", (0.8, 2.0)), ("integers", (8, 30)),
          ("normal", (0.0, 1.0, 120)), ("integers", (0, 1)), ("integers", (-5, 5, True)),
          ("integers", (0, 2 ** 32 - 1, True)), ("integers", (3, 2 ** 40)), ("integers", (-2 ** 62, 2 ** 62)),
          ("uniform", (-3.0, 7.5)), ("integers", (8, 30)), ("random", ()), ("normal", (2.5, 0.05, 7)),
          ("integers", (0, 7)), ("normal", (-120.0, 4.0, 3)), ("integers", (1, 2 ** 63 - 1))]


def run_op(g, op, a):
    if op == "uniform":
        return float(g.uniform(*a))
    if op == "random":
        return float(g.random())
    if op == "normal":
        return [float(v) for v in g.normal(a[0], a[1], a[2])]
    endpoint = len(a) > 2 and a[2]
    return int(g.integers(a[0], a[1], endpoint=endpoint))


def rng_cases():
    out = {"crc32": [[k, zlib.crc32(k.encode("utf-8"))] for k in KEYS], "seeds": []}
    for seed in SEEDS:
        state = SeedSequence(seed).generate_state(8, np.uint32)
        raw = np.random.PCG64(seed).random_raw(12)
        g = np.random.default_rng(seed)
        script = [[op, list(a), run_op(g, op, a)] for op, a in SCRIPT]
        if seed in LONG_SEEDS:
            g = np.random.default_rng(seed)
            ARRAYS[f"rng.{seed}.normal"] = g.normal(0.0, 1.0, 8000)  # every ziggurat branch
            ARRAYS[f"rng.{seed}.uniform"] = g.uniform(0.0, 1.0, 500)
            ARRAYS[f"rng.{seed}.integers"] = np.array([g.integers(0, 1000003) for _ in range(300)], dtype=np.int64)
        out["seeds"].append({"seed": str(seed), "state": [int(v) for v in state],
                             "raw": [str(int(v)) for v in raw], "script": script})
    return out


# --- coloring.py ----------------------------------------------------------------------------------

def glow_cases():
    r = np.random.default_rng(7)
    cases = []
    d0 = np.array([[9, 9, 2, 2, 9, 9, 9, 9, 9, 9], [0] * 10], dtype=float)  # test_render_cache
    m0 = np.array([False, True])
    d1 = np.round(r.uniform(0.0, 9.0, (6, 40)), 2)
    d1[2, 5] = np.nan
    d1[3, :] = np.inf
    m1 = np.array([False, False, False, False, True, False])
    for name, d, m in (("mdatk_test", d0, m0), ("random", d1, m1)):
        for mode in coloring.GLOW_MODES + ("disco",):
            for fade in (0, 3, 15, -2):
                for cutoff in (4.0, 2.5):
                    case = {"input": name, "mode": mode, "fade_frames": fade, "cutoff": cutoff}
                    try:
                        case["glow"] = coloring.glow_intensity(d, m, cutoff, mode, fade).tolist()
                    except ValueError as exc:
                        case["error"] = str(exc)
                    cases.append(case)
    return {"inputs": {"mdatk_test": {"d": d0.tolist(), "member": m0.tolist()},
                       "random": {"d": d1.tolist(), "member": m1.tolist()}}, "cases": cases}


def color_cases():
    out = []
    chains = [[0, 0, 1], [3, 1, 1, 0, 3, 7, 2], [0] * 4, []]
    gis = [[0, -1, 1], [0, 1, 2, -1, 5, 1, 0], [-1, -1, 0, 1], []]
    for chain, gi in zip(chains, gis):
        for mode in ("chain", "groups", "proximity", "other"):
            for pal in ("tab10", "viridis", "Set2_r"):
                for gcols in (["#ff0000", "#0000ff"], ["tab:green", "C3", "0.5"]):
                    st = coloring.static_colors(mode, np.array(chain, dtype=np.int64), np.array(gi, dtype=np.int64),
                                                gcols, "#808080", pal)
                    rec = {"mode": mode, "chain": chain, "group_index": gi, "group_colors": gcols,
                           "base": "#808080", "palette": pal, "static": st.tolist()}
                    if len(chain):
                        inten = np.linspace(0, 1, len(chain)).astype(np.float32)
                        member = np.array([k % 3 == 1 for k in range(len(chain))])
                        rec["intensity"] = [float(v) for v in inten]
                        rec["member"] = member.tolist()
                        rec["frame"] = {m: np.asarray(coloring.frame_colors(st, m, "#b5b3ad", "#e34948", inten,
                                                                            member)).tolist()
                                        for m in ("proximity", "groups")}
                    out.append(rec)
    return out


ONE_LETTER = ["G5", "RA", "LG3", "0C", "U", "DT", "HIS", "XYZ", "ala", " gly ", "HIE", "HIP", "CYX", "A", "rG",
              "ADE", "GUA", "CYT", "URA", "THY", "DA5", "C3", "MG", "WAT", "", "  ", "Na+", "hıs", "LYS",
              "DG3", "RU5", "A3", "N", "PSU", "ÅLA", "0DA", "0DT", "LDG", "LDC5", "LA", "LDU"]

# --- orientation.py -------------------------------------------------------------------------------

ANGLES = [[0, 0, 0], [10, 20, 30], [-170, 85, 179], [90, 90, 0], [30, -90, 45], [180, 0, 0],
          [-180, 0, -180], [359, -361, 720], [0.0001, 89.9999, -0.0001], [45.5, -30.25, 12.125],
          [-179.999, 0, 0], [179.9995, 10, -179.9999], [1e-12, -1e-12, 0], [400, 270, -90]]
WRAP = [0.0, -0.0, 180.0, -180.0, 540.0, -540.0, 179.99999, -179.99999, -179.999, -179.9985, -179.998,
        360.0, -360.0, 1e-300, -1e-12, 720.5, -0.5, 1e6 + 0.25, float("nan"), float("inf")]
DRAGS = [([10, 20, 30], 40, -15), ([0, 0, 0], 0, 0), ([0, 0, 0], 180, 0), ([0, 89, 0], 3, 2), ([-170, 5, 170], -7, 400),
         ([12.5, -45, 100], 1, 1)]


def orientation_cases():
    out = {"rx": [], "euler": [], "wrap": [], "to_euler": [], "drag": [], "format": [], "project": []}
    for a in (0, 30, -45.5, 90, 180, 1e-9, 1234.5):
        out["rx"].append([a, orientation.rx(a).tolist(), orientation.ry(a).tolist(), orientation.rz(a).tolist()])
    for ang in ANGLES:
        r = orientation.euler_to_matrix(ang)
        out["euler"].append([ang, r.tolist()])
        out["to_euler"].append([r.tolist(), list(orientation.matrix_to_euler(r))])
    gimbal = [[[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.5, 0.0], [0.3, 0.9, -0.2], [1.0, 0.0, 0.0]],
              [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0000001, 0.0, 0.0]]]
    for r in gimbal:
        out["to_euler"].append([r, list(orientation.matrix_to_euler(np.array(r)))])
    for a in WRAP:
        out["wrap"].append([a, orientation.wrap_degrees(a)])
    for ang, dx, dy in DRAGS:
        out["drag"].append([ang, dx, dy, list(orientation.drag(ang, dx, dy))])
    out["drag"].append([[5, 6, 7], 10, 10, list(orientation.drag([5, 6, 7], 10, 10, degrees_per_pixel=0.25)), 0.25])
    for ang in ANGLES + [[-0.04, 0.05, 0.25], [0.35, -0.35, 2.25], [float("nan"), float("inf"), -float("inf")]]:
        out["format"].append([ang, orientation.format_angles(ang)])
    pts = np.array([[1.0, 2.0, 3.0], [-4.5, 0.25, 9.0], [0.0, 0.0, 0.0]])
    for ang in ANGLES[:6]:
        out["project"].append([ang, pts.tolist(), [0.5, -1.0, 2.0],
                               orientation.project(pts, ang, np.array([0.5, -1.0, 2.0])).tolist()])
    return out


# --- example.py and residues.py -------------------------------------------------------------------

def panels(*specs):
    return [PanelSpec(t, list(s)) for t, s in specs]


def preview_of(p, atoms=True):
    """A preview dict as Pipeline.structure_preview() returns it, built from the schematic duplex
    (with chain names, labels, group membership, two ions)."""
    s = example.schematic_structure(p)
    tr = s["trace"]
    n = len(tr["chain"])
    out = {"coords": np.asarray(tr["coords"], dtype=np.float64), "ions": np.array([[1.0, 2.0, 3.0], [-4.0, 0.5, 6.25]]),
           "links": np.asarray(tr["links"], dtype=np.int64), "chain": np.asarray(tr["chain"], dtype=np.int64),
           "group_index": np.asarray(tr["group_index"], dtype=np.int64),
           "group_member": np.array([np.asarray(tr["group_index"]) == k for k in range(len(p.groups))]).reshape(len(p.groups), n),
           "resid": np.asarray(tr["resid"], dtype=np.int64), "resname": np.asarray(tr["resname"], dtype="U8"),
           "chain_name": np.array(["P" if c == 0 else "Q" for c in tr["chain"]]),
           "labels": [f"{'PQ'[int(c)]}:{b}{r}" for c, b, r in zip(tr["chain"], tr["resname"], tr["resid"])]}
    if atoms:
        a = s["atoms"]
        out["atoms"] = {"coords": np.asarray(a["coords"], dtype=np.float64), "names": list(a["names"]),
                        "elements": list(a["elements"]), "resnames": list(a["resnames"]),
                        "row": np.asarray(a["row"], dtype=np.int64), "bonds": np.asarray(a["bonds"], dtype=np.int64)}
    return out


def project(**kw):
    p = Project()
    for k, v in kw.items():
        obj = p
        parts = k.split("__")
        for part in parts[:-1]:
            obj = getattr(obj, part)
        setattr(obj, parts[-1], v)
    return p


def example_projects():
    """(case, project, preview, n_frames)"""
    cases = [("default", Project(), None, 120)]
    p = project(layout__preset="large_plus_4")
    p.layout.panels = panels(("structure", []), ("timeseries", ["stacking@Aptamer+Target/count"]),
                             ("heatmap", ["mg/occupancy"]), ("contact_map", ["hbonds/map"]),
                             ("dual_axis", ["rmsd/Complex", "energy/total"]))
    cases.append(("v002", p, None, 120))
    p = project(inputs__dt_ps=2.5, inputs__time_unit="ps")
    p.layout.panels = panels(("timeseries", ["rmsd/Complex", "energy/net_force", "energy@Target+Aptamer/total",
                                             "contacts@Target+Aptamer/q_native", "external/cpp.rmsd", "nosuch/x",
                                             "rmsd", "hbonds/count", "rmsd/Complex"]),
                             ("heatmap", ["contacts/per_residue_b", "contacts/per_residue", "mg/occupancy",
                                          "hbonds@Target+Nope/per_residue"]),
                             ("contact_map", ["hbonds/map", "contacts@Target+Aptamer/map"]),
                             ("nosuchpanel", ["rmsd/Other"]), ("sequence", []), ("secondary", []))
    cases.append(("all_kinds", p, None, 30))
    p = Project()
    p.visual.representations = [RepSpec(style="cartoon")]
    p.layout.panels = panels(("structure", []))
    cases.append(("representations", p, None, 12))
    p = Project()
    p.layout.panels = panels(("structure", []), ("timeseries", ["rmsd/Complex"]))
    cases.append(("preview", p, preview_of(p), 10))
    cases.append(("preview_trace_only", copy.deepcopy(p), preview_of(p, atoms=False), 10))
    p = project(visual__structure_groups=["Aptamer"])
    p.layout.panels = panels(("structure", []), ("sequence", []))
    cases.append(("single_strand", p, None, 20))  # RN-21
    p = project(visual__structure_groups=["Aptamer"])
    p.visual.representations = [RepSpec(style="tube")]
    p.layout.panels = panels(("secondary", []))
    cases.append(("single_strand_atoms", p, None, 14))  # RN-21
    p = Project()
    p.groups[0].residues = "1-1000000000"
    p.groups[1].residues = "500-650, 700"
    p.layout.panels = panels(("structure", []), ("heatmap", ["mg/occupancy"]), ("contact_map", ["hbonds/map"]))
    cases.append(("big_groups", p, None, 6))  # RN-25
    p = Project()
    p.groups[0].residues = "B:5-8, 3, 20-18"
    p.groups[1].residues = "1-4, x"
    p.groups.append(copy.deepcopy(p.groups[0]))
    p.groups[2].name = "Target"
    p.groups[2].residues = "40-45"
    p.groups.append(copy.deepcopy(p.groups[0]))
    p.groups[3].name = "Empty"
    p.groups[3].residues = "  "
    p.visual.structure_groups = ["Target", "Empty", "Aptamer"]
    p.layout.panels = panels(("structure", []), ("heatmap", ["contacts/per_residue"]),
                             ("contact_map", ["contacts@Target+Empty/map"]))
    cases.append(("odd_groups", p, None, 13))
    p = Project()
    p.groups = []
    p.layout.panels = panels(("structure", []), ("timeseries", ["hbonds/count"]))
    cases.append(("no_groups", p, None, 8))
    p = Project()
    p.layout.panels = panels(("timeseries", ["hbonds/count", "rmsd/Complex"]), ("heatmap", ["mg/occupancy"]))
    cases.append(("one_frame", p, None, 1))  # MDATK's NaN (t / t[-1]) kept
    return cases


def dump_dataset(case, ds):
    results = []
    for name, res in ds.results.items():
        results.append({"name": res.name, "key": res.key, "params": res.params, "meta": res.meta,
                        "arrays": list(res.arrays)})
        for k, a in res.arrays.items():
            ARRAYS[f"{case}.{name}.{k}"] = np.ascontiguousarray(a)
    return {"frames": ds.frames, "time_ps": ds.time_ps, "time_unit": ds.time_unit, "results": results}


def preview_json(pv):
    if pv is None:
        return None
    out = {k: v for k, v in pv.items() if k != "atoms"}
    if "atoms" in pv:
        out["atoms"] = pv["atoms"]
    return out


SHADING_CASES = ("default", "all_kinds", "representations", "preview", "single_strand", "odd_groups", "no_groups",
                 "one_frame")
SPECS = ["", "  ", "12", "B:13-15", "A:1-3, B:20", "C:1", "1-4;7", "x", "13-12", "Q:13", "P:2-3", "1000"]
TAGS = [(0, ""), (5, "1, 12-13"), (4, ""), (0, "B:13-15"), (3, "x, 2"), (0, "  "), (-2, "5"), (7, "A:2"),
        (2, "Q:14, 1"), (0, "1-")]
SHADING_VISUALS = [("groups_fade", {}), ("chain_current", {"color_mode": "chain", "glow_mode": "current"}),
                   ("proximity_none", {"color_mode": "proximity", "glow_mode": "none"}),
                   ("proximity_ever", {"color_mode": "proximity", "glow_mode": "ever", "proximity_cutoff": 2.0,
                                       "proximity_color": "C2", "glow_color": "tab:purple"}),
                   ("groups_none", {"glow_mode": "none"}),
                   ("fade3", {"glow_mode": "fade", "fade_frames": 3, "chain_palette": "Paired"})]


def residue_cases(datasets):
    out = []
    for case, (p, ds) in datasets.items():
        if "trace" not in ds.results:
            continue
        tr = ds.results["trace"]
        rec = {"case": case, "chain_names": residues.chain_names(tr).tolist(),
               "resnames": residues.resnames(tr).tolist(), "rows": [], "tags": [], "shading": []}
        for spec in SPECS:
            try:
                rec["rows"].append([spec, residues.residue_rows(tr, spec).tolist()])
            except SelectionError as exc:
                rec["rows"].append([spec, {"error": str(exc)}])
        for every, spec in TAGS:
            v = copy.deepcopy(p.visual)
            v.label_every, v.label_residues = every, spec
            rc.reset_log()
            residues._WARNED.clear()
            rows = residues.tag_rows(tr, v).tolist()
            rec["tags"].append([every, spec, rows, list(rc.LOG_LINES)])
        for name, changes in SHADING_VISUALS if case in SHADING_CASES else ():
            q = copy.deepcopy(p)
            for k, v in changes.items():
                setattr(q.visual, k, v)
            ctx = types.SimpleNamespace(dataset=ds, project=q)
            for need in (None, True, False) if name == SHADING_VISUALS[0][0] else (None,):
                key = f"shading.{case}.{name}.{need}"
                sh = residues.residue_shading(ctx, tr, need)
                ARRAYS[key + ".static"] = np.asarray(sh.static, dtype=np.float64)
                ARRAYS[key + ".member"] = np.asarray(sh.member, dtype=bool)
                ARRAYS[key + ".glow"] = np.asarray(sh.glow, dtype=np.float32)
                ARRAYS[key + ".color_i"] = np.asarray(sh.color_i, dtype=np.float32)
                rec["shading"].append({"name": name, "visual": changes, "need": need, "key": key,
                                       "glow_rgb": sh.glow_rgb.tolist()})
        out.append(rec)
    return out


def remapped_case():
    """residue_shading with a proximity result whose rows differ from the trace's."""
    p = Project()
    p.layout.panels = panels(("structure", []))
    ds = example.example_dataset(p, n_frames=9)
    prox = ds.results["proximity"]
    order = np.array([5, 3, 99, 0, 1, 2, 4, 7, 8, 6, 3] + list(range(9, 24)), dtype=np.int64)
    d = np.asarray(prox.arrays["min_distance"])
    prox.arrays["resindex"] = order
    prox.arrays["min_distance"] = (d[np.clip(order, 0, 23)] - 2.5).astype(np.float32)
    prox.arrays["member"] = np.array([k % 5 == 0 for k in range(len(order))])
    ctx = types.SimpleNamespace(dataset=ds, project=p)
    sh = residues.residue_shading(ctx, ds.results["trace"])
    ARRAYS["remap.member"] = np.asarray(sh.member, bool)
    ARRAYS["remap.glow"] = np.asarray(sh.glow, np.float32)
    ARRAYS["remap.color_i"] = np.asarray(sh.color_i, np.float32)
    ARRAYS["remap.resindex"] = order
    ARRAYS["remap.min_distance"] = prox.arrays["min_distance"]
    ARRAYS["remap.prox_member"] = prox.arrays["member"]
    missing = error(lambda: residues.residue_shading(
        types.SimpleNamespace(dataset=types.SimpleNamespace(n_frames=3, results={}), project=p),
        ds.results["trace"], True))
    return {"project": rc.project_dict(p), "n_frames": 9, "missing": missing}


def main():
    os.makedirs(OUT, exist_ok=True)
    examples, datasets = [], {}
    for case, p, pv, n in example_projects():
        ds = example.example_dataset(p, pv, n_frames=n)
        datasets[case] = (p, ds)
        examples.append({"case": case, "project": rc.project_dict(p), "n_frames": n, "preview": preview_json(pv),
                         "dataset": dump_dataset(case, ds)})
    group_resids = []
    for case, (p, _ds) in datasets.items():
        for g in [g.name for g in p.groups] + ["Nope"]:
            for fb in (12, 3):
                group_resids.append([case, g, fb, example.group_resids(p, g, fb)])
    cases = {"rng": rng_cases(), "glow": glow_cases(), "colors": color_cases(),
             "one_letter": [[n, residues.one_letter(n)] for n in ONE_LETTER], "orientation": orientation_cases(),
             "group_resids": group_resids, "residues": residue_cases(datasets), "remap": remapped_case(),
             "schematic": {"trace": rc.plain(example.schematic_trace(datasets["default"][0], 1)),
                           "structure": rc.plain(example.schematic_structure(datasets["default"][0]))}}
    rc.write_json(os.path.join(OUT, "shading_cases.json"), cases)
    rc.write_json(os.path.join(OUT, "examples.json"), {"cases": examples})
    rc.write_npz(os.path.join(OUT, "arrays.npz"), ARRAYS)


main()
