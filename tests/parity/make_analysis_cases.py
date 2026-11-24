"""What MDATK's analyses (mdatk/analysis/*.py) compute (tests/data/analysis), for the C++ analysis
module of milestone 3 (docs/specs/m3-plan.md, sections 4.1 / 4.2, chunk W0.6).

    /opt/mdbridge-venv/bin/python -I -B make_analysis_cases.py MDATK_CHECKOUT OUT_DIR [--native]

Needs MDAnalysis 2.11.0, numpy and scipy (MDATK's dependencies); MDATK's tests/conftest.py,
test_analyses.py, test_basepairs.py and test_pdbtools.py are imported for their geometry helpers (a
stand-in for pytest is installed when pytest is missing). Nothing is written into the checkout.
MDAnalysis runs with MDBridge's independent tables (independent_tables.apply(), and the
residue-scoped ion rule on every topology file; tests/data/README.md). Reads
tests/data/topology/demo_dry.parm7, tests/data/io/trajectories/demo_dry.nc and its DCD / XTC / TRR
copies (written by earlier generators). Writes into OUT_DIR:

  micro_cases.json      the micro systems of MDATK's tests (atoms, bonds and float32 coordinates of
                        conftest.make_universe) and the analyses run on them
  demo_dry_cases.json   the demo project (mdatk.demo.demo_project, as `mdatk demo` writes it) on
                        demo_dry.parm7 + demo_dry.nc: every analysis, all 24 frames

Every case is one run of the analysis runner, built like Pipeline.analyze builds it: a project
(the base project with `--set`-style overrides, mdatk.cli.apply_override), the analysis instances
(each with instance_project's project and its own AnalysisContext, all sharing one Reporter),
frame_indices(project.inputs) and run_analyses(...) in one pass over a fresh Universe. Recorded:

  case        name, system, overrides, instances, the project sections the analyses read
              (groups, analyses, visual, inputs) as tagged JSON
  frames      analysed frame indices; times (ts.time, float.hex)
  log         the Reporter's log lines in order ("Preparing ...", "Analysing ...", "WARNING: ...")
  progress    [fraction (hex), message] of every progress call
  python_warnings
              "Category: message" of every Python warning raised during the run (in order)
  setup_frame per instance, the trajectory frame current when its setup() started (rmsd seeks
              in setup; atoms guesses bonds from the current frame in setup)
  error       {type, message} when instance_project, a setup, a frame or a finalize raises
              (then no results); paths in messages as {data} / {out}
  results     per instance: params (tagged JSON), meta (tagged JSON, with outputs and labels),
              output_order (the add_series / add_matrix / add_pairmap order, which sorted keys
              lose), key (cache_key for the case's fixed fingerprint, below) and every array:
              {dtype (numpy dtype.str), shape, sha256 (of the C-order little-endian bytes),
              values (flattened: int / bool as JSON, floats as float.hex(), str as str)}
  margins     per instance, the decision margins (below)
  fingerprint the input fingerprint the cache keys use: Pipeline.input_fingerprint's shape with
              fixed paths ({data}/..., {micro}/name), the real file sizes and mtime_ns 0

Decision margins. Wherever a result depends on a comparison with a threshold, the generator
evaluates the same quantity again, with MDATK's own functions on the same inputs (pairs_within with
cutoff + 1 Å, MDATK's formulas for the angles), and records, per analysed frame, the smallest
|quantity - threshold| over everything compared (float.hex; "inf" when nothing was within 1 Å of
the threshold; null for frames where the decision is not made). A C++ result that differs in a
count may be a last-bit flip only where that margin is tiny. Recorded quantities:

  contacts   cutoff (atom-pair distance)
  hbonds     distance_cutoff (donor-acceptor distance); angle_cos (cos D-H...A vs
             cos(radians(angle_cutoff)), for donor-acceptor pairs within the cutoff + 1 Å, only
             when the angle is used); h_guess_distance (H-donor distance vs 1.25, first frame,
             only when hydrogens are assigned from distances) and h_guess_tie (difference of the
             two nearest donor distances of one hydrogen, the "strict <" tie rule)
  mg         inner_cutoff / outer_cutoff (Mg-ligand distances), geometry_cos (min(||c| - cos90|,
             |c - cos180|) over the ligand pairs of one ion, geometry on), occupancy_inner /
             occupancy_outer (Mg-occupancy atom distances)
  stacking   planarity (ring RMS vs planarity_cutoff), distance_cutoff (centroid distances),
             angle_cos (|n_i.n_j| vs cos(radians(angle_cutoff))), offset (smaller lateral offset
             vs offset_cutoff), for centroid pairs within the cutoff + 1 Å
  energy     cutoff (pair distances); abs_elec / abs_vdw (sum of |term| of the pairs used) and
             abs_force (sum of |force vector| of the pairs) for relative tolerances
  proximity  radius
  trace      link_max (|trace-point distance - 9.0| of the link candidates, first frame)
  atoms      setup: bond_max / bond_h (|d - 1.95| / |d - 1.25| of the distance-guessed bonds, only
             when the topology has none, from the frame current at setup)
  basepairs  cutoff (every finite rule H-bond distance), c1_window (C1'-C1' distance of rule pairs
             with enough bonds within the cutoff + 1 Å vs 8.5 / 12.5), score_grid (distance, in Å, of
             each candidate's score from the rounding boundary of round(score / 0.05)), score_dev_tie
             (smallest |dev_a - dev_b| of candidates whose rounded scores are equal); detection
             frames only
  thresholds the cosines MDATK compares with, as numpy gives them and as libm's cos(x * pi / 180)
             (math.cos(math.radians(x))) gives them (numpy may use SIMD kernels)

Floating point and the CPU. numpy's SIMD kernels and OpenBLAS / LAPACK's kernels are chosen per
CPU, and their results differ in the last bits (SVD of kabsch and of the ring planes, BLAS
products, cos). So that the committed values do not depend on the machine that regenerates them,
the generator pins them before numpy is imported: OPENBLAS_CORETYPE=Nehalem (the generic SSE
kernels every x86-64 CPU runs), one BLAS / OpenMP thread, and on x86-64 numpy's dispatch limited to
its baseline (NPY_ENABLE_CPU_FEATURES=X86_V2). Values that depend on these kernels (rmsd with a
fit selection, aligned trace / atoms coordinates, ring normals and planarity) are therefore those
of the generic kernels; the C++ tests compare them within the tolerances of the plan (4.3) and use
the margins for counts. --native skips the pinning (to see what the machine's own kernels change;
its output is not the committed data).
"""

import os
import platform
import sys

NATIVE = "--native" in sys.argv
if not NATIVE:  # before numpy is imported (read when OpenBLAS / numpy load)
    os.environ["OPENBLAS_CORETYPE"] = "Nehalem"
    for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_var] = "1"
    if platform.machine().lower() in ("x86_64", "amd64"):
        os.environ["NPY_ENABLE_CPU_FEATURES"] = "X86_V2"

import hashlib  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import math  # noqa: E402
import re  # noqa: E402
import types  # noqa: E402
import warnings  # noqa: E402

ARGS = [a for a in sys.argv[1:] if a != "--native"]
if len(ARGS) != 2:
    sys.exit("usage: make_analysis_cases.py MDATK_CHECKOUT OUT_DIR [--native]")
CHECKOUT = os.path.abspath(ARGS[0])
OUT = os.path.abspath(ARGS[1])
sys.path.insert(0, CHECKOUT)

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)  # ctx.warn() also logs; the Reporter's lines are recorded

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
sys.path.insert(0, HERE)
import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)

from mdatk import cli  # noqa: E402
from mdatk import spatial as S  # noqa: E402
from mdatk.analysis import ANALYSES, AnalysisContext, run_analyses  # noqa: E402
from mdatk.analysis import stacking as stacking_mod  # noqa: E402
from mdatk.analysis import structure as structure_mod  # noqa: E402
from mdatk.analysis.base import analysis_class, instance_project, parse_instance  # noqa: E402
from mdatk.analysis.energy import COULOMB  # noqa: E402
from mdatk.config import Project  # noqa: E402
from mdatk.demo import demo_project  # noqa: E402
from mdatk.io import loader  # noqa: E402
from mdatk.progress import Reporter  # noqa: E402
from mdatk.secstruct import C1_WINDOW, SCORE_TIE  # noqa: E402

DELTA = 1.0  # margins: pairs within threshold + DELTA are evaluated


# --- MDATK's test helpers ---------------------------------------------------------------------------

def _pytest_stand_in():
    """MDATK's test modules import pytest for decorators only; without pytest a stand-in will do."""
    try:
        import pytest  # noqa: F401
        return
    except ImportError:
        pass

    def decorator(*args, **kwargs):
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return lambda f: f

    class Mark:
        def __getattr__(self, name):
            return decorator

    module = types.ModuleType("pytest")
    module.fixture = decorator
    module.mark = Mark()
    module.approx = lambda *a, **k: None
    module.raises = lambda *a, **k: None
    module.importorskip = lambda *a, **k: None
    module.skip = lambda *a, **k: None
    sys.modules["pytest"] = module


_pytest_stand_in()
sys.path.insert(0, os.path.join(CHECKOUT, "tests"))
import conftest  # noqa: E402
import test_analyses as ta  # noqa: E402
import test_basepairs as tb  # noqa: E402

from mdatk.render.orientation import euler_to_matrix, rx  # noqa: E402


# --- JSON helpers -----------------------------------------------------------------------------------

def fhex(x):
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x.hex()


def plain(o):
    """numpy scalars / arrays -> Python values (for meta and params)."""
    if isinstance(o, np.ndarray):
        return [plain(v) for v in o.tolist()]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, dict):
        return {str(k): plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [plain(v) for v in o]
    return o


def tagged(o):
    """Strict JSON: non-finite floats as {"$float": "nan" | "inf" | "-inf"} (python_cases.py)."""
    o = plain(o)
    if isinstance(o, float) and not math.isfinite(o):
        return {"$float": "nan" if math.isnan(o) else ("inf" if o > 0 else "-inf")}
    if isinstance(o, dict):
        return {str(k): tagged(v) for k, v in o.items()}
    if isinstance(o, list):
        return [tagged(v) for v in o]
    return o


def array_json(a):
    a = np.ascontiguousarray(np.asarray(a))
    kind = a.dtype.kind
    flat = a.reshape(-1)
    if kind == "b":
        values = [bool(v) for v in flat]
    elif kind in "iu":
        values = [int(v) for v in flat]
    elif kind == "f":
        values = [fhex(v) for v in flat]
    elif kind == "U":
        values = [str(v) for v in flat]
    else:
        raise TypeError(f"array dtype {a.dtype}")
    data = a.astype(a.dtype.newbyteorder("<"), copy=False).tobytes()
    return {"dtype": a.dtype.str, "shape": list(a.shape), "sha256": hashlib.sha256(data).hexdigest(),
            "values": values}


def placeholders(text):
    return str(text).replace(OUT, "{out}").replace(DATA, "{data}")


def error(exc):
    return {"type": type(exc).__name__, "message": placeholders(str(exc))}


def write_json(name, data):
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
        json.dump(data, fh, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
        fh.write("\n")


def remove_offsets():
    """MDAnalysis caches XTC / TRR frame offsets next to the files it reads: remove them."""
    for folder in (DATA, OUT):
        for root, _dirs, names in os.walk(folder):
            for name in names:
                if name.startswith(".") and re.search(r"_offsets\.(lock|npz)$", name):
                    os.remove(os.path.join(root, name))


# --- decision margins -------------------------------------------------------------------------------

def min_abs(values, threshold):
    values = np.asarray(values, dtype=np.float64)
    return float(np.min(np.abs(values - threshold))) if values.size else math.inf


def distance_margin(a, b, cutoff, box):
    _, _, d = S.pairs_within(a, b, cutoff + DELTA, box)
    return min_abs(d, cutoff)


def cos_thresholds(degrees):
    return {"angle": fhex(degrees), "numpy": fhex(np.cos(np.radians(degrees))),
            "libm": fhex(math.cos(math.radians(degrees)))}


class Margins:
    """Per-frame margins of one analysis instance."""

    def __init__(self, n_frames):
        self.n = n_frames
        self.frames = {}
        self.setup = {}
        self.thresholds = {}

    def put(self, k, name, value):
        self.frames.setdefault(name, [None] * self.n)[k] = fhex(value)

    def json(self):
        out = {"frames": self.frames}
        if self.setup:
            out["setup"] = {k: fhex(v) for k, v in self.setup.items()}
        if self.thresholds:
            out["thresholds"] = self.thresholds
        return out


def contacts_margins(obj, m, k, ts, state):
    m.put(k, "cutoff", distance_margin(obj.a.positions, obj.b.positions, obj.cutoff, obj.ctx.box(ts)))


def hbonds_margins(obj, m, k, ts, state):
    p = obj.project.analyses.hbonds
    pos = obj.ab.positions.astype(np.float64)
    box = obj.ctx.box(ts)
    if state.get("guessing"):
        heavy_idx = np.flatnonzero(obj.cand_donor)
        h_idx = np.flatnonzero(obj.is_h)
        i, j, dw = S.pairs_within(pos[h_idx], pos[heavy_idx], 1.25 + DELTA, box)
        # MDATK compares the plain (not minimum-image) norm for the nearest donor
        d = np.linalg.norm(pos[h_idx[i]] - pos[heavy_idx[j]], axis=1) if len(i) else np.zeros(0)
        m.put(k, "h_guess_distance", min_abs(dw, 1.25))
        tie = math.inf
        inside = dw <= 1.25
        for h in np.unique(i[inside]):
            ds = np.sort(d[inside][i[inside] == h])
            if len(ds) > 1:
                tie = min(tie, float(ds[1] - ds[0]))
        m.put(k, "h_guess_tie", tie)
    dist, cos = math.inf, math.inf
    for s in (0, 1):
        don, acc = obj.donor_idx[s], obj.acc_idx[1 - s]
        i, j, d = S.pairs_within(pos[don], pos[acc], obj.dcut + DELTA, box)
        dist = min(dist, min_abs(d, obj.dcut))
        if obj.use_angle and len(i):
            d_atoms, a_atoms = don[i], acc[j]
            n_h = obj.h_start[d_atoms + 1] - obj.h_start[d_atoms]
            rep_d = np.repeat(d_atoms, n_h)
            rep_a = np.repeat(a_atoms, n_h)
            offs = np.repeat(obj.h_start[d_atoms], n_h)
            within = np.arange(n_h.sum()) - np.repeat(np.cumsum(n_h) - n_h, n_h)
            hyd = obj.h_flat[offs + within]
            v1 = S.minimum_image_vectors(pos[rep_d] - pos[hyd], box)
            v2 = S.minimum_image_vectors(pos[rep_a] - pos[hyd], box)
            c = (v1 * v2).sum(1) / np.maximum(np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1), 1e-12)
            cos = min(cos, min_abs(c, obj.cos_cut))
    m.put(k, "distance_cutoff", dist)
    if obj.use_angle:
        m.put(k, "angle_cos", cos)
    m.thresholds["angle_cutoff"] = cos_thresholds(p.angle_cutoff)


def mg_margins(obj, m, k, ts, state):
    p = obj.p
    m.thresholds["cos90"] = cos_thresholds(90.0 - p.angle_tolerance)
    m.thresholds["cos180"] = cos_thresholds(180.0 - p.angle_tolerance)
    if obj.mg.n_atoms == 0:
        return
    box = obj.ctx.box(ts)
    mpos = obj.mg.positions.astype(np.float64)
    if obj.lig.n_atoms:
        lpos = obj.lig.positions.astype(np.float64)
        i, j, d = S.pairs_within(mpos, lpos, p.outer_cutoff + DELTA, box)
        m.put(k, "inner_cutoff", min_abs(d, p.inner_cutoff))
        m.put(k, "outer_cutoff", min_abs(d, p.outer_cutoff))
        if p.geometry:
            shell = (p.inner_cutoff if p.bridge_shell == "inner" else p.outer_cutoff) + DELTA
            keep = d <= shell
            bi, bj = i[keep], j[keep]
            g = math.inf
            for ion in np.unique(bi):
                ligs = bj[bi == ion]
                if len(ligs) < 2:
                    continue
                p1, p2 = np.triu_indices(len(ligs), 1)
                l1, l2 = ligs[p1], ligs[p2]
                v1 = S.minimum_image_vectors(lpos[l1] - mpos[ion], box)
                v2 = S.minimum_image_vectors(lpos[l2] - mpos[ion], box)
                c = (v1 * v2).sum(1) / np.maximum(np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1), 1e-12)
                g = min(g, float(np.min(np.minimum(np.abs(np.abs(c) - obj.cos90), np.abs(c - obj.cos180)))))
            m.put(k, "geometry_cos", g)
    if obj.occ_atoms.n_atoms:
        opos = obj.occ_atoms.positions.astype(np.float64)
        _, _, d = S.pairs_within(mpos, opos, p.outer_cutoff + DELTA, box)
        m.put(k, "occupancy_inner", min_abs(d, p.inner_cutoff))
        m.put(k, "occupancy_outer", min_abs(d, p.outer_cutoff))


def stacking_margins(obj, m, k, ts, state):
    p = obj.p
    m.thresholds["angle_cutoff"] = cos_thresholds(p.angle_cutoff)
    if obj.n_rings < 2:
        return
    box = obj.ctx.box(ts)
    pos = obj.atoms.positions.astype(np.float64)
    cent = np.zeros((obj.n_rings, 3))
    norm = np.zeros((obj.n_rings, 3))
    planar = np.zeros(obj.n_rings)
    for ring_ids, idx in obj.by_size.values():
        c, nvec, pl = stacking_mod.ring_geometry(pos[idx])
        cent[ring_ids], norm[ring_ids], planar[ring_ids] = c, nvec, pl
    m.put(k, "planarity", min_abs(planar, p.planarity_cutoff))
    cut = p.distance_cutoff + DELTA
    if p.scope == "inter":
        i, j, d = S.pairs_within(cent[obj.rows], cent[obj.cols], cut, box)
        i, j = obj.rows[i], obj.cols[j]
    else:
        i, j, d = S.pairs_within(cent, cent, cut, box)
        keep = (i < j) & (obj.res_ix[i] != obj.res_ix[j])
        i, j, d = i[keep], j[keep], d[keep]
    m.put(k, "distance_cutoff", min_abs(d, p.distance_cutoff))
    if len(i):
        v = S.minimum_image_vectors(cent[j] - cent[i], box)
        ni, nj = norm[i], norm[j]
        cosang = np.abs((ni * nj).sum(1))
        proj_i = (v * ni).sum(1)
        proj_j = (v * nj).sum(1)
        off_i = np.sqrt(np.maximum((v * v).sum(1) - proj_i ** 2, 0))
        off_j = np.sqrt(np.maximum((v * v).sum(1) - proj_j ** 2, 0))
        m.put(k, "angle_cos", min_abs(cosang, obj.cos_cut))
        m.put(k, "offset", min_abs(np.minimum(off_i, off_j), p.offset_cutoff))
    else:
        m.put(k, "angle_cos", math.inf)
        m.put(k, "offset", math.inf)


def energy_margins(obj, m, k, ts, state):
    p = obj.p
    box = obj.ctx.box(ts)
    pa = obj.a.positions.astype(np.float64)
    pb = obj.b.positions.astype(np.float64)
    m.put(k, "cutoff", distance_margin(pa, pb, p.cutoff, box))
    i, j, d = S.pairs_within(pa, pb, p.cutoff, box)
    if len(obj.excl) and len(i):
        keep = ~np.isin(i * obj.b.n_atoms + j, obj.excl)
        i, j, d = i[keep], j[keep], d[keep]
    if len(i) == 0:
        for name in ("abs_elec", "abs_vdw", "abs_force"):
            m.put(k, name, 0.0)
        return
    d = np.maximum(d, 0.5)
    qq = COULOMB * obj.qa[i] * obj.qb[j] / p.dielectric
    if p.distance_dependent:
        e_el, de_el = qq / (d * d), -2.0 * qq / d ** 3
    else:
        e_el, de_el = qq / d, -qq / (d * d)
    if obj.lj is not None:
        acoef, bcoef = obj.lj.pair_coefficients(obj.a.ix[i], obj.b.ix[j])
    else:
        rmin = obj.ra[i] + obj.rb[j]
        eps = np.sqrt(obj.ea[i] * obj.eb[j])
        r6 = rmin ** 6
        acoef, bcoef = eps * r6 * r6, 2.0 * eps * r6
    inv6 = 1.0 / d ** 6
    e_vdw = acoef * inv6 * inv6 - bcoef * inv6
    de_vdw = (-12.0 * acoef * inv6 * inv6 + 6.0 * bcoef * inv6) / d
    rvec = S.minimum_image_vectors(pa[i] - pb[j], box)
    fmag = -(de_el + de_vdw) / d
    m.put(k, "abs_elec", float(np.abs(e_el).sum()))
    m.put(k, "abs_vdw", float(np.abs(e_vdw).sum()))
    m.put(k, "abs_force", float(np.linalg.norm(rvec * fmag[:, None], axis=1).sum()))


def proximity_margins(obj, m, k, ts, state):
    if obj.struct.n_atoms:
        m.put(k, "radius", distance_margin(obj.struct.positions, obj.set_atoms.positions, obj.radius,
                                           obj.ctx.box(ts)))


def trace_margins(obj, m, k, ts, state):
    if not state.get("first"):
        return
    xyz = obj.tracer.positions()
    n = len(xyz)
    if n < 2:
        m.put(k, "link_max", math.inf)
        return
    i = np.arange(n - 1)
    j = i + 1
    cand = ((obj.chain[i] == obj.chain[j]) & ((obj.residues.resindices[j] - obj.residues.resindices[i]) == 1)
            & obj.is_point[i] & obj.is_point[j])
    m.put(k, "link_max", min_abs(np.linalg.norm(xyz[j] - xyz[i], axis=1)[cand], structure_mod.LINK_MAX))


def atoms_setup_margins(obj, m):
    atoms = obj.atoms
    try:
        has_bonds = any(True for b in atoms.bonds
                        if all(int(x) in set(atoms.ix.tolist()) for x in b.indices))
    except Exception:  # noqa: BLE001 - the topology has no bonds (as atom_bonds sees it)
        has_bonds = False
    if has_bonds or atoms.n_atoms < 2:
        return
    xyz = atoms.positions.astype(np.float64)
    i, j, d = S.pairs_within(xyz, xyz, 1.95 + DELTA, None)
    keep = i < j
    m.setup["bond_max"] = min_abs(d[keep], 1.95)
    m.setup["bond_h"] = min_abs(d[keep], 1.25)


def basepairs_margins(obj, m, k, ts, state):
    if obj.engine != "builtin" or k % obj.stride:
        return
    det = obj.detector
    xyz = obj.atoms.positions.astype(np.float64)
    cut, window, grid, tie = math.inf, math.inf, math.inf, math.inf
    lo, hi = C1_WINDOW
    for _kind, need, _rows_p, ip, c1_p, _rows_q, iq, c1_q in det.rules:
        p = xyz[np.maximum(ip, 0)]
        q = xyz[np.maximum(iq, 0)]
        d = np.linalg.norm(p[:, None, :, :] - q[None, :, :, :], axis=-1)
        d[(ip < 0)[:, None, :] | (iq < 0)[None, :, :]] = np.inf
        finite = d[np.isfinite(d)]
        cut = min(cut, min_abs(finite, det.cutoff))
        a, b = np.nonzero((d <= det.cutoff + DELTA).sum(axis=-1) >= need)
        ca, cb = c1_p[a], c1_q[b]
        both = (ca >= 0) & (cb >= 0)
        dc = np.linalg.norm(xyz[np.maximum(ca, 0)] - xyz[np.maximum(cb, 0)], axis=-1)[both]
        window = min(window, min_abs(dc, lo), min_abs(dc, hi))
    cand = det.candidates(xyz)
    if len(cand):
        x = cand[:, 3] / SCORE_TIE
        grid = float(np.min(np.abs(x - np.floor(x) - 0.5))) * SCORE_TIE
        rounded = np.round(x) * SCORE_TIE
        for r in np.unique(rounded):
            dev = np.sort(cand[rounded == r, 4])
            if len(dev) > 1:
                tie = min(tie, float(np.min(np.diff(dev))))
    m.put(k, "cutoff", cut)
    m.put(k, "c1_window", window)
    m.put(k, "score_grid", grid)
    m.put(k, "score_dev_tie", tie)


MARGINS = {"contacts": contacts_margins, "hbonds": hbonds_margins, "mg": mg_margins,
           "stacking": stacking_margins, "energy": energy_margins, "proximity": proximity_margins,
           "trace": trace_margins, "basepairs": basepairs_margins}


def instrument(obj, margins, setup_frames):
    """Wrap setup / process of one analysis object to record margins and the setup frame."""
    base = obj.name
    setup, process = obj.setup, obj.process

    def wrapped_setup(frame_indices):
        setup_frames[obj.instance] = int(obj.ctx.u.trajectory.ts.frame)
        setup(frame_indices)
        m = Margins(len(frame_indices))
        margins[obj.instance] = m
        if base == "atoms":
            atoms_setup_margins(obj, m)

    def wrapped_process(k, ts):
        state = {"guessing": bool(getattr(obj, "_finalize_donors_pending", False)),
                 "first": getattr(obj, "ref", 0) is None}
        process(k, ts)
        fn = MARGINS.get(base)
        if fn is not None:
            fn(obj, margins[obj.instance], k, ts, state)

    obj.setup = wrapped_setup
    obj.process = wrapped_process


# --- running a case ---------------------------------------------------------------------------------

ALL = list(ANALYSES)  # registry order: contacts, energy, hbonds, mg, rmsd, stacking, trace, proximity, atoms, basepairs


class System:
    """How to build the Universe of a case (fresh for every case) and its fixed fingerprint."""

    def __init__(self, name, make, topology_path="", fp_topology=None, fp_trajectories=()):
        self.name = name
        self.make = make
        self.topology_path = topology_path
        self.fp_topology = fp_topology
        self.fp_trajectories = list(fp_trajectories)


def file_system(name, topology, trajectories=(), dt=None):
    """A system loaded like MDATK loads files (load_universe), with the residue-scoped ion rule."""
    def make():
        u = loader.load_universe(topology, list(trajectories), dt)
        fmt = loader._TOPOLOGY_FORMATS.get(os.path.splitext(topology)[1].lower())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            independent_tables.residue_rule_universe(u, topology, **({"topology_format": fmt} if fmt else {}))
        return u

    def fp(path):
        return {"path": placeholders(path), "size": os.path.getsize(path), "mtime_ns": 0}

    return System(name, make, topology, fp(topology), [fp(t) for t in trajectories])


def project_sections(p):
    d = p.to_dict()
    return tagged({"groups": d["groups"], "analyses": d["analyses"], "visual": d["visual"],
                   "inputs": {k: d["inputs"].get(k) for k in ("start", "stop", "step", "dt_ps")}})


def run_case(name, system, base, overrides=(), instances=None, note=""):
    """One run of the analysis runner (as Pipeline.analyze builds it); returns the case record."""
    instances = list(instances or ALL)
    project = base()
    for assignment in overrides:
        cli.apply_override(project, assignment)
    case = {"name": name, "system": system.name, "overrides": list(overrides), "instances": instances,
            "project": project_sections(project)}
    if note:
        case["note"] = note
    fingerprint = {"topology": system.fp_topology, "trajectories": system.fp_trajectories,
                   "frames": {"start": project.inputs.start, "stop": project.inputs.stop,
                              "step": project.inputs.step}, "dt_ps": project.inputs.dt_ps}
    case["fingerprint"] = tagged(fingerprint)
    logs, progress = [], []
    reporter = Reporter(log_cb=logs.append, progress_cb=lambda f, msg: progress.append([fhex(f), msg]))
    margins, setup_frames = {}, {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            u = system.make()
            idx = loader.frame_indices(u.trajectory.n_frames, project.inputs.start, project.inputs.stop,
                                       project.inputs.step)
            projects = {n: instance_project(project, n) for n in instances}
            keys = {n: analysis_class(n).cache_key(projects[n], fingerprint) for n in instances}
            objs = [analysis_class(n)(AnalysisContext(u, projects[n], system.topology_path, reporter), n)
                    for n in instances]
            for obj in objs:
                instrument(obj, margins, setup_frames)
            results, times = run_analyses(objs[0].ctx, objs, idx, reporter)
        except Exception as exc:  # noqa: BLE001 - MDATK's error is the expected behaviour
            case["error"] = error(exc)
            results = None
    case["log"] = [placeholders(line) for line in logs]
    case["progress"] = progress
    case["python_warnings"] = [f"{w.category.__name__}: {placeholders(w.message)}" for w in caught
                               if "MDAnalysis" not in str(w.filename)]
    case["setup_frame"] = setup_frames
    if results is None:
        return case
    case["frames"] = [int(f) for f in idx]
    case["times"] = [fhex(t) for t in times]
    case["results"] = []
    for res in results:
        outputs = res.meta.get("outputs", {})
        case["results"].append({
            "instance": res.name,
            "analysis": parse_instance(res.name)[0],
            "key": keys[res.name],
            "params": tagged(res.params),
            "meta": tagged(res.meta),
            "output_order": {kind: list(outputs.get(kind, {})) for kind in ("series", "matrices", "pairmaps")},
            "arrays": {k: array_json(v) for k, v in res.arrays.items()},
            "margins": margins[res.name].json(),
        })
    return case


# --- micro systems (MDATK's tests/test_analyses.py and test_basepairs.py) ---------------------------

MICRO_SYSTEMS = {}


def micro(name, atoms, coords, bonds=(), dt=1.0):
    """A conftest.make_universe system: recorded as data for the C++ side to rebuild."""
    coords32 = np.asarray(coords, dtype=np.float32)
    MICRO_SYSTEMS[name] = {
        "atoms": [{"name": a["name"], "resid": int(a["resid"]), "resname": a["resname"],
                   "element": a["element"], "mass": fhex(a.get("mass", 12.0)),
                   "charge": fhex(a.get("charge", 0.0))} for a in atoms],
        "segid": "SYS", "dt": fhex(dt), "bonds": [[int(x), int(y)] for x, y in bonds],
        "n_frames": int(coords32.shape[0]),
        "coords": [[fhex(v) for v in frame.reshape(-1)] for frame in coords32],
    }

    def make():
        return conftest.make_universe(atoms, coords32, bonds=bonds, dt=dt)

    return System(name, make, "", {"path": "{micro}/" + name, "size": 0, "mtime_ns": 0}, [])


def micro_base(groups, **extra):
    def base():
        p = conftest.micro_project(groups)
        for key, value in extra.items():
            cli.apply_override(p, f"{key}={json.dumps(value)}")
        return p
    return base


def micro_from_frames(name, frames):
    """test_basepairs.universe's atoms from [(resname, {atom: xyz})] frames."""
    atoms = [dict(name=n, resid=r + 1, resname=resname, element=n.lstrip("0123456789")[0])
             for r, (resname, res) in enumerate(frames[0]) for n in res]
    coords = [[xyz for _, res in f for xyz in res.values()] for f in frames]
    return micro(name, atoms, coords)


def micro_cases():
    atom = ta.atom
    cases = []
    two = {"Aptamer": "1", "Target": "2"}

    # contacts: counts and native fraction
    s = micro("contacts", [atom("C1", 1, "A", "C"), atom("O1", 1, "A", "O"),
                           atom("N1", 2, "U", "N"), atom("H1", 2, "U", "H")],
              [[[0, 0, 0], [1, 0, 0], [4, 0, 0], [4.5, 0, 0]], [[0, 0, 0], [1, 0, 0], [10, 0, 0], [10.5, 0, 0]]])
    base = micro_base(two)
    cases.append(run_case("contacts_native", s, base, ["analyses.contacts.cutoff=3.5"], ["contacts"]))
    cases.append(run_case("contacts_all_atoms", s, base, ["analyses.contacts.cutoff=3.5",
                                                          "analyses.contacts.heavy_only=false"], ["contacts"]))
    cases.append(run_case("contacts_same_group", s, base, ["analyses.contacts.cutoff=5.0"],
                          ["contacts@Aptamer+Aptamer", "contacts@Target+Aptamer"]))
    cases.append(run_case("contacts_micro_all", s, base, [], ALL))

    # H-bonds: linear / bent / far, with and without bonds
    hb_atoms = [atom("N3", 1, "U", "N"), atom("H3", 1, "U", "H"), atom("O6", 2, "G", "O")]
    hb_frames = [[[0, 0, 0], [1.0, 0, 0], [2.9, 0, 0]], [[0, 0, 0], [1.0, 0, 0], [0, 2.9, 0]],
                 [[0, 0, 0], [1.0, 0, 0], [4.0, 0, 0]]]
    for bonded in (True, False):
        s = micro("hbonds_bonds" if bonded else "hbonds_nobonds", hb_atoms, hb_frames,
                  bonds=[(0, 1)] if bonded else ())
        cases.append(run_case(f"{s.name}_angle", s, base, [], ["hbonds", "hbonds@Target+Aptamer"]))
        cases.append(run_case(f"{s.name}_distance_only", s, base, ["analyses.hbonds.use_angle=false"],
                              ["hbonds"]))

    # Mg: bridges, octahedral geometry, shells, occupancy
    mg = np.array([5.0, 5.0, 5.0])
    ang60 = 2.1 * np.array([np.cos(np.radians(60)), np.sin(np.radians(60)), 0])
    mg_frames = [[mg + [2.1, 0, 0], mg + o2, mg] for o2 in ([0, 2.1, 0], ang60, [0, 3.5, 0], [0, 9.0, 0])]
    s = micro("mg", [atom("OP1", 1, "A", "O"), atom("OP2", 2, "U", "O"), atom("MG", 3, "MG", "MG")], mg_frames)
    mg_sel = 'analyses.mg.mg_selection=resname MG'
    cases.append(run_case("mg_default", s, base, [mg_sel], ["mg"]))
    cases.append(run_case("mg_geometry", s, base, [mg_sel, "analyses.mg.geometry=true"], ["mg"]))
    cases.append(run_case("mg_outer_shell", s, base, [mg_sel, "analyses.mg.bridge_shell=outer"], ["mg"]))
    cases.append(run_case("mg_default_selection", s, base, [], ["mg"]))
    cases.append(run_case("mg_no_ions", s, base, ["analyses.mg.mg_selection=resname XX"], ["mg"]))

    # stacking: stacked, offset, tilted rings
    st_atoms = [atom(n, resid, "U", n[0]) for resid in (1, 2) for n in ["N1", "C2", "N3", "C4", "C5", "C6"]]
    st_frames = []
    for center, rot in (([0, 0, 3.5], None), ([2.8, 0, 3.5], None), ([0, 0, 3.5], rx(60))):
        _, r1 = ta._ring([0, 0, 0])
        _, r2 = ta._ring(center, rot)
        st_frames.append(np.vstack([r1, r2]))
    s = micro("stacking", st_atoms, st_frames)
    sel = "analyses.stacking.selection=all"
    cases.append(run_case("stacking_default", s, base, [sel], ["stacking"]))
    cases.append(run_case("stacking_offset_3", s, base, [sel, "analyses.stacking.offset_cutoff=3.0"], ["stacking"]))
    cases.append(run_case("stacking_inter", s, base, [sel, "analyses.stacking.scope=inter"], ["stacking"]))
    cases.append(run_case("stacking_none", s, base, ["analyses.stacking.selection=resid 1"], ["stacking"]))

    # energy: analytic Coulomb + LJ (generic parameters)
    s = micro("energy_ions", [atom("NA", 1, "NA", "NA", charge=1.0), atom("CL", 2, "CL", "CL", charge=-1.0)],
              [[[0, 0, 0], [5.0, 0, 0]]])
    cases.append(run_case("energy_analytic", s, base, [], ["energy"]))
    cases.append(run_case("energy_analytic_eps4r", s, base, ["analyses.energy.dielectric=4.0",
                                                             "analyses.energy.distance_dependent=true"], ["energy"]))
    # energy: bonded exclusions across the groups
    chain_atoms = [atom(f"C{k}", 1 if k < 3 else 2, "A", "C") for k in range(6)]
    chain_bonds = [(k, k + 1) for k in range(5)]
    s = micro("energy_chain", chain_atoms, [[[1.5 * k, 0, 0] for k in range(6)]], bonds=chain_bonds)
    for depth in (0, 1, 2, 3, 4):
        cases.append(run_case(f"energy_exclude_{depth}", s, base, [f"analyses.energy.exclude_bonded={depth}"],
                              ["energy"]))

    # RMSD: superposition, fit selection, none, masses
    rng = np.random.default_rng(3)
    x = rng.normal(scale=3.0, size=(12, 3))
    y = x @ euler_to_matrix([40, 10, -30]).T + [5, 0, 0]
    s = micro("rmsd", [atom(f"C{k}", 1 + k // 6, "A", "C") for k in range(12)], [x, y])
    target = "analyses.rmsd.targets.0.selection=group:Aptamer or group:Target"
    cases.append(run_case("rmsd_superposition", s, base, [target], ["rmsd"]))
    cases.append(run_case("rmsd_fit_selection", s, base,
                          [target, "analyses.rmsd.targets.0.fit_selection=group:Aptamer"], ["rmsd"]))
    cases.append(run_case("rmsd_no_superposition", s, base, [target, "analyses.rmsd.superposition=false"], ["rmsd"]))
    cases.append(run_case("rmsd_mass_weighted", s, base, [target, "analyses.rmsd.mass_weighted=true"], ["rmsd"]))
    cases.append(run_case("rmsd_instance", s, base, [], ["rmsd@Aptamer", "rmsd@Target"]))

    # base pairs (test_basepairs.py): duplex G1-C4, A2-U3; C4 leaves in frames 1-2
    duplex = tb.duplex_frames()
    s = micro_from_frames("duplex", duplex)
    rna = micro_base({"RNA": "1-4"})
    cases.append(run_case("duplex_stride_1", s, rna, ["analyses.basepairs.stride=1"], ["basepairs", "trace"]))
    cases.append(run_case("duplex_stride_2", s, rna, ["analyses.basepairs.stride=2"], ["basepairs"]))
    cases.append(run_case("duplex_no_wobble_cutoff", s, rna, ["analyses.basepairs.stride=1",
                                                              "analyses.basepairs.wobble=false",
                                                              "analyses.basepairs.cutoff=3.0"], ["basepairs"]))
    cases.append(run_case("duplex_all", s, micro_base({"Aptamer": "1-2", "Target": "3-4"}),
                          ["analyses.basepairs.stride=1"], ALL))
    frames = [f + [("ALA", {"N1": np.array([0.0, 0.0, 50.0]), "N3": np.array([1.0, 0.0, 50.0])})]
              for f in tb.duplex_frames()]
    s = micro_from_frames("duplex_ala", frames)
    other = micro_base({"RNA": "1-3", "Other": "5"}, **{"visual.structure_groups": ["RNA", "Other"]})
    cases.append(run_case("structure_residues_only", s, other, [], ["basepairs", "trace", "atoms"]))
    cases.append(run_case("proximity_unknown_set", s, other, [], ["proximity"]))
    g = tb.nucleotide("G")
    c = tb.partner(g, "G", "C")
    s = micro_from_frames("l_rna", [[("LG5", g), ("0C", c)]] * 2)
    cases.append(run_case("l_rna", s, micro_base({"RNA": "1-2"}), [], ["basepairs"]))
    # wobble pair and a pair beyond the cutoff (test_wobble_optional_and_cutoff)
    u_ = tb.partner(g, "G", "U")
    far = tb.partner(g, "G", "U", d=3.8)
    s = micro_from_frames("wobble", [[("G", g), ("U", u_)], [("G", g), ("U", far)]])
    for name, extra in (("wobble_default", []), ("wobble_off", ["analyses.basepairs.wobble=false"]),
                        ("wobble_cutoff_4", ["analyses.basepairs.cutoff=4.0"])):
        cases.append(run_case(name, s, micro_base({"RNA": "1-2"}), ["analyses.basepairs.stride=1"] + extra,
                              ["basepairs"]))

    # energy: bonded exclusions, the function itself (keys a_local * n_b + b_local)
    u = conftest.make_universe(chain_atoms, np.asarray([[[1.5 * k, 0, 0] for k in range(6)]], np.float32),
                               bonds=chain_bonds)
    from mdatk.analysis.energy import _bonded_exclusions
    exclusions = [{"depth": depth, "keys": [int(v) for v in _bonded_exclusions(u, u.atoms[:3], u.atoms[3:], depth)]}
                  for depth in range(6)]
    return {"systems": MICRO_SYSTEMS, "cases": cases,
            "bonded_exclusions": {"system": "energy_chain", "a": [0, 1, 2], "b": [3, 4, 5], "cases": exclusions}}


# --- the demo system --------------------------------------------------------------------------------

DEMO_TOP = os.path.join(DATA, "topology", "demo_dry.parm7")
DEMO_NC = os.path.join(DATA, "io", "trajectories", "demo_dry.nc")
DEMO_FILES = os.path.join(DATA, "io", "trajectories", "files")


def demo_base():
    """The project `mdatk demo` writes for demo_dry (build_demo(n_frames=24, solvent=False))."""
    return demo_project("demo_dry.parm7", "demo_dry.nc", "MDATK_output", 12, cleaning=False,
                        mid_time_ps=(24 // 2) * 10.0)


def demo_cases():
    system = file_system("demo_dry", DEMO_TOP, [DEMO_NC])
    cases = [run_case("default", system, demo_base, [], ALL)]
    return {"systems": {"demo_dry": {"topology": "{data}/topology/demo_dry.parm7",
                                     "trajectories": ["{data}/io/trajectories/demo_dry.nc"]}},
            "base_project": "mdatk.demo.demo_project('demo_dry.parm7', 'demo_dry.nc', 'MDATK_output', 12, "
                            "cleaning=False, mid_time_ps=120.0)",
            "cases": cases}


def main():
    os.makedirs(OUT, exist_ok=True)
    write_json("micro_cases.json", micro_cases())
    write_json("demo_dry_cases.json", demo_cases())
    remove_offsets()


if __name__ == "__main__":
    main()
