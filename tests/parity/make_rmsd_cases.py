"""Probe cases for docs/specs/m3-rmsd.md: what MDAnalysis 2.11.0's RMSD function gives for the call
MDATK's rmsd analysis makes (weights None or float64 masses, centring and superposition on, float64
(n, 3) coordinates), for the C++ tests of MDBridge's own implementation.

    /opt/mdbridge-venv/bin/python -I -B make_rmsd_cases.py MDATK_CHECKOUT OUT_DIR [--report]

The MDATK checkout is not used (the argument keeps the generators' common calling convention);
nothing is written into it. Writes OUT_DIR/rmsd_cases.json:

  cases   one probe per entry: "a" (first argument, the moving frame in MDATK), "b" (second, the
          reference), both as flat lists of 3n coordinates; "weights" (null or n values); the
          result ("result", a float, or null) or the error ("error": {"type", "message"}); the
          warnings numpy emitted ("warnings": "Category: message"); "scale" = E0 / N, the mean
          weighted squared distance of the centred inputs from their centres (see the spec; the
          tolerance of non-degenerate cases is relative to it); "degenerate": true for inputs
          whose optimal rotation is not unique (one or two effective atoms, collinear or coincident
          atoms); "group" names the family. Where MDBridge improves on MDAnalysis (AN-16: degenerate
          sets; AN-17: E0 - lambda lost to cancellation; improvements_analyses.rmsd_expected),
          "result" is MDBridge's value (the exact superposition RMSD, 0 for identical inputs or one
          atom), "improved" names the item, "mdanalysis" keeps MDAnalysis' value and "near": true
          marks an AN-17 decision close to its threshold (either value is acceptable).
  demo    the series of rmsd(frame, reference) over tests/data/io/trajectories/demo_dry.nc with
          topology/demo_dry.parm7 (24 frames), as MDATK computes it: per target the atom indices
          ([first, last] runs), the weights (null or the parm7 masses), the reference frame and the
          value for every frame (MDBridge's, as for the cases; "improved": the frame positions
          whose value is the exact one, "near": those close to the AN-17 threshold, "mdanalysis":
          MDAnalysis' values).

Floats are float.hex() strings ("nan", "inf", "-inf" for non-finite values). Everything recorded
is deterministic (MDAnalysis' sums are sequential C loops and numpy's pairwise sums are fixed);
nothing goes through LAPACK or BLAS, whose kernels are chosen per CPU and round differently (the
random rotations come from quaternions and are applied elementwise, the exact superposition RMSD
of the improved values is computed in plain floats: improvements_analyses.exact_rmsd), so the
file does not depend on the CPU.

--report also prints the spread between MDAnalysis and an exact superposition (Horn's quaternion,
the residual computed directly after rotating) per family: the evidence for the spec's tolerances.
"""

import json
import math
import os
import sys
import warnings

import numpy as np

if len(sys.argv) not in (3, 4) or (len(sys.argv) == 4 and sys.argv[3] != "--report"):
    sys.exit("usage: make_rmsd_cases.py MDATK_CHECKOUT OUT_DIR [--report]")
OUT = os.path.abspath(sys.argv[2])
REPORT = len(sys.argv) == 4
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

sys.path.insert(0, HERE)
import improvements_analyses as improved  # noqa: E402
import independent_tables  # noqa: E402

import MDAnalysis as mda  # noqa: E402
from MDAnalysis.analysis import rms  # noqa: E402

rng = np.random.default_rng(20261124)
MASSES = [1.008, 12.011, 14.007, 15.999, 30.974, 24.305, 22.99, 32.06]


def fhex(x):
    return float(x).hex()


def flat(x):
    return [fhex(v) for v in np.asarray(x, dtype=np.float64).ravel()]


def runs(indices):
    out = []
    for i in indices:
        if out and out[-1][1] + 1 == i:
            out[-1][1] = int(i)
        else:
            out.append([int(i), int(i)])
    return out


def scale_of(a, b, w):
    """E0 / N of the centred inputs, with weights normalised to mean 1 (numpy, deterministic)."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    n = len(a)
    if n == 0:
        return float("nan")
    with np.errstate(all="ignore"):
        w = np.ones(n) if w is None else np.asarray(w, dtype=np.float64)
        ac = a - (a * w[:, None]).sum(axis=0) / w.sum()
        bc = b - (b * w[:, None]).sum(axis=0) / w.sum()
        wn = w / w.mean()
        return float((np.sum(wn[:, None] * ac * ac) + np.sum(wn[:, None] * bc * bc)) / (2 * n))


def call(a, b, w):
    """rms.rmsd as MDATK calls it: result, error, warnings."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            value = rms.rmsd(a, b, weights=w, center=True, superposition=True)
            result, error = value, None
        except Exception as exc:  # noqa: BLE001
            result, error = None, {"type": type(exc).__name__, "message": str(exc)}
    msgs = [f"{c.category.__name__}: {c.message}" for c in caught]
    return result, error, msgs


CASES = []
SPREAD = {}


def case(group, a, b, w=None, degenerate=False):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    wa = None if w is None else np.asarray(w, dtype=np.float64)
    result, error, msgs = call(a, b, wa)
    if result is not None:
        assert type(result) is float, type(result)
    entry = {"name": f"{group}_{len(CASES)}", "group": group, "n": int(a.shape[0]) if a.ndim == 2 else None,
             "a": flat(a), "b": flat(b), "a_shape": list(a.shape), "b_shape": list(b.shape),
             "weights": None if wa is None else flat(wa),
             "weights_shape": None if wa is None else list(wa.shape),
             "result": None if result is None else fhex(result), "error": error, "warnings": msgs,
             "scale": fhex(scale_of(a, b, wa)) if a.shape == b.shape and a.ndim == 2 and
             (wa is None or wa.shape == (a.shape[0],)) and error is None else None,
             "degenerate": degenerate}
    if error is None:
        value, item, near = improved.rmsd_expected(a, b, wa, result)
        if item is not None:
            entry.update({"result": fhex(value), "improved": item, "mdanalysis": fhex(result), "near": near})
    CASES.append(entry)
    if REPORT and result is not None and np.isfinite(result) and a.shape[0] > 0:
        try:
            ref = improved.exact_rmsd(a, b, wa)
        except ZeroDivisionError:
            return
        s = scale_of(a, b, wa)
        eps = np.finfo(float).eps
        stat = SPREAD.setdefault(group, [0, 0.0, 0.0, 0.0])
        stat[0] += 1
        stat[1] = max(stat[1], abs(result - ref))
        stat[2] = max(stat[2], abs(result * result - ref * ref) / (eps * s) if s > 0 else 0.0)
        if ref > 1e-3:
            stat[3] = max(stat[3], abs(result - ref) / ref)


def rotation():
    """A random rotation (applied to rows by rot), from a normalised random quaternion in plain
    floats: numpy's QR and matmul go to LAPACK / OpenBLAS, whose last bits depend on the CPU."""
    q0, q1, q2, q3 = (float(v) for v in rng.normal(size=4))
    s = math.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
    q0, q1, q2, q3 = q0 / s, q1 / s, q2 / s, q3 / s
    return np.array([
        [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q1 * q3 + q0 * q2)],
        [2.0 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2.0 * (q2 * q3 - q0 * q1)],
        [2.0 * (q1 * q3 - q0 * q2), 2.0 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3],
    ])


def rot(x, r):
    """x @ r without BLAS (improvements_analyses.rotate_rows)."""
    return improved.rotate_rows(x, r)


def as_f32(x):
    """Coordinates as MDATK has them: float32 positions widened to float64."""
    return np.asarray(x, dtype=np.float32).astype(np.float64)


def weights_for(n, kind):
    if kind == "none":
        return None
    if kind == "uniform":
        return rng.uniform(1.0, 32.0, size=n)
    if kind == "masses":
        return np.array([MASSES[k] for k in rng.integers(0, len(MASSES), size=n)])
    if kind == "equal":
        return np.full(n, 12.011)
    raise ValueError(kind)


def general_cases():
    kinds = ["none", "uniform", "masses", "equal"]
    noises = [0.0, 1e-9, 1e-6, 1e-3, 0.1, 1.0, 5.0]
    count = 0
    for n in [3, 4, 5, 7, 8, 9, 16, 17, 40]:
        for noise in noises:
            for kind in ("none", kinds[1 + count % 3]):
                count += 1
                a = as_f32(rng.normal(size=(n, 3)) * rng.choice([3.0, 15.0, 40.0]) + rng.normal(size=3) * 30)
                b = as_f32(rot(a, rotation()) + rng.normal(size=3) * 20 + rng.normal(size=(n, 3)) * noise)
                case(f"general_noise_{noise:g}", a, b, weights_for(n, kind))
    # beyond numpy's 128-element pairwise-summation block (the sum of the weights)
    for n, noise in [(129, 0.0), (129, 0.5), (300, 1e-6), (300, 2.0)]:
        a = as_f32(rng.normal(size=(n, 3)) * 15 + rng.normal(size=3) * 30)
        b = as_f32(rot(a, rotation()) + rng.normal(size=3) * 20 + rng.normal(size=(n, 3)) * noise)
        case(f"general_noise_{noise:g}", a, b, weights_for(n, "masses"))
    # identical inputs (MDATK's reference frame against itself), translated and rotated copies
    for n in [3, 4, 5, 9, 24, 100]:
        for kind in ["none", "masses"]:
            a = as_f32(rng.normal(size=(n, 3)) * 12 + rng.normal(size=3) * 25)
            w = weights_for(n, kind)
            case("identical", a, a.copy(), w)
            case("translated", a, a + np.array([1.5, -2.25, 40.0]), w)
            case("rotated", a, rot(a, rotation()), w)
    # coplanar sets (unique optimal rotation, a doubly degenerate inertia only)
    for n in [3, 4, 6, 12]:
        a = rng.normal(size=(n, 3)) * 8
        a[:, 2] = 1.25
        b = rot(a, rotation()) + rng.normal(size=(n, 3)) * 0.3
        case("planar", a, b)
        b = a + rng.normal(size=(n, 3)) * 0.3
        b[:, 2] = -3.0
        case("planar_both", a, rot(b, rotation()))
    # magnitudes
    for mag in [1e-6, 1e-3, 1e3, 1e6]:
        a = rng.normal(size=(10, 3)) * mag
        case("magnitude", a, rot(a, rotation()) + rng.normal(size=(10, 3)) * mag * 0.01)


def degenerate_cases():
    for k in range(4):
        a = rng.normal(size=(1, 3)) * 10
        case("one_atom", a, a + rng.normal(size=(1, 3)), [None, [12.011]][k % 2], degenerate=True)
    case("one_atom", [[1.0, 2.0, 3.0]], [[1.0, 2.0, 3.0]], degenerate=True)
    for k in range(16):
        a = as_f32(rng.normal(size=(2, 3)) * 5)
        w = None if k % 2 == 0 else weights_for(2, "masses")
        if k < 6:
            b = a.copy()
        elif k < 11:
            b = as_f32(rot(a, rotation()) * rng.choice([1.0, 1.05]) + rng.normal(size=3))
        else:
            b = as_f32(a + rng.normal(size=(2, 3)) * 0.5)
        case("two_atoms", a, b, w, degenerate=True)
    for k in range(12):
        n = int(rng.integers(3, 12))
        d = rng.normal(size=3)
        a = np.outer(rng.normal(size=n) * 6, d) + rng.normal(size=3) * 4
        if k < 4:
            b = rot(a, rotation()) * 1.1
        elif k < 8:
            b = rot(a, rotation()) + rng.normal(size=(n, 3)) * 0.3
        else:
            b = a.copy()
        case("collinear", a, b, degenerate=True)
    for k in range(4):  # collinear in float32 (MDATK's inputs) up to rounding
        n = 4 + k
        a = as_f32(np.outer(np.arange(n) * 1.5, [0.3, -0.7, 0.2]) + 10)
        case("collinear_f32", a, as_f32(rot(a, rotation()) + 2), degenerate=True)
    for n in [2, 5]:  # coincident atoms
        a = np.tile([[3.0, -1.0, 7.5]], (n, 1))
        case("coincident", a, a.copy(), degenerate=True)
        case("coincident", a, a + rng.normal(size=(n, 3)), degenerate=True)
    # weights that leave one or two atoms
    n = 8
    a = rng.normal(size=(n, 3)) * 6
    b = rot(a, rotation()) + rng.normal(size=(n, 3)) * 0.2
    for w in ([0, 0, 0, 0, 0, 0, 0, 1.0], [0, 0, 5.0, 0, 0, 0, 0, 0], [1.0, 0, 0, 0, 0, 0, 0, 2.0]):
        case("weights_few", a, b, w, degenerate=True)


def special_cases():
    n = 6
    a = rng.normal(size=(n, 3)) * 6
    b = rot(a, rotation()) + rng.normal(size=(n, 3)) * 0.2
    for where, value in [("a", np.nan), ("b", np.nan), ("a", np.inf), ("b", -np.inf)]:
        x, y = a.copy(), b.copy()
        (x if where == "a" else y)[2, 1] = value
        case(f"nonfinite_{where}", x, y)
    case("nonfinite_weight", a, b, [1, 2, np.nan, 1, 1, 1])
    case("nonfinite_weight", a, b, [1, 2, np.inf, 1, 1, 1])
    case("negative_weights", a, b, [-1.0] * n)
    case("negative_weights", a, b, [3.0, -1.0, 2.0, 1.0, 1.0, 1.0])
    case("zero_sum_weights", a, b, [0.0] * n)
    case("zero_sum_weights", a, b, [1.0, -1.0, 2.0, -2.0, 0.5, -0.5])
    case("tiny_weights", a, b, [1e-300] * n)
    case("huge_weights", a, b, [1e300] * n)
    case("empty", np.zeros((0, 3)), np.zeros((0, 3)))
    case("empty", np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0))


def error_cases():
    a = rng.normal(size=(5, 3))
    case("error_shape", a, a[:4])
    case("error_shape", a[:4], a)
    case("error_shape", a, np.zeros((5, 2)))
    case("error_weights", a, a, [1.0] * 4)
    case("error_weights", a, a, [1.0] * 6)
    case("error_weights", a, a, np.ones((5, 1)))
    case("error_weights", a, a, [])


def demo_cases():
    independent_tables.apply()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(os.path.join(DATA, "topology", "demo_dry.parm7"),
                         os.path.join(DATA, "io", "trajectories", "demo_dry.nc"))
    targets = []
    for label, sel in [("heavy", "not name H*"), ("all", "all"), ("C4'", "name C4'"), ("P", "name P"),
                       ("strand1 heavy", "resid 1:12 and not name H*"), ("Mg", "resname MG"),
                       ("three C4'", "resid 1:3 and name C4'"), ("two C4'", "resid 1:2 and name C4'"),
                       ("one C4'", "resid 1 and name C4'"), ("one Mg", "resid 25")]:
        ag = u.select_atoms(sel)
        for weighted in (False, True):
            for ref_frame in (0, 5, 23):
                if ref_frame != 0 and label not in ("heavy", "C4'", "two C4'"):
                    continue
                w = ag.masses.astype(np.float64) if weighted else None
                u.trajectory[ref_frame]
                ref = ag.positions.astype(np.float64).copy()
                values, mda_values, better, near = [], [], [], []
                for ts in u.trajectory:
                    x = ag.positions.astype(np.float64)
                    value = rms.rmsd(x, ref, weights=w, center=True, superposition=True)
                    mda_values.append(fhex(value))
                    value, item, close = improved.rmsd_expected(x, ref, w, value)
                    if item is not None:
                        better.append(ts.frame)
                    if close:
                        near.append(ts.frame)
                    values.append(fhex(value))
                    if REPORT:
                        case_stat("demo_" + label, x, ref, w)
                targets.append({"label": label, "selection": sel, "atoms": runs(ag.indices),
                                "weights": None if w is None else flat(w), "reference_frame": ref_frame,
                                "values": values, "improved": better, "near": near,
                                "mdanalysis": mda_values,
                                "degenerate": label in ("two C4'", "one C4'", "one Mg")})
    return {"topology": "topology/demo_dry.parm7", "trajectory": "io/trajectories/demo_dry.nc",
            "n_frames": int(u.trajectory.n_frames), "targets": targets}


def case_stat(group, x, ref, w):
    value = rms.rmsd(x, ref, weights=w, center=True, superposition=True)
    if not np.isfinite(value) or len(x) < 3:
        return
    r = exact_rmsd(x, ref, w)
    s = scale_of(x, ref, w)
    eps = np.finfo(float).eps
    stat = SPREAD.setdefault(group, [0, 0.0, 0.0, 0.0])
    stat[0] += 1
    stat[1] = max(stat[1], abs(value - r))
    stat[2] = max(stat[2], abs(value * value - r * r) / (eps * s) if s > 0 else 0.0)
    if r > 1e-3:
        stat[3] = max(stat[3], abs(value - r) / r)


def main():
    general_cases()
    degenerate_cases()
    special_cases()
    error_cases()
    demo = demo_cases()
    os.makedirs(OUT, exist_ok=True)
    head = {"generator": "tests/parity/make_rmsd_cases.py", "mdanalysis": mda.__version__,
            "call": "MDAnalysis.analysis.rms.rmsd(a, b, weights=weights, center=True, superposition=True)"}
    lines = ["{"]
    for key, value in head.items():
        lines.append(f" {json.dumps(key)}: {json.dumps(value)},")
    lines.append(' "cases": [')
    lines.append(",\n".join("  " + json.dumps(c, ensure_ascii=False) for c in CASES))
    lines.append(" ],")
    lines.append(' "demo": {')
    for key in ("topology", "trajectory", "n_frames"):
        lines.append(f"  {json.dumps(key)}: {json.dumps(demo[key])},")
    lines.append('  "targets": [')
    lines.append(",\n".join("   " + json.dumps(t, ensure_ascii=False) for t in demo["targets"]))
    lines.append("  ]")
    lines.append(" }")
    lines.append("}")
    with open(os.path.join(OUT, "rmsd_cases.json"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"{len(CASES)} cases, {len(demo['targets'])} demo targets")
    if REPORT:
        print("family: count, max |mda - svd|, max |mda^2 - svd^2| / (eps * E0/N), max relative (svd > 1e-3)")
        for group, (count, dabs, dsq, drel) in sorted(SPREAD.items()):
            print(f"  {group:<28} {count:5d}  {dabs:.3e}  {dsq:8.1f}  {drel:.3e}")


main()
