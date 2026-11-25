"""What MDATK's spatial helpers give (mdatk/spatial.py), for tests/analysis/test_spatial.cpp.

    /opt/mdbridge-venv/bin/python -I -B make_spatial_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy and scipy (mdatk.spatial imports them); nothing is written into the checkout. MDATK's
functions run in-process. Writes OUT_DIR/spatial_extra_cases.json:

  pairs        pairs_within(a, b, cutoff) without a box (scipy's cKDTree path): the inputs and the
               (i, j, d) MDATK returns, in its order (cKDTree traversal; the C++ test compares the
               pairs sorted by (i, j)), or the error. Random float64 and float32-widened points,
               clustered and duplicated points, the same array twice, crafted pairs at exactly the
               cutoff and pairs where sqrt(d2) <= r and d2 <= r * r disagree, a sparse spread that
               defeats a cell grid, cutoffs 0 / negative / infinite / NaN, non-finite points
  min_image    minimum_image_vectors(vectors, box) for boxes with and without 90-degree angles
               (numpy's allclose margins), unusable boxes, half-box ties (np.round: half to even):
               MDATK's result ("out") and MDBridge's improved one ("improved", minimum images for
               every cell shape, one float32 orthogonality test; docs/DIFFERENCES.md AN-04,
               improvements_periodic.minimum_image_vectors), which the C++ test expects
  usable_box   _usable_box(box): the box, or None
  kabsch       kabsch(mobile, reference, weights): R, mobile centre, reference centre; random,
               weighted, reflected, coplanar (rank 2), fewer than 3 points
  rings        mdatk.analysis.stacking.ring_geometry of (n, k, 3) rings (k = 5, 6, 9): centroids,
               normals and planarity (np.linalg.svd, batched)
  sums         np.sum / np.mean of float64 arrays of 0 to 300 and some larger lengths (numpy's
               pairwise summation); the values are value(k) of the formula below, which the C++
               test computes itself
  tests        the inputs of MDATK's test_pairs_within_* / test_kabsch_* (tests/test_analyses.py)

Every float is float.hex() ("inf", "-inf", "nan" for the non-finite ones). kabsch and ring_geometry
use OpenBLAS / LAPACK, whose kernels are chosen per CPU, so those numbers may differ in the last bits
between machines (regenerate_all.py compares this file within a relative 1e-12 where its bytes
differ); everything else comes from scipy's C++ and numpy's own loops and is reproducible.
"""

import json
import math
import os
import sys
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_spatial_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import numpy as np  # noqa: E402

from mdatk import spatial as S  # noqa: E402
from mdatk.analysis.stacking import ring_geometry  # noqa: E402
from mdatk.render.orientation import euler_to_matrix  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import improvements_periodic  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
rng = np.random.default_rng(20261124)


def h(x):
    return float(x).hex()


def hexes(values):
    return [h(v) for v in np.asarray(values, dtype=np.float64).ravel()]


def value(k):
    """The summed values: exact IEEE operations only, so C++ computes the same doubles."""
    return (((k * 7919) % 1000) - 500) / 37.0 * (1.0 + k / 3.0)


# ---- pairs_within --------------------------------------------------------------------------------

def pairs_case(name, a, b, cutoff):
    a = np.asarray(a, dtype=np.float64).reshape(-1, 3)
    b = np.asarray(b, dtype=np.float64).reshape(-1, 3)
    case = {"name": name, "a": hexes(a), "b": hexes(b), "cutoff": h(cutoff)}
    try:
        i, j, d = S.pairs_within(a, b, cutoff)
    except Exception as exc:  # noqa: BLE001
        case["error"] = {"type": type(exc).__name__, "message": str(exc)}
        return case
    assert i.dtype == np.int64 and j.dtype == np.int64 and d.dtype == np.float64
    case["i"] = [int(x) for x in i]
    case["j"] = [int(x) for x in j]
    case["d"] = hexes(d)
    return case


def boundary_pairs(n):
    """Pairs (point, other) whose squared distance is within a few ulp of r * r: exactly at the
    cutoff, and where sqrt(d2) <= r and d2 <= r * r disagree (scipy compares squares)."""
    a, b, rs = [], [], []
    found_disagree = 0
    while len(a) < n:
        r = float(rng.uniform(1.0, 6.0))
        p = rng.uniform(-20, 20, 3)
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        q = p + u * r
        # nudge q along the line by a few ulp until both tests are evaluated on either side
        for step in range(-6, 7):
            qq = q + u * (step * 2.0 ** -50 * r)
            dv = p - qq
            d2 = ((0.0 + dv[0] * dv[0]) + dv[1] * dv[1]) + dv[2] * dv[2]
            disagree = (math.sqrt(d2) <= r) != (d2 <= r * r)
            if disagree or step in (-1, 0, 1):
                a.append(p)
                b.append(qq)
                rs.append(r)
                found_disagree += disagree
    return a, b, rs, found_disagree


def pairs_cases():
    cases = []
    for k, (na, nb, size, cutoff) in enumerate([(1, 1, 5, 3.0), (2, 3, 4, 2.5), (10, 10, 10, 4.0),
                                                 (17, 33, 12, 3.5), (60, 50, 20, 4.0), (100, 1, 15, 6.0),
                                                 (1, 100, 15, 6.0), (200, 150, 25, 3.2),
                                                 (250, 300, 40, 5.0)]):
        a = rng.uniform(0, size, (na, 3))
        b = rng.uniform(0, size, (nb, 3))
        cases.append(pairs_case(f"uniform float64 {k}", a, b, cutoff))
        cases.append(pairs_case(f"uniform float32 {k}", a.astype(np.float32), b.astype(np.float32),
                                cutoff))
    # coordinates as MD positions: float32, negative and offset
    a = (rng.normal(size=(300, 3)) * 8 - 40).astype(np.float32)
    cases.append(pairs_case("self float32", a, a, 3.0))
    cases.append(pairs_case("self cutoff 1.95", a, a, 1.95))
    b = (rng.normal(size=(120, 3)) * 8 - 40).astype(np.float32)
    cases.append(pairs_case("offset float32", a, b, 4.5))
    # clusters and duplicated points (many equal distances, zero distances)
    centres = rng.uniform(0, 30, (5, 3))
    a = np.repeat(centres, 12, axis=0) + rng.normal(size=(60, 3)) * 0.5
    a[10] = a[3]
    a[11] = a[3]
    b = np.vstack([a[::3], centres])
    cases.append(pairs_case("clusters", a, b, 2.0))
    grid = np.array([[x, y, z] for x in range(4) for y in range(4) for z in range(4)], dtype=np.float64)
    for cutoff in (1.0, math.sqrt(2.0), math.sqrt(3.0), 2.0, math.nextafter(1.0, 0.0), math.nextafter(1.0, 2.0)):
        cases.append(pairs_case(f"lattice {cutoff!r}", grid, grid * 1.0, cutoff))
    cases.append(pairs_case("lattice 0.5 spacing", grid * 0.5, grid, 0.5))
    # exact boundaries
    a, b, rs, disagree = boundary_pairs(60)
    assert disagree > 0, "no sqrt / square disagreement found"
    for k in range(0, len(a), 3):
        cases.append(pairs_case(f"boundary {k}", a[k:k + 3], b[k:k + 3], rs[k]))
    # a sparse spread (a cell grid over its extent would be huge) and far-apart sets
    a = rng.uniform(-1e6, 1e6, (40, 3))
    b = a[:20] + rng.normal(size=(20, 3))
    cases.append(pairs_case("spread", a, b, 2.5))
    cases.append(pairs_case("far apart", rng.uniform(0, 5, (10, 3)), rng.uniform(1000, 1005, (10, 3)), 4.0))
    cases.append(pairs_case("large coordinates", rng.uniform(1e15, 1e15 + 64, (20, 3)),
                            rng.uniform(1e15, 1e15 + 64, (20, 3)), 30.0))
    # special cutoffs and inputs
    a = rng.uniform(0, 6, (8, 3))
    b = rng.uniform(0, 6, (9, 3))
    for cutoff in (0.0, -0.0, -2.5, 2.5, float("inf"), float("-inf"), float("nan"), 1e-300, 1e300):
        cases.append(pairs_case(f"cutoff {cutoff!r}", a, b, cutoff))
    cases.append(pairs_case("zero cutoff duplicates", a, np.vstack([a[2], a[5], a[2]]), 0.0))
    cases.append(pairs_case("empty a", np.zeros((0, 3)), b, 2.0))
    cases.append(pairs_case("empty b", a, np.zeros((0, 3)), 2.0))
    cases.append(pairs_case("empty both", np.zeros((0, 3)), np.zeros((0, 3)), 2.0))
    for bad in (float("nan"), float("inf"), float("-inf")):
        c = a.copy()
        c[3, 1] = bad
        cases.append(pairs_case(f"a has {bad!r}", c, b, 2.0))
        cases.append(pairs_case(f"b has {bad!r}", b, c, 2.0))
        cases.append(pairs_case(f"empty b, a has {bad!r}", c, np.zeros((0, 3)), 2.0))
    return cases


# ---- minimum_image_vectors / _usable_box -------------------------------------------------------

def box_value(box):
    if box is None:
        return None
    return hexes(box)


def min_image_cases():
    cases = []
    vecs = rng.uniform(-30, 30, (12, 3))
    vecs = np.vstack([vecs, [[5.0, -5.0, 15.0], [25.0, -15.0, 0.5], [-0.0, 0.0, 10.0],
                             [np.nan, np.inf, -np.inf]]])
    boxes = [None, [10.0, 10.0, 10.0, 90.0, 90.0, 90.0], [10.0, 20.0, 30.0, 90.0, 90.0, 90.0],
             [7.3, 11.1, 9.9, 90.0, 90.0, 90.0], [10.0, 10.0, 10.0, 90.0 + 1e-6, 90.0 - 1e-6, 90.0],
             [10.0, 10.0, 10.0, 90.0009, 90.0, 90.0], [10.0, 10.0, 10.0, 90.00090001, 90.0, 90.0],
             [10.0, 10.0, 10.0, 90.0011, 90.0, 90.0], [10.0, 10.0, 10.0, 60.0, 70.0, 80.0],
             [10.0, 10.0, 10.0, np.nan, 90.0, 90.0], [10.0, 10.0, 10.0, np.inf, 90.0, 90.0],
             [0.0, 10.0, 10.0, 90.0, 90.0, 90.0], [-10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
             [np.nan, 10.0, 10.0, 90.0, 90.0, 90.0], [np.inf, 10.0, 10.0, 90.0, 90.0, 90.0],
             [10.0, 10.0, 10.0], [10.0, 10.0, 10.0, 90.0, 90.0, 90.0, 1.0],
             [1e-300, 10.0, 10.0, 90.0, 90.0, 90.0]]
    for box in boxes:
        arr = None if box is None else np.asarray(box, dtype=np.float64)
        out = S.minimum_image_vectors(vecs, arr)
        improved = improvements_periodic.minimum_image_vectors(vecs, arr) if arr is not None else vecs
        cases.append({"vectors": hexes(vecs), "box": box_value(box), "out": hexes(out),
                      "improved": hexes(improved)})
    # float32 box values (ts.dimensions) widened to float64
    box32 = np.array([23.7, 31.3, 27.1, 90.0, 90.0, 90.0], dtype=np.float32)
    cases.append({"vectors": hexes(vecs), "box": box_value(box32.astype(np.float64)),
                  "out": hexes(S.minimum_image_vectors(vecs, box32)),
                  "improved": hexes(improvements_periodic.minimum_image_vectors(vecs, box32))})
    return cases


def usable_box_cases():
    cases = []
    for box in ([10.0, 10.0, 10.0, 90.0, 90.0, 90.0], [10.0, 10.0, 10.0, 60.0, 70.0, 80.0],
                [0.0, 10.0, 10.0, 90.0, 90.0, 90.0], [10.0, -0.0, 10.0, 90.0, 90.0, 90.0],
                [10.0, 10.0, -1.0, 90.0, 90.0, 90.0], [np.nan, 10.0, 10.0, 90.0, 90.0, 90.0],
                [np.inf, 10.0, 10.0, 90.0, 90.0, 90.0], [10.0, 10.0, 10.0, 0.0, 0.0, 0.0],
                [10.0, 10.0, 10.0, np.nan, np.nan, np.nan], [10.0, 10.0, 10.0], [],
                [10.0, 10.0, 10.0, 90.0, 90.0, 90.0, 0.0], [5e-324, 1.0, 1.0, 90.0, 90.0, 90.0]):
        out = S._usable_box(np.asarray(box, dtype=np.float64))
        cases.append({"box": hexes(box), "out": None if out is None else hexes(out)})
    assert S._usable_box(None) is None
    return cases


# ---- kabsch ------------------------------------------------------------------------------------

def kabsch_case(name, mobile, reference, weights=None):
    mobile = np.asarray(mobile, dtype=np.float64).reshape(-1, 3)
    reference = np.asarray(reference, dtype=np.float64).reshape(-1, 3)
    r, mc, rc = S.kabsch(mobile, reference, weights)
    return {"name": name, "mobile": hexes(mobile), "reference": hexes(reference),
            "weights": None if weights is None else hexes(weights),
            "r": hexes(r), "mc": hexes(mc), "rc": hexes(rc)}


def random_rotation():
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def kabsch_cases():
    cases = []
    for n in (3, 4, 5, 8, 20, 50, 100):
        x = rng.normal(size=(n, 3)) * 5
        r = random_rotation()
        y = x @ r.T + rng.normal(size=3) * 10
        cases.append(kabsch_case(f"exact rotation {n}", x, y))
        cases.append(kabsch_case(f"noisy {n}", x, y + rng.normal(size=y.shape) * 0.3))
        w = rng.uniform(0.5, 16.0, n)
        cases.append(kabsch_case(f"weighted {n}", x, y + rng.normal(size=y.shape) * 0.3, w))
        cases.append(kabsch_case(f"unrelated {n}", x, rng.normal(size=(n, 3)) * 5))
        mirror = y * np.array([1.0, 1.0, -1.0])
        cases.append(kabsch_case(f"mirrored {n}", x, mirror + rng.normal(size=y.shape) * 0.1))
    for n in (3, 6):  # coplanar: rank 2
        x = np.column_stack([rng.normal(size=(n, 2)) * 4, np.zeros(n)])
        r = random_rotation()
        cases.append(kabsch_case(f"coplanar {n}", x, x @ r.T + 1.0))
        cases.append(kabsch_case(f"coplanar noisy {n}", x, x @ r.T + rng.normal(size=(n, 3)) * 0.2))
    x = rng.normal(size=(10, 3)) * 1e-3
    cases.append(kabsch_case("tiny", x, x @ random_rotation().T))
    x = rng.normal(size=(10, 3)) * 1e4 + 1e5
    cases.append(kabsch_case("large offset", x, x @ random_rotation().T - 3e5))
    x = rng.normal(size=(30, 3)).astype(np.float32).astype(np.float64) * 10
    cases.append(kabsch_case("float32 points", x, (x @ random_rotation().T).astype(np.float32)))
    w = np.zeros(10)
    w[:4] = [1.0, 2.0, 3.0, 4.0]
    x = rng.normal(size=(10, 3))
    cases.append(kabsch_case("zero weights", x, x @ random_rotation().T + 2.0, w))
    for n in (0, 1, 2):
        x = rng.normal(size=(n, 3))
        cases.append(kabsch_case(f"{n} points", x, x + 1.0))
        if n:
            cases.append(kabsch_case(f"{n} points weighted", x, x + 1.0, rng.uniform(1, 2, n)))
    x = rng.normal(size=(5, 3))
    cases.append(kabsch_case("identical", x, x))
    return cases


def ring_cases():
    cases = []
    for k in (5, 6, 9):
        n = 7
        rings = []
        for m in range(n):
            ang = np.arange(k) * 2 * np.pi / k + rng.uniform(0, 1)
            flat = np.column_stack([np.cos(ang) * 1.4, np.sin(ang) * 1.4, rng.normal(size=k) * 0.05 * m])
            rings.append(flat @ random_rotation().T + rng.uniform(-20, 20, 3))
        pos = np.array(rings)
        c, normals, planarity = ring_geometry(pos)
        cases.append({"k": k, "positions": hexes(pos), "centroids": hexes(c), "normals": hexes(normals),
                      "planarity": hexes(planarity)})
    return cases


def sum_cases():
    cases = []
    lengths = list(range(0, 301)) + [511, 512, 513, 1000, 1023, 1024, 1025, 4096, 4097, 10007, 65536]
    for n in lengths:
        x = np.array([value(k) for k in range(n)], dtype=np.float64)
        cases.append({"n": n, "sum": h(np.sum(x)), "mean": h(np.mean(x)) if n else h(float("nan"))})
    return cases


def test_inputs():
    out = {}
    rng1 = np.random.default_rng(1)
    a = rng1.uniform(0, 20, (60, 3))
    b = rng1.uniform(0, 20, (50, 3))
    out["pairs_within_matches_brute_force"] = {"a": hexes(a), "b": hexes(b), "cutoff": h(4.0)}
    rng2 = np.random.default_rng(2)
    x = rng2.normal(size=(20, 3))
    r = euler_to_matrix([30, -20, 75])
    y = x @ r.T + np.array([1.0, 2.0, 3.0])
    out["kabsch_recovers_rotation"] = {"x": hexes(x), "y": hexes(y)}
    return out


def dump(data):
    parts = []
    for key, cases in data.items():
        if isinstance(cases, list):
            rows = ",\n".join("  " + json.dumps(c, separators=(",", ":")) for c in cases)
            parts.append(f" {json.dumps(key)}: [\n{rows}\n ]")
        else:
            parts.append(f" {json.dumps(key)}: {json.dumps(cases, separators=(',', ':'))}")
    return "{\n" + ",\n".join(parts) + "\n}\n"


def main():
    data = {
        "pairs": pairs_cases(),
        "min_image": min_image_cases(),
        "usable_box": usable_box_cases(),
        "kabsch": kabsch_cases(),
        "rings": ring_cases(),
        "sums": sum_cases(),
        "tests": test_inputs(),
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "spatial_extra_cases.json"), "w", encoding="utf-8", newline="\n") as f:
        f.write(dump(data))


if __name__ == "__main__":
    main()
