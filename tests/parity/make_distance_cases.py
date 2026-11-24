"""What MDAnalysis 2.11.0's capped_distance gives for the calls MDATK's pairs_within makes with a
periodic box (mdatk/spatial.py), for MDBridge's periodic pair search (docs/specs/m3-capped-distance.md).

    /opt/mdbridge-venv/bin/python -I -B make_distance_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy, scipy (mdatk.spatial imports scipy) and MDAnalysis 2.11.0; nothing is written into the
checkout. Every call is the one MDATK makes: both point sets cast to float32, the cutoff a Python
float, the box a float64 [a, b, c, alpha, beta, gamma] with a, b, c > 0, the default (automatic)
method, return_distances=True. Each case is also run through MDATK's own pairs_within, which must
give the same pairs and distances (checked here, not written). Writes OUT_DIR/distance_cases.json:

  meta       versions and the band formula below
  cases      one object per call:
    name, group      a label and the group of the spec it exercises
    a, b             the two point sets (float32 values as JSON numbers, exact; NaN / Infinity as
                     json writes them), or
    lattice_a, lattice_b
                     {"origin": [x, y, z], "spacing": s, "counts": [nx, ny, nz]}: point k (k from 0
                     to nx*ny*nz - 1) is origin + s * (k % nx, (k // nx) % ny, k // (nx * ny)),
                     computed in float64 and rounded to float32
    box, cutoff      as MDATK passes them (JSON numbers, exact)
    method           the method the automatic rule chose: "bruteforce" or "grid"
    pairs            [[i, j], ...] in MDAnalysis' order
    distances        float.hex of each returned float64 distance, same order
    or, for lattice cases, count, pairs_sha256 (SHA-256 of the pairs as little-endian int64, row by
                     row) and distances_sha256 (of the distances as little-endian float64)
    error            {"type": ..., "message": ...} instead of pairs when MDAnalysis raises
    warnings         messages of the warnings MDAnalysis emitted (usually none)
    band             the pairs [i, j] (returned or not) whose distance lies within the band
                     |d - cutoff| <= 8 * ulp32(M) of the cutoff, where d is the brute-force
                     minimum-image distance of the pair and M the largest of the pair's absolute
                     coordinates and the three box lengths (ulp32(M) = spacing of float32 at M):
                     an implementation that matches the distances only within the spec's tolerance
                     may disagree with MDAnalysis on exactly these pairs. Absent when empty or for
                     lattice cases.

Inputs MDAnalysis handles inconsistently (non-finite coordinates with a triclinic box: garbage pairs
or a crash) are not recorded; the spec says what MDBridge does instead.
"""

import hashlib
import json
import math
import os
import sys
import warnings

if len(sys.argv) != 3:
    sys.exit("usage: make_distance_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import numpy as np  # noqa: E402
import MDAnalysis  # noqa: E402
import MDAnalysis.lib.distances as mdadist  # noqa: E402
from MDAnalysis.lib.distances import capped_distance, distance_array  # noqa: E402

from mdatk import spatial  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
rng = np.random.default_rng(20261124)
f32 = np.float32

# record which method the automatic rule picks (test tooling only)
_chosen = []
_determine = mdadist._determine_method


def _spy(*args, **kwargs):
    function = _determine(*args, **kwargs)
    _chosen.append(function.__name__)
    return function


mdadist._determine_method = _spy
METHOD_NAMES = {"_bruteforce_capped": "bruteforce", "_nsgrid_capped": "grid"}


def num(x):
    return float(x)


def points_json(p):
    return [[num(v) for v in row] for row in np.asarray(p, dtype=np.float64)]


def lattice(spec):
    nx, ny, nz = spec["counts"]
    k = np.arange(nx * ny * nz)
    idx = np.stack([k % nx, (k // nx) % ny, k // (nx * ny)], axis=1).astype(np.float64)
    return (np.asarray(spec["origin"], dtype=np.float64) + spec["spacing"] * idx).astype(f32)


def ulp32(x):
    return np.spacing(np.abs(np.asarray(x, dtype=np.float64)).astype(f32)).astype(np.float64)


def band_pairs(a, b, box, cutoff):
    if not (len(a) and len(b)) or not math.isfinite(cutoff):
        return []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d = distance_array(a.copy(), b.copy(), box=np.asarray(box, dtype=np.float64))
    with np.errstate(invalid="ignore"):
        ma = np.where(np.isfinite(a), np.abs(a), 0).max(axis=1)
        mb = np.where(np.isfinite(b), np.abs(b), 0).max(axis=1)
        m = np.maximum(np.maximum(ma[:, None], mb[None, :]), max(abs(float(x)) for x in box[:3]))
        hit = np.abs(d - cutoff) <= 8 * ulp32(m)
    return [[int(i), int(j)] for i, j in zip(*np.nonzero(hit))]


def call(a, b, cutoff, box):
    _chosen.clear()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            pairs, dist = capped_distance(np.asarray(a, np.float32), np.asarray(b, np.float32),
                                          max_cutoff=cutoff, box=np.asarray(box, dtype=np.float64),
                                          return_distances=True)
            result = (pairs, dist, None)
        except Exception as exc:  # noqa: BLE001 - recorded
            result = (None, None, {"type": type(exc).__name__, "message": str(exc)})
    method = METHOD_NAMES[_chosen[0]] if _chosen else None
    msgs = [str(w.message) for w in caught]
    return result, method, msgs


def check_mdatk(a, b, cutoff, box, pairs, dist, error):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            i, j, d = spatial.pairs_within(np.asarray(a, np.float32), np.asarray(b, np.float32), cutoff,
                                           list(box))
        except Exception as exc:  # noqa: BLE001
            assert error is not None and str(exc) == error["message"], (exc, error)
            return
    assert error is None
    assert i.dtype == np.int64 and j.dtype == np.int64 and d.dtype == np.float64
    assert np.array_equal(i, pairs[:, 0]) and np.array_equal(j, pairs[:, 1])
    assert np.array_equal(d.view(np.int64), dist.view(np.int64))


CASES = []


def add(name, group, a, b, cutoff, box, band=True):
    a = np.asarray(a, dtype=np.float64).astype(f32).reshape(-1, 3)
    b = np.asarray(b, dtype=np.float64).astype(f32).reshape(-1, 3)
    box = [float(x) for x in box]
    (pairs, dist, error), method, msgs = call(a, b, cutoff, box)
    check_mdatk(a, b, cutoff, box, pairs, dist, error)
    case = {"name": name, "group": group, "a": points_json(a), "b": points_json(b), "box": box,
            "cutoff": float(cutoff), "method": method}
    if error is not None:
        case["error"] = error
    else:
        assert pairs.dtype == np.intp and pairs.shape == (len(dist), 2) and dist.dtype == np.float64
        case["pairs"] = [[int(i), int(j)] for i, j in pairs]
        case["distances"] = [float(x).hex() for x in dist]
    if msgs:
        case["warnings"] = msgs
    if band:
        bp = band_pairs(a, b, box, cutoff)
        if bp:
            case["band"] = bp
    CASES.append(case)
    return case


def add_lattice(name, group, la, lb, cutoff, box):
    a, b = lattice(la), lattice(lb)
    (pairs, dist, error), method, msgs = call(a, b, cutoff, box)
    check_mdatk(a, b, cutoff, box, pairs, dist, error)
    case = {"name": name, "group": group, "lattice_a": la, "lattice_b": lb, "box": [float(x) for x in box],
            "cutoff": float(cutoff), "method": method}
    if error is not None:
        case["error"] = error
    else:
        case["count"] = int(len(pairs))
        case["pairs_sha256"] = hashlib.sha256(pairs.astype("<i8").tobytes()).hexdigest()
        case["distances_sha256"] = hashlib.sha256(dist.astype("<f8").tobytes()).hexdigest()
    if msgs:
        case["warnings"] = msgs
    CASES.append(case)


ORTHO = [25.0, 25.0, 25.0, 90.0, 90.0, 90.0]
TRIC = [25.0, 25.0, 25.0, 60.0, 70.0, 80.0]


def pad(points, count=10):
    """Add points around (9, 9, 9) so that the set has at least COUNT points (grid route); they are
    more than 4 A from the test points near the origin and the faces of the boxes used with it."""
    points = [list(p) for p in points]
    k = 0
    while len(points) < count:
        points.append([9.0 + 0.25 * k, 9.0, 9.0])
        k += 1
    return points


def rand_points(n, lo, hi):
    return rng.uniform(lo, hi, (n, 3))


# ---------------------------------------------------------------------------------------------
# 1. the automatic rule: sizes 9 / 10, n * m >= 1e8 (lattices, below), cutoff vs 0.3 * size
def auto_rule():
    base = rand_points(12, 0.0, 25.0)
    for n, m in ((9, 12), (12, 9), (10, 10), (1, 12), (12, 1), (9, 9)):
        add(f"sizes {n}x{m}", "auto", base[:n], base[:m], 3.0, ORTHO)
    for cutoff in (7.4999995, 7.4999998, 7.5, 7.5000002, 7.5000005, 7.500001, 7.6, 12.5, 12.6):
        add(f"ortho 25 cutoff {cutoff!r}", "auto", base[:10], base[:12], cutoff, ORTHO)
    box = [10.0, 20.0, 30.0, 90.0, 90.0, 90.0]
    for cutoff in (2.9999999, 3.0, 3.0000002, 3.0001, 5.0):
        add(f"ortho 10/20/30 cutoff {cutoff!r}", "auto", base[:10] * [0.4, 0.8, 1.2], base[:12] * [0.4, 0.8, 1.2],
            cutoff, box)
    for box in ([25.0, 25.0, 25.0, 60.0, 70.0, 80.0], [40.0, 40.0, 40.0, 90.0, 90.0, 120.0],
                [30.0, 30.0, 30.0, 109.4712, 109.4712, 109.4712], [20.0, 20.0, 20.0, 30.0, 30.0, 30.0],
                [10.0, 10.0, 10.0, 90.0, 90.0, 30.0]):
        for cutoff in (1.0, 1.25, 1.3, 1.5, 2.5, 4.0, 6.0, 7.5, 9.0):
            add(f"triclinic {box[3:]} cutoff {cutoff}", "auto", base[:10] * box[0] / 25, base[:11] * box[0] / 25,
                cutoff, box)
    # angles within float32 rounding of 90 count as orthogonal; 90.00001 does not
    for ang in (90.000001, 90.00001, 89.99999):
        add(f"angle {ang!r}", "auto", base[:10], base[:12], 7.5, [25.0, 25.0, 25.0, 90.0, 90.0, ang])


# ---------------------------------------------------------------------------------------------
# 2. boundary: points at (or one float32 step around) the cutoff, both methods
def boundary():
    for L in (25.0, 30.0, 31.7, 18.3, 77.1, 40.0, 60.0):
        box = [L, L, L, 90.0, 90.0, 90.0]
        for c in (3.5, 1.25, 4.0, 2.6, 3.3):
            cf = float(f32(c))
            x = [cf, float(np.nextafter(f32(c), f32(0))), float(np.nextafter(f32(c), f32(100)))]
            pts = [[1.0 + v, 1.0, 1.0] for v in x] + [[1.0, 1.0 + v, 1.0] for v in x] + \
                  [[1.0, 1.0, 1.0 + v] for v in x] + [[1.0 - v, 1.0, 1.0] for v in x]
            origin = [[1.0, 1.0, 1.0]]
            # brute force: one reference point
            add(f"axis bf L={L} c={c}", "boundary", origin, pts, c, box)
            # grid: ten reference points (nine far away)
            add(f"axis grid L={L} c={c}", "boundary", pad(origin), pad(pts), c, box)
            # across the periodic face
            # across the periodic face: from L - 0.5 to v - 0.5 (and beyond the far face)
            across = [[v - 0.5, 1.0, 1.0] for v in x] + [[L + v - 0.5, 1.0, 1.0] for v in x]
            add(f"face bf L={L} c={c}", "boundary", [[L - 0.5, 1.0, 1.0]], across, c, box)
            add(f"face grid L={L} c={c}", "boundary", pad([[L - 0.5, 1.0, 1.0]]), pad(across), c, box)
    # cutoffs equal to distances that are sums of squares: sqrt(d2) <= c and d2 <= c*c may differ
    found = 0
    tries = 0
    while found < 24 and tries < 200000:
        tries += 1
        p = f32(rng.uniform(0.5, 4.0)); q = f32(rng.uniform(0.5, 4.0)); r = f32(rng.uniform(0.0, 2.0))
        d2 = (float(p) * float(p) + float(q) * float(q)) + float(r) * float(r)
        c = math.sqrt(d2)
        for cand in (c, math.nextafter(c, 0.0), math.nextafter(c, math.inf)):
            if (math.sqrt(d2) <= cand) != (d2 <= cand * cand):
                pts = [[2.0 + float(p), 2.0 + float(q), 2.0 + float(r)]]
                add(f"square tie bf {found}", "boundary", [[2.0, 2.0, 2.0]], pts, cand, ORTHO)
                add(f"square tie grid {found}", "boundary", pad([[2.0, 2.0, 2.0]]), pad(pts), cand, ORTHO)
                tri = [40.0, 40.0, 40.0, 70.0, 80.0, 75.0]
                add(f"square tie tric bf {found}", "boundary", [[2.0, 2.0, 2.0]], pts, cand, tri)
                add(f"square tie tric grid {found}", "boundary", pad([[2.0, 2.0, 2.0]]), pad(pts), cand, tri)
                found += 1
                break


# ---------------------------------------------------------------------------------------------
# 3. random inputs: inside the cell, up to a box outside, far outside; orthogonal and triclinic
def random_cases():
    boxes = [ORTHO, [40.0, 33.0, 51.0, 90.0, 90.0, 90.0], [31.7, 18.3, 77.1, 90.0, 90.0, 90.0], TRIC,
             [40.0, 40.0, 40.0, 109.4712, 109.4712, 109.4712], [60.0, 60.0, 60.0, 90.0, 90.0, 120.0],
             [30.0, 35.0, 40.0, 80.0, 100.0, 95.0], [20.0, 20.0, 20.0, 40.0, 50.0, 60.0]]
    for bi, box in enumerate(boxes):
        lmin = min(box[:3])
        for spread in (1, 3, 40):
            for k in range(2):
                n = int(rng.integers(10, 30)) if k == 0 else int(rng.integers(2, 10))
                m = int(rng.integers(10, 30))
                lo = -lmin * (spread - 1) / 2
                hi = lmin * (spread + 1) / 2
                a = rand_points(n, lo, hi)
                b = rand_points(m, lo, hi)
                c = float(rng.uniform(0.15, 0.3) * lmin) if spread != 40 else float(rng.uniform(0.2, 0.3) * lmin)
                add(f"random box{bi} spread{spread} #{k} {n}x{m}", "random", a, b, round(c, 3), box)


# ---------------------------------------------------------------------------------------------
# 4. outside the cell, cutoffs beyond half a box, self arrays, duplicates
def images():
    box = ORTHO
    a = [[0.0, 0.0, 0.0]]
    b = [[24.0, 0, 0], [12.5, 0, 0], [-12.5, 0, 0], [37.5, 0, 0], [25.0, 0, 0], [50.0, 0, 0], [-25.0, 0, 0],
         [1e4, 0, 0], [1e4 + 3, 0, 0], [-1e4 - 3, 0, 0], [1e6 + 1, 0, 0], [3.0, 4.0, 0.0], [28.0, -21.0, 50.0]]
    for c in (5.0, 12.5, 13.0, 20.0, 30.0, 1e9, math.inf):
        add(f"ortho images cutoff {c}", "images", a, b, c, box)
        add(f"ortho images grid cutoff {c}", "images", pad(a), pad(b), c, box)
    tri = TRIC
    for c in (5.0, 12.0, 20.0, 40.0):
        add(f"triclinic images cutoff {c}", "images", a, b, c, tri)
    # cutoff larger than half the box: one pair per (i, j), the nearest image only
    pts = rand_points(6, 0.0, 25.0)
    for box in (ORTHO, TRIC, [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]):
        add(f"big cutoff self {box}", "images", pts, pts, 30.0, box)
    # the same array on both sides: (i, i) pairs at 0 and both orders
    pts = rand_points(14, 0.0, 25.0)
    for box in (ORTHO, TRIC):
        add(f"self grid {box[3:]}", "self", pts, pts, 6.0, box)
        add(f"self bf {box[3:]}", "self", pts[:5], pts[:5], 9.0, box)
        add(f"self bf all {box[3:]}", "self", pts[:5], pts[:5], 100.0, box)
    dup = np.array([[3.0, 4.0, 5.0]] * 4 + [[3.0, 4.0, 5.0 + 25.0], [3.0 - 25.0, 4.0, 5.0], [10.0, 10.0, 10.0]] * 2)
    for box in (ORTHO, TRIC):
        add(f"duplicates grid {box[3:]}", "self", dup, dup, 1.0, box)
        add(f"duplicates bf {box[3:]}", "self", dup[:4], dup, 1.0, box)


# ---------------------------------------------------------------------------------------------
# 5. the order of the grid method: dense random sets, several cell counts
def order_cases():
    for box, c in ((ORTHO, 1.25), (ORTHO, 3.0), (ORTHO, 6.2), ([60.0, 45.0, 52.0, 90.0, 90.0, 90.0], 4.5),
                   ([60.0, 45.0, 52.0, 90.0, 90.0, 90.0], 0.5), ([60.0, 60.0, 60.0, 90.0, 90.0, 120.0], 5.0),
                   ([40.0, 40.0, 40.0, 109.4712, 109.4712, 109.4712], 3.5), (TRIC, 2.0)):
        lmin = min(box[:3])
        a = rand_points(40, -lmin * 0.5, lmin * 1.5)
        b = rand_points(60, -lmin * 0.5, lmin * 1.5)
        # concentrate points so that many pairs exist
        b[:30] = a[:30] + rng.normal(0.0, c * 0.6, (30, 3))
        add(f"order {box} c={c}", "order", a, b, c, box)
    pts = rand_points(30, 0.0, 25.0)
    add("order self 25 c=7", "order", pts, pts, 7.0, ORTHO)  # every axis a single cell


# ---------------------------------------------------------------------------------------------
# 6. degenerate inputs reachable after MDATK's lengths > 0 filter
def degenerate():
    pts = np.array([[x, y, z] for x in range(3) for y in range(2) for z in range(2)], dtype=np.float64) * 4 + 1
    few = pts[:4]
    for c in (0.0, -1.0, -1e-9, math.nan, math.inf):
        add(f"cutoff {c} grid", "degenerate", pts, pts, c, ORTHO)
        add(f"cutoff {c} bf", "degenerate", few, pts, c, ORTHO)
        add(f"cutoff {c} triclinic grid", "degenerate", pts, pts, c, TRIC)
    # cutoff 0 with points on opposite faces
    faces = [[0.0, 1.0, 1.0], [25.0, 1.0, 1.0], [-25.0, 1.0, 1.0], [12.5, 1.0, 1.0]]
    add("cutoff 0 faces bf", "degenerate", faces, faces, 0.0, ORTHO)
    add("cutoff 0 faces grid", "degenerate", pad(faces), pad(faces), 0.0, ORTHO)
    for ang in ([0.0, 0.0, 0.0], [90.0, 90.0, 180.0], [180.0, 180.0, 180.0], [math.nan, 90.0, 90.0],
                [90.0, 90.0, math.nan], [120.0, 120.0, 120.0], [60.0, 60.0, 150.0], [90.0, 90.0, 0.5],
                [30.0, 30.0, 60.0], [-90.0, 90.0, 90.0], [270.0, 90.0, 90.0], [1e-3, 1e-3, 1e-3]):
        box = [25.0, 25.0, 25.0] + ang
        add(f"angles {ang} grid-size", "degenerate", pts, pts, 3.0, box)
        add(f"angles {ang} bf", "degenerate", few, pts, 3.0, box)
        add(f"angles {ang} cutoff 0", "degenerate", pts, pts, 0.0, box)
    for lengths in ([1e-10, 25.0, 25.0], [1e-50, 25.0, 25.0], [1e-6, 1e-6, 1e-6], [1e30, 1e30, 1e30],
                    [math.inf, 25.0, 25.0], [3.0, 25.0, 25.0], [25.0, 25.0, 0.5]):
        add(f"lengths {lengths} ortho", "degenerate", pts, pts, 3.0, lengths + [90.0, 90.0, 90.0])
        add(f"lengths {lengths} ortho bf", "degenerate", few, pts, 3.0, lengths + [90.0, 90.0, 90.0])
        add(f"lengths {lengths} triclinic", "degenerate", pts, pts, 3.0, lengths + [60.0, 70.0, 80.0])
    # non-finite coordinates: orthogonal box only (triclinic is inconsistent in MDAnalysis, see spec)
    for bad in (math.nan, math.inf, -math.inf):
        for where in ("a", "b"):
            a = pts.copy()
            b = pts.copy()
            (a if where == "a" else b)[3] = [bad, 1.0, 1.0]
            add(f"non-finite {bad} in {where} grid", "degenerate", a, b, 3.0, ORTHO)
            add(f"non-finite {bad} in {where} bf", "degenerate", a[:4], b, 3.0, ORTHO)
    for big in (1e10, -1e10, 3e38):
        a = pts.copy()
        a[3] = [big, 1.0, 1.0]
        if abs(big) < 1e20:  # beyond about 2^31 box lengths the grid's shift count overflows
            add(f"huge {big} grid", "degenerate", a, pts, 3.0, ORTHO)
        add(f"huge {big} bf", "degenerate", a[:4], pts, 3.0, ORTHO)
        if abs(big) < 1e20:
            add(f"huge {big} triclinic", "degenerate", a, pts, 3.0, TRIC)
            add(f"huge {big} triclinic bf", "degenerate", a[:4], pts, 3.0, TRIC)


# ---------------------------------------------------------------------------------------------
# 7. n * m >= 1e8: the grid whatever the cutoff
def large():
    la = {"origin": [0.5, 0.5, 0.5], "spacing": 1.0, "counts": [25, 25, 16]}
    lb = {"origin": [0.75, 0.6, 0.8], "spacing": 1.0, "counts": [25, 25, 16]}
    add_lattice("lattice 10000x10000 c=2", "large", la, lb, 2.0, ORTHO)
    add_lattice("lattice self 10000x10000 c=2", "large", la, la, 2.0, ORTHO)
    add_lattice("lattice 10000x10000 c=13 (too large)", "large", la, lb, 13.0, ORTHO)
    add_lattice("lattice 10000x10000 c=inf", "large", la, lb, math.inf, ORTHO)
    # a flat box: the cutoff exceeds 0.3 * c, but n * m >= 1e8 still sends it to the grid
    flat = [1000.0, 1000.0, 1.0, 90.0, 90.0, 90.0]
    lf = {"origin": [0.5, 0.5, 0.5], "spacing": 1.0, "counts": [100, 100, 1]}
    lg = {"origin": [0.75, 0.6, 0.5], "spacing": 1.0, "counts": [100, 100, 1]}
    add_lattice("lattice flat 10000x10000 c=0.5", "large", lf, lg, 0.5, flat)
    add_lattice("lattice flat 10000x10000 c=0.6 (too large)", "large", lf, lg, 0.6, flat)
    add_lattice("lattice 10000x10000 triclinic c=3", "large", la, lb, 3.0, TRIC)
    add_lattice("lattice 10000x10000 triclinic c=10 (too large)", "large", la, lb, 10.0, TRIC)


def main():
    auto_rule()
    boundary()
    random_cases()
    images()
    order_cases()
    degenerate()
    large()
    names = [c["name"] for c in CASES]
    dups = sorted({n for n in names if names.count(n) > 1})
    assert not dups, dups
    doc = {
        "meta": {
            "mdanalysis": MDAnalysis.__version__,
            "numpy": np.__version__,
            "call": "capped_distance(a.astype(float32), b.astype(float32), max_cutoff=cutoff, box=box, "
                    "return_distances=True)",
            "band": "|d - cutoff| <= 8 * ulp32(max(|coordinates of the pair|, a, b, c))",
        },
        "cases": CASES,
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "distance_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=None, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
