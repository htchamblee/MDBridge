"""What MDATK's rigid transforms and rotation helpers give (mdatk/density.py), for
tests/core/test_rigid_transform.cpp.

    /opt/mdbridge-venv/bin/python -I -B make_rigid_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy and scipy (mdatk.density imports them); nothing is written into the checkout. Writes
OUT_DIR/rigid_cases.json:

  from_list         RigidTransform.from_list of 12-number lists (rotations rounded to 1-12 digits,
                    perturbed, scaled, reflected, singular, huge / tiny / non-finite entries, wrong
                    lengths): the re-orthonormalised rotation and the translation, or the message
  to_list           to_list(digits) of transforms given exactly
  algebra           apply / inverse / compose / about of transforms and points
  rotvec, quaternion, angle
                    rotvec_matrix, quaternion_matrix, rotation_angle
  anchor_points     anchor_points of point lists

Floats are JSON numbers (repr, NaN / Infinity / -Infinity as json writes them). numpy computes
products with BLAS and the SVD with LAPACK, so the C++ side compares within a few 1e-16.
"""

import json
import math
import os
import sys
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_rigid_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import numpy as np  # noqa: E402

from mdatk import density as D  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
rng = np.random.default_rng(20261124)


def fhex(x):
    return float(x)  # JSON numbers: repr() reads back exactly; NaN / Infinity as json writes them


def hexes(values):
    return [fhex(v) for v in np.asarray(values, dtype=np.float64).ravel()]


def random_rotation():
    return D.quaternion_matrix(rng.normal(size=4))


def near_boundary(values):
    """Singular values within 1e-9 of MDATK's 0.01 tolerance: platform-dependent, skipped."""
    v = np.asarray(values, dtype=np.float64)
    if v.size != 12 or not np.all(np.isfinite(v)):
        return False
    sv = np.linalg.svd(v[:9].reshape(3, 3), compute_uv=False)
    return bool(np.any(np.abs(np.abs(sv - 1.0) - 0.01) < 1e-9))


def from_list_case(values):
    case = {"input": hexes(values) if len(values) else []}
    try:
        t = D.RigidTransform.from_list(values)
        case["rotation"] = hexes(t.rotation)
        case["translation"] = hexes(t.translation)
    except ValueError as exc:
        case["error"] = str(exc)
    return case


def from_list_cases():
    inputs = []
    inputs.append([1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
    inputs.append([0, -1, 0, 1, 0, 0, 0, 0, 1, 1.5, -2.0, 1e6])
    inputs.append([0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0])  # cyclic permutation
    inputs.append([0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0])  # swap: reflection
    inputs.append([-1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0])  # inversion
    inputs.append([-1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0])  # 180 degrees about z
    for digits in (1, 2, 3, 4, 5, 6, 6, 6, 6, 8, 10, 12):
        r = random_rotation()
        t = rng.normal(size=3) * 10 ** rng.uniform(-2, 4)
        inputs.append([round(float(x), digits) for x in np.concatenate([r.ravel(), t])])
    for _ in range(40):  # stored as MDATK stores them
        r = random_rotation()
        inputs.append([round(float(x), 6) for x in np.concatenate([r.ravel(), rng.normal(size=3) * 30])])
    for noise in (1e-4, 1e-3, 3e-3, 5e-3, 2e-2, 0.1):
        for _ in range(4):
            r = random_rotation() + rng.normal(size=(3, 3)) * noise
            inputs.append(list(r.ravel()) + [0.0, 1.0, 2.0])
    for scale in (0.95, 0.992, 0.995, 1.004, 1.008, 1.012, 1.05, 2.0):
        r = random_rotation() * scale
        inputs.append(list(r.ravel()) + [3.0, 2.0, 1.0])
    for _ in range(6):  # reflections
        r = random_rotation()
        r[rng.integers(3)] *= -1
        inputs.append([round(float(x), 6) for x in r.ravel()] + [0.0, 0.0, 0.0])
    for _ in range(4):  # a stretched axis
        r = random_rotation() @ np.diag([1.0, 1.0, 1.0 + rng.uniform(0.003, 0.009)])
        inputs.append(list(r.ravel()) + [1.0, 1.0, 1.0])
    inputs.append([0.0] * 12)
    inputs.append([1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0])  # singular
    inputs.append([1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
    inputs.append([1e300, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
    inputs.append([1e-300, 0, 0, 0, 1e-300, 0, 0, 0, 1e-300, 0, 0, 0])
    inputs.append([5e-324, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
    inputs.append([1.0101, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
    inputs.append([1.0099, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
    inputs.append([-1.0099, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
    inputs.append([1, 0, 0, 0, 1, 0, 0, 0, 1, 1e308, -1e308, 5e-324])
    for k in range(12):
        for bad in (float("nan"), float("inf"), float("-inf")):
            v = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
            v[k] = bad
            inputs.append(v)
    for n in (0, 1, 3, 9, 11, 13, 16):
        inputs.append([float(i % 2) for i in range(n)])
    inputs.append([float("nan")] * 11)
    cases = []
    for values in inputs:
        values = [float(x) for x in values]
        if near_boundary(values):
            continue
        cases.append(from_list_case(values))
    return cases


def to_list_cases():
    cases = []
    specials = [0.5, -0.5, 2.5e-7, -2.5e-7, 1e-7, -1e-7, 0.1234565, -0.0, 0.0, 1.0000005, 123456.7891235,
                -9.9999995e-5, 1e16, 5e-324]
    for k in range(12):
        r = random_rotation()
        t = rng.normal(size=3) * 10 ** rng.uniform(-3, 5, size=3)
        if k < 4:
            flat = np.array([specials[(3 * k + i) % len(specials)] for i in range(12)])
            r = flat[:9].reshape(3, 3)
            t = flat[9:12]
        for digits in (0, 1, 3, 6, 9, 12, 15, 20):
            tr = D.RigidTransform(np.array(r, dtype=np.float64), np.array(t, dtype=np.float64))
            cases.append({"rotation": hexes(tr.rotation), "translation": hexes(tr.translation),
                          "digits": digits, "list": hexes(tr.to_list(digits))})
    return cases


def transform():
    return D.RigidTransform(random_rotation(), rng.normal(size=3) * 10 ** rng.uniform(-1, 3))


def algebra_cases():
    cases = []
    for _ in range(40):
        a, b = transform(), transform()
        points = rng.normal(size=(int(rng.integers(1, 6)), 3)) * 10 ** rng.uniform(-1, 3)
        centre = rng.normal(size=3) * 20
        inv = a.inverse()
        comp = a.compose(b)
        about = a.about(centre)
        cases.append({
            "a": {"rotation": hexes(a.rotation), "translation": hexes(a.translation)},
            "b": {"rotation": hexes(b.rotation), "translation": hexes(b.translation)},
            "points": [hexes(p) for p in points],
            "centre": hexes(centre),
            "apply": [hexes(p) for p in a.apply(points)],
            "inverse": {"rotation": hexes(inv.rotation), "translation": hexes(inv.translation)},
            "compose": {"rotation": hexes(comp.rotation), "translation": hexes(comp.translation)},
            "about": {"rotation": hexes(about.rotation), "translation": hexes(about.translation)},
        })
    return cases


def rotation_cases():
    rotvec = []
    for w in ([0, 0, 0], [1e-13, 0, 0], [0, 0, 1e-12], [0, 2e-12, 0], [math.pi, 0, 0], [0, 0, -math.pi],
              [0.3, -0.2, 0.5], [0.15, -0.1, 0.1]):
        rotvec.append({"w": hexes(w), "matrix": hexes(D.rotvec_matrix(np.array(w, dtype=float)))})
    for _ in range(30):
        w = rng.normal(size=3) * rng.uniform(0, 4)
        rotvec.append({"w": hexes(w), "matrix": hexes(D.rotvec_matrix(w))})
    quaternion = []
    for q in ([1, 0, 0, 0], [0, 1, 0, 0], [0.3, -0.8, 0.4, 0.2], [0.9, 0.1, 0.3, -0.2], [2, 0, 0, 0],
              [0.5, 0.5, -0.5, 0.5]):
        quaternion.append({"q": hexes(q), "matrix": hexes(D.quaternion_matrix(np.array(q, dtype=float)))})
    for _ in range(30):
        q = rng.normal(size=4) * 10 ** rng.uniform(-3, 3)
        quaternion.append({"q": hexes(q), "matrix": hexes(D.quaternion_matrix(q))})
    angle = []
    mats = [np.eye(3), -np.eye(3), np.diag([1.0, -1.0, -1.0]), np.eye(3) * 1.0000001,
            np.diag([-1.0, -1.0, -1.0000001]), np.full((3, 3), np.nan)]
    mats += [random_rotation() for _ in range(30)]
    mats += [D.rotvec_matrix(np.array([0, 0, a])) for a in (1e-8, 1e-4, 0.5, 3.1415)]
    for m in mats:
        angle.append({"matrix": hexes(m), "degrees": fhex(D.rotation_angle(m))})
    return rotvec, quaternion, angle


def anchor_cases():
    cases = []
    for pts in ([[0.0005, -0.0005, 1.2345]], [[2.675, -2.675, 1e-9], [123456.78951, -0.0, 0.0]]):
        cases.append({"points": [hexes(p) for p in pts], "anchor": [hexes(p) for p in D.anchor_points(np.array(pts))]})
    pts = rng.normal(size=(20, 3)) * 50
    cases.append({"points": [hexes(p) for p in pts], "anchor": [hexes(p) for p in D.anchor_points(pts)]})
    return cases


def dump(data):
    """One case per line."""
    parts = []
    for key, cases in data.items():
        rows = ",\n".join("  " + json.dumps(c, separators=(",", ":")) for c in cases)
        parts.append(f" {json.dumps(key)}: [\n{rows}\n ]")
    return "{\n" + ",\n".join(parts) + "\n}\n"


def main():
    os.makedirs(OUT, exist_ok=True)
    rotvec, quaternion, angle = rotation_cases()
    data = {"from_list": from_list_cases(), "to_list": to_list_cases(), "algebra": algebra_cases(),
            "rotvec": rotvec, "quaternion": quaternion, "angle": angle, "anchor_points": anchor_cases()}
    with open(os.path.join(OUT, "rigid_cases.json"), "w", encoding="utf-8") as fh:
        fh.write(dump(data))
    print({k: len(v) for k, v in data.items()})


main()
