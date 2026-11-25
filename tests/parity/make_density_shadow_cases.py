"""What MDATK's density shadow (mdatk/density.py: project_density, _above_level, shadow_rgba,
anchored) and the scipy.ndimage.affine_transform call behind it compute, for
src/render/molecule/density_shadow.cpp (tests/render/test_density_shadow.cpp; milestone 4, chunk
W0-I of docs/specs/m4-plan.md). The density fit is not covered here (D9).

    /opt/mdbridge-venv/bin/python -I -B make_density_shadow_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies (numpy, scipy); MDAnalysis is not used (no topology is read). MDATK runs
in-process; nothing is written into the checkout. The maps are committed files of tests/data
(io/maps/files, written by make_map_cases.py, and render/inputs/demo_dry_map.mrc, written by
make_figure_cases.py), read with MDATK's read_map; each case names its map relative to tests/data.
Writes into OUT_DIR:

  affine.json    scipy.ndimage.affine_transform(order=1, mode="constant", prefilter=False) of
                 float32 arrays: input (shape, values), 3 x 3 matrix, offset, output shape, cval,
                 the float32 output (coordinates exactly on the last grid line, outside, NaN)
  shadows.json   _above_level(map, level) of the small maps; project_density(map, placement,
                 view, centre, level, window, max_pixels, max_depth) -> image (float32, rows
                 bottom -> top) and extent, or None, or the error; shadow_rgba of some of them
  anchored.json  anchored(values, anchor, reference).to_list(12) or the error

MDBridge's improvements (docs/DIFFERENCES.md) are applied to the expected data, test-side, in the
style of the improvements_*.py modules (improvements_render.py RN-30 ... RN-32):

  RN-30  _above_level: non-finite voxels (NaN, +-inf) carry no density (MDATK: a NaN voxel makes
         the whole shadow NaN, an infinite one makes it NaN or blank); a non-finite level gives no
         shadow (None)
  RN-31  project_density with a non-finite placement, view, centre or map grid raises
         ValueError("density map placement, view or grid is not finite") (MDATK: "cannot convert
         float NaN to integer" / OverflowError "cannot convert float infinity to integer")
  RN-32  anchored with an anchor that is not n x 3 numbers (ragged rows) or with non-finite anchor /
         reference points uses the stored placement as it is and says why (MDATK: numpy's
         "setting an array element with a sequence..." / LinAlgError "SVD did not converge", and
         the map is not drawn)

Each record keeps MDATK's own outcome under "mdatk" where an improvement changes it. Floats are
written exactly (repr of the float64 value; float32 arrays widened). numpy's and OpenBLAS' kernels
are pinned before numpy is imported, as make_analysis_cases.py pins them (the view / placement
products and kabsch's SVD go through BLAS / LAPACK).
"""

import os
import platform
import sys

if True:  # before numpy is imported (read when OpenBLAS / numpy load)
    os.environ["OPENBLAS_CORETYPE"] = "Nehalem"
    for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_var] = "1"
    if platform.machine().lower() in ("x86_64", "amd64"):
        os.environ["NPY_ENABLE_CPU_FEATURES"] = "X86_V2"

import json  # noqa: E402
import math  # noqa: E402
import warnings  # noqa: E402

if len(sys.argv) != 3:
    sys.exit("usage: make_density_shadow_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
sys.path.insert(0, CHECKOUT)
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
sys.path.insert(0, HERE)

warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402

from mdatk import density as D  # noqa: E402
from mdatk.io.mrc import read_map  # noqa: E402
from mdatk.render.orientation import euler_to_matrix  # noqa: E402

NOT_FINITE = "density map placement, view or grid is not finite"


# --- output -----------------------------------------------------------------------------------------

def num(x):
    return float(x)


def flat(a):
    return [num(v) for v in np.asarray(a, dtype=np.float64).ravel()]


def write_cases(name, header, cases):
    """One case per line (compact), so that diffs stay readable."""
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as f:
        f.write('{"note": ' + json.dumps(header) + ',\n"cases": [\n')
        for k, case in enumerate(cases):
            f.write(json.dumps(case, ensure_ascii=False, separators=(",", ":"), allow_nan=True))
            f.write(",\n" if k + 1 < len(cases) else "\n")
        f.write("]}\n")


def error_of(exc):
    return {"type": type(exc).__name__, "message": str(exc)}


# --- MDATK with the improvements (test-side) ----------------------------------------------------------

_mdatk_above_level = D._above_level
_mdatk_project = D.project_density
_mdatk_anchored = D.anchored


def improved_above_level(dmap, level):
    """RN-30: max(map - level, 0) with non-finite voxels (and any non-finite difference) as 0."""
    diff = dmap.data - np.float32(level)
    return np.where(np.isfinite(diff), np.maximum(diff, 0), 0).astype(np.float32)


def improved_project(dmap, world_to_map, view, centre, level, **kw):
    """RN-30 / RN-31 around MDATK's project_density (its _above_level replaced)."""
    if not math.isfinite(level):
        return None
    D._above_level = improved_above_level
    D._LEVEL_CACHE.clear()
    try:
        return _mdatk_project(dmap, world_to_map, view, centre, level, **kw)
    except (ValueError, OverflowError) as exc:
        if "cannot convert float" in str(exc):
            raise ValueError(NOT_FINITE) from exc
        raise
    finally:
        D._above_level = _mdatk_above_level
        D._LEVEL_CACHE.clear()


def anchor_problem(anchor, reference):
    """RN-32: why the anchor cannot be used (None when MDATK's rule applies)."""
    rows = list(anchor)
    if any(not hasattr(r, "__len__") or len(r) != len(rows[0]) for r in rows) if rows else False:
        return "density anchor ignored: its points are not all x, y, z"
    a = np.asarray(rows, dtype=np.float64)
    ref = np.asarray(reference, dtype=np.float64)
    if a.ndim != 2 or a.shape != ref.shape or len(a) < 3:
        return None  # MDATK: the placement as stored
    if not np.all(np.isfinite(a)):
        return "density anchor ignored: it holds non-finite numbers"
    if not np.all(np.isfinite(ref)):
        return "density anchor ignored: the reference frame holds non-finite coordinates"
    return None


# --- affine_transform ---------------------------------------------------------------------------------

def affine_cases():
    rng = np.random.RandomState(20261125)
    cases = []

    def case(name, data, matrix, offset, out_shape, cval=0.0):
        data = np.asarray(data, dtype=np.float32)
        matrix = np.asarray(matrix, dtype=np.float64)
        offset = np.asarray(offset, dtype=np.float64)
        out = ndimage.affine_transform(data, matrix, offset=offset, output_shape=tuple(out_shape), order=1,
                                       mode="constant", cval=cval, prefilter=False)
        assert out.dtype == np.float32
        cases.append({"name": name, "shape": list(data.shape), "data": flat(data), "matrix": flat(matrix),
                      "offset": flat(offset), "output_shape": list(out_shape), "cval": cval,
                      "output": flat(out)})

    a = rng.rand(4, 5, 6).astype(np.float32) * 3 - 0.5
    case("identity", a, np.eye(3), [0, 0, 0], a.shape)  # last grid lines exactly (edge taps)
    case("shift_half", a, np.eye(3), [-0.5, 0.25, 0.5], a.shape)
    case("shift_out", a, np.eye(3), [3.5, -4.5, 6.0], a.shape)
    case("cval", a, np.eye(3), [0.5, 0.5, 0.5], (5, 6, 7), cval=-2.5)
    case("scale_down", a, np.diag([0.5, 0.4, 0.3]), [0.1, 0.2, 0.3], (8, 9, 10))
    r = D.rotvec_matrix(np.array([0.3, -0.4, 0.5]))
    case("rotated", a, r * 1.3, [1.0, -0.5, 2.0], (6, 7, 5))
    case("rotated_neg", a, -r, [3.0, 4.0, 5.0], (5, 5, 5))
    b = rng.rand(1, 3, 2).astype(np.float32)
    case("thin", b, np.eye(3), [0, 0, 0], (1, 3, 2))
    case("thin_moved", b, [[0, 0, 0], [0, 0.5, 0], [0, 0, 0.5]], [0, 0.25, 0.25], (2, 5, 3))
    c = np.arange(8, dtype=np.float32).reshape(2, 2, 2)
    case("nan_offset", c, np.eye(3), [float("nan"), 0, 0], (2, 2, 2))
    case("tiny", c, np.full((3, 3), 0.25), [0.1, 0.2, 0.3], (3, 3, 3))
    case("empty_output", c, np.eye(3), [0, 0, 0], (0, 2, 2))
    return cases


# --- shadows --------------------------------------------------------------------------------------------

MAPS = {
    "demo": "render/inputs/demo_dry_map.mrc",
    "axes_321": "io/maps/files/axes_321le.mrc",
    "nstart_231": "io/maps/files/axes_nstart_231le.mrc",
    "permuted": "io/maps/files/mxyz_permuted.mrc",
    "large": "io/maps/files/large.map.gz",
    "constant": "io/maps/files/mode2_constant.mrc",
    "nan": "io/maps/files/mode2_nan.mrc",
    "inf": "io/maps/files/mode2_inf.mrc",
    "specials": "io/maps/files/mode2_specials_le.mrc",
    "cell_nan": "io/maps/files/cell_nan.mrc",
}


def load(key):
    return read_map(os.path.join(DATA, *MAPS[key].split("/")))


def transform(rotation, translation):
    return D.RigidTransform(np.asarray(rotation, dtype=np.float64), np.asarray(translation, dtype=np.float64))


def above_level_cases():
    cases = []
    for key, levels in (("axes_321", [0.0, 1.25, 2.4]), ("nan", [0.0, 1.5]), ("inf", [0.5]),
                        ("specials", [0.0, -1e38]), ("constant", [0.1, 0.0]), ("nan", [float("nan")]),
                        ("axes_321", [float("-inf")])):
        m = load(key)
        for level in levels:
            rec = {"kind": "above_level", "map": MAPS[key], "level": level, "values": flat(improved_above_level(m, level))}
            mdatk = _mdatk_above_level(m, level)
            D._LEVEL_CACHE.clear()
            if not np.array_equal(mdatk, improved_above_level(m, level), equal_nan=False):
                rec["mdatk"] = flat(mdatk)
            cases.append(rec)
    return cases


def shadow_case(name, key, placement, view, centre, level, window=None, max_pixels=220, max_depth=220,
                rgba=()):
    m = load(key)
    kw = {"window": window, "max_pixels": max_pixels, "max_depth": max_depth}
    rec = {"kind": "project", "name": name, "map": MAPS[key], "rotation": flat(placement.rotation),
           "translation": flat(placement.translation), "view": flat(view), "centre": flat(centre),
           "level": num(level), "window": None if window is None else [num(x) for x in window],
           "max_pixels": max_pixels, "max_depth": max_depth}
    try:
        sh = improved_project(m, placement, np.asarray(view, dtype=np.float64), np.asarray(centre, dtype=np.float64),
                              level, **kw)
    except Exception as exc:  # noqa: BLE001 - recorded
        rec["error"] = error_of(exc)
        sh = None
    else:
        if sh is None:
            rec["shadow"] = None
        else:
            assert sh.image.dtype == np.float32
            rec["shadow"] = {"rows": sh.image.shape[0], "cols": sh.image.shape[1], "image": flat(sh.image),
                             "extent": [num(x) for x in sh.extent]}
            rec["rgba"] = []
            for color, opacity, gamma in rgba:
                out = D.shadow_rgba(sh, color, opacity, gamma)
                assert out.dtype == np.float32 and out.shape == sh.image.shape + (4,)
                rec["rgba"].append({"color": list(color), "opacity": opacity, "gamma": gamma, "values": flat(out)})
    # MDATK's own outcome where the improvements change it
    try:
        D._LEVEL_CACHE.clear()
        theirs = _mdatk_project(m, placement, np.asarray(view, dtype=np.float64),
                                np.asarray(centre, dtype=np.float64), level, **kw)
        D._LEVEL_CACHE.clear()
        if "error" in rec:
            rec["mdatk"] = {"shadow": "image" if theirs is not None else None}
        elif (theirs is None) != (sh is None) or (
                theirs is not None and not np.array_equal(theirs.image, sh.image)):
            rec["mdatk"] = {"shadow": None if theirs is None else
                            ("all NaN" if np.isnan(theirs.image).all() else
                             "with NaN" if np.isnan(theirs.image).any() else "other image")}
    except Exception as exc:  # noqa: BLE001 - recorded
        D._LEVEL_CACHE.clear()
        if rec.get("error") != error_of(exc):
            rec["mdatk"] = {"error": error_of(exc)}
    return rec


def shadow_cases():
    cases = []
    ident = D.RigidTransform.identity()
    eye = np.eye(3)
    demo = load("demo")
    demo_centre = (demo.origin + demo.extent) / 2
    lvl1 = D.map_level(demo, 1.0)
    lvl0 = D.map_level(demo, 0.0)
    placed = transform(D.quaternion_matrix(np.array([0.9, 0.1, 0.3, -0.2])), [5.0, 8.0, -3.0])
    rgba = [((0.1, 0.2, 0.9), 0.5, 0.7), ((0.8901960784313725, 0.28627450980392155, 0.2823529411764706), 0.45, 0.7),
            ((1.0, 1.0, 1.0), 1.3, 1.0), ((0.0, 0.5, 0.25), -0.2, 2.2), ((0.3, 0.6, 0.9), 0.8, 0.0)]
    cases.append(shadow_case("demo_identity", "demo", ident, eye, demo_centre, lvl1, rgba=rgba))
    cases.append(shadow_case("demo_level0", "demo", ident, eye, demo_centre, lvl0, rgba=rgba[:1]))
    for angles in ([0, 0, 0], [40, -30, 10], [90, 0, 45], [-120, 75, 200]):
        v = euler_to_matrix(angles)
        cases.append(shadow_case(f"demo_view_{angles[0]}_{angles[1]}_{angles[2]}", "demo", ident, v, demo_centre,
                                 lvl1, rgba=rgba[:1]))
        world_centre = placed.inverse().apply(demo_centre[None, :])[0] + np.array([1.5, -2.0, 0.5])
        cases.append(shadow_case(f"demo_placed_{angles[0]}_{angles[1]}_{angles[2]}", "demo", placed, v,
                                 world_centre, lvl1))
    v = euler_to_matrix([40, -30, 10])
    # windows: clipped, fast quality (100), outside, NaN entries (Python's max / min ignore them)
    cases.append(shadow_case("demo_window", "demo", ident, v, demo_centre, lvl1, window=(-10.0, 12.0, -8.0, 5.0)))
    cases.append(shadow_case("demo_window_fast", "demo", ident, v, demo_centre, lvl1,
                             window=(-30.0, 30.0, -25.0, 25.0), max_pixels=100, max_depth=100))
    cases.append(shadow_case("demo_window_outside", "demo", ident, v, demo_centre, lvl1,
                             window=(100.0, 120.0, -8.0, 5.0)))
    cases.append(shadow_case("demo_window_nan", "demo", ident, v, demo_centre, lvl1,
                             window=(float("nan"), 12.0, -8.0, float("nan"))))
    cases.append(shadow_case("demo_coarse", "demo", placed, v, demo_centre, lvl1, max_pixels=12, max_depth=7))
    cases.append(shadow_case("demo_level_high", "demo", ident, v, demo_centre, 1e6))  # nothing above: None
    cases.append(shadow_case("demo_far", "demo", transform(eye, [500.0, 0, 0]), eye, demo_centre, lvl1,
                             window=(-30.0, 30.0, -30.0, 30.0)))
    # small maps: permuted axes, non-cubic voxels, NSTART origins
    for key in ("axes_321", "nstart_231", "permuted", "large"):
        m = load(key)
        c = (m.origin + m.extent) / 2
        for angles in ([0, 0, 0], [25, 50, -35]):
            cases.append(shadow_case(f"{key}_{angles[1]}", key, ident, euler_to_matrix(angles), c,
                                     D.map_level(m, 0.5), max_pixels=60 if key == "large" else 220,
                                     max_depth=60 if key == "large" else 220))
    m = load("constant")
    cases.append(shadow_case("constant_mean", "constant", ident, eye, (m.origin + m.extent) / 2,
                             D.map_level(m, 0.0)))  # nothing above the mean: None
    cases.append(shadow_case("constant_below", "constant", ident, eye, (m.origin + m.extent) / 2, 0.05))
    # RN-30: non-finite voxels / levels
    for key in ("nan", "inf"):
        m = load(key)
        c = (m.origin + m.extent) / 2
        cases.append(shadow_case(f"{key}_voxel", key, ident, euler_to_matrix([10, 20, 30]), c, 1.0, rgba=rgba[:1]))
        cases.append(shadow_case(f"{key}_map_level", key, ident, eye, c, D.map_level(m, 1.0)))
    cases.append(shadow_case("specials", "specials", ident, eye, (load("specials").origin + load("specials").extent) / 2,
                             0.0))
    # RN-31: non-finite inputs
    cases.append(shadow_case("nan_centre", "demo", ident, eye, demo_centre * float("nan"), lvl1))
    nan_view = eye.copy()
    nan_view[1, 2] = float("nan")
    cases.append(shadow_case("nan_view", "demo", ident, nan_view, demo_centre, lvl1))
    cases.append(shadow_case("inf_translation", "demo", transform(eye, [float("inf"), 0, 0]), eye, demo_centre, lvl1))
    cases.append(shadow_case("huge_translation", "demo", transform(euler_to_matrix([30, 0, 0]), [1e308, -1e308, 0]),
                             eye, demo_centre, lvl1))
    cases.append(shadow_case("nan_voxel_size", "cell_nan", ident, eye, np.zeros(3), 0.0))
    return cases


# --- anchored --------------------------------------------------------------------------------------------

def anchored_cases():
    rng = np.random.RandomState(4242)
    cases = []
    ref = rng.rand(12, 3) * 30.0 - 15.0
    anchor_exact = D.anchor_points(ref)
    values = D.RigidTransform(D.quaternion_matrix(np.array([0.4, 0.3, -0.8, 0.2])),
                              np.array([20.0, -15.0, 8.0])).to_list()
    moved = transform(D.rotvec_matrix(np.array([0.2, -0.1, 0.35])), [3.0, -1.0, 2.5])
    ref_moved = moved.apply(ref) + rng.normal(scale=0.05, size=ref.shape)
    noisy = [list(map(float, r)) for r in (np.asarray(anchor_exact) + 0.0004)]
    nan_anchor = [list(r) for r in anchor_exact]
    nan_anchor[3][1] = float("nan")
    ref_nan = ref.copy()
    ref_nan[5, 0] = float("nan")
    ragged = [list(r) for r in anchor_exact]
    ragged[2] = ragged[2][:2]
    rounded = [round(x, 4) for x in values]  # rounded rotation: re-orthonormalised by from_list
    for name, vals, anchor, reference in (
            ("same", values, anchor_exact, ref),
            ("within_tolerance", values, noisy, ref),
            ("moved", values, anchor_exact, ref_moved),
            ("moved_rounded_values", rounded, anchor_exact, ref_moved),
            ("identity_moved", D.RigidTransform.identity().to_list(), anchor_exact, ref_moved),
            ("other_count", values, anchor_exact[:-1], ref_moved),
            ("two_points", values, anchor_exact[:2], ref_moved[:2]),
            ("three_points", values, anchor_exact[:3], ref_moved[:3]),
            ("empty", values, [], ref_moved),
            ("xy_only", values, [r[:2] for r in anchor_exact], ref_moved),
            ("ragged", values, ragged, ref_moved),
            ("nan_anchor", values, nan_anchor, ref_moved),
            ("nan_reference", values, anchor_exact, ref_nan),
            ("eleven_numbers", values[:11], anchor_exact, ref_moved),
            ("reflection", [-1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0, 0, 0, 0], anchor_exact, ref_moved),
            ("not_rotation", [2.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0, 0, 0, 0], anchor_exact, ref_moved),
            ("nan_values", [float("nan")] + values[1:], anchor_exact, ref_moved)):
        rec = {"name": name, "values": [num(x) for x in vals], "anchor": [[num(x) for x in r] for r in anchor],
               "reference": [flat(r) for r in np.asarray(reference)]}
        problem = None
        try:
            problem = anchor_problem(anchor, reference)
            if problem is not None:
                rec["transform"] = D.RigidTransform.from_list(vals).to_list(12)
                rec["note"] = problem
            else:
                rec["transform"] = _mdatk_anchored(vals, anchor, reference).to_list(12)
        except Exception as exc:  # noqa: BLE001 - recorded
            rec["error"] = error_of(exc)
        if problem is not None:
            try:
                _mdatk_anchored(vals, anchor, reference)
                rec["mdatk"] = {"error": None}
            except Exception as exc:  # noqa: BLE001 - recorded
                rec["mdatk"] = {"error": error_of(exc)}
        cases.append(rec)
    return cases


def main():
    ver = {"numpy": np.__version__, "scipy": __import__("scipy").__version__}
    write_cases("affine.json", f"scipy.ndimage.affine_transform, order 1, mode constant ({ver})", affine_cases())
    write_cases("shadows.json", f"mdatk.density _above_level / project_density / shadow_rgba, RN-30, RN-31 ({ver})",
                above_level_cases() + shadow_cases())
    write_cases("anchored.json", f"mdatk.density.anchored, RN-32 ({ver})", anchored_cases())


main()
