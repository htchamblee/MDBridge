"""What MDATK's rigid-body density fit (mdatk/density.py) computes, for MDBridge's chem::density_fit
(tests/chem/test_density_fit.cpp; docs/specs/m6-plan.md 1.5, 4.1, 4.2, owner decisions D2, D10).

    python3 -I -B make_density_fit_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy and scipy; MDAnalysis for the demo frame (mdatk.demo's synthetic system read through the
pipeline, with MDBridge's independent tables patched in, independent_tables.py). MDBridge's chem
improvements are applied (improvements_chem.py: CH-11 non-finite voxels as no density). numpy's and
OpenBLAS' kernels are pinned (chem_common.pin_kernels); the fits still depend on the CPU in the last
bits (Powell's path follows every CC value), so regenerate_all.py compares this file within the chem
tolerance (FITTED) and the C++ test within D10's (0.5°, 0.2 Å, CC 1e-3). Writes:

  OUT_DIR/density_fit_cases.json
    rotations   rotation_samples(step): count, the first, second and last matrices (and all of them
                for 30° and 45°), or the error
    fft_size    _fft_size(n) for n = 1 ... 400
    problems    per model x map: the FitProblem (centre, radius, rg, envelope point count, the first
                points and values, sum of values), cc at given transforms, inside_fraction,
                at_resolution(8) cc, flips, the coarse search (spacing, origin, grid, var_floor, the
                peaks of three rotations), refine from a start, fit_model global and manual (transform,
                cc, start, candidates, evaluations, message), or the errors
    errors      FitError / ValueError cases (step below 3°, no resolution, fewer than three atoms, a
                map smaller than the coarse grid)
  OUT_DIR/maps/fit_*.mrc   the maps (target maps as tests/test_density.py makes them: simulated at the
                true pose, padded, Gaussian noise of a fixed seed; one with unrelated density, one cut to
                a partial overlap, one with a NaN voxel)
  The models are given as coordinates (float.hex) and elements.
"""

import json
import os
import sys
import tempfile

if len(sys.argv) != 3:
    sys.exit("usage: make_density_fit_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
sys.path.insert(0, CHECKOUT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_common  # noqa: E402

chem_common.pin_kernels()
import numpy as np  # noqa: E402

import improvements_chem  # noqa: E402

improvements_chem.apply_mdatk()
from mdatk import builder as B  # noqa: E402
from mdatk import density as D  # noqa: E402
from mdatk import pdbtools as pt  # noqa: E402
from mdatk.io.mrc import DensityMap, write_map  # noqa: E402

error, fhex, matrix, vector = chem_common.error, chem_common.fhex, chem_common.matrix, chem_common.vector


def transform_dict(t):
    return {"rotation": matrix(t.rotation), "translation": vector(t.translation)}


def apply(t, xyz):
    """t.apply in plain Python (no BLAS: the maps' bytes must not depend on the CPU)."""
    r, s = t.rotation.tolist(), t.translation.tolist()
    return np.array([[sum(r[i][k] * p[k] for k in range(3)) + s[i] for i in range(3)] for p in xyz.tolist()])


def target_map(xyz, elements, true, res=6.0, pad=6, noise=0.1, seed=3, spacing=2.0):
    """tests/test_density.py's _target_map (a smaller pad and grid)."""
    sim = D.simulate_map(apply(true, xyz), D.atom_weights(elements), res, spacing=spacing)
    big = np.zeros(tuple(np.array(sim.data.shape) + 2 * pad), np.float32)
    big[pad:-pad, pad:-pad, pad:-pad] = sim.data
    big += np.random.default_rng(seed).normal(scale=noise * sim.data.max(), size=big.shape).astype(np.float32)
    return DensityMap(big, sim.origin - pad * spacing, sim.voxel)


def save(name, dmap):
    os.makedirs(os.path.join(OUT, "maps"), exist_ok=True)
    write_map(os.path.join(OUT, "maps", name), dmap)
    return f"maps/{name}"


def models():
    out = {}
    text = B.build_helix("GGCAUGCC")
    xyz, el = B.model_atoms(pt.read_pdb(text))
    out["helix"] = (xyz, el)
    single = B.build_helix("GAAGAGAAGGC", partner="none")
    out["strand"] = B.model_atoms(pt.read_pdb(single))
    try:
        import independent_tables

        independent_tables.apply()
        from mdatk.config import Project
        from mdatk.demo import build_demo
        from mdatk.pipeline import Pipeline

        with tempfile.TemporaryDirectory(prefix="mdatk-dens-") as d:
            demo = build_demo(d, n_frames=3, solvent=False)
            p = Project.load(demo["preset"])
            p.output.out_dir = os.path.join(d, "out")
            xyz, elements, _f, _ref = Pipeline(p).density_frame_atoms("ref")
        out["demo"] = (np.array(xyz, dtype=np.float64), list(elements))
    except ImportError as exc:  # pragma: no cover - the reference environment has MDAnalysis
        sys.exit(f"the demo frame needs MDAnalysis: {exc}")
    return out


def rotation_cases():
    out = []
    for step in (60.0, 45.0, 30.0, 20.0, 15.0, 10.0, 3.0, 2.99, float("nan"), -5.0):
        case = {"step": fhex(step)}
        try:
            r = D.rotation_samples(step)
            case.update(count=len(r), first=matrix(r[0]), second=matrix(r[1]), last=matrix(r[-1]))
            if step in (30.0, 45.0):
                case["all"] = [matrix(m) for m in r]
        except D.FitError as exc:
            case["error"] = error(exc)
        out.append(case)
    return out


def problem_case(name, xyz, el, map_path, dmap, res, true, starts, global_step=30.0):
    case = {"model": name, "map": map_path, "resolution": fhex(res), "true": transform_dict(true)}
    try:
        prob = D.FitProblem(improvements_chem.finite_map(dmap)[0], xyz, D.atom_weights(el), res)
    except D.FitError as exc:
        case["error"] = error(exc)
        return case
    case["problem"] = {"centre": vector(prob.centre), "radius": fhex(prob.radius), "rg": fhex(prob.rg),
                       "points": len(prob.points), "first_points": [vector(p) for p in prob.points[:5]],
                       "first_values": vector(prob.s[:5]), "sum": fhex(float(prob.s.sum()))}
    transforms = [D.RigidTransform.identity(), true] + starts
    case["cc"] = [[transform_dict(t), fhex(prob.cc(t)), fhex(prob.inside_fraction(t))] for t in transforms]
    coarse8 = prob.at_resolution(D.COARSE_RESOLUTION) if res < D.COARSE_RESOLUTION else prob
    case["cc_at_8"] = [fhex(coarse8.cc(t)) for t in transforms]
    case["flips"] = [transform_dict(f) for f in prob.flips()]
    cs = D._CoarseSearch(prob, global_step)
    rots = D.rotation_samples(global_step)
    case["coarse"] = {"h": fhex(cs.h), "lo": vector(cs.lo), "n": [int(v) for v in cs.n], "var_floor": fhex(cs.var_floor),
                      "peaks": [[k, [[fhex(s), vector(c)] for s, c in cs.peaks(rots[k], 2)]] for k in (0, 1, len(rots) // 2)]}
    refined = []
    for t in starts:
        before = prob.evaluations
        fit, cc = prob.refine(t)
        refined.append({"start": transform_dict(t), "transform": transform_dict(fit), "cc": fhex(cc),
                        "evaluations": prob.evaluations - before})
    case["refine"] = refined
    fits = []
    for search, start in (("global", None),) + tuple(("manual", t) for t in starts[:1]) + (("manual", None),):
        try:
            r = D.fit_model(dmap, xyz, el, res, start=start, search=search, angular_step=global_step)
            fits.append({"search": search, "start": None if start is None else transform_dict(start),
                         "transform": transform_dict(r.transform), "cc": fhex(r.cc), "start_kind": r.start,
                         "candidates": r.candidates, "evaluations": r.evaluations, "message": r.message})
        except D.FitError as exc:
            fits.append({"search": search, "error": error(exc)})
    case["fit_model"] = fits
    return case


def near(true, centre, rotvec, shift):
    c = true.apply(np.asarray(centre)[None, :])[0]
    return D.RigidTransform(D.rotvec_matrix(np.array(rotvec)), np.array(shift, dtype=np.float64)).about(c).compose(true)


def problem_cases(ms):
    out = []
    helix_xyz, helix_el = ms["helix"]
    true = D.RigidTransform(D.quaternion_matrix([0.3, -0.8, 0.4, 0.2]), np.array([25.0, -40.0, 12.0]))
    c = helix_xyz.mean(axis=0)
    emap = target_map(helix_xyz, helix_el, true)
    starts = [near(true, c, [0.15, -0.1, 0.1], [2.0, -1.5, 1.0]), near(true, c, [-0.3, 0.2, 0.0], [0.0, 3.0, -2.0])]
    out.append(problem_case("helix", helix_xyz, helix_el, save("fit_helix.mrc", emap), emap, 6.0, true, starts))
    # unrelated density: a denser blob beside the molecule (Roseman's local normalisation ignores it)
    other = emap.data.copy()
    nz, ny, nx = other.shape
    zz, yy, xx = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    blob = np.exp(-((zz - 3) ** 2 + (yy - 3) ** 2 + (xx - 3) ** 2) / (2 * 2.5 ** 2)).astype(np.float32)
    other += 3.0 * float(emap.data.max()) * blob
    omap = DensityMap(other, emap.origin, emap.voxel)
    out.append(problem_case("helix", helix_xyz, helix_el, save("fit_unrelated.mrc", omap), omap, 6.0, true, starts[:1]))
    # partial overlap: the box cut through the molecule
    cut = DensityMap(np.ascontiguousarray(emap.data[:, :, : nx // 2 + 2]), emap.origin, emap.voxel)
    out.append(problem_case("helix", helix_xyz, helix_el, save("fit_partial.mrc", cut), cut, 6.0, true, starts[:1]))
    # a NaN voxel inside the molecule's density (CH-11)
    nan = emap.data.copy()
    nan[nz // 2, ny // 2, nx // 2] = np.nan
    nmap = DensityMap(nan, emap.origin, emap.voxel)
    out.append(problem_case("helix", helix_xyz, helix_el, save("fit_nan.mrc", nmap), nmap, 6.0, true, starts[:1]))
    # a single strand at 4 Å, 45° step
    s_xyz, s_el = ms["strand"]
    true2 = D.RigidTransform(D.rotvec_matrix(np.array([0.4, -0.9, 0.3])), np.array([30.0, -12.0, 55.0]))
    smap = target_map(s_xyz, s_el, true2, res=4.0, spacing=1.5, noise=0.05, seed=5)
    out.append(problem_case("strand", s_xyz, s_el, save("fit_strand.mrc", smap), smap, 4.0, true2,
                            [near(true2, s_xyz.mean(axis=0), [0.05, 0.05, -0.1], [1.0, 1.0, 0.5])], global_step=45.0))
    # the demo frame (tests/test_density.py's model)
    d_xyz, d_el = ms["demo"]
    dmap = target_map(d_xyz, d_el, true, res=6.0, spacing=2.0)
    out.append(problem_case("demo", d_xyz, d_el, save("fit_demo.mrc", dmap), dmap, 6.0, true,
                            [near(true, d_xyz.mean(axis=0), [0.15, -0.1, 0.1], [2.0, -1.5, 1.0])]))
    return out


def error_cases(ms):
    xyz, el = ms["helix"]
    small = DensityMap(np.ones((3, 3, 3), np.float32), np.zeros(3), np.full(3, 1.0))
    out = []
    runs = [("step below 3", lambda: D.rotation_samples(2.0)),
            ("no resolution", lambda: D.FitProblem(small, xyz, D.atom_weights(el), 0.0)),
            ("negative resolution", lambda: D.FitProblem(small, xyz, D.atom_weights(el), -1.0)),
            ("two atoms", lambda: D.FitProblem(small, xyz[:2], D.atom_weights(el[:2]), 6.0)),
            ("no atoms to simulate", lambda: D.simulate_map(np.zeros((0, 3)), np.zeros(0), 6.0)),
            ("tiny map global", lambda: D.fit_model(small, xyz, el, 6.0, search="global", angular_step=45.0)),
            ("tiny map manual", lambda: D.fit_model(small, xyz, el, 6.0, search="manual"))]
    for label, run in runs:
        try:
            r = run()
            got = {"value": {"cc": fhex(r.cc), "message": r.message, "start": r.start}} if isinstance(r, D.FitResult) \
                else {"value": "ok"}
        except (D.FitError, ValueError) as exc:
            got = {"error": error(exc)}
        out.append([label, got])
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    ms = models()
    data = {"models": {k: {"xyz": [vector(p) for p in xyz], "elements": list(el)} for k, (xyz, el) in ms.items()},
            "rotations": rotation_cases(), "fft_size": [[n, D._fft_size(n)] for n in range(1, 401)],
            "problems": problem_cases(ms), "errors": error_cases(ms)}
    with open(os.path.join(OUT, "density_fit_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
