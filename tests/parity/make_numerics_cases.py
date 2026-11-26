"""What SciPy 1.18.1 computes for the numerics MDBridge ports or vendors (src/numerics), for
tests/numerics/test_*.cpp.

    /opt/mdbridge-venv/bin/python -I -B make_numerics_cases.py OUT_DIR

Needs numpy 2.5.5 and scipy 1.18.1 (the reference environment); MDATK is not needed. Writes
OUT_DIR/numerics_cases.json:

  fft            scipy.fft.rfftn(a, workers=-1) and irfftn(c, s=shape, axes=(0, 1, 2), workers=-1) of
                 3-D float64 arrays a[i] = value(i) (formula below; the C++ test builds them itself):
                 every output value for small shapes, the SHA-256 of the little-endian bytes for all;
                 irfftn of conj(rfftn(a)) * rfftn(b) (b[i] = value(i + 7)), MDATK's correlation
                 product (mdatk/density.py _CoarseSearch.peaks), and of a half spectrum that is not
                 the transform of a real array (imaginary parts the inverse ignores)
  next_fast_len  scipy.fft.next_fast_len(n, real=True)
  gaussian       scipy.ndimage.gaussian_filter of float32 arrays d[i] = value32(i) with one sigma per
                 axis, modes "reflect" (scipy's default, MDATK's density fit) and "constant",
                 kernels longer than the axis included
  map_coords     scipy.ndimage.map_coordinates(d, coords, order=1, mode="constant", cval, prefilter=False)
                 (mdatk/density.py sample_map): points inside, on the last grid line, just outside,
                 far outside and NaN, cval 0 / 1.5 / NaN; float32 results
  powell         scipy.optimize.minimize(f, x0, method="Powell", callback, options) of the named
                 objectives below (the C++ test defines them the same way): every iterate the
                 callback sees (x, fun), the final x, fun, direc, nit, nfev, status, message; MDATK's
                 options (xtol 0.02, ftol 1e-7, maxfev 4000, direc = step * I) among them
  brent          scipy.optimize.minimize_scalar(f, method="brent", options={"xtol": ...}): x, fun,
                 nit, nfev, success (a monotone objective: no valid bracket)
  least_squares  scipy.optimize.least_squares(fun, x0, ...) of the named residual functions below:
                 dense (tr_solver "exact", m >= n and m < n) and sparse (jac_sparsity as MDATK's
                 _span_sparsity builds it, tr_solver "lsmr", x_scale, tr_options atol / btol 1e-5 as
                 MDATK's _CLOSE_LSMR), x_scale "jac", diff_step, max_nfev; every iterate the callback
                 sees (x, cost, nit, nfev) and the result (x, cost, fun, grad, optimality, nfev, njev,
                 status, message); ulp_spread: how far SciPy's own iterates move (relative, at most 1)
                 when every residual moves by about an ulp (the largest of three such runs;
                 ulp_same_count: each with as many iterations),
                 the finite differences' sensitivity to rounding that the C++ test allows for
  groups         scipy.optimize._numdiff.group_columns of sparse structures, and
                 numpy.random.RandomState(0).permutation(n) (its column order)
  lsmr           scipy.sparse.linalg.lsmr(A, b, damp, atol, btol, maxiter) of small dense matrices:
                 x, istop, itn, normr, normar, normA, condA, normx; ulp_spread as for least_squares
                 (b moved by about an ulp)
  roots          the real roots of np.roots(c) (np.isreal), sorted, of quartics and lower degrees
                 (leading and trailing zeros included)
  pairs          scipy.spatial.cKDTree(points).query_pairs(r) as a sorted list: random points on a
                 1/64 grid, a unit lattice at exactly r = 1 and sqrt(2), duplicates, one point, none

value(i) = ((i * 2654435761) mod 2^32) / 2^32 - 0.5 (exact in float64); value32(i) = float32(value(i)).

Floats are JSON numbers (repr: they read back exactly; NaN as json writes it). The Powell and
least-squares objectives use Python's math module (C libm) and arithmetic in a fixed order, never
numpy reductions, so the C++ test can evaluate them identically. least_squares uses BLAS (dot
products, norms, J.T @ f) and LAPACK (SVD, QR, np.roots' eigenvalues), so its iterates depend on the
CPU's kernels in the last bits: the generator pins OpenBLAS to its generic kernels and one thread
and numpy to its baseline (as make_analysis_cases.py does), and regenerate_all.py compares this file
within a relative 1e-12 where its bytes differ (MACHINE_DEPENDENT). The C++ tests compare the FFT
bit for bit, the filters within float32 rounding, Powell and least-squares iterates within 1e-9,
pairs as sets.
"""

import os
import platform
import sys

if len(sys.argv) != 2:
    sys.exit("usage: make_numerics_cases.py OUT_DIR")
os.environ["OPENBLAS_CORETYPE"] = "Nehalem"  # before numpy is imported
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_var] = "1"
if platform.machine().lower() in ("x86_64", "amd64"):
    os.environ["NPY_ENABLE_CPU_FEATURES"] = "X86_V2"

import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy import fft as sfft  # noqa: E402
from scipy import ndimage, optimize  # noqa: E402
from scipy.optimize._numdiff import group_columns  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.sparse.linalg import lsmr  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

warnings.filterwarnings("ignore")
OUT = os.path.abspath(sys.argv[1])
if scipy.__version__ != "1.18.1" or np.__version__ != "2.5.5":
    sys.exit(f"error: needs scipy 1.18.1 and numpy 2.5.5 (found {scipy.__version__}, {np.__version__})")


def value(i):
    return ((i * 2654435761) % 2**32) / 2**32 - 0.5


def values(n, offset=0):
    return np.array([value(i + offset) for i in range(n)], dtype=np.float64)


def floats(a):
    return [float(v) for v in np.asarray(a, dtype=np.float64).ravel()]


def sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).astype("<f8", copy=False).tobytes()).hexdigest()


def csha(c):
    return sha(np.ascontiguousarray(c).view(np.float64))


# ---- FFT --------------------------------------------------------------------------------------------

def fft_cases():
    out = []
    small = [(1, 1, 1), (1, 1, 2), (2, 3, 4), (3, 1, 5), (1, 7, 1), (5, 6, 7), (8, 9, 10), (4, 4, 3)]
    big = [(12, 10, 9), (16, 15, 8), (27, 25, 24), (11, 13, 17), (64, 48, 40), (96, 90, 80), (30, 32, 45)]
    for shape in small + big:
        n = int(np.prod(shape))
        a = values(n).reshape(shape)
        b = values(n, 7).reshape(shape)
        c = sfft.rfftn(a, workers=-1)
        r = sfft.irfftn(c, s=shape, axes=(0, 1, 2), workers=-1)
        prod = sfft.irfftn(np.conj(sfft.rfftn(b, workers=-1)) * c, s=shape, axes=(0, 1, 2), workers=-1)
        # a half spectrum with arbitrary imaginary parts in the bins a real signal has real
        hs = (values(c.size, 3) + 1j * values(c.size, 11)).reshape(c.shape)
        raw = sfft.irfftn(hs, s=shape, axes=(0, 1, 2), workers=-1)
        case = {"shape": list(shape), "rfftn_sha256": csha(c), "irfftn_sha256": sha(r),
                "product_sha256": sha(prod), "raw_sha256": sha(raw)}
        if shape in small:
            case["rfftn"] = floats(c.view(np.float64))
            case["irfftn"] = floats(r)
            case["raw"] = floats(raw)
        out.append(case)
    return out


# ---- ndimage ----------------------------------------------------------------------------------------

def gaussian_cases():
    out = []
    for shape, sigma, mode in [((5, 6, 7), (0.7, 1.3, 0.0), "reflect"), ((5, 6, 7), (2.5, 0.4, 3.1), "reflect"),
                               ((3, 9, 2), (1.0, 1.0, 1.0), "reflect"), ((1, 4, 12), (0.0, 0.8, 5.0), "reflect"),
                               ((5, 6, 7), (0.7, 1.3, 0.0), "constant"), ((4, 4, 4), (1e-16, 2.0, 0.6), "reflect")]:
        d = values(int(np.prod(shape))).astype(np.float32).reshape(shape)
        g = ndimage.gaussian_filter(d, sigma=list(sigma), mode=mode)
        out.append({"shape": list(shape), "sigma": list(sigma), "mode": mode, "out": floats(g)})
    return out


def map_coords_cases():
    out = []
    shape = (4, 5, 6)
    d = values(int(np.prod(shape))).astype(np.float32).reshape(shape)
    d[1, 2, 3] = np.nan  # NaN data propagate to the points that touch them
    pts = [(0, 0, 0), (3, 4, 5), (1.5, 2.25, 3.75), (3.0, 0.5, 5.0), (2.999999, 4, 4.5), (-1e-9, 1, 1),
           (3.0000001, 1, 1), (1, -0.5, 2), (1, 2, 6.0), (1, 2, 5.5), (10, 2, 2), (-3, -3, -3),
           (float("nan"), 1, 1), (1, float("nan"), 1), (0.5, 1.5, 2.5), (1.2, 2.4, 3.1), (0, 4, 0),
           (2.5, 3.5, 4.999), (0.75, 0.25, 0.125)]
    rng = np.random.default_rng(20261126)
    pts += [tuple(float(v) for v in p) for p in rng.uniform(-0.6, 5.6, size=(30, 3))]
    coords = np.array(pts, dtype=np.float64).T
    for cval in (0.0, 1.5, float("nan")):
        r = ndimage.map_coordinates(d, coords, order=1, mode="constant", cval=cval, prefilter=False)
        out.append({"shape": list(shape), "nan_at": [1, 2, 3], "cval": cval, "points": [list(p) for p in pts],
                    "out": floats(r)})
    return out


# ---- Powell / Brent ---------------------------------------------------------------------------------

def rosen(x):
    s = 0.0
    for i in range(len(x) - 1):
        s += 100.0 * (x[i + 1] - x[i] * x[i]) ** 2 + (1.0 - x[i]) ** 2
    return s


QUAD_A = [[4.0, 1.0, 0.5], [1.0, 3.0, -0.25], [0.5, -0.25, 2.0]]
QUAD_C = [0.5, -1.25, 2.0]


def quad3(x):
    s = 0.0
    for i in range(3):
        for j in range(3):
            s += (x[i] - QUAD_C[i]) * QUAD_A[i][j] * (x[j] - QUAD_C[j])
    return s + 0.125


GAUSS_C = [0.4, -0.3, 0.2, 1.1, -0.7, 0.5]
GAUSS_W = [0.9, 1.4, 0.6, 0.3, 0.5, 0.8]


def gauss6(x):
    """A smooth score like FitProblem.refine's: -exp(-sum w (x - c)^2) plus a small coupling."""
    s = 0.0
    for i in range(6):
        s += GAUSS_W[i] * (x[i] - GAUSS_C[i]) ** 2
    s += 0.1 * x[0] * x[3]
    return -math.exp(-s)


def gauss6_shift(x):
    """gauss6 with MDATK's max_shift guard (1.0 outside a shift of 1.2)."""
    if math.sqrt(x[3] * x[3] + x[4] * x[4] + x[5] * x[5]) > 1.2:
        return 1.0
    return gauss6(x)


POWELL_OBJECTIVES = {"rosen": rosen, "quad3": quad3, "gauss6": gauss6, "gauss6_shift": gauss6_shift}


def powell_cases():
    out = []
    eye6 = [[1.0 if i == j else 0.0 for j in range(6)] for i in range(6)]
    runs = [
        ("rosen", [-1.2, 1.0], {}),
        ("rosen", [-1.2, 1.0], {"xtol": 1e-8, "ftol": 1e-10}),
        ("rosen", [1.3, 0.7, 0.8, 1.9], {}),
        ("rosen", [-1.2, 1.0], {"maxfev": 50}),
        ("rosen", [-1.2, 1.0], {"maxiter": 3}),
        ("quad3", [0.0, 0.0, 0.0], {}),
        ("quad3", [3.0, -2.0, 1.0], {"direc": [[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]]}),
        ("gauss6", [0.0] * 6, {"xtol": 0.02, "ftol": 1e-7, "maxfev": 4000,
                               "direc": [[0.75 * v for v in row] for row in eye6]}),
        ("gauss6", [0.1, 0.2, -0.1, 0.0, 0.3, 0.1], {"xtol": 0.02, "ftol": 1e-7, "maxfev": 4000,
                                                     "direc": [[2.0 * v for v in row] for row in eye6]}),
        ("gauss6_shift", [0.0] * 6, {"xtol": 0.02, "ftol": 1e-7, "maxfev": 4000,
                                     "direc": [[0.5 * v for v in row] for row in eye6]}),
        ("gauss6", [0.0] * 6, {"xtol": 0.02, "ftol": 1e-7, "maxfev": 40,
                               "direc": [[0.75 * v for v in row] for row in eye6]}),
    ]
    for name, x0, options in runs:
        f = POWELL_OBJECTIVES[name]
        iterates = []

        def cb(intermediate_result):
            iterates.append({"x": floats(intermediate_result.x), "fun": float(intermediate_result.fun)})

        opts = dict(options)
        if "direc" in opts:
            opts["direc"] = np.array(opts["direc"], dtype=np.float64)
        res = optimize.minimize(lambda p: f(p), np.array(x0, dtype=np.float64), method="Powell", callback=cb,
                                options=opts)
        out.append({"objective": name, "x0": x0, "options": options, "iterates": iterates, "x": floats(res.x),
                    "fun": float(res.fun), "direc": [floats(r) for r in res.direc], "nit": int(res.nit),
                    "nfev": int(res.nfev), "status": int(res.status), "message": res.message})
    return out


BRENT_OBJECTIVES = {
    "parabola": lambda t: (t - 2.0) ** 2 + 1.0,
    "quartic": lambda t: t ** 4 - 3.0 * t ** 3 + 2.0,
    "double_well": lambda t: (t - 0.3) ** 2 * (t + 1.0) ** 2 + 0.1 * t,
    "decreasing": lambda t: -t,
    "flat_far": lambda t: (t - 250.0) ** 2,
}


def brent_cases():
    out = []
    for name, f in BRENT_OBJECTIVES.items():
        for xtol in (1.48e-8, 1e-3):
            res = optimize.minimize_scalar(f, method="brent", options={"xtol": xtol})
            out.append({"objective": name, "xtol": xtol, "x": float(res.x), "fun": float(res.fun),
                        "nit": int(res.nit), "nfev": int(res.nfev), "success": bool(res.success)})
    return out


# ---- least_squares ----------------------------------------------------------------------------------

def r_rosen(x):
    return np.array([10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]])


EXP_T = [i / 4.0 for i in range(12)]
EXP_Y = [2.5 * math.exp(-1.3 * t) + 0.3 + 0.01 * (((i * 37) % 11) - 5) / 5.0 for i, t in enumerate(EXP_T)]


def r_expfit(x):
    return np.array([x[0] * math.exp(-x[1] * t) + x[2] - y for t, y in zip(EXP_T, EXP_Y)])


def r_under(x):
    return np.array([x[0] + 2.0 * x[1] - 3.0 * x[2] - 1.0, x[0] * x[1] - 0.5])


ROT_P = [[1.0, 0.0, 0.0], [0.0, 1.5, 0.0], [0.0, 0.0, 2.0], [1.0, 1.0, -0.5], [-0.5, 0.25, 1.0]]


def rodrigues(v, p):
    """Rotate p by the rotation vector v (Rodrigues), in a fixed order of operations."""
    th = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if th == 0.0:
        return list(p)
    k = [v[0] / th, v[1] / th, v[2] / th]
    c, s = math.cos(th), math.sin(th)
    kxp = [k[1] * p[2] - k[2] * p[1], k[2] * p[0] - k[0] * p[2], k[0] * p[1] - k[1] * p[0]]
    kp = k[0] * p[0] + k[1] * p[1] + k[2] * p[2]
    return [p[i] * c + kxp[i] * s + k[i] * kp * (1.0 - c) for i in range(3)]


ROT_TRUE = [0.3, -0.2, 0.5]
ROT_SHIFT = [0.4, -1.0, 0.25]
ROT_Q = [[a + b for a, b in zip(rodrigues(ROT_TRUE, p), ROT_SHIFT)] for p in ROT_P]


def r_rigid(x):
    """Six rigid-body parameters placing ROT_P on ROT_Q, plus soft priors (like _pair_fit's)."""
    res = []
    for p, q in zip(ROT_P, ROT_Q):
        r = rodrigues(x[:3], p)
        for i in range(3):
            res.append(r[i] + x[3 + i] - q[i])
    for i in range(6):
        res.append(0.05 * x[i])
    return np.array(res)


CHAIN_N = 10
CHAIN_C = [0.5 + 0.1 * i for i in range(CHAIN_N)]


def r_chain(x):
    """Links between neighbouring variables, a three-variable term, and priors: a banded Jacobian."""
    n = CHAIN_N
    res = [x[i] * x[i] + x[i + 1] - CHAIN_C[i] for i in range(n - 1)]
    res += [0.5 * (x[i] - x[i + 1] * x[i + 2]) for i in range(n - 2)]
    res += [0.3 * (x[i] - 1.0) for i in range(n)]
    return np.array(res)


def chain_sparsity():
    """_span_sparsity's construction: residual k depends on the variables lo[k] .. hi[k]."""
    n = CHAIN_N
    lo = list(range(n - 1)) + list(range(n - 2)) + list(range(n))
    hi = [i + 1 for i in range(n - 1)] + [i + 2 for i in range(n - 2)] + list(range(n))
    lo, hi = np.asarray(lo, dtype=int), np.asarray(hi, dtype=int)
    width = hi - lo + 1
    rows = np.repeat(np.arange(len(lo)), width)
    cols = np.repeat(lo, width) + np.arange(width.sum()) - np.repeat(np.cumsum(width) - width, width)
    return rows, cols, (len(lo), n)


LSQ_RESIDUALS = {"rosen": r_rosen, "expfit": r_expfit, "under": r_under, "rigid": r_rigid, "chain": r_chain}


def lsq_cases():
    rows, cols, shape = chain_sparsity()
    sparsity = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=shape)
    chain_x0 = [0.2 * ((i * 7) % 5) for i in range(CHAIN_N)]
    runs = [
        ("rosen", [2.0, 2.0], {}),
        ("rosen", [-1.2, 1.0], {"x_scale": [1.0, 3.0]}),
        ("expfit", [1.0, 1.0, 0.0], {}),
        ("expfit", [1.0, 1.0, 0.0], {"x_scale": [1.0, 0.5, 0.1]}),
        ("expfit", [1.0, 1.0, 0.0], {"x_scale": "jac"}),
        ("expfit", [1.0, 1.0, 0.0], {"diff_step": 1e-6}),
        ("expfit", [1.0, 1.0, 0.0], {"max_nfev": 5}),
        ("expfit", [1.0, 1.0, 0.0], {"tr_solver": "lsmr"}),
        ("under", [0.0, 0.0, 0.0], {}),
        ("rigid", [0.0] * 6, {}),
        ("rigid", [0.0] * 6, {"x_scale": [0.5, 0.5, 0.5, 1.0, 1.0, 1.0], "ftol": 1e-6}),
        ("chain", chain_x0, {}),
        ("chain", chain_x0, {"jac_sparsity": "chain"}),
        ("chain", chain_x0, {"jac_sparsity": "chain", "ftol": 1e-6,
                             "x_scale": [3.0 if i % 2 == 0 else 0.2 for i in range(CHAIN_N)],
                             "tr_options": {"atol": 1e-5, "btol": 1e-5}, "max_nfev": 40}),
        ("chain", chain_x0, {"jac_sparsity": "chain", "tr_options": {"regularize": False}}),
    ]
    out = []
    for name, x0, options in runs:
        fun = LSQ_RESIDUALS[name]
        iterates = []

        def cb(intermediate_result):
            iterates.append({"x": floats(intermediate_result.x), "cost": float(intermediate_result.cost),
                             "nit": int(intermediate_result.nit), "nfev": int(intermediate_result.nfev)})

        kw = dict(options)
        if kw.get("jac_sparsity") == "chain":
            kw["jac_sparsity"] = sparsity
        res = optimize.least_squares(fun, np.array(x0, dtype=np.float64), callback=cb, **kw)
        # SciPy's own sensitivity to rounding: the same run with every residual moved by about an ulp
        # (times 1 + 2^-52, 1 - 2^-53, 1 + 2^-51). The finite-difference Jacobian divides residual
        # differences by steps of sqrt(eps) max(1, |x|), so an ulp in a residual or an iterate moves
        # the next Jacobian by about 1e-8 relative, and the iterates after it by as much; no port that
        # rounds differently from LAPACK and BLAS can do better
        main = iterates
        spread, same_count = 0.0, True
        for factor in (1.0 + 2.0 ** -52, 1.0 - 2.0 ** -53, 1.0 + 2.0 ** -51):
            iterates = []
            optimize.least_squares(lambda x, k=factor: fun(x) * k, np.array(x0, dtype=np.float64), callback=cb,
                                   **kw)
            for a, b in zip(main, iterates):
                for u, v in zip(a["x"] + [a["cost"]], b["x"] + [b["cost"]]):
                    spread = max(spread, abs(u - v) / max(1.0, abs(v)))
            same_count = same_count and len(main) == len(iterates)
        iterates = main
        out.append({"residuals": name, "x0": x0, "options": options, "iterates": iterates, "x": floats(res.x),
                    "cost": float(res.cost), "fun": floats(res.fun), "grad": floats(res.grad),
                    "optimality": float(res.optimality), "nfev": int(res.nfev), "njev": int(res.njev),
                    "status": int(res.status), "message": res.message, "ulp_spread": spread,
                    "ulp_same_count": same_count})
    return out, {"rows": [int(v) for v in rows], "cols": [int(v) for v in cols], "shape": list(shape)}


def group_cases():
    perms = []
    for n in (1, 2, 3, 5, 10, 37, 100, 1000):
        perms.append({"n": n, "permutation": [int(v) for v in np.random.RandomState(0).permutation(n)]})
    groups = []
    rng = np.random.default_rng(20261126)
    structures = [chain_sparsity()]
    for m, n, density in [(6, 4, 0.4), (20, 15, 0.15), (40, 30, 0.1), (5, 5, 1.0), (8, 12, 0.05)]:
        mask = rng.uniform(size=(m, n)) < density
        r, c = np.nonzero(mask)
        structures.append((r, c, (m, n)))
    for r, c, shape in structures:
        a = csr_matrix((np.ones(len(r), dtype=np.int8), (r, c)), shape=shape)
        groups.append({"shape": list(shape), "rows": [int(v) for v in r], "cols": [int(v) for v in c],
                       "groups": [int(v) for v in group_columns(a)]})
    return perms, groups


def lsmr_cases():
    out = []
    rng = np.random.default_rng(7)
    for m, n, damp, atol, btol, maxiter in [(8, 5, 0.0, 1e-6, 1e-6, None), (8, 5, 0.5, 1e-6, 1e-6, None),
                                            (12, 12, 0.0, 1e-5, 1e-5, None), (5, 8, 0.0, 1e-6, 1e-6, None),
                                            (10, 6, 0.1, 1e-10, 1e-10, 3), (6, 6, 0.0, 1e-8, 1e-8, None)]:
        a = rng.integers(-8, 9, size=(m, n)) / 8.0
        b = rng.integers(-16, 17, size=m) / 16.0
        res = lsmr(a, b, damp=damp, atol=atol, btol=btol, maxiter=maxiter)
        # as for least_squares: how far x moves when b moves by about an ulp (Krylov iterations
        # without reorthogonalisation amplify rounding)
        spread = 0.0
        for factor in (1.0 + 2.0 ** -52, 1.0 - 2.0 ** -53, 1.0 + 2.0 ** -51):
            other = lsmr(a, b * factor, damp=damp, atol=atol, btol=btol, maxiter=maxiter)[0]
            spread = max(spread, float(np.max(np.abs(other - res[0]) / np.maximum(1.0, np.abs(res[0])))))
        out.append({"ulp_spread": spread, "a": [floats(row) for row in a], "b": floats(b), "damp": damp, "atol": atol, "btol": btol,
                    "maxiter": maxiter, "x": floats(res[0]), "istop": int(res[1]), "itn": int(res[2]),
                    "normr": float(res[3]), "normar": float(res[4]), "normA": float(res[5]),
                    "condA": float(res[6]), "normx": float(res[7])})
    return out


def roots_cases():
    out = []
    rng = np.random.default_rng(11)
    polys = [[1.0, -10.0, 35.0, -50.0, 24.0], [1.0, 0.0, -2.0, 0.0, 1.0], [0.0, 0.0, 1.0, -3.0, 2.0],
             [2.0, -3.0, 0.0, 0.0], [1.0, 0.0, 1.0, 0.0, 1.0], [0.0, 3.0, -1.0], [5.0], [0.0, 0.0, 2.0, 0.0],
             [1.0, -1.0, -1.0, 1.0, 0.0], [1e-3, 1.0, -2.0, 0.5, 0.1]]
    polys += [list(rng.integers(-9, 10, size=5).astype(float)) for _ in range(20)]
    polys += [list(rng.normal(size=5)) for _ in range(20)]
    for c in polys:
        t = np.roots(np.array(c, dtype=np.float64))
        real = sorted(float(v) for v in np.real(t[np.isreal(t)]))
        # a multiple root: LAPACK splits it by about sqrt(eps) (two real or a complex pair), which
        # the C++ test does not try to reproduce
        multiple = any(abs(t[i] - t[j]) < 1e-6 for i in range(len(t)) for j in range(i + 1, len(t)))
        out.append({"c": [float(v) for v in c], "real": real, "multiple": bool(multiple)})
    return out


def pair_cases():
    out = []
    rng = np.random.default_rng(3)
    sets = []
    pts = np.round(rng.uniform(0, 10, size=(200, 3)) * 64) / 64
    for r in (0.5, 1.0, 2.5, 3.1):
        sets.append((pts, r))
    lattice = np.array([[i, j, k] for i in range(4) for j in range(4) for k in range(3)], dtype=np.float64)
    sets += [(lattice, 1.0), (lattice, math.sqrt(2.0)), (lattice, 0.999999)]
    dup = np.array([[1.0, 2.0, 3.0]] * 4 + [[1.0, 2.0, 3.5]], dtype=np.float64)
    sets += [(dup, 0.0), (dup, 0.5), (np.array([[0.0, 0.0, 0.0]]), 1.0), (np.zeros((0, 3)), 1.0)]
    spread = np.array([[0.0, 0.0, 0.0], [1e6, 0.0, 0.0], [1e6 + 0.5, 0.0, 0.0], [0.0, -2e5, 3.0]])
    sets += [(spread, 1.0)]
    for p, r in sets:
        found = sorted((int(a), int(b)) for a, b in cKDTree(p).query_pairs(r)) if len(p) else []
        out.append({"points": floats(p), "r": r, "pairs": [list(x) for x in found]})
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    lsq, chain = lsq_cases()
    perms, groups = group_cases()
    data = {
        "versions": {"scipy": scipy.__version__, "numpy": np.__version__},
        "fft": fft_cases(),
        "next_fast_len": [{"n": n, "out": int(sfft.next_fast_len(n, real=True))}
                          for n in (1, 2, 5, 6, 7, 11, 13, 17, 31, 97, 121, 127, 1000, 1001, 4099)],
        "gaussian": gaussian_cases(),
        "map_coords": map_coords_cases(),
        "powell": powell_cases(),
        "brent": brent_cases(),
        "least_squares": lsq,
        "chain_sparsity": chain,
        "permutations": perms,
        "groups": groups,
        "lsmr": lsmr_cases(),
        "roots": roots_cases(),
        "pairs": pair_cases(),
    }
    parts = []
    for key, val in data.items():
        if isinstance(val, list):
            rows = ",\n".join("  " + json.dumps(c, separators=(",", ":")) for c in val)
            parts.append(f" {json.dumps(key)}: [\n{rows}\n ]")
        else:
            parts.append(f" {json.dumps(key)}: {json.dumps(val, separators=(',', ':'))}")
    with open(os.path.join(OUT, "numerics_cases.json"), "w", encoding="utf-8", newline="\n") as f:
        f.write("{\n" + ",\n".join(parts) + "\n}\n")


if __name__ == "__main__":
    main()
