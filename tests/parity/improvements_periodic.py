"""MDBridge's periodic geometry (docs/IMPROVEMENTS.md AN-01 to AN-05, AN-11, AN-12), written
independently in numpy for the parity generators: the expected values of the improved behaviour.

MDATK computes periodic distances with MDAnalysis' capped_distance and corrects difference vectors
for orthorhombic cells only; MDBridge replaces both with one float64 minimum image for every cell
shape (docs/DIFFERENCES.md). The generators keep recording MDATK / MDAnalysis and use this module
for what MDBridge must give instead:

  AN-01 / AN-02 / AN-03  pairs_within_periodic: the minimum-image distance of every pair, kept iff
                         d^2 <= cutoff^2, each (i, j) once, sorted by (i, j), no errors
  AN-04                  minimum_image_vectors: the same image for any cell; one orthogonality test
                         (the three float32 angles exactly 90)
  AN-05                  energy net force: r_hat times the force at the clamped distance
  AN-11                  base pairs: minimum-image donor-acceptor and C1'-C1' distances with use_pbc
  AN-12                  stacking: ring atoms at the minimum image of the ring's first atom with
                         use_pbc

apply_mdatk() patches MDATK's modules in-process (spatial helpers in every analysis module, the
energy, stacking and base-pair steps above) and gives MDATK's analyses MDBridge's versions (one table
for every improvement: analysis_versions.py; periodic geometry changes contacts, energy, hbonds, mg,
stacking, proximity and basepairs), so that a generator records MDBridge's expected results and cache
keys. improvements_analyses.apply_mdatk() calls it first, so that its items (AN-20's finite-point
filter, the energy exclusions) wrap the periodic functions.

Every float operation is the one src/analysis/spatial.cpp does, in the same order (numpy float64,
np.rint = nearbyint, half to even), so the values agree bit for bit; the box vectors of triclinic
cells use math.cos / math.sin (libm), as the C++ code does.

    import improvements_periodic
    improvements_periodic.apply_mdatk()   # after MDATK is importable
"""

import math

import numpy as np


MAX_REDUCTIONS = 64
_RAD = 3.14159265358979323846 / 180.0


def _usable(box):
    """MDATK's _usable_box: a float64 array of shape (6,) with lengths > 0 (NaN passes), else None."""
    if box is None:
        return None
    box = np.asarray(box, dtype=np.float64)
    if box.shape != (6,) or np.any(box[:3] <= 0):
        return None
    return box


class Cell:
    """PeriodicCell of src/analysis/spatial.cpp."""

    def __init__(self, orthogonal, periodic, vectors):
        self.orthogonal = orthogonal
        self.periodic = periodic
        self.vectors = vectors  # 3 x 3 float64, rows a, b, c

    @staticmethod
    def of(box):
        box = _usable(box)
        if box is None:
            return None
        length = [float(np.float32(x)) for x in box[:3]]
        angle = [float(np.float32(x)) for x in box[3:]]
        vec = np.zeros((3, 3))
        if all(x == 90.0 for x in angle):
            periodic = [math.isfinite(x) and x > 0 for x in length]
            for k in range(3):
                vec[k, k] = length[k] if periodic[k] else 0.0
            return Cell(True, periodic, vec)
        for k in range(3):
            if not (math.isfinite(length[k]) and length[k] > 0 and angle[k] > 0 and angle[k] < 180):
                return None
        cos = lambda d: 0.0 if d == 90.0 else math.cos(d * _RAD)  # noqa: E731
        sin = lambda d: 1.0 if d == 90.0 else math.sin(d * _RAD)  # noqa: E731
        a, b, c = length
        ca, cb, cg = cos(angle[0]), cos(angle[1]), cos(angle[2])
        sg = sin(angle[2])
        bx = b * cg
        by = b * sg
        cx = c * cb
        cy = c * (ca - cb * cg) / sg
        rad = c * c - cx * cx - cy * cy
        cz = math.sqrt(rad) if rad >= 0 else math.nan
        if not (by > 0 and cz > 0 and math.isfinite(cz)):
            return None
        vec[:] = [[a, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]]
        return Cell(False, [True, True, True], vec)

    def minimum_image(self, v):
        """Rows of v (n x 3 float64) mapped to their minimum images (a new array)."""
        v = np.array(v, dtype=np.float64, copy=True).reshape(-1, 3)
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            if self.orthogonal:
                for k in range(3):
                    if not self.periodic[k]:
                        continue
                    length = self.vectors[k, k]
                    half = length / 2
                    for _ in range(MAX_REDUCTIONS):
                        m = np.abs(v[:, k]) > half
                        if not m.any():
                            break
                        v[m, k] = v[m, k] - length * np.rint(v[m, k] / length)
                return v
            for k in (2, 1, 0):
                u = self.vectors[k]
                half = u[k] / 2
                for _ in range(MAX_REDUCTIONS):
                    m = np.abs(v[:, k]) > half
                    if not m.any():
                        break
                    n = np.rint(v[m, k] / u[k])
                    for c in range(k + 1):
                        v[m, c] = v[m, c] - n * u[c]
            va, vb, vc = self.vectors
            best = v.copy()
            bs = _sq(v)
            for i in (-1.0, 0.0, 1.0):
                for j in (-1.0, 0.0, 1.0):
                    for k in (-1.0, 0.0, 1.0):
                        if i == 0 and j == 0 and k == 0:
                            continue
                        w = np.empty_like(v)
                        w[:, 0] = ((v[:, 0] + i * va[0]) + j * vb[0]) + k * vc[0]
                        w[:, 1] = (v[:, 1] + j * vb[1]) + k * vc[1]
                        w[:, 2] = v[:, 2] + k * vc[2]
                        s = _sq(w)
                        upd = s < bs
                        best[upd] = w[upd]
                        bs[upd] = s[upd]
            return best


def _sq(v):
    return (v[:, 0] * v[:, 0] + v[:, 1] * v[:, 1]) + v[:, 2] * v[:, 2]


def minimum_image_vectors(vectors, box=None):
    """AN-04: MDATK's minimum_image_vectors, for every cell with periodic geometry."""
    cell = Cell.of(box)
    vectors = np.asarray(vectors, dtype=np.float64)
    if cell is None:
        return vectors
    return cell.minimum_image(vectors).reshape(vectors.shape)


def pairs_within_periodic(a, b, cutoff, box):
    """AN-01 / AN-02 / AN-03: (i, j, d) of every pair within the cutoff under the minimum image,
    sorted by (i, j). a and b are rounded to float32 (as MDATK passes them) and widened."""
    a = np.asarray(a, dtype=np.float32).astype(np.float64).reshape(-1, 3)
    b = np.asarray(b, dtype=np.float32).astype(np.float64).reshape(-1, 3)
    cutoff = float(cutoff)
    bound = cutoff * cutoff
    cell = Cell.of(box)
    ok_a = np.isfinite(a).all(axis=1)
    ok_b = np.isfinite(b).all(axis=1)
    out_i, out_j, out_d = [], [], []
    if len(a) and len(b) and not math.isnan(bound):
        chunk = max(1, 1_000_000 // max(1, len(b)))
        for s in range(0, len(a), chunk):
            diff = (a[s:s + chunk, None, :] - b[None, :, :]).reshape(-1, 3)
            if cell is not None:
                diff = cell.minimum_image(diff)
            sq = _sq(diff).reshape(-1, len(b))
            keep = (sq <= bound) & ok_a[s:s + chunk, None] & ok_b[None, :]
            i, j = np.nonzero(keep)
            out_i.append(i + s)
            out_j.append(j)
            out_d.append(np.sqrt(sq[i, j]))
    if not out_i:
        return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0, np.float64)
    return (np.concatenate(out_i).astype(np.int64), np.concatenate(out_j).astype(np.int64),
            np.concatenate(out_d).astype(np.float64))


# --- MDATK patches ----------------------------------------------------------------------------------

def apply_mdatk():
    """Patch MDATK's analyses in-process to MDBridge's improved periodic behaviour (idempotent)."""
    from mdatk import secstruct, spatial
    from mdatk.analysis import ANALYSES, basepairs, contacts, energy, hbonds, mg, stacking, structure

    if getattr(spatial, "_mdbridge_periodic", False):
        return
    spatial._mdbridge_periodic = True
    original_pairs = spatial.pairs_within
    # MDATK's own functions, for generators that record them as the reference behaviour
    spatial.mdatk_pairs_within = original_pairs
    spatial.mdatk_minimum_image_vectors = spatial.minimum_image_vectors

    def pairs_within(a, b, cutoff, box=None):
        if len(a) == 0 or len(b) == 0 or _usable(box) is None:
            return original_pairs(a, b, cutoff, box)
        return pairs_within_periodic(a, b, cutoff, box)

    for module in (spatial, contacts, energy, hbonds, mg, stacking, structure):
        if hasattr(module, "pairs_within"):
            module.pairs_within = pairs_within
        if hasattr(module, "minimum_image_vectors"):
            module.minimum_image_vectors = minimum_image_vectors

    import analysis_versions  # the versions of every improvement (stacking 2 among them)
    analysis_versions.apply_mdatk()

    _patch_energy(energy)
    _patch_stacking(stacking)
    _patch_basepairs(basepairs, secstruct)


def _pow(x, exponent):
    """x ** exponent elementwise through libm's pow, as numpy's generic power loop and energy.cpp's
    std::pow compute it. numpy's own `**` dispatches to SVML's pow on AVX-512 CPUs, whose results
    differ in the last bits from libm's, which made the generated energies depend on the machine
    (CI run 45551493216)."""
    return np.array([math.pow(v, exponent) for v in np.asarray(x, dtype=np.float64).tolist()],
                    dtype=np.float64)


def _patch_energy(energy):
    """AN-05: the net force is r_hat times the force at the clamped distance (0 for r = 0)."""
    from mdatk.analysis import ANALYSES

    cls = ANALYSES["energy"]
    COULOMB = energy.COULOMB

    def process(self, k, ts):
        p = self.p
        box = self.ctx.box(ts)
        pa = self.a.positions.astype(np.float64)
        pb = self.b.positions.astype(np.float64)
        i, j, d = energy.pairs_within(pa, pb, p.cutoff, box)
        if len(self.excl) and len(i):
            keep = ~np.isin(i * self.b.n_atoms + j, self.excl)
            i, j, d = i[keep], j[keep], d[keep]
        self.n_pairs[k] = len(i)
        if len(i) == 0:
            return
        r = d
        d = np.maximum(d, 0.5)
        qq = COULOMB * self.qa[i] * self.qb[j] / p.dielectric
        if p.distance_dependent:
            e_el = qq / (d * d)
            de_el = -2.0 * qq / _pow(d, 3.0)
        else:
            e_el = qq / d
            de_el = -qq / (d * d)
        if self.lj is not None:
            acoef, bcoef = self.lj.pair_coefficients(self.a.ix[i], self.b.ix[j])
        else:
            rmin = self.ra[i] + self.rb[j]
            eps = np.sqrt(self.ea[i] * self.eb[j])
            r6 = _pow(rmin, 6.0)
            acoef = eps * r6 * r6
            bcoef = 2.0 * eps * r6
        inv6 = 1.0 / _pow(d, 6.0)
        e_vdw = acoef * inv6 * inv6 - bcoef * inv6
        de_vdw = (-12.0 * acoef * inv6 * inv6 + 6.0 * bcoef * inv6) / d
        self.e_el[k] = e_el.sum()
        self.e_vdw[k] = e_vdw.sum()
        rvec = energy.minimum_image_vectors(pa[i] - pb[j], box)
        with np.errstate(divide="ignore", invalid="ignore"):
            fmag = np.where(r > 0, -(de_el + de_vdw) / r, 0.0)
        total = (rvec * fmag[:, None]).sum(axis=0).tolist()
        # np.linalg.norm of a 3-vector is sqrt(x . x) through BLAS' ddot (CPU-dependent kernels):
        # the sum in energy.cpp's order instead
        self.force[k] = math.sqrt((total[0] * total[0] + total[1] * total[1]) + total[2] * total[2])

    cls.process = process


_STATE = {"cell": None}


def _patch_stacking(stacking):
    """AN-12: with use_pbc, ring atoms at the minimum image of the ring's first atom."""
    from mdatk.analysis import ANALYSES

    cls = ANALYSES["stacking"]
    original_process = cls.process
    original_geometry = stacking.ring_geometry

    def ring_geometry(positions):
        cell = _STATE["cell"]
        if cell is not None and positions.shape[1] > 1:
            positions = np.array(positions, dtype=np.float64, copy=True)
            x0 = positions[:, :1, :]
            rest = (positions[:, 1:, :] - x0).reshape(-1, 3)
            positions[:, 1:, :] = x0 + cell.minimum_image(rest).reshape(positions.shape[0], -1, 3)
        return original_geometry(positions)

    def process(self, k, ts):
        _STATE["cell"] = Cell.of(self.ctx.box(ts))
        try:
            return original_process(self, k, ts)
        finally:
            _STATE["cell"] = None

    stacking.ring_geometry = ring_geometry
    cls.process = process


def _patch_basepairs(basepairs, secstruct):
    """AN-11: with use_pbc, every distance of the builtin detector is a minimum-image distance."""
    from mdatk.analysis import ANALYSES

    cls = ANALYSES["basepairs"]
    original_process = cls.process
    detector = secstruct.PairDetector
    original_candidates = detector.candidates
    lo, hi = secstruct.C1_WINDOW

    def norm(diff, cell):
        shape = diff.shape
        v = cell.minimum_image(diff.reshape(-1, 3)).reshape(shape)
        return np.linalg.norm(v, axis=-1)

    def candidates(self, xyz):
        cell = _STATE["cell"]
        if cell is None:
            return original_candidates(self, xyz)
        # MDATK's PairDetector.candidates (mdatk/secstruct.py) with minimum-image distances
        xyz = np.asarray(xyz, dtype=np.float64)
        out = []
        for kind, need, rows_p, ip, c1_p, rows_q, iq, c1_q in self.rules:
            q = xyz[np.maximum(iq, 0)]
            q_missing = iq < 0
            chunk = max(1, 2_000_000 // max(1, iq.size))
            for s in range(0, len(ip), chunk):
                pi = ip[s:s + chunk]
                p = xyz[np.maximum(pi, 0)]
                d = norm(p[:, None, :, :] - q[None, :, :, :], cell)
                d[(pi < 0)[:, None, :] | q_missing[None, :, :]] = np.inf
                a, b = np.nonzero((d <= self.cutoff).sum(axis=-1) >= need)
                ca, cb = c1_p[s + a], c1_q[b]
                both = (ca >= 0) & (cb >= 0)
                dc = norm(xyz[np.maximum(ca, 0)] - xyz[np.maximum(cb, 0)], cell)
                ok = ~both | ((dc >= lo) & (dc <= hi))
                a, b = a[ok], b[ok]
                if len(a) == 0:
                    continue
                score = np.sort(d[a, b], axis=-1)[:, :2].sum(axis=-1)
                dev = np.where(both[ok], np.abs(dc[ok] - secstruct.C1_IDEAL), secstruct.C1_IDEAL)
                out.append(np.column_stack([rows_p[s + a], rows_q[b], np.full(len(a), kind), score, dev]))
        return np.concatenate(out) if out else np.zeros((0, 5))

    def process(self, k, ts):
        _STATE["cell"] = Cell.of(self.ctx.box(ts))
        try:
            return original_process(self, k, ts)
        finally:
            _STATE["cell"] = None

    detector.candidates = candidates
    cls.process = process
