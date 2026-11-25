"""Test-side adjustments for MDBridge's deliberate improvements of the analyses (docs/IMPROVEMENTS.md
section 4 and SEL-17; docs/DIFFERENCES.md "Analyses"). Pure Python / numpy on top of MDATK and
MDAnalysis' public API: the generators record what MDATK computes with these items fixed, so the
committed expectations are MDBridge's improved behaviour. Nothing here is MDAnalysis code.

Used by make_rmsd_cases.py (rmsd_expected) and make_analysis_cases.py / make_pipeline_cases.py
(apply_mdatk, before they import from mdatk.analysis).

Items:
  AN-16  RMSD of degenerate sets (fewer than 3 atoms, collinear or coincident atoms, by the weighted
         inner-product matrix): the exact superposition RMSD (SVD, residual summed after rotating;
         one atom: 0) instead of the QCP iteration's value.
  AN-17  RMSD where E0 - lambda < 4096 eps E0 (rounding of the iteration and cancellation, e.g. the
         reference frame against itself): the exact residual; identical inputs give exactly 0.
  AN-18  rmsd fit selections of fewer than 3 or collinear atoms: translation only, one warning.
  AN-06  energy: use_parm7 with a non-parm7 topology: a warning.
  AN-07  energy: exclude_bonded > 0 and a topology without bonds: a warning, and the exclusions
         follow the distance-guessed bonds (AN-14's guesser) of the first analysed frame; the warning
         cites the covalent radii of the guess (CIT-03).
  AN-08  energy: exclusions follow any path of at most exclude_bonded bonds.
  AN-25  energy without parm7 parameters: Lennard-Jones parameters of Amber's residue templates (OpenMM's
         amber14 ff14SB / OL15 / OL3 / TIP3P conversion) for residues that match one (energy_parameters.py);
  AN-26  the other atoms: Amber atom types by element, hydrogens by their bonded atom, CT / OH / OS / OW by
         neighbours, carbon's for unknown elements (warning);
  AN-27  without topology charges: the charges of completely matched templates (warning for the rest);
  CIT-04 one log line naming the parameters' sources and references; metadata lj_atoms, charge_atoms,
         lj_references, lj_default_elements; lj_table (MDATK 0.1.1's key) without parm7 parameters.
  AN-09  H-bonds: a guessed hydrogen goes to the donor nearest by the pair search's distance.
  AN-10  H-bonds: nitrogens carrying a guessed hydrogen are not acceptors.
  AN-13  atoms: bonds guessed from the first analysed frame.
  AN-14  atoms: bonds guessed from covalent radii (sum + 0.4 A), consecutive residues only within one
         chain, never from a one-atom residue.
  CIT-05 atoms: when the bonds are guessed, a warning citing the covalent radii (as CIT-03's).
  AN-15  proximity: members are 0 also when no structure atoms remain.
  AN-19  use_pbc without a usable cell: one warning per run.
  AN-20  non-finite coordinates: such atoms are in no pair; one warning per run of frames.
  AN-21  structural view: unknown residue groups are named in a warning.
  SEL-17 the analyses' selections are evaluated on the first analysed frame whatever the setup
         order, and distance keywords in them honour analyses.use_pbc.

From MDATK 0.1.0 (same owner; its code ported):
  H-bonds by distance only (no hydrogens, or use_angle off): each donor-acceptor atom pair counts
         once, not once per direction (hbonds.py process).
  Contacts of a group with itself: per_res_a / per_res_b are each residue's row plus its column;
         sequence neighbours (|resindex difference| <= 1) are left out only within one chain
         (structure.chain_ids); the atom_pairs label says "atom contacts" without heavy_only.
  Structure trace: the cache key's groups are every residue group, not only the shown ones.
"""

import math
import os

import numpy as np

EPS = float(np.finfo(float).eps)
COLLINEAR = 1e-10
CANCELLATION = 4096.0  # AN-17 threshold in eps * E0 (rmsd_qcp.cpp kCancellation)
NEAR = 4.0  # the AN-17 decision is not predicted within this factor of the threshold

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


# --- RMSD (AN-16, AN-17) ------------------------------------------------------------------------

def _centred(x, w):
    return x - (x * w[:, None]).sum(axis=0) / w.sum()


def _collinear(x, w):
    """AN-16: q <= 1e-10 p^2 for the weighted inner-product matrix of the centred set."""
    if len(x) < 3:
        return True
    m = (w[:, None, None] * (x[:, :, None] * x[:, None, :])).sum(axis=0)
    p = m[0, 0] + m[1, 1] + m[2, 2]
    q = (m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]) + (m[0, 0] * m[2, 2] - m[0, 2] * m[2, 0]) + \
        (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
    return bool(q <= COLLINEAR * p * p)


def rotate_rows(x, r):
    """x @ r for rows x (..., 3) and a 3 x 3 matrix r, as elementwise products and sums in a fixed
    order. numpy's matmul goes to OpenBLAS, whose kernels (chosen per CPU) round differently, which
    made the generated data depend on the machine (CI run 45551493216)."""
    x = np.asarray(x, dtype=np.float64)
    r = np.asarray(r, dtype=np.float64)
    return np.stack([(x[..., 0] * r[0, k] + x[..., 1] * r[1, k]) + x[..., 2] * r[2, k] for k in range(3)],
                    axis=-1)


def _top_eigenvector(k):
    """The unit eigenvector of the largest eigenvalue of the symmetric 4 x 4 matrix k: cyclic Jacobi
    rotations in plain Python floats (ties to the lowest index), the operations and order of
    topEigenvector in src/analysis/rmsd_qcp.cpp."""
    k = [[float(v) for v in row] for row in k]
    v = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    for _sweep in range(64):
        off = 0.0
        for p in range(4):
            for q in range(p + 1, 4):
                off += k[p][q] * k[p][q]
        if off == 0.0:
            break
        for p in range(4):
            for q in range(p + 1, 4):
                if k[p][q] == 0.0:
                    continue
                theta = (k[q][q] - k[p][p]) / (2.0 * k[p][q])
                t = (1.0 if theta >= 0.0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for r in range(4):  # k <- J^T k J
                    krp, krq = k[r][p], k[r][q]
                    k[r][p] = c * krp - s * krq
                    k[r][q] = s * krp + c * krq
                for r in range(4):
                    kpr, kqr = k[p][r], k[q][r]
                    k[p][r] = c * kpr - s * kqr
                    k[q][r] = s * kpr + c * kqr
                k[p][q] = 0.0
                k[q][p] = 0.0
                for r in range(4):  # v <- v J
                    vrp, vrq = v[r][p], v[r][q]
                    v[r][p] = c * vrp - s * vrq
                    v[r][q] = s * vrp + c * vrq
    best = 0
    for i in range(1, 4):
        if k[i][i] > k[best][best]:
            best = i
    q = [v[r][best] for r in range(4)]
    norm = math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3])
    return [x / norm for x in q]


def exact_rmsd(a, b, w=None):
    """The exact superposition RMSD: the optimal rotation from the top eigenvector of Horn's key
    matrix, the residual summed after rotating; one atom 0, identical inputs exactly 0. Plain
    Python floats in the order of exactRmsd in src/analysis/rmsd_qcp.cpp (no LAPACK / BLAS, whose
    last bits depend on the CPU), so the value is the same on every machine."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    n = len(a)
    if n == 1 or (a.shape == b.shape and np.array_equal(a, b)):
        return 0.0
    w = np.ones(n) if w is None else np.asarray(w, dtype=np.float64)
    wn = (w / w.mean()).tolist()
    ac = _centred(a, w).tolist()
    bc = _centred(b, w).tolist()
    m = [[0.0] * 3 for _ in range(3)]  # sum w a b^T
    for i in range(n):
        for r in range(3):
            for c in range(3):
                m[r][c] += wn[i] * (ac[i][r] * bc[i][c])
    k = [
        [m[0][0] + m[1][1] + m[2][2], m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0]],
        [m[1][2] - m[2][1], m[0][0] - m[1][1] - m[2][2], m[0][1] + m[1][0], m[2][0] + m[0][2]],
        [m[2][0] - m[0][2], m[0][1] + m[1][0], -m[0][0] + m[1][1] - m[2][2], m[1][2] + m[2][1]],
        [m[0][1] - m[1][0], m[2][0] + m[0][2], m[1][2] + m[2][1], -m[0][0] - m[1][1] + m[2][2]],
    ]
    q0, q1, q2, q3 = _top_eigenvector(k)
    rot = [
        [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q1 * q3 + q0 * q2)],
        [2.0 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2.0 * (q2 * q3 - q0 * q1)],
        [2.0 * (q1 * q3 - q0 * q2), 2.0 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3],
    ]
    total = 0.0
    for i in range(n):
        ai, bi = ac[i], bc[i]
        d2 = 0.0
        for c in range(3):
            d = ((rot[c][0] * ai[0] + rot[c][1] * ai[1]) + rot[c][2] * ai[2]) - bi[c]
            d2 += d * d
        total += wn[i] * d2
    return math.sqrt(total / n)


def rmsd_expected(a, b, w, mda):
    """MDBridge's value for rms.rmsd(a, b, weights=w, center=True, superposition=True) whose
    MDAnalysis value is `mda`: (value, item, near). item is None where MDBridge keeps the QCP value
    bit for bit, else "AN-16" / "AN-17"; near marks AN-17 decisions within a factor NEAR of the
    threshold (either value is MDBridge's)."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if mda is None or a.ndim != 2 or a.shape != b.shape or len(a) == 0 or \
            not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return mda, None, False
    n = len(a)
    wa = np.ones(n) if w is None else np.asarray(w, dtype=np.float64)
    if not np.all(np.isfinite(wa)) or np.any(wa < 0) or wa.sum() == 0:
        return mda, None, False
    wn = wa / wa.mean()
    ac, bc = _centred(a, wa), _centred(b, wa)
    if _collinear(ac, wn) or _collinear(bc, wn):
        return exact_rmsd(a, b, w), "AN-16", False
    if not math.isfinite(mda):  # the iteration failed (NaN): AN-16 as well
        return exact_rmsd(a, b, w), "AN-16", False
    e0 = (np.sum(wn * np.sum(ac * ac, axis=1)) + np.sum(wn * np.sum(bc * bc, axis=1))) / 2
    gap = mda * mda * n / 2  # E0 - lambda as QCP had it
    threshold = CANCELLATION * EPS * e0
    if gap < threshold / NEAR:
        return exact_rmsd(a, b, w), "AN-17", False
    if gap < threshold * NEAR:
        return exact_rmsd(a, b, w), "AN-17", True
    return mda, None, False


# --- MDATK with the analysis improvements ---------------------------------------------------------

def _covalent_radii():
    path = os.path.join(ROOT, "data", "sources", "covalent_radii", "cordero2008.tsv")
    out = {}
    with open(path, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            out[fields[1].upper()] = float(fields[2])
    return out


CARBON_RADIUS = 0.76
BOND_TOLERANCE = 0.4
# energy: metadata "lj_table" without parm7 parameters (MDATK 0.1.1's key; src/analysis/energy.cpp kLjTable)
LJ_TABLE = ("Amber force-field residue templates (OpenMM amber14: ff14SB, OL15, OL3, TIP3P with Joung-Cheatham and "
            "Li-Merz ions), then Amber parm10 atom types, Joung-Cheatham and Li-Merz (CM) ions for TIP3P "
            "(docs/GENERIC_LJ.md)")
CORDERO_CITATION = "(covalent radii: Cordero et al., Dalton Trans. 2008, 2832-2838, doi:10.1039/b801115j)"


def guess_bonds(atoms):
    """AN-14: bonds among `atoms` (local index pairs, sorted) from their current positions: within
    r_i + r_j + 0.4 A (covalent radii, carbon's for unknown elements), never H-H, never an atom of a
    one-atom residue, within a residue or between consecutive residues of one segment and chain ID."""
    from mdatk.selections import atom_elements
    from mdatk.spatial import pairs_within

    if atoms.n_atoms < 2:
        return np.zeros((0, 2), dtype=np.int64)
    radii = _covalent_radii()
    el = atom_elements(atoms)
    r = np.array([radii.get(str(e), CARBON_RADIUS) for e in el])
    u = atoms.universe
    size = np.bincount(u.atoms.resindices, minlength=len(u.residues))
    xyz = atoms.positions.astype(np.float64)
    i, j, d = pairs_within(xyz, xyz, 2 * float(r.max()) + BOND_TOLERANCE, None)
    keep = i < j
    i, j, d = i[keep], j[keep], d[keep]
    res = atoms.resindices
    seg = atoms.segindices
    chain = np.asarray(atoms.chainIDs, dtype=object) if hasattr(atoms, "chainIDs") else None
    same_chain = seg[i] == seg[j]
    if chain is not None:
        same_chain &= chain[i] == chain[j]
    near = (res[i] == res[j]) | ((np.abs(res[i] - res[j]) == 1) & same_chain)
    hh = (el[i] == "H") & (el[j] == "H")
    single = (size[res[i]] < 2) | (size[res[j]] < 2)
    ok = near & ~hh & ~single & (d <= r[i] + r[j] + BOND_TOLERANCE)
    out = np.stack([i[ok], j[ok]], axis=1).astype(np.int64)
    return out[np.lexsort((out[:, 1], out[:, 0]))] if len(out) else out.reshape(0, 2)


def bonded_exclusions(u, a, b, depth, bonds=None):
    """AN-08: keys a_local * n_b + b_local of A-B pairs joined by <= depth bonds through any atoms
    (`bonds`: index pairs, default the topology's)."""
    from collections import deque

    if depth <= 0:
        return np.empty(0, dtype=np.int64)
    if bonds is None:
        try:
            bonds = u.bonds.indices
        except Exception:  # noqa: BLE001 - no bond information
            return np.empty(0, dtype=np.int64)
    if len(bonds) == 0:
        return np.empty(0, dtype=np.int64)
    adj = {}
    for x, y in bonds:
        adj.setdefault(int(x), []).append(int(y))
        adj.setdefault(int(y), []).append(int(x))
    b_local = {int(ix): n for n, ix in enumerate(b.ix)}
    keys = set()
    for n, s in enumerate(a.ix):
        s = int(s)
        seen = {s: 0}
        q = deque([s])
        while q:
            cur = q.popleft()
            if seen[cur] >= depth:
                continue
            for nb in adj.get(cur, ()):
                if nb in seen:
                    continue
                seen[nb] = seen[cur] + 1
                q.append(nb)
                if nb in b_local:
                    keys.add(n * b.n_atoms + b_local[nb])
    return np.array(sorted(keys), dtype=np.int64)


def _usable_cell(dims):
    if dims is None:
        return False
    dims = np.asarray(dims, dtype=np.float32)
    if dims.shape != (6,) or not np.all(np.isfinite(dims)):
        return False
    return bool(np.all(dims[:3] > 0) and np.all(dims[3:] > 0) and np.all(dims[3:] < 180))


def apply_mdatk():
    """Patch MDATK's analyses (in this process) to the improved behaviour of the items above."""
    from MDAnalysis.analysis import rms as mda_rms

    # the selection improvements (SEL-01 ... SEL-15, improvements_selections.py) and the periodic
    # geometry (AN-01 ... AN-05, AN-11, AN-12, improvements_periodic.py) first: the analyses select with
    # MDBridge's selections, and the items below wrap the periodic functions, whatever the generator
    # applies (both are idempotent)
    import improvements_periodic
    import improvements_selections
    import mdatk.selections as mdatk_selections

    improvements_selections.patch_mdatk(mdatk_selections)
    improvements_periodic.apply_mdatk()
    import improvements_core_io  # chains from bonds only when they hold the residues together (TOP-01)

    improvements_core_io.apply_chains()
    import improvements_nucleotides  # NA-01: L-DNA names, CCD atom names, backbone guard

    improvements_nucleotides.apply_mdatk()

    import mdatk.analysis as analysis_pkg
    import mdatk.analysis.base as base
    import mdatk.analysis.energy as energy
    import mdatk.analysis.hbonds as hbonds
    import mdatk.analysis.rmsd as rmsd
    import mdatk.analysis.structure as structure
    import mdatk.pipeline as pipeline
    from mdatk.spatial import kabsch as mdatk_kabsch
    from mdatk.spatial import pairs_within

    improvements_selections.patch_mdatk_analysis_base(base)

    # SEL-17: distance keywords are periodic only with use_pbc (MDBridge's selections otherwise,
    # with their messages: SEL-19, SEL-22)
    def select(self, text, allow_empty=False):
        return improvements_selections.select_mdbridge(mdatk_selections, self.u, text, self.project.groups,
                                                       allow_empty, periodic=bool(self.use_pbc))

    base.AnalysisContext.select = select

    # SEL-17 / AN-13: every setup on the first analysed frame; AN-19 / AN-20: frame warnings
    def run_analyses(ctx, analyses, frame_indices, reporter=None):
        reporter = reporter or ctx.reporter
        frame_indices = np.asarray(frame_indices, dtype=np.int64)
        n = len(frame_indices)
        traj = ctx.u.trajectory
        for a in analyses:
            reporter.log(f"Preparing {a.label} ...")
            if n and traj.ts.frame != frame_indices[0]:
                traj[int(frame_indices[0])]
            a.setup(frame_indices)
        times = np.empty(n, dtype=np.float64)
        every = max(1, n // 200)
        reporter.log(f"Analysing {n} frame(s): {', '.join(a.label for a in analyses)}")
        state = {"pbc": False, "range": None}

        def flush():
            if state["range"] is None:
                return
            first, last, count, worst = state["range"]
            atoms = "1 atom has" if worst == 1 else f"{worst} atoms have"
            if count == 1:
                ctx.warn(f"frame {first}: {atoms} non-finite coordinates; they are left out of distance searches")
            else:
                ctx.warn(f"frames {first}-{last}: up to {atoms} non-finite coordinates; they are left out of "
                         "distance searches")
            state["range"] = None

        dt_ps = ctx.project.inputs.dt_ps  # LD-05 (MDATK 0.1.0): the override is the time axis, frame x dt
        for k, ts in enumerate(traj[frame_indices]):
            reporter.check_cancel()
            times[k] = frame_indices[k] * float(dt_ps) if dt_ps else ts.time
            if ctx.use_pbc and not state["pbc"] and not _usable_cell(ts.dimensions):
                state["pbc"] = True
                ctx.warn(f"use_pbc is on but trajectory frame {ts.frame} has no usable periodic cell; "
                         "distances are not periodic")
            bad = int(np.count_nonzero(~np.all(np.isfinite(ts.positions), axis=1)))
            if bad:
                r = state["range"]
                state["range"] = (ts.frame, ts.frame, 1, bad) if r is None else (r[0], ts.frame, r[2] + 1,
                                                                                 max(r[3], bad))
            else:
                flush()
            for a in analyses:
                a.process(k, ts)
            if k % every == 0 or k == n - 1:
                reporter.progress((k + 1) / n, f"Analysing frame {k + 1}/{n}")
        flush()
        return [a.finalize() for a in analyses], times

    base.run_analyses = run_analyses
    analysis_pkg.run_analyses = run_analyses
    pipeline.run_analyses = run_analyses

    # AN-20: non-finite points are in no pair (open path; the periodic path already leaves them out)
    import mdatk.spatial as spatial
    original_pairs = spatial.pairs_within

    def pairs_within_finite(a, b, cutoff, box=None):
        a = np.asarray(a, dtype=np.float64)
        b = np.asarray(b, dtype=np.float64)
        if spatial._usable_box(box) is not None or len(a) == 0 or len(b) == 0:
            return original_pairs(a, b, cutoff, box)
        fa = np.flatnonzero(np.all(np.isfinite(a), axis=1))
        fb = np.flatnonzero(np.all(np.isfinite(b), axis=1))
        if len(fa) == len(a) and len(fb) == len(b):
            return original_pairs(a, b, cutoff, box)
        if len(fa) == 0 or len(fb) == 0:
            return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0, np.float64)
        i, j, d = original_pairs(a[fa], b[fb], cutoff, box)
        return fa[i], fb[j], d

    for module in (spatial, energy, hbonds, structure):
        module.pairs_within = pairs_within_finite
    for name in ("contacts", "mg", "stacking", "basepairs"):
        module = getattr(analysis_pkg, name, None)
        if module is not None and hasattr(module, "pairs_within"):
            module.pairs_within = pairs_within_finite

    # AN-14 / AN-13 (the setup now runs on the first analysed frame)
    guessed = {"bonds": False}

    def atom_bonds(atoms, max_len=None):
        index = {int(ix): k for k, ix in enumerate(atoms.ix)}
        pairs = []
        try:
            for bond in atoms.bonds:
                i, j = (int(x) for x in bond.indices)
                if i in index and j in index:
                    pairs.append((index[i], index[j]))
        except Exception:  # noqa: BLE001 - topology without bonds
            pairs = []
        if pairs:
            return np.array(sorted({(min(p), max(p)) for p in pairs}), dtype=np.int64).reshape(-1, 2)
        guessed["bonds"] = True
        return guess_bonds(atoms)

    structure.atom_bonds = atom_bonds

    # CIT-05: the atoms analysis says when its bonds are guessed, citing the radii (its last warning
    # of the setup, after the trace points')
    atoms_setup = structure.AtomsAnalysis.setup

    def atoms_setup_warning(self, frame_indices):
        guessed["bonds"] = False
        atoms_setup(self, frame_indices)
        self._bonds_guessed = guessed["bonds"]
        if guessed["bonds"]:
            self.ctx.warn("structural view: the topology has no bonds among the structure atoms; bonds are "
                          f"guessed from distances on the first analysed frame {CORDERO_CITATION}")

    structure.AtomsAnalysis.setup = atoms_setup_warning

    # CIT-06: the atoms result records what it cites (Cordero et al.'s radii when bonds were guessed),
    # after MDATK's meta keys, before label / analysis
    atoms_finalize = structure.AtomsAnalysis.finalize

    def atoms_finalize_cited(self):
        res = atoms_finalize(self)
        meta = {}
        for key, value in res.meta.items():
            if key == "label":
                meta["citation_keys"] = ["cordero2008"] if self._bonds_guessed else []
            meta[key] = value
        res.meta = meta
        return res

    structure.AtomsAnalysis.finalize = atoms_finalize_cited

    # AN-15
    def proximity_process(self, k, ts):
        self.dist[self.member, k] = 0.0
        if self.struct.n_atoms == 0:
            return
        i, _, d = structure.pairs_within(self.struct.positions, self.set_atoms.positions, self.radius,
                                         self.ctx.box(ts))
        if len(i):
            best = np.full(len(self.residues), np.inf)
            np.minimum.at(best, self.row[i], d)
            self.dist[:, k] = best
        self.dist[self.member, k] = 0.0

    structure.ProximityAnalysis.process = proximity_process

    # AN-21
    trace_setup = structure.TraceAnalysis.setup

    def trace_setup_warning(self, frame_indices):
        known = set(self.project.group_names())
        for name in self.project.visual.structure_groups:
            if name not in known:
                self.ctx.warn(f"structural view: unknown residue group {name!r} ignored")
        trace_setup(self, frame_indices)

    structure.TraceAnalysis.setup = trace_setup_warning

    # AN-06 / AN-07 / AN-08; AN-25 - AN-27 and CIT-04 (parameters without parm7)
    import energy_parameters as eparams

    energy._bonded_exclusions = bonded_exclusions
    energy_process = energy.EnergyAnalysis.process
    energy_finalize = energy.EnergyAnalysis.finalize

    def count(n):
        return f"{n} atom" if n == 1 else f"{n} atoms"

    def energy_setup_improved(self, frame_indices):
        energy.Analysis.setup(self, frame_indices)
        p = self.project.analyses.energy
        self.p = p
        a = self.ctx.group_atoms(p.group_a)
        b = energy.without(self.ctx.group_atoms(p.group_b), a)
        if b.n_atoms == 0:
            raise ValueError(f"energy: group {p.group_b!r} lies entirely inside {p.group_a!r}")
        self.a, self.b = a, b
        u = self.ctx.u
        warn = self.ctx.warn
        topology_charges = hasattr(u.atoms, "charges")
        if topology_charges:
            self.qa = a.charges.astype(np.float64)
            self.qb = b.charges.astype(np.float64)
            self.charge_source = "topology"
        path = self.ctx.topology_path
        if p.use_parm7 and path and not energy.is_parm7(path):
            warn(f"energy: use_parm7 is set but the topology {path} is not a parm7 file; "
                 "Lennard-Jones parameters come from force-field residue templates and Amber atom types")
        self.lj = None
        if p.use_parm7 and path and energy.is_parm7(path):
            try:
                self.lj = energy.LJTable.from_parm7(path)
                if len(self.lj.type_index) != u.atoms.n_atoms:
                    raise ValueError("atom count mismatch")
                self.lj_source = "parm7"
            except Exception as exc:  # noqa: BLE001 - as MDATK
                warn(f"energy: could not read LJ parameters from parm7 ({exc}); Lennard-Jones parameters come from "
                     "force-field residue templates and Amber atom types")
                self.lj = None
        na, nb = a.n_atoms, b.n_atoms
        self.lj_atoms = {"parm7": 0, "template": 0, "typed": 0}
        self.charge_atoms = {"topology": na + nb if topology_charges else 0, "template": 0, "none": 0}
        self.lj_references = []
        self.lj_default_elements = []
        self.lj_guessed_bonds = False
        charged = 0
        if self.lj is not None:
            self.lj_atoms["parm7"] = na + nb
        else:
            indices = [int(x) for x in a.ix] + [int(x) for x in b.ix]
            par = eparams.energy_parameters(u, indices, mdatk_selections.atom_elements, guess_bonds)
            self.ra = np.array(par["rstar"][:na], dtype=np.float64)
            self.rb = np.array(par["rstar"][na:], dtype=np.float64)
            self.ea = np.array(par["epsilon"][:na], dtype=np.float64)
            self.eb = np.array(par["epsilon"][na:], dtype=np.float64)
            self.lj_atoms["template"] = par["layer"].count("template")
            self.lj_atoms["typed"] = par["layer"].count("typed")
            t, y = self.lj_atoms["template"], self.lj_atoms["typed"]
            self.lj_source = "typed" if t == 0 else "templates" if y == 0 else "templates+typed"
            self.lj_references = list(par["references"])
            self.lj_guessed_bonds = bool(par["guessed_bonds"])
            self.lj_default_elements = list(par["default_elements"])
            if not topology_charges:
                self.qa = np.array(par["charge"][:na], dtype=np.float64)
                self.qb = np.array(par["charge"][na:], dtype=np.float64)
                charged = sum(1 for c in par["template_charge"] if c)
                self.charge_atoms = {"topology": 0, "template": charged, "none": na + nb - charged}
            if not topology_charges and charged == 0:
                warn("energy: topology has no partial charges; electrostatics will be zero")
            elif not topology_charges and charged < na + nb:
                missing = na + nb - charged
                rc = par["chargeless_residues"]
                warn(f"energy: the topology has no partial charges; {charged} of {na + nb} atoms take the charges "
                     "of the force-field residue template they match completely (hydrogens included); the other "
                     f"{count(missing)} ({rc} {'residue' if rc == 1 else 'residues'}) have no charge, so their "
                     "electrostatic terms are zero")
            if par["default_elements"]:
                listed = ", ".join(repr(e) for e in par["default_elements"])
                warn(f"energy: no Lennard-Jones parameters for element(s) {listed} ({count(par['default_atoms'])}); "
                     "carbon's (Amber type C) are used")
            line = "energy: Lennard-Jones parameters without parm7: "
            if t > 0:
                line += f"{count(t)} from Amber force-field residue templates ({', '.join(par['force_fields'])})"
            if y > 0:
                line += (", " if t > 0 else "") + f"{count(y)} from Amber atom types by element"
            if topology_charges:
                line += "; partial charges from the topology"
            elif charged > 0:
                line += f"; partial charges of {count(charged)} from the templates"
            else:
                line += "; no partial charges"
            line += ". Cite: " + "; ".join(eparams.short_citation(k) for k in par["references"])
            self.ctx.reporter.log(line)
        if not topology_charges:
            if charged == 0:
                if self.lj is not None:
                    warn("energy: topology has no partial charges; electrostatics will be zero")
                self.qa = np.zeros(na)
                self.qb = np.zeros(nb)
                self.charge_atoms = {"topology": 0, "template": 0, "none": na + nb}
            self.charge_source = ("none" if charged == 0 else
                                  "templates" if self.charge_atoms["none"] == 0 else "partial")
        self.excl = energy._bonded_exclusions(u, a, b, p.exclude_bonded)
        n = self.n_frames
        self.e_el = np.zeros(n)
        self.e_vdw = np.zeros(n)
        self.force = np.zeros(n)
        self.n_pairs = np.zeros(n, dtype=np.int64)
        try:
            has_bonds = len(self.ctx.u.bonds) > 0
        except Exception:  # noqa: BLE001
            has_bonds = False
        self._guess_exclusions = p.exclude_bonded > 0 and not has_bonds
        self._guess_exclusions_cited = self._guess_exclusions  # CIT-06 (process() resets the flag)
        if self._guess_exclusions:
            warn("energy: the topology has no bonds; bonded pairs are excluded along bonds guessed from "
                 f"distances on the first analysed frame {CORDERO_CITATION}")

    def energy_finalize_improved(self):
        res = energy_finalize(self)
        meta = {}
        for key, value in res.meta.items():
            if key == "label":  # MDBridge's keys before label / analysis, as newResult appends them
                if self.lj_source != "parm7":  # MDATK 0.1.1's key (energy.cpp kLjTable)
                    meta["lj_table"] = LJ_TABLE
                meta.update(lj_atoms=dict(self.lj_atoms), charge_atoms=dict(self.charge_atoms),
                            lj_references=list(self.lj_references),
                            lj_default_elements=list(self.lj_default_elements),
                            citation_keys=list(self.lj_references) +
                            (["cordero2008"] if self._guess_exclusions_cited or self.lj_guessed_bonds else []))
            meta[key] = value
        res.meta = meta
        return res

    def energy_process_improved(self, k, ts):
        if getattr(self, "_guess_exclusions", False):
            self._guess_exclusions = False
            u = self.ctx.u
            bonds = guess_bonds(u.atoms)  # u.atoms: local index = atom index
            self.excl = bonded_exclusions(u, self.a, self.b, self.p.exclude_bonded, bonds)
        energy_process(self, k, ts)

    energy.EnergyAnalysis.setup = energy_setup_improved
    energy.EnergyAnalysis.process = energy_process_improved
    energy.EnergyAnalysis.finalize = energy_finalize_improved

    # AN-09 / AN-10
    def guess_hydrogens(self, pos, box):
        heavy_idx = np.flatnonzero(self.cand_donor)
        h_idx = np.flatnonzero(self.is_h)
        i, j, d = pairs_within(pos[h_idx], pos[heavy_idx], 1.25, box)
        order = np.lexsort((j, i))  # MDBridge's pair order (i, j): ties go to the lowest donor
        best = {}
        for hi, hj, dist in zip(i[order], j[order], d[order]):
            if hi not in best or dist < best[hi][1]:
                best[hi] = (hj, dist)
        is_n = mdatk_selections.atom_elements(self.ab) == "N"  # SEL-13's, whenever patched
        for hi, (hj, _) in best.items():
            donor = int(heavy_idx[hj])
            self.h_of.setdefault(donor, []).append(int(h_idx[hi]))
            if is_n[donor]:
                self.cand_acc[donor] = False
        self._build_donor_lists()

    hbonds.HBondAnalysis._guess_hydrogens = guess_hydrogens

    # MDATK 0.1.0: H-bonds by distance only count each donor-acceptor atom pair once
    def hbonds_process(self, k, ts):
        pos = self.ab.positions.astype(np.float64)
        box = self.ctx.box(ts)
        if self._finalize_donors_pending:
            self._guess_hydrogens(pos, box)
            self._finalize_donors_pending = False
        rows_a, rows_b = [], []
        total = 0
        counted = None  # distance only: the atom pairs counted with the donor in group A
        n_atoms = len(pos)
        for s in (0, 1):
            don = self.donor_idx[s]
            acc = self.acc_idx[1 - s]
            i, j, _ = hbonds.pairs_within(pos[don], pos[acc], self.dcut, box)
            if len(i) == 0:
                continue
            d_atoms = don[i]
            a_atoms = acc[j]
            if not self.use_angle:
                if s == 0:
                    counted = d_atoms * n_atoms + a_atoms
                elif counted is not None:
                    new = ~np.isin(a_atoms * n_atoms + d_atoms, counted)
                    d_atoms, a_atoms = d_atoms[new], a_atoms[new]
            if self.use_angle:
                n_h = self.h_start[d_atoms + 1] - self.h_start[d_atoms]
                rep_d = np.repeat(d_atoms, n_h)
                rep_a = np.repeat(a_atoms, n_h)
                offs = np.repeat(self.h_start[d_atoms], n_h)
                within = np.arange(n_h.sum()) - np.repeat(np.cumsum(n_h) - n_h, n_h)
                hyd = self.h_flat[offs + within]
                v1 = hbonds.minimum_image_vectors(pos[rep_d] - pos[hyd], box)
                v2 = hbonds.minimum_image_vectors(pos[rep_a] - pos[hyd], box)
                cosang = (v1 * v2).sum(1) / np.maximum(
                    np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1), 1e-12)
                ok = cosang <= self.cos_cut  # angle >= cutoff
                d_atoms, a_atoms = rep_d[ok], rep_a[ok]
            total += len(d_atoms)
            if s == 0:
                rows_a.append(self.row[d_atoms])
                rows_b.append(self.row[a_atoms])
            else:
                rows_a.append(self.row[a_atoms])
                rows_b.append(self.row[d_atoms])
        self.count[k] = total
        if total:
            ra = np.concatenate(rows_a).astype(np.int64)
            rb = np.concatenate(rows_b).astype(np.int64)
            uniq, cnt = np.unique(ra * self.n_b + rb, return_counts=True)
            self.csr.add(uniq // self.n_b, uniq % self.n_b, cnt)
        else:
            self.csr.add([], [], [])

    hbonds.HBondAnalysis.process = hbonds_process

    # MDATK 0.1.0: contacts of a group with itself
    import mdatk.analysis.contacts as contacts

    contacts_setup = contacts.ContactAnalysis.setup

    def contacts_setup_chains(self, frame_indices):
        contacts_setup(self, frame_indices)
        self.chain = structure.chain_ids(self.res_a) if self.same else None

    def contacts_process(self, k, ts):
        i, j, _ = contacts.pairs_within(self.a.positions, self.b.positions, self.cutoff, self.ctx.box(ts))
        if self.same and len(i):
            ri = self.a.resindices[i]
            rj = self.b.resindices[j]
            bonded = (np.abs(ri - rj) <= 1) & (self.chain[self.row_a[i]] == self.chain[self.row_b[j]])
            keep = (self.a.ix[i] < self.b.ix[j]) & ~bonded
            i, j = i[keep], j[keep]
        self.atom_pairs[k] = len(i)
        if len(i) == 0:
            self.csr.add([], [], [])
            return
        keys = self.row_a[i].astype(np.int64) * self.n_b + self.row_b[j]
        uniq, cnt = np.unique(keys, return_counts=True)
        self.residue_pairs[k] = len(uniq)
        self.csr.add(uniq // self.n_b, uniq % self.n_b, cnt)

    contacts_finalize = contacts.ContactAnalysis.finalize

    def contacts_finalize_improved(self):
        res = contacts_finalize(self)
        if self.same:  # the map holds each pair once (i < j): a residue's contacts are its row plus its column
            res.arrays["per_res_a"] = res.arrays["per_res_b"] = res.arrays["per_res_a"] + res.arrays["per_res_b"]
        if not self.project.analyses.contacts.heavy_only:
            series = res.meta["outputs"]["series"]["atom_pairs"]
            series["label"] = series["label"].replace("heavy-atom contacts", "atom contacts")
        return res

    base_output_label = contacts.ContactAnalysis.output_label.__func__

    def contacts_output_label(cls, project, out):
        label = base_output_label(cls, project, out)
        if out == "atom_pairs" and not project.analyses.contacts.heavy_only:
            label = label.replace("heavy-atom contacts", "atom contacts")  # hydrogens are counted too
        return label

    if not getattr(contacts.ContactAnalysis, "_mdbridge_self_group", False):
        contacts.ContactAnalysis.setup = contacts_setup_chains
        contacts.ContactAnalysis.process = contacts_process
        contacts.ContactAnalysis.finalize = contacts_finalize_improved
        contacts.ContactAnalysis.output_label = classmethod(contacts_output_label)
        contacts.ContactAnalysis._mdbridge_self_group = True

    # MDATK 0.1.0: the structure trace's cache key covers every residue group
    structure.TraceAnalysis.selection_texts = classmethod(
        lambda cls, project: project.group_names() + [project.analyses.mg.mg_selection])

    # AN-16 / AN-17 / AN-18
    class _Rms:
        @staticmethod
        def rmsd(a, b, weights=None, center=True, superposition=True):
            value = mda_rms.rmsd(a, b, weights=weights, center=center, superposition=superposition)
            return rmsd_expected(a, b, weights, value)[0]

    rmsd.rms = _Rms
    translation_only = set()

    def kabsch(mobile, reference, weights=None):
        r, mc, rc = mdatk_kabsch(mobile, reference, weights)
        if id(reference) in translation_only:
            r = np.eye(3)
        return r, mc, rc

    rmsd.kabsch = kabsch
    rmsd_setup = rmsd.RMSDAnalysis.setup

    def rmsd_setup_improved(self, frame_indices):
        rmsd_setup(self, frame_indices)
        for t, spec in zip(self.targets, self.project.analyses.rmsd.targets):
            if t["fit"] is None or not self.superposition:
                continue
            x = t["ref_fit"]
            w = np.ones(len(x)) if t["wf"] is None else t["wf"]
            collinear = len(x) < 3 or (w.sum() > 0 and _collinear(x - (x * w[:, None]).sum(axis=0) / w.sum(), w))
            if collinear:
                translation_only.add(id(t["ref_fit"]))
                n = len(x)
                what = (f"{n} atom" if n == 1 else f"{n} atoms") if n < 3 else f"{n} collinear atoms"
                self.ctx.warn(f"RMSD target {spec.label!r}: fit selection has {what}; "
                              "superposition is translation only")

    rmsd.RMSDAnalysis.setup = rmsd_setup_improved

    # the versions of the analyses whose results these items change (energy, hbonds, rmsd,
    # proximity, atoms): MDBridge's, from the one table of them
    import analysis_versions
    analysis_versions.apply_mdatk()
