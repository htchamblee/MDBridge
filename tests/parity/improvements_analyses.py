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
         follow the distance-guessed bonds (AN-14's guesser) of the first analysed frame.
  AN-08  energy: exclusions follow any path of at most exclude_bonded bonds.
  AN-09  H-bonds: a guessed hydrogen goes to the donor nearest by the pair search's distance.
  AN-10  H-bonds: nitrogens carrying a guessed hydrogen are not acceptors.
  AN-13  atoms: bonds guessed from the first analysed frame.
  AN-14  atoms: bonds guessed from covalent radii (sum + 0.4 A), consecutive residues only within one
         chain, never from a one-atom residue.
  AN-15  proximity: members are 0 also when no structure atoms remain.
  AN-19  use_pbc without a usable cell: one warning per run.
  AN-20  non-finite coordinates: such atoms are in no pair; one warning per run of frames.
  AN-21  structural view: unknown residue groups are named in a warning.
  SEL-17 the analyses' selections are evaluated on the first analysed frame whatever the setup
         order, and distance keywords in them honour analyses.use_pbc.
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


def exact_rmsd(a, b, w=None):
    """The exact superposition RMSD: SVD rotation, residual summed after rotating; one atom 0,
    identical inputs exactly 0."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    n = len(a)
    if n == 1 or (a.shape == b.shape and np.array_equal(a, b)):
        return 0.0
    w = np.ones(n) if w is None else np.asarray(w, dtype=np.float64)
    wn = w / w.mean()
    ac = _centred(a, w)
    bc = _centred(b, w)
    h = (wn[:, None] * ac).T @ bc
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt)) or 1.0
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    diff = ac @ r - bc
    return float(np.sqrt(np.sum(wn * np.sum(diff * diff, axis=1)) / n))


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

    import mdatk.analysis as analysis_pkg
    import mdatk.analysis.base as base
    import mdatk.analysis.energy as energy
    import mdatk.analysis.hbonds as hbonds
    import mdatk.analysis.rmsd as rmsd
    import mdatk.analysis.structure as structure
    import mdatk.pipeline as pipeline
    from mdatk.selections import SelectionError, atom_elements, resolve_selection
    from mdatk.spatial import kabsch as mdatk_kabsch
    from mdatk.spatial import pairs_within

    # SEL-17: distance keywords are periodic only with use_pbc
    def select(self, text, allow_empty=False):
        resolved = resolve_selection(text, self.project.groups, self.u)
        try:
            ag = self.u.select_atoms(resolved, periodic=bool(self.use_pbc))
        except Exception as exc:  # noqa: BLE001 - MDAnalysis raises several exception types
            raise SelectionError(f"invalid selection {text!r} -> {resolved!r}: {exc}") from exc
        if not allow_empty and ag.n_atoms == 0:
            raise SelectionError(f"selection {text!r} matched no atoms")
        return ag

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

        for k, ts in enumerate(traj[frame_indices]):
            reporter.check_cancel()
            times[k] = ts.time
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
        return guess_bonds(atoms)

    structure.atom_bonds = atom_bonds
    structure.AtomsAnalysis.version = structure.AtomsAnalysis.version + 1

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
    structure.ProximityAnalysis.version = structure.ProximityAnalysis.version + 1

    # AN-21
    trace_setup = structure.TraceAnalysis.setup

    def trace_setup_warning(self, frame_indices):
        known = set(self.project.group_names())
        for name in self.project.visual.structure_groups:
            if name not in known:
                self.ctx.warn(f"structural view: unknown residue group {name!r} ignored")
        trace_setup(self, frame_indices)

    structure.TraceAnalysis.setup = trace_setup_warning

    # AN-06 / AN-07 / AN-08
    energy._bonded_exclusions = bonded_exclusions
    energy_setup = energy.EnergyAnalysis.setup
    energy_process = energy.EnergyAnalysis.process

    def energy_setup_improved(self, frame_indices):
        p = self.project.analyses.energy
        path = self.ctx.topology_path
        real_warn = self.ctx.warn
        notes = []
        self.ctx.warn = notes.append  # MDATK's own warnings (charges, parm7) come first
        try:
            energy_setup(self, frame_indices)
        finally:
            del self.ctx.warn
        for message in notes:
            real_warn(message)
        if p.use_parm7 and path and not energy.is_parm7(path):
            real_warn(f"energy: use_parm7 is set but the topology {path} is not a parm7 file; "
                      "using generic element parameters")
        try:
            has_bonds = len(self.ctx.u.bonds) > 0
        except Exception:  # noqa: BLE001
            has_bonds = False
        self._guess_exclusions = p.exclude_bonded > 0 and not has_bonds
        if self._guess_exclusions:
            real_warn("energy: the topology has no bonds; bonded pairs are excluded along bonds guessed from "
                      "distances on the first analysed frame")

    def energy_process_improved(self, k, ts):
        if getattr(self, "_guess_exclusions", False):
            self._guess_exclusions = False
            u = self.ctx.u
            bonds = guess_bonds(u.atoms)  # u.atoms: local index = atom index
            self.excl = bonded_exclusions(u, self.a, self.b, self.p.exclude_bonded, bonds)
        energy_process(self, k, ts)

    energy.EnergyAnalysis.setup = energy_setup_improved
    energy.EnergyAnalysis.process = energy_process_improved
    energy.EnergyAnalysis.version = energy.EnergyAnalysis.version + 1

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
        is_n = atom_elements(self.ab) == "N"
        for hi, (hj, _) in best.items():
            donor = int(heavy_idx[hj])
            self.h_of.setdefault(donor, []).append(int(h_idx[hi]))
            if is_n[donor]:
                self.cand_acc[donor] = False
        self._build_donor_lists()

    hbonds.HBondAnalysis._guess_hydrogens = guess_hydrogens
    hbonds.HBondAnalysis.version = hbonds.HBondAnalysis.version + 1

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
    rmsd.RMSDAnalysis.version = rmsd.RMSDAnalysis.version + 1
