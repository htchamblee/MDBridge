# Optimal-superposition RMSD as MDATK calls it — behaviour specification (MDAnalysis 2.11.0)

Clean-room behaviour description of MDAnalysis 2.11.0's RMSD function
(`MDAnalysis.analysis.rms.rmsd`) for the one call MDATK v0.0.4.1 makes (the rmsd analysis,
targets without a fit selection, superposition on). Derived from running MDAnalysis 2.11.0 with
NumPy 2.5.5 on x86-64 and from the published method it implements. "f64" = IEEE binary64.
Probe data: `tests/data/analysis/rmsd_cases.json`, written by `tests/parity/make_rmsd_cases.py`
(`--report` prints the comparison with an exact SVD superposition quoted in §7).

## 1. The call

MDATK calls the function once per frame and target as

    rmsd(a, b, weights=w, center=True, superposition=True)

* `a`: the target's atom positions in the current frame, f32 positions widened to f64, shape (n, 3);
* `b`: the same atoms' positions in the reference frame, same shape and type;
* `w`: either absent (None) or the atoms' masses as f64, length n (`mass_weighted`);
* n ≥ 1 (MDATK skips targets that match no atoms before calling).

The result is a Python `float` (an f64), also when it is NaN. No other value is returned and the
inputs are not modified.

Only this combination is specified. (With `superposition` off the function computes a plain
weighted RMSD; MDATK never calls it that way: it computes its non-superposed and fit-selection
RMSDs itself.)

## 2. Input checks and errors, in order

1. If the shapes of `a` and `b` differ: `ValueError` with the message
   `a and b must have same shape`. This is checked first (also before any weight check).
2. Centring (§3) uses the weights as given. A weights array whose shape is not (n,) — wrong
   length, empty, two-dimensional (n, 1) — raises `ValueError` with the message
   `Shape of weights must be consistent with shape of a along specified axis.`
3. If the weights sum to exactly zero (all zeros, or mixed signs that cancel, or an empty weight
   array with empty coordinates): `ZeroDivisionError` with the message
   `Weights sum to zero, can't be normalized`.

The function has a second length check of its own (message
`weights must have same length as a and b`) which centring always pre-empts in this call: it is
unreachable here. All of these are unreachable from MDATK (masses of the same atom group, at least
one atom); MDATK cannot produce a zero mass sum unless every atom of the target has mass 0
(MDBridge's tables give virtual sites mass 0), in which case the analysis fails with the
`ZeroDivisionError` above.

Empty inputs (n = 0, no weights) give NaN with four numpy warnings, the pair
`RuntimeWarning: Mean of empty slice`, `RuntimeWarning: invalid value encountered in divide`
once per structure (they are only logged by Python's warning machinery; MDATK never reaches this).

Shapes other than (n, 3) with equal shapes are not specified (undefined in MDAnalysis).

## 3. Centring and weights

Both structures are centred with the **raw** weights (before normalisation); then the weights are
normalised for the superposition.

Unweighted: for each coordinate column k, the centre is c_k = (Σ_i x_ik) / n, where the sum is
accumulated sequentially in atom order (plain left-to-right additions in f64) and divided once by n.
Each coordinate becomes x_ik − c_k.

Weighted: c_k = (Σ_i (x_ik · w_i)) / S, where each product x_ik · w_i is rounded to f64 and the
products are summed sequentially in atom order; S = Σ_i w_i is summed with **numpy's pairwise
summation** (the same algorithm as numpy's 1-D `sum`: eight interleaved accumulators for blocks of
up to 128 elements, recursive halving above; sequential below eight elements). Each coordinate
becomes x_ik − c_k.

Normalised weights: m = S' / n with S' again numpy's pairwise sum of the raw weights (the same
value as S), and w'_i = w_i / m (one f64 division per weight). These w'_i (mean 1) are used in §4.

Consequences (observed):

* Negative weights are accepted (only a zero sum is an error); results are finite.
* Uniform scaling of all weights (1e-300 … 1e300) changes nothing beyond the last bits.
* Weights that are all equal do **not** give the unweighted result bit for bit (the weighted
  inner products are evaluated differently, §4.1): in the probes 10 of 24 equal-weight cases
  differ from the unweighted call in the last bits. An implementation must not route equal
  weights to the unweighted path.
* A NaN or infinite weight, or a NaN / infinite coordinate, makes the result NaN (no exception).
  An infinite coordinate additionally emits one `RuntimeWarning: invalid value encountered in
  subtract` (from centring), an infinite weight three `RuntimeWarning: invalid value encountered in
  divide`; NaN coordinates or a NaN weight emit nothing. MDATK does not inspect warnings.

## 4. The superposition RMSD (QCP)

The function returns the minimum, over proper rotations R (no reflections), of

    sqrt( Σ_i w'_i |R·a_i − b_i|² / n )      (w'_i = 1 when unweighted)

of the centred coordinates, computed with the quaternion characteristic polynomial (QCP) method:
D. L. Theobald, "Rapid calculation of RMSD using a quaternion-based characteristic polynomial",
Acta Cryst. A 61 (2005) 478–480, and P. Liu, D. K. Agrafiotis, D. L. Theobald, "Fast determination
of the optimal rotational matrix for macromolecular superpositions", J. Comput. Chem. 31 (2010)
1561–1563. Theobald's reference C implementation (`qcprot.c`, from the authors' site) is
BSD-licensed; MDAnalysis' numbers are those of that reference algorithm's RMSD-only path evaluated
in f64, with the settings below. The implementer may follow the papers and this section; porting
the BSD C code itself would be a permissive-licence port that keeps its notice (LICENSING.md), an
owner decision. Only the RMSD is needed (no rotation matrix).

All arithmetic is f64, each operation rounded (no fused multiply-add: contraction changes the
results of degenerate inputs, §6, and the last bits elsewhere), sqrt correctly rounded.
Below, x, y, z are the components of the centred `b` (second argument, the reference in MDATK)
and of the centred `a` (first argument), atom i, written b_i = (bx, by, bz), a_i = (ax, ay, az).

### 4.1 Inner products

Accumulated sequentially in atom order, each starting at 0:

Unweighted:

* G_b = Σ_i ((bx·bx + by·by) + bz·bz)
* G_a = Σ_i ((ax·ax + ay·ay) + az·az)
* the nine sums S_pq = Σ_i b_p · a_q for p, q ∈ {x, y, z}
  (S_xy = Σ bx·ay, S_yx = Σ by·ax, and so on: the first index is the reference's component).

Weighted (w'_i as in §3): first u = (w'_i·bx, w'_i·by, w'_i·bz) (three rounded products), then

* G_b = Σ_i ((ux·bx + uy·by) + uz·bz)
* G_a = Σ_i w'_i · ((ax·ax + ay·ay) + az·az)
* S_pq = Σ_i u_p · a_q.

E0 = (G_b + G_a) · 0.5 (sum first, then the multiplication).

### 4.2 Characteristic polynomial

The largest eigenvalue λ of Theobald's symmetric 4×4 key matrix built from the S_pq is the largest
root of P(λ) = λ⁴ + C2·λ² + C1·λ + C0. With the squares Sxx² = S_xx·S_xx etc.:

* C2 = −2 · (Sxx² + Syy² + Szz² + Sxy² + Syx² + Sxz² + Szx² + Syz² + Szy²) — summed left to right
  in exactly this order, then multiplied by −2.0.
* C1 = 8 · (Sxx·Syz·Szy + Syy·Szx·Sxz + Szz·Sxy·Syx − Sxx·Syy·Szz − Syz·Szx·Sxy − Szy·Syx·Sxz)
  — each triple product left to right, the six terms added / subtracted left to right, then
  multiplied by 8.0.
* C0 = det of the key matrix, evaluated in Theobald's (2005) grouped form. With
  * D1 = Sxy² + Sxz² − Syx² − Szx² (left to right),
  * E = Syy² + Szz² − Sxx² + Syz² + Szy² (left to right), F = 2·(Syz·Szy − Syy·Szz),
  * the sums and differences p_xz = Sxz + Szx, p_yz = Syz + Szy, p_xy = Sxy + Syx,
    m_yz = Syz − Szy, m_xz = Sxz − Szx, m_xy = Sxy − Syx, p_xx = Sxx + Syy, m_xx = Sxx − Syy,

  C0 is the sum, left to right, of these six products:
  1. D1 · D1
  2. (E + F) · (E − F)
  3. (−p_xz·m_yz + m_xy·(m_xx − Szz)) · (−m_xz·p_yz + m_xy·(m_xx + Szz))
  4. (−p_xz·p_yz − p_xy·(p_xx − Szz)) · (−m_xz·m_yz − p_xy·(p_xx + Szz))
  5. (p_xy·p_yz + p_xz·(m_xx + Szz)) · (−m_xy·m_yz + p_xz·(p_xx + Szz))
  6. (p_xy·m_yz + m_xz·(m_xx − Szz)) · (−m_xy·p_yz + m_xz·(p_xx − Szz))

  where −u·v means (−u)·v and each bracket is evaluated left to right.

### 4.3 Newton iteration

λ starts at E0 (an upper bound of the root). Each step, with the current λ:

* x2 = λ·λ; β = (x2 + C2)·λ; α = β + C1;
* δ = (α·λ + C0) / ((2·x2)·λ + β + α)   — numerator P(λ), denominator P′(λ), both in this
  Horner-like grouping; the denominator's additions left to right;
* λ ← λ − δ.

After each step the iteration stops if |λ_new − λ_old| < |1e-14 · λ_new|. At most **50** steps are
made; if the test never succeeds, the 50th iterate is used **silently** (no warning, no error, no
fallback). A division 0/0 or x/0 inside the iteration is not trapped: it produces NaN or ±inf,
which then propagates (a NaN λ never satisfies the test, so all 50 steps run).

The tolerance is 1e-14 relative. (Some releases of the reference C code use a looser 1e-11; the
results here need 1e-14.) Typical step counts in the probes: 3–4 for general structures, 1 for
identical structures, 4–9 for coplanar sets, 17–50 or "not converged" for two atoms and collinear
sets (§6).

### 4.4 Result

rmsd = sqrt( | 2·(E0 − λ) / n | ): the difference E0 − λ, times 2.0, divided by n (the atom count,
converted to f64; with weights it is still n, because the normalised weights sum to n), absolute
value (it guards tiny negative values from rounding), square root. Returned as a Python float.

## 5. Identical and near-identical structures

For b = a (MDATK's reference frame against itself), translated copies and rigidly rotated copies,
the exact answer is 0 but the function returns 0 **or a small positive value**: E0 − λ is a
difference of two numbers of size E0 that agree to rounding, so the result is of order
sqrt(k·ε·E0/n) with ε = 2⁻⁵² and small k. Observed:

| Inputs | Results |
|---|---|
| random identical sets, 3 to 3000 atoms, 450 probes | 0.0 in 2–38 % (fewer for more atoms), otherwise 1e-7 … 5e-6 Å (max 5.0e-6 at n = 3) |
| demo heavy atoms (514), reference frame vs itself, ref 0 / 5 / 23 | 7.89e-7 / 2.38e-7 / 5.32e-7 Å; mass-weighted 5.32e-7 / 3.37e-7 / 4.76e-7 Å |
| demo C4′ (24 atoms), ref 0 / 5 / 23 | 2.75e-7 / 3.37e-7 / 1.95e-7 Å; mass-weighted 1.95e-7 / 0.0 / 2.75e-7 Å |
| demo strand-1 heavy atoms, ref 0 | 0.0 (both) |

The iteration stops after one step there. Near-identical inputs (noise 1e-9 … 1e-6 Å) are
dominated by the same rounding floor: results differ from the exact RMSD by up to ~1e-6 Å.

## 6. Degenerate inputs: one atom, two atoms, collinear and coincident sets

When the optimal rotation is not unique, the top eigenvalue of the key matrix is (at least)
double, P′ vanishes at the root, Newton's method loses its quadratic convergence, and the results
are governed by rounding and the 50-step limit:

| Inputs | Observed |
|---|---|
| one atom (n = 1), or weights that leave one nonzero weight | the centred coordinates are usually exactly 0, all sums are 0, the first step is 0/0: **NaN**. When centring leaves rounding residue (weighted centring, or coordinates whose mean is inexact) the result is a tiny number or 0.0 instead (e.g. 1.78e-15; demo: the single C4′ target with masses gives 0.0 in 6 of 24 frames, NaN otherwise) |
| coincident atoms (all positions equal) | identical copy: NaN; against a different structure: a finite value off the exact RMSD (by 4e-7 Å in the probes) |
| two atoms (always collinear) | 1–50 steps or not converged (≈40 % of probes hit the limit); about 20 % NaN; otherwise errors up to 10 Å against the exact RMSD (identical pairs: 0.0, 4e-4 … 9e-4 Å, or NaN); demo two-C4′ target against its own reference frame: 5.52e-4 / 5.54e-4 / 4.51e-4 Å |
| three or more collinear atoms | 17–50 steps or not converged (≈20 %); errors up to ~2e-3 Å (identical copies: 1e-4 … 1e-3 Å instead of 0) |
| coplanar atoms | not degenerate (4–9 steps), errors ≤ 1e-12 relative |

These results are **chaotic in the input bits**: in a probe, moving one coordinate of a collinear
set by one ulp changed the result by up to 35 Å. No tolerance makes an independent method agree
with them; only an evaluation in exactly the order of §3–§4 reproduces them (verified bit for bit,
below). MDATK reaches them with targets of one or two atoms (for example the "Backbone" target
`name C4'` of a one- or two-residue group) and in principle with collinear atom sets.

## 7. Tolerances for an independent implementation

A faithful f64 evaluation of §3–§4 (sequential sums, the stated groupings, no FMA, numpy's
pairwise sum for the weight sums) reproduced MDAnalysis **bit for bit** on all 232 non-error
cases of `rmsd_cases.json` and on 3 000 further random probes (1 to 3 000 atoms, weighted and
unweighted, general, planar, collinear, identical, two-atom). Exact equality is therefore the
target; the C++ tests should compare exactly (NaN = NaN) and fall back to the tolerance below only
for non-degenerate cases if a platform cannot be made exact.

Spread of MDAnalysis against an exact superposition (numpy SVD, residual summed directly after
rotating; `make_rmsd_cases.py --report`), and of an independent QCP variant (C0 as a numpy 4×4
determinant, numpy sums) — the error an implementation that does not follow §4's order sees.
"Ratio" is |r_MDA² − r_other²| / (ε · E0/n), E0/n being the case's `scale`:

| Family | max \|r_MDA − r_SVD\| (Å) | max ratio vs SVD | max ratio vs other QCP |
|---|---|---|---|
| general, noise 0 … 1e-6 Å | 3.3e-6 | 24.8 | 14 |
| general, noise 1e-3 … 5 Å | 4.5e-10 | 4.9 | 63 |
| identical / translated / rotated copies | 1.4e-6 | 24.4 | 14 |
| demo series (heavy, all, C4′, P, Mg, three C4′) | 7.9e-7 | 18.1 | — |
| planar | 2.2e-13 | 5.4 | — |
| collinear / two atoms / coincident / one atom | up to 6 Å, NaN | 10⁸ … 10¹⁶ | NaN, 10⁷ Å |

Recommended comparison for **non-degenerate** cases (`degenerate` false):

    | r_MDA² − r² | ≤ 1024 · ε · E0/n        (ε = 2⁻⁵², E0/n = `scale`)

i.e. at least 16 times the largest spread observed. For the demo (E0/n ≈ 155–166 Å²) this is
3.6e-11 Å² on the squared RMSD: up to 6e-6 Å absolute near 0, 6e-11 Å (relative 2e-10) at
0.3 Å. An absolute tolerance of 1e-6 Å near 0 (the milestone plan's guess) is **not** enough: an
exact method returns 0 where MDAnalysis returns up to 5e-6 Å. A relative tolerance alone fails for
the same reason. **Degenerate** cases (`degenerate` true): exact only (NaN must stay NaN); a
platform that cannot reproduce them should report them as known differences rather than widen the
tolerance.

## 8. The probe data (`tests/data/analysis/rmsd_cases.json`)

* `cases`: `a`, `b` (flat lists of 3n coordinates, `a_shape` / `b_shape` give the shapes),
  `weights` (null or a list, `weights_shape`), `result` (null on error) or `error`
  `{type, message}`, `warnings` (`"Category: message"` in emission order), `scale` (E0/n of the
  centred inputs with normalised weights, computed with numpy; null on errors), `degenerate`,
  `group` (the family). Families: `general_noise_*` (n = 3 … 300, none / uniform / mass /
  equal weights; inputs rounded through f32 like MDATK's), `identical`, `translated`, `rotated`,
  `planar`, `planar_both`, `magnitude` (coordinates 1e-6 … 1e6), `one_atom`, `two_atoms`,
  `collinear`, `collinear_f32`, `coincident`, `weights_few`, `nonfinite_*`, `negative_weights`,
  `zero_sum_weights`, `tiny_weights`, `huge_weights`, `empty`, `error_shape`, `error_weights`.
* `demo`: `topology/demo_dry.parm7` + `io/trajectories/demo_dry.nc` (24 frames). Per target:
  `selection` (for information), `atoms` (`[first, last]` index runs — the C++ test uses these,
  not the selection), `weights` (null or the parm7 masses as f64), `reference_frame` (0, 5 or 23)
  and `values` (one per frame, in frame order): exactly MDATK's loop — reference positions read at
  the reference frame, then rmsd(frame positions, reference positions) for every frame.
* All floats are `float.hex()` strings; non-finite values are `nan`, `inf`, `-inf`.
