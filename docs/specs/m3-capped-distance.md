# Periodic pair search (`capped_distance`) — behaviour specification (MDAnalysis 2.11.0)

> **Reference behaviour only.** MDBridge no longer reproduces this function: its periodic search is
> an exact float64 minimum image for every cell shape, sorted by (i, j), without these errors
> (docs/DIFFERENCES.md, "Periodic geometry"; AN-01 to AN-04). This spec records what MDATK and
> MDAnalysis do, which distance_cases.json still lists beside MDBridge's `exact` expectations.

Clean-room behaviour description of MDAnalysis 2.11.0's public function `capped_distance`, as MDATK
v0.0.4.1 calls it, derived by reading and running MDAnalysis 2.11.0 (NumPy 2.5.5, x86-64, Linux).
"f32"/"f64" are IEEE binary32/binary64; fl32(x) and fl64(x) mean x rounded to nearest (ties to
even) in that format; ulp32(x) is the spacing of f32 numbers at |x|. Every table below comes from
actual runs; most rows are cases of `tests/data/analysis/distance_cases.json` (generator
`tests/parity/make_distance_cases.py`), cited by their `name`.

## 1. The call and the inputs MDATK can produce

MDATK's periodic pair search (used by contacts, H-bonds and their hydrogen guess at 1.25 Å, Mg
(two searches), stacking (centroids against centroids, the same array on both sides), energy and
proximity, when periodic boundaries are switched on and the frame has a usable cell) makes exactly
one call:

* **reference** A: n points, **configuration** B: m points; both converted to f32 by MDATK before
  the call (atom positions are already f32; centroids are f64 and are rounded to f32).
* **cutoff** c: a Python float (f64), the analysis parameter as given. Only the contact cutoff is
  validated (> 0); the others can be 0, negative, very large, infinite or NaN.
* **cell**: six f64 numbers [a, b, c, α, β, γ] (lengths in Å, angles in degrees), the frame's cell;
  in practice f32 values widened to f64 (trajectory cells are stored in f32).
* no lower cutoff, automatic method choice, distances requested.

Inputs that never reach the call: n = 0 or m = 0 (MDATK returns three empty arrays itself), a
cell that is absent, not six numbers, or has any length ≤ 0 (MDATK then uses its non-periodic
search instead). Everything else is reachable: angles of any value (0, 180, negative, NaN),
lengths that are tiny, huge or infinite, NaN / infinite coordinates, duplicated points, and the
same points on both sides.

**What the function returns** (MDATK converts and passes them on unchanged, in the same order):

* `pairs`: an n_pairs × 2 array of 64-bit signed integers; row k is (i, j) with i an index into A
  and j an index into B. With no pair: shape (0, 2).
* `distances`: n_pairs f64 values, `distances[k]` belonging to row k. With no pair: shape (0,).
* MDATK returns `pairs[:,0]` and `pairs[:,1]` as int64 arrays and the distances as f64.

Each index pair (i, j) occurs **at most once** (one image per pair — the nearest one; §4.4); the
order of the rows depends on the method (§5). The inputs are never modified.

**Errors.** Three `ValueError`s are possible, all from the grid method (§3), never from brute
force; MDATK does not catch them, so the analysis fails with this message:

| condition (grid method only) | message (verbatim) |
|---|---|
| c < 0 | `Cutoff must be positive` |
| c larger than the grid's limit for the cell (§3.2) | `Cutoff <C> too large for box (max <M>)` |
| a cell length that is exactly 0 after rounding to f32 | `Any of the box dimensions cannot be 0` |

`<C>` is Python's `str()` of the f64 cutoff (`13.0`, `1.3`, `inf`); `<M>` is the f32 limit
printed as a Python float (`12.5`, `1.25`, `8.945474624633789`). The third message is
unreachable from MDATK (a positive f32 length stays positive) and is listed for completeness.
A NaN cutoff raises nothing.

**Warnings.** For a cell whose angles admit no third vector (the radicand of §2's c_z is negative,
e.g. [25, 25, 25, 60, 60, 150] or [25, 25, 25, 30, 30, 60]), NumPy emits the RuntimeWarning
`invalid value encountered in sqrt`, once per call when n < 10 or m < 10 and twice otherwise
(the method rule builds the vectors too). Python prints a
warning at most once per code location per process; it does not change the result. MDBridge need
not reproduce it.

## 2. The cell as single-precision box vectors

All six cell numbers are first rounded to f32. The cell is **orthogonal** iff the three rounded
angles are exactly 90 (so 90.000001 counts as 90, while 90.00001 → f32 90.00000763 and 89.99999
do not; cases `angle 90.000001`, `angle 90.00001`, `angle 89.99999`). Orthogonal cells use the
three rounded lengths La, Lb, Lc directly.

Otherwise the cell is turned into three box vectors in lower-triangular form, computed in f64 from
the f32-rounded numbers and then each component rounded to f32:

* a = (a, 0, 0)
* b = (b·cos γ, b·sin γ, 0)
* c = (c_x, c_y, c_z) with c_x = c·cos β, c_y = c·(cos α − cos β·cos γ)/sin γ,
  c_z = sqrt(c² − c_x² − c_y²)

Angles are converted to radians by multiplying by π/180 (f64); an angle of exactly 90 uses
cos = 0 (and sin γ = 1) instead of the trigonometric value. The cell is **invalid** — all nine
components 0 — unless all six numbers are > 0, each angle is < 180, and the computed c_z is > 0
(NaN fails every test). Below, a_x, b_x, b_y, c_x, c_y, c_z are the f32 components.
(cos and sin are NumPy's f64 functions; their last bit can differ from the C library's, which
can at most change a component's f32 rounding at a tie.)

## 3. Method choice (automatic)

The two methods give the same pairs in exact arithmetic but differ in rounding, in the boundary
rule (§4.2), in the order of the pairs (§5) and in their errors (§1), so the choice must be
reproduced exactly. It is made per call:

1. n < 10 **or** m < 10 → **brute force**.
2. else n·m ≥ 100 000 000 → **grid** (no other check: this is how a cutoff larger than the grid
   allows reaches it, cases `lattice 10000x10000 c=13 (too large)`, `... c=inf`).
3. else a size vector S: for an orthogonal cell S = (La, Lb, Lc); otherwise, per Cartesian axis,
   the largest minus the smallest component among the three f32 box vectors (the origin is not
   included; an invalid cell gives S = 0). If fl32(c) > fl32(fl32(0.3)·S_k) on **any** axis
   (comparison and product in f32, strictly greater) → **brute force**, else → **grid**. A NaN
   cutoff never compares greater (grid if n, m ≥ 10); +inf always does.

This is the same rule as for distance selections (`m2-selection-grid-method.md` §3, with its table
of size vectors).

| cell | cutoff | n × m | method (case) |
|---|---|---|---|
| [25,25,25,90,90,90] | 3 | 9 × 12, 12 × 9, 9 × 9, 1 × 12 | brute force (`sizes ...`) |
| same | 3 | 10 × 10 | grid |
| same | 7.4999995 … 7.5000005 | 10 × 12 | grid: fl32(0.3)·25 rounds to f32 7.5000005 |
| same | 7.500001, 7.6, 12.5 | 10 × 12 | brute force |
| [10,20,30,90,90,90] | 3.0 / 3.0000002 | 10 × 12 | grid / brute force |
| [25,25,25,60,70,80] | 6 / 7.5 | 10 × 11 | grid / brute force |
| [40,40,40,90,90,120] | 9 | 10 × 11 | grid |
| [30,30,30,109.4712×3] | 6 / 7.5 | 10 × 11 | grid / brute force |
| [20,20,20,30,30,30] | ≥ 1 | 10 × 11 | brute force (S_x small) |
| [10,10,10,90,90,30] | 1.25 / 1.3 | 10 × 11 | grid / grid → `Cutoff 1.3 too large for box (max 1.25)` |
| invalid cell (angles 0, 180, NaN, −90, 270 …) | 3 | 12 × 12 | brute force (S = 0) |
| same | 0 | 12 × 12 | grid |
| any | −1 | ≥ 10 each | grid → `Cutoff must be positive` |
| any | −1 | 4 × 12 | brute force → no pairs |

### 3.1 Grid limit (the "too large" error)

With the f32 box-vector diagonal a_x, b_y, c_z and the f32 angles, let
s(θ) = fl32(sin(fl32(θ·π/180))) (the product rounded to f32 before the f64 sine).
The limit M is half the smallest of the f32 products b_y·s(α), c_z·s(α), a_x·s(β), c_z·s(β),
a_x·s(γ), b_y·s(γ), all in f32; for an orthogonal cell it is half the shortest length. The grid
raises the error when c > M. Observed: 12.5 for a 25 Å cube, 1.25 for [10,10,10,90,90,30],
8.945474624633789 for [25,25,25,60,70,80], 0.5 for [1000,1000,1,90,90,90].

## 4. The pair set and the distances

Common to both methods: coordinates are f32; the per-component difference of two points is
computed **in single precision** (fl32 of the difference of the two f32 values) and only then
widened to f64; all further arithmetic is f64 with each operation rounded; square roots are
correctly rounded f64.

### 4.1 Orthogonal cell

**Brute force.** For each reference point i and configuration point j, per axis k with length L
(f32) and r = fl32(fl64(1/L)) (the reciprocal computed in f64, then rounded to f32):

* δ = fl32(B_j,k − A_i,k);
* if L > 2^−23: s = δ·r, t = s − round(s) with round = nearest integer, halves away from zero,
  δ' = L·t (all f64); if L ≤ 2^−23 the axis is not periodic: δ' = δ;
* d = sqrt((δ'_x² + δ'_y²) + δ'_z²).

The pair is kept iff d ≤ c (f64 comparison); the returned distance is d.

Consequence: because r is not exactly 1/L, every distance is scaled by L·r ≠ 1 (relative error up
to 2^−24), *including pairs that need no image*, and the error grows with the raw separation for
points several cells apart. Whether a point at exactly the cutoff is kept therefore depends on the
cell length:

| L (f32) | r (hex) | L·r | f32 difference exactly 3.5 / 1.25 / 4.0 along one axis → d |
|---|---|---|---|
| 25 | 0x1.47ae14p-5 | 0.99999998 | 3.4999999217689037 / 1.2499999720603228 / 3.999999910593033 (kept at c = 3.5 / 1.25 / 4) |
| 31.7 | 0x1.026c36p-5 | 0.99999997 | 3.49999989517773 / 1.249999962563475 / 3.99999988020312 (kept) |
| 30 | 0x1.111112p-5 | 1.00000005 | 3.5000001825392246 / 1.2500000651925802 / 4.000000208616257 (**dropped**) |
| 18.3 | 0x1.bfa67ap-5 | 1.00000002 | 3.5000000561796085 / … (dropped) |
| 40, 60, 77.1 | | > 1 | dropped |

(cases `axis bf L=... c=...`; their grid twins `axis grid ...` keep all of these at exactly 3.5,
1.25, 4.0.) More examples, reference point (0, 0, 0), L = 25 (`ortho images cutoff 5.0` vs
`ortho images grid cutoff 5.0`):

| configuration point | brute force d | grid d |
|---|---|---|
| (24, 0, 0) | 1.000000536441803 | 1.0 |
| (25, 0, 0), (−25, 0, 0) | 5.587935447692871e-07 | 0.0 |
| (50, 0, 0) | 1.1175870895385742e-06 | 0.0 |
| (10000, 0, 0) | 0.00022351741790771484 | 0.0 |
| (10003, 0, 0), (−10003, 0, 0) | 2.999776415526867 | 3.0 |
| (1000001, 0, 0) | 0.9776482358574867 | 1.0 |
| (28, −21, 50) | 4.999999888241291 | 5.0 |
| (12.5, 0, 0) | 12.499999720603228 (c = 12.5) | — |

**Grid.** First every point of both sets is **wrapped** into [0, L) per axis, in single
precision (each "+", "−", "·" below rounded to f32 unless marked f64):

1. If x < 0: let y = x + L. If y ≥ 0, x becomes y. Otherwise (more than one cell away)
   x becomes x − n·L with n = ⌊x·(1/L)⌋ (the reciprocal and the product in f64, n an integer,
   n·L rounded to f32), and if the result is still < 0, L is added once more; the wrap ends here.
2. Then, with v the value so far (x itself if it was ≥ 0, else y): if v ≥ L, let w = v − L; if
   w < L, x becomes w; otherwise x becomes v − n·L with n = ⌊v·(1/L)⌋ (as above), and if the
   result is still ≥ L, L is subtracted once more.

(So a tiny negative x can become exactly L in step 1 — e.g. −1e−7 + 10 = 10 in f32 — and step 2
then maps it to 0.) For |x| beyond about 2^31·L the integer n overflows and the result is
undefined (not recorded). Then, for wrapped reference point p and configuration point q, per axis:

* δ = fl32(p_k − q_k) widened to f64 (note: reference minus configuration);
* if δ > L/2: δ − L; else if δ ≤ −L/2: δ + L (f64; L/2 is exact);
* d² = (δ_x² + δ_y²) + δ_z².

The pair is kept iff **d² ≤ fl64(c·c)**; the returned distance is sqrt(d²).

A spec-writer reimplementation of exactly these rules (with §5's order) reproduced the pairs, the
order and every distance bit for bit in all 267 orthogonal cases of the probe file with finite
inputs (the cases with non-finite or > 1e6 lengths, non-finite cutoffs and the 3e38 coordinate
were not part of that check).

### 4.2 The boundary rule and its precision

* Both methods are **inclusive** at the cutoff, but they test different quantities: brute force
  compares the rounded distance fl64(sqrt(d²)) ≤ c, the grid the squared distance d² ≤ fl64(c²).
  The two disagree when d² lies within an ulp of c². Example (triclinic cell
  [40,40,40,70,80,75], where both methods compute the same d², §4.3): reference (2, 2, 2),
  configuration (4.941999912261963, 2.587775945663452, 3.607692003250122), c = 3.403750523079556
  (0x1.b3ae18de2f4eap+1): brute force keeps the pair with d = c exactly, the grid drops it
  (`square tie tric bf 8` / `square tie tric grid 8`; also 16 and 17).
* For orthogonal cells the two methods also compute different d (§4.1), so near the cutoff their
  pair sets differ more often (`square tie bf/grid 0, 6, 8, 10, 14, 16, 17`; `axis bf/grid ...`).
* The comparison is in f64, with the f64 cutoff; nothing about it is f32 except the inputs and
  the single-precision differences and wraps.

### 4.3 Triclinic cell (both methods)

Both methods compute **the same d²** for a pair (verified on every triclinic grid case of the
probe file: identical distance bits), and differ only in the test of §4.2 and the order.

1. Both point sets are wrapped into the primary cell, in single precision, so that the fractional
   coordinates (s_a, s_b, s_c) of each point satisfy 0 ≤ s < 1 (half-open, as computed): first by
   whole c vectors until 0 ≤ z < c_z, then by whole b vectors until z·(c_y/c_z) ≤ y <
   z·(c_y/c_z) + b_y (the bound in f64), then by whole a vectors until the x bound of the same kind
   holds. A point less than one cell away is moved by adding or subtracting the f32 components of
   one vector (f32 sums); a point further away is first moved by an integer multiple of each vector
   (the integer from its fractional coordinate computed in f64, the multiple and the subtraction in
   f32) and then corrected by single vectors in the same order.
2. δ = fl32 difference of the two wrapped points per component, widened to f64.
3. Among the 27 vectors δ + i·a + j·b + k·c (i, j, k ∈ {−1, 0, 1}; components summed in f64),
   the one with the smallest squared length is used — the first one strictly smaller than all
   before it, enumerating i, then j, then k from −1 to 1 (k fastest), so ties are irrelevant for
   d. d² is its squared length ((x² + y²) + z²).

For every valid cell probed (angles down to 0.5° and up to 170°, 20–80 Å, points up to 20 cells
outside), the resulting d equals the true minimum-image distance (over all lattice translations,
with the f32 box vectors, computed exactly from the f32 inputs) to within a few ulp32 of the
coordinates' magnitude; no pair was ever lost or added beyond that rounding (§6).

**Invalid cell** (all components zero, §2): nothing is wrapped, all 27 candidates are the same
vector, so d is the plain (non-periodic) distance |δ|. With c > 0 the method is always brute
force (S = 0). With c = 0 and n, m ≥ 10 the grid runs and treats the cell as **orthogonal with the
given lengths** (cases `angles [0.0, 0.0, 0.0] cutoff 0` etc.: only coincident points are
reported). A **flat** cell (valid but with c_z ≈ 0, e.g. angles 120/120/120) gives meaningless
but deterministic pairs (`angles [120.0, 120.0, 120.0] ...`).

### 4.4 Large cutoffs and points outside the cell

* Each (i, j) is reported **once**, with its nearest-image distance, even when the cutoff exceeds
  half a cell or the whole cell (brute force; the grid cannot run there). Cutoff 30 or 1e9 or inf
  in a 25 Å cube: all n·m pairs, each once (`ortho images cutoff ...`, `big cutoff self ...`,
  `cutoff inf ...`). So the result never contains several images of the same pair.
* Points outside the primary cell, even far outside, are handled by the formulas above: by
  construction of the brute-force orthogonal formula (with its growing error, table in §4.1), by
  wrapping otherwise.

### 4.5 Degenerate inputs

| input | result (cases) |
|---|---|
| A and B the same points | every (i, i) at distance 0, and both (i, j) and (j, i) (`self ...`) |
| duplicated points | each coincident pair at 0 (brute force orthogonal: points one cell apart at 5.587935447692871e-07 for L = 25) (`duplicates ...`) |
| c = 0 | only pairs with d = 0 (brute force: δ' = 0 exactly; grid: coincident after wrapping, e.g. x = 0 and x = 25 in a 25 Å cell, `cutoff 0 faces grid`; brute force keeps only exact coincidences there, `cutoff 0 faces bf`) |
| c < 0 | brute force: no pairs; grid: `Cutoff must be positive` (also −1e−9) |
| c = NaN | no pairs, no error |
| c = +inf | brute force: all pairs (each once); grid route (n·m ≥ 1e8): `Cutoff inf too large for box (max …)` |
| NaN or ±inf coordinate, orthogonal cell | that point is in no pair, both methods (`non-finite ...`) |
| NaN or ±inf coordinate, triclinic cell | **inconsistent in MDAnalysis**: brute force reports the point at distance 0 from *every* point; the grid reports spurious and even duplicated pairs, or crashes the process (non-finite configuration point). Not recorded; MDBridge should treat such a point like the orthogonal case (in no pair) |
| huge finite coordinate (1e10) | grid: correct (wrapped); brute force orthogonal: per §4.1 (1.4825820922851562 instead of 1.0, `huge 1e10 bf`); triclinic: correct |
| angles 0, 180, NaN, negative, ≥ 180, or with negative radicand | invalid cell (§4.3): plain distances |
| length ≤ 2^−23 (orthogonal) | that axis is not periodic in brute force (`lengths [1e-10, ...]`) |
| infinite length (orthogonal) | brute force: every distance NaN → no pairs; grid: that axis not periodic |

## 5. Order of the returned pairs

**Brute force:** row-major — i ascending, and for each i, j ascending.

**Grid:** i ascending (each reference point is handled completely before the next). For one
reference point, the configuration points are visited cell by cell of a grid laid over the cell:

* Grid size per axis: with the f32 box components, the perpendicular widths are
  w_x = V / sqrt((b_y·c_z)² + (b_x·c_z)² + (b_x·c_y − b_y·c_x)²), w_y = b_y·c_z / sqrt(c_y² + c_z²),
  w_z = c_z, V = a_x·b_y·c_z (f64; orthogonal: the three lengths; invalid: 0). Let
  h = fl32(fl32(max(c, 1)) + 0.001) (the sum in f64, rounded to f32). Then
  N_k = min(⌊w_k / h⌋, 1290); an axis with N_k ≤ 3 gets N_k = 1 and is **not wrapped** in the
  enumeration below; an axis with N_k ≥ 4 is.
* Cell of a wrapped point (x, y, z) (§4.1 / §4.3 wrapping), in f64 with integer parts truncated
  toward zero and then reduced modulo N (remainder with the sign of the dividend):
  k_z = ⌊z / (c_z/N_z)⌋, y' = y − z·(c_y/c_z), k_y = ⌊y' / (b_y/N_y)⌋,
  k_x = ⌊(x − y'·(b_x/b_y) − z·(c_x/c_z)) / (a_x/N_x)⌋ (orthogonal: the tilt terms are 0).
* For a reference point in cell (k_x, k_y, k_z), the cells (k_x + d_x, k_y + d_y, k_z + d_z) are
  visited with d_x, d_y, d_z ∈ (−1, 0, +1), d_x outermost and d_z innermost. On a wrapped axis −1
  becomes N − 1 and N becomes 0; on an unwrapped axis (N = 1) only offset 0 exists.
* Within a cell, configuration points come in **descending** index order.
* A visited configuration point is output when it passes the test of §4.1 / §4.3.

Examples: `axis grid L=25.0 c=3.5` (all points in one cell: j descending, 11 … 0); `order ...`
cases (several cell counts, orthogonal and triclinic); `order self 25 c=7` (25 / 7.001 < 4 cells
per axis, so every axis has one cell: for each i, j descending). The rule was checked against
every grid case of the probe file with a valid cell and finite inputs (orthogonal and triclinic,
the latter with wrapped coordinates as MDAnalysis computes them).

## 6. Precision summary and the tolerance an implementation needs

* **Orthogonal cells: exact reproduction is possible** — pairs, order and every distance bit —
  with §3, §4.1 and §5; tests can compare exactly.
* **Triclinic cells:** reproducing the wrap of §4.3 step 1 to the last bit is what makes results
  exact. An implementation that wraps differently (or computes the minimum image exactly in f64)
  differs from MDAnalysis because MDAnalysis wraps in single precision: observed
  |d_MDBridge − d_MDAnalysis| ≤ 3.1·ulp32(M) for usual cells and up to 7·ulp32(M) for a 30° cell,
  where M is the largest of the pair's absolute coordinates and the three cell lengths (≈ 1e−5 Å
  for coordinates below 64 Å). **Recommended tolerance: |Δd| ≤ 8·ulp32(M)**; pairs whose
  distance lies within that band of the cutoff may then be in or out (listed per case in
  `band`); all other pairs must match exactly.
* **Order:** brute force row-major; grid as §5. A different wrap can move a point near a grid-cell
  boundary to the neighbouring cell and change the order of that reference point's pairs;
  order-sensitive results (f64 sums over pairs, first-wins ties) then differ in the last bits.
* **The two methods are not interchangeable:** for orthogonal cells their distances differ by up
  to ~|raw separation|·2^−24 + 1 ulp32 (≈ 4e−6 Å observed for points within two 25 Å cells) and
  their boundary tests differ (§4.2); for triclinic cells only the boundary test differs. The
  method choice (§3) is therefore part of the behaviour.
* Distances are f64 values computed from f32 data: they are not f32 numbers, and they are not the
  exact distances of the f32 inputs.

## 7. The probe file

`tests/data/analysis/distance_cases.json` (format in the generator's docstring): 495 cases
with the inputs (f32 values as exact JSON numbers), the cell and cutoff as MDATK passes them, the
method the rule chose, the pairs in MDAnalysis' order and the distances as `float.hex`, or the
error; lattice cases (n·m ≥ 1e8) give point-generation parameters, the pair count and SHA-256
hashes of the pairs (little-endian int64) and distances (little-endian f64). Groups: `auto`
(method rule), `boundary` (points at and one f32 step around the cutoff along axes and across
faces, sqrt/square ties), `random` (inside, up to one cell outside and 20 cells outside;
orthogonal and triclinic cells), `images` (outside points, cutoffs beyond half a cell), `self`
(same array, duplicates), `order` (dense sets, several grid sizes), `degenerate` (§4.5) and
`large`. Each case was also run through MDATK's own pair search, which returned the same pairs
and distances.
