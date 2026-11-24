# Distance selections: grid search vs brute force

Clean-room behaviour specification derived by observing MDAnalysis 2.11.0
`select_atoms`. The choice of search method is normally invisible (both
methods return the same atoms), but it decides **which errors can occur** and
a few edge results, so a faithful port must reproduce the choice. Method
choices quoted in the tables were observed by instrumenting MDAnalysis at run
time; atoms/errors are the user-visible results.

## 1. Which keywords perform a pair search, and with which point sets

| keyword | "reference" set (size n) | "configuration" set (size m) | searches |
|---|---|---|---|
| `around R sel` | atoms matched by `sel` | atoms of the searched group **not** matched by `sel` | one, cutoff R |
| `isolayer R1 R2 sel` | atoms matched by `sel` | group atoms not matched by `sel` | two independent searches: cutoff R2 (outer) and cutoff R1 (inner); method chosen separately for each |
| `sphzone R sel` | **one point**: centre of geometry of `sel` | all atoms of the group | one, cutoff R |
| `sphlayer R1 R2 sel` | one point (centre of geometry) | all atoms of the group | one, cutoff R2, keeps distances > R1 |
| `point X Y Z R` | one point (X,Y,Z) | all atoms of the group | one, cutoff R |
| `cyzone`, `cylayer` | — | — | no pair search (own geometric test, see §6) |

If either set is empty, nothing is searched and the result is empty (for
`around`/`isolayer`: when the inner selection or its complement is empty;
for the sphere keywords: when the inner selection is empty). An empty group
gives an empty result without evaluating anything.

Distances are inclusive at the outer cutoff (distance ≤ R) and exclusive at
an inner radius (distance > R1) for both methods.

## 2. The periodic box that is used

* `select_atoms` has a `periodic` option, default **on**.
* With `periodic` on, the box is the current frame's unit cell
  [a, b, c, α, β, γ] if the universe has one; if it has none (MDAnalysis
  `None`, i.e. all three lengths 0) no box is used.
* With `periodic` off, no box is used.
* All box numbers are handled in single precision.

## 3. The decision rule (made for every search separately)

Inputs: n, m (set sizes from §1), cutoff R (the user's number, converted to
a double), and the box (or none).

1. If n < 10 **or** m < 10 → **brute force**. Therefore `sphzone`,
   `sphlayer` and `point` (n = 1) are *always* brute force, whatever the
   system size.
2. Else if n · m ≥ 100 000 000 (1e8) → **grid search**, with no further
   check (this is the only route by which a cutoff larger than the box can
   reach the grid search; see §4).
3. Else compute a per-axis "size" vector S = (Sx, Sy, Sz):
   * **no box**: for each Cartesian axis, (largest coordinate) − (smallest
     coordinate) over the union of both point sets;
   * **box with α = β = γ = 90 exactly** (compared after conversion to single
     precision; 90.0001 is *not* orthogonal): S = (a, b, c);
   * **any other box**: build the three box vectors in the standard
     lower-triangular form (first vector along x, second in the xy-plane;
     an invalid box — any length ≤ 0, any angle ≤ 0 or ≥ 180 — gives three
     zero vectors), and for each Cartesian axis take (largest component) −
     (smallest component) **among the three vectors only** (the origin is
     not included).
   Then, if R > 0.3 · S on **any** axis → **brute force**; otherwise → **grid
   search**. The product 0.3·S and the comparison are carried out in single
   precision: R is rounded to float32, 0.3 is rounded to float32, and the
   float32 product is compared with "strictly greater than". Example: for an
   orthogonal box of 100 Å the threshold is float32(30.000001907…), so
   R = 30.000002 (rounds to the same float32) still uses the grid, while
   R = 30.0001 uses brute force.

Practical consequences of rule 3:

* A system that is flat or linear without a box (extent 0 on some axis) uses
  brute force for every positive cutoff, but the grid for R = 0 and for any
  negative R.
* A negative cutoff never exceeds 0.3·S (S ≥ 0), so with n, m ≥ 10 it always
  goes to the grid search (→ error, §4); with n or m < 10 it goes to brute
  force and simply selects nothing.
* An invalid triclinic box (e.g. angles 0,0,0, or γ = 180) gives S = 0 →
  brute force for every positive cutoff.

### Size vectors for some boxes (float32), threshold = 0.3·S

| box [a, b, c, α, β, γ] | S | 0.3·S (brute force if R exceeds any) |
|---|---|---|
| [10, 10, 10, 90, 90, 90] | (10, 10, 10) | (3, 3, 3) |
| [100, 100, 100, 90, 90, 90] | (100, 100, 100) | (30.000001907…, …) |
| [10, 20, 30, 90, 90, 90] | (10, 20, 30) | (3, 6, 9) |
| [10, 10, 10, 60, 60, 60] | (5.0, 8.66025, 8.16497) | (1.5, 2.59808, 2.44949) |
| [10, 10, 10, 90, 90, 120] | (15.0, 8.66025, 10.0) | (4.5, 2.59808, 3.0) |
| [10, 10, 10, 90, 90, 30] | (10.0, 5.0, 10.0) | (3.0, 1.5, 3.0) |
| [10, 10, 10, 90, 90, 45] | (10.0, 7.07107, 10.0) | (3.0, 2.12132, 3.0) |
| [20, 20, 20, 30, 30, 30] | (2.67949, 10.0, 8.85782) | (0.80385, 3.0, 2.65735) |
| [10, 10, 10, 90, 90, 90.0001] (treated as triclinic) | (10.0000172, 10, 10) | (3.0000052, 3, 3) |
| [10, 10, 10, 0, 0, 0] or [10, 10, 10, 90, 90, 180] (invalid) | (0, 0, 0) | (0, 0, 0) |

## 4. Errors and limits that only the grid search has

When the grid search runs it validates its inputs; brute force has none of
these checks. All are ValueError, raised at evaluation time (not wrapped in
"Selection failed"):

| condition | message |
|---|---|
| cutoff < 0 | `Cutoff must be positive` |
| cutoff larger than the box allows (box present) | `Cutoff <R> too large for box (max <M>)` |
| a box length of exactly 0 in the (pseudo-)box | `Any of the box dimensions cannot be 0` |

* `<R>` is the cutoff printed as a Python float (`6.0`, `5.0001`, `1.3`);
  `<M>` is the limit computed in single precision and printed as a double
  (`5.0`, `1.25`, `3.535533905029297`).
* **Maximum allowed cutoff with a box**: with box vectors in lower-triangular
  form, let ax, by, cz be the diagonal entries (x of the first vector, y of
  the second, z of the third). The limit is half of the smallest of
  by·sin α, cz·sin α, ax·sin β, cz·sin β, ax·sin γ, by·sin γ. For an
  orthogonal box this is half the shortest edge. Because the decision rule
  sends R > 0.3·edge to brute force, an orthogonal box can only hit this
  error through rule 2 (n·m ≥ 1e8); a strongly skewed triclinic box can hit
  it with small systems (example below).
* **No box**: the grid search builds its own non-periodic cubic pseudo-box of
  edge max(largest extent, 2R) + 2.2R, so "too large" cannot occur; but a
  negative R gives `Cutoff must be positive`, or, if the pseudo-box edge comes
  out exactly 0 (e.g. extent 2.2 and R = −1), `Any of the box dimensions
  cannot be 0` (the zero-edge check comes first).
* **Cutoff exactly 0 without a box** (grid route, n, m ≥ 10): the pseudo-box
  edge then equals the extent, and atoms lying on the maximum face are
  reported at distance 0 from atoms on the minimum face, so spurious atoms
  are selected (see examples). Brute force with R = 0 selects only atoms
  coinciding exactly with a reference atom.
* Internal grid cell counts are capped at 1290 per axis; this has no
  user-visible effect on results.

## 5. Examples

`index a:b` selections are used to control n; the rest of the system forms m.

### 5.1 20 atoms on a line along x, spacing 1 Å, no box (extent (19, 0, 0))

| selection | n, m | method | result |
|---|---|---|---|
| `around 1.5 index 0:9` | 10, 10 | brute force (Sy = 0) | 1 atom |
| `around 5.7 index 0:9` | 10, 10 | brute force | 5 atoms |
| `around 1.5 index 0:8` | 9, 11 | brute force (n < 10) | 1 atom |
| `around 1.5 index 0:18` | 19, 1 | brute force (m < 10) | 1 atom |
| `sphzone 3 index 0:9` | 1, 20 | brute force | 6 atoms |
| `point 0 0 0 3` | 1, 20 | brute force | 4 atoms |
| `isolayer 1 3 index 0:9` | 10, 10 | brute force (both searches) | 2 atoms |
| `sphlayer 1 3 index 0:9` | 1, 20 | brute force | 4 atoms |
| `around 0 index 0:9` | 10, 10 | grid | 1 atom (index 19 — spurious, see §4) |
| `around -1 index 0:9` | 10, 10 | grid | ValueError `Cutoff must be positive` |

### 5.2 Same line, box [100, 100, 100, 90, 90, 90]

| selection | periodic on | periodic off |
|---|---|---|
| `around 1.5 index 0:9` | grid → 1 atom | brute force → 1 atom |
| `around 29.9 index 0:9` | grid → 10 atoms | brute force → 10 atoms |
| `around 30 index 0:9` | grid → 10 atoms | brute force → 10 atoms |
| `around 30.000002 index 0:9` | grid → 10 atoms (float32 tie) | brute force → 10 atoms |
| `around 30.0001 index 0:9` | brute force → 10 atoms | brute force → 10 atoms |
| `around -1 index 0:9` | grid → `Cutoff must be positive` | grid → `Cutoff must be positive` |

### 5.3 27 atoms on a 3×3×3 grid, spacing 2 Å (extent (4, 4, 4)), no box

| selection | n, m | method | result |
|---|---|---|---|
| `around 1.2 index 0:9` | 10, 17 | grid (1.2 ≤ float32(0.3·4)) | 0 atoms |
| `around 1.2000001 index 0:9` | 10, 17 | grid (same float32) | 0 atoms |
| `around 1.200001 index 0:9` | 10, 17 | brute force | 0 atoms |
| `around 1.21 index 0:9` | 10, 17 | brute force | 0 atoms |
| `around 1.2 index 0:8` | 9, 18 | brute force | 0 atoms |
| `around 1.2 index 0:16` | 17, 10 | grid | 0 atoms |
| `around 1.2 index 0:17` | 18, 9 | brute force | 0 atoms |
| `isolayer 1 3 index 0:9` | 10, 17 | outer R=3: brute force; inner R=1: grid | 11 atoms |
| `sphzone -1 index 0:9` | 1, 27 | brute force | 0 atoms |
| `around -1 index 0:8` | 9, 18 | brute force | 0 atoms |
| `around -1 index 0:9` | 10, 17 | grid | ValueError `Cutoff must be positive` |
| `isolayer -1 3 index 0:9` | 10, 17 | outer brute force, inner grid | ValueError `Cutoff must be positive` |
| `isolayer 1 -3 index 0:9` | 10, 17 | outer grid | ValueError `Cutoff must be positive` |
| `around 0 index 0:9` | 10, 17 | grid | 12 atoms (spurious, §4; brute force would give 0) |

With spacing 1.1 Å instead (extent 2.2): `around -1 index 0:9` → ValueError
`Any of the box dimensions cannot be 0`; `around -0.5 index 0:9`,
`around -3 index 0:9`, `around -1e-9 index 0:9` → `Cutoff must be positive`.

### 5.4 Same 27 atoms with various boxes (periodic on)

| box | `around 0.5` | `around 1.2` | `around 1.3` / `1.4` / `1.5` | `around 2.9` | `around 3` | `around 3.01` | `around 5` |
|---|---|---|---|---|---|---|---|
| [10,10,10,90,90,90] | grid, 0 | grid, 0 | — | grid, 11 | grid, 11 | brute, 11 | brute, 17 |
| [4,4,4,90,90,90] | grid, 12 | grid, 12 | — | brute, 17 | brute, 17 | brute, 17 | brute, 17 |
| [10,20,30,90,90,90] | grid, 0 | grid, 0 | — | grid, 11 | grid, 11 | brute, 11 | brute, 17 |
| [10,10,10,60,60,60] | grid, 0 | grid, 0 | — | brute, 11 | brute, 11 | brute, 11 | brute, 17 |
| [10,10,10,90,90,120] | grid, 0 | grid, 0 | — | brute, 11 | brute, 11 | brute, 11 | brute, 17 |
| [10,10,10,90,90,90.0001] | grid, 0 | grid, 0 | — | grid, 11 | grid, 11 | brute, 11 | brute, 17 |
| [10,10,10,90,90,30] | — | grid, 2 atoms | grid → `Cutoff 1.3 too large for box (max 1.25)` / `Cutoff 1.4 too large for box (max 1.25)` / `Cutoff 1.5 too large for box (max 1.25)` | — | brute, 12 | — | — |
| [10,10,10,90,90,45] | — | grid, 0 | grid, 0 / 0 / 0 | — | brute, 11 | — | — |
| [20,20,20,30,30,30] | — | brute, 0 | brute, 0 / 0 / 0 | — | brute, 11 | — | — |
| [10,10,10,0,0,0] (invalid) | brute, 0 | brute, 0 | brute | brute, 11 | brute, 11 | brute, 11 | brute, 17 |
| [0,0,0,90,90,90] (= no cell) | grid, 0 | grid, 0 | — | brute, 11 | brute, 11 | brute, 11 | brute, 17 |

(`around R index 0:9` in every cell; n = 10, m = 17. "—" = not run.
For [10,10,10,90,90,30] the size rule allows the grid up to R = 1.5 but the
grid's own limit is half of min(by·sin γ) = 2.5/2 = 1.25. `around 0` on an
invalid box goes to the grid and returns 0 atoms.)

### 5.5 Rule 2 (n·m ≥ 1e8): 20 000 random atoms in a [10,10,10,90,90,90] box

| selection | n, m | n·m | method | result |
|---|---|---|---|---|
| `around 6 index 0:9998` | 9 999, 10 001 | 99 999 999 | brute force (6 > 3) | 10 001 atoms |
| `around 6 index 0:9999` | 10 000, 10 000 | 1e8 | grid | ValueError `Cutoff 6.0 too large for box (max 5.0)` |
| `around 5 index 0:9999` | 10 000, 10 000 | 1e8 | grid | 10 000 atoms |
| `around 5.0001 index 0:9999` | 10 000, 10 000 | 1e8 | grid | ValueError `Cutoff 5.0001 too large for box (max 5.0)` |
| `around 3.5 index 0:9999` | 10 000, 10 000 | 1e8 | grid | 10 000 atoms |
| `around 6 index 0:9999`, periodic off | 10 000, 10 000 | 1e8 | grid (no box, pseudo-box) | 10 000 atoms |
| same atoms, box [10,10,10,60,60,60]: `around 4 index 0:9999` | | 1e8 | grid | ValueError `Cutoff 4.0 too large for box (max 3.535533905029297)` |
| same, `around 3 index 0:9999` | | 1e8 | grid | 10 000 atoms |

## 6. Related limits of the cylinder keywords (no method choice)

`cyzone R ZMAX ZMIN sel` and `cylayer R1 R2 ZMAX ZMIN sel` never call the
pair search. With `periodic` on **and** a unit cell present they refuse
cylinders that do not fit the cell, raising NotImplementedError at
evaluation (numbers formatted with three decimals; checks in this order:
x, y, then z):

* `The diameter of the cylinder selection (<2R>) is larger than the unit cell's x dimension (<a>). Can only do selections where it is smaller or equal.`
* `The diameter of the cylinder selection (<2R>) is larger than the unit cell's y dimension (<b>). Can only do selections where it is smaller or equal.`
* `The total length of the cylinder selection in z (<ZMAX − ZMIN>) is larger than the unit cell's z dimension (<c>). Can only do selections where it is smaller or equal.`

(R is the outer radius.) Examples on the fixture PDB with cell
[10, 8, 6, 90, 90, 90]:

| selection | periodic on | periodic off |
|---|---|---|
| `cyzone 3 2 -2 name N9` | 33 atoms (images folded in) | 5 atoms |
| `cyzone 4.5 2 -2 name N9` | `The diameter of the cylinder selection (9.000) is larger than the unit cell's y dimension (8.000). Can only do selections where it is smaller or equal.` | 7 atoms |
| `cyzone 5 2 -2 name N9` | `… (10.000) is larger than the unit cell's y dimension (8.000). …` (x passes: 10 ≤ 10) | 9 atoms |
| `cylayer 1 4.01 2 -2 name N9` | `The diameter of the cylinder selection (8.020) is larger than the unit cell's y dimension (8.000). Can only do selections where it is smaller or equal.` | 6 atoms |
| `cyzone 3.9 4 -3 name N9` | `The total length of the cylinder selection in z (7.000) is larger than the unit cell's z dimension (6.000). Can only do selections where it is smaller or equal.` | 7 atoms |
| `cyzone 3 3 -3 name N9` | 33 atoms (6 ≤ 6 allowed) | 5 atoms |
| `cyzone 3 3 -3.5 name N9` | `The total length of the cylinder selection in z (6.500) is larger than the unit cell's z dimension (6.000). Can only do selections where it is smaller or equal.` | 5 atoms |
