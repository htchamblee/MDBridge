# Box vectors → (a, b, c, α, β, γ) — behaviour specification (MDAnalysis 2.11.0)

Clean-room behaviour description of the conversion MDAnalysis applies when a
reader obtains three box (cell) vectors and turns them into the six unit-cell
"dimensions". Derived from running MDAnalysis 2.11.0 (NumPy 2.5.5, OpenBLAS
0.3.34, x86-64). "f32"/"f64" = IEEE binary32/binary64.

## 1. Where it is used

Readers that store a cell as three vectors use this conversion directly for
the per-frame box: GROMACS XTC, TRR, TPR, TNG, GRO (triclinic case),
DCD (vector-form unit cells), TRZ, DL_POLY CONFIG/HISTORY, DMS, H5MD,
FHI-aims, IMD, LAMMPS dump / data. (The AMBER readers do **not** use it:
NetCDF gives lengths+angles directly; mdcrd assumes 90° angles; inpcrd has
no box.)

The result is stored into the frame's f32 6-vector. Separately, the frame
reports "no box" (None) whenever the three stored lengths are all exactly
zero, so an invalid result (all zeros, §3) appears to the user as "no box".

## 2. The conversion

Inputs: three vectors **x**, **y**, **z** (the first, second, third box
vectors), each with three components, in any numeric type. Each component
is first converted exactly to f64 (f32 inputs are widened exactly).

All arithmetic below is in f64:

* lengths: a = |x|, b = |y|, c = |z|, where |v| = sqrt(v·v);
* angles in degrees:
  * α = angle between y and z = acos( (y·z) / (b·c) ) × (180/π)
  * β = angle between x and z = acos( (x·z) / (a·c) ) × (180/π)
  * γ = angle between x and y = acos( (x·y) / (a·b) ) × (180/π)

Numerical details that reproduce the observed bits:

* Dot product of 3-vectors u·v: on the reference machine it equals the
  fused chain fma(u3, v3, fma(u2, v2, u1·v1)) (OpenBLAS kernel). A plain
  left-to-right sum (u1v1 + u2v2) + u3v3 gave **identical f32 results** on
  600 000 random and GROMACS-style triples; the difference can only matter
  at f32 rounding ties.
* sqrt is the correctly rounded f64 square root.
* The cosine is one f64 division (dot) / (product of the two lengths, one
  f64 multiplication).
* acos is a faithful f64 arccos (NumPy's; differs from glibc `acos` by at
  most 1 ulp in ~9% of arguments, never changing the f32 result in 10^6
  tests). acos of an argument outside [−1, 1], or NaN, gives NaN (no error,
  no warning).
* Conversion to degrees is one multiplication by the f64 constant 180/π
  (0x1.ca5dc1a63c1f8p+5 = 57.29577951308232), not "×180 then ÷π".
* Division by zero (zero-length vector) yields NaN/±inf silently.

Then the six f64 values (a, b, c, α, β, γ) are each rounded to nearest f32.

## 3. Validity rule ("no box" / all zeros)

The result is returned **as is** iff **both**:

1. all six **f32-rounded** values are strictly greater than 0 (NaN fails;
   +inf passes); and
2. each of the three **f64** (unrounded) angles is strictly less than 180.

Otherwise the result is six f32 zeros.

Consequences (all verified, see table):

* Any zero-length vector → zero length → invalid. A length that is positive
  in f64 but rounds to 0 in f32 (|x| = 1e-46) → invalid; a length that
  rounds to an f32 subnormal (1e-40, 1e-45) → valid, kept.
* Parallel vectors: cos = 1 → angle 0 → invalid. Due to rounding the cos
  may come out slightly below 1 and then the angle is a tiny positive value
  → **valid** (e.g. x=(0.1,0.1,0.1), y=(0.3,0.3,0.3) gives γ = 8.5377366e-7°).
* Antiparallel: cos = −1 gives exactly 180 → invalid; cos < −1 → NaN →
  invalid.
* An angle that is < 180 in f64 but rounds to 180.0 in f32 is **valid** and
  returned as 180.0 (rule 2 uses the f64 angle, rule 1 the f32 value).
* NaN anywhere in an input that reaches a length or cosine → invalid.
  ±inf component → its length is inf and dot products involving it
  become NaN (inf×0) or inf/inf → NaN cosines → invalid.
* Lengths that overflow f32 (1e39) or whose square overflows f64 (1e200,
  length becomes +inf, cosines with an orthogonal vector are 0/inf = 0 →
  90°) → **valid**, returned with a +inf length.
* Negative components, left-handed sets, coplanar (zero volume) vector sets
  are all valid; only lengths and pairwise angles matter.
* Nothing checks the vectors' orientation/triangular form or the volume.

## 4. Example table

Inputs are written as Python floats (exact f64 values; float32 inputs shown
after widening). Results: f32 values (shortest repr) and the raw f32 bit
patterns, element order a b c α β γ.

| # | case | x | y | z | result (f32) | bits |
|---|---|---|---|---|---|---|
| 1 | orthorhombic | (10.0, 0.0, 0.0) | (0.0, 20.0, 0.0) | (0.0, 0.0, 30.0) | [10., 20., 30., 90., 90., 90.] | `41200000 41a00000 41f00000 42b40000 42b40000 42b40000` |
| 2 | unit cube | (1.0, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | [1., 1., 1., 90., 90., 90.] | `3f800000 3f800000 3f800000 42b40000 42b40000 42b40000` |
| 3 | f32 orthorhombic | (50.12300109863281, 0.0, 0.0) | (0.0, 60.45600128173828, 0.0) | (0.0, 0.0, 70.78900146484375) | [50.123, 60.456, 70.789, 90., 90., 90.] | `42487df4 4271d2f2 428d93f8 42b40000 42b40000 42b40000` |
| 4 | hexagonal | (10.0, 0.0, 0.0) | (-5.0, 8.660254037844386, 0.0) | (0.0, 0.0, 12.0) | [10., 10., 12., 90., 90., 120.] | `41200000 41200000 41400000 42b40000 42b40000 42f00000` |
| 5 | hexagonal, f32 inputs | (10.0, 0.0, 0.0) | (-5.0, 8.66025447845459, 0.0) | (0.0, 0.0, 12.0) | [10., 10., 12., 90., 90., 120.] | `41200000 41200000 41400000 42b40000 42b40000 42f00000` |
| 6 | rhombic dodecahedron (xy-square) | (7.0, 0.0, 0.0) | (0.0, 7.0, 0.0) | (3.5, 3.5, 4.949747468305833) | [7., 7., 7., 60., 60., 90.] | `40e00000 40e00000 40e00000 42700000 42700000 42b40000` |
| 7 | truncated octahedron | (7.0, 0.0, 0.0) | (2.3333333333333335, 6.599663257598877, 0.0) | (-2.3333333333333335, 3.2998316287994385, 5.715476036071777) | [7., 7., 7., 70.52878, 109.47122, 70.52878] | `40e00000 40e00000 40e00000 428d0ebc 42daf144 428d0ebc` |
| 8 | general triclinic | (30.1, 0.0, 0.0) | (5.2, 40.3, 0.0) | (-3.3, 7.7, 50.5) | [30.1, 40.634098, 51.190136, 81.898155, 93.69617, 82.64762] | `41f0cccd 42228951 424cc2b3 42a3cbdb 42bb6470 42a54b95` |
| 9 | arbitrary orientation | (3.0, 4.0, 0.0) | (0.0, 5.0, 12.0) | (8.0, 0.0, 15.0) | [5., 13., 17., 35.464058, 73.59931, 72.07979] | `40a00000 41500000 41880000 420ddb32 429332d9 429028da` |
| 10 | negative first vector | (-10.0, 0.0, 0.0) | (5.0, 10.0, 0.0) | (0.0, 0.0, 10.0) | [10., 11.18034, 10., 90., 90., 116.56505] | `41200000 4132e2ac 41200000 42b40000 42b40000 42e9214e` |
| 11 | all negative | (-10.0, 0.0, 0.0) | (0.0, -20.0, 0.0) | (0.0, 0.0, -30.0) | [10., 20., 30., 90., 90., 90.] | `41200000 41a00000 41f00000 42b40000 42b40000 42b40000` |
| 12 | 120° | (1.0, 0.0, 0.0) | (-0.5, 0.8660254037844386, 0.0) | (0.0, 0.0, 1.0) | [1., 1., 1., 90., 90., 120.] | `3f800000 3f800000 3f800000 42b40000 42b40000 42f00000` |
| 13 | 60° | (1.0, 0.0, 0.0) | (0.5, 0.8660254037844386, 0.0) | (0.0, 0.0, 1.0) | [1., 1., 1., 90., 90., 60.] | `3f800000 3f800000 3f800000 42b40000 42b40000 42700000` |
| 14 | z zero vector | (10.0, 0.0, 0.0) | (0.0, 20.0, 0.0) | (0.0, 0.0, 0.0) | all 0 | `00000000 ×6` |
| 15 | all zero | (0,0,0) | (0,0,0) | (0,0,0) | all 0 | `00000000 ×6` |
| 16 | x zero vector | (0.0, 0.0, 0.0) | (0.0, 20.0, 0.0) | (0.0, 0.0, 30.0) | all 0 | `00000000 ×6` |
| 17 | x ∥ y | (1.0, 0.0, 0.0) | (2.0, 0.0, 0.0) | (0.0, 0.0, 1.0) | all 0 | `00000000 ×6` |
| 18 | x antiparallel y | (1.0, 0.0, 0.0) | (-2.0, 0.0, 0.0) | (0.0, 0.0, 1.0) | all 0 | `00000000 ×6` |
| 19 | parallel, irrational lengths | (0.1, 0.1, 0.1) | (0.3, 0.3, 0.3) | (0.0, 0.0, 1.0) | [0.17320508, 0.51961523, 1., 54.73561, 54.73561, 8.5377366e-07] | `3e315cac 3f050581 3f800000 425af144 425af144 35652ee1` |
| 20 | nearly parallel (1e-7) | (1.0, 0.0, 0.0) | (1.0, 1e-07, 0.0) | (0.0, 0.0, 1.0) | [1., 1., 1., 90., 90., 5.663294e-06] | `3f800000 3f800000 3f800000 42b40000 42b40000 36be0753` |
| 21 | nearly parallel (1e-9) | (1.0, 0.0, 0.0) | (1.0, 1e-09, 0.0) | (0.0, 0.0, 1.0) | all 0 (cos rounds to 1 → angle 0) | `00000000 ×6` |
| 22 | γ ≈ 179.99999° | (1.0, 0.0, 0.0) | (-1.0, 1.7453292519943474e-07, 0.0) | (0.0, 0.0, 1.0) | [1., 1., 1., 90., 90., 179.99998] | `3f800000 3f800000 3f800000 42b40000 42b40000 4333ffff` |
| 23 | γ < 180 in f64, 180 in f32 | (1.0, 0.0, 0.0) | (-1.0, 1e-07, 0.0) | (0.0, 0.0, 1.0) | [1., 1., 1., 90., 90., 180.] (valid) | `3f800000 3f800000 3f800000 42b40000 42b40000 43340000` |
| 24 | nearly antiparallel (1e-9) | (1.0, 0.0, 0.0) | (-1.0, 1e-09, 0.0) | (0.0, 0.0, 1.0) | all 0 (exactly 180 in f64) | `00000000 ×6` |
| 25 | NaN component | (nan, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | all 0 | `00000000 ×6` |
| 26 | NaN off-diagonal | (1.0, 0.0, 0.0) | (nan, 1.0, 0.0) | (0.0, 0.0, 1.0) | all 0 | `00000000 ×6` |
| 27 | +inf component | (inf, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | all 0 | `00000000 ×6` |
| 28 | 1e39 (f32 overflow) | (1e+39, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | [inf, 1., 1., 90., 90., 90.] | `7f800000 3f800000 3f800000 42b40000 42b40000 42b40000` |
| 29 | 1e200 (square overflows f64) | (1e+200, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | [inf, 1., 1., 90., 90., 90.] | `7f800000 3f800000 3f800000 42b40000 42b40000 42b40000` |
| 30 | 1e-40 (f32 subnormal) | (1e-40, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | [1e-40, 1., 1., 90., 90., 90.] | `000116c2 3f800000 3f800000 42b40000 42b40000 42b40000` |
| 31 | 1e-46 (rounds to f32 0) | (1e-46, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | all 0 | `00000000 ×6` |
| 32 | 1e-200 (square underflows) | (1e-200, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | all 0 | `00000000 ×6` |
| 33 | 1e-45 (smallest f32 subnormal) | (1e-45, 0.0, 0.0) | (0.0, 1.0, 0.0) | (0.0, 0.0, 1.0) | [1e-45, 1., 1., 90., 90., 90.] | `00000001 3f800000 3f800000 42b40000 42b40000 42b40000` |
| 34 | coplanar (zero volume) | (1.0, 0.0, 0.0) | (0.0, 1.0, 0.0) | (1.0, 1.0, 0.0) | [1., 1., 1.4142135, 45., 45., 90.] | `3f800000 3f800000 3fb504f3 42340000 42340000 42b40000` |
| 35 | left-handed | (10.0, 0.0, 0.0) | (0.0, 0.0, 20.0) | (0.0, 30.0, 0.0) | [10., 20., 30., 90., 90., 90.] | `41200000 41a00000 41f00000 42b40000 42b40000 42b40000` |
| 36 | integers | (3, 0, 0) | (1, 3, 0) | (1, 1, 3) | [3., 3.1622777, 3.3166249, 67.5805, 72.4516, 71.56505] | `40400000 404a62c2 40544395 42872937 4290e738 428f214e` |

Row 7 inputs are the usual GROMACS truncated-octahedron vectors for edge 7
with the irrational components stored as f32 (7/3 kept in f64).

## 5. Cross-check performed

An independent re-implementation of §2–§3 (written only for verification,
not part of this spec) matched MDAnalysis bit-for-bit on 100 000 random
triples (components uniform in [−50, 50], 30 % pre-rounded to f32) using
the fused dot product, and on 300 000 further random triples plus 300 000
random lower-triangular f32 boxes using the plain left-to-right dot
product.
