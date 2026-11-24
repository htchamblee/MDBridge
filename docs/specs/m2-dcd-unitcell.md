# DCD unit cell: stored six numbers → [a, b, c, α, β, γ]

Clean-room behaviour specification derived by observing MDAnalysis 2.11.0
reading DCD files (CHARMM, NAMD and X-PLOR flavours; the separate
LAMMPS-flavoured DCD reader is out of scope). Test files were produced with
MDAnalysis' own low-level DCD writer (raw six numbers) and with hand-built
byte streams (endianness, missing blocks), then read back.

## 1. Where the six numbers come from

* The header's first record is 84 bytes and begins with the characters
  `CORD`. Counting byte offsets from the first byte of that 84-byte payload
  (i.e. `C` is offset 0):
  * the 32-bit integer at offset **80** (the "CHARMM version" slot) decides
    the flavour: non-zero → CHARMM-style file (CHARMM and NAMD both write
    non-zero here; MDAnalysis' writer puts 24); zero → X-PLOR-style file.
    Only zero / non-zero matters (any non-zero value, even negative, counts).
  * for CHARMM-style files only, the 32-bit integer at offset **44** decides
    whether frames carry a unit-cell record: non-zero → yes.
* An X-PLOR-style file, or a CHARMM-style file whose offset-44 flag is zero,
  never yields a unit cell: every frame's dimensions are "no unit cell"
  (MDAnalysis reports `None`).
* When present, each frame starts with a Fortran-style record: a 4-byte
  length, the payload, the same 4-byte length. If the length is 48, the
  payload is read as **six 64-bit IEEE doubles** in the file's byte order
  (big- and little-endian files both work; endianness is detected from the
  header). If the length is anything else the payload is skipped and that
  frame gets the default described below. (Frame counting assumes a 48-byte
  payload, so files with other lengths are not reliably readable; a one-frame
  file with a 24-byte record is rejected with
  `opened empty file. No frames are saved`.)
* Before each frame is read, the six-number buffer is reset to the default
  stored values `[0, 90, 0, 90, 90, 0]`, so a frame without a usable record
  never inherits the previous frame's cell. The default converts (by the
  rules below) to lengths 0,0,0 → "no unit cell".
* The decision of §3 is made **independently for every frame** (one file may
  mix conventions; see the multi-frame example in §6).

## 2. Meaning of the six stored slots

Call the stored numbers s1…s6 in file order. Their nominal meaning (X-PLOR /
NAMD / VMD ordering) is

| slot | s1 | s2 | s3 | s4 | s5 | s6 |
|---|---|---|---|---|---|---|
| nominal content | A (length) | γ term | B (length) | β term | α term | C (length) |

where the "angle terms" are either cosines of the angles or the angles in
degrees, and in the new-CHARMM convention all six slots are instead the
lower triangle of a symmetric box matrix (§3, case 3):

```
        | s1  s2  s4 |
  H  =  | s2  s3  s5 |
        | s4  s5  s6 |
```

The output order is always [a, b, c, α, β, γ]; reading s1,s3,s6 as a,b,c and
s5,s4,s2 as α,β,γ is the base re-ordering.

## 3. Decision rule (per frame)

Let a = s1, b = s3, c = s6, and the three angle terms tα = s5, tβ = s4,
tγ = s2. All tests are on the stored doubles.

1. **Cosine convention** — if each of tα, tβ, tγ lies in the closed interval
   [−1, +1] (a NaN fails this test), the angle terms are cosines. The lengths
   are taken as stored *without any sign or range check* (negative or NaN or
   infinite lengths pass through unchanged). Each angle is

   angle = 90 − arcsin(t) · 180/π   (degrees)

   evaluated in double precision. This equals arccos(t) mathematically but
   yields exactly 90 for t = 0 (also for −0), exactly 0 for t = 1, exactly
   180 for t = −1, and exactly 60/120 for ±0.5.
2. **Box-matrix (new CHARMM) convention** — otherwise, if **any** of the six
   values (lengths or angle terms) is negative, or any of tα, tβ, tγ is
   greater than 180, the six values are the matrix H above. With row vectors
   e1 = (s1, s2, s4), e2 = (s2, s3, s5), e3 = (s4, s5, s6):

   a = |e1|, b = |e2|, c = |e3|,
   α = arccos(e2·e3 / (b·c)), β = arccos(e1·e3 / (a·c)),
   γ = arccos(e1·e2 / (a·b)), in degrees,

   computed in double precision, then rounded to single precision. The result
   is then validated: if any of the six single-precision values is not > 0
   (NaN counts as failing), or any of the double-precision angles is not
   < 180, all six become 0 (→ "no unit cell"). (A double angle a hair below
   180 can still round to 180.0 in single precision and is kept.)
3. **Degree convention** (NAMD ≤ 2.5 style) — otherwise the values are used
   as they are: a, b, c as stored, α = tα, β = tβ, γ = tγ already in degrees.
   Note this branch is also taken for angle terms that are NaN, or slightly
   outside [−1,1] but non-negative (1.0000001 is "1.0000001 degrees"), and
   any *other* angle term that is e.g. 0 is then 0 **degrees**, not a cosine.

Note the ordering consequences: a value of exactly −1 or +1 is a cosine; a
negative length with valid cosines stays a negative length (cosine branch
wins), whereas a negative length with degree-style angles triggers the
box-matrix branch; an angle term of exactly 180 is degrees (not > 180).

## 4. Precision and final normalisation

* Inputs are 64-bit doubles; all conversions are done in double precision;
  the six results are stored in **single precision** (IEEE float32,
  round-to-nearest). Hence 1e40 → inf, 1e-30 → 1.0000000031710769e-30,
  12.3456789 → 12.34567928314209, a cosine of 0.25 → 75.52249145507812.
* After conversion, if a, b and c are **all exactly zero**, the frame has
  **no unit cell** (MDAnalysis reports `None`) regardless of the angles.
  Otherwise all six numbers are reported as computed (including zero or
  negative single lengths, zero angles, NaN, inf).
* Length unit is Å in the file and in the result (no scaling).

## 5. Cosine → angle examples

| stored cosine c | resulting angle (single precision) | 90 − asin(c)·180/π evaluated in double |
|---|---|---|
| 0.0 | 90.0 | 90.0 |
| -0.0 | 90.0 | 90.0 |
| 1e-08 | 90.0 | 89.9999994270422 |
| -1e-08 | 90.0 | 90.0000005729578 |
| 0.1 | 84.26082611083984 | 84.26082952273322 |
| 0.25 | 75.52249145507812 | 75.52248781407008 |
| 0.5 | 60.0 | 60.0 |
| 0.49999999999999994 | 60.0 | 60.0 |
| 0.7071067811865476 | 45.0 | 44.99999999999999 |
| 0.8660254037844386 | 30.0 | 30.000000000000007 |
| -0.3333333333333333 | 109.47122192382812 | 109.47122063449069 |
| -0.5 | 120.0 | 120.0 |
| 0.9999999 | 0.02562345191836357 | 0.025623451769789085 |
| 0.99999999999 | 0.0002562345180194825 | 0.00025623452621914566 |
| 1.0 | 0.0 | 0.0 |
| -1.0 | 180.0 | 180.0 |

## 6. Stored six numbers → dimensions (CHARMM-style file with cell records)

"branch" names the rule of §3 that applies.

| stored s1..s6 = [A, γ-slot, B, β-slot, α-slot, C] | branch | dimensions read back [a, b, c, α, β, γ] |
|---|---|---|
| [10, 0, 20, 0, 0, 30] | cosine | [10.0, 20.0, 30.0, 90.0, 90.0, 90.0] |
| [10, 0.5, 20, 0, 0, 30] | cosine | [10.0, 20.0, 30.0, 90.0, 90.0, 60.0] |
| [10, 0, 20, 0.5, 0, 30] | cosine | [10.0, 20.0, 30.0, 90.0, 60.0, 90.0] |
| [10, 0, 20, 0, 0.5, 30] | cosine | [10.0, 20.0, 30.0, 60.0, 90.0, 90.0] |
| [10, -0.5, 20, 0.2, -0.3, 30] | cosine | [10.0, 20.0, 30.0, 107.45760345458984, 78.46304321289062, 120.0] |
| [10, 1, 20, 1, 1, 30] | cosine | [10.0, 20.0, 30.0, 0.0, 0.0, 0.0] |
| [10, -1, 20, -1, -1, 30] | cosine | [10.0, 20.0, 30.0, 180.0, 180.0, 180.0] |
| [10, 1.0000001, 20, 0, 0, 30] | degrees | [10.0, 20.0, 30.0, 0.0, 0.0, 1.0000001192092896] |
| [10, 90, 20, 90, 90, 30] | degrees | [10.0, 20.0, 30.0, 90.0, 90.0, 90.0] |
| [10, 120, 20, 90, 90, 30] | degrees | [10.0, 20.0, 30.0, 90.0, 90.0, 120.0] |
| [10, 60, 20, 70, 80, 30] | degrees | [10.0, 20.0, 30.0, 80.0, 70.0, 60.0] |
| [10, 0.9, 20, 90, 90, 30] | degrees | [10.0, 20.0, 30.0, 90.0, 90.0, 0.8999999761581421] |
| [10, 180, 20, 90, 90, 30] | degrees | [10.0, 20.0, 30.0, 90.0, 90.0, 180.0] |
| [10, 180.0001, 20, 90, 90, 30] | box matrix | [201.4945068359375, 202.23757934570312, 130.7669677734375, 38.48900604248047, 41.28325653076172, 70.65288543701172] |
| [10, 200, 20, 90, 90, 30] | box matrix | [219.54498291015625, 220.22715759277344, 130.7669677734375, 38.620941162109375, 41.203800201416016, 73.04473114013672] |
| [10, 0, 20, 0, 0, -30] | cosine | [10.0, 20.0, -30.0, 90.0, 90.0, 90.0] |
| [-10, 0, 20, 0, 0, 30] | cosine | [-10.0, 20.0, 30.0, 90.0, 90.0, 90.0] |
| [-10, 90, 20, 90, 90, 30] | box matrix | [127.67145538330078, 128.8409881591797, 130.7669677734375, 41.594932556152344, 53.630916595458984, 56.82930374145508] |
| [10, 90, -20, 90, 90, 30] | box matrix | [127.67145538330078, 128.8409881591797, 130.7669677734375, 57.711448669433594, 45.5087890625, 64.0422592163086] |
| [10, 0, 20, 0, 0, 0] | cosine | [10.0, 20.0, 0.0, 90.0, 90.0, 90.0] |
| [0, 0, 0, 0, 0, 0] | cosine | None |
| [0, 90, 0, 90, 90, 0] | degrees | None |
| [0, 0.5, 0, 0.5, 0.5, 0] | cosine | None |
| [10, 0, 0, 0, 0, 0] | cosine | [10.0, 0.0, 0.0, 90.0, 90.0, 90.0] |
| [10, -1.5, 20, 90, 90, 30] | box matrix | [90.56627655029297, 92.20764923095703, 130.7669677734375, 68.77638244628906, 72.98776245117188, 15.298059463500977] |
| [10, -5, 20, 0, 0, 30] | box matrix | [11.180339813232422, 20.615528106689453, 30.0, 90.0, 90.0, 130.60128784179688] |
| [10, 2, 20, -3, 1, 30] | box matrix | [10.630146026611328, 20.124610900878906, 30.16620635986328, 85.84368896484375, 111.59098052978516, 74.54714965820312] |
| [30, -5, 30, 2, -4, 30] | box matrix | [30.479501724243164, 30.675724029541016, 30.331501007080078, 105.58631896972656, 81.28990173339844, 109.23348999023438] |
| [10, 0, 10, 0, 0, 10] | cosine | [10.0, 10.0, 10.0, 90.0, 90.0, 90.0] |
| [10, 5, 10, 0, 0, 10] | degrees | [10.0, 10.0, 10.0, 0.0, 0.0, 5.0] |
| [10, -5, 10, 0, 0, 10] | box matrix | [11.180339813232422, 11.180339813232422, 10.0, 90.0, 90.0, 143.13009643554688] |
| [10, -10, 10, 0, 0, 10] | box matrix | [14.142135620117188, 14.142135620117188, 10.0, 90.0, 90.0, 180.0] |
| [10, -11, 10, 0, 0, 10] | box matrix | [14.866068840026855, 14.866068840026855, 10.0, 90.0, 90.0, 174.54737854003906] |
| [10, 0, -10, 0, 0, 10] | cosine | [10.0, -10.0, 10.0, 90.0, 90.0, 90.0] |
| [10, 0, 20, 0, 0, nan] | cosine | [10.0, 20.0, nan, 90.0, 90.0, 90.0] |
| [10, nan, 20, 0, 0, 30] | degrees | [10.0, 20.0, 30.0, 0.0, 0.0, nan] |
| [10, inf, 20, 90, 90, 30] | box matrix | None |
| [10, 0, 20, 0, 0, inf] | cosine | [10.0, 20.0, inf, 90.0, 90.0, 90.0] |
| [10, 0.5, 20, 90, 90, 30] | degrees | [10.0, 20.0, 30.0, 90.0, 90.0, 0.5] |
| [10, -0.0, 20, -0.0, -0.0, 30] | cosine | [10.0, 20.0, 30.0, 90.0, 90.0, 90.0] |
| [10, 0.5000000000000001, 20, 0.3420201433256688, 0.17364817766693041, 30] | cosine | [10.0, 20.0, 30.0, 80.0, 70.0, 60.0] |
| [1e-30, 0, 1e-30, 0, 0, 1e-30] | cosine | [1.0000000031710769e-30, 1.0000000031710769e-30, 1.0000000031710769e-30, 90.0, 90.0, 90.0] |
| [1e+40, 0, 20, 0, 0, 30] | cosine | [inf, 20.0, 30.0, 90.0, 90.0, 90.0] |
| [10, -1e-09, 20, 90, 90, 30] | box matrix | [90.55384826660156, 92.19544219970703, 130.7669677734375, 68.08353424072266, 72.30097961425781, 14.018765449523926] |
| [10, 0, 20, 0, 0, -1e-09] | cosine | [10.0, 20.0, -9.999999717180685e-10, 90.0, 90.0, 90.0] |
| [-10, -10, -10, 0, 0, 10] | box matrix | [14.142135620117188, 14.142135620117188, 10.0, 90.0, 90.0, 8.537736562175269e-07] |
| [0, 0, 0, 0, 0, -5] | cosine | [0.0, 0.0, -5.0, 90.0, 90.0, 90.0] |
| [0, 0, 20, 0, 0, 0] | cosine | [0.0, 20.0, 0.0, 90.0, 90.0, 90.0] |
| [10, 190, 20, 90, 90, 30] | box matrix | [210.47564697265625, 211.18711853027344, 130.7669677734375, 38.542442321777344, 41.2281608581543, 71.9128646850586] |
| [10, 90, 20, 190, 90, 30] | box matrix | [210.47564697265625, 128.8409881591797, 212.3675994873047, 37.86809539794922, 69.4365463256836, 43.10142517089844] |
| [10, 90, 20, 90, 181, 30] | box matrix | [127.67145538330078, 203.12803649902344, 204.35508728027344, 65.5970687866211, 40.3280143737793, 42.92437744140625] |
| [10, 90, 20, 90, 90, 1000] | degrees | [10.0, 20.0, 1000.0, 90.0, 90.0, 90.0] |
| [10, 1.0000000000000002, 20, 0, 0, 30] | degrees | [10.0, 20.0, 30.0, 0.0, 0.0, 1.0] |
| [10, -1.0000000000000002, 20, 0, 0, 30] | box matrix | [10.049875259399414, 20.02498435974121, 30.0, 90.0, 90.0, 98.572998046875] |

Other files:

| file | dimensions of each frame |
|---|---|
| X-PLOR style (offset-80 integer = 0), any frame content | None |
| CHARMM style, offset-44 flag = 0 (no cell record) | None |
| big-endian CHARMM file, stored [10, 0.5, 20, 0, 0, 30] | [10.0, 20.0, 30.0, 90.0, 90.0, 60.0] |
| 4 frames stored [10,0,20,0,0,30], [10,90,20,90,90,30], [10,−5,10,0,0,10], [0,0,0,0,0,0] | [10,20,30,90,90,90]; [10,20,30,90,90,90]; [11.180339813232422, 11.180339813232422, 10.0, 90.0, 90.0, 143.13009643554688]; None |

## 7. For reference: what MDAnalysis' writer stores

The writer always produces a CHARMM-style file (version 24) with cell records
in the cosine convention: stored = [a, cos γ, b, cos β, cos α, c], where each
cosine is computed as sin((90 − angle)·π/180) from the single-precision
dimensions. If the frame has no unit cell it warns
`No dimensions set for current frame, zeroed unitcell will be written` and
writes all-zero dimensions, which become stored [0, 1, 0, 1, 1, 0] and read
back as no unit cell.

| dimensions written [a, b, c, α, β, γ] | six stored doubles | dimensions read back | warning |
|---|---|---|---|
| [10, 20, 30, 90, 90, 90] | [10.0, 0.0, 20.0, 0.0, 0.0, 30.0] | [10.0, 20.0, 30.0, 90.0, 90.0, 90.0] |  |
| [10, 20, 30, 60, 70, 80] | [10.0, 0.1736481785774231, 20.0, 0.3420201241970062, 0.5, 30.0] | [10.0, 20.0, 30.0, 60.0, 70.0, 80.0] |  |
| [10, 20, 30, 80, 70, 60] | [10.0, 0.5, 20.0, 0.3420201241970062, 0.1736481785774231, 30.0] | [10.0, 20.0, 30.0, 80.0, 70.0, 60.0] |  |
| [50, 50, 50, 109.4712206, 109.4712206, 109.4712206] | [50.0, -0.3333333432674408, 50.0, -0.3333333432674408, -0.3333333432674408, 50.0] | [50.0, 50.0, 50.0, 109.47122192382812, 109.47122192382812, 109.47122192382812] |  |
| [10, 10, 10, 60, 60, 90] | [10.0, 0.0, 10.0, 0.5, 0.5, 10.0] | [10.0, 10.0, 10.0, 60.0, 60.0, 90.0] |  |
| [10, 10, 10, 90, 90, 120] | [10.0, -0.5, 10.0, 0.0, 0.0, 10.0] | [10.0, 10.0, 10.0, 90.0, 90.0, 120.0] |  |
| [10, 20, 30, 0, 0, 0] | [10.0, 1.0, 20.0, 1.0, 1.0, 30.0] | [10.0, 20.0, 30.0, 0.0, 0.0, 0.0] |  |
| [10, 20, 30, 180, 180, 180] | [10.0, -1.0, 20.0, -1.0, -1.0, 30.0] | [10.0, 20.0, 30.0, 180.0, 180.0, 180.0] |  |
| [10, 20, 30, 90.000001, 90, 90] | [10.0, 0.0, 20.0, 0.0, 0.0, 30.0] | [10.0, 20.0, 30.0, 90.0, 90.0, 90.0] |  |
| [12.3456789, 20, 30, 90, 90, 90] | [12.34567928314209, 0.0, 20.0, 0.0, 0.0, 30.0] | [12.34567928314209, 20.0, 30.0, 90.0, 90.0, 90.0] |  |
| None | [0.0, 1.0, 0.0, 1.0, 1.0, 0.0] | None | `No dimensions set for current frame, zeroed unitcell will be written` |
