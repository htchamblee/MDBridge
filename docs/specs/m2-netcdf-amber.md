# AMBER NetCDF trajectory reading — behaviour specification (MDAnalysis 2.11.0)

Clean-room behaviour description. Everything below was observed by running
MDAnalysis 2.11.0 (NumPy 2.5.5, SciPy 1.18.1, no `netCDF4` module installed)
on crafted files. It describes *what* happens, not how MDAnalysis does it.

Notation: `<path>` stands for the file name **exactly as the caller supplied
it** (a relative path stays relative; it is not made absolute). "f32" means
IEEE binary32, "f64" binary64. "Warning" means a Python `UserWarning`
(it is also logged); a warning never stops reading.

---

## 1. Which files reach this reader

The reader is selected by file extension, case-insensitively: `nc`, `ncdf`
(also `x.nc.gz` is routed here, but then fails as "not a valid NetCDF 3
file" because compressed input is not decompressed). `netcdf`, `rst7`,
`ncrst` are **not** recognised extensions (unknown-format error from the
format lookup, not from this reader).

## 2. Container format support

The file is parsed as a netCDF *classic* file (public netCDF-3 spec). Only
the behaviour that matters for the reader is listed.

| file content | outcome |
|---|---|
| does not start with the 3 bytes `CDF` (e.g. HDF5/NetCDF-4 files that start with `\x89HDF\r\n\x1a\n`, plain text) | `TypeError: Error: <path> is not a valid NetCDF 3 file` |
| empty file (0 bytes) | `ValueError: cannot mmap an empty file` |
| file does not exist | `FileNotFoundError: [Errno 2] No such file or directory: '<path>'` |
| `CDF` + version byte 1 (classic, 32-bit offsets) | header parses fine, then rejected by the AMBER checks (step 4 below) |
| `CDF` + version byte 2 (64-bit offset) | accepted |
| `CDF` + version byte 5 (CDF-5) or 3 | `IndexError: list index out of range` (raised while parsing the header, i.e. before any AMBER check) |
| `CDF` + version byte 0 | header parsed as if 64-bit offsets; later rejected at step 4 |
| truncated inside the header | a low-level error, e.g. `IndexError: index 0 is out of bounds for axis 0 with size 0` |
| record data shorter than `numrecs × recsize` (truncated data section) | `ValueError: When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.` |
| header `numrecs` larger than the records present | e.g. `ValueError: cannot reshape array of size 2 into shape (3,)` |
| header `numrecs` smaller than the records present | only `numrecs` records exist (extra data ignored) |
| header `numrecs` = 0xFFFFFFFF (streaming) | treated as no usable frames: `IndexError: frame index must be 0 <= frame < 0` |

Header-parsing facts the implementer needs:

* A dimension whose stored length is 0 is the record (unlimited) dimension;
  its length is "unknown" (see frame count, §6). Any dimension stored with
  length 0 is treated this way — so a `spatial` dimension of length 0 would
  compare unequal to 3 (step 5 fails with the `TypeError`).
* Names (dimension, variable, attribute) are decoded as Latin-1.
* Attribute values: a character (NC_CHAR) attribute becomes a byte string
  with **trailing NUL bytes removed** (other trailing bytes such as spaces
  are kept). A numeric attribute with exactly one element becomes a scalar
  of that numeric type (NC_FLOAT → f32 scalar, NC_DOUBLE → f64 scalar,
  NC_SHORT → int16, NC_INT → int32, …). A numeric attribute with several
  elements stays an array.
* Text checks below decode attribute bytes as UTF-8; invalid UTF-8 would
  raise a `UnicodeDecodeError` (not otherwise handled).
* Variable data is big-endian, typed NC_BYTE/SHORT/INT/FLOAT/DOUBLE (int8,
  int16, int32, f32, f64). Record variables are interleaved per record as in
  the classic format.
* `_FillValue`, `missing_value`, `add_offset`, `valid_*` attributes are
  **ignored** (no masking, no offset). E.g. coordinates equal to the
  default fill 9.96921e36 are returned as 9.969209968386869e+36.
* The memory-mapping option (`mmap=True/False/None`) changes only how bytes
  are fetched; observed values and errors are identical (except the empty
  file message above, which is the mmap one; default is mmap for a path).

## 3. Open-time checks, in the exact order performed

The first failing check raises; warnings from earlier checks have already
been emitted. After all checks pass, frame 0 is read (§5) — so any per-frame
error for frame 0 is also an open-time error.

1. **Parse the container** (§2 errors).
2. **`Conventions` global attribute.**
   * Absent, or present but not a text attribute (e.g. an integer) →
     `ValueError: NCDF trajectory <path> is missing Conventions`
   * Text value accepted iff the token `AMBER` (exact, case-sensitive)
     appears either among the pieces obtained by splitting the value at
     every comma (pieces are *not* trimmed), or among the pieces obtained by
     splitting at runs of whitespace (leading/trailing whitespace ignored).
   * Otherwise →
     `TypeError: NCDF trajectory <path> does not conform to AMBER specifications, http://ambermd.org/netcdf/nctraj.xhtml ('AMBER' must be one of the token in attribute Conventions)`
3. **`ConventionVersion` global attribute.**
   * Absent or non-text → `ValueError: NCDF trajectory <path> is missing ConventionVersion`
   * Text but not exactly `1.0` (byte-exact after NUL stripping; no
     trimming) → Warning
     `NCDF trajectory format is <value> but the reader implements format 1.0`
     where `<value>` is the decoded text verbatim (spaces preserved).
4. **Format version byte must be 2** (64-bit offset). Otherwise →
   `TypeError: NCDF trajectory <path> does not conform to AMBER specifications, as detailed in https://ambermd.org/netcdf/nctraj.xhtml (NetCDF file does not use 64 bit offsets [version_byte = 2])`
   (Note the different URL scheme/wording from step 2.)
5. **`spatial` dimension.**
   * Absent → `ValueError: NCDF trajectory does not contain spatial dimension` (no path in this message)
   * Length ≠ 3 (including unlimited) → `TypeError: Incorrect spatial value for NCDF trajectory file`
6. **`program` / `programVersion` global attributes.** If either is absent
   (any type counts as present) → Warning
   ``NCDF trajectory <path> may not fully adhere to AMBER standards as either the `program` or `programVersion` attributes are missing``
7. **`atom` dimension → atom count.**
   * Absent → `ValueError: NCDF trajectory <path> does not contain atom information` (raised even if the caller supplied an atom count).
   * If the caller supplied an atom count (when loaded into a Universe, the
     topology's atom count is always supplied) and it differs →
     `ValueError: Supplied n_atoms (<given>) != natom from ncdf (<file>). Note: n_atoms can be None and then the ncdf value is used!`
     This is the message a user sees for a topology/trajectory atom-count
     mismatch (the Universe-level mismatch message is never reached).
8. **`frame` dimension → frame count.**
   * Absent → `ValueError: NCDF trajectory <path> does not contain frame information`
   * Fixed length → that length is the frame count.
   * Unlimited → frame count = leading extent of the `coordinates` variable
     (= the record count). If `coordinates` does not exist in this case, the
     same "does not contain frame information" `ValueError` is raised.
9. **`title` global attribute** → remarks (raw bytes, unmodified); absent →
   empty string. Never an error.
10. **`time` variable.**
    * Absent → Warning
      ``NCDF trajectory does not contain `time` information; `time` will be set as an increasing index``
      and the file is treated as having no time (§7).
    * Present: its `units` attribute must equal exactly `picosecond`.
      Missing `units` → `AttributeError: 'netcdf_variable' object has no attribute 'units'`;
      numeric `units` → `AttributeError: 'numpy.int32' object has no attribute 'decode'` (type name varies with the attribute type);
      any other text →
      `NotImplementedError: NETCDFReader currently assumes that the trajectory was written in units of <found> instead of picosecond`
11. **`coordinates` variable `units`** must be exactly `angstrom` (same three
    failure forms, `<expected>` = `angstrom`). If the `coordinates` variable
    does not exist at this point (only possible with a fixed-length `frame`
    dimension) → `KeyError: 'coordinates'`.
12. **`scale_factor` attributes**, visiting variables in file order:
    * Variable has no `scale_factor` → nothing.
    * Variable is one of `time`, `cell_lengths`, `cell_angles`,
      `coordinates`, `velocities`, `forces`: the attribute must be a
      floating scalar (NC_FLOAT or NC_DOUBLE, single element). Otherwise
      `TypeError: <repr> is not a float`, where `<repr>` is the value as
      printed: `2` for an int attribute, `[1. 2.]` for a 2-element double
      array, `b'2.0'` for a text attribute.
    * Any other variable carrying `scale_factor` →
      `NotImplementedError: scale_factors for variable <name> are not implemented`
    * Because the visit follows file order, whichever offending variable
      comes first decides which of these two errors is raised.
13. **`velocities` variable** (optional): if present, `units` must be
    exactly `angstrom/picosecond` (same failure forms).
14. **`forces` variable** (optional): if present, `units` must be exactly
    `kilocalorie/mole/angstrom`.
15. **Periodic box**: the file is periodic iff a `cell_lengths` variable
    exists. If so, `cell_lengths` units must be `angstrom`, then
    `cell_angles` units must be `degree` (singular; `degrees`, `Degree`,
    `radian` are rejected with the `NotImplementedError` above). If
    `cell_lengths` exists but `cell_angles` does not → `KeyError: 'cell_angles'`.
    A `cell_angles` variable without `cell_lengths` is ignored (no box).
16. **Read frame 0** (§5).

Units comparisons are byte-exact (no trimming; `picosecond ` with a trailing
space is rejected and the message then shows the space). The units
attribute is only examined for these variables; other variables (e.g.
`cell_spatial` label variables, `foo`) are not checked.

### Ordering examples (observed)

| file defects | error raised |
|---|---|
| Conventions `X` and no `spatial` | the Conventions `TypeError` |
| no `spatial`, no `atom`, atom count mismatch | `ValueError: NCDF trajectory does not contain spatial dimension` |
| atom count mismatch and no `frame` dimension | the `Supplied n_atoms (7) != natom from ncdf (3)…` `ValueError` |
| missing `program` and atom-count mismatch | Warning about program, then the n_atoms `ValueError` |
| ConventionVersion `2` and CDF-1 | Warning `…format is 2…`, then the 64-bit-offset `TypeError` |
| time units `ns` and coordinates units `nm` | the `ns` error (time checked first) |
| time units `ns` and an int `scale_factor` on time | the `ns` error (units before scale factors) |
| coordinates units `nm` and int `scale_factor` on coordinates | the `nm` error |
| int `scale_factor` on coordinates and bad velocity units | `TypeError: 2 is not a float` |
| velocities units `angstrom/ps` and forces units `kcal/mol/A` | the velocities error |

### Conventions examples

| `Conventions` value (bytes) | result |
|---|---|
| `AMBER` | OK |
| `AMBER,CF-1.0` | OK (comma split) |
| `CF-1.0,AMBER` | OK |
| `CF AMBER` | OK (whitespace split) |
| `CF, AMBER` | OK (whitespace split gives `CF,` and `AMBER`) |
| ` AMBER` | OK (whitespace split) |
| `AMBER\tCF` | OK |
| `AMBER\0\0` | OK (NULs stripped) |
| `CF,AMBER ` | TypeError (comma piece is `AMBER ` with a space; whitespace split gives one token `CF,AMBER`) |
| `amber`, `AMBERX`, `CF`, empty | TypeError |
| int32 attribute `3` | `ValueError … is missing Conventions` |
| absent | `ValueError … is missing Conventions` |

### ConventionVersion examples

| value | result |
|---|---|
| `1.0` | silent |
| `1.00`, `2.0`, ` 1.0`, `1.0 `, empty | Warning, e.g. `NCDF trajectory format is 1.0  but the reader implements format 1.0` (two spaces when value is `1.0 `) |
| float attribute 1.0, or absent | `ValueError … is missing ConventionVersion` |

## 4. Atom count, frame count, remarks

* Atom count = `atom` dimension length.
* Frame count per step 8. A file with zero records (unlimited frame,
  numrecs 0) fails at open while reading frame 0:
  `IndexError: frame index must be 0 <= frame < 0`.
* The reader's standalone "how many atoms" query (used for formats where
  topology needs it) opens the file and returns the `atom` dimension; no
  AMBER checks are made by that query (a missing `atom` dimension gives
  `KeyError: 'atom'`).

## 5. Per-frame read

Random access by 0-based integer frame index. Public indexing converts
negative indices (−1 = last) and range-checks first:
`IndexError: Index <i> exceeds length of trajectory (<n>).` where `<i>` is
the index *after* adding n for negatives (e.g. index −3 on 2 frames reports
`Index -1 exceeds…`). A float index →
`TypeError: Trajectories must be an indexed using an integer, slice or list of indices`.
A boolean index (`True`) passes the public check but then gives
`TypeError: frame must be a positive integer or zero`. (The internal frame
read rejects any index that is not the platform default integer type with
that same message, and out-of-range with
`IndexError: frame index must be 0 <= frame < <n>`.) After the file is
closed: `OSError: Trajectory is closed`.

For frame k the following are produced (all stored arrays are f32):

1. **Positions**: `coordinates[k]` (shape atom × 3), multiplied by its scale
   factor if one applies (see §5.1), then stored as f32.
2. **Time** (only when a `time` variable exists): `time[k]` × scale factor
   if applicable. The value is kept **as a scalar of the resulting numeric
   type**, not converted to f32/f64: an f32 `time` variable gives an f32
   time (`10.0` as f32), an f64 variable gives f64, an int32 variable gives
   int32; an f32 variable with an f64 scale factor gives f64. A user
   `time_offset` is added on access (f32 10.0 + offset 100 → f32 110.0).
3. **Velocities** (if the variable exists): `velocities[k]` × scale factor,
   stored as f32. Native unit Å/ps equals the library unit, no conversion.
4. **Forces** (if the variable exists): `forces[k]` × scale factor, stored
   as f32; then, if unit conversion is on (default), multiplied **in f32
   arithmetic** by f32(4.184) = 0x4085E354 (4.184000015258789), i.e. kcal →
   kJ, result rounded to f32. Example: stored force 112.908646 (f32) →
   472.4098 (f32); a f64 product rounded to f32 would give 472.40976, which
   is **not** what MDAnalysis returns. With `convert_units=False` forces stay
   112.908646.
5. **Box** (if periodic): a 6-vector `[cell_lengths[k]×sf, cell_angles[k]×sf]`
   built in f64 then stored as f32. It is taken verbatim — no validation,
   no conversion of angles, negative/NaN/out-of-range values are kept
   (lengths `[-1,40,50]` with angles `[0,-5,400]` → `[-1,40,50,0,-5,400]`;
   NaN kept). If the `cell_*` variables are not per-frame but 1-D over
   `cell_spatial`, indexing element k of them is used for all three entries
   (observed: lengths `[10,20,30]` give box `[10,10,10,…]` at frame 0 and
   `[20,20,20,…]` at frame 1) — an artefact, not a feature.
   **Zero-box defect:** the library reports "no box" whenever the three
   stored lengths are all exactly zero (−0.0 counts as zero). In that case,
   with unit conversion on (the default), reading that frame raises
   `TypeError: 'NoneType' object is not subscriptable` (at open time if it
   is frame 0; later only when that frame is visited). With
   `convert_units=False` the frame reads fine and the box is reported as
   absent. Partially-zero lengths (e.g. `[0,0,5]`) are kept as is.
6. Non-periodic files report no box (all-zero cell internally, "None" to the user).
7. Frame label = k.

### 5.1 Scale factors and dtype consequences

* A scale factor applies to a variable iff the variable has a floating
  `scale_factor` and it is **not** close to 1, where "close" means
  |s − 1| ≤ 1e-9 · max(|s|, 1) (relative tolerance 1e-9, no absolute
  tolerance). If close, the raw values are used unmultiplied, in their raw
  type.
* The multiplication happens in the NumPy-promoted type of (variable
  element type, scale-factor scalar type):
  f32 × f32 → f32; f32 × f64 → f64; f64 × any float → f64;
  int16/int8 × f32 → f32; int32 × f32 → f64; int × f64 → f64.
  The product is then rounded to f32 for positions/velocities/forces, to
  f64 then f32 for the box, and kept as-is for time. So f32 data with an f64
  factor is computed in double and rounded once to f32; f32 data with an f32
  factor is rounded after an f32 multiply — these can differ in the last
  bit.
* Raw values with no scale factor are converted directly to f32 (f64
  coordinates `0.1+2^-30` → f32 0.1).

| variable / data type | scale_factor (type) | stored raw | observed value |
|---|---|---|---|
| coordinates f32 | 2.0 (f32) | 0.12345679 | 0.2469135820865631 |
| coordinates f32 | 2.0 (f64) | 0.12345679 | 0.2469135820865631 |
| coordinates f32 | 1.0 (f32) | — | unscaled |
| coordinates f32 | 1+1e-10 (f64) | — | unscaled (close to 1) |
| coordinates f32 | 0.1 (f32) | f32(1.2234567) | 0.12234567850828171 (f32 multiply) |
| coordinates f32 | 0.1 (f64) | f32(1.2234567) | 0.12234567105770111 (f64 multiply, one rounding) |
| coordinates f32 | 0.1 (f64) | f32(2.3234568) | 0.23234567046165466 (vs 0.23234568536281586 with f32 factor) |
| coordinates f64 | 0.1 (f32) | 1.2234567890 | 0.12234567850828171 |
| coordinates int16 | none | 1, 2 | 1.0, 2.0 |
| coordinates int16 | 0.01 (f32) | 1, 2 | 0.009999999776482582, 0.019999999552965164 |
| time f32 | 2.0 (f32) | 10.0 | time 20.0 (f32) |
| time f32 | 2.0 (f64) | 10.0 | time 20.0 (f64) |
| time f64 | 1+5e-10 (f64) | 10.0 | 10.0 (not applied) |
| time f64 | 1+2e-9 (f64) | 10.0 | 10.00000002 (applied) |
| time f64 | 1−1.1e-9 (f64) | 10.0 | 9.999999989 |
| cell_lengths f64 | 0.5 (f64) | 30,40,50 | box lengths 15,20,25 |
| cell_angles f64 | 0.5 (f32) | 90,100,110 | angles 45,50,55 |
| velocities f32 | 20.455 (f64) | 1.5 | 30.682498931884766 |
| forces f32 | 2.0 (f64) | 112.908646 | 944.8196 (scaled, then ×f32(4.184) in f32) |
| any data var | int32 2 | — | `TypeError: 2 is not a float` |
| `foo` (non-AMBER var) | 2.0 | — | `NotImplementedError: scale_factors for variable foo are not implemented` |

## 6. Frame count, dt and time

* Frame count: §3 step 8.
* **dt** (time between frames) is computed lazily on first need as
  `time[1] − time[0]` using the **raw, unscaled** `time` values in their raw
  type (so a scale factor on `time` is ignored for dt — e.g. raw times 10,
  12.5 with scale 2 report times 20, 25 but dt 2.5). Type follows the raw
  variable: f32 times 0.1, 0.3 → dt f32 0.20000002; f64 → 0.19999999999999998;
  int32 3, 7 → int32 4.
* If there is no `time` variable or only one frame, dt falls back to 1.0
  (Python float) and **each** access of dt (including each computed time)
  emits the Warning `Reader has no dt information, set to 1.0 ps`.
* A user-supplied `dt` keyword overrides dt but does **not** override times
  read from the file (dt=0.5 with file times 10, 12.5 → times 10.0, 12.5,
  dt 0.5).
* With no `time` variable, time of frame k = k × dt (+ time_offset), as a
  Python float (0.0, 1.0, …).
* Time unit conversion is a no-op (picoseconds are native).

## 7. Iteration

Sequential iteration yields frames 0 … n−1 then repositions to frame 0
(after a full loop the current frame label is 0). Slicing / fancy indexing
use random access.

## 8. Baseline example

File: CDF-2, `Conventions=AMBER`, `ConventionVersion=1.0`, program attrs
present, dims frame(unlimited, 2 records), spatial 3, atom 3, cell_spatial 3,
cell_angular 3; `time` f32 [10, 12.5]; coordinates f32 = 1.1·i + 0.123456789;
cell_lengths f64 [30,40,50]; cell_angles f64 [90,100,110].

| quantity | value |
|---|---|
| frame count / atoms | 2 / 3 |
| periodic | yes |
| dt | f32 2.5 |
| frame 0 time | f32 10.0 |
| frame 0 atom 0 | [0.12345679104328156, 1.2234567403793335, 2.3234567642211914] |
| frame 0 box | f32 [30, 40, 50, 90, 100, 110] |
| frame 1 atom 0 | [10.023456573486328, 11.123456954956055, 12.223456382751465] |
| remarks | empty |
| warnings | none |
