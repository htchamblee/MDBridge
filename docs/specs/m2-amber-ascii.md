# AMBER ASCII trajectory (mdcrd) and restart (inpcrd) reading — behaviour specification (MDAnalysis 2.11.0)

Clean-room behaviour description, derived from running MDAnalysis 2.11.0 on
crafted files. "f32"/"f64" = IEEE binary32/binary64. "Warning" = Python
`UserWarning`. `<path>` = file name as supplied by the caller.

---

## Part A — formatted trajectory (mdcrd / trj / crdbox)

### A.1 Format selection

Extensions (case-insensitive) `mdcrd`, `trj`, `crdbox` select this reader;
`.mdcrd.gz` / `.mdcrd.bz2` are transparently decompressed. **`.crd` is not
this reader** — it selects the CHARMM CRD reader. Lines are read in text
mode with universal newlines: `\n`, `\r\n` and lone `\r` all end a line and
behave identically.

### A.2 Atom count is mandatory

The file does not store the atom count. It must be supplied (a Universe
always supplies the topology's count). Without it:
`ValueError: AMBER TRJ reader REQUIRES the n_atoms keyword`.
There is **no explicit check** that the supplied count matches the file;
a wrong count shows up only as parse/shape errors (A.7) or silently
mis-framed data.

### A.3 Line layout

* Line 1: title, ignored except for one check: if the title, after removing
  trailing whitespace (incl. the newline), is longer than 80 characters →
  `OSError: Header of AMBER formatted trajectory has more than 80 chars. This is probably not a AMBER trajectory.`
  (80 chars followed by trailing spaces is fine.)
* Each frame: the 3N coordinates x1 y1 z1 x2 … in consecutive lines of
  **10 fixed-width fields of 8 characters** (columns 1–8, 9–16, …, 73–80).
  Number of coordinate lines per frame L = ceil(3N / 10). The last line of a
  frame carries r = 3N mod 10 fields (if r = 0 the last line is a full line).
* Optional box line after each frame's coordinates (periodic files): three
  8-character fields a, b, c in columns 1–24.
* Characters beyond column 80 (or beyond the fields a line needs) are
  ignored, so an 11th value or trailing text on a 10-value line is ignored.

### A.4 Field conversion and precision

A field is the exact 8-character slice of the line (shorter if the line
ends earlier; the slice can then include the newline character or be
empty). It is converted with Python-`float` semantics:

* surrounding whitespace allowed (space, tab, `\v`, `\f`, newline);
* optional sign, digits, optional `.`, optional exponent with `e`/`E`;
  `5.`, `.5000`, `+1.5`, `-0.` OK;
* `nan`, `inf`, `Infinity` (any case, optional sign) OK;
* single underscores between digits are accepted (`1_0.5` → 10.5;
  `1__0` rejected);
* Unicode decimal digits are accepted (`١٢٣` → 123);
* rejected: empty or all-blank field, Fortran `D`/`d` exponents
  (`1.0D+01`), hex (`0x1`), embedded comma/space (`1,5`, `- 1.0`),
  trailing garbage (`1.0e5x`), overflow asterisks `********`.

Values run together are fine because the slicing is fixed-width
(`-100.001-200.501-300.001` → −100.001, −200.501, −300.001); free-format
lines that do not respect the 8-column grid fail.

Precision: the decimal text is converted to the nearest f64, then stored
rounded to nearest f32 (two roundings). Examples:
`   0.001` → 0.0010000000474974513; `   0.501` → 0.5009999871253967;
`  10.001` → 10.00100040435791; `12345678` → 12345678.0;
`  1.5E-1` → 0.15000000596046448. Box values the same way (box
`30.1234` → f32 30.1234).

### A.5 Line parsing rule (important for malformed input)

For each of the L coordinate lines of a frame:

1. Try to read all 10 fields. If all 10 convert, the line contributes **10
   values** — even if it is the last line of the frame, which needs only r.
2. If any of the 10 fails, read only the first r fields (r = 3N mod 10). If
   those r convert, the line contributes r values; if any of them fails, a
   `ValueError` is raised (message below). When r = 0 this fallback reads
   zero fields and the line contributes **nothing** (no error).

After L lines the collected values must number exactly 3N; otherwise
`ValueError: cannot reshape array of size <count> into shape (<N>,3)`.

Field error message (verbatim pattern):
`ValueError: FixedcolumnEntry(<start>,<stop>,'F'): Failed to read&convert <slice>`
where `<start>`,`<stop>` are the 0-based column bounds of the failing field
among the first r fields (e.g. `(0,8)`, `(16,24)`) and `<slice>` is the raw
text of that field, including any spaces or newline character it contains
(e.g. eight spaces, `\n`, or the 8 characters `30.0 40.`).

### A.6 Box detection (done once, at open)

* If N = 1: Warning `Trajectory contains a single atom: assuming periodic=False`;
  the file is treated as having no box lines (a file that does have them
  then reads the box lines as extra 1-atom frames: 3 frames + 3 box lines
  → 6 frames).
* Otherwise: frame 0's L coordinate lines are consumed (with the A.5 rule),
  then the **next line** is examined: count how many of its 10 8-char
  fields convert as numbers (any positions, not necessarily the first
  three). If exactly 3 convert, the file is periodic, else not.
  * If periodic, the first three fields of that line are parsed as the box
    (a, b, c); if one of them fails (e.g. the three numbers sit in columns
    25–48) the open fails with the field `ValueError` (A.5 pattern, e.g.
    `FixedcolumnEntry(0,8,'F'): Failed to read&convert ` + 8 spaces).
  * If there is **no** next line (single-frame file without box): the
    open fails with a bare `StopIteration` exception (empty message). A
    single frame followed by one blank line is accepted as a non-periodic
    1-frame file (the blank line has 0 numeric fields). A file with only
    the title line raises `EOFError` (empty message). A first frame cut
    short by end of file gives the shape error (e.g. one 10-value line for
    N=4: `ValueError: cannot reshape array of size 10 into shape (4,3)`).
    A completely empty (0-byte) file raises
    `EOFError: Compressed file ended before the end-of-stream marker was reached`.
* Box lines with 2, 4 or 6 numbers ⇒ non-periodic; those lines are then
  consumed as coordinate lines of following frames (garbage frames, or
  shape errors).
* After detection the reader restarts and reads frame 0 normally.

| line after first frame (N=4) | periodic | consequence |
|---|---|---|
| `  30.000  40.000  50.000` | yes | box 30,40,50,90,90,90 |
| `  30.000  40.000  50.000 abcdefg` | yes (text field does not convert) | box 30,40,50 |
| `   0.000   0.000   0.000` | yes | box all-zero ⇒ reported as "no box" (None) |
| `  -1.000  40.000  50.000` | yes | box −1,40,50,90,90,90 (not validated) |
| ` 30.1234 40.1     1e1` | yes | box 30.1234, 40.1, 10 |
| three numbers in columns 25–48 | yes | open fails: field `ValueError` on `(0,8)` |
| `30.0 40.0 50.0` (free format) | no (fields don't convert) | that line becomes frame 1 data ⇒ `ValueError … Failed to read&convert 30.0 40.` when frame 1 is read |
| 2, 4 or 6 numbers | no | the line is read as coordinates; e.g. 6-number line + next line give frame 1 atom 0 = (30, 40, 10.001), iteration ends with `cannot reshape array of size 4 into shape (4,3)` |
| absent (EOF) | — | `StopIteration` at open |

### A.7 Per-frame read

* Read L coordinate lines (A.5); if periodic, read one more line and parse
  its first three fields as a, b, c; box = (a, b, c, 90, 90, 90) stored as
  f32. Later box lines are parsed this way regardless of how many values
  they have (a 6-value box line just uses its first three).
* Box values are not validated. An all-zero (a=b=c=0) box is reported as
  "no box".
* Non-periodic files: no box.
* Velocities/forces: never present.
* **Time**: no time in the file. Time of frame k = k × dt (Python float).
  dt defaults to 1.0 and **every** access of dt or time then emits the
  Warning `Reader has no dt information, set to 1.0 ps`. A caller-supplied
  `dt` removes the warning (dt=2 → frame 2 time 4.0).

### A.8 Frame count

Computed by counting lines: let T = number of lines after the title line
(every line counts, blank lines included; a final line without newline
counts; a trailing `\n` does not create an extra line). Lines per frame
P = L, plus 1 if periodic. Frame count = floor(T / P). A partial trailing
frame is therefore not counted. Frame k starts at line 1 + k·P (after the
title).

| file (N=4, L=2) | frame count |
|---|---|
| 3 frames, no box | 3 |
| 3 frames, box lines | 3 (P = 3) |
| 2 frames, no trailing newline | 2 |
| 2 frames + one blank line at end | 2 |
| 2 frames + two blank lines at end | 3 (blank lines form a "frame") |
| 2 frames + 1 extra full line | 2 |
| periodic, 2 frames but last box line missing | 1 |
| N=1, 3 frames each followed by a box line | 6 |

### A.9 Random access vs. sequential iteration

* Indexing frame k (−n ≤ k < n, negatives from the end) seeks to its line
  and reads it; out of range →
  `IndexError: Index <k> exceeds length of trajectory (<n>).`
* Sequential iteration does **not** stop at the frame count: it keeps
  reading until the file is exhausted.
  * A partial trailing frame is read and fails:
    `ValueError: cannot reshape array of size 10 into shape (4,3)`.
  * Trailing blank line(s): the blank line fails the field parse:
    `ValueError: FixedcolumnEntry(0,8,'F'): Failed to read&convert ` + newline.
    (Indexing works for the counted frames; only iteration errors.)
  * Periodic file whose final box line is missing: iteration silently stops
    after the last complete frame (no error).
  * Clean end of file ends iteration normally; the reader then rewinds to
    frame 0.

### A.10 Malformed-input examples (N = 4 unless noted)

| input | result |
|---|---|
| title of 81 chars | `OSError: Header of AMBER formatted trajectory has more than 80 chars. This is probably not a AMBER trajectory.` |
| atom count not supplied | `ValueError: AMBER TRJ reader REQUIRES the n_atoms keyword` |
| supplied N=5 for a 4-atom file | `ValueError: FixedcolumnEntry(16,24,'F'): Failed to read&convert ` + newline (second line has 2 values, needs r=5) |
| supplied N=3 for a 4-atom file | `ValueError: cannot reshape array of size 10 into shape (3,3)` (full 10-value line taken) |
| blank 5th field on line 1 | line falls back to r=2 values ⇒ `cannot reshape array of size 4 into shape (4,3)` |
| `********` in field 1 | `ValueError: FixedcolumnEntry(0,8,'F'): Failed to read&convert ********` |
| last line has 10 values instead of 2 | `cannot reshape array of size 20 into shape (4,3)` |
| last line has 3 values (needs 2) | OK, the 3rd value is ignored |
| blank line right after title | `ValueError: FixedcolumnEntry(0,8,'F'): Failed to read&convert ` + newline |
| unpadded `0.001 0.501 …` | `ValueError: FixedcolumnEntry(0,8,'F'): Failed to read&convert 0.001 0.` |
| N=10 (r=0), 2nd line short (2 values) | short line contributes nothing ⇒ `cannot reshape array of size 20 into shape (10,3)` |
| `     nan`, `    -inf` fields | accepted as NaN / −inf |

### A.11 Normal examples

N=4, values v_i = 100·frame + 0.5·i + 0.001 written `%8.3f`, 3 frames:

| case | periodic | frames | frame 0 atom 0 | frame 2 atom 0 | box | time f0/f1/f2 |
|---|---|---|---|---|---|---|
| no box | no | 3 | [0.0010000000474974513, 0.5009999871253967, 1.0010000467300415] | [200.00100708007812, 200.50100708007812, 201.00100708007812] | none | 0.0 / 1.0 / 2.0 (+warnings) |
| box 30 40 50 | yes | 3 | same | same | [30,40,50,90,90,90] f32 | same |

Other observed sizes: N=2 (6 values, one line), N=3 (9 values, one line),
N=10 (30 values, three full lines) all read correctly with and without box.
A single-frame file **with** a box line reads fine (1 frame); a
single-frame file without one fails at open with `StopIteration`.

---

## Part B — restart / inpcrd (inpcrd, restrt)

### B.1 Format selection

Extensions `inpcrd` and `restrt` (case-insensitive). `rst7`, `rst`, `ncrst`
are **not** recognised. The file is opened as plain text (no
decompression), universal newlines.

### B.2 Layout and parsing

* Line 1: title = the line with leading and trailing whitespace removed.
* Line 2: split on whitespace. Token 1 = atom count, parsed as a base-10
  integer (`2.0` → `ValueError: invalid literal for int() with base 10: '2.0'`).
  If a second token exists it is the time, parsed as a float
  (`abc` → `ValueError: could not convert string to float: 'abc'`);
  further tokens are ignored. Empty line 2 or empty file →
  `IndexError: list index out of range`.
  Negative count → `OverflowError: can't convert negative value to uint64_t`.
  Count 0 is accepted (no atoms).
* Coordinates: for each pair of atoms one line with **6 fields of 12
  characters** (columns 1–12, …, 61–72): x,y,z of atom 2p+1 then of atom
  2p+2. If N is odd, one final line with 3 fields. Each field is the raw
  12-character slice converted with Python-float semantics (same acceptance
  rules as A.4: whitespace, `nan`, `inf`, `e`/`E` exponent OK; `D`
  exponent, blank field rejected). Characters beyond column 72 are ignored.
  Values may touch (`-123.4567890-234.5678901-345.6789012` OK).
* Errors are reported with Python's float message, the slice shown with
  Python repr quoting:
  * short line (field slice is just the newline): `ValueError: could not convert string to float: '\n'`
  * missing line (EOF): `ValueError: could not convert string to float: ''`
  * all-blank field: `ValueError: could not convert string to float: '            '`
  * free format `1.0 2.0 3.0 4.0…`: `ValueError: could not convert string to float: '1.0 2.0 3.0 '`
  * `  1.5D+00   `: `ValueError: could not convert string to float: '  1.5D+00   '`
* Nothing after the coordinates is read: **velocity lines and the box line
  are ignored**; the restart reader never reports a box (always "no box")
  and never reports velocities.
* Precision: decimal → f64 → f32. `   1.2345678` → 1.2345677614212036,
  `   0.1000000` → 0.10000000149011612, ` 123456.1234` → 123456.125,
  `-123.4567890` → −123.456787109375.

### B.3 Frame, time, atom count

* Exactly one frame (frame 0).
* Time: if line 2 had a time token, time = that value as a Python float
  (`1.0000000E+01` → 10.0); otherwise time = 0 × dt = 0.0 and the Warning
  `Reader has no dt information, set to 1.0 ps` is emitted on access.
* An atom count supplied by the caller is **ignored** by the reader (no
  error). When loaded into a Universe whose topology has a different count,
  the Universe raises (multi-line message, format list printed as a Python
  list):

```
ValueError: The topology and ['INPCRD', 'RESTRT'] trajectory files don't have the same number of atoms!
Topology number of atoms 5
Trajectory: <path> Number of atoms 2
```

* The standalone atom-count query reads line 2's first token the same way.

### B.4 Examples

| file | result |
|---|---|
| `TITLE line   ` / `    3  1.0000000E+01` / 6 values / 3 values / box line | 3 atoms, title `TITLE line`, positions [[1.1,2.2,3.3],[4.4,5.5,6.6],[7.7,8.8,9.9]] as f32 (1.100000023841858 …), time 10.0, no box |
| `T` / `    2` / one line of 6 values | 2 atoms, time 0.0 (+ dt warning) |
| `T` / `    2  5.0  7.0` / … | time 5.0 (third token ignored) |
| `T` / `    2  5.0` / coords / velocities line / box line | coords read, velocities and box ignored, time 5.0 |
| line 2 = `    2`, coordinate line has 5 fields | `ValueError: could not convert string to float: '\n'` |
| line 2 = `    3`, only one coordinate line | `ValueError: could not convert string to float: ''` |
| CRLF line endings | identical to LF |
