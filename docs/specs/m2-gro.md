# GROMACS GRO topology parser and coordinate reader — behaviour specification

Reference behaviour: MDAnalysis 2.11.0, format `GRO`. All examples were
produced by running MDAnalysis on crafted files. Columns are 1-based in
prose. `<FILE>` is the path as given.

A GRO file is read twice, independently:

* the **topology pass**, which gives atoms, names, ids and residues, and
* the **coordinate pass**, which gives positions, velocities and the box,
  for the **first frame only** (n_frames is always 1).

After that, types and masses are guessed from the atom names when a Universe
is built.

Files may be gzip/bzip2 compressed. LF and CRLF line ends both work. A
completely empty file fails with EOFError
`Compressed file ended before the end-of-stream marker was reached`.

---------------------------------------------------------------------------

## 1. File layout

| line | content |
|---|---|
| 1 | title. Ignored by both passes (the text is not kept). |
| 2 | atom count N |
| 3 … N+2 | one line per atom |
| N+3 | box line |

Lines after the box line are ignored. That includes later frames.

### 1.1 Atom count (line 2)

The whole line, including its terminator, is read as an integer using Python
integer rules: surrounding whitespace is allowed, as is a sign. Otherwise it
fails with ValueError `invalid literal for int() with base 10: '<line including \n>'`.
In a Universe this is wrapped as:

    Failed to construct topology from file <FILE> with parser <class 'MDAnalysis.topology.GROParser.GROParser'>.
    Error: invalid literal for int() with base 10: 'two\n'

| line 2 | result |
|---|---|
| `2` | 2 |
| `   1   ` | 1 |
| `two` | error `'two\n'` |
| `1 2` | error `'1 2\n'` |
| `1.0` | error `'1.0\n'` |

A file with only a title line fails with StopIteration (empty message).

---------------------------------------------------------------------------

## 2. Topology pass

### 2.1 Fields of an atom line (fixed columns)

| attribute | columns | treatment |
|---|---|---|
| resid | 1–5 | integer (Python rules; blanks around are fine; negatives allowed) |
| resname | 6–10 | whitespace stripped |
| name | 11–15 | whitespace stripped |
| id (atom number) | 16–20 | integer |

Coordinates are not read here. Tabs are not treated as separators. A line
that uses tabs puts text in the wrong columns and normally fails.

If the resid or id is not an integer (blank, `*****`, letters, a tab-shifted
line), parsing fails with an IOError (OSError) whose message has the
offending line appended, including its newline:

    Couldn't read the following line of the .gro file:
    <line>

In a Universe, IOError/OSError is wrapped as:

    Failed to load from the topology file <FILE> with parser <class 'MDAnalysis.topology.GROParser.GROParser'>.
    Error: Couldn't read the following line of the .gro file:
    <line>

| atom line | resid | resname | name | id |
|---|---|---|---|---|
| `    1SOL     OW    1   0.100   0.200   0.300` | 1 | SOL | OW | 1 |
| `    1SOL  O W    1 …` | 1 | SOL | `O W` | 1 |
| `    1SOL HW1      2 …` | 1 | `SOL H` | `W1` | 2 (fixed columns) |
| `   -5SOL     OW    0 …` | −5 | SOL | OW | 0 |
| `    1SOL     OW   -1 …` | 1 | SOL | OW | −1 |
| `  x  SOL …` | error: Couldn't read … | | | |
| `     SOL     OW    2 …` (blank resid) | error | | | |
| `    1SOL     OW***** …` | error | | | |

### 2.2 How many lines are read

Exactly N lines after line 2 are read as atom lines. If the file has fewer
lines than that:

* the box line is read as an atom line instead. Its resid columns normally
  fail, for example with
  `Couldn't read the following line of the .gro file:\n   1.00000   1.00000   1.00000\n`;
* if the file ends even earlier (no box line), the topology is left with
  unfilled atoms, and the result is IndexError
  `index 0 is out of bounds for axis 0 with size 0`.

N = 0 also fails with IndexError `index 0 is out of bounds for axis 0 with size 0`.

### 2.3 Missing names

After reading, if any atom's name is empty, parsing fails with IOError:

    Missing atom name on line: <L>

L is the 1-based file line number of the **first** atom with an empty name,
which is its atom index + 3. Example: atoms 2 and 3 unnamed gives
`Missing atom name on line: 4`. In a Universe this is wrapped in the same
"Failed to load from the topology file …" text as above.

### 2.4 Resid wrapping (resid 0 marks a wrap)

GRO resids have 5 digits. Writers wrap 99999 → 0, so MDAnalysis treats a
**block of consecutive atoms with resid 0** as the start of a new
100000-block:

* Find every maximal run of consecutive atoms whose resid is exactly 0.
* For each run that does **not** start at the very first atom, add 100000 to
  the resid of every atom from the run's first atom to the end of the file.
  Each such run adds again, so the increases accumulate.
* A run of zeros at the very start of the file is left alone (resid 0 is
  legal there).

The atom ids in columns 16–20 are never adjusted.

| resid column sequence | resulting per-atom resids | residues |
|---|---|---|
| 99999, 0, 1, 2 | 99999, 100000, 100001, 100002 | 4 |
| 99998, 99999, 0, 0, 1 | 99998, 99999, 100000, 100000, 100001 | 4 |
| 0, 1, 2 | 0, 1, 2 | 3 |
| 0, 0, 1, 0, 1 | 0, 0, 1, 100000, 100001 | 4 |
| 5, 0, 7, 0, 9 | 5, 100000, 100007, 200000, 200009 | 5 |
| 1, 0, 0, 0, 2 | 1, 100000, 100000, 100000, 100002 | 3 |
| 99999, 0, 99999, 0 | 99999, 100000, 199999, 200000 | 4 |

### 2.5 Residues and segments

* Residues are **consecutive runs** of atoms with the same (resid after
  wrapping, resname). Resid 1, 2, 1 gives 3 residues. Resid 1,1,1 with
  resnames SOL, NA, SOL gives 3 residues.
* resnums equal resids.
* There is a single segment named `SYSTEM`.
* No other attributes come from the file: no charges, elements, chain IDs or
  bonds.

### 2.6 Guessed attributes (Universe)

When a Universe is built, `types` are guessed from the names and `masses`
from the types. This uses exactly the name→type and type→mass rules given in
`pdb-topology.md` §9.3–§9.4. Examples: OW → O / 15.999, HW1 → H / 1.008,
NA → NA / 22.98977, CL → CL / 35.45, LI → LI / 0.0.

---------------------------------------------------------------------------

## 3. Coordinate pass

### 3.1 Field width ("precision") detection

The coordinate field width *w* is taken from the **first atom line only**,
and is used for every atom line:

* find the first `.` at or after **column 26** of the first atom line. Call
  its column (1-based) *c*;
* w = c − 25.

For files written as GROMACS writes them (each field is `%(d+5).df`), this is
d + 5: 8 for 3 decimals, 10 for 5 decimals, 7 for 2 decimals. The search
starts at column 26, so with standard widths the dot it finds is the one in
the **y** field. Edge cases:

| first atom line, columns 21+ | w | effect |
|---|---|---|
| `   0.100   0.200   0.300` | 8 | normal |
| `   0.12345   1.50000   2.25000` | 10 | normal |
| `   0.12   1.50   2.25` | 7 | normal |
| `1234.567   2.000   3.000` | 8 | normal (the x dot is at column 25 and is skipped) |
| `12345.6   0.200   0.300` (x dot at column 26) | 1 | fields are single characters: positions read as `1`,`2`,`3` (→ 10, 20, 30 Å) |
| `   1   2   3` (no dot) | 0 | every field is empty → ValueError `could not convert string to float: ''` |
| line ends at column 20 | 0 | same error |

### 3.2 Positions

For every atom line i (1 … N), x, y and z are the fields at columns
21…20+w, 21+w…20+2w and 21+2w…20+3w. Each is parsed as a real (Python
float syntax, surrounding blanks allowed) and rounded to float32. A field
that is blank or invalid fails with ValueError
`could not convert string to float: '<field text>'` (for example
`'        '` or `'   0.1x0'`). This error is not wrapped.

### 3.3 Velocities

Velocities are the next three fields, using the **same width w**: columns
21+3w … 20+6w. They are read for every atom line:

* If all three parse, they are stored.
* If parsing fails partway, the fields **before** the failing one are kept,
  and the failing field and the ones after it stay 0. A "missing velocity"
  flag is set. This also happens if the line simply has no velocity fields.
* After all atoms are read:
  * if **every** velocity component of every atom is exactly 0, the frame has
    **no velocities**, and there is no warning even if some lines lacked
    them;
  * otherwise the frame has velocities. Missing ones are 0. If the
    missing-velocity flag is set, the warning
    `Not all velocities were present.  Unset velocities set to zero.` is
    issued (two spaces).

| atom lines (columns 21+) | velocities (Å/ps) | warning |
|---|---|---|
| `   0.100   0.200   0.300` only | none | — |
| `… 0.300  0.0000  0.0000  0.0000` | none (all zero) | — |
| atom 1 `… 0.1000 -0.2000  0.3000`, atom 2 `… 1.0000  abcdef  3.0000` (w=8) | (1, −2, 3), (10, 0, 0) | yes |
| atom 1 no velocity, atom 2 `  1.0000  2.0000  3.0000` | (0, 0, 0), (10, 20, 30) | yes |
| `… 0.300  1.0000  2.0000` (two of three) | (10, 20, 0) | yes |

### 3.4 Unit conversion: nm → Å (float32 arithmetic)

Lengths are stored in nm and reported in Å. The conversion is done **in
single precision**: the float32 value parsed from the text is multiplied by
float32 10, and the product is rounded to float32. This applies to positions,
velocities (nm/ps → Å/ps) and the three box lengths. It is not the same as
rounding (text × 10) to float32 directly:

| text | float32(text) | reported (float32 × 10 in float32) | bits | float32(10 × text) for contrast |
|---|---|---|---|---|
| `0.001` | 0.0010000000474974513 | 0.010000000707805157 | 3c23d70b | 0.009999999776482582 |
| `0.003` | 0.003000000026077032 | 0.030000001192092896 | 3cf5c290 | 0.029999999329447746 |
| `3.333` | 3.3329999446868896 | 33.32999801635742 | 420551eb | 33.33000183105469 |
| `12.345` | 12.345000267028809 | 123.45000457763672 | 42f6e667 | 123.44999694824219 |
| `-0.001` | | -0.010000000707805157 | bc23d70b | -0.009999999776482582 |
| `1.234` | 1.2339999675750732 | 12.34000015258789 | 414570a4 | same |
| `0.007` | | 0.07000000029802322 | 3d8f5c29 | same |
| `2.345` | | 23.450000762939453 | 41bb999a | same |
| `999.999` | 999.9990234375 | 9999.990234375 | 461c3ff6 | same |
| `-99.999` | | -999.989990234375 | c479ff5c | same |
| `-7.777` | | -77.7699966430664 | c29b8a3d | same |
| `0.12345` | | 1.2345000505447388 | 3f9e0419 | same |
| `0.1` | 0.10000000149011612 | 1.0 | 3f800000 | same |
| `1.00001` | | 10.000100135803223 | 41200069 | same |

All of these were confirmed on positions, velocities and box lengths in real
GRO files. For example, the box line `   3.33300  12.34500   0.00100` gives
lengths 33.32999801635742, 123.45000457763672, 0.010000000707805157.

---------------------------------------------------------------------------

## 4. Box line

The box line is the line straight after the N-th atom line. If the file ends
before it, the result is UnboundLocalError
`cannot access local variable 'unitcell' where it is not associated with a value`.

### 4.1 Tokenising

1. Split the line on whitespace. If every token is a valid real (Python
   float syntax: `1`, `1e0` and `nan` are fine), those are the values,
   rounded to float32.
2. Otherwise (any token invalid), the values are **every match of "one or
   more digits, a dot, then exactly five digits"** found in the line, in
   order. This needs no separators, but **signs are lost**, and values
   without exactly 5 decimals are not found.
3. If the number of values is neither 3 nor 9 → ValueError
   `GRO unitcell has neither 3 nor 9 entries.` An empty line gives 0 values
   and so this error.

| box line | values |
|---|---|
| `   1.00000   2.00000   3.00000` | 1, 2, 3 |
| `   1 2 3` | 1, 2, 3 |
| `1e0 2e0 3e0` | 1, 2, 3 |
| `   1.00000   2.00000   3.00000   ` | 1, 2, 3 |
| `   1.00000   2.00000   3.00000xyz` | split fails → pattern: 1, 2, 3 |
| `1.000002.000003.00000` | pattern: 1, 2, 3 |
| `10.0000020.0000030.00000` | pattern: 10, 20, 30 |
| `-1.00000-2.00000-3.00000` | pattern: 1, 2, 3 (signs lost) |
| `-1 -2 -3` | −1, −2, −3 (split works, signs kept) |
| `   1.234   2.345   3.456  junk` | pattern finds nothing → `GRO unitcell has neither 3 nor 9 entries.` |
| `   1.00000` | error (1 value) |
| `   1.0 2.0 3.0 4.0` | error (4 values) |
| (empty line) | error (0 values) |

### 4.2 Three values (a, b, c in nm)

* If all three are close to 0 (each |v| ≤ 1e-8), the warning is
  ``Empty box [0., 0., 0.] found - treating as missing unit cell. Dimensions set to `None`.``
  and there is **no** unit cell.
* Otherwise the cell is (10·a, 10·b, 10·c, 90, 90, 90). The lengths use the
  float32 × 10 rule from §3.4, and there is no validation, so negative and
  nan values are kept. If, after conversion, all three lengths are exactly 0,
  the cell is reported as none.

| values | dimensions (Å, degrees) | warning |
|---|---|---|
| 1 2 3 | 10, 20, 30, 90, 90, 90 | |
| 0 0 0 (`0.00000` or `0.0`) | none | Empty box … |
| 1 0 0 | 10, 0, 0, 90, 90, 90 | |
| 0.00001 0 0 | 9.999999747378752e-05, 0, 0, 90, 90, 90 | |
| 1.23456 2.34567 3.45678 | 12.345600128173828, 23.45669937133789, 34.56779861450195, 90, 90, 90 | |
| −1 −2 −3 | −10, −20, −30, 90, 90, 90 | |
| nan 1 1 | nan, 10, 10, 90, 90, 90 | |

### 4.3 Nine values (GROMACS triclinic order)

The nine values are, in order:
v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y (in nm). They define three box
vectors:

* **a** = (v1x, v1y, v1z) — values 1, 4, 5
* **b** = (v2x, v2y, v2z) — values 6, 2, 7
* **c** = (v3x, v3y, v3z) — values 8, 9, 3

From the float32 values, working in double precision:

* the lengths are |a|, |b| and |c|;
* α is the angle between b and c, β between a and c, and γ between a and b,
  each in degrees (arccos of the normalised dot product);
* the six results are rounded to float32;
* **validity**: if any of the six is not > 0, or any angle is not < 180 (nan
  also fails), the cell is treated as all zeros and is reported as **none**,
  **with no warning**.

The three lengths are then converted to Å with the float32 × 10 rule. Angles
are not scaled.

| nine values | dimensions |
|---|---|
| 5 5 5 0 0 0 0 0 0 | 50, 50, 50, 90, 90, 90 |
| −5 5 5 0 0 0 0 0 0 | 50, 50, 50, 90, 90, 90 (length of (−5,0,0) is 5) |
| 5 5 5 0 0 2.5 0 2.5 2.5 | 50, 55.90169906616211, 61.23724365234375, 56.78908920288086, 65.90515899658203, 63.4349479675293 |
| 5 4 3 0 0 1 0 1 1.5 | 50, 41.231056213378906, 35.0, 60.98286056518555, 73.39845275878906, 75.96375274658203 |
| 3.333 3.333 3.333 0 0 1.6665 0 1.6665 1.6665 | 33.32999801635742, 37.26407241821289, 40.82074737548828, 56.78908920288086, 65.90515899658203, 63.4349479675293 |
| 1.23456 2.34567 3.45678 0 0 0 0 0 0 | 12.345600128173828, 23.45669937133789, 34.56779861450195, 90, 90, 90 |
| 5 5 0 0 0 0 0 0 0 | none (zero length, silently) |
| all nine 0 | none (silently, no "Empty box" warning) |

---------------------------------------------------------------------------

## 5. Other reader facts

* Only the first frame is read, even if the file holds more.
* The time step is not known: the reported step is 1.0 ps, with the warning
  `Reader has no dt information, set to 1.0 ps` when it is asked for, and
  frame 0 has time 0.0.
* Positions, velocities and the box are float32.

---------------------------------------------------------------------------

## 6. Message catalogue

| trigger | type | message |
|---|---|---|
| line 2 not an integer | ValueError | `invalid literal for int() with base 10: '<line with \n>'` (in a Universe: wrapped as `Failed to construct topology from file <FILE> with parser <class 'MDAnalysis.topology.GROParser.GROParser'>.\nError: …`) |
| bad resid or id field (or box line read as an atom) | IOError/OSError | `Couldn't read the following line of the .gro file:\n<line with \n>` (in a Universe: `Failed to load from the topology file <FILE> with parser <class 'MDAnalysis.topology.GROParser.GROParser'>.\nError: Couldn't read …`) |
| empty atom name | IOError/OSError | `Missing atom name on line: <atom index + 3>` (wrapped the same way) |
| N = 0, or the file ends before N atom lines and no box line | IndexError | `index 0 is out of bounds for axis 0 with size 0` |
| only a title line | StopIteration | *(empty)* |
| empty file | EOFError | `Compressed file ended before the end-of-stream marker was reached` |
| bad or blank coordinate field, or no `.` found for the width | ValueError | `could not convert string to float: '<field>'` |
| no box line | UnboundLocalError | `cannot access local variable 'unitcell' where it is not associated with a value` |
| box count ≠ 3 or 9 | ValueError | `GRO unitcell has neither 3 nor 9 entries.` |
| warning: 3-value zero box | UserWarning | ``Empty box [0., 0., 0.] found - treating as missing unit cell. Dimensions set to `None`.`` |
| warning: some velocities missing | UserWarning | `Not all velocities were present.  Unset velocities set to zero.` |
