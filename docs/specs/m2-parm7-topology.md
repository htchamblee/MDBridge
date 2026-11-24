# Amber parm7 / prmtop topology — behaviour specification

Reference behaviour: MDAnalysis 2.11.0 (topology formats `TOP`, `PRMTOP`, `PARM7`).
Every example below was produced by running MDAnalysis on crafted files. In
messages, `<FILE>` is the file name exactly as the caller passed it (not
normalised).

This document describes observable behaviour only.

---------------------------------------------------------------------------

## 1. Overview

The parser reads the file line by line, from top to bottom, in one pass. It:

1. checks two header lines,
2. skips forward to the `%FLAG POINTERS` section and reads the counts,
3. walks the remaining `%FLAG` sections. Sections it knows are read with a
   line count worked out from the POINTERS counts. Sections it does not know
   are skipped.
4. builds the topology from what it collected.

Files may be plain text or gzip/bzip2 compressed. Line terminators LF, CRLF
and lone CR are all accepted, and give identical results. A missing final
newline is fine.

---------------------------------------------------------------------------

## 2. Header requirements

### 2.1 Line 1

Line 1 must start with the three characters `%VE`, matched case-sensitively
at column 1. Anything after them is ignored, so `%VE`, `%VERSION ...` and
`%VERSIO x` are all accepted.

If it does not start with `%VE`, parsing fails with ValueError:

    <FILE> is not a valid TOP file. %VE Missing in header

| line 1 | result |
|---|---|
| `%VERSION  VERSION_STAMP = V0001.000 ...` | ok |
| `%VE` | ok |
| `%VERSIO x` | ok |
| `%ve...` | error (%VE Missing) |
| ` %VERSION` (leading space) | error (%VE Missing) |

### 2.2 Line 2

Line 2 is split on runs of whitespace, and its **second token** is examined.
The first token can be anything, and tokens after the second are ignored.

| line 2 | result |
|---|---|
| `%FLAG TITLE` | ok |
| `xx TITLE` | ok (first token not checked) |
| `%FLAG TITLE junk` | ok |
| `%FLAG CTITLE` | ValueError `<FILE> is detected as a Chamber-style TOP file. At this time MDAnalysis does not support such topologies` |
| `%FLAG NAME` | ValueError `<FILE> is not a valid TOP file. 'TITLE' missing in header` |
| `%FLAG title` | ValueError `... 'TITLE' missing in header` (case-sensitive) |
| `%FLAG` (only one token) | IndexError `list index out of range` |
| empty line | IndexError `list index out of range` |

If the file has fewer than two lines, a StopIteration with an empty message
is raised. This happens wherever the parser needs a line and the file has
ended (see §8).

### 2.3 Locating POINTERS

Starting from line 1, the parser skips lines until it finds one that
**starts with** `%FLAG POINTERS`. This is a prefix test, so
`%FLAG POINTERS_X` also matches. Lines between the title and POINTERS are
ignored, including the TITLE contents and any unknown sections. If no such
line exists, the file ends and a StopIteration (empty message) is raised.

### 2.4 Reading POINTERS

* The line straight after the `%FLAG POINTERS` line is discarded without
  being examined. It is normally `%FORMAT(10I8)`, but it is not checked.
* The next **exactly four lines** are read. Each one is split on whitespace,
  and every token is converted to an integer. The format is not used here,
  so field widths do not matter, but numbers that touch each other become one
  token. All the integers from the four lines are joined into one list, the
  "pointer list". It is indexed from 0.
* Any token that is not an integer fails with ValueError
  `invalid literal for int() with base 10: '<token>'`.
  * If POINTERS has only 3 data lines, the 4th line read is the next
    `%FLAG ...` line, and the error is `invalid literal for int() with base 10: '%FLAG'`.
  * If the discarded FORMAT line is missing, the first data line is
    discarded instead, and the same `'%FLAG'` error follows.
* The line after those four must contain `%FLAG` (see §3.1). If POINTERS has
  5 data lines, that line has no `%FLAG`, and the result is IndexError
  `list index out of range`. Note that this is **not** the custom "%FLAG
  section not found" message.
* Only these pointer-list entries are used (index: Amber name):

| index | Amber name | used for |
|---|---|---|
| 0 | NATOM | number of values in the per-atom sections |
| 2 | NBONH | bonds-with-hydrogen count |
| 3 | MBONA | bonds-without-hydrogen count |
| 4 | NTHETH | angles-with-hydrogen count |
| 5 | MTHETA | angles-without-hydrogen count |
| 6 | NPHIH | dihedrals-with-hydrogen count |
| 7 | MPHIA | dihedrals-without-hydrogen count |
| 11 | NRES | number of values in the per-residue sections |

* A pointer entry is only looked up when a section that needs it is met. If
  the list is too short for that entry (for example only 4 values in total),
  the result is IndexError `list index out of range`.

The atom and residue counts of the final topology do **not** come from
pointers 0 and 11 directly. They are the numbers of values actually
collected from ATOM_NAME and RESIDUE_LABEL (§6).

---------------------------------------------------------------------------

## 3. Walking the sections

### 3.1 The line after a section must be a `%FLAG` line

After POINTERS, and after every *known* section has been read (§4), the
parser reads one more line and takes the section name from it:

* The line must contain the text `%FLAG` somewhere. It does not have to be
  at the start: `xx%FLAG CHARGE` works. The name is everything after the
  first `%FLAG` (up to any later `%FLAG`), with leading and trailing
  whitespace removed. So `%FLAG   CHARGE   ` gives `CHARGE`, and
  `%FLAG CHARGE extra` gives `CHARGE extra`. The second is an unknown name,
  so that section is skipped.
* If the line has no `%FLAG`, parsing fails with IndexError (note the
  trailing space):

      %FLAG section not found, formatting error for PARM7 file <FILE> 

* Special case: if the section just read had a line count of 0 (§4.1), and
  the line read after it is empty or whitespace-only, that line is ignored
  and the next line is used instead. Only **one** such blank line is
  forgiven. A blank line after a section that has data, or two blank lines
  after an empty section, gives the "%FLAG section not found" error.
* If the file ends at this point, the walk finishes normally.

### 3.2 Unknown sections

A section whose name is not in the table in §4 is skipped. The parser reads
lines until it reaches one that **starts with** `%FLAG`, and uses that line
as the next section header. The contents of skipped sections are never
checked, so blank lines and garbage inside them are harmless. If the file
ends while skipping, the walk finishes normally.

Quirk: if the **first** section after POINTERS is an unknown section, parsing
fails with UnboundLocalError:

    cannot access local variable 'numlines' where it is not associated with a value

This needs another `%FLAG` line after it. If the unknown section runs to the
end of the file, the walk simply ends, and the missing-section KeyError
`'name'` follows (§4.2). An unknown section anywhere later is fine.

### 3.3 Repeated sections

If a known section appears twice, the later one replaces the earlier one
completely. For example, two RESIDUE_CHAINID sections: only the second
counts.

### 3.4 Section header lines (`%COMMENT` / `%FORMAT`)

For a known section, the lines after `%FLAG` are examined in order:

* a line whose first 8 characters are `%COMMENT` is skipped (any number of
  them);
* a line whose first 7 characters are `%FORMAT` gives the format. The format
  text is everything from character 9 onward (that is, after `%FORMAT` plus
  one more character, normally `(`), up to the first `)`. If there is no `)`
  it runs to the end of the line, including the newline;
* any other line fails with ValueError. The message ends with the offending
  line, including its line terminator:

      Invalid header line. Does not begin with either %FLAG, %FORMAT nor %COMMENT:
      <the line, including its newline>

  Examples that trigger it: a missing FORMAT line (the first data line is
  quoted), `%COMMEN hello`, and lower-case `%format(...)`.

| FORMAT line | format text | result |
|---|---|---|
| `%FORMAT(5E16.8)` | `5E16.8` | ok |
| `%FORMAT(5e16.8)` | `5e16.8` | ok (descriptors are case-insensitive) |
| `%FORMAT(5E16.8` | `5E16.8` + newline | ok |
| `%FORMAT (5E16.8)` | `(5E16.8` | ValueError `unrecognized FORTRAN format '(5E16.8'` |
| `%FORMAT(5D16.8)` | | ValueError `unrecognized FORTRAN format '5D16.8'` |
| `%FORMAT(1P,5E16.8)` | | ValueError `unrecognized FORTRAN format '1P'` |

### 3.5 Format descriptors → column layout

The format text is split on commas. Each item is read, ignoring case, as:
optional repeat count (digits), one letter from `I F E A X`, optional width
(digits), and optional `.decimals`. Only the start of each item has to match
this shape; trailing characters are ignored. An item that does not match even
with an implied repeat of 1 fails with ValueError
`unrecognized FORTRAN format '<item as written>'`.

* `I`, `F`, `E` and `A` items with repeat r and width w add r consecutive
  fields of w columns each.
* An `X` item skips `r` columns. Its width is ignored and taken as 1.
* Fields are laid out left to right from column 1.

The letter only fixes the field positions. How each field's text is converted
depends on the **section**, not on the letter (§4).

### 3.6 Reading the data lines of a known section

* The number of data lines is fixed **before** any are read:
  `ceil(count × values_per_record / per_line)`. Here `count` comes from the
  POINTERS entry, and `per_line` is a constant for each section (table in
  §4). It does **not** depend on the `%FORMAT`.
* Exactly that many lines are read, whatever they contain. A `%FLAG` line
  read this way is treated as data.
* On each line, the text in each field's columns is taken. A field past the
  end of the line gives empty text. Leading and trailing whitespace is
  removed (any whitespace, including tabs).
  * **Empty fields are dropped silently.** They do not produce a value and do
    not hold a place. So a short last line is fine, and so is a blank name
    field, but a blank field in the middle shifts all later values.
  * Text past the last field of the format is ignored. Fields defined by the
    format beyond the usual count *are* read: a `20a4` line with 24
    characters of names gives 6 names.
* All values of a section are joined in reading order.

Consequences (all checked against MDAnalysis):

| situation | effect |
|---|---|
| CHARGE declared `3E16.8` with 4 atoms | 1 line read (ceil(4/5)), 3 charges collected; the parser succeeds, but building a Universe fails (§7) |
| same, but the 4th charge is on a 2nd line | that 2nd line is where a `%FLAG` is expected → "%FLAG section not found" |
| CHARGE declared `5E20.8` with 20-wide fields | read correctly (field widths follow the format) |
| `%FORMAT(2X,5E16.8)` on standard 16-wide data | fields are misaligned; e.g. ValueError `could not convert string to float: '1.82223000E+01 -'` |
| NATOM=24 but only one ATOM_NAME line | 2 lines read; the 2nd is the `%FLAG CHARGE` line, read as names; the next line (`%FORMAT…`) is not `%FLAG` → "%FLAG section not found" |
| NATOM=6, CHARGE has 4 values on one line | 2 lines read; the 2nd is `%FLAG ATOMIC_NUMBER`, so ValueError `could not convert string to float: '%FLAG ATOMIC_NUM'` (the first 16 columns) |
| ATOMIC_NUMBER field `     6.0` | ValueError `invalid literal for int() with base 10: '6.0'` |
| ATOM_TYPE_INDEX `3000000000` written across 8-wide fields | split by columns into `3000000` and `000`, i.e. 3000000 and 0 |

### 3.7 Number syntax

* Integer fields accept an optional sign, leading zeros and decimal digits
  (`+6` → 6, `06` → 6). Python's integer syntax is used, so `1_0` would be
  10. Anything else fails with ValueError
  `invalid literal for int() with base 10: '<field text>'`. The quoted text
  is the stripped field.
* Real fields use Python float syntax: `1.0E+00`, `1e3`, `-2.5`, `inf`,
  `-inf` and `nan` (any case), and underscores between digits (`1_0` → 10).
  Fortran `D` exponents are **not** accepted: `1.00000000D+00` fails with
  ValueError `could not convert string to float: '1.00000000D+00'`. Hex is
  not accepted (`0x10` → same error with `'0x10'`).

---------------------------------------------------------------------------

## 4. Known sections

| %FLAG name | count source (pointer index) | values per record | per_line (for the line count) | conversion of each field | becomes |
|---|---|---|---|---|---|
| ATOM_NAME | 0 (NATOM) | 1 | 20 | stripped text | atom names |
| CHARGE | 0 | 1 | 5 | real → see §5.2 | charges |
| ATOMIC_NUMBER | 0 | 1 | 10 | integer → element symbol (§5.3) | elements |
| MASS | 0 | 1 | 5 | real (64-bit) | masses |
| ATOM_TYPE_INDEX | 0 | 1 | 10 | integer (stored as 32-bit) | type_indices |
| AMBER_ATOM_TYPE | 0 | 1 | 20 | stripped text | types |
| RESIDUE_LABEL | 11 (NRES) | 1 | 20 | stripped text | residue names |
| RESIDUE_POINTER | 11 | 1 | 10 | integer minus 1 | residue membership (§5.6) |
| BONDS_INC_HYDROGEN | 2 | 3 | 10 | integer floor-divided by 3 | bonds (2nd group) |
| BONDS_WITHOUT_HYDROGEN | 3 | 3 | 10 | integer floor-divided by 3 | bonds (1st group) |
| ANGLES_INC_HYDROGEN | 4 | 4 | 10 | integer floor-divided by 3 | angles (2nd group) |
| ANGLES_WITHOUT_HYDROGEN | 5 | 4 | 10 | integer floor-divided by 3 | angles (1st group) |
| DIHEDRALS_INC_HYDROGEN | 6 | 5 | 10 | integer floor-divided by 3 | dihedrals/impropers (2nd group) |
| DIHEDRALS_WITHOUT_HYDROGEN | 7 | 5 | 10 | integer floor-divided by 3 | dihedrals/impropers (1st group) |
| RESIDUE_CHAINID | 11 | 1 | 20 | stripped text | chain IDs / segments (§5.7) |

Line count = ceil(pointer × values-per-record / per_line). For example, 4
atoms give 1 CHARGE line; 6 atoms give 2 CHARGE lines; 3 bonds (9 integers)
give 1 line; 4 bonds (12 integers) give 2 lines.

### 4.1 Empty sections

If the count is 0, the section still needs its `%FORMAT` line, and no data
lines are read. A single blank line after it is allowed (§3.1). Writers such
as ParmEd emit one, but it is optional.

### 4.2 Required vs optional sections

Once the walk finishes, these must have been seen. If one is missing, a
KeyError is raised whose message is the quoted internal key. The checks run
in this order, and the first missing one is reported:

| missing section | KeyError message |
|---|---|
| ATOM_NAME | `'name'` |
| RESIDUE_POINTER | `'respoint'` |
| RESIDUE_LABEL | `'resname'` |
| BONDS_WITHOUT_HYDROGEN | `'bonda'` |
| BONDS_INC_HYDROGEN | `'bondh'` |
| ANGLES_WITHOUT_HYDROGEN | `'anga'` |
| ANGLES_INC_HYDROGEN | `'angh'` |
| DIHEDRALS_WITHOUT_HYDROGEN | `'diha'` |
| DIHEDRALS_INC_HYDROGEN | `'dihh'` |

Optional: CHARGE, MASS, AMBER_ATOM_TYPE and ATOM_TYPE_INDEX. If one is
absent, that attribute is simply not created; there is no warning. ATOMIC_NUMBER
is also optional, but if it is absent a warning is issued (§5.3). RESIDUE_CHAINID
is optional too (§5.7).

---------------------------------------------------------------------------

## 5. Attributes

### 5.1 Atom count, ids, names, types

* **n_atoms** is the number of non-empty ATOM_NAME fields collected. Blank
  name fields are dropped (§3.6), so a blank name reduces n_atoms. This
  usually shows up later as a length mismatch on another attribute (§7).
* **names**: the stripped field text. Spaces inside are kept (`C A`), and
  non-ASCII is kept (`Nä`). A tab inside a field is kept, because only the
  ends are stripped.
* **ids**: 1, 2, …, n_atoms. They are always generated and never read from
  the file.
* **types**: AMBER_ATOM_TYPE text, read in the same way as names.
* **type_indices**: ATOM_TYPE_INDEX integers, stored as signed 32-bit, so
  negative values are kept as they are.

### 5.2 Charges — exact arithmetic

For each CHARGE field:

1. parse the text as an IEEE double;
2. round it to IEEE **single** precision (float32), round-to-nearest-even.
   Values too large become ±inf, and very small ones become 0 or a
   subnormal;
3. divide, **in single precision**, by the single-precision value of
   18.2223. That value is 18.222299575805664 exactly (bits 0x4191C78F). The
   quotient is rounded to single precision;
4. the stored and reported value is that single-precision number, widened
   exactly to double.

This is *not* the same as dividing in double precision and then rounding to
single. They differ for many inputs:

| CHARGE text | float32(x) | reported charge (float32 quotient, shown as double) | float32 bits | (for contrast) float32(x / 18.2223 in double) |
|---|---|---|---|---|
| `1.0` | 1.0 | 0.05487781763076782 | 3d60c790 | 0.054877813905477524 |
| `-1.0` | -1.0 | -0.05487781763076782 | bd60c790 | -0.054877813905477524 |
| `2.0` | 2.0 | 0.10975563526153564 | 3de0c790 | 0.10975562781095505 |
| `0.5` | 0.5 | 0.02743890881538391 | 3ce0c790 | 0.027438906952738762 |
| `3.0` | 3.0 | 0.16463345289230347 | 3e2895ac | 0.16463343799114227 |
| `7.0` | 7.0 | 0.38414472341537476 | 3ec4ae9e | 0.38414469361305237 |
| `10.0` | 10.0 | 0.5487781763076782 | 3f0c7cba | 0.5487781167030334 |
| `-9.0` | -9.0 | -0.493900328874588 | befce081 | (same) |
| `18.2223` / `1.82223000E+01` | 18.222299575805664 | 1.0 | 3f800000 | (same) |
| `-18.2223` | | -1.0 | bf800000 | (same) |
| `3.32181810E+00` | 3.3218181133270264 | 0.1822941154241562 | 3e3aab4f | (same) |
| `-7.61156190E+00` | -7.6115617752075195 | -0.4177058935165405 | bed5dd8c | (same) |
| `8.24183910E+00` | 8.241839408874512 | 0.45229414105415344 | 3ee79319 | 0.45229411125183105 |
| `-1.49004819E+01` | -14.900482177734375 | -0.8177059292793274 | bf51552d | -0.8177058696746826 |
| `2.18356749E+01` | 21.835674285888672 | 1.1982941627502441 | 3f9961b4 | (same) |
| `12.345678901234` | 12.34567928314209 | 0.6775039434432983 | 3f2d70e6 | 0.6775038838386536 |
| `-3.70815013E+00` | | -0.20349517464637756 | | |
| `5.06425158E+00` | | 0.2779150605201721 | | |
| `0.1` | 0.10000000149011612 | 0.00548778148368001 | 3bb3d2d9 | (same) |
| `-0.4157` | | -0.022812707349658012 | bcbae1b7 | (same) |
| `0.0` | | 0.0 | 00000000 | |
| `-0.00000000E+00` | -0.0 | -0.0 (sign kept) | 80000000 | |
| `3.4028235e+38` | 3.4028234663852886e+38 | 1.8673952407202682e+37 | 7d60c78f | 1.8673951139552082e+37 |
| `1e39` | inf | inf | 7f800000 | |
| `1e-45` | 1.4e-45 (subnormal) | 0.0 | 00000000 | |
| `1e-50` | 0.0 | 0.0 | | |
| `inf` / `-inf` / `nan` | | inf / -inf / nan | | |

Charges are reported as 64-bit numbers, but they always hold one of these
float32 values.

### 5.3 Elements from ATOMIC_NUMBER

Each integer Z is mapped as follows:

* 1 ≤ Z ≤ 118: the standard symbol, with the first letter upper case and any
  second or third letter lower case: 1→`H`, 2→`He`, 6→`C`, 11→`Na`,
  17→`Cl`, 30→`Zn`, 118→`Og`.
* Z ≤ 0 (0, -1, -5, …): the empty string `""`.
* Z ≥ 119: parsing fails with KeyError, whose message is the number:
  `119` (Z = 200 gives `200`).

After the walk:

* if ATOMIC_NUMBER was absent, there is no elements attribute and a warning
  (UserWarning, also logged) is issued:

      ATOMIC_NUMBER record not found, elements attribute will not be populated. If needed these can be guessed using universe.guess_TopologyAttrs(to_guess=['elements']).

* if it was present and at least one element is `""`, **one** warning is
  issued however many atoms are affected:

      Unknown ATOMIC_NUMBER value found for some atoms, these have been given an empty element record. If needed these can be guessed using universe.guess_TopologyAttrs(to_guess=['elements']).

Elements are never guessed by this parser.

| ATOMIC_NUMBER values | elements | warning |
|---|---|---|
| 7 6 6 8 | N C C O | none |
| 0 -1 6 1 | "" "" C H | unknown-value warning (once) |
| -5 0 0 0 | "" "" "" "" | unknown-value warning (once) |
| 118 1 2 3 | Og H He Li | none |
| 119 1 2 3 | — | KeyError `119` |
| +6 06 6 8 | C C C O | none |

### 5.4 Masses

Each MASS field is parsed as a double and stored as given, with no
rounding to single precision: `14.01` → 14.01, `-2.5` → -2.5, `1e3` → 1000.0,
`0.000001` → 1e-06, `2.82245899E+01` → 28.2245899.

### 5.5 Residue names, resids, resnums

* **resnames**: the stripped RESIDUE_LABEL texts. **n_residues** is the number
  of non-empty labels collected. This can differ from NRES, and from the
  number of residue pointers.
* **resids** and **resnums** are both 1, 2, …, n_residues. Nothing about them
  is read from the file.

### 5.6 Residue membership from RESIDUE_POINTER

Each pointer p is the 1-based index of the first atom of a residue. It is
used as q = p − 1. Membership is assigned like this:

* Every atom starts in residue 0.
* For residue k = 0, 1, …, (number of pointers − 1) in order, the atoms
  from position q_k up to, but not including, q_{k+1} are put in residue k.
  For the last pointer, the end is n_atoms. A later residue overwrites
  an earlier assignment.
* The range bounds follow "slice" rules:
  * a negative bound b means n_atoms + b, and is 0 if that is still
    negative;
  * a bound greater than n_atoms is cut to n_atoms;
  * if, after this, start ≥ end, the residue gets no atoms.

The pointer count is whatever was collected, which may differ from the number
of labels. No validation is done. Pointers that are unsorted, zero, negative
or out of range are all accepted silently. If there are more pointers than
labels, atoms can end up with residue indices ≥ n_residues. The parser does
not complain, and building a Universe also succeeds. Examples (6 atoms):

| RESIDUE_POINTER values | atom → residue index |
|---|---|
| 1 3 5 | 0 0 1 1 2 2 |
| 1 | 0 0 0 0 0 0 |
| 2 4 | 0 0 0 1 1 1 (atom 0 was never covered, so it stays in 0) |
| 1 1 4 | 1 1 1 2 2 2 (residue 0 is empty) |
| 1 5 3 | 0 0 2 2 2 2 (residue 1 is empty, residue 2 covers 2..5) |
| 3 1 | 1 1 1 1 1 1 |
| 0 3 | 0 0 1 1 1 1 (first range is 5..2, empty) |
| -1 3 | 0 0 1 1 1 1 |
| 1 -2 | 0 0 0 1 1 1 (bound −3 means 3) |
| 1 7 | 0 0 0 0 0 0 (second range 6..6 is empty) |
| 1 10 | 0 0 0 0 0 0 |
| 6 | 0 0 0 0 0 0 |
| 1 3 5 with only 2 labels | 0 0 1 1 2 2, but n_residues = 2 |

### 5.7 Segments and chain IDs from RESIDUE_CHAINID

* The section is read as stripped text, NRES values at 20 per line. The
  field widths come from its own `%FORMAT`, so `20a8` with 8-wide fields
  works. Blank fields are dropped (§3.6).
* If the number of collected chain IDs **equals n_residues**:
  * residues are grouped into segments by **runs of consecutive residues with
    equal chain ID**. The comparison is exact and case-sensitive. Each run
    becomes one segment, in file order, and the segment's segid is the chain
    ID text. The same ID appearing in two separate runs gives two segments
    with the same segid;
  * a per-atom `chainIDs` attribute is created. Each atom gets the chain ID
    of the residue it belongs to, using the membership from §5.6;
  * n_segments is the number of runs.
* Otherwise, if the section was present (count mismatch, including the case
  where blank fields dropped some), a warning is issued (UserWarning, also
  logged):

      Number of residues (<n_res>) does not match number of %RESIDUE_CHAINID (<count>). Skipping section.

  and the section is ignored, as if it were absent.
* If the section is absent or ignored, there is one segment with segid
  `SYSTEM`, every residue belongs to it, and there is **no** chainIDs
  attribute.

Examples (4 residues of 2 atoms each):

| RESIDUE_CHAINID line (20a4) | segids | per-atom chainIDs | residue→segment | warning |
|---|---|---|---|---|
| `A   A   B   B   ` | A, B | A A A A B B B B | 0 0 1 1 | — |
| `A   B   A   B   ` | A, B, A, B | A A B B A A B B | 0 1 2 3 | — |
| `A   B   C   D   ` | A, B, C, D | A A B B C C D D | 0 1 2 3 | — |
| `a   A   a   a   ` | a, A, a | a a A A a a a a | 0 1 2 2 | — |
| `AB  CDE ABCDEFGH` | AB, CDE, ABCD, EFGH | AB AB CDE CDE ABCD ABCD EFGH EFGH | 0 1 2 3 | — |
| `A B A B A B A B ` | `A B` | `A B` ×8 | 0 0 0 0 | — |
| `A       B   B   ` (blank 2nd field) | SYSTEM | (none) | 0 0 0 0 | `Number of residues (4) does not match number of %RESIDUE_CHAINID (3). Skipping section.` |
| 16 spaces | SYSTEM | (none) | 0 0 0 0 | `... (4) ... (0). Skipping section.` |
| `A   B   ` | SYSTEM | (none) | 0 0 0 0 | `... (4) ... (2) ...` |
| `A   B   C   D   E   ` | SYSTEM | (none) | 0 0 0 0 | `... (4) ... (5) ...` |
| `SYSTEMX ` | SYSTEM | (none) | | `... (4) ... (2) ...` (fields `SYST`, `EMX`) |
| `A       B       C       D       ` with `%FORMAT(20a8)` | A, B, C, D | A A B B C C D D | 0 1 2 3 | — |
| two sections, `A A B B` then `X Y X Y` | X, Y, X, Y | X X Y Y X X Y Y | 0 1 2 3 | — |
| section is last in file, FORMAT line present, no data line | — | — | — | StopIteration (empty) |

The section may appear anywhere after POINTERS; position does not matter.

### 5.8 Bonds, angles, dihedrals, impropers

Every integer in these sections is floor-divided by 3. Floor division rounds
toward −∞: 4→1, 5→1, −3→−1, −4→−2, −7→−3, −10→−4. The values are then cut
into consecutive chunks of 3 (bonds), 4 (angles) or 5 (dihedrals). The
**last** integer of each chunk (the parameter index) is dropped, and the rest
form the tuple. Nothing checks whether the integers were multiples of 3.

The pointer count only decides how many lines are read (§3.6). Every value
present on those lines is used:

* POINTERS says 1 bond, but the line holds 6 integers: 2 bonds are produced.
* The line holds only `3 6` (no parameter index): a bond (1, 2) is produced.
  Precisely: each tuple is made of the *first N−1 values of its chunk*, where
  N is 3, 4 or 5. The chunk normally ends with the parameter index. If the
  final chunk is short, it may give fewer values than that:
  * a final chunk that still holds N−1 values works (bond `3 6` → (1,2);
    dihedral `0 3 6 -9` → improper (0,1,2,3));
  * bonds/angles with fewer: ValueError
    `bonds must be an iterable of tuples with 2 atom indices` /
    `angles must be an iterable of tuples with 3 atom indices`
    (e.g. bond values `3 6 1 9`; angle values `0 3 6 1 9 12`);
  * dihedrals with fewer than 4 values in the final chunk: IndexError
    `tuple index out of range`.
  * `3 6 1 9 12` (one full chunk, then a 2-value chunk) gives (1,2), (3,4).

**Orientation rule** (applies to bonds, angles, dihedrals and impropers): if
a tuple's first index is greater than its last index, the whole tuple is
reversed. Equal ends are left as they are.

**Bonds.** The list is all BONDS_WITHOUT_HYDROGEN tuples in file order,
followed by all BONDS_INC_HYDROGEN tuples in file order, each oriented. The
raw list is not de-duplicated, and self-bonds (i, i) are kept. Index values
are not range-checked by the parser. A negative index such as −1 is kept.
An index ≥ n_atoms is accepted by the parser, but building a Universe then
fails with IndexError `index 100 is out of bounds for axis 0 with size 12`
(numbers vary).

When a Universe's bond list is queried, it is the set of distinct bonds,
sorted ascending by (first, second).

| WITHOUT_H raw | INC_H raw | parser bond list | universe bond list |
|---|---|---|---|
| `30 33 1  6 3 2  0 3 1` | `3 0 1  9 6 1` | (10,11) (1,2) (0,1) (0,1) (2,3) | (0,1) (1,2) (2,3) (10,11) |
| `5 2 1` | `4 8 1` | (0,1) (1,2) | |
| — | `3 3 1` | (1,1) | (1,1) |
| — | `-3 6 1` | (−1,2) | |
| — | `-4 6 1` | (−2,2) | |
| `117 0 1` | | (0,39) | |

**Angles.** WITHOUT_H tuples then INC_H tuples, each oriented, with no
de-duplication. Example: INC_H `6 3 0 1  0 3 6 1  0 3 0 1` gives
(0,1,2) (0,1,2) (0,1,0).

**Dihedrals and impropers.** Take all WITHOUT_H 4-tuples followed by all
INC_H 4-tuples, and classify each in that order:

* 4th index negative → **improper**. The 3rd and 4th indices are replaced by
  their absolute values. Impropers keep this order and are not
  de-duplicated.
* else 3rd index negative → **dihedral**, with the 3rd index replaced by its
  absolute value (an Amber "no 1-4" flag).
* else → **dihedral** as is.

The dihedral collection is then de-duplicated (exact tuple equality), sorted
ascending lexicographically, and only **after that** oriented. As a result,
the final dihedral list can contain repeats and need not be sorted. "−0" is
0 and counts as non-negative. The absolute value is applied after floor
division, so −7 → −3 → 3, and −10 → −4 → 4.

| WITHOUT_H dihedral records (raw) | INC_H (raw) | dihedrals | impropers |
|---|---|---|---|
| `0 3 6 9 1`, `0 3 -6 9 1`, `0 3 6 -9 1`, `0 3 -6 -9 1`, `9 6 3 0 1`, `0 3 6 9 2` | `0 3 6 9 1`, `30 3 6 9 1` | (0,1,2,3) (0,1,2,3) (3,2,1,10) | (0,1,2,3) (0,1,2,3) |
| `0 3 -7 -10 1` | | — | (0,1,3,4) |
| `9 6 3 0 1`, `9 6 3 -0 1`, `12 3 6 -9 1` | | (0,1,2,3) | (3,2,1,4) |
| `30 33 36 39 1`, `30 33 -36 39 1`, `42 45 48 -51 1` | `0 3 6 9 1` | (0,1,2,3) (10,11,12,13) | (14,15,16,17) |

The bonds, angles, dihedrals and impropers attributes are always created.
They may be empty.

---------------------------------------------------------------------------

## 6. Counts summary

* n_atoms = number of ATOM_NAME values collected.
* n_residues = number of RESIDUE_LABEL values collected. NRES = 0 with an
  empty label section gives 0 residues. Atoms are then still put in residue
  0, and a Universe can still be built.
* n_segments = 1 (`SYSTEM`) or the number of chain-ID runs.
* NATOM = 0 with all sections empty is accepted (0 atoms, 0 residues, 1
  segment `SYSTEM`).

---------------------------------------------------------------------------

## 7. Errors when building a Universe from the parsed topology

The parser itself does not compare attribute lengths. When a Universe is
built, each per-atom attribute is checked against n_atoms, in the order the
attributes were created: file order of the sections, then bonds, angles,
dihedrals, impropers, ids, resids, resnums, segids, chainIDs. The first
mismatch raises ValueError:

    Length of <attr> does not match number of atoms.
    Expect: <n_atoms> Have: <length>

`<attr>` is one of `names`, `charges`, `elements`, `masses`, `types` or
`type_indices`. Examples:

| situation | message |
|---|---|
| MASS has 3 values for 4 atoms | `Length of masses does not match number of atoms.\nExpect: 4 Have: 3` |
| CHARGE format 3E16.8, 4 atoms (3 charges read) | `Length of charges does not match number of atoms.\nExpect: 4 Have: 3` |
| ATOM_NAME has 3 names, CHARGE (next in file) has 4 | `Length of charges does not match number of atoms.\nExpect: 3 Have: 4` |
| a blank ATOM_NAME field | same as above (names shrank) |
| AMBER_ATOM_TYPE has 2 values | `Length of types does not match number of atoms.\nExpect: 4 Have: 2` |
| a blank mass field (40 atoms) | `Length of masses does not match number of atoms.\nExpect: 40 Have: 39` |

Error wrapping at Universe level: a ValueError raised **by the parser** (not
the length check above) is re-raised as ValueError with this text:

    Failed to construct topology from file <FILE> with parser <class 'MDAnalysis.topology.TOPParser.TOPParser'>.
    Error: <original message>

IndexError, KeyError, StopIteration and UnboundLocalError pass through
unchanged.

---------------------------------------------------------------------------

## 8. Error catalogue (parser level)

| trigger | exception | exact message |
|---|---|---|
| line 1 not starting `%VE` | ValueError | `<FILE> is not a valid TOP file. %VE Missing in header` |
| 2nd token of line 2 is `CTITLE` | ValueError | `<FILE> is detected as a Chamber-style TOP file. At this time MDAnalysis does not support such topologies` |
| 2nd token of line 2 not `TITLE` | ValueError | `<FILE> is not a valid TOP file. 'TITLE' missing in header` |
| line 2 has < 2 tokens | IndexError | `list index out of range` |
| file ends anywhere a line is needed: before line 2, while looking for POINTERS, inside POINTERS, inside a known section's header or data | StopIteration | *(empty)* |
| non-integer POINTERS token (including a `%FLAG` line read as POINTERS data) | ValueError | `invalid literal for int() with base 10: '<token>'` |
| POINTERS has more than 4 data lines | IndexError | `list index out of range` |
| needed POINTERS entry beyond the list | IndexError | `list index out of range` |
| line after a known section (or after POINTERS) has no `%FLAG` | IndexError | `%FLAG section not found, formatting error for PARM7 file <FILE> ` (trailing space) |
| first section after POINTERS is unknown | UnboundLocalError | `cannot access local variable 'numlines' where it is not associated with a value` |
| line after `%FLAG` (and comments) not `%FORMAT`/`%COMMENT` | ValueError | `Invalid header line. Does not begin with either %FLAG, %FORMAT nor %COMMENT:\n<line incl. newline>` |
| bad format item | ValueError | `unrecognized FORTRAN format '<item>'` |
| bad integer field | ValueError | `invalid literal for int() with base 10: '<stripped field>'` |
| bad real field | ValueError | `could not convert string to float: '<stripped field>'` |
| atomic number ≥ 119 | KeyError | `<Z>` (e.g. `119`) |
| a required section missing | KeyError | `'name'`, `'respoint'`, `'resname'`, `'bonda'`, `'bondh'`, `'anga'`, `'angh'`, `'diha'`, `'dihh'` (first missing in this order) |
| final bond/angle chunk with too few values | ValueError | `bonds must be an iterable of tuples with 2 atom indices` / `angles must be an iterable of tuples with 3 atom indices` |
| final dihedral chunk with < 4 values | IndexError | `tuple index out of range` |

Warnings (UserWarning, also logged; parsing continues):

* `ATOMIC_NUMBER record not found, elements attribute will not be populated. If needed these can be guessed using universe.guess_TopologyAttrs(to_guess=['elements']).`
* `Unknown ATOMIC_NUMBER value found for some atoms, these have been given an empty element record. If needed these can be guessed using universe.guess_TopologyAttrs(to_guess=['elements']).`
* `Number of residues (<n>) does not match number of %RESIDUE_CHAINID (<m>). Skipping section.`

(Universe construction from a topology-only file also warns
`No coordinate reader found for <FILE>. Skipping this file.`; that warning
does not come from this parser.)
