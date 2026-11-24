# PDB topology parser and PDB coordinate reader — behaviour specification

Reference behaviour: MDAnalysis 2.11.0, formats `PDB` / `ENT`. All examples
were produced by running MDAnalysis on crafted files. Columns are 1-based
in prose, as in the PDB standard. `<FILE>` is the path as given.

Loading a PDB as a "Universe" has three parts:

1. **topology pass** (§1–§8): atoms, residues, segments, per-atom attributes,
   bonds;
2. **default guessing** (§9): types and masses are filled in when they are
   not provided;
3. **coordinate pass** (§10–§12): frames, positions, unit cell, per-frame
   occupancy/tempfactor, and HEADER/TITLE/COMPND/REMARK metadata.

The two passes read the file **independently, using different line-recognition
rules**. When they disagree on how many atoms there are, the coordinate pass
raises an error (§12).

---------------------------------------------------------------------------

## 1. Topology pass: which lines are atoms

The topology pass reads the file as text, one line at a time. For **each line,
leading and trailing whitespace is removed first**. All column positions
below apply to that stripped line. A line indented by k spaces therefore has
every field shifted k columns to the left.

* Empty lines (after stripping) are skipped.
* The first line that **starts with `END`** stops the pass. This includes
  `END`, `ENDMDL` and `ENDX`. Only the first model of a multi-model file
  becomes the topology. ATOM records after an `END…` line are ignored here,
  but the coordinate pass still sees them (§12).
* A line that **starts with `ATOM` or `HETATM`** is an atom record. This is a
  prefix test and case-sensitive, so `ATOMS…` and `ATOM100000…` count, and
  `atom` does not.
* All other lines are ignored, including TER, ANISOU, HEADER and so on.

`record_types` is columns 1–6 with whitespace stripped:

| line start | record_type |
|---|---|
| `ATOM      1 …` | `ATOM` |
| `HETATM    2 …` | `HETATM` |
| `ATOMS     3 …` | `ATOMS` |
| `ATOM100000 …` | `ATOM10` (and the serial field holds `0000 `, giving 0) |
| `HETATM10000 …` | `HETATM` (serial 10000) |

If there are no atom records at all, the topology pass fails with IndexError
`index 0 is out of bounds for axis 0 with size 0`.

---------------------------------------------------------------------------

## 2. Per-atom fields (columns of the stripped line)

A field that lies past the end of the line is empty text.

| attribute | columns | treatment |
|---|---|---|
| ids (serial) | 7–11 | §3 |
| names | 13–16 | whitespace stripped (`CA`, `1HB1`; `C A` keeps its inner space; all blank gives `""`) |
| altLocs | 17 | stripped (blank → `""`) |
| resnames | 18–21 | stripped; 4 characters are kept (`ABCD`) |
| chainIDs | 22 | stripped (blank → `""`) |
| resids | 23–26 | §4 |
| icodes | 27 | stripped |
| occupancies | 55–60 | real number; if the text is not a valid real (including blank or missing) → **0.0**; stored as float32 |
| tempfactors | 61–66 | real number; if not valid → **1.0**; stored as float32 |
| segids | 73–76 | stripped (§5) |
| element column | 77–78 | stripped (§6) |
| formal charge | 79–80 | stripped (§7) |

Coordinates (31–54) are **not** read by the topology pass.

Real-number syntax for occupancy and tempfactor follows Python float rules:
surrounding blanks are allowed, as are `1e-1`, `nan` and `-inf` (any case),
and `1_0.5` (= 10.5). The value is parsed as a double and rounded to
float32.

| occ/tempfactor text | occupancy | tempfactor |
|---|---|---|
| `  1.00` | 1.0 | 1.0 |
| 6 blanks | 0.0 | 1.0 |
| ` 0.5  ` | 0.5 | 0.5 |
| `abcdef` | 0.0 | 1.0 |
| `  1e-1` | 0.10000000149011612 | same |
| `   nan` | nan | nan |
| `  -inf` | -inf | -inf |
| ` 1_0.5` | 10.5 | 10.5 |
| `999999` | 999999.0 | 999999.0 |
| `1.2345` | 1.2345000505447388 | same |
| line ends before column 55 | 0.0 | 1.0 |

---------------------------------------------------------------------------

## 3. Atom serials (ids)

Columns 7–11, which is 5 characters if the line is long enough, are decoded
as follows:

1. If the text is a valid integer (Python rules: blanks around it, an
   optional sign and leading zeros are allowed), that integer is used. For
   example, `    1` → 1, `-9999` → −9999, ` -1  ` → −1, `1    ` → 1,
   `    0` → 0.
2. Otherwise, if the text is exactly 5 characters long, it is decoded as
   **hybrid-36**:
   * first character `-`, a blank or a digit → it should have been a decimal
     integer. It is not, so treat it as undecodable (step 3). Example:
     `  1 2`, `0000a`;
   * first character an upper-case letter A–Z → read all 5 characters as a
     base-36 number with digits `0–9A–Z` (values 0–35). The serial is that
     number − 10·36⁴ + 100000;
   * first character a lower-case letter a–z → read as base-36 with digits
     `0–9a–z`. The serial is that number + 16·36⁴ + 100000;
   * any other first character (e.g. `*`) → undecodable.
   * If any **later** character is not in the chosen alphabet (mixed case,
     blank, `-`, …), the parse **aborts with an error**, not a fallback:
     KeyError whose message is the quoted offending character, e.g. `'b'` or
     `' '`. If this happens on the very first atom record, the error is
     instead UnboundLocalError
     `cannot access local variable 'serial' where it is not associated with a value`.
3. If the text is undecodable (including text shorter than 5 characters that
   is not an integer), the atom gets a **made-up serial**. These come from a
   counter that starts at 100000 and goes up by 1 each time it is used, over
   the whole file. A "wrapped serials" flag is also set.

| columns 7–11 | serial |
|---|---|
| `99999` | 99999 |
| `A0000` | 100000 |
| `A0001` | 100001 |
| `ZZZZZ` | 43770015 |
| `a0000` | 43770016 |
| `zzzzz` | 87440031 |
| `*****` | 100000, then 100001 for the next such atom, … |
| `  1 2` | made-up (100000…) |
| `0000a` | made-up |
| `Ab000` | error: KeyError `'b'` (UnboundLocalError if first atom) |
| `A-000` | error: KeyError `'-'` (UnboundLocalError if first atom) |
| `A 000` | error: KeyError `' '` |
| `12` (line ends at column 8) | 12 |

If the wrapped-serials flag was set, then after the atoms are read:

* warning `Serial numbers went over 100,000.  Higher serials have been guessed`
  (note the two spaces);
* CONECT records are **not** parsed, and there is no bonds attribute. The
  warning `Invalid atom serials were present, bonds will not be parsed` is
  issued **twice in a row**. This happens whether or not the file has any
  CONECT records.

Serials are stored as 32-bit integers. Duplicate serials are allowed.

---------------------------------------------------------------------------

## 4. Residue numbers (resids) and wrapping

Columns 23–26 are parsed as an integer (Python rules: blanks and sign
allowed). Hybrid-36 is not supported.

* If parsing fails (blank, letters, `A000`, …): the resid is **1**, and the
  warning `PDB file is missing resid information.  Defaulted to '1'` (two
  spaces) is issued **for each such atom**. Such atoms do not change the
  "previous resid" used for wrapping.
* Otherwise, wrapping applies. Keep a "previous resid", initially 0. While
  (this resid − previous resid) < −5000, add 10000 to this resid. The result
  then becomes the new previous resid.

| resid column sequence | resulting resids |
|---|---|
| 9998, 9999, 0, 1, 9999, 0 | 9998, 9999, 10000, 10001, 9999, 10000 |
| 1, 9999, 1 | 1, 9999, 10001 |
| 9000, 3999, 3000, 4001 | 9000, 13999, 13000, 14001 |
| −501, 500 | −501, 500 |
| blank, 5, blank, 2 | 1, 5, 1, 2 (two warnings) |
| `A000`, 5 | 1, 5 (one warning) |
| 9999, blank, 0 | 9999, 1, 10000 |

`resnums` are equal to resids.

---------------------------------------------------------------------------

## 5. Residues, segments, chain IDs

### 5.1 Which per-atom value becomes the segid

* Normally the per-atom segid is columns 73–76, stripped.
* If **every** atom's segid is blank, the per-atom segid is the atom's chain
  ID instead.
* With the option `force_chainids_to_segids=True`, the per-atom segid is
  always the chain ID, whatever columns 73–76 contain.

### 5.2 Residues

Atoms are grouped into residues by **consecutive runs**. A new residue starts
at every atom where at least one of these differs from the previous atom:
resid (after wrapping), resname, insertion code, or per-atom segid (after
§5.1). The chain ID itself is **not** part of this test, except when it is
being used as the segid.

Each residue takes resid, resname, icode and segid from its atoms. Because
these are the run criteria, all atoms in a residue share them.

### 5.3 Segments

* If at least one residue segid is non-empty, **or** the force option is set,
  segments are the **consecutive runs of residues with equal segid**. Each run
  is one segment named by that segid. The same segid appearing in two
  separate runs gives two segments. A blank segid gives a segment named `""`.
* Otherwise (all segids blank and no force option), there is one segment
  named `SYSTEM`.

### 5.4 chainIDs

There is always a per-atom `chainIDs` attribute: column 22, stripped. It is
not affected by segid selection.

| atoms (resid, chain, segid col) | residues (resid/segid) | segments | note |
|---|---|---|---|
| (1,A,""),(1,A,""),(2,A,""),(1,B,""),(1,A,"") | 1/A, 2/A, 1/B, 1/A | A, B, A | segid from chain |
| (1,A,""),(1,B,"") | 1/A, 1/B | A, B | chain differs ⇒ segid differs ⇒ new residue |
| (1,A,PROT),(1,B,PROT) | one residue 1/PROT | PROT | chain change ignored; chainIDs still A, B |
| (1,A,PROT),(2,A,""),(3,A,PROT) | 3 residues | PROT, "", PROT | partial segids used as-is |
| same as above with chains A,B,B and force option | 3 residues | A, B | |
| (1," ",""),(2," ","") | 2 | SYSTEM | |
| same with force option | 2 | `""` (one segment) | |
| (1,A,""),(2," ",""),(3,A,"") | 3 | A, "", A | |
| icodes " ","A","A"," " at resid 1 | 3 residues: 1/"", 1/A, 1/"" | | |
| resnames ALA, GLY, ALA at resid 1 | 3 residues | | |
| segid columns `ABCD`, ` BC ` | | ABCD, BC | |

---------------------------------------------------------------------------

## 6. Elements and types from columns 77–78

Let *E* be the stripped element-column text of each atom.

* If **every** atom's E is empty, there is no `elements` attribute and no
  `types` attribute from this column. A warning is issued:

      Element information is missing, elements attribute will not be populated. If needed these can be guessed using universe.guess_TopologyAttrs(context='default', to_guess=['elements']).

* Otherwise:
  * `types` is set to E for every atom, **exactly as written**. Case is
    kept, invalid values are kept, and blanks stay `""`.
  * `elements`: E is "capitalised" (first character upper case, the rest
    lower case). If the result is one of the 118 element symbols (H … Og), it
    is used. Otherwise the element is `""`, and for **each such atom** this
    warning is issued, with `<E>` being the raw stripped text, possibly
    empty:

        Unknown element <E> found for some atoms. These have been given an empty element record. If needed they can be guessed using universe.guess_TopologyAttrs(context='default', to_guess=['elements']).

    A blank E in a column that is partly filled gives the message with two
    spaces: `Unknown element  found for some atoms. …`.

| cols 77–78 | type | element | warning |
|---|---|---|---|
| ` C` | `C` | `C` | |
| `C ` | `C` | `C` | |
| `CL` | `CL` | `Cl` | |
| `cl` | `cl` | `Cl` | |
| `cL` | `cL` | `Cl` | |
| `ZN` / `Zn` | as written | `Zn` | |
| ` h` | `h` | `H` | |
| `NA` | `NA` | `Na` | |
| `OG` | `OG` | `Og` | |
| `HE` | `HE` | `He` | |
| ` D` | `D` | `""` | `Unknown element D …` |
| ` X` | `X` | `""` | `Unknown element X …` |
| `1H` | `1H` | `""` | `Unknown element 1H …` |
| `Xx` | `Xx` | `""` | yes |
| `**` | `**` | `""` | yes |
| `Q ` | `Q` | `""` | yes |
| blank (others not blank) | `""` | `""` | `Unknown element  found …` |

---------------------------------------------------------------------------

## 7. Formal charges from columns 79–80

Let *Q* be the stripped text.

* If every atom's Q is empty, there is no `formalcharges` attribute and no
  warning.
* Otherwise each atom is converted. An empty Q is 0. `0` is 0. A Q containing
  `+` or `-` is **reversed** and read as an integer, so the sign must come
  **after** the digits. Anything else is invalid.
* At the first invalid entry, in atom order, conversion stops. The attribute
  is **not created at all**, and one warning is issued naming that entry:

      Unknown entry <Q> encountered in formal charge field. This likely indicates that the PDB file is not fully standard compliant. The formalcharges attribute will not be populated.

| Q | value |
|---|---|
| (blank) | 0 |
| `0` (as `0 ` or ` 0`) | 0 |
| `1+` | 1 |
| `2-` | −2 |
| `9+` | 9 |
| `0+` | 0 |
| `0-` | 0 |
| `+1` | invalid → `Unknown entry +1 …` |
| `-2` | invalid |
| `+` / `-` | invalid |
| `+-` / `--` | invalid |
| `10` | invalid (no sign) |
| `1` | invalid (no sign) |

Formal charges are integers.

---------------------------------------------------------------------------

## 8. CONECT records → bonds

CONECT lines are found by scanning the **whole file**, including the part
after `END`/`ENDMDL`. A line counts if its **raw** first six characters are
exactly `CONECT`, which is case-sensitive and not indented. Each one is
handled like this:

1. The line is stripped of surrounding whitespace.
2. Columns 7–11 must hold an integer (the base atom). If not, parsing fails
   with ValueError `invalid literal for int() with base 10: '<5 chars>'`
   (for example `'    a'`). As a Universe error, this is wrapped as
   `Failed to construct topology from file <FILE> with parser <class 'MDAnalysis.topology.PDBParser.PDBParser'>.\nError: invalid literal for int() with base 10: '    a'`.
   Hybrid-36 is **not** accepted here, so `CONECTA0000…` fails the same way
   with `'A0000'`.
3. Let *L* be the length of the stripped line from column 12 onward, and
   *n* = floor(L / 5).
   * n = 0 (L between 0 and 4) → warning
     `Found CONECT record with single entry, ignoring this`, and the line
     gives no bonds.
   * L not divisible by n → the whole CONECT parse is **abandoned**. There is
     no bonds attribute, even from earlier good lines, and the warning is
     `CONECT records was corrupt, bonds will not be parsed`.
     (The test really is divisibility by n, not by 5. L = 12 with n = 2 is
     accepted and reads 2 fields; L = 13 is rejected; L = 6…9 is accepted
     with 1 field.)
   * otherwise n partner fields of 5 characters each, from column 12, are
     read as integers. Any extra trailing characters are ignored. A field
     that is not an integer, including an all-blank field in the middle,
     fails with ValueError `invalid literal for int() with base 10: '     '`.
4. Each base→partner pair is looked up by serial. If a serial is duplicated,
   the **last** atom with that serial is the one used. If either serial is
   not an atom, the pair is dropped with the warning
   `PDB file contained CONECT record to TER entry. These are not included in bonds.`,
   once per dropped pair.

The bond set is **symmetrised**. A bond listed in only one direction exists.
A bond listed in both directions is one bond. Self-bonds (`CONECT 1 1`) are
kept as (i, i). The Universe bond list holds the distinct pairs (i ≤ j),
sorted ascending.

If no bond survives, there is **no bonds attribute** at all (not an empty
one).

| atoms with serials 1–4; CONECT lines | bonds (atom indices) | warnings |
|---|---|---|
| `CONECT    1    2` | (0,1) | |
| `CONECT    1    2    3` + `CONECT    2    1` | (0,1) (0,2) | |
| `CONECT    1    1` | (0,0) | |
| `CONECT    1` | none | single entry |
| `CONECT    1 222` (L=4) | none | single entry |
| `CONECT    1    2 x` (L=7) | (0,1) | |
| `CONECT    1    2    3 x` (L=12) | (0,1) (0,2) | |
| `CONECT    1    2    3xxx` (L=13) | none (all bonds dropped) | corrupt |
| `CONECT    1    2    3    4x` (L=16) | none | corrupt |
| `CONECT    1    9` + `CONECT    9    1` | none | TER-entry warning ×2 |
| `CONECT    1    5    2` | (0,1) | TER-entry ×1 |
| `CONECT    1         2` | — | ValueError `'     '` |
| `END` then `CONECT    3    4` | (2,3) | |
| `conect    1    2` | none | |
| ` CONECT    1    2` (indented) | none | |
| `CONECT    1    2    3    4    1    2    3    4    1    2    3    4` | (0,0) (0,1) (0,2) (0,3) | |
| atoms with serials 1,1,2; `CONECT    1    2` | (1,2) | |

---------------------------------------------------------------------------

## 9. Guessing types and masses (when building a Universe)

By default, building a Universe guesses **types** and **masses**, in that
order. Elements are **not** guessed by default.

### 9.1 Types

* If the topology already has `types` (§6: the element column had at least
  one non-blank entry), they are kept as they are. Atoms whose element column
  was blank keep type `""`; they are not guessed.
* Otherwise every atom's type is guessed from its **name** (§9.3). Guessed
  types are upper case.

### 9.2 Masses

* If an `elements` attribute exists, each mass is looked up from the
  element (§9.4). Elements that are `""` (invalid or blank) get 0.0.
* Otherwise each mass is looked up from the type.
* A lookup that fails gives 0.0, with a PendingDeprecationWarning
  `Unknown masses are set to 0.0 for current version, this will be deprecated in version 3.0.0 and replaced by Masse's no_value_label (np.nan)`.
  This warning class is hidden by default in Python.

Masses are 64-bit.

| element column | name | → type | → element | → mass |
|---|---|---|---|---|
| all blank in file | CA | C | (none) | 12.011 |
| all blank in file | SOD | NA | (none) | 22.98977 |
| all blank in file | 1HB | H | (none) | 1.008 |
| all blank in file | LI | LI | (none) | 0.0 (upper-case `LI` is not in the mass table) |
| blank (others filled) | CA | `""` | `""` | 0.0 |
| `NA` | SOD | NA | Na | 22.989768 (from `Na`) |
| ` X` | 1HB | X | `""` | 0.0 |
| `CA` | CA | CA | Ca | 40.08 |
| `LI` | LI | LI | Li | 6.941 |
| `cl` | CL | cl | Cl | 35.45 |

### 9.3 Guessing an element/type from an atom name

These rules are listed in priority order; the first one that applies gives
the result. "Upper-casing" means upper-casing letters.

Data:

* **Special-name table** (key → result). Keys are matched against the
  upper-cased name:

| name | → | name | → | name | → |
|---|---|---|---|---|---|
| BR | BR | LI+ | LI | CS | CS |
| CAL | CA | QL | LI | CS+ | CS |
| C0 | CA | MG | MG | QC | CE |
| CA2+ | CA | MG2+ | MG | RB | RB |
| CES | CS | K | K | QR | RB |
| CLA | CL | POT | K | BC | C |
| CLAL | CL | K+ | K | AC | C |
| CL | CL | QK | K | MW | DUMMY |
| CL- | CL | SOD | NA | FE | FE |
| IOD | I | NA | NA | FE2 | FE |
| LIT | LI | NA+ | NA | ZN | ZN |
| LI | LI | QN | NA | CU | CU |

* **Short element list**: H, LI, BE, B, C, N, O, F, NA, MG, AL, P, SI, S, CL,
  K.

Rules:

1. The empty name gives `""`.
2. If the upper-cased name is a special-name key → the table value.
3. Form the **stem**: remove every `*`, `+` and `-`, then remove every digit,
   then upper-case. If the stem is a special-name key → the table value.
4. Otherwise examine the stem, and then successively shorter versions of it
   made by dropping its **last** character. For each version V, take the
   first of these that holds:
   * V is in the short element list → V;
   * V without its last character is in the list → that;
   * V without its first character is in the list → that;
   * V has 2 or fewer characters → the first character of V.

   If none hold, drop V's last character and continue.
5. If the stem is empty (the name was only `*`, `+`, `-` and digits), the
   result is the name with `*`, `+` and `-` removed. Digits and case are kept,
   so `123` → `123`, `1` → `1`, `*` → `""`.

Example table (name → guessed type → mass from that type):

| name | type | mass | | name | type | mass |
|---|---|---|---|---|---|---|
| CA | C | 12.011 | | CAL | CA | 40.08 |
| C | C | 12.011 | | CAL2 | CA | 40.08 |
| N | N | 14.007 | | CA2+ | CA | 40.08 |
| O | O | 15.999 | | CES | CS | 132.9 |
| H | H | 1.008 | | LIT | LI | 0.0 |
| HA | H | 1.008 | | LI | LI | 0.0 |
| HB1 | H | 1.008 | | LI+ | LI | 0.0 |
| 1HB | H | 1.008 | | QL | LI | 0.0 |
| 2HG1 | H | 1.008 | | IOD | I | 126.9045 |
| 3HH2 | H | 1.008 | | FE / FE2 | FE | 55.847 |
| CB / CG / CD1 / CZ / CT | C | 12.011 | | FEZ | F | 18.998 |
| OG1 / OXT / OT1 / OH / OC1 / OHH | O | 15.999 | | CU | CU | 63.546 |
| SG / SD | S | 32.06 | | QN | NA | 22.98977 |
| NZ / NE / NE2 / N1 | N | 14.007 | | QK | K | 39.102 |
| P | P | 30.974 | | QC | CE | 0.0 |
| OP1 / O1P | O | 15.999 | | QR | RB | 85.4678 |
| C1' | C | 12.011 | | BC / AC | C | 12.011 |
| H5'' / HO2' | H | 1.008 | | MW / MW1 | DUMMY | 0.0 |
| MG / MG2+ | MG | 24.305 | | BR / Br | BR | 79.904 |
| ZN | ZN | 65.37 | | OW | O | 15.999 |
| NA / Na / na / NA+ / SOD / SOD+ / nA | NA | 22.98977 | | HW1 / HW2 / H1 / HN / HT1 / HXT / HE / HE2 | H | 1.008 |
| CL / Cl / cl / CL- / CLA / CLX | CL | 35.45 | | OH2 | O | 15.999 |
| K / K+ / POT / KX | K | 39.102 | | SE | S | 32.06 |
| SI | SI | 0.0 | | SOL | S | 32.06 |
| AL | AL | 0.0 | | BE | BE | 0.0 |
| B | B | 10.811 | | F | F | 18.998 |
| S | S | 32.06 | | XX / X | X | 0.0 |
| Q / QZ | Q | 0.0 | | Z / ZZZ | Z | 0.0 |
| 1 | 1 | 0.0 | | 123 | 123 | 0.0 |
| `*` / `+` / `-` | `""` | 0.0 | | `""` | `""` | 0.0 |
| H* | H | 1.008 | | C* / C- / *CA / ++C | C | 12.011 |
| AB / BA / BAC / ABCD / YB | B | 10.811 | | ZNA / NAX | NA | 22.98977 |
| LP | P | 30.974 | | DUM / DA | D | 0.0 |
| CH3 | C | 12.011 | | 1H / h | H | 1.008 |
| c / ca | C | 12.011 | | o | O | 15.999 |
| mg | MG | 24.305 | | Nn | N | 14.007 |

### 9.4 Mass table and lookup

A mass is looked up by **exact key first, then by the upper-cased key**. If
neither is in the table, the mass is 0.0. Table keys have mixed case
(data):

```
Ac 227.028  Al 26.981539  Am 243  Sb 121.757  Ar 39.948  As 74.92159  At 210
Ba 137.327  Bk 247  Be 9.012182  Bi 208.98037  Bh 264  B 10.811  BR 79.90400
Cd 112.411  CA 40.08000  Cf 251  C 12.01100  Ce 140.11600  CS 132.90000
CL 35.45000  Cr 51.9961  Co 58.9332  CU 63.54600  Cm 247  Db 262  Dy 162.5
Es 252  Er 167.26  Eu 151.965  Fm 257  F 18.99800  Fr 223  Gd 157.25
Ga 69.723  Ge 72.61  Au 196.96654  Hf 178.49  Hs 265  HE 4.00260
Ho 164.93032  H 1.00800  In 114.82  I 126.90450  Ir 192.22  FE 55.84700
Kr 83.8  La 138.9055  Lr 262  Pb 207.2  Li 6.941  Lu 174.967  MG 24.30500
Mn 54.93805  Mt 266  Md 258  Hg 200.59  Mo 95.94  N 14.00700  NA 22.98977
Nd 144.24  NE 20.17970  Np 237.048  Ni 58.6934  Nb 92.90638  No 259
Os 190.2  O 15.99900  Pd 106.42  P 30.97400  Pt 195.08  Pu 244  Po 209
K 39.10200  Pr 140.90765  Pm 145  Pa 231.0359  Ra 226.025  Rn 222
Re 186.207  Rh 102.9055  RB 85.46780  Ru 101.07  Rf 261  Sm 150.36
Sc 44.95591  Sg 263  Se 78.96  Si 28.0855  Ag 107.8682  Na 22.989768
Sr 87.62  S 32.06000  Ta 180.9479  Tc 98  Te 127.6  Tb 158.92534
Tl 204.3833  Th 232.0381  Tm 168.93421  Sn 118.71  Ti 47.88  W 183.85
U 238.0289  V 50.9415  Xe 131.29  Yb 173.04  Y 88.90585  ZN 65.37000
Zr 91.224  DUMMY 0.0
```

Results of this lookup for capitalised element symbols, which is how they
appear when masses come from the `elements` attribute:

H 1.008, He 4.0026, Li 6.941, Be 9.012182, B 10.811, C 12.011, N 14.007,
O 15.999, F 18.998, Ne 20.1797, Na 22.989768, Mg 24.305, Al 26.981539,
Si 28.0855, P 30.974, S 32.06, Cl 35.45, Ar 39.948, K 39.102, Ca 40.08,
Sc 44.95591, Ti 47.88, V 50.9415, Cr 51.9961, Mn 54.93805, Fe 55.847,
Co 58.9332, Ni 58.6934, Cu 63.546, Zn 65.37, Ga 69.723, Ge 72.61,
As 74.92159, Se 78.96, Br 79.904, Kr 83.8, Rb 85.4678, Sr 87.62,
Y 88.90585, Zr 91.224, Nb 92.90638, Mo 95.94, Tc 98, Ru 101.07,
Rh 102.9055, Pd 106.42, Ag 107.8682, Cd 112.411, In 114.82, Sn 118.71,
Sb 121.757, Te 127.6, I 126.9045, Xe 131.29, Cs 132.9, Ba 137.327,
La 138.9055, Ce 140.116, Pr 140.90765, Nd 144.24, Pm 145, Sm 150.36,
Eu 151.965, Gd 157.25, Tb 158.92534, Dy 162.5, Ho 164.93032, Er 167.26,
Tm 168.93421, Yb 173.04, Lu 174.967, Hf 178.49, Ta 180.9479, W 183.85,
Re 186.207, Os 190.2, Ir 192.22, Pt 195.08, Au 196.96654, Hg 200.59,
Tl 204.3833, Pb 207.2, Bi 208.98037, Po 209, At 210, Rn 222, Fr 223,
Ra 226.025, Ac 227.028, Th 232.0381, Pa 231.0359, U 238.0289, Np 237.048,
Pu 244, Am 243, Cm 247, Bk 247, Cf 251, Es 252, Fm 257, Md 258, No 259,
Lr 262, Rf 261, Db 262, Sg 263, Bh 264, Hs 265, Mt 266; Ds, Rg, Cn, Nh,
Fl, Mc, Lv, Ts, Og → 0.0.

For **upper-case** symbols, as produced by name guessing or by an upper-case
type, only these are non-zero: H 1.008, HE 4.0026, B 10.811, C 12.011,
N 14.007, O 15.999, F 18.998, NE 20.1797, NA 22.98977, MG 24.305, P 30.974,
S 32.06, CL 35.45, K 39.102, CA 40.08, V 50.9415, FE 55.847, CU 63.546,
ZN 65.37, BR 79.904, RB 85.4678, Y 88.90585, I 126.9045, CS 132.9,
W 183.85, U 238.0289, DUMMY 0.0. Every other upper-case symbol (LI, BE, AL,
SI, AR, …, CE, …) gives **0.0**.

---------------------------------------------------------------------------

## 10. Coordinate pass: frames

The coordinate pass reads the file as raw bytes.

* **MODEL detection**: every line whose first 5 bytes are `MODEL` (so
  `MODELX` also counts) marks the start of a frame. The frame begins just
  **after** that line.
* **CRYST detection**: every line whose first 5 bytes are `CRYST` is noted at
  the **start** of the line.
* If there are no MODEL lines, there is one frame starting at the beginning
  of the file.
* If the number of CRYST lines equals the number of frames, frame k starts at
  whichever comes first: the k-th MODEL position or the k-th CRYST position.
  This allows CRYST1 either before or after its MODEL line. Otherwise frames
  start at the MODEL positions only.
* Frame k runs from its start up to the start of frame k+1, or to the end of
  the file. Anything before the first frame start is never read as
  coordinates. That includes ATOM lines and a lone CRYST1 that come before
  the first MODEL when the counts differ.
* n_frames is the number of frames. The time step is not known; it reports
  1.0 ps, with the warning `Reader has no dt information, set to 1.0 ps`
  when the step is queried. Frame k has time k × 1.0.

Inside a frame, the text is decoded (UTF-8) and split into lines (LF, CRLF and
CR are all accepted). Lines are **not** stripped:

* a line whose first 6 characters are exactly `ATOM␣␣` (ATOM plus two
  spaces) or `HETATM` is an atom line, in order;
* a line whose first 6 characters are `CRYST1` is a unit cell (§11);
* every other line is ignored. END/ENDMDL do **not** stop the frame here.

For each atom line:

* x, y, z are columns 31–38, 39–46 and 47–54. Each is parsed as a real
  number (Python float syntax, surrounding blanks allowed) and rounded to
  float32. A field that is blank or invalid fails with ValueError
  `could not convert string to float: '<the 8-character field>'`, for
  example `'        '` or `'   abc  '`. Values are in Å; no unit conversion
  is done.
* the per-frame occupancy (columns 55–60) defaults to 0.0 if it is not a
  valid number. The per-frame tempfactor (61–66) defaults to 1.0. Both are
  kept as 64-bit values in the frame's data. The `tempfactor` entry exists
  **only if** at least one atom in that frame had a valid tempfactor; the
  `occupancy` entry always exists.

| columns 31–54 | position |
|---|---|
| `   1.000   2.000   3.000` | (1, 2, 3) |
| `1.0e+002   2.000   3.000` | (100, 2, 3) |
| `     nan     inf    -inf` | (nan, inf, −inf) |
| ` 1_0.000   2.000   3.000` | (10, 2, 3) |
| `1234.5678.0000002.000000` | (1234.5670166015625, 8.0, 2.0) (fixed columns) |
| `12345678.12345678.1234567` | (12345678.0, 0.12345670163631439, 8.123456001281738) |
| `   1.000           3.000` | ValueError `could not convert string to float: '        '` |
| `   abc     2.000   3.000` | ValueError `could not convert string to float: '   abc  '` |

---------------------------------------------------------------------------

## 11. CRYST1 handling

For each `CRYST1` line in a frame, in order, with the last one winning:

* the six fields a, b, c (columns 7–15, 16–24, 25–33) and α, β, γ (34–40,
  41–47, 48–54) are parsed as reals and rounded to float32;
* if any field is not a valid real (blank, cut short, misaligned): warning
  `Failed to read CRYST1 record, possibly invalid PDB file, got:\n<the line without its terminator>`,
  and the cell becomes "none";
* else if all six values are close to (1, 1, 1, 90, 90, 90) — closeness is
  defined below — the warning is
  `1 A^3 CRYST1 record, this is usually a placeholder. Unit cell dimensions will be set to None.`
  and the cell becomes "none";
* else the cell is the six float32 values, with no further checks: negative
  lengths, zero angles and nan are all kept.
* Closeness: each float32 value v (widened to double) must satisfy
  |v − r| ≤ 1e-8 + 1e-5·|r|, where r is the matching reference value.
* A cell whose three lengths are all exactly 0 is reported as "none",
  silently.

These warnings are issued each time the frame is read. Building a Universe
reads frame 0 more than once, so a warning can appear several times.

| CRYST1 a b c α β γ | dimensions | warning |
|---|---|---|
| 10 20 30 90 90 90 | (10, 20, 30, 90, 90, 90) | |
| 10 20 30 60 70 80 | (10, 20, 30, 60, 70, 80) | |
| 1 1 1 90 90 90 | none | placeholder |
| 1.000009 1 1 90 90 90 | none | placeholder |
| 1.00001 1 1 90 90 90 | (1.0000100135803223, 1, 1, 90, 90, 90) | — (just outside the tolerance after float32 rounding) |
| 1.00002 1 1 90 90 90 | (1.0000200271606445, …) | |
| 0.999991 1 1 90 90 90 | none | placeholder |
| 0.99999 1 1 90 90 90 | (0.9999899864196777, …) | |
| 1 1 1 90.0008 90 90 | none | placeholder |
| 1 1 1 90.0009 90 90 | (…, 90.00090026855469, …) | |
| 1 1 1 89.9991 90 90 | (…, 89.99909973144531, …) | |
| 1 1 1 90 90 90.01 | (…, 90.01000213623047) | |
| 0 0 0 0 0 0 | none | (none) |
| −1 1 1 90 90 90 | (−1, 1, 1, 90, 90, 90) | |
| nan 1 1 90 90 90 | (nan, 1, 1, 90, 90, 90) | |
| 1 1 1 0 0 0 | (1, 1, 1, 0, 0, 0) | |
| `1e1 2e1 3e1 9e1 9e1 9e1` (in their columns) | (10, 20, 30, 90, 90, 90) | |
| `1_0` 1 1 90 90 90 | (10, 1, 1, 90, 90, 90) | |
| only a, b, c present (line ends at col 33) | none | `Failed to read CRYST1 record, possibly invalid PDB file, got:\nCRYST1   10.000   20.000   30.000` |
| a field overflowing its columns (`12.3456789` in a 9-wide slot) | none | Failed to read … |

**Cell persistence between frames:** a frame without its own CRYST1 line keeps
the cell of the **previously read** frame. Examples:

| layout | frame 0 cell | frame 1 cell |
|---|---|---|
| CRYST1, MODEL 1, …, MODEL 2, … (1 CRYST1, 2 models) | none (the CRYST1 is before the first frame start) | none |
| MODEL 1, CRYST1(10…), …, MODEL 2, CRYST1(15…), … | 10,20,30 | 15,20,30 |
| CRYST1(10…), MODEL 1, …, CRYST1(15…), MODEL 2, … | 10,… | 15,… |
| MODEL 1, CRYST1(10…), …, MODEL 2, … (no CRYST1) | 10,20,30 | 10,20,30 (carried over) |

---------------------------------------------------------------------------

## 12. Atom-count consistency between the passes

The topology defines n_atoms (§1). In each frame, the coordinate pass counts
its atom lines (§10):

* fewer than n_atoms → ValueError:

      Inconsistency in file '<FILE>': The number of atoms (<count>) in trajectory frame <k> differs from the number of atoms (<n_atoms>) in the corresponding topology.
      Trajectories with varying numbers of atoms are currently not supported.

* more than n_atoms → IndexError
  `index <n_atoms> is out of bounds for axis 0 with size <n_atoms>`.

Frame 0 is read when the Universe is built, so these errors can come from
building it. Common triggers:

| situation | error |
|---|---|
| model 2 has fewer atoms than model 1 | ValueError (frame 1) when frame 1 is read |
| model 2 has more atoms | IndexError `index 1 is out of bounds for axis 0 with size 1` |
| ATOM lines after `END` (single model) | IndexError at load (the topology stopped at END) |
| `MODEL 1`, atom, `MODEL 2`, atom, no ENDMDL/END | ValueError frame 0: 1 vs 2 (the topology read both atoms) |
| an indented atom line | ValueError (the topology counts it after stripping; the coordinate pass does not) |
| `ATOM 123456 …` (6-digit serial) | ValueError (prefix `ATOM 1` is not `ATOM␣␣`) |
| ATOM lines before the first MODEL | ValueError |

Asking for a frame index ≥ n_frames through normal indexing raises
IndexError `Index 5 exceeds length of trajectory (2).` (numbers vary).

---------------------------------------------------------------------------

## 13. Metadata (HEADER, TITLE, COMPND, REMARK)

These are collected over the **whole file** when it is opened. The line is
taken including its terminator, and stripping removes surrounding whitespace
including CR/LF.

| record | test on line start | text kept | storage |
|---|---|---|---|
| HEADER | first 6 bytes `HEADER` | columns 11–66, stripped | single string; **the last HEADER line wins**; `""` if none |
| TITLE | first 5 bytes `TITLE` | columns 9–80, stripped | list, one entry per line, in order |
| COMPND | first 6 bytes `COMPND` | columns 8–80, stripped | list |
| REMARK | first 6 bytes `REMARK` | column 7 to end of line, stripped | list (an empty REMARK gives `""`) |

Examples:

| line | stored |
|---|---|
| `HEADER    PROTEIN                                 01-JAN-00   1ABC` then `HEADER    SECOND` | header = `SECOND` |
| `TITLE     FIRST LINE OF TITLE   ` | `FIRST LINE OF TITLE` |
| `TITLE    2 CONTINUATION` | `2 CONTINUATION` |
| `TITLEX extra` | `xtra` (it starts with `TITLE`; text from column 9) |
| `TITLE     ` + 70×`A` + `BBBBBB` | 70×`A` (cut at column 80) |
| `COMPND    MOL_ID: 1;` | `MOL_ID: 1;` |
| `COMPND   2 MOLECULE: X;` | `2 MOLECULE: X;` |
| `REMARK   2 RESOLUTION. 2.00 ANGSTROMS.` | `2 RESOLUTION. 2.00 ANGSTROMS.` |
| `REMARK` | `""` |
| `REMARKS hello` | `S hello` |
| CRLF file: `HEADER    PROTEIN  ` | `PROTEIN` |

---------------------------------------------------------------------------

## 14. Message catalogue

Warnings (UserWarning unless noted):

* `Serial numbers went over 100,000.  Higher serials have been guessed`
* `Invalid atom serials were present, bonds will not be parsed` (issued twice)
* `PDB file is missing resid information.  Defaulted to '1'` (per atom)
* `Element information is missing, elements attribute will not be populated. If needed these can be guessed using universe.guess_TopologyAttrs(context='default', to_guess=['elements']).`
* `Unknown element <E> found for some atoms. These have been given an empty element record. If needed they can be guessed using universe.guess_TopologyAttrs(context='default', to_guess=['elements']).` (per atom)
* `Unknown entry <Q> encountered in formal charge field. This likely indicates that the PDB file is not fully standard compliant. The formalcharges attribute will not be populated.`
* `Found CONECT record with single entry, ignoring this`
* `CONECT records was corrupt, bonds will not be parsed`
* `PDB file contained CONECT record to TER entry. These are not included in bonds.`
* `Failed to read CRYST1 record, possibly invalid PDB file, got:\n<line>`
* `1 A^3 CRYST1 record, this is usually a placeholder. Unit cell dimensions will be set to None.`
* `Reader has no dt information, set to 1.0 ps`
* PendingDeprecationWarning `Unknown masses are set to 0.0 for current version, this will be deprecated in version 3.0.0 and replaced by Masse's no_value_label (np.nan)`

Errors:

| trigger | error |
|---|---|
| no atom records | IndexError `index 0 is out of bounds for axis 0 with size 0` |
| malformed hybrid-36 serial | KeyError `'<char>'`, or UnboundLocalError `cannot access local variable 'serial' where it is not associated with a value` if it is the first atom |
| non-integer CONECT field | ValueError `invalid literal for int() with base 10: '<5 chars>'`. In a Universe it is wrapped: `Failed to construct topology from file <FILE> with parser <class 'MDAnalysis.topology.PDBParser.PDBParser'>.\nError: <message>` |
| bad coordinate field | ValueError `could not convert string to float: '<8 chars>'` (not wrapped) |
| frame has too few atoms | ValueError `Inconsistency in file '<FILE>': The number of atoms (<c>) in trajectory frame <k> differs from the number of atoms (<n>) in the corresponding topology.\nTrajectories with varying numbers of atoms are currently not supported.` |
| frame has too many atoms | IndexError `index <n> is out of bounds for axis 0 with size <n>` |
