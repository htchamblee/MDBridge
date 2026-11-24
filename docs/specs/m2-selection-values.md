# Selection language: value tokens, ranges, `prop`, and parser behaviour

Behaviour specification (clean-room) derived by observing MDAnalysis 2.11.0
(`select_atoms`). Everything below is stated as observable behaviour: what
string goes in, which atoms or which error come out. Error/warning texts are
quoted verbatim; they are part of the behaviour.

Conventions used in this document

* "token" = one whitespace-separated piece of the selection string after the
  tokenisation rule in §0.
* "END" = the end-of-input marker that the parser places after the last token.
  When it shows up in a message it is printed as the four characters `None`.
* In message texts, `<x>` stands for the offending text inserted verbatim.
* Exceptions are named with their MDAnalysis/Python class (SelectionError,
  ValueError, TypeError, IndexError, AttributeError, NoDataError, …). A port
  that has no such classes should still preserve the *message text* and the
  distinction "selection error raised while reading the string" vs. "error
  raised while evaluating against a topology" (§4.6).

## Example fixture

All tables were produced by running MDAnalysis on
`tests/data/topology/fixture_chains.pdb` (60 atoms, read as-is; no box, no
bonds, no charges). Results are listed as `N atoms: <atom serials>` (the PDB
serial numbers, which equal index+1). Useful facts about the fixture:

| serials | resname | resid (icode) | chain/segid | atom names |
|---|---|---|---|---|
| 1–8 | G | 1 | A | P OP1 C4' HO2' N9 N1 H1 O6 |
| 9–13 | C | 2 | A | P C4' N1 N3 N4 |
| 14–17 | A | 3 | A | P C4' N1 N6 |
| 18–19 | A | 3 (A) | A | P C4' |
| 20–23 | U | 4 | A | P C4* N3 H3 |
| 24–26 | G | 10 | B | P C4' N1 |
| 27–29 | C | 11 | B | P C4' N3 |
| 30–31 | C | 11 (A) | B | P C4' |
| 32–33 | U | 12 | B | P C4' |
| 34–40 | ALA | 5 | C | N CA C O CB H 1HB |
| 41–48 | PHE | 6 | C | N CA CG CD1 CD2 CE1 CE2 CZ |
| 49 / 50 / 51 | MG / Mg / MG2 | 201 / 202 / 203 | D | MG / Mg / MG |
| 52 / 53 / 54 | Na+ / CL / K | 204 / 205 / 206 | D | Na+ / CL / K |
| 55–57, 58–60 | HOH | 301, 302 | W | O H1 H2 |

x coordinate of atom *s* is 1.2·*s* Å; y cycles 0.3,0.6,0.9,1.2,0.0; z cycles
0.2,…,1.2,0.0. Masses are the element defaults (P 30.974, O 15.999,
C 12.011, H 1.008, N 14.007, Mg 24.305, Na 22.99, Cl 35.45, K 39.098);
tempfactor 0 and occupancy 1 for all atoms. Types: element symbols, upper-case
(`MG`), elements: `Mg`.

Some tables also use a second variant of the same topology with the insertion
code attribute removed ("no icodes").

---

## 0. Tokenisation (needed to understand everything else)

1. Every `(` and every `)` in the string is surrounded by spaces, then the
   string is split on runs of whitespace (any Unicode whitespace, including
   tab, newline and non-breaking space). Leading/trailing whitespace yields no
   tokens. Consequently a parenthesis is *always* its own token, even inside
   a would-be value: `name (P)` is the four tokens `name`, `(`, `P`, `)`, and
   `name \(` is `name`, `\`, `(`.
2. END is appended after the last token.
3. Keywords are case-sensitive. `NAME P` → `Unknown selection token: 'NAME'`;
   `AND` is not an operator (it is an ordinary value).
4. The empty string (exactly `""`) is special: no parsing happens, an empty
   result is returned and the warning
   `Empty string to select atoms, empty group returned.` is issued. A string
   consisting only of whitespace is *not* special: it fails with
   `Unknown selection token: 'None'`.
5. Several selection strings may be passed at once; all of them are parsed
   first (so a parse error in a later string is reported before any
   evaluation error of an earlier one), then each is evaluated and the
   results are concatenated in order (duplicates are kept; e.g. `'name P'`,
   `'resid 1'` gives 17 atoms with serial 1 appearing twice).
6. The result of one selection string is the set of matching atoms, unique
   and sorted by index.

---

## 1. (a) Collecting the value tokens after a value keyword

Applies to the string, integer, float, `resid`, boolean and single-character
value keywords (list in §1.3). It does **not** apply to `prop`, `around`,
`sphzone`, `sphlayer`, `isolayer`, `cyzone`, `cylayer`, `point`, `atom`,
`group`, `same`, `smarts` (those read a fixed number of raw tokens, see §4.4).

### 1.1 Rule

Starting right after the keyword, tokens are taken one by one, in order, as
long as the next token is **not a keyword**. Taking stops at the first
keyword (which is left for the parser). For each taken token, **one** leading
backslash, if present, is removed; nothing else is altered (no case folding,
no unquoting).

A token "is a keyword" exactly when it is one of:

* a selection keyword (full list in §1.3, case-sensitive; `chainID`/`chainid`
  and `altLoc`/`altloc` are both keywords),
* `and`, `or`,
* `(` or `)`,
* END.

Notes that follow from the rule:

* The keyword test is done on the token *before* the backslash is removed, so
  `\protein`, `\and`, `\resid`, `\all` are values (`protein`, `and`, `resid`,
  `all`). `\\P` loses only one backslash and is the value `\P`. A lone `\`
  becomes the empty value. A trailing backslash (`P\`) is kept.
* Words that only *look* like syntax are values: `AND`, `OR`, `NOT`, `to`,
  `:`, `-`, `foo`, `as`. So `name P foo` silently selects names `P` and `foo`.
* `not` and `global` and the distance keywords are keywords, so they stop
  value collection.
* Since parentheses are always tokens and always keywords, a value can never
  contain `(` or `)`.

If **zero** tokens were collected (the very next token is a keyword or END),
the keyword fails during parsing. The message depends on the keyword family
(note the colon difference):

| keyword family | message (inside `Selection failed: '…'`) |
|---|---|
| string keywords, `chirality` | `Unexpected token '<next token>'` |
| `resid`, integer and float keywords | `Unexpected token: '<next token>'` |
| boolean (`aromaticity`) | no error; behaves as if `true` had been given |

The full user-visible text is e.g. `Selection failed: 'Unexpected token 'and''`
or `Selection failed: 'Unexpected token: 'None''` (END prints as `None`).

### 1.2 What the collected values mean per family

* **String keywords**: each value is a case-sensitive shell-style pattern
  (`*` any run, `?` one char, `[seq]`, `[!seq]`); an atom matches if its
  attribute matches any value. (`resname \M*` → pattern `M*`.)
* **Boolean (`aromaticity`)**: each value must equal `true` or `false`
  ignoring case; otherwise `Selection failed: ''<value>' is an invalid value
  for boolean selection. Use 'True' or 'False''`. No values → `true`.
  Several values are OR-ed (`aromaticity true false` selects everything).
* **Single character (`chirality`)**: exact, case-sensitive equality with any
  value.
* **Integer, float and `resid`**: values are first joined into range pieces
  (§2), then each piece is converted (§2.3–§2.5).

### 1.3 Keyword inventory (as registered in 2.11.0)

* Logical: `and`, `or` (binary); `not`, `global` (prefix).
* Operand-taking: `byres`, `bonded`, `same … as`, `around`, `sphzone`,
  `sphlayer`, `isolayer`, `cyzone`, `cylayer`.
* Fixed-argument: `point`, `atom`, `group`, `prop`, `smarts`.
* No-argument: `all`, `protein`, `nucleic`, `water`, `backbone`,
  `nucleicbackbone`, `nucleicbase`, `nucleicsugar`, `aromatic`.
* String values: `name`, `type`, `element`, `record_type`, `chainID`,
  `chainid`, `altLoc`, `altloc`, `resname`, `icode`, `moltype`, `segid`.
* Integer values: `index`, `resindex`, `segindex`, `id`, `bynum` (compares
  index+1), `formalcharge`, `nbindex`, `resnum`, `molnum`, `model`.
* Float values: `radius`, `tempfactor`, `bfactor` (alias of tempfactor; each
  evaluation emits `The bfactor topology attribute is only provided as an
  alias to the tempfactor attribute. It will be removed in 3.0. Please use
  the tempfactor attribute instead.`), `mass`, `charge`, `occupancy`,
  `gbscreen`, `solventradius`, `rmin`, `epsilon`, `rmin14`, `epsilon14`.
* Special integer-with-insertion-code values: `resid`.
* Boolean: `aromaticity`. Single character: `chirality`.

The keyword set is global (it does not depend on which attributes the loaded
topology has); a keyword whose attribute is missing parses fine and fails at
evaluation with `This Universe does not contain <singular> information`
(NoDataError), e.g. `charge 0` → `This Universe does not contain charge
information`.

### 1.4 Examples (value collection)

| Selection string | MDAnalysis result |
|---|---|
| `name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name P C4'` | 17 atoms: 1,3,9,10,14,15,18-20,24,25,27,28,30-33 |
| `name P C4' and resid 1` | 2 atoms: 1,3 |
| `name P foo` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name P and C4'` | SelectionError: `Unknown selection token: 'C4''` |
| `name P name C4'` | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `resid 1 resname G` | SelectionError: `Unexpected token at end of selection string: 'resname'` |
| `protein P` | SelectionError: `Unexpected token at end of selection string: 'P'` |
| `name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `name and resid 1` | SelectionError: `Selection failed: 'Unexpected token 'and''` |
| `name protein` | SelectionError: `Selection failed: 'Unexpected token 'protein''` |
| `name \protein` | 0 atoms: (none) |
| `name resid` | SelectionError: `Selection failed: 'Unexpected token 'resid''` |
| `name \resid` | 0 atoms: (none) |
| `name all` | SelectionError: `Selection failed: 'Unexpected token 'all''` |
| `name \all` | 0 atoms: (none) |
| `name group` | SelectionError: `Selection failed: 'Unexpected token 'group''` |
| `name not` | SelectionError: `Selection failed: 'Unexpected token 'not''` |
| `name \not` | 0 atoms: (none) |
| `name \P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name \\P` | 0 atoms: (none) |
| `name P\` | 0 atoms: (none) |
| `name AND` | 0 atoms: (none) |
| `name And or name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name OR` | 0 atoms: (none) |
| `NAME P` | SelectionError: `Unknown selection token: 'NAME'` |
| `Name P` | SelectionError: `Unknown selection token: 'Name'` |
| `name to` | 0 atoms: (none) |
| `name :` | 0 atoms: (none) |
| `name -` | 0 atoms: (none) |
| `name H1 \and` | 3 atoms: 7,56,59 |
| `name P)` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `name P ( C4' )` | SelectionError: `Unexpected token at end of selection string: '('` |
| `name (P)` | SelectionError: `Selection failed: 'Unexpected token '(''` |
| `name \(` | SelectionError: `Unexpected token at end of selection string: '('` |
| `chainID A` | 23 atoms: 1-23 |
| `chainid A` | 23 atoms: 1-23 |
| `chainID B C` | 25 atoms: 24-48 |
| `resname Na+` | 1 atoms: 52 |
| `resname MG*` | 2 atoms: 49,51 |
| `resname \M*` | 3 atoms: 49-51 |
| `name C?'` | 8 atoms: 3,10,15,19,25,28,31,33 |
| `name [CN]?` | 16 atoms: 5,6,11-13,16,17,22,26,29,35,38,42,43,48,53 |
| `name [!C]*` | 40 atoms: 1,2,4-9,11-14,16-18,20,22-24,26,27,29,30,32,34,37,39-41,49-52,54-60 |
| `name Mg` | 1 atoms: 50 |
| `name mg` | 0 atoms: (none) |
| `element Mg` | 3 atoms: 49-51 |
| `element MG` | 0 atoms: (none) |
| `type MG` | 3 atoms: 49-51 |

---

## 2. (b) Range pieces for `resid` and the numeric keywords

### 2.1 Joining collected tokens into pieces

The collected values (after backslash removal, §1.1) are scanned left to
right and grouped into *pieces*; a piece is a string that may contain single
spaces. A token is a **separator token** when it is exactly `to`, `:` or `-`
(case-sensitive; `TO`, `To` are not). For each token, in order:

1. **Separator token.** It glues the current last piece, the separator and the
   *next* token together as `«last» «sep» «next»` (single spaces), consuming
   that next token whatever it is (even another separator). If there is no
   current piece, or no next token, this fails immediately with a
   SelectionError (not wrapped in "Selection failed") whose text is

   `Invalid expression given: ` + *L* + ` ` + *sep* + ` ` + *R*

   where *L* = all pieces formed so far joined by single spaces (empty if
   none) and *R* = all not-yet-consumed tokens joined by single spaces (empty
   if none; when there was no previous piece, *R* still includes the token
   that would have been the right-hand side). Note the resulting doubled or
   trailing spaces (shown with ␣ = one space character):

   | input | exact message |
   |---|---|
   | `resid : 3` | `Invalid␣expression␣given:␣␣:␣3` |
   | `resid : 3 4` | `Invalid␣expression␣given:␣␣:␣3␣4` |
   | `resid 1 :` | `Invalid␣expression␣given:␣1␣:␣` |
   | `resid 5 1 :` | `Invalid␣expression␣given:␣5␣1␣:␣` |
   | `resid 1 : 3 :` | `Invalid␣expression␣given:␣1␣:␣3␣:␣` |
   | `resid -` | `Invalid␣expression␣given:␣␣-␣` |
   | `resid to` | `Invalid␣expression␣given:␣␣to␣` |

   Only collected value tokens count: in `resid 1 : and name P` the `and`
   stops collection, so the `:` is dangling → same message as `resid 1 :`.
2. **Continuation token.** Otherwise, if a current piece exists and **any** of
   the following holds, the token is appended to the last piece with one
   space: the token begins with `--`; the token begins with `to` (so `to3`,
   `toto`, `to:3` qualify); the token begins with `:`; the last piece ends
   with `to`, `:` or `-`.
3. Otherwise the token starts a new piece.

Consequences: `1 -3` is two pieces (`1` and `-3`, a negative number), while
`1 --3` is one piece `1 --3` (range 1 to −3); `1- 3`, `1 :3`, `1: 3` are one
piece each; `1to 3` is one piece (`1to` ends with `to`) that later fails to
convert. `1 : 3 : 5` → `1 : 3 : 5` (one piece). Gluing happens before
conversion, so failures to convert report the whole piece
(`Failed to parse value: 1 : : 3`).

### 2.2 Range text syntax (shared idea)

After joining, each piece is either a single value or a range
"*lower* SEP *upper*" where SEP is one of: optional whitespace, `:`, optional
whitespace; optional whitespace, `-`, optional whitespace; or a space, `to`,
a space (exactly the letters `to` with whitespace on **both** sides; `1to3`,
`1 to3` do not contain a valid separator). The range text is matched as a
**prefix only**: anything after the upper bound is ignored (`1:3:5` ≡ `1:3`,
`1-3-5` ≡ `1-3`, `1 to 3 to 5` ≡ `1 to 3`, `1:3abc` ≡ `1:3` for integers,
`1:3.9` ≡ `1:3` for integers). A range is inclusive at both ends; if
lower > upper it selects nothing (except the `resid` + insertion-code case in
§2.4).

### 2.3 Integer keywords (`index`, `bynum`, `id`, `resnum`, `resindex`, …)

For each piece:

* First, the whole piece is tried as an integer literal with general integer
  syntax: optional surrounding whitespace, optional `+`/`-`, decimal digits,
  single underscores allowed between digits, leading zeros allowed
  (`+1`, `01`, `1_0` = 10 are accepted; `1.0`, `1e1`, `0x10` are not). If it
  converts, it is a single value (equality).
* Otherwise it must match, as a prefix, *int* SEP *int* where *int* is an
  optional `-` followed by one or more decimal digits (no `+`, no `_`). Then
  it is a range.
* Otherwise: `Selection failed: 'Failed to parse number: <piece>'`.

`bynum` compares against index+1, all others against the raw attribute.

### 2.4 `resid` (integer plus optional insertion code)

For each piece:

* **Single:** the entire piece is an optional `-`, one or more decimal digits,
  and at most one trailing "word character" (letter, digit or underscore,
  Unicode-aware) which is the insertion code (empty if absent). `3A` → resid 3,
  icode `A`; `+1`, `1.5`, `1AB`, `1,2`, `0x1` are rejected; `01` → 1.
* **Range:** otherwise, as a prefix: *resid-with-optional-icode* SEP
  *resid-with-optional-icode* (each side: optional `-`, digits, optional one
  word character). Trailing text after the upper side is ignored
  (`1:3abc` = `1` … `3a`; `1:3.5` = `1` … `3`).
* Otherwise: `Selection failed: 'Failed to parse value: <piece>'` (note
  "value" here vs "number" for the other numeric keywords).

Evaluation depends on whether the topology has insertion codes:

* **No insertion-code attribute:** if any parsed value or bound carries a
  non-empty insertion code, evaluation fails with ValueError
  `Selection specified icodes, while the topology doesn't have any.`
  Otherwise singles are equality on resid and ranges are inclusive
  `lower ≤ resid ≤ upper` (reversed range → nothing).
* **With insertion codes** (codes compare as strings; empty string sorts
  before every letter; `_` sorts after upper-case letters):
  * single `N` (no code) matches only resid N with **empty** code (so
    `resid 3` does *not* include 3A); `NX` matches resid N with code X.
  * range with equal resids `N[a]`…`N[b]`: resid N and a ≤ code ≤ b.
  * range with different resids `L[a]`…`U[b]`: (resid L and code ≥ a) or
    (L < resid < U, any code) or (resid U and code ≤ b). This is applied even
    when L > U, so reversed ranges are **not** empty: `resid 3:1` selects
    resid 3 (all codes) plus resid 1 with empty code, and `resid 1 --3`
    (range 1 to −3) selects resid 1.

### 2.5 Float keywords (`mass`, `tempfactor`, `occupancy`, `charge`, …)

For each piece:

* First, the whole piece is tried as a float literal with general float
  syntax: optional sign (`+`/`-`), digits with optional single underscores,
  optional fraction, optional exponent with `e` or `E`, or case-insensitive
  `inf`/`infinity`/`nan`. Accepted: `12`, `1_2.011`, `inf`, `nan`, `1E1`.
  A single value selects with tolerance (|v − value| ≤ 1e-8 + 1e-5·|value|)
  and **emits a warning at evaluation time, once per single value**:

  ```
  Using float equality to select atoms is not recommended because of inherent limitations in representing numbers on computers (see https://docs.python.org/3.8/tutorial/floatingpoint.html for more). Instead, we recommend using a range to select, e.g. 'mass 11.011 to 13.011'. If you still want to compare floats, use the `atol` and `rtol` keywords to modify the tolerance for `np.isclose`.
  ```

  where `mass` is the keyword as typed and the two numbers are value−1 and
  value+1 printed as Python prints floats (shortest round-trip repr, e.g.
  `11.0 to 13.0` for `12`, `0.008000000000000007 to 2.008` for `1.008`,
  `-14.0 to -12.0` for `-13`, `inf to inf`, `nan to nan`).
* Otherwise the piece is matched as a prefix against
  *f* SEP *f*, where *f* is: optional `-`, zero or more digits, optional `.`,
  zero or more digits, optional (`e`, optional `+`/`-`, one or more digits) —
  lower-case `e` only, no `+` sign on the mantissa, and *f* may be empty.
  If this does not match: `Selection failed: 'Failed to parse number: <piece>'`.
  If it matches, both captured texts are converted as floats; because *f* can
  be empty or just `-`/`.`, conversion can fail with
  `Selection failed: 'could not convert string to float: <repr>'`
  (examples: `:13`, `12:`, `-inf:inf`, `12:+13` all give `''`; `-:13` gives
  `'-'`; `.:13` gives `'.'`; `e1:13` gives `'e1'`). Ranges are inclusive
  and produce no warning.

`<repr>` is the Python representation of the text: wrapped in single quotes,
except that a text containing `'` and no `"` is wrapped in double quotes;
a backslash is shown doubled.

### 2.6 Examples: `resid` joining and parsing (fixture has insertion codes)

| Selection string | MDAnalysis result |
|---|---|
| `resid 1` | 8 atoms: 1-8 |
| `resid 1 2` | 13 atoms: 1-13 |
| `resid 1:3` | 17 atoms: 1-17 |
| `resid 1-3` | 17 atoms: 1-17 |
| `resid 1 to 3` | 17 atoms: 1-17 |
| `resid 1 : 3` | 17 atoms: 1-17 |
| `resid 1 - 3` | 17 atoms: 1-17 |
| `resid 1- 3` | 17 atoms: 1-17 |
| `resid 1 :3` | 17 atoms: 1-17 |
| `resid 1: 3` | 17 atoms: 1-17 |
| `resid 1 -3` | 8 atoms: 1-8 |
| `resid 1 --3` | 8 atoms: 1-8 |
| `resid 1--3` | 8 atoms: 1-8 |
| `resid 1 -- 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 -- 3'` |
| `resid 1to3` | SelectionError: `Selection failed: 'Failed to parse value: 1to3'` |
| `resid 1 to3` | SelectionError: `Selection failed: 'Failed to parse value: 1 to3'` |
| `resid 1to 3` | SelectionError: `Selection failed: 'Failed to parse value: 1to 3'` |
| `resid 1 TO 3` | SelectionError: `Selection failed: 'Failed to parse value: TO'` |
| `resid 1 To 3` | SelectionError: `Selection failed: 'Failed to parse value: To'` |
| `resid 1 toto 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 toto 3'` |
| `resid 1 to:3` | SelectionError: `Selection failed: 'Failed to parse value: 1 to:3'` |
| `resid 1 :to 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 :to 3'` |
| `resid 1:to 3` | SelectionError: `Selection failed: 'Failed to parse value: 1:to 3'` |
| `resid 1 \to 3` | 17 atoms: 1-17 |
| `resid 1 \: 3` | 17 atoms: 1-17 |
| `resid 1 : 3 10 to 11 201` | 24 atoms: 1-17,24-29,49 |
| `resid 1:3 10:11` | 23 atoms: 1-17,24-29 |
| `resid 1 2- 3` | 17 atoms: 1-17 |
| `resid 1 2 - 3 4` | 21 atoms: 1-17,20-23 |
| `resid 1 to 3 to 5` | 17 atoms: 1-17 |
| `resid 1 : 3 : 5` | 17 atoms: 1-17 |
| `resid 1:3:5` | 17 atoms: 1-17 |
| `resid 1-3-5` | 17 atoms: 1-17 |
| `resid 1:3 to 5` | 17 atoms: 1-17 |
| `resid 1 : : 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 : : 3'` |
| `resid 1 to to 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 to to 3'` |
| `resid 1 - - 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 - - 3'` |
| `resid 1 --- 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 --- 3'` |
| `resid 1 :: 3` | SelectionError: `Selection failed: 'Failed to parse value: 1 :: 3'` |
| `resid 1::3` | SelectionError: `Selection failed: 'Failed to parse value: 1::3'` |
| `resid 1: :3` | SelectionError: `Selection failed: 'Failed to parse value: 1: :3'` |
| `resid 1:` | SelectionError: `Selection failed: 'Failed to parse value: 1:'` |
| `resid :3` | SelectionError: `Selection failed: 'Failed to parse value: :3'` |
| `resid 1-` | SelectionError: `Selection failed: 'Failed to parse value: 1-'` |
| `resid 1 :` | SelectionError: `Invalid expression given: 1 : ` |
| `resid 1 to` | SelectionError: `Invalid expression given: 1 to ` |
| `resid 1 -` | SelectionError: `Invalid expression given: 1 - ` |
| `resid 5 1 :` | SelectionError: `Invalid expression given: 5 1 : ` |
| `resid 1 2 3 :` | SelectionError: `Invalid expression given: 1 2 3 : ` |
| `resid 1 : 3 :` | SelectionError: `Invalid expression given: 1 : 3 : ` |
| `resid : 3` | SelectionError: `Invalid expression given:  : 3` |
| `resid : 3 4` | SelectionError: `Invalid expression given:  : 3 4` |
| `resid to 3` | SelectionError: `Invalid expression given:  to 3` |
| `resid - 3` | SelectionError: `Invalid expression given:  - 3` |
| `resid :` | SelectionError: `Invalid expression given:  : ` |
| `resid to` | SelectionError: `Invalid expression given:  to ` |
| `resid -` | SelectionError: `Invalid expression given:  - ` |
| `resid 1 : and name P` | SelectionError: `Invalid expression given: 1 : ` |
| `resid 1 to and name P` | SelectionError: `Invalid expression given: 1 to ` |
| `resid 1 : 3 xyz` | SelectionError: `Selection failed: 'Failed to parse value: xyz'` |
| `resid 3 A` | SelectionError: `Selection failed: 'Failed to parse value: A'` |
| `resid -3` | 0 atoms: (none) |
| `resid -3:-1` | 0 atoms: (none) |
| `resid -3--1` | 0 atoms: (none) |
| `resid -3 - -1` | 0 atoms: (none) |
| `resid -3 to -1` | 0 atoms: (none) |
| `resid 1 : -3` | 8 atoms: 1-8 |
| `resid -3 : 1` | 8 atoms: 1-8 |
| `resid 3:1` | 14 atoms: 1-8,14-19 |
| `resid 01` | 8 atoms: 1-8 |
| `resid +1` | SelectionError: `Selection failed: 'Failed to parse value: +1'` |
| `resid 0x1` | SelectionError: `Selection failed: 'Failed to parse value: 0x1'` |
| `resid 1,2` | SelectionError: `Selection failed: 'Failed to parse value: 1,2'` |
| `resid 1.5` | SelectionError: `Selection failed: 'Failed to parse value: 1.5'` |
| `resid 1e2` | SelectionError: `Selection failed: 'Failed to parse value: 1e2'` |
| `resid abc` | SelectionError: `Selection failed: 'Failed to parse value: abc'` |
| `resid 1:abc` | SelectionError: `Selection failed: 'Failed to parse value: 1:abc'` |
| `resid 1AB` | SelectionError: `Selection failed: 'Failed to parse value: 1AB'` |
| `resid 1:3abc` | 19 atoms: 1-19 |
| `resid 1:3.5` | 17 atoms: 1-17 |
| `resid 1 to 3abc` | 19 atoms: 1-19 |
| `resid 1 \and` | SelectionError: `Selection failed: 'Failed to parse value: and'` |
| `resid \` | SelectionError: `Selection failed: 'Failed to parse value: '` |
| `resid` | SelectionError: `Selection failed: 'Unexpected token: 'None''` |
| `resid and name P` | SelectionError: `Selection failed: 'Unexpected token: 'and''` |

### 2.7 Examples: `resid` with and without insertion codes

| Selection string | with insertion codes (fixture as read) | same topology after removing the insertion-code attribute |
|---|---|---|
| `resid 3` | 4 atoms: 14-17 | 6 atoms: 14-19 |
| `resid 3A` | 2 atoms: 18,19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3a` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3B` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3 3A` | 6 atoms: 14-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3:3A` | 6 atoms: 14-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3-3A` | 6 atoms: 14-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3A:3A` | 2 atoms: 18,19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3A-3` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3A:4` | 6 atoms: 18-23 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3A to 4` | 6 atoms: 18-23 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 2-3A` | 11 atoms: 9-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 2A:4` | 10 atoms: 14-23 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 2:4A` | 15 atoms: 9-23 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11A` | 2 atoms: 30,31 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 10:11A` | 8 atoms: 24-31 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11A:12` | 4 atoms: 30-33 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 1x:3y` | 11 atoms: 9-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 1_` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 1:3_` | 19 atoms: 1-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 1:3abc` | 19 atoms: 1-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 1 --3` | 8 atoms: 1-8 | 0 atoms: (none) |
| `resid 3:1` | 14 atoms: 1-8,14-19 | 0 atoms: (none) |
| `resid 1 : -3` | 8 atoms: 1-8 | 0 atoms: (none) |
| `resid -3:1` | 8 atoms: 1-8 | 8 atoms: 1-8 |

### 2.8 Examples: integer keywords

| Selection string | MDAnalysis result |
|---|---|
| `index 0` | 1 atoms: 1 |
| `index 0 1 2` | 3 atoms: 1-3 |
| `index 0:2` | 3 atoms: 1-3 |
| `index 0-2` | 3 atoms: 1-3 |
| `index 0 to 2` | 3 atoms: 1-3 |
| `index 0 : 2` | 3 atoms: 1-3 |
| `index 0 - 2` | 3 atoms: 1-3 |
| `index 0 -2` | 1 atoms: 1 |
| `index 0 --2` | 0 atoms: (none) |
| `index 0--2` | 0 atoms: (none) |
| `index 1--3` | 0 atoms: (none) |
| `index -3--1` | 0 atoms: (none) |
| `index -3 - -1` | 0 atoms: (none) |
| `index 2:0` | 0 atoms: (none) |
| `index -1` | 0 atoms: (none) |
| `index +1` | 1 atoms: 2 |
| `index 01` | 1 atoms: 2 |
| `index 1_0` | 1 atoms: 11 |
| `index 1.0` | SelectionError: `Selection failed: 'Failed to parse number: 1.0'` |
| `index 1e1` | SelectionError: `Selection failed: 'Failed to parse number: 1e1'` |
| `index 0x10` | SelectionError: `Selection failed: 'Failed to parse number: 0x10'` |
| `index a` | SelectionError: `Selection failed: 'Failed to parse number: a'` |
| `index 1a` | SelectionError: `Selection failed: 'Failed to parse number: 1a'` |
| `index 1:a` | SelectionError: `Selection failed: 'Failed to parse number: 1:a'` |
| `index 1:3abc` | 3 atoms: 2-4 |
| `index 1:3:5` | 3 atoms: 2-4 |
| `index 1:3.9` | 3 atoms: 2-4 |
| `index 1 to 3 to 5` | 3 atoms: 2-4 |
| `index 0:` | SelectionError: `Selection failed: 'Failed to parse number: 0:'` |
| `index :2` | SelectionError: `Selection failed: 'Failed to parse number: :2'` |
| `index 0 :` | SelectionError: `Invalid expression given: 0 : ` |
| `index : 2` | SelectionError: `Invalid expression given:  : 2` |
| `index -` | SelectionError: `Invalid expression given:  - ` |
| `index 1 -- 3` | SelectionError: `Selection failed: 'Failed to parse number: 1 -- 3'` |
| `index 1 to3` | SelectionError: `Selection failed: 'Failed to parse number: 1 to3'` |
| `index 1to3` | SelectionError: `Selection failed: 'Failed to parse number: 1to3'` |
| `index 1 TO 3` | SelectionError: `Selection failed: 'Failed to parse number: TO'` |
| `bynum 1` | 1 atoms: 1 |
| `bynum 1:3` | 3 atoms: 1-3 |
| `bynum 0` | 0 atoms: (none) |
| `bynum 61` | 0 atoms: (none) |
| `bynum 58 to 70` | 3 atoms: 58-60 |
| `id 1 - 3` | 3 atoms: 1-3 |
| `resnum 1:3` | 19 atoms: 1-19 |
| `resnum 1 3` | 14 atoms: 1-8,14-19 |
| `resnum 3A` | SelectionError: `Selection failed: 'Failed to parse number: 3A'` |

### 2.9 Examples: float keywords

| Selection string | MDAnalysis result |
|---|---|
| `mass 12.011` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 (+float-equality warning) |
| `mass 30.974` | 9 atoms: 1,9,14,18,20,24,27,30,32 (+float-equality warning) |
| `mass 12.011 1.008` | 28 atoms: 3,4,7,10,15,19,21,23,25,28,31,33,35,36,38-40,42-48,56,57,59,60 (+2× float-equality warning) |
| `mass 1_2.011` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 (+float-equality warning) |
| `mass 12 to 13` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass 12-13` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass 12:13` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass 12 - 13` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass 12 :13` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass 12 -13` | 0 atoms: (none) (+2× float-equality warning) |
| `mass 12--13` | 0 atoms: (none) |
| `mass -1--13` | 0 atoms: (none) |
| `mass 13:12` | 0 atoms: (none) |
| `mass 1e1-2e1` | 36 atoms: 2,3,5,6,8,10-13,15-17,19,21,22,25,26,28,29,31,33-38,41-48,55,58 |
| `mass 1e1 to 2e1` | 36 atoms: 2,3,5,6,8,10-13,15-17,19,21,22,25,26,28,29,31,33-38,41-48,55,58 |
| `mass 1e+1:2e+1` | 36 atoms: 2,3,5,6,8,10-13,15-17,19,21,22,25,26,28,29,31,33-38,41-48,55,58 |
| `mass 1.0e1:1.3e1` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass 1E1-2E1` | SelectionError: `Selection failed: 'Failed to parse number: 1E1-2E1'` |
| `mass 10.:13.` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass .5:13` | 28 atoms: 3,4,7,10,15,19,21,23,25,28,31,33,35,36,38-40,42-48,56,57,59,60 |
| `mass -1:13` | 28 atoms: 3,4,7,10,15,19,21,23,25,28,31,33,35,36,38-40,42-48,56,57,59,60 |
| `mass 12:13abc` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 |
| `mass +12:13` | SelectionError: `Selection failed: 'Failed to parse number: +12:13'` |
| `mass 12:+13` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `mass -inf:inf` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `mass inf` | 0 atoms: (none) (+float-equality warning) |
| `mass nan` | 0 atoms: (none) (+float-equality warning) |
| `mass :13` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `mass 12:` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `mass -:13` | SelectionError: `Selection failed: 'could not convert string to float: '-''` |
| `mass e1:13` | SelectionError: `Selection failed: 'could not convert string to float: 'e1''` |
| `mass .:13` | SelectionError: `Selection failed: 'could not convert string to float: '.''` |
| `mass 12 :` | SelectionError: `Invalid expression given: 12 : ` |
| `mass : 13` | SelectionError: `Invalid expression given:  : 13` |
| `mass 12 to13` | SelectionError: `Selection failed: 'Failed to parse number: 12 to13'` |
| `mass abc` | SelectionError: `Selection failed: 'Failed to parse number: abc'` |
| `mass 30.974 to 30.974` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `tempfactor 0` | 60 atoms: 1-60 (+float-equality warning) |
| `tempfactor 0:0` | 60 atoms: 1-60 |
| `occupancy 0.5:1.5` | 60 atoms: 1-60 |
| `charge 0` | NoDataError: `This Universe does not contain charge information` |
| `formalcharge 0` | NoDataError: `This Universe does not contain formalcharge information` |
| `radius 1` | NoDataError: `This Universe does not contain radius information` |

---

## 3. (c) The `prop` keyword

### 3.1 Accepted forms

`prop [abs] PROPERTY OP VALUE`, where OP is one of `<`, `>`, `<=`, `>=`,
`==`, `!=`, with or without whitespace on either side of OP
(`prop x<5`, `prop x< 5`, `prop x <5`, `prop x < 5` are equivalent), and
the reversed order `prop [abs] VALUE OP PROPERTY` (e.g. `prop 5 > x`).
`5 > prop x` is **not** accepted (`Unknown selection token: '5'`).

* `abs` must be a separate, lower-case token directly after `prop`
  (`prop absx`, `prop ABS x` do not work); it makes the comparison use the
  absolute value of the attribute.
* `prop` does **not** use the value-collection rule of §1: it reads raw
  tokens (keywords, parentheses and END can be swallowed), and backslashes
  are not removed.
* Comparison semantics: `<`, `>`, `<=`, `>=`, `!=` are exact comparisons;
  `==` is the tolerance test |v − value| ≤ 1e-8 + 1e-5·|value|. When the
  attribute is stored in single precision (coordinates x, y, z) the value is
  rounded to single precision before comparing, so `prop x < 6.0000001`
  excludes an atom at x = 6.0 and `prop x <= 5.9999999` includes it;
  double-precision attributes such as mass compare exactly.
* x, y, z are the Cartesian coordinate components; all other properties are
  the per-atom value of the named attribute (residue/segment attributes are
  broadcast to atoms; `index`, `resindex`, `segindex` are 0-based).

### 3.2 How the three parts are recovered (observable rules)

The six operator symbols are always considered in this priority order:
`<=`, `>=`, `==`, `!=`, `<`, `>`.

1. **abs**: if the first token after `prop` is exactly `abs`, absolute mode is
   on and the property text is the following token; otherwise the first
   token is the property text.
2. **Operator embedded in the property token**: going through the symbols in
   priority order, the first symbol that occurs **exactly once** in the
   property token splits it: the text before becomes the property name, and
   the symbol plus the text after becomes the *operator text*. A symbol that
   occurs twice or more is skipped (the next symbol is tried). If no symbol
   splits the token, the *operator text* is the next raw token.
3. **Value embedded in the operator text**: find the first symbol (priority
   order) that occurs anywhere in the operator text.
   * None occurs → the operator text is taken verbatim as the operator and the
     value is the next raw token.
   * It occurs more than once → error `Selection failed: 'too many values to
     unpack (expected 2)'`.
   * It occurs once and something follows it → operator = that symbol,
     value = the text after it. Whatever preceded the symbol is silently
     discarded (`prop x ><5` means `x < 5`).
   * It occurs once and nothing follows it → the whole operator text
     (including anything before the symbol) is the operator, and the value is
     the next raw token.
4. **Reversed operands**: if the value text is exactly one of the registered
   property names (§3.3, which includes END), property and value are swapped
   and the operator is replaced: `<`→`>=`, `>`→`<=`, `<=`→`>`, `>=`→`<`,
   `==`→`==`, `!=`→`!=`. **This replacement is the logical negation, not the
   mirror**, so the reversed form behaves inclusively/exclusively "wrongly":
   `prop 6 > x` means x ≤ 6, `prop 6 < x` means x ≥ 6, `prop 6 >= x` means
   x < 6, `prop 6 <= x` means x > 6. If the operator text at this point is
   not one of the six symbols, the error is
   `Unknown selection token: 'prop'` (SelectionError, no "Selection failed").
5. **Operator check**: if the operator is not one of the six symbols:
   `Selection failed: 'Invalid operator : '<op>' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''`
   (this check happens after the value token was consumed and after step 4,
   but before the value is converted).
6. **Value**: converted with general float syntax (as in §2.5); failure →
   `Selection failed: 'could not convert string to float: <repr>'`.
7. **Property** (checked only at evaluation, i.e. after the whole string
   parsed successfully):
   * not a registered property name → SelectionError
     `Expected one of [<list>]` (full list in §3.3);
   * registered but absent from the topology → SelectionError
     `This Universe does not contain <plural> information` where `<plural>`
     is the attribute's plural name (§3.3), e.g. `charges`, `formalcharges`,
     `radii`, `bonds`, `topologyattrs`;
   * registered and present but text-valued (name, resname, chainID,
     element, segid, icode, …): `<`, `>`, `<=`, `>=` raise TypeError
     `'<' not supported between instances of 'str' and 'float'` (with the
     operator used), `==` raises TypeError
     `unsupported operand type(s) for -: 'str' and 'float'`, `!=` selects
     every atom, and `abs` raises TypeError `bad operand type for abs(): 'str'`.
   * Boolean attributes compare as 0/1 (`prop aromaticity == 1` works).

**END inside `prop`** (because `prop` reads raw tokens):

* property token is END (`prop`, `prop abs`) → AttributeError
  `'NoneType' object has no attribute 'split'`.
* operator text would be END (`prop x`, `prop x5`, `prop x<<5`, `prop x<5<6`,
  `prop x==5==`) → TypeError `argument of type 'NoneType' is not iterable`.
* value would be END → step 4 swaps (END is a registered name): with a valid
  operator the former property text becomes the value
  (`prop x <` → `Selection failed: 'could not convert string to float: 'x''`;
  `prop <` → `… float: '''`); with an invalid operator →
  `Unknown selection token: 'prop'` (`prop x 5`, `prop x !5`). If the former
  property text *is* numeric (`prop 5 <`), the selection constructs and the
  parser then fails with IndexError `deque index out of range`.

### 3.3 Registered property names and their plural attribute names

Valid property names (exactly this list, in this order, is printed by the
"Expected one of" error; note the bare `None`):

```
Expected one of ['x', 'y', 'z', 'topologyattr', 'index', 'resindex', 'segindex', 'atomattr', 'id', None, 'name', 'type', 'element', 'radius', 'record_type', 'chainID', 'tempfactor', 'mass', 'charge', 'formalcharge', 'occupancy', 'altLoc', 'gbscreen', 'solventradius', 'nbindex', 'rmin', 'epsilon', 'rmin14', 'epsilon14', 'aromaticity', 'chirality', 'residueattr', 'resid', 'resname', 'resnum', 'icode', 'moltype', 'molnum', 'segmentattr', 'segid', 'bonds', 'ureybradleys', 'angles', 'dihedrals', 'impropers', 'cmaps', 'atomnum', 'atomiccharge', 'chargegroup', 'model', 'type_index']
```

Plural names used in "does not contain … information": x/y/z → positions,
topologyattr → topologyattrs, index → indices, resindex → resindices,
segindex → segindices, atomattr → atomattrs, id → ids, name → names,
type → types, element → elements, radius → radii, record_type → record_types,
chainID → chainIDs, tempfactor → tempfactors, mass → masses,
charge → charges, formalcharge → formalcharges, occupancy → occupancies,
altLoc → altLocs, gbscreen → gbscreens, solventradius → solventradii,
nbindex → nbindices, rmin → rmins, epsilon → epsilons, rmin14 → rmin14s,
epsilon14 → epsilon14s, aromaticity → aromaticities, chirality → chiralities,
residueattr → residueattrs, resid → resids, resname → resnames,
resnum → resnums, icode → icodes, moltype → moltypes, molnum → molnums,
segmentattr → segmentattrs, segid → segids, bonds → bonds,
ureybradleys → ureybradleys, angles → angles, dihedrals → dihedrals,
impropers → impropers, cmaps → cmaps, atomnum → atomnums,
atomiccharge → atomiccharges, chargegroup → chargegroups, model → models,
type_index → type_indices. Plural spellings themselves (`names`, `masses`)
are **not** valid property names. `bfactor` is not a property name.

### 3.4 Examples: valid `prop` forms

| Selection string | MDAnalysis result |
|---|---|
| `prop x < 5` | 4 atoms: 1-4 |
| `prop x< 5` | 4 atoms: 1-4 |
| `prop x <5` | 4 atoms: 1-4 |
| `prop x<5` | 4 atoms: 1-4 |
| `prop x > 70` | 2 atoms: 59,60 |
| `prop x>70` | 2 atoms: 59,60 |
| `prop x >= 69.6` | 3 atoms: 58-60 |
| `prop x>=69.6` | 3 atoms: 58-60 |
| `prop x <= 2.4` | 2 atoms: 1,2 |
| `prop x<=2.4` | 2 atoms: 1,2 |
| `prop x == 6` | 1 atoms: 5 |
| `prop x==6` | 1 atoms: 5 |
| `prop x == 6.0000001` | 1 atoms: 5 |
| `prop x == 6.00005` | 1 atoms: 5 |
| `prop x == 6.00007` | 0 atoms: (none) |
| `prop x != 1.2` | 59 atoms: 2-60 |
| `prop x!=1.2` | 59 atoms: 2-60 |
| `prop x!= 5` | 60 atoms: 1-60 |
| `prop x !=5` | 60 atoms: 1-60 |
| `prop x < 6.0000001` | 4 atoms: 1-4 |
| `prop x <= 5.9999999` | 5 atoms: 1-5 |
| `prop x != 6.0000001` | 59 atoms: 1-4,6-60 |
| `prop y == 0` | 12 atoms: 5,10,15,20,25,30,35,40,45,50,55,60 |
| `prop z > 1.1` | 8 atoms: 6,13,20,27,34,41,48,55 |
| `prop abs z > 1.1` | 8 atoms: 6,13,20,27,34,41,48,55 |
| `prop abs x < 3` | 2 atoms: 1,2 |
| `prop abs x<3` | 2 atoms: 1,2 |
| `prop abs y>=1.2` | 12 atoms: 4,9,14,19,24,29,34,39,44,49,54,59 |
| `prop absx < 5` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop ABS x < 5` | SelectionError: `Selection failed: 'Invalid operator : 'x' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop mass > 30` | 11 atoms: 1,9,14,18,20,24,27,30,32,53,54 |
| `prop mass>30` | 11 atoms: 1,9,14,18,20,24,27,30,32,53,54 |
| `prop mass == 30.974` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `prop mass < 30.9740001` | 58 atoms: 1-52,55-60 |
| `prop mass <= 30.97399999` | 49 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-52,55-60 |
| `prop tempfactor == 0` | 60 atoms: 1-60 |
| `prop occupancy < 1` | 0 atoms: (none) |
| `prop resid < 3` | 13 atoms: 1-13 |
| `prop resnum < 3` | 13 atoms: 1-13 |
| `prop index < 3` | 3 atoms: 1-3 |
| `prop id < 3` | 2 atoms: 1,2 |
| `prop resindex == 1` | 5 atoms: 9-13 |
| `prop segindex == 1` | 10 atoms: 24-33 |
| `prop x < 1e1` | 8 atoms: 1-8 |
| `prop x < 1E1` | 8 atoms: 1-8 |
| `prop x < .5e1` | 4 atoms: 1-4 |
| `prop x < 1_0` | 8 atoms: 1-8 |
| `prop x < inf` | 60 atoms: 1-60 |
| `prop x < Infinity` | 60 atoms: 1-60 |
| `prop x < nan` | 0 atoms: (none) |
| `prop x < -inf` | 0 atoms: (none) |
| `prop x < -5` | 0 atoms: (none) |
| `prop x<-5` | 0 atoms: (none) |
| `prop x>-1` | 60 atoms: 1-60 |
| `prop x < 5 and name P` | 1 atoms: 1 |
| `(prop x < 5)` | 4 atoms: 1-4 |

### 3.5 Examples: reversed operand order

| Selection string | MDAnalysis result |
|---|---|
| `prop 5 > x` | 4 atoms: 1-4 |
| `prop 6 > x` | 5 atoms: 1-5 |
| `prop x < 6` | 4 atoms: 1-4 |
| `prop 6 < x` | 56 atoms: 5-60 |
| `prop x > 6` | 55 atoms: 6-60 |
| `prop 6 >= x` | 4 atoms: 1-4 |
| `prop 6 <= x` | 55 atoms: 6-60 |
| `prop x <= 6` | 5 atoms: 1-5 |
| `prop 70 < x` | 2 atoms: 59,60 |
| `prop 70<x` | 2 atoms: 59,60 |
| `prop 6 == x` | 1 atoms: 5 |
| `prop 1.2 != x` | 59 atoms: 2-60 |
| `prop 30 < mass` | 11 atoms: 1,9,14,18,20,24,27,30,32,53,54 |
| `prop 30<mass` | 11 atoms: 1,9,14,18,20,24,27,30,32,53,54 |
| `prop 2 > resid` | 13 atoms: 1-13 |
| `prop abs 1.1 < z` | 8 atoms: 6,13,20,27,34,41,48,55 |
| `5 > prop x` | SelectionError: `Unknown selection token: '5'` |
| `prop x < y` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop y < x` | SelectionError: `Selection failed: 'could not convert string to float: 'y''` |
| `prop x < x` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop mass<x` | SelectionError: `Selection failed: 'could not convert string to float: 'mass''` |
| `prop x<mass` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop 5 < 6` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop 5<6` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |

### 3.6 Examples: malformed `prop` input (every distinct message)

| Selection string | MDAnalysis result |
|---|---|
| `prop` | AttributeError: `'NoneType' object has no attribute 'split'` |
| `prop abs` | AttributeError: `'NoneType' object has no attribute 'split'` |
| `prop x` | TypeError: `argument of type 'NoneType' is not iterable` |
| `prop abs x` | TypeError: `argument of type 'NoneType' is not iterable` |
| `prop x5` | TypeError: `argument of type 'NoneType' is not iterable` |
| `prop 5` | TypeError: `argument of type 'NoneType' is not iterable` |
| `prop x<<5` | TypeError: `argument of type 'NoneType' is not iterable` |
| `prop x<5<6` | TypeError: `argument of type 'NoneType' is not iterable` |
| `prop x==5==` | TypeError: `argument of type 'NoneType' is not iterable` |
| `prop x <` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop abs x <` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop x!=` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop x ==` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop <` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `prop <=` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `prop ==` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `prop x 5` | SelectionError: `Unknown selection token: 'prop'` |
| `prop x !5` | SelectionError: `Unknown selection token: 'prop'` |
| `prop 5 = mass` | SelectionError: `Unknown selection token: 'prop'` |
| `prop x = 5` | SelectionError: `Selection failed: 'Invalid operator : '=' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x = abc` | SelectionError: `Selection failed: 'Invalid operator : '=' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x => 5` | SelectionError: `Selection failed: 'Invalid operator : '=>' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x =< 5` | SelectionError: `Selection failed: 'Invalid operator : '=<' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x ! 5` | SelectionError: `Selection failed: 'Invalid operator : '!' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x - 5` | SelectionError: `Selection failed: 'Invalid operator : '-' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x : 5` | SelectionError: `Selection failed: 'Invalid operator : ':' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x to 5` | SelectionError: `Selection failed: 'Invalid operator : 'to' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x 5 <` | SelectionError: `Selection failed: 'Invalid operator : '5' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x \< 5` | SelectionError: `Selection failed: 'Invalid operator : '\<' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop abs abs x < 1` | SelectionError: `Selection failed: 'Invalid operator : 'x' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop ( x < 5 )` | SelectionError: `Selection failed: 'Invalid operator : 'x' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop abs(x) < 5` | SelectionError: `Selection failed: 'Invalid operator : 'x' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` |
| `prop x <> 5` | SelectionError: `Selection failed: 'could not convert string to float: '>''` |
| `prop x ><5` | 4 atoms: 1-4 |
| `prop x === 5` | SelectionError: `Selection failed: 'could not convert string to float: '=''` |
| `prop x << 5` | SelectionError: `Selection failed: 'too many values to unpack (expected 2)'` |
| `prop x <<5` | SelectionError: `Selection failed: 'too many values to unpack (expected 2)'` |
| `prop x <5<6` | SelectionError: `Selection failed: 'too many values to unpack (expected 2)'` |
| `prop x ==5==6` | SelectionError: `Selection failed: 'too many values to unpack (expected 2)'` |
| `prop x < < 5` | SelectionError: `Selection failed: 'could not convert string to float: '<''` |
| `prop x > = 5` | SelectionError: `Selection failed: 'could not convert string to float: '=''` |
| `prop x >=< 5` | SelectionError: `Selection failed: 'could not convert string to float: '<''` |
| `prop x == == 5` | SelectionError: `Selection failed: 'could not convert string to float: '==''` |
| `prop x <= <= 5` | SelectionError: `Selection failed: 'could not convert string to float: '<=''` |
| `prop x != != 5` | SelectionError: `Selection failed: 'could not convert string to float: '!=''` |
| `prop x< =5` | SelectionError: `Selection failed: 'could not convert string to float: '=5''` |
| `prop x <= =5` | SelectionError: `Selection failed: 'could not convert string to float: '=5''` |
| `prop x <=5=` | SelectionError: `Selection failed: 'could not convert string to float: '5=''` |
| `prop x<5>6` | SelectionError: `Selection failed: 'could not convert string to float: '5>6''` |
| `prop x< 5<6` | SelectionError: `Selection failed: 'could not convert string to float: '5<6''` |
| `prop x < abc` | SelectionError: `Selection failed: 'could not convert string to float: 'abc''` |
| `prop x < 5abc` | SelectionError: `Selection failed: 'could not convert string to float: '5abc''` |
| `prop x < \5` | SelectionError: `Selection failed: 'could not convert string to float: '\\5''` |
| `prop x < 5'` | SelectionError: `Selection failed: 'could not convert string to float: "5'"'` |
| `prop x < 5"` | SelectionError: `Selection failed: 'could not convert string to float: '5"''` |
| `prop x < (5)` | SelectionError: `Selection failed: 'could not convert string to float: '(''` |
| `prop x < and name P` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `prop x < and` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `prop name<P` | SelectionError: `Selection failed: 'could not convert string to float: 'P''` |
| `prop x<y` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `prop < x 5` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `prop x5 <` | SelectionError: `Selection failed: 'could not convert string to float: 'x5''` |
| `prop x < 5 6` | SelectionError: `Unexpected token at end of selection string: '6'` |
| `prop x <5 <6` | SelectionError: `Unexpected token at end of selection string: '<6'` |
| `prop <5 x` | SelectionError: `Unexpected token at end of selection string: 'x'` |
| `prop x<5)` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `prop foo > 1` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop X < 5` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop \x < 5` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop names == 1` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop masses > 30` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop abs < 5` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop and < 5` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop None > 1` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop charge > 0` | SelectionError: `This Universe does not contain charges information` |
| `prop formalcharge == 0` | SelectionError: `This Universe does not contain formalcharges information` |
| `prop radius > 0` | SelectionError: `This Universe does not contain radii information` |
| `prop bonds > 1` | SelectionError: `This Universe does not contain bonds information` |
| `prop topologyattr > 1` | SelectionError: `This Universe does not contain topologyattrs information` |
| `prop name == 1` | TypeError: `unsupported operand type(s) for -: 'str' and 'float'` |
| `prop name != 1` | 60 atoms: 1-60 |
| `prop name < 5` | TypeError: `'<' not supported between instances of 'str' and 'float'` |
| `prop resname > 1` | TypeError: `'>' not supported between instances of 'str' and 'float'` |
| `prop chainID == 1` | TypeError: `unsupported operand type(s) for -: 'str' and 'float'` |
| `prop chainID != 1` | 60 atoms: 1-60 |
| `prop abs name > 1` | TypeError: `bad operand type for abs(): 'str'` |

Additional END cases:

| Selection string | MDAnalysis result |
|---|---|
| `prop 5 <` | IndexError: `deque index out of range` |
| `prop 5 ==` | IndexError: `deque index out of range` |
| `prop abs 5 <` | IndexError: `deque index out of range` |
| `prop 5<` | IndexError: `deque index out of range` |
| `prop -1 <` | IndexError: `deque index out of range` |
| `prop x>` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |

(`…PROPERTY-LIST…` abbreviates the list printed in §3.3.)

---

## 4. (d) Parser structure, precedence and errors

### 4.1 Grammar in words

* An *expression* is a sequence of *terms* joined by `and`/`or`.
* `and` and `or` have the **same** precedence and associate to the **left**:
  `A or B and C` = `(A or B) and C`; `A and B or C` = `(A and B) or C`.
* A *term* is one of: `( expression )`; a prefix operator `not`/`global`
  followed by a *term-level* operand; an operand-taking keyword
  (`byres`, `bonded`, `same <prop> as`, `around R`, `sphzone R`,
  `sphlayer R1 R2`, `isolayer R1 R2`, `cyzone R ZMAX ZMIN`,
  `cylayer R1 R2 ZMAX ZMIN`) followed by its operand; or a simple selection
  (value keyword with values, no-argument keyword, fixed-argument keyword).
* Operand binding strengths:
  * `not` and `global` take the tightest operand: a single term. They never
    swallow a following `and`/`or`. `not A and B` = `(not A) and B`.
  * `byres`, `bonded`, `same … as`, `around`, `sphzone`, `sphlayer`,
    `isolayer`, `cyzone`, `cylayer` take the **loosest** operand: everything
    up to the end of input or the closing parenthesis of the enclosing group,
    including all `and`/`or`. `byres A or B` = `byres (A or B)`;
    `around 2 A and B` = `around 2 (A and B)`. This holds even when the
    keyword itself sits in a tight position: `not byres A and B` =
    `not (byres (A and B))`, `X and around 1.3 A or B` =
    `X and (around 1.3 (A or B))`. Parentheses limit the swallowing:
    `(around 1.3 A) or B`.
  * The right-hand side of `and`/`or` is a single term (with the swallowing
    exception just described), which produces the left associativity.
* After a complete expression the next token must be END.

### 4.2 Evaluation order (visible only through which error wins)

All strings are parsed before evaluation. During evaluation, for `and` the
**right** operand is evaluated before the left one; for `or` the **left**
operand is evaluated first. So `charge 0 and prop foo > 1` reports the
`prop` error while `prop foo > 1 and charge 0` reports the missing charges.

### 4.3 Parse-time error messages

| situation | message (SelectionError unless stated) |
|---|---|
| a term is expected but the token is not a keyword, or is `and`/`or`/`)`, or END | `Unknown selection token: '<token>'` (END → `'None'`) |
| expression complete but tokens remain | `Unexpected token at end of selection string: '<token>'` |
| `(` without matching `)` | `Unexpected token: '<token>' Expected: ')'` |
| `same <prop>` not followed by `as` | `Unexpected token: '<token>' Expected: 'as'` |
| a keyword's own parsing fails with a "value" problem | `Selection failed: '<inner message>'` |
| a keyword's own parsing hits an unknown-name lookup (only in `prop`, §3.2) | `Unknown selection token: '<keyword>'` |

"Value problems" that get wrapped into `Selection failed: '…'` include: no
values after a value keyword; `Failed to parse value: …`;
`Failed to parse number: …`; boolean value errors; float/int conversion
failures of distance/`atom` arguments
(`could not convert string to float: 'abc'`,
`invalid literal for int() with base 10: 'x'`); `Failed to find group: <name>`;
`Unknown same property : <prop>Choose one of : dict_keys([...])` (no space
before "Choose"; the list is
`['fragment', 'x', 'y', 'z', 'residue', 'segment', 'name', 'type', 'resname', 'resid', 'segid', 'mass', 'charge', 'radius', 'bfactor', 'resnum', 'topologyattr', 'index', 'resindex', 'segindex', 'atomattr', 'id', None, 'element', 'record_type', 'chainID', 'tempfactor', 'formalcharge', 'occupancy', 'altLoc', 'gbscreen', 'solventradius', 'nbindex', 'rmin', 'epsilon', 'rmin14', 'epsilon14', 'aromaticity', 'chirality', 'residueattr', 'icode', 'moltype', 'molnum', 'segmentattr', 'bonds', 'ureybradleys', 'angles', 'dihedrals', 'impropers', 'cmaps', 'atomnum', 'atomiccharge', 'chargegroup', 'model', 'type_index']`);
`Invalid operator : …`; `too many values to unpack (expected 2)`.

Errors that are **not** wrapped (raised as-is):

* `Invalid expression given: …` (range joining, §2.1) — SelectionError.
* Nested errors keep the innermost message: `byres name` →
  `Selection failed: 'Unexpected token 'None''`.
* A missing argument to a fixed-argument keyword when END is consumed:
  `around` / `sphlayer 1` / `cyzone 1 2` / `point 1 2` / `point 6 0 1` →
  TypeError `float() argument must be a string or a real number, not 'NoneType'`;
  `atom A` → TypeError `int() argument must be a string, a bytes-like object or a real number, not 'NoneType'`;
  `atom` → IndexError `pop from an empty deque`; `same` → IndexError
  `deque index out of range`; and whenever a keyword successfully swallows
  END as its last argument (`atom A 1`, `prop 5 <`), the parser then fails
  with IndexError `deque index out of range`.
* `group updating` → TypeError
  `The 'updating' keyword is reserved and cannot be used as a selection group name.`

### 4.4 Fixed-argument keywords (raw tokens, no keyword stop, no unescaping)

`around R`, `sphzone R`: one float. `sphlayer R1 R2`, `isolayer R1 R2`: two
floats. `cyzone R ZMAX ZMIN`: three floats. `cylayer R1 R2 ZMAX ZMIN`: four.
`point X Y Z R`: four floats, no operand. `atom SEGID RESID NAME`: three raw
tokens, RESID must be an integer literal. `group NAME`: one raw token.
`same PROP as`: one raw token, which must be in the "same" list above, then
literally `as`. Floats use general float syntax (`around -1 …` parses;
see the grid-method spec for what negative cutoffs do).

### 4.5 Examples: precedence

| Selection string | MDAnalysis result |
|---|---|
| `resname G and name P or name N1` | 6 atoms: 1,6,11,16,24,26 |
| `name N1 or name P and resname G` | 4 atoms: 1,6,24,26 |
| `resname G and (name P or name N1)` | 4 atoms: 1,6,24,26 |
| `name P or name N1 and resid 1` | 2 atoms: 1,6 |
| `name P or (name N1 and resid 1)` | 10 atoms: 1,6,9,14,18,20,24,27,30,32 |
| `not name P and resname G` | 9 atoms: 2-8,25,26 |
| `not (name P and resname G)` | 58 atoms: 2-23,25-60 |
| `name P and not resname G` | 7 atoms: 9,14,18,20,27,30,32 |
| `name P and not not resname G` | 2 atoms: 1,24 |
| `not not not name P` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `not resid 1 or resid 2` | 52 atoms: 9-60 |
| `not (resid 1 or resid 2)` | 47 atoms: 14-60 |
| `resid 1 and not name P or name H3` | 8 atoms: 2-8,23 |
| `global resid 1 or resid 2` | 13 atoms: 1-13 |
| `name P and global resname G` | 2 atoms: 1,24 |
| `byres name N9 or name H3` | 12 atoms: 1-8,20-23 |
| `byres (name N9 or name H3)` | 12 atoms: 1-8,20-23 |
| `byres name N9 and resid 2` | 0 atoms: (none) |
| `byres name N9 and not name P` | 8 atoms: 1-8 |
| `not byres name N9 and resname G` | 52 atoms: 9-60 |
| `resid 1 and byres name N1 or name H3` | 8 atoms: 1-8 |
| `byres byres name N9` | 8 atoms: 1-8 |
| `same resid as name N9 or resid 2` | 13 atoms: 1-13 |
| `not same resid as name N9 or resid 2` | 47 atoms: 14-60 |
| `name P and same resid as name N9 and name C4'` | 0 atoms: (none) |
| `around 1.3 name N9 or name H1` | 4 atoms: 6,8,57,58 |
| `(around 1.3 name N9) or name H1` | 4 atoms: 6,7,56,59 |
| `resid 1 and around 1.3 name N9 or name H1` | 2 atoms: 6,8 |
| `not around 1.3 name N9` | 59 atoms: 1-5,7-60 |
| `sphzone 2 name N9 and name H1` | 0 atoms: (none) |
| `point 6 0 1 1.5 and name N9` | 1 atoms: 5 |
| `all and name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `(name P)` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `((name P))` | 9 atoms: 1,9,14,18,20,24,27,30,32 |

### 4.6 Examples: parse and evaluation errors

| Selection string | MDAnalysis result |
|---|---|
| `foo` | SelectionError: `Unknown selection token: 'foo'` |
| `foo bar` | SelectionError: `Unknown selection token: 'foo'` |
| `name P and foo` | SelectionError: `Unknown selection token: 'foo'` |
| `NAME P` | SelectionError: `Unknown selection token: 'NAME'` |
| `and name P` | SelectionError: `Unknown selection token: 'and'` |
| `or` | SelectionError: `Unknown selection token: 'or'` |
| `name P and` | SelectionError: `Unknown selection token: 'None'` |
| `name P or` | SelectionError: `Unknown selection token: 'None'` |
| `not` | SelectionError: `Unknown selection token: 'None'` |
| `byres` | SelectionError: `Unknown selection token: 'None'` |
| `same resid as` | SelectionError: `Unknown selection token: 'None'` |
| `name P and and name N1` | SelectionError: `Unknown selection token: 'and'` |
| `name P or or name N1` | SelectionError: `Unknown selection token: 'or'` |
| `not and name P` | SelectionError: `Unknown selection token: 'and'` |
| `name P not name N1` | SelectionError: `Unexpected token at end of selection string: 'not'` |
| `name P name P` | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `(` | SelectionError: `Unknown selection token: 'None'` |
| `)` | SelectionError: `Unknown selection token: ')'` |
| `()` | SelectionError: `Unknown selection token: ')'` |
| `( )` | SelectionError: `Unknown selection token: ')'` |
| `not ( )` | SelectionError: `Unknown selection token: ')'` |
| `name P and ( )` | SelectionError: `Unknown selection token: ')'` |
| `name P and )` | SelectionError: `Unknown selection token: ')'` |
| `(name P` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `((name P)` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `resid 1 and (name P` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `name P and (resid 1 or` | SelectionError: `Unknown selection token: 'None'` |
| `name P )` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `(name P))` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `name P (` | SelectionError: `Unexpected token at end of selection string: '('` |
| `same resid name N9` | SelectionError: `Unexpected token: 'name' Expected: 'as'` |
| `same resid of name N9` | SelectionError: `Unexpected token: 'of' Expected: 'as'` |
| `same resid` | SelectionError: `Unexpected token: 'None' Expected: 'as'` |
| `same foo as name P` | SelectionError: `Selection failed: 'Unknown same property : fooChoose one of : dict_keys([ …SAME-LIST… ])'` |
| `same` | IndexError: `deque index out of range` |
| `around name N9` | SelectionError: `Selection failed: 'could not convert string to float: 'name''` |
| `around abc name N9` | SelectionError: `Selection failed: 'could not convert string to float: 'abc''` |
| `around 2` | SelectionError: `Unknown selection token: 'None'` |
| `around` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` |
| `sphlayer 1` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` |
| `cyzone 1 2` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` |
| `point 1 2` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` |
| `point 6 0 1` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` |
| `atom` | IndexError: `pop from an empty deque` |
| `atom A` | TypeError: `int() argument must be a string, a bytes-like object or a real number, not 'NoneType'` |
| `atom A 1` | IndexError: `deque index out of range` |
| `atom A 1A P` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: '1A''` |
| `atom A x P` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: 'x''` |
| `group` | SelectionError: `Selection failed: 'Failed to find group: None'` |
| `group foo` | SelectionError: `Selection failed: 'Failed to find group: foo'` |
| `group updating` | TypeError: `The 'updating' keyword is reserved and cannot be used as a selection group name.` |
| `not prop x 5` | SelectionError: `Unknown selection token: 'prop'` |
| `byres name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `name P and resid 1 :` | SelectionError: `Invalid expression given: 1 : ` |
| `resid 1 : and foo` | SelectionError: `Invalid expression given: 1 : ` |
| `foo and resid 1 :` | SelectionError: `Unknown selection token: 'foo'` |
| `resid 1 : 3 and` | SelectionError: `Unknown selection token: 'None'` |
| `charge 0 and prop foo > 1` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `prop foo > 1 and charge 0` | NoDataError: `This Universe does not contain charge information` |
| `charge 0 or prop foo > 1` | NoDataError: `This Universe does not contain charge information` |
| `prop foo > 1 or charge 0` | SelectionError: `Expected one of [ …PROPERTY-LIST… ]` |
| `bonded name P` | NoDataError: `This Universe does not contain bonds information` |
| `same fragment as name P` | NoDataError: `AtomGroup.fragments not available; this requires bonds` |
