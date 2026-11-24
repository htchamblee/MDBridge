# Selection language: parser supplement (tokens, grammar, keyword catalogue, evaluation)

Clean-room behaviour specification derived by observing MDAnalysis 2.11.0
(`select_atoms`, NumPy 2.5.5, Python 3.13, no RDKit installed). It is a
**supplement** to `selection-values.md` (cited below as **[SV §n]**). The two
files together describe the whole selection-string language for a single,
static (non-updating) selection. Where [SV] already states a rule it is only
referenced here, not repeated.

## 0. Conventions and fixtures

* Result notation is the one used in [SV]: `N atoms: <list>` where the list
  gives **index + 1** of every selected atom (for the PDB fixture this equals
  the serial number), ranges compressed (`1-8`). When a result is not in
  ascending order, or contains duplicates, the list is printed verbatim in
  result order without compression.
* `(+float-equality warning)` = the warning text given in [SV §2.5];
  `(+bfactor-alias warning)` = the alias warning text in [SV §1.3] (note that
  one evaluation of `bfactor` emits it **twice**); `(+empty-string warning)` =
  `Empty string to select atoms, empty group returned.` Other warnings are
  quoted verbatim. `2×` means the warning was emitted twice.
* Errors are written `ExceptionClass: \`message\``. Inside selection strings and
  messages, a character that is whitespace (other than the plain space) or
  non-printable is written `<U+XXXX>`; inside a message `\n` stands for one
  newline character.
* In the first column, several selection strings passed to **one** call are
  written one after the other separated by ` , `; `(empty string)` is the
  zero-length string.
* "the fixture" = `tests/data/topology/fixture_chains.pdb` exactly as described
  in [SV, "Example fixture"]. "the fixture without icodes" = the same topology
  after removing the insertion-code attribute.
* **SYN** is a 12-atom synthetic system built in memory (no box; y = z = 0;
  no elements, chain IDs, occupancies, radii, record types, alternate
  locations; ids = index + 1):

| index+1 | name | type | resname | resid | icode | resnum | resindex | segid | segindex | mass | charge | formal charge | aromaticity | chirality | tempfactor | x |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | N | N | ALA | 1 | (empty) | 1 | 0 | PRO | 0 | 14.007 | −0.3 | 0 | false | (empty) | 0.1 | 0.0 |
| 2 | CA | C | ALA | 1 | (empty) | 1 | 0 | PRO | 0 | 12.011 | 0.1 | 0 | true | S | 0.2 | 1.5 |
| 3 | C | C | ALA | 1 | (empty) | 1 | 0 | PRO | 0 | 12.011 | 0.5 | 0 | true | (empty) | 0.3 | 3.0 |
| 4 | O | O | ALA | 1 | (empty) | 1 | 0 | PRO | 0 | 15.999 | −0.5 | 0 | false | (empty) | 0.4 | 4.5 |
| 5 | N | N | GLY | 2 | (empty) | 2 | 1 | PRO | 0 | 14.007 | −0.3 | 0 | false | (empty) | 0.5 | 6.0 |
| 6 | CA | C | GLY | 2 | (empty) | 2 | 1 | PRO | 0 | 12.011 | 0.1 | 0 | true | R | 0.6 | 7.5 |
| 7 | C | C | GLY | 2 | (empty) | 2 | 1 | PRO | 0 | 12.011 | 0.5 | 0 | true | (empty) | 0.7 | 9.0 |
| 8 | O | O | GLY | 2 | (empty) | 2 | 1 | PRO | 0 | 15.999 | −0.5 | 0 | false | (empty) | 0.8 | 10.5 |
| 9 | OW | OW | SOL | 3 | (empty) | 3 | 2 | SOL | 1 | 15.999 | −0.834 | 0 | false | (empty) | 0.9 | 12.0 |
| 10 | HW1 | HW | SOL | 3 | (empty) | 3 | 2 | SOL | 1 | 1.008 | 0.417 | 0 | true | (empty) | 1.0 | 13.5 |
| 11 | HW2 | HW | SOL | 3 | (empty) | 3 | 2 | SOL | 1 | 1.008 | 0.417 | 0 | true | (empty) | 1.1 | 15.0 |
| 12 | NA | NA | NA | 3 | B | 4 | 3 | SOL | 1 | 22.99 | 1.0 | 1 | false | (empty) | 1.2 | 16.5 |

  Bonds (by index+1): 1–2, 2–3, 3–4, 3–5, 5–6, 6–7, 7–8, 9–10, 9–11 (atom 12
  has no bond). All float attributes are held in double precision; positions
  in single precision. Further small systems are described where they are
  used.

---

## 1. The selection call as a whole

### 1.1 Inputs

A selection call receives: one or more selection strings (the first is
mandatory); a *periodic* switch (default on; only affects the distance
keywords, see `selection-grid-method.md`); two float-equality tolerances
*atol* (default 1e-8) and *rtol* (default 1e-5); a *sorted* switch (default
on); and optionally any number of **named groups** (name → atom group) that
the `group` keyword can refer to. The call is made on a *searched group*: the
whole system, or any atom group of it (possibly unordered, possibly with
repeated atoms).

### 1.2 Order of processing (observable through which error wins)

1. If the **first** string is empty (zero characters), the call returns an
   empty group and issues the empty-string warning. Nothing else is looked
   at: later strings are not parsed and named groups are not checked
   (`'' , 'foo'` returns empty; an invalid named group is not reported).
2. Every named group must be an atom group; otherwise TypeError
   `Passed groups must be AtomGroups. You provided <type name> for group '<name>'`
   (`<type name>` is the Python class name of the value, e.g. `str`, `list`).
   This is checked before any parsing.
3. All strings are parsed, in order. The first parse error aborts the call.
   An empty string in second or later position is **not** special: it is
   parsed and fails with `Unknown selection token: 'None'`.
4. Each parsed string is evaluated against the searched group, in order; the
   first evaluation error aborts the call.
5. The results are concatenated in string order. Duplicates across strings
   are kept; each individual string's result has no duplicates.

With *sorted* on (the default), each string's result is sorted by atom index.
With *sorted* off, see §6.11.

The keyword `global` is meaningful only when the searched group is not the
whole system (§6.1). On the whole system it changes nothing.

| Selection string | MDAnalysis result |
|---|---|
| `name P` , `resid 1` | 17 atoms: 1,9,14,18,20,24,27,30,32,1,2,3,4,5,6,7,8 |
| `resid 1` , `name P` | 17 atoms: 1,2,3,4,5,6,7,8,1,9,14,18,20,24,27,30,32 |
| `name P` , `name P` | 18 atoms: 1,9,14,18,20,24,27,30,32,1,9,14,18,20,24,27,30,32 |
| (empty string) , `name P` | 0 atoms: (none) (+empty-string warning) |
| (empty string) , `foo` | 0 atoms: (none) (+empty-string warning) |
| `name P` , (empty string) | SelectionError: `Unknown selection token: 'None'` |
| `name P` , ` ` | SelectionError: `Unknown selection token: 'None'` |
| `name P` , `foo` | SelectionError: `Unknown selection token: 'foo'` |
| `charge 0` , `foo` | SelectionError: `Unknown selection token: 'foo'` |
| `charge 0` , `name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `charge 0` , `resid 1 :` | SelectionError: `Invalid expression given: 1 : ` |
| `name P` , `charge 0` | NoDataError: `This Universe does not contain charge information` |
| `charge 0` , `name P` | NoDataError: `This Universe does not contain charge information` |
| `name P and` , `resid 1` | SelectionError: `Unknown selection token: 'None'` |
| `index 5` , `index 1` , `index 5` | 3 atoms: 6,2,6 |
| `none` | SelectionError: `Unknown selection token: 'none'` |
| `resid 1 or` , `resid 2` | SelectionError: `Unknown selection token: 'None'` |

### 1.3 Strings that contain no tokens

| Selection string | MDAnalysis result |
|---|---|
| (empty string) | 0 atoms: (none) (+empty-string warning) |
| ` ` | SelectionError: `Unknown selection token: 'None'` |
| `<U+0009>` | SelectionError: `Unknown selection token: 'None'` |
| `<U+000A>` | SelectionError: `Unknown selection token: 'None'` |
| `()` | SelectionError: `Unknown selection token: ')'` |
| `( )` | SelectionError: `Unknown selection token: ')'` |

| call | MDAnalysis result |
|---|---|
| first string empty, grp = the string 'x' | 0 atoms: (none) (+empty-string warning) |
| `name P`, grp = the string 'x' | TypeError: `Passed groups must be AtomGroups. You provided str for group 'grp'` |
| `foo`, grp = the string 'x' | TypeError: `Passed groups must be AtomGroups. You provided str for group 'grp'` |

---

## 2. Tokenisation (supplement to [SV §0])

* The only characters that separate tokens are whitespace characters, and
  the only characters that form tokens on their own are `(` and `)`.
  Whitespace is exactly the set of characters that Python treats as
  whitespace: U+0009–U+000D, U+001C–U+001F, U+0020, U+0085, U+00A0, U+1680,
  U+2000–U+200A, U+2028, U+2029, U+202F, U+205F, U+3000. Zero-width
  characters such as U+200B and U+FEFF are **not** whitespace and stay inside
  tokens. Commas, semicolons, quotes, `!`, `|`, `&`, `;`, `=` are ordinary
  token characters (no quoting mechanism exists; `"P"` is a three-character
  value).
* Every token is compared verbatim (case-sensitive) against the keyword
  list. The four-character token `None` is an ordinary token: as a value it
  is the text `None`; as a would-be keyword it gives
  `Unknown selection token: 'None'`, which is textually identical to the
  message produced at end of input.
* Parentheses are recognised anywhere, also glued to words: `(name P)or(name N1)`
  is `(` `name` `P` `)` `or` `(` `name` `N1` `)`. Words glued to `and`/`or`
  are not split: `andresid` is one token.
* Decimal digits in numbers (resid values, integer/float keywords, numeric
  arguments of `around`, `point`, `atom` …) may be any Unicode decimal digit
  (general category Nd: Arabic-Indic `١`, full-width `１`, …), mixed freely
  with ASCII digits. Superscripts (`²`) and dash-like characters (`–`, `—`)
  are not digits / separators. An insertion code may be any single Unicode
  letter, digit or underscore.

| Selection string | MDAnalysis result |
|---|---|
| `name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `   name    P   ` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+0009>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+000A>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `<U+0009>name P<U+000A>` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+000D>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+000B>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+000C>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+001C>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+001F>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+0085>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+00A0>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+1680>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+2003>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+2028>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+202F>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+3000>P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+200B>P` | SelectionError: `Unknown selection token: 'name<U+200B>P'` |
| `name P<U+200B>` | 0 atoms: (none) |
| `name<U+FEFF>P` | SelectionError: `Unknown selection token: 'name<U+FEFF>P'` |
| `(name P)` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `( name P )` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `(name P)or(name N1)` | 13 atoms: 1,6,9,11,14,16,18,20,24,26,27,30,32 |
| `not(name P)` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `name P and(resid 1)` | 1 atoms: 1 |
| `name P and(resid 1)or(name H3)` | 2 atoms: 1,23 |
| `byres(name N9)` | 8 atoms: 1-8 |
| `around 1.3(name N9)` | 1 atoms: 6 |
| `same resid as(name N9)` | 8 atoms: 1-8 |
| `name P)or(name N1` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `name P andresid 1` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `resid 1and name P` | SelectionError: `Selection failed: 'Failed to parse value: 1and'` |
| `name P,N1` | 0 atoms: (none) |
| `name P;` | 0 atoms: (none) |
| `name "P"` | 0 atoms: (none) |
| `name 'P'` | 0 atoms: (none) |
| `name P\ and resid 1` | 0 atoms: (none) |
| `resid 1 and(not(name P))` | 7 atoms: 2-8 |
| `name P and((resid 1))` | 1 atoms: 1 |
| `(((name P)))` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `name<U+0009><U+0009>(P)` | SelectionError: `Selection failed: 'Unexpected token '(''` |
| `resid 1:2)` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `(resid 1:2` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `resid (1:2)` | SelectionError: `Selection failed: 'Unexpected token: '(''` |
| `resid 1:(2)` | SelectionError: `Selection failed: 'Failed to parse value: 1:'` |

| Selection string | MDAnalysis result |
|---|---|
| `resid ١` | 8 atoms: 1-8 |
| `resid 1١` | 3 atoms: 27-29 |
| `index ١` | 1 atoms: 2 |
| `index １` | 1 atoms: 2 |
| `resid １:３` | 17 atoms: 1-17 |
| `mass ١٢.011` | 19 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48 (+float-equality warning) |
| `around ١.3 name N9` | 1 atoms: 6 |
| `atom A ١ P` | 1 atoms: 1 |
| `resid ²` | SelectionError: `Selection failed: 'Failed to parse value: ²'` |
| `index ²` | SelectionError: `Selection failed: 'Failed to parse number: ²'` |
| `resid 1é` | 0 atoms: (none) |
| `resid 1 – 3` | SelectionError: `Selection failed: 'Failed to parse value: –'` |
| `resid 1—3` | SelectionError: `Selection failed: 'Failed to parse value: 1—3'` |
| `name None` | 0 atoms: (none) |
| `None` | SelectionError: `Unknown selection token: 'None'` |
| `name<U+200B>P` | SelectionError: `Unknown selection token: 'name<U+200B>P'` |
| `<U+00A0><U+3000>` | SelectionError: `Unknown selection token: 'None'` |

---

## 3. Grammar, binding and associativity (supplement to [SV §4.1])

### 3.1 Productions in words

The parser reads tokens left to right with one token of look-ahead and
never backtracks.

* An **expression** is: one *primary*, followed by zero or more pairs
  "*binary operator*, *primary*", where the binary operators are exactly the
  tokens `and` and `or`. The pairs are combined strictly left to right
  (both operators have the same strength): `A or B and C or D` means
  `((A or B) and C) or D`. An expression ends at the first token, after a
  complete primary, that is not `and`/`or` (that token is left for whoever
  asked for the expression: a closing `)`, the end check, or nothing).
* A **primary** is one of:
  1. `(` *expression* `)`. If the token after the inner expression is not
     `)`: SelectionError `Unexpected token: '<token>' Expected: ')'`.
  2. `not` *primary* or `global` *primary* (the operand is exactly one
     primary, so these never extend over a following `and`/`or`).
  3. A **swallowing keyword** — `byres`, `bonded`, `same P as`, `around R`,
     `sphzone R`, `sphlayer R1 R2`, `isolayer R1 R2`, `cyzone R ZMAX ZMIN`,
     `cylayer R1 R2 ZMAX ZMIN` — whose operand is a whole *expression*
     (primary plus any number of `and`/`or` pairs). It therefore extends to
     the matching `)` or to the end of input, through every `and`/`or`,
     even when the keyword itself is the operand of `not`, `global` or the
     right-hand side of `and`/`or`.
  4. Any other keyword with its arguments (catalogue in §4); these take no
     operand.
  Any other token at a primary position (a non-keyword word, `and`, `or`,
  `)`, END) → SelectionError `Unknown selection token: '<token>'`
  (END prints `None`).
* The **whole string** is one expression followed by END; otherwise
  SelectionError `Unexpected token at end of selection string: '<token>'`.

Consequences worth stating explicitly:

* The right operand of `and`/`or` is a single primary. Since a swallowing
  keyword is a primary that itself swallows the rest, `X and byres A or B`
  is `X and (byres (A or B))`, whereas `byres A or B and X` is
  `byres ((A or B) and X)`.
* `not not A` = `not (not A)`; `not global A` = `not (global A)`;
  `global not A` = `global (not A)`.
* `not byres A or B` = `not (byres (A or B))` but `not (byres A) or B` =
  `(not (byres A)) or B`.
* Swallowing keywords nest: `byres around 1.3 A or B` =
  `byres (around 1.3 (A or B))`.
* `point X Y Z R`, `prop …`, `atom …`, `group …` and the value keywords are
  ordinary primaries; an `and`/`or` after them belongs to the enclosing
  expression.
* Evaluation order (which error is reported first) is in [SV §4.2]: for
  `and` the right operand is evaluated before the left, for `or` the left
  first; a unary keyword evaluates its operand first, then its own
  attribute access.

### 3.2 Grouping verified on the fixture

Each row was run twice (as written and fully parenthesised) and gave the
identical result.

| Selection string | equivalent fully parenthesised form | result (identical for both strings) |
|---|---|---|
| `name P or name N1 and resid 1` | `(name P or name N1) and resid 1` | 2 atoms: 1,6 |
| `resid 1 and name P or name N1` | `(resid 1 and name P) or name N1` | 5 atoms: 1,6,11,16,26 |
| `name P and resid 1 or name N1 and resid 10` | `((name P and resid 1) or name N1) and resid 10` | 1 atoms: 26 |
| `name P and resid 1 or (name N1 and resid 10)` | `(name P and resid 1) or (name N1 and resid 10)` | 2 atoms: 1,26 |
| `name P or name N1 or name H3` | `(name P or name N1) or name H3` | 14 atoms: 1,6,9,11,14,16,18,20,23,24,26,27,30,32 |
| `name P or resid 1 and resid 2 or name H3` | `((name P or resid 1) and resid 2) or name H3` | 2 atoms: 9,23 |
| `name P or (resid 1 and resid 2) or name H3` | `(name P or (resid 1 and resid 2)) or name H3` | 10 atoms: 1,9,14,18,20,23,24,27,30,32 |
| `not name P or name N1` | `(not name P) or name N1` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `not name P and not name N1` | `(not name P) and (not name N1)` | 47 atoms: 2-5,7,8,10,12,13,15,17,19,21-23,25,28,29,31,33-60 |
| `not (name P or name N1)` | `not (name P or name N1)` | 47 atoms: 2-5,7,8,10,12,13,15,17,19,21-23,25,28,29,31,33-60 |
| `not resid 1 and name P or name H3` | `((not resid 1) and name P) or name H3` | 9 atoms: 9,14,18,20,23,24,27,30,32 |
| `resid 1 and not name P and not name N1` | `(resid 1 and (not name P)) and (not name N1)` | 6 atoms: 2-5,7,8 |
| `not not name P` | `not (not name P)` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `not not not resid 1 or resid 2` | `(not (not (not resid 1))) or resid 2` | 52 atoms: 9-60 |
| `global name P and resid 1` | `(global name P) and resid 1` | 1 atoms: 1 |
| `not global name P` | `not (global name P)` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `global not name P or name N1` | `(global (not name P)) or name N1` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `byres name N9 or name H3 and resid 4` | `byres ((name N9 or name H3) and resid 4)` | 4 atoms: 20-23 |
| `resid 1 or byres name H3 and name P` | `resid 1 or (byres (name H3 and name P))` | 8 atoms: 1-8 |
| `(byres name H3) and name P` | `(byres name H3) and name P` | 1 atoms: 20 |
| `not byres name H3 or name N9` | `not (byres (name H3 or name N9))` | 48 atoms: 9-19,24-60 |
| `not (byres name H3) or name N9` | `(not (byres name H3)) or name N9` | 56 atoms: 1-19,24-60 |
| `name P and byres name N9 or name H3` | `name P and (byres (name N9 or name H3))` | 2 atoms: 1,20 |
| `name P and byres name N9 or name H3 and not resid 1` | `name P and (byres ((name N9 or name H3) and (not resid 1)))` | 1 atoms: 20 |
| `byres name N9 and byres name P or name H3` | `byres (name N9 and (byres (name P or name H3)))` | 8 atoms: 1-8 |
| `same resid as name N9 or name H3` | `same resid as (name N9 or name H3)` | 12 atoms: 1-8,20-23 |
| `same resid as name N9 and name P` | `same resid as (name N9 and name P)` | 0 atoms: (none) |
| `name P or same resid as name N9 and name P` | `name P or (same resid as (name N9 and name P))` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `around 1.3 name N9 or name H1` | `around 1.3 (name N9 or name H1)` | 4 atoms: 6,8,57,58 |
| `not around 1.3 name N9 or name H1` | `not (around 1.3 (name N9 or name H1))` | 56 atoms: 1-5,7,9-56,59,60 |
| `around 1.3 name N9 and not name O6` | `around 1.3 (name N9 and (not name O6))` | 1 atoms: 6 |
| `sphzone 2 name N9 or name H1` | `sphzone 2 (name N9 or name H1)` | 3 atoms: 31-33 |
| `point 6 0 1 1.5 or name P` | `(point 6 0 1 1.5) or name P` | 11 atoms: 1,5,6,9,14,18,20,24,27,30,32 |
| `point 6 0 1 1.5 and name N9 or name P` | `((point 6 0 1 1.5) and name N9) or name P` | 10 atoms: 1,5,9,14,18,20,24,27,30,32 |
| `name P and point 6 0 1 1.5 or name N9` | `(name P and (point 6 0 1 1.5)) or name N9` | 1 atoms: 5 |
| `prop x < 5 or name N9` | `(prop x < 5) or name N9` | 5 atoms: 1-5 |
| `name P or prop x < 5 and resid 1` | `(name P or (prop x < 5)) and resid 1` | 4 atoms: 1-4 |
| `all and not all or name P` | `(all and (not all)) or name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `resid 1 and (name P or name N1) and not (name P)` | `(resid 1 and (name P or name N1)) and (not (name P))` | 1 atoms: 6 |
| `around 1.3 around 1.3 name N9` | `around 1.3 (around 1.3 name N9)` | 1 atoms: 5 |
| `byres around 1.3 name N9 or name H1` | `byres (around 1.3 (name N9 or name H1))` | 14 atoms: 1-8,55-60 |

Pairs that differ only by parentheses:

| Selection string | MDAnalysis result |
|---|---|
| `name P or name N1 and resid 1` | 2 atoms: 1,6 |
| `name P or (name N1 and resid 1)` | 10 atoms: 1,6,9,14,18,20,24,27,30,32 |
| `not name P or name N1` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `not (name P or name N1)` | 47 atoms: 2-5,7,8,10,12,13,15,17,19,21-23,25,28,29,31,33-60 |
| `byres name H3 or name N9` | 12 atoms: 1-8,20-23 |
| `(byres name H3) or name N9` | 5 atoms: 5,20-23 |
| `around 1.3 name N9 or name H1` | 4 atoms: 6,8,57,58 |
| `(around 1.3 name N9) or name H1` | 4 atoms: 6,7,56,59 |
| `same resid as name N9 or name H3` | 12 atoms: 1-8,20-23 |
| `(same resid as name N9) or name H3` | 9 atoms: 1-8,23 |
| `name P and byres name N9 or name H3` | 2 atoms: 1,20 |
| `(name P and byres name N9) or name H3` | 2 atoms: 1,23 |

### 3.3 Nesting depth

There is no explicit limit, but MDAnalysis hits Python's recursion limit
and raises RecursionError `maximum recursion depth exceeded` for very deep
nesting. Observed thresholds (largest count that still works, on the
fixture): 329 consecutive `not` before `name P` (also 329 consecutive
`byres`); 494 nested parenthesis pairs; a chain of 495 `name P` terms joined
by `or` (or by `and`). A port need not reproduce these thresholds exactly,
but should not crash on inputs within them.

---

## 4. Keyword catalogue

### 4.1 What each keyword consumes

"value list" = the value-collection rule of [SV §1.1] (stop at a keyword /
`and` / `or` / parenthesis / END, strip one leading backslash).
"raw token" = the next token taken unconditionally, whatever it is (keyword,
`and`, parenthesis, or END), with no backslash removal. "operand" = a whole
expression (§3.1, swallowing). For every raw numeric token the conversion is
general float syntax ([SV §2.5]: sign, `_` between digits, exponent `e`/`E`,
`inf`, `nan`, Unicode digits) or, for `atom`, general integer syntax
([SV §2.3]).

| keyword(s) | consumes | operand | evaluation (details in §6) |
|---|---|---|---|
| `and`, `or` | — | binary | intersection / union |
| `not` | — | one primary | searched group minus operand |
| `global` | — | one primary | operand evaluated on the whole system; result not restricted to the searched group |
| `all` | nothing | none | the searched group |
| `protein`, `nucleic`, `water`, `backbone`, `nucleicbackbone`, `nucleicbase`, `nucleicsugar` | nothing | none | name/resname lists (§6.6) |
| `aromatic` | nothing | none | aromaticity true |
| `byres` | nothing | operand | whole residues (§6.7) |
| `bonded` | nothing | operand | bond partners (§6.7) |
| `same` | raw token P (must be a known property name), then the literal token `as` | operand | equal property (§6.7) |
| `around`, `sphzone` | 1 raw float | operand | §6.10 |
| `sphlayer`, `isolayer` | 2 raw floats (inner, outer) | operand | §6.10 |
| `cyzone` | 3 raw floats (R, ZMAX, ZMIN) | operand | §6.10 |
| `cylayer` | 4 raw floats (R1, R2, ZMAX, ZMIN) | operand | §6.10 |
| `point` | 4 raw floats (X, Y, Z, R) | none | §6.10 |
| `atom` | 3 raw tokens: segid text, resid (general integer syntax), name text | none | §6.8 |
| `group` | 1 raw token: group name | none | §6.8 |
| `prop` | raw tokens as in [SV §3] | none | [SV §3] |
| `smarts` | special raw collection (§4.4) | none | needs RDKit (§4.4) |
| string keywords `name`, `type`, `element`, `record_type`, `chainID`, `chainid`, `altLoc`, `altloc`, `resname`, `icode`, `moltype`, `segid` | value list, ≥ 1 value | none | glob match [SV §1.2], §6.5 |
| integer keywords `index`, `resindex`, `segindex`, `id`, `bynum`, `formalcharge`, `nbindex`, `resnum`, `molnum`, `model` | value list, ≥ 1 value, joined into ranges | none | [SV §2.3] |
| float keywords `radius`, `tempfactor`, `bfactor`, `mass`, `charge`, `occupancy`, `gbscreen`, `solventradius`, `rmin`, `epsilon`, `rmin14`, `epsilon14` | value list, ≥ 1 value, joined into ranges | none | [SV §2.5], §6.2 |
| `resid` | value list, ≥ 1 value, joined into ranges | none | [SV §2.4], §6.3 |
| `aromaticity` | value list, 0 or more (none = `true`) | none | [SV §1.2] |
| `chirality` | value list, ≥ 1 value | none | exact single-character text equality |

No other keywords exist (in particular there is no `none`, `hydrogen`,
`heavy`, `within`, `to`, `segname`, `resi`, `chain`, `fullgroup`).
`and`/`or`/`not` and all keywords are lower-case except `chainID`/`altLoc`,
which also exist in lower case.

### 4.2 End of input, `)` and `and` right after a keyword

When a raw token is required but the input is exhausted, the "token" is
END. What happens then depends on the conversion applied to it:

* float conversion of END → TypeError
  `float() argument must be a string or a real number, not 'NoneType'`
  (not wrapped);
* integer conversion of END (`atom` resid) → TypeError
  `int() argument must be a string, a bytes-like object or a real number, not 'NoneType'`;
* `atom` with nothing at all after it → IndexError `pop from an empty deque`;
* `same` with nothing after it → IndexError `deque index out of range`;
* `group` with nothing after it → the name is END →
  `Selection failed: 'Failed to find group: None'`;
* a keyword that **successfully** consumed END as its last raw token (e.g.
  `atom A 1` where END becomes the atom name) leaves the parser with no
  END marker → IndexError `deque index out of range`.

A `)` or `and`/`or` used where a raw float is needed is converted like any
text: `Selection failed: 'could not convert string to float: ')''`. A `)`
or `and` where a primary is needed (operand of `not`/`global`/swallowing
keyword, right side of `and`/`or`) → `Unknown selection token: ')'` /
`'and'`. A `)` or `and` where a value list starts → the "no values" message of
[SV §1.1].

The following matrix shows six probes for every keyword (searched group =
whole fixture). `K` stands for the keyword in the first column.

| keyword | `K` | `(K)` | `K and name P` | `K 1` | `(K 1)` | `K 1 and name P` |
|---|---|---|---|---|---|---|
| `all` | 60 atoms: 1-60 | 60 atoms: 1-60 | 9 atoms: 1,9,14,18,20,24,27,30,32 | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `not` | SelectionError: `Unknown selection token: 'None'` | SelectionError: `Unknown selection token: ')'` | SelectionError: `Unknown selection token: 'and'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` |
| `global` | SelectionError: `Unknown selection token: 'None'` | SelectionError: `Unknown selection token: ')'` | SelectionError: `Unknown selection token: 'and'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` |
| `byres` | SelectionError: `Unknown selection token: 'None'` | SelectionError: `Unknown selection token: ')'` | SelectionError: `Unknown selection token: 'and'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` |
| `around` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` | SelectionError: `Unknown selection token: 'None'` | SelectionError: `Unknown selection token: ')'` | SelectionError: `Unknown selection token: 'and'` |
| `sphlayer` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `isolayer` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `sphzone` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` | SelectionError: `Unknown selection token: 'None'` | SelectionError: `Unknown selection token: ')'` | SelectionError: `Unknown selection token: 'and'` |
| `cyzone` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `cylayer` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `point` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'could not convert string to float: ')''` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `atom` | IndexError: `pop from an empty deque` | TypeError: `int() argument must be a string, a bytes-like object or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: 'name''` | TypeError: `int() argument must be a string, a bytes-like object or a real number, not 'NoneType'` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: ')''` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: 'and''` |
| `bonded` | SelectionError: `Unknown selection token: 'None'` | SelectionError: `Unknown selection token: ')'` | SelectionError: `Unknown selection token: 'and'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` | SelectionError: `Unknown selection token: '1'` |
| `group` | SelectionError: `Selection failed: 'Failed to find group: None'` | SelectionError: `Selection failed: 'Failed to find group: )'` | SelectionError: `Selection failed: 'Failed to find group: and'` | SelectionError: `Selection failed: 'Failed to find group: 1'` | SelectionError: `Selection failed: 'Failed to find group: 1'` | SelectionError: `Selection failed: 'Failed to find group: 1'` |
| `aromatic` | NoDataError: `This Universe does not contain aromaticity information` | NoDataError: `This Universe does not contain aromaticity information` | NoDataError: `This Universe does not contain aromaticity information` | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `smarts` | ImportError: `RDKit is required …` | ImportError: `RDKit is required …` | ImportError: `RDKit is required …` | ImportError: `RDKit is required …` | ImportError: `RDKit is required …` | ImportError: `RDKit is required …` |
| `resid` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 8 atoms: 1-8 | 8 atoms: 1-8 | 1 atoms: 1 |
| `bynum` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 1 atoms: 1 | 1 atoms: 1 | 1 atoms: 1 |
| `protein` | 15 atoms: 34-48 | 15 atoms: 34-48 | 0 atoms: (none) | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `nucleic` | 33 atoms: 1-33 | 33 atoms: 1-33 | 9 atoms: 1,9,14,18,20,24,27,30,32 | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `water` | 6 atoms: 55-60 | 6 atoms: 55-60 | 0 atoms: (none) | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `backbone` | 6 atoms: 34-37,41,42 | 6 atoms: 34-37,41,42 | 0 atoms: (none) | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `nucleicbackbone` | 9 atoms: 1,9,14,18,20,24,27,30,32 | 9 atoms: 1,9,14,18,20,24,27,30,32 | 9 atoms: 1,9,14,18,20,24,27,30,32 | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `nucleicbase` | 11 atoms: 5,6,8,11-13,16,17,22,26,29 | 11 atoms: 5,6,8,11-13,16,17,22,26,29 | 0 atoms: (none) | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `nucleicsugar` | 8 atoms: 3,10,15,19,25,28,31,33 | 8 atoms: 3,10,15,19,25,28,31,33 | 0 atoms: (none) | SelectionError: `Unexpected token at end of selection string: '1'` | SelectionError: `Unexpected token: '1' Expected: ')'` | SelectionError: `Unexpected token at end of selection string: '1'` |
| `prop` | AttributeError: `'NoneType' object has no attribute 'split'` | TypeError: `argument of type 'NoneType' is not iterable` | SelectionError: `Selection failed: 'Invalid operator : 'name' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])''` | TypeError: `argument of type 'NoneType' is not iterable` | SelectionError: `Unknown selection token: 'prop'` | SelectionError: `Unknown selection token: 'prop'` |
| `same` | IndexError: `deque index out of range` | SelectionError: `Selection failed: 'Unknown same property : )Choose one of : dict_keys([ …SAME-LIST… ])'` | SelectionError: `Selection failed: 'Unknown same property : andChoose one of : dict_keys([ …SAME-LIST… ])'` | SelectionError: `Selection failed: 'Unknown same property : 1Choose one of : dict_keys([ …SAME-LIST… ])'` | SelectionError: `Selection failed: 'Unknown same property : 1Choose one of : dict_keys([ …SAME-LIST… ])'` | SelectionError: `Selection failed: 'Unknown same property : 1Choose one of : dict_keys([ …SAME-LIST… ])'` |
| `index` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 1 atoms: 2 | 1 atoms: 2 | 0 atoms: (none) |
| `resindex` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 5 atoms: 9-13 | 5 atoms: 9-13 | 1 atoms: 9 |
| `segindex` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 10 atoms: 24-33 | 10 atoms: 24-33 | 4 atoms: 24,27,30,32 |
| `id` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 1 atoms: 1 | 1 atoms: 1 | 1 atoms: 1 |
| `name` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `type` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `element` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `radius` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain radius information` | NoDataError: `This Universe does not contain radius information` | NoDataError: `This Universe does not contain radius information` |
| `record_type` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `chainID` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `chainid` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `tempfactor` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 0 atoms: (none) (+float-equality warning) | 0 atoms: (none) (+float-equality warning) | 0 atoms: (none) (+float-equality warning) |
| `bfactor` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 0 atoms: (none) (+2× bfactor-alias warning, float-equality warning) | 0 atoms: (none) (+2× bfactor-alias warning, float-equality warning) | 0 atoms: (none) (+2× bfactor-alias warning, float-equality warning) |
| `mass` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 0 atoms: (none) (+float-equality warning) | 0 atoms: (none) (+float-equality warning) | 0 atoms: (none) (+float-equality warning) |
| `charge` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain charge information` | NoDataError: `This Universe does not contain charge information` | NoDataError: `This Universe does not contain charge information` |
| `formalcharge` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain formalcharge information` | NoDataError: `This Universe does not contain formalcharge information` | NoDataError: `This Universe does not contain formalcharge information` |
| `occupancy` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 60 atoms: 1-60 (+float-equality warning) | 60 atoms: 1-60 (+float-equality warning) | 9 atoms: 1,9,14,18,20,24,27,30,32 (+float-equality warning) |
| `altLoc` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `altloc` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `gbscreen` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain gbscreen information` | NoDataError: `This Universe does not contain gbscreen information` | NoDataError: `This Universe does not contain gbscreen information` |
| `solventradius` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain solventradius information` | NoDataError: `This Universe does not contain solventradius information` | NoDataError: `This Universe does not contain solventradius information` |
| `nbindex` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain nbindex information` | NoDataError: `This Universe does not contain nbindex information` | NoDataError: `This Universe does not contain nbindex information` |
| `rmin` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain rmin information` | NoDataError: `This Universe does not contain rmin information` | NoDataError: `This Universe does not contain rmin information` |
| `epsilon` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain epsilon information` | NoDataError: `This Universe does not contain epsilon information` | NoDataError: `This Universe does not contain epsilon information` |
| `rmin14` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain rmin14 information` | NoDataError: `This Universe does not contain rmin14 information` | NoDataError: `This Universe does not contain rmin14 information` |
| `epsilon14` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain epsilon14 information` | NoDataError: `This Universe does not contain epsilon14 information` | NoDataError: `This Universe does not contain epsilon14 information` |
| `aromaticity` | NoDataError: `This Universe does not contain aromaticity information` | NoDataError: `This Universe does not contain aromaticity information` | NoDataError: `This Universe does not contain aromaticity information` | SelectionError: `Selection failed: ''1' is an invalid value for boolean selection. Use 'True' or 'False''` | SelectionError: `Selection failed: ''1' is an invalid value for boolean selection. Use 'True' or 'False''` | SelectionError: `Selection failed: ''1' is an invalid value for boolean selection. Use 'True' or 'False''` |
| `chirality` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | NoDataError: `This Universe does not contain chirality information` | NoDataError: `This Universe does not contain chirality information` | NoDataError: `This Universe does not contain chirality information` |
| `resname` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `resnum` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | 8 atoms: 1-8 | 8 atoms: 1-8 | 1 atoms: 1 |
| `icode` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `moltype` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | AttributeError: `'Topology' object has no attribute 'moltypes'` | AttributeError: `'Topology' object has no attribute 'moltypes'` | AttributeError: `'Topology' object has no attribute 'moltypes'` |
| `molnum` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain molnum information` | NoDataError: `This Universe does not contain molnum information` | NoDataError: `This Universe does not contain molnum information` |
| `segid` | SelectionError: `Selection failed: 'Unexpected token 'None''` | SelectionError: `Selection failed: 'Unexpected token ')''` | SelectionError: `Selection failed: 'Unexpected token 'and''` | 0 atoms: (none) | 0 atoms: (none) | 0 atoms: (none) |
| `model` | SelectionError: `Selection failed: 'Unexpected token: 'None''` | SelectionError: `Selection failed: 'Unexpected token: ')''` | SelectionError: `Selection failed: 'Unexpected token: 'and''` | NoDataError: `This Universe does not contain model information` | NoDataError: `This Universe does not contain model information` | NoDataError: `This Universe does not contain model information` |

### 4.3 More malformed input

| Selection string | MDAnalysis result |
|---|---|
| `name P and` | SelectionError: `Unknown selection token: 'None'` |
| `name P or` | SelectionError: `Unknown selection token: 'None'` |
| `name P and or name N1` | SelectionError: `Unknown selection token: 'or'` |
| `and` | SelectionError: `Unknown selection token: 'and'` |
| `or` | SelectionError: `Unknown selection token: 'or'` |
| `not` | SelectionError: `Unknown selection token: 'None'` |
| `not not` | SelectionError: `Unknown selection token: 'None'` |
| `not or name P` | SelectionError: `Unknown selection token: 'or'` |
| `name P and not` | SelectionError: `Unknown selection token: 'None'` |
| `global` | SelectionError: `Unknown selection token: 'None'` |
| `global and name P` | SelectionError: `Unknown selection token: 'and'` |
| `(` | SelectionError: `Unknown selection token: 'None'` |
| `((` | SelectionError: `Unknown selection token: 'None'` |
| `)` | SelectionError: `Unknown selection token: ')'` |
| `))` | SelectionError: `Unknown selection token: ')'` |
| `( (name P) ` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `(name P))` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `(name P) (name N1)` | SelectionError: `Unexpected token at end of selection string: '('` |
| `(name P) name N1` | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `(name P) or` | SelectionError: `Unknown selection token: 'None'` |
| `(name P) and ()` | SelectionError: `Unknown selection token: ')'` |
| `name P and (` | SelectionError: `Unknown selection token: 'None'` |
| `name P and ( resid 1` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `name P and ( resid 1 ))` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `not ( name P` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `byres ( name P` | SelectionError: `Unexpected token: 'None' Expected: ')'` |
| `byres name P )` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `not name P )` | SelectionError: `Unexpected token at end of selection string: ')'` |
| `name P or )` | SelectionError: `Unknown selection token: ')'` |
| `( or name P )` | SelectionError: `Unknown selection token: 'or'` |
| `( name P or )` | SelectionError: `Unknown selection token: ')'` |
| `( and )` | SelectionError: `Unknown selection token: 'and'` |
| `( not )` | SelectionError: `Unknown selection token: ')'` |
| `( byres )` | SelectionError: `Unknown selection token: ')'` |
| `( global )` | SelectionError: `Unknown selection token: ')'` |
| `( around 1 )` | SelectionError: `Unknown selection token: ')'` |
| `( same resid as )` | SelectionError: `Unknown selection token: ')'` |
| `( same resid )` | SelectionError: `Unexpected token: ')' Expected: 'as'` |
| `name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `not name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `not ( name )` | SelectionError: `Selection failed: 'Unexpected token ')''` |
| `byres name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `byres resid` | SelectionError: `Selection failed: 'Unexpected token: 'None''` |
| `byres resid 1 :` | SelectionError: `Invalid expression given: 1 : ` |
| `global resid :` | SelectionError: `Invalid expression given:  : ` |
| `same resid as name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `around 1 name` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `around 1 resid` | SelectionError: `Selection failed: 'Unexpected token: 'None''` |
| `name P resid 1` | SelectionError: `Unexpected token at end of selection string: 'resid'` |
| `resid 1 name P` | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `all all` | SelectionError: `Unexpected token at end of selection string: 'all'` |
| `all name P` | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `protein and` | SelectionError: `Unknown selection token: 'None'` |
| `protein or protein` | 15 atoms: 34-48 |
| `protein protein` | SelectionError: `Unexpected token at end of selection string: 'protein'` |
| `not protein P` | SelectionError: `Unexpected token at end of selection string: 'P'` |
| `point 1 2 3 4 5` | SelectionError: `Unexpected token at end of selection string: '5'` |
| `point 1 2 3 4 name P` | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `point 1 2 3 x` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `point a 2 3 4` | SelectionError: `Selection failed: 'could not convert string to float: 'a''` |
| `point 1 2 3` | TypeError: `float() argument must be a string or a real number, not 'NoneType'` |
| `point 1 2 3 and` | SelectionError: `Selection failed: 'could not convert string to float: 'and''` |
| `point ( 1 2 3 4 )` | SelectionError: `Selection failed: 'could not convert string to float: '(''` |
| `around 1 2 name P` | SelectionError: `Unknown selection token: '2'` |
| `around -1 name P` | 0 atoms: (none) |
| `around 1e0 name P` | 0 atoms: (none) |
| `around inf name P` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `around nan name P` | 0 atoms: (none) |
| `around 1_0 name P` | 31 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-40 |
| `around +1 name P` | 0 atoms: (none) |
| `around \1 name P` | SelectionError: `Selection failed: 'could not convert string to float: '\\1''` |
| `sphlayer 1 name P` | SelectionError: `Selection failed: 'could not convert string to float: 'name''` |
| `sphlayer 1 2` | SelectionError: `Unknown selection token: 'None'` |
| `isolayer 1 2` | SelectionError: `Unknown selection token: 'None'` |
| `cyzone 1 2 3` | SelectionError: `Unknown selection token: 'None'` |
| `cylayer 1 2 3 4` | SelectionError: `Unknown selection token: 'None'` |
| `cylayer 1 2 3 x all` | SelectionError: `Selection failed: 'could not convert string to float: 'x''` |
| `atom A 1 P Q` | SelectionError: `Unexpected token at end of selection string: 'Q'` |
| `atom A 1.0 P` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: '1.0''` |
| `atom A +1 P` | 1 atoms: 1 |
| `atom A 1_0 P` | 0 atoms: (none) |
| `atom A 01 P` | 1 atoms: 1 |
| `atom A -0 P` | 0 atoms: (none) |
| `atom A \1 P` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: '\\1''` |
| `atom A 1 and` | 0 atoms: (none) |
| `atom A 1 (` | 0 atoms: (none) |
| `atom and 1 P` | 0 atoms: (none) |
| `atom A 1 \P` | 0 atoms: (none) |
| `same resid as` | SelectionError: `Unknown selection token: 'None'` |
| `same resid` | SelectionError: `Unexpected token: 'None' Expected: 'as'` |
| `same` | IndexError: `deque index out of range` |
| `same as name P` | SelectionError: `Selection failed: 'Unknown same property : asChoose one of : dict_keys([ …SAME-LIST… ])'` |
| `same resid of name P` | SelectionError: `Unexpected token: 'of' Expected: 'as'` |
| `same resid AS name P` | SelectionError: `Unexpected token: 'AS' Expected: 'as'` |
| `same RESID as name P` | SelectionError: `Selection failed: 'Unknown same property : RESIDChoose one of : dict_keys([ …SAME-LIST… ])'` |
| `same \resid as name P` | SelectionError: `Selection failed: 'Unknown same property : \residChoose one of : dict_keys([ …SAME-LIST… ])'` |
| `same None as name P` | SelectionError: `Selection failed: 'Unknown same property : NoneChoose one of : dict_keys([ …SAME-LIST… ])'` |
| `same residue as name P and` | SelectionError: `Unknown selection token: 'None'` |
| `same fragments as all` | SelectionError: `Selection failed: 'Unknown same property : fragmentsChoose one of : dict_keys([ …SAME-LIST… ])'` |
| `same residues as all` | SelectionError: `Selection failed: 'Unknown same property : residuesChoose one of : dict_keys([ …SAME-LIST… ])'` |
| `aromaticity maybe` | SelectionError: `Selection failed: ''maybe' is an invalid value for boolean selection. Use 'True' or 'False''` |
| `aromaticity True TRUE tRuE` | NoDataError: `This Universe does not contain aromaticity information` |
| `aromaticity 1` | SelectionError: `Selection failed: ''1' is an invalid value for boolean selection. Use 'True' or 'False''` |
| `chirality` | SelectionError: `Selection failed: 'Unexpected token 'None''` |
| `chirality and` | SelectionError: `Selection failed: 'Unexpected token 'and''` |
| `resid 1 resid 2` | SelectionError: `Unexpected token at end of selection string: 'resid'` |
| `resid 1 or 2` | SelectionError: `Unknown selection token: '2'` |
| `name P or N1` | SelectionError: `Unknown selection token: 'N1'` |
| `index 0 or 1` | SelectionError: `Unknown selection token: '1'` |
| `foo` | SelectionError: `Unknown selection token: 'foo'` |
| `P` | SelectionError: `Unknown selection token: 'P'` |
| `1` | SelectionError: `Unknown selection token: '1'` |
| `\name P` | SelectionError: `Unknown selection token: '\name'` |
| `NOT name P` | SelectionError: `Unknown selection token: 'NOT'` |
| `Not name P` | SelectionError: `Unknown selection token: 'Not'` |
| `name P AND resid 1` | SelectionError: `Unexpected token at end of selection string: 'resid'` |
| `name P && resid 1` | SelectionError: `Unexpected token at end of selection string: 'resid'` |
| `name P \| name N1` | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `!name P` | SelectionError: `Unknown selection token: '!name'` |
| `name P and !name N1` | SelectionError: `Unknown selection token: '!name'` |
| `name P and NOT resid 1` | SelectionError: `Unknown selection token: 'NOT'` |

### 4.4 `smarts`

`smarts` collects raw tokens into a pattern:

* collection stops at the first token that is a keyword, `and`, `or` or
  END — **except** that `(` and `)` are collected; a `)` is collected only
  while it closes a `(` collected by this same `smarts` (a `)` that would make
  the count of collected `)` exceed the count of collected `(` stops the
  collection and is left for the parser);
* no backslash removal;
* the collected tokens are concatenated **without** separators (the
  whitespace and the spaces inserted around parentheses are lost:
  `C ( C ) O` → `C(C)O`); zero tokens give an empty pattern, which is not a
  parse error.

The reference installation has no RDKit, so every `smarts` primary that is
actually evaluated raises ImportError
`RDKit is required for SMARTS-based atom selection but it's not installed. Try installing it with \nconda install -c conda-forge rdkit`
(even on an empty searched group). Pattern validity is therefore never
checked in the reference. Parsing behaviour:

| Selection string | pattern text handed to the SMARTS engine | what the parser sees next / outcome |
|---|---|---|
| `smarts` | (empty) | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `smarts C` | `C` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `smarts C ( C ) O and name P` | `C(C)O` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `(smarts C(C))` | `C(C)` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `smarts [#6](=O)O` | `[#6](=O)O` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `smarts C)` | – | SelectionError: `Unexpected token at end of selection string: ')'` |
| `smarts C ( ( O ) ) ) and all` | – | SelectionError: `Unexpected token at end of selection string: ')'` |
| `smarts \name P` | `\nameP` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `smarts c1ccccc1 name P` | – | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `smarts C not O` | – | SelectionError: `Unexpected token at end of selection string: 'not'` |
| `smarts C ( and )` | – | SelectionError: `Unknown selection token: ')'` |
| `smarts ( C` | `(C` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `smarts ( name P )` | – | SelectionError: `Unexpected token at end of selection string: 'name'` |
| `smarts 1 and name P` | `1` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |
| `not smarts C or name P` | `C` | parses; evaluation raises ImportError (RDKit missing in the reference installation) |

---

## 5. Error classification (supplement to [SV §4.3])

### 5.1 The wrapping rule

While a keyword reads its own arguments, any failure of the *value* kind
(Python ValueError family: no values after a value keyword, failed number /
value parsing, invalid boolean, failed float()/int() conversion of a raw
argument, unknown group name, unknown `same` property, invalid `prop`
operator, "too many values to unpack") is re-raised as SelectionError with
the text `Selection failed: '<original message>'`. Exactly one level of
wrapping happens, at the keyword where the failure occurred; enclosing
keywords (`not`, `byres`, parentheses, …) do not wrap again.

Not wrapped (raised as-is with their own class): every SelectionError
produced by the parser itself (`Unknown selection token`,
`Unexpected token …`, `Invalid expression given: …`, `prop`'s
`Unknown selection token: 'prop'`); TypeError, IndexError and
AttributeError from consuming END (§4.2) or from `prop` ([SV §3.6]);
`group updating` TypeError; RecursionError.

### 5.2 Evaluation-time errors (never wrapped)

| situation | class: message |
|---|---|
| numeric / boolean / single-character keyword, attribute absent | NoDataError: `This Universe does not contain <singular> information` (`bfactor` reports `tempfactor`; `aromatic` reports `aromaticity`) |
| string keyword or macro, attribute absent | AttributeError: `'Topology' object has no attribute '<plural>'` (plural = `names`, `types`, `elements`, `record_types`, `chainIDs`, `altLocs`, `resnames`, `icodes`, `moltypes`, `segids`; macros need `resnames`, and the name-based macros check `names` first) |
| `same P as` with attribute absent | NoDataError `This Universe does not contain <P> information`; `same fragment as` without bonds: NoDataError `AtomGroup.fragments not available; this requires bonds` |
| `bonded` with no bond attribute | NoDataError `This Universe does not contain bonds information` |
| `atom` with resid / segid absent | NoDataError `This Universe does not contain resid information` / `… segid information` (only reached if the earlier filters left atoms, §6.8) |
| `prop` | [SV §3.2 item 7] |
| `resid` with insertion codes on a topology without them | ValueError `Selection specified icodes, while the topology doesn't have any.` |
| `smarts` | ImportError (§4.4) |
| cylinder too large for the box | NotImplementedError, see `selection-grid-method.md` §6 |
| grid-search failures | see `selection-grid-method.md` §4 |
| `same bonds as …` on a system with bonds | IndexError `too many indices for array: array is 1-dimensional, but 2 were indexed` |

NoDataError is both a ValueError and an AttributeError in MDAnalysis, but
since it only arises during evaluation it is never wrapped.

**Empty searched group.** For string keywords, the single-character keyword
(`chirality`), and every distance keyword (`around`, `sphzone`, `sphlayer`,
`isolayer`, `cyzone`, `cylayer`, `point`), evaluation on an empty searched
group returns empty *without* touching the attribute, so a missing attribute
is not reported. The macros (`protein`, `backbone`, …), `atom`, `resid`, the
integer and float keywords, `aromaticity`/`aromatic`, `bonded`, `prop` and
`smarts` still report it. `same …` and `byres` with an empty operand result
return empty without touching their own attribute. An operand's error always
comes before the enclosing keyword's own attribute check (`bonded charge 0`,
`byres charge 0`, `same fragment as charge 0`, `around 1 charge 0` all report
the missing charges on the fixture).

| Selection string | whole fixture (60 atoms) | empty group taken from the fixture |
|---|---|---|
| `name X` | 0 atoms: (none) | 0 atoms: (none) |
| `type X` | 0 atoms: (none) | 0 atoms: (none) |
| `element X` | 0 atoms: (none) | 0 atoms: (none) |
| `record_type ATOM` | 48 atoms: 1-48 | 0 atoms: (none) |
| `chainID A` | 23 atoms: 1-23 | 0 atoms: (none) |
| `altLoc A` | 0 atoms: (none) | 0 atoms: (none) |
| `resname X` | 0 atoms: (none) | 0 atoms: (none) |
| `icode A` | 4 atoms: 18,19,30,31 | 0 atoms: (none) |
| `moltype X` | AttributeError: `'Topology' object has no attribute 'moltypes'` | 0 atoms: (none) |
| `segid A` | 23 atoms: 1-23 | 0 atoms: (none) |
| `index 0` | 1 atoms: 1 | 0 atoms: (none) |
| `resindex 0` | 8 atoms: 1-8 | 0 atoms: (none) |
| `segindex 0` | 23 atoms: 1-23 | 0 atoms: (none) |
| `id 1` | 1 atoms: 1 | 0 atoms: (none) |
| `bynum 1` | 1 atoms: 1 | 0 atoms: (none) |
| `formalcharge 0` | NoDataError: `This Universe does not contain formalcharge information` | NoDataError: `This Universe does not contain formalcharge information` |
| `nbindex 0` | NoDataError: `This Universe does not contain nbindex information` | NoDataError: `This Universe does not contain nbindex information` |
| `resnum 1` | 8 atoms: 1-8 | 0 atoms: (none) |
| `molnum 0` | NoDataError: `This Universe does not contain molnum information` | NoDataError: `This Universe does not contain molnum information` |
| `model 0` | NoDataError: `This Universe does not contain model information` | NoDataError: `This Universe does not contain model information` |
| `radius 1` | NoDataError: `This Universe does not contain radius information` | NoDataError: `This Universe does not contain radius information` |
| `tempfactor 0` | 60 atoms: 1-60 (+float-equality warning) | 0 atoms: (none) (+float-equality warning) |
| `bfactor 0` | 60 atoms: 1-60 (+2× bfactor-alias warning, float-equality warning) | 0 atoms: (none) (+2× bfactor-alias warning, float-equality warning) |
| `mass 1` | 0 atoms: (none) (+float-equality warning) | 0 atoms: (none) (+float-equality warning) |
| `charge 0` | NoDataError: `This Universe does not contain charge information` | NoDataError: `This Universe does not contain charge information` |
| `occupancy 1` | 60 atoms: 1-60 (+float-equality warning) | 0 atoms: (none) (+float-equality warning) |
| `gbscreen 0` | NoDataError: `This Universe does not contain gbscreen information` | NoDataError: `This Universe does not contain gbscreen information` |
| `solventradius 0` | NoDataError: `This Universe does not contain solventradius information` | NoDataError: `This Universe does not contain solventradius information` |
| `rmin 0` | NoDataError: `This Universe does not contain rmin information` | NoDataError: `This Universe does not contain rmin information` |
| `epsilon 0` | NoDataError: `This Universe does not contain epsilon information` | NoDataError: `This Universe does not contain epsilon information` |
| `rmin14 0` | NoDataError: `This Universe does not contain rmin14 information` | NoDataError: `This Universe does not contain rmin14 information` |
| `epsilon14 0` | NoDataError: `This Universe does not contain epsilon14 information` | NoDataError: `This Universe does not contain epsilon14 information` |
| `resid 1` | 8 atoms: 1-8 | 0 atoms: (none) |
| `aromaticity` | NoDataError: `This Universe does not contain aromaticity information` | NoDataError: `This Universe does not contain aromaticity information` |
| `aromatic` | NoDataError: `This Universe does not contain aromaticity information` | NoDataError: `This Universe does not contain aromaticity information` |
| `chirality R` | NoDataError: `This Universe does not contain chirality information` | 0 atoms: (none) |
| `protein` | 15 atoms: 34-48 | 0 atoms: (none) |
| `nucleic` | 33 atoms: 1-33 | 0 atoms: (none) |
| `water` | 6 atoms: 55-60 | 0 atoms: (none) |
| `backbone` | 6 atoms: 34-37,41,42 | 0 atoms: (none) |
| `nucleicbackbone` | 9 atoms: 1,9,14,18,20,24,27,30,32 | 0 atoms: (none) |
| `nucleicbase` | 11 atoms: 5,6,8,11-13,16,17,22,26,29 | 0 atoms: (none) |
| `nucleicsugar` | 8 atoms: 3,10,15,19,25,28,31,33 | 0 atoms: (none) |
| `all` | 60 atoms: 1-60 | 0 atoms: (none) |
| `bonded all` | NoDataError: `This Universe does not contain bonds information` | NoDataError: `This Universe does not contain bonds information` |
| `same fragment as all` | NoDataError: `AtomGroup.fragments not available; this requires bonds` | 0 atoms: (none) |
| `same charge as all` | NoDataError: `This Universe does not contain charge information` | 0 atoms: (none) |
| `same radius as all` | NoDataError: `This Universe does not contain radius information` | 0 atoms: (none) |
| `same bfactor as all` | 60 atoms: 1-60 (+2× bfactor-alias warning) | 0 atoms: (none) |
| `same x as all` | 60 atoms: 1-60 | 0 atoms: (none) |
| `byres all` | 60 atoms: 1-60 | 0 atoms: (none) |
| `around 1 all` | 0 atoms: (none) | 0 atoms: (none) |
| `sphzone 1 all` | 2 atoms: 30,31 | 0 atoms: (none) |
| `point 0 0 0 1` | 0 atoms: (none) | 0 atoms: (none) |
| `atom A 1 P` | 1 atoms: 1 | 0 atoms: (none) |
| `prop charge > 0` | SelectionError: `This Universe does not contain charges information` | SelectionError: `This Universe does not contain charges information` |
| `prop x > 0` | 60 atoms: 1-60 | 0 atoms: (none) |
| `smarts C` | ImportError: `RDKit is required for SMARTS-based atom selection but it's not installed. Try installing it with \nconda install -c conda-forge rdkit` | ImportError: `RDKit is required for SMARTS-based atom selection but it's not installed. Try installing it with \nconda install -c conda-forge rdkit` |
| `not all` | 0 atoms: (none) | 0 atoms: (none) |
| `global all` | 60 atoms: 1-60 | 60 atoms: 1-60 |
| `cyzone 1 1 -1 all` | 2 atoms: 30,31 | 0 atoms: (none) |

A system with no attributes at all except the built-in index, residue and
segment numbering (3 atoms, 1 residue, 1 segment):

| Selection string | MDAnalysis result |
|---|---|
| `name X` | AttributeError: `'Topology' object has no attribute 'names'` |
| `type X` | AttributeError: `'Topology' object has no attribute 'types'` |
| `element X` | AttributeError: `'Topology' object has no attribute 'elements'` |
| `record_type ATOM` | AttributeError: `'Topology' object has no attribute 'record_types'` |
| `chainID A` | AttributeError: `'Topology' object has no attribute 'chainIDs'` |
| `chainid A` | AttributeError: `'Topology' object has no attribute 'chainIDs'` |
| `altLoc A` | AttributeError: `'Topology' object has no attribute 'altLocs'` |
| `resname X` | AttributeError: `'Topology' object has no attribute 'resnames'` |
| `icode A` | AttributeError: `'Topology' object has no attribute 'icodes'` |
| `moltype X` | AttributeError: `'Topology' object has no attribute 'moltypes'` |
| `segid A` | AttributeError: `'Topology' object has no attribute 'segids'` |
| `index 0` | 1 atoms: 1 |
| `resindex 0` | 3 atoms: 1-3 |
| `segindex 0` | 3 atoms: 1-3 |
| `id 1` | NoDataError: `This Universe does not contain id information` |
| `bynum 1` | 1 atoms: 1 |
| `resnum 1` | NoDataError: `This Universe does not contain resnum information` |
| `resid 1` | NoDataError: `This Universe does not contain resid information` |
| `mass 1` | NoDataError: `This Universe does not contain mass information` |
| `mass 0:1` | NoDataError: `This Universe does not contain mass information` |
| `tempfactor 0` | NoDataError: `This Universe does not contain tempfactor information` |
| `bfactor 0` | NoDataError: `This Universe does not contain tempfactor information` |
| `protein` | AttributeError: `'Topology' object has no attribute 'resnames'` |
| `nucleic` | AttributeError: `'Topology' object has no attribute 'resnames'` |
| `water` | AttributeError: `'Topology' object has no attribute 'resnames'` |
| `backbone` | AttributeError: `'Topology' object has no attribute 'names'` |
| `nucleicbackbone` | AttributeError: `'Topology' object has no attribute 'names'` |
| `nucleicbase` | AttributeError: `'Topology' object has no attribute 'names'` |
| `nucleicsugar` | AttributeError: `'Topology' object has no attribute 'names'` |
| `atom A 1 P` | NoDataError: `This Universe does not contain name information` |
| `same name as all` | NoDataError: `This Universe does not contain name information` |
| `same resid as all` | NoDataError: `This Universe does not contain resid information` |
| `same resnum as all` | NoDataError: `This Universe does not contain resnum information` |
| `same segid as all` | NoDataError: `This Universe does not contain segid information` |
| `same residue as all` | 3 atoms: 1-3 |
| `same segment as all` | 3 atoms: 1-3 |
| `same type as all` | NoDataError: `This Universe does not contain type information` |
| `same resname as all` | NoDataError: `This Universe does not contain resname information` |
| `same mass as all` | NoDataError: `This Universe does not contain mass information` |
| `prop mass > 0` | SelectionError: `This Universe does not contain masses information` |
| `prop name == 1` | SelectionError: `This Universe does not contain names information` |
| `prop id > 0` | SelectionError: `This Universe does not contain ids information` |

The same attribute-less system, empty searched group:

| Selection string | MDAnalysis result |
|---|---|
| `protein` | AttributeError: `'Topology' object has no attribute 'resnames'` |
| `water` | AttributeError: `'Topology' object has no attribute 'resnames'` |
| `backbone` | AttributeError: `'Topology' object has no attribute 'names'` |
| `nucleicsugar` | AttributeError: `'Topology' object has no attribute 'names'` |
| `name X` | 0 atoms: (none) |
| `resname X` | 0 atoms: (none) |
| `chirality R` | 0 atoms: (none) |
| `atom A 1 P` | NoDataError: `This Universe does not contain name information` |
| `resid 1` | NoDataError: `This Universe does not contain resid information` |
| `mass 1` | NoDataError: `This Universe does not contain mass information` |
| `id 1` | NoDataError: `This Universe does not contain id information` |
| `index 0` | 0 atoms: (none) |
| `same name as all` | 0 atoms: (none) |
| `byres all` | 0 atoms: (none) |
| `around 1 all` | 0 atoms: (none) |
| `point 0 0 0 1` | 0 atoms: (none) |
| `prop mass > 1` | SelectionError: `This Universe does not contain masses information` |
| `aromatic` | NoDataError: `This Universe does not contain aromaticity information` |

Empty-operand short cuts (whole fixture):

| Selection string | MDAnalysis result |
|---|---|
| `same charge as not all` | 0 atoms: (none) |
| `same charge as name XX` | 0 atoms: (none) |
| `same fragment as not all` | 0 atoms: (none) |
| `same radius as name XX` | 0 atoms: (none) |
| `bonded not all` | NoDataError: `This Universe does not contain bonds information` |
| `byres not all and charge 0` | NoDataError: `This Universe does not contain charge information` |
| `not all and charge 0` | NoDataError: `This Universe does not contain charge information` |
| `charge 0 and not all` | NoDataError: `This Universe does not contain charge information` |

---

## 6. Evaluation semantics

### 6.1 Searched group, `and`, `or`, `not`, `global`

* Every keyword selects only from the searched group, **except**: `global X`
  evaluates X on the whole system and returns that result unrestricted;
  `bonded` returns bond partners anywhere in the system (§6.7); `or` returns
  the union of its operands' results (so a `global` operand of `or` can bring
  in outside atoms).
* `A and B` = atoms in both results. `A or B` = atoms in either. `not A` =
  atoms of the searched group not in A's result.
* `global` does not change what the surrounding keywords search: in
  `around 2 global X` the reference atoms come from the whole system but the
  candidate atoms are still the searched group; in `byres global X` the
  residues come from X but the returned atoms are those of the searched group.
* `index`, `resindex`, `segindex`, `bynum` always refer to the system-wide
  numbering, not to positions inside the searched group.

| Selection string | called on the whole fixture (60 atoms) | called on the group `resid 1 2` (serials 1–13) | called on the group `chainID C` (serials 34–48) |
|---|---|---|---|
| `all` | 60 atoms: 1-60 | 13 atoms: 1-13 | 15 atoms: 34-48 |
| `name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 | 2 atoms: 1,9 | 0 atoms: (none) |
| `global name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 | 9 atoms: 1,9,14,18,20,24,27,30,32 | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `not name P` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 | 11 atoms: 2-8,10-13 | 15 atoms: 34-48 |
| `global not name P` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 |
| `not global name P` | 51 atoms: 2-8,10-13,15-17,19,21-23,25,26,28,29,31,33-60 | 11 atoms: 2-8,10-13 | 15 atoms: 34-48 |
| `name P or global name N1` | 13 atoms: 1,6,9,11,14,16,18,20,24,26,27,30,32 | 6 atoms: 1,6,9,11,16,26 | 4 atoms: 6,11,16,26 |
| `global name P and resid 1` | 1 atoms: 1 | 1 atoms: 1 | 0 atoms: (none) |
| `global (name P or name N1) and not resid 1` | 11 atoms: 9,11,14,16,18,20,24,26,27,30,32 | 2 atoms: 9,11 | 0 atoms: (none) |
| `global global name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 | 9 atoms: 1,9,14,18,20,24,27,30,32 | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `index 0 20` | 2 atoms: 1,21 | 1 atoms: 1 | 0 atoms: (none) |
| `bynum 1:100` | 60 atoms: 1-60 | 13 atoms: 1-13 | 15 atoms: 34-48 |
| `resid 1 to 100` | 48 atoms: 1-48 | 13 atoms: 1-13 | 15 atoms: 34-48 |
| `byres name N3` | 12 atoms: 9-13,20-23,27-29 | 5 atoms: 9-13 | 0 atoms: (none) |
| `byres global name N3` | 12 atoms: 9-13,20-23,27-29 | 5 atoms: 9-13 | 0 atoms: (none) |
| `byres global resid 3` | 4 atoms: 14-17 | 0 atoms: (none) | 0 atoms: (none) |
| `same resname as global name P` | 33 atoms: 1-33 | 13 atoms: 1-13 | 0 atoms: (none) |
| `same resid as global resid 3` | 6 atoms: 14-19 | 0 atoms: (none) | 0 atoms: (none) |
| `around 1.3 name N9` | 1 atoms: 6 | 1 atoms: 6 | 0 atoms: (none) |
| `around 1.3 global name H1` | 3 atoms: 8,57,58 | 1 atoms: 8 | 0 atoms: (none) |
| `around 2 global resid 3` | 2 atoms: 13,18 | 1 atoms: 13 | 0 atoms: (none) |
| `sphzone 2 global name N9` | 3 atoms: 4-6 | 3 atoms: 4-6 | 0 atoms: (none) |
| `point 6 0 1 1.5` | 2 atoms: 5,6 | 2 atoms: 5,6 | 0 atoms: (none) |
| `point 18 0 1 1.5` | 2 atoms: 15,16 | 0 atoms: (none) | 0 atoms: (none) |
| `atom A 3 P` | 2 atoms: 14,18 | 0 atoms: (none) | 0 atoms: (none) |
| `global atom A 3 P` | 2 atoms: 14,18 | 2 atoms: 14,18 | 2 atoms: 14,18 |
| `prop x < 3` | 2 atoms: 1,2 | 2 atoms: 1,2 | 0 atoms: (none) |
| `global prop x < 3` | 2 atoms: 1,2 | 2 atoms: 1,2 | 2 atoms: 1,2 |
| `nucleic` | 33 atoms: 1-33 | 13 atoms: 1-13 | 0 atoms: (none) |
| `global protein` | 15 atoms: 34-48 | 15 atoms: 34-48 | 15 atoms: 34-48 |
| `protein` | 15 atoms: 34-48 | 0 atoms: (none) | 15 atoms: 34-48 |

### 6.2 Numeric and float keywords (supplement to [SV §2.3, §2.5])

* Integer keywords compare exactly; `bynum` compares index + 1; `id` compares
  the stored id (which need not be index + 1). Ranges inclusive; reversed
  range selects nothing.
* Float keywords: ranges inclusive, compared in double precision against the
  stored value (a reader may store single-precision-rounded numbers in it;
  e.g. parm7 charges, so a typed `-0.4177` will not equal a stored
  −0.4177058935165405). A single value V selects atoms whose value v satisfies
  **v = V exactly, or (V finite and |v − V| ≤ atol + rtol·|V|)** with the
  call's atol/rtol (defaults 1e-8, 1e-5). The relative term uses the
  *typed* value. With negative tolerances only exact equality matches. NaN
  never matches. The float-equality warning is emitted once per single value
  per evaluation, even when the value is repeated.
* atol/rtol also apply to `prop … == …`, but **not** to `same x/y/z as`
  (fixed 1e-8 / 1e-5) and not to ranges.

| Selection string | MDAnalysis result |
|---|---|
| `charge -0.3` | 2 atoms: 1,5 (+float-equality warning) |
| `charge -0.30000001` | 2 atoms: 1,5 (+float-equality warning) |
| `charge -0.300003` | 2 atoms: 1,5 (+float-equality warning) |
| `charge -0.3000031` | 0 atoms: (none) (+float-equality warning) |
| `charge -0.29999` | 0 atoms: (none) (+float-equality warning) |
| `charge 0.417` | 2 atoms: 10,11 (+float-equality warning) |
| `charge 0.4170041` | 2 atoms: 10,11 (+float-equality warning) |
| `charge 0.4170042` | 0 atoms: (none) (+float-equality warning) |
| `charge 0` | 0 atoms: (none) (+float-equality warning) |
| `charge -0.5:0` | 4 atoms: 1,4,5,8 |
| `charge 0:-0.5` | 0 atoms: (none) |
| `charge -0.5 - 0` | 4 atoms: 1,4,5,8 |
| `charge -0.5-0` | 4 atoms: 1,4,5,8 |
| `charge -0.5 -0` | 2 atoms: 4,8 (+2× float-equality warning) |
| `charge -0.5--0.3` | 4 atoms: 1,4,5,8 |
| `charge -1:-0.5` | 3 atoms: 4,8,9 |
| `charge -0.834 1` | 2 atoms: 9,12 (+2× float-equality warning) |
| `charge -0.834 0.4:2` | 6 atoms: 3,7,9-12 (+float-equality warning) |
| `charge 0.1 0.1 0.1` | 2 atoms: 2,6 (+3× float-equality warning) |
| `charge 1 to 1` | 1 atoms: 12 |
| `charge 1.0` | 1 atoms: 12 (+float-equality warning) |
| `charge 1e0` | 1 atoms: 12 (+float-equality warning) |
| `charge +1` | 1 atoms: 12 (+float-equality warning) |
| `charge 1.` | 1 atoms: 12 (+float-equality warning) |
| `charge .417` | 2 atoms: 10,11 (+float-equality warning) |
| `charge 0.4 to 0.42 -0.3` | 4 atoms: 1,5,10,11 (+float-equality warning) |
| `formalcharge 1` | 1 atoms: 12 |
| `formalcharge 0` | 11 atoms: 1-11 |
| `formalcharge -1:0` | 11 atoms: 1-11 |
| `formalcharge 1:0` | 0 atoms: (none) |
| `formalcharge 0 1` | 12 atoms: 1-12 |
| `formalcharge 1.0` | SelectionError: `Selection failed: 'Failed to parse number: 1.0'` |
| `formalcharge +1` | 1 atoms: 12 |
| `mass 12.011` | 4 atoms: 2,3,6,7 (+float-equality warning) |
| `mass 12.0110001` | 4 atoms: 2,3,6,7 (+float-equality warning) |
| `mass 12.0111` | 4 atoms: 2,3,6,7 (+float-equality warning) |
| `mass 12.0112` | 0 atoms: (none) (+float-equality warning) |
| `mass 22.99` | 1 atoms: 12 (+float-equality warning) |
| `mass 1:15` | 8 atoms: 1-3,5-7,10,11 |
| `mass 15.999 to 22.99` | 4 atoms: 4,8,9,12 |
| `mass 1.008 22.99` | 3 atoms: 10-12 (+2× float-equality warning) |
| `tempfactor 0.3` | 1 atoms: 3 (+float-equality warning) |
| `tempfactor 0.1:0.3` | 3 atoms: 1-3 |
| `tempfactor 0.30000001` | 1 atoms: 3 (+float-equality warning) |
| `tempfactor 0.3:0.3` | 1 atoms: 3 |
| `tempfactor 1.2:inf` | SelectionError: `Selection failed: 'could not convert string to float: '''` |
| `occupancy 1` | NoDataError: `This Universe does not contain occupancy information` |
| `index 0 11` | 2 atoms: 1,12 |
| `index 11:0` | 0 atoms: (none) |
| `index 3-5 8` | 4 atoms: 4-6,9 |
| `index 3 - 5 8` | 4 atoms: 4-6,9 |
| `index -1` | 0 atoms: (none) |
| `index 12` | 0 atoms: (none) |
| `bynum 1 12` | 2 atoms: 1,12 |
| `bynum 0:1` | 1 atoms: 1 |
| `bynum 12:13` | 1 atoms: 12 |
| `id 1 12` | 2 atoms: 1,12 |
| `id 0:1` | 1 atoms: 1 |
| `resindex 0` | 4 atoms: 1-4 |
| `resindex 3` | 1 atoms: 12 |
| `resindex 1:2` | 7 atoms: 5-11 |
| `segindex 1` | 4 atoms: 9-12 |
| `segindex 0:1` | 12 atoms: 1-12 |
| `resnum 3` | 3 atoms: 9-11 |
| `resnum 4` | 1 atoms: 12 |
| `resnum 3:4` | 4 atoms: 9-12 |
| `resnum 3B` | SelectionError: `Selection failed: 'Failed to parse number: 3B'` |
| `resid 3` | 3 atoms: 9-11 |
| `resid 3B` | 1 atoms: 12 |
| `resid 3 3B` | 4 atoms: 9-12 |
| `resid 3:3B` | 4 atoms: 9-12 |
| `resid 3A:3C` | 1 atoms: 12 |
| `resid 3B:3` | 0 atoms: (none) |
| `resid 2:3` | 7 atoms: 5-11 |
| `resid 2:3A` | 7 atoms: 5-11 |
| `resid 2:3B` | 8 atoms: 5-12 |
| `resid 1 3` | 7 atoms: 1-4,9-11 |
| `resid 3:1` | 8 atoms: 1-4,9-12 |
| `resid 4:2` | 4 atoms: 5-8 |
| `resid 1 1 1` | 4 atoms: 1-4 |
| `resid 2 1:1` | 8 atoms: 1-8 |
| `resid -5:1` | 4 atoms: 1-4 |
| `resid 3b` | 0 atoms: (none) |

With non-default tolerances:

| Selection string | keyword arguments | MDAnalysis result |
|---|---|---|
| `charge -0.3` | atol=0, rtol=0 | 2 atoms: 1,5 (+float-equality warning) |
| `charge -0.30000001` | atol=0, rtol=0 | 0 atoms: (none) (+float-equality warning) |
| `charge -0.3` | atol=0.2 | 4 atoms: 1,4,5,8 (+float-equality warning) |
| `charge -0.3` | atol=0.21 | 4 atoms: 1,4,5,8 (+float-equality warning) |
| `charge 0` | atol=0.1 | 2 atoms: 2,6 (+float-equality warning) |
| `charge 0` | atol=0.1000001 | 2 atoms: 2,6 (+float-equality warning) |
| `charge 0.4` | atol=0, rtol=0.05 | 2 atoms: 10,11 (+float-equality warning) |
| `charge 0.4` | atol=0, rtol=0.03 | 0 atoms: (none) (+float-equality warning) |
| `mass 12` | atol=0.011 | 4 atoms: 2,3,6,7 (+float-equality warning) |
| `mass 12` | atol=0.0111 | 4 atoms: 2,3,6,7 (+float-equality warning) |
| `mass 12` | rtol=0.001 | 4 atoms: 2,3,6,7 (+float-equality warning) |
| `mass 12:12.011` | atol=1 | 4 atoms: 2,3,6,7 |
| `charge -0.5:0` | atol=1, rtol=1 | 4 atoms: 1,4,5,8 |
| `prop charge == 0.4` | atol=0.02 | 2 atoms: 10,11 |
| `prop charge == 0.4` | atol=0.01 | 0 atoms: (none) |
| `prop mass == 12` | atol=0.02 | 4 atoms: 2,3,6,7 |
| `same x as index 0` | atol=10 | 1 atoms: 1 |
| `charge -0.3` | atol=-1 | 2 atoms: 1,5 (+float-equality warning) |

SYN with ids replaced by 115, 105, …, 5 (index 0 has id 115, index 11 has
id 5):

| Selection string | MDAnalysis result |
|---|---|
| `id 5` | 1 atoms: 12 |
| `id 115` | 1 atoms: 1 |
| `id 5:25` | 3 atoms: 10-12 |
| `bynum 1` | 1 atoms: 1 |
| `bynum 12` | 1 atoms: 12 |
| `index 0` | 1 atoms: 1 |
| `index 11` | 1 atoms: 12 |
| `prop id < 20` | 2 atoms: 11,12 |
| `id 25 to 5` | 0 atoms: (none) |

### 6.3 `resid` and `resnum` (supplement to [SV §2.4, §2.7])

`resnum` is a plain integer keyword: it ignores insertion codes and does not
accept them (`resnum 3B` → `Failed to parse number`). `icode` is a string
keyword over the insertion-code text; the empty insertion code can be
matched with a lone backslash value (`icode \`, [SV §1.1]). On a topology
without insertion codes, `icode` raises AttributeError (§5.2) while `resid`
with a code raises the ValueError of [SV §2.4].

| Selection string | fixture as read (has insertion codes) | fixture, insertion-code attribute removed |
|---|---|---|
| `resid 10A` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 10A-12B` | 7 atoms: 27-33 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 10:12` | 10 atoms: 24-33 | 10 atoms: 24-33 |
| `resid 10 - 12` | 10 atoms: 24-33 | 10 atoms: 24-33 |
| `resid 10A:12` | 7 atoms: 27-33 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11A` | 2 atoms: 30,31 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11 11A` | 5 atoms: 27-31 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11A-11A` | 2 atoms: 30,31 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11-11A` | 5 atoms: 27-31 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11B` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11:11` | 3 atoms: 27-29 | 5 atoms: 27-31 |
| `resid 11A:11` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11Z` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 10Z:11` | 3 atoms: 27-29 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11Z:12` | 2 atoms: 32,33 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11_` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 11a` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3A-4` | 6 atoms: 18-23 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 3:3Z` | 6 atoms: 14-19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 2Z:3` | 4 atoms: 14-17 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 10A 3A` | 2 atoms: 18,19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 12:10` | 5 atoms: 24-26,32,33 | 0 atoms: (none) |
| `resid 12A:10` | 3 atoms: 24-26 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 4 : 1` | 12 atoms: 1-8,20-23 | 0 atoms: (none) |
| `resid 0:2` | 13 atoms: 1-13 | 13 atoms: 1-13 |
| `resid 201:206` | 6 atoms: 49-54 | 6 atoms: 49-54 |
| `resid 300 to 400` | 6 atoms: 55-60 | 6 atoms: 55-60 |
| `resid 3 and icode A` | 0 atoms: (none) | AttributeError: `'Topology' object has no attribute 'icodes'` |
| `icode A` | 4 atoms: 18,19,30,31 | AttributeError: `'Topology' object has no attribute 'icodes'` |
| `icode \` | 56 atoms: 1-17,20-29,32-60 | AttributeError: `'Topology' object has no attribute 'icodes'` |
| `not icode \` | 4 atoms: 18,19,30,31 | AttributeError: `'Topology' object has no attribute 'icodes'` |
| `resnum 3` | 6 atoms: 14-19 | 6 atoms: 14-19 |
| `resnum 3:3` | 6 atoms: 14-19 | 6 atoms: 14-19 |
| `resnum 11` | 5 atoms: 27-31 | 5 atoms: 27-31 |
| `resid 11 and not icode A` | 3 atoms: 27-29 | AttributeError: `'Topology' object has no attribute 'icodes'` |
| `resid 3A and resnum 3` | 2 atoms: 18,19 | ValueError: `Selection specified icodes, while the topology doesn't have any.` |
| `resid 1 and resid 1A` | 0 atoms: (none) | ValueError: `Selection specified icodes, while the topology doesn't have any.` |

### 6.4 Missing attributes

Covered by §5.2 and its tables.

### 6.5 String keywords (supplement to [SV §1.2])

Matching is shell-style, case-sensitive, over the whole attribute text:
`*` (any run, also empty), `?` (exactly one character), `[set]`,
`[!set]` (negation only with `!`; `^` is an ordinary set member), ranges
`a-z` inside sets (a reversed range such as `[Z-A]` matches nothing), `]` as
the first set member is literal, and an unclosed `[` is a literal `[`.
Backslash has no escaping role inside a pattern (it is a literal character;
only the single leading backslash of the token is removed, [SV §1.1]).
Residue-level attributes (`resname`, `icode`, `moltype`) and segment-level
ones (`segid`) are matched per residue / segment and broadcast to its atoms.

| Selection string | MDAnalysis result |
|---|---|
| `name *` | 60 atoms: 1-60 |
| `name ?` | 17 atoms: 1,9,14,18,20,24,27,30,32,34,36,37,39,41,54,55,58 |
| `name ??` | 26 atoms: 5-8,11-13,16,17,22,23,26,29,35,38,42,43,48-51,53,56,57,59,60 |
| `name ???` | 16 atoms: 2,3,10,15,19,21,25,28,31,33,40,44-47,52 |
| `name [` | 0 atoms: (none) |
| `name [P` | 0 atoms: (none) |
| `name P]` | 0 atoms: (none) |
| `name []]` | 0 atoms: (none) |
| `name [!]]` | 17 atoms: 1,9,14,18,20,24,27,30,32,34,36,37,39,41,54,55,58 |
| `name [!P]` | 8 atoms: 34,36,37,39,41,54,55,58 |
| `name [A-Z]?` | 26 atoms: 5-8,11-13,16,17,22,23,26,29,35,38,42,43,48-51,53,56,57,59,60 |
| `name [a-z]*` | 0 atoms: (none) |
| `name *[0-9]` | 22 atoms: 2,5-8,11-13,16,17,22,23,26,29,44-47,56,57,59,60 |
| `name C4[*']` | 9 atoms: 3,10,15,19,21,25,28,31,33 |
| `name C4\*` | 0 atoms: (none) |
| `name C4[*]` | 1 atoms: 21 |
| `name *\\*` | 0 atoms: (none) |
| `name H*` | 8 atoms: 4,7,23,39,56,57,59,60 |
| `name *H*` | 9 atoms: 4,7,23,39,40,56,57,59,60 |
| `name *1` | 10 atoms: 2,6,7,11,16,26,44,46,56,59 |
| `name [^C]*` | 20 atoms: 3,10,15,19,21,25,28,31,33,35,36,38,42-48,53 |
| `name [!C-Z]*` | 1 atoms: 40 |
| `name [Z-A]*` | 0 atoms: (none) |
| `name ** *` | 60 atoms: 1-60 |
| `name O[P2]*` | 1 atoms: 2 |
| `resname [GC]` | 21 atoms: 1-13,24-31 |
| `resname ?` | 34 atoms: 1-33,54 |
| `segid [AB]` | 33 atoms: 1-33 |
| `element M?` | 3 atoms: 49-51 |
| `type M*` | 3 atoms: 49-51 |
| `chainID W` | 6 atoms: 55-60 |
| `record_type HETATM` | 12 atoms: 49-60 |
| `record_type ATOM HETATM` | 60 atoms: 1-60 |
| `altLoc *` | 60 atoms: 1-60 |
| `icode *` | 60 atoms: 1-60 |
| `icode ?` | 4 atoms: 18,19,30,31 |
| `icode A B` | 4 atoms: 18,19,30,31 |

### 6.6 Macros

A macro selects atoms of the searched group whose **residue name** is in a
fixed list (exact, case-sensitive) and, for the name-based macros, whose
**atom name** is also in a fixed list (exact, no patterns).

* `protein` — resname in: ACE, ALA, ALAD, ARG, ARGN, ASF, ASH, ASN, ASN1, ASP, ASPH, CALA, CARG, CASF, CASN, CASP, CCYS, CCYX, CGLN, CGLU, CGLY, CHID, CHIE, CHIP, CILE, CLEU, CLYS, CME, CMET, CPHE, CPRO, CSER, CTHR, CTRP, CTYR, CVAL, CYM, CYS, CYS1, CYS2, CYSH, CYX, DAB, GLH, GLN, GLU, GLUH, GLY, HID, HIE, HIP, HIS, HIS1, HIS2, HISA, HISB, HISD, HISE, HISH, HSD, HSE, HSP, HYP, ILE, LEU, LYN, LYS, LYSH, MET, MSE, NALA, NARG, NASN, NASP, NCYS, NCYX, NGLN, NGLU, NGLY, NHID, NHIE, NHIP, NILE, NLEU, NLYS, NME, NMET, NPHE, NPRO, NSER, NTHR, NTRP, NTYR, NVAL, ORN, PGLU, PHE, PRO, QLN, SER, THR, TRP, TYR, VAL (104 names).
* `backbone` — protein resname (list above) **and** atom name in
  {C, CA, N, O}. (Terminal oxygens such as OT1, OXT are not included.)
* `nucleic` — resname in: A, ADE, C, CYT, DA, DA3, DA5, DC, DC3, DC5, DG,
  DG3, DG5, DT, DT3, DT5, G, GUA, RA, RA3, RA5, RC, RC3, RC5, RG, RG3, RG5,
  RU, RU3, RU5, T, THY, U, URA (34 names).
* `nucleicbackbone` — nucleic resname and atom name in {C3', C5', O3', O5', P}.
* `nucleicbase` — nucleic resname and atom name in {C2, C4, C5, C5M, C6, C8,
  N1, N2, N3, N4, N6, N7, N9, O2, O4, O6}.
* `nucleicsugar` — nucleic resname and atom name in {C1', C2', C3', C4', O4'}
  (note: O2' and O3' are **not** sugar atoms; C3' is both backbone and sugar).
* `water` — resname in: H2O, HHO, HOH, OH2, OHH, SOL, T3P, T4P, T5P, TIP,
  TIP2, TIP3, TIP4, WAT (14 names).
* `all` — the searched group.

Residue-name test on a synthetic system (one residue per row, three atoms
named N, CA, P in each):

| resname | atoms (serial: name) | `protein` | `backbone` | `nucleic` | `nucleicbackbone` | `water` |
|---|---|---|---|---|---|---|
| `ALA` | 1:N 2:CA 3:P | 1,2,3 | 1,2 | – | – | – |
| `ala` | 4:N 5:CA 6:P | – | – | – | – | – |
| `HIS` | 7:N 8:CA 9:P | 7,8,9 | 7,8 | – | – | – |
| `HSD` | 10:N 11:CA 12:P | 10,11,12 | 10,11 | – | – | – |
| `HIE` | 13:N 14:CA 15:P | 13,14,15 | 13,14 | – | – | – |
| `ACE` | 16:N 17:CA 18:P | 16,17,18 | 16,17 | – | – | – |
| `NME` | 19:N 20:CA 21:P | 19,20,21 | 19,20 | – | – | – |
| `CME` | 22:N 23:CA 24:P | 22,23,24 | 22,23 | – | – | – |
| `ASF` | 25:N 26:CA 27:P | 25,26,27 | 25,26 | – | – | – |
| `ALAD` | 28:N 29:CA 30:P | 28,29,30 | 28,29 | – | – | – |
| `MSE` | 31:N 32:CA 33:P | 31,32,33 | 31,32 | – | – | – |
| `NALA` | 34:N 35:CA 36:P | 34,35,36 | 34,35 | – | – | – |
| `CALA` | 37:N 38:CA 39:P | 37,38,39 | 37,38 | – | – | – |
| `HOH` | 40:N 41:CA 42:P | – | – | – | – | 40,41,42 |
| `WAT` | 43:N 44:CA 45:P | – | – | – | – | 43,44,45 |
| `SOL` | 46:N 47:CA 48:P | – | – | – | – | 46,47,48 |
| `TIP3` | 49:N 50:CA 51:P | – | – | – | – | 49,50,51 |
| `TIP4` | 52:N 53:CA 54:P | – | – | – | – | 52,53,54 |
| `TIP5` | 55:N 56:CA 57:P | – | – | – | – | – |
| `T3P` | 58:N 59:CA 60:P | – | – | – | – | 58,59,60 |
| `H2O` | 61:N 62:CA 63:P | – | – | – | – | 61,62,63 |
| `DA` | 64:N 65:CA 66:P | – | – | 64,65,66 | 66 | – |
| `DA5` | 67:N 68:CA 69:P | – | – | 67,68,69 | 69 | – |
| `RU3` | 70:N 71:CA 72:P | – | – | 70,71,72 | 72 | – |
| `A` | 73:N 74:CA 75:P | – | – | 73,74,75 | 75 | – |
| `G` | 76:N 77:CA 78:P | – | – | 76,77,78 | 78 | – |
| `URA` | 79:N 80:CA 81:P | – | – | 79,80,81 | 81 | – |
| `ADE` | 82:N 83:CA 84:P | – | – | 82,83,84 | 84 | – |
| `DT` | 85:N 86:CA 87:P | – | – | 85,86,87 | 87 | – |
| `LIG` | 88:N 89:CA 90:P | – | – | – | – | – |
| `Na+` | 91:N 92:CA 93:P | – | – | – | – | – |
| `PRO` | 94:N 95:CA 96:P | 94,95,96 | 94,95 | – | – | – |
| `GLY` | 97:N 98:CA 99:P | 97,98,99 | 97,98 | – | – | – |

Atom-name test of the nucleic macros (residue `G`, and the same names in a
residue `LIG`):

| atom name | `nucleicbackbone` | `nucleicbase` | `nucleicsugar` | (same name in resname `LIG`: any macro) |
|---|---|---|---|---|
| `P` | yes | – | – | – |
| `C5'` | yes | – | – | – |
| `C3'` | yes | – | yes | – |
| `O3'` | yes | – | – | – |
| `O5'` | yes | – | – | – |
| `C1'` | – | – | yes | – |
| `C2'` | – | – | yes | – |
| `C4'` | – | – | yes | – |
| `O4'` | – | – | yes | – |
| `O2'` | – | – | – | – |
| `N9` | – | yes | – | – |
| `N7` | – | yes | – | – |
| `C8` | – | yes | – | – |
| `C5` | – | yes | – | – |
| `C4` | – | yes | – | – |
| `N3` | – | yes | – | – |
| `C2` | – | yes | – | – |
| `N1` | – | yes | – | – |
| `C6` | – | yes | – | – |
| `O6` | – | yes | – | – |
| `N2` | – | yes | – | – |
| `N6` | – | yes | – | – |
| `O2` | – | yes | – | – |
| `N4` | – | yes | – | – |
| `O4` | – | yes | – | – |
| `C5M` | – | yes | – | – |
| `C7` | – | – | – | – |
| `H1` | – | – | – | – |
| `p` | – | – | – | – |
| `C4*` | – | – | – | – |

| atom name | in `backbone` when resname is `GLY` | when resname is `LIG` |
|---|---|---|
| `N` | yes | – |
| `CA` | yes | – |
| `C` | yes | – |
| `O` | yes | – |
| `OT1` | – | – |
| `OXT` | – | – |
| `CB` | – | – |
| `H` | – | – |
| `HA` | – | – |
| `n` | – | – |
| `ca` | – | – |

### 6.7 `byres`, `same … as`, `bonded`

* `byres X`: the residues (by residue identity, i.e. resindex — not by resid
  number) containing at least one atom of X's result; returns all atoms of
  the **searched group** in those residues.
* `same P as X`: if X's result is empty → empty (no attribute access).
  Otherwise returns the atoms of the searched group whose P equals P of at
  least one atom of X's result:
  * `residue` / `resindex` → same residue; `segment` / `segindex` → same
    segment; `fragment` → same bonded fragment (connected component of the
    whole system's bond graph; an atom without bonds is its own fragment;
    needs a bond attribute, which may be empty).
  * `x`, `y`, `z` → coordinate equal within |a − b| ≤ 1e-8 + 1e-5·|b| (b = the
    reference atom's coordinate; always these tolerances).
  * any other accepted property → exact equality of the attribute value (no
    float tolerance, e.g. `same mass as`, `same charge as`); `bfactor` reads
    the tempfactor values (and emits the alias warning twice).
  * Accepted property names = the "Unknown same property" list in
    [SV §4.3]. Plural or capitalised spellings are rejected; the entry `None`
    in that list cannot be typed (the text `None` is rejected).
* `bonded X`: consider the bonds that have at least one atom in the
  searched group. If the system has no bond attribute → NoDataError; if it
  has one but no such bond exists → warning `Bonded selection has 0 bonds`
  and empty result. Otherwise returns every atom that is bonded (through one
  of those bonds) to an atom of X's result — **from the whole system**, not
  restricted to the searched group. Atoms of X are included only if bonded to
  another atom of X.

SYN (whole system):

| Selection string | MDAnalysis result |
|---|---|
| `byres name CA` | 8 atoms: 1-8 |
| `byres index 11` | 1 atoms: 12 |
| `byres resid 3` | 3 atoms: 9-11 |
| `byres (name OW or name NA)` | 4 atoms: 9-12 |
| `byres not all` | 0 atoms: (none) |
| `byres all` | 12 atoms: 1-12 |
| `same residue as name OW` | 3 atoms: 9-11 |
| `same residue as index 0 11` | 5 atoms: 1-4,12 |
| `same segment as name NA` | 4 atoms: 9-12 |
| `same segment as resid 1` | 8 atoms: 1-8 |
| `same resid as name NA` | 4 atoms: 9-12 |
| `same resnum as name NA` | 1 atoms: 12 |
| `same name as index 0` | 2 atoms: 1,5 |
| `same type as name HW1` | 2 atoms: 10,11 |
| `same resname as index 11` | 1 atoms: 12 |
| `same segid as index 0` | 8 atoms: 1-8 |
| `same mass as index 1` | 4 atoms: 2,3,6,7 |
| `same mass as name O` | 3 atoms: 4,8,9 |
| `same charge as index 0` | 2 atoms: 1,5 |
| `same charge as name OW` | 1 atoms: 9 |
| `same x as index 3` | 1 atoms: 4 |
| `same y as index 0` | 12 atoms: 1-12 |
| `same z as all` | 12 atoms: 1-12 |
| `same fragment as index 0` | 8 atoms: 1-8 |
| `same fragment as index 11` | 1 atoms: 12 |
| `same fragment as index 9 11` | 4 atoms: 9-12 |
| `same bfactor as index 0` | 1 atoms: 1 (+2× bfactor-alias warning) |
| `same radius as all` | NoDataError: `This Universe does not contain radius information` |
| `same resid as none_such` | SelectionError: `Unknown selection token: 'none_such'` |
| `same index as all` | 12 atoms: 1-12 |
| `same element as all` | NoDataError: `This Universe does not contain element information` |
| `same tempfactor as all` | 12 atoms: 1-12 |
| `same chainID as all` | NoDataError: `This Universe does not contain chainID information` |
| `same residue as not all` | 0 atoms: (none) |
| `same fragment as not all` | 0 atoms: (none) |
| `bonded index 0` | 1 atoms: 2 |
| `bonded index 2` | 3 atoms: 2,4,5 |
| `bonded index 11` | 0 atoms: (none) |
| `bonded name C` | 5 atoms: 2,4-6,8 |
| `bonded all` | 11 atoms: 1-11 |
| `bonded not all` | 0 atoms: (none) |
| `bonded bonded index 0` | 2 atoms: 1,3 |
| `bonded resid 1` | 5 atoms: 1-5 |
| `not bonded index 2` | 9 atoms: 1,3,6-12 |

| Selection string | MDAnalysis result |
|---|---|
| `same aromaticity as index 1` | 6 atoms: 2,3,6,7,10,11 |
| `same formalcharge as index 11` | 1 atoms: 12 |
| `same chirality as index 1` | 1 atoms: 2 |
| `same chirality as index 0` | 10 atoms: 1,3-5,7-12 |
| `same bonds as index 0` | IndexError: `too many indices for array: array is 1-dimensional, but 2 were indexed` |
| `same index as index 0 1` | 2 atoms: 1,2 |
| `same id as index 0` | 1 atoms: 1 |

Fixture (additional accepted `same` properties):

| Selection string | MDAnalysis result |
|---|---|
| `same index as name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `same id as name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `same element as name MG` | 3 atoms: 49-51 |
| `same chainID as name N9` | 23 atoms: 1-23 |
| `same icode as resid 3A` | 4 atoms: 18,19,30,31 |
| `same altLoc as all` | 60 atoms: 1-60 |
| `same record_type as name MG` | 12 atoms: 49-60 |
| `same occupancy as all` | 60 atoms: 1-60 |
| `same tempfactor as name P` | 60 atoms: 1-60 |
| `same resindex as name N9` | 8 atoms: 1-8 |
| `same segindex as name N9` | 23 atoms: 1-23 |
| `same segment as name N9` | 23 atoms: 1-23 |
| `same bonds as all` | NoDataError: `This Universe does not contain bonds information` |
| `same topologyattr as all` | NoDataError: `This Universe does not contain topologyattr information` |
| `same atomattr as all` | NoDataError: `This Universe does not contain atomattr information` |
| `same model as all` | NoDataError: `This Universe does not contain model information` |
| `same type_index as all` | NoDataError: `This Universe does not contain type_index information` |
| `same x as name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `same z as name P` | 52 atoms: 1-4,6-11,13-18,20-25,27-32,34-39,41-46,48-53,55-60 |
| `same y as name H1` | 36 atoms: 1,2,4,6,7,9,11,12,14,16,17,19,21,22,24,26,27,29,31,32,34,36,37,39,41,42,44,46,47,49,51,52,54,56,57,59 |

SYN with a bond attribute that holds zero bonds:

| Selection string | MDAnalysis result |
|---|---|
| `bonded all` | 0 atoms: (none) (+warning `Bonded selection has 0 bonds`) |
| `bonded index 0` | 0 atoms: (none) (+warning `Bonded selection has 0 bonds`) |
| `same fragment as index 0` | 1 atoms: 1 |
| `same fragment as all` | 12 atoms: 1-12 |

SYN, searched group = atoms 1–3 (index+1):

| Selection string | MDAnalysis result |
|---|---|
| `bonded name C` | 3 atoms: 2,4,5 |
| `bonded index 2` | 3 atoms: 2,4,5 |
| `bonded global name N` | 2 atoms: 2,3 |
| `bonded global index 4` | 1 atoms: 3 |
| `same fragment as name N` | 3 atoms: 1-3 |
| `same fragment as global name OW` | 0 atoms: (none) |
| `byres index 0` | 3 atoms: 1-3 |
| `same residue as global index 4` | 0 atoms: (none) |
| `same resname as global resid 2` | 0 atoms: (none) |

SYN, empty searched group:

| Selection string | MDAnalysis result |
|---|---|
| `bonded all` | 0 atoms: (none) (+warning `Bonded selection has 0 bonds`) |
| `same fragment as all` | 0 atoms: (none) |
| `charge 0` | 0 atoms: (none) (+float-equality warning) |
| `name N` | 0 atoms: (none) |

### 6.8 `atom`, `group`

* `atom SEG RES NAME`: atoms of the searched group whose name equals NAME
  exactly (no glob), whose resid equals RES (insertion code ignored, so
  `atom A 3 P` matches resid 3 and 3A), and whose segid equals SEG exactly.
  The three tests are applied in the order name, resid, segid, and a test is
  only performed (so its attribute only required) if the previous ones left
  at least one atom. No backslash removal on any of the three tokens.
* `group NAME`: atoms of the searched group that are also in the named
  group passed with the call (comparison by system index); `global group
  NAME` returns the named group's atoms unrestricted. The name `updating` is
  reserved: TypeError
  `The 'updating' keyword is reserved and cannot be used as a selection group name.`
  (raised at parse time, before the lookup). Unknown name →
  `Selection failed: 'Failed to find group: <name>'`.

| Selection string | MDAnalysis result |
|---|---|
| `atom A 1 P` | 1 atoms: 1 |
| `atom A 1 C4'` | 1 atoms: 3 |
| `atom A 3 P` | 2 atoms: 14,18 |
| `atom A 3A P` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: '3A''` |
| `atom B 11 P` | 2 atoms: 27,30 |
| `atom A 1 C*` | 0 atoms: (none) |
| `atom a 1 P` | 0 atoms: (none) |
| `atom A 01 P` | 1 atoms: 1 |
| `atom A +1 P` | 1 atoms: 1 |
| `atom A 1_0 P` | 0 atoms: (none) |
| `atom D 201 MG` | 1 atoms: 49 |
| `atom W 301 O or atom W 302 O` | 2 atoms: 55,58 |
| `atom A 1 \P` | 0 atoms: (none) |
| `atom \A 1 P` | 0 atoms: (none) |
| `atom A 1 and` | 0 atoms: (none) |
| `not atom A 1 P` | 59 atoms: 2-60 |

| Selection string | 3 atoms with names A,B,C only (no resid, no segid) | same plus resid 1 for the single residue (no segid) |
|---|---|---|
| `atom X 1 Q` | 0 atoms: (none) | 0 atoms: (none) |
| `atom X 1 A` | NoDataError: `This Universe does not contain resid information` | NoDataError: `This Universe does not contain segid information` |
| `atom X 2 A` | NoDataError: `This Universe does not contain resid information` | 0 atoms: (none) |
| `atom X x A` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: 'x''` | SelectionError: `Selection failed: 'invalid literal for int() with base 10: 'x''` |

| call | MDAnalysis result |
|---|---|
| whole fixture, `group grp`, grp = `name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| whole fixture, `group grp and resid 1`, grp = `name P` | 1 atoms: 1 |
| whole fixture, `group a or group b`, a = `name P`, b = `resid 1` | 16 atoms: 1-9,14,18,20,24,27,30,32 |
| whole fixture, `group c`, a = `name P` | SelectionError: `Selection failed: 'Failed to find group: c'` |
| whole fixture, `group \a`, a = `name P` (no unescaping) | SelectionError: `Selection failed: 'Failed to find group: \a'` |
| group `resid 1 2`, `group grp`, grp = `name P` | 2 atoms: 1,9 |
| group `resid 1 2`, `global group grp`, grp = `name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| whole fixture, `group updating` | TypeError: `The 'updating' keyword is reserved and cannot be used as a selection group name.` |
| whole fixture, `name P`, grp = the string 'x' | TypeError: `Passed groups must be AtomGroups. You provided str for group 'grp'` |
| whole fixture, `name P`, grp = a Residue list | TypeError: `Passed groups must be AtomGroups. You provided list for group 'grp'` |

### 6.9 `aromaticity`, `aromatic`, `chirality`

| Selection string | MDAnalysis result |
|---|---|
| `aromatic` | 6 atoms: 2,3,6,7,10,11 |
| `aromaticity` | 6 atoms: 2,3,6,7,10,11 |
| `aromaticity true` | 6 atoms: 2,3,6,7,10,11 |
| `aromaticity True` | 6 atoms: 2,3,6,7,10,11 |
| `aromaticity FALSE` | 6 atoms: 1,4,5,8,9,12 |
| `aromaticity true false` | 12 atoms: 1-12 |
| `aromaticity yes` | SelectionError: `Selection failed: ''yes' is an invalid value for boolean selection. Use 'True' or 'False''` |
| `aromaticity \true` | 6 atoms: 2,3,6,7,10,11 |
| `aromaticity and name C` | 2 atoms: 3,7 |
| `aromaticity false and name C` | 0 atoms: (none) |
| `not aromatic` | 6 atoms: 1,4,5,8,9,12 |
| `chirality S` | 1 atoms: 2 |
| `chirality s` | 0 atoms: (none) |
| `chirality R S` | 2 atoms: 2,6 |
| `chirality *` | 0 atoms: (none) |
| `chirality \` | 10 atoms: 1,3-5,7-12 |
| `chirality SR` | 0 atoms: (none) |
| `not chirality \` | 2 atoms: 2,6 |

### 6.10 Geometric keywords

Box handling, periodic images and the search-method choice are in
`selection-grid-method.md`. Without a box (or with *periodic* off),
distances are Euclidean:

* `around R X`: atoms of the searched group **not** in X's result whose
  distance to at least one atom of X's result is ≤ R.
* `sphzone R X`: atoms of the searched group within distance ≤ R of the
  centre of geometry (unweighted mean position) of X's result (X's own atoms
  can be selected).
* `sphlayer R1 R2 X`: same centre; R1 < distance ≤ R2.
* `isolayer R1 R2 X`: atoms of the searched group not in X whose distance
  to the nearest atom of X is > R1 and ≤ R2 (i.e. within R2 of some X atom and
  not within R1 of any).
* `point X Y Z R`: atoms of the searched group within ≤ R of (X, Y, Z)
  (coordinates rounded to single precision).
* `cyzone R ZMAX ZMIN X`: with d = atom position − centre of geometry of X:
  ZMIN < d_z < ZMAX (strict) and d_x² + d_y² < R² (strict).
  `cylayer R1 R2 ZMAX ZMIN X` additionally requires d_x² + d_y² > R1²
  (strict).
* Empty X result → empty. Empty searched group → empty.

SYN (atoms on the x axis 1.5 Å apart, no box):

| Selection string | MDAnalysis result |
|---|---|
| `around 1.5 index 0` | 1 atoms: 2 |
| `around 1.4999 index 0` | 0 atoms: (none) |
| `around 3 index 0` | 2 atoms: 2,3 |
| `around 1.5 index 0 1` | 1 atoms: 3 |
| `around 0 index 0` | 0 atoms: (none) |
| `around -1 index 0` | 0 atoms: (none) |
| `around 100 index 0` | 11 atoms: 2-12 |
| `around 1.5 not all` | 0 atoms: (none) |
| `around 1.5 all` | 0 atoms: (none) |
| `sphzone 1.5 index 0` | 2 atoms: 1,2 |
| `sphzone 0 index 0` | 1 atoms: 1 |
| `sphzone 0.75 index 0 1` | 2 atoms: 1,2 |
| `sphzone 0.7499 index 0 1` | 0 atoms: (none) |
| `sphzone 3 index 0 2` | 4 atoms: 1-4 |
| `sphzone 1 not all` | 0 atoms: (none) |
| `sphlayer 1.5 3 index 0` | 1 atoms: 3 |
| `sphlayer 0 1.5 index 0` | 1 atoms: 2 |
| `sphlayer 1.4999 3 index 0` | 2 atoms: 2,3 |
| `sphlayer 3 1.5 index 0` | 0 atoms: (none) |
| `sphlayer -1 0 index 0` | 1 atoms: 1 |
| `isolayer 1.5 3 index 0` | 1 atoms: 3 |
| `isolayer 0 1.5 index 0` | 1 atoms: 2 |
| `isolayer 1.4999 3 index 5` | 4 atoms: 4,5,7,8 |
| `isolayer 1.5 3 index 5 6` | 2 atoms: 4,9 |
| `point 0 0 0 1.5` | 2 atoms: 1,2 |
| `point 0 0 0 1.4999` | 1 atoms: 1 |
| `point 3 0 0 0` | 1 atoms: 3 |
| `point 3 0 0 -1` | 0 atoms: (none) |
| `point 1 1 1 1.8` | 2 atoms: 1,2 |
| `point 0.75 0 0 0.75` | 2 atoms: 1,2 |
| `cyzone 1.5 1 -1 index 0` | 1 atoms: 1 |
| `cyzone 1.5001 1 -1 index 0` | 2 atoms: 1,2 |
| `cyzone 3.1 1 -1 index 2` | 5 atoms: 1-5 |
| `cyzone 100 0 -1 index 0` | 0 atoms: (none) |
| `cyzone 100 1 0 index 0` | 0 atoms: (none) |
| `cyzone 100 0.1 -0.1 index 0` | 12 atoms: 1-12 |
| `cyzone 2 -1 1 index 0` | 0 atoms: (none) |
| `cylayer 1.5 3.1 1 -1 index 2` | 2 atoms: 1,5 |
| `cylayer 1.4999 3.1 1 -1 index 2` | 4 atoms: 1,2,4,5 |
| `cylayer 0 1 1 -1 index 0` | 0 atoms: (none) |

Fixture (no box):

| Selection string | MDAnalysis result |
|---|---|
| `around 1.3 name N9` | 1 atoms: 6 |
| `around 1.3 name N9 or name H1` | 4 atoms: 6,8,57,58 |
| `sphzone 2 name N9` | 3 atoms: 4-6 |
| `sphlayer 1 2 name N9` | 2 atoms: 4,6 |
| `isolayer 1 2 name N9` | 2 atoms: 4,6 |
| `point 6 0 1 1.5` | 2 atoms: 5,6 |
| `cyzone 2 1 -1 name N9` | 3 atoms: 4-6 |
| `cylayer 1 2 1 -1 name N9` | 2 atoms: 4,6 |

### 6.11 Result order with *sorted* off

With *sorted* off each string's result is de-duplicated keeping the first
occurrence, without sorting: keywords that filter the searched group keep
the searched group's order; `or` yields ascending index order; `global`
yields the whole-system order. Searched group below = the atoms with
index+1 = 11, 4, 6, 4, 1, 59 in that order (atom 4 twice).

| Selection string | MDAnalysis result |
|---|---|
| `all` | 5 atoms: 11,4,6,1,59 |
| `name N*` | 2 atoms: 11,6 |
| `not name N1` | 3 atoms: 4,1,59 |
| `index 3 or index 10` | 2 atoms: 4,11 |
| `index 10 0` | 2 atoms: 11,1 |
| `index 10 and index 0 10` | 1 atoms: 11 |
| `index 0 10 and index 10` | 1 atoms: 11 |
| `global name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |
| `byres name N1` | 4 atoms: 11,4,6,1 |
| `resid 1 or resid 2` | 4 atoms: 1,4,6,11 |
| `around 2 index 0` | 0 atoms: (none) |
| `same resid as index 10` | 1 atoms: 11 |

The same group with *sorted* on (default):

| Selection string | MDAnalysis result |
|---|---|
| `all` | 5 atoms: 1,4,6,11,59 |
| `name N*` | 2 atoms: 6,11 |
| `global name P` | 9 atoms: 1,9,14,18,20,24,27,30,32 |

Several strings, sorted on (each string sorted, then concatenated):

| Selection string | MDAnalysis result |
|---|---|
| `index 5` , `index 1` | 2 atoms: 6,2 |
| `index 5 1` , `index 1` | 3 atoms: 2,6,2 |
| `index 5 1` , `index 1` | 3 atoms: 2,6,2 |
