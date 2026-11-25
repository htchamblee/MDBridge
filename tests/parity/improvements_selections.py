"""MDBridge's improvements to the selection behaviour of MDATK v0.0.4.1 / MDAnalysis 2.11.0
(docs/IMPROVEMENTS.md section 1, docs/DIFFERENCES.md "Selections"), applied on the test side.

The generators record what MDAnalysis and MDATK do; this module turns that into what MDBridge does
where MDBridge deliberately differs. It is pure Python (plus numpy and MDAnalysis' public
select_atoms) and contains no MDAnalysis code. Items:

SEL-01  prop with reversed operands mirrors the operator (`prop 3 > x` means x < 3).
SEL-02  prop rejects malformed operators: "invalid operator '><' in prop (use <, >, <=, >=, ==, !=)".
SEL-03  keywords whose arguments run out name what they expect ("around: expected a distance and a
        selection, e.g. 'around 3.5 resname MG'", "atom: expected segid resid name", ...).
SEL-04  the end of the selection is called that, not 'None' ("Unexpected end of selection after
        'and'", "Unexpected end of selection (expected ')')"); a blank selection is "empty selection".
SEL-05  atom refuses keywords as its arguments (a backslash escape such as \\and still works).
SEL-06  no grid-search errors; a negative distance or radius is an error ("around: distance must not
        be negative").
SEL-08  readable 'same' property errors; upper-case keywords get a hint.
SEL-09  float range bounds accept the single-value syntax (sign, E, inf); a bad bound gives
        "Failed to parse number: <piece>".
SEL-10  a reversed resid range selects nothing, with or without insertion codes.
SEL-13  atom_elements falls back to the name guess per atom, not per group.
SEL-14  group: tokens with characters a group name cannot hold, and groups with invalid names, give
        "group name 'x y' may only contain letters, digits, '_', '.' and '-'".
SEL-15  the first group of a name wins in group: expansion (as everywhere else).

How the atom selections are modelled: `Model` parses a selection string with MDBridge's grammar
(src/core/atom_selection.cpp) to find its first parse error, if any. Errors whose text MDBridge
keeps are compared with MDAnalysis' text (a mismatch stops the generator: the model is wrong).
Without a parse error, the clauses whose meaning changed are rewritten into strings MDAnalysis
evaluates the same way MDBridge does (prop in canonical order, explicit index lists), and
MDAnalysis evaluates the result.
"""

import re

import numpy as np

# ---- keywords (MDBridge's table, src/core/atom_selection.cpp keywords()) ------------------------

KEYWORDS = {
    "not": "not", "global": "global", "all": "all",
    "protein": "macro", "nucleic": "macro", "water": "macro", "backbone": "macro",
    "nucleicbackbone": "macro", "nucleicbase": "macro", "nucleicsugar": "macro",
    "aromatic": "aromatic", "byres": "byres", "bonded": "bonded", "same": "same",
    "around": "around", "sphzone": "sphzone", "sphlayer": "sphlayer", "isolayer": "isolayer",
    "cyzone": "cyzone", "cylayer": "cylayer", "point": "point", "atom": "atom", "group": "group",
    "prop": "prop", "smarts": "smarts",
    **{k: "text" for k in ("name", "type", "element", "record_type", "chainID", "chainid", "altLoc",
                           "altloc", "resname", "icode", "moltype", "segid")},
    **{k: "integer" for k in ("index", "resindex", "segindex", "id", "bynum", "formalcharge",
                              "nbindex", "resnum", "molnum", "model")},
    **{k: "real" for k in ("radius", "tempfactor", "bfactor", "mass", "charge", "occupancy",
                           "gbscreen", "solventradius", "rmin", "epsilon", "rmin14", "epsilon14")},
    "resid": "resid", "aromaticity": "aromaticity", "chirality": "chirality",
}

# same <property> as: the accepted properties, in the order the message lists them
SAME_PROPERTIES = [
    "fragment", "x", "y", "z", "residue", "segment", "name", "type", "resname", "resid", "segid",
    "mass", "charge", "radius", "bfactor", "resnum", "topologyattr", "index", "resindex",
    "segindex", "atomattr", "id", "element", "record_type", "chainID", "tempfactor",
    "formalcharge", "occupancy", "altLoc", "gbscreen", "solventradius", "nbindex", "rmin",
    "epsilon", "rmin14", "epsilon14", "aromaticity", "chirality", "residueattr", "icode",
    "moltype", "molnum", "segmentattr", "bonds", "ureybradleys", "angles", "dihedrals",
    "impropers", "cmaps", "atomnum", "atomiccharge", "chargegroup", "model", "type_index"]

PROP_NAMES = {
    "x", "y", "z", "topologyattr", "index", "resindex", "segindex", "atomattr", "id", "name",
    "type", "element", "radius", "record_type", "chainID", "tempfactor", "mass", "charge",
    "formalcharge", "occupancy", "altLoc", "gbscreen", "solventradius", "nbindex", "rmin",
    "epsilon", "rmin14", "epsilon14", "aromaticity", "chirality", "residueattr", "resid",
    "resname", "resnum", "icode", "moltype", "molnum", "segmentattr", "segid", "bonds",
    "ureybradleys", "angles", "dihedrals", "impropers", "cmaps", "atomnum", "atomiccharge",
    "chargegroup", "model", "type_index"}

OPERATORS = ("<", ">", "<=", ">=", "==", "!=")
MIRRORED = {"<": ">", ">": "<", "<=": ">=", ">=": "<=", "==": "==", "!=": "!="}
OPERATOR_CHARACTERS = "<>=!"

SIGNATURES = {  # SEL-03
    "around": "around: expected a distance and a selection, e.g. 'around 3.5 resname MG'",
    "isolayer": "isolayer: expected an inner and an outer distance and a selection, "
                "e.g. 'isolayer 2 3.5 resname MG'",
    "sphzone": "sphzone: expected a radius and a selection, e.g. 'sphzone 6 resname MG'",
    "sphlayer": "sphlayer: expected an inner and an outer radius and a selection, "
                "e.g. 'sphlayer 3 6 resname MG'",
    "point": "point: expected x y z radius",
    "cyzone": "cyzone: expected externalRadius zMax zMin and a selection",
    "cylayer": "cylayer: expected innerRadius externalRadius zMax zMin and a selection",
    "atom": "atom: expected segid resid name",
    "same": "same: expected 'same <property> as <selection>'",
    "prop": "prop: expected 'prop [abs] <property> <op> <value>'",
}
NEGATIVE = {"around": "distance", "isolayer": "distance", "sphzone": "radius", "sphlayer": "radius",
            "point": "radius"}  # SEL-06

# float keywords MDBridge evaluates (the others have no data in its topologies)
FLOAT_ATTRIBUTES = {"mass": "masses", "charge": "charges", "occupancy": "occupancies",
                    "tempfactor": "tempfactors", "bfactor": "tempfactors"}

INT_RANGE = re.compile(r"(-?\d+)\s*(?:[:-]| to )\s*(-?\d+)")
RESID_SINGLE = re.compile(r"(-?\d+)(\w?)\Z")
RESID_RANGE = re.compile(r"(-?\d+)(\w?)\s*(?:[:-]| to )\s*(-?\d+)(\w?)")


class SelError(Exception):
    """A selection error with MDBridge's text; `new`: the text differs from MDAnalysis'."""

    def __init__(self, text, new=False):
        super().__init__(text)
        self.text, self.new = text, new


class _ValueErr(Exception):
    """A ValueError inside a keyword's arguments (reported as "Selection failed: '...'")."""

    def __init__(self, text, new=False):
        super().__init__(text)
        self.text, self.new = text, new


def tokenize(text):
    words, word = [], ""
    for c in text:
        if c.isspace():
            if word:
                words.append(word)
            word = ""
        elif c in "()":
            if word:
                words.append(word)
            word = ""
            words.append(c)
        else:
            word += c
    if word:
        words.append(word)
    return words


def py_float(text):
    try:
        return float(text)
    except ValueError as exc:
        raise _ValueErr(str(exc))


def py_int(text):
    try:
        return int(text)
    except ValueError as exc:
        raise _ValueErr(str(exc))


def float_range(piece):
    """SEL-09: (lower, upper) of a range piece, None when it is not one."""
    for k, c in enumerate(piece):
        width = 1 if c in ":-" else 4 if piece.startswith(" to ", k) else 0
        if not width:
            continue
        try:
            return float(piece[:k]), float(piece[k + width:])
        except ValueError:
            continue
    return None


class Model:
    def __init__(self, sel, universe):
        self.words = tokenize(sel)
        self.tokens = self.words + [None]
        self.pos = 0
        self.u = universe
        self.rewrites = []  # (start, end, replacement tokens)

    # -- tokens
    def front(self):
        if self.pos >= len(self.tokens):
            raise SelError("deque index out of range")
        return self.tokens[self.pos]

    def pop(self):
        if self.pos >= len(self.tokens):
            raise SelError("pop from an empty deque")
        self.pos += 1
        return self.tokens[self.pos - 1]

    def end_error(self):
        if not self.words:
            return SelError("empty selection", True)
        return SelError(f"Unexpected end of selection after '{self.words[-1]}'", True)

    @staticmethod
    def is_keyword(token):
        return token is None or token in ("and", "or", "(", ")") or token in KEYWORDS

    # -- grammar
    def parse(self):
        self.expression()
        rest = self.front()
        if rest is not None:
            raise SelError(f"Unexpected token at end of selection string: '{rest}'")

    def rewritten(self):
        if not self.rewrites:
            return None
        out, pos = [], 0
        for start, end, tokens in sorted(self.rewrites):
            out += self.words[pos:start] + tokens
            pos = end
        out += self.words[pos:]
        return " ".join(out)

    def expression(self):
        self.primary()
        while self.front() in ("and", "or"):
            self.pop()
            self.primary()

    def primary(self):
        token = self.pop()
        if token is None:
            raise self.end_error()
        if token == "(":
            self.expression()
            close = self.front()
            if close is None:
                raise SelError("Unexpected end of selection (expected ')')", True)
            if close != ")":
                raise SelError(f"Unexpected token: '{close}' Expected: ')'")
            self.pop()
            return
        kind = KEYWORDS.get(token)
        if kind is None:
            message = f"Unknown selection token: '{token}'"
            lower = token.lower()
            if lower != token and (lower in KEYWORDS or lower in ("and", "or")):
                raise SelError(f"{message} (keywords are lower case: did you mean '{lower}'?)", True)
            raise SelError(message)
        try:
            self.keyword(token, kind, self.pos - 1)
        except _ValueErr as exc:
            raise SelError(f"Selection failed: '{exc.text}'", exc.new)

    def number(self, keyword, check_negative=False):
        token = self.pop()
        if token is None:
            raise SelError(SIGNATURES[keyword], True)
        value = py_float(token)
        if check_negative and value < 0:
            raise SelError(f"{keyword}: {NEGATIVE[keyword]} must not be negative", True)
        return value

    def selection_of(self, keyword):
        if self.front() is None:
            raise SelError(SIGNATURES[keyword], True)
        self.expression()

    def keyword(self, word, kind, start):
        if kind in ("not", "global"):
            self.primary()
        elif kind in ("all", "macro", "aromatic"):
            pass
        elif kind in ("byres", "bonded"):
            self.expression()
        elif kind == "same":
            prop = self.pop()
            if prop is None:
                raise SelError(SIGNATURES["same"], True)
            if prop not in SAME_PROPERTIES:
                raise SelError(f"Unknown 'same' property '{prop}' (choose one of: "
                               f"{', '.join(SAME_PROPERTIES)})", True)
            as_ = self.front()
            if as_ is None:
                raise SelError(SIGNATURES["same"], True)
            if as_ != "as":
                raise SelError(f"Unexpected token: '{as_}' Expected: 'as'")
            self.pop()
            self.selection_of("same")
        elif kind in ("around", "isolayer", "sphzone", "sphlayer"):
            if kind in ("isolayer", "sphlayer"):
                self.number(kind, True)
            self.number(kind, True)
            self.selection_of(kind)
        elif kind == "point":
            for _ in range(3):
                self.number(kind)
            self.number(kind, True)
        elif kind in ("cyzone", "cylayer"):
            for _ in range(3 if kind == "cyzone" else 4):
                self.number(kind)
            self.selection_of(kind)
        elif kind == "smarts":
            depth = 0
            while True:
                t = self.front()
                if t == "(":
                    depth += 1
                elif t == ")":
                    if depth == 0:
                        break
                    depth -= 1
                elif self.is_keyword(t):
                    break
                self.pop()
        elif kind == "atom":
            args, escaped = [], False
            for i in range(3):
                t = self.pop()
                if self.is_keyword(t):
                    raise SelError(SIGNATURES["atom"], True)
                stripped = t.startswith("\\")
                if stripped:
                    t, escaped = t[1:], True
                if i == 1:
                    try:
                        py_int(t)
                    except _ValueErr as exc:
                        exc.new = stripped
                        raise
                args.append(t)
            if escaped:
                assert all(args), args
                self.rewrites.append((start, self.pos, ["atom", *args]))
        elif kind == "group":
            name = self.pop()
            if name is None:
                raise self.end_error()
            if name == "updating":
                raise SelError("The 'updating' keyword is reserved and cannot be used as a "
                               "selection group name.")
            raise _ValueErr(f"Failed to find group: {name}")
        elif kind == "prop":
            self.prop(start)
        elif kind in ("text", "chirality"):
            self.required_values()
        elif kind == "aromaticity":
            for value in self.grab():
                if value.lower() not in ("true", "false"):
                    raise _ValueErr(f"'{value}' is an invalid value for boolean selection. "
                                    "Use 'True' or 'False'")
        elif kind == "integer":
            for piece in self.range_values():
                try:
                    int(piece)
                except ValueError:
                    if not INT_RANGE.match(piece):
                        raise _ValueErr(f"Failed to parse number: {piece}")
        elif kind == "real":
            self.reals(word, start)
        elif kind == "resid":
            self.resids(start)
        else:
            raise AssertionError(kind)

    def grab(self):
        values = []
        while not self.is_keyword(self.front()):
            value = self.pop()
            values.append(value[1:] if value.startswith("\\") else value)
        return values

    def required_values(self):
        values = self.grab()
        if not values:
            if self.front() is None:
                raise self.end_error()
            raise _ValueErr(f"Unexpected token '{self.front()}'")
        return values

    def range_values(self):
        tokens = self.grab()
        if not tokens:
            if self.front() is None:
                raise self.end_error()
            raise _ValueErr(f"Unexpected token: '{self.front()}'")
        pieces = []
        i = 0
        while i < len(tokens):
            token = tokens[i]
            if token in ("to", ":", "-"):
                if not pieces or i + 1 == len(tokens):
                    raise SelError("Invalid expression given: " + " ".join(pieces) + " " + token + " " +
                                   " ".join(tokens[i + 1:]))
                pieces[-1] += " " + token + " " + tokens[i + 1]
                i += 2
                continue
            if pieces and (token.startswith(("--", "to", ":")) or pieces[-1].endswith(("to", ":", "-"))):
                pieces[-1] += " " + token
            else:
                pieces.append(token)
            i += 1
        return pieces

    def reals(self, keyword, start):
        singles, ranges = [], []
        for piece in self.range_values():
            try:
                singles.append(float(piece))
                continue
            except ValueError:
                pass
            bounds = float_range(piece)
            if bounds is None:
                raise _ValueErr(f"Failed to parse number: {piece}", True)
            ranges.append(bounds)
        if not ranges:
            return
        # SEL-09: MDAnalysis' range syntax is narrower; the clause becomes the atoms it selects
        attr = FLOAT_ATTRIBUTES.get(keyword)
        if attr is None or not hasattr(self.u.atoms, attr):
            replacement = [keyword, "0:1"]  # the attribute is missing: the same error
        else:
            v = np.asarray(getattr(self.u.atoms, attr), dtype=np.float64)
            hit = np.zeros(len(v), dtype=bool)
            for x in singles:
                hit |= np.isclose(v, x)
            for lo, hi in ranges:
                hit |= (v >= lo) & (v <= hi)
            replacement = index_list(np.nonzero(hit)[0])
        self.rewrites.append((start, self.pos, replacement))

    def resids(self, start):
        pieces = self.range_values()
        keep, reversed_ = [], False
        for piece in pieces:
            if RESID_SINGLE.match(piece):
                keep.append(piece)
                continue
            m = RESID_RANGE.match(piece)
            if not m:
                raise _ValueErr(f"Failed to parse value: {piece}")
            if int(m.group(1)) > int(m.group(3)):
                reversed_ = True  # SEL-10: selects nothing
            else:
                keep.append(piece)
        if reversed_ and hasattr(self.u.atoms, "icodes"):
            replacement = ["(", "resid", *keep, ")"] if keep else ["(", "not", "all", ")"]
            self.rewrites.append((start, self.pos, replacement))

    def prop(self, start):
        absolute = False
        token = self.pop()
        if token == "abs":
            absolute = True
            token = self.pop()
        if token is None:
            raise SelError(SIGNATURES["prop"], True)
        at = next((k for k, c in enumerate(token) if c in OPERATOR_CHARACTERS), None)
        if at is not None:
            prop, rest = token[:at], token[at:]
            if not prop:
                raise SelError(SIGNATURES["prop"], True)
        else:
            prop, rest = token, self.pop()
            if rest is None:
                raise SelError(SIGNATURES["prop"], True)
        run = len(rest) - len(rest.lstrip(OPERATOR_CHARACTERS))
        op = rest[:run] if run else rest
        if op not in OPERATORS:
            raise SelError(f"invalid operator '{op}' in prop (use <, >, <=, >=, ==, !=)", True)
        value = rest[run:]
        if not value:
            value = self.pop()
            if value is None:
                raise SelError(SIGNATURES["prop"], True)
        if value in PROP_NAMES:  # reversed operands: mirror the operator (SEL-01)
            prop, value, op = value, prop, MIRRORED[op]
        py_float(value)
        self.rewrites.append((start, self.pos, ["prop", *(["abs"] if absolute else []), prop, op, value]))


def index_list(indices):
    if len(indices) == 0:
        return ["(", "not", "all", ")"]
    return ["(", "index", *[str(int(i)) for i in indices], ")"]


def model(universe, sel):
    """('error', text, new) or ('ok', rewritten string or None)."""
    m = Model(sel, universe)
    try:
        m.parse()
    except SelError as exc:
        return "error", exc.text, exc.new
    return "ok", m.rewritten()


GRID_ERRORS = ("too large for box", "Cutoff must be positive", "box dimensions cannot be 0")
INTERNAL = ("'None'", "NoneType", "deque", "dict_keys", "unpack")


def select_atoms(universe, sel, periodic=True):
    """universe.select_atoms(sel, periodic=periodic) as MDBridge evaluates it; raises SelError on
    errors."""
    if sel == "":
        return universe.atoms[[]]
    verdict = model(universe, sel)
    if verdict[0] == "error":
        raise SelError(verdict[1], verdict[2])
    return universe.select_atoms(verdict[1] if verdict[1] is not None else sel, periodic=periodic)


def improve_case(universe, sel, case):
    """What MDBridge gives for `sel`, from MDAnalysis' answer `case` ({"atoms": runs} or {"error":
    text}, as make_atom_selection_cases records it): {"atoms": runs} or {"error": text}. Checks the
    model against MDAnalysis wherever MDBridge keeps MDAnalysis' answer."""
    if sel == "":
        return case
    verdict = model(universe, sel)
    if verdict[0] == "error":
        _, text, new = verdict
        if not new and case.get("error") != text:
            raise AssertionError(f"selection model disagrees with MDAnalysis on {sel!r}: "
                                 f"{text!r} vs {case.get('error', case.get('atoms'))!r}")
        return {"error": text}
    rewritten = verdict[1]
    if rewritten is None:
        error = case.get("error", "")
        if any(s in error for s in INTERNAL + GRID_ERRORS):
            raise AssertionError(f"selection model accepts {sel!r}, MDAnalysis: {error!r}")
        return case
    try:
        return {"atoms": _runs(universe.select_atoms(rewritten).indices)}
    except Exception as exc:  # noqa: BLE001 - what the generators record
        if any(s in str(exc) for s in GRID_ERRORS):
            raise AssertionError(f"grid-search error for {sel!r}: {exc}")
        return {"error": str(exc)}


def _runs(indices):
    out = []
    for i in indices:
        i = int(i)
        if out and out[-1][1] == i - 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return out


# ---- MDATK's selection layer (mdatk/selections.py) --------------------------------------------

PLAIN_NAME = re.compile(r"[A-Za-z0-9_.\-]+")
GROUP_WORD = re.compile(r"[^\s()]*")


def invalid_group_name(name):
    return f"group name {name!r} may only contain letters, digits, '_', '.' and '-'"


def patch_mdatk(S):
    """Replaces the functions of mdatk.selections (module S) whose behaviour MDBridge improves:
    expand_groups (SEL-14, SEL-15), select (the atom selections above) and atom_elements (SEL-13).
    Call before importing MDATK modules that import these names."""

    def expand_groups(text, groups, universe=None, _depth=0):
        if _depth > 8:
            raise S.SelectionError("group definitions are nested too deeply (circular reference?)")
        names = list(dict.fromkeys(g.name for g in groups))
        out, pos = [], 0
        while True:
            at = text.find("group:", pos)
            if at < 0:
                out.append(text[pos:])
                return "".join(out)
            rest = text[at + 6:]
            bad = [g.name for g in groups
                   if g.name and not PLAIN_NAME.fullmatch(g.name) and rest.startswith(g.name)
                   and (len(rest) == len(g.name) or rest[len(g.name)].isspace() or rest[len(g.name)] in "()")]
            if bad:
                raise S.SelectionError(invalid_group_name(max(bad, key=len)))
            name = GROUP_WORD.match(rest).group(0)
            if not name:
                out.append(text[pos:at + 6])
                pos = at + 6
                continue
            if not PLAIN_NAME.fullmatch(name):
                raise S.SelectionError(invalid_group_name(name))
            group = next((g for g in groups if g.name == name), None)  # the first one wins
            if group is None:
                raise S.SelectionError(
                    f"unknown residue group {name!r} (defined: {', '.join(names) or 'none'})")
            out.append(text[pos:at])
            out.append(expand_groups(S.group_selection(group, universe), groups, universe, _depth + 1))
            pos = at + 6 + len(name)

    def select(universe, text, groups=(), allow_empty=False):
        resolved = S.resolve_selection(text, groups, universe)
        try:
            ag = select_atoms(universe, resolved)
        except Exception as exc:  # noqa: BLE001
            raise S.SelectionError(f"invalid selection {text!r} -> {resolved!r}: {exc}") from exc
        if not allow_empty and ag.n_atoms == 0:
            raise S.SelectionError(f"selection {text!r} matched no atoms")
        return ag

    def atom_elements(atoms):
        elements = None
        if hasattr(atoms, "elements"):
            elements = np.asarray(atoms.elements, dtype=object)
        out = []
        for i, a in enumerate(atoms):
            e = elements[i] if elements is not None else None
            if isinstance(e, str) and e.strip():
                out.append(e.strip().upper())
            else:
                out.append(S.guess_element(a.name, a.resname, a.residue.atoms.n_atoms == 1))
        return np.array(out, dtype=object)

    if getattr(S, "_mdbridge_selections", False):
        return
    replaced = {name: getattr(S, name) for name in ("expand_groups", "select", "atom_elements")}
    S.expand_groups = expand_groups
    S.select = select
    S.atom_elements = atom_elements
    S._mdbridge_selections = True
    # MDATK modules imported before this call hold the old functions (`from ..selections import
    # atom_elements`): give them the new ones too, so the order of the patches does not matter
    import sys
    for module_name, module in list(sys.modules.items()):
        if module_name.startswith("mdatk.") and module is not S and module is not None:
            for name, old in replaced.items():
                if getattr(module, name, None) is old:
                    setattr(module, name, getattr(S, name))


def patch_mdatk_analysis_base(base):
    """SEL-15 in the cache keys: referenced_groups of mdatk.analysis.base (module `base`) takes the
    first group of a repeated name, the one the selections use (MDATK: the last)."""
    import copy

    if getattr(base, "_mdbridge_first_groups", False):
        return
    base._mdbridge_first_groups = True
    original = base.referenced_groups

    def referenced_groups(texts, project):
        first, seen = [], set()
        for g in project.groups:
            if g.name not in seen:
                seen.add(g.name)
                first.append(g)
        proxy = copy.copy(project)
        proxy.groups = first
        return original(texts, proxy)

    base.referenced_groups = referenced_groups
