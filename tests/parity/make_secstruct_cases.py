"""What MDATK's secondary-structure code (mdatk/secstruct.py) gives, for src/chem/secstruct.cpp
(tests/chem/test_secstruct.cpp).

    /opt/mdbridge-venv/bin/python -I -B make_secstruct_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy and MDAnalysis (MDATK's tests/conftest.py imports it lazily; the independent tables are
applied first, as in every generator that runs MDATK in-process). The geometry helpers of MDATK's
tests/test_basepairs.py (nucleotide, partner, detector_input) build the base-pair geometries; they
import pytest, which is replaced by a minimal stand-in when it is not installed (only its decorators
run at import time). Nothing is written into the checkout; the output is reproducible (fixed seeds).
Writes OUT_DIR/secstruct_cases.json:

  constants       OPENERS, CLOSERS, UNPAIRED, CHAIN_SEPARATOR, PAIR_RULES, C1_WINDOW, C1_IDEAL,
                  SCORE_TIE, sorted BASE_PAIR_ATOMS, KIND_LABELS
  residue_names   base_letter and canonical_name of every name of up to three characters over an
                  alphabet of base letters, prefixes, terminals and others (both cases), aliases,
                  padded names, non-nucleotides
  lw_kind         lw_kind of Leontis-Westhof codes
  parse           dotbracket_to_pairs of MDATK's test texts, every bracket level, separators
                  (Unicode whitespace included), random balanced and random unbalanced texts:
                  pairs, n, chain breaks, or the ValueError message
  write           pairs_to_dotbracket of pair lists (MDATK's tests, random nested and crossing
                  sets up to 30 levels and beyond, reversed / duplicated / invalid pairs, chain
                  breaks in and out of range): the text or the message
  stems           _stems of pair lists
  crosses         _crosses of pair couples
  unique_pairs    unique_pairs of pair lists
  detector        PairDetector(resnames, atom_index, cutoff, wobble) on frames: the geometries of
                  MDATK's tests/test_basepairs.py, perturbed copies, dense random clouds, missing
                  atoms, C1* names, NaN / infinite coordinates, distances at the cutoff and C1'
                  windows exactly, score near-ties, L-RNA names and every residue-name scheme.
                  Per frame: candidates(xyz) as rows [i, j, kind, score, dev] and detect(xyz)
  select          detect() on crafted candidate arrays (what MDATK's tie test does by replacing
                  candidates): scores on and around the 0.05 grid, equal deviations
  lattice         a large detector case rebuilt by the C++ test from a formula (coordinates are
                  multiples of 0.25 A, see lattice_xyz): candidate count, SHA-256 of the candidate
                  array (float64, little-endian, row-major) and of detect's int64 array

Coordinates the C++ test receives are float.hex() strings (flattened x y z per row): the geometries
built by the test helpers go through numpy's LAPACK SVD (CPU-dependent in the last bits), so they
are rounded to float32 first, as MDBridge's frames are, which makes the file reproducible; the
perturbations added afterwards are float64. Floats in results are float.hex() strings; "nan",
"inf" and "-inf" for non-finite values.
"""

import hashlib
import importlib.util
import itertools
import json
import math
import os
import random
import sys
import types

if len(sys.argv) != 3:
    sys.exit("usage: make_secstruct_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
sys.path.insert(0, CHECKOUT)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)

import numpy as np  # noqa: E402

from mdatk import secstruct as S  # noqa: E402

import improvements_nucleotides  # noqa: E402

improvements_nucleotides.apply_mdatk()  # NA-01 (docs/DIFFERENCES.md): L-DNA names, CCD atom names, backbone guard


def _pytest_stand_in():
    """Enough of pytest for importing MDATK's test modules (decorators only)."""
    try:
        import pytest  # noqa: F401
        return
    except ImportError:
        pass
    mod = types.ModuleType("pytest")

    def fixture(*args, **kwargs):
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return lambda f: f

    mod.fixture = fixture
    mod.mark = types.SimpleNamespace(parametrize=lambda *a, **k: (lambda f: f))
    mod.importorskip = lambda name: __import__(name)
    sys.modules["pytest"] = mod


def _load_test(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(CHECKOUT, "tests", name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


_pytest_stand_in()
sys.path.insert(0, os.path.join(CHECKOUT, "tests"))  # test_basepairs imports conftest, test_pdbtools
_load_test("conftest")
_load_test("test_pdbtools")
TB = _load_test("test_basepairs")

OUT = os.path.abspath(sys.argv[2])
rng = np.random.default_rng(20261124)
prng = random.Random(20261124)


def fhex(x):
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x.hex()


def error(fn, *args):
    try:
        return {"result": fn(*args)}
    except ValueError as exc:
        return {"error": str(exc)}


def pairs_json(pairs):
    return [[int(i), int(j)] for i, j in pairs]


# --- constants and names ---------------------------------------------------------------------------

def constants():
    return {
        "openers": S.OPENERS, "closers": S.CLOSERS, "unpaired": S.UNPAIRED,
        "chain_separator": S.CHAIN_SEPARATOR,
        "pair_rules": [{"base": b, "partners": list(p), "kind": k, "atoms": [list(a) for a in atoms],
                        "required": need} for b, p, k, atoms, need in S.PAIR_RULES],
        "c1_window": [fhex(v) for v in S.C1_WINDOW], "c1_ideal": fhex(S.C1_IDEAL),
        "score_tie": fhex(S.SCORE_TIE), "base_pair_atoms": sorted(S.BASE_PAIR_ATOMS),
        "kind_labels": {str(k): v for k, v in S.KIND_LABELS.items()},
        "kinds": {"WATSON_CRICK": S.WATSON_CRICK, "WOBBLE": S.WOBBLE, "OTHER": S.OTHER},
    }


def residue_names():
    alphabet = "ACGTUNRLD053X"
    names = [""]
    for length in (1, 2, 3):
        names += ["".join(t) for t in itertools.product(alphabet, repeat=length)]
    names += [n.lower() for n in names if n.lower() != n][::7]
    names += ["ADE", "CYT", "GUA", "URA", "URI", "THY", "ade", "Thy", " G ", "\tDA5\n", "　LG3　",
              "G​", "RG5 ", "LIG", "ALA", "MG", "LA", "HOH", "DT3", "0U5", "DTN", "RAN", "GUAN", "ADE5",
              "  ", "É", "ǵ", "D", "R", "L", "0"]
    out = []
    for name in names:
        out.append({"name": name, "base": S.base_letter(name), "canonical": S.canonical_name(name)})
    return out


def lw_kinds():
    codes = ["WCc", "GUc", "WHt", "WWt", "SWc", "HSt", "XXX", "", "wcc", "WCc ", "GUt", "WC"]
    return [{"code": c, "kind": S.lw_kind(c)} for c in codes]


# --- dot-bracket ---------------------------------------------------------------------------------------

def parse_case(text):
    try:
        pairs, n, breaks = S.dotbracket_to_pairs(text)
        return {"text": text, "pairs": pairs_json(pairs), "n": n, "breaks": breaks}
    except ValueError as exc:
        return {"text": text, "error": str(exc)}


def random_pairs(n, n_helices, crossing=True, max_len=6):
    """Random helices on n residues; with crossing=False only nested ones."""
    free = set(range(n))
    pairs = []
    for _ in range(n_helices):
        length = prng.randint(1, max_len)
        i = prng.randrange(n)
        j = prng.randrange(n)
        i, j = min(i, j), max(i, j)
        stem = [(i + k, j - k) for k in range(length) if i + k < j - k - 0]
        stem = [(a, b) for a, b in stem if a in free and b in free and a != b]
        if not crossing and any(S._crosses(p, q) for p in stem for q in pairs):
            continue
        for a, b in stem:
            free.discard(a)
            free.discard(b)
            pairs.append((a, b))
    return pairs


def parse_cases():
    texts = ["((..))", "((..[[..))..]]{<Aa>}", "(((&)))", "((( )))", "(((  &  )))", "&(((&)))&", "..)",
             "([.)", "((.", "(#)", "((((....))))", "((..[[..))..]]", "(((&)))..((..))", "(([[{{<<..))]]}}>>",
             ".....", "(&)", "((..[[..))..{{..]]..}}", "AB..((..ab..))", "[{..((..]}..))", "", "&", "&&", " ",
             "&.", ".&", "..&&..", ". & . & .", "(\t)", "(\n)", "(\r\n)", "(\x0b\x0c\x1c\x1d\x1e\x1f\x85)",
             "(     　)", "(​)", "(é)", "(\x00)", "(')", "(\")", "(\\)",
             "(]", "[)", ")(", "a", "A", "Aa", "aA", "(((", ")))", "((.))))", "((((.))", "(.[.{.<.A.)", "(.Z.z.)",
             "-,:_.", "(-,:_)", "(.)　(.)", "\U0001F600", "(\U0001F600)", "(.)&&(.)&", "&&(.)", " &(.)",
             "x", "(x)", "(1)", "<>", "{}", "[]", "()" * 40, "(" * 30 + ")" * 30, "<<<.>>>", "((.((.)).))",
             "(((([[[[))))]]]]", "ABCDEFGHIJKLMNOPQRSTUVWXYZ" + "zyxwvutsrqponmlkjihgfedcba",
             "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz", "([{<" + "ABCDEFGHIJKLMNOPQRSTUVWXYZ" +
             ")]}>" + "abcdefghijklmnopqrstuvwxyz", ".((...))&.[[..]].&..{..}..", "(.(.(.(.).).).)"]
    out = [parse_case(t) for t in texts]
    # random balanced texts from random structures (pseudoknots, chain breaks)
    for _ in range(150):
        n = prng.randint(1, 120)
        pairs = random_pairs(n, prng.randint(0, 12), crossing=prng.random() < 0.7)
        breaks = sorted(prng.sample(range(-2, n + 3), prng.randint(0, 4)))
        try:
            text = S.pairs_to_dotbracket(pairs, n, breaks)
        except ValueError:
            continue
        if prng.random() < 0.3:  # other separators and unpaired symbols
            text = "".join(prng.choice(" \t&") if c == "&" else (prng.choice(S.UNPAIRED) if c == "." else c)
                           for c in text)
        out.append(parse_case(text))
    # random soups: mostly errors
    soup = "().[]{}<>AaBbZz.-,:_& \t#x　"
    for _ in range(150):
        out.append(parse_case("".join(prng.choice(soup) for _ in range(prng.randint(1, 30)))))
    # one character changed in a balanced text
    for _ in range(60):
        n = prng.randint(4, 60)
        text = S.pairs_to_dotbracket(random_pairs(n, 4, crossing=False), n)
        k = prng.randrange(len(text))
        out.append(parse_case(text[:k] + prng.choice("()[]<>Aa.#") + text[k + 1:]))
    return out


def write_case(pairs, n, breaks=None):
    args = (pairs, n) if breaks is None else (pairs, n, breaks)
    case = {"pairs": pairs_json(pairs), "n": n}
    if breaks is not None:
        case["breaks"] = list(breaks)
    try:
        case["text"] = S.pairs_to_dotbracket(*args)
    except ValueError as exc:
        case["error"] = str(exc)
    return case


def write_cases():
    out = []
    for text in ["((((....))))", "((..[[..))..]]", "(((&)))..((..))", "(([[{{<<..))]]}}>>", ".....", "(&)",
                 "((..[[..))..{{..]]..}}", "AB..((..ab..))"]:
        pairs, n, breaks = S.dotbracket_to_pairs(text)
        out.append(write_case(pairs, n, breaks))
    out.append(write_case([(0, 7), (1, 6), (3, 13), (4, 12), (5, 11), (2, 14)], 15))
    out.append(write_case([(5, 0)], 6))
    out.append(write_case([(0, 5), (0, 5)], 6))
    out.append(write_case([(0, 5), (5, 0)], 6))
    out.append(write_case([(0, 5), (0, 4)], 6))
    out.append(write_case([(0, 5), (1, 5)], 6))
    out.append(write_case([(0, 5), (2, 3), (3, 4)], 6))
    out.append(write_case([(0, 6)], 6))
    out.append(write_case([(-1, 3)], 6))
    out.append(write_case([(3, 3)], 6))
    out.append(write_case([(2, 1)], 2))
    out.append(write_case([], 4, [0, 2, 4]))
    out.append(write_case([], 4, [2, 2, 3, -1, 9]))
    out.append(write_case([], 0))
    out.append(write_case([], 0, [0, 1]))
    out.append(write_case([], -3))
    out.append(write_case([(0, 1)], -3))
    out.append(write_case([(0, 1)], 2, [1]))
    out.append(write_case([(0, 2), (1, 3)], 4))
    out.append(write_case([(1, 3), (0, 2)], 4))
    out.append(write_case([(0, 3), (1, 2), (4, 7), (5, 6)], 8))
    out.append(write_case([(0, 9), (1, 8), (3, 6), (2, 7)], 10))  # one stem out of order
    out.append(write_case([(0, 4), (1, 3), (2, 6)], 7))
    # mutually crossing single pairs: 30 levels fit, 31 do not
    for count in (26, 29, 30, 31, 32):
        out.append(write_case([(k, count + k) for k in range(count)], 2 * count))
    # equal-length stems: order by first pair
    out.append(write_case([(0, 4), (2, 6), (1, 5), (3, 7)], 8))
    out.append(write_case([(2, 9), (3, 8), (0, 5), (1, 4)], 10))
    for _ in range(200):
        n = prng.randint(1, 150)
        pairs = random_pairs(n, prng.randint(0, 25), crossing=prng.random() < 0.8, max_len=prng.randint(1, 8))
        prng.shuffle(pairs)
        pairs = [(j, i) if prng.random() < 0.2 else (i, j) for i, j in pairs]
        if pairs and prng.random() < 0.15:
            pairs.append(prng.choice(pairs))
        if prng.random() < 0.05 and pairs:
            i, j = prng.choice(pairs)
            pairs.append((i, prng.randrange(-1, n + 1)))
        breaks = sorted(prng.sample(range(-1, n + 2), prng.randint(0, 5))) if prng.random() < 0.5 else None
        out.append(write_case(pairs, n, breaks))
    # many crossing helices of random lengths (deep pseudoknot levels)
    for count in (10, 20, 28, 31, 40):
        pairs = []
        for k in range(count):
            pairs.append((k, count + k))
        n = 2 * count + 3
        prng.shuffle(pairs)
        out.append(write_case(pairs, n))
    return out


def stems_cases():
    lists = [[], [(0, 5), (1, 4)], [(0, 9), (1, 8), (3, 6), (2, 7)], [(0, 5), (2, 3)], [(0, 3), (1, 2), (4, 7)],
             [(0, 4), (2, 6), (1, 5), (3, 7)], [(1, 2), (0, 3)]]
    for _ in range(40):
        n = prng.randint(2, 80)
        lists.append(random_pairs(n, prng.randint(1, 10)))
    return [{"pairs": pairs_json(p), "stems": [pairs_json(s) for s in S._stems(p)]} for p in lists]


def crosses_cases():
    out = []
    for a in itertools.combinations(range(5), 2):
        for b in itertools.permutations(range(5), 2):
            out.append({"a": list(a), "b": list(b), "crosses": S._crosses(a, tuple(b))})
    return out


def unique_cases():
    lists = [[(0, 5), (1, 5), (2, 4)], [], [(0, 1), (0, 1)], [(1, 0), (0, 2)], [(3, 3), (3, 4)],
             [(0, 1), (2, 3), (1, 2), (4, 0)]]
    for _ in range(60):
        n = prng.randint(2, 30)
        lists.append([(prng.randrange(n), prng.randrange(n)) for _ in range(prng.randint(1, 25))])
    out = []
    for p in lists:
        keep, dropped = S.unique_pairs(p)
        out.append({"pairs": pairs_json(p), "keep": pairs_json(keep), "dropped": pairs_json(dropped)})
    return out


# --- detector ------------------------------------------------------------------------------------------

def f32(res):
    """Coordinates rounded to float32 (see the module docstring)."""
    return {k: np.asarray(v, dtype=np.float32).astype(np.float64) for k, v in res.items()}


def nucleotide(base, shift=(0.0, 0.0, 0.0)):
    return f32(TB.nucleotide(base, shift))


def partner(first, fbase, sbase, d=2.9, atoms=None):
    return f32(TB.partner(first, fbase, sbase, d=d, atoms=atoms))


def moved(res, shift):
    return {k: v + np.asarray(shift, dtype=np.float64) for k, v in res.items()}


KEEP_ATOMS = S.BASE_PAIR_ATOMS | {"C1'", "C1*", "P", "LA"}


def detector_case(name, residues_frames, cutoff=3.5, wobble=True, resnames=None, index=None, full=False):
    """residues_frames: list of frames, each [(resname, {atom: xyz}), ...] with the same atoms;
    or (resnames, index, [xyz, ...]) given directly. Unless ``full``, atoms the detector never reads
    (all but the H-bond atoms, C1', C1* and P) are left out to keep the file small."""
    if resnames is None:
        if not full:
            residues_frames = [[(n, {a: v for a, v in r.items() if a in KEEP_ATOMS}) for n, r in fr]
                               for fr in residues_frames]
        resnames, index, _ = TB.detector_input(residues_frames[0])
        frames = [TB.detector_input(f)[2] for f in residues_frames]
    else:
        frames = residues_frames
    det = S.PairDetector(resnames, index, cutoff=cutoff, wobble=wobble)
    out_frames = []
    for xyz in frames:
        xyz = np.asarray(xyz, dtype=np.float64).reshape(-1, 3)
        cand = det.candidates(xyz)
        pairs = det.detect(xyz)
        out_frames.append({
            "xyz": [fhex(v) for v in xyz.ravel()],
            "candidates": [[int(r[0]), int(r[1]), int(r[2]), fhex(r[3]), fhex(r[4])] for r in cand],
            "pairs": pairs.tolist(),
        })
    return {"name": name, "resnames": list(resnames),
            "atoms": [[[a, int(r)] for a, r in idx.items()] for idx in index],
            "cutoff": fhex(cutoff), "wobble": wobble, "frames": out_frames}


def detector_cases():
    out = []
    PA = TB.PAIR_ATOMS
    # test_canonical_pairs_detected
    for p, q, names in [("G", "C", ("G", "C")), ("A", "U", ("A5", "U3")), ("G", "U", ("RG", "RU")),
                        ("G", "C", ("0G", "LC5")), ("A", "U", ("LA", "0U")), ("A", "U", ("DA", "DT")),
                        ("G", "C", ("DG5", "DC3"))]:
        first = nucleotide(p)
        out.append(detector_case(f"canonical {names}", [[(names[0], first), (names[1], partner(first, p, q))]],
                                 full=True))
    # NA-01: L-DNA names, the CCD's own atom names of 0DA (N1A, N6A ...), the backbone guard
    own = {std: name for name, std in improvements_nucleotides.STANDARD_NAMES["0DA"].items()}
    a0 = nucleotide("A")
    t0 = partner(a0, "A", "U")
    out.append(detector_case("L-DNA 0DA with CCD atom names", [[("0DA", {own.get(k, k): v for k, v in a0.items()}),
                                                                ("0DT", t0)]], full=True))
    out.append(detector_case("L-DNA LDA with a backbone", [[("LDA", a0), ("LDT", t0)]]))
    out.append(detector_case("LDA without a backbone (a detergent)",
                             [[("LDA", {k: v for k, v in a0.items() if "'" not in k and k != "P"}), ("U", t0)]]))
    out.append(detector_case("LA without a backbone (a lanthanum ion)",
                             [[("LA", {k: v for k, v in a0.items() if "'" not in k and k != "P"}), ("U", t0)]]))
    # test_wobble_optional_and_cutoff
    g = nucleotide("G")
    u = partner(g, "G", "U")
    out.append(detector_case("wobble off", [[("G", g), ("U", u)]], wobble=False))
    out.append(detector_case("wobble on", [[("G", g), ("U", u)]]))
    far = partner(g, "G", "U", d=3.8)
    for cutoff in (3.5, 4.0, 3.0, math.inf):  # inf: missing atoms (distance inf) count as bonds
        out.append(detector_case(f"wobble far cutoff {cutoff}", [[("G", g), ("U", far)]], cutoff=cutoff))
    # test_stacked_neighbours_are_not_pairs
    c = partner(g, "G", "C")
    c2 = dict(c)
    c2["C1'"] = g["C1'"] + np.array([5.5, 0.0, 0.0])
    out.append(detector_case("stacked", [[("G", g), ("C", c2)]]))
    c3 = dict(c2)
    del c3["C1'"]
    out.append(detector_case("stacked without C1'", [[("G", g), ("C", c3)]]))
    c4 = {("C1*" if k == "C1'" else k): v for k, v in c.items()}
    out.append(detector_case("C1* name", [[("G", g), ("C", c4)]]))
    g5 = {("C1*" if k == "C1'" else k): v for k, v in g.items()}
    g5["C1'"] = g["C1'"] + np.array([0.0, 30.0, 0.0])  # C1' wins over C1*
    out.append(detector_case("C1' before C1*", [[("G", g5), ("C", c)]]))
    # C1'-C1' exactly at the window edges and one ulp beyond
    for dist in (8.5, 12.5, np.nextafter(8.5, 0), np.nextafter(12.5, 20), 10.5, 8.0, 13.0):
        ga = dict(g)
        cb = dict(c)
        ga["C1'"] = np.array([0.0, 0.0, 0.0])
        cb["C1'"] = np.array([float(dist), 0.0, 0.0])
        out.append(detector_case(f"C1' window {float(dist)!r}", [[("G", ga), ("C", cb)]]))
    # H-bond distances exactly at the cutoff (inclusive)
    names, index, xyz = TB.detector_input([("G", g), ("C", c)])
    d = sorted(float(np.linalg.norm(xyz[index[0][a]] - xyz[index[1][b]])) for a, b in PA[("G", "C")])
    for cut in (d[1], np.nextafter(d[1], 0), d[0], d[2], np.nextafter(d[2], 0)):
        out.append(detector_case(f"cutoff at distance {float(cut)!r}", [[("G", g), ("C", c)]], cutoff=float(cut)))
    # test_gc_needs_two_of_three_bonds
    one_off = dict(c)
    one_off["O2"] = c["O2"] + 5.0 * np.array([0.0, 0.0, 1.0])
    two_off = dict(one_off)
    two_off["N4"] = c["N4"] + 5.0 * np.array([0.0, 0.0, 1.0])
    out.append(detector_case("gc two of three", [[("G", g), ("C", c)], [("G", g), ("C", one_off)],
                                                 [("G", g), ("C", two_off)]]))
    # missing atoms
    for drop in (["O2"], ["O2", "N4"], ["N3"], ["N3", "O2", "N4"]):
        cm = {k: v for k, v in c.items() if k not in drop}
        out.append(detector_case(f"C without {drop}", [[("G", g), ("C", cm)]]))
        out.append(detector_case(f"C without {drop} cutoff inf", [[("G", g), ("C", cm)]], cutoff=math.inf))
    for drop in (["N2"], ["N1", "O6"], ["C1'"]):
        gm = {k: v for k, v in g.items() if k not in drop}
        out.append(detector_case(f"G without {drop}", [[("G", gm), ("C", c)]]))
    a = nucleotide("A")
    ua = partner(a, "A", "U")
    for drop in (["N6"], ["O4"]):
        out.append(detector_case(f"A-U without {drop}", [[("A", {k: v for k, v in a.items() if k not in drop}),
                                                          ("U", {k: v for k, v in ua.items() if k not in drop})]]))
    # test_conflicts_resolved_by_distance
    close = partner(g, "G", "C", d=2.9)
    loose = partner(g, "G", "C", d=3.3)
    out.append(detector_case("conflict", [[("C", loose), ("G", g), ("C", close)]]))
    g2 = nucleotide("G", shift=(0.0, 0.0, 30.0))
    out.append(detector_case("conflict loser pairs elsewhere",
                             [[("C", loose), ("G", g), ("C", close), ("G", g2), ("C", partner(g2, "G", "C"))]]))
    # near-ties: a second partner at distances around the first's
    for dl in [2.80 + 0.005 * k for k in range(41)]:
        other = partner(g, "G", "C", d=dl)
        out.append(detector_case(f"tie sweep {dl:.3f}", [[("C", other), ("G", g), ("C", close)]]))
    # test_non_nucleotides_never_pair
    for name in ("LIG", "ALA", "MG", "X"):
        out.append(detector_case(f"non-nucleotide {name}", [[("G", g), (name, c)]]))
    out.append(detector_case("lanthanum", [[("A", nucleotide("A")), ("LA", {"LA": np.zeros(3)})]]))
    # duplex frames (analysis tests), with an ALA and L-RNA names
    duplex = [[(n, f32(r)) for n, r in fr] for fr in TB.duplex_frames()]
    out.append(detector_case("duplex", duplex, full=True))
    ala = {"N1": np.array([0.0, 0.0, 50.0]), "N3": np.array([1.0, 0.0, 50.0])}
    out.append(detector_case("duplex + ALA", [f + [("ALA", ala)] for f in duplex]))
    out.append(detector_case("L-RNA", [[("LG5", g), ("0C", c)]] * 2))
    # every residue-name scheme of the four pairs
    schemes = {
        "G": ["G", "G5", "G3", "GN", "RG", "RG5", "RG3", "RGN", "GUA", "LG", "LG5", "LG3", "LGN", "0G", "0G5",
              "0G3", "DG", "DG5", "DG3", "DGN", "g", " rg ", "lg3"],
        "C": ["C", "C5", "C3", "CN", "RC", "RC3", "CYT", "LC", "LC5", "0C", "0C3", "DC", "DC5", "DC3", "c"],
        "A": ["A", "A5", "A3", "AN", "RA", "RA5", "ADE", "LA", "LA3", "0A", "0A5", "DA", "DA3", "DAN", "a"],
        "U": ["U", "U5", "U3", "UN", "RU", "RU3", "URA", "URI", "LU", "LU5", "0U", "0U3", "u", "DU", "RT",
              "T", "DT", "DT5", "DT3", "THY", "LT", "0T", "t"],
    }
    for fb, sb in (("G", "C"), ("A", "U"), ("G", "U")):
        first = nucleotide(fb)
        second = partner(first, fb, sb)
        names2 = schemes[sb] + (schemes["C"] if sb == "C" else [])
        residues = []
        for k, n1 in enumerate(schemes[fb]):
            n2 = names2[k % len(names2)]
            shift = np.array([0.0, 0.0, 40.0 * k])
            residues.append((n1, moved(first, shift)))
            residues.append((n2, moved(second, shift)))
        out.append(detector_case(f"schemes {fb}-{sb}", [residues]))
        residues = []
        for k, n2 in enumerate(schemes["U"] if sb == "U" else names2):
            shift = np.array([40.0 * k, 0.0, 0.0])
            residues.append((schemes[fb][k % len(schemes[fb])], moved(first, shift)))
            residues.append((n2, moved(second, shift)))
        out.append(detector_case(f"schemes {fb}-{sb} partners", [residues]))
    # non-finite coordinates
    cn = dict(c)
    cn["O2"] = np.array([np.nan, 0.0, 0.0])
    out.append(detector_case("NaN H-bond atom", [[("G", g), ("C", cn)]]))
    cn = dict(c)
    cn["C1'"] = np.array([np.inf, 0.0, 0.0])
    out.append(detector_case("infinite C1'", [[("G", g), ("C", cn)]]))
    gn = dict(g)
    gn["N1"] = np.array([np.nan, np.nan, np.nan])
    gn["O6"] = np.array([-np.inf, 1.0, 1.0])
    out.append(detector_case("NaN and inf G", [[("G", gn), ("C", c)]]))
    # perturbed helices: several stacked pairs, noisy frames (float64 perturbations)
    out.append(detector_case("perturbed duplex", perturbed_frames(duplex[0], [0.05, 0.1, 0.2, 0.3, 0.5, 0.8])))
    helix = helix_residues(12)
    out.append(detector_case("perturbed helix", perturbed_frames(helix, [0.0, 0.1, 0.25, 0.4, 0.6, 1.0, 1.5])))
    out.append(detector_case("perturbed helix wobble off", perturbed_frames(helix, [0.2, 0.5]), wobble=False))
    out.append(detector_case("perturbed helix cutoff 3.0", perturbed_frames(helix, [0.2, 0.5]), cutoff=3.0))
    out.append(detector_case("perturbed helix cutoff 4.5", perturbed_frames(helix, [0.2, 0.5]), cutoff=4.5))
    # dense random clouds: many competing candidates
    for k, (n_res, box) in enumerate([(12, 6.0), (30, 9.0), (60, 12.0), (40, 7.0), (50, 6.0)]):
        out.append(random_cloud_case(f"random cloud {k}", n_res, box, n_frames=4))
    return out


def helix_residues(n_pairs):
    """n_pairs paired residues stacked roughly along z: G-C, A-U, G-U and reversed in turn."""
    kinds = [("G", "C"), ("A", "U"), ("G", "U"), ("C", "G"), ("U", "A"), ("U", "G")]
    first_strand, second_strand = [], []
    for k in range(n_pairs):
        fb, sb = kinds[k % len(kinds)]
        p, q = (fb, sb) if (fb, sb) in TB.PAIR_ATOMS else (sb, fb)
        x = nucleotide(p)
        y = partner(x, p, q, d=2.75 + 0.05 * (k % 5))
        angle = math.radians(33.0 * k)
        rot = np.array([[math.cos(angle), -math.sin(angle), 0.0], [math.sin(angle), math.cos(angle), 0.0],
                        [0.0, 0.0, 1.0]])
        shift = np.array([0.0, 0.0, 3.4 * k])
        x = f32({n: rot @ v + shift for n, v in x.items()})
        y = f32({n: rot @ v + shift for n, v in y.items()})
        a, b = ((p, x), (q, y)) if (fb, sb) == (p, q) else ((q, y), (p, x))
        first_strand.append(a)
        second_strand.append(b)
    return first_strand + second_strand[::-1]


def perturbed_frames(residues, sigmas):
    frames = []
    for sigma in sigmas:
        frames.append([(n, {a: v + rng.normal(scale=sigma, size=3) for a, v in r.items()}) for n, r in residues])
    return frames


def random_cloud_case(name, n_res, box, n_frames):
    names = []
    index = []
    rows = 0
    pool = ["G", "C", "A", "U", "DT", "LG", "0C", "GUA", "ALA", "DA", "RU"]
    for r in range(n_res):
        resname = prng.choice(pool)
        base = S.base_letter(resname) or "G"
        atoms = [a for a in ["C1'"] + [x for x, *_ in TB.BASES.get(base if base != "T" else "U", [])]
                 if a in S.BASE_PAIR_ATOMS or a == "C1'"]
        atoms = [a for a in atoms if prng.random() > 0.08]  # some atoms missing
        if prng.random() < 0.05:
            atoms = ["C1*" if a == "C1'" else a for a in atoms]
        names.append(resname)
        index.append({a: rows + k for k, a in enumerate(atoms)})
        rows += len(atoms)
    frames = []
    for _ in range(n_frames):
        # clustered points: residue centres spread over the box, atoms within ~4 A of them; C1' 5-6 A away
        xyz = np.zeros((rows, 3))
        for r in range(n_res):
            centre = rng.uniform(0.0, box, size=3)
            for a, row in index[r].items():
                spread = 5.5 if a in ("C1'", "C1*") else 1.8
                xyz[row] = centre + rng.normal(scale=spread, size=3)
        frames.append(xyz)
    return detector_case(name, frames, resnames=names, index=index)


def select_cases():
    """detect() with candidates replaced (MDATK's tie test)."""
    WC = S.WATSON_CRICK
    lists = [
        ("mdatk tie test", 3, [[0, 1, WC, 5.63, 1.3], [0, 2, WC, 5.63, 0.0]]),
        ("equal dev keeps order", 3, [[0, 1, WC, 5.63, 0.5], [0, 2, WC, 5.63, 0.5]]),
        ("equal dev keeps order 2", 3, [[0, 2, WC, 5.63, 0.5], [0, 1, WC, 5.63, 0.5]]),
        ("reversed rows", 4, [[3, 0, 2, 5.7, 0.1], [1, 2, 1, 5.6, 0.2]]),
        ("empty", 3, []),
    ]
    out = []
    for name, n, cand in lists:
        out.append(select_case(name, n, np.array(cand, dtype=np.float64).reshape(-1, 5)))
    # scores on and next to grid points and half-way points of 0.05
    grid = []
    for k in range(100, 140):
        for base in (k * 0.05, (k + 0.5) * 0.05, k / 20, (k + 0.5) / 20, k * S.SCORE_TIE + 0.025):
            grid += [base, np.nextafter(base, 0), np.nextafter(base, 10)]
    for t in range(40):
        n = 12
        m = 30
        cand = np.zeros((m, 5))
        for r in range(m):
            i = prng.randrange(n)
            j = prng.randrange(n)
            while j == i:
                j = prng.randrange(n)
            cand[r] = [i, j, prng.choice([1, 2]), prng.choice(grid), prng.choice([0.0, 0.25, 0.5, 1.0, 10.5,
                                                                                    prng.random() * 2])]
        out.append(select_case(f"grid {t}", n, cand))
    return out


def select_case(name, n, cand):
    det = S.PairDetector(["G"] * n, [{} for _ in range(n)])
    det.candidates = lambda xyz: cand
    pairs = det.detect(np.zeros((0, 3)))
    return {"name": name, "n": n,
            "candidates": [[int(r[0]), int(r[1]), int(r[2]), fhex(r[3]), fhex(r[4])] for r in cand],
            "rounded": [fhex(v) for v in np.round(cand[:, 3] / S.SCORE_TIE) * S.SCORE_TIE],
            "pairs": pairs.tolist()}


# --- lattice (rebuilt in C++) ------------------------------------------------------------------------

LATTICE_NAMES = ["G", "C", "C", "C", "A", "U", "U", "DT", "GUA", "LC5"]  # residue r: LATTICE_NAMES[r % 10]


def lattice_atoms(resname):
    base = S.base_letter(resname)
    atoms = {"G": ["N1", "N2", "O6", "C1'"], "C": ["N3", "O2", "N4", "C1'"], "A": ["N6", "N1", "C1'"],
             "U": ["O4", "N3", "O2", "C1'"], "T": ["O4", "N3", "O2", "C1'"]}[base]
    return atoms


def lattice_xyz(n_res):
    """Residue r, atom slot a (in lattice_atoms order), row = running count:
    x = ((37 r + 11 a) mod 53) / 4, y = ((17 r + 5 a) mod 47) / 4, z = ((29 r + 3 a) mod 41) / 4."""
    names, index, rows = [], [], []
    for r in range(n_res):
        name = LATTICE_NAMES[r % len(LATTICE_NAMES)]
        names.append(name)
        idx = {}
        for a, atom in enumerate(lattice_atoms(name)):
            idx[atom] = len(rows)
            rows.append([((37 * r + 11 * a) % 53) / 4, ((17 * r + 5 * a) % 47) / 4, ((29 * r + 3 * a) % 41) / 4])
        index.append(idx)
    return names, index, np.array(rows, dtype=np.float64)


def lattice_case(n_res, cutoff):
    names, index, xyz = lattice_xyz(n_res)
    det = S.PairDetector(names, index, cutoff=cutoff)
    cand = det.candidates(xyz)
    pairs = det.detect(xyz)
    return {"n_residues": n_res, "cutoff": fhex(cutoff), "n_candidates": len(cand), "n_pairs": len(pairs),
            "candidates_sha256": hashlib.sha256(np.ascontiguousarray(cand, dtype="<f8").tobytes()).hexdigest(),
            "pairs_sha256": hashlib.sha256(np.ascontiguousarray(pairs, dtype="<i8").tobytes()).hexdigest(),
            "first_candidates": [[int(r[0]), int(r[1]), int(r[2]), fhex(r[3]), fhex(r[4])] for r in cand[:20]],
            "first_pairs": pairs[:20].tolist()}


def main():
    data = {
        "constants": constants(),
        "residue_names": residue_names(),
        "lw_kind": lw_kinds(),
        "parse": parse_cases(),
        "write": write_cases(),
        "stems": stems_cases(),
        "crosses": crosses_cases(),
        "unique_pairs": unique_cases(),
        "detector": detector_cases(),
        "select": select_cases(),
        "lattice": [lattice_case(3000, 3.5), lattice_case(400, 2.5)],
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "secstruct_cases.json"), "w", encoding="utf-8", newline="\n") as fh:
        # one line per case: small diffs, compact file
        fh.write("{\n")
        for n, (key, value) in enumerate(data.items()):
            fh.write(json.dumps(key) + ": ")
            if isinstance(value, list):
                fh.write("[\n" + ",\n".join(dumps(v) for v in value) + "\n]")
            else:
                fh.write(dumps(value))
            fh.write(",\n" if n + 1 < len(data) else "\n")
        fh.write("}\n")


def dumps(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


if __name__ == "__main__":
    main()
