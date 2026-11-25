"""What MDATK's selection layer computes (tests/data/selections/selections_cases.json).

    /opt/mdbridge-venv/bin/python -I -B make_selections_cases.py MDATK_CHECKOUT OUT_DIR

Every function of mdatk/selections.py and resolve_export_selection / frames_from_specs of
mdatk/io/pdb_export.py on edge cases: Unicode whitespace and digits, chains, numbers beyond 64 bits,
nested / cyclic / duplicate residue groups, frame and time lists. Results are JSON values; errors
are {"error": str(exc)}. Python ints beyond 64 bits are written as {"$int": "<digits>"}. Systems
are the fixtures of tests/data/topology (by file name) with the frame of the dump. Needs MDATK's
dependencies (numpy, MDAnalysis 2.11.0); nothing is written into the checkout. Where MDBridge
deliberately differs (SEL-13, SEL-14, SEL-15, SEL-19 ... SEL-22 and the atom selections of
docs/IMPROVEMENTS.md), improvements_selections.patch_mdatk replaces MDATK's function with MDBridge's
behaviour.
"""

import json
import logging
import os
import sys
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_selections_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import improvements_selections  # noqa: E402
import independent_tables  # noqa: E402

from mdatk import selections as S  # noqa: E402

# MDBridge's deliberate differences (improvements_selections.py, docs/DIFFERENCES.md): group
# expansion, select() and atom_elements as MDBridge has them, before pdb_export binds select
improvements_selections.patch_mdatk(S)
import improvements_core_io  # noqa: E402

improvements_core_io.apply_export()  # the export subset of commas or unknown groups (EXP-01, MDATK 0.0.10)

from mdatk.config import ResidueGroup  # noqa: E402
from mdatk.io.pdb_export import frames_from_specs, resolve_export_selection  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)

TOPOLOGY = os.path.join(HERE, "..", "data", "topology")
FIXTURES = ["fixture_chains.pdb", "fixture_noelements.pdb", "fixture_segids.pdb", "fixture_chains_ba.pdb",
            "demo_dry.parm7"]


def plain(value):
    if isinstance(value, (bool, type(None), str)):
        return value
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        value = int(value)
        return value if -2**63 <= value < 2**63 else {"$int": str(value)}
    if isinstance(value, (float, np.floating)):
        return float(value)
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    raise TypeError(type(value))


def call(fn, *args, **kwargs):
    try:
        return {"result": plain(fn(*args, **kwargs))}
    except Exception as exc:
        return {"error": str(exc), "type": type(exc).__name__}


def groups_json(groups):
    return [{"name": g.name, "residues": g.residues, "selection": g.selection} for g in groups]


RANGES = ["", "   ", "1-40, 45 ,50 - 60", "A:1-4;B:7", "3:1", "1-x", "1", "1:3", "1 : 3", "1 -3", "A : 1 - 3",
          "A1:5", "a:5:6", "AB12:7-9", "1A:5", "_A:5", "A:", ":5", "5-", "-5", "1--3", "1-3-5", "1,,2;;3 4",
          "\u0661-\u0663", "A:\u0662", "007-010", "99999999999999999999", "1-99999999999999999999",
          "99999999999999999999-1", "1\u2003-\u20033", "1\u00a0,2", "\t1\n-\n2\t", "1\u3000\u30002",
          "A:1-3, B:4-6, 7", "\u00e9:1", "A:1 - B:3", "0", "0-0", "12-12", "A:3:1", "1-2:3", "1-\u0663",
          "１-３", "2147483648", "9223372036854775807", "9223372036854775808", "1;2,3 4\t5"]

ELEMENT_NAMES = [("C4'", "A", False), ("HO2'", "G", False), ("MG", "MG", True), ("CA", "GLY", False),
                 ("Na+", "Na+", True), ("CA", "CA", True), ("CL", "CL", True), ("1HB", "ALA", False),
                 ("123", "X", False), ("", "", False), ("", "", True), (" mg ", " mg2 ", True),
                 ("MG2+", "XYZ", True), ("X", "K+", True), ("Zn", "ZN2", True), ("ZN", "ZN", False),
                 ("C*", "", False), ("'", "", False), ("1'", "", False), ("o1p", "g", False), ("N\u00e9", "", False),
                 ("_N", "", False), ("1_N", "", False), ("\u0661N", "", False),
                 ("CLA", "CLA", True), ("SOD", "SOD", True), ("POT", "", True), ("FE3", "HEM", True),
                 ("BR-", "", True), ("Br", "BR", True), ("LI+", "LI", True), ("CS", "CS", True),
                 ("RB+", "", True), ("MN", "MN2", True), ("OW", "HOH", False), ("H1", "HOH", True),
                 ("\tNA\n", "", True), ("\u2003P\u2003", "", False), ("+1C", "", False)]

FRAME_SPECS = ["", "  ", "0, 10, 20-50:10, last, -1", "100", "0", "last", "end", "first", "start", "LAST",
               "End", "-1", "-100", "-101", "5-2", "2-5", "2-5:2", "2-5:0", "2-last", "2-last:3", "0-last:7",
               "1-", "-1-3", "a", "1.5", "1e2", "+5", "1_0", "\u0663", "\u0661-\u0663", "1-\u0663:\u0661",
               "3;4,5\t6\u20037", "0,0,0", "99999999999999999999", "-99999999999999999999", "1, 200, -300",
               "2-5:", "2-5:x", "2-:3", "lastly", "0-99", "50-200", "10-20:3, 15", " 1 ,, 2 ", "0x10",
               "2-5:99999999999999999999", "07", "-0", "1 - 3"]

TIME_SPECS = ["", "1.5ns, 2000ps, 3", "5 parsecs", "1ns", "1NS", "1 ns", ".5", "5.", "5.ns", "-1ps",
              "+2fs", "1e3", "1e3ps", "1E-1ns", "1e", "1ens", "us", "3us", "3\u00b5s", "\u0663ns", "1_0ps",
              "1.5.5", "1,2;3", "nan", "inf", "1ms", "2Ps", "1e400", "0.0001fs", "12 34", "1.2e+3NS"]


def system(fixture):
    return mda.Universe(os.path.join(TOPOLOGY, fixture))


def main():
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    out = {"mdanalysis": mda.__version__}

    out["parse_residue_ranges"] = [{"text": t, **call(S.parse_residue_ranges, t)} for t in RANGES]
    out["ranges_to_selection"] = []
    for fixture in (None, "fixture_chains.pdb", "demo_dry.parm7"):
        u = system(fixture) if fixture else None
        for t in RANGES:
            out["ranges_to_selection"].append({"text": t, "fixture": fixture, **call(S.ranges_to_selection, t, u)})

    group_sets = {
        "basic": [ResidueGroup("A", "1-5"), ResidueGroup("B", "", "group:A and name P"),
                  ResidueGroup("Loop1", "", "group:Loop2"), ResidueGroup("Loop2", "", "group:Loop1")],
        "chains": [ResidueGroup("Aptamer", "A:1-4"), ResidueGroup("Target", "B:10-12"),
                   ResidueGroup("Prot", "C:5-6, 201"), ResidueGroup("Both", "", " group:Aptamer or group:Target ")],
        "duplicates": [ResidueGroup("X", "1-2"), ResidueGroup("Y", "3"), ResidueGroup("X", "", "resid 10"),
                       ResidueGroup("Empty", "", "  "), ResidueGroup("Bad", "1-x"), ResidueGroup("Sel", "1", "name P")],
        "deep": [ResidueGroup(f"G{i}", "", f"group:G{i + 1}") for i in range(8)] + [ResidueGroup("G8", "1")],
        "deeper": [ResidueGroup(f"G{i}", "", f"group:G{i + 1}") for i in range(9)] + [ResidueGroup("G9", "1")],
        "odd": [ResidueGroup("a.b+c-d'e_f", "2"), ResidueGroup("g", "", "group:a.b+c-d'e_f or mygroup:g2"),
                ResidueGroup("g2", "3"), ResidueGroup("sel:x", "4"), ResidueGroup(" pad ", "5"),
                ResidueGroup("\u00e9", "6")],
        "none": [],
        # SEL-19 ... SEL-22: hints, suggestions, circles, empty groups
        "typos": [ResidueGroup("Aptamer", "A:1-4"), ResidueGroup("Target", "B:10-12"),
                  ResidueGroup("Uses", "", "group:Tagret and name P"), ResidueGroup("Far", "900-910"),
                  ResidueGroup("FarSel", "", "resid 900:910"), ResidueGroup("Circle1", "", "group:Circle2 or name P"),
                  ResidueGroup("Circle2", "", "resid 1 and group:Circle3"), ResidueGroup("Circle3", "", "group:Circle1"),
                  ResidueGroup("protein", "5")],
    }
    texts = ["group:A", "group:B or resid 9", "A", "B", "sel:name P", "sel: group:A ", "group:Nope", "group:Loop1",
             "", "   ", "group:", "group:A)", "(group:A)and(group:B)", "xgroup:Ay", "group:A.B", "group:Aptamer",
             "Both", "Target", "group:Both and name P", "X", "group:X", "group:Empty", "group:Bad", "Sel",
             "G0", "group:G0", "group:G1", "group:G8", "a.b+c-d'e_f", "g", "group:g", "sel:x", "sel:",
             "sel:sel:x", " pad ", "pad", "\u00e9", "group:\u00e9", "\u2003A\u2003", "group:Agroup:B",
             "group:A group:A", "group:G9", "Loop2", "group:Tagret", "group:target", "Uses", "group:Circle1",
             "Circle3", "group:Loop2 or group:Nope"]
    out["groups"] = {}
    for name, groups in group_sets.items():
        cases = []
        for fixture in (None, "fixture_chains.pdb", "demo_dry.parm7"):
            u = system(fixture) if fixture else None
            for t in texts:
                cases.append({"text": t, "fixture": fixture,
                              "expand": call(S.expand_groups, t, groups, u),
                              "resolve": call(S.resolve_selection, t, groups, u)})
            for g in groups:
                cases.append({"group": g.name, "fixture": fixture, "group_selection": call(S.group_selection, g, u)})
        out["groups"][name] = {"groups": groups_json(groups), "cases": cases}

    select_texts = texts + ["name P", "name Q", "resid 1:3", "foo", "around 3 name P", "group:A and name P",
                            "sel:resid 200:206", "not all", "name", "resid 3A", "protein", "byres group:A",
                            # SEL-19: hints
                            "resnam G", "Target and name P", "name P and Aptamer", "group Target", "group Nope",
                            "group target", "(all x)", "(name P", "RESID 1", "namee P", "foo and name P",
                            "Proteins", "sel:Target",
                            # SEL-22: why nothing matched
                            "group:Far", "group:FarSel", "group:Far or group:FarSel", "Far", "sel:resid 900",
                            "name P and resid 999", "name P or resid 999", "name Q and resid 1:3",
                            "name P and name O1P", "(name P) and (name Q or resid 999)", "group:Aptamer and name Q",
                            "resid 999 and around 3 name P", "name P andname Q", "name Q and\u2003resid 1"]
    out["select"] = []
    for fixture in FIXTURES:
        u = system(fixture)
        for name, groups in group_sets.items():
            if name in ("deep", "deeper") and fixture != "fixture_chains.pdb":
                continue
            for t in select_texts:
                for allow_empty in (False, True):
                    r = call(lambda: S.select(u, t, groups, allow_empty).indices)
                    out["select"].append({"fixture": fixture, "groups": name, "text": t, "allow_empty": allow_empty, **r})
        for t in ["", "  ", "Aptamer, Target", "Aptamer,Target,", "Aptamer, Nope", "Aptamer", " Prot ,Both ",
                  "name P", "group:Target", ",", "Aptamer,,Target"]:
            r = call(lambda: resolve_export_selection(u, t, group_sets["chains"]).indices)
            out["select"].append({"fixture": fixture, "groups": "chains", "text": t, "export": True, **r})

    # core::suggestName / similarity (MDATK 0.1.x selections.suggest: difflib)
    import difflib
    words = ["resname", "resid", "name", "and", "or", "Target", "target", "TARGET", "Aptamer", "", "chainID",
             "chainid", "altLoc", "nucleicbackbone", "\u00e9t\u00e9", "\u0130x", "abcabcabc", "bcab", "x" * 199]
    names = ["resnam", "RESID", "nmae", "adn", "o", "Tagret", "tARGET", "aptamr", "", "chainjd", "CHAINID",
             "altlok", "nucleicbackbon", "\u00c9T\u00c9", "ix", "cab", "abc", "zzz", "x" * 199, "x" * 200,
             "Proteins", "nam", "rs", "bcabca"]
    out["suggest"] = [{"name": n, "candidates": words,
                       "result": improvements_selections.suggest(n, words)} for n in names]
    out["suggest"] += [{"name": n, "candidates": c, "result": improvements_selections.suggest(n, c)}
                       for n, c in [("ab", ["ba", "ab_"]), ("ab", ["ba", "aB"]), ("xy", ["xa", "ay"]),
                                    ("Ab", ["aB", "AB"]), ("q", [])]]
    out["similarity"] = [{"a": a, "b": b, "ratio": difflib.SequenceMatcher(None, a, b).ratio()}
                         for a in words for b in names if len(b) < 200]

    out["guess_element"] = [{"args": [n, r, s], **call(S.guess_element, n, r, s)} for n, r, s in ELEMENT_NAMES]
    out["atoms"] = []
    for fixture in FIXTURES:
        u = system(fixture)
        subsets = {"all": u.atoms, "first": u.atoms[:7], "ions": u.atoms[u.atoms.resnames != "XXX"][-8:],
                   "empty": u.atoms[[]], "reversed": u.atoms[::-3]}
        for label, ag in subsets.items():
            out["atoms"].append({
                "fixture": fixture, "subset": label, "indices": [int(i) for i in ag.indices],
                "elements": call(S.atom_elements, ag), "heavy_mask": call(S.heavy_mask, ag),
                "heavy": call(lambda: S.heavy(ag).indices),
                "filter_NO": call(lambda: S.element_filter(ag, ["n", " O ", ""]).indices),
                "filter_MG": call(lambda: S.element_filter(ag, ["Mg", "na", "cl", "K"]).indices),
                "filter_none": call(lambda: S.element_filter(ag, []).indices),
            })
    out["residue_labels"] = []
    for fixture in FIXTURES:
        u = system(fixture)
        for label, res in {"all": u.residues, "first": u.residues[:3], "last": u.residues[-4:],
                           "empty": u.residues[[]], "mixed": u.residues[[0, -1]]}.items():
            out["residue_labels"].append({"fixture": fixture, "residues": [int(r) for r in res.ix],
                                          **call(S.residue_labels, res)})
    out["split_csv"] = [{"text": t, **call(S.split_csv, t)} for t in
                        ["", "a, ,b", " a ,b,, c ", ",", "\u2003a\u2003,b", "a;b", "a,\nb\t", None]]
    out["parse_frame_spec"] = [{"text": t, "n": n, **call(S.parse_frame_spec, t, n)}
                               for t in FRAME_SPECS for n in (100, 1)]
    out["parse_time_spec"] = [{"text": t, "unit": unit, **call(S.parse_time_spec, t, unit)}
                              for t in TIME_SPECS for unit in ("ps", "ns", "NS", "frame", "fs", "min")]
    out["time_unit_factor"] = [{"unit": unit, **call(S.time_unit_factor, unit)}
                               for unit in ("ps", "ns", "us", "fs", "PS", "frame", "", "ms")]

    # frames_from_specs on in-memory trajectories (times as MDAnalysis reports them)
    out["frames_from_specs"] = []
    base = system("fixture_chains.pdb")
    for n, dt in ((10, 2.0), (1, 5.0), (7, 0.1), (25, 1.0)):
        u = system("fixture_chains.pdb")
        u.load_new(np.repeat(base.atoms.positions[None], n, axis=0), format=MemoryReader, dt=dt)
        times = [float(ts.time) for ts in u.trajectory]
        for frames, times_text, unit, dt_ps in [
                ("", "", "ps", None), ("0, last", "", "ns", None), ("", "4", "ps", None), ("1", "3ps", "ps", None),
                ("", "0.004ns, 7", "ps", None), ("", "100", "ps", None), ("", "-50", "ps", None),
                ("", "3", "frame", None), ("", "3", "ns", None), ("", "5", "ps", 2.5), ("", "5", "ps", 0.0),
                ("2", "1000fs", "ps", 0.5), ("", "0.05", "ps", None), ("", "x", "ps", None), ("99", "", "ps", None),
                ("", "1e9", "ps", None), ("", "0.15, 0.25", "ps", None)]:
            records = []

            class Capture(logging.Handler):
                def emit(self, record):
                    records.append(record.getMessage())

            handler = Capture()
            logging.getLogger("mdatk").addHandler(handler)
            logging.getLogger("mdatk").setLevel(logging.INFO)
            r = call(frames_from_specs, u, frames, times_text, unit, dt_ps)
            logging.getLogger("mdatk").removeHandler(handler)
            out["frames_from_specs"].append({"n": n, "dt": dt, "times": times, "frames": frames,
                                             "times_text": times_text, "unit": unit, "dt_ps": dt_ps,
                                             "warnings": [m for m in records], **r})
    # int() / float() of str, which the selection parsers use (mdbridge/core/python.hpp)
    numbers = ["0", "1", "-1", "+1", " 7 ", "\u20037\u3000", "007", "1_000", "1__000", "_1", "1_", "+_1",
               "\u0661\u0662", "\uff11\uff12", "\U0001d7ce", "1.5", "1e3", ".5", "5.", ".", "e3", "1e", "1e+",
               "1E-3", "1_0.5", "1._5", "1e1_0", "inf", "-Infinity", "+INF", "infinit", "nan", "-NaN", "nAn ",
               "0x10", "", " ", "-", "+", "--1", "1 2", "\x00", "1\x00", "\u00e9", "9223372036854775807",
               "-9223372036854775808", "9223372036854775808", "-9223372036854775809", "1" * 30, "-" + "9" * 25,
               "1e309", "-1e309", "1e-400", "4.9e-324", "2.4703282292062328e-324", "2.4703282292062327e-324",
               "1.7976931348623157e308", "1.7976931348623159e308", "0.1", "123456789012345678901234567890.5",
               "x" * 250, "\u00e9" * 210, "1\n", "\t-5\r\n", "\x1c3\x1f", "\x853", "١.٥", "1٫5", "１e２"]
    out["python_numbers"] = []
    for text in numbers:
        case = {"text": text}
        try:
            case["int"] = plain(int(text))
        except ValueError as exc:
            case["int_error"] = str(exc)
        try:
            case["float"] = repr(float(text))
        except ValueError as exc:
            case["float_error"] = str(exc)
        out["python_numbers"].append(case)
    with open(os.path.join(out_dir, "selections_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=0)
        fh.write("\n")


if __name__ == "__main__":
    main()
