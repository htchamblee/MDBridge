#!/usr/bin/env python3
"""Test cases of the step-type schemas (mdatk/workflows/schema.py, base.py, registry.py and the
settings comments of projectfile.py), computed by MDATK v0.0.4.1.

    python3 -I -B tests/parity/make_step_type_cases.py MDATK_CHECKOUT OUT_DIR

Run with Python 3.13 (float() parsing and the Unicode data depend on the version). Nothing is
written into the checkout; synthetic step types are registered in this process only.
Writes into OUT_DIR (tests/data/step-types), JSON, UTF-8:

* ``coerce_cases.json``        Option.coerce of many inputs for an option of every shape used by the
                               step types (and synthetic ones): {"inputs_<set>": [value], "options":
                               [{"option", "inputs": set, "extra": [value], "results": [{"out"} |
                               {"error"[, "exception"]}]}]}, one result per input of the set, then per
                               extra input; "repr": true means the message ends with repr(input),
                               "lossy" marks an integer result beyond 64 bits
* ``digit_cases.json``         float() of every Unicode decimal digit and around every whitespace
                               character (through an int option without limits)
* ``normalise_cases.json``     Tool.normalise of settings per step type: the result, the issues and
                               the log lines with an issue list, the result and log lines without one
* ``requirements_cases.json``  schema.requirements_met
* ``setting_comments.json``    projectfile._settings_doc for every key of every step type, of unknown
                               keys and step types, and of synthetic step types (search order)
* ``option_text_cases.json``   projectfile._option_text of synthetic labels, tooltips and choices
* ``registry_cases.json``      get_tool's KeyError text, register()'s ValueErrors, tools_by_category,
                               Section.template's KeyError, guess_kind and kind_label

Values are tagged JSON: {"t": "none" | "bool" | "int" | "float" | "str" | "datetime" | "date" |
"time" | "array" | "table", ...}; int "v" is a decimal string, float "hex" is float.hex() ("nan",
"inf", "-inf"), str "v" the text, date / time "iso" its isoformat(), table "v" a list of [key, value]
pairs.
Step types and options are given in the form of tests/data/step-types/step_types.json.
"""

from __future__ import annotations

import datetime as dt
import decimal
import json
import logging
import os
import sys
import types
import unicodedata

sys.dont_write_bytecode = True
if len(sys.argv) != 3:
    sys.exit("usage: make_step_type_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools"))

from gen_step_types import enc, option as dump_option, tool as dump_tool  # noqa: E402
from mdatk import projectfile as PF  # noqa: E402
from mdatk.workflows import registry, schema  # noqa: E402
from mdatk.workflows.base import Tool  # noqa: E402
from mdatk.workflows.schema import Option, Section  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
os.makedirs(OUT, exist_ok=True)
BUILTIN = list(registry.TOOLS.values())


def dump(name, obj):
    """JSON with one item of every top-level list per line."""
    def line(v):
        return json.dumps(v, ensure_ascii=False)

    def items(values, pad):
        return "[\n" + ",\n".join(pad + line(v) for v in values) + "\n" + pad[:-1] + "]"

    if isinstance(obj, list):
        text = items(obj, " ")
    else:
        text = "{\n" + ",\n".join(f" {line(k)}: " + (items(v, "  ") if isinstance(v, list) else line(v))
                                   for k, v in obj.items()) + "\n}"
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text + "\n")


class Capture(logging.Handler):
    def __init__(self):
        super().__init__(logging.DEBUG)
        self.lines = []

    def emit(self, record):
        self.lines.append(record.getMessage())


CAPTURE = Capture()
_base_log = logging.getLogger("mdatk.workflows.base")
_base_log.addHandler(CAPTURE)
_base_log.propagate = False
_base_log.setLevel(logging.DEBUG)


def fits_int64(v):
    return not isinstance(v, int) or isinstance(v, bool) or -2**63 <= v < 2**63


# ------------------------------------------------------------------------------------ coerce

def subnormal_strings():
    decimal.getcontext().prec = 2000
    tiny = decimal.Decimal(2) ** -1074
    out = []
    for k in (1, 2, 3, 2**52 - 1, 2**52):
        half = tiny * k + tiny / 2  # halfway between k and k + 1 units of the last place
        out += [str(half), str(half + decimal.Decimal("1e-1500")), str(half - decimal.Decimal("1e-1500"))]
    out += [str(tiny / 2), str(tiny / 2 + decimal.Decimal("1e-1400")), str(tiny * 7)]
    return out


NUMBERS = [None, True, False, 0, 1, -1, 2, 7, 100, 100000, 100001, -5, 2**53 + 1, 2**62, 2**63 - 1,
           -2**63, 0.0, -0.0, 1.0, 0.5, 1.5, 2.0, -1.0, 1e-07, 0.001, 1e300, -1e300, 1e16,
           float("nan"), float("inf"), float("-inf"), 9.5, 99999.99, 100000.5]
STRINGS = ["", " ", "\t\n", " ", "\x1c", "\x85", "0", "1", "5", " 5 ", "+5", "-5", "5.0", "5.5",
           "05", "007", "1e3", "1E3", "1e+3", "1e-3", "1e", "1e+", "e5", ".5", "5.", ".", "-.5",
           "+.5e-3", "1_000", "1__000", "_1", "1_", "1_.5", "1._5", "1.5_0", "1e1_0", "1e_10",
           "1_e5", "-_1", "0x10", "0b1", "0o7", "1j", "١٢", "١٢.٥",
           "٣_٤", " 5 ", "\x1c5", "5\x00", "\x7f5", "５", "\U0001D7D7",
           "inf", "-inf", "+inf", "Infinity", "-INFINITY", "iNfInItY", "infinit", "infinityy",
           "nan", "-nan", "+NaN", "nan(1)", "nana", "in f", "1 2", "--1", "+-1", "1e400", "-1e400",
           "1e-400", "-1e-400", "4.9e-324", "5e-324", "1e-310", "2.2250738585072011e-308",
           "2.2250738585072014e-308", "1.7976931348623157e308", "1.7976931348623158e308",
           "1.7976931348623159e308", "0." + "0" * 400 + "1", "1" + "0" * 400, "1" * 900,
           "0.5" + "0" * 900 + "1", "9007199254740993", "9007199254740993.0000000001",
           "9223372036854775807", "9223372036854775808", "1e18", "1e19", "-1e19", "123456789e-9",
           "0e999999999999999999999", "1e-999999999999999999999", "1e999999999999999999999",
           "true", "True", "TRUE", "yes", "Yes", "YES", "on", "ON", "off", "OFF", "no", "No",
           "false", "False", "maybe", "t", "y", " true", "true ", "K", "İ",
           "pmemd", "pmemd.MPI", "pmemd.cuda", " pmemd", "Pmemd", "Ⅷ", "²", "5²",
           "é", "\U0001F600", "a'b\"c", "\\", "it's", "line\nbreak"] + subnormal_strings()
OTHERS = [[], [1], ["a", 2.5], {}, {"a": 1}, dt.datetime(2026, 1, 2, 3, 4, 5),
          dt.datetime(2026, 1, 2, 3, 4, 5, 250000, tzinfo=dt.timezone(dt.timedelta(hours=2))),
          dt.date(2026, 1, 2), dt.time(1, 2, 3)]
SHORT = [None, True, False, 0, 1, 2, -0.0, 1.0, 1.5, float("nan"), "", " ", "a", "x y", "é",
         "it's", "1", "true", "True", "On", "pmemd.MPI", "26", "24", [], ["a"], {"k": "v"},
         dt.date(2026, 1, 2)]


def signature(o):
    return (o.kind, o.optional, o.minimum, type(o.minimum).__name__, o.maximum, type(o.maximum).__name__,
            tuple(c[0] for c in o.choices))


def synthetic_options():
    return [
        Option("plain_int", "Plain int", "int", 0),
        Option("plain_float", "Plain float", "float", 0.0),
        Option("opt_int", "Optional int", "int", None, optional=True, minimum=1),
        Option("neg_float", "Negative float", "float", -1.0, minimum=-273.15, maximum=-0.5),
        Option("half_int", "Int with a float minimum", "int", 1, minimum=0.5, maximum=1e7),
        Option("tiny_float", "Tiny float", "float", 1e-5, minimum=1e-05, maximum=123456789.0),
        Option("big_int", "Big int", "int", 10, minimum=-10**7, maximum=10**16),
        Option("zero_float", "Zero bounds", "float", 0.0, minimum=-0.0, maximum=0.0),
        Option("odd_choice", "Odd choice", "choice", "", choices=[("", "none"), ("1", "one"), ("True", "T"),
                                                                    ("1.0", "f"), ("é", "e")]),
        Option("plain_bool", "Bool", "bool", False),
        Option("plain_text", "Text", "text", ""),
    ]


def coerce_result(o, value, with_input=False):
    """{"out": value} or {"error": message, "exception": type}; an error that ends with the input's
    repr gives only the text before it and "repr": true (keeps the file small)."""
    case = {"in": enc(value)} if with_input else {}
    try:
        out = o.coerce(value)
    except (ValueError, OverflowError) as exc:
        text = str(exc)
        if text.endswith(repr(value)) and len(repr(value)) > 12:
            case.update(error=text[:-len(repr(value))], repr=True)
        else:
            case["error"] = text
        if not isinstance(exc, ValueError):
            case["exception"] = type(exc).__name__
        return case
    case["out"] = enc(out) if fits_int64(out) else {"t": "int", "v": str(out)}
    if not fits_int64(out):
        case["lossy"] = True
    return case


def coerce_cases():
    seen, chosen = set(), []
    for t in BUILTIN:
        for s in t.sections:
            for o in s.options:
                sig = signature(o)
                if sig not in seen:
                    seen.add(sig)
                    chosen.append(o)
    chosen += synthetic_options()
    sets = {"numeric": NUMBERS + STRINGS + OTHERS, "short": SHORT}
    out, count = [], 0
    for o in chosen:
        name = "numeric" if o.kind in ("bool", "int", "float") else "short"
        extra = [c[0] for c in o.choices] + [" " + c[0] for c in o.choices[:2]]
        results = [coerce_result(o, v) for v in sets[name] + extra]
        count += len(results)
        out.append({"option": dump_option(o), "inputs": name, "extra": [enc(v) for v in extra],
                    "results": results})
    dump("coerce_cases.json", {**{f"inputs_{k}": [enc(v) for v in vs] for k, vs in sets.items()},
                               "options": out})
    return count


def digit_cases():
    o = Option("n", "N", "int", 0)
    out = []
    for code in range(0x110000):
        c = chr(code)
        if unicodedata.category(c) == "Nd" or c.isdecimal():
            out.append(coerce_result(o, c, True))
        if c.isspace() or code in (0x1F, 0x7F, 0x200B, 0xFEFF):
            out.append(coerce_result(o, c + "5" + c, True))
            out.append(coerce_result(o, "4" + c + "5", True))
    dump("digit_cases.json", out)
    return len(out)


# ------------------------------------------------------------------------------------ normalise

BAD = {"bool": "maybe", "int": "x", "float": "1,5", "choice": "nope", "radio": "nope"}


def alternative(o):
    if o.kind == "bool":
        return "yes" if not o.default else "off"
    if o.kind in ("int", "float"):
        base = o.minimum if o.minimum is not None else 3
        return str(base) if o.kind == "float" else (int(base) if base == int(base) else "7")
    if o.kind in ("choice", "radio"):
        return o.choices[-1][0]
    return 42


def normalise_case(t, settings, label):
    CAPTURE.lines = []
    issues = []
    result = t.normalise(settings, issues)
    with_issues = CAPTURE.lines
    CAPTURE.lines = []
    result_without = t.normalise(settings)
    without = CAPTURE.lines
    assert result == result_without
    return {"tool": t.name, "label": label, "settings": enc(settings), "result": enc(result),
            "issues": [[i.level, i.message, list(i.keys)] for i in issues],
            "issue_texts": [str(i) for i in issues], "log_with_issues": with_issues,
            "log_without_issues": without}


def normalise_cases():
    out = []
    for t in BUILTIN:
        flat = [o for s in t.sections if not s.repeat for o in s.options]
        out.append(normalise_case(t, {}, "empty"))
        out.append(normalise_case(t, {o.key: alternative(o) for o in flat}, "alternatives"))
        bad = {o.key: BAD[o.kind] for o in flat if o.kind in BAD}
        bad.update({"zzz": 1, "aaa": "x", "ÿ": True, "Zed": None})
        out.append(normalise_case(t, bad, "bad values and unknown keys"))
        nulls = {o.key: None for o in flat}
        out.append(normalise_case(t, nulls, "None everywhere"))
        mixed = {o.key: [1] if k % 2 else {"a": 1} for k, o in enumerate(flat)}
        out.append(normalise_case(t, mixed, "lists and tables"))
        for s in t.sections:
            if not s.repeat:
                continue
            out.append(normalise_case(t, {s.key: None}, f"{s.key}: None"))
            out.append(normalise_case(t, {s.key: "stage"}, f"{s.key}: not a list"))
            out.append(normalise_case(t, {s.key: {"name": "x"}}, f"{s.key}: a table"))
            out.append(normalise_case(t, {s.key: []}, f"{s.key}: empty"))
            records = [{"name": "first"}, "text", 5, {o.key: BAD.get(o.kind, 7) for o in s.options},
                       {o.key: alternative(o) for o in s.options}, {"unknown": 1}, None, [1, 2]]
            out.append(normalise_case(t, {s.key: records, "zz": 0}, f"{s.key}: mixed records"))
            defaults = t.defaults()[s.key]
            out.append(normalise_case(t, {s.key: defaults}, f"{s.key}: defaults"))
    dump("normalise_cases.json", out)
    return len(out)


# ------------------------------------------------------------------------------------ requirements

def requirements_cases():
    values = [None, True, False, 0, 1, 2, 0.0, 1.0, -0.0, float("nan"), "", "a", "pmemd.MPI", "1", [],
              [1], ["a"], {}, {"a": 1}, dt.time(0, 0), dt.date(2026, 1, 1), 2**53, 2**53 + 1,
              float(2**53), 9.5]
    wanted = [True, False, "a", "pmemd.MPI", 1, 1.0, 0, ("a", "b"), ("pmemd.MPI",), (1, "x"), (True,),
              (0.0,), (), [1], ("a", [1]), None, 2**53 + 1, float(2**53), ["a"]]
    out = []
    for w in wanted:
        for v in values:
            for present in (True, False):
                req = {"k": w}
                vals = {"k": v, "other": 1} if present else {"other": v}
                out.append({"requires": enc(req), "values": enc(vals),
                            "met": schema.requirements_met(req, vals)})
    for w1 in (True, "a", ("a",)):
        for w2 in (False, 1, (2, 3)):
            for v1 in ("a", "", None):
                for v2 in (0, 1, 2, None):
                    req = {"x": w1, "y": w2}
                    vals = {"y": v2, "x": v1}
                    out.append({"requires": enc(req), "values": enc(vals),
                                "met": schema.requirements_met(req, vals)})
    dump("requirements_cases.json", out)
    return len(out)


# ------------------------------------------------------------------------------------ comments

def settings_doc(tool, key):
    for root in ("nodes", "scratch"):
        path = ("workflow", root, "settings", key)
        stack = [{}, {"nodes": []}, {"tool": tool}, {}]
        text = PF._settings_doc(path, 3, stack)
        if root == "nodes":
            first = text
        else:
            assert text == first
    return first


SEARCH_TOOLS = [
    # a repeated section matches before the options of later sections; options of repeated
    # sections are found too; a plain section's key is no setting
    ("parity-search-a", [
        Section("plain", "Plain", [Option("x", "X first", "text", "", "flat x"),
                                   Option("plain2", "", "text", "", "")], tooltip="plain section"),
        Section("rec", "Records", [Option("x", "X in records", "text", "", "record x"),
                                   Option("inner", "Inner", "bool", False, "only in records",
                                          choices=())], repeat=True, tooltip="rec tooltip:"),
        Section("other", "Other", [Option("rec", "Rec option", "text", "", "shadowed by the section")],
                repeat=True, tooltip=""),
        Section("late", "Late", [Option("y", "Y", "choice", "a", "pick", choices=[("a", "A"), ("b", "B")])]),
    ]),
    ("parity-search-b", [
        Section("first", "First", [Option("same", "", "text", "", ": colon tooltip ::")], repeat=True,
                tooltip="  spaced  "),
        Section("same", "Same", [Option("z", "Z", "text", "", "")], repeat=True, tooltip="section same"),
        Section("flat", "Flat", [Option("solo", "Solo", "text", "", "flat option")]),
        Section("rr", "Records", [Option("flat", "Flat option", "text", "", "named like a plain section")],
                repeat=True),
    ]),
]


def search_tools():
    made = []
    for name, sections in SEARCH_TOOLS:
        cls = type("Synthetic", (Tool,), {"name": name, "label": name, "category": "Parity", "tag": "x",
                                          "sections": sections})
        registry.register(cls)
        made.append(registry.TOOLS[name])
    return made


def setting_comments(synthetic):
    out = []
    for t in BUILTIN + synthetic:
        keys = []
        for s in t.sections:
            keys.append(s.key)
            keys += [o.key for o in s.options]
        keys += ["", "nope", "Name", " name"]
        for key in dict.fromkeys(keys):
            out.append({"tool": t.name, "key": key, "comment": settings_doc(t.name, key)})
    for tool in ("", "unknown", "Amber-build"):
        out.append({"tool": tool, "key": "engine", "comment": settings_doc(tool, "engine")})
    dump("setting_comments.json", {"synthetic": [dump_tool(t) for t in synthetic], "cases": out})
    return len(out)


def option_text_cases():
    labels = ["", "Label", " Label ", "Label:", ":Label", "  ", "a  b", " L"]
    tips = ["", "Tip.", " Tip. ", "Tip:", "ends with colon:", ": starts", "two  spaces", "trailing  ",
            " lead", "x: y", ":: ::", "\n", "tab\t", ":", " : "]
    choices = [[], [("a", "A")], [("a", "A"), ("b b", "B")], [("", "Empty")], [(" ", "Space")],
               [(":", "colon")]]
    out = []
    for label in labels:
        for tip in tips:
            for ch in choices:
                opt = types.SimpleNamespace(label=label, tooltip=tip, choices=ch)
                out.append({"label": label, "tooltip": tip, "choices": [list(c) for c in ch],
                            "text": PF._option_text(opt)})
    for tip in tips:  # a repeated section: no label, no choices
        sec = Section("s", "Title", [], tooltip=tip, repeat=True)
        out.append({"section": True, "tooltip": tip, "text": PF._option_text(sec)})
    dump("option_text_cases.json", out)
    return len(out)


# ------------------------------------------------------------------------------------ registry

def message(fn):
    try:
        fn()
    except (KeyError, ValueError) as exc:
        return exc.args[0]
    return None


REGISTER = [
    ("", [Section("a", "A", [Option("k", "K", "text")])]),
    ("amber-md", [Section("a", "A", [Option("k", "K", "text")])]),
    ("parity-dupes", [Section("a", "A", [Option("k", "K", "text"), Option("b", "B", "text")]),
                      Section("b", "B", [Option("k", "K", "text"), Option("b", "B2", "text")]),
                      Section("z", "Z", [Option("k", "K", "text")], repeat=True)]),
    ("parity-dupe-section", [Section("a", "A", [Option("stages", "S", "text")]),
                             Section("stages", "S", [Option("x", "X", "text")], repeat=True),
                             Section("rep", "R", [Option("x", "X", "text")], repeat=True),
                             Section("rep", "R", [Option("y", "Y", "text")], repeat=True)]),
    ("parity-ok", [Section("a", "A", [Option("x", "X", "text")]),
                   Section("r", "R", [Option("x", "X", "text")], repeat=True)]),
    ("parity-ok", [Section("a", "A", [Option("x", "X", "text")])]),
]


def registry_cases():
    names = ["x", "", "it's", "Ünï", "amber", "Amber-build", " amber-build"]
    get_tool = [{"name": n, "error": message(lambda n=n: registry.get_tool(n))} for n in names]
    by_category = [[cat, [t.name for t in tools]] for cat, tools in registry.tools_by_category()]
    templates = []
    for t in BUILTIN:
        for s in t.sections:
            if s.repeat:
                for ident in ["nope", ""] + [x[0] for x in s.templates]:
                    err = message(lambda s=s, ident=ident: s.template(ident))
                    templates.append({"tool": t.name, "section": s.key, "id": ident, "error": err,
                                      "record": None if err else enc(s.template(ident))})
    register = []
    for name, sections in REGISTER:
        cls = type("Synthetic", (Tool,), {"name": name, "label": name, "category": "Parity", "tag": "x",
                                          "sections": sections})
        register.append({"type": dump_tool(cls()), "error": message(lambda cls=cls: registry.register(cls))})
    paths = ["a.pdb", "A.PDB", "x.fa", "x.FASTA", "x.dbn", "x.ss", "x.map", "x.ccp4", "x.cif", "x.smi",
             "x.out", "x.prmtop", "x.parm7", "x.rst7", "x.ncrst", "x.inpcrd", "x.crd", "x.nc", "x.mdcrd",
             "x.RST", "x.rst", "x.dist", "x.top", "x.gro", "x.xtc", "x.trr", "x.trafl", "x.png", "x.tga",
             "x.mp4", "x.txt", "x.dat", "x", "", ".pdb", "pdb", "x.pdb.gz", "x.pdbx", "dir.pdb/x",
             "x.K", "É.PDB", "x.Mrc", "a.b.c.ccp4", "x.secstruct", "x.mdcrd.nc"]
    kinds = list(schema.FILE_KINDS) + ["", "PDB", "unknown"]
    dump("registry_cases.json", {
        "get_tool": get_tool, "tools_by_category": by_category, "templates": templates,
        "register": register,
        "guess_kind": [[p, schema.guess_kind(p)] for p in paths],
        "kind_label": [[k, schema.kind_label(k)] for k in kinds],
    })
    return len(get_tool) + len(templates) + len(register) + len(paths) + len(kinds)


def main():
    counts = {
        "coerce_cases.json": coerce_cases(),
        "digit_cases.json": digit_cases(),
        "normalise_cases.json": normalise_cases(),
        "requirements_cases.json": requirements_cases(),
        "option_text_cases.json": option_text_cases(),
        "registry_cases.json": registry_cases(),
    }
    counts["setting_comments.json"] = setting_comments(search_tools())
    for name, n in counts.items():
        print(f"{name}: {n} cases")


main()
