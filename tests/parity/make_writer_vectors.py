#!/usr/bin/env python3
"""Test vectors for the C++ port of mdatk/tomlio.py (and the textwrap / round() behaviour it relies on).

    python3 -I -B make_writer_vectors.py MDATK_CHECKOUT OUT_DIR

Writes into OUT_DIR (JSON, UTF-8, ensure_ascii=False):

* ``textwrap_vectors.json``  [{"text", "width", "lines"}]: ``textwrap.wrap(text, width)`` with default options.
* ``comment_vectors.json``   [{"text", "indent", "lines"}]: the lines ``tomlio.Writer.comment(text, indent)`` adds.
* ``value_vectors.json``     [{"input", "output"}] for ``tomlio.value``; [{"input", "error"}] for TypeErrors.
* ``key_vectors.json``       [{"key", "output"}] for ``tomlio.key`` (keys given as UTF-8 hex too).
* ``dumps_vectors.json``     [{"name", "data", "docs", "header", "output"}]: ``tomlio.dumps(data, doc, header)``
                             where doc(path, stack) = docs.get("/".join(path)) (missing: None).
* ``round2_vectors.json``    [{"hex", "repr_in", "round2"}]: ``repr(round(x, 2))`` (LogEntry.duration_s).
* ``entry_vectors.json``     [{"entry", "with_doc", "output"}]: ``projectfile._entry_text`` (one [[log]] entry).
                             A long one-character string is given as {"t": "str", "repeat": c, "count": n}; long
                             outputs as output_head / output_tail / output_len_codepoints.
* ``unicode_classes.json``   code-point ranges of Python 3.13's ``re`` ``\\w`` and ``\\d`` (str patterns),
                             ``str.isspace`` and the line boundaries of ``str.splitlines``.

Tagged values (``data`` / ``input``): {"t": "none" | "bool" | "int" | "float" | "str" | "datetime" | "date" |
"time" | "array" | "table", ...}.  int: "v" is a decimal string.  float: "hex" is ``float.hex()`` (exact; "nan",
"inf", "-inf" for those), "repr" is Python's repr.  str: "utf8" is the hex of the bytes a C++ std::string holds
(UTF-8; a lone surrogate U+DCxx from an undecodable file name is the single byte xx, as os.fsencode gives it),
"v" the text when it is valid Unicode.  table: "v" is a list of [key, value] pairs in order.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import random
import re
import sys
import textwrap

sys.dont_write_bytecode = True
if len(sys.argv) != 3:
    sys.exit("usage: make_writer_vectors.py MDATK_CHECKOUT OUT_DIR")
MDATK = sys.argv[1]
sys.path.insert(0, MDATK)

from mdatk import projectfile as PF  # noqa: E402
from mdatk import tomlio  # noqa: E402

HERE = os.path.abspath(sys.argv[2])
os.makedirs(HERE, exist_ok=True)


def dump(name, obj):
    with open(os.path.join(HERE, name), "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, ensure_ascii=False)
        fh.write("\n")


def has_surrogate(s):
    return any(0xD800 <= ord(c) <= 0xDFFF for c in s)


def enc(v):
    if v is None:
        return {"t": "none"}
    if isinstance(v, bool):
        return {"t": "bool", "v": v}
    if isinstance(v, int):
        return {"t": "int", "v": str(v)}
    if isinstance(v, float):
        return {"t": "float", "hex": v.hex(), "repr": repr(v)}
    if isinstance(v, str):
        out = {"t": "str", "utf8": v.encode("utf-8", "surrogateescape").hex()}
        if not has_surrogate(v):
            out["v"] = v
        return out
    if isinstance(v, dt.datetime):
        return {"t": "datetime", "iso": v.isoformat()}
    if isinstance(v, dt.date):
        return {"t": "date", "iso": v.isoformat()}
    if isinstance(v, dt.time):
        return {"t": "time", "iso": v.isoformat()}
    if isinstance(v, (list, tuple)):
        return {"t": "array", "v": [enc(x) for x in v]}
    if isinstance(v, dict):
        return {"t": "table", "v": [[k, enc(x)] for k, x in v.items()]}
    raise TypeError(type(v))


def text_out(s):
    """Output text: plain JSON string plus its UTF-8 hex when it holds characters JSON cannot show plainly."""
    return s


# ---------------------------------------------------------------------------- textwrap
CRAFTED = [
    "", " ", "   ", "a", "hello world", "  leading spaces kept on the first line", "trailing spaces   ",
    "two  spaces inside", "tab\there", "\ttab first", "a\tb\tc", "x" * 120, "x" * 98, "x" * 99,
    "word " * 30, "a-b", "ab-cd", "ab-cd-ef", "a-b-c", "a-b-c-d", "abc-def-ghi", "well-known", "e-mail",
    "co-op", "non-trivial-thing", "x-1", "1-x", "ab-1", "ab--cd", "ab---cd", "ab -- cd", "word--word",
    "a.--b", "--option", "-b option", "pre-  post", "Mg2+-bridged", "5′-terminal", "été-hiver",
    "Ångström-scale", "٣٤-ab", "a_b-c_d", "x y", " ", "a   b",
    "Look, goof-ball -- use the -b option!", "Hello there -- you goof-ball, use the -b option!",
    "/home/user/a-very-long-path/with-hyphens-and-more-hyphens/and/even/more/segments/than/fit/on/one/line.pdb",
    "https://example.org/some/really/long/url/that/goes/on/and/on/and/on/and/on/and/on/and/on/and/on/x.html",
    "-" * 120, "a" + "-" * 120, "ab-" * 50, "a-" * 60, "ab-" * 33 + "c",
    "Sequence / secondary structure: residues shown, e.g. 'A:1-20, B:5-30' (empty: all residues of the "
    "structure view).",
    "RNA force field: leaprc file for RNA. 'OL3 + CP' adds the phosphate oxygen vdW correction of Steinbrecher "
    "et al. (JCTC 2012) from the frcmod below; NBfix19 / 0BPh assume it. Choices: ol3cp, ol3, ljbb, roc, yil, shaw.",
    "x" * 97 + " -y", "x" * 96 + "-yy", "x" * 96 + "--yy", "y " + "x" * 200, "yy-" + "x" * 200,
    "a" * 50 + "-" + "b" * 60, "a" * 97 + "-" + "b" * 10, "ab-cd " * 20, "é" * 120,
    " " * 120 + "word", "word" + " " * 120, " " * 99 + "x", "x " * 49 + "y", "x " * 49 + "yz",
]


def textwrap_vectors():
    out = [{"text": t, "width": 98, "lines": textwrap.wrap(t, 98)} for t in CRAFTED]
    small = ["ab-cd-ef gh", "abcdef", "a-b-c-d-e", "aa-bb-cc-dd", "--x", "x--y z", "ab---cd", "a  b",
             "  ab cd", "ab\tcd", "abc-defg", "éé-éé", "12-34", "ab-12", "a_-_b"]
    for t in small:
        for w in range(1, 9):
            out.append({"text": t, "width": w, "lines": textwrap.wrap(t, w)})
    rng = random.Random(20261015)
    alphabet = list("abcxyzAB019-- -  \t.,!?'\"&_()/:;") + ["é", "Å", "٣", " ", "′",
                                                           "²", "中", "\U0001F600"]
    for _ in range(1500):
        n = rng.randint(0, 60)
        t = "".join(rng.choice(alphabet) for _ in range(n))
        w = rng.choice([1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 98])
        out.append({"text": t, "width": w, "lines": textwrap.wrap(t, w)})
    for _ in range(300):
        words = []
        for _ in range(rng.randint(1, 40)):
            words.append("".join(rng.choice("abcdefghij-") for _ in range(rng.randint(1, 14))))
        t = rng.choice([" ", "  ", " \t"]).join(words)
        out.append({"text": t, "width": 98, "lines": textwrap.wrap(t, 98)})
    dump("textwrap_vectors.json", out)
    return len(out)


# ---------------------------------------------------------------------------- comments
def comment_lines(text, indent=""):
    w = tomlio.Writer()
    w.comment(text, indent)
    return w.lines


def comment_vectors():
    texts = ["", None, "One line.", "x" * 150, "para one\npara two", "para one\n\npara three", "\n", "\n\n",
             "trailing newline\n", "\nleading newline", "crlf\r\nline", "cr\rline", "vt\x0bline", "ff\x0cline",
             "fs\x1cline", "nel\x85line", "ls line", "ps line", "tab\tinside", "  indented",
             "trailing nbsp ", "trailing spaces   ", "Ends with colon:", "a" * 97 + " b",
             "a" * 98 + " b", "a" * 99, "Contact map: 'instant' (evolves with the movie), 'frame' (one fixed frame) "
             "or 'frequency' (fraction of frames in contact)."]
    out = []
    for t in texts:
        for indent in ("", "    "):
            out.append({"text": t, "indent": indent, "lines": comment_lines(t, indent)})
    out.append({"text": 0, "indent": "", "lines": comment_lines(0), "note": "falsy non-string: nothing"})
    dump("comment_vectors.json", out)
    return len(out)


# ---------------------------------------------------------------------------- values
def value_vectors():
    ins = [
        "", "plain", 'a "b" \\ c\ttab', "line 1\nline \"2\"\n", "a\n", "\n", "\r", "a\r\nb",
        "\x00\x01\x08\x0b\x0c\x1f\x7f\x80\x9f", "é Å ≤ \U0001F600", "caf\udce9.pdb", "\udc80x\udcff",
        '"""', 'say """hi"""\nnext', "\\n literal", " ", "tab\tand\nnewline", "\x1b[0m",
        0, -3, 42, 2 ** 63 - 1, -(2 ** 63), True, False,
        0.0, -0.0, 1.0, -1.5, 0.1, 2.5, 1 / 3, 1e-7, 1e-5, 0.0001, 1e15, 1e16, 123456789012345680.0, 1e22,
        5e-324, 1.7976931348623157e308, float("nan"), float("inf"), float("-inf"), 12.345678, 300.0,
        dt.datetime(2026, 10, 15, 14, 32, 5), dt.datetime(2026, 10, 15, 14, 32, 5, 123456),
        dt.datetime(2026, 10, 15, 14, 32, 5, tzinfo=dt.timezone(dt.timedelta(hours=2))),
        dt.datetime(2026, 10, 15, 14, 32, 5, 500, tzinfo=dt.timezone.utc), dt.date(2026, 10, 15),
        dt.time(14, 32, 5), dt.time(14, 32, 5, 7),
        [], [1, 2.5, "x"], list(range(60)), [None], [1, None, 2], [[1, 2], [3]], [list(range(40))],
        ["a" * 92], ["a" * 93], ["é" * 92], ["é" * 93], [0.1] * 31, [0.1] * 32,
        ["multi\nline", "x"], [{"a": 1}, 2], [{"a": None}, 2], [{}, 1], [{"a": 1}], (1, 2),
        [True, False], [dt.date(2026, 1, 2)], [["x" * 50, "y" * 50]],
        {}, {"a": 1}, {"a": None}, {"weird key": 1, "b": [1, 2]}, {"a": {"b": 1}}, {"": 0},
    ]
    out = [{"input": enc(v), "output": tomlio.value(v)} for v in ins]
    errors = []
    for v in (object(), b"bytes", {1, 2}):
        try:
            tomlio.value(v)
        except TypeError as exc:
            errors.append({"python_type": type(v).__name__, "error": str(exc)})
    dump("value_vectors.json", {"values": out, "type_errors": errors})
    return len(out)


def key_vectors():
    keys = ["a", "A-b_9", "", "a.b", "a b", "é", "abc\n", "a\nb", 'q"', "1", "-", "_", "tab\t",
            "caf\udce9", "x" * 120, "with/slash", "plus+"]
    out = [{"key": enc(k), "output": tomlio.key(k)} for k in keys]
    dump("key_vectors.json", out)
    return len(out)


# ---------------------------------------------------------------------------- documents
def dumps_vectors():
    cases = []

    def add(name, data, docs=None, header=()):
        docs = docs or {}
        out = tomlio.dumps(data, lambda p, s: docs.get("/".join(p)), list(header))
        cases.append({"name": name, "data": enc(data), "docs": docs, "header": list(header), "output": out})

    add("empty", {})
    add("empty with header", {}, header=["Title line."])
    add("header lines", {"a": 1}, header=["One.", "", "Trailing spaces   ", "x" * 150])
    add("order: scalars, None, empty tables, tables, arrays of tables",
        {"t": {"x": 1}, "a": {}, "arr": [{"y": 1}, {"y": None}], "b": 1, "e": None, "f": [], "z": "last"},
        {"t": "doc t", "a": "doc a", "arr": "doc arr", "arr/y": "doc y", "b": "doc b", "e": "doc e", "f": "doc f",
         "t/x": "doc x"})
    add("no docs at all", {"a": 1, "t": {"b": {"c": [1, 2]}}, "arr": [{"k": "v"}, {"k": "w"}]})
    add("empty doc string means no comment", {"a": 1, "t": {"b": 2}}, {"a": "", "t": ""})
    add("multi-line and long docs", {"a": 1}, {"a": "first paragraph\nsecond paragraph that is long " + "w " * 60})
    add("quoted keys in headers", {"a b": {"c.d": {"k": 1}}, "é": {"x": None}})
    add("array of empty tables", {"arr": [{}, {}], "after": {"x": 1}})
    add("nested arrays of tables", {"n": [{"id": 1, "s": {"o": 1, "rec": [{"r": 1}, {"r": 2}]}},
                                          {"id": 2, "s": {"rec": [{"r": 3}]}}]},
        {"n": "nodes", "n/s": "settings", "n/s/rec": "records", "n/id": "id"})
    add("mixed list stays inline", {"m": [{"a": 1}, 2, "x"]})
    add("long arrays", {"l": list(range(60)), "s": ["x" * 40, "y" * 40, "z" * 40]})
    add("multi-line strings", {"m": "a\nb", "n": "ends\n", "o": "only\r"})
    add("table that only holds tables", {"outer": {"inner": {"k": 1}}}, {"outer": "o", "outer/inner": "i"})
    add("datetime", {"saved": dt.datetime(2026, 10, 15, 14, 32, 5)})
    dump("dumps_vectors.json", cases)
    return len(cases)


# ---------------------------------------------------------------------------- round(x, 2)
def round2_vectors():
    xs = [0.0, 0.001, 0.004, 0.005, 0.0050000001, 0.015, 0.025, 0.035, 0.045, 0.125, 0.375, 0.625, 0.875,
          1.005, 1.015, 2.675, 1.23456, 2.555000000000007, 0.13500000000000512, 7.0, 12345.678, 1e-9,
          1e17, 2.0 ** 53, 0.995, 9.995, 99.995, 1234567.125, 0.1 + 0.2, 3600.0, 86399.999]
    rng = random.Random(7)
    xs += [rng.uniform(0, 10) for _ in range(200)] + [rng.uniform(0, 100000) for _ in range(100)]
    xs += [k / 1000 for k in range(0, 2001, 5)]
    out = [{"hex": x.hex(), "repr_in": repr(x), "round2": repr(round(float(x), 2))} for x in xs]
    dump("round2_vectors.json", out)
    return len(out)


# ---------------------------------------------------------------------------- [[log]] entries
def entry_vectors():
    t = dt.datetime(2026, 10, 15, 14, 32, 5)
    entries = [
        PF.LogEntry(step="02", title="tleap build", time=t, duration_s=1.23456, command="bash x.sh",
                    folder="/w", files=["a.parm7", "b.rst7"], summary="exit status 0", citations=["k1"],
                    output="hello\nworld"),
        PF.LogEntry(step="03", time=t),
        PF.LogEntry(step="", title="", status="", time=t, duration_s=0.0),
        PF.LogEntry(step="04", time=t, output="x" * (PF.LOG_OUTPUT_CHARS + 5)),
        PF.LogEntry(step="05", time=t, output="é" * PF.LOG_OUTPUT_CHARS),
        PF.LogEntry(step="06", time=t, output="é" * (PF.LOG_OUTPUT_CHARS + 1)),
        PF.LogEntry(step="07", time=t, files=["caf\udce9.pdb"], summary="multi\nline\n"),
    ]
    out = []
    for e in entries:
        for with_doc in (True, False):
            fields = {}
            for k, v in vars(e).items():
                if isinstance(v, str) and len(v) > 1000 and len(set(v)) == 1:  # "x" * n: say so instead
                    fields[k] = {"t": "str", "repeat": v[0], "count": len(v)}
                else:
                    fields[k] = enc(v)
            text = PF._entry_text(e, with_doc)
            item = {"entry": fields, "with_doc": with_doc}
            if len(text) > 4000:  # keep the file small: give the long output's head and tail
                item["output_head"] = text[:600]
                item["output_tail"] = text[-200:]
                item["output_len_codepoints"] = len(text)
            else:
                item["output"] = text
            out.append(item)
    dump("entry_vectors.json", out)
    return len(out)


# ---------------------------------------------------------------------------- Unicode classes
def ranges(pred):
    out = []
    start = None
    for cp in range(0x110000):
        ok = (not 0xD800 <= cp <= 0xDFFF) and pred(chr(cp))
        if ok and start is None:
            start = cp
        elif not ok and start is not None:
            out.append([start, cp - 1])
            start = None
    if start is not None:
        out.append([start, 0x10FFFF])
    return out


def unicode_classes():
    w = re.compile(r"\w")
    d = re.compile(r"\d")
    lines = [c for c in map(chr, range(0x110000)) if not 0xD800 <= ord(c) <= 0xDFFF
             and len(("a" + c + "b").splitlines()) == 2]
    data = {
        "python": sys.version.split()[0],
        "unicodedata": __import__("unicodedata").unidata_version,
        "re_word": ranges(lambda c: w.match(c) is not None),
        "re_digit": ranges(lambda c: d.match(c) is not None),
        "str_isspace": ranges(str.isspace),
        "splitlines_boundaries": [ord(c) for c in lines] + ["\\r\\n counts as one boundary"],
        "textwrap_whitespace": [ord(c) for c in "\t\n\x0b\x0c\r "],
        "textwrap_word_punct": [ord(c) for c in "!\"'&.,?"] + ["plus every \\w character"],
    }
    dump("unicode_classes.json", data)
    return {k: len(v) for k, v in data.items() if isinstance(v, list)}


def main():
    print("textwrap", textwrap_vectors())
    print("comment", comment_vectors())
    print("value", value_vectors())
    print("key", key_vectors())
    print("dumps", dumps_vectors())
    print("round2", round2_vectors())
    print("entries", entry_vectors())
    print("unicode", unicode_classes())


if __name__ == "__main__":
    main()
