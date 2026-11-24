"""What MDATK's mmCIF code (mdatk/dictionaries.py) gives, for tests/io/test_cif.cpp.

    python3 -I -B make_cif_cases.py MDATK_CHECKOUT OUT_DIR

Imports MDATK (a checkout at v0.0.4.1, with pytest installed: the synthetic acedrg / eLBOW
dictionaries of its tests/test_dictionaries.py are used as fixtures; nothing is written into the
checkout) and writes OUT_DIR/cif_cases.json:

  parse        texts: MDATK's test fixtures, the synthetic dictionaries and what finish_dictionary
               makes of them, and crafted texts (quoting rules, text fields, loops, pairs, several
               blocks, comments, keywords in any case, line breaks and Unicode whitespace, every
               error). For each: list(_cif_tokens) or its error; parse_cif or its error; write_cif
               of the result and parse_cif of that; block / get / col / column lookups
  write        write_cif of documents built directly (every quoting branch, ragged and multi-row
               pair items, line breaks in cells, empty items, non-ASCII widths)
  cif_value    cif_value of many values
  lower        str.lower() of strings around final sigma and other special cases, and the SHA-256
               of lower() over every code point alone and in sigma contexts (core::pythonLower)
  dictionary   is_dictionary_file of files given as byte parts (hex, repeat count) and of a
               folder and a missing file
  find         find_dictionary scenarios: files (with ages in seconds before a fixed time),
               folders, and calls with sources, code and exclude, run in the scenario's root
               folder; paths in results and messages have the root written {root}

Texts that are not valid UTF-8 are given as "hex" (MDATK decodes with errors="replace"). The output
is deterministic.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import random
import shutil
import sys
import tempfile

if len(sys.argv) != 3:
    sys.exit("usage: make_cif_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
sys.path.insert(0, CHECKOUT)

from mdatk import dictionaries as D  # noqa: E402

RNG = random.Random(20261124)


def load_tests():
    spec = importlib.util.spec_from_file_location(
        "mdatk_test_dictionaries", os.path.join(CHECKOUT, "tests", "test_dictionaries.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def plain(doc):
    return {"header": list(doc.header),
            "blocks": [{"name": b.name, "items": [{"tags": list(i.tags), "rows": [list(r) for r in i.rows],
                                                   "pairs": i.pairs} for i in b.items]} for b in doc.blocks]}


def from_plain(data):
    doc = D.CifDocument(list(data["header"]))
    for b in data["blocks"]:
        doc.blocks.append(D.CifBlock(b["name"], [D.CifLoop(i["tags"], [list(r) for r in i["rows"]], i["pairs"])
                                                 for i in b["items"]]))
    return doc


def attempt(fn):
    try:
        return {"ok": fn()}
    except D.DictionaryError as e:
        return {"error": str(e)}


MDATK_FIXTURE = """# header comment
data_comp_list
_chem_comp.id   XYZ
_chem_comp.name 'name with "quotes" inside'
_chem_comp.group non-polymer
data_comp_XYZ
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
XYZ "O5'" O
XYZ 'C1' C # trailing comment
XYZ N1'x N
loop_
_x.text
;
line one
line two
;
"""

CRAFTED = [
    "",
    "\n\n\n",
    "# only a comment\n#\n  # indented\n",
    "data_a\n",
    "data_\n_a.b 1\n",
    "DATA_Upper\n_A.B 1\n_a.c 2\n_B.x 3\n",
    "Data_x _a.b 1 _a.c 2 _a.b 3\n",
    "data_a\n_a.b 'it''s'\n_a.c \"say \"hi\" now\"\n_a.d 'x'y' z\n",
    "data_a\n_a.b a'b\n_a.c a\"b\"\n_a.d ''\n_a.e \"\"\n_a.f '' _a.g \"\"x\"\n",
    "data_a\n_a.b\t'tab'\t_a.c 'nbsp' _a.d　'ideo'　\n",
    "data_a\n_a.b 'unterminated\n",
    "data_a\n_a.b 'end'x\n",
    "data_a\n_a.b \"q\"\"\n",
    "data_a _a.b 1 # comment _a.c 2\n_a.d a#b\n#_a.e 3\n",
    "data_a\nloop_\n_a.b\n_a.c\n1 2 3\n",
    "data_a\nloop_\n_a.b\n_a.c\n1 2\n3 4\n_b.x 5\n",
    "data_a\nloop_\n1 2\n",
    "data_a\nloop_\n",
    "data_a\nloop_\ndata_b\n",
    "data_a\nloop_\nloop_\n_a.b\n",
    "data_a\nloop_\n_a.x\n_b.y\n1 2\n",
    "data_a\nloop_\n_a.x\n1\n2\nloop_\n_b.y\n3\n",
    "data_a\nloop_\n_a.x\n1\n_a.y\n2\n",
    "data_a\nloop_\n_a.x\n_a.y\n1 2 _a.z 3\n",
    "data_a\nLOOP_\n_a.x\n1\nLoop_x 2\n",
    "data_a\n_a.b\n_a.c 1\n",
    "data_a\n_a.b\nloop_\n_a.c\n1\n",
    "data_a\n_a.b\n",
    "data_a\n_a.b\ndata_b\n1\n",
    "data_a\n_a.b\ndata_b\n_a.c 1\n",
    "data_a\n1\n",
    "data_a\n\"it's\"\n",
    "data_a\n;\ntext\n;\n",
    "data_a\n_a.b\n;line one\n  line two\n;\n_a.c 2\n",
    "data_a\n_a.b\n;\n\nfirst blank kept\n;\n",
    "data_a\n_a.b\n;   \nwhitespace first line dropped\n;  trailing words ignored\n_a.c 3\n",
    "data_a\n_a.b\n; x\n ;not the end\n;\n",
    "data_a\n_a.b\n;\nnever closed\n",
    "data_a\n_a.b\n;\n;\n",
    "data_a\n_a.b\n;only line\n;",
    "data_a\r\n_a.b\r\n;\r\ncrlf text\r\n;\r\n_a.c 'x'\r\n",
    "data_a\r_a.b 1\r_a.c 2",
    "data_a\n_a.b 1\x0b_a.c 2\x0c_a.d 3\x1c_a.e 4\x1d_a.f 5\x1e_a.g 6\x85_a.h 7 _a.i 8 _a.j 9\n",
    "data_a\n_a.b 'x\x0by'\n",
    "data_a\n_a.b v _a.c w\n",
    "_x.y 1\nloop_\n_p.q\n1\nglobal_\n_lib.name mon\ndata_a\n_a.b 1\n",
    "# h1\n# h2\nvalue_before\ndata_a\n# dropped\n_a.b 1\n",
    "global_ save_x stop_ GLOBAL_ Save_ STOP_\ndata_a\nsave_frame\n_a.b 1\nsave_\nstop_\n",
    "data_a\n_a.b global_x\n_a.c stop_x\n_a.d save_\n",
    "data_a\n_a.b 1\ndata_a\n_a.b 2\ndata_A\n_a.b 3\n",
    "data_a\n_a.b 1\n_b.c 2\n_a.d 3\n_a.b 4\n",
    "data_a\n_foo 1\n_foo 2\n_bar 3\n",
    "data_a\n_A.X 1\n_a.x 2\n_a.Y 3\n",
    "data_a\n_Äb.c 1\n_äB.d 2\n_ΣΑ.x 1\n_ΑΣ.y 2\n_İ.z 3\n_Ka.x 4\n_ka.y 5\n",
    "data_été\n_té.ü 'café ü'\nloop_\n_m.n\nÅngström\n☃\n",
    "data_a\n_a.b \"O5'\"\n_a.c 'O5\"'\n_a.d \"a' b\"\n",
    "data_a\nloop_\n_a.x\n_a.y\n;\nmulti\nline\n;\n2\n3\n;\nend\n;\n",
    "data_a\n   loop_\n  _a.x\n\t1\n",
    "data_a\n_a.b 1 _a.b 2\n",
    "data_a\n_a.b x y\n",
    "data_a\nloop_\n_a.b\n_a.c\n",
]

BYTE_TEXTS = [
    b"data_a\n_a.b caf\xe9\n_a.c '\xff\xfe'\n",
    b"data_\xe2\x82\n_a.\xc3 1\n",
]


def parse_case(name, text):
    case = {"name": name, "text": text}
    try:
        case["tokens"] = [[k, v] for k, v in D._cif_tokens(text)]
    except D.DictionaryError as e:
        case["tokens_error"] = str(e)
    try:
        doc = D.parse_cif(text)
    except D.DictionaryError as e:
        case["error"] = str(e)
        return case
    case["doc"] = plain(doc)
    written = D.write_cif(doc)
    case["written"] = written
    case["reparsed"] = attempt(lambda: plain(D.parse_cif(written)))
    lookups = []
    for b in doc.blocks + [D.CifBlock("nope")]:
        for q in (b.name, b.name.upper(), b.name.swapcase()):
            found = doc.block(q)
            lookups.append({"block": q, "index": next((k for k, x in enumerate(doc.blocks) if x is found), -1)})
        real = next((k for k, x in enumerate(doc.blocks) if x is b), -1)
        if real < 0:
            continue
        for item in b.items + [D.CifLoop(["_nope.x"])]:
            for q in (item.category, item.category.upper(), item.category.swapcase()):
                found = b.get(q)
                lookups.append({"in": real, "get": q,
                                "index": next((k for k, x in enumerate(b.items) if x is found), -1)})
            if item not in b.items:
                continue
            ii = b.items.index(item)
            for tag in item.tags + ["_x.nope", "nope"]:
                for q in (tag.split(".", 1)[-1], tag.split(".", 1)[-1].upper(), tag):
                    entry = {"in": real, "item": ii, "col": q, "index": item.col(q)}
                    try:
                        entry["column"] = item.column(q)
                    except IndexError:
                        entry["column_error"] = True
                    lookups.append(entry)
    case["lookups"] = lookups
    return case


def parse_cases(tests):
    out = [parse_case("mdatk-fixture", MDATK_FIXTURE)]
    out.append(parse_case("mdatk-values-for", "data_a\nloop_\n_a.b\n_a.c\n1 2 3\n"))
    out.append(parse_case("mdatk-unterminated", "data_a\n_a.b 'open\n"))
    for code in ("0A5", "0C3", "0G", "0U5"):
        for style in ("acedrg", "elbow"):
            text, _new, _values = tests.synthetic_dictionary(D.COMPONENTS[code], style, seed=len(code) + len(style))
            out.append(parse_case(f"synthetic-{code}-{style}", text))
            out.append(parse_case(f"finished-{code}-{style}", D.finish_dictionary(text, code).text))
    for k, text in enumerate(CRAFTED):
        out.append(parse_case(f"crafted-{k}", text))
    for k, raw in enumerate(BYTE_TEXTS):
        case = parse_case(f"bytes-{k}", raw.decode("utf-8", "replace"))
        del case["text"]
        case["hex"] = raw.hex()
        out.append(case)
    # random token soup
    pool = ["data_a", "data_b", "loop_", "_a.x", "_a.y", "_b.z", "_A.X", "1", "'q q'", "\"d'q\"", "#c", "\n;t\n;\n",
            "save_", "global_", "x", "\n", "\t", "''"]
    for k in range(60):
        words = [RNG.choice(pool) for _ in range(RNG.randint(1, 25))]
        text = "data_r\n" + " ".join(words) + "\n" if k % 2 else " ".join(words)
        out.append(parse_case(f"random-{k}", text))
    return out


VALUES = ["", "a", "A1", "O5'", "O5\"", "a'b\"c", "with space", "tab\there", "　", "x y", "line\nbreak",
          "\n", "cr\r", "_under", "#hash", "$dollar", "'quote", "\"dq", "[br", "]br", ";semi", "a;b", "a_b",
          "data_x", "DATA_X", "Loop_", "loop", "save_", "Global_x", "stop_", "stops", "data", ".", "?", "-",
          "1.000", "été", "☃", "KLOOP_", "İdata_", "trailing ", " leading", "a#b",
          "'", "\"", "''", "\"'", "x\x0by", "x\x85y", " ", "a b"]


def write_cases():
    docs = []
    docs.append({"header": [], "blocks": []})
    docs.append({"header": ["# only header"], "blocks": []})
    docs.append({"header": ["# h", "#"], "blocks": [{"name": "a", "items": []}, {"name": "b", "items": []}]})
    docs.append({"header": [], "blocks": [{"name": "values", "items": [
        {"tags": [f"_v.t{k}" for k in range(len(VALUES))], "rows": [VALUES], "pairs": True}]}]})
    docs.append({"header": ["#"], "blocks": [{"name": "loopvalues", "items": [
        {"tags": ["_l.k", "_l.v", "_l.w"], "rows": [[str(k), v, "x" * (k % 5)] for k, v in enumerate(VALUES)],
         "pairs": False}]}]})
    docs.append({"header": [], "blocks": [{"name": "odd", "items": [
        {"tags": [], "rows": [], "pairs": False},
        {"tags": [], "rows": [["1"]], "pairs": True},
        {"tags": ["_p.a", "_p.bb"], "rows": [["1", "2"], ["3", "4"]], "pairs": True},
        {"tags": ["_p.a", "_p.bb"], "rows": [], "pairs": True},
        {"tags": ["_p.a", "_p.bbb", "_p.c"], "rows": [["1"]], "pairs": True},
        {"tags": ["_p.a"], "rows": [["1", "2", "3"]], "pairs": True},
        {"tags": ["_l.a", "_l.b", "_l.c"], "rows": [["1", "2", "3"], ["44"], ["5", "6", "7", "8"]], "pairs": False},
        {"tags": ["_l.a", "_l.b"], "rows": [], "pairs": False},
        {"tags": ["_l.a", "_l.b"], "rows": [["x\ny", "z"], ["long value", "w"]], "pairs": False},
        {"tags": ["_l.a"], "rows": [["x\ny"], ["p\nq"]], "pairs": False},
        {"tags": ["_l.a", "_l.b"], "rows": [["", "it's \"q\""], ["éé", "☃"]], "pairs": False},
        {"tags": ["_éé.ü", "_e.long_tag_name"], "rows": [["1", "2"]], "pairs": True},
        {"tags": ["_t.a"], "rows": [["v  "]], "pairs": True},
        {"tags": ["_t.a"], "rows": [["multi\nline"]], "pairs": True},
    ]}, {"name": "", "items": [{"tags": ["_x.y"], "rows": [["1"]], "pairs": False}]}]})
    for k in range(25):
        items = []
        for _ in range(RNG.randint(0, 4)):
            ntags = RNG.randint(0, 4)
            tags = [f"_c{RNG.randint(0, 3)}.{RNG.choice(['a', 'bb', 'ccc', 'é'])}" for _ in range(ntags)]
            pairs = RNG.random() < 0.4
            nrows = 1 if pairs and RNG.random() < 0.8 else RNG.randint(0, 4)
            rows = [[RNG.choice(VALUES) for _ in range(ntags if RNG.random() < 0.85 else RNG.randint(0, 5))]
                    for _ in range(nrows)]
            items.append({"tags": tags, "rows": rows, "pairs": pairs})
        docs.append({"header": RNG.choice([[], ["#"], ["# x"], ["# x", "#"]]),
                     "blocks": [{"name": f"r{k}", "items": items}]})
    out = []
    for d in docs:
        written = D.write_cif(from_plain(d))
        out.append({"doc": d, "written": written, "reparsed": attempt(lambda: plain(D.parse_cif(written)))})
    return out


def lower_cases():
    strings = ["", "ABC", "Σ", "AΣ", "ΣA", "AΣA", "AΣ ", "AΣ.", "ÁΣ",
               "AΣ́", "AΣ́B", "A'Σ", "AΣ'B", "1Σ", "_ΣΑ", "_ΑΣ",
               "ΑΣΣ", "ΣΣΣ", "İ", "Iİi", "K", "ß", "ẞ",
               "İ", "ǅ", "Ⅰ", "Ⓐ", "\U00010400", "\U0001e900", "ͅΣ",
               "ʰΣ", "A­Σ", "Σ­", "AΣʰ"]
    pool = "AaZzΣσΆ'.­ 1_İKÅʰͅⅠ\U00010400"
    for _ in range(300):
        strings.append("".join(RNG.choice(pool) for _ in range(RNG.randint(1, 8))))
    chars = [chr(c) for c in range(0x110000) if not 0xD800 <= c <= 0xDFFF]
    single = hashlib.sha256("\n".join(c.lower() for c in chars).encode("utf-8")).hexdigest()
    context = hashlib.sha256("\n".join((c + "Σ").lower() + "|" + ("AΣ" + c).lower() + "|"
                                       + ("A" + c + "Σ").lower() for c in chars).encode("utf-8")).hexdigest()
    return {"strings": [[s, s.lower()] for s in strings], "single_sha256": single, "context_sha256": context}


def dictionary_cases(root):
    marker = b"_chem_comp_atom.x 1\n_chem_comp_bond.y 2\n"
    specs = [
        [[marker.hex(), 1]],
        [[b"_chem_comp_atom.x\n".hex(), 1]],
        [[b"_chem_comp_bond.x\n_chem_comp_atom".hex(), 1]],
        [[b"_CHEM_COMP_ATOM._CHEM_COMP_BOND.".hex(), 1]],
        [[b"".hex(), 1]],
        [[b"\xff\xfe".hex(), 3], [marker.hex(), 1]],
        [[b"x".hex(), 1_999_970], [marker.hex(), 1]],
        [[b"x".hex(), 1_999_960], [marker.hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"x".hex(), 1_999_970], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"x".hex(), 1_999_967], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"x".hex(), 1_999_966], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"\r\n".hex(), 1_999_966], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"\r\n".hex(), 1_999_967], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"\r".hex(), 1_999_966], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], ["\U0001f600".encode().hex(), 1_999_966],
         [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], ["\U0001f600".encode().hex(), 1_999_967],
         [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"\xe2\x82".hex(), 1_999_966], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"\xf0\x9f\x98".hex(), 1_999_967], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"\xc0".hex(), 1_999_966], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"\xc0".hex(), 1_999_967], [b"_chem_comp_bond.".hex(), 1]],
        [[b"_chem_comp_atom.\n".hex(), 1], [b"x".hex(), 1_999_965], [b"\r\n_chem_comp_bond.".hex(), 1]],
    ]
    out = []
    for k, parts in enumerate(specs):
        path = os.path.join(root, f"d{k}.cif")
        with open(path, "wb") as fh:
            for unit, count in parts:
                fh.write(bytes.fromhex(unit) * count)
        out.append({"parts": parts, "result": D.is_dictionary_file(path)})
        os.remove(path)
    os.mkdir(os.path.join(root, "folder.cif"))
    out.append({"folder": True, "result": D.is_dictionary_file(os.path.join(root, "folder.cif"))})
    out.append({"missing": True, "result": D.is_dictionary_file(os.path.join(root, "missing.cif"))})
    os.rmdir(os.path.join(root, "folder.cif"))
    return out


DICT = "data_x\n_chem_comp_atom.a 1\n_chem_comp_bond.b 2\n"
NOT_DICT = "data_x\n_a.b 1\n"

SCENARIOS = [
    {   # MDATK's test_find_dictionary
        "dirs": ["work_0A5/0A5_TMP"],
        "files": [["work_0A5/0A5_TMP/inter.cif", DICT, 1], ["work_0A5/notes.cif", NOT_DICT, 2],
                  ["work_0A5/0A5.pdb", "END\n", 3], ["work_0A5/0A5.cif", DICT, 4]],
        "calls": [
            {"sources": ["{root}/work_0A5"], "code": "0A5"},
            {"sources": ["{root}/work_0A5/*.cif"]},
            {"sources": ["{root}/work_0A5"], "exclude": ["{root}/work_0A5/0A5.cif"]},
            {"sources": ["{root}/missing"]},
            {"sources": ["work_0A5", "work_0A5/*/*.cif"]},
            {"sources": ["work_0A5/0A5_TMP/inter.cif", "work_0A5/0A5.cif"]},
            {"sources": ["work_0A5/*/*.cif", "work_0A5/notes.cif"]},
            {"sources": ["work_0A5/*/*.cif"], "exclude": ["./work_0A5/0A5_TMP/../0A5_TMP/inter.cif"]},
        ],
    },
    {   # code, tmp and age order; ties by path
        "dirs": ["out", "tmpdir", "TMP2"],
        "files": [["out/a.cif", DICT, 10], ["out/0g5_out.cif", DICT, 50], ["out/new.cif", DICT, 1],
                  ["out/tmp_new.cif", DICT, 0], ["tmpdir/x.cif", DICT, 0], ["TMP2/y.cif", DICT, 0],
                  ["out/same1.cif", DICT, 5], ["out/same0.cif", DICT, 5], ["out/0G5.CIF", DICT, 0]],
        "calls": [
            {"sources": ["out"]}, {"sources": ["out"], "code": "0G5"}, {"sources": ["out"], "code": "0g5_OUT"},
            {"sources": ["out"], "code": "zz"}, {"sources": ["tmpdir", "TMP2", "out/a.cif"]},
            {"sources": ["tmpdir", "TMP2"]}, {"sources": ["out/same1.cif", "out/same0.cif"]},
            {"sources": ["out/*.CIF"]}, {"sources": ["out/*.cif"], "code": "0G5"},
            {"sources": ["out/same1.cif", "./out/same1.cif"]},
            {"sources": ["out"], "exclude": ["out/new.cif", "out/tmp_new.cif", "{root}/out/0g5_out.cif"]},
        ],
    },
    {   # hidden names, patterns, folders named *.cif, odd names
        "dirs": ["h", "h/.hid", "h/sub1", "h/sub2", "h/dir.cif", "we[ir]d*", ".dot"],
        "files": [["h/.hidden.cif", DICT, 0], ["h/vis.cif", DICT, 9], ["h/.hid/in.cif", DICT, 0],
                  ["h/sub1/x1.cif", DICT, 3], ["h/sub2/x2.cif", DICT, 2], ["h/sub2/xy.cif", DICT, 1],
                  ["we[ir]d*/w.cif", DICT, 0], ["we[ir]d*/.w.cif", DICT, 0], [".dot/d.cif", DICT, 0],
                  ["h/a.cif", DICT, 7], ["h/b.cif", DICT, 6], ["h/-.cif", DICT, 8], ["h/é.cif", DICT, 8],
                  ["h/[x].cif", DICT, 5]],
        "calls": [
            {"sources": ["h"]}, {"sources": ["h/.*.cif"]}, {"sources": ["h/*/*.cif"]}, {"sources": ["h/sub?/x?.cif"]},
            {"sources": ["h/[ab].cif"]}, {"sources": ["h/[!a].cif"]}, {"sources": ["h/[!ab].cif"]},
            {"sources": ["h/[]-a].cif"]}, {"sources": ["h/[a-].cif"]}, {"sources": ["h/[z-a].cif"]},
            {"sources": ["h/[c-b-a].cif"]}, {"sources": ["h/[à-ÿ].cif"]}, {"sources": ["h/[[]x].cif"]},
            {"sources": ["h/[x].cif"]}, {"sources": ["h/[.cif"]}, {"sources": ["h/dir.cif"]},
            {"sources": ["h/*.cif/"]}, {"sources": ["we[ir]d*"]}, {"sources": ["we[ir]d*/w.cif"]},
            {"sources": ["we[[]ir]d[*]/*.cif"]}, {"sources": ["*/w.cif"]}, {"sources": [".*/*.cif"]},
            {"sources": ["*"]}, {"sources": ["h//sub1//x1.cif"]}, {"sources": ["h//sub1//*.cif"]},
            {"sources": ["{root}/h/sub*/*.cif"]}, {"sources": []}, {"sources": ["", "nothing*", "h/*.pdb"]},
            {"sources": ["h/sub1/"]}, {"sources": ["h/**/*.cif"]},
        ],
    },
]


def find_cases(root):
    out = []
    old_cwd = os.getcwd()
    base = 1_700_000_000.25
    for k, sc in enumerate(SCENARIOS):
        here = os.path.join(root, f"s{k}")
        os.mkdir(here)
        for d in sc["dirs"]:
            os.makedirs(os.path.join(here, d), exist_ok=True)
        for path, text, age in sc["files"]:
            full = os.path.join(here, path)
            with open(full, "w", encoding="utf-8") as fh:
                fh.write(text)
            os.utime(full, (base - age, base - age))
        os.chdir(here)
        calls = []
        try:
            for call in sc["calls"]:
                sources = [s.replace("{root}", here) for s in call["sources"]]
                exclude = [s.replace("{root}", here) for s in call.get("exclude", [])]
                entry = dict(call)
                try:
                    entry["result"] = D.find_dictionary(sources, call.get("code", ""), exclude).replace(here, "{root}")
                except D.DictionaryError as e:
                    entry["error"] = str(e).replace(here, "{root}")
                except IndexError as e:
                    entry["python_error"] = f"IndexError: {e}"
                calls.append(entry)
        finally:
            os.chdir(old_cwd)
        out.append({"dirs": sc["dirs"], "files": sc["files"], "calls": calls})
    return out


def main() -> int:
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    tests = load_tests()
    root = tempfile.mkdtemp(prefix="mdbridge-cif-")
    try:
        data = {
            "parse": parse_cases(tests),
            "write": write_cases(),
            "cif_value": [[v, D.cif_value(v)] for v in VALUES],
            "lower": lower_cases(),
            "dictionary": dictionary_cases(root),
            "find": find_cases(root),
        }
    finally:
        shutil.rmtree(root)
    path = os.path.join(out_dir, "cif_cases.json")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
