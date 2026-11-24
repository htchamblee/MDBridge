"""What MDATK's parm7 module and MDAnalysis' TOPParser make of parm7 files (tests/data/io/parm7).

    /opt/mdbridge-venv/bin/python -I -B make_parm7_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDAnalysis 2.11.0 and numpy; nothing is written into the checkout. Writes

  read_cases.json       mdatk.io.parm7.read_sections of crafted files (bytes as hex): formats, D
                        exponents, blank fields, CRLF / CR line ends, Unicode, bad values, repeated
                        flags, `wanted` subsets; and of tests/data/topology/demo_dry.parm7
  writer_cases.json     _format_values and the bytes write_parm7 writes (or its error and the bytes
                        written before it)
  lj_cases.json         LJTable.from_parm7 of demo_dry.parm7 and pair_coefficients of atom pairs;
                        the "missing parm7 sections" errors
  topology/*.parm7      variants of demo_dry.parm7 (line ends, RESIDUE_CHAINID sections, atomic
                        numbers, formats, comments, unknown sections) and of a small system, with
  topology/*.topology.json
                        their MDAnalysis topologies (the dump of tests/parity/dump_topologies.py)
  topology_errors.json  the files TOPParser rejects, with str(exception)

Floats are float.hex() strings ("nan", "-nan", "inf", "-inf" for the others). Errors are
{"type": class name, "message": str(exc)} with the file's path replaced by "{path}".
"""

import json
import math
import os
import random
import shutil
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_parm7_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
from MDAnalysis.topology.TOPParser import TOPParser  # noqa: E402

from mdatk.io import parm7  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO = os.path.join(HERE, "..", "data", "topology", "demo_dry.parm7")
sys.path.insert(0, HERE)
import dump_topologies  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
TMP = tempfile.mkdtemp(prefix="parm7-cases-")


def fhex(x):
    x = float(x)
    if math.isnan(x):
        return "-nan" if math.copysign(1.0, x) < 0 else "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x.hex()


def error(exc, path):
    return {"type": type(exc).__name__, "message": str(exc).replace(path, "{path}")}


def sections_json(sections):
    out = []
    for name, values in sections.items():
        if isinstance(values, list):
            out.append([name, "a", values])
        elif values.dtype == np.int64:
            out.append([name, "i", [int(v) for v in values]])
        else:
            out.append([name, "f", [fhex(v) for v in values]])
    return out


# ---------------------------------------------------------------- read_sections

def P(*lines, nl="\n", end=True):
    text = nl.join(lines) + (nl if end else "")
    return text.encode("utf-8")


HEAD = "%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/26  00:00:00"
READ = [
    ("basic", P(HEAD, "%FLAG TITLE", "%FORMAT(20a4)", "demo".ljust(80), "%FLAG POINTERS", "%FORMAT(10I8)",
                "       3       2       0", "%FLAG ATOM_NAME", "%FORMAT(20a4)", "N   CA  C   ",
                "%FLAG CHARGE", "%FORMAT(5E16.8)", "  1.00000000E+00 -2.50000000E-01  3.00000000E+02")),
    ("d_exponents", P("%FLAG X", "%FORMAT(5E16.8)", "  1.00000000D+00 -2.50000000d-01  3.0000000D+002",
                      "%FLAG Y", "%FORMAT(3F10.4)", "    1.2500   -2.0000 1.5D0     ")),
    ("format_spacing", P("%FLAG A", "%FORMAT( 10 I 8 )", "       1       2", "%FLAG B", "%FORMAT  junk (3e12.4) x",
                         "  1.0000E+00  2.0000E+00", "%FLAG C", "%FORMAT(2a10)", "abcdefghijklmnopqrstuvwxyz")),
    ("format_dot_on_int", P("%FLAG A", "%FORMAT(5I6.3)", "     1     2     3")),
    ("format_lower_kinds", P("%FLAG A", "%FORMAT(10i8)", "      12", "%FLAG B", "%FORMAT(5f8.2)", "    1.25",
                             "%FLAG C", "%FORMAT(20A4)", "AB  CD  ")),
    ("format_bad", P("%FLAG A", "%FORMAT(10Q8)", "1")),
    ("format_bad_spaces", P("%FLAG A", "%FORMAT   (  )   ", "1")),
    ("format_empty", P("%FLAG A", "%FORMAT", "1")),
    ("format_before_flag", P("%FORMAT(garbage)", "%FORMAT(10I8)", "       5", "%FLAG A", "%FORMAT(10I8)", "       7")),
    ("format_twice", P("%FLAG A", "%FORMAT(10I8)", "       1       2", "%FORMAT(5E16.8)", "  1.00000000E+00")),
    ("format_after_data", P("%FLAG A", "%FORMAT(2a8)", "abcdefgh12345678", "%FORMAT(4a4)", "wxyz")),
    ("format_width_zero", P("%FLAG A", "%FORMAT(10I0)", "1")),
    ("format_width_zero_empty", P("%FLAG A", "%FORMAT(10I0)")),
    ("format_width_zero_unwanted", P("%FLAG A", "%FORMAT(10I0)", "1", "%FLAG B", "%FORMAT(10I8)", "       1"),
     ["B"]),
    ("format_unicode_digits", P("%FLAG A", "%FORMAT(١٠I٨)", "       1       2")),
    ("format_unicode_space", P("%FLAG A", "%FORMAT( 10I8　)", "       4")),
    ("format_second_paren", P("%FLAG A", "%FORMAT(x) (2I3)", "  1  2")),
    ("no_flag", P("hello", "%FORMAT(10I8)", "       1")),
    ("empty_file", b""),
    ("only_flag", P("%FLAG A")),
    ("flag_no_newline", P("%FLAG A", "%FORMAT(10I8)", "       1       2", end=False)),
    ("flag_names", P("%FLAG", "%FORMAT(10I8)", "       1", "%FLAGX", "       2", "%FLAG  spaced name  ",
                     "%FORMAT(10I8)", "       3", "%FLAG\tTAB", "%FORMAT(10I8)", "       4")),
    ("default_kind_a80", P("%FLAG A", "first line", "", "x" * 100)),
    ("version_comment", P("%VERSION x", "%FLAG A", "%COMMENT some note", "%FORMAT(10I8)", "%COMMENT 12",
                          "       1", "%VERSION       2", "       3")),
    ("other_percent", P("%FLAG A", "%FORMAT(20a4)", "%FOO", "%COMMENTARY")),
    ("data_before_flag", P("junk line", "       1", "%FLAG A", "%FORMAT(10I8)", "       9")),
    ("repeated_flag", P("%FLAG A", "%FORMAT(10I8)", "       1", "%FLAG B", "%FORMAT(10I8)", "       2",
                        "%FLAG A", "%FORMAT(5E16.8)", "  1.00000000E+00")),
    ("blank_fields_int", P("%FLAG A", "%FORMAT(10I8)", "       1               3        ", "",
                           "                       4")),
    ("blank_fields_a", P("%FLAG A", "%FORMAT(20a4)", "N       CA  ", "    ", "  B ")),
    ("blank_fields_float", P("%FLAG A", "%FORMAT(5E16.8)", "  1.00000000E+00                  2.0E+00")),
    ("touching_floats", P("%FLAG A", "%FORMAT(5E15.8)", "-7.61156190E+00 3.32181810E+00-7.61156190E+00")),
    ("long_lines", P("%FLAG A", "%FORMAT(3I4)", "   1   2   3   4   5", "%FLAG B", "%FORMAT(2a3)", "abcdefgh")),
    ("short_last_chunk", P("%FLAG A", "%FORMAT(10I8)", "       1      2", "%FLAG B", "%FORMAT(20a4)", "AB  C")),
    ("tabs", P("%FLAG A", "%FORMAT(10I8)", "\t1\t\t2", "%FLAG B", "%FORMAT(20a4)", "\tA\t B")),
    ("crlf", P(HEAD, "%FLAG A", "%FORMAT(10I8)", "       1       2", "%FLAG B", "%FORMAT(20a4)", "N   CA  ",
               nl="\r\n")),
    ("cr", P(HEAD, "%FLAG A", "%FORMAT(10I8)", "       1       2", "%FLAG B", "%FORMAT(20a4)", "N   CA  ",
             nl="\r")),
    ("mixed_newlines", b"%FLAG A\r\n%FORMAT(10I8)\r       1\n       2\r\n%FLAG B\r%FORMAT(20a4)\nO   \r"),
    ("other_line_breaks", P("%FLAG A", "%FORMAT(20a4)", "A\x0bB\x0cC\x1cD\x85E F G   H")),
    ("unicode_names", P("%FLAG A", "%FORMAT(20a4)", "étéCα  中文字符\U0001F600xyzw",
                        "%FLAG B", "%FORMAT(10I8)", "١٢".rjust(8) + "    １２")),
    ("invalid_utf8", b"%FLAG A\n%FORMAT(20a4)\nab\xffc\xe2\x82d   \xc3(  \xed\xa0\x80xx\xf0\x9f\x98\n"
                     b"\xc0\xaf\xe0\x80\x80\xf4\x90\x80\x80\xf8\x88\x80\x80\x80\xc2\n%FLAG \xfe\n%FORMAT(10I8)\n"),
    ("invalid_utf8_int", b"%FLAG B\n%FORMAT(10I8)\n\xff      1\n"),
    ("bom", b"\xef\xbb\xbf%FLAG A\n%FORMAT(10I8)\n       1\n%FLAG B\n%FORMAT(10I8)\n       2\n"),
    ("bad_int", P("%FLAG A", "%FORMAT(10I8)", "     12x")),
    ("bad_int_float", P("%FLAG A", "%FORMAT(10I8)", "     1.5")),
    ("int_underscore", P("%FLAG A", "%FORMAT(10I8)", "     1_0    -2_0   +0003  -00000")),
    ("int_bad_underscore", P("%FLAG A", "%FORMAT(10I8)", "    1__0")),
    ("int_too_large", P("%FLAG A", "%FORMAT(1I30)", "     99999999999999999999999")),
    ("int_limits", P("%FLAG A", "%FORMAT(2I24)", "     9223372036854775807    -9223372036854775808")),
    ("bad_float", P("%FLAG A", "%FORMAT(5E16.8)", "  1.0000000E+0x")),
    ("float_specials", P("%FLAG A", "%FORMAT(5E16.8)", "             nan            -inf        Infinity"
                         "           1e999          -1e999",
                         "".join(v.rjust(16) for v in ["1D-400", "1_0.5e1_0", "-0.0d0", " 1.5 ", ".5"]),
                         "".join(v.rjust(16) for v in ["4.9E-324", "2.4e-324", "1.7976931348623159e308", "nAn",
                                                       "+iNf"]))),
    ("float_bad_underscore", P("%FLAG A", "%FORMAT(5E16.8)", "            1__0")),
    ("tab_separated", P("%FLAG A", "%FORMAT(2I4)", "\t 1\t 2")),
    ("float_bad_d", P("%FLAG A", "%FORMAT(5E16.8)", "             1DD")),
    ("float_d_in_word", P("%FLAG A", "%FORMAT(2E16.8)", "             nDn     inDinity")),
    ("unwanted_bad_values", P("%FLAG A", "%FORMAT(10I8)", "     xyz", "%FLAG B", "%FORMAT(10I8)", "       5"), ["B"]),
    ("wanted_missing", P("%FLAG A", "%FORMAT(10I8)", "       1"), ["Z", "A"]),
    ("wanted_empty", P("%FLAG A", "%FORMAT(10I8)", "       1"), []),
    ("empty_sections", P("%FLAG A", "%FORMAT(10I8)", "", "%FLAG B", "%FORMAT(20a4)", "", "%FLAG C",
                         "%FORMAT(5E16.8)")),
    ("whitespace_only_a", P("%FLAG A", "%FORMAT(20a4)", "        ")),
]

random.seed(7)
for k in range(12):
    lines = []
    for s in range(random.randint(1, 4)):
        kind = random.choice("aief")
        width = random.randint(1, 12)
        lines.append(f"%FLAG S{s}")
        lines.append(f"%FORMAT({random.randint(1, 20)}{kind.upper() if random.random() < 0.5 else kind}{width}"
                     f"{'.3' if random.random() < 0.3 else ''})")
        for _ in range(random.randint(0, 3)):
            cells = []
            for _ in range(random.randint(0, 6)):
                if kind == "a":
                    v = "".join(random.choice("AB C'*") for _ in range(random.randint(0, width)))
                elif kind == "i":
                    v = str(random.randint(-10 ** (width - 1) + 1, 10 ** max(width - 1, 0)))
                else:
                    v = f"{random.uniform(-100, 100):.{max(width - 8, 0)}E}"
                    if random.random() < 0.3:
                        v = v.replace("E", "D")
                    if len(v) > width:
                        v = random.choice(["1", "0", "-2"])[:width]
                cells.append(v.rjust(width) if random.random() < 0.8 else " " * width)
            lines.append("".join(cells))
    READ.append((f"random_{k}", P(*lines, nl=random.choice(["\n", "\r\n"]))))


def read_cases():
    cases = []
    for entry in READ:
        name, data = entry[0], entry[1]
        wanted = entry[2] if len(entry) > 2 else None
        path = os.path.join(TMP, name + ".parm7")
        with open(path, "wb") as fh:
            fh.write(data)
        case = {"name": name, "bytes": data.hex(), "wanted": wanted}
        try:
            case["sections"] = sections_json(parm7.read_sections(path, wanted))
        except Exception as exc:  # noqa: BLE001
            case["error"] = error(exc, path)
        cases.append(case)
    for wanted in (None, ["POINTERS", "LENNARD_JONES_ACOEF", "AMBER_ATOM_TYPE", "NOPE"]):
        case = {"name": "demo_dry", "file": "topology/demo_dry.parm7", "wanted": wanted}
        case["sections"] = sections_json(parm7.read_sections(DEMO, wanted))
        cases.append(case)
    return cases


# ---------------------------------------------------------------- writer

def tag(values):
    """Values of one section: {"strings": [...]}, {"ints": ["..."]} or {"floats": [hex]}."""
    if all(isinstance(v, str) for v in values) and values:
        return {"strings": values}
    if all(isinstance(v, int) for v in values) and values:
        return {"ints": [str(v) for v in values]}
    if all(isinstance(v, float) for v in values):
        return {"floats": [fhex(v) for v in values]}
    raise TypeError(values)


FLOATS = [0.0, -0.0, 1.0, -1.5, 18.2223, 1e300, -1e-300, 5e-324, 2.2250738585072014e-308, 123456789.123456789,
          9.999999995, 9.999999994999, 0.1, 1.7976931348623157e308, 1e16, -12345.678901234567, 1.25e-5,
          float("nan"), -float("nan"), float("inf"), -float("inf"), 0.5e-8, 99999999.5, 1.000000005]
INTS = [0, 1, -1, 7, 12345678, -1234567, 99999999, -99999999, 123456789, 2 ** 63 - 1, -2 ** 63]
STRINGS = ["N", "CA", "H5''", "HO2'", "C1'", "", " ", "  X", "ABCDEFG", "a b", "\t", "x" * 100, "Na+"]

WRITE = [
    ("all_kinds", "MDATK", [("POINTERS", "i", INTS), ("ATOM_NAME", "a", STRINGS), ("CHARGE", "e", FLOATS),
                            ("TITLE2", "s", STRINGS)]),
    ("empty_sections", "", [("A", "a", []), ("B", "i", []), ("C", "e", []), ("D", "s", [])]),
    ("long_title", "T" * 120, [("A", "i", list(range(25)))]),
    ("short_title", "x", [("A", "i", list(range(10))), ("B", "i", list(range(11))), ("C", "e", FLOATS[:5]),
                          ("D", "e", FLOATS[:6]), ("E", "a", STRINGS[:3] * 7), ("F", "a", ["AB"] * 20)]),
    ("title_with_spaces", "  padded  title  ", []),
    ("unicode_title", "café demo", [("A", "i", [1])]),
    ("unicode_title_run", "abéè中cdé", [("A", "i", [1])]),
    ("unicode_title_after_80", "x" * 80 + "é", [("A", "i", [1])]),
    ("unicode_flag", "MDATK", [("A", "i", [1]), ("ÉT", "i", [1])]),
    ("unicode_value", "MDATK", [("A", "i", [1]), ("B", "a", ["N", "Cα", "O"])]),
    ("unicode_value_cut", "MDATK", [("B", "a", ["ABCDé", "\U0001F600"])]),
    ("unicode_s", "MDATK", [("B", "s", ["ok", "x" * 79 + "é" + "y"])]),
    ("unicode_s_cut", "MDATK", [("B", "s", ["x" * 80 + "é"])]),
    ("ints_as_a", "MDATK", [("A", "a", [1, -2, 123456])]),
    ("floats_as_a", "MDATK", [("A", "a", [1.5, -0.0, 1e300, float("nan"), float("inf"), 0.1])]),
    ("floats_as_s", "MDATK", [("A", "s", [1.5, 1e-7, 123456789012345678.0])]),
    ("ints_as_e", "MDATK", [("A", "e", [1, -2, 2 ** 63 - 1])]),
    ("floats_as_i", "MDATK", [("A", "i", [1.9, -1.9, 0.0, 1e18, -0.5])]),
    ("nan_as_i", "MDATK", [("A", "i", [1.0]), ("B", "i", [2.0, float("nan")])]),
    ("inf_as_i", "MDATK", [("A", "i", [float("-inf")])]),
    ("strings_as_i", "MDATK", [("A", "i", ["12", " -3 ", "1_000", "١٢"])]),
    ("bad_string_as_i", "MDATK", [("A", "i", ["12", "x"])]),
    ("strings_as_e", "MDATK", [("A", "e", ["1.5", " nan", "-inf", "1e5", "1_0.5"])]),
    ("bad_string_as_e", "MDATK", [("A", "e", ["1.5D0"])]),
    ("ints_as_s", "MDATK", [("A", "s", [1, 22, -333])]),
]


def writer_cases():
    cases = []
    for name, title, sections in WRITE:
        path = os.path.join(TMP, name + ".out.parm7")
        case = {"name": name, "title": title,
                "sections": [{"flag": f, "kind": k, "values": tag(v)} for f, k, v in sections]}
        try:
            parm7.write_parm7(path, sections, title=title)
        except Exception as exc:  # noqa: BLE001
            case["error"] = error(exc, path)
        with open(path, "rb") as fh:
            case["bytes"] = fh.read().hex()
        cases.append(case)
    fmt = []
    for kind, values in [("a", STRINGS), ("i", INTS), ("e", FLOATS), ("s", STRINGS), ("a", []), ("i", []),
                         ("e", []), ("s", []), ("a", ["ééééé"]), ("i", [1.5]),
                         ("e", [3]), ("a", [12345]), ("s", [1.0])]:
        fmt.append({"kind": kind, "values": tag(values), "lines": parm7._format_values(kind, values)})
    rng = random.Random(11)
    for _ in range(300):
        x = rng.choice([rng.uniform(-1e3, 1e3), 10 ** rng.uniform(-320, 308) * rng.choice([1, -1]),
                        rng.randint(-10 ** 9, 10 ** 9) / 10 ** rng.randint(0, 9),
                        float.fromhex(f"0x1.{rng.getrandbits(52):013x}p{rng.randint(-1074, 1023)}")])
        fmt.append({"kind": "e", "values": tag([x]), "lines": parm7._format_values("e", [x])})
    return cases, fmt


# ---------------------------------------------------------------- LJ table

def lj_cases():
    table = parm7.LJTable.from_parm7(DEMO)
    out = {
        "file": "topology/demo_dry.parm7",
        "type_index": [int(v) for v in table.type_index],
        "ntypes": table.ntypes,
        "nb_index": [int(v) for v in table.nb_index],
        "acoef": [fhex(v) for v in table.acoef],
        "bcoef": [fhex(v) for v in table.bcoef],
    }
    n = len(table.type_index)
    rng = np.random.default_rng(5)
    i = np.concatenate([np.arange(n), rng.integers(-n, n, 500), [0, 5, 30, 100]])
    j = np.concatenate([np.arange(n)[::-1], rng.integers(-n, n, 500), [1, 7, 64, 200]])
    a, b = table.pair_coefficients(i, j)
    out["pairs"] = {"i": [int(v) for v in i], "j": [int(v) for v in j], "a": [fhex(v) for v in a],
                    "b": [fhex(v) for v in b]}
    # MDATK's own test (tests/test_analyses.py::test_parm7_lj_table_roundtrip) reads the types MDAnalysis
    # gives; LJ_TYPES of mdatk/demo.py are the (rmin/2, eps) of each type
    from mdatk.demo import LJ_TYPES
    u = mda.Universe(DEMO)
    out["types"] = [str(t) for t in u.atoms.types]
    out["lj_types"] = {k: [fhex(r), fhex(e)] for k, (r, e) in LJ_TYPES.items()}
    errors = []
    text = DEMO_TEXT  # the small system
    flags = ["POINTERS", "ATOM_TYPE_INDEX", "NONBONDED_PARM_INDEX", "LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF"]
    for drop in ([flags[2]], [flags[4], flags[0]], flags, [flags[1], flags[3]]):
        t = text
        for f in drop:
            t = t.replace(f"%FLAG {f}\n", f"%FLAG X{f}\n")
        path = os.path.join(TMP, "lj_missing.parm7")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
        try:
            parm7.LJTable.from_parm7(path)
        except ValueError as exc:
            errors.append({"bytes": t.encode().hex(), "error": error(exc, path)})
    out["errors"] = errors
    return out


# ---------------------------------------------------------------- MDAnalysis topologies

def split_sections(text):
    head, secs, cur = [], [], None
    for line in text.splitlines():
        if line.startswith("%FLAG"):
            cur = [line[6:].strip(), None, []]
            secs.append(cur)
        elif line.startswith("%FORMAT") and cur is not None and cur[1] is None:
            cur[1] = line
        elif cur is None:
            head.append(line)
        else:
            cur[2].append(line)
    return head, secs


def join_sections(head, secs, nl="\n"):
    out = list(head)
    for name, fmt, lines in secs:
        out.append("%FLAG " + name)
        if fmt is not None:
            out.append(fmt)
        out.extend(lines)
    return nl.join(out) + nl


def get(secs, name):
    return next(s for s in secs if s[0] == name)


def index(secs, name):
    return [s[0] for s in secs].index(name)


def chain_section(values, fmt="%FORMAT(20a4)", width=4, per=20):
    cells = [f"{v:<{width}s}" for v in values]
    return ["RESIDUE_CHAINID", fmt, ["".join(cells[k:k + per]) for k in range(0, len(cells), per)] or [""]]


def small_system():
    """A 40-atom, 8-residue system with every section TOPParser needs, written by write_parm7."""
    rng = random.Random(3)
    natom, nres = 40, 8
    names = [rng.choice(["N", "CA", "C", "O", "H", "HO2'", "P", "OP1", "MG", "C1'"]) for _ in range(natom)]
    types = [rng.choice(["N", "CT", "C", "O", "H", "HO", "P", "O2", "MG"]) for _ in range(natom)]
    znum = [rng.choice([1, 6, 7, 8, 12, 15, 11, 17, 30]) for _ in range(natom)]
    charges = [rng.uniform(-1.2, 1.2) * parm7.AMBER_CHARGE_FACTOR for _ in range(natom)]
    charges[:4] = [0.0, -0.0, 18.2223, -18.2223]
    masses = [rng.uniform(1.0, 40.0) for _ in range(natom)]
    respoint = [1, 6, 11, 12, 20, 27, 33, 38]
    resnames = ["G5", "A", "C", "MG", "U", "G", "NA", "CL"]
    bonds_h = [(k, k + 1) for k in range(0, 10, 2)]
    bonds_a = [(k, k + 3) for k in range(10, 36, 3)] + [(39, 0)]
    angles_h = [(0, 1, 2), (3, 4, 5)]
    angles_a = [(10, 11, 12)]
    dihed_h = [(0, 1, 2, 3)]
    dihed_a = [(10, 11, 12, 13), (10, 11, -12, 13), (14, 15, 16, -17)]
    pointers = [0] * 31
    pointers[0], pointers[1], pointers[11] = natom, 3, nres
    pointers[2], pointers[3], pointers[12] = len(bonds_h), len(bonds_a), len(bonds_a)
    pointers[4], pointers[5] = len(angles_h), len(angles_a)
    pointers[6], pointers[7] = len(dihed_h), len(dihed_a)

    def bonded(items):
        out = []
        for item in items:
            out.extend(3 * v if v >= 0 else -3 * -v for v in item)
            out.append(1)
        return out

    sections = [
        ("POINTERS", "i", pointers), ("ATOM_NAME", "a", names), ("CHARGE", "e", charges),
        ("ATOMIC_NUMBER", "i", znum), ("MASS", "e", masses),
        ("ATOM_TYPE_INDEX", "i", [rng.randint(1, 3) for _ in range(natom)]),
        ("NONBONDED_PARM_INDEX", "i", [1, 2, 4, 2, 3, 5, 4, 5, 6]),
        ("RESIDUE_LABEL", "a", resnames), ("RESIDUE_POINTER", "i", respoint),
        ("BONDS_INC_HYDROGEN", "i", bonded(bonds_h)), ("BONDS_WITHOUT_HYDROGEN", "i", bonded(bonds_a)),
        ("ANGLES_INC_HYDROGEN", "i", bonded(angles_h)), ("ANGLES_WITHOUT_HYDROGEN", "i", bonded(angles_a)),
        ("DIHEDRALS_INC_HYDROGEN", "i", bonded(dihed_h)), ("DIHEDRALS_WITHOUT_HYDROGEN", "i", bonded(dihed_a)),
        ("LENNARD_JONES_ACOEF", "e", [1.0e6, 2.0e5, 3.0e4, 4.0e3, 5.0e2, 6.0e1]),
        ("LENNARD_JONES_BCOEF", "e", [1.0e3, 2.0e2, 3.0e1, 4.0, 5.0e-1, 6.0e-2]),
        ("AMBER_ATOM_TYPE", "a", types),
    ]
    path = os.path.join(TMP, "small.parm7")
    parm7.write_parm7(path, sections, title="small test system")
    with open(path, encoding="ascii") as fh:
        return fh.read()


DEMO_TEXT = small_system()
NRES = 8


def demo_variant(fn, nl="\n"):
    head, secs = split_sections(DEMO_TEXT)
    fn(head, secs)
    return join_sections(head, secs, nl)


def after(name, section):
    def fn(head, secs):
        secs.insert(index(secs, name) + 1, section)
    return fn


def set_atomic_numbers(values):
    def fn(head, secs):
        sec = get(secs, "ATOMIC_NUMBER")
        first = sec[2][0].split()
        first[:len(values)] = [str(v) for v in values]
        sec[2][0] = "".join(f"{int(x):8d}" for x in first)
    return fn


def reformat(name, fmt, per, width, render):
    def fn(head, secs):
        sec = get(secs, name)
        vals = " ".join(sec[2]).split()
        cells = [render(v, width) for v in vals]
        sec[1] = fmt
        sec[2] = ["".join(cells[k:k + per]) for k in range(0, len(cells), per)]
    return fn


def drop(name):
    def fn(head, secs):
        secs.remove(get(secs, name))
    return fn


def edit_line(name, k, new):
    def fn(head, secs):
        get(secs, name)[2][k] = new
    return fn


def chain(*fns):
    def fn(head, secs):
        for f in fns:
            f(head, secs)
    return fn


def touching(v, w):
    return f"{float(v):.8E}".rjust(w)[-w:]


def unicode_names(head, secs):
    sec = get(secs, "ATOM_NAME")
    line = sec[2][0]
    sec[2][0] = "Été" + "α" + line[4:]


def comment_before_format(head, secs):
    sec = get(secs, "CHARGE")
    sec[1] = "%COMMENT charges in units of e * 18.2223\n%COMMENT second comment\n" + sec[1]


def empty_bond_section(head, secs):
    # BONDS_INC_HYDROGEN moved to the end of BONDS_WITHOUT_HYDROGEN; NBONH = 0
    hyd = get(secs, "BONDS_INC_HYDROGEN")
    heavy = get(secs, "BONDS_WITHOUT_HYDROGEN")
    values = " ".join(heavy[2]).split() + " ".join(hyd[2]).split()
    cells = [f"{int(v):8d}" for v in values]
    heavy[2] = ["".join(cells[k:k + 10]) for k in range(0, len(cells), 10)]
    hyd[2] = [""]
    pointers = get(secs, "POINTERS")
    vals = " ".join(pointers[2]).split()
    nbonh, mbona = int(vals[2]), int(vals[3])
    vals[2], vals[3] = "0", str(nbonh + mbona)
    vals[12] = str(int(vals[12]) + nbonh)  # NBONA
    cells = [f"{int(v):8d}" for v in vals]
    pointers[2] = ["".join(cells[k:k + 10]) for k in range(0, len(cells), 10)]


def residue_pointer_offset(head, secs):
    sec = get(secs, "RESIDUE_POINTER")
    vals = " ".join(sec[2]).split()
    vals[0] = "3"
    cells = [f"{int(v):8d}" for v in vals]
    sec[2] = ["".join(cells[k:k + 10]) for k in range(0, len(cells), 10)]


E = (lambda v, w: f"{float(v):{w}.8E}")
I8 = (lambda v, w: f"{int(v):{w}d}")
TOPOLOGY_VARIANTS = [
    ("small", demo_variant(lambda h, s: None)),
    ("small_crlf", demo_variant(lambda h, s: None, nl="\r\n")),
    ("small_cr", demo_variant(lambda h, s: None, nl="\r")),
    ("chain_all_a", demo_variant(after("RESIDUE_POINTER", chain_section(["A"] * NRES)))),
    ("chain_runs", demo_variant(after("RESIDUE_POINTER", chain_section(["A"] * 2 + ["B"] * 3 + ["A"] * 3)))),
    ("chain_each", demo_variant(after("RESIDUE_POINTER", chain_section([chr(65 + k % 26) for k in range(NRES)])))),
    ("chain_blank", demo_variant(after("RESIDUE_POINTER", chain_section(["A"] * 2 + [""] + ["B"] * 5)))),
    ("chain_short", demo_variant(after("RESIDUE_POINTER", chain_section(["A"] * (NRES - 1))))),
    ("chain_long", demo_variant(after("RESIDUE_POINTER", chain_section(["A"] * (NRES + 1))))),
    ("chain_last", demo_variant(lambda h, s: s.append(chain_section(["X"] * 4 + ["Y"] * 4)))),
    ("chain_before_pointers", demo_variant(lambda h, s: s.insert(1, chain_section(["X"] * 4 + ["Y"] * 4)))),
    ("chain_long_names", demo_variant(after("RESIDUE_POINTER", chain_section(["AB", "C D"] + ["LONG"] * 6)))),
    ("chain_wide_format", demo_variant(after("RESIDUE_POINTER", chain_section(
        ["A1"] * 5 + ["B2"] * 3, fmt="%FORMAT(20a8)", width=8)))),
    ("chain_comment", demo_variant(after("RESIDUE_POINTER", ["RESIDUE_CHAINID", "%COMMENT chains\n%FORMAT(20a4)",
                                                             chain_section(["P"] * 3 + ["Q"] * 5)[2]]))),
    ("chain_twice", demo_variant(chain(after("RESIDUE_POINTER", chain_section(["A"] * NRES)),
                                       after("MASS", chain_section(["B"] * 4 + ["C"] * 4))))),
    ("resptr_zero", demo_variant(edit_line("RESIDUE_POINTER", 0, "       0       6      11      12"
                                                                 "      20      27      33      38"))),
    ("resptr_unsorted", demo_variant(edit_line("RESIDUE_POINTER", 0, "       1      11       6      12"
                                                                     "      20      27      33      38"))),
    ("mass_blank_field", demo_variant(edit_line("MASS", 0, " " * 16 + get(split_sections(DEMO_TEXT)[1],
                                                                          "MASS")[2][0][16:]))),
    ("title_one_token", DEMO_TEXT.replace("%FLAG TITLE", "%FLAG", 1)),
    ("no_pointers", DEMO_TEXT.replace("%FLAG POINTERS", "%FLAG POINTERZ", 1)),
    ("pointer_lines_short", demo_variant(lambda h, s: get(s, "POINTERS")[2].pop())),
    ("not_flag_after_section", demo_variant(lambda h, s: get(s, "MASS")[2].append("trailing junk"))),
    ("empty_line_after_section", demo_variant(lambda h, s: get(s, "MASS")[2].append(""))),
    ("atomic_number_zero", demo_variant(set_atomic_numbers([0, -1, 12, 1, 200 - 200, 118, 54]))),
    ("no_atomic_number", demo_variant(drop("ATOMIC_NUMBER"))),
    ("charge_touching", demo_variant(reformat("CHARGE", "%FORMAT(5E15.8)", 5, 15, touching))),
    ("charge_wide", demo_variant(reformat("CHARGE", "%FORMAT(5E20.12)", 5, 20, lambda v, w: f"{float(v):{w}.12E}"))),
    ("mass_d_free", demo_variant(reformat("MASS", "%FORMAT(5F16.4)", 5, 16, lambda v, w: f"{float(v):{w}.4f}"))),
    ("comment_before_format", demo_variant(comment_before_format)),
    ("unknown_first", demo_variant(after("POINTERS", ["MY_SECTION", "%FORMAT(10I8)", ["       1"]]))),
    ("unknown_between", demo_variant(after("MASS", ["MY_SECTION", "%FORMAT(20a4)", ["abc", "def"]]))),
    ("empty_bond_section", demo_variant(empty_bond_section)),
    ("residue_pointer_offset", demo_variant(residue_pointer_offset)),
    ("unicode_names", demo_variant(unicode_names)),
    ("names_trailing_text", demo_variant(edit_line("ATOM_NAME", 0, get(split_sections(DEMO_TEXT)[1],
                                                                         "ATOM_NAME")[2][0] + "EXTRA"))),
    ("charge_three_per_line", demo_variant(reformat("CHARGE", "%FORMAT(3E24.16)", 3, 24, E))),
    ("bond_twelve_per_line", demo_variant(reformat("BONDS_WITHOUT_HYDROGEN", "%FORMAT(12I6)", 12, 6, I8))),
    ("atomic_number_119", demo_variant(set_atomic_numbers([119]))),
    ("no_atom_name", demo_variant(drop("ATOM_NAME"))),
    ("no_angles", demo_variant(drop("ANGLES_WITHOUT_HYDROGEN"))),
    ("no_version", DEMO_TEXT.split("\n", 1)[1]),
    ("no_title", DEMO_TEXT.replace("%FLAG TITLE", "%FLAG TITEL", 1)),
    ("ctitle", DEMO_TEXT.replace("%FLAG TITLE", "%FLAG CTITLE", 1)),
    ("bad_format_line", demo_variant(lambda h, s: get(s, "MASS").__setitem__(1, "%FROMAT(5E16.8)"))),
    ("bad_int", demo_variant(edit_line("RESIDUE_POINTER", 0, "       1      1x"))),
]


def topology_cases():
    folder = os.path.join(OUT, "topology")
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    errors = []
    for name, text in TOPOLOGY_VARIANTS:
        path = os.path.join(folder, name + ".parm7")
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        try:
            TOPParser(path).parse()
        except Exception as exc:  # noqa: BLE001
            errors.append({"file": name + ".parm7", "stage": "parser", "error": error(exc, path)})
            continue
        try:
            data = dump_topologies.dump(path)
        except Exception as exc:  # noqa: BLE001  (Universe checks the attribute lengths)
            errors.append({"file": name + ".parm7", "stage": "universe", "error": error(exc, path)})
            continue
        with open(path + ".topology.json", "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=1, ensure_ascii=False)
            fh.write("\n")
    return errors


def is_parm7_cases():
    paths = ["a.parm7", "a.PARM7", "a.Prmtop", "a.top", "x/y.TOP", "a.parm7.bak", ".parm7", "dir.parm7/file",
             "a.", "a", "a.pdb", "a.parm", "/abs/path/sys.prmtop", "a..top", "..top", "a/.top", "a/b.c/d",
             "a.top/", "a.İ", "é.parm7", "x.parm7 ", "x.top\n", "a.tOp"]
    return [[p, os.path.splitext(p)[1].lower() in (".parm7", ".prmtop", ".top")] for p in paths]


def main():
    os.makedirs(OUT, exist_ok=True)
    writer, fmt = writer_cases()
    files = {
        "read_cases.json": read_cases(),
        "writer_cases.json": {"write": writer, "format_values": fmt},
        "lj_cases.json": lj_cases(),
        "topology_errors.json": topology_cases(),
        "is_parm7_cases.json": is_parm7_cases(),
    }
    for name, data in files.items():
        with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=None if name == "lj_cases.json" else 1, ensure_ascii=False)
            fh.write("\n")
    shutil.rmtree(TMP)


if __name__ == "__main__":
    main()
