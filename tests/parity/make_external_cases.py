"""What MDATK's external-series reader (mdatk/io/external.py) and numpy.interp compute
(tests/data/io/external).

    /opt/mdbridge-venv/bin/python -I -B make_external_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy; nothing is written into the checkout. Writes

  table_cases.json      read_table of crafted files (bytes as hex) with and without a delimiter:
                        CSV, tabs, whitespace, per-line detection, comments (# @ ; %), xvg files,
                        header fallbacks, ragged rows, Unicode whitespace and digits, CR / CRLF,
                        invalid UTF-8, errors
  column_cases.json     _column_index (names, indices, Unicode digits, errors) and _sanitize
  interp_cases.json     numpy.interp(x, xp, fp, left, right) of random inputs with duplicates, NaN
                        and infinities in every argument, and its errors
  load_cases.json       load_external of files and specs on datasets: the arrays and repr(meta) of
                        the result, or the error

Floats are repr() strings ("nan", "inf", "-inf"), which read back exactly. Errors are str(exc)
with the temporary folder replaced by "{dir}".
"""

import json
import math
import os
import random
import shutil
import sys
import tempfile

if len(sys.argv) != 3:
    sys.exit("usage: make_external_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import numpy as np  # noqa: E402

from mdatk.config import ExternalSeries  # noqa: E402
from mdatk.data import Dataset  # noqa: E402
from mdatk.io import external  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import improvements_core_io  # noqa: E402

improvements_core_io.apply_external()  # two series of one key keep both (EXT-01, MDATK 0.1.0)

TMP = tempfile.mkdtemp(prefix="external-cases-")


def num(x):
    x = float(x)
    if math.isnan(x):
        return "nan"
    return repr(x)


def err(exc):
    return str(exc).replace(TMP, "{dir}")


# ---------------------------------------------------------------- read_table

def B(text):
    return text.encode("utf-8")


TABLES = [
    ("mdatk_test.dat", B("#Frame  rmsd  other\n1 0.0 5\n2 1.0 6\n3 2.0 7\n"), ""),
    ("mdatk_test.csv", B("time_ns,energy\n0,1\n0.02,3\n"), ""),
    ("tabs.tsv", B("a\tb\tc\n1\t2\t3\n4\t5\t6\n"), ""),
    ("tabs_spaces.tsv", B("a\t b \t c\n 1 \t2\t3\n"), ""),
    ("per_line.txt", B("x y\n1,2\n3\t4\n5 6\n"), ""),
    ("explicit_semicolon.txt", B("a;b\n1;2\n3; 4\n"), ";"),
    ("explicit_space.txt", B("a b\n1 2\n3  4\n"), " "),
    ("explicit_multi.txt", B("a::b\n1::2\n"), "::"),
    ("explicit_comma_in_ws.txt", B("1 2\n3 4\n"), ","),
    ("comments.dat", B("# first\n@ title\n; semi\n% percent\n#  col_a   col_b\n1 2\n# after data\n3 4\n"), ""),
    ("xvg.xvg", B("# GROMACS\n@    title \"RMSD\"\n@    xaxis  label \"Time (ps)\"\n@TYPE xy\n0.0 0.1\n"
                  "10.0 0.2\n&\n"), ""),
    ("xvg_no_amp.xvg", B("@ s0 legend \"rmsd\"\n0 1\n1 2\n"), ""),
    ("comment_csv.csv", B("# t, value\n0,1\n1,2\n"), ""),
    ("comment_width_mismatch.dat", B("# a b c\n1 2\n"), ""),
    ("header_and_comment.dat", B("# p q\nx y\n1 2\n"), ""),
    ("header_width_mismatch.dat", B("x y z\n1 2\n3 4\n"), ""),
    ("two_headers.dat", B("x y\nu v\n1 2\n"), ""),
    ("header_after_data.dat", B("1 2\nx y\n"), ""),
    ("ragged.dat", B("1 2 3\n4 5\n6 7 8 9\n"), ""),
    ("empty.dat", b"", ""),
    ("blank_lines.dat", B("\n   \n\t\n"), ""),
    ("only_comments.dat", B("# a b\n@ c\n"), ""),
    ("only_header.dat", B("a b\n"), ""),
    ("unicode_ws.dat", B("t v\n1 2\n3　4\n 5 6 \n"), ""),
    ("ctrl_ws.dat", B("1\x1c2\n3\x0c4\n5\x0b6\n"), ""),
    ("numbers.dat", B("1_000 nan\ninf -Infinity\n 1e5 +1.5\n١٢ .5\n"), ""),
    ("bad_underscore.dat", B("1 2\n1__0 3\n"), ""),
    ("crlf.csv", B("a,b\r\n1,2\r\n3,4\r\n"), ""),
    ("cr.csv", B("a,b\r1,2\r3,4\r"), ""),
    ("invalid_utf8.dat", b"# t \xff\xfe val\n1 2\n\xc3( 3\n", ""),
    ("invalid_utf8_header.dat", b"t \xe2\x82 v\n1 2\n", ""),
    ("hash_only.dat", B("#\n1 2\n"), ""),
    ("empty_tokens.csv", B("a,b\n1,,2\n"), ""),
    ("empty_token_first.csv", B("1,,2\n3,4,5\n"), ""),
    ("trailing_comma.csv", B("1,2,\n"), ""),
    ("comma_spaces.csv", B("a , b\n1 , 2\n"), ""),
    ("leading_ws.dat", B("   # c d\n   1 2\n"), ""),
    ("comment_strip_chars.dat", B("#@;% a b\n1 2\n"), ""),
    ("comment_inner_hash.dat", B("# a #b\n1 2\n"), ""),
    ("comment_tabs.dat", B("#\ta\tb\n1\t2\n"), ""),
    ("comment_after_header.dat", B("x y\n# p q\n1 2\n"), ""),
    ("numeric_header.dat", B("1 2\n3 4\n"), ""),
    ("bad_float_e.dat", B("1 2\n3 1e\n"), ""),
    ("unicode_header.dat", B("temps énergie ΔG\n0 1 2\n"), ""),
    ("one_column.dat", B("v\n1\n2\n"), ""),
    ("long_header.dat", B("a b c d\n1 2\n"), ""),
    ("bom.csv", b"\xef\xbb\xbfa,b\n1,2\n", ""),
    ("nul_byte.dat", B("a\x00b c\n1 2\n"), ""),
    ("vertical_tab_line.dat", B("a b\n1 2\x0b3 4\n"), ""),
]

random.seed(4)
for k in range(10):
    lines = []
    width = random.randint(1, 5)
    if random.random() < 0.5:
        lines.append(random.choice(["#", "@", ";", "%", ""]) + " ".join(f"c{n}" for n in range(width)))
    for _ in range(random.randint(1, 6)):
        sep = random.choice([" ", "  ", "\t", ",", ", "])
        cells = [random.choice([str(random.randint(-9, 99)), f"{random.uniform(-1e3, 1e3):.4g}", "1e-3", "nan"])
                 for _ in range(width + random.randint(-1, 1) if width > 1 else width)]
        lines.append(sep.join(cells))
        if random.random() < 0.2:
            lines.append("# mid")
    TABLES.append((f"random_{k}.dat", B(random.choice(["\n", "\r\n"]).join(lines) + "\n"), ""))


def table_cases():
    out = []
    for name, data, delimiter in TABLES:
        path = os.path.join(TMP, name)
        with open(path, "wb") as fh:
            fh.write(data)
        case = {"name": name, "bytes": data.hex(), "delimiter": delimiter}
        try:
            header, rows = external.read_table(path, delimiter)
            case["header"] = header
            case["rows"] = [[num(v) for v in row] for row in rows]
            case["shape"] = list(rows.shape)
        except ValueError as exc:
            case["error"] = err(exc)
        out.append(case)
    return out


# ---------------------------------------------------------------- columns, sanitize

def column_cases():
    headers = [["Frame", "rmsd", "1", "x"], ["a", "a", "b"], ["col0", "col1"], ["é", "2"], []]
    specs = ["", "rmsd", "Frame", "1", "0", "3", "4", "٣", "１", "01", "-1", " 1", "1 ", "1\n", "x",
             "99999999999999999999999", "²", "a", "b", "col1", "é", "2", "frame", "1.0", "+1", "٠"]
    out = []
    for header in headers:
        for spec in specs:
            for default in (0, -1):
                try:
                    out.append({"header": header, "spec": spec, "default": default,
                                "index": external._column_index(header, spec, default)})
                except ValueError as exc:
                    out.append({"header": header, "spec": spec, "default": default, "error": str(exc)})
    texts = ["cpp.rmsd", "cpp: rmsd", "", "___", "a b c", "énergie", "xΔG", "--x++", "a..b", "_a_",
             "\U0001F600", "tab\tsep", "a/b\\c", "ok-name_1.2+3", "  ", "é", "aééb", "x\n", "\x00"]
    sanitize = [[t, external._sanitize(t)] for t in texts]
    return {"column_index": out, "sanitize": sanitize}


# ---------------------------------------------------------------- numpy.interp

def interp_cases():
    rng = random.Random(9)
    pool = [0.0, -0.0, 1.0, 2.0, 2.5, 3.0, -1.0, 1e300, -1e300, 5e-324, float("inf"), float("-inf"), float("nan")]
    cases = []

    def value(kind):
        r = rng.random()
        if kind == "grid":
            return float(rng.randint(-3, 12)) / rng.choice([1, 2, 4])
        if r < 0.08:
            return rng.choice(pool)
        if r < 0.5:
            return round(rng.uniform(-5, 15), rng.randint(0, 3))
        return float(rng.randint(-5, 15))

    for n in range(2000):
        nxp = rng.choice([1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 20, 40])
        nx = rng.choice([0, 1, 2, 3, 5, 8, 13, 30, 60])
        kind = "grid" if rng.random() < 0.5 else "mixed"
        xp = sorted(value(kind) for _ in range(nxp))
        if rng.random() < 0.15:
            rng.shuffle(xp)  # unsorted sample points: numpy's result is still defined
        if rng.random() < 0.2 and nxp > 1:
            k = rng.randrange(nxp - 1)
            xp[k + 1] = xp[k]  # duplicate
        fp = [value("mixed") for _ in range(nxp)]
        if rng.random() < 0.5:
            x = [value(kind) for _ in range(nx)]
        else:
            x = sorted(value(kind) for _ in range(nx))
        if rng.random() < 0.2 and nx:
            x[rng.randrange(nx)] = rng.choice(xp)
        left = float("nan") if rng.random() < 0.5 else value("mixed")
        right = float("nan") if rng.random() < 0.5 else value("mixed")
        res = np.interp(np.array(x, dtype=np.float64), np.array(xp, dtype=np.float64),
                        np.array(fp, dtype=np.float64), left=left, right=right)
        cases.append([[num(v) for v in x], [num(v) for v in xp], [num(v) for v in fp], num(left), num(right),
                      [num(v) for v in res]])
    # longer inputs: the guessed binary search over many sample points
    for n in range(16):
        nxp = rng.choice([50, 100, 300])
        xp = sorted(round(rng.uniform(0, 1000), 1) for _ in range(nxp))
        fp = [round(rng.uniform(-10, 10), 3) for _ in range(nxp)]
        x = [round(rng.uniform(-50, 1050), 2) for _ in range(rng.choice([20, 200]))]
        if n % 2:
            x.sort()
        res = np.interp(x, xp, fp, left=np.nan, right=np.nan)
        cases.append([[num(v) for v in x], [num(v) for v in xp], [num(v) for v in fp], "nan", "nan",
                      [num(v) for v in res]])
    errors = []
    for x, xp, fp in [([1.0], [], []), ([1.0], [1.0, 2.0], [1.0]), ([], [], []), ([1.0], [1.0], [1.0, 2.0]),
                      ([], [], [1.0]), ([1.0], [], [1.0]), ([], [1.0], [])]:
        try:
            np.interp(x, xp, fp, left=np.nan, right=np.nan)
            errors.append([x, xp, fp, None])
        except ValueError as exc:
            errors.append([x, xp, fp, str(exc)])
    return {"cases": cases, "errors": errors}


# ---------------------------------------------------------------- load_external

FILES = {
    "rmsd.dat": "#Frame  rmsd  other\n1 0.0 5\n2 1.0 6\n3 2.0 7\n",
    "e.csv": "time_ns,energy\n0,1\n0.02,3\n",
    "unsorted.dat": "# t a b\n30 3 30\n10 1 10\nnan 9 90\n20 2 20\n0 0 0\n20 5 50\n",
    "dupes.dat": "t v v w\n0 1 2 3\n10 4 5 6\n",
    "noheader.dat": "0 1 2\n100 3 4\n200 5 6\n",
    "a.b.dat": "t,x y\n0,1\n1000,2\n",
    "odd name!.dat": "t v\n0 1\n5 2\n",
    "fs.dat": "t v\n0 0\n1e6 1\n",
    "frames.dat": "f v\n0 0\n1 1\n2 4\n3 9\n",
    "spaces.dat": "time value (kcal/mol)\n0 1\n",
    "bad.dat": "t v\n0 x\n",
}

DATASETS = {
    "mdatk": ([0, 1, 2], [0.0, 10.0, 20.0]),
    "wide": ([0, 1, 2, 3, 5], [-5.0, 0.0, 10.0, 25.0, 1000.0]),
    "empty": ([], []),
}


def S(path, **kw):
    return dict(path=path, **kw)


LOADS = [
    ("mdatk_test_1", "mdatk", [S("rmsd.dat", label="cpp", x_column="Frame", x_unit="frame1", columns=["rmsd"])]),
    ("mdatk_test_2", "mdatk", [S("e.csv", x_unit="ns")]),
    ("all_columns", "wide", [S("rmsd.dat", x_unit="frame1")]),
    ("skip_empty_path", "wide", [S("", label="none"), S("rmsd.dat", x_unit="frame", columns=["2", "rmsd"])]),
    ("unsorted_nan", "wide", [S("unsorted.dat")]),
    ("duplicate_names", "wide", [S("dupes.dat"), S("dupes.dat", columns=["v", "2", "w"], label="d2")]),
    ("fallback_header", "wide", [S("noheader.dat", x_unit="ps")]),
    ("last_column", "wide", [S("noheader.dat", columns=[""])]),
    ("label_from_name", "wide", [S("a.b.dat", delimiter=",", x_unit="us")]),
    ("odd_label", "wide", [S("odd name!.dat", label="  my label/é "), S("odd name!.dat", label="___")]),
    ("fs_unit", "wide", [S("fs.dat", x_unit="fs")]),
    ("empty_unit", "wide", [S("fs.dat", x_unit="")]),
    ("frames", "wide", [S("frames.dat", x_unit="frame"), S("frames.dat", x_unit="frame1", label="f1")]),
    ("x_by_index", "wide", [S("rmsd.dat", x_column="2", x_unit="frame")]),
    ("x_is_column", "wide", [S("rmsd.dat", columns=["Frame", "0"])]),
    ("same_key_twice", "wide", [S("rmsd.dat", label="x"), S("rmsd.dat", label="x")]),
    ("spaces_header", "wide", [S("spaces.dat")]),
    ("empty_dataset", "empty", [S("rmsd.dat")]),
    ("bad_unit", "wide", [S("rmsd.dat", x_unit="ms")]),
    ("bad_unit_frame_ok", "wide", [S("rmsd.dat", x_unit="frames")]),
    ("missing_column", "wide", [S("rmsd.dat", columns=["nope"])]),
    ("missing_x_column", "wide", [S("rmsd.dat", x_column="time")]),
    ("index_out_of_range", "wide", [S("rmsd.dat", columns=["3"])]),
    ("bad_file", "wide", [S("rmsd.dat"), S("bad.dat")]),
    ("missing_file", "wide", [S("no_such_file.dat")]),
    ("no_specs", "wide", []),
]


def load_cases():
    for name, text in FILES.items():
        with open(os.path.join(TMP, name), "w", encoding="utf-8") as fh:
            fh.write(text)
    out = []
    for name, ds_name, specs in LOADS:
        frames, times = DATASETS[ds_name]
        ds = Dataset(np.array(frames, dtype=np.int64), np.array(times, dtype=np.float64))
        objs = [ExternalSeries(**dict(s, path=os.path.join(TMP, s["path"]) if s["path"] else "")) for s in specs]
        case = {"name": name, "frames": frames, "time_ps": [num(t) for t in times], "specs": specs}
        try:
            res = external.load_external(objs, ds)
            case["arrays"] = [[k, [num(v) for v in a]] for k, a in res.arrays.items()]
            case["meta_repr"] = repr(res.meta)
        except (ValueError, OSError) as exc:
            case["error"] = err(exc)
            case["error_type"] = type(exc).__name__
        out.append(case)
    return {"files": FILES, "cases": out}


def main():
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    files = {
        "table_cases.json": (table_cases(), 1),
        "column_cases.json": (column_cases(), 1),
        "interp_cases.json": (interp_cases(), None),
        "load_cases.json": (load_cases(), 1),
    }
    for name, (data, indent) in files.items():
        with open(os.path.join(out_dir, name), "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=indent, ensure_ascii=False)
            fh.write("\n")
    shutil.rmtree(TMP)


if __name__ == "__main__":
    main()
