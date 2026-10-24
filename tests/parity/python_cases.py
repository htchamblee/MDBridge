"""Write the test vectors of MDBridge's Python layer (core: python, json, hash, naming).

    python3 -I -B tests/parity/python_cases.py MDATK_CHECKOUT OUT_DIR

Run with Python 3.13 (MDATK's reference interpreter; json's messages and the Unicode data differ in
other versions). -I / -B: the MDATK checkout is put on sys.path explicitly and gets no byte code.
Writes into OUT_DIR (tests/data/python):

The vector files of the porting specs, rebuilt byte for byte (generators of the specs: section 4-6
and 8 of gen_config_model.py, round2_vectors of make_vectors.py, naming_cases of
dump_tool_metadata.py):

* ``python_repr_cases.json``   repr() / str() / type().__name__ of every TOML value type
* ``stable_hash_cases.json``   config.stable_hash: input (as repr), the JSON blob, sha256, hashes of
                               length 0 / 8 / 16 / 64 / 100; file_fingerprint of a file and a missing file
* ``resolve_paths_cases.json`` Project.resolve_paths with HOME=/home/tester and base /base/proj
* ``dirs_cases.json``          out_dir / cache_dir / workflow_dir, sanitize_project, group()
* ``round2_vectors.json``      repr(round(x, 2)) of 732 doubles (given as float.hex())
* ``naming_cases.json``        mdatk.workflows.naming

New here:

* ``stable_hash_inputs.json``  the inputs of stable_hash_cases.json as data: tagged JSON, see below
* ``posixpath_cases.json``     posixpath functions over awkward paths; expanduser with several HOME values
* ``string_repr_cases.json``   repr() and json.dumps() of characters on both sides of every printable /
                               non-printable boundary, and of undecodable bytes (surrogateescape), as UTF-8 hex
* ``json_loads_cases.json``    json.loads of valid and invalid texts: repr() of the result or the message
* ``round_cases.json``         repr(round(x, n)) for n from -309 to 324 (float.hex() inputs)
* ``sha256_cases.json``        hashlib.sha256 of "a" * n (n = 0..200, block boundaries) and of longer texts
* ``naming_extra_cases.json``  naming.parse_name / sanitize_* / new_stamp of further edge cases

Tagged JSON (stable_hash_inputs.json): {"$float": "nan" | "inf" | "-inf"}, {"$datetime": isoformat},
{"$date": isoformat}, {"$int": decimal} for an integer beyond 64 bits; tuples are lists and int keys
strings (what a core::Value holds); a dataclass is given as its _to_plain() form and marked
"dataclass": true (_to_plain keeps its None fields).
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import hashlib
import json
import math
import os
import posixpath
import random
import shutil
import struct
import sys
import tempfile
import unicodedata

if len(sys.argv) != 3:
    sys.exit("usage: python3 -I -B tests/parity/python_cases.py MDATK_CHECKOUT OUT_DIR")
MDATK, OUT = (os.path.abspath(a) for a in sys.argv[1:3])
sys.dont_write_bytecode = True
sys.path.insert(0, MDATK)

from mdatk import config as C  # noqa: E402
from mdatk import tomlio  # noqa: E402
from mdatk.workflows import naming  # noqa: E402
from mdatk.workflows.naming import sanitize_project  # noqa: E402

os.makedirs(OUT, exist_ok=True)


# ------------------------------------------------------------------------------------------ helpers


def _finite(o):
    """Strict JSON: a non-finite float becomes {"$float": "nan" | "inf" | "-inf"}."""
    if isinstance(o, float) and not math.isfinite(o):
        return {"$float": "nan" if math.isnan(o) else ("inf" if o > 0 else "-inf")}
    if isinstance(o, dict):
        return {k: _finite(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_finite(v) for v in o]
    return o


def dump2(name, data):
    """The format of gen_config_model.py."""
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
        json.dump(_finite(data), fh, indent=2, ensure_ascii=False, allow_nan=False)
        fh.write("\n")
    print("wrote", name)


def dump1(name, data):
    """The format of make_vectors.py and dump_tool_metadata.py."""
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print("wrote", name)


# ------------------------------------------------------------ gen_config_model.py, section 4


def resolve_paths_cases():
    os.environ["HOME"] = "/home/tester"
    base = "/base/proj"
    p = C.Project()
    p.inputs.topology = "top.parm7"
    p.inputs.trajectories = ["a.nc", "../up/b.nc", "/abs//x/../c.nc", "~/d.nc", "", "./e.nc", "sub//f.nc"]
    p.inputs.external = [C.ExternalSeries("data/rmsd.dat"), C.ExternalSeries(""), C.ExternalSeries("~")]
    p.cleaning.reference_file = "ref.pdb"
    p.cleaning.output_dir = "cleaned/"
    p.cleaning.cpptraj_path = "bin/cpptraj"
    p.output.out_dir = "MDATK_output"
    p.output.cache_dir = ""
    p.density.map_path = "~nosuchuser_xyz/map.mrc"
    p.workflow.out_dir = "wf/./x/.."
    p.output.video.filename = "videos/v.mp4"   # not resolved
    p.output.video.ffmpeg_path = "bin/ffmpeg"  # not resolved
    p.cleaning.conda_exe = "bin/conda"         # not resolved
    p.workflow.nodes = [
        C.PlanNode("n1", "input-file", 0, settings={"path": "maps/map.mrc", "kind": "map"},
                   files={"structure": " a.pdb ; ; /abs/b.pdb;~/c.pdb ", "empty": "", "blank": " ; "}),
        C.PlanNode("n2", "input-file", 0, settings={"path": "", "kind": "pdb"}),
        C.PlanNode("n3", "input-file", 0, settings={"path": 5, "kind": "pdb"}),
        C.PlanNode("n4", "amber-build", 1, settings={"path": "not/a/path/setting.pdb"},
                   files={"structure": "inputs/model.pdb"}),
        C.PlanNode("n5", "no-such-tool", 1, settings={"path": "x.pdb"}, files={"x": "y.pdb"}),
    ]
    p.workflow.scratch = [C.PlanNode("wizard-input-file", "input-file", 1, settings={"path": "s.pdb"},
                                     files={"p": "q.pdb"})]
    before = p.to_dict()
    p.resolve_paths(base)
    after = p.to_dict()

    def diff(a, b, prefix=""):
        out = []
        if isinstance(a, dict) and isinstance(b, dict):
            for k in a:
                out += diff(a[k], b.get(k), f"{prefix}.{k}" if prefix else k)
        elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
            for i, (x, y) in enumerate(zip(a, b)):
                out += diff(x, y, f"{prefix}[{i}]")
        elif a != b:
            out.append({"path": prefix, "before": a, "after": b})
        return out

    cpp = []
    for cp in ["cpptraj", "bin/cpptraj", "/opt/amber/bin/cpptraj", "~/amber/bin/cpptraj", "~", "", "./cpptraj"]:
        q = C.Project()
        q.cleaning.cpptraj_path = cp
        q.resolve_paths(base)
        cpp.append({"before": cp, "after": q.cleaning.cpptraj_path})
    fixes = []
    for raw in ["", "a", "a/", "./a", "a/./b", "a//b", "a/../b", "../a", "../../../../x", "/abs/../x", "//net/x",
                "///x", "~", "~/", "~/x/../y", "~root/x", "~nosuchuser_xyz/x", "a\\b", " a ", "caf\u00e9.pdb", "."]:
        q = C.Project()
        q.inputs.topology = raw
        q.resolve_paths(base)
        fixes.append({"before": raw, "after": q.inputs.topology})
    fixes_rel_base = []
    for raw in ["a", "../a"]:
        q = C.Project()
        q.inputs.topology = raw
        q.resolve_paths("relbase/sub")
        fixes_rel_base.append({"base_dir": "relbase/sub", "before": raw, "after": q.inputs.topology})
    dump2("resolve_paths_cases.json", {
        "how": "HOME=/home/tester (root's home from the passwd database: " + os.path.expanduser("~root") + "); "
               "base_dir = " + base + "; Project.resolve_paths(base_dir)",
        "base_dir": base,
        "home": "/home/tester",
        "project_before": before,
        "changes": diff(before, after),
        "unchanged_fields_note": "every path not listed under changes is untouched, e.g. output.video.filename, "
                                 "output.video.ffmpeg_path, cleaning.conda_exe, settings of tools without "
                                 "path_settings, non-string settings",
        "cpptraj_path": cpp,
        "fix": fixes,
        "fix_relative_base": fixes_rel_base,
    })


# ------------------------------------------------------------ gen_config_model.py, section 5


def dirs_cases():
    tmp = tempfile.mkdtemp(prefix="mdatk-dirs-")
    old = os.getcwd()
    os.chdir(tmp)
    cwd = os.getcwd()
    dirs = []
    for out_dir, cache_dir, wf_dir, name in [
        ("MDATK_output", "", "", "project"),
        ("", "", "", "project"),
        ("/abs/out/", "", "", "apt 1"),
        ("/abs/out/../o2", "cache/./x", "", "  ..My Plan!!  "),
        ("rel//out", "/abs/cache/", "wf", ""),
        ("rel", "", "/abs/wf/", "caf\u00e9-\u00fc.v1_"),
        ("rel", "", "", "---"),
        ("rel", "", "", "a..b"),
    ]:
        q = C.Project()
        q.output.out_dir, q.output.cache_dir, q.workflow.out_dir, q.workflow.project_name = (
            out_dir, cache_dir, wf_dir, name)

        def sub(s):
            return s.replace(cwd, "{cwd}")
        dirs.append({"output.out_dir": out_dir, "output.cache_dir": cache_dir, "workflow.out_dir": wf_dir,
                     "workflow.project_name": name, "out_dir": sub(q.out_dir), "cache_dir": sub(q.cache_dir),
                     "workflow_dir": sub(q.workflow_dir)})
    os.chdir(old)
    shutil.rmtree(tmp)
    san = [{"name": n, "sanitized": sanitize_project(n)} for n in
           ["project", "", "   ", "apt 1", "  ..My Plan!!  ", "caf\u00e9", "a/b\\c", "-._x._-", "x\u00a0y",
            "\t a\tb \n", "\u3000z\u3000", "A_b.c-d", "\u00e9\u00e9", "..", "a  b", "a\u200bb"]]
    groups_case = []
    q = C.Project()
    for name in ["Aptamer", "Target", "aptamer", "Nope"]:
        try:
            g = q.group(name)
            groups_case.append({"groups": q.group_names(), "name": name, "found": C._to_plain(g)})
        except KeyError as exc:
            groups_case.append({"groups": q.group_names(), "name": name, "error": {"type": "KeyError",
                                                                                 "args0": exc.args[0],
                                                                                 "str": str(exc)}})
    q.groups = []
    try:
        q.group("X")
    except KeyError as exc:
        groups_case.append({"groups": [], "name": "X", "error": {"type": "KeyError", "args0": exc.args[0],
                                                                 "str": str(exc)}})
    G = C.ResidueGroup
    q.groups = [G("A"), G("A", "2")]
    groups_case.append({"groups": q.group_names(), "name": "A", "found": C._to_plain(q.group("A"))})
    dump2("dirs_cases.json", {
        "how": "working directory = {cwd}; Project() with the four fields set; the out_dir / cache_dir / "
               "workflow_dir properties; sanitize_project(); Project.group()",
        "dirs": dirs, "sanitize_project": san, "group": groups_case,
    })


# ------------------------------------------------------------ gen_config_model.py, section 6


HASH_INPUTS = [
    ("empty_dict", {}),
    ("simple", {"b": 1, "a": [1, 2.5, "x", True, None]}),
    ("none_values_in_dicts_dropped", {"a": None, "b": {"c": None, "d": 0}}),
    ("floats", {"f": [0.1, 1.0, 1e-7, 1e16, 1e15, 123456789.123, -0.0, 2.5e-5, 1e22, 0.0001]}),
    ("nonfinite", {"f": [float("nan"), float("inf"), float("-inf")]}),
    ("unicode_escaped", {"s": "caf\u00e9 \u00c5 \u2013 \U0001F600 \x00 \x1f \x7f \" \\ / \n\t\r\b\f"}),
    ("key_order_codepoints", {"b": 1, "B": 2, "a": 3, "_": 4, "\u00e9": 5, "z": 6, "10": 7, "9": 8, "\uff01": 9,
                              "\U0001F600": 10}),
    ("tuple_to_list", {"t": (1, 2)}),
    ("int_keys_to_str", {1: "a", 2: "b"}),
    ("datetime_default_str", {"t": dt.datetime(2026, 10, 15, 14, 32, 5), "d": dt.date(2026, 1, 2),
                              "u": dt.datetime(2026, 10, 15, 14, 32, 5, 250000, tzinfo=dt.timezone.utc)}),
    ("big_int", {"i": 2 ** 63, "j": -2 ** 63, "k": 12345678901234567890123}),
    ("bools", {"t": True, "f": False}),
    ("nested_lists", [[1, [2, [3]]], []]),
    ("plain_string", "abc"),
    ("plain_number", 3),
    ("group", C.ResidueGroup("Aptamer", "1-12")),
    ("default_project", C.Project()),
    ("default_project_dict", C.Project().to_dict()),
    ("cleaning_params_vars", {k: v for k, v in vars(C.CleaningParams()).items()
                              if k not in ("enabled", "cpptraj_path", "output_dir", "conda_env", "conda_exe")}),
]


def stable_hash_cases():
    def blob_of(obj):
        return json.dumps(C._to_plain(obj), sort_keys=True, separators=(",", ":"), default=str)

    hcases = []
    for name, obj in HASH_INPUTS:
        blob = blob_of(obj)
        hcases.append({"name": name,
                       "input_repr": repr(obj) if not dataclasses.is_dataclass(obj) else f"{type(obj).__name__}(...)",
                       "blob": blob, "sha256": hashlib.sha256(blob.encode("utf-8")).hexdigest(),
                       "stable_hash_16": C.stable_hash(obj), "stable_hash_8": C.stable_hash(obj, 8),
                       "stable_hash_64": C.stable_hash(obj, 64), "stable_hash_100": C.stable_hash(obj, 100),
                       "stable_hash_0": C.stable_hash(obj, 0)})
    fd, fpath = tempfile.mkstemp(prefix="fp-")
    os.write(fd, b"hello")
    os.close(fd)
    os.utime(fpath, ns=(1_700_000_000_123_456_789, 1_700_000_000_123_456_789))
    fpr = C.file_fingerprint(fpath)
    fingerprint = {"note": "file_fingerprint of a 5-byte file whose mtime was set with os.utime(ns=...)",
                   "mtime_set_ns": 1_700_000_000_123_456_789,
                   "result": {**fpr, "path": "{tmpfile}"},
                   "key_order": list(fpr.keys())}
    os.unlink(fpath)
    try:
        C.file_fingerprint(fpath)
    except Exception as exc:  # noqa: BLE001
        fingerprint["missing_file_error"] = {"type": type(exc).__name__, "errno": getattr(exc, "errno", None),
                                             "message": str(exc).replace(fpath, "{tmpfile}")}
    dump2("stable_hash_cases.json", {
        "how": "blob = json.dumps(_to_plain(obj), sort_keys=True, separators=(',', ':'), default=str); "
               "stable_hash(obj, n) = sha256(blob.encode('utf-8')).hexdigest()[:n]",
        "cases": hcases, "file_fingerprint": fingerprint,
    })


def tagged(o):
    """A core::Value-shaped copy of o as strict JSON (see the module docstring)."""
    if isinstance(o, float) and not math.isfinite(o):
        return {"$float": "nan" if math.isnan(o) else ("inf" if o > 0 else "-inf")}
    if isinstance(o, dt.datetime):
        return {"$datetime": o.isoformat()}
    if isinstance(o, dt.date):
        return {"$date": o.isoformat()}
    if isinstance(o, int) and not isinstance(o, bool) and not -2 ** 63 <= o < 2 ** 63:
        return {"$int": str(o)}
    if isinstance(o, dict):
        return {str(k): tagged(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [tagged(v) for v in o]
    return o


def stable_hash_inputs():
    cases = []
    for name, obj in HASH_INPUTS:
        is_dc = dataclasses.is_dataclass(obj)
        cases.append({"name": name, "dataclass": is_dc, "input": tagged(C._to_plain(obj) if is_dc else obj)})
    dump1("stable_hash_inputs.json", {
        "how": "the inputs of stable_hash_cases.json (same order and names) as tagged JSON: "
               "{\"$float\": ...}, {\"$datetime\": isoformat}, {\"$date\": isoformat}, {\"$int\": decimal} "
               "beyond 64 bits; tuples as lists, keys as str(); a dataclass as its _to_plain() form (None fields "
               "kept), marked dataclass=true",
        "cases": cases,
    })


# ------------------------------------------------------------ gen_config_model.py, section 8


def python_repr_cases():
    rtoml = ('s1 = "far"\ns2 = "it\'s"\ns3 = "it\'s \\"q\\""\ns4 = "tab\\tnl\\nbs\\\\ bell\\u0007 del\\u007f"\n'
             's5 = "\u00e9 \u00c5 \u2013 \u00a0 \u200b \U0001F600 \ud7ff"\ni1 = 1\ni2 = -0\ni3 = 12345678901234567890\n'
             'f1 = 1.0\nf2 = 0.1\nf3 = 1e-7\nf4 = 1e16\nf5 = 1e15\nf6 = nan\nf7 = -inf\nf8 = -0.0\nf9 = 2.5e-5\n'
             'b1 = true\nb2 = false\na1 = []\na2 = [1, "x", 2.5, true, [], {}]\nt1 = {}\n'
             't2 = { a = 1, "b c" = "d" }\n'
             'd1 = 2026-10-15T14:32:05\nd2 = 2026-10-15T14:32:05.250\nd3 = 2026-10-15T14:32:05Z\n'
             'd4 = 2026-10-15T14:32:05.000001-05:30\nd5 = 2026-01-02\nd6 = 07:08:09\nd7 = 07:08:00\n'
             'd8 = 07:08:09.5\nd9 = 00:00:00\nd10 = 2026-10-15T14:32:00\nd11 = 2026-10-15T00:00:00+00:00\n')
    rdata = tomlio.loads(rtoml)
    rcases = [{"key": k, "repr": repr(v), "str": str(v), "type_name": type(v).__name__} for k, v in rdata.items()]
    dump2("python_repr_cases.json", {
        "how": "values = tomllib.loads(input_toml); for each: repr(v) (the {value!r} of error messages), str(v) "
               "(what a str field stores), type(v).__name__ (the 'got <type>' of 'expected an object')",
        "input_toml": rtoml, "cases": rcases,
    })


# ------------------------------------------------------------ make_vectors.py: round(x, 2)


def round2_vectors():
    xs = [0.0, 0.001, 0.004, 0.005, 0.0050000001, 0.015, 0.025, 0.035, 0.045, 0.125, 0.375, 0.625, 0.875,
          1.005, 1.015, 2.675, 1.23456, 2.555000000000007, 0.13500000000000512, 7.0, 12345.678, 1e-9,
          1e17, 2.0 ** 53, 0.995, 9.995, 99.995, 1234567.125, 0.1 + 0.2, 3600.0, 86399.999]
    rng = random.Random(7)
    xs += [rng.uniform(0, 10) for _ in range(200)] + [rng.uniform(0, 100000) for _ in range(100)]
    xs += [k / 1000 for k in range(0, 2001, 5)]
    out = [{"hex": x.hex(), "repr_in": repr(x), "round2": repr(round(float(x), 2))} for x in xs]
    dump1("round2_vectors.json", out)


# ------------------------------------------------------------ dump_tool_metadata.py: naming


def naming_cases():
    names = ["", "   ", "my apt 1", "apt1", "  ..My Plan!!  ", "caf\xe9", "a/b\\c", "--x--", "_a_", ".hidden",
             "a__b", "a..b", "a b\tc\nd", "x" * 3, "\xdcn\xefc\xf8d\xe9", "\u9805\u76ee", "a-._b", "\xa0nbsp\xa0",
             None]
    tags = ["", "eq_free", "a b/c", "amber-build", "--", "_", "stage 1", "min.1", "\xdcn\xef", None, "eq-free"]
    out = {
        "sanitize_project": [[n, naming.sanitize_project(n)] for n in names],
        "sanitize_tag": [[t, naming.sanitize_tag(t)] for t in tags],
        "step_label": [[k, naming.step_label(k)] for k in (0, 1, 9, 10, 99, 100, 1234, -1)],
        "step_folder": [[k, t, naming.step_folder(k, t)] for k, t in ((3, "amber-build"), (12, "eq_free"), (0, ""))],
        "file_name": [[a, b, c, d, e, naming.file_name(a, b, c, d, e)] for a, b, c, d, e in (
            ("apt1", 3, "tleap", "20261015-1432", "parm7"), ("apt1", 3, "tleap", "", "parm7"),
            ("my apt", 1, "eq_free", "20261015-1432", ".rst7"), ("", 0, "plan", "", "sh"), ("a", 2, "b", "", ""),
            ("a", 2, "b", "s", "..x"))],
        "parse_name": [[n, naming.parse_name(n)] for n in (
            "my_apt_03_eq-free_20261015-1432.rst7", "apt1_03_tleap.parm7", "apt1_03_tleap", "apt1_3_tleap.parm7",
            "apt1_003_x_20261015-1432", "x_01_a_b.pdb", "apt1_03_tleap.tar.gz", "_01_a", "a_01_b_20261015-14321.x")],
        "new_stamp_example": naming.new_stamp(dt.datetime(2026, 10, 15, 14, 32, 59)),
    }
    dump1("naming_cases.json", out)


# ------------------------------------------------------------ new: posixpath


PATHS = ["", ".", "..", "...", "/", "//", "///", "////", "a", "a/", "a//", "/a", "//a", "///a", "////a",
         "a/b", "a//b", "a/./b", "a/../b", "../a", "../../a", "a/..", "a/../..", "a/b/../../..", "/..", "/../a",
         "//..", "//../a", "///../a", "./", "./.", "./a/", "/.", "/./a", "~", "~/", "~/x", "~x/y", ".hidden",
         "a/.hidden", "..hidden", "a.b.c", "a/.b", "a.b/c", "a/b.", "a..", "..a", "...a", "/a/b.c/", "a.b/",
         "foo.tar.gz", "dir.d/file", "dir.d/.file", "dir.d/..file", "a\\b", "a\\b.c\\d", " a ", " /a", "a /",
         "caf\u00e9/\u00fc.pdb", "\U0001F600.x", "a/b/", "a//b//c//", "/a/./b/../c/", "x/../../y/./z/..", "a:b",
         "c:/x", "./..", "../..", "../../", "/a/b/c/../../../../d"]
HOMES = ["/home/tester", "/home/u/", "/", "", "//", "/a//"]


def posixpath_cases():
    old = os.getcwd()
    cwd_cases = []
    for cwd in ["/"]:
        os.chdir(cwd)
        cwd_cases += [[p, os.getcwd(), posixpath.abspath(p)] for p in PATHS]
    os.chdir(old)
    home_before = os.environ.get("HOME")
    expand = []
    for home in HOMES:
        os.environ["HOME"] = home
        expand += [[home, p, posixpath.expanduser(p)] for p in
                   ["~", "~/", "~/x", "~//x", "~/x/../y", "x/~", "/~", "", "~nosuchuser_xyz", "~nosuchuser_xyz/x",
                    "~root", "~root/x", "~root//x"]]
    if home_before is None:
        del os.environ["HOME"]
    else:
        os.environ["HOME"] = home_before
    bases = ["", "/", "a", "a/", "/a", "//", "a//"]
    dump1("posixpath_cases.json", {
        "how": "posixpath (Python " + sys.version.split()[0] + ") of each path; abspath with the working directory "
               "given; expanduser with HOME set to each value (~root: root's passwd entry, here "
               + posixpath.expanduser("~root") + "; ~nosuchuser_xyz: no such user)",
        "normpath": [[p, posixpath.normpath(p)] for p in PATHS],
        "isabs": [[p, posixpath.isabs(p)] for p in PATHS],
        "dirname": [[p, posixpath.dirname(p)] for p in PATHS],
        "basename": [[p, posixpath.basename(p)] for p in PATHS],
        "splitext": [[p, list(posixpath.splitext(p))] for p in PATHS],
        "join": [[b, p, posixpath.join(b, p)] for b in bases for p in PATHS],
        "abspath": cwd_cases + [[p, cwd, posixpath.normpath(posixpath.join(cwd, p))]
                                for cwd in ("/base/proj", "//net/share", "/x/") for p in PATHS],
        "expanduser": expand,
    })


# ------------------------------------------------------------ new: repr() of strings


def string_repr_cases():
    """Latin-1, both sides of every change of str.isprintable() (so every range of the table) and texts
    that mix quotes, escapes and undecodable bytes."""
    picks = set(range(0, 0x100)) | {0x10FFFF}
    for code in range(1, 0x110000):
        if chr(code).isprintable() != chr(code - 1).isprintable():
            picks |= {code - 1, code}
    picks = sorted(c for c in picks if not 0xD800 <= c <= 0xDFFF)
    texts = [chr(c).encode("utf-8") for c in picks]
    texts += [s.encode("utf-8") for s in [
        "", "plain", "it's", 'say "hi"', "it's \"both\"", "back\\slash", "tab\tnl\ncr\r", "\x00\x1f\x7f\x80\x9f\xa0",
        "\u00ad soft hyphen", "\u200b\u200e\u2028\u2029\u3000", "\ue000 private", "\U000e0001 tag", "\U0010ffff",
        "caf\u00e9 \u00c5 \u2013 \U0001F600", "\ufeff BOM", "\ufffd replacement", "'", '"', "'\"", "\\'"]]
    # bytes that are not UTF-8 reach MDBridge raw; Python holds them as lone surrogates (surrogateescape)
    texts += [b"caf\xe9.pdb", b"\xff", b"\xe2\x82", b"\xe2\x82A", b"\xed\xa0\x80", b"\xc0\xaf", b"\xf4\x90\x80\x80",
              b"ok\x80'\""]
    cases = []
    for raw in texts:
        s = raw.decode("utf-8", "surrogateescape")
        cases.append([raw.hex(), repr(s), json.dumps(s)])
    dump1("string_repr_cases.json", {
        "how": "[utf8 hex, repr(s), json.dumps(s)] with s = bytes.fromhex(hex).decode('utf-8', 'surrogateescape'), "
               "Python " + sys.version.split()[0] + " (Unicode " + unicodedata.unidata_version + ")",
        "cases": cases,
    })


# ------------------------------------------------------------ new: json.loads


JSON_TEXTS = [
    "", " ", "x", "[", "[1,", "[1,]", '{"a":1,}', "{", '{"a"', '{"a" 1}', '{"a":}', "{1:2}", "[1 2]",
    '{"a":1 "b":2}', '"abc', '"a\\x"', '"a\\u12"', '"a\\u12g4"', '"\x01"', "01", "1.", "-", "-a", "1e", "1e+",
    ".5", "+1", "NaN", "-Infinity", "Infinity", "nan", "tru", "true false", "[] x", '"\\ud800"', '"\\udc00\\ud800"',
    '"\\ud83d\\ude00"', '"\\ud83d\\ude00', '"\\ud83d\\u0041"', '"\\ud83d\\uzzzz"', '"\\ud83d\\"', '"\\ud83dx"',
    "1\n2", "[\n1,\n]", "\ufeff1", '{"a":1,"a":2,"b":3}', '{"b":1,"a":2,"b":{"c":3}}', "-0", "-0.0", "1E5",
    "1e-400", "-1e-400", "1e400", "-1e400", "123456789012345678901234567890", "9223372036854775807",
    "9223372036854775808", "-9223372036854775808", "-9223372036854775809", '"\\/"', '"\u00e9"',
    '  {"a" : [ 1 , 2 ] }  ', '"\\', '"\\u', "[1,,2]", "{,}", '{"a":1,,}', "[,]", '"a\nb"', '"\t"', '"\x7f"',
    "-Infinityx", "Infinit", "nul", "\u00e9", '["\u00e9", x]', '["\U0001F600",\n "\u00e9" x]',
    '{"\u00e9\u00e9": 1, "k": [true, false, null]}', "0.1", "5e-324", "2.4703282292062328e-324",
    "2.4703282292062327e-324", "1.7976931348623157e308", "1.7976931348623158e308", "1.7976931348623159e308",
    "2.2250738585072011e-308", "0.1000000000000000055511151231257827021181583404541015625", "1e22", "1e23",
    "-1.5e-07", "1.0", "100", "[[[[]]]]", '{"a":{"b":{"c":[{}]}}}', '"\\b\\f\\n\\r\\t\\"\\\\"', '"\\u0000x"',
    '"\\u00e9\\u00C9"', "\t\n\r [1]\r\n", " \u00a0[1]", "[1]\u00a0", "null", "true", "false", '"\\uD83D\\uDE00"',
    '{"a":1}{', "[1]]", "{}}", '{"a"::1}', '{"a":1,"b"}', "[-]", "[-01]", "[0.]", "[1.e5]", "[1e5.]", "[.e]",
    "0e0", "0E-0", "-0e+00", "1.5E+3", "[1,2,3,", '"unterminated \\u00e9', '{"\\ud800": "\\udfff"}',
]


def json_loads_cases():
    rng = random.Random(13)
    alphabet = ['{', '}', '[', ']', ',', ':', '"', '\\', 'u', '0', '1', '9', '-', '+', '.', 'e', 'E', ' ', '\n',
                't', 'r', 'n', 'N', 'a', 'I', '\u00e9', 'x', 'd', '8', 'D']
    seeds = [t for t in JSON_TEXTS if t]
    texts = list(JSON_TEXTS)
    while len(texts) < len(JSON_TEXTS) + 300:  # mutations of the hand-written texts
        t = list(rng.choice(seeds))
        for _ in range(rng.randint(1, 3)):
            k = rng.randrange(len(t) + 1)
            op = rng.randrange(3)
            if op == 0:
                t.insert(k, rng.choice(alphabet))
            elif op == 1 and k < len(t):
                del t[k]
            elif k < len(t):
                t[k] = rng.choice(alphabet)
        texts.append("".join(t))

    texts += [random_document(rng) for _ in range(150)]

    def big_int(v):
        if isinstance(v, int) and not isinstance(v, bool):
            return not -2 ** 63 <= v < 2 ** 63
        if isinstance(v, list):
            return any(big_int(x) for x in v)
        if isinstance(v, dict):
            return any(big_int(x) for x in v.values())
        return False

    def replaced(v):  # what MDBridge holds: U+FFFD for a lone surrogate
        if isinstance(v, str):
            return "".join("\ufffd" if 0xD800 <= ord(c) <= 0xDFFF else c for c in v)
        if isinstance(v, list):
            return [replaced(x) for x in v]
        if isinstance(v, dict):
            out = {}
            for k, x in v.items():
                out[replaced(k)] = replaced(x)
            return out
        return v

    cases = []
    for text in texts:
        case = {"text": text}
        try:
            v = json.loads(text)
            case["repr"] = repr(v)
            # MDBridge's documented deviations (core/json.hpp)
            if big_int(v):
                case["lossy"] = "int64"
            elif repr(replaced(v)) != repr(v):
                case["lossy"] = "surrogate"
                case["repr_mdbridge"] = repr(replaced(v))
        except json.JSONDecodeError as exc:
            case["error"] = str(exc)
        except RecursionError as exc:
            case["error"] = str(exc)
        cases.append(case)
    dump1("json_loads_cases.json", {
        "how": "json.loads(text) (Python " + sys.version.split()[0] + "): repr() of the result, or str() of the "
               "JSONDecodeError; lossy: \"int64\" (an integer outside 64 bits) or \"surrogate\" (a lone surrogate, "
               "repr_mdbridge with U+FFFD in its place)",
        "cases": cases,
    })


def random_float_text(rng):
    """A float in one of several spellings (also many digits, and integral ones that json reads as int)."""
    bits = struct.unpack("<d", struct.pack("<Q", rng.getrandbits(64)))[0]
    x = rng.choice([rng.uniform(-1e3, 1e3), rng.uniform(0, 1e-300), rng.uniform(0, 1e300), bits,
                    rng.randint(-10 ** 6, 10 ** 6) / 8])
    if not math.isfinite(x):
        x = 1.5
    return rng.choice(["%r", "%.17e", "%.30g", "%.3E"]) % x


def random_document(rng):
    """A random JSON text written in one of several layouts."""
    def value(d):
        kind = rng.randrange(9 if d < 4 else 6)
        if kind == 0:
            return rng.choice([None, True, False])
        if kind == 1:
            return rng.choice([0, -1, 7, rng.randint(-2 ** 63, 2 ** 63 - 1), rng.randint(-1000, 1000)])
        if kind == 2:
            return RawFloat(random_float_text(rng))
        if kind in (3, 4, 5):
            return "".join(rng.choice(["a", "\"", "\\", "/", "\n", "\t", "\x01", "\u00e9", "\u2028",
                                       "\U0001f600", " ", "z", "\x7f", "\ufffd"]) for _ in range(rng.randrange(6)))
        if kind in (6, 7):
            return [value(d + 1) for _ in range(rng.randrange(4))]
        return {rng.choice(["a", "b", "c", "\u00e9", "k k", ""]): value(d + 1) for _ in range(rng.randrange(4))}

    return dumps_raw(value(0), rng.choice([None, 0, 2]), rng.random() < 0.5, rng.random() < 0.5)


class RawFloat(str):
    """A float written with exactly this text."""


def dumps_raw(v, indent, ascii_only, spaced, level=0):
    pad = "" if indent is None else "\n" + " " * (indent * (level + 1))
    end = "" if indent is None else "\n" + " " * (indent * level)
    sep = ", " if spaced and indent is None else ","
    if isinstance(v, RawFloat):
        return str(v)
    if isinstance(v, list):
        if not v:
            return "[]"
        return "[" + sep.join(pad + dumps_raw(x, indent, ascii_only, spaced, level + 1) for x in v) + end + "]"
    if isinstance(v, dict):
        if not v:
            return "{}"
        colon = ": " if spaced else ":"
        return "{" + sep.join(pad + json.dumps(k, ensure_ascii=ascii_only) + colon +
                              dumps_raw(x, indent, ascii_only, spaced, level + 1) for k, x in v.items()) + end + "}"
    return json.dumps(v, ensure_ascii=ascii_only)


# ------------------------------------------------------------ new: round(x, n), sha256, naming


def round_cases():
    xs = [0.0, -0.0, 0.5, 1.5, 2.5, -2.5, 0.125, 0.375, 2.675, 1.005, 5.0, 15.0, 25.0, 149.5, 150.0, 250.0,
          -150.0, 0.7, 49.999, 50.0, 51.0, 12345.678, -12345.678, 123.456, -123.456, 5e22, 1e22, 1e300, -1e300,
          5e-324, 1e-300, 2.0 ** 53 + 2, 1.7976931348623157e308, 0.1 + 0.2, 999999.5, 9.5, 0.05, 0.15, 0.25]
    digits = [-309, -308, -307, -300, -23, -22, -10, -5, -3, -2, -1, 0, 1, 2, 3, 5, 10, 15, 17, 20, 100, 300,
              323, 324]
    out = []
    for x in xs:
        for n in digits:
            try:
                out.append([x.hex(), n, repr(round(x, n))])
            except OverflowError as exc:
                out.append([x.hex(), n, "OverflowError: " + str(exc)])
    dump1("round_cases.json", {"how": "[float.hex(x), n, repr(round(x, n)) or the OverflowError]", "cases": out})


def sha256_cases():
    texts = [("a" * n).encode() for n in range(201)]
    texts += [b"abc", b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
              "Mg²⁺ 😀".encode(), bytes(range(256)) * 5]
    dump1("sha256_cases.json", {"how": "[bytes as hex, hashlib.sha256(bytes).hexdigest()]",
                                "cases": [[t.hex(), hashlib.sha256(t).hexdigest()] for t in texts]})


def naming_extra_cases():
    names = ["apt1_03_tleap\n", "apt1_03_tleap\n\n", "a\nb_03_x", "apt_\u0660\u0663_x", "a_\u0663\u0663_t_\u0661" * 1,
             "a_01_b_02_c", "a_01_b_02_c.d", "a_01_b_20261015-1432_20261015-1432", "\u00e9_12_x.y", "a_01_-",
             "a_01_b.", "a_01_b._", "a_1234567_b_20261015-1432.tar", "a_01_b_2026101\u0665-1432.x", "_01_a_01_b",
             "a__01_b", "a_01__b", "a_01_b_c_02_d.e", "x_99_y.Z9", "x_99_y.\u00e9", "\n_01_a", "a_01_b\r"]
    texts = ["", " \u00a0", "\u2003x\u2003", "a\u0085", "\x1c_x_\x1f", "a\udce9b", "\U0001F600tag", "1.2.3",
             "...", "_-_", "A-Z_a.z", "\u0663", "tab\tsep"]
    stamps = [dt.datetime(999, 1, 2, 3, 4), dt.datetime(5, 12, 31, 23, 59, 59), dt.datetime(9999, 1, 1, 0, 0),
              dt.datetime(2026, 10, 15, 14, 32, 59, 999999)]
    dump1("naming_extra_cases.json", {
        "how": "mdatk.workflows.naming on further edge cases (Unicode digits, newlines, greedy matching, "
               "unpadded %Y of strftime on Linux); texts as UTF-8 hex (surrogateescape)",
        "parse_name": [[os.fsencode(n).hex(), naming.parse_name(n)] for n in names],
        "sanitize_project": [[os.fsencode(t).hex(), naming.sanitize_project(t)] for t in texts],
        "sanitize_tag": [[os.fsencode(t).hex(), naming.sanitize_tag(t)] for t in texts],
        "new_stamp": [[d.isoformat(), naming.new_stamp(d)] for d in stamps],
    })


def main() -> int:
    resolve_paths_cases()
    dirs_cases()
    stable_hash_cases()
    stable_hash_inputs()
    python_repr_cases()
    round2_vectors()
    naming_cases()
    posixpath_cases()
    string_repr_cases()
    json_loads_cases()
    round_cases()
    sha256_cases()
    naming_extra_cases()
    return 0


if __name__ == "__main__":
    sys.exit(main())
