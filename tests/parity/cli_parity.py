#!/usr/bin/env python3
"""The milestone-2 commands of the command line, MDATK v0.0.4.1 against MDBridge.

    python3 -I -B tests/parity/cli_parity.py --mdatk MDATK_CHECKOUT --mdbridge EXE [--only NAME ...]
    python3 -I -B tests/parity/cli_parity.py --mdatk MDATK_CHECKOUT --write tests/data/cli/cli_cases.json
    python3 -I -B tests/parity/cli_parity.py --list

Every case (CASES below) runs `mdatk ARGS` and `EXE ARGS` in two fresh work folders holding the
same inputs (copies of files of tests/data, and project files written by the case), then compares
the exit status, stdout, stderr and every file the command left in the folder (byte for byte).
`init` / `init-params`, `info`, `export-pdb` and `density info` are covered: their results, files
(project files, PDB files), messages and errors.

MDATK runs in-process in one child interpreter (this script with --child), from MDATK_CHECKOUT
(nothing is written into it), with MDBridge's independent data tables patched into MDAnalysis
(independent_tables.py, tests/data/README.md): MDBridge deliberately uses its own residue,
element and mass tables, so the reference is MDAnalysis' behaviour with those tables. The child
needs Python 3.13 with MDATK's dependencies at the reference versions (MDAnalysis 2.11.0;
requirements-regenerate.txt); EXE is the `mdbridge` program, the AppImage, or the test build's
`mdbridge_commands_cli` (the same commands without the Qt front end).

Normalised on both sides before comparing (the only differences allowed):
  * the work folder's path -> {work} (each side runs in its own folder);
  * `saved = <time>` of project files -> {time}; run time stamps (2026-10-15_143205) -> {stamp};
  * MDATK's stderr loses the line MDAnalysis logs when netCDF4 is missing (an environment notice,
    not MDATK's output; MDBridge reads NetCDF itself);
  * MDAnalysis' hidden XTC / TRR offset caches (.NAME_offsets.npz / .lock) are not files of the run;
  * MDATK messages that MDBridge words with its own name (NAMING) are adapted.
Cases marked `exit_only` compare the exit status alone: usage errors (argparse and CLI11 word them
differently; both exit with 2) and documented deviations (DEVIATIONS: map errors that
src/io/include/mdbridge/io/mrc.hpp words differently on purpose). Where MDATK crashes (an exception
its main() does not catch: a Python traceback, exit status 1), MDBridge must exit with 1 and print
"error: <message>"; stdout and files are compared as usual.

Help texts are not compared (argparse and CLI11 lay them out differently; MDBridge's `--help` lists
only the commands it has).

--write runs MDATK alone and writes the normalised results as the expected data of
tests/cli/test_commands.cpp (tests/data/cli/cli_cases.json); tests/parity/regenerate_all.py runs it.
Files up to 16 KiB are stored as text, larger ones as their SHA-256 and size. Standard library only
(the child imports MDATK). Exit status 1 on any difference.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
LIMIT = 16384

DEMO = {"demo_dry.parm7": "topology/demo_dry.parm7", "demo_dry.nc": "io/trajectories/demo_dry.nc"}
DEMO_RUN = ["--top", "demo_dry.parm7", "--traj", "demo_dry.nc"]
NO_STAMP = ["-s", "output.timestamp_files=false"]
MAPS = "io/maps/files"

# Maps whose error MDBridge words differently on purpose (src/io/include/mdbridge/io/mrc.hpp: a
# broken gzip member raises MapError "{path}: <gzip's message>", Python gzip.BadGzipFile without
# the path): exit status only. (Maps on which MDATK crashes are compared as described above.)
DEVIATIONS = {
    "gzip_bad_crc.map.gz", "gzip_bad_isize.map.gz", "gzip_bad_magic.map.gz", "gzip_bad_magic_quote.map.gz",
    "gzip_bad_method.map.gz", "gzip_one_byte_tail.map.gz",
}

PROJECT = """format = "mdatk-project"
notes = "relative paths, resolved against this file's folder"

[inputs]
topology = "demo_dry.parm7"
trajectories = ["demo_dry.nc"]
dt_ps = 5.0
time_unit = "ps"

[[groups]]
name = "Aptamer"
residues = "1-4, 21-24"
color = "#123456"

[[groups]]
name = "Target"
selection = "resid 5:20 and not name H*"

[output]
out_dir = "results"
timestamp_files = false

[output.pdb]
frames = "1, 3"
times = "100"
selection = "Aptamer"
prefix = "snap"
"""

GROUPS = json.dumps([
    {"name": "A", "residues": "1-3"},
    {"name": "bad name", "selection": "resname ZZZ"},
    {"name": "A", "residues": "20-22"},
    {"name": "Ions", "selection": "resname MG"},
    {"name": "Broken", "selection": "name ("},
    {"name": "Chain", "residues": "X:5-6"},
    {"name": "Empty", "residues": ""},
    {"name": "Ünï", "residues": "7"},
])


def case(name, argv, inputs=None, setup=None, exit_only=False):
    """inputs: work-folder path -> file of tests/data copied there; setup: path -> text written there."""
    return {"name": name, "argv": argv, "inputs": dict(inputs or {}), "setup": dict(setup or {}),
            "exit_only": exit_only}


def map_cases():
    out = []
    for name in sorted(os.listdir(os.path.join(DATA, MAPS))):
        out.append(case(f"density_info_{name}", ["density", "info", f"maps/{name}"],
                        {f"maps/{name}": f"{MAPS}/{name}"}, exit_only=name in DEVIATIONS))
    return out


CASES = [
    # ---- init / init-params
    case("init_default", ["init", "new.mdatk.toml"]),
    case("init_alias_options", ["init-params", "p.mdatk.toml", *DEMO_RUN, "--frames", "2:20:3", "-o", "out",
                                "--cache-dir", "cache", "--cpptraj", "/opt/amber/bin/cpptraj", "--conda-env",
                                "AmberTools26", "--clean", "-s", "analyses.contacts.cutoff=5", "-s",
                                "groups.0.residues=1-6", "-s", "notes=5", "-s", 'output.charts.formats=["png","svg"]',
                                "-s", "inputs.dt_ps=2.5", "-s", "groups.-1.name=Last", "-s",
                                'layout.panels.3.options={"mode": "mean"}', "-s", "visual.rotation.1=90",
                                "-s", "cleaning.enabled=false", "--no-clean", "--clean"], DEMO),
    case("init_relative", ["init", "rel.mdatk.toml", "--top", "sub/../demo_dry.parm7", "--traj", "a.nc", "b.nc",
                           "--frames", ":7", "-s", "inputs.trajectories=[\"x.nc\"]", "-s", "output.out_dir=o"]),
    case("init_from_project", ["init", "copy.mdatk.toml", "-p", "base.mdatk.toml", "--no-clean"], DEMO,
         {"base.mdatk.toml": PROJECT}),
    case("init_from_params_log", ["init", "replay.mdatk.toml", "--params", "MDATK_params_log.toml"],
         {"MDATK_params_log.toml": "projects/MDATK_params_log.toml"}),
    case("init_exists", ["init", "exists.mdatk.toml"], setup={"exists.mdatk.toml": "keep me\n"}),
    case("init_overwrite", ["init", "base.mdatk.toml", "--overwrite", "-s", "notes=rewritten"],
         setup={"base.mdatk.toml": PROJECT}),
    case("init_json", ["init", "preset.json"]),
    case("init_no_folder", ["init", "missing/new.mdatk.toml"]),
    case("init_set_no_equals", ["init", "x.mdatk.toml", "-s", "nonsense"]),
    case("init_set_empty_key", ["init", "x.mdatk.toml", "-s", " . =3"]),
    case("init_set_unknown", ["init", "x.mdatk.toml", "-s", "analyses.bogus=1"]),
    case("init_set_free_form", ["init", "x.mdatk.toml", "-s", "layout.panels.3.options.mode=x"]),
    case("init_set_bool", ["init", "x.mdatk.toml", "-s", "cleaning.enabled=1"]),
    case("init_set_index", ["init", "x.mdatk.toml", "-s", "groups.7.name=x"]),
    case("init_set_index_text", ["init", "x.mdatk.toml", "-s", "groups.first.name=x"]),
    case("init_set_index_dashes", ["init", "x.mdatk.toml", "-s", "groups.--1.name=x"]),
    case("init_set_type", ["init", "x.mdatk.toml", "-s", "analyses.contacts.cutoff=far"]),
    case("init_set_table", ["init", "x.mdatk.toml", "-s", 'groups={"name": "x"}']),
    case("init_set_property", ["init", "x.mdatk.toml", "-s", "out_dir=5"]),
    case("init_set_method", ["init", "x.mdatk.toml", "-s", "validate.x=5"]),
    case("init_set_section", ["init", "x.mdatk.toml", "-s", 'inputs={"topology": "t.parm7", "bogus": 1}']),
    case("init_frames_bad", ["init", "x.mdatk.toml", "--frames", "1:x"]),
    case("init_missing_output", ["init"], exit_only=True),
    case("init_unknown_option", ["init", "x.mdatk.toml", "--bogus"], exit_only=True),

    # ---- info
    case("info_demo", ["info", *DEMO_RUN], DEMO),
    case("info_no_topology", ["info"]),
    case("info_groups", ["info", *DEMO_RUN, "-s", f"groups={GROUPS}"], DEMO),
    case("info_project", ["info", "-p", "base.mdatk.toml"], DEMO, {"base.mdatk.toml": PROJECT}),
    case("info_two_trajectories", ["info", "--top", "demo_dry.parm7", "--traj", "demo_dry.nc", "demo_dry.xtc",
                                   "-s", "inputs.dt_ps=4"],
         {**DEMO, "demo_dry.xtc": "io/trajectories/files/demo_dry.xtc"}),
    case("info_dcd", ["info", "--top", "demo_dry.parm7", "--traj", "demo_dry.dcd", "--frames", "::0"],
         {**DEMO, "demo_dry.dcd": "io/trajectories/files/demo_dry.dcd"}),
    case("info_pdb", ["info", "--top", "elements.pdb", "-s", 'groups=[{"name": "All", "selection": "all"}]'],
         {"elements.pdb": "io/loader/files/elements.pdb"}),
    case("info_gro", ["info", "--top", "precise.gro"], {"precise.gro": "io/loader/files/precise.gro"}),
    case("info_missing_topology", ["info", "--top", "missing.parm7"]),
    case("info_missing_trajectory", ["info", "--top", "demo_dry.parm7", "--traj", "missing.nc"], DEMO),
    case("info_parm7_alone", ["info", "--top", "demo_dry.parm7"], DEMO),
    case("info_unknown_option", ["info", "--bogus"], exit_only=True),

    # ---- export-pdb
    case("export_demo", ["export-pdb", *DEMO_RUN, "--frame-list", "0, 10, 20-23:2, last", "--selection", "Aptamer",
                         *NO_STAMP], DEMO),
    case("export_times", ["export-pdb", *DEMO_RUN, "--times", "55ps, 0.1ns, 10ns", "--selection", "resname MG",
                          "-o", "pdbs", *NO_STAMP, "-s", "output.pdb.prefix=ion"], DEMO),
    case("export_stamp", ["export-pdb", *DEMO_RUN, "--frame-list", "3", "--selection", "Target"], DEMO),
    case("export_all_atoms", ["export-pdb", *DEMO_RUN, "--frame-list", "last", *NO_STAMP], DEMO),
    case("export_project", ["export-pdb", "-p", "base.mdatk.toml"], DEMO, {"base.mdatk.toml": PROJECT}),
    case("export_project_override", ["export-pdb", "-p", "base.mdatk.toml", "--frame-list", "", "--times",
                                     "0.01ns", "--selection", "Aptamer, Target", "-q"], DEMO,
         {"base.mdatk.toml": PROJECT}),
    case("export_verbose_layout", ["export-pdb", *DEMO_RUN, "--frame-list", "2", "-v", *NO_STAMP, "-s",
                                   "layout.panels.1.type=structure", "--selection", "resid 3"], DEMO),
    case("export_no_frames", ["export-pdb", *DEMO_RUN, *NO_STAMP], DEMO),
    case("export_invalid", ["export-pdb", "--top", "demo_dry.parm7", "--frame-list", "0"], DEMO),
    case("export_bad_selection", ["export-pdb", *DEMO_RUN, "--frame-list", "0", "--selection", "name ("], DEMO),
    case("export_empty_selection", ["export-pdb", *DEMO_RUN, "--frame-list", "0", "--selection", "resname ZZZ"],
         DEMO),
    case("export_bad_frames", ["export-pdb", *DEMO_RUN, "--frame-list", "0, 99"], DEMO),
    case("export_pdb_trajectory", ["export-pdb", "--top", "chains.pdb", "--traj", "chains.pdb", "--frame-list",
                                   "0", *NO_STAMP, "-s", 'groups=[{"name": "All", "selection": "all"}]',
                                   "-s", "analyses.contacts.enabled=false", "-s", "analyses.hbonds.enabled=false",
                                   "-s", "analyses.mg.enabled=false", "-s", "analyses.stacking.enabled=false",
                                   "-s", "analyses.energy.enabled=false"],
         {"chains.pdb": "io/pdb/export/files/chains.pdb"}),

    # ---- density info
    *map_cases(),
    case("density_missing", ["density", "info", "missing.mrc"]),
    case("density_folder", ["density", "info", "."]),
    case("density_no_action", ["density"], exit_only=True),
]

STAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}_\d{6}")
SAVED_RE = re.compile(r"^saved = \d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?$", re.M)
OFFSETS_RE = re.compile(r"^\..+_offsets\.(npz|lock)$")
NETCDF_NOTICE = "WARNING: netCDF4 is not available. Writing AMBER ncdf files will be slow.\n"
TRACEBACK = "Traceback (most recent call last):"
# MDATK's words that MDBridge's libraries replace by its own name (applied to MDATK's output)
NAMING = [("JSON presets are no longer supported; MDATK 0.0.4.1 uses TOML project files",
           "JSON presets are no longer supported; MDBridge uses TOML project files")]


# ---------------------------------------------------------------------------------------------
# work folders

def prepare(case_, folder):
    os.makedirs(folder)
    for dest, source in case_["inputs"].items():
        target = os.path.join(folder, dest)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        shutil.copyfile(os.path.join(DATA, source), target)
    for dest, text in case_["setup"].items():
        with open(os.path.join(folder, dest), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)


def normalise(text, folder):
    for path in {os.path.realpath(folder), os.path.abspath(folder)}:
        text = text.replace(path, "{work}")
    return STAMP_RE.sub("{stamp}", text)


def collect(case_, folder):
    """Files left in the folder (inputs included when the run changed them): relative path -> bytes."""
    out = {}
    for here, _dirs, names in os.walk(folder):
        for name in names:
            if OFFSETS_RE.match(name):
                continue
            full = os.path.join(here, name)
            rel = os.path.relpath(full, folder).replace(os.sep, "/")
            with open(full, "rb") as fh:
                data = fh.read()
            if rel in case_["inputs"]:
                continue  # copies of test data: never written by these commands
            if rel in case_["setup"] and data == case_["setup"][rel].encode("utf-8"):
                continue
            out[normalise(rel, folder)] = data
    return dict(sorted(out.items()))


def file_entry(data, folder):
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return {"sha256": hashlib.sha256(data).hexdigest(), "size": len(data)}
    text = SAVED_RE.sub("saved = {time}", normalise(text, folder))
    if len(text.encode("utf-8")) <= LIMIT:
        return {"text": text}
    raw = text.encode("utf-8")
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size": len(raw)}


def result(case_, folder, code, stdout, stderr):
    files = collect(case_, folder)
    texts = [stdout, stderr, *(data.decode("utf-8", "replace") for data in files.values())]
    return {"exit": code, "stdout": normalise(stdout, folder), "stderr": normalise(stderr, folder),
            "written": {rel: file_entry(data, folder) for rel, data in files.items()},
            # results that hold absolute paths or OS error texts, or read a project file's paths:
            # defined on Linux only (the C++ test skips them on Windows)
            "posix_only": any("{work}" in normalise(t, folder) for t in texts) or "[Errno" in stderr
            or any(a in ("-p", "--project", "--params") for a in case_["argv"])}


# ---------------------------------------------------------------------------------------------
# MDATK, in-process in a child interpreter

def child(mdatk, jobs_file, results_file):
    import logging
    import warnings

    warnings.filterwarnings("ignore")
    sys.path.insert(0, os.path.abspath(mdatk))
    sys.path.insert(0, HERE)
    import independent_tables

    independent_tables.apply()
    from mdatk import cli

    with open(jobs_file, encoding="utf-8") as fh:
        jobs = json.load(fh)
    results = []
    start = os.getcwd()
    for argv, folder in jobs:
        out, err = io.StringIO(), io.StringIO()
        root = logging.getLogger()
        for handler in list(root.handlers):  # basicConfig() in cli.main binds the current sys.stderr
            root.removeHandler(handler)
        os.chdir(folder)
        old = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = out, err
        try:
            code = cli.main(argv)
        except SystemExit as exc:  # argparse
            code = exc.code if isinstance(exc.code, int) else 1
        except BaseException:  # noqa: BLE001 - an uncaught exception: Python prints it, exit 1
            traceback.print_exc()
            code = 1
        finally:
            sys.stdout, sys.stderr = old
            os.chdir(start)
        results.append([code, out.getvalue(), err.getvalue().replace(NETCDF_NOTICE, "")])
    with open(results_file, "w", encoding="utf-8") as fh:
        json.dump(results, fh)


def run_mdatk(cases, root, mdatk, python):
    jobs = []
    for c in cases:
        folder = os.path.join(root, "mdatk", c["name"])
        prepare(c, folder)
        jobs.append([c["argv"], folder])
    jobs_file, results_file = os.path.join(root, "jobs.json"), os.path.join(root, "results.json")
    with open(jobs_file, "w", encoding="utf-8") as fh:
        json.dump(jobs, fh)
    subprocess.run([python, "-I", "-B", os.path.abspath(__file__), "--child", mdatk, jobs_file, results_file],
                   check=True)
    with open(results_file, encoding="utf-8") as fh:
        outputs = json.load(fh)
    results = {}
    for c, (_argv, folder), (code, stdout, stderr) in zip(cases, jobs, outputs):
        for old, new in NAMING:
            stdout, stderr = stdout.replace(old, new), stderr.replace(old, new)
        res = result(c, folder, code, stdout, stderr)
        if TRACEBACK in stderr:
            # MDATK crashed (an exception its main() does not catch: MemoryError, zlib.error,
            # AttributeError ...): Python prints the traceback and exits with 1, MDBridge reports
            # "error: <message>"; the exit status, stdout and files are compared
            res["stderr"] = None
            res["mdatk_crash"] = stderr.strip().splitlines()[-1]
        results[c["name"]] = res
    return results


def run_mdbridge(cases, root, exe):
    env = dict(os.environ, APPIMAGE_EXTRACT_AND_RUN="1")
    out = {}
    for c in cases:
        folder = os.path.join(root, "mdbridge", c["name"])
        prepare(c, folder)
        run = subprocess.run([*exe, *c["argv"]], cwd=folder, capture_output=True, env=env, timeout=300)
        decode = lambda b: b.decode("utf-8", "surrogateescape")  # noqa: E731
        out[c["name"]] = result(c, folder, run.returncode, decode(run.stdout), decode(run.stderr))
    return out


# ---------------------------------------------------------------------------------------------

def differences(c, expected, actual):
    if expected["exit"] != actual["exit"]:
        problems = [f"exit status: MDATK {expected['exit']}, MDBridge {actual['exit']}"]
        if c["exit_only"]:
            return problems
    else:
        problems = []
    if c["exit_only"]:
        return []
    if expected["stderr"] is None and not actual["stderr"].startswith("error: "):
        problems.append(f"MDATK crashed ({expected['mdatk_crash']}), MDBridge printed no error: "
                        f"{actual['stderr']!r}")
    for part in ("stdout", "stderr"):
        if expected[part] is not None and expected[part] != actual[part]:
            problems.append(f"{part} differs\n{diff(expected[part], actual[part])}")
    for rel in sorted(set(expected["written"]) | set(actual["written"])):
        e, a = expected["written"].get(rel), actual["written"].get(rel)
        if e is None or a is None:
            problems.append(f"{rel}: written by {'MDBridge' if e is None else 'MDATK'} only")
        elif e != a:
            if "text" in e and "text" in a:
                problems.append(f"{rel} differs\n{diff(e['text'], a['text'])}")
            else:
                problems.append(f"{rel} differs: {e} / {a}")
    return problems


def diff(expected, actual):
    lines = difflib.unified_diff(expected.splitlines(), actual.splitlines(), "MDATK", "MDBridge", lineterm="", n=2)
    return "\n".join(list(lines)[:40])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mdatk", help="checkout of MDATK v0.0.4.1")
    parser.add_argument("--python", default=sys.executable, help="interpreter that runs MDATK (default: this one)")
    parser.add_argument("--mdbridge", help="the mdbridge program, AppImage or mdbridge_commands_cli to compare")
    parser.add_argument("--write", help="write MDATK's normalised results to this JSON file instead")
    parser.add_argument("--only", nargs="+", help="run these cases only")
    parser.add_argument("--keep", help="run in this (new) folder and keep it")
    parser.add_argument("--list", action="store_true", help="list the cases")
    parser.add_argument("--child", nargs=3, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.child:
        child(*args.child)
        return 0
    cases = [c for c in CASES if not args.only or c["name"] in args.only]
    if args.list:
        for c in cases:
            print(f"{c['name']:40s} {' '.join(c['argv'])}")
        return 0
    if not args.mdatk or bool(args.mdbridge) == bool(args.write):
        parser.error("give --mdatk and one of --mdbridge / --write")
    root = os.path.abspath(args.keep) if args.keep else tempfile.mkdtemp(prefix="cli_parity_")
    os.makedirs(root, exist_ok=True)
    try:
        expected = run_mdatk(cases, root, args.mdatk, args.python)
        if args.write:
            data = {"cases": [{**c, **expected[c["name"]]} for c in cases]}
            os.makedirs(os.path.dirname(os.path.abspath(args.write)), exist_ok=True)
            with open(args.write, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(data, fh, indent=1, ensure_ascii=False)
                fh.write("\n")
            print(f"{len(cases)} cases written to {args.write}")
            return 0
        actual = run_mdbridge(cases, root, [args.mdbridge])
        failed = 0
        for c in cases:
            problems = differences(c, expected[c["name"]], actual[c["name"]])
            print(f"{'ok' if not problems else 'FAIL':4s} {c['name']}")
            for problem in problems:
                print("     " + problem.replace("\n", "\n     "))
            failed += bool(problems)
        print(f"{len(cases) - failed} of {len(cases)} cases give the same results in MDATK and MDBridge")
        return 1 if failed else 0
    finally:
        if not args.keep:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
