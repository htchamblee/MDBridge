#!/usr/bin/env python3
"""Regenerate every committed reference file and compare it with the committed one, byte for byte.

    python3 tests/parity/regenerate_all.py --mdatk MDATK_CHECKOUT --out FRESH_DIR \\
        [--python PY] [--pytest-python PY] [--only STEP ...] [--list]

Runs every generator of tests/data/README.md (tests/parity/*.py and tools/gen_*.py) in the
documented order into FRESH_DIR, a new folder tree that mirrors the repository: FRESH_DIR/tests/data
for the test data, FRESH_DIR/src/core and src/io for the tables of tools/gen_*.py. The generators are
copied to FRESH_DIR/tests/parity and FRESH_DIR/tools first and run from there, so the ones that read
earlier outputs (the topology fixtures, the loader's files, demo_dry.parm7 ...) read the fresh ones,
never the committed ones. The committed public-source extracts of data/sources/ (inputs, made by
tools/extract_*.py from downloads) are copied to FRESH_DIR/data/sources: tools/gen_data_tables.py
writes src/core/include/mdbridge/core/detail/data_tables_independent.inc from them, and
tests/parity/independent_tables.py loads the same tables through it to patch MDAnalysis in the
generators that depend on them. Then every
committed file is compared with its fresh copy; each difference is printed with a unified-diff
excerpt (a hex excerpt for binary files). Exit status 1 on any difference, missing or unexpected
file, or generator failure.

PY (default: this interpreter) runs the generators; it needs Python 3.13 (textwrap, json's messages
and the Unicode tables depend on the version), MDAnalysis 2.11.0, numpy and scipy (MDATK's
dependencies, at the versions of the reference environment). --pytest-python (default: PY) runs
the two generators that import MDATK's tests (make_pdb_cases.py, make_cif_cases.py) and needs pytest
too. MDATK_CHECKOUT is a checkout of MDATK v0.0.4.1; nothing is written into it. FRESH_DIR must not
exist or be empty. Standard library only.

--only runs the named steps alone (see --list), on a tree seeded with the committed files the
steps do not write; their outputs are compared. For trying out a change to one generator.

Special cases handled here (tests/data/README.md):

* make_projects.py writes time stamps: it runs with --fixed-clock, the clock of the committed run.
  The hand-made files of projects/ (no generator) are copied in, then reference.py records what
  MDATK reads from all of them (run in projects/, names relative: messages name the file).
* MDATK's logged() leaves lock files (.NAME.mdatk.toml.lock, ignored by git) next to the project
  files make_projects.py writes; they are removed before reference.py runs.
* make_project_file_golden.py writes manifest.txt (sizes that depend on the temporary folder),
  which is deleted; make_writer_vectors.py writes two files besides the six the tests use, which
  are deleted (round2_vectors.json is python/round2_vectors.json's twin, unicode_classes.json
  was folded into src/core/unicode_tables.inc).
* Expected data is MDAnalysis 2.11.0 with MDBridge's independent tables patched in
  (independent_tables.py, applied by the generators that depend on the tables); it needs only the
  copied data/sources/ and tools/gen_data_tables.py, which run first.
* gen_step_types.py writes tests/data/step-types/step_types.json besides its table; the step-type
  cases run after it. The loader cases run before the export cases (which use the loader's files),
  the topology fixtures and their dumps before everything that reads them.
* MDAnalysis caches XTC / TRR frame offsets in hidden files next to the trajectories
  (.NAME_offsets.npz / .lock); they are removed after every step and never compared.
* Version stamps: a few files record Python's full version (3.13.N) in a header; the patch number
  is masked on both sides (PATCH_STAMPED below).
* rigid/rigid_cases.json depends on the CPU's OpenBLAS / LAPACK kernels (MACHINE_DEPENDENT below):
  where its bytes differ, it is compared as JSON with numbers within a relative 1e-12.
  Everything else must match exactly.

Committed files without a generator (NOT_GENERATED: the hand-made project files, the golden
documents of tests/core/test_text.cpp and test_topology.cpp) are listed and not compared.
"""

from __future__ import annotations

import argparse
import difflib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = "tests/data"

# committed files no generator writes (inputs and hand-made goldens); copied where a step reads them
HANDMADE_PROJECTS = ["handwritten.mdatk.toml", "datetimes.mdatk.toml", "empty_params_log.toml",
                     "missing_log.mdatk.toml", "other_format.toml", "broken.mdatk.toml", "preset.json"]
NOT_GENERATED = (["tests/data/README.md",
                  "tests/data/text/tomlio_golden.toml",       # by hand from tomlio.dumps; repeated in test_text.cpp
                  "tests/data/topology/golden_model.json"]    # MDATK's chain numbering / Reporter events
                 + [f"{DATA}/projects/{name}" for name in HANDMADE_PROJECTS])

# files whose header records Python's full version: "3.13.N" compares as "3.13.x"
PATCH_STAMPED = {
    "src/core/case_tables.inc", "src/core/unicode_tables.inc",
    "tests/data/text/textwrap_cases.toml", "tests/data/python/json_loads_cases.json",
    "tests/data/python/posixpath_cases.json", "tests/data/python/string_repr_cases.json",
    "tests/data/project-file/golden/params_log.toml",   # MDATK records platform.python_version()
}
PATCH_RE = re.compile(rb"(?<![0-9.])3\.13\.[0-9]+(?![0-9.])")

# the tables of tools/gen_*.py
DATA_TABLES = "src/core/include/mdbridge/core/detail/data_tables_independent.inc"
TABLES = ["src/core/unicode_tables.inc", "src/core/case_tables.inc", "src/core/citation_data.inc",
          "src/core/step_types.inc", DATA_TABLES, "src/core/selection_tables.inc", "src/io/pdb_tables.inc"]
SOURCES = "data/sources"  # committed inputs of tools/gen_data_tables.py, copied into the fresh tree

# Outputs that legitimately depend on the CPU: compared byte for byte, and where the bytes differ, as
# JSON documents with the same structure, strings and case count whose numbers agree within
# FLOAT_TOLERANCE (relative, at least 1 in scale). Nothing else gets this leeway.
MACHINE_DEPENDENT = {
    # numpy computes the rotations, products and SVDs of mdatk/density.py with OpenBLAS / LAPACK,
    # whose kernels are chosen per CPU (SkylakeX on the reference machine, Haswell / Zen on runners
    # without AVX-512): results differ in the last bits (relative 2e-16 on such a runner); the C++
    # test compares them within a few 1e-15 (tests/data/README.md)
    "tests/data/rigid/rigid_cases.json",
}
FLOAT_TOLERANCE = 1e-12

# files a generator writes that are deliberately not kept (tests/data/README.md)
DISCARDED = ["tests/data/project-file/golden/manifest.txt",
             "tests/data/project-file/vectors/round2_vectors.json",
             "tests/data/project-file/vectors/unicode_classes.json"]

OFFSETS_RE = re.compile(r"^\..+_offsets\.(npz|lock)$")  # MDAnalysis' hidden XTC / TRR offset caches


class Step:
    def __init__(self, name, script, args, outputs, *, flags=("-I", "-B"), pytest=False, cwd=None,
                 mdatk_path=False, before=None, after=None):
        self.name = name
        self.script = script          # relative to the root of the tree
        self.args = args              # function(ctx) -> list of arguments
        self.outputs = outputs        # paths (files or folders) relative to the root it writes
        self.flags = list(flags)
        self.pytest = pytest
        self.cwd = cwd                # function(ctx) -> folder; default: an empty folder of its own
        self.mdatk_path = mdatk_path  # MDATK importable through PYTHONPATH (scripts without -I)
        self.before = before
        self.after = after


def d(ctx, rel=""):
    return os.path.join(ctx["root"], DATA, rel) if rel else os.path.join(ctx["root"], DATA)


def copy_handmade_projects(ctx):
    folder = d(ctx, "projects")
    for name in os.listdir(folder):  # the lock files of MDATK's logged() (.NAME.mdatk.toml.lock)
        if name.startswith(".") and name.endswith(".lock"):
            os.remove(os.path.join(folder, name))
    for name in HANDMADE_PROJECTS:
        shutil.copyfile(os.path.join(ROOT, DATA, "projects", name), os.path.join(d(ctx, "projects"), name))


def remove(*rel):
    def run(ctx):
        for path in rel:
            full = os.path.join(ctx["root"], path)
            if os.path.exists(full):
                os.remove(full)
    return run


def topology_files(ctx):
    folder = d(ctx, "topology")
    return [os.path.join(folder, n) for n in sorted(os.listdir(folder)) if n.endswith(".pdb")] + \
        [os.path.join(folder, "demo_dry.parm7")]


def project_files(ctx):
    return sorted(os.listdir(d(ctx, "projects")))


M = lambda ctx: ctx["mdatk"]  # noqa: E731


def case_step(name, script, out, **kw):
    """A generator called as SCRIPT MDATK_CHECKOUT OUT_DIR."""
    return Step(name, f"tests/parity/{script}", lambda ctx: [M(ctx), d(ctx, out)], [f"{DATA}/{out}"], **kw)


# In the order of tests/data/README.md, dependencies first.
STEPS = [
    Step("unicode_tables", "tools/gen_unicode_tables.py", lambda ctx: [], ["src/core/unicode_tables.inc"],
         flags=["-I"]),
    Step("case_tables", "tools/gen_case_tables.py", lambda ctx: [], ["src/core/case_tables.inc"], flags=["-I"]),
    Step("citation_data", "tools/gen_citation_data.py", lambda ctx: [M(ctx)], ["src/core/citation_data.inc"]),
    Step("step_types", "tools/gen_step_types.py", lambda ctx: [M(ctx)],
         ["src/core/step_types.inc", f"{DATA}/step-types/step_types.json"]),
    Step("data_tables", "tools/gen_data_tables.py", lambda ctx: [], [DATA_TABLES], flags=["-I", "-B"]),
    Step("selection_tables", "tools/gen_selection_tables.py", lambda ctx: [], ["src/core/selection_tables.inc"],
         flags=["-I"]),
    Step("pdb_tables", "tools/gen_pdb_tables.py", lambda ctx: [], ["src/io/pdb_tables.inc"], flags=["-I"]),
    # project files: MDATK's own files (fixed clock), the hand-made ones, and what MDATK reads from them
    Step("projects", "tests/parity/make_projects.py", lambda ctx: ["--fixed-clock", d(ctx, "projects")],
         [f"{DATA}/projects/{n}" for n in ("default.mdatk.toml", "full.mdatk.toml", "split.mdatk.toml",
                                            "split.log.toml", "MDATK_params_log.toml")],
         flags=["-B"], mdatk_path=True, after=copy_handmade_projects),
    Step("expected", "tests/parity/reference.py", lambda ctx: ["--out", d(ctx, "expected"), *project_files(ctx)],
         [f"{DATA}/expected"], flags=["-B"], mdatk_path=True, cwd=lambda ctx: d(ctx, "projects")),
    Step("textwrap", "tests/parity/textwrap_cases.py",
         lambda ctx: [M(ctx), d(ctx, "text/textwrap_cases.toml")], [f"{DATA}/text/textwrap_cases.toml"]),
    case_step("python", "python_cases.py", "python"),
    case_step("project_file_golden", "make_project_file_golden.py", "project-file/golden",
              after=remove(f"{DATA}/project-file/golden/manifest.txt")),
    case_step("writer_vectors", "make_writer_vectors.py", "project-file/vectors",
              after=remove(f"{DATA}/project-file/vectors/round2_vectors.json",
                           f"{DATA}/project-file/vectors/unicode_classes.json")),
    case_step("workflow_goldens", "make_workflow_goldens.py", "project-file/workflow"),
    case_step("model", "make_model_cases.py", "project-model"),
    Step("step_type_cases", "tests/parity/make_step_type_cases.py", lambda ctx: [M(ctx), d(ctx, "step-types")],
         [f"{DATA}/step-types/{n}" for n in ("coerce_cases.json", "digit_cases.json", "normalise_cases.json",
                                             "requirements_cases.json", "setting_comments.json",
                                             "option_text_cases.json", "registry_cases.json")]),
    Step("selection_fixtures", "tests/parity/make_selection_fixtures.py", lambda ctx: [M(ctx), d(ctx, "topology")],
         [f"{DATA}/topology/{n}" for n in ("fixture_chains.pdb", "fixture_noelements.pdb", "fixture_segids.pdb",
                                           "fixture_chains_ba.pdb", "demo_dry.parm7")]),
    Step("topology_dumps", "tests/parity/dump_topologies.py",
         lambda ctx: [d(ctx, "topology"), *topology_files(ctx)],
         [f"{DATA}/topology/{n}.topology.json" for n in ("fixture_chains.pdb", "fixture_noelements.pdb",
                                                         "fixture_segids.pdb", "fixture_chains_ba.pdb",
                                                         "demo_dry.parm7")], flags=["-I", "-B"]),
    Step("atom_selection", "tests/parity/make_atom_selection_cases.py", lambda ctx: [M(ctx), d(ctx, "selections")],
         [f"{DATA}/selections/atom_selection_cases.json", f"{DATA}/selections/atom_selection_distance_cases.json"]),
    Step("selections", "tests/parity/make_selections_cases.py", lambda ctx: [M(ctx), d(ctx, "selections")],
         [f"{DATA}/selections/selections_cases.json"]),
    case_step("parm7", "make_parm7_cases.py", "io/parm7"),
    case_step("external", "make_external_cases.py", "io/external"),
    case_step("dataset", "make_dataset_cases.py", "dataset"),
    case_step("trajectories", "make_trajectory_cases.py", "io/trajectories"),
    Step("pdb", "tests/parity/make_pdb_cases.py", lambda ctx: [M(ctx), d(ctx, "io/pdb")],
         [f"{DATA}/io/pdb/pdb_cases.json"], pytest=True),
    case_step("loader", "make_loader_cases.py", "io/loader"),
    case_step("export", "make_export_cases.py", "io/pdb/export"),          # uses the loader's files
    case_step("maps", "make_map_cases.py", "io/maps"),
    # the commands of the command line on files the steps above wrote (demo system, loader, export
    # and map files, the parameter log of projects/)
    Step("cli", "tests/parity/cli_parity.py", lambda ctx: ["--mdatk", M(ctx), "--write", d(ctx, "cli/cli_cases.json")],
         [f"{DATA}/cli/cli_cases.json"]),
    case_step("rigid", "make_rigid_cases.py", "rigid"),
    Step("distances", "tests/parity/make_distance_cases.py", lambda ctx: [M(ctx), d(ctx, "analysis")],
         [f"{DATA}/analysis/distance_cases.json"]),
    case_step("citations", "make_citation_cases.py", "citations"),
    case_step("cif", "make_cif_cases.py", "io/cif", pytest=True),
]


def remove_offsets(folder: str) -> None:
    for here, _dirs, files in os.walk(folder):
        for name in files:
            if OFFSETS_RE.match(name):
                os.remove(os.path.join(here, name))


def committed_files(prefixes: list[str]) -> list[str]:
    """Committed files (git ls-files; a plain walk outside a checkout) under the given paths."""
    try:
        run = subprocess.run(["git", "-C", ROOT, "ls-files", "-z", "--", *prefixes], capture_output=True,
                             check=True)
        files = [p for p in run.stdout.decode("utf-8", "surrogateescape").split("\0") if p]
    except (OSError, subprocess.CalledProcessError):
        files = []
        for prefix in prefixes:
            full = os.path.join(ROOT, prefix)
            if os.path.isfile(full):
                files.append(prefix)
            for here, _dirs, names in os.walk(full):
                files += [os.path.relpath(os.path.join(here, n), ROOT).replace(os.sep, "/") for n in names]
    return sorted(f for f in set(files) if not OFFSETS_RE.match(os.path.basename(f)))


def fresh_files(root: str, prefixes: list[str]) -> list[str]:
    out = []
    for prefix in prefixes:
        full = os.path.join(root, prefix)
        if os.path.isfile(full):
            out.append(prefix)
        for here, _dirs, names in os.walk(full):
            out += [os.path.relpath(os.path.join(here, n), root).replace(os.sep, "/") for n in names]
    return sorted(set(out))


def under(path: str, prefixes) -> bool:
    return any(path == p or path.startswith(p.rstrip("/") + "/") for p in prefixes)


def read(path: str) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def normalised(rel: str, data: bytes) -> bytes:
    return PATCH_RE.sub(b"3.13.x", data) if rel in PATCH_STAMPED else data


WIDTH = 240  # characters of a diff line shown


def numerically_equal(old: bytes, new: bytes) -> tuple[bool, str]:
    """JSON documents equal but for float rounding (MACHINE_DEPENDENT); the reason when not."""
    try:
        a, b = json.loads(old), json.loads(new)
    except ValueError as exc:
        return False, f"not JSON: {exc}"
    worst = [0.0, ""]

    def walk(x, y, path):
        if isinstance(x, bool) or isinstance(y, bool) or not isinstance(x, (int, float)) \
                or not isinstance(y, (int, float)):
            if type(x) is not type(y):
                return f"{path}: {x!r} / {y!r}"
            if isinstance(x, list):
                if len(x) != len(y):
                    return f"{path}: {len(x)} / {len(y)} entries"
                return next((r for r in (walk(p, q, f"{path}[{i}]") for i, (p, q) in enumerate(zip(x, y))) if r),
                            None)
            if isinstance(x, dict):
                if list(x) != list(y):
                    return f"{path}: keys {list(x)} / {list(y)}"
                return next((r for r in (walk(x[k], y[k], f"{path}.{k}") for k in x) if r), None)
            return None if x == y else f"{path}: {x!r} / {y!r}"
        if math.isnan(x) or math.isnan(y) or math.isinf(x) or math.isinf(y):
            return None if repr(x) == repr(y) else f"{path}: {x!r} / {y!r}"
        diff = abs(x - y) / max(1.0, abs(x), abs(y))
        if diff > worst[0]:
            worst[:] = [diff, f"{path}: {x!r} / {y!r}"]
        return None if diff <= FLOAT_TOLERANCE else f"{path}: {x!r} / {y!r} (relative {diff:.1e})"

    problem = walk(a, b, "$")
    return (False, problem) if problem else (True, f"largest relative difference {worst[0]:.1e} ({worst[1]})")


def lines_of(text: str) -> list[str]:
    """Lines split at "\n" only (str.splitlines also splits at \f, \x1c, U+2028 ... inside JSON strings)."""
    lines = [line + "\n" for line in text.split("\n")]
    lines[-1] = lines[-1][:-1]
    return lines if lines[-1] else lines[:-1]


def excerpt(rel: str, old: bytes, new: bytes, max_lines: int) -> list[str]:
    try:
        a, b = old.decode("utf-8"), new.decode("utf-8")
    except UnicodeDecodeError:
        a = b = None
    if a is not None and "\0" not in a and "\0" not in b:
        diff = list(difflib.unified_diff(lines_of(a), lines_of(b), f"committed/{rel}", f"fresh/{rel}", n=2))
        # long lines (one-line JSON documents): show each changed line around its first difference
        removed = [i for i, line in enumerate(diff) if line.startswith("-") and not line.startswith("---")]
        added = [i for i, line in enumerate(diff) if line.startswith("+") and not line.startswith("+++")]
        focus = {}
        for i, j in zip(removed, added):
            x, y = diff[i], diff[j]
            k = next((n for n, (p, q) in enumerate(zip(x[1:], y[1:])) if p != q), min(len(x), len(y)) - 1)
            focus[i] = focus[j] = k
        lines = []
        for i, line in enumerate(diff):
            if len(lines) >= max_lines:
                lines.append("    ... (diff truncated)\n")
                break
            body = line.rstrip("\n")
            if len(body) > WIDTH:
                k = max(0, focus.get(i, 0) - WIDTH // 3)
                body = body[0] + ("... " if k else "") + body[1 + k:1 + k + WIDTH] + " ..."
            lines.append("    " + body + "\n")
            if not line.endswith("\n"):
                lines.append("    \\ No newline at end of file\n")
        if not lines:  # identical lines: line ends or a trailing newline differ
            lines.append("    (the texts differ only in line ends)\n")
        return lines
    first = next((i for i, (x, y) in enumerate(zip(old, new)) if x != y), min(len(old), len(new)))
    start = max(0, first - 16)
    return [f"    binary: committed {len(old)} bytes, fresh {len(new)} bytes, first difference at byte {first}\n",
            f"    committed[{start}:{start + 48}] {old[start:start + 48].hex(' ')}\n",
            f"    fresh    [{start}:{start + 48}] {new[start:start + 48].hex(' ')}\n"]


def prepare(root: str, steps: list[Step], only: bool) -> None:
    """The fresh tree: the generators, the source extracts, empty output folders and (with --only)
    the committed inputs."""
    for sub in ("tools", "tests/parity"):
        os.makedirs(os.path.join(root, sub))
        for name in os.listdir(os.path.join(ROOT, sub)):
            if name.endswith(".py"):
                shutil.copyfile(os.path.join(ROOT, sub, name), os.path.join(root, sub, name))
    for rel in committed_files([SOURCES]):
        os.makedirs(os.path.dirname(os.path.join(root, rel)), exist_ok=True)
        shutil.copyfile(os.path.join(ROOT, rel), os.path.join(root, rel))
    for sub in ("src/core", "src/io", DATA):
        os.makedirs(os.path.join(root, sub), exist_ok=True)
    for step in STEPS:
        for out in step.outputs:
            folder = out if not os.path.splitext(out)[1] else os.path.dirname(out)
            os.makedirs(os.path.join(root, folder), exist_ok=True)
    if only:
        written = [o for s in steps for o in s.outputs]
        for rel in committed_files([DATA]):
            if not under(rel, written) and rel not in NOT_GENERATED[:1]:
                os.makedirs(os.path.dirname(os.path.join(root, rel)), exist_ok=True)
                shutil.copyfile(os.path.join(ROOT, rel), os.path.join(root, rel))


def run_step(step: Step, ctx: dict, logs: str) -> tuple[bool, float]:
    root = ctx["root"]
    python = ctx["pytest_python"] if step.pytest else ctx["python"]
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env["LC_ALL"] = "C.UTF-8"
    env["PYTHONIOENCODING"] = "utf-8"
    if step.mdatk_path:
        env["PYTHONPATH"] = ctx["mdatk"]
    if step.cwd:
        cwd = step.cwd(ctx)
    else:
        cwd = os.path.join(root, "run", step.name)
        os.makedirs(cwd, exist_ok=True)
    if step.before:
        step.before(ctx)
    cmd = [python, *step.flags, os.path.join(root, step.script), *step.args(ctx)]
    log_path = os.path.join(logs, step.name + ".log")
    start = time.monotonic()
    with open(log_path, "wb") as log:
        log.write(("$ " + " ".join(cmd) + f"\n(in {cwd})\n").encode())
        log.flush()
        code = subprocess.run(cmd, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    ok = code == 0
    if ok and step.after:
        step.after(ctx)
    remove_offsets(os.path.join(root, DATA))
    elapsed = time.monotonic() - start
    print(f"{'ok  ' if ok else 'FAIL'} {step.name:<20} {elapsed:7.1f} s", flush=True)
    if not ok:
        with open(log_path, "rb") as fh:
            tail = fh.read().decode("utf-8", "replace").splitlines()[-25:]
        print(f"     exit status {code}; last lines of {log_path}:")
        for line in tail:
            print("     | " + line)
    return ok, elapsed


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mdatk", help="MDATK checkout at v0.0.4.1 (read only)")
    parser.add_argument("--out", help="new folder for the fresh tree (must not exist or be empty)")
    parser.add_argument("--python", default=sys.executable,
                        help="interpreter of the generators: Python 3.13 with MDAnalysis 2.11.0, numpy, scipy")
    parser.add_argument("--pytest-python", help="interpreter with pytest too, for make_pdb_cases.py and "
                                                "make_cif_cases.py (default: --python)")
    parser.add_argument("--only", nargs="+", metavar="STEP", help="run these steps only (see --list)")
    parser.add_argument("--max-diff-lines", type=int, default=40, help="diff lines shown per file")
    parser.add_argument("--list", action="store_true", help="list the steps and their outputs")
    args = parser.parse_args(argv)

    if args.list:
        for step in STEPS:
            print(f"{step.name:<20} {step.script}  ->  {', '.join(step.outputs)}")
        print("not generated (not compared):", ", ".join(NOT_GENERATED))
        return 0
    if not args.mdatk or not args.out:
        parser.error("--mdatk and --out are required")
    mdatk = os.path.abspath(args.mdatk)
    if not os.path.isfile(os.path.join(mdatk, "mdatk", "projectfile.py")):
        parser.error(f"{mdatk} is not an MDATK checkout")
    root = os.path.abspath(args.out)
    if os.path.exists(root) and os.listdir(root):
        parser.error(f"{root} is not empty")
    steps = STEPS
    if args.only:
        names = {s.name for s in STEPS}
        unknown = [n for n in args.only if n not in names]
        if unknown:
            parser.error(f"unknown steps {unknown}; --list shows them")
        steps = [s for s in STEPS if s.name in args.only]

    os.makedirs(root, exist_ok=True)
    logs = os.path.join(root, "logs")
    os.makedirs(logs)
    prepare(root, steps, bool(args.only))
    ctx = {"root": root, "mdatk": mdatk, "python": args.python, "pytest_python": args.pytest_python or args.python}
    for key in ("python", "pytest_python"):
        version = subprocess.run([ctx[key], "-I", "-c", "import sys; print(sys.version.split()[0])"],
                                 capture_output=True, text=True)
        print(f"{key.replace('_', ' ')}: {ctx[key]} (Python {version.stdout.strip() or '?'})")
    print(f"MDATK: {mdatk}\nfresh tree: {root}\n")

    total = time.monotonic()
    failed = []
    for step in steps:
        ok, _ = run_step(step, ctx, logs)
        if not ok:
            failed.append(step.name)
    total = time.monotonic() - total

    # compare
    prefixes = sorted({o for s in steps for o in s.outputs})
    committed = [f for f in committed_files(prefixes if args.only else [DATA, *TABLES])
                 if f not in NOT_GENERATED]
    generated_by_steps = [o for s in STEPS for o in s.outputs]
    unaccounted = [f for f in committed if not under(f, generated_by_steps)]
    committed = [f for f in committed if under(f, prefixes)]
    fresh = [f for f in fresh_files(root, prefixes) if f not in DISCARDED]
    missing = sorted(set(committed) - set(fresh))
    extra = sorted(set(fresh) - set(committed))
    differing, tolerated = [], []
    for rel in sorted(set(committed) & set(fresh)):
        old, new = read(os.path.join(ROOT, rel)), read(os.path.join(root, rel))
        if old == new or normalised(rel, old) == normalised(rel, new):
            continue
        if rel in MACHINE_DEPENDENT:
            ok, why = numerically_equal(old, new)
            if ok:
                tolerated.append((rel, why))
                continue
            print(f"(machine-dependent {rel}: beyond the tolerance at {why})")
        differing.append((rel, old, new))

    print(f"\n{len(steps)} generators in {total:.0f} s; {len(committed)} committed files compared")
    if not args.only:
        print("not generated, not compared: " + ", ".join(NOT_GENERATED))
    for rel, why in tolerated:
        print(f"TOLERATED    {rel}: not byte-identical on this CPU (machine-dependent), numbers within "
              f"{FLOAT_TOLERANCE:g}; {why}")
    for rel in unaccounted:
        print(f"NO GENERATOR {rel} (committed, but written by no step; list it in NOT_GENERATED)")
    for rel in missing:
        print(f"MISSING      {rel} (committed, not written)")
    for rel in extra:
        print(f"UNEXPECTED   {rel} (written, not committed)")
    for rel, old, new in differing:
        print(f"DIFFERS      {rel}")
        sys.stdout.writelines(excerpt(rel, old, new, args.max_diff_lines))
    if failed:
        print("failed generators: " + ", ".join(failed) + f" (logs in {logs})")
    bad = failed or unaccounted or missing or extra or differing
    print("RESULT: " + ("differences found" if bad else "the committed files are what the generators write"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
