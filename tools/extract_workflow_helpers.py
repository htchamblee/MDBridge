"""Extract the texts of MDATK 0.0.4.1's six small embedded helper scripts into src/workflow/helpers/.

    python3 -I -B tools/extract_workflow_helpers.py MDATK_CHECKOUT

MDATK writes a helper script into a generated workflow with `embed` (mdatk/workflows/scripts/
__init__.py): the module's source (`inspect.getsource`) without its docstring, its `from __future__`
line and the blank lines before the code, between a header (shebang, a new docstring, extra import
lines) and a tail (the settings as Python constants and the `sys.exit(<call>)` line). MDBridge
writes the same helpers byte for byte (owner decision D4 of docs/specs/m7-plan.md, section 4.6), so
this tool stores that middle part, the body, once per helper:

    src/workflow/helpers/<name>.py    the body: the module's lines as embed keeps them, joined with
                                      "\\n", and one final "\\n"
    src/workflow/helpers/helpers.sha256   SHA-256 of each body (sha256sum format)

<name> is the tag of the file the step writes (`<project>_<step>_<name>_<stamp>.py`): helix
(helix_dist.py, amber-restraints), prep (prep_pdb.py), cp-lib (cp_retype.py), tip4p2005
(tip4p2005.py; amber-build), convert (amber2gmx.py), posre (gmx_posre.py; gromacs).

Checks (a failure writes nothing):

* for each helper and several settings, calls, docstrings and import lines, MDATK's own `embed`
  output equals the script rebuilt from the stored body with embed's layout (`assemble` below,
  the layout the C++ side reproduces);
* the published numbers in the helpers equal data/sources/nbfix: the CP Lennard-Jones values and
  masses (cp_types.tsv) in cp_retype.py, the TIP4P/2005 constants (tip4p2005.tsv) in tip4p2005.py.

The two large embedded copies (pdbtools.py, dictionaries.py) are not extracted: owner decision D3
replaces them with calls to `mdbridge pdb clean` / `dict finish`. Standard library only.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "src", "workflow", "helpers")
SOURCES = os.path.join(ROOT, "data", "sources", "nbfix")

# written name -> MDATK module (mdatk/workflows/scripts/<module>.py), settings and call of a sample embed
HELPERS = {
    "helix": ("helix_dist", {"PAIRS": [(1, 12, "GC"), (2, 11, "AU")], "BONDS": [("O6", "N4"), ("N1", "N3")],
                             "WINDOW": (2.3, 2.6, 3.1, 3.6, 0.0, 32.0), "OFFSET": 0},
              "main(sys.argv[1:], PAIRS, BONDS, WINDOW, OFFSET)"),
    "prep": ("prep_pdb", {}, "main(sys.argv[1:], strip_h=True, strip_5p=True, gap_ter=True, strip_water=True)"),
    "cp-lib": ("cp_retype", {}, "main(sys.argv[1:])"),
    "tip4p2005": ("tip4p2005", {}, "main(sys.argv[1:])"),
    "convert": ("amber2gmx", {"ALLOW_NBFIX_LOSS": False}, "main(sys.argv[1:], ALLOW_NBFIX_LOSS)"),
    "posre": ("gmx_posre", {"ITP_PATTERN": "demo_04_posre-{mol}_20261127-1200.itp", "FORCE": 1000.0},
              "main(sys.argv[1:], ITP_PATTERN, FORCE)"),
}
DOCS = ["Prepare an RNA PDB for tleap - written by MDATK 0.0.4.1.\n\nUsage: python3 THIS_SCRIPT IN OUT",
        "  one line with spaces around  \n", "Two\nlines"]
IMPORTS = [(), ("import os", "")]


def body_of(source: str) -> list[str]:
    """The lines embed keeps of a module's source (mdatk/workflows/scripts/__init__.py:15-29)."""
    lines = source.splitlines()
    k = 0
    if lines and lines[0].startswith('"""'):
        k = 1
        if not (lines[0].count('"""') >= 2 and len(lines[0]) > 3):
            while k < len(lines) and '"""' not in lines[k]:
                k += 1
            k += 1
    rest = [ln for ln in lines[k:] if not ln.startswith("from __future__")]
    while rest and not rest[0].strip():
        rest.pop(0)
    return rest


def assemble(body: str, constants: dict, call: str, doc: str, imports=()) -> str:
    """A helper script from a stored body, in embed's layout (constants as Python repr)."""
    out = ["#!/usr/bin/env python3", f'"""{doc.strip()}', '"""', ""]
    out += list(imports)
    out += body[:-1].split("\n")  # the stored text ends with one "\n"
    out += ["", "", "# --- settings written by MDATK ---"]
    out += [f"{key} = {value!r}" for key, value in constants.items()]
    out += ["", "", 'if __name__ == "__main__":', f"    sys.exit({call})", ""]
    return "\n".join(out)


def load(path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def constants_of(source: str, names) -> dict:
    tree = ast.parse(source)
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name) and node.targets[0].id in names:
            out[node.targets[0].id] = ast.literal_eval(node.value)
    return out


def tsv(name):
    with open(os.path.join(SOURCES, name), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_numbers(sources: dict) -> list[str]:
    """Differences between the helpers' published constants and data/sources/nbfix."""
    problems = []
    cp = constants_of(sources["cp-lib"], ("CP_NONBOND", "CP_MASS", "SUBST"))
    for row in tsv("cp_types.tsv"):
        t = row["type"]
        want = (float(row["rmin_half_angstrom"]), float(row["epsilon_kcal_per_mol"]))
        if cp["CP_NONBOND"].get(t) != want:
            problems.append(f"cp_retype.py CP_NONBOND[{t}] = {cp['CP_NONBOND'].get(t)}, source {want}")
        mass = cp["CP_MASS"].get(t, "").split()
        if mass != [row["mass"], row["polarizability"]]:
            problems.append(f"cp_retype.py CP_MASS[{t}] = {cp['CP_MASS'].get(t)!r}, source {row['mass']} "
                            f"{row['polarizability']}")
        if cp["SUBST"].get(row["replaces"]) != t:
            problems.append(f"cp_retype.py SUBST does not map {row['replaces']} to {t}")
    names = {"epsilon_over_k": "EPS_K", "sigma": "SIGMA", "q_h": "Q_H", "d_om": "D_OM", "r_oh": "R_OH",
             "hoh_angle": "HOH"}
    w = constants_of(sources["tip4p2005"], tuple(names.values()))
    for row in tsv("tip4p2005.tsv"):
        ours = w.get(names[row["parameter"]])
        if ours != float(row["value"]):
            problems.append(f"tip4p2005.py {names[row['parameter']]} = {ours}, source {row['value']}")
    return problems


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        raise SystemExit(__doc__.split("\n\n")[1])
    scripts = os.path.join(os.path.abspath(argv[0]), "mdatk", "workflows", "scripts")
    embed = load(os.path.join(scripts, "__init__.py"), "_mdatk_scripts").embed
    bodies, sources, problems = {}, {}, []
    for name, (module_name, constants, call) in HELPERS.items():
        path = os.path.join(scripts, module_name + ".py")
        module = load(path, f"_mdatk_scripts_{module_name}")
        with open(path, encoding="utf-8") as fh:
            sources[name] = fh.read()
        body = "\n".join(body_of(sources[name])) + "\n"
        bodies[name] = body
        cases = 0
        for doc in DOCS:
            for imports in IMPORTS:
                for consts in (constants, {}):
                    theirs = embed(module, consts, call, doc, imports)
                    if assemble(body, consts, call, doc, imports) != theirs:
                        problems.append(f"{name}: the rebuilt script differs from MDATK's embed output")
                    cases += 1
        print(f"{name:<10} {module_name}.py: {body.count(chr(10))} lines, {cases} embed cases equal")
    problems += check_numbers(sources)
    for p in problems:
        print("DIFFERENCE: " + p)
    if problems:
        print("nothing written")
        return 1
    os.makedirs(OUT, exist_ok=True)
    sums = []
    for name, body in bodies.items():
        data = body.encode("utf-8")
        with open(os.path.join(OUT, name + ".py"), "wb") as fh:
            fh.write(data)
        sums.append(f"{hashlib.sha256(data).hexdigest()}  {name}.py\n")
    with open(os.path.join(OUT, "helpers.sha256"), "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(sums)
    print(f"wrote {len(bodies)} helper bodies and helpers.sha256 to {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
