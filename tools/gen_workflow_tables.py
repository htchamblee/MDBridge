"""Write the workflow generators' tables (src/workflow/tables.inc) and check them against MDATK 0.0.4.1.

    python3 -I -B tools/gen_workflow_tables.py MDATK_CHECKOUT [--mg-extracts DIR --out FILE]

The names and numbers come from the public-source extracts in data/sources only:

* data/sources/amber/leap_files.tsv: the AmberTools leaprc, solvent box, 12-6-4 frcmod and ParmEd
  water-model names of every force field and water model, the four-site models, the residue and ion
  names of the strip masks (tools/extract_leap_files.py, checked against AmberTools 26);
* data/sources/nbfix/nbfix_pairs.tsv: the NBfix19 and 0BPh pair terms;
* data/sources/mg/mg_sets.tsv, mg_rna_pairs.tsv: the microMg / nanoMg sets of Grotz, Cruz-Leon &
  Schwierz (tools/extract_mg_params.py; numbers recorded as facts, data/sources/mg/README.md);
  without them the Mg2+ tables are left out (`kHaveMgSets = false`) and the run says so. --mg-extracts reads them from DIR instead (made by
  tools/extract_mg_params.py --out DIR), to check them against MDATK; with --out the tables go to
  FILE, not into the repository.

MDATK's own texts (labels, comment lines of the ParmEd input, the chloride models' names, the
citation texts and resources of the Mg2+ sets, the hydrogen-bond atom table of the helix
restraints) are read from MDATK_CHECKOUT: mdatk/workflows/amber.py and rna.py through their syntax
trees, mg_params.py (standard library only) imported from its file. Nothing is written into it.

Then every table is compared with MDATK's: names, numbers (as doubles, exactly), line texts and
order. A difference that is not in DOCUMENTED (each a docs/DIFFERENCES.md entry: MDBridge keeps
the source) fails the run, and so does a DOCUMENTED difference that no longer occurs.
Deterministic; standard library only.
"""

from __future__ import annotations

import argparse
import ast
import csv
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCES = os.path.join(ROOT, "data", "sources")
INC = os.path.join(ROOT, "src", "workflow", "tables.inc")
WIDTH = 110

# Differences from MDATK 0.0.4.1 that MDBridge makes on purpose (docs/DIFFERENCES.md): key -> entry id.
DOCUMENTED: dict[str, str] = {}

MG_REPOSITORIES = {  # tools/extract_mg_params.py REPOSITORIES
    "Magnesium-FFs": "e71379d1a88b54221018af53ca54d3c18066ba05",
    "optimizedMgFFs": "1df53e9c94fc14a97864ac35c6951085bb60d1c0",
    "MgFF_OPC": "cf8148a3131c198b19f3a2896b68bad6f6f62f5b",
}


# --- reading ---------------------------------------------------------------------------------------

def tsv(*parts):
    with open(os.path.join(SOURCES, *parts), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def module_literals(path, names):
    """{name: value} of module-level assignments with literal values."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), path)
    out = {}
    for node in tree.body:
        target = node.targets[0] if isinstance(node, ast.Assign) else \
            node.target if isinstance(node, ast.AnnAssign) else None
        if isinstance(target, ast.Name) and target.id in names and node.value is not None:
            out[target.id] = ast.literal_eval(node.value)
    missing = set(names) - set(out)
    if missing:
        raise SystemExit(f"{path}: no literal {sorted(missing)}")
    return out, tree


def module_function(path, tree, name):
    """A module-level function of MDATK's file, compiled alone (it may use only typing names)."""
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            mod = ast.Module(body=[ast.parse("from typing import Dict, List, Tuple").body[0], node],
                             type_ignores=[])
            scope: dict = {}
            exec(compile(mod, path, "exec"), scope)  # noqa: S102 - MDATK's own code, read only
            return scope[name]
    raise SystemExit(f"{path}: no function {name}")


def import_file(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # dataclasses look their module up
    spec.loader.exec_module(module)
    return module


# --- C++ text ----------------------------------------------------------------------------------------

def cpp_string(text):
    """A C++ string literal (non-ASCII as octal UTF-8 bytes, which no following digit can extend)."""
    out = []
    for ch in text:
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif " " <= ch <= "~":
            out.append(ch)
        else:
            out.append("".join(f"\\{b:03o}" for b in ch.encode("utf-8")))
    return '"' + "".join(out) + '"'


def cpp_double(x):
    text = repr(float(x))
    if text in ("inf", "-inf", "nan"):
        raise SystemExit(f"no literal for {text}")
    return text


def initialiser(fields, indent=4):
    """`{a, b, ...},` wrapped at WIDTH."""
    pad = " " * indent
    line = pad + "{" + ", ".join(fields) + "},"
    if len(line) <= WIDTH:
        return [line]
    out, current = [], pad + "{"
    for k, f in enumerate(fields):
        piece = f + (", " if k + 1 < len(fields) else "},")
        if len(current) + len(piece.rstrip()) > WIDTH and current.strip() != "{":
            out.append(current.rstrip())
            current = pad + "    "
        current += piece
    out.append(current.rstrip())
    return out


# --- the tables ----------------------------------------------------------------------------------------

def leap_tables(leap_rows):
    """{(table, id): {field: value}} in the TSV's order."""
    out = {}
    for row in leap_rows:
        out.setdefault((row["table"], row["id"]), {}).setdefault(row["field"], []).append(row["value"])
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mdatk", help="MDATK checkout at v0.0.4.1 (read only)")
    parser.add_argument("--mg-extracts", help="folder with mg_sets.tsv / mg_rna_pairs.tsv (default data/sources/mg)")
    parser.add_argument("--out", default=INC, help="file to write (default src/workflow/tables.inc)")
    args = parser.parse_args(argv)
    mdatk = os.path.abspath(args.mdatk)
    wf = os.path.join(mdatk, "mdatk", "workflows")
    amber_py = os.path.join(wf, "amber.py")
    lit, amber_tree = module_literals(amber_py, ["FORCE_FIELDS", "WATERS", "FOUR_SITE_WATERS", "DEFAULT_STRIP",
                                                 "SOLVENT_MASK", "NBFIX19_LINES", "BPH0_LINES", "LRNA_NOTE"])
    lrna_pair_lines = module_function(amber_py, amber_tree, "lrna_pair_lines")
    rna, _ = module_literals(os.path.join(wf, "rna.py"), ["HBOND_ATOMS", "CANONICAL", "WOBBLE"])
    mg = import_file(os.path.join(wf, "mg_params.py"), "_mdatk_mg_params")

    differences: list[tuple[str, str]] = []  # (key, text)

    def compare(key, ours, theirs):
        if ours != theirs:
            differences.append((key, f"MDBridge {ours!r}, MDATK {theirs!r}"))

    leap = leap_tables(tsv("amber", "leap_files.tsv"))
    out = ["// Tables of the workflow generators (docs/specs/m7-plan.md 3.3). Generated by",
           "// tools/gen_workflow_tables.py from data/sources/amber/leap_files.tsv (AmberTools 26 / ParmEd 4.3.1",
           "// names), data/sources/nbfix/nbfix_pairs.tsv and data/sources/mg (Grotz, Cruz-Leon & Schwierz), with",
           "// MDATK 0.0.4.1's own texts (labels, comment lines, citations); checked equal to MDATK's tables.",
           "// Do not edit.",
           "//",
           "// Included inside a namespace that defines (std::string_view members, all literals UTF-8):",
           "//   ForceField {id, label, leaprc, short_name}",
           "//   Water {id, label, leaprc, box, parmed_1264, frcmod_1264, four_site (bool)}",
           "//   LjPair {modification, type_i, type_j, rmin, epsilon}   (numbers as the ParmEd lines write them)",
           "//   HbondAtoms {pair, atom_i, atom_j}",
           "//   MgModel {id, name}; MgCitation {water, text}; MgResource {model, water, resource}",
           "//   MgSet {model, water, ion_sigma_nm, ion_epsilon_kj, chloride, chloride_sigma_nm, chloride_epsilon_kj,",
           "//          mg_cl_sigma_nm, mg_cl_epsilon_kj, source, first_pair, pair_count}  (numbers double,",
           "//          first_pair / pair_count std::size_t into kMgRnaPairs)",
           "//   MgPair {rna_type, sigma_nm, epsilon_kj}",
           ""]

    # force fields
    out += ["// RNA force fields: {id, label, leaprc, short name used in generated files} (amber.py FORCE_FIELDS)",
            "inline constexpr ForceField kForceFields[] = {"]
    ours_ff = {}
    for fid, (label, _rc, short) in lit["FORCE_FIELDS"].items():
        rc = leap.get(("forcefield", fid), {}).get("leaprc", [None])[0]
        if rc is None:
            raise SystemExit(f"leap_files.tsv: no force field {fid}")
        ours_ff[fid] = (label, rc, short)
        out += initialiser([cpp_string(fid), cpp_string(label), cpp_string(rc), cpp_string(short)])
    out += ["};", ""]
    compare("FORCE_FIELDS", ours_ff, lit["FORCE_FIELDS"])
    extra = [i for (t, i) in leap if t == "forcefield" and i not in lit["FORCE_FIELDS"]]
    if extra:
        differences.append(("FORCE_FIELDS ids", f"leap_files.tsv has {extra}, MDATK does not"))

    # waters
    out += ["// Water models: {id, label, leaprc, solvent box unit, ParmEd add12_6_4 water model, 12-6-4 ion",
            "// frcmod, four-site (an extra point)} (amber.py WATERS, FOUR_SITE_WATERS); '' = none",
            "inline constexpr Water kWaters[] = {"]
    ours_w, four = {}, []
    for wid, (label, *_rest) in lit["WATERS"].items():
        if wid == "none":
            ours_w[wid] = (label, "", "", "", "")
            out += initialiser([cpp_string(wid), cpp_string(label)] + [cpp_string("")] * 4 + ["false"])
            continue
        f = leap.get(("water", wid))
        if f is None:
            raise SystemExit(f"leap_files.tsv: no water {wid}")
        vals = tuple(f[k][0] for k in ("leaprc", "box", "parmed_1264", "frcmod_1264"))
        ours_w[wid] = (label,) + vals
        if f["four_site"][0] == "1":
            four.append(wid)
        out += initialiser([cpp_string(wid), cpp_string(label)] + [cpp_string(v) for v in vals]
                           + ["true" if f["four_site"][0] == "1" else "false"])
    out += ["};", ""]
    compare("WATERS", ours_w, lit["WATERS"])
    compare("FOUR_SITE_WATERS", tuple(four), lit["FOUR_SITE_WATERS"])

    # strip masks
    names = leap[("strip", "water")]["alias"] + leap[("strip", "ion")]["unit"]
    strip = ":" + ",".join(names)
    compare("DEFAULT_STRIP", strip, lit["DEFAULT_STRIP"])
    compare("SOLVENT_MASK", "!" + strip, lit["SOLVENT_MASK"])
    out += ["// Residues stripped from dry outputs, and the solute mask (amber.py DEFAULT_STRIP, SOLVENT_MASK)",
            f"inline constexpr std::string_view kDefaultStrip = {cpp_string(strip)};",
            f"inline constexpr std::string_view kSolventMask = {cpp_string('!' + strip)};", ""]

    # NBfix19 / 0BPh
    pairs = tsv("nbfix", "nbfix_pairs.tsv")
    out += ["// NBfix19 and 0BPh pair terms: {modification, type i, type j, Rmin_ij (Angstrom), epsilon_ij",
            "// (kcal/mol)} (data/sources/nbfix/nbfix_pairs.tsv)",
            "inline constexpr LjPair kNbfixPairs[] = {"]
    for p in pairs:
        out += initialiser([cpp_string(p[k]) for k in ("modification", "type_i", "type_j", "rmin_angstrom",
                                                       "epsilon_kcal_per_mol")])
    out += ["};", ""]

    def pair_line(p, type_j=None):
        return f"changeLJPair @%{p['type_i']} @%{type_j or p['type_j']} {p['rmin_angstrom']} " \
               f"{p['epsilon_kcal_per_mol']}"

    def rebuild(name, lines, modification):
        """MDATK's lines with every changeLJPair line rebuilt from the TSV (in its order)."""
        rows = [p for p in pairs if p["modification"] == modification]
        built, k = [], 0
        for line in lines:
            if line.startswith("changeLJPair"):
                if k >= len(rows):
                    raise SystemExit(f"{name}: more changeLJPair lines than {modification} pairs")
                built.append(pair_line(rows[k]))
                k += 1
            else:  # comments and addLJType: MDATK's text
                built.append(line)
        if k != len(rows):
            differences.append((name, f"{len(rows) - k} {modification} pairs without a line"))
        compare(name, built, lines)
        return built

    for cname, mod, what in (("kNbfix19Lines", "nbfix19", "NBFIX19_LINES"),
                             ("kBph0Lines", "0bph", "BPH0_LINES")):
        built = rebuild(what, lit[what], mod)
        out += [f"// ParmEd input lines of {mod} for OL3 (amber.py {what}; MDATK's comment lines)",
                f"inline constexpr std::string_view {cname}[] = {{"]
        out += [f"    {cpp_string(line)}," for line in built] + ["};", ""]
    out += ["// note written when L-RNA pairs are left out (amber.py LRNA_NOTE)",
            "inline constexpr std::string_view kLrnaNote[] = {"]
    out += [f"    {cpp_string(line)}," for line in lit["LRNA_NOTE"]] + ["};", ""]
    # L-RNA: the OS pairs again for the library's bridging type (amber.py lrna_pair_lines)
    for nb in (False, True):
        for bp in (False, True):
            for bridging in ("LD", "LE"):
                ours = ["#", f"# L-RNA pairs: bridging phosphate oxygens of the L-RNA library have the type {bridging}"]
                ours += [pair_line(p, bridging) for p in pairs if p["type_j"] == "OS" and
                         ((p["modification"] == "nbfix19" and nb) or (p["modification"] == "0bph" and bp))]
                compare(f"lrna_pair_lines({nb}, {bp}, {bridging})", ours, lrna_pair_lines(nb, bp, bridging))

    # hydrogen-bond atoms of the helix restraints (MDATK's own table)
    out += ["// Hydrogen-bonded atom pairs of the helix-preserving restraints, per pair type <base i><base j>",
            "// (rna.py HBOND_ATOMS; CANONICAL / WOBBLE: the pair types restrained without / with wobble)",
            "inline constexpr HbondAtoms kHbondAtoms[] = {"]
    for ptype, atoms in rna["HBOND_ATOMS"].items():
        for a, b in atoms:
            out += initialiser([cpp_string(ptype), cpp_string(a), cpp_string(b)])
    out += ["};",
            "inline constexpr std::string_view kCanonicalPairs[] = {"
            + ", ".join(cpp_string(p) for p in rna["CANONICAL"]) + "};",
            "inline constexpr std::string_view kWobblePairs[] = {"
            + ", ".join(cpp_string(p) for p in rna["WOBBLE"]) + "};", ""]

    # Mg2+ sets
    out += mg_tables(mg, differences, args.mg_extracts or os.path.join(SOURCES, "mg"))

    # differences
    status = 0
    seen = set()
    for key, text in differences:
        entry = DOCUMENTED.get(key)
        seen.add(key)
        print(f"{'documented ' + entry if entry else 'DIFFERENCE'}: {key}: {text}")
        if not entry:
            status = 1
    for key in sorted(set(DOCUMENTED) - seen):
        print(f"STALE: {key} ({DOCUMENTED[key]}) no longer differs")
        status = 1
    if status:
        print("tables.inc not written")
        return status
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out).rstrip("\n") + "\n")
    print(f"wrote {args.out}; the tables equal MDATK 0.0.4.1's")
    return 0


def mg_tables(mg, differences, folder):
    sets_path = os.path.join(folder, "mg_sets.tsv")
    if not os.path.exists(sets_path):
        print("note: data/sources/mg has no extracts (licence question open, data/sources/mg/README.md): "
              "the Mg2+ tables are left out")
        return ["// microMg / nanoMg sets: not generated, data/sources/mg holds no extracts (licence question",
                "// open, data/sources/mg/README.md)",
                "inline constexpr bool kHaveMgSets = false;"]
    sets = tsv(folder, "mg_sets.tsv")
    rna_pairs = tsv(folder, "mg_rna_pairs.tsv")
    out = ["// microMg / nanoMg (mg_params.py): models, citation per water model, citation resource per set",
           "inline constexpr bool kHaveMgSets = true;",
           f"inline constexpr double kKjPerKcal = {cpp_double(mg.KJ_PER_KCAL)};",
           "inline constexpr MgModel kMgModels[] = {"]
    out += initialiser_list([(cpp_string(k), cpp_string(v)) for k, v in mg.MODELS.items()])
    out += ["};", "inline constexpr MgCitation kMgCitations[] = {"]
    out += initialiser_list([(cpp_string(k), cpp_string(v)) for k, v in mg.CITATIONS.items()])
    out += ["};", "inline constexpr MgResource kMgResources[] = {"]
    out += initialiser_list([(cpp_string(m), cpp_string(w), cpp_string(r)) for (m, w), r in mg.RESOURCES.items()])
    out += ["};"]
    if mg.UNAVAILABLE:
        differences.append(("UNAVAILABLE", f"MDATK lists {mg.UNAVAILABLE}; no table for it here"))

    # the sets in MDATK's order; numbers from the extracts
    by_key = {(r["model"], r["water"]): r for r in sets}
    if set(by_key) != set(mg.MG_SETS):
        differences.append(("MG_SETS keys", f"extracts {sorted(by_key)}, MDATK {sorted(mg.MG_SETS)}"))
    set_lines, pair_lines, first = [], [], 0
    for key, theirs in mg.MG_SETS.items():
        row = by_key.get(key)
        if row is None:
            continue
        repo, name = row["file"].split("/", 1)
        source = f"https://github.com/bio-phys/{repo}/blob/{MG_REPOSITORIES[repo]}/{name}"
        label = theirs.chloride  # MDATK's name of the chloride model, with the files' type in parentheses
        if not label.endswith(f"({row['chloride_type']})"):
            differences.append((f"MG_SETS {key} chloride", f"type {row['chloride_type']}, MDATK {label!r}"))
        ion = (float(row["ion_sigma_nm"]), float(row["ion_epsilon_kj_per_mol"]))
        cl = (float(row["chloride_sigma_nm"]), float(row["chloride_epsilon_kj_per_mol"]))
        mgcl = (float(row["mg_cl_sigma_nm"]), float(row["mg_cl_epsilon_kj_per_mol"]))
        mine = [(r["rna_type"], (float(r["sigma_nm"]), float(r["epsilon_kj_per_mol"]))) for r in rna_pairs
                if (r["model"], r["water"]) == key]
        mine.sort(key=lambda t: t[0])  # MDATK's order: by type name (code points)
        for what, a, b in (("ion", ion, theirs.ion), ("chloride", cl, theirs.chloride_lj),
                           ("Mg-Cl", mgcl, theirs.mg_cl), ("source", source, theirs.source)):
            if a != b:
                differences.append((f"MG_SETS {key} {what}", f"source {a!r}, MDATK {b!r}"))
        if [t for t, _ in mine] != list(theirs.rna):
            differences.append((f"MG_SETS {key} RNA types", f"source {[t for t, _ in mine]}, MDATK {list(theirs.rna)}"))
        for t, value in mine:
            if t in theirs.rna and theirs.rna[t] != value:
                differences.append((f"MG_SETS {key} RNA pair {t}", f"source {value!r}, MDATK {theirs.rna[t]!r}"))
        set_lines += initialiser([cpp_string(key[0]), cpp_string(key[1]), cpp_double(ion[0]), cpp_double(ion[1]),
                                  cpp_string(label), cpp_double(cl[0]), cpp_double(cl[1]), cpp_double(mgcl[0]),
                                  cpp_double(mgcl[1]), cpp_string(source), str(first), str(len(mine))])
        for t, (s, e) in mine:
            pair_lines += initialiser([cpp_string(t), cpp_double(s), cpp_double(e)])
        first += len(mine)
    out += ["", "// the sets' Mg2+-RNA pairs (sigma nm, epsilon kJ/mol; GROMACS [ nonbond_params ]), per set by type name",
            "inline constexpr MgPair kMgRnaPairs[] = {"] + pair_lines + ["};", "",
            "// the sets (mg_params.py MG_SETS order): ion, chloride model, Mg2+-Cl- pair (sigma nm, epsilon kJ/mol),",
            "// the authors' parameter file, the range of kMgRnaPairs",
            "inline constexpr MgSet kMgSets[] = {"] + set_lines + ["};"]
    return out


def initialiser_list(rows):
    out = []
    for fields in rows:
        out += initialiser(list(fields))
    return out


if __name__ == "__main__":
    sys.exit(main())
