"""Write the energy analysis' Lennard-Jones and charge tables as C++ data.

    python3 -I -B tools/gen_energy_lj.py

Reads data/sources/forcefields (see its README):
  openmm_amber14_templates.tsv            residue templates of OpenMM's conversion of Amber ff14SB,
                                          OL15, OL3 and TIP3P water / ions (MIT): atom names,
                                          elements, Amber classes, charges, sigma / epsilon
  openmm_pdbnames_aliases.tsv             residue-name aliases, generic atom-name aliases (MIT)
  openmm_pdbnames_residue_atom_aliases.tsv  residue-specific atom-name aliases (MIT)
  amber_lj_types.tsv                      R* / epsilon of the Amber atom types of the typed layer
                                          (numbers read from parm10.dat and two ion frcmod files)
and writes src/analysis/energy_lj_tables.inc, included by src/analysis/energy_parameters.cpp.
sigma (nm) and epsilon (kJ/mol) become R* = sigma * 10 * 2^(1/6) / 2 (Angstrom) and
epsilon / 4.184 (kcal/mol), Amber's units. tests/parity/improvements_analyses.py imports load()
from here, so MDATK's patched energy analysis uses the same numbers. Deterministic; the build needs
no Python.
"""

import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCES = os.path.join(ROOT, "data", "sources", "forcefields")
INC = os.path.join(ROOT, "src", "analysis", "energy_lj_tables.inc")

KJ_PER_KCAL = 4.184

# force-field file -> (label, family, references of its templates). Owner policy: cite who
# originated the numbers used. The energy analysis uses only the templates' charges and
# Lennard-Jones parameters; for these protein and nucleic residues both are ff94's (Cornell et al.
# 1995): every Lennard-Jones class of the three files has the R* / epsilon of a parm94.dat type
# (checked against parm94.dat, docs/DIFFERENCES.md CIT-04), and ff14SB, OL15 and OL3 changed
# torsions only, so their own papers are not cited for this use.
FILES = {
    "amber14_protein.ff14SB.xml": ("ff14SB", "protein", ["cornell1995"]),
    "amber14_DNA.OL15.xml": ("OL15", "nucleic", ["cornell1995"]),
    "amber14_RNA.OL3.xml": ("OL3", "nucleic", ["cornell1995"]),
    "amber14_tip3p.xml": ("TIP3P", "solvent", []),  # per residue below
}
TIP3P_WATER_REFERENCES = ["jorgensen1983"]
# tip3p.xml's ions by charge: frcmod.ionsjc_tip3p (+1, -1), frcmod.ions234lm_126_tip3p (+2: the
# 12-6 CM set; +3, +4: the IOD set of Li et al. 2015), as its Info block lists them
ION_REFERENCES = {1: ["joung2008"], 2: ["li2013"], 3: ["li2015"], 4: ["li2015"]}

# typed layer: Amber type -> reference. parm94 (Cornell et al. 1995) defines C, CT, N, O, OH, OS,
# P, S, F, I, H, HO, HS, H1; parm99 (Wang et al. 2000) adds Cl and Br; OW is TIP3P's water oxygen
# (Jorgensen et al. 1983); the ions are Joung & Cheatham's (2008) and Li et al.'s (2013).
# EP (extra point: no Lennard-Jones term) cites nobody.
TYPE_REFERENCES = {"Cl": "wang2000", "Br": "wang2000", "OW": "jorgensen1983", "EP": ""}
FILE_REFERENCES = {"parm10.dat": "cornell1995", "frcmod.ionsjc_tip3p": "joung2008",
                   "frcmod.ions234lm_126_tip3p": "li2013"}

# short citations for the energy analysis' log line, by reference key (the full entries are in the
# citation table, src/core/citation_data.inc)
SHORT_CITATIONS = {
    "cornell1995": "Cornell et al. 1995, doi:10.1021/ja00124a002",
    "wang2000": "Wang et al. 2000, doi:10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F",
    "jorgensen1983": "Jorgensen et al. 1983, doi:10.1063/1.445869",
    "joung2008": "Joung and Cheatham 2008, doi:10.1021/jp8001614",
    "li2013": "Li et al. 2013, doi:10.1021/ct400146w",
    "li2015": "Li et al. 2015, doi:10.1021/jp505875v",
}


def _tsv(name):
    with open(os.path.join(SOURCES, name), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rstar_from_sigma(sigma_nm):
    """R* = Rmin / 2 in Angstrom from OpenMM's sigma in nm (Rmin = 2^(1/6) sigma)."""
    return float(sigma_nm) * 10.0 * 2.0 ** (1.0 / 6.0) / 2.0


def kcal_from_kj(epsilon_kj):
    return float(epsilon_kj) / KJ_PER_KCAL


def load():
    """The tables as plain Python data (the C++ table is written from this):

    types      {amber type: (rstar, epsilon, reference key)} of the typed layer
    templates  [{name, forcefield, family, base, variant, references, generic, specific,
                 atoms: [(name, element, class, charge, rstar, epsilon)]}] in file order
    variants   {base template name: [variant template names]}
    residue_aliases {residue name that is not a template: [template names]} (pdbNames.xml)
    atom_aliases {alias group: {alias: atom name}} ('Protein', 'Nucleic', 'HOH' and pdbNames
               residue names)
    reference_keys  every reference key, first use (the bit order of the C++ table's masks)
    element_symbols chemical element symbols, upper case, sorted (CIAAW)
    virtual_sites   [(residue, atom)] of OpenMM's water virtual sites (no element), upper case
    short      SHORT_CITATIONS
    """
    types = {}
    for row in _tsv("amber_lj_types.tsv"):
        ref = TYPE_REFERENCES.get(row["type"], FILE_REFERENCES[row["file"]])
        types[row["type"]] = (float(row["rstar_angstrom"]), float(row["epsilon_kcal_per_mol"]), ref)

    rows = _tsv("openmm_amber14_templates.tsv")
    order, by_name = [], {}
    for row in rows:
        key = row["residue"]
        if key not in by_name:
            label, family, refs = FILES[row["file"]]
            by_name[key] = {"name": key, "forcefield": label, "family": family, "references": list(refs),
                            "atoms": []}
            order.append(key)
        elif by_name[key]["forcefield"] != FILES[row["file"]][0]:
            raise SystemExit(f"template {key} in two files")
        by_name[key]["atoms"].append((row["atom"], row["element"].upper(), row["class"], float(row["charge"]),
                                      rstar_from_sigma(row["sigma_nm"]), kcal_from_kj(row["epsilon_kj_per_mol"])))
    names = set(order)
    # terminal variants: Amber's N/C-prefixed amino acids (manual 13.4.2) and 5/3/N-suffixed
    # nucleotides (13.4.3) whose stem is a template of the same file
    for t in by_name.values():
        t["base"], t["variant"] = t["name"], ""
        n = t["name"]
        if t["family"] == "protein" and len(n) == 4 and n[0] in "NC" and n[1:] in names \
                and by_name[n[1:]]["forcefield"] == t["forcefield"]:
            t["base"], t["variant"] = n[1:], n[0]
        if t["family"] == "nucleic" and n[-1] in "53N" and n[:-1] in names \
                and by_name[n[:-1]]["forcefield"] == t["forcefield"]:
            t["base"], t["variant"] = n[:-1], n[-1]
        if t["family"] == "solvent":
            t["family"] = "water" if n == "HOH" else "ion"
            if n == "HOH":
                t["references"] = list(TIP3P_WATER_REFERENCES)
            else:
                if len(t["atoms"]) != 1:
                    raise SystemExit(f"ion template {n} with {len(t['atoms'])} atoms")
                t["references"] = list(ION_REFERENCES[abs(round(t["atoms"][0][3]))])
    variants = {}
    for n in order:
        t = by_name[n]
        if t["variant"]:
            variants.setdefault(t["base"], []).append(n)

    # residue-name aliases of pdbNames.xml: classes of a name and its alternatives; a residue name
    # that is no template maps to the templates of its class (HIS -> HID, HIE, HIP; WAT -> HOH)
    classes = {}
    for row in _tsv("openmm_pdbnames_aliases.tsv"):
        if row["kind"] == "residue":
            classes.setdefault(row["residue"], {row["residue"]}).add(row["alias"])
    residue_aliases = {}
    for members in classes.values():
        targets = sorted(m for m in members if m in names and not by_name[m]["variant"])
        for m in sorted(members):
            if m not in names and targets:
                residue_aliases.setdefault(m, [])
                residue_aliases[m] = sorted(set(residue_aliases[m]) | set(targets))
    # atom-name aliases: the generic groups and each pdbNames residue's own
    atom_aliases = {}
    for row in _tsv("openmm_pdbnames_aliases.tsv"):
        if row["kind"] == "atom":
            atom_aliases.setdefault(row["residue"], {})[row["alias"]] = row["atom"]
    for row in _tsv("openmm_pdbnames_residue_atom_aliases.tsv"):
        atom_aliases.setdefault(row["residue"], {})[row["alias"]] = row["atom"]
    for t in by_name.values():
        t["generic"] = {"protein": "Protein", "nucleic": "Nucleic", "water": "HOH"}.get(t["family"], "")
        # the pdbNames residue entry of the template's stem (its name or one of its alternatives)
        t["specific"] = ""
        for name, members in sorted(classes.items()):
            if t["base"] in members and name in atom_aliases and name != "HOH":
                t["specific"] = name
        if not t["specific"] and t["base"] in atom_aliases and t["base"] != "HOH":
            t["specific"] = t["base"]
    templates = [by_name[n] for n in order]
    keys = []  # every reference key, first use (templates in file order, then the types)
    for t in templates:
        for k in t["references"]:
            if k not in keys:
                keys.append(k)
    for t in types.values():
        if t[2] and t[2] not in keys:
            keys.append(t[2])
    # chemical element symbols (upper case; CIAAW) and the virtual sites of OpenMM's water models
    # (atoms whose type has no element), by residue: the extra points of the typed layer
    symbols = set()
    for name in ("abridged_standard_atomic_weights_2024.tsv", "radioactive_most_stable_isotopes.tsv"):
        with open(os.path.join(ROOT, "data", "sources", "ciaaw", name), encoding="utf-8", newline="") as fh:
            symbols |= {row["symbol"].upper() for row in csv.DictReader(fh, delimiter="\t")}
    virtual_sites = sorted({(row["residue"].upper(), row["atom"].upper()) for row in _tsv("openmm_residue_atoms.tsv")
                            if row["group"] == "solvent" and row["element"] == ""})
    return {"types": types, "templates": templates, "variants": variants, "residue_aliases": residue_aliases,
            "atom_aliases": atom_aliases, "reference_keys": keys, "short": dict(SHORT_CITATIONS),
            "element_symbols": sorted(symbols), "virtual_sites": virtual_sites}


def _s(text):
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def main() -> int:
    data = load()
    out = ["// Generated by tools/gen_energy_lj.py from data/sources/forcefields (OpenMM 8.1.2's amber14",
           "// force-field XML and pdbNames.xml, MIT; Amber parm10.dat / frcmod ion numbers). Do not edit.",
           "// Included by src/analysis/energy_parameters.cpp inside its anonymous namespace.",
           "",
           "// typed layer: {Amber type, R* (Angstrom), epsilon (kcal/mol), reference key}",
           "constexpr AmberType kAmberTypes[] = {"]
    for name, (r, e, ref) in data["types"].items():
        out.append(f"    {{{_s(name)}, {r!r}, {e!r}, {_s(ref)}}},")
    out += ["};", "",
            "// template atoms: {name, element (upper case), Amber class, charge (e), R* (Angstrom),",
            "// epsilon (kcal/mol)}, by template in file order",
            "constexpr TemplateAtom kTemplateAtoms[] = {"]
    offsets = []
    count = 0
    for t in data["templates"]:
        offsets.append(count)
        for name, element, cls, q, r, e in t["atoms"]:
            out.append(f"    {{{_s(name)}, {_s(element)}, {_s(cls)}, {q!r}, {r!r}, {e!r}}},")
            count += 1
    out += ["};", "",
            "// references of the templates, by template (indices into kReferenceKeys)"]
    keys = data["reference_keys"]
    out.append("constexpr std::string_view kReferenceKeys[] = {")
    for k in keys:
        out.append(f"    {_s(k)},")
    out += ["};", "// {reference key, short citation of the log line}", "constexpr ShortCitation kShortCitations[] = {"]
    for k in keys:
        out.append(f"    {{{_s(k)}, {_s(data['short'][k])}}},")
    out += ["};", "",
            "// templates: {name, force field, family, base, variant, first atom, atom count, generic alias",
            "// group, specific alias group, reference mask (bit k: kReferenceKeys[k])}",
            "constexpr Template kTemplates[] = {"]
    for t, first in zip(data["templates"], offsets):
        mask = 0
        for k in t["references"]:
            mask |= 1 << keys.index(k)
        out.append(f"    {{{_s(t['name'])}, {_s(t['forcefield'])}, {_s(t['family'])}, {_s(t['base'])}, "
                   f"{_s(t['variant'])}, {first}, {len(t['atoms'])}, {_s(t['generic'])}, {_s(t['specific'])}, "
                   f"{mask}u}},")
    out += ["};", "",
            "// residue names that are no template: {name, template} (pdbNames.xml residue aliases)",
            "constexpr NamePair kResidueAliases[] = {"]
    for name in sorted(data["residue_aliases"]):
        for target in data["residue_aliases"][name]:
            out.append(f"    {{{_s(name)}, {_s(target)}}},")
    out += ["};", "",
            "// atom-name aliases: {group, alias, atom name} (pdbNames.xml), sorted",
            "constexpr AtomAlias kAtomAliases[] = {"]
    for group in sorted(data["atom_aliases"]):
        for alias in sorted(data["atom_aliases"][group]):
            out.append(f"    {{{_s(group)}, {_s(alias)}, {_s(data['atom_aliases'][group][alias])}}},")
    out.append("};")
    out += ["", "// chemical element symbols, upper case, sorted (data/sources/ciaaw)",
            "constexpr std::string_view kElementSymbols[] = {"]
    for sym in data["element_symbols"]:
        out.append(f"    {_s(sym)},")
    out += ["};", "", "// water virtual sites of OpenMM's force fields: {residue, atom}, upper case, sorted",
            "constexpr NamePair kVirtualSites[] = {"]
    for res, atom in data["virtual_sites"]:
        out.append(f"    {{{_s(res)}, {_s(atom)}}},")
    out.append("};")
    os.makedirs(os.path.dirname(INC), exist_ok=True)
    with open(INC, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
