"""Write src/core/include/mdbridge/core/detail/data_tables_independent.inc from the public-source
extracts in data/sources/.

    python3 -I tools/gen_data_tables.py

Every table is derived only from the extracts committed under data/sources/ (see
data/sources/ATTRIBUTION.md for their origin and licences):

* ccd/             wwPDB Chemical Component Dictionary (CC0 1.0)
* pdb_format/      the PDB format v3.3 and IUPAC-IUB residue names (facts)
* ciaaw/           CIAAW/IUPAC abridged standard atomic weights 2024 and the most stable
                   isotopes of radioactive elements (NUBASE2020, as listed by CIAAW) (facts)
* forcefields/     residue and atom names of the Amber, CHARMM and GROMACS force fields: the Amber 2024
                   and GROMACS manuals (names, facts) and OpenMM 8.1.2's force-field XML files (MIT)

Collision rule. MDBridge guesses elements from atom names alone (MDAnalysis' procedure: the whole
name, then the name without '*', '+', '-' and digits, is looked up in the name -> element table;
then one- and two-letter prefixes are tried against the guessable elements), and its PDB writer
starts a name in column 13 when its first two letters are an ion name. So that ordinary
biomolecules come out right, every candidate entry is checked against the reference names S: every
atom name (and PDB v2.3 alternative spelling) of the CCD standard residues (amino acids,
nucleotides, water) and every atom name of the CCD peptide-linking components whose id is a
non-standard name of the protein list (ccd/protein_list_residue_atoms.tsv: HIP, HSE, HYP, LYN;
components that only share the id with a force-field residue, such as the drug "ASH", are not
used), each with its CCD element:

* a force-field name N is dropped from the name -> element table if N is the upper-cased name or
  the stripped name (no '*', '+', '-', digits; upper case) of a name in S with another element
  (for example CA, CD, CE, HG, ND: carbon alpha / delta / epsilon, hydrogen gamma, nitrogen delta);
* a two-letter element X is dropped from the guessable set if the stripped name of a name in S
  with another element starts with X (for example Ca, Cd: CA, CD1);
* an ion name I is dropped from the PDB writer's list if it is the first two letters of a name in
  S with another element (for example CA, HG, HO: CA, HG1, HO2').

check_collisions() verifies the written tables against these rules and the run fails if one is
left. The dropped entries are listed in comments of the output.

The output is committed, so a build needs no Python. It is included (inside a namespace) by
mdbridge/core/detail/data_tables.hpp, which the selection macros, the PDB / GRO topology guessing,
the parm7 element column and the PDB writer use.
"""

from __future__ import annotations

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "sources")
OUTPUT = os.path.join(ROOT, "src", "core", "include", "mdbridge", "core", "detail",
                      "data_tables_independent.inc")


def read_tsv(*path: str) -> list[dict[str, str]]:
    """Rows of a tab-separated extract; lines starting with '#' are comments, the first other line
    is the header."""
    with open(os.path.join(SRC, *path), encoding="utf-8", newline="") as fh:
        lines = [ln for ln in fh.read().splitlines() if ln and not ln.startswith("#")]
    if path[-1] in ("standard_residues.tsv", "amber24_manual_names.tsv", "gromacs_manual_names.tsv"):
        header = {"standard_residues.tsv": ["residue", "class", "source"],
                  "amber24_manual_names.tsv": ["name", "class", "section"],
                  "gromacs_manual_names.tsv": ["name", "class", "meaning"]}[path[-1]]
    else:
        header, lines = lines[0].split("\t"), lines[1:]
    return [dict(zip(header, ln.split("\t"))) for ln in lines]


# --- residue classes ------------------------------------------------------------------------------

STANDARD = read_tsv("pdb_format", "standard_residues.tsv")
CCD_STD = read_tsv("ccd", "standard_residue_atoms.tsv")
OMM = read_tsv("forcefields", "openmm_residue_atoms.tsv")
ALIASES = read_tsv("forcefields", "openmm_pdbnames_aliases.tsv")
AMBER = read_tsv("forcefields", "amber24_manual_names.tsv")
GROMACS = read_tsv("forcefields", "gromacs_manual_names.tsv")


def residues_of(rows, cls, key="class", name="residue") -> set[str]:
    return {r[name] for r in rows if r[key] == cls}


def water_residues_in_openmm() -> set[str]:
    """Residues of the OpenMM solvent files made of one O, two H and only virtual sites."""
    atoms: dict[tuple[str, str], list[str]] = {}
    for r in OMM:
        if r["group"] == "solvent":
            atoms.setdefault((r["file"], r["residue"]), []).append(r["element"])
    return {res for (_, res), els in atoms.items()
            if sorted(e for e in els if e) == ["H", "H", "O"]}


def protein_residues() -> set[str]:
    out = residues_of(STANDARD, "protein")
    out |= {r["residue"] for r in OMM if r["group"] == "protein"}       # Amber ff14SB + CHARMM36
    out |= {r["name"] for r in AMBER if r["class"] == "protein" and "*" not in r["name"]}
    out |= residues_of(GROMACS, "protein", name="name")
    standard = residues_of(STANDARD, "protein")
    # pdbNames.xml aliases of the standard amino acids (protonation states); aliases of caps or
    # of non-standard groups are not taken.
    out |= {a["alias"] for a in ALIASES
            if a["kind"] == "residue" and a["type"] == "Protein" and a["residue"] in standard}
    return out


def nucleic_residues() -> set[str]:
    out = residues_of(STANDARD, "nucleic")
    # Amber OL15/OL3 + CHARMM36 residues that carry a base (N1 or N9); not caps such as OHE.
    out |= {r["residue"] for r in OMM if r["group"] == "nucleic" and r["atom"] in ("N1", "N9")}
    out |= {a["alias"] for a in ALIASES if a["kind"] == "residue" and a["type"] == "Nucleic"}
    return out


def water_residues() -> set[str]:
    out = residues_of(STANDARD, "water")
    out |= {r["name"] for r in AMBER if r["class"] == "water"}
    out |= residues_of(GROMACS, "water", name="name")
    out |= water_residues_in_openmm()
    out |= {a["alias"] for a in ALIASES if a["kind"] == "residue" and a["residue"] == "HOH"}
    return out


# --- atom-name sets -------------------------------------------------------------------------------

def protein_backbone() -> set[str]:
    """Heavy atoms the CCD flags as backbone in every one of the 20 IUPAC amino acids, without the
    C-terminal leaving atom (OXT)."""
    iupac20 = {r["residue"] for r in STANDARD if r["class"] == "protein" and r["residue"] not in ("SEC", "PYL")}
    sets = []
    for res in sorted(iupac20):
        sets.append({r["atom_id"] for r in CCD_STD if r["comp_id"] == res and r["element"] != "H"
                     and r["backbone_flag"] == "Y" and r["leaving_flag"] == "N"})
    return set.intersection(*sets)


def old_spelling(atom_id: str, alt: str) -> bool:
    """True if alt is the PDB v2.3 spelling of the v3 atom_id: "*" for "'", or OnP for OPn."""
    if alt == atom_id.replace("'", "*"):
        return True
    m = re.fullmatch(r"OP(\d)", atom_id)
    return bool(m) and alt == f"O{m.group(1)}P"


def nucleic_heavy_atom_names() -> set[str]:
    nucleic = {r["residue"] for r in STANDARD if r["class"] == "nucleic"}
    names = set()
    for r in CCD_STD:
        if r["comp_id"] in nucleic and r["element"] != "H":
            names.add(r["atom_id"])
            if old_spelling(r["atom_id"], r["alt_atom_id"]):
                names.add(r["alt_atom_id"])
    names |= {r["atom"] for r in OMM if r["group"] == "nucleic"}
    # pdbNames.xml spellings of heavy nucleic atoms (hydrogen names start with H).
    names |= {a["alias"] for a in ALIASES if a["kind"] == "atom" and a["residue"] == "Nucleic"
              and not a["atom"].startswith("H")}
    return names


def primed(name: str) -> bool:
    return name.endswith("'") or name.endswith("*")


def unprime(name: str) -> str:
    return name[:-1] if primed(name) else name


def phosphate(name: str) -> bool:
    return name == "P" or re.fullmatch(r"OP\d|O\dP", name) is not None


# Sugar-phosphate backbone: the atoms along the IUPAC-IUB backbone torsions alpha..zeta
# (P-O5'-C5'-C4'-C3'-O3') and the phosphate oxygens.
BACKBONE_PATH = {"P", "O5", "C5", "C4", "C3", "O3"}
# Furanose: ring atoms C1'..C4', O4' and the 2'-oxygen of ribose.
SUGAR = {"C1", "C2", "C3", "C4", "O4", "O2"}


def nucleic_sets() -> tuple[set[str], set[str], set[str]]:
    backbone, sugar, base = set(), set(), set()
    for n in nucleic_heavy_atom_names():
        if phosphate(n) or (primed(n) and unprime(n) in BACKBONE_PATH):
            backbone.add(n)
        if primed(n) and unprime(n) in SUGAR:
            sugar.add(n)
        if not primed(n) and not phosphate(n):
            base.add(n)
    return backbone, sugar, base


# --- elements and masses --------------------------------------------------------------------------

def elements() -> list[tuple[int, str, str, str, bool]]:
    """(Z, symbol, name, mass literal, is a standard atomic weight)."""
    weights = read_tsv("ciaaw", "abridged_standard_atomic_weights_2024.tsv")
    stable: dict[str, str] = {}
    for r in read_tsv("ciaaw", "radioactive_most_stable_isotopes.tsv"):
        stable.setdefault(r["symbol"], r["mass_number"])  # the first one CIAAW lists
    out = []
    for r in weights:
        if r["abridged_standard_atomic_weight"]:
            out.append((int(r["Z"]), r["symbol"], r["element"], r["abridged_standard_atomic_weight"], True))
        else:
            out.append((int(r["Z"]), r["symbol"], r["element"], stable[r["symbol"]] + ".0", False))
    assert [e[0] for e in out] == list(range(1, 119))
    return out


SYMBOLS = {e[1].upper(): e[1] for e in elements()}


def starts_with_element(name: str, element: str) -> bool:
    return bool(element) and name.upper().startswith(element.upper())


def forcefield_name_elements() -> list[tuple[str, str, str, str]]:
    """(residue, atom name, element or "" for a virtual site, origin): every atom of the one-atom
    (ion) residues of the Amber and CHARMM36 ion libraries, water atoms whose name does not start
    with their element (virtual sites), and the water atom spellings of pdbNames.xml."""
    out = set()
    sizes: dict[tuple[str, str], int] = {}
    for r in OMM:
        sizes[(r["file"], r["residue"])] = sizes.get((r["file"], r["residue"]), 0) + 1
    for r in OMM:
        if r["group"] != "solvent":
            continue
        origin = "CHARMM36" if r["file"].startswith("charmm36") else "Amber"
        if sizes[(r["file"], r["residue"])] == 1:
            out.add((r["residue"], r["atom"], r["element"], origin + " ion"))
        elif not starts_with_element(r["atom"], r["element"]):
            out.add((r["residue"], r["atom"], r["element"], origin + " water virtual site"))
    for a in ALIASES:
        if a["kind"] == "atom" and a["residue"] == "HOH":
            out.add(("*water*", a["alias"], a["atom"][0], "OpenMM pdbNames.xml"))
    return sorted(out)


def guessable_elements() -> tuple[list[str], list[str]]:
    """Elements a one- or two-letter atom-name prefix may be guessed as: the elements of the
    standard residues (CCD) and of water, and those of the CHARMM36 ion set."""
    els = {r["element"].upper() for r in CCD_STD}
    els |= {r["element"].upper() for r in OMM if r["file"] == "charmm36_water.xml" and r["element"]}
    els = sorted(SYMBOLS[e] for e in els if e in SYMBOLS)
    return [e for e in els if len(e) == 1], [e for e in els if len(e) == 2]


# --- reference names and the collision rule -------------------------------------------------------

def reference_names() -> dict[str, set[str]]:
    """S: atom name -> CCD elements (upper case) of the standard residues and of the peptide-linking
    CCD components named in the protein list. For the standard residues the alternative spelling
    (alt_atom_id, the PDB v2.3 name of the remediation: 1HB, C4*) is a reference name too; for
    other components alt_atom_id is a legacy name with no common convention (homoserine's N is
    "NA") and is not used."""
    out: dict[str, set[str]] = {}
    for r in CCD_STD:
        for name in {r["atom_id"], r["alt_atom_id"]} - {"", "?", "."}:
            out.setdefault(name, set()).add(r["element"].upper())
    for r in read_tsv("ccd", "protein_list_residue_atoms.tsv"):
        if "PEPTIDE" in r["type"].upper() and r["comp_id"] in protein_residues():
            out.setdefault(r["atom_id"], set()).add(r["element"].upper())
    return out


def stripped(name: str) -> str:
    """The name without '*', '+', '-' and digits, upper case (the guesser's second lookup key)."""
    return re.sub(r"[0-9]", "", re.sub(r"[*+-]", "", name)).upper()


def lookup_keys(name: str) -> set[str]:
    return {name.upper(), stripped(name)}


def collides(element: str, elements_of_name: set[str]) -> bool:
    return any(e != element.upper() for e in elements_of_name)


def name_table_collisions(name: str, element: str, ref: dict[str, set[str]]) -> list[str]:
    return sorted(s for s, els in ref.items() if name.upper() in lookup_keys(s) and collides(element, els))


def two_letter_collisions(symbol: str, ref: dict[str, set[str]]) -> list[str]:
    return sorted(s for s, els in ref.items() if stripped(s).startswith(symbol.upper()) and collides(symbol, els))


def ion_collisions(ion: str, ref: dict[str, set[str]]) -> list[str]:
    return sorted(s for s, els in ref.items() if s[:2].upper() == ion.upper() and collides(ion, els))


def name_elements() -> tuple[list[tuple[str, str, str]], dict[str, list[str]]]:
    """The name -> element table, keyed by upper-case name: [(name, element, origins)], and the
    dropped names with the reference names they collide with."""
    merged: dict[str, tuple[str, set[str]]] = {}
    for _res, name, element, origin in forcefield_name_elements():
        key = name.upper()
        if key in merged and merged[key][0] != element:
            raise SystemExit(f"atom name {key} maps to {merged[key][0]!r} and {element!r}")
        merged.setdefault(key, (element, set()))[1].add(origin)
    ref = reference_names()
    kept, dropped = [], {}
    for key, (element, origins) in sorted(merged.items()):
        hits = name_table_collisions(key, element, ref)
        if hits:
            dropped[key] = hits
        else:
            kept.append((key, element, ", ".join(sorted(origins))))
    return kept, dropped


def guessable_two_letter() -> tuple[list[str], dict[str, list[str]]]:
    ref = reference_names()
    _one, two = guessable_elements()
    dropped = {x: two_letter_collisions(x, ref) for x in two if two_letter_collisions(x, ref)}
    return [x for x in two if x not in dropped], dropped


def check_collisions(names, two, ions) -> None:
    """Fails if a written table still collides with a reference name (the rules above)."""
    ref = reference_names()
    problems = [f"name table {n} -> {e!r}: {name_table_collisions(n, e, ref)}"
                for n, e, _ in names if name_table_collisions(n, e, ref)]
    problems += [f"two-letter element {x}: {two_letter_collisions(x, ref)}" for x in two
                 if two_letter_collisions(x, ref)]
    problems += [f"ion name {i}: {ion_collisions(i, ref)}" for i in ions if ion_collisions(i, ref)]
    if problems:
        raise SystemExit("collisions with CCD reference names:\n  " + "\n  ".join(problems))


# --- PDB writer -----------------------------------------------------------------------------------

def monatomic_ion_names() -> list[str]:
    """Atom names of one-atom CCD components whose two-letter element symbol is the atom name."""
    rows = read_tsv("ccd", "monatomic_components.tsv")
    return sorted({r["atom_id"] for r in rows
                   if len(r["element"]) == 2 and r["atom_id"].upper() == r["element"].upper()})


def writer_ion_names() -> tuple[list[str], dict[str, list[str]]]:
    ref = reference_names()
    ions = monatomic_ion_names()
    dropped = {i: ion_collisions(i, ref) for i in ions if ion_collisions(i, ref)}
    return [i for i in ions if i not in dropped], dropped


def dropped_comment(what: str, dropped: dict[str, list[str]]) -> list[str]:
    lines = [f"// Dropped by the collision rule ({what}):"]
    for key, hits in sorted(dropped.items()):
        shown = " ".join(hits[:8]) + (" ..." if len(hits) > 8 else "")
        lines.append(f"//   {key}: {shown}")
    return lines


def hg_cl_residues() -> tuple[list[str], list[str]]:
    rows = read_tsv("ccd", "hg_cl_atom_names.tsv")
    hg = sorted({r["comp_id"] for r in rows if r["atom_id"] == "HG"})
    cl = sorted({r["comp_id"] for r in rows if r["atom_id"] == "CL"})
    return hg, cl


# --- output ---------------------------------------------------------------------------------------

def cstr(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def string_array(name: str, values, comment: str, per_line: int = 8) -> list[str]:
    values = sorted(values)
    lines = [f"// {comment}", f"inline constexpr std::string_view {name}[] = {{"]
    for i in range(0, len(values), per_line):
        lines.append("    " + " ".join(cstr(v) + "," for v in values[i:i + per_line]))
    lines.append("};")
    lines.append("")
    return lines


def main() -> int:
    backbone, sugar, base = nucleic_sets()
    one, _ = guessable_elements()
    two, two_dropped = guessable_two_letter()
    names, names_dropped = name_elements()
    ions, ions_dropped = writer_ion_names()
    check_collisions(names, two, ions)
    hg, cl = hg_cl_residues()
    out = [
        "// Residue, atom-name, element and mass tables rebuilt from public sources only:",
        "// wwPDB Chemical Component Dictionary (CC0), PDB format v3.3, IUPAC-IUB amino-acid symbols,",
        "// CIAAW abridged standard atomic weights 2024, Amber 2024 and GROMACS manuals, OpenMM 8.1.2",
        "// force-field XML (MIT). See data/sources/ATTRIBUTION.md.",
        "// Generated by tools/gen_data_tables.py; do not edit. Needs <string_view>; included inside",
        "// a namespace (see mdbridge/core/detail/data_tables.hpp).",
        "",
        "// ---- Selection macros ----",
        "",
    ]
    out += string_array("kIndProteinResidues", protein_residues(),
                        "Residue names of protein: IUPAC/PDB amino acids, Amber ff14SB (incl. N*/C* termini "
                        "and caps), CHARMM36, GROMACS variants.")
    out += string_array("kIndNucleicResidues", nucleic_residues(),
                        "Residue names of nucleic acids: PDB, Amber OL15/OL3 (5/3/N ends), CHARMM36, "
                        "OpenMM pdbNames aliases.")
    out += string_array("kIndWaterResidues", water_residues(),
                        "Residue names of water: PDB HOH, Amber solvents.lib, GROMACS SOL, CHARMM36, "
                        "OpenMM pdbNames aliases.")
    out += string_array("kIndProteinBackboneAtoms", protein_backbone(),
                        "Protein backbone: CCD backbone heavy atoms common to the 20 amino acids, OXT excluded.")
    out += string_array("kIndNucleicBackboneAtoms", backbone,
                        "Nucleic sugar-phosphate backbone: P, phosphate oxygens, O5' C5' C4' C3' O3' "
                        "(PDB v3, v2.3 '*'/OnP and force-field spellings).")
    out += string_array("kIndNucleicSugarAtoms", sugar,
                        "Nucleic sugar: C1' C2' C3' C4' O4' and O2' in every spelling.")
    out += string_array("kIndNucleicBaseAtoms", base,
                        "Nucleic base heavy atoms of the standard nucleotides (CCD, Amber, CHARMM names).")
    out += [
        "// ---- Elements and masses ----",
        "",
        "struct IndElement {",
        "    int z;",
        "    std::string_view symbol;",
        "    double mass;               // CIAAW 2024 abridged standard atomic weight, or the mass number",
        "    bool standardWeight;       // false: mass number of the most stable isotope (radioactive)",
        "};",
        "",
        "inline constexpr IndElement kIndElements[] = {",
    ]
    for z, sym, name, mass, std in elements():
        out.append(f"    {{{z}, {cstr(sym)}, {mass}, {'true' if std else 'false'}}},  // {name}")
    out += ["};", "",
            "struct IndNameElement {",
            "    std::string_view name;     // upper case",
            "    std::string_view element;  // \"\" = virtual site / lone pair (no element, mass 0)",
            "};",
            "",
            "// Atom name -> element, by name alone (upper case, sorted): force-field ion atom names, water",
            "// virtual sites, and water atom spellings (OW, HW1, HW2, OH2).",
            *dropped_comment("CCD reference names of another element", names_dropped),
            "inline constexpr IndNameElement kIndForceFieldNameElements[] = {"]
    for name, el, origin in names:
        out.append(f"    {{{cstr(name)}, {cstr(el)}}},  // {origin}")
    out += ["};", ""]
    out += string_array("kIndGuessOneLetterElements", one,
                        "One-letter elements an atom-name prefix may be guessed as (standard residues, water, ions).")
    out += dropped_comment("CCD reference names of another element start with them", two_dropped)
    out += string_array("kIndGuessTwoLetterElements", two,
                        "Two-letter elements an atom-name prefix may be guessed as (standard residues, CHARMM36 ions).")
    out += ["// ---- PDB writer atom-name placement ----", ""]
    out += dropped_comment("first two letters of CCD reference names of another element", ions_dropped)
    out += string_array("kIndMonatomicIonNames", ions,
                        "Atom names of one-atom CCD components with a two-letter element (written from column 13).")
    out += string_array("kIndMercuryHgResidues", hg,
                        "CCD components in which the atom named HG is mercury.")
    out += string_array("kIndChlorineClResidues", cl,
                        "CCD components in which the atom named CL is chlorine.", per_line=10)
    text = "\n".join(out).rstrip("\n") + "\n"
    with open(OUTPUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print(f"wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
