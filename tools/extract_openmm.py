"""Cut data/sources/forcefields/openmm_*.tsv out of OpenMM 8.1.2's force-field XML files.

    python3 -I tools/extract_openmm.py DIR

DIR holds these files of https://github.com/openmm/openmm/tree/8.1.2/wrappers/python/openmm/app/data
(MIT licence), with '/' in their path replaced by '_': amber14_protein.ff14SB.xml,
amber14_DNA.OL15.xml, amber14_RNA.OL3.xml, amber14_tip3p.xml, amber14_tip4pew.xml, amber14_opc.xml,
tip5p.xml, charmm36.xml, charmm36_water.xml, charmm36_tip4pew.xml, charmm36_tip5p.xml,
charmm36_spce.xml, pdbNames.xml. The amber14 files are OpenMM's conversion of AmberTools' ff14SB,
OL15, OL3 and water/ion libraries; the charmm36 files are its conversion of the CHARMM36 toppar
(c36, Aug 2015). Only residue names (protein), heavy-atom names (nucleic) and every atom name
(solvent and ions) are kept, each with the element of the atom's type ('' for virtual sites).

openmm_amber14_templates.tsv: every residue template of the four Amber files the energy analysis
matches residues against (protein.ff14SB, DNA.OL15, RNA.OL3, tip3p: TIP3P water and the ions of
frcmod.ionsjc_tip3p and frcmod.ions234lm_126_tip3p), one row per atom in the file's order, with its
atom type, Amber atom class, element, partial charge and its type's NonbondedForce sigma (nm) and
epsilon (kJ/mol), copied as the file writes them. openmm_pdbnames_residue_atom_aliases.tsv: the
residue-specific atom-name aliases of pdbNames.xml (the generic Protein / Nucleic / HOH ones are in
openmm_pdbnames_aliases.tsv). tools/gen_energy_lj.py turns them into src/analysis/energy_lj_tables.inc.
"""

from __future__ import annotations

import os
import sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "forcefields")

# (file, residue filter) — None keeps every residue of the file.
CHARMM_PROTEIN = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD", "HSE", "HSP", "ILE",
                  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "CYM"}
CHARMM_NUCLEIC = {"ADE", "CYT", "GUA", "THY", "URA"}
FILES = [
    ("amber14_protein.ff14SB.xml", "protein", None),
    ("amber14_DNA.OL15.xml", "nucleic", None),
    ("amber14_RNA.OL3.xml", "nucleic", None),
    ("amber14_tip3p.xml", "solvent", None),
    ("amber14_tip4pew.xml", "solvent", None),
    ("amber14_opc.xml", "solvent", None),
    ("tip5p.xml", "solvent", None),
    ("charmm36.xml", "protein", CHARMM_PROTEIN),
    ("charmm36.xml", "nucleic", CHARMM_NUCLEIC),
    ("charmm36_water.xml", "solvent", None),
    ("charmm36_tip4pew.xml", "solvent", None),
    ("charmm36_tip5p.xml", "solvent", None),
    ("charmm36_spce.xml", "solvent", None),
]
# the residue templates of the energy analysis (openmm_amber14_templates.tsv)
TEMPLATE_FILES = ["amber14_protein.ff14SB.xml", "amber14_DNA.OL15.xml", "amber14_RNA.OL3.xml", "amber14_tip3p.xml"]


def main() -> int:
    src = sys.argv[1]
    rows = set()
    for name, group, keep in FILES:
        root = ET.parse(os.path.join(src, name)).getroot()
        elements = {t.get("name"): t.get("element") or "" for t in root.iter("Type")}
        for res in root.iter("Residue"):
            rname = res.get("name")
            if keep is not None and rname not in keep:
                continue
            if group == "protein":
                rows.add((name, group, rname, "-", "-"))  # only the residue name is needed
                continue
            for atom in res.findall("Atom"):
                element = elements.get(atom.get("type"), "")
                if group == "nucleic" and element == "H":
                    continue  # nucleic atom sets are heavy-atom sets
                rows.add((name, group, rname, atom.get("name"), element))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "openmm_residue_atoms.tsv"), "w", encoding="utf-8") as fh:
        fh.write("file\tgroup\tresidue\tatom\telement\n")
        for r in sorted(rows):
            fh.write("\t".join(r) + "\n")
    aliases = set()
    root = root_pdbnames = ET.parse(os.path.join(src, "pdbNames.xml")).getroot()
    for res in root.iter("Residue"):
        rtype = res.get("type", "")
        for key, value in res.attrib.items():
            if key.startswith("alt"):
                aliases.add(("residue", res.get("name"), "", value, rtype))
        if res.get("name") == "HOH" or rtype == "" and res.get("name") in ("Protein", "Nucleic"):
            for atom in res.findall("Atom"):
                for key, value in atom.attrib.items():
                    if key.startswith("alt") or key == "atl2":  # 'atl2' is a typo in pdbNames.xml
                        aliases.add(("atom", res.get("name"), atom.get("name"), value, rtype))
    with open(os.path.join(OUT, "openmm_pdbnames_aliases.tsv"), "w", encoding="utf-8") as fh:
        fh.write("kind\tresidue\tatom\talias\ttype\n")
        for r in sorted(aliases):
            fh.write("\t".join(r) + "\n")
    templates = []
    for name in TEMPLATE_FILES:
        root = ET.parse(os.path.join(src, name)).getroot()
        types = {t.get("name"): t for t in root.iter("Type")}
        nonbonded = {a.get("type"): a for a in root.find("NonbondedForce").findall("Atom")}
        for res in root.iter("Residue"):
            for atom in res.findall("Atom"):
                t = types[atom.get("type")]
                lj = nonbonded[atom.get("type")]
                templates.append((name, res.get("name"), atom.get("name"), atom.get("type"), t.get("class"),
                                  t.get("element") or "", atom.get("charge"), lj.get("sigma"), lj.get("epsilon")))
    with open(os.path.join(OUT, "openmm_amber14_templates.tsv"), "w", encoding="utf-8") as fh:
        fh.write("file\tresidue\tatom\ttype\tclass\telement\tcharge\tsigma_nm\tepsilon_kj_per_mol\n")
        for r in templates:
            fh.write("\t".join(r) + "\n")
    specific = set()
    for res in root_pdbnames.iter("Residue"):
        if res.get("name") in ("All", "Protein", "Nucleic", "HOH"):
            continue
        for atom in res.findall("Atom"):
            for key, value in atom.attrib.items():
                if key.startswith("alt") or key == "atl2":
                    specific.add((res.get("name"), atom.get("name"), value))
    with open(os.path.join(OUT, "openmm_pdbnames_residue_atom_aliases.tsv"), "w", encoding="utf-8") as fh:
        fh.write("residue\tatom\talias\n")
        for r in sorted(specific):
            fh.write("\t".join(r) + "\n")
    print(len(rows), len(aliases), len(templates), len(specific))
    return 0


if __name__ == "__main__":
    sys.exit(main())
