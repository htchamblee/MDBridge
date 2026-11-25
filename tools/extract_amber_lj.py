"""Cut data/sources/forcefields/amber_lj_types.tsv out of three AmberTools parameter files.

    python3 -I tools/extract_amber_lj.py DIR

DIR holds parm10.dat, frcmod.ionsjc_tip3p and frcmod.ions234lm_126_tip3p, the copies named in
data/sources/forcefields/README.md (their SHA-256 is checked here). Only numbers are kept: for each
Amber atom type the energy analysis' typed parameters use (src/analysis/energy_parameters.cpp, the
"typed" layer), its R* (Rmin/2, Angstrom) and epsilon (kcal/mol) as the file writes them, with the
file and line they were read from. Nothing else of the files (comments, other sections) is copied.
"""

from __future__ import annotations

import hashlib
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "forcefields", "amber_lj_types.tsv")

SHA256 = {
    "parm10.dat": "1ca292d5679e35450ea111a8d1b9442b5ece8518e824a8e43ed290dfff1a2ae6",
    "frcmod.ionsjc_tip3p": "eb7a1a59bf356581f89534287dd59d1bebbb2e5de7176d7c7dd69cb1d4928a55",
    "frcmod.ions234lm_126_tip3p": "f8bdb3df35ebdd9047b2e63786b4147265d5f73eb9b5197be3cb71027a2e25f3",
}

# (file, Amber atom types read from it)
TYPES = [
    # parm10.dat, section MOD4 RE: heavy-atom types, hydrogen types and the TIP3P water oxygen
    ("parm10.dat", ["C", "CT", "N", "O", "OH", "OS", "OW", "P", "S", "F", "Cl", "Br", "I",
                    "H", "HO", "HS", "H1"]),
    # Joung & Cheatham 2008, TIP3P set
    ("frcmod.ionsjc_tip3p", ["Li+", "Na+", "K+", "Rb+", "Cs+", "F-", "Cl-", "Br-", "I-"]),
    # Li et al. 2013, 12-6 CM set for TIP3P
    ("frcmod.ions234lm_126_tip3p", ["Mg2+", "Ca2+", "Sr2+", "Ba2+", "Mn2+", "Fe2+", "Co2+", "Ni2+",
                                    "Cu2+", "Zn2+", "Cd2+"]),
]


def nonbonded_rows(path, name):
    """{type: (rstar, epsilon, line number)} of the Lennard-Jones section: parm10.dat's 'MOD4 RE'
    block (until its blank line), a frcmod file's NONBON block (until a blank line or the end)."""
    with open(path, encoding="ascii") as fh:
        lines = fh.read().split("\n")
    start = None
    for k, line in enumerate(lines):
        if (name == "parm10.dat" and line.startswith("MOD4") and line.split()[1:2] == ["RE"]) or \
                (name != "parm10.dat" and line.strip() == "NONBON"):
            start = k + 1
            break
    if start is None:
        raise SystemExit(f"{name}: no Lennard-Jones section")
    out = {}
    for k in range(start, len(lines)):
        fields = lines[k].split()
        if not fields:
            break
        if fields[0] in out:
            raise SystemExit(f"{name}: type {fields[0]} twice")
        out[fields[0]] = (fields[1], fields[2], k + 1)
    return out


def main() -> int:
    src = sys.argv[1]
    rows = []
    for name, types in TYPES:
        path = os.path.join(src, name)
        with open(path, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        if digest != SHA256[name]:
            raise SystemExit(f"{name}: SHA-256 {digest}, expected {SHA256[name]} (another copy?)")
        table = nonbonded_rows(path, name)
        for t in types:
            rstar, epsilon, line = table[t]
            float(rstar), float(epsilon)  # numbers
            rows.append((name, t, rstar, epsilon, str(line)))
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("file\ttype\trstar_angstrom\tepsilon_kcal_per_mol\tline\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
