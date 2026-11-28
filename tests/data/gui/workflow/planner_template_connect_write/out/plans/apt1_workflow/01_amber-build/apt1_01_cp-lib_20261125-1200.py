#!/usr/bin/env python3
"""OL3 + CP library and parameters - written by MDBridge <VERSION>.

Usage: python3 THIS_SCRIPT AMBERHOME OUTPUT_DIR
"""

import itertools
import os
import sys
from typing import Dict, List, Tuple

RETYPE = {("OP1", "O2"): "OP", ("OP2", "O2"): "OP", ("O5'", "OS"): "OR", ("O3'", "OS"): "OR"}
SUBST = {"O2": "OP", "OS": "OR"}
# Steinbrecher, Latzer & Case 2012, Table 1 (Rmin/2 in Å, epsilon in kcal/mol)
CP_NONBOND = {"OP": (1.7493, 0.2100), "OR": (1.7718, 0.1700)}
CP_MASS = {"OP": "16.00         0.434", "OR": "16.00         0.465"}


def retype_lib(text: str) -> Tuple[str, Dict[str, int]]:
    """RNA.lib text with the CP types; returns (text, counts per new type)."""
    out, counts, in_atoms = [], {"OP": 0, "OR": 0}, False
    for line in text.splitlines():
        if line.startswith("!"):
            in_atoms = line.split()[0].endswith((".atoms", ".atomspertinfo"))
            out.append(line)
            continue
        if in_atoms:
            parts = line.split('"')
            # ' "NAME" "TYPE" typex resx flags seq elmnt chg'
            if len(parts) >= 5:
                name, atype = parts[1], parts[3]
                new = RETYPE.get((name, atype))
                if new:
                    parts[3] = new
                    counts[new] += 1
                    line = '"'.join(parts)
        out.append(line)
    return "\n".join(out) + "\n", counts


def _parm_sections(text: str) -> Dict[str, List[str]]:
    """BOND / ANGLE / DIHE / IMPROPER lines of an Amber parm*.dat file."""
    lines = text.splitlines()[1:]
    blocks: List[List[str]] = [[]]
    for line in lines:
        if line.strip():
            blocks[-1].append(line)
        elif blocks[-1]:
            blocks.append([])
    # MASS block, then the hydrophilic-types line starts the BOND block
    mass, bond = blocks[0], blocks[1][1:]
    return {"MASS": mass, "BOND": bond, "ANGLE": blocks[2], "DIHE": blocks[3], "IMPROPER": blocks[4]}


def _key_width(section: str) -> int:
    return {"BOND": 5, "ANGLE": 8, "DIHE": 11, "IMPROPER": 11}[section]


def _variants(types: List[str]) -> List[List[str]]:
    idx = [i for i, t in enumerate(types) if t in SUBST]
    out = []
    for k in range(1, len(idx) + 1):
        for chosen in itertools.combinations(idx, k):
            out.append([SUBST[t] if i in chosen else t for i, t in enumerate(types)])
    return out


def cp_terms(parm_text: str) -> Dict[str, List[str]]:
    """Bonded terms of ``parm_text`` repeated for OP / OR, keeping each multi-term dihedral together."""
    sections = _parm_sections(parm_text)
    out: Dict[str, List[str]] = {}
    for section in ("BOND", "ANGLE", "DIHE", "IMPROPER"):
        width = _key_width(section)
        groups: List[Tuple[List[str], List[str]]] = []  # (types, lines)
        for line in sections[section]:
            types = [t.strip() for t in line[:width].split("-")]
            if groups and groups[-1][0] == types and section == "DIHE":
                groups[-1][1].append(line)
            else:
                groups.append((types, [line]))
        seen, lines = set(), []
        for types, glines in groups:
            for var in _variants(types):
                key = tuple(var)
                if key in seen or tuple(reversed(var)) in seen:
                    continue
                seen.add(key)
                head = "-".join(f"{t:<2}" for t in var)
                lines += [head + gl[width:] for gl in glines]
        out[section] = lines
    return out


def frcmod_text(parm_text: str) -> str:
    terms = cp_terms(parm_text)
    lines = ["CP phosphate oxygen types OP / OR (Steinbrecher, Latzer & Case, JCTC 8, 4405 (2012)) "
             "with the parm10.dat bonded terms of O2 / OS", "MASS"]
    lines += [f"{t}  {m}" for t, m in CP_MASS.items()] + [""]
    for section in ("BOND", "ANGLE", "DIHE", "IMPROPER"):
        lines += [section if section != "IMPROPER" else "IMPROPER"] + terms[section] + [""]
    lines += ["NONB"] + [f"  {t}       {r:.4f}  {e:.4f}" for t, (r, e) in CP_NONBOND.items()] + ["", ""]
    return "\n".join(lines)


def pol_table(text: str) -> str:
    """lj_1264_pol.dat with OP / OR lines that reuse the O2 / OS polarizabilities."""
    values = {ln.split()[0]: ln.split()[1] for ln in text.splitlines() if len(ln.split()) >= 2}
    extra = [f"{new:<6}{values[old]}    CP type, value of {old}" for old, new in SUBST.items() if old in values]
    return text.rstrip("\n") + "\n" + "\n".join(extra) + "\n"


def main(argv):
    """argv: [AMBERHOME, output directory]; writes RNA_CP.lib and frcmod.cp."""
    if len(argv) != 2 or argv[0] in ("-h", "--help"):
        print("usage: python3 THIS_SCRIPT AMBERHOME OUTPUT_DIR", file=sys.stderr)
        return 2
    leap = os.path.join(argv[0], "dat", "leap")
    try:
        with open(os.path.join(leap, "lib", "RNA.lib")) as fh:
            lib, counts = retype_lib(fh.read())
        with open(os.path.join(leap, "parm", "parm10.dat")) as fh:
            frcmod = frcmod_text(fh.read())
        if not counts["OP"] or not counts["OR"]:
            raise ValueError(f"RNA.lib has no O2 / OS phosphate oxygens to retype ({counts})")
        os.makedirs(argv[1], exist_ok=True)
        with open(os.path.join(argv[1], "RNA_CP.lib"), "w") as fh:
            fh.write(lib)
        with open(os.path.join(argv[1], "frcmod.cp"), "w") as fh:
            fh.write(frcmod)
        pol = os.path.join(leap, "parm", "lj_1264_pol.dat")
        if os.path.isfile(pol):
            with open(pol) as fh:
                pol_text = pol_table(fh.read())
            with open(os.path.join(argv[1], "lj_1264_pol_cp.dat"), "w") as fh:
                fh.write(pol_text)
    except (OSError, ValueError, IndexError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"{argv[1]}: RNA_CP.lib ({counts['OP']} OP, {counts['OR']} OR atoms) and frcmod.cp written",
          file=sys.stderr)
    return 0


# --- settings written by MDBridge ---


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
