#!/usr/bin/env python3
"""Prepare an RNA PDB for tleap - written by MDATK <VERSION>.

Usage: python3 THIS_SCRIPT INPUT.pdb OUTPUT.pdb
"""

import sys

PHOSPHATE = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}
WATERS = {"HOH", "WAT", "H2O", "TIP", "TIP3", "SOL", "DOD"}
RESNAMES = {"RA": "A", "RC": "C", "RG": "G", "RU": "U", "ADE": "A", "CYT": "C", "GUA": "G",
            "URA": "U", "URI": "U", "A5": "A", "C5": "C", "G5": "G", "U5": "U", "A3": "A",
            "C3": "C", "G3": "G", "U3": "U"}


def is_hydrogen(line):
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper() in ("H", "D")
    name = line[12:16].strip().lstrip("0123456789")
    return name.upper().startswith(("H", "D"))


def prepare(lines, strip_h=True, strip_5p=True, gap_ter=True, strip_water=True, rename=True):
    out = []
    prev_key = None
    prev_chain = None
    prev_num = None
    pending_ter = False
    first_residue = False
    stats = {"atoms": 0, "hydrogens": 0, "phosphates": 0, "waters": 0, "chains": 0}
    for raw in lines:
        line = raw.rstrip("\n")
        rec = line[:6]
        if rec.startswith("ENDMDL"):
            break
        if rec.startswith("TER"):
            pending_ter = True
            continue
        if rec not in ("ATOM  ", "HETATM"):
            continue
        line = line.ljust(54)
        resname = line[17:20].strip()
        if strip_water and resname.upper() in WATERS:
            stats["waters"] += 1
            continue
        chain = line[21]
        try:
            num = int(line[22:26])
        except ValueError:
            num = prev_num
        key = (chain, line[22:27])
        if key != prev_key:
            gap = (gap_ter and prev_num is not None and num is not None
                   and num != prev_num and num != prev_num + 1)
            new_chain = prev_key is None or chain != prev_chain or pending_ter or gap
            if new_chain:
                if out and not out[-1].startswith("TER"):
                    out.append("TER")
                stats["chains"] += 1
            first_residue = new_chain
            pending_ter = False
            prev_key, prev_chain, prev_num = key, chain, num
        if strip_h and is_hydrogen(line):
            stats["hydrogens"] += 1
            continue
        name = line[12:16].strip()
        if strip_5p and first_residue and name.upper() in PHOSPHATE:
            stats["phosphates"] += 1
            continue
        if rename and resname.upper() in RESNAMES:
            line = line[:17] + f"{RESNAMES[resname.upper()]:>3s}" + line[20:]
        out.append(line.rstrip())
        stats["atoms"] += 1
    if out and not out[-1].startswith("TER"):
        out.append("TER")
    out.append("END")
    return out, stats


def main(argv, strip_h=True, strip_5p=True, gap_ter=True, strip_water=True, rename=True):
    if len(argv) != 2 or argv[0] in ("-h", "--help"):
        print("usage: python3 THIS_SCRIPT INPUT.pdb OUTPUT.pdb", file=sys.stderr)
        return 2
    try:
        with open(argv[0], "r", encoding="utf-8", errors="replace") as fh:
            lines, stats = prepare(fh, strip_h, strip_5p, gap_ter, strip_water, rename)
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if not stats["atoms"]:
        print(f"ERROR: no atoms found in {argv[0]}", file=sys.stderr)
        return 1
    with open(argv[1], "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"{argv[1]}: {stats['atoms']} atoms in {stats['chains']} chain(s); removed "
          f"{stats['hydrogens']} hydrogens, {stats['phosphates']} 5'-phosphate atoms, "
          f"{stats['waters']} water atoms", file=sys.stderr)
    return 0


# --- settings written by MDATK ---


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:], strip_h=True, strip_5p=True, gap_ter=True, strip_water=True))
