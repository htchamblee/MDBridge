"""Cut data/sources/helix/arnott_arna.tsv out of AmberTools' nucgen.dat (the numbers only).

    python3 -I tools/extract_arnott.py nucgen.dat

nucgen.dat is the copy named in data/sources/helix/README.md (ParmEd's test/files/parm/nucgen.dat,
AmberTools' dat/leap/parm/nucgen.dat); its SHA-256 is checked here. Only the numbers of the
`$ARNA` block ("ARNOTT ARNA MONOMER UNIT") are kept, as published facts: the Arnott A-RNA fibre
monomer (Arnott, Hukins & Dover 1972) in cylindrical coordinates, one row per line of the block with
the line number, the file's atom name, the PDB v3 name, r (Angstrom), phi (degrees) and z (Angstrom)
as the file writes them. Nothing else of the file (other blocks, comments, layout) is copied.

The block lists the backbone and then the bases without naming them; the rows are assigned in the
block's order (two placeholder hydrogens, the backbone, then A, G, C, U), each base's run being
checked to consist of the atoms of that base's standard residue in the CCD
(data/sources/ccd/ribonucleotide_atoms.tsv) and to hold all of its heavy atoms. PDB v3 names:
OA / OB are OP1 / OP2 and O1' is O4' (the old names); tools/gen_chem_tables.py checks OA / OB
against the chirality at P of the CCD's adenosine monophosphate.
"""

from __future__ import annotations

import csv
import hashlib
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "helix", "arnott_arna.tsv")
CCD = os.path.join(ROOT, "data", "sources", "ccd", "ribonucleotide_atoms.tsv")
SHA256 = "2442dcb7769b00cec209c6a81710e9b26fad7e4a23e2b81389a23b56bad812bc"
BLOCK = "$ARNA"
OLD_NAMES = {"OA": "OP1", "OB": "OP2", "O1'": "O4'"}
BASES = ("A", "G", "C", "U")  # the block's order
BACKBONE = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'")


def main() -> int:
    src = sys.argv[1]
    with open(src, "rb") as fh:
        data = fh.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != SHA256:
        raise SystemExit(f"{src}: SHA-256 {digest}, expected {SHA256} (another copy?)")
    lines = data.decode("ascii").split("\n")
    start = next(k for k, ln in enumerate(lines) if ln.split()[:1] == [BLOCK])
    block = []
    for k in range(start + 1, len(lines)):
        fields = lines[k].split()
        if not fields:
            break
        if len(fields) != 4:
            raise SystemExit(f"line {k + 1}: expected name r phi z")
        for v in fields[1:]:
            float(v)
        block.append((k + 1, *fields))
    with open(CCD, encoding="utf-8", newline="") as fh:
        ccd = list(csv.DictReader(fh, delimiter="\t"))
    names = {b: {r["atom_id"]: r["element"] for r in ccd if r["comp_id"] == b} for b in BASES}

    rows, k = [], 0
    while k < len(block) and block[k][1] == "H":  # placeholders on P's position before the backbone
        rows.append((block[k][0], "placeholder", block[k][1], "", *block[k][2:]))
        k += 1
    backbone = block[k:k + len(BACKBONE)]
    if sorted(OLD_NAMES.get(r[1], r[1]) for r in backbone) != sorted(BACKBONE):
        raise SystemExit("the backbone rows are not the backbone atoms")
    rows += [(r[0], "backbone", r[1], OLD_NAMES.get(r[1], r[1]), *r[2:]) for r in backbone]
    k += len(BACKBONE)
    for base in BASES:
        run = []
        while k < len(block) and block[k][1] in names[base] and block[k][1] not in {r[1] for r in run} \
                and block[k][1] not in BACKBONE:
            run.append(block[k])
            k += 1
        heavy = {n for n, el in names[base].items() if el != "H" and n not in BACKBONE and n != "OP3"}
        if not heavy <= {r[1] for r in run}:
            raise SystemExit(f"base {base}: rows {[r[1] for r in run]} lack {sorted(heavy - {r[1] for r in run})}")
        rows += [(r[0], base, r[1], r[1], *r[2:]) for r in run]
    if k != len(block):
        raise SystemExit(f"line {block[k][0]}: row not assigned to a base")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("line\tset\tnucgen_name\tpdb_name\tr_angstrom\tphi_degrees\tz_angstrom\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")
    print(os.path.relpath(OUT, ROOT), len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
