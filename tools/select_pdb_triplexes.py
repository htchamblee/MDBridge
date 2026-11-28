"""Select the stacked base triples of data/sources/pdb/triplex_levels.tsv from the wwPDB mmCIF files of
the RNA triple helices the structure builder's Hoogsteen duplex and pyrimidine-motif triplex take their
geometry from (docs/DIFFERENCES.md CH-18; data/sources/pdb/README.md).

    python3 -I tools/select_pdb_triplexes.py 3p22.cif.gz 4plx.cif.gz 6svs.cif.gz 7jnh.cif.gz
    python3 -I tools/extract_pdb_triplexes.py 3p22.cif.gz 4plx.cif.gz 6svs.cif.gz 7jnh.cif.gz

The files are the ones tools/extract_pdb_triplexes.py names (their SHA-256 is checked). A triple is a
standard purine (A, G; model 1, no alternate location or the first) with its Watson-Crick partner and a
pyrimidine on its Hoogsteen edge, by donor-acceptor distances under 3.7 A:

* U·A-U: A N1 - U N3 and A N6 - U O4 (Watson-Crick), A N7 - U N3 and A N6 - U O4 (Hoogsteen);
* C+·G-C: G N1 - C N3 and G O6 - C N4 (Watson-Crick), G N7 - C N3 and G O6 - C N4 (Hoogsteen).

Two triples are consecutive (a step) when the purine, Watson-Crick and Hoogsteen residues of the first
are linked O3'-P (under 2.0 A) to the next triple's purine, Hoogsteen residues (O3' of the first, P of
the second) and, for the antiparallel Watson-Crick strand, from the second triple's to the first's. Kept
are the triples that take part in a step, whose three residues are standard A / C / G / U residues with
every atom of the nucleoside the tables use (C1', base, sugar). The table lists them per entry in the
order of the purine's author chain and number, numbered from 1.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "pdb", "triplex_levels.tsv")
CUTOFF = 3.7   # A: donor-acceptor distance of a base-base H-bond of a triple
LINK = 2.0     # A: O3'-P of a link
CRITERIA = {  # purine, Watson-Crick, Hoogsteen residue: (purine atom, partner atom) per H-bond
    ("A", "U", "U"): ((("N1", "N3"), ("N6", "O4")), (("N7", "N3"), ("N6", "O4"))),
    ("G", "C", "C"): ((("N1", "N3"), ("O6", "N4")), (("N7", "N3"), ("O6", "N4"))),
}
NUCLEOSIDE = {  # the atoms the tables use (tools/gen_chem_tables.py nucleoside_names)
    "A": ("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"),
    "G": ("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"),
    "C": ("N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"),
    "U": ("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"),
}
SUGAR = ("C1'", "C2'", "O2'", "C3'", "O3'", "C4'", "O4'", "C5'")
LABEL = {("A", "U", "U"): "U·A-U", ("G", "C", "C"): "C+·G-C"}


def _tool(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "tools", name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def residues(atoms):
    """(author chain, author number) -> {"comp", "atoms": name -> xyz} of the standard RNA residues."""
    out = {}
    for a in atoms:
        if a["group_PDB"] != "ATOM" or a.get("pdbx_PDB_model_num", "1") != "1":
            continue
        if a["label_alt_id"] not in (".", "A") or a["label_comp_id"] not in NUCLEOSIDE:
            continue
        key = (a["auth_asym_id"], int(a["auth_seq_id"]))
        res = out.setdefault(key, {"comp": a["label_comp_id"], "atoms": {}})
        res["atoms"].setdefault(a["label_atom_id"], tuple(float(a[c]) for c in ("Cartn_x", "Cartn_y", "Cartn_z")))
    return out


def complete(res):
    return all(n in res["atoms"] for n in SUGAR + NUCLEOSIDE[res["comp"]])


def bonded(res, a, b, x, y):
    ra, rb = res[a]["atoms"], res[b]["atoms"]
    return x in ra and y in rb and math.dist(ra[x], rb[y]) < CUTOFF


def triples(res):
    """[(purine, Watson-Crick, Hoogsteen residue keys)] of the triples meeting CRITERIA: per purine the
    one whose H-bond distances (the six of CRITERIA) add up to the least (a stacked neighbour of the
    partner can also come within the cutoff)."""
    out = []
    keys = sorted(res)
    for p in keys:
        if res[p]["comp"] not in "AG":
            continue
        partners = [k for k in keys if k != p and res[k]["comp"] in "CU"]
        best = None
        for w in partners:
            for h in partners:
                types = (res[p]["comp"], res[w]["comp"], res[h]["comp"])
                if w == h or types not in CRITERIA:
                    continue
                wc, hoog = CRITERIA[types]
                if all(bonded(res, p, w, x, y) for x, y in wc) and all(bonded(res, p, h, x, y) for x, y in hoog):
                    total = sum(distance(res, p, w, x, y) for x, y in wc) + sum(distance(res, p, h, x, y)
                                                                                 for x, y in hoog)
                    if best is None or total < best[0]:
                        best = (total, (p, w, h))
        if best:
            out.append(best[1])
    return out


def distance(res, a, b, x, y):
    return math.dist(res[a]["atoms"][x], res[b]["atoms"][y])


def linked(res, a, b):
    ra, rb = res[a]["atoms"], res[b]["atoms"]
    return "O3'" in ra and "P" in rb and math.dist(ra["O3'"], rb["P"]) < LINK


def selected(res):
    trios = [t for t in triples(res) if all(complete(res[k]) for k in t)]
    keep = set()
    for a in trios:
        for b in trios:
            if linked(res, a[0], b[0]) and linked(res, b[1], a[1]) and linked(res, a[2], b[2]):
                keep |= {a, b}
    return sorted(keep)


def main() -> int:
    extract = _tool("extract_pdb_triplexes")
    entry_tool = _tool("extract_pdb_entry")
    rows = []
    for src in sys.argv[1:]:
        entry = os.path.basename(src).split(".")[0].upper()
        with open(src, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        if digest != extract.SHA256.get(entry):
            raise SystemExit(f"{src}: SHA-256 {digest}, expected {extract.SHA256.get(entry)} (another release?)")
        res = residues(entry_tool.read_cif(src)[1])
        for level, (p, w, h) in enumerate(selected(res), start=1):
            label = LABEL[(res[p]["comp"], res[w]["comp"], res[h]["comp"])]
            rows.append((entry, str(level), label, p[0], str(p[1]), w[0], str(w[1]), h[0], str(h[1])))
        print(entry, sum(r[0] == entry for r in rows), "triples")
    order = {e: i for i, e in enumerate(extract.SHA256)}
    rows.sort(key=lambda r: (order[r[0]], int(r[1])))
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("entry\tlevel\ttriple\tpurine_chain\tpurine_seq\twc_chain\twc_seq\thoogsteen_chain\thoogsteen_seq\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
