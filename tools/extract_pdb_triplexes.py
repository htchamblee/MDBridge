"""Cut data/sources/pdb/<id>_triplex_atoms.tsv and <id>_entry.tsv out of the wwPDB mmCIF files of the
RNA triple helices the structure builder's Hoogsteen duplex and pyrimidine-motif triplex take their
geometry from (docs/DIFFERENCES.md CH-18).

    for id in 3p22 4plx 6svs 7jnh; do
        curl -O https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/${id:1:2}/$id.cif.gz
    done
    python3 -I tools/extract_pdb_triplexes.py 3p22.cif.gz 4plx.cif.gz 6svs.cif.gz 7jnh.cif.gz

PDB entries are released by the wwPDB under CC0 1.0. The base triples used are listed in
data/sources/pdb/triplex_levels.tsv (see data/sources/pdb/README.md for how they were chosen); this
tool keeps the ATOM records of the residues of those triples, as the files write them (the fields of
tools/extract_pdb_entry.py), and the entries' id, title, method, resolution, deposition date, latest
revision and primary citation. The SHA-256 of every file is checked, so the extract is the one the
README describes.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "pdb")
SHA256 = {  # *.cif.gz of 2026-11-27
    "3P22": "22aca08798b54da2619dcc588b7b2211f98c23815a2dcd8b376c1b01a8ff885c",
    "4PLX": "ab735099ead9988d39993742e91e24013e1dd283728974b8ae4934f98946d4f3",
    "6SVS": "57b79d3b391209907832d42412fd120f60674429dd10a9a9ec8664fe26d6a65d",
    "7JNH": "f0494fb55ac4f815d6e8264e0e2b57f651f35d6d9240edc77feadecac4f31abc",
}


def _entry_tool():
    spec = importlib.util.spec_from_file_location("extract_pdb_entry",
                                                  os.path.join(ROOT, "tools", "extract_pdb_entry.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def single_items(path, token, unquote):
    """key -> value of the single-valued items, also those whose value is on the next line."""
    items: dict[str, str] = {}
    pending = None
    in_loop = False
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if pending is not None:
                if line.startswith(";"):
                    value = line[1:].strip()
                    for more in fh:  # a text field: up to the closing ';'
                        if more.startswith(";"):
                            break
                        value += " " + more.strip()
                    items[pending] = value.strip()
                else:
                    items[pending] = unquote(token.findall(line)[0]) if line.strip() else ""
                pending = None
                continue
            if line.startswith("#"):
                in_loop = False
                continue
            if line.startswith("loop_"):
                in_loop = True
                continue
            if line.startswith("_") and not in_loop:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    items[parts[0]] = unquote(parts[1].strip())
                else:
                    pending = parts[0]
    return items


def levels():
    rows = []
    with open(os.path.join(OUT, "triplex_levels.tsv"), encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    return rows


def main() -> int:
    tool = _entry_tool()
    token, unquote = tool._tokens()
    wanted: dict[str, set[tuple[str, str]]] = {}
    for row in levels():
        for strand in ("purine", "wc", "hoogsteen"):
            wanted.setdefault(row["entry"], set()).add((row[strand + "_chain"], row[strand + "_seq"]))
    for src in sys.argv[1:]:
        entry = os.path.basename(src).split(".")[0].upper()
        with open(src, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        if digest != SHA256.get(entry):
            raise SystemExit(f"{src}: SHA-256 {digest}, expected {SHA256.get(entry)} (another release of the file?)")
        _items, atoms, revisions = tool.read_cif(src)
        items = single_items(src, token, unquote)
        keep = wanted[entry]
        rows = [tuple(a[k] for k in tool.ATOM_FIELDS) for a in atoms
                if a["group_PDB"] == "ATOM" and (a["auth_asym_id"], a["auth_seq_id"]) in keep]
        found = {(r[5], r[6]) for r in rows}
        if found != keep:
            raise SystemExit(f"{src}: residues {sorted(keep - found)} of triplex_levels.tsv have no ATOM records")
        name = entry.lower()
        with open(os.path.join(OUT, f"{name}_triplex_atoms.tsv"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(tool.ATOM_FIELDS) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        with open(os.path.join(OUT, f"{name}_entry.tsv"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("key\tvalue\n")
            for key in tool.ENTRY_FIELDS:
                fh.write(f"{key}\t{items.get(key, '?')}\n")
            fh.write(f"_pdbx_audit_revision_history.revision_date (latest)\t{max(revisions)}\n")
        print(f"{name}_triplex_atoms.tsv", len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
