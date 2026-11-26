"""Cut data/sources/pdb/4jrd_atoms.tsv and 4jrd_entry.tsv out of the wwPDB mmCIF file of PDB entry 4JRD.

    curl -O https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/jr/4jrd.cif.gz
    python3 -I tools/extract_pdb_entry.py 4jrd.cif.gz

PDB entries are released by the wwPDB under CC0 1.0. 4JRD is the parallel duplex of poly(rA)
(Safaee et al., Angew. Chem. Int. Ed. 52 (2013) 10370-10373), the source of the structure builder's
4JRD nucleotide and screw (tools/gen_chem_tables.py; MDATK builder.py `_JRD_NUCLEOTIDE`,
`JRD_TWIST`, `JRD_RISE`). The SHA-256 of the file is checked (data/sources/pdb/README.md), so the
extract is the one the README describes. Kept:

* 4jrd_atoms.tsv: every ATOM record of `_atom_site` (the two RNA strands; the waters and ammonium
  ions, HETATM, are not used) as the file writes them: id, element, label atom / residue names,
  alternate location, author chain and residue number, x, y, z, occupancy, B, model number.
* 4jrd_entry.tsv: the entry's id, title, method, resolution, deposition date, latest revision and
  the primary citation (journal, volume, pages, year, DOI, PubMed id), as key / value rows.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "pdb")
SHA256 = "fede84f06b2771793f447a598756aaab40c0f35a2611fdf8ede2bb016acbe408"  # 4jrd.cif.gz of 2026-11-26

ATOM_FIELDS = ("id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "auth_asym_id",
               "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
               "pdbx_PDB_model_num")
ENTRY_FIELDS = ("_entry.id", "_struct.title", "_exptl.method", "_refine.ls_d_res_high",
                "_pdbx_database_status.recvd_initial_deposition_date", "_citation.title",
                "_citation.journal_abbrev", "_citation.journal_volume", "_citation.page_first",
                "_citation.page_last", "_citation.year", "_citation.pdbx_database_id_DOI",
                "_citation.pdbx_database_id_PubMed")


def _tokens():
    spec = importlib.util.spec_from_file_location("extract_ccd", os.path.join(ROOT, "tools", "extract_ccd.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.TOKEN, module.unquote


def read_cif(path):
    """(key -> value of the single-valued items, atom_site rows, revision dates) of a small mmCIF file
    whose values fit on their lines (true of the items read here)."""
    token, unquote = _tokens()
    items: dict[str, str] = {}
    atoms: list[dict[str, str]] = []
    revisions: list[str] = []
    headers: list[str] = []
    in_loop = reading = False
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                in_loop = reading = False
                continue
            if line.startswith("loop_"):
                in_loop, reading, headers = True, False, []
                continue
            if in_loop and not reading and line.startswith("_"):
                headers.append(line.split()[0])
                continue
            if in_loop:
                reading = True
                toks = [unquote(t) for t in token.findall(line)]
                if len(toks) != len(headers):
                    continue
                row = dict(zip(headers, toks))
                if headers[0].startswith("_atom_site."):
                    atoms.append({h.split(".", 1)[1]: v for h, v in row.items()})
                elif headers[0].startswith("_pdbx_audit_revision_history."):
                    revisions.append(row["_pdbx_audit_revision_history.revision_date"])
                continue
            if line.startswith("_"):
                parts = line.split(None, 1)
                if len(parts) == 2:
                    items[parts[0]] = unquote(parts[1].strip())
    return items, atoms, revisions


def main() -> int:
    src = sys.argv[1]
    with open(src, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    if digest != SHA256:
        raise SystemExit(f"{src}: SHA-256 {digest}, expected {SHA256} (another release of the file?)")
    items, atoms, revisions = read_cif(src)
    rows = [tuple(a[k] for k in ATOM_FIELDS) for a in atoms if a["group_PDB"] == "ATOM"]
    if not rows:
        raise SystemExit(f"{src}: no ATOM records")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "4jrd_atoms.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(ATOM_FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(os.path.join(OUT, "4jrd_entry.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("key\tvalue\n")
        for key in ENTRY_FIELDS:
            fh.write(f"{key}\t{items[key]}\n")
        fh.write(f"_pdbx_audit_revision_history.revision_date (latest)\t{max(revisions)}\n")
    print("4jrd_atoms.tsv", len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
