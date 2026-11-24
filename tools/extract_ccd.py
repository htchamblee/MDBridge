"""Cut the small extracts under data/sources/ccd/ out of the wwPDB Chemical Component Dictionary.

    curl -O https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
    python3 -I tools/extract_ccd.py components.cif.gz

The CCD is released by the wwPDB under CC0 1.0. Only the fields that tools/gen_data_tables.py
needs are kept:

* standard_residue_atoms.tsv: comp_id, atom_id, alt_atom_id (the PDB v2.3 spelling), element and
  the CCD backbone and leaving-atom flags of every atom of the standard residues listed in
  data/sources/pdb_format/standard_residues.tsv.
* monatomic_components.tsv: every released component made of exactly one atom
  (comp_id, atom_id, element, formal charge, name).
* hg_cl_atom_names.tsv: every released component with an atom named HG whose element is Hg, or
  an atom named CL whose element is Cl (comp_id, atom_id, element).
"""

from __future__ import annotations

import gzip
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "ccd")
TOKEN = re.compile(r"""'(?:[^']|'(?!\s|$))*'(?=\s|$)|"(?:[^"]|"(?!\s|$))*"(?=\s|$)|\S+""")


def unquote(tok: str) -> str:
    if len(tok) >= 2 and tok[0] == tok[-1] and tok[0] in "'\"":
        return tok[1:-1]
    return tok


def components(path: str):
    """Yield (comp dict, list of atom dicts) per data block."""
    comp: dict[str, str] = {}
    atoms: list[dict[str, str]] = []
    headers: list[str] = []
    in_loop = False
    reading_rows = False
    pending_key = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("data_"):
                if comp:
                    yield comp, atoms
                comp, atoms, headers = {}, [], []
                in_loop = reading_rows = False
                pending_key = None
                continue
            if line.startswith("#"):
                in_loop = reading_rows = False
                continue
            if line.startswith("loop_"):
                in_loop, reading_rows, headers = True, False, []
                continue
            if in_loop and not reading_rows and line.startswith("_"):
                headers.append(line.split()[0])
                continue
            if in_loop:
                reading_rows = True
                if headers and headers[0].startswith("_chem_comp_atom."):
                    toks = [unquote(t) for t in TOKEN.findall(line)]
                    if len(toks) == len(headers):
                        atoms.append({h.split(".", 1)[1]: t for h, t in zip(headers, toks)})
                continue
            if line.startswith("_chem_comp_atom."):
                # A one-atom component lists its atom as key-value pairs, not as a loop.
                parts = line.split(None, 1)
                if not atoms or parts[0] == "_chem_comp_atom.comp_id":
                    atoms.append({})
                atoms[-1][parts[0].split(".", 1)[1]] = unquote(parts[1].strip()) if len(parts) == 2 else ""
                continue
            if line.startswith("_chem_comp."):
                parts = line.split(None, 1)
                key = parts[0].split(".", 1)[1]
                if len(parts) == 2:
                    comp[key] = unquote(parts[1].strip())
                else:
                    pending_key = key
                continue
            if pending_key is not None:
                if line.startswith(";"):
                    if pending_key not in comp:
                        comp[pending_key] = line[1:]
                    else:
                        pending_key = None
                else:
                    comp[pending_key] = unquote(line.strip())
                    pending_key = None
    if comp:
        yield comp, atoms


def main() -> int:
    src = sys.argv[1]
    with open(os.path.join(ROOT, "data", "sources", "pdb_format", "standard_residues.tsv")) as fh:
        standard = {ln.split("\t")[0] for ln in fh if ln.strip() and not ln.startswith("#")}
    std_rows, mono_rows, hgcl_rows = [], [], []
    for comp, atoms in components(src):
        cid = comp.get("id", "")
        if cid in standard:
            for a in atoms:
                std_rows.append((cid, a["atom_id"], a["alt_atom_id"], a["type_symbol"],
                                 a.get("pdbx_backbone_atom_flag", "?"),
                                 a.get("pdbx_leaving_atom_flag", "?")))
        if comp.get("pdbx_release_status") != "REL":
            continue
        if len(atoms) == 1:
            a = atoms[0]
            mono_rows.append((cid, a["atom_id"], a["type_symbol"], a["charge"], comp.get("name", "")))
        for a in atoms:
            if (a["atom_id"], a["type_symbol"].upper()) in (("HG", "HG"), ("CL", "CL")):
                hgcl_rows.append((cid, a["atom_id"], a["type_symbol"]))
    os.makedirs(OUT, exist_ok=True)
    for name, header, rows in (
        ("standard_residue_atoms.tsv", "comp_id\tatom_id\talt_atom_id\telement\tbackbone_flag\tleaving_flag", std_rows),
        ("monatomic_components.tsv", "comp_id\tatom_id\telement\tcharge\tname", mono_rows),
        ("hg_cl_atom_names.tsv", "comp_id\tatom_id\telement", hgcl_rows),
    ):
        with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
            fh.write(header + "\n")
            for r in sorted(set(rows)):
                fh.write("\t".join(r) + "\n")
        print(name, len(set(rows)))
    missing = standard - {r[0] for r in std_rows}
    if missing:
        print("standard residues not in the CCD:", sorted(missing))
    return 0


if __name__ == "__main__":
    sys.exit(main())
