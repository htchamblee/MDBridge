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
* l_deoxyribonucleotide_atoms.tsv: every atom of the released L-DNA linking components whose one-letter
  code is A, C, G or T and whose atoms and bonds match those of the standard DNA residue of that code
  (DA, DC, DG, DT; the mirror images 0DA, 0DC, 0DG, 0DT): comp_id, base, atom_id, the matching
  atom of the standard residue (standard_atom_id; a graph match that keeps names where the bonds
  allow, see match_atoms) and element.
* protein_list_residue_atoms.tsv: comp_id, atom_id, alt_atom_id, element and the component's name
  and type of every atom of the released components whose id is a non-standard residue name of
  tools/gen_data_tables.py's protein list (force-field names such as ACE, NME, HYP that the CCD
  also defines). The list is computed from the committed force-field and PDB-format extracts
  (gen_data_tables.py is imported, so the other extracts must already exist).
"""

from __future__ import annotations

import gzip
import importlib.util
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
    """Yield (comp dict, list of atom dicts, list of bond dicts) per data block."""
    comp: dict[str, str] = {}
    atoms: list[dict[str, str]] = []
    bonds: list[dict[str, str]] = []
    headers: list[str] = []
    in_loop = False
    reading_rows = False
    pending_key = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("data_"):
                if comp:
                    yield comp, atoms, bonds
                comp, atoms, bonds, headers = {}, [], [], []
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
                if headers and headers[0].startswith(("_chem_comp_atom.", "_chem_comp_bond.")):
                    toks = [unquote(t) for t in TOKEN.findall(line)]
                    if len(toks) == len(headers):
                        rows = atoms if headers[0].startswith("_chem_comp_atom.") else bonds
                        rows.append({h.split(".", 1)[1]: t for h, t in zip(headers, toks)})
                continue
            if line.startswith("_chem_comp_bond."):
                # A one-bond component lists its bond as key-value pairs, not as a loop.
                parts = line.split(None, 1)
                if not bonds or parts[0] == "_chem_comp_bond.comp_id":
                    bonds.append({})
                bonds[-1][parts[0].split(".", 1)[1]] = unquote(parts[1].strip()) if len(parts) == 2 else ""
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
        yield comp, atoms, bonds


def atom_graph(atoms, bonds):
    """(atom ids in order, element per atom, neighbour set per atom) of a component."""
    order = [a["atom_id"] for a in atoms]
    element = {a["atom_id"]: a["type_symbol"].upper() for a in atoms}
    neighbours: dict[str, set[str]] = {name: set() for name in order}
    for b in bonds:
        one, two = b.get("atom_id_1"), b.get("atom_id_2")
        if one in neighbours and two in neighbours:
            neighbours[one].add(two)
            neighbours[two].add(one)
    return order, element, neighbours


def match_atoms(own, standard):
    """The atoms of a component matched to those of the standard component it mirrors: a map atom id ->
    standard atom id that keeps elements and bonds (a graph isomorphism), or None when there is none.
    Atoms are matched in breadth-first order from the first; each tries the standard atom of the same name
    first, then the others in the standard component's order, so atoms keep their names wherever the
    bonds allow and symmetric atoms (H61 / H62, OP1 / OP2) pair in order."""
    order, element, neighbours = own
    s_order, s_element, s_neighbours = standard
    if len(order) != len(s_order) or sorted(element.values()) != sorted(s_element.values()):
        return None
    sequence, seen = [], set()
    for start in order:
        queue = [start] if start not in seen else []
        seen.update(queue)
        while queue:
            atom = queue.pop(0)
            sequence.append(atom)
            for n in sorted(neighbours[atom], key=order.index):
                if n not in seen:
                    seen.add(n)
                    queue.append(n)
    mapping: dict[str, str] = {}
    used: set[str] = set()

    def extend(k: int) -> bool:
        if k == len(sequence):
            return True
        atom = sequence[k]
        candidates = ([atom] if atom in s_element else []) + [c for c in s_order if c != atom]
        for c in candidates:
            if c in used or s_element[c] != element[atom] or len(s_neighbours[c]) != len(neighbours[atom]):
                continue
            if any((m in neighbours[atom]) != (mapping[m] in s_neighbours[c]) for m in mapping):
                continue
            mapping[atom] = c
            used.add(c)
            if extend(k + 1):
                return True
            del mapping[atom]
            used.discard(c)
        return False

    return dict(mapping) if extend(0) else None


def protein_list() -> set[str]:
    spec = importlib.util.spec_from_file_location(
        "gen_data_tables", os.path.join(ROOT, "tools", "gen_data_tables.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.protein_residues()


def main() -> int:
    src = sys.argv[1]
    with open(os.path.join(ROOT, "data", "sources", "pdb_format", "standard_residues.tsv")) as fh:
        standard = {ln.split("\t")[0] for ln in fh if ln.strip() and not ln.startswith("#")}
    extra = protein_list() - standard
    std_rows, mono_rows, hgcl_rows, extra_rows, ldna_rows = [], [], [], [], []
    d_dna, l_dna = {}, []
    for comp, atoms, bonds in components(src):
        cid = comp.get("id", "")
        if cid in ("DA", "DC", "DG", "DT"):
            d_dna[cid] = atom_graph(atoms, bonds)
        if (comp.get("type", "").upper() == "L-DNA LINKING" and comp.get("pdbx_release_status") == "REL"
                and comp.get("one_letter_code") in ("A", "C", "G", "T")):
            l_dna.append((cid, comp["one_letter_code"], atoms, atom_graph(atoms, bonds)))
        if cid in standard:
            for a in atoms:
                std_rows.append((cid, a["atom_id"], a["alt_atom_id"], a["type_symbol"],
                                 a.get("pdbx_backbone_atom_flag", "?"),
                                 a.get("pdbx_leaving_atom_flag", "?")))
        if comp.get("pdbx_release_status") != "REL":
            continue
        if cid in extra:
            for a in atoms:
                extra_rows.append((cid, a["atom_id"], a["alt_atom_id"], a["type_symbol"],
                                   comp.get("name", ""), comp.get("type", "")))
        if len(atoms) == 1:
            a = atoms[0]
            mono_rows.append((cid, a["atom_id"], a["type_symbol"], a["charge"], comp.get("name", "")))
        for a in atoms:
            if (a["atom_id"], a["type_symbol"].upper()) in (("HG", "HG"), ("CL", "CL")):
                hgcl_rows.append((cid, a["atom_id"], a["type_symbol"]))
    for cid, letter, atoms, graph in l_dna:  # the mirror images of DA DC DG DT
        mapping = match_atoms(graph, d_dna["D" + letter]) if "D" + letter in d_dna else None
        if mapping is None:
            continue
        for a in atoms:
            ldna_rows.append((cid, letter, a["atom_id"], mapping[a["atom_id"]], a["type_symbol"]))
    os.makedirs(OUT, exist_ok=True)
    for name, header, rows in (
        ("standard_residue_atoms.tsv", "comp_id\tatom_id\talt_atom_id\telement\tbackbone_flag\tleaving_flag", std_rows),
        ("monatomic_components.tsv", "comp_id\tatom_id\telement\tcharge\tname", mono_rows),
        ("hg_cl_atom_names.tsv", "comp_id\tatom_id\telement", hgcl_rows),
        ("protein_list_residue_atoms.tsv", "comp_id\tatom_id\talt_atom_id\telement\tname\ttype", extra_rows),
        ("l_deoxyribonucleotide_atoms.tsv", "comp_id\tbase\tatom_id\tstandard_atom_id\telement", ldna_rows),
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
