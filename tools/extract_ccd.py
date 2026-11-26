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
* ribonucleotide_components.tsv, ribonucleotide_atoms.tsv, ribonucleotide_bonds.tsv: the standard
  ribonucleotides A, C, G, U and their mirror images, the L-ribonucleotides 0A, 0C, 0G, 0U (the
  reference of MDATK's restraint dictionaries; tools/gen_chem_tables.py). Components: comp_id, name,
  type, formula, formal charge, one-letter code, modification date and InChI. Atoms: every atom with
  element, formal charge, aromatic / leaving flags, CIP label, ordinal, the model and ideal coordinates and
  standard_atom_id / standard_ordinal, the matching atom of the standard residue of the same base
  (for A ... U the atom itself). The match is the bond-graph match of match_atoms, which keeps the
  CCD's names where they are the standard ones; the two hydrogens of a CH2 or NH2 group that the
  CCD names its own way (0A's HN6 / HN6A, H5' / H5'A), which the bonds cannot tell apart, are named
  by their position in the ideal coordinates (name_by_geometry: the sign of the chiral volume at
  the CH2, mirrored for an L component; cis or trans to the ring atom for an NH2): HN6 / HN6A / H5' /
  H5'A are the mirror images of H62 / H61 / H5'' / H5'. (The CCD's model and ideal coordinates do not
  always agree on such pairs; tools/gen_chem_tables.py reports where.) Bonds: atom ids, order,
  aromatic flag and ordinal.
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
                elif headers and headers[0].startswith("_pdbx_chem_comp_descriptor."):
                    toks = [unquote(t) for t in TOKEN.findall(line)]
                    if len(toks) == len(headers):
                        row = {h.split(".", 1)[1]: t for h, t in zip(headers, toks)}
                        if row.get("type") == "InChI":
                            comp["_inchi"] = row.get("descriptor", "")
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


RIBO_STANDARD = ("A", "C", "G", "U")
RIBO_L = {"0A": "A", "0C": "C", "0G": "G", "0U": "U"}


def _sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def _volume(c, a, b, d):
    """Signed volume (a - c) . ((b - c) x (d - c))."""
    u, v, w = _sub(a, c), _sub(b, c), _sub(d, c)
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0]))


def _cis(h, x, a, b) -> bool:
    """True when h and b are on the same side of the x-a bond (|torsion h-x-a-b| < 90 degrees)."""
    xa = _sub(a, x)
    n = sum(t * t for t in xa) ** 0.5
    xa = [t / n for t in xa]

    def perp(p, q):
        d = _sub(p, q)
        k = sum(s * t for s, t in zip(d, xa))
        return [s - k * t for s, t in zip(d, xa)]
    return sum(s * t for s, t in zip(perp(h, x), perp(b, a))) > 0


def coordinates(atoms, kind):
    """atom id -> [x, y, z] of the model ("model_Cartn_") or ideal ("pdbx_model_Cartn_", "_ideal")
    coordinates, or None when an atom has none."""
    pre, suf = ("model_Cartn_", "") if kind == "model" else ("pdbx_model_Cartn_", "_ideal")
    try:
        return {a["atom_id"]: [float(a[pre + k + suf]) for k in "xyz"] for a in atoms}
    except (KeyError, ValueError):
        return None


def name_by_geometry(mapping, own, standard, own_xyz, std_xyz, mirror):
    """Re-pair the hydrogens of each CH2 / NH2 group of a graph match whose names are not the standard
    residue's (0A's HN6 / HN6A, H5' / H5'A) by position in the ideal coordinates: a hydrogen of a CH2 with
    heavy neighbours A, B (standard order) gets the standard name whose signed volume V(X; A, B, H) has the
    same sign (the opposite sign when the component is the mirror image); a hydrogen of an NH2 bonded to
    A gets the standard name that is cis (trans) to the first other heavy neighbour of A. Raises when a
    group's two positions do not differ in that sign."""
    order, element, neighbours = own
    s_order, s_element, s_neighbours = standard
    inverse = {v: k for k, v in mapping.items()}
    out = dict(mapping)
    for x in order:
        hs = sorted((n for n in neighbours[x] if element[n] == "H"), key=order.index)
        if len(hs) != 2 or all(h in s_element for h in hs):
            continue  # a pair that has the standard names keeps them
        sx = mapping[x]
        heavy = sorted((n for n in s_neighbours[sx] if s_element[n] != "H"), key=s_order.index)
        std_h = sorted((n for n in s_neighbours[sx] if s_element[n] == "H"), key=s_order.index)
        if len(heavy) == 2:  # sp3 CH2
            a, b = heavy

            def sig_std(h):
                return _volume(std_xyz[sx], std_xyz[a], std_xyz[b], std_xyz[h]) > 0

            def sig_own(h):
                v = _volume(own_xyz[x], own_xyz[inverse[a]], own_xyz[inverse[b]], own_xyz[h]) > 0
                return v != mirror
        elif len(heavy) == 1:  # NH2 on a ring atom
            a = heavy[0]
            b = sorted((n for n in s_neighbours[a] if n != sx and s_element[n] != "H"), key=s_order.index)[0]

            def sig_std(h):
                return _cis(std_xyz[h], std_xyz[sx], std_xyz[a], std_xyz[b])

            def sig_own(h):
                return _cis(own_xyz[h], own_xyz[x], own_xyz[inverse[a]], own_xyz[inverse[b]])
        else:
            continue
        if sig_std(std_h[0]) == sig_std(std_h[1]) or sig_own(hs[0]) == sig_own(hs[1]):
            raise SystemExit(f"{x}: the two hydrogens are not told apart by their position")
        for h in hs:
            out[h] = next(s for s in std_h if sig_std(s) == sig_own(h))
    return out


def ribonucleotide_rows(blocks):
    """(component rows, atom rows, bond rows) of A C G U and 0A 0C 0G 0U; blocks: comp_id -> (comp, atoms,
    bonds). The L components are checked to be the mirror images of the standard ones (every CIP label
    R <-> S, and the four sugar centres' signed volumes inverted)."""
    comp_rows, atom_rows, bond_rows = [], [], []
    for cid in RIBO_STANDARD + tuple(RIBO_L):
        if cid not in blocks:
            raise SystemExit(f"component {cid} not in the CCD")
        comp, atoms, bonds = blocks[cid]
        comp_rows.append((cid, comp.get("name", ""), comp.get("type", ""), comp.get("formula", ""),
                          comp.get("pdbx_formal_charge", ""), comp.get("one_letter_code", ""),
                          comp.get("pdbx_modified_date", ""), comp.get("_inchi", "")))
        if cid in RIBO_L:
            _, std_atoms, std_bonds = blocks[RIBO_L[cid]]
            own, standard = atom_graph(atoms, bonds), atom_graph(std_atoms, std_bonds)
            mapping = match_atoms(own, standard)
            if mapping is None:
                raise SystemExit(f"{cid}: atoms and bonds do not match {RIBO_L[cid]}")
            own_xyz, std_xyz = coordinates(atoms, "ideal"), coordinates(std_atoms, "ideal")
            mapping = name_by_geometry(mapping, own, standard, own_xyz, std_xyz, mirror=True)

            std_by_id = {a["atom_id"]: a for a in std_atoms}
            for a in atoms:
                s = std_by_id[mapping[a["atom_id"]]]
                pair = {a.get("pdbx_stereo_config", "N"), s.get("pdbx_stereo_config", "N")}
                if pair not in ({"N"}, {"R", "S"}):
                    raise SystemExit(f"{cid} {a['atom_id']}: not the mirror image of {s['atom_id']}")
            for centre, (p, q, r) in (("C4'", ("O4'", "C3'", "C5'")), ("C3'", ("O3'", "C4'", "C2'")),
                                      ("C2'", ("O2'", "C1'", "C3'")), ("C1'", ("O4'", "C2'", "C3'"))):
                inv = {v: k for k, v in mapping.items()}
                v_own = _volume(*(own_xyz[inv[n]] for n in (centre, p, q, r)))
                v_std = _volume(*(std_xyz[n] for n in (centre, p, q, r)))
                if v_own * v_std >= 0:
                    raise SystemExit(f"{cid} {centre}: the ideal coordinates are not the mirror image")
            ordinal = {a["atom_id"]: a["pdbx_ordinal"] for a in std_atoms}
        else:
            mapping = {a["atom_id"]: a["atom_id"] for a in atoms}
            ordinal = {a["atom_id"]: a["pdbx_ordinal"] for a in atoms}
        for a in atoms:
            xyz = [a[k] for k in ("model_Cartn_x", "model_Cartn_y", "model_Cartn_z", "pdbx_model_Cartn_x_ideal",
                                  "pdbx_model_Cartn_y_ideal", "pdbx_model_Cartn_z_ideal")]
            atom_rows.append((cid, a["atom_id"], a["type_symbol"], a["charge"], a["pdbx_aromatic_flag"],
                              a["pdbx_leaving_atom_flag"], a["pdbx_stereo_config"], a["pdbx_ordinal"],
                              *xyz, mapping[a["atom_id"]], ordinal[mapping[a["atom_id"]]]))
        for b in bonds:
            bond_rows.append((cid, b["atom_id_1"], b["atom_id_2"], b["value_order"], b["pdbx_aromatic_flag"],
                              b["pdbx_ordinal"]))
    return comp_rows, atom_rows, bond_rows


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
    ribo = {}
    for comp, atoms, bonds in components(src):
        cid = comp.get("id", "")
        if cid in RIBO_STANDARD or cid in RIBO_L:
            ribo[cid] = (comp, atoms, bonds)
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
    ribo_comps, ribo_atoms, ribo_bonds = ribonucleotide_rows(ribo)
    os.makedirs(OUT, exist_ok=True)
    for name, header, rows in (  # in the components' own order (atom and bond ordinals)
        ("ribonucleotide_components.tsv",
         "comp_id\tname\ttype\tformula\tformal_charge\tone_letter_code\tmodified_date\tinchi", ribo_comps),
        ("ribonucleotide_atoms.tsv",
         "comp_id\tatom_id\telement\tcharge\taromatic_flag\tleaving_flag\tstereo_config\tordinal"
         "\tmodel_x\tmodel_y\tmodel_z\tideal_x\tideal_y\tideal_z\tstandard_atom_id\tstandard_ordinal",
         ribo_atoms),
        ("ribonucleotide_bonds.tsv",
         "comp_id\tatom_id_1\tatom_id_2\tvalue_order\taromatic_flag\tordinal", ribo_bonds),
    ):
        with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(header + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        print(name, len(rows))
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
