"""Write the chemistry reference tables of src/chem from the public-source extracts in data/sources.

    python3 -I -B tools/gen_chem_tables.py MDATK_CHECKOUT

Writes, deterministically and only from the committed extracts:

* src/chem/helix_tables.inc: the Arnott A-RNA fibre monomer (data/sources/helix/arnott_arna.tsv,
  Arnott, Hukins & Dover 1972 as tabulated in AmberTools' nucgen.dat `$ARNA`), heavy atoms in the
  order of the CCD's standard residues, and the 4JRD nucleotide (chain A residue 10 of PDB 4JRD,
  data/sources/pdb/4jrd_atoms.tsv) in the frame of the screw that maps the duplex one level on,
  least squares over residues A3-A11 / B2-B10 (every atom present in both residues of a step),
  with its twist and rise. Frame: z along the screw axis, strand 1 (chain A) running 5' -> 3'
  towards -z, x towards the C1' atom of A10 (phi 0, z 0 there). Derived numbers are rounded to
  1e-6 (Angstrom, degrees), far below their meaning, so that the last bits of libm do not reach
  the file.
* src/chem/component_tables.inc: the L-ribonucleotides 0A 0C 0G 0U of the wwPDB CCD
  (data/sources/ccd/ribonucleotide_*.tsv) in the form MDATK's restraint dictionaries use: heavy
  atoms, bonds and hydrogens under the names of the standard residue (H61 / H62, H5' / H5'' ...,
  the CCD's 0A names HN6 / HN6A / H5'A mapped by tools/extract_ccd.py), in the standard residue's
  order, with the phosphate deprotonated (each phosphate oxygen that carries a hydrogen, OP2 and
  OP3, loses it and is charged -1); the positions that tell the two hydrogens of a CH2 / NH2
  group apart (from the standard residues' ideal coordinates; a CH2 mirrored for the L form); and
  the signed chiral volumes of the L-ribose centres (from the components' ideal coordinates).

Then the tables are compared with MDATK 0.0.4.1's embedded ones (builder.py `_BACKBONE`, `_BASES`,
`_H_RULES`, `_JRD_NUCLEOTIDE`, `JRD_TWIST`, `JRD_RISE`; dictionaries.py `COMPONENTS`,
`L_CHIRALITY`), read from MDATK_CHECKOUT (builder.py through its syntax tree; dictionaries.py,
standard library only, imported from its file). Every difference is printed; one that is not in
DOCUMENTED (each a docs/DIFFERENCES.md entry: MDBridge keeps the primary source) fails the run,
and so does a DOCUMENTED difference that no longer occurs. The CCD's own internal disagreements
(model against ideal coordinates) are printed as notes.

Pure Python (no numpy): the superpositions are Horn's quaternion method with a Jacobi
eigensolver, so the output does not depend on a BLAS.
"""

from __future__ import annotations

import ast
import csv
import importlib.util
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCES = os.path.join(ROOT, "data", "sources")
HELIX_INC = os.path.join(ROOT, "src", "chem", "helix_tables.inc")
COMPONENT_INC = os.path.join(ROOT, "src", "chem", "component_tables.inc")

BASES = ("A", "C", "G", "U")
L_CODES = {"0A": "A", "0C": "C", "0G": "G", "0U": "U"}
GLYCOSIDIC_N = {"A": "N9", "G": "N9", "C": "N1", "U": "N1"}
# MDATK's chirality descriptors (atom-name conventions; the signs are computed): (centre, three
# neighbours, "N" = the glycosidic nitrogen); the fourth neighbour is the centre's hydrogen
CHIRAL_CENTRES = (("C4'", ("O4'", "C3'", "C5'")), ("C3'", ("O3'", "C4'", "C2'")),
                  ("C2'", ("O2'", "C1'", "C3'")), ("C1'", ("O4'", "N", "C2'")))
JRD_REFERENCE = ("A", 10)  # chain, residue of the 4JRD nucleotide
JRD_SCREW = (("A", 3, 11), ("B", 2, 10))  # steps (chain, first, last residue) of the screw fit
ROUND = 6

# Differences from MDATK 0.0.4.1 that MDBridge makes on purpose (docs/DIFFERENCES.md): key -> entry id.
DOCUMENTED = {
    "hpos G H21": "CHD-01",
    "hpos G H22": "CHD-01",
    "hpos 0G H21": "CHD-01",
    "hpos 0G H22": "CHD-01",
    "builder H21": "CHD-01",
    "builder H22": "CHD-01",
    "jrd frame": "CHD-02",
    "jrd values": "CHD-02",
    "jrd twist": "CHD-02",
    "jrd rise": "CHD-02",
}


# --- small vector algebra -------------------------------------------------------------------------

def sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def add(a, b):
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def scale(a, s):
    return [a[0] * s, a[1] * s, a[2] * s]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    return scale(a, 1.0 / norm(a))


def volume(c, a, b, d):
    """Signed volume (a - c) . ((b - c) x (d - c))."""
    return dot(sub(a, c), cross(sub(b, c), sub(d, c)))


def torsion(a, b, c, d):
    b0, b1, b2 = sub(a, b), unit(sub(c, b)), sub(d, c)
    v = sub(b0, scale(b1, dot(b0, b1)))
    w = sub(b2, scale(b1, dot(b2, b1)))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))


def mat_vec(m, v):
    return [dot(m[0], v), dot(m[1], v), dot(m[2], v)]


def jacobi_max_eigenvector(n):
    """Eigenvector of the largest eigenvalue of a symmetric 4 x 4 matrix (cyclic Jacobi rotations)."""
    a = [row[:] for row in n]
    v = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    for _ in range(100):
        off = sum(a[i][j] ** 2 for i in range(4) for j in range(4) if i != j)
        if off < 1e-30:
            break
        for p in range(3):
            for q in range(p + 1, 4):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(4):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(4):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(4):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    best = max(range(4), key=lambda i: a[i][i])
    return [v[k][best] for k in range(4)]


def superpose(p, q):
    """Least-squares rigid motion x -> R x + t taking the points p onto q (Horn 1987): (quaternion
    (w, x, y, z), R, t)."""
    n = len(p)
    pc = scale([sum(x[i] for x in p) for i in range(3)], 1.0 / n)
    qc = scale([sum(x[i] for x in q) for i in range(3)], 1.0 / n)
    s = [[0.0] * 3 for _ in range(3)]
    for a, b in zip(p, q):
        da, db = sub(a, pc), sub(b, qc)
        for i in range(3):
            for j in range(3):
                s[i][j] += da[i] * db[j]
    (sxx, sxy, sxz), (syx, syy, syz), (szx, szy, szz) = s
    m = [[sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
         [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
         [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
         [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz]]
    w, x, y, z = jacobi_max_eigenvector(m)
    if w < 0:
        w, x, y, z = -w, -x, -y, -z
    r = [[w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
         [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
         [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z]]
    return (w, x, y, z), r, sub(qc, mat_vec(r, pc))


def rounded(v):
    return float(f"{v:.{ROUND}f}") + 0.0  # + 0.0: no "-0.0"


# --- sources --------------------------------------------------------------------------------------

def read_tsv(*parts):
    with open(os.path.join(SOURCES, *parts), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ccd():
    atoms, bonds = read_tsv("ccd", "ribonucleotide_atoms.tsv"), read_tsv("ccd", "ribonucleotide_bonds.tsv")
    comps = {r["comp_id"]: r for r in read_tsv("ccd", "ribonucleotide_components.tsv")}
    return comps, atoms, bonds


def xyz(row, kind):
    return [float(row[f"{kind}_{k}"]) for k in "xyz"]


def standard_order():
    """comp -> [atom ids in CCD ordinal order] of A C G U."""
    _, atoms, _ = ccd()
    out = {}
    for b in BASES:
        rows = sorted((r for r in atoms if r["comp_id"] == b), key=lambda r: int(r["ordinal"]))
        out[b] = [r["atom_id"] for r in rows]
    return out


def arnott():
    """(backbone, {base: atoms}) of the Arnott monomer: (PDB name, r, phi, z), heavy atoms in the CCD
    standard residues' order (A's backbone without OP3)."""
    rows = read_tsv("helix", "arnott_arna.tsv")
    order = standard_order()
    _, atoms, _ = ccd()
    element = {(r["comp_id"], r["atom_id"]): r["element"] for r in atoms}

    def pick(sets, comp):
        got = {r["pdb_name"]: (r["pdb_name"], float(r["r_angstrom"]), float(r["phi_degrees"]),
                               float(r["z_angstrom"])) for r in rows if r["set"] in sets}
        return [got[n] for n in order[comp] if n in got and element[(comp, n)] != "H"]
    backbone = pick({"backbone"}, "A")
    bases = {b: pick({b}, b) for b in BASES}
    return backbone, bases


def jrd_atoms():
    """(chain, residue) -> {atom name: [x, y, z]} of the ATOM records of 4JRD (first altloc)."""
    out = {}
    for r in read_tsv("pdb", "4jrd_atoms.tsv"):
        if r["label_alt_id"] not in (".", "A") or r["pdbx_PDB_model_num"] != "1":
            continue
        out.setdefault((r["auth_asym_id"], int(r["auth_seq_id"])), {}).setdefault(
            r["label_atom_id"], [float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])])
    return out


def jrd_screw(res):
    """(origin, unit axis, rotation angle about it in degrees, rise) of the least-squares screw that maps
    each residue of JRD_SCREW onto the next."""
    p, q = [], []
    for chain, first, last in JRD_SCREW:
        for k in range(first, last):
            a, b = res[(chain, k)], res[(chain, k + 1)]
            for name in a:
                if name in b:
                    p.append(a[name])
                    q.append(b[name])
    (w, x, y, z), r, t = superpose(p, q)
    v = [x, y, z]
    angle = 2.0 * math.atan2(norm(v), w)
    axis = unit(v)
    rise = dot(axis, t)
    t_perp = sub(t, scale(axis, rise))
    origin = scale(add(t_perp, scale(cross(axis, t_perp), 1.0 / math.tan(angle / 2.0))), 0.5)
    return origin, axis, math.degrees(angle), rise


def jrd_nucleotide():
    """(twist, rise, [(name, r, phi, z)]) of the 4JRD nucleotide in its helix frame (see the module doc);
    also the frame (origin, x, y, z axes in 4JRD coordinates) for the checks."""
    res = jrd_atoms()
    origin, axis, angle, rise = jrd_screw(res)
    first, last = res[(JRD_SCREW[0][0], JRD_SCREW[0][1])]["C1'"], res[(JRD_SCREW[0][0], JRD_SCREW[0][2])]["C1'"]
    ez = axis if dot(axis, sub(last, first)) < 0 else scale(axis, -1.0)  # strand 1 5'->3' towards -z
    # the screw about +z: a turn of -twist per level (MDATK's _Screw convention, right-handed)
    twist = -angle if ez == axis else angle
    rise_z = -dot(ez, scale(axis, rise))  # moving -rise along z per level
    if twist <= 0 or rise_z <= 0:
        raise SystemExit(f"4JRD screw: twist {twist}, rise {rise_z} (expected a right-handed helix)")
    ref = res[JRD_REFERENCE]
    c1 = sub(ref["C1'"], origin)
    ex = unit(sub(c1, scale(ez, dot(c1, ez))))
    ey = cross(ez, ex)
    z0 = dot(c1, ez)
    order = standard_order()["A"]
    out = []
    for name in order:
        if name not in ref:
            continue
        d = sub(ref[name], origin)
        out.append((name, rounded(math.hypot(dot(d, ex), dot(d, ey))),
                    rounded(math.degrees(math.atan2(dot(d, ey), dot(d, ex)))), rounded(dot(d, ez) - z0)))
    return rounded(twist), rounded(rise_z), out, (origin, ex, ey, ez)


def duplex_rmsd(nucleotide, twist, rise, to_pdb):
    """RMSD (A) of the rigid duplex the nucleotide makes with its screw (one level on: -twist about z,
    -rise along z; the partner strand by the two-fold rotation about z) from 4JRD's residues A3-A11 /
    B2-B10 (the reference nucleotide is A10; B k pairs with A k+1); to_pdb maps frame -> 4JRD coordinates."""
    res = jrd_atoms()
    sq = []
    for k in range(JRD_SCREW[0][1], JRD_SCREW[0][2] + 1):
        level = k - JRD_REFERENCE[1]
        t = math.radians(-twist * level)
        for name, (x, y, z) in nucleotide.items():
            a = [math.cos(t) * x - math.sin(t) * y, math.sin(t) * x + math.cos(t) * y, z - rise * level]
            b = [-a[0], -a[1], a[2]]
            sq.append(dot(sub(to_pdb(a), res[("A", k)][name]), sub(to_pdb(a), res[("A", k)][name])))
            if name in res[("B", k - 1)]:
                d = sub(to_pdb(b), res[("B", k - 1)][name])
                sq.append(dot(d, d))
    return math.sqrt(sum(sq) / len(sq))


def cylindrical(r, phi, z):
    p = math.radians(phi)
    return [r * math.cos(p), r * math.sin(p), z]


# --- the dictionary components --------------------------------------------------------------------

def components():
    """{code: dict(name, base, atoms [(name, element, charge)], bonds [(a, b, order)], hydrogens
    [(heavy, (H ...))], hpos [(H, X, A, B, kind)], chirality [(centre, (a, b, c), sign)])} of 0A 0C 0G 0U
    in the form of MDATK's restraint dictionaries (module doc)."""
    comps, atoms, bonds = ccd()
    std_bond_rank = {}
    for b in bonds:
        if b["comp_id"] in BASES:
            std_bond_rank[(b["comp_id"], b["atom_id_1"], b["atom_id_2"])] = int(b["ordinal"])
    rules = hydrogen_rules()
    out = {}
    for code, base in L_CODES.items():
        rows = [r for r in atoms if r["comp_id"] == code]
        name_of = {r["atom_id"]: r["standard_atom_id"] for r in rows}
        rank = {r["standard_atom_id"]: int(r["standard_ordinal"]) for r in rows}
        element = {r["standard_atom_id"]: r["element"] for r in rows}
        charge = {r["standard_atom_id"]: int(r["charge"]) for r in rows}
        neighbours = {n: [] for n in element}
        own_bonds = []
        for b in bonds:
            if b["comp_id"] != code:
                continue
            a1, a2 = name_of[b["atom_id_1"]], name_of[b["atom_id_2"]]
            order = {"SING": 1, "DOUB": 2}[b["value_order"]]
            neighbours[a1].append(a2)
            neighbours[a2].append(a1)
            if (base, a1, a2) not in std_bond_rank:
                a1, a2 = a2, a1  # the standard residue's orientation
            own_bonds.append((std_bond_rank[(base, a1, a2)], a1, a2, order))
        std_orders = {(b["atom_id_1"], b["atom_id_2"]): b["value_order"] for b in bonds if b["comp_id"] == base}
        for _, a1, a2, order in own_bonds:
            if {"SING": 1, "DOUB": 2}[std_orders[(a1, a2)]] != order:
                raise SystemExit(f"{code} {a1}-{a2}: bond order differs from {base}")
        # deprotonated phosphate: a phosphate oxygen's hydrogen goes, the oxygen is charged -1
        removed = set()
        for o in element:
            if element[o] == "O" and any(element[n] == "P" for n in neighbours[o]):
                hs = [n for n in neighbours[o] if element[n] == "H"]
                if hs:
                    removed.update(hs)
                    charge[o] -= len(hs)
        heavy = sorted((n for n in element if element[n] != "H"), key=rank.get)
        hydrogens = sorted((n for n in element if element[n] == "H" and n not in removed), key=rank.get)
        parent = {h: next(n for n in neighbours[h] if element[n] != "H") for h in hydrogens}
        groups = []
        for h in hydrogens:
            if groups and groups[-1][0] == parent[h]:
                groups[-1][1].append(h)
            elif any(g[0] == parent[h] for g in groups):
                raise SystemExit(f"{code}: the hydrogens of {parent[h]} are not consecutive")
            else:
                groups.append((parent[h], [h]))
        ideal = {r["standard_atom_id"]: xyz(r, "ideal") for r in rows}
        chirality = []
        for centre, nbrs in CHIRAL_CENTRES:
            names = [GLYCOSIDIC_N[base] if n == "N" else n for n in nbrs]
            v = volume(ideal[centre], *(ideal[n] for n in names))
            chirality.append((centre, tuple(nbrs), 1 if v > 0 else -1))
        out[code] = dict(
            name=comps[code]["name"], base=base,
            atoms=[(n, element[n], charge[n]) for n in heavy],
            bonds=[(a1, a2, order) for _, a1, a2, order in sorted(own_bonds)
                   if element[a1] != "H" and element[a2] != "H"],
            hydrogens=[(x, tuple(hs)) for x, hs in groups],
            hpos=rules[code], chirality=chirality)
    return out


def hydrogen_rules():
    """{comp: [(H, X, A, B, kind)]} for A C G U and 0A 0C 0G 0U: how the two hydrogens of each CH2 and NH2
    group of the standard residue are told apart, from the standard residues' ideal coordinates. CH2 with
    heavy neighbours A, B (standard order): kind '+' / '-', the sign of the chiral volume V(X; A, B, H),
    inverted for the L components (mirror images). NH2 on X bonded to A: B is A's first ring nitrogen
    (standard order) other than X; kind 'cis' / 'trans' of the torsion H-X-A-B."""
    _, atoms, bonds = ccd()
    out = {}
    for base in BASES:
        rows = {r["atom_id"]: r for r in atoms if r["comp_id"] == base}
        rank = {n: int(r["ordinal"]) for n, r in rows.items()}
        el = {n: r["element"] for n, r in rows.items()}
        nb = {n: [] for n in rows}
        for b in bonds:
            if b["comp_id"] == base:
                nb[b["atom_id_1"]].append(b["atom_id_2"])
                nb[b["atom_id_2"]].append(b["atom_id_1"])
        ideal = {n: xyz(r, "ideal") for n, r in rows.items()}
        d_rules, l_rules = [], []
        for x in sorted(rows, key=rank.get):
            hs = sorted((n for n in nb[x] if el[n] == "H"), key=rank.get)
            heavy = sorted((n for n in nb[x] if el[n] != "H"), key=rank.get)
            if len(hs) != 2:
                continue
            for h in hs:
                if el[x] == "C" and len(heavy) == 2:
                    a, b = heavy
                    sign = volume(ideal[x], ideal[a], ideal[b], ideal[h]) > 0
                    d_rules.append((h, x, a, b, "+" if sign else "-"))
                    l_rules.append((h, x, a, b, "-" if sign else "+"))
                elif el[x] == "N" and len(heavy) == 1:
                    a = heavy[0]
                    b = sorted((n for n in nb[a] if n != x and el[n] == "N"), key=rank.get)[0]
                    kind = "cis" if abs(torsion(ideal[h], ideal[x], ideal[a], ideal[b])) < 90 else "trans"
                    d_rules.append((h, x, a, b, kind))
                    l_rules.append((h, x, a, b, kind))
        out[base] = d_rules
        out[next(c for c, b in L_CODES.items() if b == base)] = l_rules
    return out


def ccd_notes():
    """Where the CCD's model coordinates, or a component's own ideal coordinates, contradict the rules of
    hydrogen_rules() (informational)."""
    _, atoms, _ = ccd()
    rules = hydrogen_rules()
    notes = []
    for comp, rs in rules.items():
        rows = [r for r in atoms if r["comp_id"] == comp]
        by_std = {r["standard_atom_id"]: r for r in rows}
        for kind_xyz in ("model", "ideal"):
            try:
                pos = {n: xyz(r, kind_xyz) for n, r in by_std.items()}
            except ValueError:
                continue
            for h, x, a, b, kind in rs:
                if kind in "+-":
                    got = "+" if volume(pos[x], pos[a], pos[b], pos[h]) > 0 else "-"
                else:
                    got = "cis" if abs(torsion(pos[h], pos[x], pos[a], pos[b])) < 90 else "trans"
                if got != kind:
                    notes.append(f"CCD {comp} {by_std[h]['atom_id']} ({h}): its {kind_xyz} coordinates "
                                 f"give {got}, the standard residue's ideal ones {kind}")
    return notes


# --- C++ output -----------------------------------------------------------------------------------

def cpp_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def cpp_num(v):
    return repr(float(v))


def helix_inc():
    backbone, bases = arnott()
    twist, rise, jrd, _ = jrd_nucleotide()
    lines = [
        "// Generated by tools/gen_chem_tables.py from data/sources/helix/arnott_arna.tsv (the Arnott A-RNA",
        "// fibre monomer, Arnott, Hukins & Dover 1972, as tabulated in AmberTools' nucgen.dat $ARNA) and",
        "// data/sources/pdb/4jrd_atoms.tsv (PDB 4JRD, CC0). Do not edit.",
        "// Included inside a namespace that defines CylindricalAtom {name, r (Angstrom), phi (degrees),",
        "// z (Angstrom)} and CylindricalBaseAtom {base, name, r, phi, z}.",
        "",
        "// Arnott A-RNA monomer, helix frame (D sugar, strand 1, level 0): backbone heavy atoms",
        "inline constexpr CylindricalAtom kArnottBackbone[] = {",
    ]
    lines += [f"    {{{cpp_str(n)}, {cpp_num(r)}, {cpp_num(p)}, {cpp_num(z)}}}," for n, r, p, z in backbone]
    lines += ["};", "", "// base heavy atoms, per base in the order A C G U, atoms in the CCD residue's order",
              "inline constexpr CylindricalBaseAtom kArnottBases[] = {"]
    for b in BASES:
        lines += [f"    {{'{b}', {cpp_str(n)}, {cpp_num(r)}, {cpp_num(p)}, {cpp_num(z)}}},"
                  for n, r, p, z in bases[b]]
    lines += ["};", "",
              f"// PDB 4JRD chain {JRD_REFERENCE[0]} residue {JRD_REFERENCE[1]} in the frame of the duplex's screw "
              "(least squares over",
              "// A3-A11 / B2-B10): z along the axis, strand 1 5'->3' towards -z, x through C1' (phi 0, z 0);",
              "// one level on: turn by -kJrdTwist about z, move -kJrdRise along z",
              f"inline constexpr double kJrdTwist = {cpp_num(twist)};  // degrees per level",
              f"inline constexpr double kJrdRise = {cpp_num(rise)};  // Angstrom per level",
              "inline constexpr CylindricalAtom kJrdNucleotide[] = {"]
    lines += [f"    {{{cpp_str(n)}, {cpp_num(r)}, {cpp_num(p)}, {cpp_num(z)}}}," for n, r, p, z in jrd]
    lines += ["};"]
    return "\n".join(lines) + "\n"


def component_inc():
    comps = components()
    rules = hydrogen_rules()
    lines = [
        "// Generated by tools/gen_chem_tables.py from data/sources/ccd/ribonucleotide_*.tsv (wwPDB Chemical",
        "// Component Dictionary, CC0). Do not edit.",
        "// Included inside a namespace that defines ComponentInfo {code, base, name}, ComponentAtom {code,",
        "// name, element, charge}, ComponentBond {code, atom 1, atom 2, order}, ComponentHydrogen {code,",
        "// heavy atom, hydrogen}, HydrogenPosition {residue, H, X, A, B, kind} and ChiralCentre {centre, a, b,",
        "// c, sign}.",
        "",
        "// the CCD's L-ribonucleotides (restraint-dictionary form: standard atom names, phosphate",
        "// deprotonated)",
        "inline constexpr ComponentInfo kComponents[] = {",
    ]
    lines += [f"    {{{cpp_str(c)}, '{d['base']}', {cpp_str(d['name'])}}}," for c, d in comps.items()]
    lines += ["};", "", "// heavy atoms in the standard residue's order: {code, name, element, formal charge}",
              "inline constexpr ComponentAtom kComponentAtoms[] = {"]
    for c, d in comps.items():
        lines += [f"    {{{cpp_str(c)}, {cpp_str(n)}, {cpp_str(e)}, {q}}}," for n, e, q in d["atoms"]]
    lines += ["};", "", "// bonds between heavy atoms in the standard residue's order: {code, atom 1, atom 2, "
              "order (1, 2)}", "inline constexpr ComponentBond kComponentBonds[] = {"]
    for c, d in comps.items():
        lines += [f"    {{{cpp_str(c)}, {cpp_str(a)}, {cpp_str(b)}, {o}}}," for a, b, o in d["bonds"]]
    lines += ["};", "", "// hydrogens by heavy atom, in the standard residue's order: {code, heavy atom, hydrogen}",
              "inline constexpr ComponentHydrogen kComponentHydrogens[] = {"]
    for c, d in comps.items():
        lines += [f"    {{{cpp_str(c)}, {cpp_str(x)}, {cpp_str(h)}}}," for x, hs in d["hydrogens"] for h in hs]
    lines += ["};", "",
              "// the two hydrogens of a CH2 / NH2 group: {residue (A C G U: D; 0A ...: L), H, X, A, B, kind}; kind",
              "// \"+\" / \"-\": sign of the chiral volume V(X; A, B, H); \"cis\" / \"trans\": torsion H-X-A-B",
              "inline constexpr HydrogenPosition kHydrogenPositions[] = {"]
    for comp, rs in rules.items():
        lines += [f"    {{{cpp_str(comp)}, {cpp_str(h)}, {cpp_str(x)}, {cpp_str(a)}, {cpp_str(b)}, {cpp_str(k)}}},"
                  for h, x, a, b, k in rs]
    lines += ["};", "",
              "// L-ribose: sign of the chiral volume V(centre; a, b, c) (\"N\": N9 of a purine, N1 of a pyrimidine);",
              "// the centre's fourth neighbour is its hydrogen. D-ribose has the opposite signs.",
              "inline constexpr ChiralCentre kLChirality[] = {"]
    first = next(iter(comps.values()))["chirality"]
    for c, d in comps.items():
        if d["chirality"] != first:
            raise SystemExit(f"{c}: chirality signs differ between the components")
    lines += [f"    {{{cpp_str(c)}, {cpp_str(a)}, {cpp_str(b)}, {cpp_str(e)}, {s}}}," for c, (a, b, e), s in first]
    lines += ["};"]
    return "\n".join(lines) + "\n"


# --- MDATK's embedded tables ----------------------------------------------------------------------

def mdatk_constants(path, names):
    tree = ast.parse(open(path, encoding="utf-8").read())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in names:
                out[node.targets[0].id] = ast.literal_eval(node.value)
    return out


def mdatk_dictionaries(path):
    spec = importlib.util.spec_from_file_location("mdatk_dictionaries", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["mdatk_dictionaries"] = module  # dataclasses look the module up
    spec.loader.exec_module(module)
    return module


def compare(mdatk):
    """[(key, message)] of every difference between MDATK's tables and these."""
    diffs = []
    builder = mdatk_constants(os.path.join(mdatk, "mdatk", "builder.py"),
                              {"_BACKBONE", "_BASES", "_H_RULES", "_JRD_NUCLEOTIDE", "JRD_TWIST", "JRD_RISE"})
    backbone, bases = arnott()
    if tuple(backbone) != tuple(builder["_BACKBONE"]):
        diffs.append(("arnott backbone", f"builder.py _BACKBONE {builder['_BACKBONE']} != {backbone}"))
    for b in BASES:
        if tuple(bases[b]) != tuple(builder["_BASES"][b]):
            diffs.append((f"arnott {b}", f"builder.py _BASES[{b}] {builder['_BASES'][b]} != {bases[b]}"))
    # OA / OB -> OP1 / OP2: the same chirality at P as every phosphate of PDB 4JRD (an RNA polymer
    # named by the PDB's rule; the CCD's ideal monomer coordinates do not follow it, its model ones do)
    _, atoms, _ = ccd()
    arn = {n: cylindrical(r, p, z) for n, r, p, z in backbone}
    v_arn = volume(arn["P"], arn["O5'"], arn["OP1"], arn["OP2"])
    v_pdb = [volume(r["P"], r["O5'"], r["OP1"], r["OP2"]) for r in jrd_atoms().values() if "P" in r]
    if not all(v * v_arn > 0 for v in v_pdb):
        diffs.append(("arnott OP1/OP2", f"OA / OB as OP1 / OP2 give V(P; O5', OP1, OP2) {v_arn:.2f}, "
                                        f"4JRD {min(v_pdb):.2f} .. {max(v_pdb):.2f}"))
    print(f"Arnott: OA / OB as OP1 / OP2 give V(P; O5', OP1, OP2) {v_arn:.2f}; the {len(v_pdb)} phosphates "
          f"of 4JRD {min(v_pdb):.2f} .. {max(v_pdb):.2f}")
    # the builder's hydrogen rules against the standard residues' positions
    rules = hydrogen_rules()
    std_res = {b: {r["atom_id"]: r for r in atoms if r["comp_id"] == b} for b in BASES}
    for base in BASES:
        for h, x, a, b, kind in rules[base]:
            for name, parent, rule, refs in builder["_H_RULES"]["sugar"] + builder["_H_RULES"][base]:
                if name != h:
                    continue
                if rule == "nh2":  # placed parallel to refs[1] -> refs[0]: cis to the other ring neighbour
                    other = next(n for n in (r_ for r_ in nb_heavy(std_res[base], parent=refs[0], bonds_of=base))
                                 if n not in (parent, refs[1]))
                    placed = "cis" if other == b else "trans"
                    if placed != kind:
                        diffs.append((f"builder {h}", f"builder.py _H_RULES places {h} ({base}) {placed} to {b} "
                                                       f"(torsion {h}-{x}-{a}-{b}); the CCD's {base}: {kind}"))
                elif rule == "ch2" and kind != ("+" if h == "H5'" else "-"):
                    diffs.append((f"builder {h}", f"builder.py names {h} by a positive volume; the CCD: {kind}"))
    # 4JRD
    twist, rise, jrd, (origin, ex, ey, ez) = jrd_nucleotide()
    md = {n: (r, p, z) for n, r, p, z in builder["_JRD_NUCLEOTIDE"]}
    res = jrd_atoms()[JRD_REFERENCE]
    if [n for n, *_ in jrd] != list(md):
        diffs.append(("jrd names", f"_JRD_NUCLEOTIDE atoms {list(md)} != {[n for n, *_ in jrd]}"))
    else:
        p = [cylindrical(*md[n]) for n, *_ in jrd]
        q = [res[n] for n, *_ in jrd]
        _, r_m, t_m = superpose(p, q)
        worst = max(norm(sub(add(mat_vec(r_m, a), t_m), b)) for a, b in zip(p, q))
        print(f"4JRD: MDATK's _JRD_NUCLEOTIDE is chain A residue 10 to {worst:.4f} A (largest deviation after "
              "superposition; its numbers are rounded to 0.001 A / 0.01 deg)")
        if worst > 0.0011:
            diffs.append(("jrd shape", f"_JRD_NUCLEOTIDE is not 4JRD A10 (deviation {worst:.4f} A)"))
        md_z = [r_m[0][2], r_m[1][2], r_m[2][2]]  # MDATK's z axis in 4JRD coordinates
        ours = {n: cylindrical(r, p_, z) for n, r, p_, z in jrd}
        theirs = {n: cylindrical(*md[n]) for n in md}

        def ours_to_pdb(v):
            return add(origin, add(scale(ex, v[0]), add(scale(ey, v[1]), scale(ez, v[2] + z_ref))))
        z_ref = dot(sub(res["C1'"], origin), ez)
        print(f"4JRD: the rigid duplex of the nucleotide and its screw against A3-A11 / B2-B10: RMSD "
              f"{duplex_rmsd(ours, twist, rise, ours_to_pdb):.3f} A in the screw frame, "
              f"{duplex_rmsd(theirs, builder['JRD_TWIST'], builder['JRD_RISE'], lambda v: add(mat_vec(r_m, v), t_m)):.3f} A "
              "in MDATK's")
        tilt = math.degrees(math.acos(max(-1.0, min(1.0, dot(md_z, ez)))))
        if tilt > 0.05:
            diffs.append(("jrd frame", f"MDATK's helix frame: z is {tilt:.2f} deg from the axis of the 4JRD screw "
                                       f"(least squares, A3-A11 / B2-B10)"))
        if any(abs(md[n][0] - r) > 0.0005 or abs(((md[n][1] - p_ + 180) % 360) - 180) > 0.005
               or abs(md[n][2] - z) > 0.0005 for n, r, p_, z in jrd):
            diffs.append(("jrd values", "_JRD_NUCLEOTIDE's cylindrical coordinates differ from the screw frame's "
                                        "(and are rounded)"))
    if builder["JRD_TWIST"] != twist:
        diffs.append(("jrd twist", f"JRD_TWIST {builder['JRD_TWIST']} = the fitted {twist} rounded to "
                                   f"{round(twist, 2)}" if round(twist, 2) == builder["JRD_TWIST"]
                      else f"JRD_TWIST {builder['JRD_TWIST']} != {twist}"))
    if builder["JRD_RISE"] != rise:
        diffs.append(("jrd rise", f"JRD_RISE {builder['JRD_RISE']} = the fitted {rise} rounded to {round(rise, 2)}"
                      if round(rise, 2) == builder["JRD_RISE"] else f"JRD_RISE {builder['JRD_RISE']} != {rise}"))
    # the dictionary components
    dic = mdatk_dictionaries(os.path.join(mdatk, "mdatk", "dictionaries.py"))
    for code, d in components().items():
        c = dic.COMPONENTS[code]
        for field, mine, theirs in (("name", d["name"], c.name), ("base", d["base"], c.base),
                                    ("atoms", tuple(d["atoms"]), c.atoms), ("bonds", tuple(d["bonds"]), c.bonds),
                                    ("hydrogens", tuple(d["hydrogens"]), c.hydrogens)):
            if mine != theirs:
                diffs.append((f"{field} {code}", f"COMPONENTS[{code}].{field}: MDATK {theirs}, CCD {mine}"))
        theirs = {h: (x, a, b, k) for h, x, a, b, k in c.hpos}
        for h, x, a, b, k in d["hpos"]:
            if theirs.get(h) != (x, a, b, k):
                diffs.append((f"hpos {code} {h}", f"COMPONENTS[{code}].hpos {h}: MDATK {theirs.get(h)}, "
                                                  f"CCD {(x, a, b, k)}"))
        if len(theirs) != len(d["hpos"]):
            diffs.append((f"hpos {code}", f"COMPONENTS[{code}].hpos has {len(theirs)} rules, the CCD {len(d['hpos'])}"))
        if tuple((ce, tuple(n), s) for ce, n, s in dic.L_CHIRALITY) != tuple(d["chirality"]):
            diffs.append((f"chirality {code}", f"L_CHIRALITY {dic.L_CHIRALITY} != {d['chirality']}"))
    # the D residues' rules as the dictionaries would mirror them (hpos of the standard residues)
    for base in BASES:
        for h, x, a, b, k in rules[base]:
            if h.startswith("H2") and base == "G":
                lc = next(cc for cc, bb in L_CODES.items() if bb == base)
                theirs = {r[0]: r for r in dic.COMPONENTS[lc].hpos}
                if theirs[h][4] != k:
                    diffs.append((f"hpos G {h}", f"the CCD's G: {h} {k} to {b}; MDATK's 0G rule {theirs[h][4]}"))
    return diffs


def nb_heavy(rows, parent, bonds_of):
    """Heavy neighbours of `parent` in the standard residue `bonds_of`."""
    _, _, bonds = ccd()
    el = {n: r["element"] for n, r in rows.items()}
    out = []
    for b in bonds:
        if b["comp_id"] != bonds_of:
            continue
        if b["atom_id_1"] == parent and el[b["atom_id_2"]] != "H":
            out.append(b["atom_id_2"])
        elif b["atom_id_2"] == parent and el[b["atom_id_1"]] != "H":
            out.append(b["atom_id_1"])
    return out


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__.split("\n\n")[1])
    mdatk = sys.argv[1]
    for path, text in ((HELIX_INC, helix_inc()), (COMPONENT_INC, component_inc())):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    for note in ccd_notes():
        print("note:", note)
    diffs = compare(mdatk)
    failed = False
    seen = set()
    for key, message in diffs:
        seen.add(key)
        entry = DOCUMENTED.get(key)
        print(f"{'documented (' + entry + ')' if entry else 'UNDOCUMENTED'} difference from MDATK 0.0.4.1: {message}")
        failed |= entry is None
    for key in sorted(set(DOCUMENTED) - seen):
        print(f"documented difference no longer found: {key} ({DOCUMENTED[key]})")
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
