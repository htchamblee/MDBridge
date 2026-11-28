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
* src/chem/triplex_tables.inc: the experimental geometry of the Hoogsteen duplex and the
  pyrimidine-motif triplex (docs/DIFFERENCES.md CH-18): mean U·A-U and C+·G-C triples, the screws,
  their spreads, from the RNA triple helices of PDB 3P22, 4PLX, 6SVS, 7JNH (data/sources/pdb/
  triplex_levels.tsv, *_triplex_atoms.tsv); triplex_tables() (see its section) also gives them to
  tests/parity/improvements_chem.py.
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
TRIPLEX_INC = os.path.join(ROOT, "src", "chem", "triplex_tables.inc")

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


# --- RNA triple helices (CH-18) -------------------------------------------------------------------
#
# The experimental reference geometry of the builder's Hoogsteen duplex ('hoogsteen': the purine
# strand and the Hoogsteen strand) and pyrimidine-motif triplex ('th': purine, Watson-Crick and
# Hoogsteen strands): the stacked U·A-U and C+·G-C base triples of the RNA crystal structures in
# TRIPLEX_ENTRIES (data/sources/pdb/triplex_levels.tsv, the ATOM records in <id>_triplex_atoms.tsv).
#
# * Levels: each triple's three nucleosides (C1', the base heavy atoms, the sugar ring with O2',
#   O3', C5'). The mean level of each type (U·A-U, C+·G-C) is their generalised Procrustes mean
#   (unweighted superpositions onto the mean, repeated until it no longer changes); the C+·G-C mean
#   is then superposed onto the U·A-U mean by the glycosidic frames (C1', N9 C4 C8 / N1 C2 C6) of
#   the three strands, which puts both in one level frame.
# * Steps (C1' and base atoms, as the report measures them): two triples whose purine, Watson-Crick and Hoogsteen residues are linked O3'-P (< 2 A)
#   to each other's 3' / 5' / 3' neighbour. The step of the level frame is the mean (quaternions
#   summed, translations averaged) of mean(k)^-1 o mean(k + 1), each mean superposed onto its
#   experimental level; its screw defines the helix frame: z along the axis (one level on: a turn by
#   -twist about z and a move of -rise along z, MDATK's _Screw), x through the C1' of the U·A-U
#   purine, z 0 there. The duplex screw ('hoogsteen') is the same mean over the same steps with the
#   purine and Hoogsteen strands only, written in the same frame.
# * Spreads (the restraints of the reference fit): per step the screw's axis offset (two components
#   in the plane z 0), axis tilts (two angles), twist and rise about the mean screw: their standard
#   deviations; the separation of the base planes of consecutive levels (all base heavy atoms of the
#   kind's strands) and its standard deviation; the donor-acceptor distances of the U·A-U triples'
#   Watson-Crick and Hoogsteen pairs; the RMS rotation (degrees) and C1' shift (A) per component of
#   a U·A-U partner nucleoside about its mean pose relative to the purine nucleoside.
TRIPLEX_ENTRIES = ("3P22", "4PLX", "6SVS", "7JNH")
TRIPLEX_BASE_ATOMS = {"A": ("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"),
                      "G": ("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"),
                      "C": ("N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"),
                      "U": ("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6")}
TRIPLEX_SUGAR = ("C2'", "O2'", "C3'", "O3'", "C4'", "O4'", "C5'")
TRIPLEX_GLYCOSIDIC = {"A": ("C1'", "N9", "C4", "C8"), "G": ("C1'", "N9", "C4", "C8"),
                      "C": ("C1'", "N1", "C2", "C6"), "U": ("C1'", "N1", "C2", "C6")}
TRIPLEX_HBONDS = (("wc", "N1", "N3"), ("wc", "N6", "O4"), ("hoogsteen", "N7", "N3"), ("hoogsteen", "N6", "O4"))
TRIPLEX_KINDS = (("triplex-hoogsteen", (0, 1, 2)), ("hoogsteen", (0, 2)))
TRIPLEX_STRANDS = ("purine", "wc", "hoogsteen")
# the H-bonds of a triple that tools/select_pdb_triplexes.py selects by (each under 3.7 A): purine atom,
# partner atom, which pair
TRIPLEX_CRITERIA = {
    "AUU": (("N1", "N3", "Watson-Crick"), ("N6", "O4", "Watson-Crick"), ("N7", "N3", "Hoogsteen"),
            ("N6", "O4", "Hoogsteen")),
    "GCC": (("N1", "N3", "Watson-Crick"), ("O6", "N4", "Watson-Crick"), ("N7", "N3", "Hoogsteen"),
            ("O6", "N4", "Hoogsteen")),
}


def nucleoside_names(base):
    return ("C1'",) + TRIPLEX_BASE_ATOMS[base] + TRIPLEX_SUGAR


def triplex_residues(entry):
    out = {}
    for r in read_tsv("pdb", f"{entry.lower()}_triplex_atoms.tsv"):
        if r["label_alt_id"] not in (".", "A") or r["pdbx_PDB_model_num"] != "1":
            continue
        res = out.setdefault((r["auth_asym_id"], int(r["auth_seq_id"])), {"comp": r["label_comp_id"], "atoms": {}})
        res["atoms"].setdefault(r["label_atom_id"],
                                [float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])])
    return out


def _linked(a, b):
    """O3' of a bonded to P of b."""
    return "O3'" in a["atoms"] and "P" in b["atoms"] and norm(sub(a["atoms"]["O3'"], b["atoms"]["P"])) < 2.0


def triplex_data():
    """(levels, steps): levels [(entry, type 'AUU' / 'GCC', [3 residues])]; steps [(index k, index k+1)]."""
    levels, keys = [], []
    for entry in TRIPLEX_ENTRIES:
        res = triplex_residues(entry)
        for row in read_tsv("pdb", "triplex_levels.tsv"):
            if row["entry"] != entry:
                continue
            trio = [res[(row[f"{s}_chain"], int(row[f"{s}_seq"]))] for s in TRIPLEX_STRANDS]
            kind = "".join(r["comp"] for r in trio)
            if kind not in ("AUU", "GCC"):
                raise SystemExit(f"{entry} level {row['level']}: triple {kind} (expected U·A-U or C+·G-C)")
            missing = [(s, n) for s, r in zip(TRIPLEX_STRANDS, trio) for n in nucleoside_names(r["comp"])
                       if n not in r["atoms"]]
            if missing:
                raise SystemExit(f"{entry} level {row['level']}: atoms {missing} missing")
            p, w, t = (r["atoms"] for r in trio)
            for x, y, d in TRIPLEX_CRITERIA[kind]:
                partner = w if d.startswith("Watson-Crick") else t
                if norm(sub(p[x], partner[y])) >= 3.7:
                    raise SystemExit(f"{entry} level {row['level']}: {d} {x}-{y} {norm(sub(p[x], partner[y])):.2f} A")
            # the copies of one molecule in the crystal (NCS / crystallographic copies) share their
            # Watson-Crick and Hoogsteen residue numbers: one site
            site = (entry, int(row["wc_seq"]), int(row["hoogsteen_seq"]))
            levels.append((entry, kind, trio, site))
            keys.append((entry, row))
    steps = []
    for i, (ea, _ka, a, _sa) in enumerate(levels):
        for j, (eb, _kb, b, _sb) in enumerate(levels):
            if ea == eb and _linked(a[0], b[0]) and _linked(b[1], a[1]) and _linked(a[2], b[2]):
                steps.append((i, j))
    return levels, steps


def level_points(trio, strands=(0, 1, 2)):
    return [r["atoms"][n] for s in strands for r in (trio[s],) for n in nucleoside_names(r["comp"])]


def moved(r, t, points):
    return [add(mat_vec(r, x), t) for x in points]


def mat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def transpose(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def procrustes_mean(sets, weights):
    """Weighted generalised Procrustes mean of point sets (each superposed onto the mean, until it is
    stable)."""
    mean = [list(x) for x in sets[0]]
    total = sum(weights)
    for _ in range(200):
        fitted = [moved(*superpose(points, mean)[1:], points) for points in sets]
        new = [scale([sum(w * f[i][c] for f, w in zip(fitted, weights)) for c in range(3)], 1.0 / total)
               for i in range(len(mean))]
        new = moved(*superpose(new, mean)[1:], new)
        change = max(norm(sub(a, b)) for a, b in zip(new, mean))
        mean = new
        if change < 1e-12:
            break
    return mean


def quat_matrix(q):
    w, x, y, z = q
    n = math.sqrt(w * w + x * x + y * y + z * z)
    w, x, y, z = w / n, x / n, y / n, z / n
    return [[w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z]]


def matrix_quat(r):
    w = math.sqrt(max(0.0, 1.0 + r[0][0] + r[1][1] + r[2][2])) / 2.0
    return [w, (r[2][1] - r[1][2]) / (4 * w), (r[0][2] - r[2][0]) / (4 * w), (r[1][0] - r[0][1]) / (4 * w)]


def screw_of(r, t):
    """(angle degrees, unit axis (right-hand rule), translation along it, point on the axis) of x -> R x + t."""
    w, x, y, z = matrix_quat(r)
    v = [x, y, z]
    angle = 2.0 * math.atan2(norm(v), w)
    axis = unit(v)
    d = dot(axis, t)
    t_perp = sub(t, scale(axis, d))
    point = scale(add(t_perp, scale(cross(axis, t_perp), 1.0 / math.tan(angle / 2.0))), 0.5)
    return math.degrees(angle), axis, d, point


def symmetric_min_eigenvector(m):
    """Unit eigenvector of the smallest eigenvalue of a symmetric 3 x 3 matrix (cyclic Jacobi)."""
    a = [row[:] for row in m]
    v = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(100):
        if sum(a[i][j] ** 2 for i in range(3) for j in range(3) if i != j) < 1e-30:
            break
        for p_ in range(2):
            for q in range(p_ + 1, 3):
                if abs(a[p_][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p_][p_]) / (2.0 * a[p_][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(3):
                    akp, akq = a[k][p_], a[k][q]
                    a[k][p_], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(3):
                    apk, aqk = a[p_][k], a[q][k]
                    a[p_][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(3):
                    vkp, vkq = v[k][p_], v[k][q]
                    v[k][p_], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    best = min(range(3), key=lambda i: a[i][i])
    return unit([v[k][best] for k in range(3)])


def plane(points):
    c = scale([sum(x[i] for x in points) for i in range(3)], 1.0 / len(points))
    m = [[sum((x[i] - c[i]) * (x[j] - c[j]) for x in points) for j in range(3)] for i in range(3)]
    return c, symmetric_min_eigenvector(m)


def plane_separation(a, b):
    ca, na = plane(a)
    cb, nb = plane(b)
    if dot(na, nb) < 0:
        nb = scale(nb, -1.0)
    n = unit(add(na, nb))
    return abs(dot(sub(cb, ca), n))


def mean_sd(values):
    m = sum(values) / len(values)
    return m, math.sqrt(sum((v - m) ** 2 for v in values) / (len(values) - 1))


def site_means(values, sites):
    """The mean of the values of each site (the copies of one molecule counted once), in the order of
    the sites' first values."""
    groups = {}
    for v, s in zip(values, sites):
        groups.setdefault(s, []).append(v)
    return [sum(g) / len(g) for g in groups.values()]


def site_weights(sites):
    """1 / (number of values of its site) per value: each site weighs one."""
    count = {}
    for s in sites:
        count[s] = count.get(s, 0) + 1
    return [1.0 / count[s] for s in sites]


def perp_basis(a):
    u1 = unit(cross(a, [1.0, 0.0, 0.0]))
    return u1, cross(a, u1)


def step_parameters(r, t, origin, axis):
    """(offset 1, offset 2, tilt 1, tilt 2, twist, rise) of the step x -> R x + t about the screw
    (origin, axis): the step's axis a = cos t1 cos t2 axis + sin t2 u1 - sin t1 cos t2 u2 meets the
    plane through origin normal to axis at origin + o1 u1 + o2 u2 (u1 = axis x (1, 0, 0) normalised,
    u2 = axis x u1); its twist and rise as _Screw has them (a turn of -twist about a)."""
    angle, u, d, point = screw_of(r, t)
    a = scale(u, -1.0)
    u1, u2 = perp_basis(axis)
    s = dot(sub(origin, point), axis) / dot(a, axis)
    q = add(point, scale(a, s))
    t2 = math.degrees(math.asin(max(-1.0, min(1.0, dot(a, u1)))))
    t1 = -math.degrees(math.atan2(dot(a, u2), dot(a, axis)))
    return [dot(sub(q, origin), u1), dot(sub(q, origin), u2), t1, t2, angle, d]


def subset(points, trio, strands):
    out, k = [], 0
    for s in range(3):
        n = len(nucleoside_names(trio[s]["comp"]))
        if s in strands:
            out += points[k:k + n]
        k += n
    return out


def matched_points(a, b, strands):
    """The C1' and base atoms that two levels' residues of each strand both have (a U·A-U and a C+·G-C level:
    C1', N9 ... C4 of the purine, N1 ... C6 but O4 / N4 of the pyrimidines), as two point lists."""
    p, q = [], []
    for s in strands:
        ra, rb = a[s], b[s]
        for n in ("C1'",) + TRIPLEX_BASE_ATOMS[ra["comp"]]:
            if n in TRIPLEX_BASE_ATOMS[rb["comp"]] or n == "C1'":
                p.append(ra["atoms"][n])
                q.append(rb["atoms"][n])
    return p, q


def bases_only(points, trio, strands):
    """The C1' and base atoms of the points of a level's strands (level_points order)."""
    flags = [n not in TRIPLEX_SUGAR for s in strands for n in nucleoside_names(trio[s]["comp"])]
    return [x for x, keep in zip(points, flags) if keep]


def triplex_tables():
    """The CH-18 tables (see the section comment): {"levels": {type: [[(name, x, y, z)] per strand]},
    "screws": {kind: {"origin", "axis", "twist", "rise", "sd" (6), "stack", "stack_sd", "steps", "sites",
    "incl", "incl_sd", "twist_range", "rise_range", "incl_range" ([min, max]), "stack_max" (per strand)}},
    "hbonds": [(rule, atom1, atom2, mean, sd)], "pose_sd": [(strand, rotation deg, shift A)],
    "counts": {type: levels}, "sites": {type: independent levels}}, rounded to 1e-6. The copies of one
    molecule in a crystal count once: the means are weighted by 1 / copies, the standard deviations are
    those of the sites' means (the sites: triplex_data)."""
    levels, steps = triplex_data()
    by_type = {}
    for i, (_e, kind, _trio, _site) in enumerate(levels):
        by_type.setdefault(kind, []).append(i)
    means = {kind: procrustes_mean([level_points(levels[i][2]) for i in idx],
                                   site_weights([levels[i][3] for i in idx])) for kind, idx in by_type.items()}
    step_sites = [levels[i][3] for i, _j in steps]  # a step's site: its first level's
    step_weights = site_weights(step_sites)
    first = {kind: levels[idx[0]][2] for kind, idx in by_type.items()}

    def glyco(kind):
        out = []
        trio = first[kind]
        for s in range(3):
            names = nucleoside_names(trio[s]["comp"])
            off = sum(len(nucleoside_names(trio[x]["comp"])) for x in range(s))
            out += [means[kind][off + names.index(n)] for n in TRIPLEX_GLYCOSIDIC[trio[s]["comp"]]]
        return out
    for kind in means:
        if kind != "AUU":
            _q, r, t = superpose(glyco(kind), glyco("AUU"))
            means[kind] = moved(r, t, means[kind])

    def step_moves(strands):
        out = []
        for i, j in steps:
            (_e, ka, a, _sa), (_f, kb, b, _sb) = levels[i], levels[j]
            _q, ra, ta = superpose(bases_only(subset(means[ka], first[ka], strands), first[ka], strands),
                                   bases_only(level_points(a, strands), a, strands))
            _q, rb, tb = superpose(bases_only(subset(means[kb], first[kb], strands), first[kb], strands),
                                   bases_only(level_points(b, strands), b, strands))
            rat = transpose(ra)
            out.append((mat_mul(rat, rb), mat_vec(rat, sub(tb, ta))))
        return out

    def mean_move(moves, centre):
        """The mean move about centre, each site weighing one: quaternions summed (signs aligned with
        the first), the moved centres averaged (a mean of translations would depend on the origin)."""
        q0 = matrix_quat(moves[0][0])
        acc = [0.0] * 4
        for (r, _t), w in zip(moves, step_weights):
            q = matrix_quat(r)
            sign = 1.0 if sum(a * b for a, b in zip(q, q0)) >= 0 else -1.0
            acc = [a + w * sign * b for a, b in zip(acc, q)]
        rm = quat_matrix(acc)
        c = scale([sum(w * add(mat_vec(r, centre), t)[k] for (r, t), w in zip(moves, step_weights))
                   for k in range(3)], 1.0 / sum(step_weights))
        return rm, sub(c, mat_vec(rm, centre))

    def centroid(points):
        return scale([sum(x[k] for x in points) for k in range(3)], 1.0 / len(points))

    # the helix frame from the triplex screw
    tri_moves = step_moves((0, 1, 2))
    tri_mean = mean_move(tri_moves, centroid(means["AUU"]))
    angle, u, d, point = screw_of(*tri_mean)
    ez = scale(u, -1.0)
    c1 = sub(means["AUU"][0], point)
    ex = unit(sub(c1, scale(ez, dot(c1, ez))))
    ey = cross(ez, ex)
    z0 = dot(c1, ez)
    frame = [ex, ey, ez]

    def to_frame(x):
        v = mat_vec(frame, sub(x, point))
        return [v[0], v[1], v[2] - z0]
    frame_means = {kind: [to_frame(x) for x in pts] for kind, pts in means.items()}

    def frame_move(r, t):
        """x -> R x + t in frame coordinates y = F (x - p) - z0 e: y -> F R F^T y + F (R p + t - p) +
        (F R F^T - 1) z0 e."""
        rf = mat_mul(mat_mul(frame, r), transpose(frame))
        e = [0.0, 0.0, z0]
        shift = add(mat_vec(frame, sub(add(mat_vec(r, point), t), point)), sub(mat_vec(rf, e), e))
        return rf, shift

    out = {"levels": {}, "screws": {}, "hbonds": [], "pose_sd": [], "counts": {}, "sites": {}}
    for kind, pts in frame_means.items():
        trio = first[kind]
        rows, k = [], 0
        for s in range(3):
            names = nucleoside_names(trio[s]["comp"])
            rows.append([(n, *[rounded(v) for v in pts[k + i]]) for i, n in enumerate(names)])
            k += len(names)
        out["levels"][kind] = rows
        out["counts"][kind] = len(by_type[kind])
        out["sites"][kind] = len({levels[i][3] for i in by_type[kind]})
    for name, strands in TRIPLEX_KINDS:
        raw = step_moves(strands)
        moves = [frame_move(r, t) for r, t in raw]
        angle, u, d, point_f = screw_of(*frame_move(*mean_move(raw, centroid(subset(means["AUU"], first["AUU"],
                                                                                    strands)))))
        axis = scale(u, -1.0)
        s = -point_f[2] / axis[2]
        origin = add(point_f, scale(axis, s))
        params = [step_parameters(r, t, origin, axis) for r, t in moves]
        sds = [mean_sd(site_means([p[c] for p in params], step_sites))[1] for c in range(6)]
        seps, incl, twists, rises, stacking = [], [], [], [], [[] for _ in strands]
        for (i, j), (r, t) in zip(steps, raw):
            pa = [x for st in strands for x in (levels[i][2][st]["atoms"][n] for n in TRIPLEX_BASE_ATOMS[levels[i][2][st]["comp"]])]
            pb = [x for st in strands for x in (levels[j][2][st]["atoms"][n] for n in TRIPLEX_BASE_ATOMS[levels[j][2][st]["comp"]])]
            seps.append(plane_separation(pa, pb))
            # inclination: the level's base plane against the step's axis (the steps of the experimental
            # levels themselves: the move of level i onto level j by their C1' and base atoms)
            _q, rr, tt = superpose(*matched_points(levels[i][2], levels[j][2], strands))
            twist_step, u_step, rise_step, _p = screw_of(rr, tt)
            incl.append(math.degrees(math.acos(min(1.0, abs(dot(plane(pa)[1], u_step))))))
            twists.append(twist_step)
            rises.append(abs(rise_step))
            for k, st in enumerate(strands):  # each strand's stacking: its bases' centres from the other's plane
                ca, na = plane([levels[i][2][st]["atoms"][n] for n in TRIPLEX_BASE_ATOMS[levels[i][2][st]["comp"]]])
                cb, nb = plane([levels[j][2][st]["atoms"][n] for n in TRIPLEX_BASE_ATOMS[levels[j][2][st]["comp"]]])
                v = sub(cb, ca)
                stacking[k].append(0.5 * (abs(dot(v, na)) + abs(dot(v, nb))))
        stack, stack_sd = mean_sd(site_means(seps, step_sites))
        incl_mean, incl_sd = mean_sd(site_means(incl, step_sites))
        out["screws"][name] = {"origin": [rounded(v) for v in origin], "axis": [rounded(v) for v in axis],
                               "twist": rounded(angle), "rise": rounded(d), "sd": [rounded(v) for v in sds],
                               "stack": rounded(stack), "stack_sd": rounded(stack_sd), "steps": len(moves),
                               "sites": len(set(step_sites)), "incl": rounded(incl_mean), "incl_sd": rounded(incl_sd),
                               # D19: the extremes of the experimental steps (every step, not the sites' means)
                               "twist_range": [rounded(min(twists)), rounded(max(twists))],
                               "rise_range": [rounded(min(rises)), rounded(max(rises))],
                               "incl_range": [rounded(min(incl)), rounded(max(incl))],
                               "stack_max": [rounded(max(v)) for v in stacking]}
    auu = [levels[i][2] for i in by_type["AUU"]]
    auu_sites = [levels[i][3] for i in by_type["AUU"]]
    for rule, x, y in TRIPLEX_HBONDS:
        s = 1 if rule == "wc" else 2
        m, sd = mean_sd(site_means([norm(sub(t[0]["atoms"][x], t[s]["atoms"][y])) for t in auu], auu_sites))
        out["hbonds"].append((rule, x, y, rounded(m), rounded(sd)))
    mean_level = means["AUU"]
    n0 = len(nucleoside_names("A"))
    nu = len(nucleoside_names("U"))
    for s in (1, 2):
        off = n0 + (s - 1) * nu
        ref0, refs = mean_level[:n0], mean_level[off:off + nu]
        c1 = refs[0]
        rot2, shift2 = [], []
        for t, w in zip(auu, site_weights(auu_sites)):
            pts = level_points(t)
            _q, r, tt = superpose(pts[:n0], ref0)  # the experimental level onto the mean by its purine
            partner = moved(r, tt, pts[off:off + nu])
            (w, qx, qy, qz), rr, t2 = superpose(refs, partner)
            ang = 2.0 * math.atan2(math.sqrt(qx * qx + qy * qy + qz * qz), w)
            vec = scale(unit([qx, qy, qz]), math.degrees(ang)) if ang > 1e-12 else [0.0, 0.0, 0.0]
            rot2 += [w * v * v for v in vec]
            shift2 += [w * v * v for v in sub(add(mat_vec(rr, c1), t2), c1)]
        total = 3 * len(set(auu_sites))  # three components per site
        out["pose_sd"].append((s, rounded(math.sqrt(sum(rot2) / total)), rounded(math.sqrt(sum(shift2) / total))))
    return out


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


def triplex_inc():
    t = triplex_tables()
    lines = [
        "// Generated by tools/gen_chem_tables.py from data/sources/pdb/triplex_levels.tsv and",
        "// data/sources/pdb/{3p22,4plx,6svs,7jnh}_triplex_atoms.tsv (PDB 3P22, 4PLX, 6SVS, 7JNH, CC0): the",
        "// experimental geometry of the builder's Hoogsteen duplex and pyrimidine-motif triplex (docs/",
        "// DIFFERENCES.md CH-18; tools/gen_chem_tables.py triplex_tables). Do not edit.",
        "// Included inside a namespace that defines TriplexAtom, TriplexScrew, TriplexHbond and",
        "// TriplexPoseSpread (mdbridge/chem/reference_tables.hpp).",
        "",
        f"// mean levels ({t['counts']['AUU']} U·A-U and {t['counts']['GCC']} C+·G-C triples) in the helix frame of the "
        "triplex screw:",
        "// strand 0 purine, 1 Watson-Crick, 2 Hoogsteen; C1', base heavy atoms, C2' O2' C3' O3' C4' O4' C5'",
        "inline constexpr TriplexAtom kTriplexLevels[] = {"]
    for kind in ("AUU", "GCC"):
        for strand, rows in enumerate(t["levels"][kind]):
            lines += [f"    {{{cpp_str(kind)}, {strand}, {cpp_str(n)}, {cpp_num(x)}, {cpp_num(y)}, {cpp_num(z)}}},"
                      for n, x, y, z in rows]
    lines += ["};", "",
              "// per kind: the screw (origin, axis; a turn of -twist about the axis and a move of -rise along it",
              "// per level), the standard deviations of the steps' axis offsets (2, Angstrom), axis tilts (2,",
              "// degrees), twist and rise about it, the base-plane separation of consecutive levels and its",
              "// standard deviation (over the independent steps: the copies of one molecule's step averaged), the",
              "// number of steps and of independent steps, the inclination of the levels' base planes to the",
              "// steps' axes (degrees) and its standard deviation; the experimental steps' extremes of twist,",
              "// rise and inclination, and per strand the largest stacking of consecutive bases",
              "inline constexpr TriplexScrew kTriplexScrews[] = {"]
    for kind, sc in t["screws"].items():
        nums = ", ".join(cpp_num(v) for v in sc["origin"] + sc["axis"] + [sc["twist"], sc["rise"]])
        sds = ", ".join(cpp_num(v) for v in sc["sd"])
        lines.append(f"    {{{cpp_str(kind)}, {nums}, {{{sds}}}, {cpp_num(sc['stack'])}, {cpp_num(sc['stack_sd'])}, "
                     f"{sc['steps']}, {sc['sites']}, {cpp_num(sc['incl'])}, {cpp_num(sc['incl_sd'])}, "
                     + ", ".join("{" + ", ".join(cpp_num(v) for v in sc[key]) + "}"
                                 for key in ("twist_range", "rise_range", "incl_range"))
                     + ", {" + ", ".join(cpp_num(v) for v in (sc["stack_max"] + [0.0] * 3)[:3]) + "}},")
    sites = t["sites"]["AUU"]
    lines += ["};", "", "// donor-acceptor distances (mean, standard deviation over the independent U·A-U triples, Angstrom) of",
              "// the U·A-U triples' pairs, the number of independent triples",
              "inline constexpr TriplexHbond kTriplexHbonds[] = {"]
    lines += [f"    {{{cpp_str(r)}, {cpp_str(a)}, {cpp_str(b)}, {cpp_num(m)}, {cpp_num(sd)}, {sites}}},"
              for r, a, b, m, sd in t["hbonds"]]
    lines += ["};", "", "// RMS per component of a U·A-U partner's rotation (degrees) and C1' shift (Angstrom) about its",
              "// mean pose relative to the purine nucleoside (each independent triple weighing one), by table",
              "// strand; the number of independent triples",
              "inline constexpr TriplexPoseSpread kTriplexPoseSpreads[] = {"]
    lines += [f"    {{{st}, {cpp_num(rot)}, {cpp_num(sh)}, {sites}}}," for st, rot, sh in t["pose_sd"]]
    lines += ["};"]
    return "\n".join(lines) + "\n"


def protonated_cytidine():
    """(code, base, proton, proton's heavy atom) of the CCD's N3-protonated cytidine (data/sources/ccd/
    protonated_cytidine_*.tsv): the builder's C+ (docs/DIFFERENCES.md CH-21). Checked: a released RNA
    linking component whose parent (mon_nstd_parent_comp_id) is C, whose atoms are C's (names and
    elements) and one hydrogen more, bonded to one heavy atom only."""
    comps = {r["comp_id"]: r for r in read_tsv("ccd", "protonated_cytidine_components.tsv")}
    code = "CH"
    comp = comps[code]
    if comp["parent"] != "C" or comp["release_status"] != "REL" or comp["type"] != "RNA LINKING":
        raise SystemExit(f"{code}: not a released RNA-linking component of parent C ({comp})")
    own = {r["atom_id"]: r["element"] for r in read_tsv("ccd", "protonated_cytidine_atoms.tsv") if r["comp_id"] == code}
    _comps, atoms, bonds = ccd()
    std = {r["atom_id"]: r["element"] for r in atoms if r["comp_id"] == "C"}
    extra = [n for n in own if n not in std]
    if any(own[n] != e for n, e in std.items() if n in own) or set(std) - set(own) or len(extra) != 1 \
            or own[extra[0]] != "H":
        raise SystemExit(f"{code}: atoms are not those of C and one hydrogen ({sorted(set(own) ^ set(std))})")
    proton = extra[0]
    partners = [b["atom_id_2"] if b["atom_id_1"] == proton else b["atom_id_1"]
                for b in read_tsv("ccd", "protonated_cytidine_bonds.tsv")
                if b["comp_id"] == code and proton in (b["atom_id_1"], b["atom_id_2"])]
    if len(partners) != 1 or own[partners[0]] != "N":
        raise SystemExit(f"{code}: {proton} bonded to {partners}")
    return code, comp["one_letter_code"], proton, partners[0]


def component_inc():
    comps = components()
    rules = hydrogen_rules()
    lines = [
        "// Generated by tools/gen_chem_tables.py from data/sources/ccd/ribonucleotide_*.tsv and",
        "// protonated_cytidine_*.tsv (wwPDB Chemical Component Dictionary, CC0). Do not edit.",
        "// Included inside a namespace that defines ComponentInfo {code, base, name}, ComponentAtom {code,",
        "// name, element, charge}, ComponentBond {code, atom 1, atom 2, order}, ComponentHydrogen {code,",
        "// heavy atom, hydrogen}, HydrogenPosition {residue, H, X, A, B, kind}, ChiralCentre {centre, a, b,",
        "// c, sign} and ProtonatedBase {code, base, proton, heavy atom}.",
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
    code, base, proton, heavy = protonated_cytidine()
    lines += ["};", "",
              "// the CCD's N3-protonated cytidine (the builder's C+): its code, parent base, the extra hydrogen and",
              "// the atom it is bonded to",
              f"inline constexpr ProtonatedBase kProtonatedCytidine = {{{cpp_str(code)}, '{base}', {cpp_str(proton)}, "
              f"{cpp_str(heavy)}}};"]
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
    for path, text in ((HELIX_INC, helix_inc()), (COMPONENT_INC, component_inc()), (TRIPLEX_INC, triplex_inc())):
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
