"""What MDATK's structure builder (mdatk/builder.py) builds, for MDBridge's chem::helix, helix_fit and
mutate (tests/chem/test_builder.cpp; docs/specs/m6-plan.md 1.3, 4.1, 4.2).

    python3 -I -B make_builder_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy and scipy (least_squares, cKDTree, the density fit). MDBridge's chem improvements are
applied first (improvements_chem.py: CH-01 / CH-02 numbering, CH-08 mutation specs, CH-09 / D11
hydrogens of mutated bases, CH-10 / CH-11 in the fit, D13 "REMARK 951 MDBridge IDEAL ..."). numpy's
and OpenBLAS' kernels are pinned (chem_common.pin_kernels). Writes:

  OUT_DIR/builder_cases.json
    sequences   check_helix(sequence, partner, pairing, triplex, third) messages over sequences,
                partners, pairings and triplexes; parse_sequence, resolve_partner, resolve_third,
                helix_kind, pair_types, pairing_kinds, triplex_kinds, _parse_ids results or errors
    templates   template(base, hydrogens=True): atom -> xyz (float.hex), the A-form monomer
    helices     Watson-Crick build_helix over sequence, partner, chirality, chains, start, scheme,
                terminals, phosphate5, hydrogens: the PDB text and notes, or the error
    mutations   mutate(text, spec, compensatory, wobble) of built helices (D, L, hydrogens, wobble
                pairs, amber names, altlocs, a DNA strand) and parse_mutations specs: text and notes
                or the error
  OUT_DIR/builder_fitted_cases.json (least-squares and density fits; compared within the tolerance
  of docs/specs/m6-plan.md 4.1 / D10)
    reference   per helix kind of _HELIX_START: the start, the refined parameters least_squares
                reaches (the "refined starts": W2-A may start from them), cost, nfev, screw twist /
                rise, per-strand torsions, closest contacts (_helix_fit)
    pair_fits   _pair_fit(rule, b1, b2): rotation and shift of every canonical pair (float.hex)
    rebuilt     build_helix with other pairings and triplexes (MDATK's tests' builds, mismatches,
                L, phosphate5): text, notes and the heavy-atom coordinates as numbers
    fits        fit_to_map(model, map, resolution, search, step): text, cc, notes or the error, on the
                maps of OUT_DIR/maps (simulated from moved models; one with a NaN voxel, CH-11)
  OUT_DIR/maps/*.mrc  the fit cases' maps (mdatk.io.mrc.write_map)

Floats are float.hex() strings unless noted; errors are {"type": class name, "message": str(exc)}.
"""

import json
import os
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_builder_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
sys.path.insert(0, CHECKOUT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_common  # noqa: E402

chem_common.pin_kernels()
import numpy as np  # noqa: E402

import improvements_chem  # noqa: E402

improvements_chem.apply_mdatk()
from mdatk import builder as B  # noqa: E402
from mdatk import density as D  # noqa: E402
from mdatk import pdbtools as pt  # noqa: E402
from mdatk.io.mrc import write_map  # noqa: E402

error, fhex, matrix, vector = chem_common.error, chem_common.fhex, chem_common.matrix, chem_common.vector
ERRORS = (B.BuildError, pt.PDBError, ValueError)


def attempt(fn, *args, **kwargs):
    try:
        return {"value": fn(*args, **kwargs)}
    except ERRORS as exc:
        return {"error": error(exc)}


# --- sequences ---------------------------------------------------------------------------------------

SEQUENCES = ["GGCAUC", "5'-GGCAUC-3'", "ggcauc", "GG CA-UC", "GGTAC", "AAGAGAAGGA", "AUGCAU", "GAGGAUCA", "A",
             "", "  ", "GGXCU", "5'GGC3'", "AAAAAAAA", "CUCUCU", "GGGG", "UUUU", "’5-GGC-3’"]
PARTNERS = ["auto", "none", "", "  ", "GAUGCC", "AAGCU", "GAUGUC", "GGTAC", "CCAC", "UUUU", "UCUCUUCUU", "X", None]
PAIRINGS = ["wc", "reverse-wc", "hoogsteen", "reverse-hoogsteen", "hh", "", None, "WC", "sheared"]
TRIPLEXES = [None, "none", "hoogsteen", "reverse-hoogsteen", "No", "parallel"]


def sequence_cases():
    out = {"check_helix": [], "parse_sequence": [], "resolve_partner": [], "resolve_third": [], "helix_kind": [],
           "pair_types": [], "pairing_kinds": [], "triplex_kinds": [], "parse_ids": []}
    for seq in SEQUENCES:
        for partner in PARTNERS:
            for pairing in PAIRINGS[:6]:
                out["check_helix"].append([seq, partner, pairing, None, "auto",
                                           B.check_helix(seq, partner, pairing)])
        for triplex in TRIPLEXES:
            for third in ("auto", "UUCUCUUCCU", "AGG", "GGC"):
                out["check_helix"].append([seq, "auto", "wc", triplex, third,
                                           B.check_helix(seq, "auto", "wc", triplex, third)])
        out["check_helix"].append([seq, "auto", "reverse-wc", "hoogsteen", "auto",
                                   B.check_helix(seq, "auto", "reverse-wc", "hoogsteen")])
        out["parse_sequence"].append([seq, attempt(B.parse_sequence, seq)])
        for partner in PARTNERS:
            for pairing in ("wc", "hoogsteen", "reverse-hoogsteen", "hh"):
                out["resolve_partner"].append([seq, partner, pairing, attempt(B.resolve_partner, seq, partner, pairing)])
        for third in ("auto", "AUTO", "UUC", "GGCAUC", ""):
            for motif in ("hoogsteen", "reverse-hoogsteen"):
                clean = "".join(c for c in seq.upper() if c in "ACGU")
                out["resolve_third"].append([clean, third, motif, attempt(B.resolve_third, clean, third, motif)])
    for pairing in PAIRINGS:
        for triplex in TRIPLEXES:
            out["helix_kind"].append([pairing, triplex, attempt(B.helix_kind, pairing, triplex)])
    for seq, other in (("GGCAUC", "GAUGCC"), ("GGUAGCC", "GGUUAUC"), ("ACGU", "ACGU"), ("AAGA", "UUUU")):
        out["pair_types"].append([seq, other, B.pair_types(seq, other)])
        for pairing in B.PAIRINGS:
            out["pairing_kinds"].append([seq, other, pairing, [list(k) for k in B.pairing_kinds(seq, other, pairing)]])
        for motif in B.TRIPLEX_MOTIFS:
            out["triplex_kinds"].append([seq, other, motif, [list(k) for k in B.triplex_kinds(seq, other, motif)]])
    for value in ("A,B", "A", "A,B,C", "P, Q", "A,B,C,D", "", "1,1", "5", "-5,1", "1,x", "1,2,3", "7,8"):
        for count in (2, 3):
            kind = int if value[:1].lstrip("-").isdigit() or value == "1,x" else str
            got = attempt(B._parse_ids, value, count, kind)
            out["parse_ids"].append([value, count, kind.__name__, {k: list(v) if k == "value" else v
                                                                   for k, v in got.items()}])
    return out


def template_cases():
    return {base: {n: vector(x) for n, x in B.template(base, hydrogens=True).items()} for base in "ACGU"}


# --- Watson-Crick helices -------------------------------------------------------------------------

def helix(seq, **kw):
    notes = []
    case = {"sequence": seq, "options": kw}
    try:
        case["text"] = B.build_helix(seq, notes=notes, **kw)
        case["notes"] = notes
    except ERRORS as exc:
        case["error"] = error(exc)
    return case


def helix_cases():
    out = [helix("GGCAUGCC"), helix("GGCAUGCC", chirality="L"), helix("GGACU", partner="none", chains="C", start=10),
           helix("GGCAU", partner="AAGCU"), helix("GGGG", partner="CCAC"), helix("GGUAGCC", partner="GGUUAUC"),
           helix("GCA", scheme="amber"), helix("GCA", chirality="L", scheme="amber"),
           helix("GCA", chirality="L", scheme="ccd", terminals=False), helix("GCA", chirality="L", scheme="pdb"),
           helix("GCA", phosphate5=True), helix("GCA", phosphate5=True, scheme="amber"),
           helix("GCAU", hydrogens=True), helix("GCAU", hydrogens=True, chirality="L", scheme="ccd"),
           helix("AAA", partner="none", hydrogens=True), helix("AAA", partner="none", chirality="L", hydrogens=True),
           helix("GGCAUC", chirality="L", scheme="ccd", hydrogens=True), helix("A"), helix("A", partner="none"),
           helix("GGCAU", start="-5,1"), helix("GGCAU", start="9998,1"), helix("GGCAU", start="1,9997"),
           helix("GGCAU", partner="none", start=9996), helix("GGCAU", start="-1002,1"), helix("GGCAU", start="-999,-995"),
           helix("GGCAU", chains="X,Y"), helix("GGCAU", chains="A,A"), helix("GGCAU", chains="AB,C"),
           helix("GGCAU", chains=" ,B"), helix("GGCAU", start="1,x"), helix("GGCAU", start="1,2,3"),
           helix("GGCAU", partner=""), helix("GGACU", partner="AGU"), helix("GGXCU"), helix("  "),
           helix("GGCAU", chirality="X"), helix("GGCAU", chirality="l"), helix("GGCAU", scheme="charmm"),
           helix("GGTAC", partner="GGTAC"),
           helix("GGCAUCGAUCGGAUCCGAUGCAAGCUUGCAU", hydrogens=True, phosphate5=True, scheme="ccd")]
    for scheme in (None, "amber", "ccd", "rosetta", "pdb"):
        for chirality in ("D", "L"):
            for terminals in (True, False):
                for phosphate5 in (False, True):
                    out.append(helix("GAUC", chirality=chirality, scheme=scheme, terminals=terminals,
                                     phosphate5=phosphate5, hydrogens=scheme == "amber", partner="none" if
                                     scheme == "pdb" else "auto"))
    return out


# --- mutations -------------------------------------------------------------------------------------

def mutation_inputs():
    src = B.build_helix("GGCAUCG")
    out = {"helix": src, "helix_h": B.build_helix("GGCAUCG", hydrogens=True),
           "wobble": B.build_helix("GGCAUCG", partner="CGAUGUC"),
           "l_amber": B.build_helix("GAUC", chirality="L", scheme="amber"),
           "l_ccd_h": B.build_helix("GAUC", chirality="L", scheme="ccd", hydrogens=True),
           "blank": B.build_helix("GGCAUCGAUCGGAUCC", partner="none", chains=" ", start=1)}
    alt, dna = [], []
    for ln in src.splitlines():
        if ln.startswith("ATOM") and ln[21] == "A" and ln[22:26].strip() == "3":
            alt.append(ln[:16] + "A" + ln[17:54] + "  0.60" + ln[60:])
            alt.append(ln[:16] + "B" + ln[17:54] + "  0.40" + ln[60:])
        else:
            alt.append(ln)
        if ln.startswith("ATOM") and ln[21] == "B":
            if ln[12:16].strip() == "O2'":
                continue
            name = ln[17:20].strip().rstrip("35")
            ln = ln[:17] + f"{'D' + name:>3}" + ln[20:]
        dna.append(ln)
    out["altlocs"] = "\n".join(alt) + "\n"
    out["dna_partner"] = "\n".join(dna) + "\n"
    out["no_frame"] = "\n".join(ln for ln in src.splitlines() if not (ln[21:26] == "A   2" and ln[12:16].strip() == "C8"))
    out["annotated"] = ("REMARK 950 MDBridge BASE-PAIR ANNOTATION\nREMARK 950 PAIR x\n"
                        + src.replace("END\n", "CONECT    1    2\nEND\n"))
    return out


MUTATIONS = [
    ("helix", "A:2U", {}), ("helix", "A:2U, B:6A", {}), ("helix", "A:2A, B:6A", {}), ("helix", "A:2U", {"compensatory": False}),
    ("helix", "A:99G", {}), ("helix", "A5", {}), ("helix", "A:5G, A:5C", {}), ("helix", "", {}), ("helix", "A:3G", {}),
    ("helix", "A:3A", {}), ("helix", "a:3G", {}), ("helix", "A:3>c; B:1->U", {}), ("helix", "A:-2G", {}),
    ("helix_h", "A:1C, A:4G", {}), ("helix_h", "B:3U", {"compensatory": False}),
    ("wobble", "A:6U", {"wobble": True}), ("wobble", "A:6U", {}), ("wobble", "A:2A", {"wobble": True}),
    ("l_amber", "A:1C", {}), ("l_ccd_h", "A:1C", {}), ("l_ccd_h", "B:4A", {"wobble": True}),
    ("blank", "15G", {}), ("blank", ":15G", {}), ("blank", "2A", {}), ("blank", "1-2A", {}), ("blank", "16T", {}),
    ("altlocs", "A:3G", {"compensatory": False}), ("altlocs", "A:3G", {}),
    ("dna_partner", "A:2C", {}), ("dna_partner", "B:2C", {}), ("no_frame", "A:2C", {}), ("annotated", "A:1C", {}),
]
SPECS = ["A:5g; B:-2>C :7U", "15G", "A15G", ":15G", "A:15G", "-3G", "A:-3G", "1A:5G", "A:5AG", "A:5A>G", "A:5->U",
         "A:5:G", "A:5T", "A:5", "5", "A:5X", "A:5G,B:5G", "A:5G A:5G", " ", "AB:5G", "é:5G", "A:5Ag"]


def mutation_cases():
    inputs = mutation_inputs()
    out = {"inputs": inputs, "mutate": [], "parse_mutations": []}
    for name, spec, kw in MUTATIONS:
        case = {"input": name, "spec": spec, "options": kw}
        try:
            text, notes = B.mutate(inputs[name], spec, **kw)
            case.update(text=text, notes=notes)
        except ERRORS as exc:
            case["error"] = error(exc)
        out["mutate"].append(case)
    for spec in SPECS:
        got = attempt(B.parse_mutations, spec)
        if "value" in got:
            got["value"] = [[m.chain, m.resseq, m.icode, m.base] for m in got["value"]]
        out["parse_mutations"].append([spec, got])
    return out


# --- fitted outputs --------------------------------------------------------------------------------

def reference_cases():
    from scipy.optimize import least_squares

    out = []
    for kind in B._HELIX_START:
        start = np.array(B._HELIX_START[kind], dtype=np.float64)
        sol = least_squares(lambda p, k=kind, s=start: B._fit_residuals(k, p, s), start, x_scale=B._fit_scale(kind))
        fit = B._helix_fit(kind)
        out.append({"kind": kind, "start": list(B._HELIX_START[kind]), "x": vector(sol.x), "cost": fhex(sol.cost),
                    "nfev": int(sol.nfev), "status": int(sol.status), "twist": fhex(fit.screw.twist),
                    "rise": fhex(fit.screw.rise), "torsions": [vector(q) for q in fit.torsions],
                    "contact": fhex(fit.contact), "contact_h": fhex(fit.contact_h)})
    return out


def pair_fit_cases():
    out = []
    for (rule, b1, b2) in sorted(B._PAIR_HBONDS):
        rot, shift = B._pair_fit(rule, b1, b2)
        out.append({"rule": rule, "b1": b1, "b2": b2, "rotation": matrix(rot), "shift": vector(shift)})
    return out


REBUILT = [("reverse-wc", None, "GAGGAUCA", {"hydrogens": True}), ("reverse-wc", None, "AUGCAU", {"hydrogens": True}),
           ("hoogsteen", None, "GAGGAUCA", {"hydrogens": True}), ("reverse-hoogsteen", None, "GAAGAG", {"hydrogens": True}),
           ("reverse-hoogsteen", None, "AUGCAU", {"hydrogens": True}), ("hh", None, "GAGGAUCA", {"hydrogens": True}),
           ("wc", "hoogsteen", "GAGGAUCA", {"hydrogens": True}), ("wc", "reverse-hoogsteen", "AUGCAU", {"hydrogens": True}),
           ("reverse-wc", None, "GAUCAGUAGC", {}), ("hoogsteen", None, "AAGAGAAGGA", {}),
           ("reverse-hoogsteen", None, "AAGAAGGAAG", {}), ("hh", None, "AAAAAAAA", {}),
           ("hh", None, "AAAAAAAAAA", {"phosphate5": True}), ("wc", "hoogsteen", "AAGAGAAGGA", {}),
           ("wc", "reverse-hoogsteen", "AAGAGAAGGA", {}), ("hoogsteen", None, "AAGA", {"partner": "UUUU"}),
           ("wc", "hoogsteen", "AAUAG", {}), ("hoogsteen", None, "AAGAG", {"chirality": "L", "scheme": "amber"}),
           ("wc", "hoogsteen", "AAGAG", {"chains": "P,Q", "start": "5,10"}), ("reverse-wc", None, "A", {}),
           ("hh", None, "AA", {"partner": "none"}), ("wc", "hoogsteen", "AAGAG", {"third": "UUC"}),
           ("wc", "hoogsteen", "AAGAG", {"partner": "none"}), ("reverse-wc", "hoogsteen", "AAGAG", {})]


def heavy_coordinates(text):
    """[[chain, resseq, atom, x, y, z], ...] of the heavy atoms (numbers, for the tolerance comparison)."""
    out = []
    for r in pt.read_pdb(text).residues():
        for a in r.atoms:
            if not a.name.lstrip("0123456789").startswith("H"):
                out.append([r.chain, r.resseq, a.name, *a.xyz])
    return out


def rebuilt_cases():
    out = []
    for pairing, motif, seq, kw in REBUILT:
        notes = []
        case = {"pairing": pairing, "triplex": motif, "sequence": seq, "options": kw}
        try:
            case["text"] = B.build_helix(seq, pairing=pairing, triplex=motif, notes=notes, **kw)
            case["notes"] = notes
            case["heavy"] = heavy_coordinates(case["text"])
        except ERRORS as exc:
            case["error"] = error(exc)
        out.append(case)
    return out


def write_moved_map(name, text, rotvec, shift, resolution, nan_voxel=False):
    xyz, el = B.model_atoms(pt.read_pdb(text))
    rot = D.rotvec_matrix(np.array(rotvec, dtype=np.float64)).tolist()
    # x' = R x + t in plain Python (no BLAS product: the map's bytes must not depend on the CPU)
    moved = np.array([[sum(rot[i][k] * p[k] for k in range(3)) + shift[i] for i in range(3)] for p in xyz.tolist()])
    dmap = D.simulate_map(moved, D.atom_weights(el), resolution)
    if nan_voxel:
        data = dmap.data.copy()
        data[tuple(s // 2 for s in data.shape)] = np.nan
        dmap = type(dmap)(data, dmap.origin, dmap.voxel)
    os.makedirs(os.path.join(OUT, "maps"), exist_ok=True)
    write_map(os.path.join(OUT, "maps", name), dmap)
    return f"maps/{name}"


def fit_cases():
    out = []
    model = B.build_helix("GAAGAGAAGGC")
    single = B.build_helix("GGCAU", partner="none")
    maps = {"global": write_moved_map("helix_moved.mrc", model, [0.4, -0.9, 0.3], [30.0, -12.0, 55.0], 4.0),
            "local": write_moved_map("strand_moved.mrc", single, [0.05, -0.04, 0.03], [1.0, -0.6, 0.4], 4.0),
            "nan": write_moved_map("strand_nan.mrc", single, [0.05, -0.04, 0.03], [1.0, -0.6, 0.4], 4.0, True)}
    near = single
    def chimerax(ln):  # residue name fused with the chain id and a colon (tests/test_pdbtools.py)
        fused = f"{ln[17:20].strip()}{ln[21]}:"
        return ln[:22 - len(fused)] + fused + ln[22:]

    damaged = "\n".join(chimerax(ln) if ln.startswith("ATOM") else ln for ln in single.splitlines()) + "\n"
    cryst = "CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1\n" + single
    far = B.transform_pdb_text(single, np.eye(3), np.array([200.0, 0.0, 0.0]))
    runs = [("model", "global", 4.0, "global", 30.0), ("single", "local", 4.0, "local", 30.0),
            ("single", "nan", 4.0, "local", 30.0), ("damaged", "local", 4.0, "local", 30.0),
            ("cryst", "local", 4.0, "local", 30.0), ("far", "local", 4.0, "local", 30.0),
            ("single", "local", float("nan"), "local", 30.0), ("single", "local", 0.0, "local", 30.0),
            ("single", "local", -2.0, "local", 30.0), ("single", "local", 4.0, "exhaustive", 30.0),
            ("single", "local", 4.0, "global", 2.0), ("single", "missing", 4.0, "local", 30.0),
            ("single", "nan", 4.0, "global", 45.0)]
    texts = {"model": model, "single": near, "damaged": damaged, "cryst": cryst, "far": far}
    here = os.getcwd()
    os.chdir(OUT)  # map paths relative to OUT_DIR (messages name them so)
    for which, map_key, resolution, search, step in runs:
        path = maps.get(map_key, "maps/missing.mrc")
        case = {"model": which, "map": path, "resolution": fhex(resolution), "search": search, "step": fhex(step)}
        try:
            text, cc, notes = B.fit_to_map(texts[which], path, resolution, search, step)
            case.update(text=text, cc=cc, notes=notes)
        except (OSError, *ERRORS) as exc:
            case["error"] = error(exc)
        out.append(case)
    os.chdir(here)
    return {"models": texts, "cases": out}


def main():
    os.makedirs(OUT, exist_ok=True)
    data = {"sequences": sequence_cases(), "templates": template_cases(), "helices": helix_cases(),
            "mutations": mutation_cases()}
    with open(os.path.join(OUT, "builder_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")
    fitted = {"reference": reference_cases(), "pair_fits": pair_fit_cases(), "rebuilt": rebuilt_cases(),
              "fits": fit_cases()}
    with open(os.path.join(OUT, "builder_fitted_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(fitted, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
