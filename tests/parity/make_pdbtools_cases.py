"""What MDATK's PDB tools (mdatk/pdbtools.py) make of structures, for MDBridge's chem::pdbtools
(tests/chem/test_pdbtools.cpp; docs/specs/m6-plan.md 1.1, 4.2).

    python3 -I -B make_pdbtools_cases.py MDATK_CHECKOUT OUT_DIR

Needs pytest (MDATK's tests/test_pdbtools.py builds the fixtures: make_pdb, atom_line, numbered and
the HYBRID / CPPTRAJ / REFERENCE / IONS / GAPPED structures). Standard library otherwise (the PDB
tools are pure Python). MDBridge's chem improvements are applied first (improvements_chem.py:
CH-01, CH-03 ... CH-07, CH-15 and NA-01's L-DNA names). The output is reproducible (fixed seeds).
Writes OUT_DIR/pdbtools_cases.json:

  inputs        PDB texts the cases name by index
  scheme_name   scheme_name(base, chirality, terminal, scheme) of every base x hand x terminal x
                scheme (and an unknown scheme): the name or the error
  chiral        chiral_volume of point quadruples and sugar_chirality(coords, base) of MDATK's test
                nucleotides (D, mirrored, distorted, flat, incomplete; volumes as float.hex)
  ranges        parse_ranges of texts and tuple lists: [[chain, start, end], ...] or the error
  operations    one function alone on a read input (read_pdb, then assign_chirality where the
                function needs the hands): assign_chirality, convert_names, to_lrna, set_records,
                restore_chains, split_chains, combine_chains, order_chains, strip_5prime_phosphates,
                renumber_residues, seqres_records, fasta, summary: notes (or the error), the residues
                after it ([chain, resseq, icode, name, chirality]) and write_pdb of the structure
  clean         clean(text, **options) over MDATK's test calls and a grid of option combinations
                (fixed seed): text, notes, fasta or the error

Errors are {"type": class name, "message": str(exc)}.
"""

import importlib.util
import json
import os
import random
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_pdbtools_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
sys.path.insert(0, CHECKOUT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import improvements_chem  # noqa: E402

improvements_chem.apply_mdatk()
from mdatk import pdbtools as pt  # noqa: E402

spec = importlib.util.spec_from_file_location("mdatk_test_pdbtools", os.path.join(CHECKOUT, "tests", "test_pdbtools.py"))
tp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tp)


def error(exc):
    return {"type": type(exc).__name__, "message": str(exc)}


def residues(s):
    """[chain, resseq, icode, name, chirality, record of the first atom] per residue, by TER segment."""
    return [[[r.chain, r.resseq, r.icode, r.name, r.chirality, r.atoms[0].record if r.atoms else ""] for r in seg]
            for seg in s.segments]


# --- inputs ---------------------------------------------------------------------------------------

def line(serial, name, resname, chain, resseq, xyz, record="ATOM"):
    """tests/test_pdbtools.py's atom_line, with a 4-character residue name in columns 18-21."""
    if len(resname) <= 3:
        return tp.atom_line(serial, name, resname, chain, resseq, xyz, record)
    ln = tp.atom_line(serial, name, "XXX", chain, resseq, xyz, record)
    return ln[:17] + resname[:4] + ln[21:]


def long_names(text, names):
    """make_pdb text with the placeholder residue names XXX, XXY ... replaced by 4-character ones."""
    out = []
    for ln in text.splitlines():
        key = ln[17:21]
        if ln.startswith(("ATOM", "HETATM", "TER")) and key.strip() in names:
            ln = ln[:17] + names[key.strip()] + ln[21:]
        out.append(ln)
    return "\n".join(out) + "\n"


def linked(names, chain="A", first=1, records=None, hydrogens=True, deoxy=(), mirror=()):
    """Residues bonded O3'-P (1.6 Å), as a real strand: (name, base) per residue."""
    lines, serial = [], 0
    shift = [0.0, 0.0, 0.0]
    for k, (name, base) in enumerate(names):
        atoms = tp.residue_atoms(base, k in mirror)
        if not hydrogens:
            atoms = [a for a in atoms if not a[0].startswith("H")]
        if k in deoxy:
            atoms = [a for a in atoms if a[0] not in ("O2'", "HO2'")]
        pos = {n: (x, y, z) for n, x, y, z in atoms}
        if k:
            shift = [prev_o3[i] + (1.6 if i == 0 else 0.0) - pos["P"][i] for i in range(3)]
        for n, x, y, z in atoms:
            serial += 1
            xyz = (x + shift[0], y + shift[1], z + shift[2])
            lines.append(line(serial, n, name, chain, first + k, xyz, (records or {}).get(name, "ATOM")))
        prev_o3 = [pos["O3'"][i] + shift[i] for i in range(3)]
    lines.append("TER")
    return "\n".join(lines) + "\nEND\n"


def inputs():
    mk = tp.make_pdb
    out = [
        mk(tp.HYBRID),                                                                   # 0
        mk(tp.CPPTRAJ, extra=tp.IONS),                                                   # 1
        mk(tp.REFERENCE),                                                                # 2
        tp.numbered(tp.GAPPED),                                                          # 3
        mk([("", [("G5", 1, "G", True), ("C", 2, "C", True), ("C3", 3, "C", True)]),     # 4 mirrored plain names
            ("", [("A5", 4, "A", False), ("U3", 5, "U", False), ("G", 6, "G", False)])]),
        mk([("A", [("G5", 1, "G", False), ("C", 2, "C", False), ("A", 3, "A", False), ("U3", 4, "U", False)]),
            ("B", [("G5", 1, "G", False), ("C3", 2, "C", False)])]),                     # 5 to_lrna
        mk([("A", [("G", 1, "G", False)]), ("A", [("C", 2, "C", False)]),
            ("B", [("U", 3, "U", True)]), ("A", [("A", 4, "A", False)])], extra=tp.IONS),  # 6 repeated chain
        mk([("C", [("G", 1, "G", False)]), ("A", [("C", 2, "C", False)]), ("B", [("A", 3, "A", False)])],
           extra=tp.IONS),                                                               # 7 order
        mk([("A", [("G", 1, "G", False), ("C", 1, "C", False)])], first_serial=10,
           extra=["CONECT   10   11", "CONECT   11   10   12"]),                          # 8 duplicate numbers
        mk([("A", [("G", 1, "G", False)])]).replace("END\n", "") + mk([("B", [("C", 1, "C", False)])]).replace(
            "END\n", "") + "CONECT    1    2\nCONECT    2    1\nEND\n",                     # 9 repeated serials (CH-04)
        linked([("G", "G"), ("C", "C"), ("U", "U")]),                                    # 10 bonded strand, H
        linked([("PSU", "U"), ("G", "G"), ("C", "C"), ("GTP", "G")]),                    # 11 CH-06
        linked([("ADE", "A"), ("CYT", "C"), ("GUA", "G")], deoxy=(0, 1, 2)),             # 12 CH-05 DNA names
        linked([("A", "A"), ("C", "C")], deoxy=(1,)),                                    # 13 CH-05 mixed
        linked([("LDA5", "A"), ("LDC", "C"), ("LDG3", "G")], deoxy=(0, 1, 2), mirror=(0, 1, 2)),  # 14 CH-07
        long_names(mk([("A", [("G", 1, "G", False), ("XXX", 2, "C", False), ("XXY", 3, "C", True)])]),
                   {"XXX": "LONG", "XXY": "LDC3"}),                                     # 15 CH-07
        mk([("A", [("G", 9998, "G", False), ("C", 9999, "C", False)]), ("B", [("A", -999, "A", False)])]),
        mk([("", [("0A", 1, "A", False), ("G", 2, "G", True), ("U", 3, "U", False)])],
           header=["HEADER    RNA", "REMARK 950 MDATK BASE-PAIR ANNOTATION", "REMARK 950 PAIR  x",
                   "SEQRES   1 A    1    G", "CRYST1   10.000"]),                       # 17 stale REMARK 950
        mk([("A", [(b, k + 1, b, False) for k, b in enumerate("GCAUGCAUGCAUGCAU")])]),   # 18 long chain
        mk([("A", [("A", 1, "A", True)]), ("B", [("A", 1, "A", False)])],
           extra=[tp.atom_line(900, "MG", "MG", "C", 1, (9, 9, 9), "HETATM", "MG"), "TER",
                  tp.atom_line(901, "CA", "ALA", "D", 1, (8, 8, 8)), tp.atom_line(902, "CB", "ALA", "D", 1, (8, 8, 7)),
                  tp.atom_line(903, "CA", "HIE", "D", 2, (7, 8, 8)), "TER"]),           # 19 protein, ions
    ]
    # 20: a lone 5'-phosphate residue and a P-only first residue (strip_5p)
    lone = [tp.atom_line(1, "P", "PO4", "A", 0, (2.3, 1.8, -0.5), "HETATM"),
            tp.atom_line(2, "OP1", "PO4", "A", 0, (1.9, 3.1, -1.1), "HETATM")]
    out.append("\n".join(lone + mk([("A", [("G", 1, "G", False), ("C", 2, "C", False)])]).splitlines()) + "\n")
    # 21: residues without sugar atoms / flat sugars (undetermined chirality)
    flat = mk([("A", [("G", 1, "G", False), ("C", 2, "C", False)])]).splitlines()
    flat = [ln for ln in flat if not (ln[17:20].strip() == "G" and ln[12:16].strip() in ("C3'", "O4'"))]
    out.append("\n".join(flat) + "\n")
    # 22: ChimeraX-damaged lines and a second model
    good = mk([("C", [("G", 1, "G", False), ("0A", 2, "A", True)])]).splitlines()
    damaged = [tp.chimerax_damage(ln) if ln.startswith("ATOM") else ln for ln in good]
    out.append("\n".join(["MODEL        1"] + damaged[:-1] + ["ENDMDL", "MODEL        2", good[0], "ENDMDL", "END"]) + "\n")
    return out


# --- single functions ------------------------------------------------------------------------------

def scheme_cases():
    out = []
    for scheme in pt.SCHEMES + ("charmm",):
        for base in "ACGU":
            for hand in "DL":
                for term in ("", "5", "3", "N"):
                    try:
                        out.append([base, hand, term, scheme, pt.scheme_name(base, hand, term, scheme)])
                    except pt.PDBError as exc:
                        out.append([base, hand, term, scheme, error(exc)])
    return out


def chiral_cases():
    rng = random.Random(3)
    out = {"volume": [], "sugar": []}
    for _ in range(40):
        pts = [[round(rng.uniform(-3, 3), 3) for _ in range(3)] for _ in range(4)]
        out["volume"].append([pts, pt.chiral_volume(*pts).hex()])
    for base in "ACGU":
        for mirror in (False, True):
            coords = tp.coords(base, mirror)
            variants = [("plain", coords), ("no-base-arg", coords)]
            for centre in ("C4'", "C3'", "C1'"):  # one centre flattened onto its neighbours' plane
                c = dict(coords)
                c[centre] = tuple(sum(coords[n][i] for n in ("O4'", "C2'", "C5'") if n in coords) / 3 for i in range(3))
                variants.append((f"flat {centre}", c))
            for drop in ("C4'", "O3'", "N9", "N1", "C2'"):
                variants.append((f"without {drop}", {k: v for k, v in coords.items() if k != drop}))
            mixed = dict(coords)
            mixed["C1'"] = tuple(-v if i == 0 else v for i, v in enumerate(tp.coords(base, not mirror)["C1'"]))
            variants.append(("moved C1'", mixed))
            jitter = {k: tuple(v + rng.uniform(-0.3, 0.3) for v in xyz) for k, xyz in coords.items()}
            variants.append(("jitter", jitter))
            for label, c in variants:
                call, vols = pt.sugar_chirality(c, None if label == "no-base-arg" else base)
                out["sugar"].append({"base": base, "label": label, "mirror": mirror,
                                     "coords": {k: [float(x).hex() for x in v] for k, v in c.items()},
                                     "base_arg": label != "no-base-arg", "call": call,
                                     "volumes": [[k, v.hex()] for k, v in vols.items()]})
    out["sugar"].append({"base": "A", "label": "empty", "mirror": False, "coords": {}, "base_arg": True,
                         "call": pt.sugar_chirality({}, "A")[0], "volumes": []})
    return out


RANGE_TEXTS = ["", " ", "C:45-73", "C:45", "C", ":45-73", ":45", "C:45-73, A:1-5", "C:45-73;A:1 B", "A:-5--1",
               "A:-5-3", "A:4-1", ":7-2", "AB:1", "A:x", "A:1-", "A:1-2-3", "A::1", "A:1:2", "é:3", "*", "A:+1",
               "A:01-02", "  C:45-73 ,, A:1  ", "A\t:1", "A:1\nB:2"]
RANGE_LISTS = [[["A", 1, 5]], [["A", None, None]], [["B", 3]], [["", 2, 2], ["A", 5, 4]], ["A:1-2", "B"], [["A"]]]


def range_cases():
    out = []
    for value in RANGE_TEXTS + RANGE_LISTS:
        arg = [tuple(x) if isinstance(x, list) else x for x in value] if isinstance(value, list) else value
        try:
            out.append({"spec": value, "ranges": [list(r) for r in pt.parse_ranges(arg)]})
        except pt.PDBError as exc:
            out.append({"spec": value, "error": error(exc)})
    return out


VALUES = ("seqres_records", "fasta", "summary")  # the others return notes
WRITTEN = ("set_records", "strip_5prime_phosphates", "renumber_residues")  # write_pdb recorded


def op_case(texts, k, op, args, run, chirality=True):
    """read_pdb, assign_chirality (detection) unless ``chirality`` is False, then the function."""
    case = {"input": k, "op": op, "args": args}
    try:
        s = pt.read_pdb(texts[k])
        if chirality:
            pt.assign_chirality(s)
        value = run(s)
    except (pt.PDBError, ValueError) as exc:
        case["error"] = error(exc)
        return case
    case["value" if op in VALUES else "notes"] = value
    case["residues"] = residues(s)
    if op in WRITTEN:
        try:
            case["write"] = pt.write_pdb(s)
        except pt.PDBError as exc:
            case["write_error"] = error(exc)
    return case


def operation_cases(texts):
    out = []
    n = len(texts)
    for k in range(n):
        for detect in (True, False):
            for ranges in (None, "B", "A:1-2", ":2", "Z:1", "A:2-1"):
                out.append(op_case(texts, k, "assign_chirality", {"detect": detect, "l_ranges": ranges},
                                   lambda s, d=detect, r=ranges: pt.assign_chirality(s, detect=d, l_ranges=r),
                                   chirality=False))
        for scheme in pt.SCHEMES + ("charmm",):
            for terminals in (True, False):
                for only in (None, "B", "A:2-3"):
                    out.append(op_case(texts, k, "convert_names",
                                       {"scheme": scheme, "only": only, "terminals": terminals},
                                       lambda s, sc=scheme, o=only, t=terminals: pt.convert_names(s, sc, o, t)))
        for policy in pt.RECORD_POLICIES + ("all",):
            out.append(op_case(texts, k, "set_records", {"policy": policy}, lambda s, p=policy: pt.set_records(s, p)))
        out.append(op_case(texts, k, "strip_5prime_phosphates", {}, pt.strip_5prime_phosphates))
        for start, per_chain, gaps in ((1, False, False), (1, True, False), (5, False, True), (-3, True, True),
                                       (9990, False, False), (-2000, False, False)):
            out.append(op_case(texts, k, "renumber_residues", {"start": start, "per_chain": per_chain, "keep_gaps": gaps},
                               lambda s, a=start, b=per_chain, c=gaps: pt.renumber_residues(s, a, b, c)))
        out.append(op_case(texts, k, "seqres_records", {}, pt.seqres_records))
        out.append(op_case(texts, k, "fasta", {"title": "t"}, lambda s: pt.fasta(s, "t")))
        out.append(op_case(texts, k, "summary", {}, pt.summary))
        out.append(op_case(texts, k, "split_chains", {"at": None}, pt.split_chains))
        for at in ("A:2", "A:2, A:3", "B:1", "A", "Z:4", ":2"):
            out.append(op_case(texts, k, "split_chains", {"at": at}, lambda s, a=at: pt.split_chains(s, a)))
        for chains, into in ((["A", "B"], None), (["B", "A"], "Z"), (["A", "Q"], None), ([], None), (["A"], "XY"),
                             (["", "A"], " C ")):
            out.append(op_case(texts, k, "combine_chains", {"chains": chains, "into": into},
                               lambda s, c=chains, i=into: pt.combine_chains(s, c, i)))
        for order in (None, ["B"], ["C", "A", "B"], []):
            out.append(op_case(texts, k, "order_chains", {"order": order}, lambda s, o=order: pt.order_chains(s, o)))
        out.append(op_case(texts, k, "restore_chains", {"assign": True}, lambda s: pt.restore_chains(s)))
        for ids in (["X", "Y"], ["X"], ["X", "YZ"], [" P ", "Q", "R"]):
            out.append(op_case(texts, k, "restore_chains", {"chain_ids": ids},
                               lambda s, i=ids: pt.restore_chains(s, chain_ids=i)))
        for ref in (2, 0):
            for copy_numbers in (False, True):
                out.append(op_case(texts, k, "restore_chains", {"reference": ref, "copy_numbers": copy_numbers},
                                   lambda s, r=ref, c=copy_numbers: pt.restore_chains(s, texts[r], copy_numbers=c)))
        for ranges in (("A:3-4, B", "B:1-2", "A:9", "A:4-3") if k in (0, 1, 4, 5, 11, 14) else ()):
            for scheme, terminals in (("amber", True), ("amber", False), ("ccd", True), ("rosetta", True)):
                case = {"input": k, "op": "to_lrna", "args": {"ranges": ranges, "scheme": scheme, "terminals": terminals}}
                try:
                    r = pt.to_lrna(texts[k], ranges, scheme=scheme, terminals=terminals)
                    case.update(text=r.text, notes=r.notes)
                except pt.PDBError as exc:
                    case["error"] = error(exc)
                out.append(case)
    return out


# --- clean ----------------------------------------------------------------------------------------

POOL = {
    "scheme": [None, "amber", "ccd", "rosetta", "pdb"],
    "records": [None, "keep", "rosetta", "pdb"],
    "strip_conect": [True, False],
    "detect_l": [True, False],
    "l_ranges": [None, "B", "A:1-2", "Z:1"],
    "chains": [None, "assign", "ids", "reference", "reference+numbers"],
    "split_segments": [False, True],
    "split_at": [None, "A:2"],
    "combine": [None, ["A", "B"]],
    "combine_into": [None, "Q"],
    "sort_chains": [False, True, ["B"]],
    "renumber": [None, (1, False, False), (1, True, False), (10, False, True)],
    "seqres": [False, True],
    "make_fasta": [False, True],
    "terminal_names": [True, False],
    "strip_5p": [False, True],
}


def clean_options(choice, texts, k):
    kw = {}
    for key in ("scheme", "records", "strip_conect", "detect_l", "l_ranges", "split_segments", "split_at", "combine",
                "combine_into", "sort_chains", "seqres", "make_fasta", "terminal_names", "strip_5p"):
        if key in choice:
            kw[key] = choice[key]
    chains = choice.get("chains")
    if chains == "assign":
        kw["assign_chains"] = True
    elif chains == "ids":
        kw["chain_ids"] = ["X", "Y"]
    elif chains in ("reference", "reference+numbers"):
        kw["chain_reference"] = texts[2]
        kw["copy_numbers"] = chains == "reference+numbers"
    if choice.get("renumber"):
        start, per_chain, gaps = choice["renumber"]
        kw.update(renumber=True, start=start, per_chain=per_chain, keep_gaps=gaps)
    return kw


def clean_case(texts, k, kw, label):
    shown = dict(kw)
    if "chain_reference" in shown:
        shown["chain_reference"] = texts.index(kw["chain_reference"])  # the input index
    case = {"input": k, "label": label, "options": shown}
    try:
        r = pt.clean(texts[k], **kw)
        case.update(text=r.text, notes=r.notes, fasta=r.fasta)
    except (pt.PDBError, ValueError) as exc:
        case["error"] = error(exc)
    return case


def clean_cases(texts):
    out = []
    mdatk_tests = [  # tests/test_pdbtools.py's calls
        (0, {"scheme": "amber"}), (0, {"scheme": "ccd"}), (0, {"scheme": "rosetta"}), (0, {"scheme": "pdb"}),
        (4, {"scheme": "amber"}), (4, {"scheme": "amber", "detect_l": False}), (7, {"sort_chains": True}),
        (7, {"sort_chains": ["B"]}), (1, {"chain_reference": texts[2]}), (1, {"assign_chains": True}),
        (3, {"renumber": True}), (3, {"renumber": True, "keep_gaps": True}),
        (3, {"renumber": True, "per_chain": True, "start": 5}), (8, {"strip_conect": False}),
        (9, {"strip_conect": False}), (10, {"scheme": "amber", "strip_5p": True}), (20, {"strip_5p": True}),
        (11, {"scheme": "amber"}), (11, {"scheme": "ccd", "seqres": True}), (12, {"scheme": "amber"}),
        (13, {"scheme": "amber", "make_fasta": True}), (14, {"seqres": True, "make_fasta": True}),
        (15, {"seqres": True}), (16, {"renumber": True, "start": -2000}), (16, {"renumber": True, "start": 9990}),
        (16, {"scheme": "amber"}), (17, {"scheme": "ccd", "seqres": True, "make_fasta": True, "title": "x"}),
        (18, {"seqres": True}), (19, {"seqres": True, "make_fasta": True, "scheme": "pdb"}),
        (0, {"combine": ["A", "B"], "combine_into": "XY"}), (0, {"l_ranges": "A:4-1"}),
        (0, {"combine": ["A", "B"], "renumber": True}), (6, {"split_segments": True}),
        (6, {"split_at": "A:2", "assign_chains": True}), (22, {}), (22, {"scheme": "amber", "seqres": True}),
        (1, {"chain_reference": texts[0]}), (1, {"chain_ids": ["A"]}), (21, {"scheme": "amber"}),
    ]
    for k, kw in mdatk_tests:
        out.append(clean_case(texts, k, kw, "test"))
    rng = random.Random(17)
    for k in range(len(texts)):
        for _ in range(6):
            choice = {key: rng.choice(values) for key, values in POOL.items() if rng.random() < 0.35}
            out.append(clean_case(texts, k, clean_options(choice, texts, k), "grid"))
    return out


def main():
    texts = inputs()
    data = {"inputs": texts, "scheme_name": scheme_cases(), "chiral": chiral_cases(), "ranges": range_cases(),
            "operations": operation_cases(texts), "clean": clean_cases(texts)}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "pdbtools_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
