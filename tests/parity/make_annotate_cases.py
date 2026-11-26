"""What MDATK's base-pair annotation (mdatk/annotate.py) writes, for MDBridge's chem::annotate
(tests/chem/test_annotate.cpp; docs/specs/m6-plan.md 1.2, 4.2).

    python3 -I -B make_annotate_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy (the builtin pair detector, and the structure builder that makes the helices).
MDBridge's chem improvements are applied first (improvements_chem.py: CH-14 'auto' / 'barnaba' as
the builtin engine with the warning, D13 "REMARK 950 MDBridge BASE-PAIR ANNOTATION", the
'mdbridge pdb clean' hint, NA-01). No Barnaba is involved (DIFFERENCES "No Barnaba"). The output is
reproducible (no randomness). Writes OUT_DIR/annotate_cases.json:

  inputs   PDB texts (Watson-Crick helices built by builder.build_helix: D, L, a single strand,
           mismatches, G·U wobbles, several strand groups, blank chain ids, ions, two models,
           ChimeraX-damaged lines, a 70-nucleotide duplex for the 60-column wrap, numbering gaps,
           an earlier annotation and MASTER, DNA, no nucleotides)
  cases    annotate(inputs[k], pairs, stacking, dotbracket, engine, cutoff, wobble): records, notes,
           engine and the text (or the error); "insert": true where the text is
           _insert_records(input, records) (the text is then not repeated)
  structure annotate of a Structure (read_pdb of an input): the text written through write_pdb

Errors are {"type": class name, "message": str(exc)}.
"""

import json
import os
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_annotate_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
sys.path.insert(0, CHECKOUT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_common  # noqa: E402

chem_common.pin_kernels()  # before numpy (the helices' products; chem_common.py)
import improvements_chem  # noqa: E402

improvements_chem.apply_mdatk()
from mdatk import annotate as an  # noqa: E402
from mdatk import builder as B  # noqa: E402
from mdatk import pdbtools as pt  # noqa: E402


error = chem_common.error


def shifted(text, dx=0.0, chains=None):
    """ATOM/HETATM/TER lines moved by dx along x, chain ids renamed (old -> new)."""
    out = []
    for ln in text.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            x = float(ln[30:38]) + dx
            ln = ln[:30] + f"{x:8.3f}" + ln[38:]
        if ln.startswith(("ATOM", "HETATM", "TER")) and chains and len(ln) > 21 and ln[21] in chains:
            ln = ln[:21] + chains[ln[21]] + ln[22:]
        out.append(ln)
    return out


def body(text):
    return [ln for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM", "TER"))]


def inputs():
    d = B.build_helix("GGGAUGCC")
    out = [d,                                                          # 0 D duplex
           B.build_helix("GGGAUGCC", chirality="L"),                   # 1 L duplex
           B.build_helix("GGCAUC", partner="none"),                    # 2 single strand
           B.build_helix("GGCAUC", partner="GAUGUC"),                  # 3 mismatches
           B.build_helix("GGUAGCC", partner="GGUUAUC"),                # 4 G·U wobbles
           B.build_helix("GGGAUGCC", scheme="amber", hydrogens=True)]  # 5 amber names, hydrogens
    two = shifted(B.build_helix("GCGC"), 40.0, {"A": "C", "B": "D"})
    out.append("\n".join(body(d) + [ln for ln in two if ln.startswith(("ATOM", "TER"))] + ["END"]) + "\n")  # 6
    blank = shifted(d, 0.0, {"A": " ", "B": " "})
    out.append("\n".join(blank) + "\n")                                # 7 blank chain ids (cpptraj)
    ions = body(d) + ["HETATM 9001 MG    MG M   1      30.000  30.000  30.000  1.00  0.00          MG",
                      "TER", "HETATM 9002  O   HOH W   1      31.000  30.000  30.000  1.00  0.00           O",
                      "TER", "END"]
    out.append("\n".join(ions) + "\n")                                 # 8 ions and water
    out.append("\n".join(["HEADER    RNA", "MODEL        1"] + body(d) + ["ENDMDL", "MODEL        2"]
                         + shifted("\n".join(body(d)), 1.0) + ["ENDMDL", "END"]) + "\n")   # 9 two models
    damaged = []
    for ln in d.splitlines():
        if ln.startswith("ATOM") and ln[21] == "B":
            fused = f"{ln[17:20].strip()}{ln[21]}:"
            ln = ln[:22 - len(fused)] + fused + ln[22:]
        damaged.append(ln)
    out.append("\n".join(damaged) + "\n")                              # 10 ChimeraX lines
    long = "GGCAUCGAUCGGAUCCGAUGCAAGCUUGCAUGCCAUGGCAUCGAUCCGAUCGAUGCAUCGAUCGAUCGGC"
    out.append(B.build_helix(long))                                    # 11 70 bp: the 60-column wrap
    gapped = pt.read_pdb(d)
    for k, r in enumerate(gapped.segments[0]):
        r.resseq = 1 + k + (5 if k >= 4 else 0)
    out.append(pt.write_pdb(gapped))                                   # 12 numbering gaps
    out.append("\n".join(["HEADER    RNA", "REMARK 950 OLD ANNOTATION LINE", "REMARK   2 kept", "MASTER    1 2 3"]
                         + body(d) + ["CONECT    1    2", "END"]) + "\n")  # 13 earlier annotation, MASTER
    dna = pt.read_pdb(B.build_helix("GCAU"))
    for r in dna.residues():
        r.name = "D" + r.name.replace("U", "T")
        r.atoms = [a for a in r.atoms if a.name != "O2'"]
    out.append(pt.write_pdb(dna))                                      # 14 DNA names
    out.append("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")  # 15
    out.append("REMARK only\n")                                        # 16 no atoms
    l_single = B.build_helix("GAC", partner="none", chirality="L", scheme="amber")
    out.append(l_single)                                               # 17 L names (amber)
    return out


OPTIONS = [
    {},
    {"engine": "builtin"},
    {"engine": "barnaba"},
    {"engine": "auto", "stacking": False},
    {"pairs": False},
    {"dotbracket": False},
    {"pairs": False, "stacking": False, "dotbracket": False},
    {"engine": "barnaba", "stacking": False, "dotbracket": False},
    {"wobble": False},
    {"cutoff": 2.5},
    {"cutoff": 4.25, "wobble": False},
    {"engine": "dssr"},
]


def annotate_case(texts, k, options):
    case = {"input": k, "options": options}
    try:
        res = an.annotate(texts[k], **options)
    except (pt.PDBError, ValueError) as exc:
        case["error"] = error(exc)
        return case
    case.update(records=res.records, notes=res.notes, engine=res.engine)
    if res.text == an._insert_records(texts[k], res.records):
        case["insert"] = True
    else:
        case["text"] = res.text
    return case


def main():
    texts = inputs()
    cases = [annotate_case(texts, k, o) for k in range(len(texts)) for o in OPTIONS]
    structure = []
    for k in (0, 2, 6, 9, 13):
        res = an.annotate(pt.read_pdb(texts[k]))
        structure.append({"input": k, "text": res.text, "notes": res.notes, "records": res.records})
    data = {"inputs": texts, "cases": cases, "structure": structure}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "annotate_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
