"""What MDATK's PDB tools make of PDB text (tests/data/io/pdb/pdb_cases.json).

    python3 -I -B make_pdb_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy and pytest (MDATK's tests/test_pdbtools.py builds the fixtures); nothing is written
into the checkout. The output is reproducible (fixed random seeds). Writes pdb_cases.json:

  read        read_pdb of PDB text: the structure (header, segments, atoms with xyz, coords(),
              label, trailer, CONECT, MODEL, notes, chain_ids) or the error, and write_pdb of it
              with CONECT stripped / kept and with SEQRES records (the chains and three-character
              names MDATK's hygiene passes, given as data), or the error of writing it
              ("write_error"); MDATK's test fixtures, ChimeraX damage, models, ANISOU /
              SIGATM / SIGUIJ, MASTER, records after END, CONECT variants, serial and residue
              number wraps, lower-case records, Unicode, line ends, broken fields, random mutations
  files       read_pdb of a file (bytes as hex): strict UTF-8, universal newlines
  residue_names parse_residue_name of every name of up to three characters over an alphabet (with
              MDBridge's L-DNA names, NA-01, and the CCD's CH, CH-22: improvements_nucleotides.py)
  insert      annotate._insert_records
  transform   builder.transform_pdb_text (ANISOU included), or None / the error
  tail        builder._tail, or the error

Floats are float.hex() strings; errors are {"type": class name, "message": str(exc)}. MDBridge's
chem improvements are applied (improvements_chem.py, which includes NA-01's names): CH-01 (a residue
number below -999 is a write error) and CH-07 (SEQRES names have three characters) change what
write_pdb gives here.
"""

import importlib.util
import itertools
import json
import math
import os
import random
import sys
import tempfile

if len(sys.argv) != 3:
    sys.exit("usage: make_pdb_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
sys.path.insert(0, CHECKOUT)

import numpy as np  # noqa: E402

from mdatk import annotate, builder  # noqa: E402
from mdatk import pdbtools as pt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import improvements_chem  # noqa: E402  (CH-01, CH-07 in write_pdb; NA-01's names)
import improvements_nucleotides  # noqa: E402  (NA-01: the L-DNA names of parse_residue_name)

improvements_chem.apply_mdatk()

spec = importlib.util.spec_from_file_location("mdatk_test_pdbtools",
                                              os.path.join(CHECKOUT, "tests", "test_pdbtools.py"))
tp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tp)

OUT = os.path.abspath(sys.argv[2])


def fhex(x):
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x.hex()


def error(exc):
    return {"type": type(exc).__name__, "message": str(exc)}


def atom_dict(a):
    """[record, name, alt, tail, element, charge, serial, anisou]"""
    return [a.record, a.name, a.alt, a.tail, a.element, a.charge, a.serial, a.anisou]


def structure_dict(s):
    return {
        "header": s.header,
        # residues: [name, chain, resseq, icode, label, coords, atoms]
        "segments": [[[r.name, r.chain, r.resseq, r.icode, r.label,
                       [[k] + [fhex(v) for v in xyz] for k, xyz in r.coords().items()],
                       [atom_dict(a) for a in r.atoms]] for r in seg] for seg in s.segments],
        "trailer": s.trailer,
        "conect": s.conect,
        "model": s.model,
        "notes": s.notes,
        "chain_ids": s.chain_ids(),
    }


def read_case(text):
    case = {"text": text}
    try:
        s = pt.read_pdb(text)
    except Exception as exc:  # noqa: BLE001
        case["error"] = error(exc)
        return case
    case["structure"] = structure_dict(s)
    try:
        case["write"] = pt.write_pdb(s)
    except pt.PDBError as exc:  # CH-01
        case["write_error"] = error(exc)
        return case
    conect = pt.write_pdb(s, strip_conect=False)
    if conect != case["write"]:
        case["write_conect"] = conect  # else equal to "write"
    case["seqres_chains"] = [[c, names] for c, names in pt._mdbridge_seqres_names(s)[0].items()]  # CH-07
    case["write_seqres"] = pt.write_pdb(s, strip_conect=False, seqres=True)
    return case


# --- inputs ---------------------------------------------------------------------------------------

def fixtures():
    out = []
    segs = [("C", [("G", 1, "G", True), ("0A", 2, "A", True), ("LG5", 3, "G", True)]),
            ("B", [("C", 4, "C", False), ("U", 5, "U", False)])]
    good = tp.make_pdb(segs)
    out.append(good)
    damaged = "\n".join(tp.chimerax_damage(ln) if ln.startswith(("ATOM", "TER")) else ln
                        for ln in good.splitlines())
    out.append(damaged)
    line = good.splitlines()[0]
    out.append(line[:17] + "  LG5C:" + line[22:])
    out.append("ATOM      1  P     0GC:  7   -2.314   1.777   -0.478 1.00 0.00 P")
    out.append(line[:17] + "G   C" + line[22:])
    out.append("\n".join([
        "ATOM      7 C1'    G A   1       1.000   2.000   3.000  1.00  0.00",
        "ATOM      8 H5''   G A   1       1.000   2.000   3.000  1.00  0.00",
        "HETATM    9 MG    MG A   2       4.000   5.000   6.000  1.00  0.00",
        "HETATM   10 Na+  Na+ A   3       4.000   5.000   6.000  1.00  0.00",
        "ATOM     11  CA  ALA B   1       0.000   0.000   0.000  1.00  0.00           C",
    ]))
    body = tp.make_pdb([("A", [("G", 1, "G", False)])]).splitlines()[:-1]
    second = [ln.replace("  1.00  0.00", "  1.00  9.99") for ln in body]
    anisou = "ANISOU    1  P     G A   1     1000   2000   3000    100    200    300       P"
    out.append("\n".join(["HEADER    RNA", "REMARK   1 kept", "SEQRES   1 A    1    G", "CRYST1   10.000",
                          "MODEL        1", body[0], anisou] + body[1:] + ["ENDMDL", "MODEL        2"]
                         + second + ["ENDMDL", "MASTER        0", "END"]))
    out.append("REMARK only\nEND\n")
    water = [tp.atom_line(90, "O", "HOH", "W", 1, (9, 9, 9), "HETATM"),
             tp.atom_line(91, "H1", "HOH", "W", 1, (9, 9, 9), "HETATM"),
             tp.atom_line(92, "MG", "MG", "W", 2, (8, 8, 8), "HETATM", "MG"), "TER"]
    out.append(tp.make_pdb([("A", [("G", 1, "G", False), ("0A", 2, "A", True), ("DG", 3, "G", False)])],
                           records={"0A": "HETATM", "G": "HETATM"}, extra=water))
    out.append(tp.make_pdb([("A", [("G", 1, "G", False)])], first_serial=101,
                           extra=["CONECT  101  102  103", "CONECT  102  101"]))
    out.append(tp.make_pdb(tp.CPPTRAJ, extra=tp.IONS))
    out.append(tp.make_pdb(tp.REFERENCE))
    out.append(tp.make_pdb([("", tp.CPPTRAJ[0][1] + tp.CPPTRAJ[1][1])]))
    out.append(tp.numbered(tp.GAPPED))
    seq = "GCAUGCAUGCAUGCA"
    out.append(tp.make_pdb([("A", [(b, k + 1, b, False) for k, b in enumerate(seq)]),
                            ("B", [("0G", 1, "G", True), ("0C", 2, "C", True)])],
                           header=["REMARK   1 note", "SEQRES   1 Z    1    A", "CRYST1   10.000"],
                           extra=tp.IONS))
    out.append(tp.make_pdb([("A", [("G", 1, "G", False), ("C", 2, "C", False), ("LA", 3, "A", True),
                                   ("LU", 4, "U", True)]),
                            ("B", [("G", 1, "G", True), ("C", 2, "C", True)]),
                            ("C", [("DG", 1, "G", False), ("PSU", 2, "U", False)])], extra=tp.IONS))
    return out


A = "ATOM      1  P     G A   1       1.000   2.000   3.000  1.00  0.00           P"
B = "ATOM      2  OP1   G A   1       1.500   2.500   3.500  1.00 10.00           O"
C = "HETATM    3 MG    MG B   2      -4.000  15.250   0.125  0.50  5.00          MG2+"


def crafted():
    out = [
        "",
        "\n\n",
        A,
        A + "\n" + B + "\nTER\n" + C + "\nEND\n",
        "atom      1  P     G A   1       1.000   2.000   3.000  1.00  0.00           P\n"
        "hetatm    2 MG    MG B   2       1.000   2.000   3.000\nter\nend",
        "  ATOM    1  P     G A   1       1.000   2.000   3.000",  # leading spaces
        A.replace("ATOM  ", "ATOM\t "),
        A + "\r\n" + B + "\r\n",
        A + "\r" + B + "\r",
        A + "\x0c" + B + "\x1c" + C + " END",
        A + "\x85" + B,
        "HEADER    TEST\nTITLE     x\n" + A + "\nANISOU    1  P     G A   1     1   2   3\n"
        "SIGATM    1  P     G A   1       0.010   0.010   0.010\nSIGUIJ    1\n" + B + "\nEND\n",
        "ANISOU    1  P     G A   1     1000   2000   3000    100    200    300       P\n" + A,
        A + "\nTER\nANISOU    1  P     G A   1     1000\n" + B,
        A + "\nMASTER    1\nMASTER    2\n" + B,
        A + "\nEND\n" + B + "\n" + C + "\nREMARK after\n",
        A + "\nENDMDL\n" + B + "\nMODEL 3\n" + C,
        "MODEL        1\n" + A + "\nENDMDL\nMODEL        2\n" + B + "\nENDMDL\nMODEL 3\n" + C + "\nEND",
        "MODEL        1\n" + A + "\nEND\n",
        A + "\nCONECT    1    2\nCONECT    2    1    3\nCONECT    3\nCONECT    x    1\n"
            "CONECT    1    9\nCONECT    9    1\nCONECT\nconect    1    2    3    1    2    3\n" + B + "\n" + C,
        A + "\n" + A + "\n" + B + "\nCONECT    1    2\n",  # duplicate serial: CONECT dropped
        A.replace("    1  P", "*****  P") + "\n" + B + "\nCONECT    2    1",
        A.replace("   1    ", "12345   ") + "\n" + B.replace("   1    ", "99999   "),
        A.replace("G A   1", "G A9999") + "\n" + B.replace("G A   1", "G A-999"),
        A.replace("G A   1", "G A-1234"),
        A.replace("G A   1", "G A   x"),
        A.replace("G A   1", "G A    "),
        A.replace("G A   1", "G A ١٢"),
        A.replace("1.000", "1.0x0"),
        A.replace("1.000   2.000   3.000", "                        "),
        A.replace("1.000   2.000", "  1e0   2_0  "),
        A.replace("ATOM      1", "ATOM  1_000"),
        A.replace("ATOM      1", "ATOM   -001"),
        A + "\n" + A.replace("G A   1 ", "G A   1B") + "\n" + A.replace(" G A", "GG A"),
        A.replace(" G A   1", "LONG   1") + "\n" + A.replace(" G A   1", "ABCDA  1"),
        A.replace("  P   ", " P1234") + "\n" + A.replace("  P   ", " H5'' "),
        A.replace("  P   ", "  Na+ ").replace(" G A", "Na+ A").replace("           P", ""),
        A.replace("  P   ", "  CL  ").replace(" G A", " CL A").replace("           P", ""),
        A.replace("  P   ", "  1HB ").replace("           P", ""),
        A.replace("  P   ", " 12   ").replace("           P", ""),
        A.replace("  P   ", "  éX  ").replace("           P", ""),
        A.replace("  P   ", "  ßa  ").replace("           P", ""),
        A.replace("  P   ", "  ſ   ").replace("           P", ""),
        A + "\n" + A.replace(" A   1", "     1") + "\n" + A.replace(" A   1", " a   1"),
        A.replace(" G A", "αβγ A") + " é" + "\n" + "MASTEr  x\nmaſter\n" + B,
        "REMARK 950 a\nMASTER\nMODRES x\nSEQRES   1 Q\nHELIX  1\n" + A + "\nTER\nCRYST1  1\n",
        "SEQRES x\nREMARK 1\nseqres y\nscale1 x\n" + A,
        "﻿" + A,
        A + "  1.00  0.00   extra columns beyond eighty ......",
        A[:54],
        A[:30],
        A + "\n" + "TER       2      G A   1\n" + B.replace(" A   1", " B   2"),
        A + "\n" + B.replace(" A   1", " B   1") + "\n" + C.replace(" B   2", " A   2"),
        # ChimeraX-damaged variants
        "ATOM      1  P     GA:   1       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P    AGA:   1       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P    BGA:   1       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P   A GA:   1       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P    A GA:   1       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P     GA:  1A   1.000 2.000 3.000 0.5 7.25 P",
        "ATOM      1  P     GA:  1A   1.000 2.000 3.000 0.5 P",
        "ATOM      1  P     GA:  -12  1.000 2.000 3.000",
        "ATOM      1  P     GA:  12 1.000 2.000",
        "ATOM      1  P     GA:  x  1.000 2.000 3.000",
        "ATOM      1  P     GA:  12b 1.0004 -2.0006 3.12345 0.25 abc",
        "ATOM      1  P     GA:  12 1.000 2.000 3.000 1.5 2.5 XYZ",
        "ATOM      1  P     GA:  12 1.000 2.000 3.000 1.5 2.5 é",
        "ATOM      1  P     GA:  ١٢ ١.٠٠٠ 2.000 3.000",
        "ATOM      1  P     GA :  1       1.000   2.000   3.000",
        "ATOM      1  P     G  A:  1       1.000   2.000   3.000",
        "ATOM      1  P  LONGNAMEA:  1   1.000   2.000   3.000",
        "ATOM      1  P                 :  1   1.000   2.000   3.000",
        "ATOM      1  P     LG5A:   3       1.000   2.000   3.000  1.00  0.00           P  \n"
        "ANISOU    1  P     LG5A:   3     1000   2000   3000    100    200    300       P",
        "ATOM      1  P       G5A:   3       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P         GA:   3       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P     GA:9999A       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      1  P     GA:12345       1.000   2.000   3.000  1.00  0.00           P",
    ]
    return out


def mutations(base, n, seed):
    rng = random.Random(seed)
    alphabet = " -.0123456789ABCabcATOMHETATMTERENDMODELCONECTANISOU:\t*'+é"
    out = []
    lines = base.splitlines()
    for _ in range(n):
        text = list(rng.choice(lines))
        mutated = [ln for ln in lines]
        k = rng.randrange(len(mutated))
        chars = list(mutated[k])
        for _ in range(rng.randint(1, 4)):
            pos = rng.randrange(max(1, len(chars) + 1))
            op = rng.random()
            if op < 0.5 and pos < len(chars):
                chars[pos] = rng.choice(alphabet)
            elif op < 0.75:
                chars.insert(pos, rng.choice(alphabet))
            elif pos < len(chars):
                del chars[pos]
        mutated[k] = "".join(chars)
        del text
        out.append("\n".join(mutated))
    return out


def read_cases():
    texts = fixtures() + crafted()
    base = "\n".join([
        "HEADER    X", "MODEL        1", A, "ANISOU    1  P     G A   1     1000   2000   3000"
        "    100    200    300       P", B, "TER", C, "ENDMDL", "CONECT    1    2", "END"])
    texts += mutations(base, 120, 7)
    damaged = "\n".join([
        "ATOM      1  P     GA:   1       1.000   2.000   3.000  1.00  0.00           P",
        "ATOM      2  OP1  LG5B:  12A      1.500   2.500   3.500  1.00 10.00           O",
        "ATOM      3  C1'A GC:   3 1.000 2.000 3.000 0.5 7.25 C"])
    texts += mutations(damaged, 60, 11)
    return [read_case(t) for t in texts]


def file_cases(tmp):
    raws = [
        A.encode() + b"\r\n" + B.encode() + b"\rEND\r\n",
        A.encode() + b"\n" + "HETATM    3 Mé    MG B   2  ".encode() + b"\n",
        A.encode() + b"\nREMARK \xe9t\xe9\n",
        A.encode() + b"\nREMARK \xe2\x82\n",
        A.encode() + b"\nREMARK \xe2\x82",
        A.encode() + b"\nREMARK \xff",
        A.encode() + b"\nREMARK \xf0\x9f\x98a",
        b"\xef\xbb\xbf" + A.encode(),
        b"",
        b"REMARK only\n",
    ]
    out = []
    for k, raw in enumerate(raws):
        path = os.path.join(tmp, f"case{k}.pdb")
        with open(path, "wb") as fh:
            fh.write(raw)
        case = {"bytes": raw.hex()}
        try:
            case["structure"] = structure_dict(pt.read_pdb(path))
        except Exception as exc:  # noqa: BLE001
            case["error"] = error(exc)
        out.append(case)
    return out


def residue_name_cases():
    alphabet = "ACGTU RL0D53NXalgu"
    names = {""}
    for n in (1, 2, 3):
        names.update("".join(p) for p in itertools.product(alphabet, repeat=n))
    names.update(["ADE", "cyt", " GUA ", "URA", "URI", "THY", "thy", "RAN", "LA5 ", "0A3", "DAN", "DT5",
                  "LT", "0T", "RT", "ADEN", "ıA", "ɐ", "A5\n", "\tG\t", "AAAA", "ﬀ",
                  "LDA5", "LDT3", "LDGN", "ldc", " 0dt ", "LDU", "0DU", "LDA6", "0DAN",
                  "CH", " ch ", "CH5", "CH3", "CHN", "RCH", "0CH", "DCH", "HC"])  # CH-22
    out = []
    for name in sorted(names):
        r = improvements_nucleotides.parse_residue_name(name)  # NA-01: also L-DNA names
        out.append([name, None if r is None else [r.base, r.sugar, r.chirality, r.terminal]])
    return out


def insert_cases():
    texts = [
        "",
        "HEADER    X\nREMARK 950 old\nREMARK 950 older\nREMARK   1 kept\n" + A + "\nMASTER   1\nEND\n",
        "HEADER    X\r\n" + A + "\r\nEND",
        "REMARK 1\nmodel 1\n" + A,
        "  ATOM\n" + A,
        "REMARK 950X\nREMARK 95\nmaster\n",
        "HEADER\x0cTITLE " + A,
        "HEADER    no atoms\n",
    ]
    records = [[], ["REMARK 950 ENGINE builtin"], ["REMARK 950 A", "REMARK 950 B é"]]
    return [{"text": t, "records": r, "out": annotate._insert_records(t, r)}
            for t in texts for r in records]


RANDOM_ROTATIONS = [
    [["-0x1.830eb884349f2p-1", "-0x1.77af99835e17dp-3", "-0x1.41ba80a7f9385p-1"],
     ["0x1.95dc8132d8befp-2", "-0x1.c8d3f3949c735p-1", "-0x1.bb1ae33e8ee0fp-3"],
     ["-0x1.0abd276142039p-1", "-0x1.a6855e2fa18ffp-2", "0x1.7e931d708c2adp-1"]],
    [["-0x1.460693f4a361ep-1", "-0x1.4795f1a0189fdp-2", "0x1.6732aadb9c9d8p-1"],
     ["0x1.7e36f88c2d640p-2", "-0x1.d91c8bae8b216p-1", "-0x1.523a837796c43p-4"],
     ["0x1.5970cc65b670cp-1", "0x1.ac9aef5076f5dp-3", "0x1.6a6643c4d6258p-1"]],
    [["-0x1.596c4a8c8736ap-1", "-0x1.06256cb3dd60ap-4", "0x1.7880512bf1e36p-1"],
     ["0x1.2fc1620ed1cb8p-1", "-0x1.478d0a38351c5p-1", "0x1.f45918d062f65p-2"],
     ["0x1.c1b59017d4b9ap-2", "0x1.882612f27de52p-1", "0x1.e0d92bf1cd43dp-2"]],
    [["-0x1.859d84c42fb30p-3", "0x1.4ffd242fdb6c3p-1", "-0x1.75db079c5da2bp-1"],
     ["-0x1.e935230c6c0cbp-1", "-0x1.2db09f7e88de3p-2", "-0x1.037b9fae7b7d6p-6"],
     ["-0x1.cddd6d30e97abp-3", "0x1.63ab9305a4819p-1", "0x1.5dba94061e992p-1"]],
]


def transform_cases():
    rng = np.random.default_rng(5)
    texts = [
        "HEADER    X\n" + A + "\n" + B + "\r\n" + C + "\rTER\nEND",
        A + "\nANISOU    1  P     G A   1     1000   2000   3000    100    200    300       P\r\n" + B,
        A + "\nANISOU    1  P     G A   1     1000   2000   3000    100    200    x\n",
        A + "\nANISOU    1  P     G A   1     1000   2000   3000    100    200\n",
        A + "\nANISOU    1  P     G A   1      -15     -3      7      1     -1      2\n",
        "atom      1" + A[11:] + "\n" + A.replace("ATOM  ", "ATOM"),
        A.replace("1.000", "1.0x0"),
        A[:40],
        A + "\x0c" + B + " " + C + "\x85",
        A.replace("   1.000   2.000   3.000", "9999.0009999.0009999.000"),
        A.replace("   1.000   2.000   3.000", "-999.000-999.000-999.000"),
        A.replace("   1.000", "     nan"),
        "",
        "REMARK only\n",
    ]
    rotations = [np.eye(3), np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])]
    # Four random rotations: the Q of np.linalg.qr(rng.normal(size=(3, 3))) with det(Q) made +1, as
    # computed once with OpenBLAS's SkylakeX kernels. LAPACK's QR differs in the last bits between
    # CPU kernels, so the matrices are written out (exactly) to keep the cases machine-independent;
    # the normal draws are still taken so the translations below come from the same stream.
    for _ in range(4):
        rng.normal(size=(3, 3))
    rotations += [np.array([[float.fromhex(v) for v in row] for row in r]) for r in RANDOM_ROTATIONS]
    translations = [np.zeros(3), np.array([10.0, -20.0, 0.0005]), rng.normal(size=3) * 50]
    out = []
    for t in texts:
        for r in rotations:
            for tr in translations:
                case = {"text": t, "rotation": [[fhex(v) for v in row] for row in r],
                        "translation": [fhex(v) for v in tr]}
                try:
                    case["out"] = builder.transform_pdb_text(t, r, tr)
                except Exception as exc:  # noqa: BLE001
                    case["error"] = error(exc)
                out.append(case)
    big = np.array([[5000.0, 0, 0], [0, 1, 0], [0, 0, 1]])
    for t in texts[:2]:
        case = {"text": t, "rotation": [[fhex(v) for v in row] for row in big],
                "translation": [fhex(0.0)] * 3}
        try:
            case["out"] = builder.transform_pdb_text(t, big, np.zeros(3))
        except Exception as exc:  # noqa: BLE001
            case["error"] = error(exc)
        out.append(case)
    return out


def tail_cases():
    values = [(1.0, 2.0, 3.0), (-999.9994, 9999.9994, 0.0), (-999.9995, 0, 0), (0, 9999.9995, 0),
              (0, -9999.9995, 0), (0.0005, 0.0015, 0.0025), (-0.0004, 1e-9, -1e-9),
              (float("nan"), 1e6, 0), (1e6, float("nan"), 0), (float("inf"), 0, 0), (0, float("-inf"), 0),
              (123.4565, -12.3455, 1.0005)]
    out = []
    for xyz in values:
        for occ, bf in ((1.0, 0.0), (0.125, 99.995), (100.0, -12.5)):
            case = {"xyz": [fhex(v) for v in xyz], "occupancy": fhex(occ), "bfactor": fhex(bf)}
            try:
                case["out"] = builder._tail(xyz, occ, bf)
            except Exception as exc:  # noqa: BLE001
                case["error"] = error(exc)
            out.append(case)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pdb-cases-") as tmp:
        files = file_cases(tmp)
    data = {"read": read_cases(), "files": files, "residue_names": residue_name_cases(),
            "insert": insert_cases(), "transform": transform_cases(), "tail": tail_cases()}
    with open(os.path.join(OUT, "pdb_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
