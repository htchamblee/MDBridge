"""Cut the microMg / nanoMg Mg2+ parameter sets out of the authors' GROMACS parameter files.

    python3 -I tools/extract_mg_params.py DIR [--out OUT_DIR]

DIR holds the three repositories of Grotz, Cruz-Leon & Schwierz (github.com/bio-phys) at the
commits MDATK 0.0.4.1 cites in mdatk/workflows/mg_params.py (each set's `source`), as clones or
extracted archives named after the repository:

    DIR/Magnesium-FFs   e71379d1a88b54221018af53ca54d3c18066ba05  (TIP3P)
    DIR/optimizedMgFFs  1df53e9c94fc14a97864ac35c6951085bb60d1c0  (SPC/E, TIP3P-FB, TIP4P-Ew,
                                                                    TIP4P-D, TIP4P/2005)
    DIR/MgFF_OPC        cf8148a3131c198b19f3a2896b68bad6f6f62f5b  (OPC)

The SHA-256 of every file read is checked (SHA256 below). Only numbers are kept, as the files write
them, each with the file and line it was read from: per (model, water) the ion's sigma / epsilon
(`[ atomtypes ]`), the chloride model the set was optimised with, its sigma / epsilon, the
Mg2+-Cl- pair and every Mg2+-RNA pair (`[ nonbond_params ]`). Nothing else of the files (comments,
the copied RNA atom types, molecule definitions) is written. Units: sigma in nm, epsilon in kJ/mol.

Writes OUT_DIR/mg_sets.tsv and OUT_DIR/mg_rna_pairs.tsv (default OUT_DIR: data/sources/mg).
The repositories state no licence; the numbers are recorded as facts (owner's decision of
2026-11-28, data/sources/mg/README.md). tools/gen_workflow_tables.py builds the Mg2+ tables from them.

Checks (printed; a failure stops the run): every Mg2+-RNA pair of a set is the Lorentz-Berthelot
pair of the ion and the RNA type (the file's own `[ atomtypes ]` row of that type) scaled by one
lambda_sigma and one lambda_epsilon for the whole set, to the digits the file prints; those two
factors are printed next to the ones the papers' Table 1 prints (PAPER_LAMBDA, CC BY 4.0 articles).
Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "mg")

REPOSITORIES = {
    "Magnesium-FFs": "e71379d1a88b54221018af53ca54d3c18066ba05",
    "optimizedMgFFs": "1df53e9c94fc14a97864ac35c6951085bb60d1c0",
    "MgFF_OPC": "cf8148a3131c198b19f3a2896b68bad6f6f62f5b",
}
# water model (MDATK's id) -> (repository, parameter file); the files' notes name the water model
FILES = {
    "tip3p": ("Magnesium-FFs", "ff_Mg.itp"),
    "spce": ("optimizedMgFFs", "ff_Mg_spce.itp"),
    "opc": ("MgFF_OPC", "ff_Mg_opc.itp"),
    "fb3": ("optimizedMgFFs", "ff_Mg_tip3pfb.itp"),
    "tip4pew": ("optimizedMgFFs", "ff_Mg_tip4pew.itp"),
    "tip4pd": ("optimizedMgFFs", "ff_Mg_tip4pd.itp"),
    "tip4p2005": ("optimizedMgFFs", "ff_Mg_tip4p2005.itp"),
}
SHA256 = {
    "Magnesium-FFs/ff_Mg.itp": "08e7349831324e2c34e6917f5a7eecefaf8c301a7c4056c4218d7cbb05999ddb",
    "optimizedMgFFs/ff_Mg_spce.itp": "7188362ebd24fe902b0caf7692878ec4b6674ea45834e141384bcbeadfa18f27",
    "MgFF_OPC/ff_Mg_opc.itp": "8f84f1f759b8d1cb3458241f807f163d236dd8a5416181fac2d3afaafdc55d91",
    "optimizedMgFFs/ff_Mg_tip3pfb.itp": "52a7c973b6d5b96951cb5ce9ef78c30488e3e7a4371324d59479f81c70158667",
    "optimizedMgFFs/ff_Mg_tip4pew.itp": "1a88d2f74aed0b3889ef294cb2030fc4fbfd2c29c07a39c6def8bf5dfa17718d",
    "optimizedMgFFs/ff_Mg_tip4pd.itp": "71cde72bd6ff0db2e3a4b7d92f09b27653b8ec1e900f24792e2f2f4cbb54c7bf",
    "optimizedMgFFs/ff_Mg_tip4p2005.itp": "a8b29ee0d8ed35a073efe04c2b6716f4ea2b60f9749b5d0b08878eba27a4cd23",
}
MODELS = (("micro", "mMg"), ("nano", "nMg"))  # model -> the files' ion type
CHLORIDES = ("CXY", "Cl")  # the files' chloride types (CXY: optimised with the ion; Cl: Smith-Dang)
# Table 1 of the papers (lambda_sigma RNA, lambda_epsilon RNA), as printed: Grotz, Cruz-Leon &
# Schwierz, JCTC 17, 2530 (2021) (TIP3P); Grotz & Schwierz, JCTC 18, 526 (2022) (the others); the
# OPC paper (J. Chem. Phys. 156, 114501 (2022)) transfers the SPC/E set unchanged
PAPER_LAMBDA = {
    ("micro", "tip3p"): ("1.1375", "0.3200"), ("nano", "tip3p"): ("1.1435", "0.2500"),
    ("micro", "spce"): ("1.1019", "0.4856"), ("nano", "spce"): ("1.1107", "0.3300"),
    ("micro", "opc"): ("1.1019", "0.4856"), ("nano", "opc"): ("1.1107", "0.3300"),
    ("micro", "fb3"): ("1.0957", "0.4913"), ("nano", "fb3"): ("1.1002", "0.4172"),
    ("micro", "tip4p2005"): ("1.1484", "0.2648"), ("nano", "tip4p2005"): ("1.1345", "0.2923"),
    ("micro", "tip4pew"): ("1.1456", "0.2778"), ("nano", "tip4pew"): ("1.1489", "0.2371"),
    ("micro", "tip4pd"): ("1.3255", "0.3264"), ("nano", "tip4pd"): ("1.1231", "0.2834"),
}


def read_itp(path):
    """{section: [(line number, fields)]} of a GROMACS include file (comments and blanks dropped; the
    header before the first section is skipped: ff_Mg_opc.itp has a comment line there that lost its
    ';')."""
    sections, current = {}, None
    with open(path, encoding="utf-8") as fh:
        for number, line in enumerate(fh.read().split("\n"), start=1):
            text = line.split(";", 1)[0].strip()
            if not text:
                continue
            if text.startswith("[") and text.endswith("]"):
                current = text[1:-1].strip()
                sections.setdefault(current, [])
                continue
            if current is None:
                continue
            sections[current].append((number, text.split()))
    return sections


def printed_digits(text):
    """Decimal places of a number as printed ('0.25133119375' -> 11)."""
    mantissa = text.lower().split("e")[0]
    return len(mantissa.split(".")[1]) if "." in mantissa else 0


def extract(src):
    sets, pairs, notes = [], [], []
    for water, (repo, name) in FILES.items():
        rel = f"{repo}/{name}"
        path = os.path.join(src, repo, name)
        with open(path, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        if digest != SHA256[rel]:
            raise SystemExit(f"{rel}: SHA-256 {digest}, expected {SHA256[rel]} (another commit?)")
        sec = read_itp(path)
        atomtypes = {}
        for number, f in sec.get("atomtypes", []):
            if len(f) != 7 or f[0] in atomtypes:
                raise SystemExit(f"{rel}:{number}: unexpected [ atomtypes ] row")
            atomtypes[f[0]] = (f[5], f[6], number)
        nonbond = {}
        for number, f in sec.get("nonbond_params", []):
            if len(f) != 5 or f[2] != "1.000" or (f[0], f[1]) in nonbond:
                raise SystemExit(f"{rel}:{number}: unexpected [ nonbond_params ] row")
            nonbond[(f[0], f[1])] = (f[3], f[4], number)
        chloride = [c for c in CHLORIDES if c in atomtypes]
        if len(chloride) != 1:
            raise SystemExit(f"{rel}: chloride type not found")
        cl = chloride[0]
        rna_types = sorted(t for t in atomtypes if t not in (cl, "mMg", "nMg"))
        for model, ion in MODELS:
            isig, ieps, iline = atomtypes[ion]
            csig, ceps, cline = atomtypes[cl]
            msig, meps, mline = nonbond[(ion, cl)]
            sets.append((model, water, rel, ion, isig, ieps, str(iline), cl, csig, ceps, str(cline),
                         msig, meps, str(mline)))
            # the set's pairs: one lambda_sigma / lambda_epsilon for every RNA type
            lam_s, lam_e = set(), set()
            for t in rna_types:
                psig, peps, pline = nonbond[(ion, t)]
                pairs.append((model, water, rel, t, psig, peps, str(pline)))
                tsig, teps, _ = atomtypes[t]
                lb_s = (float(isig) + float(tsig)) / 2.0
                lb_e = math.sqrt(float(ieps) * float(teps))
                lam_s.add(round(float(psig) / lb_s, 4))
                if lb_e > 0.0:
                    lam_e.add(round(float(peps) / lb_e, 4))
                elif float(peps) != 0.0:
                    raise SystemExit(f"{rel}:{pline}: epsilon of a zero-epsilon type")
            if len(lam_s) != 1 or len(lam_e) != 1:
                raise SystemExit(f"{rel}: {model} pairs scale the combination rule unevenly "
                                 f"(lambda_sigma {sorted(lam_s)}, lambda_epsilon {sorted(lam_e)})")
            ls, le = lam_s.pop(), lam_e.pop()
            # every pair to the digits printed
            for t in rna_types:
                psig, peps, pline = nonbond[(ion, t)]
                tsig, teps, _ = atomtypes[t]
                want_s = ls * (float(isig) + float(tsig)) / 2.0
                want_e = le * math.sqrt(float(ieps) * float(teps))
                tol_s = 0.6 * 10.0 ** -printed_digits(psig)
                tol_e = 0.6 * 10.0 ** -printed_digits(peps)
                if abs(want_s - float(psig)) > tol_s or abs(want_e - float(peps)) > tol_e:
                    raise SystemExit(f"{rel}:{pline}: {t} pair is not lambda x Lorentz-Berthelot "
                                     f"({want_s!r}, {want_e!r})")
            paper = PAPER_LAMBDA[(model, water)]
            same = (f"{ls:.4f}", f"{le:.4f}") == paper
            notes.append(f"{model}Mg {water}: lambda_sigma {ls:.4f}, lambda_epsilon {le:.4f} "
                         f"(paper Table 1: {paper[0]}, {paper[1]}){'' if same else '  DIFFERS'}")
    return sets, pairs, notes


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dir", help="folder with the three repositories at the pinned commits")
    parser.add_argument("--out", default=OUT, help="folder for mg_sets.tsv and mg_rna_pairs.tsv")
    args = parser.parse_args(argv)
    sets, pairs, notes = extract(args.dir)
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "mg_sets.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("model\twater\tfile\tion_type\tion_sigma_nm\tion_epsilon_kj_per_mol\tion_line\t"
                 "chloride_type\tchloride_sigma_nm\tchloride_epsilon_kj_per_mol\tchloride_line\t"
                 "mg_cl_sigma_nm\tmg_cl_epsilon_kj_per_mol\tmg_cl_line\n")
        for row in sets:
            fh.write("\t".join(row) + "\n")
    with open(os.path.join(args.out, "mg_rna_pairs.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("model\twater\tfile\trna_type\tsigma_nm\tepsilon_kj_per_mol\tline\n")
        for row in pairs:
            fh.write("\t".join(row) + "\n")
    for note in notes:
        print(note)
    print(f"{len(sets)} sets, {len(pairs)} Mg2+-RNA pairs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
