"""What MDATK's restraint-dictionary code (mdatk/dictionaries.py) gives, for MDBridge's chem::smiles
and chem::dictionaries (tests/chem/test_dictionaries.cpp; docs/specs/m6-plan.md 1.4, 4.2).

    python3 -I -B make_dict_cases.py MDATK_CHECKOUT OUT_DIR

Needs pytest (MDATK's tests/test_dictionaries.py builds the synthetic acedrg / eLBOW dictionaries:
synthetic_dictionary) and numpy (the structure builder places the atoms of the dictionaries that
carry coordinates). The dictionaries code is pure Python; its outputs are exact. The output is
reproducible (fixed seeds). Writes OUT_DIR/dict_cases.json:

  components    COMPONENTS: code, parent, base, terminal, name, smiles, atoms, bonds, hydrogens,
                hpos, leaving, formula, residue_atoms, chirality (the table W0-B's CCD extract must equal)
  smiles        parse_smiles of MDATK's test strings, the components' SMILES and a fuzzed set (fixed
                seed): elements, charges, hcount, adjacency, stereo, formula, or the error
  match         match_graph of permuted component graphs (and one with a bond removed): the mapping
  check         check_smiles(code, smiles) of every code x the 12 SMILES and MDATK's test SMILES:
                notes or the error; identify_smiles of each SMILES
  finish        finish_dictionary of synthetic dictionaries (every code, both styles, with and
                without coordinates, mirrored coordinates, flipped / missing chirality rows,
                protonation, two components in one file, CCD-style input, options group / name /
                source_id / allow_protonation / keep_names, wrong codes): text, notes, mapping,
                source_id, program, or the error
  OUT_DIR/dictionaries/0A5_acedrg.cif, 0U3_elbow.cif: two synthetic dictionaries, the inputs of
  cli_parity.py's `dict finish` cases

Errors are {"type": class name, "message": str(exc)}.
"""

import importlib.util
import json
import os
import random
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_dict_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
sys.path.insert(0, CHECKOUT)
sys.path.insert(0, os.path.join(CHECKOUT, "tests"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_common  # noqa: E402

chem_common.pin_kernels()
import improvements_chem  # noqa: E402

improvements_chem.apply_mdatk()
from mdatk import builder as B  # noqa: E402
from mdatk import dictionaries as D  # noqa: E402
from mdatk import pdbtools as pt  # noqa: E402

spec = importlib.util.spec_from_file_location("mdatk_test_dictionaries",
                                              os.path.join(CHECKOUT, "tests", "test_dictionaries.py"))
tests = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tests)
error = chem_common.error


def graph_dict(g):
    return {"elements": g.elements, "charges": g.charges, "hcount": g.hcount,
            "adj": [sorted([j, o] for j, o in a.items()) for a in g.adj],
            "stereo": sorted([k, list(v[0]), v[1]] for k, v in g.stereo.items()), "formula": g.formula()}


def component_cases():
    out = []
    for code, c in D.COMPONENTS.items():
        out.append({"code": code, "parent": c.parent, "base": c.base, "terminal": c.terminal, "name": c.name,
                    "smiles": c.smiles, "atoms": [list(a) for a in c.atoms], "bonds": [list(b) for b in c.bonds],
                    "hydrogens": [[x, list(h)] for x, h in c.hydrogens], "hpos": [list(h) for h in c.hpos],
                    "leaving": list(c.leaving), "description": c.description, "formula": c.formula(),
                    "residue_atoms": c.residue_atoms(),
                    "chirality": [[centre, list(order), sign] for centre, order, sign in c.chirality()]})
    return out


TEST_SMILES = ["c1ccncc1", "[NH4+].[Cl-]", "[H][C@](F)(Cl)Br", "C1CC%10CC1.C%10", "C1CC", "C(C", "C)", "[Xx", "",
               "Q", tests.D_ADENOSINE, tests.L_ADENOSINE, "Nc1ncnc2n(cnc12)C3OC(CO)C(O)C3O",
               "OC[C@@H]1O[C@H](N2C=[NH+]C3=C2N=C(N)NC3=O)[C@@H](O)[C@H]1O", "C#N", "C$C", "O=C=O", "[13CH4]",
               "[Fe+2]", "[O-][N+](=O)C", "C/C=C\\C", "CC(=O)O.[Na]", "[C@@H](F)(Cl)Br", "[nH]1cccc1",
               "c1cccc1", "C%1CC%1", "C1CC1C2CC2", "[CH3:1]C", "*C", "BrCCl", "S(=O)(=O)O", "P(=O)(O)(O)O",
               "[2H]C", "C((C))", "C=1CC1", "C1CC=1", "[NH3+][C@@H](C)C(=O)[O-]", "c1ccc2ccccc2c1", "C.C",
               "[se]1cccc1", "B(O)(O)O", "I", "[I-]", "[H+]", "[H][H]", "C[C@@TH1](F)Cl", "[C@TH2H](F)(Cl)Br"]


def fuzzed(n, seed):
    rng = random.Random(seed)
    alphabet = "CNOPcno()[]=#123@+-H%0."
    sources = [c.smiles for c in D.COMPONENTS.values()] + TEST_SMILES[:4]
    out = []
    for _ in range(n):
        s = list(rng.choice(sources))
        for _ in range(rng.randint(1, 3)):
            pos = rng.randrange(len(s) + 1)
            op = rng.random()
            if op < 0.4 and pos < len(s):
                s[pos] = rng.choice(alphabet)
            elif op < 0.7:
                s.insert(pos, rng.choice(alphabet))
            elif pos < len(s):
                del s[pos]
        out.append("".join(s))
    return out


def smiles_cases():
    out = []
    for s in TEST_SMILES + [c.smiles for c in D.COMPONENTS.values()] + fuzzed(150, 9):
        try:
            out.append({"smiles": s, "graph": graph_dict(D.parse_smiles(s))})
        except D.DictionaryError as exc:
            out.append({"smiles": s, "error": error(exc)})
    return out


def match_cases():
    out = []
    rng = random.Random(3)
    for code in ("0G3", "0A5", "0C", "0U3"):
        ref = D.COMPONENTS[code].graph()
        for k in range(3):
            perm = list(range(len(ref)))
            rng.shuffle(perm)
            tgt = D.MolGraph()
            inverse = {old: new for new, old in enumerate(perm)}
            for old in perm:
                tgt.add_atom(ref.elements[old], ref.charges[old], ref.hcount[old])
            for i in range(len(ref)):
                for j, order in ref.adj[i].items():
                    if i < j:
                        tgt.add_bond(inverse[i], inverse[j], order)
            if k == 2:  # a bond removed: no match
                a, b = next((i, j) for i in range(len(tgt)) for j in tgt.adj[i])
                del tgt.adj[a][b], tgt.adj[b][a]
            out.append({"code": code, "perm": perm, "removed": k == 2, "match": D.match_graph(ref, tgt)})
    return out


def check_cases():
    smiles = [c.smiles for c in D.COMPONENTS.values()] + [tests.D_ADENOSINE, tests.L_ADENOSINE,
                                                         "Nc1ncnc2n(cnc12)C3OC(CO)C(O)C3O", "CCO", "C1CC"]
    out = {"check": [], "identify": []}
    for s in smiles:
        for code in list(D.COMPONENTS) + ["0a5", "XY1"]:
            try:
                out["check"].append([code, s, D.check_smiles(code, s)])
            except D.DictionaryError as exc:
                out["check"].append([code, s, error(exc)])
        try:
            out["identify"].append([s, D.identify_smiles(s)])
        except D.DictionaryError as exc:
            out["identify"].append([s, error(exc)])
    return out


# --- finish ------------------------------------------------------------------------------------------

def coordinates(comp, mirror=False):
    """Positions for the component's atoms from an L helix of the structure builder (CCD names; atoms the
    builder does not place are put 1 Å from their bonded neighbour along a fixed direction)."""
    text = B.build_helix(comp.base * 3, partner="none", chirality="L", scheme="ccd", hydrogens=True, phosphate5=True)
    res = list(pt.read_pdb(text).residues())[1]
    pos = {a.name: a.xyz for a in res.atoms}
    names = [a[0] for a in comp.atoms] + [h for _x, hs in comp.hydrogens for h in hs]
    parent = {h: x for x, hs in comp.hydrogens for h in hs}
    for a, b, _o in comp.bonds:
        parent.setdefault(b, a)
        parent.setdefault(a, b)
    out = {}
    for k, n in enumerate(names):
        if n in pos:
            out[n] = pos[n]
    for k, n in enumerate(names):
        if n not in out:
            base = out.get(parent.get(n), (0.0, 0.0, 0.0))
            out[n] = (base[0] + 0.6, base[1] + 0.5 + 0.01 * k, base[2] - 0.4)
    if mirror:
        out = {n: (-p[0], p[1], p[2]) for n, p in out.items()}
    return {n: tuple(round(v, 4) for v in p) for n, p in out.items()}


def finish_case(label, text, code, **kw):
    case = {"label": label, "input": text, "code": code, "options": kw}
    try:
        res = D.finish_dictionary(text, code, **kw)
        case.update(text=res.text, notes=res.notes, mapping=res.mapping, source_id=res.source_id, program=res.program)
    except D.DictionaryError as exc:
        case["error"] = error(exc)
    return case


def finish_cases():
    out = []
    syn = tests.synthetic_dictionary
    for code, comp in D.COMPONENTS.items():
        for style in ("acedrg", "elbow"):
            text, _n, _v = syn(comp, style, seed=len(code) + len(style))
            out.append(finish_case(f"{code} {style}", text, code))
        for mirror in (False, True):
            for seed in (1, 2):
                text, _n, _v = syn(comp, "acedrg", seed=seed, coords=coordinates(comp, mirror))
                out.append(finish_case(f"{code} coords mirror={mirror} seed={seed}", text, code))
        text, _n, _v = syn(comp, "acedrg", coords=coordinates(comp, True), chir=False)
        out.append(finish_case(f"{code} mirrored coords without chirality rows", text, code))
    comp = D.COMPONENTS["0U5"]
    text, _n, _v = syn(comp, "acedrg")
    out.append(finish_case("group and name", text, "0U5", group="non-polymer", name="my uridine"))
    other, _n, _v = syn(D.COMPONENTS["0G5"], "acedrg")
    both = text + other.replace("data_comp_list", "data_comp_list_2")
    out += [finish_case("two components", both, "0U5"), finish_case("source id", both, "0U5", source_id="0U5"),
            finish_case("missing source id", both, "0U5", source_id="ZZZ"), finish_case("no reference", text, "XY1"),
            finish_case("bad code", text, "bad code"), finish_case("keep names", text, "0U5", keep_names=True),
            finish_case("lower-case code", text, "0u5"), finish_case("group RNA", text, "0U5", group="RNA")]
    a5, _n, _v = syn(D.COMPONENTS["0A5"], "acedrg")
    out.append(finish_case("wrong compound", a5, "0G5"))
    flipped, _n, _v = syn(D.COMPONENTS["0A5"], "elbow", flip_chir=True)
    out.append(finish_case("flipped chirality", flipped, "0A5"))
    out.append(finish_case("not a dictionary", "data_x\n_a.b 1\n", "0A5"))
    nochir, _n, _v = syn(D.COMPONENTS["0A5"], "acedrg", chir=False)
    out.append(finish_case("no chirality rows", nochir, "0A5"))
    prot, _n, _v = syn(D.COMPONENTS["0C3"], "elbow", extra_h=[("HX1", "OP2"), ("HX2", "OP3")])
    out.append(finish_case("phosphate protonation", prot, "0C3"))
    n7, _n, _v = syn(D.COMPONENTS["0G5"], "acedrg", extra_h=[("HX1", "N7")])
    out += [finish_case("base protonation", n7, "0G5"),
            finish_case("base protonation allowed", n7, "0G5", allow_protonation=True)]
    g, _n, _v = syn(D.COMPONENTS["0G"], "acedrg")
    first = D.finish_dictionary(g, "0G")
    out += [finish_case("CCD-style input", first.text, "0G3"), finish_case("finished again", first.text, "0G")]
    prefixed = "# monomer library\n" + g.replace("data_comp_0G", "data_comp_0g")
    out.append(finish_case("monomer-library header", prefixed, "0g3"))
    return out


def main():
    data = {"components": component_cases(), "smiles": smiles_cases(), "match": match_cases(),
            "check": check_cases(), "finish": finish_cases()}
    os.makedirs(os.path.join(OUT, "dictionaries"), exist_ok=True)
    for code, style in (("0A5", "acedrg"), ("0U3", "elbow")):  # inputs of cli_parity.py's dict finish cases
        text, _n, _v = tests.synthetic_dictionary(D.COMPONENTS[code], style, seed=len(code) + len(style))
        with open(os.path.join(OUT, "dictionaries", f"{code}_{style}.cif"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    with open(os.path.join(OUT, "dict_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
