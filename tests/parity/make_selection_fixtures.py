"""Topology fixtures of the selection and reader tests (tests/data/topology).

    python3 -I -B make_selection_fixtures.py MDATK_CHECKOUT OUT_DIR

writes fixture_chains.pdb, fixture_noelements.pdb, fixture_segids.pdb (60 atoms: two RNA chains with
insertion codes and old-style C4*, a peptide with digit-first 1HB, ions, water), fixture_chains_ba.pdb (chains B and A of the first, B first) and demo_dry.parm7, the
dry demo system of MDATK's tests (mdatk.demo.build_demo(n_frames=24, solvent=False); needs MDATK's
dependencies).
The PDB columns follow the wwPDB v3.3 ATOM/HETATM layout; only the element and segid columns
differ between the three variants.
"""
import os, shutil, sys, tempfile, warnings
warnings.filterwarnings("ignore")
if len(sys.argv) != 3:
    sys.exit("usage: make_selection_fixtures.py MDATK_CHECKOUT OUT_DIR")
MDATK = sys.argv[1]
DATA = os.path.abspath(sys.argv[2])

# (record, chain, resid, icode, resname, segid, [(name, element), ...])
RESIDUES = [
    ("ATOM", "A", 1, "", "G", "RNA1", [("P", "P"), ("OP1", "O"), ("C4'", "C"), ("HO2'", "H"), ("N9", "N"),
                                       ("N1", "N"), ("H1", "H"), ("O6", "O")]),
    ("ATOM", "A", 2, "", "C", "RNA1", [("P", "P"), ("C4'", "C"), ("N1", "N"), ("N3", "N"), ("N4", "N")]),
    ("ATOM", "A", 3, "", "A", "RNA1", [("P", "P"), ("C4'", "C"), ("N1", "N"), ("N6", "N")]),
    ("ATOM", "A", 3, "A", "A", "RNA1", [("P", "P"), ("C4'", "C")]),
    ("ATOM", "A", 4, "", "U", "RNA1", [("P", "P"), ("C4*", "C"), ("N3", "N"), ("H3", "H")]),
    ("ATOM", "B", 10, "", "G", "RNA2", [("P", "P"), ("C4'", "C"), ("N1", "N")]),
    ("ATOM", "B", 11, "", "C", "RNA2", [("P", "P"), ("C4'", "C"), ("N3", "N")]),
    ("ATOM", "B", 11, "A", "C", "RNA2", [("P", "P"), ("C4'", "C")]),
    ("ATOM", "B", 12, "", "U", "RNA2", [("P", "P"), ("C4'", "C")]),
    ("ATOM", "C", 5, "", "ALA", "PROT", [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"),
                                         ("H", "H"), ("1HB", "H")]),
    ("ATOM", "C", 6, "", "PHE", "PROT", [("N", "N"), ("CA", "C"), ("CG", "C"), ("CD1", "C"), ("CD2", "C"),
                                         ("CE1", "C"), ("CE2", "C"), ("CZ", "C")]),
    ("HETATM", "D", 201, "", "MG", "ION", [("MG", "MG")]),
    ("HETATM", "D", 202, "", "Mg", "ION", [("Mg", "MG")]),
    ("HETATM", "D", 203, "", "MG2", "ION", [("MG", "MG")]),
    ("HETATM", "D", 204, "", "Na+", "ION", [("Na+", "NA")]),
    ("HETATM", "D", 205, "", "CL", "ION", [("CL", "CL")]),
    ("HETATM", "D", 206, "", "K", "ION", [("K", "K")]),
    ("HETATM", "W", 301, "", "HOH", "WAT", [("O", "O"), ("H1", "H"), ("H2", "H")]),
    ("HETATM", "W", 302, "", "HOH", "WAT", [("O", "O"), ("H1", "H"), ("H2", "H")]),
]


def atom_name_field(name: str, element: str) -> str:
    if len(name) >= 4 or name[0].isdigit() or len(element) == 2:
        return f"{name:<4s}"
    return f" {name:<3s}"


def write_pdb(path: str, with_elements: bool, with_segids: bool) -> None:
    lines = ["REMARK   MDBridge selection fixture (see make_fixtures.py)\n"]
    k = 0
    for record, chain, resid, icode, resname, segid, atoms in RESIDUES:
        for name, element in atoms:
            k += 1
            x, y, z = 1.2 * k, 0.3 * (k % 5), 0.2 * (k % 7)
            el = element if with_elements else ""
            seg = segid if with_segids else ""
            lines.append(f"{record:<6s}{k:5d} {atom_name_field(name, element)} {resname:>3s} {chain:1s}{resid:4d}"
                         f"{icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}      {seg:<4s}{el:>2s}\n")
    lines.append("END\n")
    with open(path, "w") as fh:
        fh.write("".join(lines))


def main() -> None:
    os.makedirs(DATA, exist_ok=True)
    write_pdb(os.path.join(DATA, "fixture_chains.pdb"), with_elements=True, with_segids=False)
    write_pdb(os.path.join(DATA, "fixture_noelements.pdb"), with_elements=False, with_segids=False)
    write_pdb(os.path.join(DATA, "fixture_segids.pdb"), with_elements=True, with_segids=True)
    # the atoms of chain B written before those of chain A (chain numbering follows the names)
    lines = open(os.path.join(DATA, "fixture_chains.pdb")).read().splitlines(True)
    atoms = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    with open(os.path.join(DATA, "fixture_chains_ba.pdb"), "w") as fh:
        fh.write("".join([a for a in atoms if a[21] == "B"] + [a for a in atoms if a[21] == "A"]) + "END\n")
    sys.path.insert(0, MDATK)
    from mdatk.demo import build_demo
    tmp = tempfile.mkdtemp()
    paths = build_demo(tmp, n_frames=24, solvent=False)  # the tests' demo_dry fixture
    shutil.copy(paths["topology"], os.path.join(DATA, "demo_dry.parm7"))
    shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
