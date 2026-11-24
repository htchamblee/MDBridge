"""What MDAnalysis' select_atoms gives for selection strings (tests/data/selections).

    /opt/mdbridge-venv/bin/python -I -B make_atom_selection_cases.py MDATK_CHECKOUT OUT_DIR

MDATK passes its selections to MDAnalysis (2.11.0), so MDAnalysis is the reference of
src/core/atom_selection.cpp. Writes

  atom_selection_cases.json   per fixture of tests/data/topology (the topology and first frame the
                              C++ tests rebuild from the *.topology.json dumps): thousands of
                              selection strings with the selected atoms or the error text
  atom_selection_distance_cases.json
                              frames with random positions (no box, orthorhombic, triclinic, atoms
                              far outside the cell) and the distance keywords evaluated on them

Selected atoms are stored as [first, last] runs of indices. "error" is str(exception) (what MDATK
puts in its messages); "unsupported" marks selections MDBridge rejects on purpose (smarts, cyzone,
cylayer with coordinates, type_index), where MDAnalysis gives a result or an ImportError.
MDATK_CHECKOUT is accepted for symmetry with the other generators; MDATK itself is not imported.
"""

import json
import os
import random
import sys
import warnings

import MDAnalysis as mda
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)

TOPOLOGY = os.path.join(HERE, "..", "data", "topology")
FIXTURES = ["fixture_chains.pdb", "fixture_noelements.pdb", "fixture_segids.pdb",
            "fixture_chains_ba.pdb", "demo_dry.parm7"]


def runs(indices):
    out = []
    for i in indices:
        i = int(i)
        if out and out[-1][1] == i - 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return out


def unsupported(sel, u):
    words = sel.replace("(", " ").replace(")", " ").split()
    if "smarts" in words or "type_index" in sel:
        return True
    connectivity = ("bonds", "angles", "dihedrals", "impropers", "ureybradleys", "cmaps")
    if any(w in sel for w in connectivity) and hasattr(u, "bonds"):  # prop / same on bond objects
        return True
    has_coords = hasattr(u, "trajectory") and u.trajectory is not None and len(u.trajectory) > 0
    return has_coords and ("cyzone" in words or "cylayer" in words)


def evaluate(u, sel):
    case = {"sel": sel}
    try:
        case["atoms"] = runs(u.select_atoms(sel).indices)
    except Exception as exc:  # what MDATK's select() wraps
        case["error"] = str(exc)
    if unsupported(sel, u):
        case["unsupported"] = True
    return case


GLOBS = ["*", "?", "??", "C*", "*'", "C4\\'", "\\*", "N?", "[CN]*", "[!C]*", "[A-C]?", "N[0-9]", "[]]*",
         "[!]", "[!]]*", "[z-a]", "[!z-a]", "[a-]*", "[-a]*", "[C-N-P]*", "[^C]*", "[[]*", "[\\]*", "H*",
         "*1", "O?'", "[c-a-z]*", "[A-Z][A-Z]", "*[0-9]*", "[CN]?'", "C4[']", "c*", "M*", "[MN][Ag]*",
         "*+", "[!A-M]*", "*G", "[P]", "[&&]*", "[~|]*", "[!-]*", "[a-c-e]", "[C-]", "[--C]*", "*?*",
         "**", "[*]", "[?]*", "P[", "[", "]", "[!", "é*", "1HB", "1*", "MG", "Mg", "Na+", "HOH"]

INT_VALUES = ["0", "1", "3", "-1", "1:3", "1-3", "1 to 3", "1 : 3", "1 :3", "1: 3", "1 - 3", "3-1",
              "-3--1", "-3:-1", "1-3abc", "١٢", "1_0", "+2", "007", "x", "1.5", "1e3",
              "99999999999999999999", "-99999999999999999999:5", "2:99999999999999999999",
              "1 2 5", "1 3-5 7", "3 to", "to 3", ": 3", "1 -", "--5", "1 --5", "1 to5", "5:", "1::3",
              "1:-1", "0- -1", "1 to 3 5 : 7", "1-", "-", ":", "to", "1 to", "2 3 to 4 -6",
              "0x10", "1__0", "_1", "٣:٥", "4 5 to 2", "4 tomorrow", "1 to - 3", "200:206",
              "10:12", "773", "1:773", "772", "24:26", "1to3", "5-7 9:11 13 to 15"]

RESID_VALUES = ["3A", "3A:4", "1:3A", "3A-3B", "11A", "11:11A", "3a", "3 A", "4:3", "3A:3", "1-3xyz",
                "3A to 4", "10-12", "200:206", "-1:3", "11A:12", "3:3", "3A:3A", "3B:3A", "2:11",
                "٣A", "3_", "3é", "1A-99999999999999999999", "11 11A", "10 to 11A"]

FLOAT_VALUES = ["1", "1.0", "1.008", "12.011", "0", "0.0", "1:13", "1-13", "1 to 13", "-1:1",
                "-inf:inf", "inf", "-inf", "nan", "1e3", "1e-9", ".5:2", "1e1:2e1", "-0.5", "abc",
                "1.0:abc", "14.007 15.999", "12.0110001", "12.0111", "-0.4157", "0.5973", ":3", "1-",
                "1e999", "1_0", "١.5", "-.5:.5", "1.:2.", "2-1", "0.0 1.0", "1E3", "1e+1:2e-1",
                "-1--0.1", "1 - 13", "24.305", "22.99", "0:1e999", "-1e999:0", ". : .", "0.8:", "-",
                "1.5-", "e:1", "1e:3", "0.00:1.00", "1.000000001", "0.99999"]

PROP_ATTRS = ["x", "y", "z", "mass", "charge", "occupancy", "tempfactor", "index", "resindex",
              "segindex", "id", "resid", "resnum", "formalcharge", "name", "type", "element",
              "resname", "segid", "chainID", "icode", "altLoc", "record_type", "radius", "bonds",
              "topologyattr", "moltype", "type_index", "bynum", "q", "None", "angles", "cmaps"]

SAME_PROPS = ["fragment", "x", "y", "z", "residue", "segment", "name", "type", "resname", "resid",
              "segid", "mass", "charge", "radius", "bfactor", "resnum", "index", "resindex",
              "segindex", "id", "element", "record_type", "chainID", "tempfactor", "formalcharge",
              "occupancy", "altLoc", "icode", "moltype", "molnum", "bonds", "topologyattr",
              "type_index", "chirality", "foo", "None", "Name", "as"]

HAND = r"""
all
not all
all and not all
name P
name P C4'
name P and resid 1
not name P and resid 1
resname G and name P or name N1
name P or name N1 and resname G
(name P or name N1) and resname G
name P or (name N1 and resname G)
not (name P or name N1)
not not name P
not name P or name N1
byres name P and resid 1
byres (name P) and resid 1
byres name P or name N1
global name P
global name P and resid 1
not global name P
name P and not global resid 1
byres byres name P
byres not name P
(((name P)))
( ( name P ) )
(name P)and(resid 1)
name P)and(resid 1
name P )
( name P
(name P
name P and
and
or
and name P
name P or or name N1
not
(
)
()
( )
name
resid
mass
bynum
index
resname
segid
name and resid 1
name ( P )
name P (
not and name P
name \and
name \not P
name \\P
name \(
resname \protein
protein
nucleic
water
backbone
nucleicbackbone
nucleicbase
nucleicsugar
protein and backbone
nucleic and not nucleicbackbone
nucleicsugar or nucleicbase
protein A
backbone and name CA
water and name O
not water
nucleic or protein or water
NAME P
Name P
AND
name P AND resid 1
name P And resid 1
chainID A
chainid A
chainID A B
altLoc A
altloc A
altLoc ""
record_type ATOM
record_type HETATM
icode A
icode ""
moltype x
moltype *
segid RNA1
segid ""
segid *
segid RNA*
element C
element MG
element Mg
type C
type MG
type *
same residue as name P
same segment as index 3
same resname as name P
same name as resid 1
same type as index 0
same mass as name P
same charge as name P
same resid as index 18
same resnum as index 18
same chainID as index 30
same element as index 3
same x as name P
same y as index 3
same z as index 3
same fragment as name P
same fragment as index 0
same residue as not all
same name as not all
same chainID as not all
same
same residue
same residue name P
same as name P
same residue as
same residue as and
same foo as name P
same None as name P
same Name as name P
bonded name P
bonded not all
bonded all
bonded name P and name O*
bonded bonded name P
name H* and bonded name O*
around 3 name P
around 3 not all
around abc name P
around
around 3
sphzone 3 name P
sphzone 3 not all
sphlayer 1 3 name P
sphlayer 3 1 name P
isolayer 1 3 name P
point 1 1 1 3
point 1 1 1
point 1
point
point a b c d
cyzone 1 1 1 name P
cyzone 1 1 1 not all
cylayer 1 2 1 1 name P
cylayer 1 2 1 1 not all
cyzone 1 1 name P
smarts C
smarts C(C)C
smarts [#6] and name P
smarts ( C ) )
smarts
aromatic
aromaticity
aromaticity true
aromaticity False
aromaticity yes
aromatic true
chirality R
chirality
chirality R S
group a
group updating
group
atom A 1 P
atom A 1 P and name P
atom A x P
atom A
atom A 1
atom RNA1 1 P
atom "" 1 P
atom  1 P
atom SYSTEM 1 P
atom A 99999999999999999999 P
atom A 1 Q
prop x > 1
prop x
prop
prop x <
prop x << 3
prop x <> 3
prop x = 3
prop 3 = x
prop 3 > x
prop 3 <= x
prop 3 == x
prop 3 != x
prop abs x < 3
prop abs
prop abs x
prop abs 3 > x
prop x>=1.2
prop x>= 1.2
prop x >=1.2
prop 1.2<=x
prop 1.2 <= x
prop x <= abc
prop x == 1.2
prop x != 1.2
prop x == 1.2000001
prop q < 3
prop x a<5
prop x <5 6
prop y >= 0.3
prop z < 0
prop abs z <= 0.6
prop x < 1e999
prop x > -1e999
prop x == 1e300
prop x != nan
prop x < nan
prop mass >= 12
prop mass == 12.011
prop charge > 0
prop index >= 3.5
prop index == 3
prop resid > 3
prop abs resid > 3
prop name > 1
prop name == 1
prop name != 1
prop abs name > 1
prop chainID > 1
prop bonds > 1
prop radius > 1
prop topologyattr > 1
prop icode > 1
prop None > 1
prop x>None
prop x > y
prop x == x
resid 3A
resid 1A-3
resid 3 :4
resid 3- 4
resid 3 to4
resid 3to 4
resid 3--1
resid 1:3A
index 1-3abc
index 3--1
index -1--3
mass 1--3
mass -inf:inf
mass :3
mass 1 to 3
index 3 :
index : 3
index to
resnum 1e3
index 1_0
mass 1e999
mass nan
name P	and	resid 1
name P
and resid 1
name P and resid 1
name P and resid 1
name P and resid 1
name P and resid 1
name P and resid 1
name P and resid ١
resid ٣
name P and resid 1
"""

SPECIAL = ["", " ", "\t", "\n", "　", "name P\u0085and resid 1", "name P\x1cand resid 1",
           "name P​and resid 1", "name P", "name P\x00", "name é", "resname G or resname C"]


def string_cases(u, rng):
    sels = []
    values = set()
    for attr in ("names", "types", "resnames", "segids"):
        try:
            values.update(str(v) for v in np.unique(getattr(u.atoms, attr)))
        except Exception:
            pass
    values = sorted(values)
    if len(values) > 40:
        values = rng.sample(values, 40)
    keywords = ["name", "type", "element", "altLoc", "altloc", "record_type", "chainID", "chainid",
                "resname", "icode", "moltype", "segid"]
    for kw in keywords:
        for v in values + GLOBS:
            sels.append(f"{kw} {v}")
        sels.append(f"{kw} " + " ".join(rng.sample(values, min(3, len(values)))))
        sels.append(f"{kw} {rng.choice(GLOBS)} {rng.choice(values)}")
    return sels


def numeric_cases():
    sels = []
    for kw in ("index", "bynum", "resindex", "segindex", "id", "resnum", "resid", "formalcharge", "model",
               "molnum", "nbindex"):
        for v in INT_VALUES:
            sels.append(f"{kw} {v}")
    for v in RESID_VALUES:
        sels.append(f"resid {v}")
        sels.append(f"resnum {v}")
    for kw in ("mass", "charge", "occupancy", "tempfactor", "bfactor", "radius", "epsilon", "epsilon14",
               "gbscreen", "rmin", "rmin14", "solventradius"):
        for v in FLOAT_VALUES:
            sels.append(f"{kw} {v}")
    return sels


def prop_cases(rng):
    sels = []
    ops = [">", "<", ">=", "<=", "==", "!="]
    for attr in PROP_ATTRS:
        for op in ops:
            v = rng.choice(["0", "1.2", "3", "-1", "12.011", "5.5", "0.5973"])
            sels.append(f"prop {attr} {op} {v}")
            sels.append(f"prop abs {attr} {op} -{v}")
            sels.append(f"prop {v} {op} {attr}")
            sels.append(f"prop {attr}{op}{v}")
    return sels


def same_cases():
    sels = []
    for p in SAME_PROPS:
        sels.append(f"same {p} as name P")
        sels.append(f"same {p} as index 3 18 40")
        sels.append(f"same {p} as not all")
    return sels


def random_cases(u, rng, count):
    atoms = ["name P", "name C4'", "name N*", "resname G", "resname A C", "resid 1:3", "resid 3A",
             "resid 11", "index 0:9", "bynum 20-30", "mass 12.011", "mass 1:13", "segid RNA1",
             "chainID A", "chainID B C", "element C", "type O*", "protein", "nucleic", "water",
             "backbone", "nucleicbackbone", "resnum 5", "resindex 3:5", "segindex 0", "id 5:15",
             "all", "name H*", "charge -1:0", "occupancy 1", "tempfactor 0", "record_type HETATM",
             "icode A", "altLoc A", "same residue as name N1", "bonded name P", "same resname as index 4",
             "prop mass > 10", "prop x < 20", "around 2 name P", "sphzone 4 resid 1", "point 10 1 1 5"]
    def expr(depth):
        r = rng.random()
        if depth <= 0 or r < 0.3:
            return rng.choice(atoms)
        if r < 0.45:
            return "not " + expr(depth - 1)
        if r < 0.55:
            return "(" + expr(depth - 1) + ")"
        if r < 0.6:
            return rng.choice(["byres ", "global ", "same segid as ", "bonded "]) + expr(depth - 1)
        return expr(depth - 1) + rng.choice([" and ", " or ", " and not "]) + expr(depth - 1)
    return [expr(4) for _ in range(count)]


FUZZ_TOKENS = ["name", "resname", "resid", "index", "bynum", "mass", "charge", "prop", "same", "as",
               "and", "or", "not", "global", "byres", "bonded", "around", "sphzone", "sphlayer", "point",
               "atom", "group", "protein", "nucleic", "all", "(", ")", "P", "G", "C*", "1", "3A", "1:3",
               "-", ":", "to", "2.5", "x", "<", ">=", "==", "residue", "fragment", "segid", "RNA1", "A",
               "\\and", "chainID", "icode", "type", "element", "abs", "x<3", "3>x", "None", "1-3", "-1",
               "smarts", "cyzone", "aromaticity", "true", "chirality", "isolayer", "updating", "id",
               "resnum", "tempfactor", "bfactor", "occupancy", "record_type", "altLoc", "moltype"]


def fuzz_cases(rng, count):
    out = []
    for _ in range(count):
        length = rng.randint(1, 9)
        out.append(" ".join(rng.choice(FUZZ_TOKENS) for _ in range(length)))
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: make_atom_selection_cases.py MDATK_CHECKOUT OUT_DIR")
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    hand = [line for line in HAND.strip("\n").split("\n")] + SPECIAL
    out = {"mdanalysis": mda.__version__, "fixtures": {}}
    for fixture in FIXTURES:
        rng = random.Random(fixture)
        u = mda.Universe(os.path.join(TOPOLOGY, fixture))
        sels = hand + string_cases(u, rng) + numeric_cases() + prop_cases(rng) + same_cases()
        sels += random_cases(u, rng, 600)
        sels += fuzz_cases(rng, 2500)
        seen = set()
        cases = []
        for sel in sels:
            if sel not in seen:
                seen.add(sel)
                cases.append(evaluate(u, sel))
        out["fixtures"][fixture] = cases
    with open(os.path.join(out_dir, "atom_selection_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")
    write_distance_cases(out_dir)


DISTANCE = ["around {r} name P", "around {r} resid 1", "around {r} index {i}", "around {r} not name P",
            "around {r} name C*", "around 0 name P", "around -1 name P", "around {r} all",
            "around {r} not all", "sphzone {r} name P", "sphzone {r} index {i}", "sphzone {r} resid 1:3",
            "sphzone 0 index {i}", "sphzone {r} not all", "sphlayer {r2} {r} name P", "sphlayer 0 {r} index {i}",
            "sphlayer {r} {r2} name P", "isolayer {r2} {r} name P", "isolayer {r2} {r} resid 1",
            "point {x} {y} {z} {r}", "point {x} {y} {z} 0", "point -100 -100 -100 {r}",
            "point {x} {y} {z} 1e9", "byres around {r} index {i}", "around {r} name P and name N*",
            "not around {r} name P", "around {r} around {r2} index {i}", "around {r} (name P or name N1)",
            "same x as index {i}", "same y as name P", "same z as resid 1", "prop x < {x}", "prop y >= {y}",
            "prop abs z <= {z}", "prop x == {px}", "prop {x} > x", "prop z != {pz}",
            "cyzone {r} 2 -2 name P", "cylayer 1 {r} 2 -2 name P", "around {r} index {i} {j} {k}",
            "sphzone {r} index {i} {j} {k}", "isolayer 1 {r} index {i} {j}", "around {big} name P",
            "around {r} prop x < {x}", "global around {r} name P"]


def write_distance_cases(out_dir):
    rng = np.random.default_rng(2026)
    pyrng = random.Random(7)
    boxes = {
        "nobox": None,
        "ortho": [18.0, 21.0, 24.0, 90.0, 90.0, 90.0],
        "ortho_small": [9.0, 10.0, 12.0, 90.0, 90.0, 90.0],
        "triclinic": [20.0, 22.0, 25.0, 70.0, 80.0, 100.0],
        "triclinic_skew": [15.0, 16.0, 18.0, 60.0, 75.0, 65.0],
        "zero": [0.0, 0.0, 0.0, 90.0, 90.0, 90.0],
        "invalid_angles": [10.0, 10.0, 10.0, 10.0, 10.0, 170.0],
    }
    out = {"mdanalysis": mda.__version__, "frames": []}
    for fixture in ("fixture_chains.pdb", "demo_dry.parm7"):
        top = mda.Universe(os.path.join(TOPOLOGY, fixture))
        n = top.atoms.n_atoms
        for name, box in boxes.items():
            for spread in ("inside", "outside"):
                if fixture == "demo_dry.parm7" and spread == "outside" and name not in ("ortho", "triclinic"):
                    continue
                lo, hi = (0.0, 20.0) if spread == "inside" else (-45.0, 70.0)
                pos = rng.uniform(lo, hi, size=(n, 3)).astype(np.float32)
                u = mda.Universe(os.path.join(TOPOLOGY, fixture))
                u.load_new(pos[None, :, :], format=MemoryReader,
                           dimensions=None if box is None else np.array(box, dtype=np.float32))
                if box is None:
                    u.trajectory.ts.dimensions = None
                cases = []
                for template in DISTANCE:
                    i, j, k = (pyrng.randrange(n) for _ in range(3))
                    x, y, z = (round(pyrng.uniform(lo, hi), 3) for _ in range(3))
                    sel = template.format(r=round(pyrng.uniform(1.5, 6.0), 2), r2=round(pyrng.uniform(0.5, 1.5), 2),
                                          i=i, j=j, k=k, x=x, y=y, z=z, big=1000,
                                          px=float(pos[i, 0]), pz=float(pos[k, 2]))
                    cases.append(evaluate(u, sel))
                dims = u.trajectory.ts.dimensions
                out["frames"].append({
                    "fixture": fixture,
                    "label": f"{name}/{spread}",
                    "positions": [[float(v) for v in row] for row in pos],
                    "dimensions": None if box is None else [float(np.float32(v)) for v in box],
                    "mdanalysis_dimensions": None if dims is None else [float(v) for v in dims],
                    "cases": cases,
                })
    with open(os.path.join(out_dir, "atom_selection_distance_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
