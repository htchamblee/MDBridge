"""What MDATK's loader (mdatk/io/loader.py) and MDAnalysis 2.11.0 make of topology files and file
combinations (tests/data/io/loader).

    /opt/mdbridge-venv/bin/python -I -B make_loader_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDAnalysis 2.11.0 (and MDATK's dependencies); nothing is written into the checkout. Writes

  files/*.pdb, *.gro    crafted PDB / GRO files: elements (valid, invalid, partial, none), segids
                        and chain IDs, formal charges, resid and serial wraps, hybrid-36 serials,
                        CONECT variants, models with and without CRYST1, placeholder / broken
                        cells, HEADER / COMPND / REMARK, line ends, Unicode, GRO boxes and wraps,
                        ions.pdb / ions.gro (protein atoms and ions for the residue-scoped rule)
  files/*.topology.json the topology MDAnalysis builds from each readable file (the dump of
                        tests/parity/dump_topologies.py)
  topology_cases.json   per file: the dump name or the Universe's error, and the PDB / GRO reader:
                        n_frames, dt, header / compound / remarks and every frame (positions as
                        float32 bit patterns, dimensions, time) read in order and in a shuffled
                        order, with dt given and not
  loader_cases.json     load_universe(topology, trajectories, dt_ps) + system_summary over a
                        matrix of the demo system and crafted files, or the error; frame_indices;
                        the JSON `mdatk info --top ... --traj ...` prints

Paths in messages are written with the test data folder as "{data}" and OUT_DIR as "{out}".
Errors are {"type": class name, "message": str(exc)}.
"""

import contextlib
import io
import json
import os
import re
import random
import struct
import sys
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_loader_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402

from mdatk import cli  # noqa: E402
from mdatk.io import loader  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
sys.path.insert(0, HERE)
import dump_topologies  # noqa: E402
import improvements_cli  # noqa: E402
import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)

OUT = os.path.abspath(sys.argv[2])
FILES = os.path.join(OUT, "files")


def placeholders(text):
    return text.replace(FILES, "{out}/files").replace(OUT, "{out}").replace(DATA, "{data}")


def error(exc):
    # MDBridge's wording where it improves MDATK's (stop=end: improvements_cli.py, CLI-03)
    return {"type": type(exc).__name__, "message": improvements_cli.improve_text(placeholders(str(exc)))}


def f32bits(values):
    return "".join(struct.pack(">f", float(v)).hex() for v in np.asarray(values, dtype=np.float32).ravel())


# --- crafted files --------------------------------------------------------------------------------

# ions.pdb / ions.gro (no element column: every type is guessed): a protein fragment whose CA, CD1,
# HG1, CE atoms must stay C / H, lone ions named the CCD way (residue = atom name) and the force-field
# way, CHARMM's cadmium (CD2/CD), a CCD mercury component (EMC/HG) and K+ (in no source: name rule)
# for MDBridge's residue-scoped ion rule (independent_tables.py).
IONS = [("ILE", 1, ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "HD11"]),
        ("THR", 2, ["N", "CA", "CB", "OG1", "HG1"]),
        ("MET", 3, ["CA", "SD", "CE"]),
        ("CA", 4, ["CA"]), ("CD", 5, ["CD"]), ("HG", 6, ["HG"]), ("ZN", 7, ["ZN"]), ("MG", 8, ["MG"]),
        ("NA", 9, ["NA"]), ("CL", 10, ["CL"]), ("SOD", 11, ["SOD"]), ("CAL", 12, ["CAL"]),
        ("CLA", 13, ["CLA"]), ("POT", 14, ["POT"]), ("K+", 15, ["K+"]), ("CD2", 16, ["CD"]),
        ("EMC", 17, ["HG"]), ("ZN2", 18, ["ZN"])]


def ion_atoms():
    out = []
    for resname, resid, names in IONS:
        for name in names:
            out.append((resname, resid, name, len(out)))
    return out

def atom(serial, name, resname, chain, resid, xyz, record="ATOM", alt=" ", icode=" ", occ="  1.00",
         temp="  0.00", segid="    ", element="  ", charge="  "):
    """A PDB atom record from raw column texts (serial / resid as text when given as str)."""
    serial = serial if isinstance(serial, str) else f"{serial:5d}"
    resid = resid if isinstance(resid, str) else f"{resid:4d}"
    name = name if len(name) >= 4 else " " + name.ljust(3)
    x, y, z = xyz
    return (f"{record:<6}{serial:>5} {name[:4]}{alt}{resname:>3} {chain}{resid}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:>6}{temp:>6}      {segid:<4}{element:>2}{charge:<2}")


def chain_atoms(n, start=1, chain="A", resname="G", per_residue=3, element=True, segid="    ", offset=0.0):
    names = ["P", "C4'", "N1", "O6", "H1"]
    out = []
    for k in range(n):
        name = names[k % per_residue]
        out.append(atom(start + k, name, resname, chain, 1 + k // per_residue,
                        (1.5 * k + offset, 0.25 * k, -0.5 * k), segid=segid,
                        element=name[0] if element else "  "))
    return out


def pdb_files():
    files = {}
    seven = chain_atoms(7)
    files["seven.pdb"] = "\n".join(["HEADER    RNA FRAGMENT                            24-NOV-26   1ABC",
                                    "TITLE     SEVEN ATOMS", "COMPND    MOL_ID: 1;", "COMPND   2 MOLECULE: RNA;",
                                    "REMARK   1 crafted for MDBridge",
                                    "REMARK   2 " + "a long remark " * 8,
                                    "CRYST1   30.000   40.000   50.000  90.00  90.00 120.00 P 1           1"]
                                   + seven + ["TER", "CONECT    1    2", "CONECT    2    1    3", "END"]) + "\n"
    files["seven_noelements.pdb"] = "\n".join(chain_atoms(7, element=False)) + "\nEND\n"
    files["seven_crlf.pdb"] = "\r\n".join(["CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1"]
                                          + seven) + "\r\n"
    files["no_end.pdb"] = "\n".join(seven)
    lines = [
        atom(1, "MG", "MG", "A", 1, (1, 2, 3), "HETATM", element="MG", charge="2+"),
        atom(2, "CL", "CL", "A", 2, (1, 2, 3), "HETATM", element="Cl", charge="1-"),
        atom(3, "NA", "NA", "A", 3, (1, 2, 3), "HETATM", element="na", charge="0 "),
        atom(4, "X1", "UNK", "A", 4, (1, 2, 3), element="XX"),
        atom(5, "QQ", "UNK", "A", 4, (1, 2, 3), element=" Q"),
        atom(6, "CA", "CAL", "A", 5, (1, 2, 3), element="  "),
        atom(7, "1HB", "ALA", "A", 6, (1, 2, 3), element="  "),
        atom(8, "OW", "SOL", "A", 7, (1, 2, 3), element=" O"),
        atom(9, "FE2", "HEM", "A", 8, (1, 2, 3), element="FE"),
        atom(10, "K", "K", "A", 9, (1, 2, 3), "HETATM", element=" K", charge="1+"),
        atom(11, "ZN", "ZN", "A", 10, (1, 2, 3), "HETATM", element="Zn"),
        atom(12, "D", "DOD", "A", 11, (1, 2, 3), element=" D"),
    ]
    files["elements.pdb"] = "\n".join(lines) + "\nEND\n"
    files["bad_charges.pdb"] = "\n".join(lines[:3] + [atom(13, "O", "HOH", "A", 12, (0, 0, 0), element=" O",
                                                           charge="+2")]) + "\n"
    files["text_charges.pdb"] = "\n".join(lines[:2] + [atom(13, "O", "HOH", "A", 12, (0, 0, 0), element=" O",
                                                            charge="ab")]) + "\n"
    names = ["CA", "1HB", "Cl-", "MG2+", "OW", "SOD", "K+", "C0", "CAL", "BR", "*", "+-", "HO2'", "C4*",
             "QL", "MW", "LIT", "1", "12", "NE2", "SE", "Na+", "CLA", "HG1", "XYZW", "é", "ß"]
    files["guess.pdb"] = "\n".join(atom(k + 1, n, "RES", "A", 1 + k, (k, 0, 0), element="  ")
                                   for k, n in enumerate(names)) + "\n"
    files["guess_blank_name.pdb"] = "\n".join(
        [atom(1, "", "RES", "A", 1, (0, 0, 0)), atom(2, "C1", "RES", "A", 1, (1, 0, 0))]) + "\n"
    seg = []
    for k in range(9):
        seg.append(atom(k + 1, "C", "LIG", "AB"[k // 5], 1 + k // 2, (k, k, k), segid=["SEG1", "SEG1", "SEG2",
                                                                                       "SEG2", "", "", "SEG1",
                                                                                       "LONGER", "S"][k]))
    files["segids.pdb"] = "\n".join(seg) + "\n"
    files["chains_only.pdb"] = "\n".join(atom(k + 1, "C", "LIG", "AAB B"[k], 1, (k, 0, 0)) for k in range(5)) + "\n"
    files["no_chains.pdb"] = "\n".join(atom(k + 1, "C", "LIG", " ", k, (k, 0, 0)) for k in range(4)) + "\n"
    files["icodes.pdb"] = "\n".join([
        atom(1, "N", "ALA", "A", 5, (0, 0, 0)), atom(2, "N", "ALA", "A", 5, (0, 0, 0), icode="A", alt="B"),
        atom(3, "N", "ALA", "A", 5, (0, 0, 0), icode="A", alt="C"), atom(4, "N", "GLY", "A", 5, (0, 0, 0), icode="A"),
        atom(5, "N", "GLY", "B", 5, (0, 0, 0), icode="A")]) + "\n"
    wrap = [atom(k + 1, "C", "WRP", "A", r, (k, 0, 0)) for k, r in enumerate([9998, 9999, 0, 1, -999, 5000, "    ",
                                                                              "  x ", 7, 3])]
    files["resid_wrap.pdb"] = "\n".join(wrap) + "\n"
    files["hy36.pdb"] = "\n".join([atom("99999", "C", "A", "A", 1, (0, 0, 0)), atom("A0000", "C", "A", "A", 1, (1, 0, 0)),
                                   atom("A0001", "C", "A", "A", 1, (2, 0, 0)), atom("a0000", "C", "A", "A", 1, (3, 0, 0)),
                                   "CONECT99999A0000", "CONECTA0000A0001"]) + "\n"
    files["hy36_ok.pdb"] = "\n".join([atom("99999", "C", "A", "A", 1, (0, 0, 0)),
                                      atom("A0000", "C", "A", "A", 1, (1, 0, 0)),
                                      atom("a0001", "C", "A", "A", 1, (2, 0, 0)), "CONECT99999    2"]) + "\n"
    files["hy36_bad.pdb"] = "\n".join([atom("Aa000", "C", "A", "A", 1, (0, 0, 0))]) + "\n"
    files["serial_stars.pdb"] = "\n".join([atom("*****", "C", "A", "A", 1, (0, 0, 0)),
                                           atom("*****", "C", "A", "A", 1, (1, 0, 0)), atom(5, "C", "A", "A", 1, (2, 0, 0)),
                                           "CONECT    5    5"]) + "\n"
    base = chain_atoms(5)
    files["conect_variants.pdb"] = "\n".join(base + ["TER       6      G A   2", "CONECT    1    2    3",
                                                      "CONECT    2    1", "CONECT    3    6", "CONECT    4",
                                                      "CONECT    5    4    4",
                                                      "conect    1    5", " CONECT    1    4"]) + "\n"
    files["conect_corrupt.pdb"] = "\n".join(base + ["CONECT    1    2", "CONECT    1    2  3"]) + "\n"
    files["conect_bad_int.pdb"] = "\n".join(base + ["CONECT    1    x"]) + "\n"
    files["conect_after_end.pdb"] = "\n".join(base + ["END", "CONECT    1    3"]) + "\n"
    files["conect_duplicate_serial.pdb"] = "\n".join(base + [atom(1, "O", "HOH", "W", 9, (5, 5, 5)),
                                                             "CONECT    1    2"]) + "\n"
    files["occupancy.pdb"] = "\n".join([atom(1, "C", "A", "A", 1, (0, 0, 0), occ="  0.50", temp=" 12.25"),
                                        atom(2, "C", "A", "A", 1, (0, 0, 0), occ="      ", temp="      "),
                                        atom(3, "C", "A", "A", 1, (0, 0, 0), occ="  x   ", temp="   nan"),
                                        atom(4, "C", "A", "A", 1, (0, 0, 0), occ=" 1_0  ", temp="1e2   "),
                                        atom(5, "C", "A", "A", 1, (0, 0, 0), occ="0.1234", temp="0.3333")]) + "\n"
    files["short_lines.pdb"] = "\n".join(["ATOM      1  C   ALA A   1       1.000   2.000   3.000",
                                          "ATOM      2  N   ALA A   1       1.000   2.000   3.000  1.00",
                                          "HETATM    3 MG    MG A   2       0.500  -0.250  10.125"]) + "\n"
    files["lowercase.pdb"] = "\n".join(seven[:3] + [seven[3].lower()] + seven[4:]) + "\n"
    files["leading_space.pdb"] = "\n".join(seven[:3] + [" " + seven[3]] + seven[4:]) + "\n"
    files["unicode.pdb"] = "\n".join([atom(1, "Cé", "ÅÄÖ", "ü", 1, (0, 0, 0), element=" C"),
                                      "REMARK   3 Grüße", "COMPND    π", "HEADER    Ångström"] +
                                     seven[1:]) + "\n"
    models = []
    for m in range(3):
        models.append(f"MODEL     {m + 1:4d}")
        models.append(f"CRYST1{20.0 + m:9.3f}{21.0:9.3f}{22.0:9.3f}{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1")
        models += chain_atoms(7, offset=0.125 * m)
        models.append("ENDMDL")
    files["models_npt.pdb"] = "\n".join(models + ["END"]) + "\n"
    models = ["CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 1           1"]
    for m in range(3):
        models.append(f"MODEL     {m + 1:4d}")
        models += chain_atoms(7, offset=-0.5 * m)
        models.append("ENDMDL")
    files["models_one_cryst.pdb"] = "\n".join(models + ["END"]) + "\n"
    models = []
    for m in range(4):
        models.append(f"MODEL     {m + 1:4d}")
        if m == 1:
            models.append("CRYST1   12.000   13.000   14.000  80.00  85.00  95.00 P 1           1")
        if m == 3:
            models.append("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1")
        models += chain_atoms(7, offset=m)
        models.append("ENDMDL")
    files["models_some_cryst.pdb"] = "\n".join(models) + "\n"
    files["models_mismatch.pdb"] = "\n".join(["MODEL        1"] + chain_atoms(7) + ["ENDMDL", "MODEL        2"]
                                             + chain_atoms(6) + ["ENDMDL"]) + "\n"
    files["models_more.pdb"] = "\n".join(["MODEL        1"] + chain_atoms(5) + ["ENDMDL", "MODEL        2"]
                                         + chain_atoms(6) + ["ENDMDL"]) + "\n"
    files["cryst_bad.pdb"] = "\n".join(["CRYST1   abc      1.000    1.000  90.00  90.00  90.00 P 1"] + seven) + "\n"
    files["cryst_zero.pdb"] = "\n".join(["CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1"] + seven) + "\n"
    files["cryst_short.pdb"] = "\n".join(["CRYST1   10.000   10.000"] + seven) + "\n"
    files["cryst_two.pdb"] = "\n".join(["CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1",
                                        "CRYST1   11.000   12.000   13.000  90.00  90.00  90.00 P 1"] + seven) + "\n"
    files["bad_coordinate.pdb"] = "\n".join(seven[:2] + [seven[2][:30] + "   1.0x0" + seven[2][38:]] + seven[3:]) + "\n"
    files["header_cut_utf8.pdb"] = "HEADER    " + "x" * 55 + "é\n" + "\n".join(seven) + "\n"
    files["two_headers.pdb"] = "HEADER    first\nHEADER\nHEADER    \t second  \n" + "\n".join(seven) + "\n"
    files["empty.pdb"] = ""
    files["remarks_only.pdb"] = "REMARK only\nEND\n"
    files["endmdl_first.pdb"] = "\n".join(seven[:3] + ["ENDMDL"] + seven[3:]) + "\n"
    files["ions.pdb"] = "\n".join(
        atom(k + 1, name, resname, "A", resid, (1.5 * k, 0.5, -0.25 * k), "ATOM" if resid <= 3 else "HETATM")
        for resname, resid, name, k in ion_atoms()) + "\nEND\n"
    raw = {name: text.encode() for name, text in files.items()}
    raw["invalid_utf8.pdb"] = ("\n".join(seven[:2]) + "\n").encode() + b"REMARK \xff\n" + "\n".join(seven[2:]).encode()
    raw["latin1_name.pdb"] = "\n".join(seven).encode().replace(b"C4'", b"C\xe9'")
    return raw


def gro_files():
    def gro_line(resid, resname, name, index, xyz, precision=3, vel=None):
        width = precision + 5
        line = f"{resid % 100000:5d}{resname:<5}{name:>5}{index % 100000:5d}"
        line += "".join(f"{v:{width}.{precision}f}" for v in xyz)
        if vel:
            line += "".join(f"{v:{width + 1}.{precision + 1}f}" for v in vel)
        return line

    files = {}
    atoms = [gro_line(1 + k // 3, "SOL", ["OW", "HW1", "HW2"][k % 3], k + 1, (0.1 * k, 0.2, 0.3 - 0.05 * k))
             for k in range(7)]
    files["seven.gro"] = "\n".join(["Seven atoms", "    7"] + atoms + ["   3.00000   4.00000   5.00000"]) + "\n"
    files["triclinic.gro"] = "\n".join(["Triclinic", "7"] + atoms +
                                       ["   3.00000   4.00000   5.00000   0.00000   0.00000   1.00000   0.00000"
                                        "   1.50000   2.00000"]) + "\n"
    files["zero_box.gro"] = "\n".join(["Zero box", "7"] + atoms + ["   0.00000   0.00000   0.00000"]) + "\n"
    files["packed_box.gro"] = "\n".join(["Packed", "7"] + atoms + ["3.000004.000005.00000"]) + "\n"
    files["bad_box.gro"] = "\n".join(["Bad box", "7"] + atoms + ["   3.0   4.0"]) + "\n"
    files["precise.gro"] = "\n".join(["Precise", "7"] + [
        gro_line(1, "LIG", f"C{k}", k + 1, (0.123456 * k, -1.5, 2.25), precision=5) for k in range(7)]
        + ["   3.00000   4.00000   5.00000"]) + "\n"
    files["velocities.gro"] = "\n".join(["Velocities", "7"] + [
        gro_line(1, "LIG", f"C{k}", k + 1, (0.5 * k, 1.5, 2.5), vel=(0.1, 0.2, 0.3)) for k in range(7)]
        + ["   3.00000   4.00000   5.00000"]) + "\r\n"
    wrap = [gro_line(r, "RES", "CA", k + 1, (k, 0, 0)) for k, r in enumerate([99998, 99999, 100000, 100000, 100001,
                                                                               200000, 200001])]
    files["resid_wrap.gro"] = "\n".join(["Wrap", "7"] + wrap + ["   1.00000   1.00000   1.00000"]) + "\n"
    first_zero = [gro_line(r, "RES", "CA", k + 1, (k, 0, 0)) for k, r in enumerate([0, 0, 1, 0, 1, 0, 0])]
    files["resid_zero.gro"] = "\n".join(["Zero", "7"] + first_zero + ["   1.00000   1.00000   1.00000"]) + "\n"
    files["one_atom.gro"] = "One\n1\n" + gro_line(1, "NA", "NA", 1, (1, 2, 3)) + "\n   2.0   2.0   2.0\n"
    files["bad_line.gro"] = "\n".join(["Bad", "7"] + atoms[:3] + ["    xSOL     OW    4   0.400   0.200   0.100"]
                                      + atoms[4:] + ["   3.00000   4.00000   5.00000"]) + "\n"
    files["missing_name.gro"] = "\n".join(["Missing", "7"] + atoms[:2] + [atoms[2][:10] + "     " + atoms[2][15:]]
                                          + atoms[3:] + ["   3.00000   4.00000   5.00000"]) + "\n"
    files["too_few.gro"] = "\n".join(["Few", "7"] + atoms[:4]) + "\n"
    files["no_box.gro"] = "\n".join(["No box", "7"] + atoms) + "\n"
    files["bad_count.gro"] = "Count\n seven\n" + "\n".join(atoms) + "\n"
    files["ions.gro"] = "\n".join(["Ions", str(len(ion_atoms()))]
                                   + [gro_line(resid, resname, name, k + 1, (0.15 * k, 0.05, 0.1))
                                      for resname, resid, name, k in ion_atoms()]
                                   + ["   4.00000   4.00000   4.00000"]) + "\n"
    return {name: text.encode() for name, text in files.items()}


# --- topology and reader cases ---------------------------------------------------------------------

def reader_info(path, dt=None):
    kwargs = {"dt": dt} if dt else {}
    u = mda.Universe(path, **kwargs)
    traj = u.trajectory
    out = {"dt_given": dt, "n_frames": int(traj.n_frames), "dt": float(traj.dt).hex()}
    for key in ("header", "compound", "remarks"):
        out[key] = getattr(traj, key, None)
    order = list(range(traj.n_frames))
    shuffled = order[::-1] + order
    random.Random(3).shuffle(shuffled)
    frames = []
    for k in order + shuffled:
        try:
            ts = traj[k]
            frames.append({"frame": k, "time": float(ts.time).hex(), "positions": f32bits(ts.positions),
                           "dimensions": None if ts.dimensions is None else f32bits(ts.dimensions)})
        except Exception as exc:  # noqa: BLE001
            frames.append({"frame": k, "error": error(exc)})
            break
    out["frames"] = frames
    return out


def topology_cases(files):
    cases = []
    for name in files:
        path = os.path.join(FILES, name)
        case = {"file": name}
        try:
            data = dump_topologies.dump(path)
        except Exception as exc:  # noqa: BLE001
            case["error"] = error(exc)
            cases.append(case)
            continue
        target = os.path.join(FILES, name + ".topology.json")
        with open(target, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=1, ensure_ascii=False)
            fh.write("\n")
        case["dump"] = name + ".topology.json"
        case["readers"] = [reader_info(path), reader_info(path, 2.5)]
        cases.append(case)
    return cases


# --- loader cases ---------------------------------------------------------------------------------

def summary_case(topology, trajectories, dt):
    case = {"topology": placeholders(topology), "trajectories": [placeholders(t) for t in trajectories], "dt": dt}
    try:
        u = loader.load_universe(topology, trajectories, dt)
        summary = loader.system_summary(u)
        names, counts = np.unique(u.residues.resnames, return_counts=True)
        # numpy's argsort orders equal counts by the CPU's sort kernel: the full counts let the
        # test check MDBridge's order (equal counts by name)
        case["resname_counts"] = [[str(n), int(c)] for n, c in zip(names, counts)]
        case["summary"] = summary
        case["summary_json"] = json.dumps(summary, indent=2)
    except Exception as exc:  # noqa: BLE001
        case["error"] = error(exc)
    return case


def loader_cases():
    top = os.path.join(DATA, "topology", "demo_dry.parm7")
    nc = os.path.join(DATA, "io", "trajectories", "demo_dry.nc")
    tf = os.path.join(DATA, "io", "trajectories", "files")
    demo = [os.path.join(tf, f) for f in ("demo_dry.dcd", "demo_dry.xtc", "demo_dry.trr")]
    f = lambda name: os.path.join(FILES, name)  # noqa: E731
    seven_traj = [os.path.join(tf, n) for n in ("small.xtc", "small.trr", "ortho.dcd", "mda_box.nc", "time.rst7")]
    matrix = [
        (top, [nc], None), (top, [nc], 0), (top, [nc], 2.0), (top, [nc], 0.25), (top, [nc, nc], None),
        (top, [nc], float("nan")),
        (top, demo[:1], None), (top, demo[1:2], None), (top, demo[2:], None), (top, [nc] + demo, None),
        (top, [nc] + demo, 3.0),
        (top, [], None), ("", [nc], None), (os.path.join(DATA, "missing.parm7"), [nc], None),
        (top, [os.path.join(DATA, "missing.nc")], None), (top, [nc, os.path.join(DATA, "gone.dcd")], None),
        (os.path.join(DATA, "topology", "fixture_chains.pdb"), [], None),
        (os.path.join(DATA, "topology", "fixture_segids.pdb"), [], 4.0),
        (os.path.join(DATA, "topology", "fixture_chains_ba.pdb"), [], None),
        (f("seven.pdb"), [], None), (f("seven.gro"), [], None), (f("seven.gro"), [], 5.0),
        (f("models_npt.pdb"), [], None), (f("models_npt.pdb"), [], 0.5), (f("models_one_cryst.pdb"), [], None),
        (f("seven.pdb"), [f("models_npt.pdb")], None), (f("seven.pdb"), [f("seven.gro")], None),
        (f("seven.gro"), [f("models_some_cryst.pdb"), f("seven.pdb")], None),
        (f("seven.pdb"), [nc], None), (top, [f("seven.pdb")], None),
        (f("hy36_bad.pdb"), [], None), (f("conect_bad_int.pdb"), [], None), (f("bad_line.gro"), [], None),
        (f("invalid_utf8.pdb"), [], None), (f("empty.pdb"), [], None), (f("models_mismatch.pdb"), [], None),
        (f("bad_coordinate.pdb"), [], None), (f("header_cut_utf8.pdb"), [], None),
        (os.path.join(tf, "small.xtc"), [], None),
    ]
    for t in seven_traj:
        matrix.append((f("seven.pdb"), [t], None))
        matrix.append((f("seven.gro"), [t], 1.5))
    matrix.append((f("seven.pdb"), seven_traj[:4], None))
    matrix.append((f("seven.pdb"), seven_traj[:2] + [seven_traj[0]], None))
    matrix.append((f("seven.pdb"), [os.path.join(tf, "ten.mdcrd")], None))
    matrix.append((f("seven.pdb"), [os.path.join(tf, "box.mdcrd")], None))
    matrix.append((f("seven.pdb"), [os.path.join(tf, "time.rst7"), os.path.join(tf, "time.inpcrd")], None))
    summaries = [summary_case(*m) for m in matrix]

    frames = []
    for n in (0, 1, 5, 24):
        for start in (0, 1, 3, 5, 24, 30, -1, -3, -30):
            for stop in (None, 0, 2, 5, 24, 100, -1, -5, -100):
                for step in (1, 2, 7, 0, -1):
                    case = {"n": n, "start": start, "stop": stop, "step": step}
                    try:
                        case["frames"] = [int(v) for v in loader.frame_indices(n, start, stop, step)]
                    except Exception as exc:  # noqa: BLE001
                        case["error"] = error(exc)
                    frames.append(case)

    info = []
    for topology, trajectories in ((top, [nc]), (f("seven.pdb"), []), (top, demo[:1])):
        args = ["info", "--top", topology] + (["--traj", *trajectories] if trajectories else [])
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            cli.main(args)
        u = loader.load_universe(topology, trajectories)
        names, counts = np.unique(u.residues.resnames, return_counts=True)
        info.append({"resname_counts": [[str(n), int(c)] for n, c in zip(names, counts)],
                     "topology": placeholders(topology), "trajectories": [placeholders(t) for t in trajectories],
                     "stdout": placeholders(stdout.getvalue())})
    return {"summaries": summaries, "frame_indices": frames, "info": info}


def remove_offsets():
    """MDAnalysis caches XTC / TRR frame offsets next to the files it reads: remove them."""
    for root, _dirs, names in os.walk(DATA):
        for name in names:
            if name.startswith(".") and re.search(r"_offsets\.(lock|npz)$", name):
                os.remove(os.path.join(root, name))


def main():
    os.makedirs(FILES, exist_ok=True)
    files = {**pdb_files(), **gro_files()}
    for name, raw in files.items():
        with open(os.path.join(FILES, name), "wb") as fh:
            fh.write(raw)
    cases = topology_cases(files)
    with open(os.path.join(OUT, "topology_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(cases, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    with open(os.path.join(OUT, "loader_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(loader_cases(), fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")


if __name__ == "__main__":
    main()
    remove_offsets()
