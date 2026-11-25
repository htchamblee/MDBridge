"""The PDB files MDATK's frame export writes (mdatk/io/pdb_export.py export_frames, through
MDAnalysis 2.11.0's PDBWriter) for tests/data/io/pdb/export.

    /opt/mdbridge-venv/bin/python -I -B make_export_cases.py MDATK_CHECKOUT OUT_DIR

Run make_loader_cases.py first: some cases use its files (tests/data/io/loader/files). Nothing is
written into the checkout; the output is reproducible. Writes

  files/*               inputs: a 40-atom NetCDF trajectory for the parm7 chain-ID variants of
                        tests/data/io/parm7/topology, PDB / GRO files with chain-ID variants and
                        coordinates at and beyond the PDB limits
  export_cases.json     per case: the system (topology, trajectories, dt), the selection, frames,
                        prefix and stamp given to export_frames, and the files it wrote (the text,
                        or its SHA-256 and first lines when long) with the progress calls, or the
                        error and the files left behind; `big`: a 100001-atom system built from
                        formulas (the C++ test builds the same) and the SHA-256 of its exports

Paths are written with the test data folder as "{data}". Errors are {"type", "message"}.
"""

import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_export_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402
from MDAnalysis.coordinates.TRJ import NCDFWriter  # noqa: E402

from mdatk import selections  # noqa: E402
from mdatk.io import loader, pdb_export  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import improvements_readers_pdb  # noqa: E402
import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)

DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
OUT = os.path.abspath(sys.argv[2])
FILES = os.path.join(OUT, "files")
LOADER = os.path.join(DATA, "io", "loader", "files")
LIMIT = 16384  # files up to this size are stored as text


def placeholders(text):
    return text.replace(DATA, "{data}")


def error(exc):
    return {"type": type(exc).__name__, "message": placeholders(str(exc))}


def file_entry(path):
    with open(path, "rb") as fh:
        raw = independent_tables.residue_rule_pdb_text(fh.read())  # MDBridge's residue-scoped ion rule
    entry = {"name": os.path.basename(path)}
    if len(raw) <= LIMIT:
        entry["text"] = raw.decode()
    else:
        entry["sha256"] = hashlib.sha256(raw).hexdigest()
        entry["size"] = len(raw)
        entry["head"] = raw[:1500].decode()
    return entry


def export_case(topology, trajectories=(), dt=None, selection="", frames=(0,), prefix="frame", stamp=""):
    with tempfile.TemporaryDirectory(prefix="equivalent-") as eq:
        # the loader's PDB / GRO files as MDBridge reads them (improvements_readers_pdb.py)
        same = []
        for k, p in enumerate([topology, *trajectories]):
            if os.path.dirname(p) == LOADER:
                os.makedirs(os.path.join(eq, str(k)))
                p = improvements_readers_pdb.equivalent_file(p, os.path.join(eq, str(k)))
            same.append(p)
        return _export_case(topology, trajectories, same, dt, selection, frames, prefix, stamp)


def _export_case(topology, trajectories, same, dt, selection, frames, prefix, stamp):
    case = {"topology": placeholders(topology), "trajectories": [placeholders(t) for t in trajectories],
            "dt": dt, "selection": selection, "frames": list(frames), "prefix": prefix, "stamp": stamp}
    u = loader.load_universe(same[0], same[1:], dt)
    atoms = u.atoms if not selection else selections.select(u, selection, [])
    calls = []
    with tempfile.TemporaryDirectory(prefix="pdb-export-") as tmp:
        out_dir = os.path.join(tmp, "out")
        try:
            paths = pdb_export.export_frames(u, list(frames), atoms, out_dir, prefix,
                                             progress=lambda f, m: calls.append([f.hex(), m]), stamp=stamp)
            case["paths"] = [os.path.relpath(p, out_dir) for p in paths]
        except Exception as exc:  # noqa: BLE001
            case["error"] = error(exc)
        case["progress"] = calls
        case["files"] = [file_entry(os.path.join(out_dir, n)) for n in sorted(os.listdir(out_dir))] \
            if os.path.isdir(out_dir) else []
    return case


def pdb_atom(serial, name, resname, chain, resid, xyz, element="C", record="ATOM  ", occ="  1.00", temp="  0.00",
             segid="    ", charge="  ", alt=" ", icode=" "):
    name = name if len(name) >= 4 else " " + name.ljust(3)
    x, y, z = (v if isinstance(v, str) else f"{v:8.3f}" for v in xyz)
    return (f"{record}{serial:5d} {name}{alt}{resname:>3} {chain}{resid:4d}{icode}   {x}{y}{z}{occ}{temp}"
            f"      {segid:<4}{element:>2}{charge}")


def write(name, text):
    path = os.path.join(FILES, name)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    return path


def inputs():
    files = {}
    # a 40-atom trajectory for the parm7 variants (3 frames, with box)
    rng = np.random.default_rng(11)
    path = os.path.join(FILES, "forty.nc")
    u = mda.Universe.empty(40, trajectory=True)
    coordinates = rng.uniform(-30, 30, size=(3, 40, 3)).astype(np.float32)
    dimensions = np.array([[40.0 + k, 41.0, 42.5, 90.0, 90.0, 90.0] for k in range(3)], dtype=np.float32)
    u.load_new(coordinates, format=MemoryReader, dimensions=dimensions, dt=2.0)
    with NCDFWriter(path, 40) as w:
        for _ts in u.trajectory:
            w.write(u.atoms)
    files["forty.nc"] = path
    chains = ["A", "a", "1", "é", "*", " ", "_", "B"]
    lines = ["HEADER    CHAIN VARIANTS", "COMPND    MOL_ID: 1;",
             "REMARK   1 " + "word " * 30]
    for k, chain in enumerate(chains):
        lines.append(pdb_atom(k + 1, ["P", "MG", "CL", "HG", "FE", "CA", "NA", "C1'"][k],
                              ["G", "MG", "CPT", "CMH", "HEM", "CA", "NA", "OEC"][k], chain, k + 1,
                              (k * 1.5, -k * 0.25, 3.0), element=["P", "MG", "Cl", "HG", "Fe", "CA", "NA", "c"][k],
                              charge=["  ", "2+", "1-", "  ", "3+", "  ", "1+", "  "][k],
                              record="HETATM" if k % 3 == 1 else "ATOM  ",
                              occ=" 0.50" if k == 2 else "  1.00", temp=" 99.99" if k == 3 else "  0.00",
                              alt="B" if k == 4 else " ", icode="C" if k == 5 else " ",
                              segid="SEGMENTX" if k == 6 else "S1"))
    lines += ["CONECT    1    2    3", "CONECT    2    1", "CONECT    3    1    8", "CONECT    8    3", "END"]
    files["chains.pdb"] = write("chains.pdb", "\n".join(lines) + "\n")
    names = ["CA", "OEC", "CA1", "PD", "C14", "NI", "CO1", "FE2", "UNK", "UNL", "HG1", "CL", "ZN", "X", "", "LONG",
             "ABCDE", "C"]
    resnames = ["CA", "OEC", "OEC", "PLL", "C14", "BLM", "VNL", "PF5", "ALA", "UNL", "EMC", "UN9", "ZN2", "X",
                "ALA", "ALA", "ALA", "UNK"]
    lines = []
    for k, (n, r) in enumerate(zip(names, resnames)):
        field = n.ljust(4) if len(n) == 4 else (" " + n.ljust(3) if len(n) < 4 else n[:4])
        lines.append(f"ATOM  {k + 1:5d} {field} {r:>3} A{k + 1:4d}    {k:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00"
                     f"           C")
    files["names.pdb"] = write("names.pdb", "\n".join(lines) + "\nEND\n")
    lim = [("-999.999", "9999.999"), ("-999.999", "   0.000"), ("-1000.00", "   0.000"), ("99999.99", "   0.000"),
           ("     nan", "   0.000")]
    for k, (a, b) in enumerate(lim):
        files[f"limit{k}.pdb"] = write(f"limit{k}.pdb", pdb_atom(1, "C", "ALA", "A", 1, (a, b, "   1.000")) + "\n"
                                       + pdb_atom(2, "N", "ALA", "A", 1, (1, 2, 3)) + "\nCONECT    1    2\nEND\n")
    models = []
    for m, x in enumerate(["   1.000", "  -1.500", "-1000.00", "   2.000"]):
        models += [f"MODEL     {m + 1:4d}", pdb_atom(1, "C", "ALA", "A", 1, (x, "   0.000", "   0.000")),
                   pdb_atom(2, "N", "ALA", "A", 1, (1, 2, 3)), "ENDMDL"]
    files["limit_models.pdb"] = write("limit_models.pdb", "\n".join(models) + "\nEND\n")
    gro = ["Limits", "4",
           f"{1:5d}{'LIG':<5}{'C1':>5}{1:5d}{-99.99994:12.7f}{999.99994:12.7f}{0.0:12.7f}",
           f"{1:5d}{'LIG':<5}{'C2':>5}{2:5d}{-99.99995:12.7f}{999.99995:12.7f}{0.0:12.7f}",
           f"{1:5d}{'LIG':<5}{'C3':>5}{3:5d}{-99.99996:12.7f}{999.99990:12.7f}{0.0:12.7f}",
           f"{1:5d}{'LIG':<5}{'C4':>5}{4:5d}{0.5:12.7f}{999.99993:12.7f}{0.0:12.7f}",
           "   1.00000   1.00000   1.00000"]
    files["limit.gro"] = write("limit.gro", "\n".join(gro) + "\n")
    # single-atom GRO files at each value near the limits
    values = [-99.99995, -99.999951, -99.99994, 999.99995, 999.99991, 999.99990, 999.999905]
    for k, v in enumerate(values):
        files[f"edge{k}.gro"] = write(f"edge{k}.gro", "\n".join(
            ["Edge", "1", f"{1:5d}{'LIG':<5}{'C1':>5}{1:5d}{v:13.8f}{0.0:13.8f}{0.0:13.8f}", "   0.0   0.0   0.0"]) + "\n")
    return files


def cases(files):
    top = os.path.join(DATA, "topology", "demo_dry.parm7")
    nc = os.path.join(DATA, "io", "trajectories", "demo_dry.nc")
    tf = os.path.join(DATA, "io", "trajectories", "files")
    out = [
        export_case(top, [nc], frames=[0, 5, 23], prefix="demo"),
        export_case(top, [nc], dt=2.5, selection="resname G and name P C4'", frames=[3], prefix="sel",
                    stamp="20261124-101500"),
        export_case(top, [nc], selection="resid 1:3", frames=[23, 0, -1], prefix="neg"),
        export_case(top, [os.path.join(tf, "demo_dry.dcd")], selection="resid 2 to 4", frames=[1]),
        export_case(top, [os.path.join(tf, "demo_dry.xtc")], selection="resid 5", frames=[0, 2]),
        export_case(top, [nc], frames=[], prefix="none"),
        export_case(top, [nc], selection="name MG", frames=[7], prefix="ions"),
    ]
    seven = os.path.join(LOADER, "seven.pdb")
    seven_gro = os.path.join(LOADER, "seven.gro")
    for traj in ("box.mdcrd", "ten.mdcrd", "nobox.mdcrd", "title_unicode.mdcrd", "time.rst7", "time_box.rst7",
                 "small.xtc", "mda_box.nc", "mda_nobox.nc"):
        try:
            out.append(export_case(seven, [os.path.join(tf, traj)], frames=[0], prefix=traj.replace(".", "_")))
        except Exception as exc:  # noqa: BLE001 - files of other atom counts
            print("skipped", traj, exc, file=sys.stderr)
    out.append(export_case(seven, [os.path.join(tf, "small.xtc"), os.path.join(tf, "ortho.dcd")], frames=[1, 6]))
    out.append(export_case(seven, frames=[0], prefix="self"))
    out.append(export_case(seven_gro, frames=[0], prefix="gro"))
    out.append(export_case(seven_gro, [os.path.join(tf, "small.trr")], selection="index 0:3", frames=[2]))
    for name in ("models_npt.pdb", "models_one_cryst.pdb", "models_some_cryst.pdb", "elements.pdb", "unicode.pdb",
                 "icodes.pdb", "occupancy.pdb", "segids.pdb", "conect_variants.pdb", "conect_duplicate_serial.pdb",
                 "resid_wrap.pdb", "hy36_ok.pdb", "two_headers.pdb", "triclinic.gro", "resid_wrap.gro",
                 "guess.pdb", "lowercase.pdb", "short_lines.pdb", "ions.pdb", "ions.gro"):
        path = os.path.join(LOADER, name)
        u = mda.Universe(path)
        frames = list(range(u.trajectory.n_frames))
        out.append(export_case(path, frames=frames, prefix=name.replace(".", "_")))
    for name in ("fixture_chains.pdb", "fixture_segids.pdb", "fixture_noelements.pdb", "fixture_chains_ba.pdb"):
        out.append(export_case(os.path.join(DATA, "topology", name), prefix="fx"))
        out.append(export_case(os.path.join(DATA, "topology", name), selection="resid 2:4", prefix="fxsel"))
    out.append(export_case(files["chains.pdb"], prefix="chains"))
    out.append(export_case(files["chains.pdb"], selection="index 0 2 7", prefix="chains_sel"))
    out.append(export_case(files["names.pdb"], prefix="names"))
    for variant in ("chain_long_names", "chain_blank", "chain_each", "chain_all_a"):
        out.append(export_case(os.path.join(DATA, "io", "parm7", "topology", variant + ".parm7"), [files["forty.nc"]],
                               frames=[0, 2], prefix=variant))
    for k in range(5):
        out.append(export_case(files[f"limit{k}.pdb"], prefix=f"limit{k}"))
    out.append(export_case(files["limit_models.pdb"], frames=[0, 1, 2, 3], prefix="models"))
    out.append(export_case(files["limit_models.pdb"], selection="name N", frames=[2], prefix="models_n"))
    out.append(export_case(files["limit.gro"], prefix="limit_gro"))
    for k in range(4):
        out.append(export_case(files["limit.gro"], selection=f"index {k}", prefix=f"limit_gro{k}"))
    for k in range(7):
        out.append(export_case(files[f"edge{k}.gro"], prefix=f"edge{k}"))
    return out


# --- a big system: serials beyond 99999, resids beyond 9999 and below -999 -------------------------

def big_universe(n_atoms=100001):
    """The system the C++ test builds from the same formulas."""
    per_residue = 3
    n_res = (n_atoms + per_residue - 1) // per_residue
    resindex = np.arange(n_atoms) // per_residue
    u = mda.Universe.empty(n_atoms, n_residues=n_res, n_segments=2, atom_resindex=resindex,
                           residue_segindex=(np.arange(n_res) >= n_res // 2).astype(int), trajectory=True)
    u.add_TopologyAttr("names", [["C1", "N", "OXT"][k % 3] for k in range(n_atoms)])
    u.add_TopologyAttr("resnames", [["ALA", "GLY", "HOH", "MG"][r % 4] for r in range(n_res)])
    u.add_TopologyAttr("resids", [r * 7 - 20000 for r in range(n_res)])
    u.add_TopologyAttr("segids", ["PROT", "SOLVENT"])
    u.add_TopologyAttr("elements", [["C", "N", "O"][k % 3] for k in range(n_atoms)])
    bonds = [(k, k + 1) for k in range(0, 60, 3)] + [(n_atoms - 2, n_atoms - 1)]
    u.add_TopologyAttr("bonds", bonds)
    k = np.arange(n_atoms, dtype=np.float64)
    pos = np.stack([(k % 1000) * 0.5 - 250.0, (k // 1000) * 0.25 - 12.5,
                    ((k * 37) % 2001 - 1000) * 0.01], axis=1)
    u.atoms.positions = pos.astype(np.float32)
    u.dimensions = [500.0, 260.0, 30.0, 90.0, 90.0, 90.0]
    return u


def big_cases():
    u = big_universe()
    out = []
    for name, atoms in (("all", u.atoms), ("first", u.atoms[:60]), ("tail", u.atoms[99990:])):
        with tempfile.TemporaryDirectory(prefix="pdb-big-") as tmp:
            paths = pdb_export.export_frames(u, [0], atoms, tmp, "big")
            with open(paths[0], "rb") as fh:
                raw = independent_tables.residue_rule_pdb_text(fh.read())
        lines = raw.decode().splitlines()
        out.append({"selection": name, "sha256": hashlib.sha256(raw).hexdigest(), "size": len(raw),
                    "head": "\n".join(lines[:12]), "tail": "\n".join(lines[-12:])})
    return out


def remove_offsets():
    """MDAnalysis caches XTC / TRR frame offsets next to the files it reads: remove them."""
    for root, _dirs, names in os.walk(DATA):
        for name in names:
            if name.startswith(".") and re.search(r"_offsets\.(lock|npz)$", name):
                os.remove(os.path.join(root, name))


def main():
    if os.path.isdir(FILES):
        shutil.rmtree(FILES)
    os.makedirs(FILES)
    files = inputs()
    data = {"cases": cases(files), "big": big_cases()}
    with open(os.path.join(OUT, "export_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False)
        fh.write("\n")


if __name__ == "__main__":
    main()
    remove_offsets()
