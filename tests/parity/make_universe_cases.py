"""Probe cases for docs/specs/m3-universe-notes.md: the Universe / AtomGroup / trajectory behaviours
MDATK's analyses rely on that the milestone-2 data does not pin, as MDAnalysis 2.11.0 shows them.

    /opt/mdbridge-venv/bin/python -I -B make_universe_cases.py MDATK_CHECKOUT OUT_DIR

The MDATK checkout is not used (the argument keeps the generators' common calling convention);
nothing is written into it. MDAnalysis runs with MDBridge's independent tables patched in
(independent_tables.apply(), tests/data/README.md). Reads files of tests/data (topology/,
io/loader/files, io/trajectories) and writes:

  OUT_DIR/universe_files/*   small crafted inputs: a PDB whose CONECT records are out of order and
                             make several fragments, a PDB with a CONECT record without partners,
                             NetCDF / DCD / XTC / TRR files whose cells have zero lengths or angles
  OUT_DIR/universe_cases.json
    iteration     per system (topology + trajectory, paths relative to tests/data) and scenario:
                  the frame the trajectory was moved to first ("seek"), the int64 index array, how
                  many frames were consumed before stopping ("stop_after", null = all), the frames
                  yielded ([frame, time]) and the current frame and time afterwards, or the error
                  and when it was raised ("at": "index" when indexing the trajectory, "iterate"
                  while iterating) with the current frame afterwards. Every yielded frame's
                  positions and cell were checked against a plain read of that frame.
    dimensions    per crafted trajectory: the cell written and the cell read back per frame
                  (null for no cell)
    topology      per file: whether the atoms have charges / chainIDs (hasattr, as MDATK asks),
                  the error of reading them when not, the dtypes of masses / charges / positions,
                  the Universe's bonds in their order (or the error), fragindices of all atoms
                  (or the error), and per atom group (index lists, duplicates and unsorted orders
                  included): its bonds in order, fragindices, residues, the atoms of its residues,
                  a boolean-mask subset and a union with another group
    residues      per file: the atoms of residue groups given as index lists

Times are float.hex() of float(ts.time); cells are lists of floats (float32 values, exact in
JSON); errors are {"type", "message", "attribute_error", "value_error"} (whether the exception is
an AttributeError, which hasattr() swallows, and a ValueError). Index lists are written in full.
"""

import json
import os
import sys
import warnings

import numpy as np

if len(sys.argv) != 3:
    sys.exit("usage: make_universe_cases.py MDATK_CHECKOUT OUT_DIR")
OUT = os.path.abspath(sys.argv[2])
FILES = os.path.join(OUT, "universe_files")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

sys.path.insert(0, HERE)
import independent_tables  # noqa: E402

import MDAnalysis as mda  # noqa: E402
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402

warnings.filterwarnings("ignore")
independent_tables.apply()


def error(exc):
    return {"type": type(exc).__name__, "message": str(exc), "attribute_error": isinstance(exc, AttributeError),
            "value_error": isinstance(exc, ValueError)}


def path_of(rel):
    if rel.startswith("analysis/"):
        return os.path.join(OUT, rel[len("analysis/"):])
    return os.path.join(DATA, rel)


def universe(top, traj=None):
    return mda.Universe(path_of(top), path_of(traj)) if traj else mda.Universe(path_of(top))


def cell(ts):
    return None if ts.dimensions is None else [float(x) for x in ts.dimensions]


def remove_offsets(folder):
    for here, _dirs, names in os.walk(folder):
        for name in names:
            if name.startswith(".") and (name.endswith("_offsets.npz") or name.endswith("_offsets.lock")):
                os.remove(os.path.join(here, name))


# --- crafted files ---------------------------------------------------------------------------------

def write_files():
    os.makedirs(FILES, exist_ok=True)
    atoms = []
    for k in range(8):
        atoms.append("ATOM  %5d  C%-2d RNA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C"
                     % (k + 1, k, 1 + k // 3, k * 1.5, 0.25 * k, -0.5 * k))
    # CONECT out of order, both directions, a repeated pair: bonds 1-2 (as "7 2" ...), fragments
    # {0}, {1, 2, 6}, {3, 4, 7}, {5}
    frag = atoms + ["CONECT    7    2", "CONECT    5    4    8", "CONECT    3    2", "CONECT    2    7",
                    "CONECT    8    5", "END"]
    lonely = atoms[:4] + ["CONECT    1", "END"]
    for name, lines in (("frag_conect.pdb", frag), ("conect_lonely.pdb", lonely)):
        with open(os.path.join(FILES, name), "w", encoding="ascii", newline="\n") as fh:
            fh.write("\n".join(lines) + "\n")
    cells = {"zero_a": [0, 20, 20, 90, 90, 90], "zero_lengths": [0, 0, 0, 90, 90, 90],
             "zero_all": [0, 0, 0, 0, 0, 0], "zero_alpha": [20, 20, 20, 0, 90, 90],
             "ortho": [20, 21, 22, 90, 90, 90]}
    out = []
    for name, dims in cells.items():
        u = mda.Universe.empty(3, trajectory=True)
        xyz = np.arange(18, dtype=np.float32).reshape(2, 3, 3) * np.float32(0.5)
        u.load_new(xyz, format=MemoryReader, dimensions=np.array(dims, dtype=np.float32), dt=2.0)
        for ext in ("nc", "dcd", "xtc", "trr"):
            fname = f"cell_{name}.{ext}"
            with mda.Writer(os.path.join(FILES, fname), n_atoms=3) as w:
                for _ts in u.trajectory:
                    w.write(u.atoms)
            r = mda.Universe.empty(3)
            r.load_new(os.path.join(FILES, fname))
            out.append({"file": "analysis/universe_files/" + fname, "written": [float(x) for x in dims],
                        "read": [cell(ts) for ts in r.trajectory]})
    remove_offsets(FILES)
    return out


# --- trajectory iteration ----------------------------------------------------------------------------

SYSTEMS = [
    ("topology/demo_dry.parm7", "io/trajectories/demo_dry.nc"),
    ("topology/demo_dry.parm7", "io/trajectories/files/demo_dry.xtc"),
    ("topology/demo_dry.parm7", "io/trajectories/files/demo_dry.dcd"),
    ("topology/demo_dry.parm7", "io/trajectories/files/demo_dry.trr"),
    ("io/loader/files/models_npt.pdb", None),
    ("io/loader/files/seven.gro", None),
]


def scenarios(n, rng):
    last = n - 1
    out = [
        ("in_order", None, list(range(n)), None),
        ("reversed", None, list(range(n))[::-1], None),
        ("repeats", None, [last, last, 0, last, 0, 0], None),
        ("random", None, [int(x) for x in rng.integers(0, n, size=10)], None),
        ("mdatk_rmsd_seek", min(5, last), list(range(0, n, 2)), None),
        ("mdatk_rmsd_seek_last", last, list(range(n)), None),
        ("seek_then_empty", last, [], None),
        ("empty", None, [], None),
        ("stop_early", min(1, last), list(range(n))[::-1], 2),
        ("stop_after_one", None, [last, 0], 1),
        ("negative", None, [-1, -n, 0], None),
        ("out_of_range", min(1, last), [0, n], None),
        ("out_of_range_negative", min(1, last), [-n - 1], None),
        ("single", min(2, last), [min(2, last)], None),
    ]
    return out


def iteration_cases():
    rng = np.random.default_rng(20261124)
    out = []
    for top, traj in SYSTEMS:
        u = universe(top, traj)
        t = u.trajectory
        n = int(t.n_frames)
        plain = {}
        for k in range(n):
            ts = t[k]
            plain[k] = (ts.positions.copy(), cell(ts), float(ts.time).hex())
        for name, seek, indices, stop in scenarios(n, rng):
            t.rewind()
            if seek is not None:
                t[seek]
            entry = {"topology": top, "trajectory": traj, "n_frames": n, "scenario": name, "seek": seek,
                     "indices": indices, "stop_after": stop}
            frames = []
            try:
                it = t[np.array(indices, dtype=np.int64)]
            except Exception as exc:  # noqa: BLE001
                entry.update(error=error(exc), at="index", yielded=[],
                             after=[int(t.ts.frame), float(t.ts.time).hex()])
                out.append(entry)
                continue
            try:
                for ts in it:
                    pos, box, time = plain[int(ts.frame)]
                    assert np.array_equal(ts.positions, pos) and cell(ts) == box and float(ts.time).hex() == time
                    frames.append([int(ts.frame), float(ts.time).hex()])
                    if stop is not None and len(frames) == stop:
                        break
                entry.update(error=None, at=None)
            except Exception as exc:  # noqa: BLE001
                entry.update(error=error(exc), at="iterate")
            entry.update(yielded=frames, after=[int(t.ts.frame), float(t.ts.time).hex()])
            out.append(entry)
    remove_offsets(os.path.join(DATA, "io", "trajectories"))
    return out


# --- topology attributes ---------------------------------------------------------------------------

TOPOLOGIES = [
    ("topology/demo_dry.parm7", "io/trajectories/demo_dry.nc",
     [[30, 5, 5, 2], [0], [771], [400, 30, 31], [431, 30, 31], list(range(64, 96)), [772, 0, 771], []]),
    ("io/parm7/topology/chain_each.parm7", None, [[0, 1], [39, 0]]),
    ("topology/fixture_chains.pdb", None, [[3, 1], [0]]),
    ("io/loader/files/seven.pdb", None, [[2], [6, 0, 0, 3], [3, 4], []]),
    ("io/loader/files/conect_variants.pdb", None, [[1], [4, 2]]),
    ("analysis/universe_files/frag_conect.pdb", None, [[7, 4, 1], [0], [5, 5], [2], [6, 3], []]),
    ("analysis/universe_files/conect_lonely.pdb", None, [[0, 1], []]),
    ("io/loader/files/seven.gro", None, [[0, 1]]),
    ("io/loader/files/ions.gro", None, [[0, 1]]),
]


def attempt(fn):
    try:
        return {"value": fn(), "error": None}
    except Exception as exc:  # noqa: BLE001
        return {"value": None, "error": error(exc)}


def pairs(tg):
    return [[int(i), int(j)] for i, j in np.asarray(tg.indices).reshape(-1, 2)]


def ints(values):
    return [int(v) for v in values]


def topology_cases():
    out = []
    for top, traj, groups in TOPOLOGIES:
        u = universe(top, traj)
        atoms = u.atoms
        entry = {"topology": top, "trajectory": traj, "n_atoms": int(atoms.n_atoms),
                 "has_charges": hasattr(atoms, "charges"), "has_chainIDs": hasattr(atoms, "chainIDs"),
                 "charges": attempt(lambda: str(atoms.charges.dtype)),
                 "chainIDs": attempt(lambda: [str(c) for c in atoms.chainIDs]),
                 "masses_dtype": str(atoms.masses.dtype),
                 "positions_dtype": attempt(lambda: str(atoms.positions.dtype)),
                 "bonds": attempt(lambda: pairs(u.bonds)),
                 "bonds_len": attempt(lambda: len(u.bonds)),
                 "fragindices": attempt(lambda: ints(atoms.fragindices)),
                 "n_fragments": attempt(lambda: len(atoms.fragments)),
                 "groups": []}
        rng = np.random.default_rng(len(out))
        for idx in groups:
            ag = atoms[np.array(idx, dtype=np.int64)]
            mask = rng.integers(0, 2, size=len(idx)).astype(bool)
            other = atoms[np.array(idx[::-1][:2] + [atoms.n_atoms - 1], dtype=np.int64)] if idx else atoms[[0]]
            entry["groups"].append({
                "indices": idx,
                "bonds": attempt(lambda: pairs(ag.bonds)),
                "fragindices": attempt(lambda: ints(ag.fragindices)),
                "fragments": attempt(lambda: [ints(f.indices) for f in ag.fragments]),
                "residues": ints(ag.residues.resindices),
                "residue_atoms": ints(ag.residues.atoms.indices),
                "mask": [bool(m) for m in mask], "masked": ints(ag[mask].indices),
                "union_with": ints(other.indices), "union": ints((ag | other).indices),
            })
        out.append(entry)
    return out


def residue_cases():
    out = []
    for top, residues in [("topology/demo_dry.parm7", [[3, 1, 1], [26, 0], [2]]),
                          ("io/loader/files/seven.pdb", [[1, 0], [0, 0]])]:
        u = universe(top)
        for res in residues:
            rg = u.residues[np.array(res, dtype=np.int64)]
            out.append({"topology": top, "residues": res, "atoms": ints(rg.atoms.indices),
                        "atom_resindices": ints(rg.atoms.resindices)})
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    dimensions = write_files()
    data = {"generator": "tests/parity/make_universe_cases.py", "mdanalysis": mda.__version__,
            "iteration": iteration_cases(), "dimensions": dimensions, "topology": topology_cases(),
            "residues": residue_cases()}
    remove_offsets(FILES)
    lines = ["{"]
    keys = list(data)
    for n, key in enumerate(keys):
        comma = "," if n < len(keys) - 1 else ""
        value = data[key]
        if isinstance(value, list):
            lines.append(f" {json.dumps(key)}: [")
            lines.append(",\n".join("  " + json.dumps(v, ensure_ascii=False) for v in value))
            lines.append(" ]" + comma)
        else:
            lines.append(f" {json.dumps(key)}: {json.dumps(value)}{comma}")
    lines.append("}")
    with open(os.path.join(OUT, "universe_cases.json"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print({k: len(v) for k, v in data.items() if isinstance(v, list)})


main()
