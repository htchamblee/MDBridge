#!/usr/bin/env python3
"""The topology MDAnalysis builds from each file, as JSON (tests/data/topology/*.topology.json).

    python3 -I dump_topologies.py OUT_DIR FILE...

Needs MDAnalysis (the version is recorded in each file; the reference tables use 2.11.0). The C++
tests build mdbridge::core::Topology objects from these dumps (tests/core/topology_dump.hpp), so the
selection engine can be checked against MDAnalysis without the file readers, and the readers are
checked against MDAnalysis by comparing what they read with the dump.

Per file: atoms (names, types, resindex and every optional per-atom attribute MDAnalysis has for
the format, null when absent), residues (resids, resnums, resnames, icodes, segindex), segments
(segids), bonds (sorted [i, j] with i < j, null when the format has none) and, when the file has
coordinates, the first frame (positions as float32 values, dimensions or null).
"""

import json
import os
import sys
import warnings

import MDAnalysis as mda
import numpy as np

warnings.filterwarnings("ignore")


def attribute(group, name, convert):
    try:
        values = getattr(group, name)
    except (AttributeError, mda.exceptions.NoDataError):
        return None
    return [convert(v) for v in values]


def dump(path):
    u = mda.Universe(path)
    atoms, residues, segments = u.atoms, u.residues, u.segments
    out = {
        "file": os.path.basename(path),
        "mdanalysis": mda.__version__,
        "atoms": {
            "names": [str(v) for v in atoms.names],
            "types": attribute(atoms, "types", str),
            "resindex": [int(v) for v in atoms.resindices],
            "elements": attribute(atoms, "elements", str),
            "chainIDs": attribute(atoms, "chainIDs", str),
            "masses": attribute(atoms, "masses", float),
            "charges": attribute(atoms, "charges", float),
            "ids": attribute(atoms, "ids", int),
            "altLocs": attribute(atoms, "altLocs", str),
            "record_types": attribute(atoms, "record_types", str),
            "occupancies": attribute(atoms, "occupancies", float),
            "tempfactors": attribute(atoms, "tempfactors", float),
            "formalcharges": attribute(atoms, "formalcharges", int),
        },
        "residues": {
            "resids": [int(v) for v in residues.resids],
            "resnums": attribute(residues, "resnums", int),
            "resnames": [str(v) for v in residues.resnames],
            "icodes": attribute(residues, "icodes", str),
            "segindex": [int(v) for v in residues.segindices],
        },
        "segments": {"segids": [str(v) for v in segments.segids]},
        "bonds": None,
        "frame": None,
    }
    try:
        bonds = np.sort(np.asarray(u.bonds.indices, dtype=np.int64), axis=1)
        out["bonds"] = sorted(map(list, {tuple(int(x) for x in b) for b in bonds}))
    except (AttributeError, mda.exceptions.NoDataError):
        pass
    try:
        ts = u.trajectory.ts
        out["frame"] = {
            "positions": [[float(np.float32(x)) for x in row] for row in u.atoms.positions],
            "dimensions": None if ts.dimensions is None else [float(x) for x in ts.dimensions],
            "time": float(ts.time),
        }
    except (AttributeError, ValueError):
        pass
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: dump_topologies.py OUT_DIR FILE...")
    out_dir = sys.argv[1]
    os.makedirs(out_dir, exist_ok=True)
    for path in sys.argv[2:]:
        data = dump(path)
        target = os.path.join(out_dir, os.path.basename(path) + ".topology.json")
        with open(target, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=1, ensure_ascii=False)
            fh.write("\n")


if __name__ == "__main__":
    main()
