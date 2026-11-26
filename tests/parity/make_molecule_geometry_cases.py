"""What MDATK's molecule geometry (mdatk/render/molecule/topology.py, geometry.py, surface.py,
scene.py) computes, for src/render/molecule (tests/render/test_molecule_geometry.cpp; milestone 4,
chunk W0-F of docs/specs/m4-plan.md).

    /opt/mdbridge-venv/bin/python -I -B make_molecule_geometry_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies (numpy, scipy, MDAnalysis 2.11.0); nothing is written into the checkout.
MDATK runs in-process. The demo system's atoms come from MDATK's own `trace` and `atoms` analyses
on tests/data/topology/demo_dry.parm7 + tests/data/io/trajectories/demo_dry.nc, loaded as
make_analysis_cases.py loads them (MDBridge's independent tables and residue-scoped ion rule,
MDBridge's analysis improvements; tests/data/README.md), with `visual.align_frames = false` so that
the coordinates are the trajectory's float32 values (no BLAS superposition). The other systems are
built here from plain arrays. Writes into OUT_DIR:

  systems.json     the inputs of every Topology: per system the per-atom names, elements, resnames,
                   row, the bonds, n_res, per-residue chain and resindex, groups (names and
                   per-residue membership) and the frames used (float32 coordinates widened to
                   float64, as MDATK's tests use them; NaN where a case blanks an atom)
  topology.json    normalise_name / base_letter of name lists, and per system and backbone atom
                   the Topology roles (normalised names, elements, resnames, letters, radii,
                   element colours, control / anchor / side / C1' / C4' atoms, ring atoms and
                   glycosidic atoms of purines / pyrimidines / sugars, segments, nucleic)
  styles.json      per case (system, frame, style, residue mask, size, sugars, samples, backbone
                   atom) the style's drawn() residues and its DrawList: meshes (verts, normals,
                   tris, owner, role, atom), spheres (xyz, radius, owner, atom), sticks (a, b,
                   radius, owner, atom)
  surface.json     density_grid (shape, origin, float32 values), surface_nets and gaussian_surface
                   (verts, normals, tris, nearest atom, and each vertex's distance to it) of small
                   point sets and of demo residues
  scenes.json      MoleculeScene: specs and masks after fill_rest (residue_masks of the specs'
                   groups), which styles exist, drawn(), and frame(xyz, residue_rgb): every item's
                   kind, size, owners, atoms / roles, colours and opacity (its geometry is the
                   style's, as in styles.json); depth_range along a fixed axis
  ndimage.json     scipy.ndimage.gaussian_filter1d / gaussian_filter (mode "constant") and
                   map_coordinates (order 1, mode "nearest") on float32 arrays, as surface.py calls
                   them (the C++ port of these is in src/numerics/ndimage.cpp)

Numbers are written as repr() of the float (exact). numpy's and OpenBLAS' kernels are pinned before
numpy is imported, as make_analysis_cases.py pins them (the spline weights go through a small
matmul; np.exp makes the Gaussian kernel): the output is the same on every machine tried.
"""

import os
import platform
import sys

if True:  # before numpy is imported (read when OpenBLAS / numpy load)
    os.environ["OPENBLAS_CORETYPE"] = "Nehalem"
    for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_var] = "1"
    if platform.machine().lower() in ("x86_64", "amd64"):
        os.environ["NPY_ENABLE_CPU_FEATURES"] = "X86_V2"

import json  # noqa: E402
import logging  # noqa: E402
import math  # noqa: E402
import warnings  # noqa: E402

if len(sys.argv) != 3:
    sys.exit("usage: make_molecule_geometry_cases.py MDATK_CHECKOUT OUT_DIR")
CHECKOUT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
sys.path.insert(0, CHECKOUT)
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
sys.path.insert(0, HERE)

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

import independent_tables  # noqa: E402

independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)
import improvements_periodic  # noqa: E402

improvements_periodic.apply_mdatk()
import improvements_analyses  # noqa: E402

improvements_analyses.apply_mdatk()  # MDATK with MDBridge's analysis improvements (docs/DIFFERENCES.md)

import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402

from mdatk.analysis import AnalysisContext, run_analyses  # noqa: E402
from mdatk.analysis.base import analysis_class, instance_project  # noqa: E402
from mdatk.config import RepSpec  # noqa: E402
from mdatk.demo import demo_project  # noqa: E402
from mdatk.io import loader  # noqa: E402
from mdatk.progress import Reporter  # noqa: E402
from mdatk.render.molecule import surface as surface_mod  # noqa: E402
from mdatk.render.molecule.geometry import STYLES, Mesh, Spheres  # noqa: E402
from mdatk.render.molecule.scene import MoleculeScene, depth_range, residue_masks  # noqa: E402
from mdatk.render.molecule.topology import Topology, base_letter, normalise_name  # noqa: E402

DEMO_TOP = os.path.join(DATA, "topology", "demo_dry.parm7")
DEMO_NC = os.path.join(DATA, "io", "trajectories", "demo_dry.nc")


def num(x):
    x = float(x)
    if math.isnan(x):
        return float("nan")
    return x


def flat(a):
    """A numeric array as a flat list (floats exact, ints as ints, bools as 0 / 1)."""
    a = np.asarray(a)
    if a.dtype.kind in "iub":
        return [int(v) for v in a.ravel()]
    return [num(v) for v in a.ravel()]


def write(name, data):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(data, ensure_ascii=False, indent=1, allow_nan=True))
        f.write("\n")


def write_cases(name, header, cases):
    """One case per line (compact), so that diffs stay readable."""
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as f:
        f.write('{"note": ' + json.dumps(header) + ',\n"cases": [\n')
        for k, case in enumerate(cases):
            f.write(json.dumps(case, ensure_ascii=False, separators=(",", ":"), allow_nan=True))
            f.write(",\n" if k + 1 < len(cases) else "\n")
        f.write("]}\n")


# --- systems ----------------------------------------------------------------------------------------


class System:
    def __init__(self, name, names, elements, resnames, row, n_res, bonds, chain, resindex, groups,
                 frames):
        self.name = name
        self.names = list(names)
        self.elements = list(elements)
        self.resnames = list(resnames)
        self.row = np.asarray(row, np.int64)
        self.n_res = int(n_res)
        self.bonds = np.asarray(bonds, np.int64).reshape(-1, 2)
        self.chain = np.asarray(chain, np.int64)
        self.resindex = np.asarray(resindex, np.int64)
        self.groups = groups  # [(name, per-residue bool array)]
        self.frames = {k: np.asarray(v, np.float64) for k, v in frames.items()}

    def topology(self, backbone_atom="P"):
        return Topology(self.names, self.elements, self.resnames, self.row, self.n_res, self.bonds,
                        self.frames[0], self.chain, self.resindex, backbone_atom)

    def json(self):
        return {"names": self.names, "elements": self.elements, "resnames": self.resnames,
                "row": flat(self.row), "n_res": self.n_res, "bonds": flat(self.bonds),
                "chain": flat(self.chain), "resindex": flat(self.resindex),
                "groups": [{"name": n, "member": flat(m)} for n, m in self.groups],
                "frames": {str(k): flat(v) for k, v in sorted(self.frames.items())}}


def demo_system():
    p = demo_project(DEMO_TOP, DEMO_NC, os.path.join(OUT, "unused"), 12, cleaning=False, mid_time_ps=120.0)
    p.visual.align_frames = False  # raw float32 trajectory coordinates (no BLAS superposition)
    u = loader.load_universe(DEMO_TOP, [DEMO_NC], None)
    independent_tables.residue_rule_universe(u, DEMO_TOP, topology_format="PARM7")
    reporter = Reporter()
    objs = [analysis_class(n)(AnalysisContext(u, instance_project(p, n), DEMO_TOP, reporter), n)
            for n in ("trace", "atoms")]
    (tr, a), _times = run_analyses(objs[0].ctx, objs, np.arange(u.trajectory.n_frames), reporter)
    coords = np.asarray(a.arrays["coords"], dtype=np.float64)
    n_res = len(tr.arrays["resid"])
    member = np.asarray(tr.arrays["group_member"], bool)
    groups = [(name, member[k]) for k, name in enumerate(p.group_names())]
    frames = {k: coords[k] for k in (0, 3, 5, 17)}
    blank = coords[0].copy()
    blank[[5, 40, 41]] = np.nan  # a blown-up frame: non-finite atoms
    frames[100] = blank
    return System("demo_dry", a.meta["atom_names"], a.meta["elements"], a.meta["resnames"], a.arrays["row"],
                  n_res, a.arrays["bonds"], tr.arrays["chain"], tr.arrays["resindex"], groups, frames)


def mixed_system(demo):
    """A protein chain (with a numbering gap and a far jump), four demo nucleotides, a ligand and
    an Mg ion: protein tubes, side atoms CB / O, segment breaks, ligands and ions for fill_rest."""
    names, elements, resnames, row, xyz, bonds = [], [], [], [], [], []
    chain, resindex = [], []

    def atom(name, element, resname, r, p):
        names.append(name)
        elements.append(element)
        resnames.append(resname)
        row.append(r)
        xyz.append(p)
        return len(names) - 1

    protein = ["ALA", "GLY", "SER", "ALA", "LEU", "ALA", "GLY", "ALA"]
    prev_c = None
    for i, rn in enumerate(protein):
        ang = math.radians(100.0 * i)
        shift = 30.0 if i >= 6 else 0.0  # residues 6, 7: more than 8.5 A away (a new segment)
        ca = np.array([40.0 + 2.3 * math.cos(ang) + shift, 2.3 * math.sin(ang), 1.5 * i])
        out = np.array([math.cos(ang), math.sin(ang), 0.0])
        n = atom("N", "N", rn, i, ca + np.array([-0.5, 0.9, -0.6]))
        c_a = atom("CA", "C", rn, i, ca)
        c = atom("C", "C", rn, i, ca + np.array([0.6, -0.4, 0.9]))
        o = atom("O", "O", rn, i, ca + np.array([1.4, -0.9, 1.2]))
        bonds += [(n, c_a), (c_a, c), (c, o)]
        if rn != "GLY":
            cb = atom("CB", "C", rn, i, ca + 1.5 * out)
            bonds.append((c_a, cb))
        if prev_c is not None:
            bonds.append((prev_c, n))
        prev_c = c
        chain.append(0)
        resindex.append(i if i < 4 else i + 1)  # a numbering gap between residues 3 and 4
    base = len(protein)
    d_names, d_rows = demo.names, demo.row
    picked = [k for k in range(len(d_names)) if d_rows[k] < 4]
    local = {}
    for k in picked:
        local[k] = atom(d_names[k], demo.elements[k], demo.resnames[k], base + int(d_rows[k]), demo.frames[0][k])
    for i, j in demo.bonds:
        if int(i) in local and int(j) in local:
            bonds.append((local[int(i)], local[int(j)]))
    for r in range(4):
        chain.append(1)
        resindex.append(20 + r)
    lig = base + 4
    l1 = atom("C1", "C", "LIG", lig, np.array([35.0, 6.0, 4.0]))
    l2 = atom("O2", "O", "LIG", lig, np.array([36.2, 6.4, 4.3]))
    l3 = atom("N3", "N", "LIG", lig, np.array([34.1, 7.1, 4.4]))
    bonds += [(l1, l2), (l1, l3)]
    chain.append(2)
    resindex.append(30)
    atom("MG", "MG", "MG", lig + 1, np.array([33.0, 3.0, 2.0]))
    chain.append(3)
    resindex.append(31)
    n_res = lig + 2
    groups = [("Protein", np.arange(n_res) < base), ("RNA", (np.arange(n_res) >= base) & (np.arange(n_res) < lig)),
              ("Ligand", np.arange(n_res) >= lig)]
    xyz = np.array(xyz, np.float64)
    moved = xyz.copy()
    moved += 0.25 * np.sin(np.arange(xyz.size).reshape(xyz.shape))  # a second, distorted frame
    return System("mixed", names, elements, resnames, row, n_res, bonds, chain, resindex, groups,
                  {0: xyz, 1: moved})


def charmm_system():
    """MDATK's test_charmm_names_and_bases_told_by_their_atoms: bases told by their atoms."""
    from mdatk.render.molecule.topology import PURINE_RING, PYRIMIDINE_RING

    names, res = [], []
    for r, extra in enumerate([("O6",), ("N6",), ("N4",), ("O4",), ("O4", "C7"), ("O4", "C5M"), ()]):
        ring = PURINE_RING if r < 2 else PYRIMIDINE_RING
        for n in ring + extra:
            names.append(n)
            res.append(r)
    n = len(names)
    xyz = np.zeros((n, 3))
    return System("charmm", names, ["C"] * n, ["XXX"] * n, res, 7, np.zeros((0, 2)), np.zeros(7, int),
                  np.arange(7), [], {0: xyz})


def odd_names_system():
    """Old atom names (C1*, O1P), lower-case and blank elements, pseudouridine (C-glycosidic),
    residue names of every naming scheme, a residue without atoms and an ion."""
    rows = [
        ("PSU", ["P", "O1P", "C4*", "C1*", "C2*", "C3*", "O4*", "N1", "C2", "N3", "C4", "C5", "C6", "O4"]),
        ("RG5", ["C4'", "C1'", "C2'", "C3'", "O4'", "N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "O6"]),
        ("DT", ["P", "C4'", "C1'", "N1", "C2", "N3", "C4", "C5", "C6", "C7"]),
        ("CYT", ["P", "C4'", "C1'"]),
        ("ADE", ["P"]),
        ("XYZ", ["C1", "C2"]),
    ]
    names, elements, resnames, row, xyz = [], [], [], [], []
    for r, (rn, atoms) in enumerate(rows):
        for k, n in enumerate(atoms):
            names.append(" " + n + " " if k % 3 == 0 else n)
            el = n[0]
            elements.append(el.lower() if k % 4 == 1 else (" " if k % 5 == 2 else el))
            resnames.append(rn)
            row.append(r)
            xyz.append([3.0 * r + 0.3 * k, 0.5 * k, 0.1 * r])
    # a residue (6) without atoms, then a zinc ion (7)
    names.append("ZN")
    elements.append("ZN")
    resnames.append("ZN")
    row.append(7)
    xyz.append([0.0, 0.0, 9.0])
    n_res = 8
    return System("odd_names", names, elements, resnames, row, n_res, np.zeros((0, 2)), np.zeros(n_res, int),
                  np.arange(n_res), [], {0: np.array(xyz)})


# --- dumps ------------------------------------------------------------------------------------------


def topology_json(topo):
    return {
        "names": topo.names, "elements": topo.elements, "resname": topo.resname, "letter": topo.letter,
        "radius": flat(topo.radius), "element_rgb": flat(topo.element_rgb), "is_carbon": flat(topo.is_carbon),
        "is_hydrogen": flat(topo.is_hydrogen), "control": flat(topo.control), "anchor_atom": flat(topo.anchor_atom),
        "side": flat(topo.side), "c1": flat(topo.c1), "c4": flat(topo.c4),
        "purine_res": flat(topo.purine_res), "purine_ring": flat(topo.purine_ring), "purine_gly": flat(topo.purine_gly),
        "pyrimidine_res": flat(topo.pyrimidine_res), "pyrimidine_ring": flat(topo.pyrimidine_ring),
        "pyrimidine_gly": flat(topo.pyrimidine_gly), "sugar_res": flat(topo.sugar_res),
        "sugar_ring": flat(topo.sugar_ring), "segments": [flat(s) for s in topo.segments],
        "nucleic": flat(topo.nucleic)}


def prim_json(p):
    if isinstance(p, Mesh):
        out = {"kind": "mesh", "verts": flat(p.verts), "normals": flat(p.normals), "tris": flat(p.tris),
               "owner": flat(p.owner), "role": flat(p.role)}
        if p.atom is not None:
            out["atom"] = flat(p.atom)
        return out
    if isinstance(p, Spheres):
        return {"kind": "spheres", "xyz": flat(np.asarray(p.xyz).reshape(-1, 3)), "radius": flat(p.radius),
                "owner": flat(p.owner), "atom": flat(p.atom)}
    return {"kind": "sticks", "a": flat(np.asarray(p.a).reshape(-1, 3)), "b": flat(np.asarray(p.b).reshape(-1, 3)),
            "radius": flat(p.radius), "owner": flat(p.owner), "atom": flat(p.atom)}


def item_json(item):
    """A scene item: its primitive's kind, size, owners and atoms / roles (the geometry is that of
    the style, compared in styles.json), colours and opacity."""
    p = item.prim
    kind = "mesh" if isinstance(p, Mesh) else "spheres" if isinstance(p, Spheres) else "sticks"
    out = {"kind": kind, "owner": flat(p.owner)}
    if kind == "mesh":
        out["n_tris"] = len(p.tris)
        out["role"] = flat(p.role)
        if p.atom is not None:
            out["atom"] = flat(p.atom)
    else:
        out["atom"] = flat(p.atom)
        out["radius"] = flat(p.radius)
    out["rgb"] = flat(item.rgb)
    out["opacity"] = item.opacity
    return out


def drawlist_json(dl):
    return {"meshes": [prim_json(m) for m in dl.meshes], "spheres": [prim_json(s) for s in dl.spheres],
            "sticks": [prim_json(s) for s in dl.sticks]}


def mask_of(n_res, rows):
    m = np.zeros(n_res, bool)
    m[list(rows)] = True
    return m


def style_cases(systems):
    demo, mixed = systems["demo_dry"], systems["mixed"]
    first = range(12)
    every = range(demo.n_res)
    specs = []
    for style in sorted(STYLES):  # MDATK's test_styles_build_outward_meshes_for_selected_residues
        specs.append(("demo_dry", 3, style, first, {"sugars": True}, "P"))
    for style in ("trace", "vdw"):
        specs.append(("demo_dry", 0, style, every, {}, "P"))
    specs.append(("demo_dry", 0, "tube", every, {}, "P"))
    specs += [
        ("demo_dry", 0, "cartoon", range(4), {"sugars": False}, "P"),  # test_cartoon_rungs_are_closed_by_balls
        ("demo_dry", 0, "cartoon", range(4), {"sugars": True}, "P"),
        ("demo_dry", 0, "tube", [0, 2, 3], {}, "P"),                  # test_undrawable_residues_fall_back
        ("demo_dry", 5, "cartoon", [0, 5, 6, 7, 20], {"size": 1.3}, "P"),
        ("demo_dry", 5, "cartoon", range(12, 18), {"sugars": True}, "C4'"),
        ("demo_dry", 5, "ribbon", [1, 2, 3, 13, 14], {"samples": 5}, "C4'"),
        ("demo_dry", 17, "tube", [11, 12], {"size": 0.7}, "P"),
        ("demo_dry", 17, "trace", [0, 11, 12, 23], {"size": 1.4}, "C4'"),
        ("demo_dry", 17, "licorice", [0, 1, 23], {"size": 0.5}, "P"),
        ("demo_dry", 17, "vdw", [7], {"size": 0.8}, "P"),
        ("demo_dry", 17, "surf", [7, 8], {"size": 1.2}, "P"),
        ("demo_dry", 100, "surf", range(4), {}, "P"),    # non-finite atoms are left out of surfaces
        ("demo_dry", 100, "licorice", [0, 1, 2, 3], {}, "P"),
        ("demo_dry", 100, "vdw", [0, 1, 2], {}, "P"),
    ]
    for style in sorted(STYLES):
        specs.append(("mixed", 0, style, range(mixed.n_res), {}, "P"))
        if style != "surf":
            specs.append(("mixed", 1, style, range(mixed.n_res), {"sugars": True}, "C4'"))
    specs.append(("mixed", 1, "surf", range(8), {"size": 0.8}, "C4'"))
    specs += [
        ("mixed", 0, "tube", [2, 3, 4, 5, 9], {}, "P"),
        ("mixed", 0, "cartoon", [0, 1, 8, 9, 10, 11], {}, "P"),
        ("mixed", 1, "trace", [6, 7, 12, 13], {}, "P"),
        ("mixed", 1, "licorice", [12, 13], {}, "P"),
        ("odd_names", 0, "cartoon", range(8), {"sugars": True}, "P"),
        ("odd_names", 0, "cartoon", range(8), {"sugars": True}, "C4'"),
        ("odd_names", 0, "licorice", range(8), {}, "P"),
        ("odd_names", 0, "vdw", range(8), {}, "P"),
        ("odd_names", 0, "trace", range(8), {}, "C4'"),
    ]
    cases = []
    topos = {}
    for system, frame, style, rows, opts, backbone in specs:
        s = systems[system]
        key = (system, backbone)
        if key not in topos:
            topos[key] = s.topology(backbone)
        topo = topos[key]
        mask = mask_of(s.n_res, rows)
        size = opts.get("size", 1.0)
        kwargs = {"sugars": opts.get("sugars", False)}
        if "samples" in opts:
            kwargs["samples"] = opts["samples"]
        st = STYLES[style](topo, mask, size=size, **kwargs)
        dl = st.build(s.frames[frame])
        cases.append({"system": system, "frame": frame, "style": style, "mask": flat(mask), "size": size,
                      "sugars": kwargs["sugars"], "samples": opts.get("samples", 8), "backbone_atom": backbone,
                      "drawn": flat(st.drawn()), "drawlist": drawlist_json(dl)})
    return cases


def surface_cases(systems):
    demo = systems["demo_dry"]
    sets = [
        ("one atom", np.array([[0.3, -1.2, 2.0]]), 1.0, 1.4, 0.5),
        ("two atoms", np.array([[0.0, 0.0, 0.0], [2.63, 0.4, -0.2]]), 1.0, 1.4, 0.5),
        ("two atoms, coarse", np.array([[0.0, 0.0, 0.0], [2.63, 0.4, -0.2]]), 1.5, 1.4, 0.5),
        ("three atoms, wide", np.array([[1.0, 2.0, 3.0], [3.1, 2.2, 2.9], [2.0, 4.1, 3.3]]), 1.0, 2.1, 0.5),
        ("three atoms, level 0.3", np.array([[1.0, 2.0, 3.0], [3.1, 2.2, 2.9], [2.0, 4.1, 3.3]]), 1.0, 1.4, 0.3),
        ("demo residue 0", demo.frames[0][demo.row == 0], 1.0, 1.4, 0.5),
        ("demo residues 0-5, coarse", demo.frames[0][demo.row < 6], 1.5, 1.4, 0.5),
    ]
    cases = []
    for name, xyz, spacing, sigma, level in sets:
        case = {"name": name, "xyz": flat(xyz), "spacing": spacing, "sigma": sigma, "level": level}
        vol, lo = surface_mod.density_grid(xyz, spacing, sigma)
        case["grid_shape"] = list(vol.shape)
        case["grid_origin"] = flat(lo)
        if vol.size <= 12000:
            case["grid"] = flat(vol)
        verts, tris, coords = surface_mod.surface_nets(vol, level, lo, spacing)
        case["nets"] = {"verts": flat(verts), "tris": flat(tris), "coords": flat(coords)}
        verts, nrm, tris, near = surface_mod.gaussian_surface(xyz, spacing, sigma, level)
        dist = np.linalg.norm(verts - xyz[near], axis=1) if len(verts) else np.zeros(0)
        case["surface"] = {"verts": flat(verts), "normals": flat(nrm), "tris": flat(tris), "nearest": flat(near),
                           "nearest_distance": flat(dist)}
        cases.append(case)
    return cases


def scene_cases(systems):
    specs = [
        ("demo_dry", 5, [RepSpec("cartoon", ["Aptamer"]), RepSpec("licorice", ["Target"], "element")], "trace"),
        ("demo_dry", 5, [RepSpec("tube", ["Aptamer"], "residue", 1.2, 0.6)], "trace"),
        ("demo_dry", 5, [RepSpec("tube", ["Aptamer"], "residue", 1.2, 0.6)], "hide"),
        ("mixed", 0, [RepSpec("cartoon", [], "auto", 1.0, 1.0, True)], "trace"),
        ("mixed", 0, [RepSpec("cartoon", [], "residue", 1.0, 1.0, True)], "trace"),
        ("mixed", 0, [RepSpec("cartoon", [], "nucleotide", 1.0, 1.0, True)], "trace"),
        ("mixed", 0, [RepSpec("cartoon", [], "element", 1.0, 1.0, True)], "trace"),
        ("demo_dry", 0, [RepSpec("vdw", ["Aptamer"], "nucleotide"), RepSpec("licorice", [], "auto", 0.6)], "trace"),
        ("mixed", 0, [RepSpec("surf", ["Ligand"], "element", 1.0, 0.5), RepSpec("surf", ["Protein"], "residue")],
         "hide"),
        ("demo_dry", 0, [RepSpec("bogus"), RepSpec("ribbon", ["Nope"]), RepSpec("trace", ["Target"], "auto", 0.5)],
         "trace"),
        ("demo_dry", 0, [], "trace"),
        ("mixed", 0, [RepSpec("cartoon", ["RNA"], "auto", 1.0, 1.0, True)], "trace"),
        ("mixed", 0, [RepSpec("ribbon", [], "auto"), RepSpec("licorice", ["Ligand"], "element")], "trace"),
        ("mixed", 1, [RepSpec("licorice", ["Protein"], "auto"), RepSpec("vdw", ["Ligand"], "residue")], "trace"),
        ("mixed", 1, [RepSpec("surf", ["Protein", "Ligand"], "auto", 0.9)], "trace"),
        ("odd_names", 0, [RepSpec("cartoon", [], "auto", 1.0, 1.0, True)], "trace"),
    ]
    axis = np.array([0.3, -0.4, math.sqrt(0.75)])
    axis = axis / np.linalg.norm(axis)
    cases = []
    for system, frame, reps, rest, in specs:
        s = systems[system]
        topo = s.topology("P")
        names = [n for n, _m in s.groups]
        member = np.array([m for _n, m in s.groups], bool).reshape(-1, s.n_res)
        masks = residue_masks(reps, member, names, s.n_res)
        scene = MoleculeScene(topo, reps, masks)
        if rest == "trace":
            scene.fill_rest()
        res = np.arange(s.n_res)
        residue_rgb = np.stack([(res * 37 % 100) / 100.0, (res * 61 % 100) / 100.0, (res * 17 % 100) / 100.0], 1)
        items = scene.frame(s.frames[frame], residue_rgb)
        lo, hi = depth_range(items, axis)
        cases.append({
            "system": system, "frame": frame, "rest": rest,
            "specs": [{"style": r.style, "groups": list(r.groups), "color": r.color, "size": r.size,
                       "opacity": r.opacity, "sugars": r.sugars} for r in reps],
            "residue_rgb": flat(residue_rgb),
            "scene_specs": [{"style": r.style, "color": r.color, "size": r.size, "opacity": r.opacity,
                             "sugars": r.sugars} for r in scene.specs],
            "scene_masks": [flat(st.mask) if st is not None else None for st in scene.styles],
            "drawn": flat(scene.drawn()),
            "items": [item_json(it) for it in items],
            "depth_axis": flat(axis), "depth_range": [lo, hi]})
    return cases


def ndimage_cases():
    rng = np.random.default_rng(20261125)
    cases = {"filter1d": [], "filter3d": [], "map_coordinates": []}
    for n, sigma, truncate in ((1, 1.4, 3.0), (5, 1.4, 3.0), (20, 1.4, 3.0), (20, 0.9333333333333333, 3.0),
                               (31, 2.5, 4.0), (12, 0.2, 3.0), (40, 6.0, 3.0)):
        x = rng.random(n).astype(np.float32)
        y = ndimage.gaussian_filter1d(x, sigma, truncate=truncate, mode="constant")
        cases["filter1d"].append({"input": flat(x), "sigma": sigma, "truncate": truncate, "output": flat(y),
                                  "kernel": flat(ndimage._filters._gaussian_kernel1d(sigma, 0, int(truncate * sigma + 0.5)))})
    for shape, sigma in (((6, 7, 8), 1.4), ((3, 9, 5), 0.9333333333333333), ((1, 4, 4), 1.0)):
        x = rng.random(shape).astype(np.float32)
        y = ndimage.gaussian_filter(x, sigma, truncate=3.0, mode="constant")
        cases["filter3d"].append({"shape": list(shape), "input": flat(x), "sigma": sigma, "output": flat(y)})
    vol = rng.random((4, 5, 6)).astype(np.float32)
    pts = np.concatenate([rng.random((40, 3)) * np.array([4, 5, 6]) - 0.5,
                          np.array([[0, 0, 0], [3, 4, 5], [-2.0, 1.5, 7.2], [3.999, 4.0, 5.0], [1.5, 2.5, 3.5]])])
    out = ndimage.map_coordinates(vol, pts.T, order=1, mode="nearest")
    cases["map_coordinates"].append({"shape": list(vol.shape), "input": flat(vol), "points": flat(pts),
                                     "output": flat(out)})
    return cases


def main():
    demo = demo_system()
    systems = {"demo_dry": demo, "mixed": mixed_system(demo), "charmm": charmm_system(),
               "odd_names": odd_names_system()}
    write("systems.json", {name: s.json() for name, s in systems.items()})
    names = ["C1*", "O1P", "O2P", "O3P", " C4' ", "P", "OP1", "H5''", "C2*'", ""]
    resnames = ["G", "RA5", "C3", "DT", "PSU", "U5", "MG", "LYS", "ADE", "GUA", "CYT", "URA", "THY", "URI",
                "ADE5", "GUAX", "RU", "RC3", "DA5", "DG3", "LA", "LGX", "A3", "AN", "A5X", "C", "T", "U",
                "PU", " g ", "ra", "URIDX", "THYM", "", "A55", "TN", "X"]
    topology = {"normalise_name": [[n, normalise_name(n)] for n in names],
                "base_letter": [[r, base_letter(r)] for r in resnames],
                "topologies": []}
    for name in ("demo_dry", "mixed", "charmm", "odd_names"):
        for backbone in ("P", "C4'", "CA"):
            topology["topologies"].append({"system": name, "backbone_atom": backbone,
                                           "topology": topology_json(systems[name].topology(backbone))})
    write("topology.json", topology)
    write_cases("styles.json", "STYLES[style](topology, mask, size, sugars, samples).build(frame); "
                "make_molecule_geometry_cases.py", style_cases(systems))
    write_cases("surface.json", "mdatk/render/molecule/surface.py on small point sets; "
                "make_molecule_geometry_cases.py", surface_cases(systems))
    write_cases("scenes.json", "MoleculeScene(...).frame(xyz, residue_rgb); make_molecule_geometry_cases.py",
                scene_cases(systems))
    write("ndimage.json", ndimage_cases())


if __name__ == "__main__":
    main()
