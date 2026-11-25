"""What MDATK's 3D representations compute (mdatk/render/molecule: topology.py, geometry.py,
surface.py, scene.py, software.py), for MDBridge's molecule geometry and software renderer
(src/render molecule/ and software/; docs/specs/m4-plan.md sections 1.5 and 3).

    /opt/mdbridge-venv/bin/python -I -B make_molecule_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies; nothing is written into the checkout. MDATK runs in-process, patched by
render_common.setup(), with numpy's and OpenBLAS' kernels pinned (the renderer's float32 arithmetic
and scipy.ndimage's kernels are then the same on every x86-64 machine). Reads the input datasets
make_figure_cases.py wrote (OUT_DIR/inputs: demo_dry, protein_rna, example_example_rep; run it
first). Writes OUT_DIR/molecule/:

  molecule_cases.json
    topology   per system (demo_dry; protein_rna: RNA and protein; schematic: the example
               dataset's schematic duplex), Topology(...) as build_scene makes it (backbone_atom
               "P", and "C4'" for demo_dry_c4): n_res, names (normalised), elements, resname,
               letter, control, anchor_atom, side, c1, c4, purine_res / ring / gly, pyrimidine_res /
               ring / gly, sugar_res / ring, segments, nucleic, radius, is_carbon, is_hydrogen,
               element_rgb
    scenes     per scene (a system, representations, rep_rest, frame): specs (RepSpec fields of
               every representation, the filled-in rest included), masks / drawn (residue rows per
               representation), centre (subtracted from the atoms, as the structure panel does:
               the first frame's trace centroid), residue_rgb (the colours passed to frame(): group
               colours) and items: per Colored item its kind (mesh / spheres / sticks), opacity and
               count; the arrays are in scenes.npz as "<scene>/<k>/<field>"
    renders    per render: what was rendered (a scene's items, or synthetic items as in MDATK's
               tests/test_molecule_render.py, stored in renders.npz as "<render>/<k>/<field>"),
               the viewport (rotation, xlim, ylim, width, height, zrange), the shading (fog_rgb,
               fog_strength, glow per residue row or null, glow_rgb), quality; the image is
               renders.npz "<render>/rgba" (uint8, straight alpha, row 0 at the top), and the
               render's coverage (pixels with alpha > 0)
  scenes.npz  item arrays: mesh verts / normals (float32), tris (int32), owner (int32), role (int8),
              atom (int32, surfaces only); spheres xyz (float32), radius (float32), owner, atom
              (int32); sticks a, b, radius, owner, atom; rgb (float32) per vertex / primitive
  renders.npz the same item fields for the renders' inputs, and their images

Geometry is stored as float32 (the plan's tolerance is 1e-5 Å); every render is made from the
items exactly as stored (float32 values widened back to float64), so the C++ renderer can be
given the same input. Image parity of the software renderer: >= 99.5 % of pixels equal, all
within 2/255 (docs/specs/m4-plan.md section 3).
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import render_common as rc  # noqa: E402  (before numpy: pins its kernels)

CHECKOUT, OUT = rc.args("make_molecule_cases.py MDATK_CHECKOUT OUT_DIR")
rc.setup(CHECKOUT)
rc.add_mask(OUT, "{data}/render")

import dataclasses  # noqa: E402
import shutil  # noqa: E402

import numpy as np  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

import improvements_render  # noqa: E402
from mdatk.config import RepSpec  # noqa: E402
from mdatk.render.coloring import static_colors  # noqa: E402
from mdatk.render.molecule.geometry import Mesh, Spheres, Sticks  # noqa: E402
from mdatk.render.molecule.scene import Colored, Shading, build_scene  # noqa: E402
from mdatk.render.molecule.software import SoftwareRenderer  # noqa: E402
from mdatk.render.molecule.stage import Viewport  # noqa: E402
from mdatk.render.molecule.surface import gaussian_surface  # noqa: E402
from mdatk.render.molecule.topology import Topology  # noqa: E402
from mdatk.render.orientation import euler_to_matrix  # noqa: E402
from mdatk.render.theme import LIGHT  # noqa: E402

INPUTS = os.path.join(OUT, "inputs")
MOLECULE = os.path.join(OUT, "molecule")
WHITE = np.ones(3)
PLAIN = Shading(WHITE, 0.0)
GLOW_RGB = np.array([1.0, 0.69, 0.0])
LIGHT_SURFACE = np.array(to_rgb(LIGHT.surface))  # the fog colour of the light theme

SYSTEMS = {  # system -> (input dataset, group names and colours)
    "demo_dry": ("demo_dry", [("Aptamer", "#2a78d6"), ("Target", "#eb6834")]),
    "protein_rna": ("protein_rna", [("RNA", "#2a78d6"), ("Protein", "#eb6834")]),
    "schematic": ("example_example_rep", [("Aptamer", "#2a78d6"), ("Target", "#eb6834")]),
}


# --- items ------------------------------------------------------------------------------------------

F32 = {"verts", "normals", "xyz", "radius", "a", "b", "rgb"}
I32 = {"tris", "owner", "atom"}


def item_fields(it):
    p = it.prim
    if isinstance(p, Mesh):
        kind, fields = "mesh", {"verts": p.verts, "normals": p.normals, "tris": p.tris, "owner": p.owner,
                                "role": p.role}
        if p.atom is not None:
            fields["atom"] = p.atom
    elif isinstance(p, Spheres):
        kind, fields = "spheres", {"xyz": p.xyz, "radius": p.radius, "owner": p.owner, "atom": p.atom}
    else:
        kind, fields = "sticks", {"a": p.a, "b": p.b, "radius": p.radius, "owner": p.owner, "atom": p.atom}
    fields["rgb"] = it.rgb
    out = {}
    for k, v in fields.items():
        v = np.asarray(v)
        out[k] = v.astype(np.float32) if k in F32 else v.astype(np.int8) if k == "role" else v.astype(np.int32)
    return kind, out


def stored(items, prefix, arrays):
    """The items as stored (float32 / int32), written to `arrays`; returns (records, items rebuilt
    from the stored values)."""
    records, rebuilt = [], []
    for k, it in enumerate(items):
        kind, fields = item_fields(it)
        for name, v in fields.items():
            arrays[f"{prefix}/{k}/{name}"] = v
        f = {n: (v.astype(np.float64) if v.dtype == np.float32 else v.astype(np.int64)) for n, v in fields.items()}
        if kind == "mesh":
            prim = Mesh(f["verts"], f["normals"], f["tris"], f["owner"], fields["role"], f.get("atom"))
            count = len(f["tris"])
        elif kind == "spheres":
            prim = Spheres(f["xyz"], f["radius"], f["owner"], f["atom"])
            count = len(f["xyz"])
        else:
            prim = Sticks(f["a"], f["b"], f["radius"], f["owner"], f["atom"])
            count = len(f["a"])
        records.append({"kind": kind, "opacity": float(it.opacity), "count": count,
                        "fields": list(fields)})
        rebuilt.append(Colored(prim, f["rgb"], float(it.opacity)))
    return records, rebuilt


# --- systems ----------------------------------------------------------------------------------------

def system(name):
    inp, groups = SYSTEMS[name]
    ds = rc.load_dataset(os.path.join(INPUTS, inp))
    return ds.results["atoms"], ds.results["trace"], groups


def topology(atoms, trace, backbone_atom="P"):
    n_res = len(trace.arrays["resid"])
    return Topology(atoms.meta["atom_names"], atoms.meta["elements"], atoms.meta["resnames"], atoms.arrays["row"],
                    n_res, atoms.arrays["bonds"], np.asarray(atoms.arrays["coords"][0], dtype=np.float64),
                    trace.arrays["chain"], trace.arrays.get("resindex", np.arange(n_res)), backbone_atom)


def topology_record(t):
    keys = ["n_res", "names", "elements", "resname", "letter", "control", "anchor_atom", "side", "c1", "c4",
            "purine_res", "purine_ring", "purine_gly", "pyrimidine_res", "pyrimidine_ring", "pyrimidine_gly",
            "sugar_res", "sugar_ring", "radius", "is_carbon", "is_hydrogen", "element_rgb"]
    out = {k: getattr(t, k) for k in keys}
    out["segments"] = [s.tolist() for s in t.segments]
    out["nucleic"] = t.nucleic
    return out


def residue_rgb(trace, groups):
    return static_colors("groups", np.asarray(trace.arrays["chain"]), np.asarray(trace.arrays["group_index"]),
                         [c for _n, c in groups], "#b5b3ad", "tab10")[:, :3].astype(np.float64)


SCENES = [
    # (name, system, representations, rep_rest, frame, backbone_atom)
    ("demo_cartoon", "demo_dry", [RepSpec("cartoon")], "trace", 0, "P"),
    ("demo_cartoon_f12", "demo_dry", [RepSpec("cartoon")], "trace", 12, "P"),
    ("demo_cartoon_c4", "demo_dry", [RepSpec("cartoon")], "trace", 0, "C4'"),
    ("demo_cartoon_sugars_one", "demo_dry", [RepSpec("cartoon", ["Aptamer"], sugars=True, size=1.4)], "trace", 0, "P"),
    ("demo_ribbon", "demo_dry", [RepSpec("ribbon", color="nucleotide")], "trace", 0, "P"),
    ("demo_tube", "demo_dry", [RepSpec("tube", ["Target"], color="residue")], "hide", 0, "P"),
    ("demo_licorice", "demo_dry", [RepSpec("licorice", ["Aptamer"], color="element")], "hide", 0, "P"),
    ("demo_licorice_auto", "demo_dry", [RepSpec("licorice")], "trace", 0, "P"),
    ("demo_vdw", "demo_dry", [RepSpec("vdw", size=0.8)], "trace", 0, "P"),
    ("demo_surf", "demo_dry", [RepSpec("surf", ["Target"], color="element", opacity=0.5)], "trace", 0, "P"),
    ("demo_surf_residue", "demo_dry", [RepSpec("surf", color="residue")], "trace", 0, "P"),
    ("demo_trace", "demo_dry", [RepSpec("trace", size=1.5)], "trace", 0, "P"),
    ("demo_mixed", "demo_dry", [RepSpec("cartoon", ["Aptamer"]), RepSpec("licorice", ["Target"]),
                                RepSpec("vdw", ["Target"], opacity=0.4)], "trace", 0, "P"),
    ("protein_tube", "protein_rna", [RepSpec("tube")], "trace", 0, "P"),
    ("protein_cartoon", "protein_rna", [RepSpec("cartoon")], "trace", 1, "P"),
    ("protein_licorice_rest", "protein_rna", [RepSpec("licorice", ["Protein"])], "trace", 0, "P"),
    ("schematic_cartoon", "schematic", [RepSpec("cartoon", sugars=True)], "trace", 0, "P"),
]
SCENE_RENDERS = {  # scene -> list of (render name, rotation, quality, fog, glow every n-th residue)
    "demo_cartoon": [("demo_cartoon", [-90.0, 0.0, 0.0], "high", 0.45, 0),
                     ("demo_cartoon_glow", [-90.0, 0.0, 0.0], "high", 0.45, 5),
                     ("demo_cartoon_fast", [-90.0, 0.0, 0.0], "fast", 0.45, 0),
                     ("demo_cartoon_side", [0.0, 60.0, 30.0], "high", 0.0, 0)],
    "demo_cartoon_sugars_one": [("demo_cartoon_sugars_one", [-90.0, 0.0, 0.0], "high", 0.45, 0)],
    "demo_ribbon": [("demo_ribbon", [-90.0, 0.0, 0.0], "high", 0.45, 3)],
    "demo_tube": [("demo_tube", [-90.0, 0.0, 0.0], "high", 0.45, 0)],
    "demo_licorice": [("demo_licorice", [-90.0, 0.0, 0.0], "high", 0.45, 0)],
    "demo_vdw": [("demo_vdw", [-90.0, 0.0, 0.0], "high", 0.45, 4)],
    "demo_surf": [("demo_surf", [-90.0, 0.0, 0.0], "high", 0.45, 0)],
    "demo_surf_residue": [("demo_surf_residue", [-90.0, 0.0, 0.0], "high", 0.45, 6)],
    "demo_trace": [("demo_trace", [-90.0, 0.0, 0.0], "high", 0.45, 0)],
    "demo_mixed": [("demo_mixed", [-90.0, 0.0, 0.0], "high", 0.45, 2)],
    "protein_tube": [("protein_tube", [0.0, 0.0, 0.0], "high", 0.45, 0)],
    "protein_cartoon": [("protein_cartoon", [20.0, -30.0, 0.0], "high", 0.45, 3)],
    "schematic_cartoon": [("schematic_cartoon", [-90.0, 0.0, 0.0], "high", 0.45, 0)],
}
SIZE = 128  # pixels of the scene renders (square)


def spec_record(spec):
    return dataclasses.asdict(spec)


def scene_cases(arrays, render_arrays):
    scenes, renders = {}, {}
    for name, sys_name, specs, rest, frame, backbone in SCENES:
        atoms, trace, groups = system(sys_name)
        names = [n for n, _c in groups]
        sc = build_scene(atoms, trace, specs, names, backbone, rest)
        coords = np.asarray(trace.arrays["coords"], dtype=np.float64)
        centre = coords[0].mean(axis=0)
        xyz = np.asarray(atoms.arrays["coords"], dtype=np.float32)[frame] - centre.astype(np.float32)
        rgb = residue_rgb(trace, groups)
        items = sc.frame(xyz, rgb)
        records, rebuilt = stored(items, name, arrays)
        scenes[name] = {"system": sys_name, "frame": frame, "backbone_atom": backbone, "rep_rest": rest,
                        "specs": [spec_record(s) for s in sc.specs],
                        "masks": [st.mask if st is not None else None for st in sc.styles],
                        "drawn": [st.drawn() if st is not None else None for st in sc.styles],
                        "scene_drawn": sc.drawn(), "centre": centre, "residue_rgb": rgb,
                        "xyz_dtype": "float32", "items": records}
        radius = max(float(np.linalg.norm(coords[frame] - centre, axis=1).max()), 1.0)
        for rname, angles, quality, fog, glow_every in SCENE_RENDERS.get(name, []):
            n_res = len(trace.arrays["resid"])
            glow = (np.arange(n_res) % glow_every == 0).astype(np.float32) if glow_every else None
            shading = Shading(LIGHT_SURFACE, fog, glow, GLOW_RGB)
            lim = 1.15 * radius + 4.0
            vp = Viewport(euler_to_matrix(angles), (-lim, lim), (-lim, lim), SIZE, SIZE, (-radius, radius))
            renders[rname] = render_record(rebuilt, vp, shading, quality, render_arrays, rname,
                                           source={"scene": name, "rotation": angles})
    return scenes, renders


def render_record(items, vp, shading, quality, arrays, name, source):
    img = SoftwareRenderer().render(items, vp, shading, quality=quality)
    arrays[f"{name}/rgba"] = img
    return {"source": source,
            "viewport": {"rotation": vp.R, "xlim": list(vp.xlim), "ylim": list(vp.ylim), "width": vp.width,
                         "height": vp.height, "zrange": list(vp.zrange)},
            "shading": {"fog_rgb": shading.fog_rgb, "fog_strength": shading.fog_strength,
                        "glow": shading.glow, "glow_rgb": shading.glow_rgb},
            "quality": quality, "shape": list(img.shape), "coverage": int((img[..., 3] > 0).sum())}


# --- synthetic items (MDATK's tests/test_molecule_render.py) ----------------------------------------

def vp(R=None, lim=10.0, size=100, zrange=(-10.0, 10.0), height=None):
    return Viewport(np.eye(3) if R is None else R, (-lim, lim), (-lim, lim), size, height or size, zrange)


def spheres(xyz, r, rgb=(0.8, 0.2, 0.2), owner=None, opacity=1.0):
    xyz = np.atleast_2d(np.asarray(xyz, float))
    n = len(xyz)
    own = np.arange(n) if owner is None else np.asarray(owner)
    rgb = np.asarray(rgb, float)
    rgb = rgb if rgb.ndim == 2 else np.tile(rgb, (n, 1))
    return Colored(Spheres(xyz, np.broadcast_to(np.asarray(r, float), (n,)).copy(), own, np.full(n, -1)), rgb,
                   opacity)


def mesh(verts, tris, owner=None, rgb=(0.5, 0.5, 0.5), normals=None, opacity=1.0):
    verts = np.asarray(verts, float)
    n = len(verts)
    normals = np.tile([0, 0, 1.0], (n, 1)) if normals is None else normals
    own = np.zeros(n, np.int64) if owner is None else np.asarray(owner)
    rgb = np.asarray(rgb, float)
    rgb = rgb if rgb.ndim == 2 else np.tile(rgb, (n, 1))
    return Colored(Mesh(verts, normals, np.asarray(tris), own, np.zeros(n, np.int8)), rgb, opacity)


def synthetic():
    """(name, items, viewport, shading, quality, what) as MDATK's renderer tests draw them."""
    out = []
    out.append(("sphere", [spheres([0, 0, 0], 5.0)], vp(), PLAIN, "high", "one sphere: area and depth"))
    st = Sticks(np.array([[-6.0, 0, 0]]), np.array([[6.0, 0, 0]]), np.array([1.0]), np.array([0]), np.array([-1]))
    out.append(("stick", [Colored(st, np.array([[0.2, 0.4, 0.8]]))], vp(), PLAIN, "high",
                "one stick: width, open ends"))
    for seed in range(2):
        rng = np.random.default_rng(seed)
        atoms = rng.normal(scale=2.5, size=(12, 3))
        verts, normals, tris, _near = gaussian_surface(atoms)
        R = euler_to_matrix(rng.uniform(-180, 180, 3))
        out.append((f"closed_mesh_{seed}", [mesh(verts, tris, normals=normals)], vp(R, lim=12, size=120), PLAIN,
                    "high", "a closed surface: no pinholes"))
    rng = np.random.default_rng(3)
    xyz, rgb = rng.uniform(-6, 6, (40, 3)), rng.uniform(0.2, 0.9, (40, 3))
    out.append(("forty_spheres", [Colored(Spheres(xyz, np.full(40, 1.6), np.arange(40), np.full(40, -1)), rgb)], vp(),
                PLAIN, "high", "overlapping spheres (order independent)"))
    out.append(("beyond_depth_range", [spheres([1.0, 0, 57.0], 3.0, (0, 0, 1)), spheres([0, 0, 60.0], 3.0, (1, 0, 0))],
                vp(zrange=(-10, 10)), PLAIN, "high", "geometry beyond the panel's depth range"))
    out.append(("fog", [spheres([-5, 0, 8.0], 2.0, (0.9, 0.1, 0.1)), spheres([5, 0, -8.0], 2.0, (0.9, 0.1, 0.1))], vp(),
                Shading(WHITE, 0.45), "high", "fog toward white"))
    one = [spheres([0, 0, 0], 2.0, (0.3, 0.3, 0.8), owner=[0])]
    out.append(("glow_off", one, vp(), Shading(WHITE, 0.0, np.zeros(1, np.float32)), "high", "glow 0"))
    out.append(("glow_on", one, vp(), Shading(WHITE, 0.0, np.ones(1, np.float32), np.array([1.0, 0.7, 0.0])), "high",
                "glow halo and tint"))
    out.append(("translucent_over", [spheres([0, 0, 0], 3.0, (1, 0, 0)),
                                     spheres([0, 0, 4.0], 3.0, (0, 0, 1), opacity=0.5)],
                vp(), PLAIN, "high", "a translucent sphere over an opaque one"))
    out.append(("translucent_alone", [spheres([0, 0, 4.0], 3.0, (0, 0, 1), opacity=0.5)], vp(), PLAIN, "high",
                "straight alpha of a translucent sphere"))
    out.append(("fast", [spheres([0, 0, 0], 4.0)], vp(size=101, height=77), PLAIN, "fast", "fast quality, odd size"))
    rng = np.random.default_rng(11)
    xyz = rng.normal(scale=4.0, size=(400, 3))
    culled = [Colored(Spheres(xyz, rng.choice([1.52, 1.7, 1.8], 400), np.arange(400), np.full(400, -1)),
                      rng.uniform(0.2, 0.9, (400, 3)))]
    for angles in ([0, 0, 0], [40, -70, 10]):
        out.append((f"culling_{'_'.join(str(a) for a in angles)}", culled,
                    vp(euler_to_matrix(angles), lim=14, size=120),
                    PLAIN, "high", "400 spheres: hidden ones culled"))
    rng = np.random.default_rng(2)
    many = [spheres(rng.uniform(-8, 8, 3), 0.8, rng.uniform(0.2, 0.9, 3)) for _ in range(80)]
    out.append(("many_items", many, vp(size=120), PLAIN, "high", "80 items"))
    front = spheres([0, 0, 2.0], 3.0, (1, 0, 0), opacity=0.5)
    back = spheres([0, 0, -2.0], 3.0, (0, 0, 1), opacity=0.5)
    out.append(("translucent_front_back", [back, front], vp(), PLAIN, "high", "front-most translucent layer once"))
    stick = Sticks(np.array([[-5.0, 0, 0]]), np.array([[5.0, 0, 0]]), np.array([0.6]), np.array([0]), np.array([-1]))
    out.append(("translucent_stick_ball", [Colored(stick, np.array([[0.2, 0.4, 0.9]]), 0.5),
                                           spheres([[-5.0, 0, 0], [5.0, 0, 0]], 0.6, (0.2, 0.4, 0.9), owner=[0, 0],
                                                   opacity=0.5)],
                vp(), PLAIN, "high", "a half-stick inside a ball adds no layer"))
    out.append(("translucent_glow", [spheres([0, 0, 0], 2.0, (0.3, 0.3, 0.8), owner=[0], opacity=0.6)], vp(),
                Shading(WHITE, 0.0, np.ones(1, np.float32), np.array([1.0, 0.7, 0])), "high",
                "a translucent glowing residue keeps its halo"))
    plane = mesh([[-9, -9, 0.0], [9, -9, 0.0], [9, 9, 0.0], [-9, 9, 0.0]], [[0, 1, 2], [0, 2, 3]], rgb=(0, 0, 1.0))
    sliver = mesh([[-6, 1.25, -0.5], [6, 1.25, -0.5], [0, 1.25 + 2e-5, -3.0]], [[0, 1, 2]], rgb=(1, 0, 0.0))
    out.append(("sliver_hidden", [plane, sliver], vp(size=120), PLAIN, "high", "an edge-on sliver behind a plane"))
    x = np.repeat([-6.0, 0.0, 6.0], 2)
    y = np.tile([-4.0, 4.0], 3)
    owner = np.array([0, 0, 0, 0, 1, 1])
    out.append(("colour_boundary", [mesh(np.stack([x, y, np.zeros(6)], 1), [[0, 2, 3], [0, 3, 1], [2, 4, 5], [2, 5, 3]],
                                         owner=owner, rgb=np.where(owner[:, None] == 0, [[1.0, 0, 0]], [[0, 0, 1.0]]))],
                vp(size=120), PLAIN, "high", "a straight colour seam between residues"))
    rng = np.random.default_rng(1)
    bad = rng.normal(scale=4, size=(60, 3))
    bad[0] = np.nan
    out.append(("non_finite", [spheres(bad, 1.5)], vp(lim=14), PLAIN, "high", "a non-finite atom is skipped"))
    return out


def synthetic_cases(arrays):
    renders = {}
    for name, items, view, shading, quality, what in synthetic():
        records, rebuilt = stored(items, name, arrays)
        rec = render_record(rebuilt, view, shading, quality, arrays, name, source={"synthetic": what})
        rec["items"] = records
        renders[name] = rec
    return renders


def main():
    if os.path.isdir(MOLECULE):
        shutil.rmtree(MOLECULE)
    os.makedirs(MOLECULE)
    topo = {}
    for name in SYSTEMS:
        atoms, trace, _groups = system(name)
        topo[name] = topology_record(topology(atoms, trace))
    atoms, trace, _groups = system("demo_dry")
    topo["demo_dry_c4"] = topology_record(topology(atoms, trace, "C4'"))
    scene_arrays, render_arrays = {}, {}
    scenes, renders = scene_cases(scene_arrays, render_arrays)
    renders.update(synthetic_cases(render_arrays))
    rc.write_json(os.path.join(MOLECULE, "molecule_cases.json"),
                  improvements_render.adjust("molecule", "molecule_cases",
                                             {"improvements": improvements_render.applied(), "topology": topo,
                                              "scenes": scenes, "renders": renders}))
    rc.write_npz(os.path.join(MOLECULE, "scenes.npz"), scene_arrays)
    rc.write_npz(os.path.join(MOLECULE, "renders.npz"), render_arrays)


if __name__ == "__main__":
    main()
