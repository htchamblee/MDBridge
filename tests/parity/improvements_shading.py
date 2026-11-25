"""MDBridge's improvements to MDATK v0.0.4.1's residue shading, orientation and example dataset
(docs/specs/m4-plan.md section 4, RN-21 ... RN-25; docs/DIFFERENCES.md "Rendering: residue shading,
orientation and the example dataset"), as test-side patches of MDATK in the style of the other
improvements_*.py modules. MDBridge's code: src/render/core/shading.cpp, orientation.cpp,
example.cpp; expected data: make_shading_cases.py (tests/data/render/shading).

render_common.setup() applies them in every render generator, after improvements_render. They
change only inputs no case of make_figure_cases.py, make_chart_cases.py, make_molecule_cases.py or
make_render_cli_cases.py uses (one-strand example structures, unreadable or chain-qualified
label_residues, angles next to -180 degrees, groups of more than 120 residues), so those files do not
list them in their "improvements" field (improvements_render.applied()).
"""

from __future__ import annotations

# orientation.cpp, example.cpp; tests/parity/make_shading_cases.py) ---------------------------------

def _rn20_single_strand_example_pairs():
    """RN-21: _example_basepairs sliced a[:m - 2] with m < 2 (one strand: m = 0), so a structure
    with one strand crashed the example dataset with an IndexError. MDBridge: no pairs then."""
    import numpy as np
    from mdatk.data import AnalysisResult, CSRBuilder
    from mdatk.render import example

    def _example_basepairs(chain, n_frames):
        a = np.flatnonzero(chain == chain[0]) if len(chain) else np.zeros(0, dtype=np.int64)
        b = np.flatnonzero(chain != chain[0])[::-1] if len(chain) else np.zeros(0, dtype=np.int64)
        m = min(len(a), len(b))
        builder = CSRBuilder()
        for f in range(n_frames):
            open_end = 2 if (f // 6) % 2 else 0
            count = max(m - open_end, 0)
            i, j = a[:count], b[:count]
            builder.add(i, j, np.ones(len(i)))
        res = AnalysisResult("basepairs", "example", {}, {}, {"example": True})
        res.arrays.update(builder.arrays("map"))
        labels = [str(k) for k in range(len(chain))]
        res.meta["map_r"] = labels
        res.meta["map_c"] = labels
        res.add_pairmap("map", "map", "map_r", "map_c", "Base pairs")
        return res

    example._example_basepairs = _example_basepairs


def _rn21_wrap_degrees():
    """RN-22: wrap_degrees folded every angle within np.isclose's 0.0018 degrees of -180 to +180;
    MDBridge folds -180 itself (within 1e-9)."""
    from mdatk.render import orientation

    def wrap_degrees(a: float) -> float:
        a = (a + 180.0) % 360.0 - 180.0
        return 180.0 if abs(a + 180.0) <= 1e-9 else a

    orientation.wrap_degrees = wrap_degrees


def _rn22_rn23_tag_rows():
    """RN-23: a label_residues spec that cannot be read gave no tags without a word; MDBridge warns
    once. RN-24: a chain in label_residues ("B:12") was ignored (every chain's residue 12 got a
    tag); MDBridge tags that chain only."""
    import numpy as np
    from mdatk.render import panels, residues, seqpanels
    from mdatk.selections import SelectionError, parse_residue_ranges

    def tag_rows(tr, v):
        if "resid" not in tr.arrays or not (v.label_every > 0 or v.label_residues.strip()):
            return np.zeros(0, dtype=np.int64)
        resid = np.asarray(tr.arrays["resid"])
        wanted = np.zeros(len(resid), bool)
        if v.label_every > 0:
            wanted |= (resid % int(v.label_every)) == 0
        if v.label_residues.strip():
            try:
                chains = None
                for chain, a, b in parse_residue_ranges(v.label_residues):
                    sel = (resid >= a) & (resid <= b)
                    if chain and "chain" in tr.arrays:
                        if chains is None:
                            chains = residues.chain_names(tr)
                        sel &= chains == chain
                    wanted |= sel
            except SelectionError as exc:
                residues.warn_once(f"visual.label_residues: {exc}; no residue-number tags from it")
        return np.flatnonzero(wanted)

    residues.tag_rows = tag_rows
    panels.tag_rows = tag_rows
    seqpanels.tag_rows = tag_rows


def _rn24_large_groups_in_example():
    """RN-25: group_resids listed every number of a range before keeping 400 (a range such as
    1-1000000000 took gigabytes), and schematic_trace's sequence had 120 letters (a group of more
    than 120 residues: IndexError). MDBridge stops at 400 residues and repeats the sequence."""
    import numpy as np
    from mdatk.render import example
    from mdatk.selections import parse_residue_ranges

    def group_resids(project, name, fallback=12):
        for g in project.groups:
            if g.name == name and g.residues.strip():
                try:
                    out = []
                    for _chain, a, b in parse_residue_ranges(g.residues):
                        out.extend(range(a, min(b, a + 400 - len(out) - 1) + 1) if len(out) < 400 else [])
                    return out[:400] or list(range(1, fallback + 1))
                except ValueError:
                    break
        return list(range(1, fallback + 1))

    def schematic_trace(project, n_frames):
        names = project.visual.structure_groups or project.group_names()[:2] or ["A"]
        coords, chain, gidx, resid, links, letters = [], [], [], [], [], []
        group_pos = {g.name: k for k, g in enumerate(project.groups)}
        for s, name in enumerate(names[:2]):
            ids = example.group_resids(project, name)
            seq = ("GCAUGGCAUCGA" * (len(ids) // 12 + 1))[:len(ids)]
            if s:
                seq = "".join({"G": "C", "C": "G", "A": "U", "U": "A"}[b] for b in reversed(seq))
            for k, r in enumerate(ids):
                level = k if s == 0 else len(ids) - 1 - k
                ang = np.radians(level * 32.7 + (0 if s == 0 else 150))
                coords.append([9.0 * np.cos(ang), 9.0 * np.sin(ang), 2.8 * level])
                chain.append(s)
                gidx.append(group_pos.get(name, -1))
                resid.append(r)
                letters.append(seq[k])
                if k:
                    links.append([len(coords) - 2, len(coords) - 1])
        xyz = np.asarray(coords, dtype=np.float64)
        return {"coords": np.repeat(xyz[None], n_frames, axis=0), "chain": np.asarray(chain),
                "resname": np.asarray(letters, dtype="U8"),
                "group_index": np.asarray(gidx), "resid": np.asarray(resid),
                "links": np.asarray(links, dtype=np.int64).reshape(-1, 2),
                "ions": np.zeros((0, 3))}

    example.group_resids = group_resids
    example.schematic_trace = schematic_trace


PATCHES = {
    "RN-21": _rn20_single_strand_example_pairs,
    "RN-22": _rn21_wrap_degrees,
    "RN-23": _rn22_rn23_tag_rows,  # RN-24 in the same patch (both are tag_rows)
    "RN-25": _rn24_large_groups_in_example,
}

_applied = False


def apply_mdatk() -> None:
    """Patch the imported MDATK modules (idempotent)."""
    global _applied
    if _applied:
        return
    _applied = True
    for item in sorted(PATCHES):
        PATCHES[item]()

