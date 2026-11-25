"""MDBridge's rendering improvements over MDATK v0.0.4.1 (docs/specs/m4-plan.md section 4, items
RN-01 ... RN-17), as test-side patches of MDATK, so that the expected data of the render
generators (make_figure_cases.py, make_chart_cases.py, make_molecule_cases.py,
make_render_cli_cases.py; tests/data/render) is MDBridge's improved behaviour, in the style of the
other improvements_*.py modules.

How it is used: render_common.setup() calls ``apply_mdatk()`` once, after MDATK is importable and
before any figure is built; it applies the patch of every item in PATCHES whose function is set.
The generators also pass what they record through ``adjust(kind, case, record)`` (kind "figure",
"chart", "molecule" or "cli"), for improvements that are easier to express on the record than as a
patch of MDATK (messages, file lists).

Rules (CLAUDE.md "Improve, don't copy bugs"): an item gets its patch together with the C++ code that
implements it and its docs/DIFFERENCES.md entry (what, why, how tested), and the regenerated data
is committed with them. Until then its entry below is None and the recorded data is MDATK's own
behaviour; the cases that show the bug are recorded anyway (make_figure_cases.py records MDATK's
error, e.g. RN-04's ValueError at draw time), so the item's chunk only has to fill in the patch and
regenerate. Conventions (file names, options, output layout) are never patched here.

Status of the items (see the plan for the full text):

  RN-01  heatmap / matrix-chart column edges at the midpoints of non-uniform frame times
  RN-02  no analysed frames: "No data: no frames" panels, render / preview stop with an error
  RN-03  non-finite coordinates / series values: finite-only extents and a warning
  RN-04  vmin > vmax clamped with a warning (MDATK: ValueError at draw time)
  RN-05  matrix chart height capped at 2^15 px
  RN-06  canvas exactly figure_size(w, h) for any DPI
  RN-07  codec names mapped to the bundled encoders (video metadata, make_render_cli_cases.py)
  RN-08  ion legend label from the ions' elements
  RN-09  chain legend with chain letters (owner decision D11)
  RN-10  chart series naming a matrix key: charted as a heatmap or a warning
  RN-11  chart / preview formats validated in check()
  RN-12  colour strings validated in check()
  RN-13  chain_palette validated in check()
  RN-14  NaN-aware running mean
  RN-15  panel option types validated in check()
  RN-16  renderer "vtk" means the GPU renderer (Vulkan), never VTK
  RN-17  video positions with a stride skip the last frame: kept (convention, no patch)
  RN-21 ... RN-25  residue shading, orientation and the example dataset: improvements_shading.py

None of them is patched yet (milestone 4, wave 0: the generators and this skeleton).

RN-26 ... RN-29 (the software 3D renderer, W1-B; src/render/molecule/software.cpp) change only the
renderer's own images, which make_molecule_cases.py alone records: it calls
``apply_software_renderer()`` (SOFTWARE_RENDERER below) and lists them in molecule_cases.json.
RN-26 (anti-aliasing with the region's real neighbours) and RN-27 (fast quality at odd sizes) are
patches of MDATK's software.py; RN-28 (unusable views and primitives raise / are skipped) and RN-29
(the halo normalisation for the exact size; no collision in the recorded renders) change no recorded
image and have no patch. The structure panel's chunk moves them into PATCHES if its figure cases
record the 3D layer's pixels.

RN-30 ... RN-32 (the density shadow, W0-I) are applied by make_density_shadow_cases.py itself,
around MDATK's project_density / _above_level / anchored; no figure case uses a map with non-finite
voxels or an unusable anchor, so they are not in PATCHES (the structure panel's chunk adds them
here when its figure cases need them).
"""

from __future__ import annotations

# item -> function(mdatk modules are importable) that patches MDATK in this process, or None while
# the item is not implemented in MDBridge (the expected data is then MDATK's behaviour)
PATCHES = {
    "RN-01": None,
    "RN-02": None,
    "RN-03": None,
    "RN-04": None,
    "RN-05": None,
    "RN-06": None,
    "RN-07": None,
    "RN-08": None,
    "RN-09": None,
    "RN-10": None,
    "RN-11": None,
    "RN-12": None,
    "RN-13": None,
    "RN-14": None,
    "RN-15": None,
    "RN-16": None,
    "RN-17": None,  # a documented convention, never patched
}



def _rn26_antialias_with_neighbours():
    """composite(): the 3 x 3 alpha mean of the anti-aliasing sees the image around the molecule's
    box (MDATK: uniform_filter of the box alone, zeros outside it, a darker rectangle in a halo).
    Where the image around the box is clear the result is MDATK's to the bit."""
    import numpy as np
    from mdatk.render.molecule import software as sw
    from scipy import ndimage

    original = sw.composite

    def composite(cam, pix, rgb, depth, owner, shading, outline=0.5, outline_depth=1.6, aa=True,
                  halo_extra=None):
        out = original(cam, pix, rgb, depth, owner, shading, outline, outline_depth, False, halo_extra)
        if aa and len(pix):
            W, H = cam.W, cam.H
            py, px = pix // W, pix % W
            y0, y1 = max(int(py.min()) - 2, 0), min(int(py.max()) + 3, H)
            x0, x1 = max(int(px.min()) - 2, 0), min(int(px.max()) + 3, W)
            ey0, ey1, ex0, ex1 = max(y0 - 1, 0), min(y1 + 1, H), max(x0 - 1, 0), min(x1 + 1, W)
            a = out[ey0:ey1, ex0:ex1, 3].astype(np.float32)
            a3 = ndimage.uniform_filter(a, 3, mode="constant")
            inner = (slice(y0 - ey0, y1 - ey0), slice(x0 - ex0, x1 - ex0))
            a, a3 = a[inner], a3[inner]
            out[y0:y1, x0:x1, 3] = np.where(np.abs(a3 - a) > 1, 0.5 * (a + a3), a).astype(np.uint8)
        return out

    sw.composite = composite


def _rn27_fast_odd_sizes():
    """_Cam(scale=0.5): half the size rounded up, the limits widened to two pixels per half pixel
    (MDATK: int(round(size / 2)), a transparent last column / row at odd sizes and a view 1 px too
    small). Even sizes are unchanged."""
    from mdatk.render.molecule import software as sw

    original = sw._Cam.__init__

    def init(self, vp, items, scale=1.0):
        original(self, vp, items, scale)
        if scale != 0.5:
            return
        self.W, self.H = max(2, (vp.width + 1) // 2), max(2, (vp.height + 1) // 2)
        if 2 * self.W != vp.width:
            sx = vp.width / (vp.xlim[1] - vp.xlim[0])
            self.xlim = (vp.xlim[0], vp.xlim[0] + (2 * self.W) / sx)
        if 2 * self.H != vp.height:
            sy = vp.height / (vp.ylim[1] - vp.ylim[0])
            self.ylim = (vp.ylim[1] - (2 * self.H) / sy, vp.ylim[1])
        self.s = self.W / (self.xlim[1] - self.xlim[0])
        self.sy = self.H / (self.ylim[1] - self.ylim[0])

    sw._Cam.__init__ = init


# the software renderer's items (make_molecule_cases.py), applied by apply_software_renderer()
SOFTWARE_RENDERER = {
    "RN-26": _rn26_antialias_with_neighbours,
    "RN-27": _rn27_fast_odd_sizes,
}

_software_applied = False


def apply_software_renderer() -> None:
    """Patch MDATK's software renderer (idempotent)."""
    global _software_applied
    if _software_applied:
        return
    _software_applied = True
    for item in sorted(SOFTWARE_RENDERER):
        SOFTWARE_RENDERER[item]()


def software_renderer_applied() -> list:
    return sorted(SOFTWARE_RENDERER)


# item -> function(kind, case, record) -> record, for improvements applied to the recorded data
ADJUSTMENTS = {}

_applied = False


def applied() -> list:
    """The items whose patch is active (recorded in every generated file)."""
    return sorted(k for k, f in PATCHES.items() if f is not None)


def apply_mdatk() -> None:
    """Patch the imported MDATK modules (idempotent)."""
    global _applied
    if _applied:
        return
    _applied = True
    for item in sorted(PATCHES):
        if PATCHES[item] is not None:
            PATCHES[item]()


def adjust(kind: str, case: str, record):
    for item in sorted(ADJUSTMENTS):
        record = ADJUSTMENTS[item](kind, case, record)
    return record
