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
  RN-10  chart series naming a matrix key: charted as a matrix chart (a pair map: a warning)
  RN-11  chart formats checked before anything is drawn (and in check(), W3-B)
  RN-12  colour strings validated in check()
  RN-13  chain_palette validated in check()
  RN-14  NaN-aware running mean
  RN-15  panel option types validated in check()
  RN-16  renderer "vtk" means the GPU renderer (Vulkan), never VTK
  RN-17  video positions with a stride skip the last frame: kept (convention, no patch)
  RN-21 ... RN-25  residue shading, orientation and the example dataset: improvements_shading.py

Patched so far: RN-01, RN-04, RN-14 (the plot panels, W2-A; src/render/panels/plots.cpp); RN-01 also
for the matrix charts, RN-05, RN-10, RN-11, RN-50 and D5 (the standalone charts, W2-D;
src/render/panels/charts.cpp):

  RN-50  charts' title and colour bar in the theme's colours (MDATK styled the centre title, not
         the left one it draws, and left the colour bar black: unreadable on the dark theme)
  D5     output.charts.style: "mdatk" (or "") and "default" (Matplotlib's default look); other
         names warn and draw MDATK's chrome (owner decision D5; MDATK applied any Matplotlib style)

The others are MDATK's behaviour until their chunk fills them in.

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
    "RN-01": "_rn01_heatmap_column_edges",
    "RN-02": None,
    "RN-03": None,
    "RN-04": "_rn04_heatmap_vmin_above_vmax",
    "RN-05": "_rn05_matrix_chart_height_cap",
    "RN-06": None,
    "RN-07": None,
    "RN-08": None,
    "RN-09": None,
    "RN-10": "_rn10_chart_matrix_series",
    "RN-11": "_rn11_chart_formats",
    "RN-12": None,
    "RN-13": None,
    "RN-14": "_rn14_nan_aware_running_mean",
    "RN-15": None,
    "RN-16": None,
    "RN-17": None,  # a documented convention, never patched
    # the standalone charts (W2-D; src/render/panels/charts.cpp)
    "RN-50": "_rn50_chart_chrome_in_theme_colours",
    "D5": "_d5_chart_styles",  # owner decision D5: the styles mdatk and default
}



def _rn01_heatmap_column_edges():
    """HeatmapPanel: the frame cells span the midpoints between frame times, half a step outward at
    both ends (MDATK: every column t[1] - t[0] wide, misplaced for non-uniform times). The image
    keeps imshow's extent (first and last edge) and carries the edges as ``_mdbridge_x_edges``,
    which mpl_walker.py records as the image's "x_edges" (MDBridge draws columns of their own
    widths). Uniform times: MDATK's cells; one frame: t0 -/+ 0.5 as MDATK; a matrix whose column
    count is not the frame count keeps MDATK's cells."""
    import numpy as np
    from mdatk.render import panels

    original = panels.HeatmapPanel.setup

    def setup(self):
        original(self)
        t = np.asarray(self.ctx.time, dtype=np.float64)
        images = self.ax.get_images()
        if not images or len(t) == 0:
            return
        im = images[0]
        if im.get_array().shape[1] != len(t):
            return
        if len(t) == 1:
            edges = np.array([t[0] - 0.5, t[0] + 0.5])
        else:
            edges = np.empty(len(t) + 1)
            edges[0] = t[0] - (t[1] - t[0]) / 2
            edges[1:-1] = (t[:-1] + t[1:]) / 2
            edges[-1] = t[-1] + (t[-1] - t[-2]) / 2
        left, right, bottom, top = im.get_extent()
        im.set_extent((float(edges[0]), float(edges[-1]), bottom, top))
        im._mdbridge_x_edges = [float(e) for e in edges]

    panels.HeatmapPanel.setup = setup
    _rn01_matrix_chart_column_edges()


def _frame_edges(t):
    import numpy as np

    if len(t) == 1:
        return np.array([t[0] - 0.5, t[0] + 0.5])
    edges = np.empty(len(t) + 1)
    edges[0] = t[0] - (t[1] - t[0]) / 2
    edges[1:-1] = (t[:-1] + t[1:]) / 2
    edges[-1] = t[-1] + (t[-1] - t[-2]) / 2
    return edges


def _rn01_matrix_chart_column_edges():
    """render_matrix_chart: the same frame cells as the heatmap panel (RN-01): the image's extent
    spans the first to the last edge and carries the edges as ``_mdbridge_x_edges``."""
    import numpy as np
    from matplotlib.axes import Axes
    from mdatk.render import charts

    original = charts.render_matrix_chart

    def render_matrix_chart(dataset, project, key, out_dir, stamp=""):
        t = np.asarray(dataset.time(project.inputs.time_unit), dtype=np.float64)
        imshow = Axes.imshow

        def patched(self, X, *a, **kw):
            im = imshow(self, X, *a, **kw)
            if len(t) and np.shape(X)[1] == len(t):
                edges = _frame_edges(t)
                left, right, bottom, top = im.get_extent()
                im.set_extent((float(edges[0]), float(edges[-1]), bottom, top))
                im._mdbridge_x_edges = [float(e) for e in edges]
            return im

        Axes.imshow = patched
        try:
            return original(dataset, project, key, out_dir, stamp)
        finally:
            Axes.imshow = imshow

    charts.render_matrix_chart = render_matrix_chart


def _rn04_heatmap_vmin_above_vmax():
    """HeatmapPanel (sequential maps): vmin above vmax (the options, or vmin above the largest
    value) swapped, with a warning. MDATK did not raise (the plan's premise): its colour bar ran
    nonsingular() on the reversed norm through its callbacks, which set vmin while the old vmax was
    still in place, so vmin 5 / vmax 1 became 0.9 / 5 (an artefact of the callback order) without a
    word; Matplotlib alone raises "minvalue must be less than or equal to maxvalue" at draw time.
    A negative-valued matrix with the default vmin 0 never failed (vmax is at least 1e-9)."""
    import numpy as np
    from mdatk.render import panels

    original = panels.HeatmapPanel.setup

    def setup(self):
        try:
            m = self.ctx.dataset.matrix(self.spec.series[0]) if self.spec.series else None
        except KeyError:
            m = None
        if m is not None and m.values.size and not (m.kind == "class" and m.classes):
            values = np.ma.masked_invalid(np.asarray(m.values, dtype=np.float64))
            vmax = float(self.options.get("vmax", 0) or (values.max() if values.count() else 1.0))
            vmax = max(vmax, 1e-9)
            vmin = float(self.options.get("vmin", 0.0))
            if vmin > vmax:
                panels.log.warning("heatmap %s: vmin %s is above vmax %s; using vmin %s, vmax %s", m.key,
                                   f"{vmin:g}", f"{vmax:g}", f"{vmax:g}", f"{vmin:g}")
                self.options = {**self.options, "vmin": vmax, "vmax": vmin}
        original(self)

    panels.HeatmapPanel.setup = setup


def _rn14_nan_aware_running_mean():
    """_running_mean (the time-series panel and the series charts): the mean of the finite values
    of each window (edge padding and the window clipped to the series as MDATK), NaN only where the
    whole window is not finite (MDATK: np.convolve spread one NaN over every window it fell in)."""
    import numpy as np
    from mdatk.render import charts, panels

    def running_mean(y, window):
        if window <= 1 or len(y) < 2:
            return y
        window = min(window, len(y))
        pad_l = window // 2
        pad_r = window - 1 - pad_l
        padded = np.pad(np.asarray(y, dtype=np.float64), (pad_l, pad_r), mode="edge")
        finite = np.isfinite(padded)
        out = np.empty(len(y))
        for i in range(len(y)):
            w, f = padded[i:i + window], finite[i:i + window]
            n = int(f.sum())
            total = 0.0
            for v in w[f]:  # left to right, as MDBridge sums
                total += float(v)
            out[i] = total / n if n else np.nan
        return out

    panels._running_mean = running_mean
    charts._running_mean = running_mean


MAX_CHART_PIXELS = 32768  # RN-05: 2^15


def _rn05_matrix_chart_height_cap():
    """render_matrix_chart: a figure height of more than 2^15 pixels (0.18 in per row) is capped
    there, with a warning (MDATK: ever taller images, beyond Agg's 2^16 px limit a crash). The
    row labels are thinned to 40 already."""
    import math

    from mdatk.render import charts

    original = charts.render_matrix_chart
    figure_class = charts.Figure

    def render_matrix_chart(dataset, project, key, out_dir, stamp=""):
        rows = len(dataset.matrix(key).row_labels)

        def capped(*a, figsize=None, dpi=None, **kw):
            w, h = figsize
            if int(h * dpi) > MAX_CHART_PIXELS:
                charts.log.warning("chart %s: %d rows need %d px at %s dpi; height capped at %d px", key, rows,
                                   int(h * dpi), f"{dpi:g}", MAX_CHART_PIXELS)
                h = MAX_CHART_PIXELS / dpi
                while int(h * dpi) < MAX_CHART_PIXELS:
                    h = math.nextafter(h, math.inf)
            return figure_class(*a, figsize=(w, h), dpi=dpi, **kw)

        charts.Figure = capped
        try:
            return original(dataset, project, key, out_dir, stamp)
        finally:
            charts.Figure = figure_class

    charts.render_matrix_chart = render_matrix_chart


def _rn10_chart_matrix_series():
    """render_charts: a matrix key in output.charts.series is charted as a matrix chart (after the
    series, before the include_heatmaps matrices, each matrix once); a pair-map key is skipped with
    a warning (MDATK dropped both silently)."""
    import os

    from mdatk.progress import Reporter
    from mdatk.render import charts

    def render_charts(dataset, project, out_dir, reporter=None, stamp=""):
        reporter = reporter or Reporter()
        os.makedirs(out_dir, exist_ok=True)
        keys = charts.chart_keys(dataset, project)
        mats = []
        for k in project.output.charts.series:
            if not dataset.has(k):
                continue
            kind = dataset.kind_of(k)
            if kind == "matrices":
                if k not in mats:
                    mats.append(k)
            elif kind != "series":
                charts.log.warning("chart series %s is a pair map; charts draw series and matrices "
                                   "(a contact_map panel shows it)", k)
        if project.output.charts.include_heatmaps:
            mats += [k for k in dataset.matrix_keys() if not k.startswith(charts._INTERNAL) and k not in mats]
        jobs = [("series", k) for k in keys] + [("matrix", k) for k in mats]
        paths = []
        for n, (kind, key) in enumerate(jobs):
            reporter.check_cancel()
            if kind == "series":
                paths += charts.render_series_chart(dataset, project, key, out_dir, stamp)
            else:
                paths += charts.render_matrix_chart(dataset, project, key, out_dir, stamp)
            reporter.progress((n + 1) / max(1, len(jobs)), f"Chart {key}")
        reporter.log(f"Wrote {len(paths)} chart file(s) to {out_dir}")
        return paths

    charts.render_charts = render_charts


# RN-11: the chart formats MDBridge writes (src/render/panels/charts.cpp chartFormats())
CHART_FORMATS = ("png", "svg", "pdf", "jpg", "jpeg", "tif", "tiff", "eps", "ps")


def _rn11_chart_formats():
    """render_charts: output.charts.formats checked before anything is drawn or written: a format
    MDBridge does not write is a ValueError listing them (MDATK: Matplotlib's ValueError after the
    files of the formats before it were written, and after the analyses ran)."""
    from mdatk.render import charts

    original = charts.render_charts

    def render_charts(dataset, project, out_dir, reporter=None, stamp=""):
        bad = [f for f in (project.output.charts.formats or ["png"]) if f.lower().lstrip(".") not in CHART_FORMATS]
        if bad:
            names = ", ".join(repr(f) for f in bad)
            raise ValueError(f"chart format{'s' if len(bad) > 1 else ''} {names} "
                             f"{'are' if len(bad) > 1 else 'is'} not supported (choose from {', '.join(CHART_FORMATS)})")
        return original(dataset, project, out_dir, reporter, stamp)

    charts.render_charts = render_charts


def _rn50_chart_chrome_in_theme_colours():
    """Charts with MDATK's chrome: the title (left) in the theme's ink at 11 pt, as style_axes meant
    (it styled the centre title, which charts leave empty: the left title stayed black, 12 pt); the
    matrix charts' colour bar outline in the grid colour, its ticks muted, tick labels and label in
    ink2 (MDATK: Matplotlib's black, unreadable on the dark theme)."""
    import matplotlib.axes
    import matplotlib.figure
    from mdatk.render import charts

    style_axes = charts.style_axes

    def styled(ax, theme, grid=True, font_size=10.0):
        style_axes(ax, theme, grid=grid, font_size=font_size)
        ax._mdbridge_chart_theme = (theme, font_size)

    charts.style_axes = styled
    set_title = matplotlib.axes.Axes.set_title

    def title(self, label, *a, **kw):
        t = set_title(self, label, *a, **kw)
        th = getattr(self, "_mdbridge_chart_theme", None)
        if th is not None:
            t.set_color(th[0].ink)
            t.set_fontsize(th[1])
        return t

    matplotlib.axes.Axes.set_title = title
    colorbar = matplotlib.figure.FigureBase.colorbar

    def bar(self, mappable, cax=None, ax=None, **kw):
        cb = colorbar(self, mappable, cax=cax, ax=ax, **kw)
        th = getattr(ax, "_mdbridge_chart_theme", None)
        if th is not None:
            theme = th[0]
            cb.outline.set_edgecolor(theme.grid)
            cb.ax.tick_params(colors=theme.muted, labelcolor=theme.ink2)
            cb.ax.yaxis.label.set_color(theme.ink2)
        return cb

    matplotlib.figure.FigureBase.colorbar = bar


def _d5_chart_styles():
    """output.charts.style (owner decision D5): "mdatk" or "" is MDATK's chrome, "default" is
    Matplotlib's default look; any other name (a Matplotlib style such as ggplot, a style file, a
    typo) warns and draws MDATK's chrome (MDATK applied any Matplotlib style and warned only for
    unknown ones, with Matplotlib's message)."""
    import contextlib

    from mdatk.render import charts

    original = charts._style

    @contextlib.contextmanager
    def _style(name):
        if not name or name == "mdatk":
            yield True
            return
        if name == "default":
            with original(name) as mdatk_chrome:
                yield mdatk_chrome
            return
        charts.log.warning("unknown matplotlib style %r (MDBridge draws the styles 'mdatk' and 'default'); "
                           "using MDATK chrome", name)
        yield True

    charts._style = _style


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


# PATCHES names the functions defined above (they follow the table)
PATCHES = {k: (globals()[v] if isinstance(v, str) else v) for k, v in PATCHES.items()}

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
