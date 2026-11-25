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
  D8     output.video.ffmpeg_path is not used (the bundled FFmpeg libraries): a note
  D9     no density fit in milestone 4: "density map not drawn: no placement (run density fit)"
  RN-21 ... RN-25  residue shading, orientation and the example dataset: improvements_shading.py
  RN-40  sequence panel: an unusable per_row (negative, not a number) is a "No data" message
  RN-41  secondary panel: an unknown source is a "No data" message (MDATK: silently "dynamic")
  RN-42  sequence / secondary panels: "No data" messages keep the frame (MDATK: none once the
         secondary panel's setup had hidden the spines)

  RN-33  the density map's contour level from its finite voxels (W2-B)
  RN-50  charts' title and colour bar in the theme's colours (MDATK styled the centre title, not
         the left one it draws, and left the colour bar black: unreadable on the dark theme)
  RN-51  the structure view's centre warning names its setting: "Structure view › Orientation ›
         Centre on: ..." (MDATK 0.1.x render/panels.py _view_centre)
  D5     output.charts.style: "mdatk" (or "") and "default" (Matplotlib's default look); other
         names warn and draw MDATK's chrome (owner decision D5; MDATK applied any Matplotlib style)

Patched so far: RN-01, RN-04, RN-14 (the plot panels, W2-A; src/render/panels/plots.cpp),
RN-40 ... RN-42 (the sequence / secondary panels, W2-C; src/render/panels/seqpanels.cpp), RN-03,
RN-08, RN-09, RN-16, RN-33 (the structure view, W2-B; src/render/panels/structure.cpp) and RN-01 (for
the matrix charts too), RN-05, RN-10, RN-11, RN-50, D5 (the standalone charts, W2-D;
src/render/panels/charts.cpp) and RN-02, RN-07 (the codec warnings), RN-12, RN-13, RN-15, D8, D9
(the render, charts and preview commands, W3-B; src/app/commands/render_pipeline.cpp,
src/render/panels/checks.cpp). The others are MDATK's behaviour until their chunk fills them in.

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

# --- the structure view (W2-B; src/render/panels/structure.cpp) ------------------------------------

# MDBridge's legend label of the ions (structure.cpp ionLegendLabel): residue name of the ion labels
# ("MG25" -> "MG"), charges and trailing digits dropped, upper case
_ION_NAMES = {"MG": "Mg²⁺", "NA": "Na⁺", "SOD": "Na⁺", "K": "K⁺", "POT": "K⁺", "CA": "Ca²⁺", "CAL": "Ca²⁺",
              "ZN": "Zn²⁺", "MN": "Mn²⁺", "CL": "Cl⁻", "CLA": "Cl⁻", "LI": "Li⁺", "RB": "Rb⁺", "CS": "Cs⁺"}


def ion_legend_label(labels) -> str:
    found = set()
    for label in labels:
        key = "".join(c for c in str(label).rstrip("0123456789") if c not in "+- ").upper().rstrip("0123456789")
        found.add(_ION_NAMES.get(key, "Ions"))
    return found.pop() if len(found) == 1 else "Ions"


def _rn03_finite_structure_extent():
    """StructurePanel.setup: the bounding sphere of the finite trace / ion points (MDATK: one NaN
    point made the radius NaN, "Axis limits cannot be NaN or Inf"), and one warning naming the
    frames with non-finite points."""
    import numpy as np
    from mdatk.render import panels
    from mdatk.render.residues import warn_once

    class _Linalg:
        def __getattr__(self, name):
            return getattr(np.linalg, name)

        @staticmethod
        def norm(x, *args, **kwargs):
            r = np.linalg.norm(x, *args, **kwargs)
            return np.where(np.isfinite(r), r, -np.inf)

    class _Numpy:
        linalg = _Linalg()

        def __getattr__(self, name):
            return getattr(np, name)

    original = panels.StructurePanel.setup

    def setup(self):
        tr = self.ctx.dataset.results.get("trace")
        bad, text = set(), ""
        if tr is not None:
            bad = set()
            for key in ("coords", "ion_coords"):
                a = np.asarray(tr.arrays[key], dtype=np.float64)
                if a.size:
                    bad |= set(np.flatnonzero(~np.isfinite(a).all(axis=(1, 2))).tolist())
            if bad:
                frames = [int(self.ctx.dataset.frames[f]) for f in sorted(bad)]
                text = ", ".join(map(str, frames[:10])) + (f" and {len(frames) - 10} more" if len(frames) > 10 else "")
        saved = panels.np
        panels.np = _Numpy()
        try:
            original(self)
        finally:
            panels.np = saved
        if bad:
            warn_once(f"structure view: non-finite coordinates in frame(s) {text}; left out of the view "
                      "(those residues are not drawn there)")

    panels.StructurePanel.setup = setup


def _rn08_rn09_legend():
    """StructurePanel._legend: chains named as the trace names them ("Chain A"; the number only for
    a chain without a name; RN-09, owner decision D11) and the ions by their element (RN-08)."""
    import numpy as np
    from mdatk.render import panels
    from mdatk.render.residues import chain_names

    original = panels.StructurePanel._legend

    def _legend(self, v, tr):
        original(self, v, tr)
        if self.legend is None:
            return
        names = chain_names(tr)
        chain = np.asarray(tr.arrays["chain"])
        for text, handle in zip(self.legend.get_texts(), self.legend.legend_handles):
            label = text.get_text()
            new = label
            if v.color_mode == "chain" and label.startswith("Chain "):
                c = int(label.split()[1]) - 1
                k = int(np.flatnonzero(chain == c)[0])
                new = f"Chain {names[k] or c + 1}"
            elif label == "Mg²⁺":
                new = ion_legend_label(tr.meta.get("ion_labels") or [])
            if new != label:
                text.set_text(new)
                handle.set_label(new)

    panels.StructurePanel._legend = _legend


def _rn09_with_rn08():
    """RN-09 is applied by _rn08_rn09_legend (the same legend texts)."""


def _rn16_gpu_renderer():
    """renderer "vtk" means the GPU renderer (Vulkan); the figure cases have none: the built-in
    renderer with MDBridge's note (MDATK: "VTK is not installed (pip install vtk)")."""
    from mdatk.render.molecule import vtk_backend

    vtk_backend.available = lambda: (False, "the GPU renderer (Vulkan) is not available")


def _rn33_finite_map_level():
    """map_level from the finite voxels (MDATK: NaN for a map with one NaN voxel, so no shadow)."""
    import numpy as np
    from mdatk import density

    original = density.map_level

    def map_level(dmap, sigmas):
        d = np.asarray(dmap.data)
        finite = np.isfinite(d)
        if finite.all():
            return original(dmap, sigmas)
        if not finite.any():
            return float("nan")
        v = d[finite].astype(np.float64)
        return float(v.mean() + sigmas * v.std())

    density.map_level = map_level


# item -> function(mdatk modules are importable) that patches MDATK in this process, or None while
# the item is not implemented in MDBridge (the expected data is then MDATK's behaviour)
PATCHES = {
    "RN-01": "_rn01_heatmap_column_edges",
    "RN-02": "_rn02_no_frames",
    "RN-03": "_rn03_finite_structure_extent",
    "RN-04": "_rn04_heatmap_vmin_above_vmax",
    "RN-05": "_rn05_matrix_chart_height_cap",
    "RN-06": None,
    "RN-07": "_rn07_codec_messages",
    "RN-08": "_rn08_rn09_legend",
    "RN-09": "_rn09_with_rn08",
    "RN-10": "_rn10_chart_matrix_series",
    "RN-11": "_rn11_chart_formats",
    "RN-12": "_rn12_rn13_rn15_check",
    "RN-13": "_rn13_with_rn12",
    "RN-14": "_rn14_nan_aware_running_mean",
    "RN-15": "_rn15_with_rn12",
    "RN-16": "_rn16_gpu_renderer",
    "RN-17": None,  # a documented convention, never patched
    "RN-40": "_rn40_per_row",
    "RN-41": "_rn41_unknown_source",
    "RN-42": "_rn42_no_data_frame",
    "RN-33": "_rn33_finite_map_level",
    # the standalone charts (W2-D; src/render/panels/charts.cpp)
    "RN-50": "_rn50_chart_chrome_in_theme_colours",
    "D5": "_d5_chart_styles",  # owner decision D5: the styles mdatk and default
    # the render, charts and preview commands (W3-B; src/app/commands/render_pipeline.cpp)
    "D8": "_d8_ffmpeg_path_not_used",  # owner decision D8: the bundled FFmpeg libraries
    "D9": "_d9_no_density_fit",  # owner decision D9: no density fit in milestone 4
    # the structure view's centre warning names its setting (MDATK 0.1.x render/panels.py)
    "RN-51": "_rn51_centre_warning_names_its_setting",
}


def _rn51_centre_warning_names_its_setting():
    """StructurePanel._view_centre: "Structure view › Orientation › Centre on: <error>; using the model
    centre" (MDATK 0.1.x's readable location of visual.center_on) instead of "structure view centre:
    ..."."""
    from mdatk.render import panels

    original = panels.warn_once

    def warn_once(message):
        prefix = "structure view centre: "
        if message.startswith(prefix):
            message = "Structure view › Orientation › Centre on: " + message[len(prefix):]
        return original(message)

    panels.warn_once = warn_once


def _rn40_per_row():
    """SequencePanel.setup: per_row that int() cannot read, or a negative one, is a KeyError (the
    panel's "No data" message) once the residues and the automatic title are known (MDATK: a
    ValueError / TypeError out of build_figure, or a reversed axis with huge boxes). 0, empty and
    values int() truncates to 0 stay "automatic"."""
    from mdatk.render import seqpanels

    original = seqpanels.SequencePanel.setup

    def setup(self):
        value = self.options.get("per_row", 0)
        n = 0
        if value:
            try:
                n = int(value)
            except (TypeError, ValueError, OverflowError):
                n = None
        if n is None or n < 0:
            self._residues()
            self._bare_axes()
            self.auto_title = "Sequence" + (f" ({self.range_text})" if self.range_text else "")
            raise KeyError("option 'per_row' must be a whole number of residues per row (0: automatic), "
                           f"not {value!r}")
        original(self)

    seqpanels.SequencePanel.setup = setup


def _rn41_unknown_source():
    """SecondaryStructurePanel.setup: a source other than dynamic / first / dotbracket is a KeyError
    once the residues and the automatic title are known (MDATK: treated as "dynamic")."""
    from mdatk.render import seqpanels

    original = seqpanels.SecondaryStructurePanel.setup

    def setup(self):
        source = str(self.options.get("source", "dynamic"))
        if source not in ("dynamic", "first", "dotbracket"):
            self._residues()
            self._bare_axes()
            self.source = source
            self.auto_title = "Secondary structure" + (f" ({self.range_text})" if self.range_text else "")
            raise KeyError(f"unknown source {source!r} (option 'source': dynamic, first or dotbracket)")
        original(self)

    seqpanels.SecondaryStructurePanel.setup = setup


def _rn42_no_data_frame():
    """_ResiduePanel.show_message: the spines style_axes shows (left, bottom) stay visible, also when
    setup had hidden them (_bare_axes) before its KeyError (MDATK: a "No data" panel without a frame)."""
    from mdatk.render import panels, seqpanels

    def show_message(self, text):
        panels.Panel.show_message(self, text)
        for side in ("left", "bottom"):
            self.ax.spines[side].set_visible(True)
        for side in ("top", "right"):
            self.ax.spines[side].set_visible(False)

    seqpanels._ResiduePanel.show_message = show_message
    panels.StructurePanel.show_message = show_message  # W2-B: "No data: proximity" / "atoms"


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


# --- the render, charts and preview commands (W3-B; src/app/commands/render_pipeline.cpp,
# src/render/panels/checks.cpp) -------------------------------------------------------------------

NO_FRAMES = "no analysed frames to render"
# render/colormap.hpp bundledColormapNames (owner decision D6)
BUNDLED_COLORMAPS = ("tab10", "tab20", "tab20b", "tab20c", "Set1", "Set2", "Set3", "Pastel1", "Pastel2", "Dark2",
                     "Accent", "Paired", "viridis", "plasma", "inferno", "magma", "cividis")


def _colormap_problem(name):
    if name in BUNDLED_COLORMAPS or (name.endswith("_r") and name[:-2] in BUNDLED_COLORMAPS):
        return None
    message = (f"{name!r} is not a bundled colormap; choose from {', '.join(BUNDLED_COLORMAPS)} "
               '(or one of them reversed, with "_r")')
    lower = name.lower() if name.isascii() else "".join(c.lower() if "A" <= c <= "Z" else c for c in name)
    for candidate in BUNDLED_COLORMAPS:
        for spelled in (candidate, candidate + "_r"):
            if spelled.lower() == lower:
                return f"{message}; did you mean {spelled!r}?"
    return message


def _reads_as_number(v):
    if v is None or isinstance(v, (bool, int, float)):
        return True
    if isinstance(v, str):
        try:
            float(v)
            return True
        except ValueError:
            return False
    return False


def _reads_as_int(v):
    """(ok, value) as Panel::optionInt reads an option."""
    import math

    if v is None:
        return True, None
    if isinstance(v, bool):
        return True, int(v)
    if isinstance(v, int):
        return True, v
    if isinstance(v, float):
        return (True, math.trunc(v)) if math.isfinite(v) else (False, None)
    if isinstance(v, str):
        try:
            return True, int(v)
        except ValueError:
            return False, None
    return False, None


def render_problems(project):
    """src/render/panels/checks.cpp renderProblems: RN-11 (chart formats), RN-12 (colours), RN-13
    (chain palette), RN-15 / RN-04 (panel options), in this order, each naming its field."""
    from matplotlib.colors import to_rgba

    out = []
    formats = project.output.charts.formats
    bad = [f for f in formats if f.lower().lstrip(".") not in CHART_FORMATS]
    if bad:
        names = ", ".join(repr(f) for f in bad)
        out.append(f"output.charts.formats: chart format{'s' if len(bad) > 1 else ''} {names} "
                   f"{'are' if len(bad) > 1 else 'is'} not supported (choose from {', '.join(CHART_FORMATS)})")

    def color(field, value):
        try:
            to_rgba(value)
        except ValueError as exc:
            out.append(f"{field}: {exc}")

    for k, g in enumerate(project.groups):
        color(f"groups[{k}].color", g.color)
    v = project.visual
    for name in ("base_color", "proximity_color", "glow_color", "ion_color"):
        color(f"visual.{name}", getattr(v, name))
    if project.density.map_path:
        color("density.color", project.density.color)
    if project.output.charts.color:
        color("output.charts.color", project.output.charts.color)
    problem = _colormap_problem(v.chain_palette)
    if problem:
        out.append(f"visual.chain_palette: {problem}")
    for k, spec in enumerate(project.layout.panels):
        opts = spec.options or {}
        head = f"layout.panels[{k}].options."

        def number(key):
            if key in opts and not _reads_as_number(opts[key]):
                out.append(f"{head}{key}: expected a number, got {opts[key]!r}")

        def integer(key, minimum=None):
            if key not in opts:
                return
            ok, value = _reads_as_int(opts[key])
            if not ok:
                out.append(f"{head}{key}: expected a whole number, got {opts[key]!r}")
            elif minimum is not None and value is not None and value < minimum:
                out.append(f"{head}{key}: expected a whole number of at least {minimum}, got {opts[key]!r}")

        def pair(key):
            if key not in opts:
                return
            x = opts[key]
            ok = isinstance(x, list) and len(x) == 2 and all(
                isinstance(e, (int, float)) and not isinstance(e, bool) for e in x)
            if not ok:
                out.append(f"{head}{key}: expected two numbers, got {x!r}")

        def choice(key, choices):
            if key in opts and opts[key] not in choices:
                listed = ", ".join(choices[:-1]) + " or " + choices[-1]
                out.append(f"{head}{key}: expected {listed}, got {opts[key]!r}")

        if spec.type == "timeseries":
            integer("smooth")
            number("target")
            pair("ylim")
        elif spec.type == "heatmap":
            number("vmin")
            number("vmax")
            integer("max_ticks", 1)
        elif spec.type == "contact_map":
            choice("mode", ("instant", "frame", "frequency"))
            integer("frame")
            integer("max_ticks", 1)
    return out


def _rn12_rn13_rn15_check():
    """check(): Project.validate() also reports render_problems() (RN-11 ... RN-15: MDATK found these
    mid-render, after the analyses, as Matplotlib tracebacks); the runs from the cache (which skip
    check()) report them alone before reading the cache."""
    from mdatk import config
    from mdatk import pipeline as pipeline_module

    validate = config.Project.validate

    def validate_with_render(self):
        return validate(self) + render_problems(self)

    config.Project.validate = validate_with_render
    load = pipeline_module.Pipeline.load_cached_dataset

    def load_cached_dataset(self):
        problems = render_problems(self.project)
        if problems:
            raise pipeline_module.PipelineError("Invalid parameters:\n  - " + "\n  - ".join(problems))
        return load(self)

    pipeline_module.Pipeline.load_cached_dataset = load_cached_dataset


def _rn13_with_rn12():
    """RN-13 is applied by _rn12_rn13_rn15_check."""


def _rn15_with_rn12():
    """RN-15 is applied by _rn12_rn13_rn15_check."""


def _rn02_no_frames():
    """No analysed frames: the video, the charts and the preview stop with PipelineError("no
    analysed frames to render") before any figure is built (MDATK: IndexError in the panels, and
    the preview at position -1). The preview also checks its file's format before the analyses run
    (RN-11) and, for a name without an extension, prints the file Matplotlib writes (name.png)."""
    import os

    from mdatk import cli
    from mdatk import pipeline as pipeline_module

    Pipeline = pipeline_module.Pipeline
    video, charts = Pipeline.render_video, Pipeline.render_charts

    def render_video(self, dataset, reporter=None):
        if dataset.n_frames == 0:
            raise pipeline_module.PipelineError(NO_FRAMES)
        return video(self, dataset, reporter)

    def render_charts(self, dataset, reporter=None):
        bad = [f for f in (self.project.output.charts.formats or ["png"]) if f.lower().lstrip(".") not in CHART_FORMATS]
        if not bad and dataset.n_frames == 0:
            raise pipeline_module.PipelineError(NO_FRAMES)
        return charts(self, dataset, reporter)

    Pipeline.render_video, Pipeline.render_charts = render_video, render_charts

    def cmd_preview(args):
        project, pipe = cli._pipeline(args)
        ext = os.path.splitext(os.path.basename(args.output or "x.png"))[1]
        fmt = ext.lower().lstrip(".") if ext else "png"
        if fmt not in CHART_FORMATS:
            raise ValueError(f"preview format {fmt!r} is not supported (choose from {', '.join(CHART_FORMATS)})")
        ds = cli._dataset(pipe, args)
        if ds.n_frames == 0:
            raise pipeline_module.PipelineError(NO_FRAMES)
        pos = min(max(0, args.frame), ds.n_frames - 1)
        out = args.output or os.path.join(project.out_dir, f"MDATK_preview_{pos:05d}.png")
        if not ext:
            out += ".png"
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        print(pipe.render_preview(ds, pos, out))
        return 0

    cli.cmd_preview = cmd_preview
    parser = cli.build_parser

    def build_parser():
        p = parser()
        for action in p._subparsers._group_actions:
            for name, sub in action.choices.items():
                if name == "preview":
                    sub.set_defaults(func=cmd_preview)
        return p

    cli.build_parser = build_parser


# the encoders of MDBridge's FFmpeg build that a codec name can name (encoder.cpp chooseCodec)
MDBRIDGE_ENCODERS = ("libsvtav1", "mpeg4")


def mdbridge_codec_message(requested):
    """encoder.cpp chooseCodec's warning for output.video.codec (None: no warning)."""
    name = requested or "libx264"
    gpl = name in ("libx264", "libx265")
    if gpl or name in ("h264", "hevc", "av1"):
        return None if name == "av1" else f"video: codec {name} is not bundled{' (GPL)' if gpl else ''}; using libsvtav1 (AV1)"
    if name in MDBRIDGE_ENCODERS:
        return None
    return f"video: encoder {name} is not available; using libsvtav1"


def _rn07_codec_messages():
    """The codec names map to MDBridge's bundled encoders (RN-07, owner decision D1): MDBridge's
    warning instead of MDATK's fallback warning. MDATK's ffmpeg still encodes with what it has (the
    video's codec name is not compared, its metadata is)."""
    from mdatk.render import video

    choose = video.choose_codec

    def choose_codec(ffmpeg, requested):
        message = mdbridge_codec_message(requested)
        if message:
            video.log.warning(message)
        encoders = video.available_encoders(ffmpeg)
        for name in (requested, "libx264", "libopenh264", "mpeg4"):
            if not encoders or f" {name} " in encoders:
                return name
        return choose(ffmpeg, requested)

    video.choose_codec = choose_codec


def _d8_ffmpeg_path_not_used():
    """output.video.ffmpeg_path is not used (owner decision D8: the bundled FFmpeg libraries encode);
    another path than "ffmpeg" gets a note (MDATK: VideoError "ffmpeg not executable at P")."""
    from mdatk import pipeline as pipeline_module

    original = pipeline_module.Pipeline.render_video

    def render_video(self, dataset, reporter=None):
        vp = self.project.output.video
        path = vp.ffmpeg_path
        if dataset.n_frames and path and path != "ffmpeg":
            (reporter or self.reporter).log(f"Note: output.video.ffmpeg_path ({path}) is not used; MDBridge "
                                            "encodes with its bundled FFmpeg libraries")
        vp.ffmpeg_path = "ffmpeg"
        try:
            return original(self, dataset, reporter)
        finally:
            vp.ffmpeg_path = path

    pipeline_module.Pipeline.render_video = render_video


def _d9_no_density_fit():
    """A shown density map without a placement is not fitted (owner decision D9: no density fit in
    milestone 4): "WARNING: density map not drawn: no placement (run density fit)" (with auto_fit
    off MDATK's own warning)."""
    from mdatk import pipeline as pipeline_module

    def ensure_density_fit(self, reporter=None):
        p = self.project
        rep = reporter or self.reporter
        if p.density_active() and len(p.density.transform) != 12:
            if not p.density.auto_fit:
                rep.log("WARNING: the density map has no placement and automatic fitting is off; it is not drawn")
            else:
                rep.log("WARNING: density map not drawn: no placement (run density fit)")

    pipeline_module.Pipeline.ensure_density_fit = ensure_density_fit


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
