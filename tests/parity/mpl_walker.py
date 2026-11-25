"""The contents of a drawn Matplotlib figure as data: what make_figure_cases.py and
make_chart_cases.py record from MDATK's figures (tests/data/README.md "Rendering").

Image parity is contents only (owner decision D10, docs/specs/m4-plan.md section 3): MDBridge's
figures are compared with MDATK's by what they show - data, colours, ticks and labels, legends,
texts, layout boxes - never by pixels. ``figure(fig)`` walks a figure that has been drawn
(``canvas.draw()``: ticks, legend positions and the 3D layer are decided while drawing) and returns:

  figure     size_px [w, h], dpi, facecolor, suptitle, texts (fig.text: the watermark)
  axes       one record per axes in fig.axes order (each with its "id", "ax<k>"):
             position_px [x0, y0, x1, y1] (figure pixels, origin bottom left as in Matplotlib),
             tight_px (the box with tick labels, axis labels and title), facecolor, visible, axison,
             xlim, ylim, aspect, titles {left, center, right}, spines {side: visible, color, lw},
             xaxis / yaxis (below), and "artists": the axes' children in insertion order (the
             order Matplotlib draws artists of equal zorder in), each with an "id" "ax<k>/<n>"
  colorbars  on the colour-bar axes ("colorbar" in the axes record): vmin, vmax, the colour map
             (name, N, its lookup table in the case's .npz as "<id>.lut", over / under / bad),
             orientation, outline colour, the axis label

An axis record: label {text, color, size}, ticks [[loc, label text, label1 shown, label2 shown]]
(the ticks drawn: inside the view interval), tick style (tick colour, label colour, label size,
length, width, label rotation), grid lines (shown, colour, width), offset text, scale, inverted, the
locator (class; MaxNLocator: nbins as set and "auto_nbins", the value nbins="auto" resolved to from
the axis length and font size, steps, integer, min_n_ticks; FixedLocator: locs) and the formatter
(class; ScalarFormatter: offset, orderOfMagnitude, format; FixedFormatter: seq).

Artist records (common: kind, zorder, visible, alpha, animated, label when not "_..."):
  line        x, y, color, lw, ls (name), dashes [offset, [on, off ...]] or null, marker, ms, mfc,
              mec, mew, drawstyle, solid_capstyle
  scatter     offsets, sizes, facecolors, edgecolors, linewidths, marker ("o", "s" or the path)
  segments    segments [[[x, y] ...] ...], colors, linewidths, capstyle
  image       shape, extent, origin, interpolation, aspect, x_edges (MDBridge's column edges of a
              heatmap, RN-01: improvements_render.py), the array (2-D data: "<id>.data"
              float64 with NaN for masked cells; RGBA data: none) and what it shows
              ("<id>.rgba", uint8, row 0 at the top as stored) in the .npz, cmap name, vmin, vmax
  text        text, xy, coordinates ("data", "axes", "figure"), ha, va, color, size, weight,
              style, family, rotation, clip_on, wrap, extent_px (its drawn box)
  legend      loc (Matplotlib's code; 0 = "best"), extent_px, frame {facecolor, edgecolor, lw},
              entries [{text, color, size, handle}] (handle: the legend's proxy line / patch)
  pixel_image the 3D layer of the structure view (mdatk/render/molecule/stage.py PixelImage):
              shape and coverage (pixels with alpha > 0); its pixels are make_molecule_cases.py's
              subject (fixed viewports), not the figures'.
  other       class name only

Floats are written by render_common.dumps (Python's repr). Colours are RGBA floats 0..1.
"""

from __future__ import annotations

import numpy as np
from matplotlib import ticker
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.colors import to_rgba
from matplotlib.image import AxesImage
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text


def rgba(c):
    if c is None:
        return None
    if isinstance(c, str) and c.lower() == "none":
        return "none"
    a = np.asarray(c)
    if a.ndim == 2:
        return [list(map(float, row)) for row in a]
    return [float(v) for v in to_rgba(c)]


def box(b):
    return [float(b.x0), float(b.y0), float(b.x1), float(b.y1)]


def num(x):
    return float(x)


class Walker:
    """One figure's walk: ids for every artist and the arrays for the case's .npz."""

    def __init__(self, fig, prefix: str = ""):
        self.fig = fig
        self.renderer = fig.canvas.get_renderer()
        self.ids = {}       # id(artist) -> "ax<k>/<n>"
        self.arrays = {}    # npz member -> array
        self.prefix = prefix

    # --- identifiers ------------------------------------------------------------------------------
    def register(self):
        for k, ax in enumerate(self.fig.axes):
            self.ids[id(ax)] = f"ax{k}"
            for n, art in enumerate(ax._children):
                self.ids[id(art)] = f"ax{k}/{n}"
            if ax.get_legend() is not None:
                self.ids[id(ax.get_legend())] = f"ax{k}/legend"
            for loc, t in (("left", ax._left_title), ("center", ax.title), ("right", ax._right_title)):
                self.ids[id(t)] = f"ax{k}/title_{loc}"
        for n, t in enumerate(self.fig.texts):
            self.ids[id(t)] = f"fig/text{n}"
        if self.fig._suptitle is not None:
            self.ids[id(self.fig._suptitle)] = "fig/suptitle"

    def id_of(self, art) -> str:
        return self.ids.get(id(art), f"?{type(art).__name__}")

    def put(self, name, array):
        self.arrays[self.prefix + name] = np.ascontiguousarray(array)

    # --- figure -----------------------------------------------------------------------------------
    def figure(self, roles=None) -> dict:
        """The whole figure; roles: id(axes) -> {"panel": k, "role": ...}."""
        self.register()
        fig = self.fig
        w, h = fig.canvas.get_width_height()
        out = {"size_px": [w, h], "dpi": num(fig.dpi), "facecolor": rgba(fig.get_facecolor()),
               "suptitle": self.text(fig._suptitle) if fig._suptitle is not None else None,
               "texts": [self.text(t) for t in fig.texts], "axes": []}
        for ax in fig.axes:
            rec = self.axes(ax)
            if roles and id(ax) in roles:
                rec.update(roles[id(ax)])
            out["axes"].append(rec)
        return out

    def axes(self, ax) -> dict:
        r = self.renderer
        rec = {"id": self.id_of(ax), "position_px": box(ax.get_window_extent(r)),
               "tight_px": box(ax.get_tightbbox(r)) if ax.get_visible() else None,
               "facecolor": rgba(ax.get_facecolor()), "visible": ax.get_visible(), "axison": ax.axison,
               "xlim": [num(v) for v in ax.get_xlim()], "ylim": [num(v) for v in ax.get_ylim()],
               "aspect": ax.get_aspect() if isinstance(ax.get_aspect(), str) else num(ax.get_aspect()),
               "zorder": num(ax.get_zorder()),
               "titles": {loc: self.text(t) for loc, t in (("left", ax._left_title), ("center", ax.title),
                                                             ("right", ax._right_title)) if t.get_text()},
               "spines": {side: {"visible": s.get_visible(), "color": rgba(s.get_edgecolor()),
                                 "lw": num(s.get_linewidth())} for side, s in ax.spines.items()},
               "xaxis": self.axis(ax.xaxis), "yaxis": self.axis(ax.yaxis)}
        cb = getattr(ax, "_colorbar", None)
        if cb is not None:
            rec["colorbar"] = self.colorbar(cb)
        rec["artists"] = [self.artist(a) for a in ax._children]
        if ax.get_legend() is not None:
            rec["artists"].append(self.artist(ax.get_legend()))
        return rec

    def axis(self, axis) -> dict:
        ticks = axis._update_ticks()
        major = axis.majorTicks[0] if axis.majorTicks else None
        tick_style = None
        if major is not None:
            mark = major.tick1line if major.tick1line.get_visible() else major.tick2line
            label = major.label1 if major.label1.get_visible() else major.label2
            tick_style = {"tick1": major.tick1line.get_visible(), "tick2": major.tick2line.get_visible(),
                          "label1": major.label1.get_visible(), "label2": major.label2.get_visible(),
                          "color": rgba(mark.get_color()), "length": num(mark.get_markersize()),
                          "width": num(mark.get_markeredgewidth()), "labelcolor": rgba(label.get_color()),
                          "labelsize": num(label.get_fontsize()), "rotation": num(label.get_rotation()),
                          "gridOn": bool(major.gridline.get_visible()),
                          "grid": {"color": rgba(major.gridline.get_color()),
                                   "lw": num(major.gridline.get_linewidth()),
                                   "ls": major.gridline.get_linestyle()}}
        lab = axis.label
        return {"label": {"text": lab.get_text(), "color": rgba(lab.get_color()), "size": num(lab.get_fontsize()),
                          "extent_px": box(lab.get_window_extent(self.renderer)) if lab.get_text() else None},
                "scale": axis.get_scale(), "inverted": bool(axis.get_inverted()), "visible": axis.get_visible(),
                "ticks": [[num(t.get_loc()), (t.label1.get_text() if t.label1.get_visible() else t.label2.get_text()),
                           t.label1.get_visible(), t.label2.get_visible()] for t in ticks],
                "tick_style": tick_style,
                "offset_text": axis.get_offset_text().get_text(),
                "locator": self.locator(axis), "formatter": self.formatter(axis.get_major_formatter())}

    def locator(self, axis) -> dict:
        loc = axis.get_major_locator()
        out = {"class": type(loc).__name__}
        if isinstance(loc, ticker.MaxNLocator):
            out.update(nbins=loc._nbins, steps=[num(s) for s in loc._steps], integer=bool(loc._integer),
                       min_n_ticks=int(loc._min_n_ticks), symmetric=bool(loc._symmetric), prune=loc._prune)
            if loc._nbins == "auto":
                out["tick_space"] = int(axis.get_tick_space())
                out["auto_nbins"] = int(np.clip(axis.get_tick_space(), max(1, loc._min_n_ticks - 1), 9))
        elif isinstance(loc, ticker.FixedLocator):
            out.update(locs=[num(v) for v in loc.locs], nbins=loc.nbins)
        return out

    def formatter(self, fmt) -> dict:
        out = {"class": type(fmt).__name__}
        if isinstance(fmt, ticker.ScalarFormatter):
            out.update(offset=num(fmt.offset), orderOfMagnitude=int(fmt.orderOfMagnitude), format=fmt.format,
                       useMathText=bool(fmt._useMathText), useOffset=bool(fmt._useOffset))
        elif isinstance(fmt, ticker.FixedFormatter):
            out.update(seq=list(fmt.seq))
        return out

    def colorbar(self, cb) -> dict:
        cmap, norm = cb.cmap, cb.norm
        name = self.id_of(cb.ax)
        cmap._init() if not cmap._isinit else None
        self.put(f"{name}.lut", cmap(np.arange(cmap.N)))
        return {"mappable": self.id_of(cb.mappable), "orientation": cb.orientation,
                "vmin": num(norm.vmin), "vmax": num(norm.vmax), "norm": type(norm).__name__,
                "cmap": {"name": cmap.name, "N": int(cmap.N), "over": rgba(cmap.get_over()),
                         "under": rgba(cmap.get_under()), "bad": rgba(cmap.get_bad())},
                "boundaries": [num(v) for v in cb._boundaries] if cb._boundaries is not None else None,
                "values": [num(v) for v in cb._values] if cb._values is not None else None,
                "extend": cb.extend, "outline": {"color": rgba(cb.outline.get_edgecolor()),
                                                 "lw": num(cb.outline.get_linewidth())},
                "label": cb.ax.yaxis.label.get_text() if cb.orientation == "vertical" else
                cb.ax.xaxis.label.get_text()}

    # --- artists ----------------------------------------------------------------------------------
    def common(self, a, kind) -> dict:
        out = {"id": self.id_of(a), "kind": kind, "zorder": num(a.get_zorder()), "visible": bool(a.get_visible()),
               "alpha": None if a.get_alpha() is None else num(a.get_alpha()), "animated": bool(a.get_animated())}
        label = a.get_label()
        if label and not str(label).startswith("_"):
            out["label"] = str(label)
        return out

    def artist(self, a) -> dict:
        if isinstance(a, Legend):
            return self.legend(a)
        if isinstance(a, Line2D):
            return self.line(a)
        if isinstance(a, PathCollection):
            return self.scatter(a)
        if isinstance(a, LineCollection):
            return self.segments(a)
        if isinstance(a, AxesImage):
            return self.image(a)
        if isinstance(a, Text):
            return self.text(a)
        if type(a).__name__ == "PixelImage":
            return self.pixel_image(a)
        if isinstance(a, Patch):
            out = self.common(a, "patch")
            out.update({"class": type(a).__name__, "facecolor": rgba(a.get_facecolor()),
                        "edgecolor": rgba(a.get_edgecolor()), "lw": num(a.get_linewidth())})
            if isinstance(a, Rectangle):
                out["xywh"] = [num(a.get_x()), num(a.get_y()), num(a.get_width()), num(a.get_height())]
            return out
        out = self.common(a, "other")
        out["class"] = type(a).__name__
        return out

    def line(self, a, short=False) -> dict:
        out = self.common(a, "line")
        if not short:
            out["x"] = [num(v) for v in np.asarray(a.get_xdata(), dtype=np.float64).ravel()]
            out["y"] = [num(v) for v in np.asarray(a.get_ydata(), dtype=np.float64).ravel()]
        dashes = getattr(a, "_unscaled_dash_pattern", None)
        out.update({"color": rgba(a.get_color()), "lw": num(a.get_linewidth()), "ls": a.get_linestyle(),
                    "dashes": (None if not dashes or dashes[1] is None
                               else [num(dashes[0]), [num(v) for v in dashes[1]]]),
                    "marker": str(a.get_marker()), "ms": num(a.get_markersize()),
                    "mfc": rgba(a.get_markerfacecolor()), "mec": rgba(a.get_markeredgecolor()),
                    "mew": num(a.get_markeredgewidth()), "drawstyle": a.get_drawstyle(),
                    "solid_capstyle": a.get_solid_capstyle()})
        if not short:
            out["transform"] = self.coords(a)
        return out

    def scatter(self, a) -> dict:
        out = self.common(a, "scatter")
        offsets = np.asarray(a.get_offsets(), dtype=np.float64).reshape(-1, 2)
        out.update({"offsets": [[num(x), num(y)] for x, y in offsets],
                    "sizes": [num(s) for s in a.get_sizes()],
                    "facecolors": rgba(a.get_facecolors()), "edgecolors": rgba(a.get_edgecolors()),
                    "linewidths": [num(v) for v in a.get_linewidths()], "marker": self.marker(a)})
        return out

    @staticmethod
    def marker(a):
        paths = a.get_paths()
        if not paths:
            return None
        verts = paths[0].vertices
        for name in ("o", "s"):
            m = MarkerStyle(name)
            ref = m.get_path().transformed(m.get_transform()).vertices
            if ref.shape == verts.shape and np.allclose(ref, verts):
                return name
        return [[num(x), num(y)] for x, y in verts]

    def segments(self, a) -> dict:
        out = self.common(a, "segments")
        out.update({"segments": [[[num(x), num(y)] for x, y in np.asarray(s, dtype=np.float64).reshape(-1, 2)]
                                 for s in a.get_segments()],
                    "colors": rgba(a.get_edgecolor()), "linewidths": [num(v) for v in a.get_linewidths()],
                    "capstyle": a.get_capstyle(), "linestyles": [str(ls) for ls in a.get_linestyles()]})
        return out

    def image(self, a) -> dict:
        out = self.common(a, "image")
        arr = a.get_array()
        data = np.ma.getdata(arr)
        out.update({"shape": list(data.shape), "extent": [num(v) for v in a.get_extent()], "origin": a.origin,
                    "interpolation": a.get_interpolation(), "dtype": str(data.dtype)})
        if data.ndim == 2:
            values = np.asarray(data, dtype=np.float64).copy()
            mask = np.ma.getmaskarray(arr)
            values[mask] = np.nan
            self.put(f"{out['id']}.data", values)
            out.update({"cmap": a.get_cmap().name, "vmin": num(a.norm.vmin), "vmax": num(a.norm.vmax),
                        "masked": int(mask.sum())})
        edges = getattr(a, "_mdbridge_x_edges", None)  # columns of their own widths (RN-01)
        if edges is not None:
            out["x_edges"] = [num(v) for v in edges]
        self.put(f"{out['id']}.rgba", a.to_rgba(arr, bytes=True))
        return out

    def coords(self, a) -> str:
        t = a.get_transform()
        ax = a.axes
        if ax is not None and t == ax.transData:
            return "data"
        if ax is not None and t == ax.transAxes:
            return "axes"
        if t == self.fig.transFigure:
            return "figure"
        if ax is not None and t == ax.get_xaxis_transform():
            return "xaxis"  # x data, y axes (axvline)
        if ax is not None and t == ax.get_yaxis_transform():
            return "yaxis"  # x axes, y data (axhline)
        return "other"

    def text(self, a) -> dict:
        out = self.common(a, "text")
        x, y = a.get_position()
        shown = bool(a.get_text()) and a.get_visible()
        out.update({"text": a.get_text(), "xy": [num(x), num(y)], "coords": self.coords(a),
                    "ha": a.get_horizontalalignment(), "va": a.get_verticalalignment(),
                    "color": rgba(a.get_color()), "size": num(a.get_fontsize()), "weight": str(a.get_fontweight()),
                    "style": a.get_fontstyle(), "family": list(a.get_fontfamily()), "rotation": num(a.get_rotation()),
                    "clip_on": bool(a.get_clip_on()), "wrap": bool(a.get_wrap()),
                    "extent_px": box(a.get_window_extent(self.renderer)) if shown else None})
        return out

    def legend(self, leg) -> dict:
        out = self.common(leg, "legend")
        frame = leg.get_frame()
        entries = []
        for text, handle in zip(leg.get_texts(), leg.legend_handles):
            e = {"text": text.get_text(), "color": rgba(text.get_color()), "size": num(text.get_fontsize())}
            if isinstance(handle, Line2D):
                e["handle"] = self.line(handle, short=True)
            elif handle is not None:
                e["handle"] = {"kind": type(handle).__name__,
                               "facecolor": rgba(handle.get_facecolor()) if hasattr(handle, "get_facecolor") else None}
            entries.append(e)
        out.update({"loc": int(leg._loc), "extent_px": box(leg.get_window_extent(self.renderer)),
                    "frame": {"on": bool(leg.get_frame_on()), "facecolor": rgba(frame.get_facecolor()),
                              "edgecolor": rgba(frame.get_edgecolor()), "lw": num(frame.get_linewidth())},
                    "entries": entries})
        return out

    def pixel_image(self, a) -> dict:
        out = self.common(a, "pixel_image")
        img = a.rgba
        if img is None:
            out.update({"shape": None})
            return out
        img = np.ascontiguousarray(img)
        out.update({"shape": list(img.shape), "dtype": str(img.dtype), "coverage": int((img[..., 3] > 0).sum())})
        return out
