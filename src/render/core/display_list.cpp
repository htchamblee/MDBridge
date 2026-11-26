// A drawn figure as painting primitives (display_list.hpp), in Matplotlib's draw order: the figure
// patch, then axes and figure texts by zorder; per axes its background patch, then its artists,
// spines, axes (grid lines, ticks, tick labels, axis label, offset text), titles and legend sorted
// by zorder, ties in that order (Matplotlib 3.11.2, figure.py / axes/_base.py draw; licence notice
// in figure.hpp).

#include "mdbridge/render/display_list.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

namespace mdbridge::render {

namespace {

// An image whose columns have their own widths (Image::xEdges): resampled to one source column per
// device pixel column (nearest neighbour at the pixel centre, as the raster backends draw images)
ImageItem resampledColumns(const Image& a, ImageItem im) {
    const double x0 = im.corner0.x, x1 = im.corner1.x;
    const double e0 = a.xEdges.front(), e1 = a.xEdges.back();
    const int width = std::max(1, static_cast<int>(std::ceil(std::fabs(x1 - x0))));
    if (!(e1 != e0) || !std::isfinite(x0) || !std::isfinite(x1) || width > (1 << 16)) return im;
    std::vector<std::uint8_t> out(static_cast<std::size_t>(width) * static_cast<std::size_t>(a.rows) * 4);
    std::vector<int> source(static_cast<std::size_t>(width));
    for (int px = 0; px < width; ++px) {
        const double x = e0 + (e1 - e0) * (px + 0.5) / width;  // data x at the pixel column's centre
        // the column whose [edge i, edge i + 1) holds x (edges increasing or decreasing)
        int c = 0;
        if (e1 > e0) {
            c = static_cast<int>(std::upper_bound(a.xEdges.begin(), a.xEdges.end(), x) - a.xEdges.begin()) - 1;
        } else {
            c = static_cast<int>(std::upper_bound(a.xEdges.begin(), a.xEdges.end(), x, std::greater<>()) - a.xEdges.begin()) - 1;
        }
        source[static_cast<std::size_t>(px)] = std::clamp(c, 0, a.cols - 1);
    }
    for (int r = 0; r < a.rows; ++r) {
        for (int px = 0; px < width; ++px) {
            const std::size_t from = (static_cast<std::size_t>(r) * static_cast<std::size_t>(a.cols) +
                                      static_cast<std::size_t>(source[static_cast<std::size_t>(px)])) * 4;
            const std::size_t to = (static_cast<std::size_t>(r) * static_cast<std::size_t>(width) + static_cast<std::size_t>(px)) * 4;
            std::copy_n(im.rgba.begin() + static_cast<std::ptrdiff_t>(from), 4, out.begin() + static_cast<std::ptrdiff_t>(to));
        }
    }
    im.width = width;
    im.rgba = std::move(out);
    return im;
}


Rgba withAlpha(Rgba c, const std::optional<double>& alpha) {
    if (alpha) c.a = *alpha;
    return c;
}

class Builder {
public:
    Builder(const Figure& fig, const FontSet& fonts, DisplayLayer layer, DisplayList& out)
        : m_fig(fig), m_fonts(fonts), m_layer(layer), m_out(out), m_h(fig.heightPx()) {}

    void build();

private:
    bool wanted(bool animated) const {
        return m_layer == DisplayLayer::All || (m_layer == DisplayLayer::Static) != animated;
    }
    Point dev(Point p) const { return {p.x, m_h - p.y}; }
    Box devBox(const Box& b) const { return {b.x0, m_h - b.y1, b.x1, m_h - b.y0}; }
    double px(double pt) const { return m_fig.pointsToPixels(pt); }

    template <class T>
    void emit(T item, std::optional<Box> clip, bool animated) {
        DisplayItem d;
        d.item = std::move(item);
        d.clip = clip;
        d.animated = animated;
        m_out.items.push_back(std::move(d));
    }

    void text(const Text& t, Point anchor, std::optional<Box> clip, bool animated);
    void line(const Axes& ax, const Line& l, std::optional<Box> clip, bool animated);
    void segment(Point a, Point b, const DisplayStroke& s, std::optional<Box> clip, bool animated) {
        PathItem p;
        p.subpaths = {{dev(a), dev(b)}};
        p.stroke = s;
        emit(std::move(p), clip, animated);
    }
    void axes(std::size_t k);
    void axis(std::size_t k, AxisDirection d, std::optional<Box> clip);
    void legend(std::size_t k);
    void colorbarSolids(std::size_t k);

    const Figure& m_fig;
    const FontSet& m_fonts;
    DisplayLayer m_layer;
    DisplayList& m_out;
    double m_h;
};

void Builder::text(const Text& t, Point anchor, std::optional<Box> clip, bool animated) {
    if (!t.visible || t.text.empty() || !wanted(animated)) return;
    const Figure::TextGeometry g = m_fig.textGeometry(t, anchor, m_fonts);
    for (std::size_t i = 0; i < g.lines.size(); ++i) {
        if (g.lines[i].first.empty()) continue;
        TextItem it;
        it.text = g.lines[i].first;
        it.face = g.face;
        it.pixelSize = g.pixelSize;
        it.origin = dev(g.lines[i].second);
        it.rotation = t.rotation;
        it.color = withAlpha(t.color, t.alpha);
        it.layout = g.layouts[i];
        emit(std::move(it), t.clipOn ? clip : std::nullopt, animated);
    }
}

void Builder::line(const Axes& ax, const Line& l, std::optional<Box> clip, bool animated) {
    if (!l.visible || !wanted(animated)) return;
    const std::size_t n = std::min(l.x.size(), l.y.size());
    std::vector<Point> pts;
    pts.reserve(n);
    for (std::size_t i = 0; i < n; ++i) pts.push_back(m_fig.toPixels(ax, l.coords, {l.x[i], l.y[i]}));
    const std::optional<Box> c = l.clipOn ? clip : std::nullopt;
    if (l.drawsLine() && l.width > 0) {
        PathItem p;
        std::vector<Point> cur;
        for (const Point& q : pts) {
            if (std::isfinite(q.x) && std::isfinite(q.y)) {
                cur.push_back(dev(q));
            } else if (!cur.empty()) {
                p.subpaths.push_back(std::move(cur));
                cur.clear();
            }
        }
        if (!cur.empty()) p.subpaths.push_back(std::move(cur));
        DisplayStroke s;
        s.color = withAlpha(l.color, l.alpha);
        s.width = px(l.width);
        s.cap = l.capStyle;
        const Dashes d = l.scaledDashes();
        for (double v : d.pattern) s.dashes.push_back(px(v));
        s.dashOffset = px(d.offset);
        p.stroke = s;
        if (!p.subpaths.empty()) emit(std::move(p), c, animated);
    }
    if (l.marker != Marker::None) {
        MarkerItem m;
        m.shape = l.marker;
        for (const Point& q : pts) {
            if (std::isfinite(q.x) && std::isfinite(q.y)) m.centers.push_back(dev(q));
        }
        m.sizes = {px(l.markerSize)};
        m.faces = {withAlpha(l.markerFace.value_or(l.color), l.alpha)};
        m.edges = {withAlpha(l.markerEdge.value_or(l.color), l.alpha)};
        m.edgeWidth = px(l.markerEdgeWidth);
        if (!m.centers.empty()) emit(std::move(m), c, animated);
    }
}

void Builder::axis(std::size_t k, AxisDirection d, std::optional<Box> clip) {
    const Axes& ax = m_fig.axes[k];
    const Axis& a = d == AxisDirection::X ? ax.xaxis : ax.yaxis;
    if (!a.visible || !wanted(false)) return;
    const bool x = d == AxisDirection::X;
    const Box box = m_fig.axesBox(k);
    const TickStyle& ts = a.ticks;
    DisplayStroke grid{ts.gridColor, px(ts.gridWidth), CapStyle::Projecting, true, {}, 0.0};
    DisplayStroke tick{ts.color, px(ts.width), CapStyle::Butt, false, {}, 0.0};
    const double len = px(ts.length), pad = px(ts.length + ts.pad);
    for (const DrawnTick& t : a.drawn) {
        const Point p = m_fig.toPixels(ax, x ? Coords::XAxis : Coords::YAxis, x ? Point{t.loc, 0.0} : Point{0.0, t.loc});
        if (ts.gridOn) {
            if (x) segment({p.x, box.y0}, {p.x, box.y1}, grid, clip, false);
            else segment({box.x0, p.y}, {box.x1, p.y}, grid, clip, false);
        }
        if (ts.tick1 && len > 0) {
            if (x) segment({p.x, box.y0}, {p.x, box.y0 - len}, tick, std::nullopt, false);
            else segment({box.x0, p.y}, {box.x0 - len, p.y}, tick, std::nullopt, false);
        }
        if (ts.tick2 && len > 0) {
            if (x) segment({p.x, box.y1}, {p.x, box.y1 + len}, tick, std::nullopt, false);
            else segment({box.x1, p.y}, {box.x1 + len, p.y}, tick, std::nullopt, false);
        }
        for (int which = 1; which <= 2; ++which) {
            if ((which == 1 && !ts.label1) || (which == 2 && !ts.label2)) continue;
            Text l;
            l.text = t.label;
            l.font.size = ts.labelSize;
            l.font.family = ts.labelFamily;
            l.color = ts.labelColor;
            l.rotation = ts.rotation;
            if (x) {
                l.ha = HAlign::Center;
                l.va = which == 1 ? VAlign::Top : VAlign::Bottom;
                text(l, {p.x, which == 1 ? box.y0 - pad : box.y1 + pad}, std::nullopt, false);
            } else {
                l.ha = which == 1 ? HAlign::Right : HAlign::Left;
                l.va = VAlign::CenterBaseline;
                text(l, {which == 1 ? box.x0 - pad : box.x1 + pad, p.y}, std::nullopt, false);
            }
        }
    }
    text(a.label, m_fig.axisLabelPosition(k, d), std::nullopt, false);
    text(a.offsetText, m_fig.offsetTextPosition(k, d), std::nullopt, false);
}

void Builder::legend(std::size_t k) {
    const Axes& ax = m_fig.axes[k];
    const Legend& L = *ax.legend;
    if (!L.visible || !wanted(L.animated)) return;
    const Figure::LegendGeometry g = m_fig.legendGeometry(k);
    if (L.frameOn) {
        RoundRectItem r;
        r.rect = devBox(g.box);
        r.radius = 0.2 * px(L.fontSize);
        r.fill = L.face;
        r.stroke = DisplayStroke{L.edge, px(L.edgeWidth), CapStyle::Butt, true, {}, 0.0};
        emit(std::move(r), std::nullopt, L.animated);
    }
    for (std::size_t i = 0; i < L.entries.size() && i < g.items.size(); ++i) {
        const Figure::LegendGeometry::Item& it = g.items[i];
        const Line& h = L.entries[i].handle;
        const double y = it.handleOrigin.y + (it.handleHeight - it.handleDescent) / 2.0;
        Line proxy = h;
        proxy.coords = Coords::Pixels;
        proxy.x = {it.handleOrigin.x, it.handleOrigin.x + it.handleWidth / 2.0, it.handleOrigin.x + it.handleWidth};
        proxy.y = {y, y, y};
        proxy.animated = L.animated;
        proxy.visible = true;
        // the line over the whole handle, the marker once in the middle (numpoints 1)
        Line marker = proxy;
        proxy.marker = Marker::None;
        marker.style = LineStyle::None;
        marker.x = {proxy.x[1]};
        marker.y = {y};
        line(ax, proxy, std::nullopt, L.animated);
        line(ax, marker, std::nullopt, L.animated);
        Text t;
        t.text = L.entries[i].text;
        t.font.size = L.textSize;
        t.color = L.textColor;
        t.ha = HAlign::Left;
        t.va = VAlign::Baseline;
        text(t, it.textOrigin, std::nullopt, L.animated);
    }
}

void Builder::colorbarSolids(std::size_t k) {
    const Axes& ax = m_fig.axes[k];
    const ColorbarInfo& c = *ax.colorbar;
    if (!c.cmap || !wanted(false)) return;
    const int n = c.cmap->size();
    for (int i = 0; i < n; ++i) {
        const double b0 = c.vmin + (c.vmax - c.vmin) * i / n, b1 = c.vmin + (c.vmax - c.vmin) * (i + 1) / n;
        const Point p0 = m_fig.toPixels(ax, Coords::YAxis, {0.0, b0}), p1 = m_fig.toPixels(ax, Coords::YAxis, {1.0, b1});
        const Box r = devBox(Box{std::min(p0.x, p1.x), std::min(p0.y, p1.y), std::max(p0.x, p1.x), std::max(p0.y, p1.y)});
        emit(FillRectItem{r, c.cmap->at(i)}, devBox(m_fig.axesBox(k)), false);
    }
}

void Builder::axes(std::size_t k) {
    const Axes& ax = m_fig.axes[k];
    if (!ax.visible) return;
    const Box box = m_fig.axesBox(k);
    const std::optional<Box> clip = devBox(box);
    if (ax.patchVisible && wanted(false)) emit(FillRectItem{devBox(box), ax.face}, std::nullopt, false);
    struct Entry {
        double z;
        std::function<void()> draw;
    };
    std::vector<Entry> list;
    if (ax.colorbar) list.push_back({0.9, [this, k] { colorbarSolids(k); }});
    for (const Artist& art : ax.artists) {
        std::visit(
            [&](const auto& a) {
                using T = std::decay_t<decltype(a)>;
                const ArtistBase& base = a;
                list.push_back({base.zorder, [this, &ax, &a, clip] {
                                    if constexpr (std::is_same_v<T, Line>) {
                                        line(ax, a, clip, a.animated);
                                    } else if constexpr (std::is_same_v<T, Text>) {
                                        text(a, m_fig.toPixels(ax, a.coords, a.xy), clip, a.animated);
                                    } else if constexpr (std::is_same_v<T, Scatter>) {
                                        if (!a.visible || !wanted(a.animated) || a.offsets.empty()) return;
                                        MarkerItem m;
                                        m.shape = a.marker;
                                        for (const Point& p : a.offsets) m.centers.push_back(dev(m_fig.toPixels(ax, Coords::Data, p)));
                                        m.sizes.clear();
                                        for (double s : a.sizes) m.sizes.push_back(px(std::sqrt(std::max(s, 0.0))));
                                        if (m.sizes.empty()) m.sizes = {px(6.0)};
                                        for (const Rgba& c : a.faces) m.faces.push_back(withAlpha(c, a.alpha));
                                        for (const Rgba& c : a.edges) m.edges.push_back(withAlpha(c, a.alpha));
                                        m.edgeWidth = px(a.lineWidths.empty() ? 1.0 : a.lineWidths[0]);
                                        emit(std::move(m), a.clipOn ? clip : std::nullopt, a.animated);
                                    } else if constexpr (std::is_same_v<T, Segments>) {
                                        if (!a.visible || !wanted(a.animated)) return;
                                        for (std::size_t i = 0; i < a.segments.size(); ++i) {
                                            PathItem p;
                                            std::vector<Point> pts;
                                            for (const Point& q : a.segments[i]) pts.push_back(dev(m_fig.toPixels(ax, Coords::Data, q)));
                                            p.subpaths = {std::move(pts)};
                                            DisplayStroke s;
                                            s.color = a.colors.empty() ? Rgba{0, 0, 0, 1} : withAlpha(a.colors[i % a.colors.size()], a.alpha);
                                            s.width = px(a.widths.empty() ? 1.5 : a.widths[i % a.widths.size()]);
                                            s.cap = a.capStyle;
                                            p.stroke = s;
                                            emit(std::move(p), a.clipOn ? clip : std::nullopt, a.animated);
                                        }
                                    } else if constexpr (std::is_same_v<T, Image>) {
                                        if (!a.visible || !wanted(a.animated) || a.rows <= 0 || a.cols <= 0) return;
                                        ImageItem im;
                                        im.width = a.cols;
                                        im.height = a.rows;
                                        im.rgba = a.pixels();
                                        im.smooth = a.interpolation != "nearest";
                                        if (a.fillsAxes) {
                                            // PixelImage.draw: 1:1 at the axes box' rounded lower left corner
                                            const Box bb = m_fig.axesBox(static_cast<std::size_t>(&ax - m_fig.axes.data()));
                                            const double x0 = std::round(bb.x0), y0 = std::round(bb.y0);
                                            im.corner0 = dev({x0, y0 + a.rows});
                                            im.corner1 = dev({x0 + a.cols, y0});
                                            emit(std::move(im), clip, a.animated);
                                            return;
                                        }
                                        const double top = a.originUpper ? a.extent[3] : a.extent[2];
                                        const double bottom = a.originUpper ? a.extent[2] : a.extent[3];
                                        im.corner0 = dev(m_fig.toPixels(ax, Coords::Data, {a.extent[0], top}));
                                        im.corner1 = dev(m_fig.toPixels(ax, Coords::Data, {a.extent[1], bottom}));
                                        if (a.xEdges.size() == static_cast<std::size_t>(a.cols) + 1) {
                                            im = resampledColumns(a, std::move(im));
                                        }
                                        emit(std::move(im), a.clipOn ? clip : std::nullopt, a.animated);
                                    } else if constexpr (std::is_same_v<T, Rectangle>) {
                                        if (!a.visible || !wanted(a.animated)) return;
                                        const Point p0 = dev(m_fig.toPixels(ax, a.coords, {a.x, a.y}));
                                        const Point p1 = dev(m_fig.toPixels(ax, a.coords, {a.x + a.w, a.y + a.h}));
                                        PathItem p;
                                        p.subpaths = {{p0, {p1.x, p0.y}, p1, {p0.x, p1.y}}};
                                        p.closed = true;
                                        if (a.face) p.fill = withAlpha(*a.face, a.alpha);
                                        if (a.edge) p.stroke = DisplayStroke{withAlpha(*a.edge, a.alpha), px(a.lineWidth), CapStyle::Butt, false, {}, 0.0};
                                        emit(std::move(p), a.clipOn ? clip : std::nullopt, a.animated);
                                    }
                                }});
            },
            art);
    }
    // spines (Matplotlib's order: left, right, bottom, top), zorder 2.5
    for (SpineSide s : {kLeft, kRight, kBottom, kTop}) {
        const Spine& sp = ax.spines[s];
        if (!sp.visible || !ax.axison) continue;
        list.push_back({2.5, [this, s, box, sp] {
                            if (!wanted(false)) return;
                            DisplayStroke st{sp.color, px(sp.width), CapStyle::Projecting, false, {}, 0.0};
                            switch (s) {
                                case kLeft: segment({box.x0, box.y0}, {box.x0, box.y1}, st, std::nullopt, false); break;
                                case kRight: segment({box.x1, box.y0}, {box.x1, box.y1}, st, std::nullopt, false); break;
                                case kBottom: segment({box.x0, box.y0}, {box.x1, box.y0}, st, std::nullopt, false); break;
                                case kTop: segment({box.x0, box.y1}, {box.x1, box.y1}, st, std::nullopt, false); break;
                            }
                        }});
    }
    if (ax.colorbar && ax.colorbar->outline.visible) {
        const Spine sp = ax.colorbar->outline;
        list.push_back({2.5, [this, box, sp] {
                            if (!wanted(false)) return;
                            PathItem p;
                            p.subpaths = {{dev({box.x0, box.y0}), dev({box.x1, box.y0}), dev({box.x1, box.y1}), dev({box.x0, box.y1})}};
                            p.closed = true;
                            p.stroke = DisplayStroke{sp.color, px(sp.width), CapStyle::Butt, false, {}, 0.0};
                            emit(std::move(p), std::nullopt, false);
                        }});
    }
    if (ax.axison) {
        const double az = ax.axisBelow == AxisBelow::True ? 0.5 : 1.5;
        list.push_back({az, [this, k, clip] { axis(k, AxisDirection::X, clip); }});
        list.push_back({az, [this, k, clip] { axis(k, AxisDirection::Y, clip); }});
    }
    for (int t : {1, 0, 2}) {  // title, left title, right title
        const Text& title = ax.titles[static_cast<std::size_t>(t)];
        list.push_back({title.zorder, [this, k, t, &title] { text(title, m_fig.titlePosition(k, t), std::nullopt, title.animated); }});
    }
    if (ax.legend) list.push_back({ax.legend->zorder, [this, k] { legend(k); }});
    std::stable_sort(list.begin(), list.end(), [](const Entry& a, const Entry& b) { return a.z < b.z; });
    for (const Entry& e : list) e.draw();
}

void Builder::build() {
    m_out.width = m_fig.widthPx();
    m_out.height = m_fig.heightPx();
    m_out.dpi = m_fig.dpi();
    if (wanted(false)) emit(FillRectItem{devBox(m_fig.figureBox()), m_fig.facecolor}, std::nullopt, false);
    struct Entry {
        double z;
        std::function<void()> draw;
    };
    std::vector<Entry> list;
    for (std::size_t k = 0; k < m_fig.axes.size(); ++k) list.push_back({m_fig.axes[k].zorder, [this, k] { axes(k); }});
    if (m_fig.suptitle) {
        const Text& t = *m_fig.suptitle;
        list.push_back({t.zorder, [this, &t] { text(t, m_fig.suptitlePosition(), std::nullopt, t.animated); }});
    }
    for (const Text& t : m_fig.texts) {
        list.push_back({t.zorder, [this, &t] {
                            const Point at = t.coords == Coords::Pixels ? t.xy : Point{t.xy.x * m_fig.widthPx(), t.xy.y * m_fig.heightPx()};
                            text(t, at, std::nullopt, t.animated);
                        }});
    }
    std::stable_sort(list.begin(), list.end(), [](const Entry& a, const Entry& b) { return a.z < b.z; });
    for (const Entry& e : list) e.draw();
}

}  // namespace

DisplayList buildDisplayList(const Figure& figure, const FontSet& fonts, DisplayLayer layer) {
    DisplayList out;
    Builder(figure, fonts, layer, out).build();
    return out;
}

}  // namespace mdbridge::render
