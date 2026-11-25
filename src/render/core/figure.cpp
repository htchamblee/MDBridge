// The figure model: artists, axes, styles and the text geometry (figure.hpp). The layout and the
// rest of a draw are in figure_layout.cpp.
//
// Ported from Matplotlib 3.11.2 (text.py Text._get_layout / get_window_extent, lines.py dash
// scaling, image.py colour mapping, gridspec.py, colorbar.py make_axes, axes/_base.py twinx) under
// Matplotlib's licence (packaging/licenses/Matplotlib-LICENSE.txt):
// Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved. Changes: translated to C++
// for the artists and options MDATK uses (summary in figure.hpp).

#include "mdbridge/render/figure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace mdbridge::render {

bool Box::isFinite() const {
    const double w = width(), h = height();
    return std::isfinite(w) && std::isfinite(h) && (w != 0.0 || h != 0.0) && std::isfinite(x0) &&
           std::isfinite(y0);
}

Box unionOf(const std::vector<Box>& boxes) {
    if (boxes.empty()) {
        throw std::invalid_argument("'bboxes' cannot be empty");
    }
    Box out = boxes.front();
    for (const Box& b : boxes) {
        out.x0 = std::min({out.x0, b.x0, b.x1});
        out.x1 = std::max({out.x1, b.x0, b.x1});
        out.y0 = std::min({out.y0, b.y0, b.y1});
        out.y1 = std::max({out.y1, b.y0, b.y1});
    }
    return out;
}

Dashes Line::unscaledDashes() const {
    Dashes d;
    switch (style) {
        case LineStyle::Solid:
        case LineStyle::None:
            break;
        case LineStyle::Dashed:
            d.pattern = {3.7, 1.6};
            break;
        case LineStyle::DashDot:
            d.pattern = {6.4, 1.6, 1.0, 1.6};
            break;
        case LineStyle::Dotted:
            d.pattern = {1.0, 1.65};
            break;
        case LineStyle::Custom:
            d = dashes;
            break;
    }
    return d;
}

Dashes Line::scaledDashes() const {
    // Matplotlib's _scale_dashes: offset and pattern times the width
    Dashes d = unscaledDashes();
    d.offset *= width;
    for (double& v : d.pattern) {
        v *= width;
    }
    return d;
}

std::vector<std::uint8_t> Image::pixels() const {
    const std::size_t n = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    if (!cmap) {
        return rgba;
    }
    std::vector<std::uint8_t> out(n * 4);
    // Normalize(vmin, vmax): (v - vmin) / (vmax - vmin); vmin == vmax maps everything to 0
    const double span = vmax - vmin;
    for (std::size_t i = 0; i < n && i < values.size(); ++i) {
        const double v = values[i];
        std::array<std::uint8_t, 4> c{};
        if (!std::isfinite(v)) {
            c = cmap->bytes(std::numeric_limits<double>::quiet_NaN(), alpha);
        } else {
            const double x = span == 0.0 ? 0.0 : (v - vmin) / span;
            c = cmap->bytes(x, alpha);
        }
        std::copy(c.begin(), c.end(), out.begin() + static_cast<std::ptrdiff_t>(4 * i));
    }
    return out;
}

void styleLegend(Legend& legend, const LegendStyle& style) {
    legend.face = style.face;
    legend.edge = style.edge;
    legend.edgeWidth = style.edgeWidth;
    legend.textColor = style.textColor;
    legend.textSize = style.textSize;
}

// ---- axes ----------------------------------------------------------------------------------------

namespace {

Axis makeAxis(AxisDirection d) {
    Axis a;
    a.direction = d;
    a.label.font.size = 10.0;
    a.label.clipOn = false;
    a.offsetText.font.size = 10.0;
    a.offsetText.clipOn = false;
    if (d == AxisDirection::X) {
        a.label.va = VAlign::Top;
        a.label.ha = HAlign::Center;
        a.offsetText.va = VAlign::Top;
        a.offsetText.ha = HAlign::Right;
    } else {
        a.label.va = VAlign::Bottom;
        a.label.ha = HAlign::Center;
        a.label.rotation = 90.0;
        a.label.rotationAnchor = true;
        a.offsetText.va = VAlign::Baseline;
        a.offsetText.ha = HAlign::Left;
    }
    return a;
}

}  // namespace

Axes::Axes() : xaxis(makeAxis(AxisDirection::X)), yaxis(makeAxis(AxisDirection::Y)) {
    for (int k = 0; k < 3; ++k) {
        titles[static_cast<std::size_t>(k)].ha = k == 0 ? HAlign::Left : k == 1 ? HAlign::Center : HAlign::Right;
        titles[static_cast<std::size_t>(k)].va = VAlign::Baseline;
        titles[static_cast<std::size_t>(k)].font.size = 12.0;  // axes.titlesize "large"
        titles[static_cast<std::size_t>(k)].clipOn = true;
        titles[static_cast<std::size_t>(k)].xy = {k * 0.5, 1.0};
        titles[static_cast<std::size_t>(k)].coords = Coords::Axes;
    }
}

void Axes::setXLim(double lo, double hi) {
    xaxis.limits = setLimits(lo, hi);
    xaxis.autoscale = false;
}

void Axes::setYLim(double lo, double hi) {
    yaxis.limits = setLimits(lo, hi);
    yaxis.autoscale = false;
}

namespace {
void setTicks(Axis& axis, std::vector<double> locs, std::optional<std::vector<std::string>> labels) {
    axis.locator.kind = Locator::Kind::Fixed;
    axis.locator.locs = std::move(locs);
    axis.locator.fixedBins.reset();
    if (labels) {
        axis.formatter.kind = Formatter::Kind::Fixed;
        axis.formatter.seq = std::move(*labels);
    }
}
}  // namespace

void Axes::setXTicks(std::vector<double> locs, std::optional<std::vector<std::string>> labels) {
    setTicks(xaxis, std::move(locs), std::move(labels));
}
void Axes::setYTicks(std::vector<double> locs, std::optional<std::vector<std::string>> labels) {
    setTicks(yaxis, std::move(locs), std::move(labels));
}

Text& Axes::setTitle(std::string text, HAlign loc) {
    Text& t = titles[loc == HAlign::Left ? 0 : loc == HAlign::Center ? 1 : 2];
    t.text = std::move(text);
    return t;
}

void Axes::DataLimits::add(double v) {
    if (std::isfinite(v)) {
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
}

void Axes::updateDataLimits(const Artist& art) {
    if (const Line* l = std::get_if<Line>(&art)) {
        // update_from_path: vertices with a non-finite coordinate are skipped; blended transforms
        // (axvline, axhline) update their data direction only
        const std::size_t n = std::min(l->x.size(), l->y.size());
        for (std::size_t i = 0; i < n; ++i) {
            if (l->coords == Coords::Data && std::isfinite(l->x[i]) && std::isfinite(l->y[i])) {
                dataX.add(l->x[i]);
                dataY.add(l->y[i]);
            } else if (l->coords == Coords::XAxis) {
                dataX.add(l->x[i]);
            } else if (l->coords == Coords::YAxis) {
                dataY.add(l->y[i]);
            }
        }
    } else if (const Scatter* s = std::get_if<Scatter>(&art)) {
        for (const Point& p : s->offsets) {
            if (std::isfinite(p.x) && std::isfinite(p.y)) {
                dataX.add(p.x);
                dataY.add(p.y);
            }
        }
    } else if (const Segments* g = std::get_if<Segments>(&art)) {
        for (const auto& seg : g->segments) {
            for (const Point& p : seg) {
                if (std::isfinite(p.x) && std::isfinite(p.y)) {
                    dataX.add(p.x);
                    dataY.add(p.y);
                }
            }
        }
    } else if (const Image* im = std::get_if<Image>(&art)) {
        dataX.add(im->extent[0]);
        dataX.add(im->extent[1]);
        dataY.add(im->extent[2]);
        dataY.add(im->extent[3]);
        dataX.sticky.insert(dataX.sticky.end(), {im->extent[0], im->extent[1]});
        dataY.sticky.insert(dataY.sticky.end(), {im->extent[2], im->extent[3]});
    } else if (const Rectangle* r = std::get_if<Rectangle>(&art)) {
        if (r->coords == Coords::Data) {
            dataX.add(r->x);
            dataX.add(r->x + r->w);
            dataY.add(r->y);
            dataY.add(r->y + r->h);
        }
    }
}

std::size_t Axes::add(Artist artist) {
    updateDataLimits(artist);
    if (const Image* im = std::get_if<Image>(&artist)) {
        // AxesImage.set_extent: autoscaled axes take the extent as their limits, in its order (an
        // "upper" origin image inverts the y axis)
        if (xaxis.autoscale) xaxis.limits = {im->extent[0], im->extent[1]};
        if (yaxis.autoscale) yaxis.limits = {im->extent[2], im->extent[3]};
    }
    artists.push_back(std::move(artist));
    return artists.size() - 1;
}

void Axes::relim() {
    dataX = DataLimits{};
    dataY = DataLimits{};
    for (const Artist& a : artists) updateDataLimits(a);
}

void Axes::clear() {
    artists.clear();
    dataX = DataLimits{};
    dataY = DataLimits{};
    legend.reset();
    for (Axis* a : {&xaxis, &yaxis}) {
        Axis fresh = makeAxis(a->direction);
        fresh.labelOpposite = a->labelOpposite;
        fresh.ticks = a->ticks;  // Matplotlib's cla keeps the tick parameters
        fresh.label.font.size = a->label.font.size;
        fresh.label.color = a->label.color;
        *a = std::move(fresh);
    }
    for (Text& t : titles) {
        t.text.clear();
    }
}

void styleAxes(Axes& ax, const AxesStyle& s) {
    ax.face = s.face;
    ax.spines[kLeft] = {s.left.visible, s.left.color, s.left.width};
    ax.spines[kBottom] = {s.bottom.visible, s.bottom.color, s.bottom.width};
    ax.spines[kRight] = {s.right.visible, s.right.color, s.right.width};
    ax.spines[kTop] = {s.top.visible, s.top.color, s.top.width};
    for (Axis* a : {&ax.xaxis, &ax.yaxis}) {
        tickParams(*a, {s.tickColor, s.tickLabelColor, s.tickLabelSize, s.tickWidth, s.tickLength, std::nullopt});
        a->label.color = s.axisLabelColor;
        a->label.font.size = s.axisLabelSize;
        a->ticks.gridOn = s.grid;
        if (s.grid) {
            a->ticks.gridColor = s.gridColor;
            a->ticks.gridWidth = s.gridWidth;
        }
    }
    for (Text& t : ax.titles) {
        t.color = s.titleColor;
        t.font.size = s.titleSize;
    }
    // ax.title only (the left / right titles keep their own size until set_title(loc=) sets them)
    ax.axisBelow = s.grid ? AxisBelow::True : ax.axisBelow;
}

void tickParams(Axis& axis, const TickParams& p) {
    if (p.color) axis.ticks.color = *p.color;
    if (p.labelColor) {
        axis.ticks.labelColor = *p.labelColor;
        axis.offsetText.color = *p.labelColor;
    }
    if (p.labelSize) axis.ticks.labelSize = *p.labelSize;
    if (p.width) axis.ticks.width = *p.width;
    if (p.length) axis.ticks.length = *p.length;
    if (p.rotation) axis.ticks.rotation = *p.rotation;
}

// ---- figure --------------------------------------------------------------------------------------

Figure::Figure(int widthPx, int heightPx, double dpi) : m_width(widthPx), m_height(heightPx), m_dpi(dpi) {
    if (widthPx < 1 || heightPx < 1 || !(dpi > 0.0) || !std::isfinite(dpi)) {
        throw std::invalid_argument("figure size and dpi must be positive");
    }
}

std::size_t Figure::addSubplot(SubplotSpec spec) {
    if (spec.row0 < 0 || spec.row1 > grid.rows || spec.row0 >= spec.row1 || spec.col0 < 0 || spec.col1 > grid.cols ||
        spec.col0 >= spec.col1) {
        throw std::invalid_argument("subplot cells outside the grid");
    }
    Axes ax;
    ax.spec = spec;
    axes.push_back(std::move(ax));
    m_positioned = false;
    return axes.size() - 1;
}

std::size_t Figure::twinx(std::size_t main) {
    Axes twin;
    twin.spec = axes.at(main).spec;
    twin.twinOf = main;
    twin.patchVisible = false;
    twin.xaxis.visible = false;
    twin.xaxis.limits = axes[main].xaxis.limits;
    twin.xaxis.autoscale = axes[main].xaxis.autoscale;
    twin.yaxis.ticks.tick1 = false;
    twin.yaxis.ticks.tick2 = true;
    twin.yaxis.ticks.label1 = false;
    twin.yaxis.ticks.label2 = true;
    twin.yaxis.labelOpposite = true;
    twin.yaxis.label.va = VAlign::Top;
    twin.yaxis.offsetText.ha = HAlign::Right;
    twin.yaxis.offsetText.xy.x = 1.0;
    twin.position = axes[main].position;
    twin.active = axes[main].active;
    axes.push_back(std::move(twin));
    return axes.size() - 1;
}

std::size_t Figure::colorbar(std::size_t parent, std::optional<std::size_t> mappable, Colormap cmap, double vmin,
                             double vmax, double fraction, double pad) {
    initialPositions();
    Axes cb;
    ColorbarInfo info;
    info.parent = parent;
    info.mappable = mappable;
    info.fraction = fraction;
    info.pad = pad;
    info.cmap = std::move(cmap);
    info.vmin = vmin;
    info.vmax = vmax;
    cb.colorbar = std::move(info);
    // make_axes(location="right"): split the parent's box, shrink the parent
    Box& pb = axes.at(parent).position;
    const double w = pb.width();
    const Box pb1{pb.x0, pb.y0, pb.x0 + w * (1.0 - fraction - pad), pb.y1};
    const Box pbcb{pb.x0 + w * (1.0 - fraction), pb.y0, pb.x1, pb.y1};
    pb = pb1;
    axes[parent].active = pb1;
    axes[parent].anchor = {1.0, 0.5};
    cb.position = pbcb;
    cb.active = pbcb;
    cb.anchor = {0.0, 0.5};
    cb.boxAspect = 20.0;
    // the colour bar's axes: no spines but the outline, long axis on the right without grid
    for (Spine& s : cb.spines) {
        s.visible = false;
    }
    cb.colorbar->outline = Spine{true, Rgba{0, 0, 0, 1}, 0.8};
    cb.xaxis.locator.kind = Locator::Kind::Fixed;
    cb.xaxis.locator.locs.clear();
    cb.xaxis.limits = {0.0, 1.0};
    cb.xaxis.autoscale = false;
    cb.yaxis.ticks.tick1 = false;
    cb.yaxis.ticks.tick2 = true;
    cb.yaxis.ticks.label1 = false;
    cb.yaxis.ticks.label2 = true;
    cb.yaxis.labelOpposite = true;
    cb.yaxis.label.va = VAlign::Top;
    cb.yaxis.limits = colorbarLimits(vmin, vmax);
    cb.yaxis.autoscale = false;
    axes.push_back(std::move(cb));
    return axes.size() - 1;
}

Text& Figure::setSuptitle(std::string text, Rgba color, double size) {
    Text t;
    t.text = std::move(text);
    t.xy = {0.5, 0.98};
    t.coords = Coords::Figure;
    t.ha = HAlign::Center;
    t.va = VAlign::Top;
    t.color = color;
    t.font.size = size;  // figure.titlesize "large" unless given
    suptitle = std::move(t);
    return *suptitle;
}

Text& Figure::addText(double x, double y, std::string text) {
    Text t;
    t.text = std::move(text);
    t.xy = {x, y};
    t.coords = Coords::Figure;
    texts.push_back(std::move(t));
    return texts.back();
}

void Figure::initialPositions() {
    if (m_positioned) {
        return;
    }
    // GridSpec.get_grid_positions with the default subplot parameters (rcParams figure.subplot.*)
    const double left = 0.125, right = 0.9, bottom = 0.11, top = 0.88, wspace = 0.2, hspace = 0.2;
    const int nr = grid.rows, nc = grid.cols;
    std::vector<double> hr = grid.heightRatios.empty() ? std::vector<double>(static_cast<std::size_t>(nr), 1.0)
                                                       : grid.heightRatios;
    std::vector<double> wr = grid.widthRatios.empty() ? std::vector<double>(static_cast<std::size_t>(nc), 1.0)
                                                      : grid.widthRatios;
    const double cellH = (top - bottom) / (nr + hspace * (nr - 1));
    const double sepH = hspace * cellH;
    double sumH = 0;
    for (double r : hr) sumH += r;
    const double normH = cellH * nr / sumH;
    std::vector<double> tops(static_cast<std::size_t>(nr)), bottoms(static_cast<std::size_t>(nr));
    double acc = 0.0;
    for (int i = 0; i < nr; ++i) {
        acc += i ? sepH : 0.0;
        tops[static_cast<std::size_t>(i)] = top - acc;
        acc += hr[static_cast<std::size_t>(i)] * normH;
        bottoms[static_cast<std::size_t>(i)] = top - acc;
    }
    const double cellW = (right - left) / (nc + wspace * (nc - 1));
    const double sepW = wspace * cellW;
    double sumW = 0;
    for (double r : wr) sumW += r;
    const double normW = cellW * nc / sumW;
    std::vector<double> lefts(static_cast<std::size_t>(nc)), rights(static_cast<std::size_t>(nc));
    acc = 0.0;
    for (int i = 0; i < nc; ++i) {
        acc += i ? sepW : 0.0;
        lefts[static_cast<std::size_t>(i)] = left + acc;
        acc += wr[static_cast<std::size_t>(i)] * normW;
        rights[static_cast<std::size_t>(i)] = left + acc;
    }
    for (Axes& ax : axes) {
        if (!ax.spec) continue;
        const SubplotSpec& s = *ax.spec;
        const Box b{lefts[static_cast<std::size_t>(s.col0)], bottoms[static_cast<std::size_t>(s.row1 - 1)],
                    rights[static_cast<std::size_t>(s.col1 - 1)], tops[static_cast<std::size_t>(s.row0)]};
        ax.position = b;
        ax.active = b;
    }
    m_positioned = true;
}

Box Figure::axesBox(std::size_t a) const {
    const Box& f = axes.at(a).active;
    return {f.x0 * m_width, f.y0 * m_height, f.x1 * m_width, f.y1 * m_height};
}

Point Figure::toPixels(const Axes& ax, Coords coords, Point p) const {
    const Box& f = ax.active;
    const double bx0 = f.x0 * m_width, bx1 = f.x1 * m_width, by0 = f.y0 * m_height, by1 = f.y1 * m_height;
    auto dataX = [&](double x) {
        const Range& l = ax.xaxis.limits;
        return bx0 + (x - l.lo) / (l.hi - l.lo) * (bx1 - bx0);
    };
    auto dataY = [&](double y) {
        const Range& l = ax.yaxis.limits;
        return by0 + (y - l.lo) / (l.hi - l.lo) * (by1 - by0);
    };
    auto axesX = [&](double x) { return bx0 + x * (bx1 - bx0); };
    auto axesY = [&](double y) { return by0 + y * (by1 - by0); };
    switch (coords) {
        case Coords::Data:
            return {dataX(p.x), dataY(p.y)};
        case Coords::Axes:
            return {axesX(p.x), axesY(p.y)};
        case Coords::Figure:
            return {p.x * m_width, p.y * m_height};
        case Coords::XAxis:
            return {dataX(p.x), axesY(p.y)};
        case Coords::YAxis:
            return {axesX(p.x), dataY(p.y)};
        case Coords::Pixels:
            return p;
    }
    return p;
}

// ---- text ----------------------------------------------------------------------------------------

namespace {

Point rotatePoint(double angle, double x, double y) {
    const double a = angle * std::numbers::pi / 180.0;
    const double c = std::cos(a), s = std::sin(a);
    return {c * x - s * y, s * x + c * y};
}

std::vector<std::string> splitLines(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
        if (ch == '\n') {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    out.push_back(cur);
    return out;
}

std::vector<std::string> splitWords(const std::string& s) {
    std::vector<std::string> out{std::string()};
    for (char ch : s) {
        if (ch == ' ') {
            out.emplace_back();
        } else {
            out.back().push_back(ch);
        }
    }
    return out;
}

std::string joinWords(const std::vector<std::string>& w, std::size_t from, std::size_t to) {
    std::string out;
    for (std::size_t i = from; i < to; ++i) {
        if (i > from) out.push_back(' ');
        out += w[i];
    }
    return out;
}

// Text._get_dist_to_box
double distToBox(double rotation, double x0, double y0, const Box& fb) {
    const double rad = std::numbers::pi / 180.0;
    double h1 = 0, h2 = 0;
    if (rotation > 270) {
        const double quad = rotation - 270;
        h1 = (y0 - fb.y0) / std::cos(quad * rad);
        h2 = (fb.x1 - x0) / std::cos((90 - quad) * rad);
    } else if (rotation > 180) {
        const double quad = rotation - 180;
        h1 = (x0 - fb.x0) / std::cos(quad * rad);
        h2 = (y0 - fb.y0) / std::cos((90 - quad) * rad);
    } else if (rotation > 90) {
        const double quad = rotation - 90;
        h1 = (fb.y1 - y0) / std::cos(quad * rad);
        h2 = (x0 - fb.x0) / std::cos((90 - quad) * rad);
    } else {
        h1 = (fb.x1 - x0) / std::cos(rotation * rad);
        h2 = (fb.y1 - y0) / std::cos((90 - rotation) * rad);
    }
    return std::min(h1, h2);
}

}  // namespace

std::string Figure::wrappedText(const Text& text, Point at, const FontSet& fonts) const {
    // Text._get_wrapped_text: words joined while the line (its width rounded up) fits the distance
    // from the text's position to the figure's edge in the direction of its alignment
    if (!text.wrap) return text.text;
    const double rot = std::fmod(std::fmod(text.rotation, 360.0) + 360.0, 360.0);
    const double left = distToBox(rot, at.x, at.y, figureBox());
    const double right = distToBox(std::fmod(180.0 + rot, 360.0), at.x, at.y, figureBox());
    const double lineWidth = text.ha == HAlign::Left ? left : text.ha == HAlign::Right ? right : 2 * std::min(left, right);
    const FontFace face = chooseFace(text.font.family, text.font.style, text.font.weight);
    const double size = fontSizeInPoints(text.font.size);
    std::vector<std::string> out;
    for (const std::string& unwrapped : splitLines(text.text)) {
        std::vector<std::string> words = splitWords(unwrapped);
        std::size_t start = 0;
        while (start < words.size()) {
            const std::size_t n = words.size() - start;
            if (n == 1) {
                out.push_back(words[start++]);
                continue;
            }
            for (std::size_t i = 2; i <= n; ++i) {
                const double w = std::ceil(fonts.layout(face, joinWords(words, start, start + i), size, m_dpi).width());
                if (w > lineWidth) {
                    out.push_back(joinWords(words, start, start + i - 1));
                    start += i - 1;
                    break;
                }
                if (i == n) {
                    out.push_back(joinWords(words, start, start + i));
                    start = words.size();
                    break;
                }
            }
        }
    }
    std::string joined;
    for (std::size_t i = 0; i < out.size(); ++i) {
        if (i) joined.push_back('\n');
        joined += out[i];
    }
    return joined;
}

Figure::TextGeometry Figure::textGeometry(const Text& text, Point at, const FontSet& fonts) const {
    TextGeometry g;
    g.face = chooseFace(text.font.family, text.font.style, text.font.weight);
    const double size = fontSizeInPoints(text.font.size);
    g.pixelSize = size * m_dpi / 72.0;
    const std::vector<std::string> lines = splitLines(wrappedText(text, at, fonts));
    // a wrapped text is aligned before it is rotated (_get_wrap_line_width sets rotation_mode anchor)
    const bool rotationAnchor = text.rotationAnchor || text.wrap;
    const Font& font = fonts.font(g.face);
    const double scale = 1.0 / font.unitsPerEm() * size * m_dpi / 72.0;
    double lineGap = font.typoLineGap() * scale;
    const double minAscent = font.typoAscender() * scale;
    const double minDescent = -font.typoDescender() * scale;
    if (lines.size() == 1) {
        lineGap = 0.0;
    }
    double thisy = 0.0, baseline = 0.0, d = 0.0;
    std::vector<double> ws, ys;
    for (const std::string& line : lines) {
        double w = 0, h = 0, dd = 0;
        TextLayout run;
        if (!line.empty()) {
            run = fonts.layout(g.face, line, size, m_dpi);
            w = run.width();
            h = run.height();
            dd = run.descent();
        }
        double a = h - dd;
        a = std::max(a, minAscent) + lineGap / 2.0;
        dd = std::max(dd, minDescent) + lineGap / 2.0;
        baseline = a - thisy;
        thisy -= a;
        ws.push_back(w);
        ys.push_back(thisy);
        thisy -= dd;
        d = dd;
        g.layouts.push_back(std::move(run));
    }
    g.descent = d;
    const double width = *std::max_element(ws.begin(), ws.end());
    double xmin = 0, xmax = width, ymax = 0, ymin = ys.back() - d;
    // multialignment: Matplotlib uses the horizontal alignment
    std::vector<Point> offsetLayout;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        double x = 0;
        if (text.ha == HAlign::Center) x = width / 2 - ws[i] / 2;
        if (text.ha == HAlign::Right) x = width - ws[i];
        offsetLayout.push_back({x, ys[i]});
    }
    const std::array<Point, 4> corners{Point{xmin, ymin}, Point{xmin, ymax}, Point{xmax, ymax}, Point{xmax, ymin}};
    std::array<Point, 4> rot;
    for (std::size_t i = 0; i < 4; ++i) rot[i] = rotatePoint(text.rotation, corners[i].x, corners[i].y);
    double rxmin = rot[0].x, rxmax = rot[0].x, rymin = rot[0].y, rymax = rot[0].y;
    for (const Point& p : rot) {
        rxmin = std::min(rxmin, p.x);
        rxmax = std::max(rxmax, p.x);
        rymin = std::min(rymin, p.y);
        rymax = std::max(rymax, p.y);
    }
    const double widthRot = rxmax - rxmin, heightRot = rymax - rymin;
    double offx = 0, offy = 0;
    if (!rotationAnchor) {
        offx = text.ha == HAlign::Left ? rxmin : text.ha == HAlign::Right ? rxmax : (rxmin + rxmax) / 2;
        switch (text.va) {
            case VAlign::Bottom: offy = rymin; break;
            case VAlign::Top: offy = rymax; break;
            case VAlign::Center: offy = (rymin + rymax) / 2; break;
            case VAlign::Baseline: offy = rymin + d; break;
            case VAlign::CenterBaseline: offy = rymin + heightRot - baseline / 2; break;
        }
    } else {
        const double x0 = corners[0].x, y0 = corners[0].y, x1 = corners[2].x, y1 = corners[2].y;
        offx = text.ha == HAlign::Left ? x0 : text.ha == HAlign::Right ? x1 : (x0 + x1) / 2;
        switch (text.va) {
            case VAlign::Bottom: offy = y0; break;
            case VAlign::Top: offy = y1; break;
            case VAlign::Center: offy = (y0 + y1) / 2; break;
            case VAlign::Baseline: offy = y1 - baseline; break;
            case VAlign::CenterBaseline: offy = y1 - baseline / 2; break;
        }
        const Point r = rotatePoint(text.rotation, offx, offy);
        offx = r.x;
        offy = r.y;
    }
    g.box = Box::fromBounds(rxmin - offx + at.x, rymin - offy + at.y, widthRot, heightRot);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const Point r = rotatePoint(text.rotation, offsetLayout[i].x, offsetLayout[i].y);
        g.lines.emplace_back(lines[i], Point{r.x - offx + at.x, r.y - offy + at.y});
    }
    return g;
}

Box Figure::textExtent(const Text& text, const Axes* ax, const FontSet& fonts) const {
    Point at;
    if (ax) {
        at = toPixels(*ax, text.coords, text.xy);
    } else {
        at = text.coords == Coords::Pixels ? text.xy : Point{text.xy.x * m_width, text.xy.y * m_height};
    }
    if (!text.visible) {
        return {0, 0, 1, 1};  // Bbox.unit()
    }
    if (text.text.empty()) {
        return {at.x, at.y, at.x, at.y};
    }
    return textGeometry(text, at, fonts).box;
}

}  // namespace mdbridge::render
