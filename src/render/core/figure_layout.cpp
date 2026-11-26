// What drawing a figure decides (figure.hpp): autoscaled limits, ticks and tick labels, axis label,
// offset text and title positions, legend boxes and their "best" location, tight boxes, and the
// constrained layout of the grid.
//
// Ported from Matplotlib 3.11.2 (_constrained_layout.py, _layoutgrid.py, axis.py, spines.py,
// axes/_base.py get_tightbbox / _update_title_position / apply_aspect, legend.py, offsetbox.py
// TextArea / DrawingArea / HPacker / VPacker, transforms.py, src/_path.h path_intersects_rectangle)
// under Matplotlib's licence (packaging/licenses/Matplotlib-LICENSE.txt):
// Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved. Changes: translated to C++
// for one GridSpec, colour bars on the right and Line2D legend handles; the layout grid's
// kiwisolver system replaced by its closed-form solution (summary in figure.hpp).

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "mdbridge/render/figure.hpp"

namespace mdbridge::render {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

Box toFraction(const Box& b, double w, double h) { return {b.x0 / w, b.y0 / h, b.x1 / w, b.y1 / h}; }

// Bbox.anchored(c, container)
Box anchored(const Box& self, Point c, const Box& container) {
    const double l = container.x0, b = container.y0, w = container.width(), h = container.height();
    const double W = self.width(), H = self.height();
    const double dx = (l + c.x * (w - W)) - self.x0, dy = (b + c.y * (h - H)) - self.y0;
    return self.translated(dx, dy);
}

// Bbox.shrunk_to_aspect(box_aspect, container, fig_aspect)
Box shrunkToAspect(const Box& self, double boxAspect, const Box& container, double figAspect) {
    const double w = container.width(), h = container.height();
    double H = w * boxAspect / figAspect, W = w;
    if (H > h) {
        W = h * figAspect / boxAspect;
        H = h;
    }
    return {self.x0, self.y0, self.x0 + W, self.y0 + H};
}

// _path.h segment_intersects_rectangle / path_intersects_rectangle (filled = false), NaN vertices
// removed as PathNanRemover removes them
bool segmentIntersectsRectangle(double x1, double y1, double x2, double y2, double cx, double cy, double w, double h) {
    return std::fabs(x1 + x2 - 2.0 * cx) < std::fabs(x1 - x2) + w &&
           std::fabs(y1 + y2 - 2.0 * cy) < std::fabs(y1 - y2) + h &&
           2.0 * std::fabs((x1 - cx) * (y1 - y2) - (y1 - cy) * (x1 - x2)) <
               w * std::fabs(y1 - y2) + h * std::fabs(x1 - x2);
}

bool pathIntersectsRectangle(const std::vector<Point>& path, const Box& r) {
    std::vector<Point> pts;
    for (const Point& p : path) {
        if (std::isfinite(p.x) && std::isfinite(p.y)) pts.push_back(p);
    }
    if (pts.empty()) return false;
    const double cx = (r.x0 + r.x1) * 0.5, cy = (r.y0 + r.y1) * 0.5;
    const double w = std::fabs(r.x0 - r.x1), h = std::fabs(r.y0 - r.y1);
    double x1 = pts[0].x, y1 = pts[0].y;
    if (2.0 * std::fabs(x1 - cx) <= w && 2.0 * std::fabs(y1 - cy) <= h) return true;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        if (segmentIntersectsRectangle(x1, y1, pts[i].x, pts[i].y, cx, cy, w, h)) return true;
        x1 = pts[i].x;
        y1 = pts[i].y;
    }
    return false;
}

bool boxesOverlap(const Box& a, const Box& b) {  // count_bboxes_overlapping_bbox
    const double ax0 = std::min(a.x0, a.x1), ax1 = std::max(a.x0, a.x1), ay0 = std::min(a.y0, a.y1),
                 ay1 = std::max(a.y0, a.y1);
    const double bx0 = std::min(b.x0, b.x1), bx1 = std::max(b.x0, b.x1), by0 = std::min(b.y0, b.y1),
                 by1 = std::max(b.y0, b.y1);
    return !(bx1 <= ax0 || by1 <= ay0 || bx0 >= ax1 || by0 >= ay1);
}

Text tickLabelText(const Axis& axis, const std::string& label, bool second) {
    Text t;
    t.text = label;
    t.font.size = axis.ticks.labelSize;
    t.font.family = axis.ticks.labelFamily;
    t.color = axis.ticks.labelColor;
    t.rotation = axis.ticks.rotation;
    t.coords = Coords::Pixels;
    if (axis.direction == AxisDirection::X) {
        t.ha = HAlign::Center;  // xtick.alignment
        t.va = second ? VAlign::Bottom : VAlign::Top;
    } else {
        t.ha = second ? HAlign::Left : HAlign::Right;
        t.va = VAlign::CenterBaseline;  // ytick.alignment
    }
    return t;
}

}  // namespace

// ---- ticks and axis decorations -------------------------------------------------------------------

void Figure::computeTicks(std::size_t a, const FontSet& fonts) {
    Axes& ax = axes[a];
    const Box box = axesBox(a);
    for (Axis* axis : {&ax.xaxis, &ax.yaxis}) {
        const bool x = axis->direction == AxisDirection::X;
        const double lengthPt = (x ? box.width() : box.height()) / m_dpi * 72.0;
        axis->tickSpace = tickSpace(axis->direction, lengthPt, axis->ticks.labelSize);
        const Range lim = axis->limits;
        std::vector<double> locs;
        if (axis->locator.kind == Locator::Kind::MaxN) {
            const MaxNLocator loc(axis->locator.params);
            locs = loc.tickValues(std::min(lim.lo, lim.hi), std::max(lim.lo, lim.hi), axis->tickSpace);
        } else {
            locs = fixedLocatorTicks(axis->locator.locs, axis->locator.fixedBins);
        }
        std::vector<std::string> labels;
        axis->scalar.reset();
        if (axis->formatter.kind == Formatter::Kind::Scalar) {
            ScalarFormatter f(axis->formatter.options);
            labels = f.formatTicks(locs, lim);
            axis->scalar = f;
        } else {
            for (std::size_t i = 0; i < locs.size(); ++i) {
                labels.push_back(i < axis->formatter.seq.size() ? axis->formatter.seq[i] : std::string());
            }
        }
        axis->drawn.clear();
        const double pad = pointsToPixels(axis->ticks.length + axis->ticks.pad);
        for (std::size_t i : ticksToDraw(locs, Range{std::min(lim.lo, lim.hi), std::max(lim.lo, lim.hi)})) {
            DrawnTick t;
            t.loc = locs[i];
            t.label = labels[i];
            const Point p = toPixels(ax, x ? Coords::XAxis : Coords::YAxis, x ? Point{t.loc, 0.0} : Point{0.0, t.loc});
            Text l1 = tickLabelText(*axis, t.label, false), l2 = tickLabelText(*axis, t.label, true);
            if (x) {
                l1.xy = {p.x, box.y0 - pad};
                l2.xy = {p.x, box.y1 + pad};
            } else {
                l1.xy = {box.x0 - pad, p.y};
                l2.xy = {box.x1 + pad, p.y};
            }
            t.box1 = textExtent(l1, nullptr, fonts);
            t.box2 = textExtent(l2, nullptr, fonts);
            axis->drawn.push_back(std::move(t));
        }
        axis->offsetText.text = axis->scalar ? axis->scalar->offsetText() : std::string();
    }
}

Box Figure::spineExtent(std::size_t a, SpineSide side) const {
    const Axes& ax = axes[a];
    const Box box = axesBox(a);
    Box bb;
    switch (side) {
        case kLeft: bb = {box.x0, box.y0, box.x0, box.y1}; break;
        case kRight: bb = {box.x1, box.y0, box.x1, box.y1}; break;
        case kBottom: bb = {box.x0, box.y0, box.x1, box.y0}; break;
        case kTop: bb = {box.x0, box.y1, box.x1, box.y1}; break;
    }
    const Axis& axis = (side == kLeft || side == kRight) ? ax.yaxis : ax.xaxis;
    if (!axis.visible || axis.drawn.empty()) {
        return bb;
    }
    const double padout = pointsToPixels(axis.ticks.length);  // tickdir "out"
    Box b0 = bb;
    if (axis.ticks.tick1) {
        if (side == kLeft) b0.x0 -= padout;
        if (side == kBottom) b0.y0 -= padout;
    }
    if (axis.ticks.tick2) {
        if (side == kRight) b0.x1 += padout;
        if (side == kTop) {
            b0.y1 += padout;
            b0.y0 -= padout;  // as Matplotlib 3.11 does
        }
    }
    return unionOf({bb, b0});
}

void Figure::updateAxisLabels(std::size_t a) {
    Axes& ax = axes[a];
    Derived& d = m_derived[a];
    const Box box = axesBox(a);
    const double labelPadX = pointsToPixels(ax.xaxis.labelPad), labelPadY = pointsToPixels(ax.yaxis.labelPad);
    {  // x axis
        std::vector<Box> b1, b2;
        for (const DrawnTick& t : ax.xaxis.drawn) {
            if (ax.xaxis.ticks.label1) b1.push_back(t.box1);
            if (ax.xaxis.ticks.label2) b2.push_back(t.box2);
        }
        if (!ax.xaxis.labelOpposite) {
            std::vector<Box> all = b1;
            all.push_back(spineExtent(a, kBottom));
            d.xLabel = {box.x0 + 0.5 * box.width(), unionOf(all).y0 - labelPadX};
        } else {
            std::vector<Box> all = b2;
            all.push_back(spineExtent(a, kTop));
            d.xLabel = {box.x0 + 0.5 * box.width(), unionOf(all).y1 + labelPadX};
        }
        const double bottom = b1.empty() ? box.y0 : unionOf(b1).y0;
        d.xOffset = {box.x1, bottom - pointsToPixels(3.0)};
    }
    {  // y axis
        std::vector<Box> b1, b2;
        for (const DrawnTick& t : ax.yaxis.drawn) {
            if (ax.yaxis.ticks.label1) b1.push_back(t.box1);
            if (ax.yaxis.ticks.label2) b2.push_back(t.box2);
        }
        if (!ax.yaxis.labelOpposite) {
            std::vector<Box> all = b1;
            all.push_back(spineExtent(a, kLeft));
            d.yLabel = {unionOf(all).x0 - labelPadY, box.y0 + 0.5 * box.height()};
        } else {
            std::vector<Box> all = b2;
            all.push_back(spineExtent(a, kRight));
            d.yLabel = {unionOf(all).x1 + labelPadY, box.y0 + 0.5 * box.height()};
        }
        const double top = box.y1;  // the outline of a colour bar is the axes box
        d.yOffset = {box.x0 + ax.yaxis.offsetText.xy.x * box.width(), top + pointsToPixels(3.0)};
    }
}

Point Figure::axisLabelPosition(std::size_t a, AxisDirection d) const {
    return d == AxisDirection::X ? m_derived.at(a).xLabel : m_derived.at(a).yLabel;
}
Point Figure::offsetTextPosition(std::size_t a, AxisDirection d) const {
    return d == AxisDirection::X ? m_derived.at(a).xOffset : m_derived.at(a).yOffset;
}
Point Figure::titlePosition(std::size_t a, int which) const {
    return m_derived.at(a).titles.at(static_cast<std::size_t>(which));
}
Point Figure::suptitlePosition() const {
    return suptitle ? Point{suptitle->xy.x * m_width, suptitle->xy.y * m_height} : Point{};
}
std::optional<Box> Figure::tightBox(std::size_t a) const { return m_derived.at(a).tight; }

std::optional<Box> Figure::axisTightBox(std::size_t a, AxisDirection dir, const FontSet& fonts, bool forLayout) {
    const Axes& ax = axes[a];
    const Axis& axis = dir == AxisDirection::X ? ax.xaxis : ax.yaxis;
    if (!axis.visible) return std::nullopt;
    const Derived& d = m_derived[a];
    std::vector<Box> boxes;
    if (axis.offsetText.visible) {
        Text t = axis.offsetText;
        t.coords = Coords::Pixels;
        t.xy = dir == AxisDirection::X ? d.xOffset : d.yOffset;
        boxes.push_back(textExtent(t, nullptr, fonts));
    }
    for (const DrawnTick& t : axis.drawn) {
        if (axis.ticks.label1) boxes.push_back(t.box1);
    }
    for (const DrawnTick& t : axis.drawn) {
        if (axis.ticks.label2) boxes.push_back(t.box2);
    }
    if (axis.label.visible) {
        Text t = axis.label;
        t.coords = Coords::Pixels;
        t.xy = dir == AxisDirection::X ? d.xLabel : d.yLabel;
        Box bb = textExtent(t, nullptr, fonts);
        if (forLayout) {
            if (dir == AxisDirection::X && bb.width() > 0) {
                bb.x0 = (bb.x0 + bb.x1) / 2 - 0.5;
                bb.x1 = bb.x0 + 1.0;
            }
            if (dir == AxisDirection::Y && bb.height() > 0) {
                bb.y0 = (bb.y0 + bb.y1) / 2 - 0.5;
                bb.y1 = bb.y0 + 1.0;
            }
        }
        boxes.push_back(bb);
    }
    std::vector<Box> finite;
    for (const Box& b : boxes) {
        if (b.isFinite()) finite.push_back(b);
    }
    if (finite.empty()) return std::nullopt;
    return unionOf(finite);
}

// Axes._update_title_position
void Figure::updateTitles(std::size_t a, const FontSet& fonts) {
    Axes& ax = axes[a];
    Derived& d = m_derived[a];
    const Box box = axesBox(a);
    const double pad = pointsToPixels(6.0);  // axes.titlepad
    std::vector<std::size_t> siblings{a};
    for (std::size_t k = 0; k < axes.size(); ++k) {
        if (k == a) continue;
        const bool twin = (axes[k].twinOf && *axes[k].twinOf == a) || (ax.twinOf && (*ax.twinOf == k ||
                                                                                  (axes[k].twinOf && *axes[k].twinOf == *ax.twinOf)));
        if (twin) siblings.push_back(k);
    }
    double top = -kInf;
    for (std::size_t s : siblings) {
        const Axes& sa = axes[s];
        std::optional<Box> bb;
        if (sa.xaxis.ticks.label2 || sa.xaxis.labelOpposite) {
            bb = axisTightBox(s, AxisDirection::X, fonts, false);
        }
        if (!bb) bb = axesBox(s);
        top = std::max(top, bb->y1);
    }
    auto extent = [&](const Text& t, double y) {
        Text c = t;
        c.coords = Coords::Pixels;
        c.xy = {box.x0 + t.xy.x * box.width(), y};
        return textExtent(c, nullptr, fonts);
    };
    std::array<double, 3> ys{};
    for (std::size_t k = 0; k < 3; ++k) {
        const Text& title = ax.titles[k];
        double y = box.y1 + pad;
        double localTop = top;
        if (!title.text.empty() && !ax.yaxis.offsetText.text.empty()) {
            Text o = ax.yaxis.offsetText;
            o.coords = Coords::Pixels;
            o.xy = d.yOffset;
            const Box ob = textExtent(o, nullptr, fonts);
            const Box tb = extent(title, y);
            if (!(ob.x1 < tb.x0 || tb.x1 < ob.x0 || ob.y1 < tb.y0 || tb.y1 < ob.y0)) localTop = ob.y1;
        }
        if (localTop < 0) {
            ys[k] = y;  // the top of the axes is not on the figure: not moved
            continue;
        }
        if (extent(title, y).y0 < localTop) {
            y = localTop + pad;
            const Box e = extent(title, y);
            if (e.y0 < localTop) {
                y = (2 * localTop - e.y0) + pad;
            }
        }
        ys[k] = y;
    }
    const double ymax = std::max({ys[0], ys[1], ys[2]});
    for (std::size_t k = 0; k < 3; ++k) {
        d.titles[k] = {box.x0 + ax.titles[k].xy.x * box.width(), ymax};
    }
}

void Figure::applyAspect(Axes& ax) {
    if (!ax.boxAspect) {
        ax.active = ax.position;
        return;
    }
    const double figAspect = heightInches() / widthInches();
    const Box pb1 = shrunkToAspect(ax.position, *ax.boxAspect, ax.position, figAspect);
    ax.active = anchored(pb1, ax.anchor, ax.position);
}

std::optional<Box> Figure::computeTightBox(std::size_t a, const FontSet& fonts, bool forLayout) {
    Axes& ax = axes[a];
    if (!ax.visible) return std::nullopt;
    applyAspect(ax);
    computeTicks(a, fonts);
    updateAxisLabels(a);
    std::vector<Box> bb;
    if (ax.axison) {
        for (AxisDirection d : {AxisDirection::X, AxisDirection::Y}) {
            if (auto b = axisTightBox(a, d, fonts, forLayout)) bb.push_back(*b);
        }
    }
    updateTitles(a, fonts);
    const Box box = axesBox(a);
    bb.push_back(box);
    for (int k = 0; k < 3; ++k) {
        const Text& t = ax.titles[static_cast<std::size_t>(k)];
        if (!t.visible) continue;
        Text c = t;
        c.coords = Coords::Pixels;
        c.xy = m_derived[a].titles[static_cast<std::size_t>(k)];
        Box bt = textExtent(c, nullptr, fonts);
        if (forLayout && bt.width() > 0) {
            bt.x0 = (bt.x0 + bt.x1) / 2 - 0.5;
            bt.x1 = bt.x0 + 1.0;
        }
        bb.push_back(bt);
    }
    auto extra = [&](const Box& b) {
        if (0 < b.width() && b.width() < kInf && 0 < b.height() && b.height() < kInf) bb.push_back(b);
    };
    for (const Artist& art : ax.artists) {
        if (const Text* t = std::get_if<Text>(&art)) {
            if (t->visible && t->inLayout && !t->clipOn) extra(textExtent(*t, &ax, fonts));
        }
    }
    if (ax.axison) {
        for (SpineSide s : {kLeft, kBottom, kRight, kTop}) {
            if (ax.spines[s].visible) extra(spineExtent(a, s));
        }
        if (ax.colorbar && ax.colorbar->outline.visible) extra(box);
    }
    if (ax.legend && ax.legend->visible && ax.legend->inLayout) {
        placeLegend(a, fonts);
        extra(m_derived[a].legend.box);
    }
    if (ax.patchVisible) extra(box);
    std::vector<Box> keep;
    for (const Box& b : bb) {
        if (b.width() != 0 || b.height() != 0) keep.push_back(b);
    }
    return unionOf(keep);
}

// ---- autoscaling ----------------------------------------------------------------------------------

void Figure::autoscaleAxes(Axes& ax) {
    Axes::DataLimits lx = ax.dataX;
    // twinx shares x: the shared axes' data limits count (autoscale_view over the shared siblings)
    for (const Axes& other : axes) {
        if (&other == &ax) continue;
        const bool shared = (other.twinOf && &axes[*other.twinOf] == &ax) || (ax.twinOf && &axes[*ax.twinOf] == &other);
        if (shared) {
            lx.add(other.dataX.lo);
            lx.add(other.dataX.hi);
            lx.sticky.insert(lx.sticky.end(), other.dataX.sticky.begin(), other.dataX.sticky.end());
        }
    }
    const MaxNLocator view = MaxNLocator::autoLocator();  // view_limits ("data" mode) of any locator
    auto scale = [&](Axis& axis, const Axes::DataLimits& l) {
        if (!axis.autoscale) return;
        AutoscaleOptions o;
        o.margin = axis.margin;
        o.stickyEdges = l.sticky;
        o.inverted = axis.limits.lo > axis.limits.hi;
        axis.limits = autoscale(l.empty() ? kInf : l.lo, l.empty() ? -kInf : l.hi, view, o);
    };
    scale(ax.xaxis, lx);
    scale(ax.yaxis, ax.dataY);
}

// ---- legends --------------------------------------------------------------------------------------

Figure::LegendGeometry Figure::legendGeometry(std::size_t a) const { return m_derived.at(a).legend; }

Point Figure::legendAnchor(std::size_t a, LegendLoc loc, double W, double H, double fontSize) const {
    static const std::array<Point, 11> codes{Point{0, 0},   Point{1, 1},   Point{0, 1},   Point{0, 0},
                                             Point{1, 0},   Point{1, 0.5}, Point{0, 0.5}, Point{1, 0.5},
                                             Point{0.5, 0}, Point{0.5, 1}, Point{0.5, 0.5}};
    const int code = std::clamp(static_cast<int>(loc), 1, 10);
    const Box container = axesBox(a).padded(-0.5 * fontSize * m_dpi / 72.0);  // borderaxespad
    const Box b = anchored(Box::fromBounds(0, 0, W, H), codes[static_cast<std::size_t>(code)], container);
    return {b.x0, b.y0};
}

LegendLoc Figure::bestLegendLocation(std::size_t a, double W, double H, const FontSet& fonts) const {
    const Axes& ax = axes[a];
    const double f = ax.legend ? ax.legend->fontSize : 10.0;
    // _auto_legend_data: line paths, collection offsets, rectangle and text boxes; the badness of a
    // place counts line vertices and offsets inside it, boxes it overlaps and lines crossing it
    std::vector<std::vector<Point>> paths;
    std::vector<Point> offsets;
    std::vector<Box> boxes;
    for (const Artist& art : ax.artists) {
        if (const Line* l = std::get_if<Line>(&art)) {
            std::vector<Point> p;
            for (std::size_t i = 0; i < std::min(l->x.size(), l->y.size()); ++i) {
                p.push_back(toPixels(ax, l->coords, {l->x[i], l->y[i]}));
            }
            paths.push_back(std::move(p));
        } else if (const Rectangle* r = std::get_if<Rectangle>(&art)) {
            const Point p0 = toPixels(ax, r->coords, {r->x, r->y});
            const Point p1 = toPixels(ax, r->coords, {r->x + r->w, r->y + r->h});
            boxes.push_back({std::min(p0.x, p1.x), std::min(p0.y, p1.y), std::max(p0.x, p1.x), std::max(p0.y, p1.y)});
        } else if (const Scatter* s = std::get_if<Scatter>(&art)) {
            for (const Point& o : s->offsets) offsets.push_back(toPixels(ax, Coords::Data, o));
        } else if (const Text* t = std::get_if<Text>(&art)) {
            boxes.push_back(textExtent(*t, &ax, fonts));
        }
    }
    int bestIdx = 1;
    long bestBad = std::numeric_limits<long>::max();
    for (int idx = 1; idx <= 10; ++idx) {
        const Point p = legendAnchor(a, static_cast<LegendLoc>(idx), W, H, f);
        const Box lb = Box::fromBounds(p.x, p.y, W, H);
        long bad = 0;
        auto inside = [&lb](const Point& o) { return lb.x0 < o.x && o.x < lb.x1 && lb.y0 < o.y && o.y < lb.y1; };
        for (const auto& path : paths) {  // the lines' vertices
            for (const Point& o : path) bad += inside(o) ? 1 : 0;
        }
        for (const Point& o : offsets) bad += inside(o) ? 1 : 0;
        for (const Box& b : boxes) {
            if (boxesOverlap(lb, b)) ++bad;
        }
        for (const auto& path : paths) {
            if (pathIntersectsRectangle(path, lb)) ++bad;
        }
        if (bad < bestBad) {
            bestBad = bad;
            bestIdx = idx;
        }
        if (bad == 0) break;
    }
    return static_cast<LegendLoc>(bestIdx);
}

void Figure::placeLegend(std::size_t a, const FontSet& fonts) {
    const Axes& ax = axes[a];
    const Legend& L = *ax.legend;
    const double dpicor = m_dpi / 72.0;
    const double f = L.fontSize;
    const double wH = 2.0 * f * dpicor, hH = 0.7 * f * dpicor;  // handlelength, handleheight
    const double sepItem = 0.8 * f * dpicor;                    // handletextpad
    const double sepRow = 0.5 * f * dpicor;                      // labelspacing
    const double pad = 0.4 * f * dpicor;                         // borderpad
    struct ItemBox {
        Box box;           // the item's box (handle and text)
        double textX = 0;  // the text box' x offset in the item
        double baseShift = 0;
    };
    std::vector<ItemBox> items;
    for (const LegendEntry& e : L.entries) {
        Text t;
        t.text = e.text;
        t.font.size = L.textSize;
        t.ha = HAlign::Left;
        t.va = VAlign::Baseline;
        const TextGeometry g = textGeometry(t, {0, 0}, fonts);
        const TextLayout lp = fonts.layout(g.face, "lp", fontSizeInPoints(L.textSize), m_dpi);
        const double h_ = lp.height(), d_ = lp.descent();
        double w = g.box.width(), h = g.box.height(), lastDescent = g.descent;
        double shift = 0.0;
        if (g.lines.size() > 1) {
            const double ydNew = 0.5 * h - 0.5 * (h_ - d_);
            shift = lastDescent - ydNew;
            lastDescent = ydNew;
        } else {
            const double hd = std::max(h_ - d_, h - lastDescent);
            h = hd + lastDescent;
        }
        const Box text = Box::fromBounds(0.0, -lastDescent, w, h);
        const Box handle = Box::fromBounds(0.0, 0.0, wH, hH);
        ItemBox it;
        const double y0 = std::min(handle.y0, text.y0), y1 = std::max(handle.y1, text.y1);
        it.box = Box::fromBounds(0.0, y0, wH + sepItem + w, y1 - y0);
        it.textX = wH + sepItem;
        it.baseShift = shift;
        items.push_back(it);
    }
    // the column (VPacker, align baseline) and the legend box (VPacker, pad borderpad, align center)
    double x0 = 0, x1 = 0, height = 0;
    std::vector<double> yoff;
    if (!items.empty()) {
        x0 = kInf;
        x1 = -kInf;
        double acc = 0;
        for (std::size_t k = 0; k < items.size(); ++k) {
            x0 = std::min(x0, items[k].box.x0);
            x1 = std::max(x1, items[k].box.x1);
            yoff.push_back(acc);
            acc += items[k].box.height() + sepRow;
        }
        height = acc - sepRow;
        for (std::size_t k = 0; k < items.size(); ++k) yoff[k] = height - (yoff[k] + items[k].box.y1);
        const double yd = yoff[0];
        for (double& v : yoff) v -= yd;
        const Box col = Box::fromBounds(x0, -yd, x1 - x0, height);
        x0 = col.x0;
        x1 = col.x1;
        height = col.height();
        // legend box: the single visible child (the handle box), centred
        const Box inner = Box::fromBounds(0.0, col.y0, x1 - x0, height).padded(pad);
        const double W = inner.width(), H = inner.height();
        // placement
        Point ll;
        if (L.loc == LegendLoc::Best) {
            ll = legendAnchor(a, bestLegendLocation(a, W, H, fonts), W, H, f);
        } else {
            ll = legendAnchor(a, L.loc, W, H, f);
        }
        LegendGeometry geo;
        geo.box = Box::fromBounds(ll.x, ll.y, W, H);
        // offsets: the legend box' origin, then the column, then each item
        const Point origin{ll.x - inner.x0, ll.y - inner.y0};
        for (std::size_t k = 0; k < items.size(); ++k) {
            LegendGeometry::Item it;
            const Point io{origin.x, origin.y + yoff[k]};
            it.handleOrigin = io;
            it.handleWidth = wH;
            it.handleHeight = hH;
            it.handleDescent = 0.0;
            it.textOrigin = {io.x + items[k].textX, io.y + items[k].baseShift};
            geo.items.push_back(it);
        }
        m_derived[a].legend = geo;
    } else {
        m_derived[a].legend = LegendGeometry{};
    }
}

// ---- constrained layout ---------------------------------------------------------------------------

namespace {

struct Margins {
    double left = 0, right = 0, top = 0, bottom = 0, leftcb = 0, rightcb = 0, topcb = 0, bottomcb = 0;
};

// One LayoutGrid (Matplotlib's _layoutgrid.LayoutGrid) for a GridSpec inside the figure's grid,
// solved in closed form: columns and rows fill the parent, inner sizes in the given ratios.
struct LayoutGrid {
    int nr = 1, nc = 1;
    std::vector<double> wr, hr;
    std::vector<double> ml, mr, mlcb, mrcb;  // per column
    std::vector<double> mt, mb, mtcb, mbcb;  // per row
    std::vector<double> lefts, rights, tops, bottoms;
    // the solved margins (kiwisolver's variable values: zero before the first solve)
    std::vector<double> sml, smr, smlcb, smrcb, smt, smb, smtcb, smbcb;

    LayoutGrid(int rows, int cols, std::vector<double> heightRatios, std::vector<double> widthRatios)
        : nr(rows), nc(cols), wr(std::move(widthRatios)), hr(std::move(heightRatios)) {
        if (wr.empty()) wr.assign(static_cast<std::size_t>(nc), 1.0);
        if (hr.empty()) hr.assign(static_cast<std::size_t>(nr), 1.0);
        for (auto* v : {&ml, &mr, &mlcb, &mrcb, &lefts, &rights, &sml, &smr, &smlcb, &smrcb}) v->assign(static_cast<std::size_t>(nc), 0.0);
        for (auto* v : {&mt, &mb, &mtcb, &mbcb, &tops, &bottoms, &smt, &smb, &smtcb, &smbcb}) v->assign(static_cast<std::size_t>(nr), 0.0);
    }
    void reset() {
        for (auto* v : {&ml, &mr, &mlcb, &mrcb}) std::fill(v->begin(), v->end(), 0.0);
        for (auto* v : {&mt, &mb, &mtcb, &mbcb}) std::fill(v->begin(), v->end(), 0.0);
    }
    static void editMin(std::vector<double>& v, int cell, double size) {
        double& m = v[static_cast<std::size_t>(cell)];
        if (size > m) m = size;
    }
    void editOuterMarginMins(const Margins& m, const SubplotSpec& ss) {
        editMin(ml, ss.col0, m.left);
        editMin(mlcb, ss.col0, m.leftcb);
        editMin(mr, ss.col1 - 1, m.right);
        editMin(mrcb, ss.col1 - 1, m.rightcb);
        editMin(mt, ss.row0, m.top);
        editMin(mtcb, ss.row0, m.topcb);
        editMin(mb, ss.row1 - 1, m.bottom);
        editMin(mbcb, ss.row1 - 1, m.bottomcb);
    }
    void solve(const Box& parent) {
        double sumM = 0, sumR = 0;
        for (std::size_t i = 0; i < static_cast<std::size_t>(nc); ++i) {
            sumM += ml[i] + mlcb[i] + mr[i] + mrcb[i];
            sumR += wr[i];
        }
        const double k = (parent.width() - sumM) / sumR;
        double x = parent.x0;
        for (std::size_t i = 0; i < static_cast<std::size_t>(nc); ++i) {
            lefts[i] = x;
            x += ml[i] + mlcb[i] + wr[i] * k + mr[i] + mrcb[i];
            rights[i] = x;
        }
        sumM = sumR = 0;
        for (std::size_t i = 0; i < static_cast<std::size_t>(nr); ++i) {
            sumM += mt[i] + mtcb[i] + mb[i] + mbcb[i];
            sumR += hr[i];
        }
        const double kh = (parent.height() - sumM) / sumR;
        double y = parent.y1;
        for (std::size_t i = 0; i < static_cast<std::size_t>(nr); ++i) {
            tops[i] = y;
            y -= mt[i] + mtcb[i] + hr[i] * kh + mb[i] + mbcb[i];
            bottoms[i] = y;
        }
        sml = ml, smr = mr, smlcb = mlcb, smrcb = mrcb, smt = mt, smb = mb, smtcb = mtcb, smbcb = mbcb;
    }
    Box inner(int r0, int r1, int c0, int c1) const {
        const auto C0 = static_cast<std::size_t>(c0), C1 = static_cast<std::size_t>(c1 - 1);
        const auto R0 = static_cast<std::size_t>(r0), R1 = static_cast<std::size_t>(r1 - 1);
        return {lefts[C0] + sml[C0] + smlcb[C0], bottoms[R1] + smb[R1] + smbcb[R1], rights[C1] - smr[C1] - smrcb[C1],
                tops[R0] - smt[R0] - smtcb[R0]};
    }
    Box forColorbar(int r0, int r1, int c0, int c1) const {
        const auto C0 = static_cast<std::size_t>(c0), C1 = static_cast<std::size_t>(c1 - 1);
        const auto R0 = static_cast<std::size_t>(r0), R1 = static_cast<std::size_t>(r1 - 1);
        return {lefts[C0] + smlcb[C0], bottoms[R1] + smbcb[R1], rights[C1] - smrcb[C1], tops[R0] - smtcb[R0]};
    }
    static double maxOf(const std::vector<double>& a, const std::vector<double>& b, int from, int to) {
        double m = -kInf;
        for (int i = from; i < to; ++i) m = std::max(m, a[static_cast<std::size_t>(i)] + b[static_cast<std::size_t>(i)]);
        return m;
    }
};

}  // namespace

void Figure::constrainedLayout(const FontSet& fonts) {
    initialPositions();
    LayoutGrid lg(grid.rows, grid.cols, grid.heightRatios, grid.widthRatios);
    const double wPad = layout.wPad / widthInches(), hPad = layout.hPad / heightInches();
    const double W = m_width, H = m_height;
    std::vector<std::size_t> inLayout;
    for (std::size_t i = 0; i < axes.size(); ++i) {
        if (axes[i].spec && axes[i].visible) inLayout.push_back(i);
    }
    auto colorbarsOf = [&](std::size_t parent) {
        std::vector<std::size_t> out;
        for (std::size_t k = 0; k < axes.size(); ++k) {
            if (axes[k].colorbar && axes[k].colorbar->parent == parent) out.push_back(k);
        }
        return out;
    };
    for (int pass = 0; pass < 2; ++pass) {
        lg.reset();
        double figTop = 0.0;
        for (std::size_t i : inLayout) {
            const SubplotSpec& ss = *axes[i].spec;
            Margins m;
            // get_margin_from_padding
            const double ws = layout.wspace / 2, hs = layout.hspace / 2;
            m.leftcb = m.rightcb = wPad;
            m.bottomcb = m.topcb = hPad;
            if (ws / grid.cols > wPad) {
                if (ss.col0 > 0) m.leftcb = ws / grid.cols;
                if (ss.col1 < grid.cols) m.rightcb = ws / grid.cols;
            }
            if (hs / grid.rows > hPad) {
                if (ss.row1 < grid.rows) m.bottomcb = hs / grid.rows;
                if (ss.row0 > 0) m.topcb = hs / grid.rows;
            }
            const Box pos = axes[i].position;
            const std::optional<Box> tb = computeTightBox(i, fonts, true);
            const Box bbox = tb ? toFraction(*tb, W, H) : pos;
            m.left += pos.x0 - bbox.x0;
            m.right += bbox.x1 - pos.x1;
            m.bottom += pos.y0 - bbox.y0;
            m.top += bbox.y1 - pos.y1;
            for (std::size_t c : colorbarsOf(i)) {
                const double pad = axes[c].colorbar->pad * lg.inner(ss.row0, ss.row1, ss.col0, ss.col1).width();
                const std::optional<Box> ctb = computeTightBox(c, fonts, true);
                const Box cbbox = ctb ? toFraction(*ctb, W, H) : axes[c].position;
                m.rightcb += cbbox.width() + pad;
                if (cbbox.y0 < bbox.y0) m.bottom += bbox.y0 - cbbox.y0;
                if (cbbox.y1 > bbox.y1) m.top += cbbox.y1 - bbox.y1;
            }
            lg.editOuterMarginMins(m, ss);
        }
        if (suptitle && suptitle->visible && suptitle->inLayout) {
            suptitle->xy.y = 1.0 - hPad;
            const Box b = toFraction(textExtent(*suptitle, nullptr, fonts), W, H);
            figTop = std::max(figTop, b.height() + 2 * hPad);
        }
        // match_submerged_margins
        for (std::size_t i : inLayout) {
            const SubplotSpec& s1 = *axes[i].spec;
            if (s1.col1 - s1.col0 > 1) {
                double maxsubl = LayoutGrid::maxOf(lg.ml, lg.mlcb, s1.col0 + 1, s1.col1);
                double maxsubr = LayoutGrid::maxOf(lg.mr, lg.mrcb, s1.col0, s1.col1 - 1);
                for (std::size_t j : inLayout) {
                    const SubplotSpec& s2 = *axes[j].spec;
                    if (s2.col1 - s2.col0 > 1) {
                        maxsubl = std::max(maxsubl, LayoutGrid::maxOf(lg.ml, lg.mlcb, s2.col0 + 1, s2.col1));
                        maxsubr = std::max(maxsubr, LayoutGrid::maxOf(lg.mr, lg.mrcb, s2.col0, s2.col1 - 1));
                    }
                }
                for (int c = s1.col0 + 1; c < s1.col1; ++c) LayoutGrid::editMin(lg.ml, c, maxsubl - lg.mlcb[static_cast<std::size_t>(c)]);
                for (int c = s1.col0; c < s1.col1 - 1; ++c) LayoutGrid::editMin(lg.mr, c, maxsubr - lg.mrcb[static_cast<std::size_t>(c)]);
            }
            if (s1.row1 - s1.row0 > 1) {
                double maxsubt = LayoutGrid::maxOf(lg.mt, lg.mtcb, s1.row0 + 1, s1.row1);
                double maxsubb = LayoutGrid::maxOf(lg.mb, lg.mbcb, s1.row0, s1.row1 - 1);
                for (std::size_t j : inLayout) {
                    const SubplotSpec& s2 = *axes[j].spec;
                    if (s2.row1 - s2.row0 > 1) {
                        maxsubt = std::max(maxsubt, LayoutGrid::maxOf(lg.mt, lg.mtcb, s2.row0 + 1, s2.row1));
                        maxsubb = std::max(maxsubb, LayoutGrid::maxOf(lg.mb, lg.mbcb, s2.row0, s2.row1 - 1));
                    }
                }
                for (int r = s1.row0 + 1; r < s1.row1; ++r) LayoutGrid::editMin(lg.mt, r, maxsubt - lg.mtcb[static_cast<std::size_t>(r)]);
                for (int r = s1.row0; r < s1.row1 - 1; ++r) LayoutGrid::editMin(lg.mb, r, maxsubb - lg.mbcb[static_cast<std::size_t>(r)]);
            }
        }
        lg.solve(Box{0.0, 0.0, 1.0, 1.0 - figTop});
        bool collapsed = false;
        for (int r = 0; r < grid.rows && !collapsed; ++r) {
            for (int c = 0; c < grid.cols; ++c) {
                const Box b = lg.inner(r, r + 1, c, c + 1);
                if (b.width() <= 0 || b.height() <= 0) {
                    collapsed = true;
                    break;
                }
            }
        }
        if (collapsed || inLayout.empty()) {
            if (collapsed) {
                m_warnings.emplace_back(
                    "constrained_layout not applied because axes sizes collapsed to zero.  Try making figure larger or "
                    "Axes decorations smaller.");
            }
            continue;
        }
        for (std::size_t i : inLayout) {
            const SubplotSpec& ss = *axes[i].spec;
            const Box b = lg.inner(ss.row0, ss.row1, ss.col0, ss.col1);
            axes[i].position = b;
            axes[i].active = b;
            for (std::size_t c : colorbarsOf(i)) {
                // reposition_colorbar (location right)
                Axes& cb = axes[c];
                const ColorbarInfo& info = *cb.colorbar;
                const Box bboxparent = lg.forColorbar(ss.row0, ss.row1, ss.col0, ss.col1);
                const Box pb = b;
                Box pbcb{pb.x0, pb.y0, pb.x0 + info.fraction * pb.width(), pb.y0 + info.shrink * pb.height()};
                pbcb = anchored(pbcb, info.anchor, pb);
                const Box cbpos = cb.position;
                const std::optional<Box> ctb = computeTightBox(c, fonts, true);
                const Box cbbbox = ctb ? toFraction(*ctb, W, H) : cbpos;
                const double cbpad = info.pad * pb.width();
                const double lmargin = cbpos.x0 - cbbbox.x0;
                const double dx = bboxparent.x1 - pbcb.x0 + cbpad + lmargin;
                cb.position = pbcb.translated(dx, 0.0);
                cb.anchor = info.anchor;
                cb.boxAspect = info.aspect;
                applyAspect(cb);
            }
        }
    }
}

// ---- draw -----------------------------------------------------------------------------------------

void Figure::draw(const FontSet& fonts) {
    m_warnings.clear();
    m_derived.assign(axes.size(), Derived{});
    initialPositions();
    for (Axes& ax : axes) {
        autoscaleAxes(ax);
    }
    if (layout.constrained && !m_frozen) {
        constrainedLayout(fonts);
    }
    for (std::size_t i = 0; i < axes.size(); ++i) {
        m_derived[i].tight = computeTightBox(i, fonts, false);
        if (axes[i].legend && axes[i].legend->visible) {
            placeLegend(i, fonts);
        }
    }
}

}  // namespace mdbridge::render
