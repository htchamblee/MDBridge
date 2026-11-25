// The figure framework (src/render/core figure.*, display_list.*, figure_dump.*) against MDATK's
// figures (data/render/figures/<case>.json, tests/parity/make_figure_cases.py and mpl_walker.py):
// every case's first frame is rebuilt from its record - the grid of its layout preset, each axes in
// its slot with its role (main, twin, colour bar), the axes' chrome (spines, tick styles, locators
// and formatters, labels, titles, limits), texts, legends and the data artists - then laid out by
// MDBridge's constrained layout and written by MDBridge's record writer, and compared with what
// Matplotlib did (owner decision D10, contents only; docs/specs/m4-plan.md section 3):
//
//   layout boxes (axes, twin and colour-bar positions, tight boxes) within 1.5 % of the figure size
//   ticks and tick labels equal whenever MDBridge's automatic bin count equals Matplotlib's
//   texts, axis labels, titles and legends at the same places within the same tolerance; legend
//   locations ("best") equal
//
// What the panels put into the axes (data, colours, limits from data) is wave 2's subject (the
// panels' own tests); this test checks the framework on the chrome every panel shares.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/figure.hpp"
#include "mdbridge/render/figure_dump.hpp"
#include "mdbridge/render/theme.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using core::Array;
using core::Value;
using namespace mdbridge::render;

namespace {

const std::filesystem::path kFonts = MDBRIDGE_FONT_DIR;
const std::filesystem::path kFigures = test::kData / "render" / "figures";

const FontSet& fonts() {
    static const FontSet set = FontSet::load([](std::string_view f) { return test::readFile(kFonts / std::string(f)); });
    return set;
}

const Value& at(const Value& v, std::string_view key) {
    const Value* p = v.find(key);
    if (!p) throw std::runtime_error("missing key " + std::string(key));
    return *p;
}
const Value* opt(const Value& v, std::string_view key) { return v.find(key); }
const Array& list(const Value& v) { return *v.get<Array>(); }
std::string str(const Value& v) { return v.isString() ? *v.get<std::string>() : std::string(); }
double num(const Value& v) {
    if (const auto* d = v.get<double>()) return *d;
    if (const auto* i = v.get<std::int64_t>()) return static_cast<double>(*i);
    if (const auto* s = v.get<std::string>()) {
        if (*s == "nan" || *s == "NaN") return std::nan("");
        if (*s == "inf") return INFINITY;
        if (*s == "-inf") return -INFINITY;
    }
    return std::nan("");
}
bool flag(const Value& v) { return v.get<bool>() ? *v.get<bool>() : false; }
Rgba color(const Value& v) {
    if (!v.isArray() || list(v).size() < 4 || list(v)[0].isArray()) return {0, 0, 0, 0};
    const Array& a = list(v);
    return {num(a[0]), num(a[1]), num(a[2]), num(a[3])};
}
Box boxOf(const Value& v) {
    const Array& a = list(v);
    return {num(a[0]), num(a[1]), num(a[2]), num(a[3])};
}
std::vector<double> nums(const Value& v) {
    std::vector<double> out;
    if (v.isArray()) {
        for (const Value& x : list(v)) out.push_back(num(x));
    }
    return out;
}

HAlign ha(const std::string& s) { return s == "left" ? HAlign::Left : s == "right" ? HAlign::Right : HAlign::Center; }
VAlign va(const std::string& s) {
    if (s == "top") return VAlign::Top;
    if (s == "center") return VAlign::Center;
    if (s == "bottom") return VAlign::Bottom;
    if (s == "center_baseline") return VAlign::CenterBaseline;
    return VAlign::Baseline;
}
Coords coords(const std::string& s) {
    if (s == "axes") return Coords::Axes;
    if (s == "figure") return Coords::Figure;
    if (s == "xaxis") return Coords::XAxis;
    if (s == "yaxis") return Coords::YAxis;
    return Coords::Data;
}
Marker marker(const Value& v) {
    const std::string m = str(v);
    return m == "o" ? Marker::Circle : m == "s" ? Marker::Square : Marker::None;
}

void base(ArtistBase& a, const Value& r) {
    a.zorder = num(at(r, "zorder"));
    a.visible = flag(at(r, "visible"));
    if (const Value* al = opt(r, "alpha"); al && !al->isNull()) a.alpha = num(*al);
    a.animated = flag(at(r, "animated"));
    if (const Value* l = opt(r, "label")) a.label = str(*l);
}

Text textFrom(const Value& r) {
    Text t;
    base(t, r);
    t.text = str(at(r, "text"));
    const std::vector<double> xy = nums(at(r, "xy"));
    t.xy = {xy[0], xy[1]};
    t.coords = coords(str(at(r, "coords")));
    t.ha = ha(str(at(r, "ha")));
    t.va = va(str(at(r, "va")));
    t.color = color(at(r, "color"));
    t.font.size = num(at(r, "size"));
    t.font.weight = str(at(r, "weight")) == "bold" ? 700 : 400;
    t.font.style = str(at(r, "style"));
    t.font.family.clear();
    for (const Value& f : list(at(r, "family"))) t.font.family.push_back(str(f));
    t.rotation = num(at(r, "rotation"));
    t.clipOn = flag(at(r, "clip_on"));
    t.wrap = flag(at(r, "wrap"));
    return t;
}

Line lineFrom(const Value& r) {
    Line l;
    base(l, r);
    if (const Value* x = opt(r, "x")) l.x = nums(*x);
    if (const Value* y = opt(r, "y")) l.y = nums(*y);
    if (const Value* t = opt(r, "transform")) l.coords = coords(str(*t));
    l.color = color(at(r, "color"));
    l.width = num(at(r, "lw"));
    const std::string ls = str(at(r, "ls"));
    l.style = ls == "-" ? LineStyle::Solid : ls == "None" || ls == "" || ls == " " ? LineStyle::None : LineStyle::Custom;
    if (const Value* d = opt(r, "dashes"); d && d->isArray()) {
        l.style = LineStyle::Custom;
        l.dashes.offset = num(list(*d)[0]);
        l.dashes.pattern = nums(list(*d)[1]);
    }
    l.marker = marker(at(r, "marker"));
    l.markerSize = num(at(r, "ms"));
    l.markerFace = color(at(r, "mfc"));
    l.markerEdge = color(at(r, "mec"));
    l.markerEdgeWidth = num(at(r, "mew"));
    return l;
}

void axisFrom(Axis& axis, const Value& r) {
    axis.visible = flag(at(r, "visible"));
    const Value& lab = at(r, "label");
    axis.label.text = str(at(lab, "text"));
    axis.label.color = color(at(lab, "color"));
    axis.label.font.size = num(at(lab, "size"));
    if (const Value* ts = opt(r, "tick_style"); ts && !ts->isNull()) {
        TickStyle& s = axis.ticks;
        s.tick1 = flag(at(*ts, "tick1"));
        s.tick2 = flag(at(*ts, "tick2"));
        s.label1 = flag(at(*ts, "label1"));
        s.label2 = flag(at(*ts, "label2"));
        s.color = color(at(*ts, "color"));
        s.length = num(at(*ts, "length"));
        s.width = num(at(*ts, "width"));
        s.labelColor = color(at(*ts, "labelcolor"));
        s.labelSize = num(at(*ts, "labelsize"));
        s.rotation = num(at(*ts, "rotation"));
        s.gridOn = flag(at(*ts, "gridOn"));
        s.gridColor = color(at(at(*ts, "grid"), "color"));
        s.gridWidth = num(at(at(*ts, "grid"), "lw"));
    }
    const Value& loc = at(r, "locator");
    const std::string cls = str(at(loc, "class"));
    if (cls == "FixedLocator") {
        axis.locator.kind = Locator::Kind::Fixed;
        axis.locator.locs = nums(at(loc, "locs"));
        if (const Value* nb = opt(loc, "nbins"); nb && !nb->isNull()) axis.locator.fixedBins = static_cast<int>(num(*nb));
    } else {
        axis.locator.kind = Locator::Kind::MaxN;
        MaxNLocatorParams p;
        const Value& nb = at(loc, "nbins");
        if (nb.isString()) p.nbins = std::nullopt;
        else p.nbins = static_cast<int>(num(nb));
        p.steps = nums(at(loc, "steps"));
        p.integer = flag(at(loc, "integer"));
        p.minNTicks = static_cast<int>(num(at(loc, "min_n_ticks")));
        p.symmetric = flag(at(loc, "symmetric"));
        axis.locator.params = p;
    }
    const Value& fmt = at(r, "formatter");
    if (str(at(fmt, "class")) == "ScalarFormatter") {
        axis.formatter.kind = Formatter::Kind::Scalar;
        axis.formatter.options.useOffset = flag(at(fmt, "useOffset"));
        axis.formatter.options.useMathText = flag(at(fmt, "useMathText"));
    } else {
        // FixedFormatter / FuncFormatter of set_ticks(locs, labels): the labels by position
        axis.formatter.kind = Formatter::Kind::Fixed;
        std::map<double, std::string> byLoc;
        for (const Value& t : list(at(r, "ticks"))) byLoc[num(list(t)[0])] = str(list(t)[1]);
        if (const Value* seq = opt(fmt, "seq")) {
            for (const Value& s : list(*seq)) axis.formatter.seq.push_back(str(s));
        } else {
            for (double l : axis.locator.locs) axis.formatter.seq.push_back(byLoc.count(l) ? byLoc[l] : std::string());
        }
    }
}

// The case's first frame rebuilt in MDBridge's figure model. With `autoscaleSeries`, the y limits
// of the time-series panels (timeseries, dual_axis without an ylim option) are left to MDBridge's
// autoscaling of the recorded lines instead of being taken from the record.
Figure rebuild(const Value& rec, bool autoscaleSeries = false) {
    const std::vector<double> size = nums(at(rec, "size"));
    // figure_size(w, h): even pixel sizes (RN-06: the canvas is exactly this size)
    const int w = std::max(2, static_cast<int>(size[0]) / 2 * 2), h = std::max(2, static_cast<int>(size[1]) / 2 * 2);
    const double dpi = num(at(rec, "dpi"));
    const Value& fr = at(list(at(rec, "frames"))[0], "figure");
    Figure fig(w, h, dpi);
    fig.facecolor = color(at(fr, "facecolor"));
    const std::string presetName = str(*at(rec, "project").path({"layout", "preset"}));
    const core::LayoutPreset* preset = core::findLayoutPreset(presetName);
    if (!preset) throw std::runtime_error("unknown preset " + presetName);
    fig.grid.rows = preset->rows;
    fig.grid.cols = preset->cols;
    fig.grid.heightRatios.assign(preset->heightRatios.begin(), preset->heightRatios.end());
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    std::map<int, std::size_t> mainOf;
    const Array& axes = list(at(fr, "axes"));
    for (const Value& a : axes) {
        const int panel = static_cast<int>(num(at(a, "panel")));
        const std::string role = str(at(a, "role"));
        std::size_t k = 0;
        if (role == "main") {
            const core::LayoutSlot& s = preset->slots[static_cast<std::size_t>(panel)];
            k = fig.addSubplot({s.row0, s.row1, s.col0, s.col1});
            mainOf[panel] = k;
        } else if (role == "twin") {
            k = fig.twinx(mainOf.at(panel));
        } else {
            const Value& cb = at(a, "colorbar");
            const int n = static_cast<int>(num(at(at(cb, "cmap"), "N")));
            std::vector<Rgba> lut(static_cast<std::size_t>(std::max(n, 1)), Rgba{0.5, 0.5, 0.5, 1});
            k = fig.colorbar(mainOf.at(panel), std::nullopt, Colormap::listed(str(at(at(cb, "cmap"), "name")), lut),
                             num(at(cb, "vmin")), num(at(cb, "vmax")), 0.04, 0.02);
            fig.axes[k].colorbar->outline.color = color(at(at(cb, "outline"), "color"));
            fig.axes[k].colorbar->outline.width = num(at(at(cb, "outline"), "lw"));
        }
        Axes& ax = fig.axes[k];
        ax.panel = panel;
        ax.role = role;
        ax.face = color(at(a, "facecolor"));
        ax.visible = flag(at(a, "visible"));
        ax.axison = flag(at(a, "axison"));
        ax.zorder = num(at(a, "zorder"));
        static const char* sides[] = {"left", "bottom", "right", "top"};
        for (int s = 0; s < 4; ++s) {
            const Value& sp = at(at(a, "spines"), sides[s]);
            ax.spines[static_cast<std::size_t>(s)] = {flag(at(sp, "visible")), color(at(sp, "color")), num(at(sp, "lw"))};
        }
        axisFrom(ax.xaxis, at(a, "xaxis"));
        axisFrom(ax.yaxis, at(a, "yaxis"));
        const std::vector<double> xl = nums(at(a, "xlim")), yl = nums(at(a, "ylim"));
        ax.xaxis.limits = {xl[0], xl[1]};
        ax.yaxis.limits = {yl[0], yl[1]};
        ax.xaxis.autoscale = ax.yaxis.autoscale = false;
        if (autoscaleSeries && role != "colorbar") {
            const Value& spec = list(*at(rec, "project").path({"layout", "panels"}))[static_cast<std::size_t>(panel)];
            const std::string type = str(at(spec, "type"));
            const bool ylim = opt(at(spec, "options"), "ylim") != nullptr;
            ax.yaxis.autoscale = (type == "timeseries" || type == "dual_axis") && !ylim;
        }
        const Value& titles = at(a, "titles");
        static const char* locs[] = {"left", "center", "right"};
        for (int t = 0; t < 3; ++t) {
            if (const Value* tr = opt(titles, locs[t])) {
                Text& title = ax.titles[static_cast<std::size_t>(t)];
                title.text = str(at(*tr, "text"));
                title.font.size = num(at(*tr, "size"));
                title.color = color(at(*tr, "color"));
            }
        }
        for (const Value& art : list(at(a, "artists"))) {
            const std::string kind = str(at(art, "kind"));
            if (kind == "text") {
                ax.add(textFrom(art));
            } else if (kind == "line") {
                ax.add(lineFrom(art));
            } else if (kind == "scatter") {
                Scatter s;
                base(s, art);
                for (const Value& p : list(at(art, "offsets"))) s.offsets.push_back({num(list(p)[0]), num(list(p)[1])});
                s.sizes = nums(at(art, "sizes"));
                ax.add(s);
            } else if (kind == "segments") {
                Segments g;
                base(g, art);
                for (const Value& seg : list(at(art, "segments"))) {
                    std::vector<Point> pts;
                    for (const Value& p : list(seg)) pts.push_back({num(list(p)[0]), num(list(p)[1])});
                    g.segments.push_back(std::move(pts));
                }
                ax.add(g);
            } else if (kind == "image") {
                Image im;
                base(im, art);
                const std::vector<double> shape = nums(at(art, "shape")), ext = nums(at(art, "extent"));
                im.rows = static_cast<int>(shape[0]);
                im.cols = static_cast<int>(shape[1]);
                im.extent = {ext[0], ext[1], ext[2], ext[3]};
                im.originUpper = str(at(art, "origin")) == "upper";
                im.rgba.assign(static_cast<std::size_t>(im.rows * im.cols * 4), 128);
                ax.add(im);
            } else if (kind == "legend") {
                Legend L;
                base(L, art);
                L.loc = static_cast<LegendLoc>(static_cast<int>(num(at(art, "loc"))));
                const Value& fr2 = at(art, "frame");
                L.frameOn = flag(at(fr2, "on"));
                L.face = color(at(fr2, "facecolor"));
                L.edge = color(at(fr2, "edgecolor"));
                L.edgeWidth = num(at(fr2, "lw"));
                for (const Value& e : list(at(art, "entries"))) {
                    LegendEntry le;
                    le.text = str(at(e, "text"));
                    L.textColor = color(at(e, "color"));
                    L.textSize = num(at(e, "size"));
                    if (const Value* hd = opt(e, "handle"); hd && str(at(*hd, "kind")) == "line") le.handle = lineFrom(*hd);
                    L.entries.push_back(std::move(le));
                }
                ax.legend = std::move(L);
            } else if (kind == "patch") {
                Rectangle r;
                base(r, art);
                if (const Value* xywh = opt(art, "xywh")) {
                    const std::vector<double> v = nums(*xywh);
                    r.x = v[0], r.y = v[1], r.w = v[2], r.h = v[3];
                }
                ax.add(r);
            }
            // pixel_image (the 3D layer, wave 2) and other (colour-bar meshes) carry no layout
        }
    }
    if (const Value* st = opt(fr, "suptitle"); st && !st->isNull()) {
        const Text t = textFrom(*st);
        fig.setSuptitle(t.text, t.color, t.font.size);
    }
    for (const Value& t : list(at(fr, "texts"))) fig.texts.push_back(textFrom(t));
    return fig;
}

std::vector<std::string> caseNames() {
    std::vector<std::string> out;
    if (!std::filesystem::is_directory(kFigures)) return out;
    for (const auto& e : std::filesystem::directory_iterator(kFigures)) {
        if (e.path().extension() == ".json") out.push_back(e.path().stem().string());
    }
    std::sort(out.begin(), out.end());
    return out;
}

constexpr double kTol = 0.015;  // of the figure size, per edge

struct Checker {
    double w, h;
    double worst = 0.0;
    std::string worstWhat;
    void box(const Box& mine, const Box& theirs, const std::string& what) {
        const double d = std::max({std::fabs(mine.x0 - theirs.x0) / w, std::fabs(mine.x1 - theirs.x1) / w,
                                   std::fabs(mine.y0 - theirs.y0) / h, std::fabs(mine.y1 - theirs.y1) / h});
        if (d > worst) {
            worst = d;
            worstWhat = what;
        }
        EXPECT_LE(d, kTol) << what << ": MDBridge [" << mine.x0 << ", " << mine.y0 << ", " << mine.x1 << ", " << mine.y1
                           << "] Matplotlib [" << theirs.x0 << ", " << theirs.y0 << ", " << theirs.x1 << ", "
                           << theirs.y1 << "]";
    }
};

class FigureCase : public ::testing::TestWithParam<std::string> {};

TEST_P(FigureCase, LayoutTicksAndTextsAsMatplotlib) {
    const Value rec = test::readJson(kFigures / (GetParam() + ".json"));
    if (!opt(rec, "frames") || list(at(rec, "frames")).empty()) {
        GTEST_SKIP() << "MDATK drew no frame (" << GetParam() << ": an error case of the panels)";
    }
    Figure fig = rebuild(rec);
    fig.draw(fonts());
    fig.freeze();
    const Value mine = dumpFigure(fig, fonts());
    const Value& theirs = at(list(at(rec, "frames"))[0], "figure");
    Checker c{static_cast<double>(fig.widthPx()), static_cast<double>(fig.heightPx())};
    const Array& ma = list(at(mine, "axes"));
    const Array& ta = list(at(theirs, "axes"));
    ASSERT_EQ(ma.size(), ta.size());
    int ticksCompared = 0, ticksSkipped = 0, legendMoved = 0;
    for (std::size_t k = 0; k < ta.size(); ++k) {
        const std::string id = "ax" + std::to_string(k) + " (" + str(at(ta[k], "role")) + ")";
        c.box(boxOf(at(ma[k], "position_px")), boxOf(at(ta[k], "position_px")), id + " position");
        if (!at(ta[k], "tight_px").isNull()) {
            c.box(boxOf(at(ma[k], "tight_px")), boxOf(at(ta[k], "tight_px")), id + " tight box");
        }
        for (const char* axisName : {"xaxis", "yaxis"}) {
            const Value& mx = at(ma[k], axisName);
            const Value& tx = at(ta[k], axisName);
            const Value& tl = at(tx, "locator");
            bool comparable = true;
            if (const Value* an = opt(tl, "auto_nbins")) {
                comparable = num(*an) == num(at(at(mx, "locator"), "auto_nbins"));
            }
            if (!comparable) {
                ++ticksSkipped;
                continue;
            }
            ++ticksCompared;
            const Array& mt = list(at(mx, "ticks"));
            const Array& tt = list(at(tx, "ticks"));
            ASSERT_EQ(mt.size(), tt.size()) << id << " " << axisName << " ticks";
            for (std::size_t i = 0; i < tt.size(); ++i) {
                const double a = num(list(mt[i])[0]), b = num(list(tt[i])[0]);
                EXPECT_NEAR(a, b, 1e-9 * std::max(1.0, std::fabs(b))) << id << " " << axisName << " tick " << i;
                EXPECT_EQ(str(list(mt[i])[1]), str(list(tt[i])[1])) << id << " " << axisName << " tick " << i;
            }
            EXPECT_EQ(str(at(mx, "offset_text")), str(at(tx, "offset_text"))) << id << " " << axisName;
            const Value& tle = at(at(tx, "label"), "extent_px");
            if (!tle.isNull()) {
                c.box(boxOf(at(at(mx, "label"), "extent_px")), boxOf(tle), id + " " + axisName + " label");
            }
        }
        // titles
        const Value& tt = at(ta[k], "titles");
        for (const char* loc : {"left", "center", "right"}) {
            const Value* t = opt(tt, loc);
            if (!t) continue;
            const Value* m = opt(at(ma[k], "titles"), loc);
            ASSERT_NE(m, nullptr) << id << " title " << loc;
            EXPECT_EQ(str(at(*m, "text")), str(at(*t, "text")));
            c.box(boxOf(at(*m, "extent_px")), boxOf(at(*t, "extent_px")), id + " title " + loc);
        }
        // texts with an extent, and legends
        const Array& mArts = list(at(ma[k], "artists"));
        std::map<std::string, const Value*> mine2;
        for (const Value& a : mArts) mine2[str(at(a, "kind")) == "legend" ? "legend" : str(at(a, "id"))] = &a;
        std::size_t n = 0;
        for (const Value& a : list(at(ta[k], "artists"))) {
            const std::string kind = str(at(a, "kind"));
            if (kind == "legend") {
                const Value* m = mine2["legend"];
                ASSERT_NE(m, nullptr);
                const Box tb = boxOf(at(a, "extent_px")), mb = boxOf(at(*m, "extent_px"));
                const auto loc = static_cast<LegendLoc>(static_cast<int>(num(at(a, "loc"))));
                if (loc != LegendLoc::Best) {
                    c.box(mb, tb, id + " legend");
                    continue;
                }
                // "best": the location Matplotlib chose for its box, found again by MDBridge's search
                // with Matplotlib's box size (the port of _find_best_position)...
                int mplCode = 0;
                {  // _get_anchored_bbox in Matplotlib's own axes box
                    const Box pa = boxOf(at(ta[k], "position_px")).padded(-0.5 * 10.0 * fig.dpi() / 72.0);
                    static const double cx[] = {0, 1, 0, 0, 1, 1, 0, 1, 0.5, 0.5, 0.5};
                    static const double cy[] = {0, 1, 1, 0, 0, 0.5, 0.5, 0.5, 0, 1, 0.5};
                    for (int code = 1; code <= 10 && !mplCode; ++code) {
                        const double x = pa.x0 + cx[code] * (pa.width() - tb.width());
                        const double y = pa.y0 + cy[code] * (pa.height() - tb.height());
                        if (std::fabs(x - tb.x0) < 0.01 && std::fabs(y - tb.y0) < 0.01) mplCode = code;
                    }
                }
                ASSERT_NE(mplCode, 0) << id << " legend: not at an anchored location";
                {  // in Matplotlib's geometry: its axes box and legend box
                    Figure probe = fig;
                    const Box pa = boxOf(at(ta[k], "position_px"));
                    probe.axes[k].position = probe.axes[k].active = {pa.x0 / probe.widthPx(), pa.y0 / probe.heightPx(),
                                                                     pa.x1 / probe.widthPx(), pa.y1 / probe.heightPx()};
                    EXPECT_EQ(static_cast<int>(probe.bestLegendLocation(k, tb.width(), tb.height(), fonts())), mplCode)
                        << id << " legend location in Matplotlib's geometry";
                }
                // ... and MDBridge's own box: the same place unless its (unhinted, wider) texts make
                // another location better; the box size within the tolerance either way
                const Box mbAtMpl = Box::fromBounds(tb.x0, tb.y0, mb.width(), mb.height());
                const int mine3 = static_cast<int>(fig.bestLegendLocation(k, mb.width(), mb.height(), fonts()));
                if (mine3 == mplCode) {
                    c.box(mb, Box::fromBounds(mb.x0, mb.y0, tb.width(), tb.height()), id + " legend size");
                    const Point p = fig.legendAnchor(k, static_cast<LegendLoc>(mplCode), mb.width(), mb.height(), 10.0);
                    c.box(mb, Box::fromBounds(p.x, p.y, mb.width(), mb.height()), id + " legend place");
                    c.box(mb, tb, id + " legend");
                } else {
                    ++legendMoved;
                    c.box(mbAtMpl, tb, id + " legend size (moved: text metrics)");
                }
                continue;
            }
            if (kind == "pixel_image" || kind == "other") continue;
            const std::string mid = "ax" + std::to_string(k) + "/" + std::to_string(n++);
            if (kind == "text" && !at(a, "extent_px").isNull()) {
                const Value* m = mine2[mid];
                ASSERT_NE(m, nullptr) << mid;
                EXPECT_EQ(str(at(*m, "text")), str(at(a, "text")));
                c.box(boxOf(at(*m, "extent_px")), boxOf(at(a, "extent_px")), id + " text '" + str(at(a, "text")) + "'");
            }
        }
    }
    // figure texts (watermark) and the suptitle
    if (const Value* st = opt(theirs, "suptitle"); st && !st->isNull()) {
        c.box(boxOf(at(at(mine, "suptitle"), "extent_px")), boxOf(at(*st, "extent_px")), "suptitle");
    }
    const Array& mtexts = list(at(mine, "texts"));
    const Array& ttexts = list(at(theirs, "texts"));
    ASSERT_EQ(mtexts.size(), ttexts.size());
    for (std::size_t i = 0; i < ttexts.size(); ++i) {
        if (!at(ttexts[i], "extent_px").isNull()) {
            c.box(boxOf(at(mtexts[i], "extent_px")), boxOf(at(ttexts[i], "extent_px")), "figure text " + std::to_string(i));
        }
    }
    // the display list: something drawn for every visible axes, the size of the figure
    const DisplayList dl = buildDisplayList(fig, fonts());
    EXPECT_EQ(dl.width, fig.widthPx());
    EXPECT_EQ(dl.height, fig.heightPx());
    EXPECT_FALSE(dl.items.empty());
    std::printf("[case] %-32s worst %.4f %% of the figure (%s); axis ticks compared %d, nbins differ %d, "
                "best legends elsewhere %d\n",
                GetParam().c_str(), 100.0 * c.worst, c.worstWhat.c_str(), ticksCompared, ticksSkipped, legendMoved);
}

// A minimal axes with data lines: the time-series panels' y limits from MDBridge's autoscaling of the
// recorded lines (margins 5 %, the target line, NaN values skipped) equal Matplotlib's.
TEST_P(FigureCase, AutoscaledLinesAsMatplotlib) {
    const Value rec = test::readJson(kFigures / (GetParam() + ".json"));
    if (!opt(rec, "frames") || list(at(rec, "frames")).empty()) GTEST_SKIP();
    if (GetParam() == "ts_progressive") {
        // the record holds the progressive line as drawn at position 0 (one point); Matplotlib's limits
        // come from the whole series, added before set_data shortened it (DataLimitsAsMatplotlib)
        GTEST_SKIP() << "progressive line: limits from the data at plot time";
    }
    Figure fig = rebuild(rec, true);
    fig.draw(fonts());
    const Array& ta = list(at(at(list(at(rec, "frames"))[0], "figure"), "axes"));
    int checked = 0;
    for (std::size_t k = 0; k < ta.size(); ++k) {
        if (!fig.axes[k].yaxis.autoscale) continue;
        bool hasLines = false;
        for (const Artist& a : fig.axes[k].artists) hasLines = hasLines || std::holds_alternative<Line>(a);
        if (!hasLines) continue;  // a "No data" panel
        const std::vector<double> yl = nums(at(ta[k], "ylim"));
        EXPECT_NEAR(fig.axes[k].yaxis.limits.lo, yl[0], 1e-12 * std::max(1.0, std::fabs(yl[0]))) << "ax" << k;
        EXPECT_NEAR(fig.axes[k].yaxis.limits.hi, yl[1], 1e-12 * std::max(1.0, std::fabs(yl[1]))) << "ax" << k;
        ++checked;
    }
    if (!checked) GTEST_SKIP() << "no time-series axes";
}

INSTANTIATE_TEST_SUITE_P(Figures, FigureCase, ::testing::ValuesIn(caseNames()),
                         [](const ::testing::TestParamInfo<std::string>& i) { return i.param; });

}  // namespace

namespace {

// A one-panel figure with MDATK's chrome: a styled axes with a line, a cursor and a dot
Figure smallFigure(int w = 640, int h = 360, double dpi = 60.0) {
    Figure fig(w, h, dpi);
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    fig.facecolor = toRgba(lightTheme().page);
    const std::size_t a = fig.addSubplot({0, 1, 0, 1});
    Axes& ax = fig.axes[a];
    styleAxes(ax, axesStyle(lightTheme(), true, 11.0));
    Line l;
    l.x = {0.0, 1.0, 2.0, 3.0};
    l.y = {1.0, 3.0, 2.0, 5.0};
    l.color = toRgba(lightTheme().seriesColor(0));
    l.width = 1.6;
    l.zorder = 3;
    ax.add(l);
    Line cursor;
    cursor.coords = Coords::XAxis;
    cursor.x = {1.0, 1.0};
    cursor.y = {0.0, 1.0};
    cursor.animated = true;
    ax.add(cursor);
    ax.setXLim(0.0, 3.0);
    ax.xaxis.label.text = "Time (ns)";
    ax.setTitle("RMSD", HAlign::Left).font.size = 11.0;
    return fig;
}

}  // namespace

TEST(FigureFramework, DataLimitsAsMatplotlib) {
    // the data limits are what was added (Matplotlib: set_data does not update dataLim); relim()
    // recomputes them from the artists as they are
    Figure fig = smallFigure();
    fig.draw(fonts());
    const Range before = fig.axes[0].yaxis.limits;
    EXPECT_DOUBLE_EQ(before.lo, 1.0 - 0.05 * 4.0);
    EXPECT_DOUBLE_EQ(before.hi, 5.0 + 0.05 * 4.0);
    std::get<Line>(fig.axes[0].artists[0]).y = {1.0, 1.5, 1.2, 2.0};
    fig.draw(fonts());
    EXPECT_EQ(fig.axes[0].yaxis.limits, before);
    fig.axes[0].relim();
    fig.draw(fonts());
    EXPECT_DOUBLE_EQ(fig.axes[0].yaxis.limits.hi, 2.0 + 0.05 * 1.0);
}

TEST(FigureFramework, FrozenLayoutAndExactPixels) {
    // the canvas is exactly the requested size at any dpi (MDATK: figsize w / dpi could round the
    // canvas, RN-06); after freeze() the positions stay when the contents change
    Figure fig = smallFigure(1000, 562, 77.0);
    fig.draw(fonts());
    const Box b = fig.axesBox(0);
    fig.freeze();
    fig.axes[0].xaxis.label.text = "A much longer label that would need more room";
    fig.axes[0].setYLim(-1000.0, 100000.0);
    fig.draw(fonts());
    EXPECT_EQ(fig.axesBox(0), b);
    const DisplayList dl = buildDisplayList(fig, fonts());
    EXPECT_EQ(dl.width, 1000);
    EXPECT_EQ(dl.height, 562);
}

TEST(FigureFramework, StaticAndAnimatedLayers) {
    // the frame loop's split: the background without the animated artists, then only those
    Figure fig = smallFigure();
    fig.draw(fonts());
    const DisplayList all = buildDisplayList(fig, fonts(), DisplayLayer::All);
    const DisplayList bg = buildDisplayList(fig, fonts(), DisplayLayer::Static);
    const DisplayList fg = buildDisplayList(fig, fonts(), DisplayLayer::Animated);
    ASSERT_EQ(fg.items.size(), 1u);  // the cursor
    EXPECT_TRUE(fg.items[0].animated);
    EXPECT_TRUE(std::holds_alternative<PathItem>(fg.items[0].item));
    EXPECT_EQ(bg.items.size() + fg.items.size(), all.items.size());
    for (const DisplayItem& i : bg.items) EXPECT_FALSE(i.animated);
    // the figure background first, then the axes background
    ASSERT_GE(bg.items.size(), 2u);
    EXPECT_TRUE(std::holds_alternative<FillRectItem>(bg.items[0].item));
    EXPECT_EQ(std::get<FillRectItem>(bg.items[0].item).color, toRgba(lightTheme().page));
}

TEST(FigureFramework, NoDataPanel) {
    // show_message: a bare frame without ticks, the message centred in muted ink
    Figure fig(960, 540, 60.0);
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    const std::size_t a = fig.addSubplot({0, 1, 0, 1});
    Axes& ax = fig.axes[a];
    const Theme& th = lightTheme();
    styleAxes(ax, axesStyle(th, true, 11.0));
    ax.clear();
    ax.face = toRgba(th.surface);
    ax.setXTicks({});
    ax.setYTicks({});
    for (Spine& s : ax.spines) s.color = toRgba(th.grid);
    Text t;
    t.text = "No data: no series selected";
    t.xy = {0.5, 0.5};
    t.coords = Coords::Axes;
    t.ha = HAlign::Center;
    t.va = VAlign::Center;
    t.wrap = true;
    t.color = toRgba(th.muted);
    t.font.size = 11.0 * 0.85;
    ax.add(t);
    fig.draw(fonts());
    // Matplotlib: the axes fills the figure but for the pads (0.08 in = 4.8 px at 60 dpi)
    const Box b = fig.axesBox(a);
    EXPECT_NEAR(b.x0, 4.8, 1e-9);
    EXPECT_NEAR(b.y0, 4.8, 1e-9);
    EXPECT_NEAR(b.x1, 955.2, 1e-9);
    EXPECT_NEAR(b.y1, 535.2, 1e-9);
    const Box e = fig.textExtent(std::get<Text>(ax.artists[0]), &ax, fonts());
    EXPECT_NEAR((e.x0 + e.x1) / 2, 480.0, 1e-9);
    EXPECT_NEAR((e.y0 + e.y1) / 2, 270.0, 1e-9);
}

TEST(FigureFramework, CollapsedLayoutWarns) {
    // decorations larger than the figure: Matplotlib's warning, positions not applied
    Figure fig(60, 40, 100.0);
    const std::size_t a = fig.addSubplot({0, 1, 0, 1});
    fig.axes[a].xaxis.label.text = "x";
    fig.axes[a].xaxis.label.font.size = 60.0;
    fig.draw(fonts());
    ASSERT_FALSE(fig.warnings().empty());
    EXPECT_NE(fig.warnings()[0].find("constrained_layout not applied"), std::string::npos);
}

TEST(FigureFramework, WrappedTexts) {
    // Text(wrap=True), as the "No data" messages: words join while the line fits twice the distance
    // from the (centred) text to the nearer figure edge; one more word would not fit
    Figure fig(200, 100, 60.0);
    Text t;
    t.text = "No data: analysis 'nokey' has not been computed for this project";
    t.ha = HAlign::Center;
    t.va = VAlign::Center;
    t.wrap = true;
    t.font.size = 9.35;
    const Point at{100.0, 50.0};
    const std::string wrapped = fig.wrappedText(t, at, fonts());
    std::vector<std::string> lines;
    std::string cur;
    for (char ch : wrapped + "\n") {
        if (ch == '\n') {
            lines.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    ASSERT_GT(lines.size(), 1u);
    const FontFace face = chooseFace(t.font.family, t.font.style, t.font.weight);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        EXPECT_LE(std::ceil(fonts().layout(face, lines[i], 9.35, 60.0).width()), 200.0) << lines[i];
        if (i + 1 < lines.size()) {
            const std::string next = lines[i] + " " + lines[i + 1].substr(0, lines[i + 1].find(' '));
            EXPECT_GT(std::ceil(fonts().layout(face, next, 9.35, 60.0).width()), 200.0) << next;
        }
    }
    // the text is centred on its position, its lines stacked
    const Figure::TextGeometry g = fig.textGeometry(t, at, fonts());
    EXPECT_EQ(g.lines.size(), lines.size());
    EXPECT_NEAR((g.box.y0 + g.box.y1) / 2, 50.0, 1e-9);
    t.wrap = false;
    EXPECT_EQ(fig.wrappedText(t, at, fonts()), t.text);
}
