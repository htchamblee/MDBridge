// The record writer of drawn figures (figure_dump.hpp), mirroring tests/parity/mpl_walker.py.

#include "mdbridge/render/figure_dump.hpp"

#include <string>

namespace mdbridge::render {

namespace {

using core::Array;
using core::Table;
using core::Value;

Value rgba(const Rgba& c) { return Array{c.r, c.g, c.b, c.a}; }
Value box(const Box& b) { return Array{b.x0, b.y0, b.x1, b.y1}; }
Value nums(const std::vector<double>& v) {
    Array a;
    for (double x : v) a.emplace_back(x);
    return a;
}

const char* haName(HAlign h) { return h == HAlign::Left ? "left" : h == HAlign::Center ? "center" : "right"; }
const char* vaName(VAlign v) {
    switch (v) {
        case VAlign::Top: return "top";
        case VAlign::Center: return "center";
        case VAlign::Baseline: return "baseline";
        case VAlign::CenterBaseline: return "center_baseline";
        case VAlign::Bottom: return "bottom";
    }
    return "baseline";
}
const char* coordsName(Coords c) {
    switch (c) {
        case Coords::Data: return "data";
        case Coords::Axes: return "axes";
        case Coords::Figure: return "figure";
        case Coords::XAxis: return "xaxis";
        case Coords::YAxis: return "yaxis";
        case Coords::Pixels: return "other";
    }
    return "other";
}
const char* lsName(LineStyle s) {
    switch (s) {
        case LineStyle::Solid: return "-";
        case LineStyle::Dashed: return "--";
        case LineStyle::DashDot: return "-.";
        case LineStyle::Dotted: return ":";
        case LineStyle::None: return "None";
        case LineStyle::Custom: return "--";
    }
    return "-";
}
const char* markerName(Marker m) { return m == Marker::Circle ? "o" : m == Marker::Square ? "s" : "None"; }
const char* capName(CapStyle c) { return c == CapStyle::Butt ? "butt" : c == CapStyle::Round ? "round" : "projecting"; }

Table common(const ArtistBase& a, const std::string& id, const char* kind) {
    Table t{{"id", id}, {"kind", kind}, {"zorder", a.zorder}, {"visible", a.visible},
            {"alpha", a.alpha ? Value(*a.alpha) : Value()}, {"animated", a.animated}};
    if (!a.label.empty() && a.label[0] != '_') t.emplace_back("label", a.label);
    return t;
}

Value textRecord(const Text& t, const std::string& id, const Box* extent, Coords coords) {
    Table r = common(t, id, "text");
    const char* weight = t.font.weight >= 700 ? "bold" : "normal";
    Array family;
    for (const std::string& f : t.font.family) family.emplace_back(f);
    r.insert(r.end(), {{"text", t.text}, {"xy", Array{t.xy.x, t.xy.y}}, {"coords", coordsName(coords)},
                       {"ha", haName(t.ha)}, {"va", vaName(t.va)}, {"color", rgba(t.color)}, {"size", t.font.size},
                       {"weight", weight}, {"style", t.font.style}, {"family", family}, {"rotation", t.rotation},
                       {"clip_on", t.clipOn}, {"wrap", t.wrap},
                       {"extent_px", extent && !t.text.empty() && t.visible ? box(*extent) : Value()}});
    return r;
}

Table lineFields(const Line& l, bool full) {
    Table r;
    if (full) {
        r.emplace_back("x", nums(l.x));
        r.emplace_back("y", nums(l.y));
    }
    Value dashes;
    if (l.style != LineStyle::Solid && l.style != LineStyle::None) {
        const Dashes u = l.unscaledDashes();
        dashes = Array{u.offset, nums(u.pattern)};
    }
    const Rgba mfc = l.markerFace.value_or(l.color), mec = l.markerEdge.value_or(l.color);
    r.insert(r.end(), {{"color", rgba(l.color)}, {"lw", l.width}, {"ls", lsName(l.style)}, {"dashes", dashes},
                       {"marker", markerName(l.marker)}, {"ms", l.markerSize}, {"mfc", rgba(mfc)}, {"mec", rgba(mec)},
                       {"mew", l.markerEdgeWidth}, {"drawstyle", "default"}, {"solid_capstyle", capName(l.capStyle)}});
    if (full) r.emplace_back("transform", coordsName(l.coords));
    return r;
}

Array colors(const std::vector<Rgba>& v) {
    Array a;
    for (const Rgba& c : v) a.push_back(rgba(c));
    return a;
}

Value locatorRecord(const Axis& axis) {
    const Locator& l = axis.locator;
    if (l.kind == Locator::Kind::Fixed) {
        return Table{{"class", "FixedLocator"}, {"locs", nums(l.locs)}, {"nbins", l.fixedBins ? Value(*l.fixedBins) : Value()}};
    }
    const MaxNLocatorParams& p = l.params;
    std::vector<double> steps = p.steps.empty() ? std::vector<double>{1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10} : p.steps;
    Table t{{"class", l.isAuto() && steps == std::vector<double>{1, 2, 2.5, 5, 10} && !p.integer ? "AutoLocator" : "MaxNLocator"},
            {"nbins", p.nbins ? Value(*p.nbins) : Value("auto")}, {"steps", nums(steps)}, {"integer", p.integer},
            {"min_n_ticks", p.minNTicks}};
    if (!p.nbins) {
        t.emplace_back("tick_space", axis.tickSpace);
        t.emplace_back("auto_nbins", MaxNLocator(p).bins(axis.tickSpace));
    }
    return t;
}

Value formatterRecord(const Axis& axis) {
    if (axis.formatter.kind == Formatter::Kind::Fixed) {
        Array seq;
        for (const std::string& s : axis.formatter.seq) seq.emplace_back(s);
        return Table{{"class", "FixedFormatter"}, {"seq", seq}};
    }
    Table t{{"class", "ScalarFormatter"}};
    if (axis.scalar) {
        t.emplace_back("offset", axis.scalar->offset());
        t.emplace_back("orderOfMagnitude", axis.scalar->orderOfMagnitude());
        t.emplace_back("decimals", axis.scalar->decimals());
    }
    return t;
}

}  // namespace

core::Value dumpFigure(const Figure& fig, const FontSet& fonts) {
    Table out{{"size_px", Array{fig.widthPx(), fig.heightPx()}}, {"dpi", fig.dpi()}, {"facecolor", rgba(fig.facecolor)}};
    if (fig.suptitle) {
        const Box e = fig.textExtent(*fig.suptitle, nullptr, fonts);
        out.emplace_back("suptitle", textRecord(*fig.suptitle, "fig/suptitle", &e, Coords::Figure));
    } else {
        out.emplace_back("suptitle", Value());
    }
    Array texts;
    for (std::size_t n = 0; n < fig.texts.size(); ++n) {
        const Box e = fig.textExtent(fig.texts[n], nullptr, fonts);
        texts.push_back(textRecord(fig.texts[n], "fig/text" + std::to_string(n), &e, fig.texts[n].coords));
    }
    out.emplace_back("texts", texts);
    Array axes;
    for (std::size_t k = 0; k < fig.axes.size(); ++k) {
        const Axes& ax = fig.axes[k];
        const std::string id = "ax" + std::to_string(k);
        Table r{{"id", id}, {"position_px", box(fig.axesBox(k))},
                {"tight_px", fig.tightBox(k) ? box(*fig.tightBox(k)) : Value()}, {"facecolor", rgba(ax.face)},
                {"visible", ax.visible}, {"axison", ax.axison},
                {"xlim", Array{ax.xaxis.limits.lo, ax.xaxis.limits.hi}},
                {"ylim", Array{ax.yaxis.limits.lo, ax.yaxis.limits.hi}}, {"aspect", "auto"}, {"zorder", ax.zorder}};
        Table titles;
        static const char* locs[] = {"left", "center", "right"};
        for (int t = 0; t < 3; ++t) {
            const Text& title = ax.titles[static_cast<std::size_t>(t)];
            if (title.text.empty()) continue;
            Text placed = title;
            placed.coords = Coords::Pixels;
            placed.xy = fig.titlePosition(k, t);
            const Box e = fig.textExtent(placed, nullptr, fonts);
            Value rec = textRecord(title, id + "/title_" + locs[t], &e, Coords::Pixels);
            titles.emplace_back(locs[t], rec);
        }
        r.emplace_back("titles", titles);
        Table spines;
        static const char* sides[] = {"left", "bottom", "right", "top"};
        for (int s = 0; s < 4; ++s) {
            const Spine& sp = ax.spines[static_cast<std::size_t>(s)];
            spines.emplace_back(sides[s], Table{{"visible", sp.visible}, {"color", rgba(sp.color)}, {"lw", sp.width}});
        }
        r.emplace_back("spines", spines);
        for (AxisDirection d : {AxisDirection::X, AxisDirection::Y}) {
            const Axis& axis = d == AxisDirection::X ? ax.xaxis : ax.yaxis;
            Text label = axis.label;
            label.coords = Coords::Pixels;
            label.xy = fig.axisLabelPosition(k, d);
            Array ticks;
            for (const DrawnTick& t : axis.drawn) {
                ticks.push_back(Array{t.loc, t.label, axis.ticks.label1, axis.ticks.label2});
            }
            const TickStyle& ts = axis.ticks;
            Table style{{"tick1", ts.tick1},   {"tick2", ts.tick2},         {"label1", ts.label1},
                        {"label2", ts.label2}, {"color", rgba(ts.color)},    {"length", ts.length},
                        {"width", ts.width},   {"labelcolor", rgba(ts.labelColor)}, {"labelsize", ts.labelSize},
                        {"rotation", ts.rotation}, {"gridOn", ts.gridOn},
                        {"grid", Table{{"color", rgba(ts.gridColor)}, {"lw", ts.gridWidth}, {"ls", "-"}}}};
            Table a{{"label", Table{{"text", axis.label.text}, {"color", rgba(axis.label.color)}, {"size", axis.label.font.size},
                                    {"extent_px", axis.label.text.empty() ? Value() : box(fig.textExtent(label, nullptr, fonts))}}},
                    {"scale", "linear"}, {"inverted", axis.limits.lo > axis.limits.hi}, {"visible", axis.visible},
                    {"ticks", ticks}, {"tick_style", style}, {"offset_text", axis.offsetText.text},
                    {"locator", locatorRecord(axis)}, {"formatter", formatterRecord(axis)}};
            r.emplace_back(d == AxisDirection::X ? "xaxis" : "yaxis", a);
        }
        if (ax.panel >= 0) {
            r.emplace_back("panel", ax.panel);
            r.emplace_back("role", ax.role);
        }
        if (ax.colorbar) {
            const ColorbarInfo& c = *ax.colorbar;
            r.emplace_back("colorbar", Table{{"orientation", "vertical"}, {"vmin", c.vmin}, {"vmax", c.vmax},
                                             {"cmap", c.cmap ? Value(Table{{"name", c.cmap->name()}, {"N", c.cmap->size()}}) : Value()},
                                             {"outline", Table{{"color", rgba(c.outline.color)}, {"lw", c.outline.width}}},
                                             {"label", ax.yaxis.label.text}});
        }
        Array arts;
        for (std::size_t n = 0; n < ax.artists.size(); ++n) {
            const std::string aid = id + "/" + std::to_string(n);
            const Artist& art = ax.artists[n];
            if (const Line* l = std::get_if<Line>(&art)) {
                Table t = common(*l, aid, "line");
                for (auto& kv : lineFields(*l, true)) t.push_back(std::move(kv));
                arts.push_back(t);
            } else if (const Scatter* s = std::get_if<Scatter>(&art)) {
                Table t = common(*s, aid, "scatter");
                Array offs;
                for (const Point& p : s->offsets) offs.push_back(Array{p.x, p.y});
                t.insert(t.end(), {{"offsets", offs}, {"sizes", nums(s->sizes)}, {"facecolors", colors(s->faces)},
                                   {"edgecolors", colors(s->edges)}, {"linewidths", nums(s->lineWidths)},
                                   {"marker", markerName(s->marker)}});
                arts.push_back(t);
            } else if (const Segments* g = std::get_if<Segments>(&art)) {
                Table t = common(*g, aid, "segments");
                Array segs;
                for (const auto& seg : g->segments) {
                    Array pts;
                    for (const Point& p : seg) pts.push_back(Array{p.x, p.y});
                    segs.push_back(pts);
                }
                t.insert(t.end(), {{"segments", segs}, {"colors", colors(g->colors)}, {"linewidths", nums(g->widths)},
                                   {"capstyle", capName(g->capStyle)}});
                arts.push_back(t);
            } else if (const Image* im = std::get_if<Image>(&art); im && im->fillsAxes) {
                Table t = common(*im, aid, "pixel_image");
                if (im->rows <= 0 || im->cols <= 0) {
                    t.emplace_back("shape", Value());
                } else {
                    std::int64_t coverage = 0;
                    for (std::size_t i = 3; i < im->rgba.size(); i += 4) coverage += im->rgba[i] > 0 ? 1 : 0;
                    t.insert(t.end(), {{"shape", Array{im->rows, im->cols, 4}}, {"dtype", "uint8"}, {"coverage", coverage}});
                }
                arts.push_back(t);
            } else if (const Image* im = std::get_if<Image>(&art)) {
                Table t = common(*im, aid, "image");
                // the array's shape: rows x cols of values, or rows x cols x 4 of RGBA data
                t.insert(t.end(), {{"shape", im->cmap ? Array{im->rows, im->cols} : Array{im->rows, im->cols, 4}},
                                   {"extent", Array{im->extent[0], im->extent[1], im->extent[2], im->extent[3]}},
                                   {"origin", im->originUpper ? "upper" : "lower"}, {"interpolation", im->interpolation}});
                if (!im->xEdges.empty()) t.emplace_back("x_edges", nums(im->xEdges));
                if (im->cmap) {
                    t.insert(t.end(), {{"cmap", im->cmap->name()}, {"vmin", im->vmin}, {"vmax", im->vmax}});
                }
                arts.push_back(t);
            } else if (const Text* tx = std::get_if<Text>(&art)) {
                const Box e = fig.textExtent(*tx, &ax, fonts);
                arts.push_back(textRecord(*tx, aid, &e, tx->coords));
            } else if (const Rectangle* rc = std::get_if<Rectangle>(&art)) {
                Table t = common(*rc, aid, "patch");
                t.insert(t.end(), {{"class", "Rectangle"}, {"facecolor", rc->face ? rgba(*rc->face) : Value("none")},
                                   {"edgecolor", rc->edge ? rgba(*rc->edge) : Value("none")}, {"lw", rc->lineWidth},
                                   {"xywh", Array{rc->x, rc->y, rc->w, rc->h}}});
                arts.push_back(t);
            }
        }
        if (ax.legend) {
            const Legend& L = *ax.legend;
            Table t = common(L, id + "/legend", "legend");
            Array entries;
            for (const LegendEntry& e : L.entries) {
                Table h = common(e.handle, "?Line2D", "line");
                for (auto& kv : lineFields(e.handle, false)) h.push_back(std::move(kv));
                entries.push_back(Table{{"text", e.text}, {"color", rgba(L.textColor)}, {"size", L.textSize}, {"handle", h}});
            }
            t.insert(t.end(), {{"loc", static_cast<int>(L.loc)}, {"extent_px", box(fig.legendGeometry(k).box)},
                               {"frame", Table{{"on", L.frameOn}, {"facecolor", rgba(L.face)}, {"edgecolor", rgba(L.edge)},
                                               {"lw", L.edgeWidth}}},
                               {"entries", entries}});
            arts.push_back(t);
        }
        r.emplace_back("artists", arts);
        axes.push_back(r);
    }
    out.emplace_back("axes", axes);
    return out;
}

}  // namespace mdbridge::render
