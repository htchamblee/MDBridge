// The plot panels (panel.hpp, plots.hpp): mdatk/render/panels.py TimeSeriesPanel, DualAxisPanel,
// HeatmapPanel, ContactMapPanel and _running_mean (MDATK's own code, ported), with MDBridge's
// improvements RN-01, RN-04 and RN-14 (docs/DIFFERENCES.md "Rendering: plot panels").

#include "mdbridge/render/panels/plots.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <set>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/render/panels/panel.hpp"

namespace mdbridge::render::panels {

// ---- arithmetic ---------------------------------------------------------------------------------

std::vector<double> runningMean(std::span<const double> y, long long window) {
    const long long n = static_cast<long long>(y.size());
    if (window <= 1 || n < 2) return {y.begin(), y.end()};
    window = std::min(window, n);
    const long long padL = window / 2;
    // np.pad(y, (pad_l, pad_r), mode="edge"), then the mean of each window of the padded series
    auto at = [&](long long i) { return y[static_cast<std::size_t>(std::clamp(i - padL, 0LL, n - 1))]; };
    std::vector<double> out(static_cast<std::size_t>(n));
    for (long long i = 0; i < n; ++i) {
        double sum = 0.0;
        long long count = 0;
        for (long long k = 0; k < window; ++k) {
            const double v = at(i + k);
            if (std::isfinite(v)) {
                sum += v;
                ++count;
            }
        }
        out[static_cast<std::size_t>(i)] = count ? sum / static_cast<double>(count) : std::nan("");
    }
    return out;
}

std::vector<double> columnEdges(std::span<const double> t) {
    const std::size_t n = t.size();
    if (n == 0) return {};
    if (n == 1) return {t[0] - 0.5, t[0] + 0.5};  // dt = 1 (MDATK)
    std::vector<double> e(n + 1);
    e[0] = t[0] - (t[1] - t[0]) / 2;
    for (std::size_t i = 1; i < n; ++i) e[i] = (t[i - 1] + t[i]) / 2;
    e[n] = t[n - 1] + (t[n - 1] - t[n - 2]) / 2;
    return e;
}

Norm heatmapNorm(std::span<const double> values, std::optional<double> vminOption, std::optional<double> vmaxOption) {
    Norm n;
    n.vmin = vminOption.value_or(0.0);
    double vmax = vmaxOption.value_or(0.0);
    if (vmax == 0.0) {
        double hi = -std::numeric_limits<double>::infinity();
        bool any = false;
        for (double v : values) {
            if (std::isfinite(v)) {
                hi = std::max(hi, v);
                any = true;
            }
        }
        vmax = any ? hi : 1.0;
    }
    n.vmax = std::max(vmax, 1e-9);
    if (n.vmin > n.vmax) {
        std::swap(n.vmin, n.vmax);
        n.swapped = true;
    }
    return n;
}

namespace {

std::string fmtG(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%g", v);
    return buf;
}

// ax.plot(x, y, color=..., lw=..., label=..., zorder=...)
Line plotLine(std::vector<double> x, std::vector<double> y, Rgba color, double width, std::string label, double zorder) {
    Line l;
    l.x = std::move(x);
    l.y = std::move(y);
    l.color = color;
    l.width = width;
    l.label = std::move(label);
    l.zorder = zorder;
    return l;
}

// ax.plot([x], [y], "o", ms=6, color=c, mec=surface, mew=1.5, zorder=5): the moving dot
Line dot(double x, double y, Rgba color, Rgba surface) {
    Line d;
    d.x = {x};
    d.y = {y};
    d.color = color;
    d.style = LineStyle::None;
    d.marker = Marker::Circle;
    d.markerSize = 6.0;
    d.markerEdge = surface;
    d.markerEdgeWidth = 1.5;
    d.zorder = 5.0;
    return d;
}

// ax.axvline(x, color, lw, alpha[, zorder])
Line axvline(double x, Rgba color, double width, double alpha, double zorder = 2.0) {
    Line l;
    l.coords = Coords::XAxis;
    l.x = {x, x};
    l.y = {0.0, 1.0};
    l.color = color;
    l.width = width;
    l.alpha = alpha;
    l.zorder = zorder;
    return l;
}

// integer_axis: MaxNLocator(integer=True)
void integerAxis(Axis& axis) {
    axis.locator.kind = Locator::Kind::MaxN;
    axis.locator.params = MaxNLocatorParams{10, {}, true};
}

// ax.legend(handles, loc="best") styled by style_legend: one entry per handle (a line with a label
// that does not start with "_"), the handle drawn like the line
Legend legendOf(const std::vector<const Line*>& handles, const Theme& theme, double fontSize) {
    Legend L;
    L.loc = LegendLoc::Best;
    for (const Line* l : handles) {
        if (l->label.empty() || l->label[0] == '_') continue;
        LegendEntry e;
        e.text = l->label;
        e.handle = *l;
        e.handle.x.clear();
        e.handle.y.clear();
        e.handle.zorder = 2.0;
        e.handle.alpha.reset();
        e.handle.animated = false;
        L.entries.push_back(std::move(e));
    }
    styleLegend(L, legendStyle(theme, fontSize));
    return L;
}

std::vector<const Line*> linesOf(const Axes& a) {
    std::vector<const Line*> out;
    for (const Artist& art : a.artists) {
        if (const Line* l = std::get_if<Line>(&art)) out.push_back(l);
    }
    return out;
}

Line& lineAt(Figure& fig, std::size_t axes, std::size_t index) { return std::get<Line>(fig.axes[axes].artists[index]); }

double xEnd(const std::vector<double>& t) { return t.size() > 1 ? t.back() : t.front() + 1.0; }

// ---- timeseries ---------------------------------------------------------------------------------

class TimeSeriesPanel final : public Panel {
public:
    TimeSeriesPanel(const core::PanelSpec& spec, Context& c) : Panel("timeseries", spec, c) {}

    std::vector<ArtistRef> animatedArtists() const override {
        std::vector<ArtistRef> out{{mainAxes(), m_cursor}};
        if (m_progressive) {
            for (std::size_t l : m_lines) out.push_back({mainAxes(), l});
        }
        for (std::size_t d : m_dots) out.push_back({mainAxes(), d});
        return out;
    }

protected:
    void setup(Figure& fig) override {
        const std::vector<double>& t = ctx.time;
        const std::vector<std::string>& keys = spec().series;
        if (keys.empty()) throw NoData("no series selected");
        std::vector<core::Series> series;
        for (const std::string& k : keys) series.push_back(ctx.series(k));
        m_progressive = optionBool("progressive", false);
        const long long smooth = optionInt("smooth", 0);
        Axes& a = ax(fig);
        const Rgba surface = toRgba(theme.surface);
        for (std::size_t k = 0; k < series.size(); ++k) {
            const Rgba color = toRgba(theme.seriesColor(static_cast<std::int64_t>(k)));
            std::vector<double> y = smooth > 1 ? runningMean(series[k].values, smooth) : series[k].values;
            m_lines.push_back(a.add(plotLine(t, y, color, 1.6, series[k].label, 3.0)));
            m_dots.push_back(a.add(dot(t.at(0), y.at(0), color, surface)));
            m_values.push_back(std::move(y));
        }
        std::optional<double> target;
        if (const core::Value* v = option("target"); v && !v->isNull()) {
            target = optionNumber("target", 0.0);
        } else {
            for (const core::Series& s : series) {
                if (s.target) {
                    target = *s.target;
                    break;
                }
            }
        }
        if (target) {
            Line l;
            l.coords = Coords::YAxis;
            l.x = {0.0, 1.0};
            l.y = {*target, *target};
            l.color = toRgba(theme.ink2);
            l.width = 1.0;
            l.style = LineStyle::Custom;
            l.dashes = {0.0, {4.0, 3.0}};
            l.label = "Target (" + fmtG(*target) + ")";
            l.zorder = 2.0;
            a.add(std::move(l));
        }
        m_cursor = a.add(axvline(t.at(0), toRgba(theme.ink), 1.0, 0.7, 4.0));
        a.setXLim(t.front(), xEnd(t));
        if (const core::Value* yl = option("ylim")) {
            const core::Array* arr = yl->get<core::Array>();
            auto number = [](const core::Value& v) -> double {
                if (const double* d = v.get<double>()) return *d;
                if (const std::int64_t* i = v.get<std::int64_t>()) return static_cast<double>(*i);
                throw core::ValueError("panel option 'ylim': expected two numbers");
            };
            if (!arr || arr->size() != 2) throw core::ValueError("panel option 'ylim': expected two numbers");
            a.setYLim(number((*arr)[0]), number((*arr)[1]));
        }
        a.xaxis.label.text = ctx.timeLabel;
        std::set<std::string> units;
        for (const core::Series& s : series) units.insert(s.unit);
        if (units == std::set<std::string>{"count"}) integerAxis(a.yaxis);
        m_autoTitle = series.size() == 1 ? series[0].label : std::string();
        if (series.size() == 1) {
            a.yaxis.label.text = series[0].axisLabel();
        } else if (units.size() == 1 && *units.begin() != "" && *units.begin() != "count") {
            a.yaxis.label.text = *units.begin();
        }
        if (series.size() > 1 || target) a.legend = legendOf(linesOf(a), theme, ctx.fontSize);
    }

    void drawFrame(Figure& fig, std::size_t pos) override {
        const std::vector<double>& t = ctx.time;
        lineAt(fig, mainAxes(), m_cursor).x = {t[pos], t[pos]};
        for (std::size_t k = 0; k < m_lines.size(); ++k) {
            const std::vector<double>& y = m_values[k];
            Line& d = lineAt(fig, mainAxes(), m_dots[k]);
            d.x = {t[pos]};
            d.y = {y[pos]};
            if (m_progressive) {
                Line& l = lineAt(fig, mainAxes(), m_lines[k]);
                l.x.assign(t.begin(), t.begin() + static_cast<std::ptrdiff_t>(pos + 1));
                l.y.assign(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(pos + 1));
            }
        }
    }

private:
    bool m_progressive = false;
    std::vector<std::size_t> m_lines, m_dots;
    std::vector<std::vector<double>> m_values;
    std::size_t m_cursor = 0;
};

// ---- dual_axis ----------------------------------------------------------------------------------

class DualAxisPanel final : public Panel {
public:
    DualAxisPanel(const core::PanelSpec& spec, Context& c) : Panel("dual_axis", spec, c) {}

    std::vector<ArtistRef> animatedArtists() const override {
        return {{mainAxes(), m_cursor}, {mainAxes(), m_d1}, {*twinAxes(), m_d2}};
    }

protected:
    void setup(Figure& fig) override {
        const std::vector<double>& t = ctx.time;
        if (spec().series.size() < 2) throw NoData("dual-axis panel needs two series");
        m_s1 = ctx.series(spec().series[0]);
        m_s2 = ctx.series(spec().series[1]);
        const Rgba c1 = toRgba(theme.seriesColor(0)), c2 = toRgba(theme.seriesColor(1));
        const Rgba surface = toRgba(theme.surface);
        const std::size_t twin = fig.twinx(mainAxes());
        m_twin = twin;
        style(fig.axes[twin], false);
        fig.axes[twin].spines[kRight].visible = true;
        fig.axes[twin].spines[kRight].color = toRgba(theme.axis);
        Axes* a = &ax(fig);
        const std::size_t l1 = a->add(plotLine(t, m_s1.values, c1, 1.6, m_s1.label + " (left axis)", 2.0));
        const std::size_t l2 = fig.axes[twin].add(plotLine(t, m_s2.values, c2, 1.6, m_s2.label + " (right axis)", 2.0));
        m_d1 = a->add(dot(t.at(0), m_s1.values.at(0), c1, surface));
        m_d2 = fig.axes[twin].add(dot(t.at(0), m_s2.values.at(0), c2, surface));
        a->yaxis.label.text = m_s1.axisLabel();
        fig.axes[twin].yaxis.label.text = m_s2.axisLabel();
        if (m_s1.unit == "count") integerAxis(a->yaxis);
        if (m_s2.unit == "count") integerAxis(fig.axes[twin].yaxis);
        m_autoTitle = m_s1.label + " / " + m_s2.label;
        a->xaxis.label.text = ctx.timeLabel;
        a->setXLim(t.front(), xEnd(t));
        fig.axes[twin].xaxis.limits = a->xaxis.limits;  // the shared x axis
        fig.axes[twin].xaxis.autoscale = false;
        m_cursor = a->add(axvline(t.at(0), toRgba(theme.ink), 1.0, 0.7));
        const Line h1 = lineAt(fig, mainAxes(), l1), h2 = lineAt(fig, twin, l2);
        a->legend = legendOf({&h1, &h2}, theme, ctx.fontSize);
    }

    void drawFrame(Figure& fig, std::size_t pos) override {
        const std::vector<double>& t = ctx.time;
        lineAt(fig, mainAxes(), m_cursor).x = {t[pos], t[pos]};
        Line& d1 = lineAt(fig, mainAxes(), m_d1);
        d1.x = {t[pos]};
        d1.y = {m_s1.values[pos]};
        Line& d2 = lineAt(fig, *twinAxes(), m_d2);
        d2.x = {t[pos]};
        d2.y = {m_s2.values[pos]};
    }

private:
    core::Series m_s1, m_s2;
    std::size_t m_d1 = 0, m_d2 = 0, m_cursor = 0;
};

// ---- maps ---------------------------------------------------------------------------------------

// _set_residue_ticks(axis.set_ticks, labels, max_ticks)
void residueTicks(Axes& a, AxisDirection d, const std::vector<std::string>& labels, long long maxTicks) {
    std::vector<double> locs;
    std::vector<std::string> names;
    for (long long i : thinTicks(static_cast<long long>(labels.size()), maxTicks)) {
        locs.push_back(static_cast<double>(i));
        names.push_back(labels[static_cast<std::size_t>(i)]);
    }
    if (d == AxisDirection::X) a.setXTicks(std::move(locs), std::move(names));
    else a.setYTicks(std::move(locs), std::move(names));
}

// cb.outline.set_edgecolor(grid); cb.ax.tick_params(colors=muted, labelcolor=ink2, labelsize=0.75 x font)
void styleColorbar(Axes& cb, const Theme& theme, double fontSize) {
    cb.colorbar->outline.color = toRgba(theme.grid);
    for (Axis* axis : {&cb.xaxis, &cb.yaxis}) {
        tickParams(*axis, {toRgba(theme.muted), toRgba(theme.ink2), fontSize * 0.75, std::nullopt, std::nullopt,
                           std::nullopt});
    }
}

void colorbarLabel(Axes& cb, const std::string& text, const Theme& theme, double fontSize) {
    cb.yaxis.label.text = text;
    cb.yaxis.label.color = toRgba(theme.ink2);
    cb.yaxis.label.font.size = fontSize * 0.8;
}

// ax.grid(False)
void noGrid(Axes& a) { a.xaxis.ticks.gridOn = a.yaxis.ticks.gridOn = false; }

class HeatmapPanel final : public Panel {
public:
    HeatmapPanel(const core::PanelSpec& spec, Context& c) : Panel("heatmap", spec, c) {}

    std::vector<ArtistRef> animatedArtists() const override { return {{mainAxes(), m_cursor}}; }

protected:
    void setup(Figure& fig) override {
        if (spec().series.empty()) throw NoData("no matrix selected");
        const core::Matrix m = ctx.matrix(spec().series[0]);
        if (m.values.size() == 0) throw NoData(m.key + " is empty");
        m_autoTitle = m.label;
        const std::vector<double>& t = ctx.time;
        const std::vector<std::size_t>& shape = m.values.shape();
        Image im;
        im.rows = shape.empty() ? 0 : static_cast<int>(shape[0]);
        im.cols = shape.size() < 2 ? 1 : static_cast<int>(shape[1]);
        im.values = m.values.toDouble();
        for (double& v : im.values) {
            if (!std::isfinite(v)) v = std::nan("");  // np.ma.masked_invalid
        }
        // RN-01: the cells span the midpoints between frame times (a matrix with another number of
        // columns than frames keeps MDATK's uniform cells)
        im.xEdges = columnEdges(t);
        if (im.xEdges.size() != static_cast<std::size_t>(im.cols) + 1) {
            const double dt = t.size() > 1 ? t[1] - t[0] : 1.0;
            im.xEdges.clear();
            im.extent[0] = t.front() - dt / 2;
            im.extent[1] = t.back() + dt / 2;
        } else {
            im.extent[0] = im.xEdges.front();
            im.extent[1] = im.xEdges.back();
        }
        im.extent[2] = static_cast<double>(m.rowLabels.size()) - 0.5;
        im.extent[3] = -0.5;
        Axes& a = ax(fig);
        std::optional<std::size_t> cb;
        if (m.kind == "class" && m.classes && !m.classes->empty()) {
            const int n = static_cast<int>(m.classes->size());
            im.cmap = theme.ordinal(n);
            im.vmin = -0.5;
            im.vmax = n - 0.5;
            const std::size_t i = a.add(im);
            cb = fig.colorbar(mainAxes(), i, *im.cmap, im.vmin, im.vmax, 0.04, 0.02);
            std::vector<double> locs;
            for (int k = 0; k < n; ++k) locs.push_back(k);
            Axes& c = fig.axes[*cb];
            c.setYTicks(locs, *m.classes);
            c.yaxis.locator.fixedBins = n;  // colorbar(ticks=range(n)): FixedLocator(ticks, nbins=len(ticks))
        } else {
            im.cmap = theme.sequential().withExtremes(toRgba(theme.grid));
            const Norm norm = heatmapNorm(im.values, option("vmin") ? std::optional(optionNumber("vmin", 0.0)) : std::nullopt,
                                          option("vmax") ? std::optional(optionNumber("vmax", 0.0)) : std::nullopt);
            if (norm.swapped) {
                // RN-04: MDATK's colour bar rewrote the norm (vmin > vmax) to an artefact of its
                // callbacks; MDBridge swaps the limits and says so
                core::logWarning("heatmap " + m.key + ": vmin " + fmtG(norm.vmax) + " is above vmax " +
                                 fmtG(norm.vmin) + "; using vmin " + fmtG(norm.vmin) + ", vmax " + fmtG(norm.vmax));
            }
            im.vmin = norm.vmin;
            im.vmax = norm.vmax;
            const std::size_t i = a.add(im);
            cb = fig.colorbar(mainAxes(), i, *im.cmap, im.vmin, im.vmax, 0.04, 0.02);
            colorbarLabel(fig.axes[*cb], m.unit.empty() ? "count" : m.unit, theme, ctx.fontSize);
            if (m.kind == "count") integerAxis(fig.axes[*cb].yaxis);
        }
        m_colorbar = cb;
        styleColorbar(fig.axes[*cb], theme, ctx.fontSize);
        Axes& a2 = ax(fig);  // fig.axes grew
        noGrid(a2);
        residueTicks(a2, AxisDirection::Y, m.rowLabels, optionInt("max_ticks", 16));
        tickParams(a2.yaxis, {std::nullopt, std::nullopt, ctx.fontSize * 0.7, std::nullopt, std::nullopt, std::nullopt});
        a2.xaxis.label.text = ctx.timeLabel;
        m_cursor = a2.add(axvline(t.at(0), toRgba(theme.ink), 1.2, 0.8));
    }

    void drawFrame(Figure& fig, std::size_t pos) override {
        lineAt(fig, mainAxes(), m_cursor).x = {ctx.time[pos], ctx.time[pos]};
    }

private:
    std::size_t m_cursor = 0;
};

class ContactMapPanel final : public Panel {
public:
    ContactMapPanel(const core::PanelSpec& spec, Context& c) : Panel("contact_map", spec, c) {}

    std::vector<ArtistRef> animatedArtists() const override {
        if (m_mode == "instant") return {{mainAxes(), m_image}};
        return {};
    }

protected:
    void setup(Figure& fig) override {
        if (spec().series.empty()) throw NoData("no map selected");
        m_pm = ctx.pairmap(spec().series[0]);
        const auto [nr, nc] = m_pm.shape();
        if (nr == 0 || nc == 0) throw NoData(m_pm.key + " is empty");
        m_mode = optionString("mode", "instant");
        m_transpose = optionBool("transpose", false);
        std::vector<double> data;
        double vmax = 1.0;
        std::string label = "count";
        if (m_mode == "frequency") {
            data = m_pm.frequency();
            label = "fraction of frames";
        } else if (m_mode == "frame") {
            const long long last = static_cast<long long>(m_pm.frameCount()) - 1;
            const long long f = std::min(std::max(static_cast<long long>(optionInt("frame", 0)), 0LL), last);
            data = m_pm.dense(static_cast<std::size_t>(std::max(f, 0LL)));
            vmax = std::max(1.0, *std::max_element(data.begin(), data.end()));
        } else {
            data = m_pm.dense(0);
            vmax = m_pm.values.empty() ? 1.0 : std::max(1.0, *std::max_element(m_pm.values.begin(), m_pm.values.end()));
        }
        Image im;
        im.rows = static_cast<int>(m_transpose ? nc : nr);
        im.cols = static_cast<int>(m_transpose ? nr : nc);
        im.values = arranged(data);
        im.cmap = theme.sequential();
        im.vmin = 0.0;
        im.vmax = vmax;
        im.extent = {-0.5, im.cols - 0.5, im.rows - 0.5, -0.5};
        im.originUpper = true;
        m_image = ax(fig).add(std::move(im));
        const Image& added = std::get<Image>(ax(fig).artists[m_image]);
        const std::size_t cb = fig.colorbar(mainAxes(), m_image, *added.cmap, 0.0, vmax, 0.04, 0.02);
        m_colorbar = cb;
        colorbarLabel(fig.axes[cb], label, theme, ctx.fontSize);
        if (label == "count") integerAxis(fig.axes[cb].yaxis);
        m_autoTitle = m_pm.label;
        styleColorbar(fig.axes[cb], theme, ctx.fontSize);
        Axes& a = ax(fig);
        noGrid(a);
        const long long mt = optionInt("max_ticks", 14);
        const std::vector<std::string>& rows = m_transpose ? m_pm.colLabels : m_pm.rowLabels;
        const std::vector<std::string>& cols = m_transpose ? m_pm.rowLabels : m_pm.colLabels;
        residueTicks(a, AxisDirection::Y, rows, mt);
        residueTicks(a, AxisDirection::X, cols, mt);
        for (Axis* axis : {&a.xaxis, &a.yaxis}) {
            tickParams(*axis, {std::nullopt, std::nullopt, ctx.fontSize * 0.65, std::nullopt, std::nullopt, std::nullopt});
        }
        tickParams(a.xaxis, {std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 90.0});
    }

    void drawFrame(Figure& fig, std::size_t pos) override {
        if (m_mode != "instant") return;
        std::get<Image>(ax(fig).artists[m_image]).values = arranged(m_pm.dense(pos));
    }

private:
    // the dense map (rows x cols, row-major) as shown: transposed when asked
    std::vector<double> arranged(const std::vector<double>& dense) const {
        if (!m_transpose) return dense;
        const auto [nr, nc] = m_pm.shape();
        std::vector<double> out(dense.size());
        for (std::size_t r = 0; r < nr; ++r) {
            for (std::size_t c = 0; c < nc; ++c) out[c * nr + r] = dense[r * nc + c];
        }
        return out;
    }

    core::PairMap m_pm;
    std::string m_mode;
    bool m_transpose = false;
    std::size_t m_image = 0;
};

const PanelRegistration kPlotPanels[] = {
    {"timeseries", &Panel::make<TimeSeriesPanel>},
    {"dual_axis", &Panel::make<DualAxisPanel>},
    {"heatmap", &Panel::make<HeatmapPanel>},
    {"contact_map", &Panel::make<ContactMapPanel>},
};

}  // namespace

std::span<const PanelRegistration> plotPanels() { return kPlotPanels; }

}  // namespace mdbridge::render::panels
