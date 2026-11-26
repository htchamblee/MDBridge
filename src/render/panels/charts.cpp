// The standalone charts (charts.hpp): mdatk/render/charts.py (MDATK's own code, ported) on the
// figure model, with MDBridge's improvements RN-01, RN-05, RN-10, RN-11, RN-14, RN-50 and the
// styles of owner decision D5 (docs/DIFFERENCES.md "Rendering: charts").

#include "mdbridge/render/panels/charts.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/color.hpp"
#include "mdbridge/render/panels/plots.hpp"
#include "mdbridge/render/theme.hpp"
#include "mdbridge/render/ticks.hpp"

namespace mdbridge::render::charts {

namespace {

constexpr std::string_view kFormats[] = {"png", "svg", "pdf", "jpg", "jpeg", "tif", "tiff", "eps", "ps"};

// Internal analyses feed panels, not charts.
bool internalKey(std::string_view key) { return key.starts_with("trace/") || key.starts_with("proximity/"); }

std::string fmtG(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%g", v);
    return buf;
}

// The pixels Matplotlib's canvas has for a figure of `inches` at `dpi` (int(): truncated), at least 1.
int canvasPixels(double inches, double dpi) {
    const double px = inches * dpi;
    if (!(px >= 1.0)) return 1;
    if (px >= static_cast<double>(std::numeric_limits<int>::max())) return std::numeric_limits<int>::max();
    return static_cast<int>(px);
}

// _style(name): true for MDATK's chrome. D5: "mdatk" or "" is MDATK's chrome, "default" is
// Matplotlib's default look; other names (Matplotlib's other styles, style files, typos) warn and
// use MDATK's chrome.
bool mdatkChrome(const std::string& style) {
    if (style.empty() || style == "mdatk") return true;
    if (style == "default") return false;
    core::logWarning("unknown matplotlib style " + core::pythonRepr(core::Value(style)) +
                     " (MDBridge draws the styles 'mdatk' and 'default'); using MDATK chrome");
    return true;
}

// MaxNLocator(integer=True)
void integerAxis(Axis& axis) {
    axis.locator.kind = Locator::Kind::MaxN;
    axis.locator.params = MaxNLocatorParams{10, {}, true};
}

// The chart's figure: width x height inches at dpi, constrained layout with Matplotlib's default
// pads (MDATK sets no others for charts), one axes, the chrome of the style.
struct ChartFrame {
    Figure fig;
    std::size_t ax = 0;
    bool mdatk = true;
    const Theme* theme = nullptr;
};

ChartFrame frame(const core::Project& project, int heightPx, bool grid) {
    const core::ChartParams& cp = project.output.charts;
    const double dpi = static_cast<double>(cp.dpi);
    ChartFrame f{Figure(canvasPixels(cp.width, dpi), heightPx, dpi)};
    f.mdatk = mdatkChrome(cp.style);
    f.theme = f.mdatk ? &theme(project.visual.theme) : &lightTheme();
    f.fig.layout = LayoutParams{};
    f.fig.facecolor = f.mdatk ? toRgba(f.theme->surface) : Rgba{1.0, 1.0, 1.0, 1.0};  // rc figure.facecolor
    f.fig.grid.rows = f.fig.grid.cols = 1;
    f.ax = f.fig.addSubplot({0, 1, 0, 1});
    if (f.mdatk) styleAxes(f.fig.axes[f.ax], axesStyle(*f.theme, grid, 11.0));
    return f;
}

// ax.set_title(text, loc="left"); RN-50: in the theme's ink at 11 pt with MDATK's chrome (as
// style_axes meant; MDATK styled the centre title and left this one black, 12 pt)
void leftTitle(ChartFrame& f, std::string text) {
    Text& t = f.fig.axes[f.ax].setTitle(std::move(text), HAlign::Left);
    if (f.mdatk) {
        t.color = toRgba(f.theme->ink);
        t.font.size = 11.0;
    } else {
        t.color = {0.0, 0.0, 0.0, 1.0};
        t.font.size = 12.0;  // axes.titlesize "large"
    }
}

Line plotLine(const std::vector<double>& x, std::vector<double> y, Rgba color, double width, std::string label) {
    Line l;
    l.x = x;
    l.y = std::move(y);
    l.color = color;
    l.width = width;
    l.label = std::move(label);
    return l;
}

// ax.legend(loc="best") over the axes' labelled lines, in the order they were added
Legend legendOf(const Axes& a) {
    Legend L;
    L.loc = LegendLoc::Best;
    for (const Artist& art : a.artists) {
        const Line* l = std::get_if<Line>(&art);
        if (!l || l->label.empty() || l->label[0] == '_') continue;
        LegendEntry e;
        e.text = l->label;
        e.handle = *l;
        e.handle.x.clear();
        e.handle.y.clear();
        L.entries.push_back(std::move(e));
    }
    return L;
}

}  // namespace

std::span<const std::string_view> chartFormats() { return kFormats; }

std::string normalizeFormat(std::string_view format) {
    std::string out;
    std::size_t start = 0;
    while (start < format.size() && format[start] == '.') ++start;
    for (char c : format.substr(start)) out.push_back(c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c);
    return out;
}

std::optional<std::string> chartFormatProblem(const std::vector<std::string>& formats) {
    std::vector<std::string> bad;
    for (const std::string& f : formats) {
        const std::string n = normalizeFormat(f);
        if (std::find(std::begin(kFormats), std::end(kFormats), n) == std::end(kFormats)) {
            bad.push_back(core::pythonRepr(core::Value(f)));
        }
    }
    if (bad.empty()) return std::nullopt;
    std::string names, choices;
    for (const std::string& b : bad) names += (names.empty() ? "" : ", ") + b;
    for (std::string_view f : kFormats) choices += (choices.empty() ? "" : ", ") + std::string(f);
    const bool many = bad.size() > 1;
    return std::string("chart format") + (many ? "s " : " ") + names + (many ? " are" : " is") +
           " not supported (choose from " + choices + ")";
}

std::string safeKey(std::string_view key) {
    auto allowed = [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
               c == '+' || c == '-';
    };
    std::string out;
    for (std::size_t i = 0; i < key.size();) {
        if (allowed(key[i])) {
            out.push_back(key[i++]);
            continue;
        }
        while (i < key.size() && !allowed(key[i])) ++i;  // a run (UTF-8 bytes of one character too)
        out.push_back('_');
    }
    const std::size_t b = out.find_first_not_of('_');
    if (b == std::string::npos) return {};
    return out.substr(b, out.find_last_not_of('_') - b + 1);
}

std::string formatTitle(std::string_view tmpl, std::string_view label, std::string_view key) {
    if (tmpl.empty()) return std::string(label);
    auto replaceAll = [](std::string s, std::string_view from, std::string_view to) {
        for (std::size_t at = s.find(from); at != std::string::npos; at = s.find(from, at + to.size())) {
            s.replace(at, from.size(), to);
        }
        return s;
    };
    return replaceAll(replaceAll(std::string(tmpl), "{label}", label), "{key}", key);
}

std::vector<ChartJob> chartJobs(const core::Dataset& dataset, const core::Project& project) {
    const core::ChartParams& cp = project.output.charts;
    std::vector<ChartJob> series, matrices;
    auto hasMatrix = [&](const std::string& k) {
        return std::any_of(matrices.begin(), matrices.end(), [&](const ChartJob& j) { return j.key == k; });
    };
    if (!cp.series.empty()) {
        // chart_keys: every missing key warns first
        for (const std::string& k : cp.series) {
            if (!dataset.has(k)) core::logWarning("chart series " + k + " is not available");
        }
        for (const std::string& k : cp.series) {
            if (dataset.has(k) && dataset.kindOf(k) == "series") series.push_back({"series", k});
        }
        // RN-10: matrices are charted (once each); pair maps warn (MDATK dropped both silently)
        for (const std::string& k : cp.series) {
            if (!dataset.has(k)) continue;
            const std::string kind = dataset.kindOf(k);
            if (kind == "matrices") {
                if (!hasMatrix(k)) matrices.push_back({"matrix", k});
            } else if (kind != "series") {
                core::logWarning("chart series " + k +
                                 " is a pair map; charts draw series and matrices (a contact_map panel shows it)");
            }
        }
    } else {
        for (const std::string& k : dataset.seriesKeys()) {
            if (!internalKey(k)) series.push_back({"series", k});
        }
    }
    if (cp.include_heatmaps) {
        for (const std::string& k : dataset.matrixKeys()) {
            if (!internalKey(k) && !hasMatrix(k)) matrices.push_back({"matrix", k});
        }
    }
    series.insert(series.end(), matrices.begin(), matrices.end());
    return series;
}

Figure seriesChart(const core::Dataset& dataset, const core::Project& project, const std::string& key,
                   const FontSet& fonts) {
    const core::ChartParams& cp = project.output.charts;
    const core::Series s = dataset.series(key);
    const std::vector<double> t = dataset.time(project.inputs.time_unit);
    ChartFrame f = frame(project, canvasPixels(cp.height, static_cast<double>(cp.dpi)), true);
    Axes& a = f.fig.axes[f.ax];
    // colour: output.charts.color, palette 0 (MDATK chrome) or the next colour of Matplotlib's
    // cycle for each line drawn without one
    std::size_t cycle = 0;
    auto lineColor = [&]() -> Rgba {
        if (!cp.color.empty()) return toRgba(cp.color);
        if (f.mdatk) return toRgba(f.theme->seriesColor(0));
        const std::span<const Rgba> c = defaultColorCycle();
        return c[cycle++ % c.size()];
    };
    const std::int64_t smooth = cp.smoothing;
    bool handles = smooth > 1;
    if (smooth > 1) {
        Line raw = plotLine(t, s.values, lineColor(), cp.line_width * 0.6, "Per frame");
        raw.alpha = 0.35;
        a.add(std::move(raw));
        a.add(plotLine(t, panels::runningMean(s.values, smooth), lineColor(), cp.line_width,
                       "Running mean (" + std::to_string(smooth) + " frames)"));
    } else {
        a.add(plotLine(t, s.values, lineColor(), cp.line_width, s.label));
    }
    if (cp.show_target && s.target) {
        // ax.axhline(target, color, lw=1.0, ls=(0, (4, 3)))
        Line l;
        l.coords = Coords::YAxis;
        l.x = {0.0, 1.0};
        l.y = {*s.target, *s.target};
        l.color = f.mdatk ? toRgba(f.theme->ink2) : toRgba("0.3");
        l.width = 1.0;
        l.style = LineStyle::Custom;
        l.dashes = {0.0, {4.0, 3.0}};
        l.label = "Target (" + fmtG(*s.target) + ")";
        a.add(std::move(l));
        handles = true;
    }
    if (handles) {
        a.legend = legendOf(a);
        if (f.mdatk) styleLegend(*a.legend, legendStyle(*f.theme, 11.0));
    }
    a.xaxis.label.text = !cp.xlabel.empty() ? cp.xlabel : dataset.timeLabelFor(project.inputs.time_unit);
    a.yaxis.label.text = !cp.ylabel.empty() ? cp.ylabel : s.axisLabel();
    if (s.unit == "count") integerAxis(a.yaxis);
    leftTitle(f, formatTitle(cp.title, s.label, key));
    if (t.size() > 1) f.fig.axes[f.ax].setXLim(t.front(), t.back());
    f.fig.draw(fonts);
    return std::move(f.fig);
}

int matrixChartHeightPx(double heightInches, std::size_t rows, double dpi) {
    const double inches = std::max(heightInches, 0.18 * static_cast<double>(rows) + 1.5);
    return std::min(canvasPixels(inches, dpi), kMaxChartPixels);
}

Figure matrixChart(const core::Dataset& dataset, const core::Project& project, const std::string& key,
                   const FontSet& fonts) {
    const core::ChartParams& cp = project.output.charts;
    const core::Matrix m = dataset.matrix(key);
    const std::vector<double> t = dataset.time(project.inputs.time_unit);
    if (t.empty()) throw std::invalid_argument("chart " + key + ": no analysed frames");
    const double dpi = static_cast<double>(cp.dpi);
    const std::size_t rows = m.rowLabels.size();
    const int height = matrixChartHeightPx(cp.height, rows, dpi);
    ChartFrame f = frame(project, height, false);
    // RN-05: 0.18 in per row made images beyond any raster's limits; capped at 2^15 pixels
    const int wanted = canvasPixels(std::max(cp.height, 0.18 * static_cast<double>(rows) + 1.5), dpi);
    if (wanted > kMaxChartPixels) {
        core::logWarning("chart " + key + ": " + std::to_string(rows) + " rows need " + std::to_string(wanted) +
                         " px at " + fmtG(dpi) + " dpi; height capped at " + std::to_string(kMaxChartPixels) + " px");
    }
    Figure& fig = f.fig;
    const Theme& th = *f.theme;

    Image im;
    const std::vector<std::size_t>& shape = m.values.shape();
    im.rows = shape.empty() ? 0 : static_cast<int>(shape[0]);
    im.cols = shape.size() < 2 ? 1 : static_cast<int>(shape[1]);
    im.values = m.values.toDouble();
    for (double& v : im.values) {
        if (!std::isfinite(v)) v = std::nan("");  // np.ma.masked_invalid
    }
    // RN-01: the frame cells span the midpoints between frame times (another column count than
    // frames keeps MDATK's uniform cells)
    im.xEdges = panels::columnEdges(t);
    if (im.xEdges.size() != static_cast<std::size_t>(im.cols) + 1) {
        const double dt = t.size() > 1 ? t[1] - t[0] : 1.0;
        im.xEdges.clear();
        im.extent[0] = t.front() - dt / 2;
        im.extent[1] = t.back() + dt / 2;
    } else {
        im.extent[0] = im.xEdges.front();
        im.extent[1] = im.xEdges.back();
    }
    im.extent[2] = static_cast<double>(rows) - 0.5;
    im.extent[3] = -0.5;

    std::size_t cb = 0;
    if (m.kind == "class" && m.classes && !m.classes->empty()) {
        const int n = static_cast<int>(m.classes->size());
        im.cmap = th.ordinal(n);
        im.vmin = -0.5;
        im.vmax = n - 0.5;
        const Colormap cmap = *im.cmap;
        const std::size_t i = fig.axes[f.ax].add(std::move(im));
        cb = fig.colorbar(f.ax, i, cmap, -0.5, n - 0.5);
        std::vector<double> locs;
        for (int k = 0; k < n; ++k) locs.push_back(k);
        fig.axes[cb].setYTicks(locs, *m.classes);
        fig.axes[cb].yaxis.locator.fixedBins = n;  // colorbar(ticks=range(n))
    } else {
        im.cmap = th.sequential().withExtremes(toRgba(th.grid));
        const panels::Norm norm = panels::heatmapNorm(im.values, std::nullopt, std::nullopt);
        im.vmin = norm.vmin;
        im.vmax = norm.vmax;
        const Colormap cmap = *im.cmap;
        const std::size_t i = fig.axes[f.ax].add(std::move(im));
        cb = fig.colorbar(f.ax, i, cmap, norm.vmin, norm.vmax);
        fig.axes[cb].yaxis.label.text = m.unit.empty() ? "count" : m.unit;
        if (m.kind == "count") integerAxis(fig.axes[cb].yaxis);
    }
    if (f.mdatk) {
        // RN-50: the colour bar in the theme's colours (MDATK: Matplotlib's black)
        Axes& c = fig.axes[cb];
        c.colorbar->outline.color = toRgba(th.grid);
        for (Axis* axis : {&c.xaxis, &c.yaxis}) {
            tickParams(*axis, {toRgba(th.muted), toRgba(th.ink2), std::nullopt, std::nullopt, std::nullopt, std::nullopt});
        }
        c.yaxis.label.color = toRgba(th.ink2);
    }
    Axes& a = fig.axes[f.ax];
    std::vector<double> locs;
    std::vector<std::string> names;
    for (long long i : thinTicks(static_cast<long long>(rows), 40)) {
        locs.push_back(static_cast<double>(i));
        names.push_back(m.rowLabels[static_cast<std::size_t>(i)]);
    }
    a.setYTicks(std::move(locs), std::move(names));
    a.yaxis.ticks.labelFamily = {"monospace"};  // residue labels in the sequence font (FONT-03)
    tickParams(a.yaxis, {std::nullopt, std::nullopt, 8.0, std::nullopt, std::nullopt, std::nullopt});
    a.xaxis.label.text = !cp.xlabel.empty() ? cp.xlabel : dataset.timeLabelFor(project.inputs.time_unit);
    a.yaxis.label.text = "Residue";
    leftTitle(f, formatTitle(cp.title, m.label, key));
    fig.draw(fonts);
    return std::move(f.fig);
}

std::vector<std::string> renderCharts(const core::Dataset& dataset, const core::Project& project,
                                      const std::filesystem::path& outDir, const core::Reporter& reporter,
                                      const std::string& stamp, const FontSet& fonts, const FigureWriter& writer) {
    const core::ChartParams& cp = project.output.charts;
    const std::vector<std::string> formats = cp.formats.empty() ? std::vector<std::string>{"png"} : cp.formats;
    // RN-11: before anything is drawn or written (MDATK: Matplotlib's error at the first save in
    // that format, the files before it written)
    if (const auto problem = chartFormatProblem(formats)) throw core::ValueError(*problem);
    std::filesystem::create_directories(outDir);
    const std::vector<ChartJob> jobs = chartJobs(dataset, project);
    std::vector<std::string> paths;
    for (std::size_t n = 0; n < jobs.size(); ++n) {
        reporter.checkCancel();
        const ChartJob& job = jobs[n];
        const Figure fig = job.kind == "series" ? seriesChart(dataset, project, job.key, fonts)
                                                : matrixChart(dataset, project, job.key, fonts);
        const std::string base = (outDir / (cp.prefix + "_" + safeKey(job.key) + (stamp.empty() ? "" : "_" + stamp))).string();
        for (const std::string& f : formats) {
            const std::string format = normalizeFormat(f);
            const std::string path = base + "." + format;
            writer(fig, path, format);
            paths.push_back(path);
        }
        reporter.progress(static_cast<double>(n + 1) / static_cast<double>(std::max<std::size_t>(1, jobs.size())),
                          "Chart " + job.key);
    }
    reporter.log("Wrote " + std::to_string(paths.size()) + " chart file(s) to " + outDir.string());
    return paths;
}

}  // namespace mdbridge::render::charts
