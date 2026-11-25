// The figure's panels (src/render/panels) against MDATK's figures (data/render/figures,
// tests/parity/make_figure_cases.py): every case whose panels are all implemented is built from its
// input dataset and project as build_figure builds it, then compared with MDATK's record by its
// contents (figure_compare.hpp): the first position's whole figure (layout, chrome, data artists,
// legends, colour bars), the images and colour-bar tables of the .npz, the panels' records
// (ok, automatic title, axes, animated artists), the animated artists at the middle and last
// positions, and the log. Plus the plot panels' own arithmetic (plots.hpp: RN-01, RN-04, RN-14) and
// the frame loop's static / animated split.

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "figure_compare.hpp"
#include "figure_records.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/figure_dump.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/panels/plots.hpp"

using namespace mdbridge;
using namespace mdbridge::render;
using core::Array;
using core::Value;
using test::FigureComparison;

namespace {

const Value& at(const Value& v, std::string_view key) { return FigureComparison::at(v, key); }
std::string str(const Value& v) { return FigureComparison::str(v); }
double num(const Value& v) { return FigureComparison::num(v); }

const core::Dataset& input(const std::string& name) {
    static std::map<std::string, core::Dataset> cache;
    auto it = cache.find(name);
    if (it == cache.end()) it = cache.emplace(name, test::loadRenderInput(name)).first;
    return it->second;
}

// The cases whose panel types all have an implementation (the rest wait for their wave-2 chunk)
std::vector<std::string> caseNames() {
    std::vector<std::string> out;
    if (!std::filesystem::is_directory(test::kFigureCases)) return out;
    const std::vector<std::string_view> types = panels::registeredPanelTypes();
    for (const auto& e : std::filesystem::directory_iterator(test::kFigureCases)) {
        if (e.path().extension() != ".json") continue;
        const Value rec = test::readJson(e.path());
        bool ours = true;
        const Value* ps = rec.path({"project", "layout", "panels"});
        for (const Value& p : *ps->get<Array>()) {
            ours = ours && std::find(types.begin(), types.end(), str(at(p, "type"))) != types.end();
        }
        if (ours) out.push_back(e.path().stem().string());
    }
    std::sort(out.begin(), out.end());
    return out;
}

// the log lines ("LEVEL message") while a figure is built
struct LogCapture {
    std::vector<std::string> lines;
    core::LogSink previous;
    LogCapture() {
        previous = core::setLogSink([this](core::LogLevel level, std::string_view message) {
            if (level >= core::LogLevel::Warning) lines.push_back(std::string(core::logLevelName(level)) + " " + std::string(message));
        });
    }
    ~LogCapture() { core::setLogSink(previous); }
};

// the record of artist `id` in a dump
const Value* artistIn(const Value& dump, const std::string& id) {
    for (const Value& ax : *at(dump, "axes").get<Array>()) {
        for (const Value& a : *at(ax, "artists").get<Array>()) {
            if (str(at(a, "id")) == id) return &a;
        }
    }
    return nullptr;
}

class PanelCase : public ::testing::TestWithParam<std::string> {};

TEST_P(PanelCase, ContentsAsMdatk) {
    const Value rec = test::readJson(test::kFigureCases / (GetParam() + ".json"));
    const std::string inputName = str(at(rec, "input"));
    const core::Dataset& ds = input(inputName);
    const core::Project project = test::caseProject(rec);
    LogCapture log;
    const Value* error = FigureComparison::opt(rec, "error");
    if (error && str(at(*error, "stage")) == "build") {
        EXPECT_ANY_THROW(test::buildCase(rec, ds, project)) << str(at(*error, "message"));
        return;
    }
    panels::BuiltFigure built = test::buildCase(rec, ds, project);
    const Array& frames = *at(rec, "frames").get<Array>();
    ASSERT_FALSE(frames.empty());
    const Value mine = dumpFigure(built.figure, test::figureFonts());
    FigureComparison c;
    c.width = built.figure.widthPx();
    c.height = built.figure.heightPx();
    c.figure(mine, at(frames[0], "figure"));

    // the panels' records
    const Array& tp = *at(rec, "panels").get<Array>();
    ASSERT_EQ(built.panels.size(), tp.size());
    for (std::size_t k = 0; k < tp.size(); ++k) {
        const panels::Panel& p = *built.panels[k];
        const std::string what = "panel " + std::to_string(k);
        EXPECT_EQ(p.type(), str(at(tp[k], "type"))) << what;
        EXPECT_EQ(p.ok(), num(at(tp[k], "ok")) != 0.0) << what;
        EXPECT_EQ(p.autoTitle(), str(at(tp[k], "auto_title"))) << what;
        const Value& axes = at(tp[k], "axes");
        EXPECT_EQ("ax" + std::to_string(p.mainAxes()), str(at(axes, "main"))) << what;
        if (const Value* t = FigureComparison::opt(axes, "twin")) {
            EXPECT_EQ(p.twinAxes() ? "ax" + std::to_string(*p.twinAxes()) : "-", str(*t)) << what;
        } else {
            EXPECT_FALSE(p.twinAxes()) << what;
        }
        if (const Value* cb = FigureComparison::opt(axes, "colorbar")) {
            EXPECT_EQ(p.colorbarAxes() ? "ax" + std::to_string(*p.colorbarAxes()) : "-", str(*cb)) << what;
        } else {
            EXPECT_FALSE(p.colorbarAxes()) << what;
        }
        std::vector<std::string> ids;
        if (p.ok()) {
            for (const panels::ArtistRef& r : p.animatedArtists()) ids.push_back(r.id());
        }
        std::vector<std::string> want;
        for (const Value& a : *at(tp[k], "animated").get<Array>()) want.push_back(str(a));
        EXPECT_EQ(ids, want) << what << " animated artists";
    }

    // images and colour-bar tables
    const std::filesystem::path npzPath = test::kFigureCases / (GetParam() + ".npz");
    const auto npz = test::readStoredNpz(npzPath);
    int arrays = test::compareArrays(built.figure, npz, 0);

    // the later positions: the animated artists' changes
    for (std::size_t n = 1; n < frames.size(); ++n) {
        const std::size_t pos = static_cast<std::size_t>(num(at(frames[n], "pos")));
        built.update(pos, test::figureFonts());
        const Value now = dumpFigure(built.figure, test::figureFonts());
        for (const Value& a : *at(frames[n], "artists").get<Array>()) {
            const std::string id = str(at(a, "id"));
            const Value* m = artistIn(now, id);
            ASSERT_NE(m, nullptr) << "pos " << pos << " " << id;
            for (const auto& [key, v] : *a.get<core::Table>()) {
                if (key == "id" || key == "extent_px") continue;
                const Value* mv = m->find(key);
                ASSERT_NE(mv, nullptr) << "pos " << pos << " " << id << "." << key;
                c.value(*mv, v, "pos " + std::to_string(pos) + " " + id + "." + key, key.find("color") != std::string::npos);
            }
        }
        arrays += test::compareArrays(built.figure, npz, pos);
    }

    // the warnings MDATK logged (build_figure, the panels)
    std::vector<std::string> want;
    for (const Value& l : *at(rec, "log").get<Array>()) want.push_back(str(l));
    EXPECT_EQ(log.lines, want);
    std::printf("[case] %-28s worst box %.3f %% (%s); ticks %d (nbins differ %d); %d numbers, %d arrays\n",
                GetParam().c_str(), 100.0 * c.worstBox, c.worstWhat.c_str(), c.ticksCompared, c.ticksSkipped, c.numbers,
                arrays);
}

INSTANTIATE_TEST_SUITE_P(Figures, PanelCase, ::testing::ValuesIn(caseNames()),
                         [](const ::testing::TestParamInfo<std::string>& i) { return i.param; });

}  // namespace

TEST(PanelRegistry, UnknownTypeAsMdatk) {
    core::Dataset ds({0}, {0.0});
    core::Project p;
    panels::Context ctx(ds, p);
    Figure fig(100, 100, 72.0);
    core::PanelSpec spec;
    spec.type = "nope";
    try {
        (void)panels::createPanel(fig, {0, 1, 0, 1}, spec, ctx);
        FAIL();
    } catch (const core::ValueError& e) {
        EXPECT_STREQ(e.what(), "unknown panel type 'nope' (choose from structure, timeseries, dual_axis, heatmap, "
                               "contact_map, sequence, secondary)");
    }
}

TEST(PlotArithmetic, RunningMeanNanAware) {
    // RN-14: MDATK's convolution spread one NaN over every window it fell in; MDBridge averages the
    // finite values of each window (edge padding as MDATK: pad_l = w // 2, pad_r = w - 1 - pad_l)
    const double nan = std::nan("");
    const std::vector<double> y{1.0, 2.0, nan, 4.0, 5.0};
    const std::vector<double> m = panels::runningMean(y, 3);
    ASSERT_EQ(m.size(), 5u);
    EXPECT_DOUBLE_EQ(m[0], (1.0 + 1.0 + 2.0) / 3);
    EXPECT_DOUBLE_EQ(m[1], (1.0 + 2.0) / 2);
    EXPECT_DOUBLE_EQ(m[2], (2.0 + 4.0) / 2);
    EXPECT_DOUBLE_EQ(m[3], (4.0 + 5.0) / 2);
    EXPECT_DOUBLE_EQ(m[4], (4.0 + 5.0 + 5.0) / 3);
    // NaN only where the whole window is not finite; infinities count as missing (window 2: the
    // value and the one before it)
    const std::vector<double> holes = panels::runningMean(std::vector<double>{nan, nan, INFINITY, 3.0, 4.0}, 2);
    EXPECT_TRUE(std::isnan(holes[0]));
    EXPECT_TRUE(std::isnan(holes[1]));
    EXPECT_TRUE(std::isnan(holes[2]));
    EXPECT_DOUBLE_EQ(holes[3], 3.0);
    EXPECT_DOUBLE_EQ(holes[4], 3.5);
    // as MDATK: window <= 1 or fewer than two values: unchanged; the window clipped to the series
    EXPECT_EQ(panels::runningMean(std::vector<double>{1.0, 2.0}, 1), (std::vector<double>{1.0, 2.0}));
    EXPECT_EQ(panels::runningMean(std::vector<double>{7.0}, 5), (std::vector<double>{7.0}));
    const std::vector<double> clipped = panels::runningMean(std::vector<double>{1.0, 3.0}, 9);  // window 2
    EXPECT_DOUBLE_EQ(clipped[0], 1.0);
    EXPECT_DOUBLE_EQ(clipped[1], 2.0);
}

TEST(PlotArithmetic, ColumnEdgesAtMidpoints) {
    // RN-01: two runs of different time steps; MDATK used t[1] - t[0] for every column
    const std::vector<double> t{0.0, 10.0, 20.0, 22.0, 24.0};
    EXPECT_EQ(panels::columnEdges(t), (std::vector<double>{-5.0, 5.0, 15.0, 21.0, 23.0, 25.0}));
    // uniform times: MDATK's cells; one frame: dt = 1 as MDATK
    EXPECT_EQ(panels::columnEdges(std::vector<double>{1.0, 2.0, 3.0}), (std::vector<double>{0.5, 1.5, 2.5, 3.5}));
    EXPECT_EQ(panels::columnEdges(std::vector<double>{4.0}), (std::vector<double>{3.5, 4.5}));
    EXPECT_TRUE(panels::columnEdges(std::vector<double>{}).empty());
}

TEST(PlotArithmetic, HeatmapNormSwapsReversedLimits) {
    // RN-04: vmin above vmax swapped (MDATK: its colour bar turned vmin 5 / vmax 1 into 0.9 / 5)
    const std::vector<double> v{0.0, 2.0, 3.0, std::nan("")};
    panels::Norm n = panels::heatmapNorm(v, 5.0, 1.0);
    EXPECT_TRUE(n.swapped);
    EXPECT_EQ(n.vmin, 1.0);
    EXPECT_EQ(n.vmax, 5.0);
    // vmin above the largest value: swapped with it
    n = panels::heatmapNorm(v, 4.0, std::nullopt);
    EXPECT_TRUE(n.swapped);
    EXPECT_EQ(n.vmin, 3.0);
    EXPECT_EQ(n.vmax, 4.0);
    // as MDATK: vmax 0 or unset is the largest finite value (1 without one), at least 1e-9; a
    // negative-valued matrix with the default vmin 0 is no error (vmax 1e-9)
    n = panels::heatmapNorm(v, std::nullopt, 0.0);
    EXPECT_FALSE(n.swapped);
    EXPECT_EQ(n.vmin, 0.0);
    EXPECT_EQ(n.vmax, 3.0);
    n = panels::heatmapNorm(std::vector<double>{-3.0, -1.0}, std::nullopt, std::nullopt);
    EXPECT_FALSE(n.swapped);
    EXPECT_EQ(n.vmax, 1e-9);
    n = panels::heatmapNorm(std::vector<double>{std::nan("")}, std::nullopt, std::nullopt);
    EXPECT_EQ(n.vmax, 1.0);
}

namespace {

panels::BuiltFigure buildNamed(const std::string& name) {
    const Value rec = test::readJson(test::kFigureCases / (name + ".json"));
    static std::map<std::string, core::Dataset> inputs;
    const std::string in = str(at(rec, "input"));
    if (!inputs.count(in)) inputs.emplace(in, test::loadRenderInput(in));
    static std::map<std::string, core::Project> projects;
    projects.erase(name);
    projects.emplace(name, test::caseProject(rec));
    return test::buildCase(rec, inputs.at(in), projects.at(name));
}

}  // namespace

TEST(PlotFrameLoop, StaticBackgroundAndAnimatedArtists) {
    // render_video's blitting: the animated artists are left out of the background and drawn alone
    // per frame; the layout stays frozen while they move
    for (const char* name : {"ts_progressive", "dual_axis", "heatmap_counts", "contact_instant", "contact_frequency"}) {
        SCOPED_TRACE(name);
        panels::BuiltFigure b = buildNamed(name);
        const std::vector<panels::ArtistRef> animated = b.animatedArtists();
        panels::markAnimated(b.figure, animated);
        const DisplayList all = buildDisplayList(b.figure, test::figureFonts(), DisplayLayer::All);
        const DisplayList bg = buildDisplayList(b.figure, test::figureFonts(), DisplayLayer::Static);
        const DisplayList fg = buildDisplayList(b.figure, test::figureFonts(), DisplayLayer::Animated);
        EXPECT_EQ(bg.items.size() + fg.items.size(), all.items.size());
        EXPECT_EQ(fg.items.empty(), animated.empty());
        for (const DisplayItem& i : fg.items) EXPECT_TRUE(i.animated);
        std::vector<Box> boxes;
        for (std::size_t k = 0; k < b.figure.axes.size(); ++k) boxes.push_back(b.figure.axesBox(k));
        const std::size_t last = b.context->dataset.frameCount() - 1;
        b.update(last, test::figureFonts());
        for (std::size_t k = 0; k < b.figure.axes.size(); ++k) EXPECT_EQ(b.figure.axesBox(k), boxes[k]);
        // the background is the same at every frame
        const DisplayList bg2 = buildDisplayList(b.figure, test::figureFonts(), DisplayLayer::Static);
        EXPECT_EQ(bg2.items.size(), bg.items.size());
        const DisplayList fg2 = buildDisplayList(b.figure, test::figureFonts(), DisplayLayer::Animated);
        EXPECT_EQ(fg2.items.size(), fg.items.size());
    }
}

TEST(PlotFrameLoop, UnevenColumnsDrawnAtTheirWidths) {
    // RN-01: the heatmap of two runs (10 ps, then 2 ps steps) is drawn with columns of their own
    // widths: frames 0-11 five times as wide as frames 12-23
    panels::BuiltFigure b = buildNamed("heatmap_uneven_times");
    const DisplayList dl = buildDisplayList(b.figure, test::figureFonts());
    const ImageItem* im = nullptr;
    for (const DisplayItem& i : dl.items) {
        if (const auto* p = std::get_if<ImageItem>(&i.item)) im = p;
    }
    ASSERT_NE(im, nullptr);
    const Image& source = std::get<Image>(b.figure.axes[0].artists[0]);
    ASSERT_EQ(source.cols, 24);
    EXPECT_NEAR(im->corner1.x - im->corner0.x, im->width, 1.0);
    // the device column of each source column's start: proportional to the time edges
    const std::vector<std::uint8_t> px = source.pixels();
    auto columnOf = [&](int x) {
        for (int c = 0; c < source.cols; ++c) {
            bool same = true;
            for (int r = 0; r < source.rows && same; ++r) {
                for (int ch = 0; ch < 4; ++ch) {
                    same = same && im->rgba[static_cast<std::size_t>((r * im->width + x) * 4 + ch)] ==
                                       px[static_cast<std::size_t>((r * source.cols + c) * 4 + ch)];
                }
            }
            if (same) return c;
        }
        return -1;
    };
    const double span = source.xEdges.back() - source.xEdges.front();
    int checked = 0;
    for (int c = 0; c < source.cols; ++c) {
        // the pixel column at the centre of the source column shows it (unless a neighbour has the
        // same colours)
        const double mid = (source.xEdges[static_cast<std::size_t>(c)] + source.xEdges[static_cast<std::size_t>(c) + 1]) / 2;
        const int x = static_cast<int>((mid - source.xEdges.front()) / span * im->width);
        const int got = columnOf(x);
        if (got < 0) continue;
        bool unique = true;
        for (int o = 0; o < source.cols; ++o) {
            if (o == c) continue;
            bool same = true;
            for (int r = 0; r < source.rows && same; ++r) {
                for (int ch = 0; ch < 4; ++ch) {
                    same = same && px[static_cast<std::size_t>((r * source.cols + o) * 4 + ch)] ==
                                       px[static_cast<std::size_t>((r * source.cols + c) * 4 + ch)];
                }
            }
            unique = unique && !same;
        }
        if (!unique) continue;
        EXPECT_EQ(got, c) << "column " << c << " at pixel " << x;
        ++checked;
    }
    EXPECT_GT(checked, 10);
}
