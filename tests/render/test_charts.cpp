// The standalone charts (src/render/panels/charts.cpp) against MDATK's render_charts
// (data/render/charts, tests/parity/make_chart_cases.py, with MDBridge's improvements applied by
// improvements_render.py): every case runs renderCharts on its input dataset and project into a
// temporary folder with a writer that records each drawn figure and leaves an empty file; then
// the paths returned, the files in the folder, the log, the progress, the warnings and the error
// must be MDATK's exactly, and every chart's figure equal by its contents (figure_compare.hpp:
// layout within 1.5 %, data, colours, ticks, legends, colour bars, the images' arrays). The files'
// bytes are the Qt backend's (test_charts_qt.cpp).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "figure_compare.hpp"
#include "figure_records.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/render/figure_dump.hpp"
#include "mdbridge/render/panels/charts.hpp"

using namespace mdbridge;
using namespace mdbridge::render;
using core::Array;
using core::Value;
using test::FigureComparison;

namespace {

const std::filesystem::path kCharts = test::kRenderData / "charts";

const Value& at(const Value& v, std::string_view key) { return FigureComparison::at(v, key); }
std::string str(const Value& v) { return FigureComparison::str(v); }

std::vector<std::string> strings(const Value& v) {
    std::vector<std::string> out;
    if (const Array* a = v.get<Array>()) {
        for (const Value& x : *a) out.push_back(str(x));
    }
    return out;
}

const core::Dataset& input(const std::string& name) {
    static std::map<std::string, core::Dataset> cache;
    auto it = cache.find(name);
    if (it == cache.end()) {
        // the chart generator's own inputs (charts/inputs/<name>), else the figures' (render/inputs)
        const std::filesystem::path own = kCharts / "inputs" / name;
        it = cache.emplace(name, std::filesystem::is_directory(own) ? test::loadRenderInputFolder(own)
                                                                    : test::loadRenderInput(name))
                 .first;
    }
    return it->second;
}

std::vector<std::string> caseNames() {
    std::vector<std::string> out;
    if (!std::filesystem::is_directory(kCharts)) return out;
    for (const auto& e : std::filesystem::directory_iterator(kCharts)) {
        if (e.path().extension() == ".json") out.push_back(e.path().stem().string());
    }
    std::sort(out.begin(), out.end());
    return out;
}

// every log record ("LEVEL message", as the generator's handler records MDATK's)
struct LogCapture {
    std::vector<std::string> lines;
    core::LogSink previous;
    LogCapture() {
        previous = core::setLogSink([this](core::LogLevel level, std::string_view message) {
            if (level >= core::LogLevel::Info) lines.push_back(std::string(core::logLevelName(level)) + " " + std::string(message));
        });
    }
    ~LogCapture() { core::setLogSink(previous); }
};

struct TempDir {
    std::filesystem::path path;
    explicit TempDir(const std::string& name) {
        path = std::filesystem::temp_directory_path() / ("mdbridge-charts-" + name);
        std::filesystem::remove_all(path);
        std::filesystem::create_directories(path);
    }
    ~TempDir() {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
    }
    // a path inside as the records write it ({out}/...)
    std::string masked(std::string p) const {
        const std::string root = path.string();
        for (std::size_t at = p.find(root); at != std::string::npos; at = p.find(root, at + 5)) p.replace(at, root.size(), "{out}");
        std::replace(p.begin(), p.end(), '\\', '/');
        return p;
    }
};

// What the writer was given: each figure once (in the order written) and its files.
struct Written {
    struct Chart {
        Figure figure{1, 1, 72.0};
        std::vector<std::pair<std::string, std::string>> files;  // path, format
    };
    std::vector<Chart> charts;
};

class ChartCase : public ::testing::TestWithParam<std::string> {};

TEST_P(ChartCase, AsMdatk) {
    const Value rec = test::readJson(kCharts / (GetParam() + ".json"));
    const std::string inputName = str(at(rec, "input"));
    const core::Dataset& ds = input(inputName);
    const core::Project project = test::caseProject(rec);
    TempDir tmp(GetParam());
    const std::filesystem::path out = tmp.path / "out";
    std::vector<std::string> logs;
    std::vector<std::pair<double, std::string>> progress;
    bool next = true;  // the next file starts a chart (renderCharts reports progress after each)
    const core::Reporter reporter([&](const std::string& m) { logs.push_back(m); },
                                  [&](double f, const std::string& m) {
                                      progress.emplace_back(f, m);
                                      next = true;
                                  });
    Written written;
    const charts::FigureWriter writer = [&](const Figure& fig, const std::filesystem::path& path, const std::string& format) {
        if (next) {
            written.charts.push_back({fig, {}});
            next = false;
        }
        written.charts.back().files.emplace_back(path.string(), format);
        EXPECT_TRUE(std::ofstream(path, std::ios::binary).good()) << path;
    };
    LogCapture log;
    std::vector<std::string> paths;
    if (const Value* error = FigureComparison::opt(rec, "error")) {
        try {
            (void)charts::renderCharts(ds, project, out, reporter, str(at(rec, "stamp")), test::figureFonts(), writer);
            ADD_FAILURE() << "no error; MDATK: " << str(at(*error, "message"));
        } catch (const core::ValueError& e) {
            EXPECT_EQ(str(at(*error, "type")), "ValueError");
            EXPECT_EQ(e.what(), str(at(*error, "message")));
        }
    } else {
        for (const std::string& p : charts::renderCharts(ds, project, out, reporter, str(at(rec, "stamp")),
                                                         test::figureFonts(), writer)) {
            paths.push_back(tmp.masked(p));
        }
        EXPECT_EQ(paths, strings(at(rec, "paths")));
    }
    // the log, progress and warnings
    std::vector<std::string> mine;
    for (const std::string& l : logs) mine.push_back(tmp.masked(l));
    EXPECT_EQ(mine, strings(at(rec, "log")));
    mine.clear();
    for (const std::string& l : log.lines) mine.push_back(tmp.masked(l));
    EXPECT_EQ(mine, strings(at(rec, "messages")));
    const Array& tp = *at(rec, "progress").get<Array>();
    ASSERT_EQ(progress.size(), tp.size());
    for (std::size_t i = 0; i < tp.size(); ++i) {
        const Array& p = *tp[i].get<Array>();
        EXPECT_EQ(progress[i].first, std::strtod(str(p[0]).c_str(), nullptr)) << i;
        EXPECT_EQ(progress[i].second, str(p[1])) << i;
    }
    // the files left in the folder
    std::vector<std::string> files;
    if (std::filesystem::is_directory(out)) {
        for (const auto& e : std::filesystem::directory_iterator(out)) files.push_back(e.path().filename().string());
    }
    std::sort(files.begin(), files.end());
    EXPECT_EQ(files, strings(at(rec, "files")));

    // every chart's figure and files
    const Array& charts = *at(rec, "charts").get<Array>();
    ASSERT_EQ(written.charts.size(), charts.size());
    const std::filesystem::path npzPath = kCharts / (GetParam() + ".npz");
    const auto npz = std::filesystem::exists(npzPath) ? test::readStoredNpz(npzPath) : decltype(test::readStoredNpz(npzPath)){};
    double worst = 0.0;
    std::string worstWhat;
    int ticks = 0, skipped = 0, arrays = 0;
    for (std::size_t k = 0; k < charts.size(); ++k) {
        SCOPED_TRACE("chart " + std::to_string(k) + " " + str(at(charts[k], "key")));
        const Written::Chart& w = written.charts[k];
        const Array& tf = *at(charts[k], "files").get<Array>();
        ASSERT_EQ(w.files.size(), tf.size());
        for (std::size_t i = 0; i < tf.size(); ++i) {
            EXPECT_EQ(tmp.masked(w.files[i].first), str(at(tf[i], "path")));
            EXPECT_EQ(w.files[i].second, str(at(tf[i], "format")));
        }
        const Value dump = dumpFigure(w.figure, test::figureFonts());
        FigureComparison c;
        c.width = w.figure.widthPx();
        c.height = w.figure.heightPx();
        // at least 8 px: in the 200 px wide charts at 50 dpi a title's 7 px glyphs are hinted up to 0.3 px
        // wider each by Matplotlib (IBM Plex Serif, FONT-01; MDBridge measures the unhinted outline)
        c.minBoxPixels = 8.0;
        c.figure(dump, at(charts[k], "figure"));
        arrays += test::compareArrays(w.figure, npz, "c" + std::to_string(k) + "/");
        if (c.worstBox > worst) {
            worst = c.worstBox;
            worstWhat = c.worstWhat;
        }
        ticks += c.ticksCompared;
        skipped += c.ticksSkipped;
    }
    std::printf("[chart case] %-24s %zu chart(s); worst box %.3f %% (%s); ticks %d (nbins differ %d); %d arrays\n",
                GetParam().c_str(), charts.size(), 100.0 * worst, worstWhat.c_str(), ticks, skipped, arrays);
}

INSTANTIATE_TEST_SUITE_P(Charts, ChartCase, ::testing::ValuesIn(caseNames()),
                         [](const ::testing::TestParamInfo<std::string>& i) { return i.param; });

}  // namespace

TEST(ChartNames, SafeKeysAndTitlesAsMdatk) {
    EXPECT_EQ(charts::safeKey("rmsd/Complex"), "rmsd_Complex");
    EXPECT_EQ(charts::safeKey("rmsd@Nope/x"), "rmsd_Nope_x");
    EXPECT_EQ(charts::safeKey("/a//b c+d-e.f/"), "a_b_c+d-e.f");
    EXPECT_EQ(charts::safeKey("rmsd/\xc3\x85ngstr\xc3\xb6m"), "rmsd_ngstr_m");  // a run of non-ASCII: one "_"
    EXPECT_EQ(charts::safeKey("__"), "");
    EXPECT_EQ(charts::formatTitle("", "RMSD", "rmsd/x"), "RMSD");
    EXPECT_EQ(charts::formatTitle("{label} of {key} ({label})", "RMSD", "rmsd/x"), "RMSD of rmsd/x (RMSD)");
    EXPECT_EQ(charts::normalizeFormat("..PNG"), "png");
}

TEST(ChartFormats, CheckedBeforeAnythingIsWritten) {
    // RN-11: MDATK failed with Matplotlib's ValueError at the first save in an unknown format,
    // after the files before it were written
    EXPECT_FALSE(charts::chartFormatProblem({"png", "SVG", ".pdf", "jpg", "JPEG", "tif", "tiff", "eps", "ps"}));
    EXPECT_EQ(charts::chartFormatProblem({"png", "bogus"}).value_or(""),
              "chart format 'bogus' is not supported (choose from png, svg, pdf, jpg, jpeg, tif, tiff, eps, ps)");
    EXPECT_EQ(charts::chartFormatProblem({"gif", "it's"}).value_or(""),
              "chart formats 'gif', \"it's\" are not supported (choose from png, svg, pdf, jpg, jpeg, tif, tiff, eps, ps)");
}

TEST(ChartHeight, CappedAtTwoToTheFifteen) {
    // RN-05: 0.18 in per row; MDATK's images grew past every raster's limits
    EXPECT_EQ(charts::matrixChartHeightPx(4.5, 12, 80.0), 360);
    EXPECT_EQ(charts::matrixChartHeightPx(4.5, 24, 80.0), 465);  // 5.82 in: int(465.6)
    EXPECT_EQ(charts::matrixChartHeightPx(4.5, 700, 300.0), charts::kMaxChartPixels);
    EXPECT_EQ(charts::matrixChartHeightPx(4.5, 100000, 72.0), 32768);
}

TEST(ChartRun, CancelledBetweenCharts) {
    const std::string demo = "demo_dry";
    const core::Dataset& ds = input(demo);
    core::Project p;
    p.output.charts.series = {"rmsd/Complex", "hbonds/count"};
    p.output.charts.dpi = 20;
    TempDir tmp("cancel");
    int written = 0;
    const core::Reporter reporter;
    const charts::FigureWriter writer = [&](const Figure&, const std::filesystem::path&, const std::string&) {
        ++written;
        reporter.cancel();
    };
    EXPECT_THROW(charts::renderCharts(ds, p, tmp.path / "out", reporter, "", test::figureFonts(), writer), core::Cancelled);
    EXPECT_EQ(written, 1);
}

TEST(ChartRun, MatrixSeriesAndPairMaps) {
    // RN-10: requested matrices are charted after the series (each once, before the other
    // heatmaps); a pair map warns (MDATK dropped both silently)
    const std::string demo = "demo_dry";
    const core::Dataset& ds = input(demo);
    core::Project p;
    p.output.charts.series = {"contacts/map", "mg/occupancy", "rmsd/Complex", "mg/occupancy"};
    p.output.charts.include_heatmaps = true;
    LogCapture log;
    const std::vector<charts::ChartJob> jobs = charts::chartJobs(ds, p);
    ASSERT_GE(jobs.size(), 3u);
    EXPECT_EQ(jobs[0].kind, "series");
    EXPECT_EQ(jobs[0].key, "rmsd/Complex");
    EXPECT_EQ(jobs[1].kind, "matrix");
    EXPECT_EQ(jobs[1].key, "mg/occupancy");
    std::set<std::string> keys;
    for (const charts::ChartJob& j : jobs) EXPECT_TRUE(keys.insert(j.key).second) << j.key;
    for (const charts::ChartJob& j : jobs) EXPECT_FALSE(j.key.starts_with("proximity/")) << j.key;
    EXPECT_EQ(log.lines, (std::vector<std::string>{"WARNING chart series contacts/map is a pair map; charts draw series "
                                                   "and matrices (a contact_map panel shows it)"}));
}
