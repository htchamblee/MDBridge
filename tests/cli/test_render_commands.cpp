// The render, charts and preview commands (src/app/commands/render_pipeline.cpp) without Qt: a
// recording RenderBackend stands in for the program's QPainter writers and video loop, so what the
// commands decide is tested here: MDATK's test_cli / test_pipeline (full run, time stamps, the
// params log, the step log, re-rendering from the cache), and MDBridge's improvements (RN-02,
// RN-11 ... RN-16, D8, the density fit before the video; docs/DIFFERENCES.md "Rendering: the render, charts and preview
// commands"). The drawn files themselves are compared with MDATK's by tests/parity/cli_parity.py
// (the render cases of tests/data/cli/cli_cases.json) against the `mdbridge` program.

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "commands.hpp"
#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/render/panels/checks.hpp"
#include "render_pipeline.hpp"
#include "test_support.hpp"

namespace mdbridge {
namespace {

namespace fs = std::filesystem;
using test::kData;

const render::FontSet& fonts() {
    static const render::FontSet set = render::FontSet::load(
        [](std::string_view f) { return test::readFile(fs::path(MDBRIDGE_FONT_DIR) / std::string(f)); });
    return set;
}

// Writes "<format> <width>x<height>" for a figure and a text file for a video; counts the starts.
struct Recorder {
    int starts = 0;
    std::vector<std::string> figures;
    app::RenderBackend backend() {
        app::RenderBackend b;
        b.start = [this]() -> std::shared_ptr<void> {
            ++starts;
            return {};
        };
        b.fonts = []() -> const render::FontSet& { return fonts(); };
        b.writeFigure = [this](const render::Figure& fig, const fs::path& path, const std::string& format) {
            figures.push_back(core::pathToUtf8(path.filename()) + " " + format);
            test::writeFile(path, format + " " + std::to_string(fig.widthPx()) + "x" + std::to_string(fig.heightPx()));
        };
        b.renderVideo = [](const core::Dataset& ds, const core::Project&, const fs::path& path,
                           const core::Reporter& reporter) {
            reporter.log("Rendering video -> " + core::pathToUtf8(path));
            test::writeFile(path, "video of " + std::to_string(ds.frameCount()) + " frames");
            return path;
        };
        b.environment = [] { return core::Table{{"qt", "test"}}; };
        return b;
    }
};

struct CommandRun {
    int exit = 0;
    std::string out;
    std::string err;
};

CommandRun runIn(const fs::path& folder, const std::vector<std::string>& args, const app::RenderBackend* backend) {
    std::vector<const char*> argv{"mdbridge"};
    for (const std::string& a : args) argv.push_back(a.c_str());
    std::ostringstream out, err;
    app::Console console{out, err, false};
    const fs::path previous = fs::current_path();
    fs::current_path(folder);
    CommandRun run;
    try {
        run.exit = app::runCoreCommandLine(static_cast<int>(argv.size()), argv.data(), console, backend);
    } catch (...) {
        fs::current_path(previous);
        throw;
    }
    fs::current_path(previous);
    run.out = out.str();
    run.err = err.str();
    return run;
}

// tests/parity/cli_parity.py RENDER_PROJECT, smaller (no structure panel analyses but the trace)
const char* kProject = R"(format = "mdatk-project"

[inputs]
topology = "demo_dry.parm7"
trajectories = ["demo_dry.nc"]
time_unit = "ns"

[analyses.mg]
enabled = false

[analyses.stacking]
enabled = false

[layout]
preset = "equal_2"

[[layout.panels]]
type = "structure"

[[layout.panels]]
type = "timeseries"
series = ["rmsd/Complex"]

[output]
out_dir = "results"

[output.video]
width = 160
height = 90
dpi = 20

[output.charts]
dpi = 20
formats = ["png", "svg"]
series = ["rmsd/Complex", "energy/total"]

[output.pdb]
frames = "0, 12, last"
times = ""

[[groups]]
name = "Aptamer"
residues = "1-12"

[[groups]]
name = "Target"
residues = "13-24"
)";

class RenderCommands : public ::testing::Test {
protected:
    void SetUp() override {
        fs::copy_file(kData / "topology" / "demo_dry.parm7", work.path() / "demo_dry.parm7");
        fs::copy_file(kData / "io" / "trajectories" / "demo_dry.nc", work.path() / "demo_dry.nc");
        test::writeFile(work.path() / "demo.mdatk.toml", kProject);
        backend = recorder.backend();
    }
    CommandRun run(std::vector<std::string> args) {
        args.insert(args.begin() + 1, {"-p", "demo.mdatk.toml"});
        return runIn(work.path(), args, &backend);
    }
    test::TempDir work;
    Recorder recorder;
    app::RenderBackend backend;
};

std::string stampOf(const std::string& video) {
    std::smatch m;
    static const std::regex stamp{R"(_(\d{4}-\d\d-\d\d_\d{6})\.mp4$)"};
    return std::regex_search(video, m, stamp) ? m[1].str() : std::string();
}

TEST_F(RenderCommands, FullRunStampsEveryFileLogsTheStepAndWritesTheParamsLog) {
    const CommandRun r = run({"render", "--pdb"});
    ASSERT_EQ(r.exit, 0) << r.err;
    const core::Value out = core::parseJson(r.out, "stdout");
    const std::string video = out.find("video")->stringOr();
    const std::string stamp = stampOf(video);
    ASSERT_FALSE(stamp.empty()) << video;
    EXPECT_TRUE(fs::exists(core::pathFromUtf8(video)));
    const core::Array& charts = *out.find("charts")->get<core::Array>();
    ASSERT_EQ(charts.size(), 4u);  // 2 series x png, svg
    for (const core::Value& c : charts) {
        EXPECT_NE(c.stringOr().find("_" + stamp + "."), std::string::npos) << c.stringOr();
        EXPECT_TRUE(fs::exists(core::pathFromUtf8(c.stringOr())));
    }
    const core::Array& pdbs = *out.find("pdbs")->get<core::Array>();
    ASSERT_EQ(pdbs.size(), 3u);
    EXPECT_EQ(core::pathBasename(pdbs[1].stringOr()), "MDATK_frame_00012_" + stamp + ".pdb");
    EXPECT_EQ(recorder.starts, 1);
    // the figures: charts at 8 x 4.5 in, 20 dpi
    EXPECT_EQ(test::readFile(core::pathFromUtf8(charts[0].stringOr())), "png 160x90");
    // log lines
    EXPECT_NE(r.err.find("MDBridge "), std::string::npos);
    EXPECT_NE(r.err.find(" - output directory "), std::string::npos);
    EXPECT_NE(r.err.find("Wrote 4 chart file(s) to "), std::string::npos);
    EXPECT_TRUE(std::regex_search(r.err, std::regex(R"(\nFinished in \d+\.\d s\n$)"))) << r.err;
    // the params log: fixed name, the run's outputs and stamp, MDBridge's environment
    const std::string paramsLog = out.find("params_log")->stringOr();
    EXPECT_EQ(core::pathBasename(paramsLog), "MDATK_params_log.toml");
    const core::Value log = core::parseToml(test::readFile(core::pathFromUtf8(paramsLog)), paramsLog);
    const core::Value& runTable = *log.find("run");
    EXPECT_EQ(runTable.find("timestamp")->stringOr(), stamp);
    EXPECT_EQ(runTable.find("outputs")->find("video")->stringOr(), video);
    EXPECT_EQ(runTable.find("outputs")->find("pdbs")->get<core::Array>()->size(), 3u);
    EXPECT_EQ(runTable.find("outputs")->find("params_log"), nullptr);  // null: not written in TOML
    EXPECT_EQ(runTable.find("environment")->find("qt")->stringOr(), "test");
    EXPECT_NE(runTable.find("environment")->find("mdbridge"), nullptr);
    EXPECT_NE(runTable.find("citation_keys"), nullptr);
    // the step log of the project file
    const std::string project = test::readFile(work.path() / "demo.mdatk.toml");
    EXPECT_NE(project.find("step = \"render\""), std::string::npos);
    EXPECT_NE(project.find("status = \"ok\""), std::string::npos);
    EXPECT_NE(project.find("summary = \"24 frames analysed; 9 file(s) written to "), std::string::npos);
}

TEST_F(RenderCommands, FromCacheRendersWithoutTheTrajectory) {
    ASSERT_EQ(run({"analyze", "-q"}).exit, 0);
    fs::remove(work.path() / "demo_dry.nc");
    const CommandRun r = run({"render", "--from-cache", "--no-video", "-s", "visual.theme=\"dark\""});
    ASSERT_EQ(r.exit, 0) << r.err;
    const core::Value out = core::parseJson(r.out, "stdout");
    EXPECT_EQ(out.find("video"), nullptr);
    EXPECT_EQ(out.find("charts")->get<core::Array>()->size(), 4u);
    const CommandRun p = run({"preview", "--from-cache", "--frame", "5", "--output", "still"});
    ASSERT_EQ(p.exit, 0) << p.err;
    EXPECT_EQ(p.out, "still.png\n");  // savefig's name.png, printed as written
    EXPECT_EQ(test::readFile(work.path() / "still.png"), "png 160x90");
    const CommandRun nothing = run({"render", "--from-cache", "--no-video", "--no-charts"});
    EXPECT_EQ(nothing.out, "{}\n");
}

TEST_F(RenderCommands, ChartsOptionsAndPreviewDefaults) {
    const CommandRun c = run({"charts", "--series", " energy/total ,,", "--formats", ".PDF", "--dpi", "10"});
    ASSERT_EQ(c.exit, 0) << c.err;
    const core::Value parsed = core::parseJson(c.out, "stdout");
    const core::Array& charts = *parsed.find("charts")->get<core::Array>();
    ASSERT_EQ(charts.size(), 1u);
    EXPECT_EQ(test::readFile(core::pathFromUtf8(charts[0].stringOr())), "pdf 80x45");
    const CommandRun p = run({"preview", "--from-cache", "--frame=-4"});
    ASSERT_EQ(p.exit, 0) << p.err;
    ASSERT_FALSE(p.out.empty());
    EXPECT_EQ(core::pathBasename(p.out.substr(0, p.out.size() - 1)), "MDATK_preview_00000.png");
    const CommandRun last = run({"preview", "--from-cache", "--frame", "99", "--output", "out/x.svg"});
    ASSERT_EQ(last.exit, 0) << last.err;
    EXPECT_EQ(test::readFile(work.path() / "out" / "x.svg"), "svg 160x90");
}

TEST_F(RenderCommands, FormatsAreCheckedBeforeTheAnalyses) {  // RN-11
    const CommandRun p = run({"preview", "--output", "a.gif"});
    EXPECT_EQ(p.exit, 1);
    EXPECT_EQ(p.err, "error: preview format 'gif' is not supported (choose from png, svg, pdf, jpg, jpeg, tif, tiff, "
                     "eps, ps)\n");
    EXPECT_FALSE(fs::exists(work.path() / "results"));  // nothing analysed
    const CommandRun c = run({"charts", "--formats", "png,webp"});
    EXPECT_EQ(c.exit, 1);
    EXPECT_NE(c.err.find("  - output.charts.formats: chart format 'webp' is not supported"), std::string::npos) << c.err;
    EXPECT_FALSE(fs::exists(work.path() / "results"));
}

TEST_F(RenderCommands, InvalidSettingsNameTheirFields) {  // RN-12, RN-13, RN-15
    const CommandRun r = run({"render", "-s", "groups.1.color=\"bluish\"", "-s", "visual.chain_palette=\"jet\"", "-s",
                       "layout.panels.1.options={\"ylim\": [0, \"x\"], \"target\": []}"});
    EXPECT_EQ(r.exit, 1);
    EXPECT_EQ(r.err, "error: Invalid parameters:\n"
                     "  - groups[1].color: Invalid RGBA argument: 'bluish'\n"
                     "  - visual.chain_palette: 'jet' is not a bundled colormap; choose from tab10, tab20, tab20b, "
                     "tab20c, Set1, Set2, Set3, Pastel1, Pastel2, Dark2, Accent, Paired, viridis, plasma, inferno, "
                     "magma, cividis (or one of them reversed, with \"_r\")\n"
                     "  - layout.panels[1].options.target: expected a number, got []\n"
                     "  - layout.panels[1].options.ylim: expected two numbers, got [0, 'x']\n");
    EXPECT_EQ(recorder.starts, 0);
    // info lists them with the other problems
    const CommandRun info = run({"info", "-s", "visual.glow_color=\"\""});
    EXPECT_NE(info.out.find("visual.glow_color: Invalid RGBA argument: ''"), std::string::npos) << info.out;
}

TEST(RenderProblems, PanelOptions) {
    core::Project p;
    core::PanelSpec heat;
    heat.type = "heatmap";
    heat.options = {{"vmin", "0.5"}, {"vmax", core::Value()}, {"max_ticks", 2.7}};  // read as the panels read them
    core::PanelSpec map;
    map.type = "contact_map";
    map.options = {{"mode", "frame"}, {"frame", "3"}, {"max_ticks", -1}, {"transpose", "yes"}};
    core::PanelSpec ts;
    ts.type = "timeseries";
    ts.options = {{"smooth", true}, {"ylim", core::Array{1, 2.5}}, {"progressive", core::Array{}}};
    p.layout.panels = {heat, map, ts};
    EXPECT_EQ(render::panels::renderProblems(p),
              (std::vector<std::string>{"layout.panels[1].options.max_ticks: expected a whole number of at least 1, got -1"}));
    p.output.charts.formats = {"png", "GIF", "webp"};
    p.output.charts.color = "C3";
    p.visual.chain_palette = "viridis_r";
    EXPECT_EQ(render::panels::renderProblems(p).front(),
              "output.charts.formats: chart formats 'GIF', 'webp' are not supported (choose from png, svg, pdf, jpg, "
              "jpeg, tif, tiff, eps, ps)");
    EXPECT_EQ(render::panels::renderProblems(p).size(), 2u);
}

TEST(RenderProblems, FigureFileFormats) {
    EXPECT_EQ(render::panels::figureFileFormat("a/b.SVG"), "svg");
    EXPECT_EQ(render::panels::figureFileFormat("a.b/still"), "png");
    EXPECT_EQ(render::panels::figureFileProblem("x.tiff"), std::nullopt);
    EXPECT_EQ(render::panels::figureFileProblem("x.webp", "chart"),
              "chart format 'webp' is not supported (choose from png, svg, pdf, jpg, jpeg, tif, tiff, eps, ps)");
}

TEST_F(RenderCommands, FfmpegPathIsNotUsedAndTheDensityMapIsFitted) {  // D8, ensure_density_fit (M4 D9 retired)
    fs::copy_file(kData / "io" / "maps" / "files" / "crs_000.mrc", work.path() / "map.mrc");
    const CommandRun r = run({"render", "--no-charts", "-s", "output.video.ffmpeg_path=\"/opt/ffmpeg/bin/ffmpeg\"", "-s",
                       "density.map_path=\"map.mrc\""});
    ASSERT_EQ(r.exit, 0) << r.err;
    const std::size_t note = r.err.find("Note: output.video.ffmpeg_path (/opt/ffmpeg/bin/ffmpeg) is not used; MDBridge "
                                        "encodes with its bundled FFmpeg libraries\n");
    const std::size_t fit = r.err.find("Density fit (global search) of frame ref into map.mrc at 5 \xc3\x85\n");
    const std::size_t fitted = r.err.find("Density fit: global search, frame 0, 511 atoms, 5 \xc3\x85: CC ");
    const std::size_t video = r.err.find("Rendering video -> ");
    ASSERT_NE(note, std::string::npos) << r.err;
    ASSERT_NE(fit, std::string::npos) << r.err;
    ASSERT_NE(fitted, std::string::npos) << r.err;
    EXPECT_LT(note, fit);
    EXPECT_LT(fitted, video);
    EXPECT_EQ(r.err.find("no placement (run density fit)"), std::string::npos);  // milestone 4's D9 warning is gone
    // the fit is not saved: the next render fits again
    const CommandRun again = run({"render", "--no-charts", "-s", "density.map_path=\"map.mrc\""});
    EXPECT_NE(again.err.find("Density fit (global search)"), std::string::npos);

    const CommandRun off = run({"render", "--no-charts", "-s", "density.map_path=\"map.mrc\"", "-s", "density.auto_fit=false"});
    EXPECT_NE(off.err.find("WARNING: the density map has no placement and automatic fitting is off; it is not drawn"),
              std::string::npos);
    EXPECT_EQ(off.err.find("Density fit ("), std::string::npos);
    EXPECT_EQ(off.err.find("Note: output.video.ffmpeg_path"), std::string::npos);

    // a fit that fails is a warning, and the video is drawn without the map
    test::writeFile(work.path() / "bad.mrc", "not a map");
    const CommandRun bad = run({"render", "--no-charts", "-s", "density.map_path=\"bad.mrc\""});
    ASSERT_EQ(bad.exit, 0) << bad.err;
    const std::size_t warning = bad.err.find("WARNING: density map not drawn: density map: ");
    ASSERT_NE(warning, std::string::npos) << bad.err;
    EXPECT_LT(warning, bad.err.find("Rendering video -> "));
}

TEST_F(RenderCommands, WithoutABackend) {
    const CommandRun r = runIn(work.path(), {"charts", "-p", "demo.mdatk.toml"}, nullptr);
    EXPECT_EQ(r.exit, 1);
    EXPECT_NE(r.err.find("error: figures cannot be drawn by this program"), std::string::npos) << r.err;
    app::RenderBackend noVideo = recorder.backend();
    noVideo.renderVideo = {};
    const CommandRun v = runIn(work.path(), {"render", "-p", "demo.mdatk.toml", "--no-charts"}, &noVideo);
    EXPECT_EQ(v.exit, 1);
    EXPECT_NE(v.err.find("error: video encoding is not available in this build\n"), std::string::npos) << v.err;
    EXPECT_NE(test::readFile(work.path() / "demo.mdatk.toml").find("summary = \"VideoError: video encoding is not "
                                                                     "available in this build\""),
              std::string::npos);
}

TEST(RenderPipeline, NoAnalysedFramesToRender) {  // RN-02
    test::TempDir work;
    core::Project project;
    project.output.out_dir = core::pathToUtf8(work.path());
    Recorder recorder;
    const app::RenderBackend backend = recorder.backend();
    analysis::Pipeline analyses(project);
    app::RenderPipeline pipeline(analyses, backend);
    const core::Dataset empty({}, {});
    for (const auto& step : std::vector<std::function<void()>>{
             [&] { pipeline.renderVideo(empty); }, [&] { pipeline.renderCharts(empty); },
             [&] { pipeline.renderPreview(empty, 0, core::pathToUtf8(work.path() / "p.png")); }}) {
        try {
            step();
            ADD_FAILURE() << "no error";
        } catch (const analysis::PipelineError& e) {
            EXPECT_STREQ(e.what(), "no analysed frames to render");
        }
    }
    EXPECT_EQ(recorder.starts, 0);
    EXPECT_TRUE(recorder.figures.empty());
}

TEST(RenderPipeline, StampedNames) {
    EXPECT_EQ(app::stamped("MDATK_video.mp4", "2026-10-15_143205"), "MDATK_video_2026-10-15_143205.mp4");
    EXPECT_EQ(app::stamped("clip", "s"), "clip_s");
    EXPECT_EQ(app::stamped(".hidden", "s"), ".hidden_s");
    EXPECT_EQ(app::stamped("a.b.mp4", ""), "a.b.mp4");
}

TEST(Citations, VtkIsNeverCited) {  // RN-16
    core::Project p;
    p.visual.renderer = "vtk";
    const std::vector<std::string> ids = core::citations::pipelineResources(p);
    EXPECT_EQ(std::find(ids.begin(), ids.end(), "vtk"), ids.end());
}

}  // namespace
}  // namespace mdbridge
