// The counterpart of MDATK's tests/test_cli.py test_cli_demo_render_and_tools: the demo system
// through list, info, render, charts and preview from the cache, export-pdb and init, in-process
// with a recording RenderBackend (what the commands decide; the drawn files of `demo --run` are
// compared with MDATK's by cli_parity.py's demo_run case against the `mdbridge` program).
#include <gtest/gtest.h>

#include <filesystem>
#include <regex>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "commands.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/render/text.hpp"
#include "render_pipeline.hpp"
#include "test_support.hpp"

namespace mdbridge {
namespace {

namespace fs = std::filesystem;

const render::FontSet& demoFonts() {
    static const render::FontSet set = render::FontSet::load(
        [](std::string_view f) { return test::readFile(fs::path(MDBRIDGE_FONT_DIR) / std::string(f)); });
    return set;
}

app::RenderBackend recordingBackend() {
    app::RenderBackend b;
    b.start = []() -> std::shared_ptr<void> { return {}; };
    b.fonts = []() -> const render::FontSet& { return demoFonts(); };
    b.writeFigure = [](const render::Figure& fig, const fs::path& path, const std::string& format) {
        test::writeFile(path, format + " " + std::to_string(fig.widthPx()) + "x" + std::to_string(fig.heightPx()));
    };
    b.renderVideo = [](const core::Dataset& ds, const core::Project&, const fs::path& path, const core::Reporter&) {
        test::writeFile(path, "video of " + std::to_string(ds.frameCount()) + " frames");
        return path;
    };
    return b;
}

struct CommandResult {
    int exit = 0;
    std::string out, err;
};

CommandResult run(const fs::path& folder, const std::vector<std::string>& args, const app::RenderBackend& backend) {
    std::vector<const char*> argv{"mdbridge"};
    for (const std::string& a : args) {
        argv.push_back(a.c_str());
    }
    std::ostringstream out, err;
    app::Console console{out, err, false};
    const fs::path previous = fs::current_path();
    fs::current_path(folder);
    CommandResult r;
    try {
        r.exit = app::runCoreCommandLine(static_cast<int>(argv.size()), argv.data(), console, &backend);
    } catch (...) {
        fs::current_path(previous);
        throw;
    }
    fs::current_path(previous);
    r.out = out.str();
    r.err = err.str();
    return r;
}

const core::Value& at(const core::Value& v, std::string_view key) {
    const core::Value* f = v.find(key);
    if (f == nullptr) {
        throw std::runtime_error("missing " + std::string(key));
    }
    return *f;
}

TEST(DemoCommands, DemoRenderAndTools) {
    test::TempDir tmp;
    const app::RenderBackend backend = recordingBackend();
    const fs::path d = tmp / "demo";
    ASSERT_EQ(run(tmp.path(), {"demo", core::pathToUtf8(d), "--n-frames", "10", "--no-solvent"}, backend).exit, 0);
    const std::string preset = core::pathToUtf8(d / "demo.mdatk.toml");

    const CommandResult list = run(tmp.path(), {"list"}, backend);
    ASSERT_EQ(list.exit, 0);
    EXPECT_NE(list.out.find("contact_map"), std::string::npos);

    const CommandResult info = run(tmp.path(), {"info", "-p", preset}, backend);
    ASSERT_EQ(info.exit, 0) << info.err;
    const core::Value summary = core::parseJson(info.out, "info");
    EXPECT_EQ(at(at(at(summary, "groups"), "Aptamer"), "residues"), core::Value(std::int64_t{12}));
    EXPECT_TRUE(at(summary, "problems").get<core::Array>()->empty());

    const std::vector<std::string> args{"-p", preset, "-q", "--set", "output.video.width=320", "--set",
                                        "output.video.height=180", "--set", "output.video.dpi=40", "--set",
                                        "output.charts.dpi=40", "--set", "output.charts.series=[\"hbonds/count\"]",
                                        "-o", core::pathToUtf8(tmp / "out")};
    const auto with = [&args](std::vector<std::string> head) {
        head.insert(head.end(), args.begin(), args.end());
        return head;
    };
    const CommandResult render = run(tmp.path(), with({"render"}), backend);
    ASSERT_EQ(render.exit, 0) << render.err;
    const core::Value out = core::parseJson(render.out, "render");
    const std::string video = at(out, "video").stringOr();
    EXPECT_TRUE(fs::exists(core::pathFromUtf8(video)));
    EXPECT_EQ(at(out, "charts").get<core::Array>()->size(), 1u);
    EXPECT_FALSE(at(out, "pdbs").get<core::Array>()->empty());
    const core::Array log = core::readLog(core::pathFromUtf8(preset));  // the run is logged in the project file
    ASSERT_FALSE(log.empty());
    EXPECT_EQ(at(log.back(), "step").stringOr(), "render");
    EXPECT_EQ(at(log.back(), "status").stringOr(), "ok");
    bool listed = false;
    for (const core::Value& f : *at(log.back(), "files").get<core::Array>()) {
        listed = listed || f.stringOr() == video;
    }
    EXPECT_TRUE(listed);

    const CommandResult charts = run(tmp.path(), with({"charts", "--from-cache", "--formats", "svg"}), backend);
    ASSERT_EQ(charts.exit, 0) << charts.err;
    EXPECT_TRUE(at(core::parseJson(charts.out, "charts"), "charts").get<core::Array>()->at(0).stringOr().ends_with(".svg"));

    const CommandResult preview = run(tmp.path(), with({"preview", "--from-cache", "--frame", "3", "--output",
                                              core::pathToUtf8(tmp / "p.png")}),
                            backend);
    ASSERT_EQ(preview.exit, 0) << preview.err;
    EXPECT_TRUE(fs::exists(tmp / "p.png"));

    const CommandResult pdb = run(tmp.path(), with({"export-pdb", "--frame-list", "1", "--times", "", "--selection", "Target"}),
                        backend);
    ASSERT_EQ(pdb.exit, 0) << pdb.err;
    const core::Value exported = core::parseJson(pdb.out, "export");
    const core::Array& pdbs = *at(exported, "pdbs").get<core::Array>();
    ASSERT_EQ(pdbs.size(), 1u);
    const std::string name = core::pathToUtf8(core::pathFromUtf8(pdbs[0].stringOr()).filename());
    EXPECT_TRUE(std::regex_match(name, std::regex(R"(MDATK_frame_00001_\d{4}-\d\d-\d\d_\d{6}\.pdb)"))) << name;
    const std::string text = test::readFile(core::pathFromUtf8(pdbs[0].stringOr()));
    std::size_t atoms = 0;
    for (std::size_t pos = text.find("ATOM"); pos != std::string::npos; pos = text.find("ATOM", pos + 4)) {
        ++atoms;
    }
    EXPECT_GT(atoms, 100u);
    EXPECT_EQ(text.substr(0, text.find("ATOM")).find(" U  "), std::string::npos);

    const fs::path init = tmp / "init.mdatk.toml";
    ASSERT_EQ(run(tmp.path(), {"init", core::pathToUtf8(init), "--top", "x.parm7"}, backend).exit, 0);
    EXPECT_TRUE(core::loadProject(init).inputs.topology.ends_with("x.parm7"));
}

}  // namespace
}  // namespace mdbridge
