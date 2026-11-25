// The interface's project data without Qt (src/app/commands/project_view.cpp): the reference frame
// (MDATK's structure_preview), the 3D view's scene built from it with the project's
// representations, and the read-only figure preview (MDATK's GUI canvas: example dataset,
// watermark, size scaled from output.video). The drawing is the program's; packaging/smoke_test.py
// checks its screenshots.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>
#include <variant>
#include <vector>

#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/render/color.hpp"
#include "mdbridge/render/example.hpp"
#include "mdbridge/render/theme.hpp"
#include "project_view.hpp"
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

const char* kBase = R"(format = "mdatk-project"

[inputs]
topology = "demo_dry.parm7"
trajectories = ["demo_dry.nc"]
start = 3
step = 2

[layout]
preset = "equal_2"

[[layout.panels]]
type = "structure"

[[layout.panels]]
type = "timeseries"
series = ["rmsd/Complex"]

[output.video]
width = 1600
height = 900
dpi = 100

[[groups]]
name = "Aptamer"
residues = "1-12"
color = "#2a78d6"

[[groups]]
name = "Target"
residues = "13-24"
color = "#eb6834"
)";

class ProjectView : public ::testing::Test {
protected:
    void SetUp() override {
        fs::copy_file(kData / "topology" / "demo_dry.parm7", work.path() / "demo_dry.parm7");
        fs::copy_file(kData / "io" / "trajectories" / "demo_dry.nc", work.path() / "demo_dry.nc");
    }
    core::Project project(const std::string& extra = {}) {
        const fs::path file = work.path() / "demo.mdatk.toml";
        test::writeFile(file, std::string(kBase) + extra);
        return core::loadProject(file);
    }
    test::TempDir work;
};

template <class P>
std::size_t countOf(const std::vector<render::molecule::Colored>& items) {
    return static_cast<std::size_t>(
        std::count_if(items.begin(), items.end(), [](const render::molecule::Colored& c) { return std::holds_alternative<P>(c.prim); }));
}

TEST_F(ProjectView, ReferenceFrameIsTheFirstAnalysedFrame) {
    const core::Project p = project();
    const app::ReferenceFrame ref = app::loadReferenceFrame(p);
    EXPECT_EQ(ref.frame, 3);  // inputs.start
    EXPECT_FALSE(ref.cleaned);
    EXPECT_EQ(fs::path(ref.topology).filename(), "demo_dry.parm7");
    ASSERT_EQ(ref.dataset.frameCount(), 1U);
    EXPECT_EQ(ref.dataset.frames().front(), 3);
    ASSERT_NE(ref.dataset.result("trace"), nullptr);
    ASSERT_NE(ref.dataset.result("atoms"), nullptr);  // the 3D view draws atoms even without representations
    EXPECT_NE(ref.dataset.result("proximity"), nullptr);  // the default glow (visual.glow_mode) needs it
    EXPECT_EQ(ref.dataset.result("trace")->array("coords")->shape(), (std::vector<std::size_t>{1, 24, 3}));
}

TEST_F(ProjectView, ProximityOnlyWhenColouringOrGlowNeedIt) {
    EXPECT_EQ(app::referenceAnalyses(project("\n[visual]\nglow_mode = \"none\"\n")),
              (std::vector<std::string>{"atoms", "trace"}));
    const core::Project p = project("\n[visual]\nglow_mode = \"none\"\ncolor_mode = \"proximity\"\n");
    EXPECT_EQ(app::referenceAnalyses(p), (std::vector<std::string>{"atoms", "proximity", "trace"}));
    const app::ReferenceFrame ref = app::loadReferenceFrame(p);
    EXPECT_NE(ref.dataset.result("proximity"), nullptr);
    EXPECT_FALSE(app::structureView(p, ref.dataset).items.empty());
}

TEST_F(ProjectView, MissingInputsFail) {
    const core::Project p = project();
    fs::remove(work.path() / "demo_dry.nc");
    EXPECT_THROW(app::loadReferenceFrame(p), std::exception);
}

TEST_F(ProjectView, CleaningWithoutCleanedFilesShowsTheRawFiles) {
    const core::Project p = project("\n[cleaning]\nenabled = true\n");
    const app::ReferenceFrame ref = app::loadReferenceFrame(p);
    EXPECT_FALSE(ref.cleaned);
    ASSERT_FALSE(ref.notes.empty());
    EXPECT_NE(ref.notes.front().find("raw files"), std::string::npos);
    EXPECT_EQ(ref.dataset.frameCount(), 1U);
}

TEST_F(ProjectView, CancelledLoadThrows) {
    const core::Project p = project();
    const core::Reporter reporter;
    reporter.cancel();
    EXPECT_THROW(app::loadReferenceFrame(p, reporter), core::Cancelled);
}

TEST_F(ProjectView, ViewProjectHasOneStructurePanelAndRepresentations) {
    const core::Project p = project();
    const core::Project v = app::structureViewProject(p);
    ASSERT_EQ(v.layout.panels.size(), 1U);
    EXPECT_EQ(v.layout.panels.front().type, "structure");
    EXPECT_EQ(v.layout.preset, "single");
    ASSERT_EQ(v.visual.representations.size(), 1U);
    EXPECT_EQ(v.visual.representations.front().style, app::kDefaultViewStyle);
    EXPECT_EQ(v.visual.renderer, "software");

    const core::Project withReps = project("\n[[visual.representations]]\nstyle = \"licorice\"\n");
    EXPECT_EQ(app::structureViewProject(withReps).visual.representations, withReps.visual.representations);
}

TEST_F(ProjectView, SceneWithoutRepresentationsIsACartoon) {
    const core::Project p = project();
    const app::ReferenceFrame ref = app::loadReferenceFrame(p);
    const app::StructureView view = app::structureView(p, ref.dataset);
    EXPECT_EQ(view.residues, 24U);
    EXPECT_GT(view.atoms, 400U);
    EXPECT_GT(countOf<render::molecule::Mesh>(view.items), 0U);  // cartoon ribbons and slabs
    // a square data window around the centred structure, a symmetric depth range
    const render::molecule::Viewport& vp = view.viewport;
    EXPECT_NEAR(vp.xlim[1] - vp.xlim[0], vp.ylim[1] - vp.ylim[0], 1e-9);
    EXPECT_GT(vp.xlim[1] - vp.xlim[0], 10.0);
    EXPECT_LT(vp.xlim[1] - vp.xlim[0], 500.0);
    EXPECT_DOUBLE_EQ(vp.zrange[0], -vp.zrange[1]);
    // the background is the panel surface, the colour the depth cue fades toward
    const render::Rgba surface = render::toRgba(render::theme(p.visual.theme).surface);
    EXPECT_FLOAT_EQ(view.background[0], static_cast<float>(surface.r));
    EXPECT_FLOAT_EQ(view.background[1], static_cast<float>(surface.g));
    EXPECT_FLOAT_EQ(view.background[2], static_cast<float>(surface.b));
    EXPECT_DOUBLE_EQ(view.shading.fog_rgb[0], surface.r);
    ASSERT_TRUE(view.shading.glow.has_value());
    EXPECT_EQ(view.shading.glow->size(), 24U);
}

TEST_F(ProjectView, SceneDrawsTheProjectsRepresentations) {
    const core::Project p = project(
        "\n[visual]\nrotation = [90.0, 0.0, 0.0]\n\n[[visual.representations]]\nstyle = \"licorice\"\ngroups = [\"Aptamer\"]\n");
    const app::ReferenceFrame ref = app::loadReferenceFrame(p);
    const app::StructureView view = app::structureView(p, ref.dataset);
    EXPECT_GT(countOf<render::molecule::Sticks>(view.items), 0U);  // licorice of the aptamer
    EXPECT_GT(countOf<render::molecule::Spheres>(view.items), 0U);  // the Mg ions
    EXPECT_EQ(countOf<render::molecule::Mesh>(view.items), 0U);      // the target as a trace (rep_rest), no cartoon
    // the project's rotation: 90 degrees about x
    EXPECT_NEAR(view.viewport.rotation[0][0], 1.0, 1e-12);
    EXPECT_NEAR(std::fabs(view.viewport.rotation[1][2]), 1.0, 1e-12);
    EXPECT_NEAR(std::fabs(view.viewport.rotation[2][1]), 1.0, 1e-12);

    const core::Project hidden = project(
        "\n[visual]\nrep_rest = \"hide\"\nshow_ions = false\n\n[[visual.representations]]\nstyle = \"licorice\"\ngroups = [\"Nothing\"]\n");
    const app::ReferenceFrame ref2 = app::loadReferenceFrame(hidden);
    EXPECT_THROW(app::structureView(hidden, ref2.dataset), core::ValueError);  // nothing to draw
}

TEST_F(ProjectView, StructurePreviewFeedsTheExampleDataset) {
    const core::Project p = project();
    const app::ReferenceFrame ref = app::loadReferenceFrame(p);
    const render::StructurePreview preview = app::structurePreview(ref.dataset);
    ASSERT_EQ(preview.trace.coords.size(), 24U);
    EXPECT_EQ(preview.trace.resid.size(), 24U);
    EXPECT_EQ(preview.trace.chain.size(), 24U);
    ASSERT_TRUE(preview.trace.groupMember.has_value());
    EXPECT_EQ(preview.trace.groupMember->size(), 2U);
    EXPECT_EQ(std::count(preview.trace.groupMember->at(0).begin(), preview.trace.groupMember->at(0).end(), true), 12);
    ASSERT_TRUE(preview.atoms.has_value());
    EXPECT_EQ(preview.atoms->coords.size(), preview.atoms->row.size());
    EXPECT_EQ(preview.atoms->names.size(), preview.atoms->coords.size());

    const core::Dataset example = render::exampleDataset(p, &preview);
    EXPECT_EQ(example.frameCount(), static_cast<std::size_t>(render::kExampleFrames));
    const core::AnalysisResult* trace = example.result("trace");
    ASSERT_NE(trace, nullptr);
    EXPECT_EQ(trace->array("coords")->shape()[1], 24U);  // the real structure, not the schematic duplex
    EXPECT_THROW(app::structurePreview(core::Dataset({0}, {0.0})), core::ValueError);
}

TEST(ProjectViewPreview, SizeWatermarkAndPosition) {
    core::Project p;
    p.output.video.width = 1920;
    p.output.video.height = 1080;
    p.output.video.dpi = 120;
    const app::PreviewSize s = app::previewSize(p, 960, 900);
    EXPECT_EQ(s.width, 960);
    EXPECT_EQ(s.height, 540);
    EXPECT_DOUBLE_EQ(s.dpi, 60.0);
    const app::PreviewSize tall = app::previewSize(p, 2000, 270);
    EXPECT_EQ(tall.width, 480);
    EXPECT_EQ(tall.height, 270);
    EXPECT_EQ(app::previewPosition(120), 41U);  // int(0.35 * 119), MDATK's PREVIEW_FRAME
    EXPECT_EQ(app::previewPosition(1), 0U);
    EXPECT_EQ(app::previewPosition(0), 0U);
    EXPECT_EQ(app::previewWatermark(false), "Example data – layout preview");
    EXPECT_EQ(app::previewWatermark(true), "Example data – layout preview (structure: reference frame)");
}

TEST_F(ProjectView, PreviewFigureWithTheReferenceFrame) {
    const core::Project p = project();
    const app::ReferenceFrame ref = app::loadReferenceFrame(p);
    const render::StructurePreview preview = app::structurePreview(ref.dataset);
    const core::Dataset example = render::exampleDataset(p, &preview);
    const std::string mark = app::previewWatermark(true);
    const render::panels::BuiltFigure built = app::previewFigure(p, example, fonts(), 800, 800, mark);
    EXPECT_EQ(built.figure.widthPx(), 800);
    EXPECT_EQ(built.figure.heightPx(), 450);
    EXPECT_DOUBLE_EQ(built.figure.dpi(), 50.0);
    ASSERT_EQ(built.panels.size(), 2U);
    EXPECT_TRUE(built.panels[0]->ok());
    EXPECT_TRUE(built.panels[1]->ok());
    const auto has = [&](const std::string& text) {
        return std::any_of(built.figure.texts.begin(), built.figure.texts.end(),
                           [&](const render::Text& t) { return t.text == text; });
    };
    EXPECT_TRUE(has(mark));
}

TEST(ProjectViewPreview, SchematicPreviewWithoutAStructure) {
    core::Project p;  // the default layout, nothing loaded
    const core::Dataset example = render::exampleDataset(p);
    const render::panels::BuiltFigure built =
        app::previewFigure(p, example, fonts(), 640, 400, app::previewWatermark(false));
    EXPECT_EQ(built.figure.widthPx(), 640);
    EXPECT_EQ(built.figure.heightPx(), 360);
    EXPECT_EQ(built.panels.size(), p.layout.panels.size());
}

TEST(ProjectViewPreview, UnknownPresetIsAnError) {
    core::Project p;
    p.layout.preset = "nonsense";
    const core::Dataset example = render::exampleDataset(p);
    EXPECT_THROW(app::previewFigure(p, example, fonts(), 640, 400, {}), core::ValueError);
}

}  // namespace
}  // namespace mdbridge
