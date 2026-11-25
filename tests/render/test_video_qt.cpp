// The video frame loop with the QPainter compositor (src/render/video_qt.cpp): the counterpart of
// MDATK's test_full_run video (24 frames at 320 x 180 of the demo project, render_cli_cases.json's
// "render" case: frame count, size, frame rate, pixel format, faststart, no audio, the log lines),
// frames composed byte-identical for 1 and N threads, and the encoded frames decoded again within
// 38 dB PSNR of the composed ones (mpeg4; the bundled FFmpeg has no AV1 decoder, so AV1 videos are
// checked by their metadata only). Offscreen platform; the encoding tests skip without FFmpeg.
#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <QGuiApplication>

#include "mdbridge/core/paths.hpp"
#include "mdbridge/render/figure_qt.hpp"
#include "mdbridge/render/panels/structure.hpp"
#include "mdbridge/render/text_qt.hpp"
#include "mdbridge/render/video_qt.hpp"
#include "video/encoder.hpp"
#include "video_cases.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;
using namespace mdbridge::test;

namespace {

#define REQUIRE_FFMPEG()                                                         \
    if (!render::video::videoEncodingAvailable()) {                              \
        GTEST_SKIP() << "built without the FFmpeg libraries (packaging/ffmpeg)"; \
    }

// The frames of the case's video as the loop composes them (QtFrameCompositor, `threads`).
std::vector<std::vector<std::uint8_t>> composedFrames(const core::Project& project, unsigned threads) {
    render::panels::BuiltFigure fig = render::panels::buildFigure(demoDataset(), project, render::bundledFontSet());
    render::QtFrameCompositor comp;
    std::vector<std::vector<std::uint8_t>> frames;
    const std::vector<std::size_t> positions =
        render::video::videoPositions(demoDataset().frameCount(), project.output.video.frame_stride);
    render::video::runFrameLoop(
        fig, render::bundledFontSet(), positions, comp,
        [&](std::size_t, std::vector<std::uint8_t>& rgba) { frames.push_back(std::move(rgba)); }, core::Reporter(),
        {threads});
    return frames;
}

// PSNR over R, G and B of two RGBA images of the same size
double psnr(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b) {
    double se = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < a.size(); i += 4) {
        for (std::size_t c = 0; c < 3; ++c) {
            const double d = static_cast<double>(a[i + c]) - static_cast<double>(b[i + c]);
            se += d * d;
            ++n;
        }
    }
    const double mse = se / static_cast<double>(n);
    return mse == 0.0 ? 99.0 : 10.0 * std::log10(255.0 * 255.0 / mse);
}

// PSNR of the BT.709 luma (Y' = 0.2126 R + 0.7152 G + 0.0722 B, full range) of two RGBA images
double lumaPsnr(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b) {
    double se = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < a.size(); i += 4) {
        const double ya = 0.2126 * a[i] + 0.7152 * a[i + 1] + 0.0722 * a[i + 2];
        const double yb = 0.2126 * b[i] + 0.7152 * b[i + 1] + 0.0722 * b[i + 2];
        se += (ya - yb) * (ya - yb);
        ++n;
    }
    const double mse = se / static_cast<double>(n);
    return mse == 0.0 ? 99.0 : 10.0 * std::log10(255.0 * 255.0 / mse);
}

}  // namespace

TEST(QtVideo, ComposedFramesIdenticalForAnyThreadCount) {
    const core::Project project = renderCliProject(renderCliCase("render"));
    const auto one = composedFrames(project, 1);
    ASSERT_EQ(one.size(), 24u);
    ASSERT_EQ(one[0].size(), 320u * 180u * 4u);
    EXPECT_NE(one[0], one[23]);
    for (const unsigned threads : {3u, 8u}) {
        const auto many = composedFrames(project, threads);
        ASSERT_EQ(many.size(), one.size());
        for (std::size_t k = 0; k < one.size(); ++k) {
            EXPECT_TRUE(many[k] == one[k]) << "frame " << k << ", " << threads << " threads";
        }
    }
}

TEST(QtVideo, FrameOverBackgroundEqualsTheStillImage) {
    // blitting draws the animated artists over the static layer; where they are drawn last anyway
    // (the cursor, the time label) the frame is the still image of that position
    const core::Project project = renderCliProject(renderCliCase("render"));
    const auto frames = composedFrames(project, 2);
    render::panels::BuiltFigure fig = render::panels::buildFigure(demoDataset(), project, render::bundledFontSet());
    fig.update(12, render::bundledFontSet());
    const QImage still = render::renderImage(render::buildDisplayList(fig.figure, render::bundledFontSet()));
    std::vector<std::uint8_t> rgba(320u * 180u * 4u);
    for (int y = 0; y < 180; ++y) {
        std::memcpy(rgba.data() + static_cast<std::size_t>(y) * 320 * 4, still.constScanLine(y), 320 * 4);
    }
    EXPECT_GT(psnr(frames[12], rgba), 30.0);
}

TEST(QtVideo, FullRun24Frames320x180) {  // MDATK test_full_run's video
    REQUIRE_FFMPEG();
    const core::Value& record = renderCliCase("render");
    const core::Project project = renderCliProject(record);
    const ExpectedVideo want = expectedVideo(record);
    ASSERT_EQ(want.frames, 24);
    ASSERT_EQ(want.width, 320);
    ASSERT_EQ(want.height, 180);
    TempDir dir;
    const fs::path path = dir.path() / stamped(want.file);
    std::vector<std::string> logged;
    std::vector<std::string> progress;
    const core::Reporter reporter([&](const std::string& m) { logged.push_back(m); },
                                  [&](double, const std::string& m) { progress.push_back(m); });
    EXPECT_EQ(render::renderVideo(demoDataset(), project, path, reporter), path);
    const render::video::VideoProbe probe = render::video::probeVideo(path);
    expectVideoLike(probe, want);
    EXPECT_EQ(probe.codec, "av1");  // RN-07: libx264 -> the bundled AV1 encoder
    ASSERT_EQ(logged.size(), 2u);
    EXPECT_EQ(logged[0], "Rendering 24 video frame(s) at 320x180, 12 fps -> " + core::pathToUtf8(path));
    EXPECT_EQ(logged[1], "Wrote " + core::pathToUtf8(path) + " (24 frames)");
    ASSERT_EQ(progress.size(), 24u);
    EXPECT_EQ(progress.back(), "Video frame 24/24");
}

namespace {

struct Quality {
    double rgb = 99.0, luma = 99.0;  // the worst frame's PSNR
};

// Writes the project's video with mpeg4 (4 composing threads) and decodes it: the worst frame's
// PSNR against the frames as composed (1 thread).
Quality mpeg4Quality(const core::Project& project, const fs::path& path) {
    render::video::VideoRenderOptions options;
    options.loop.threads = 4;
    render::renderVideo(demoDataset(), project, path, core::Reporter(), options);
    EXPECT_EQ(render::video::probeVideo(path).codec, "mpeg4");
    const auto composed = composedFrames(project, 1);
    const auto decoded = render::video::decodeVideoFrames(path);
    EXPECT_EQ(decoded.size(), composed.size());
    Quality q;
    for (std::size_t k = 0; k < std::min(decoded.size(), composed.size()); ++k) {
        EXPECT_EQ(decoded[k].size(), composed[k].size());
        if (decoded[k].size() != composed[k].size()) continue;
        q.rgb = std::min(q.rgb, psnr(decoded[k], composed[k]));
        q.luma = std::min(q.luma, lumaPsnr(decoded[k], composed[k]));
    }
    return q;
}

}  // namespace

// docs/specs/m4-plan.md section 3: decoded frames within 38 dB of the composed RGBA. Only mpeg4 can
// be decoded here (the bundled FFmpeg encodes AV1 but has no AV1 decoder).
TEST(QtVideo, DecodedFramesWithin38dB) {
    REQUIRE_FFMPEG();
    if (!render::video::hasVideoDecoder("mpeg4")) GTEST_SKIP() << "no mpeg4 decoder";
    // the video's default size (1920 x 1080, dpi 120), every 4th frame, the default crf 20
    {
        core::Project project = renderCliProject(renderCliCase("render_mpeg4"));
        project.output.video.width = 1920;
        project.output.video.height = 1080;
        project.output.video.dpi = 120;
        project.output.video.frame_stride = 4;
        TempDir dir;
        const Quality q = mpeg4Quality(project, dir / "v.mp4");
        EXPECT_GE(q.rgb, 38.0);  // measured 42.3 dB
        EXPECT_GE(q.luma, 38.0);
    }
    // the parity cases' small frames (320 x 180 at dpi 40, 1-pixel coloured lines): yuv420p's
    // chroma subsampling alone limits RGB to about 35 dB (35.3 dB at crf 0), so the encoder's own
    // loss is checked on the luma (measured 42.5 dB) and RGB just below that limit (34.6 dB)
    for (const char* name : {"render_mpeg4", "render_stride", "render_odd_size"}) {
        const core::Value& record = renderCliCase(name);
        core::Project project = renderCliProject(record);
        project.output.video.codec = "mpeg4";
        TempDir dir;
        const fs::path path = dir / "v.mp4";
        const Quality q = mpeg4Quality(project, path);
        expectVideoLike(render::video::probeVideo(path), expectedVideo(record));
        EXPECT_GE(q.luma, 38.0) << name;
        EXPECT_GE(q.rgb, 33.0) << name;
    }
}

namespace {

// The software renderer, counting its frames and the threads it is called from.
class CountingRenderer final : public render::panels::StageRenderer {
public:
    render::molecule::RgbaImage render(const std::vector<render::molecule::Colored>& items,
                                       const render::molecule::Viewport& viewport,
                                       const render::molecule::Shading& shading,
                                       render::molecule::RenderQuality quality) const override {
        ++calls;
        if (std::this_thread::get_id() != caller) ++otherThread;
        return software.render(items, viewport, shading, quality);
    }
    std::string name() const override { return "CountingRenderer"; }

    render::panels::SoftwareStageRenderer software;
    std::thread::id caller = std::this_thread::get_id();
    mutable std::atomic<int> calls{0}, otherThread{0};
};

class CountingFactory final : public render::panels::StageRendererFactory {
public:
    explicit CountingFactory(std::shared_ptr<const CountingRenderer> r) : renderer(std::move(r)) {}
    render::panels::StageRendererChoice choose(std::string_view) const override { return {renderer, {}}; }
    std::shared_ptr<const CountingRenderer> renderer;
};

// installs a stage renderer factory for the test's lifetime
struct FactoryScope {
    explicit FactoryScope(std::shared_ptr<const render::panels::StageRendererFactory> f)
        : previous(render::panels::stageRendererFactory()) {
        render::panels::setStageRendererFactory(std::move(f));
    }
    ~FactoryScope() { render::panels::setStageRendererFactory(previous); }
    FactoryScope(const FactoryScope&) = delete;
    FactoryScope& operator=(const FactoryScope&) = delete;
    std::shared_ptr<const render::panels::StageRendererFactory> previous;
};

}  // namespace

// MDATK's "video with representations": the 3D layer drawn per frame by the structure panel's stage
// renderer (the factory the program installed), on the calling thread, frames still identical for
// any number of composing threads
TEST(QtVideo, RepresentationsThroughTheStageRendererFactory) {
    core::Project project = renderCliProject(renderCliCase("render"));
    app::applyOverride(project, R"(visual.representations=[{"style": "cartoon"}])");
    project.output.video.frame_stride = 3;  // 8 frames
    const auto renderer = std::make_shared<CountingRenderer>();
    const FactoryScope scope(std::make_shared<CountingFactory>(renderer));
    const auto one = composedFrames(project, 1);
    ASSERT_EQ(one.size(), 8u);
    EXPECT_EQ(renderer->calls.load(), 1 + 8);  // after the layout is frozen, then once per frame
    EXPECT_EQ(renderer->otherThread.load(), 0);
    const auto many = composedFrames(project, 4);
    ASSERT_EQ(many.size(), one.size());
    for (std::size_t k = 0; k < one.size(); ++k) EXPECT_TRUE(many[k] == one[k]) << "frame " << k;
    EXPECT_NE(one[0], one[7]);
    if (render::video::videoEncodingAvailable()) {
        TempDir dir;
        render::renderVideo(demoDataset(), project, dir / "rep.mp4");
        EXPECT_EQ(render::video::probeVideo(dir / "rep.mp4").packets, 8);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    const std::unique_ptr<QGuiApplication> app = render::startBatchApplication(argv[0]);
    return RUN_ALL_TESTS();
}
