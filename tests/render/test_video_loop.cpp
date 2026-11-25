// The video frame loop (src/render/video/frame_loop.*; MDATK's render_video) without Qt: frames
// composed by a simple test compositor. Positions and stride, frames in order and byte-identical
// for any number of threads, progress, cancellation and errors as MDATK handles them, and whole
// runs of tests/data/render/cli's video cases written through the encoder (frame count, size, frame
// rate, faststart, logs; tests/parity/make_render_cli_cases.py). The QPainter compositor and
// decoded frames: test_video_qt.cpp. The encoding tests skip in builds without FFmpeg.
#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/paths.hpp"
#include "video/encoder.hpp"
#include "video/frame_loop.hpp"
#include "video_cases.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;
using namespace mdbridge::test;
using render::DisplayList;
using render::video::FrameCompositor;

namespace {

std::uint8_t byte(double c) { return static_cast<std::uint8_t>(std::lround(std::clamp(c, 0.0, 1.0) * 255.0)); }

// A compositor without Qt: rectangles of the background filled, every vertex / marker centre / text
// origin of a frame's layer set to its colour. Content-dependent and deterministic; composing
// takes a varying time so that frames finish out of order on several threads.
class TestCompositor final : public FrameCompositor {
public:
    void setBackground(const DisplayList& bg) override {
        w = bg.width;
        h = bg.height;
        background.assign(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 4, 0);
        for (const render::DisplayItem& item : bg.items) {
            if (const auto* r = std::get_if<render::FillRectItem>(&item.item)) {
                for (int y = std::max(0, static_cast<int>(r->rect.y0)); y < std::min(h, static_cast<int>(r->rect.y1)); ++y) {
                    for (int x = std::max(0, static_cast<int>(r->rect.x0)); x < std::min(w, static_cast<int>(r->rect.x1)); ++x) {
                        set(background, x, y, r->color);
                    }
                }
            } else {
                for (const render::Point& p : points(item)) set(background, p, color(item));
            }
        }
    }
    std::vector<std::uint8_t> compose(const DisplayList& layer) const override {
        const int n = calls.fetch_add(1);
        if (failAt >= 0 && n == failAt) throw std::runtime_error("compositor failed");
        std::this_thread::sleep_for(std::chrono::microseconds((n * 7919) % 3000));
        std::vector<std::uint8_t> out = background;
        for (const render::DisplayItem& item : layer.items) {
            for (const render::Point& p : points(item)) set(out, p, color(item));
        }
        return out;
    }

    int w = 0, h = 0;
    std::vector<std::uint8_t> background;
    mutable std::atomic<int> calls{0};
    int failAt = -1;

private:
    static std::vector<render::Point> points(const render::DisplayItem& item) {
        std::vector<render::Point> out;
        if (const auto* m = std::get_if<render::MarkerItem>(&item.item)) out = m->centers;
        if (const auto* p = std::get_if<render::PathItem>(&item.item)) {
            for (const auto& sub : p->subpaths) out.insert(out.end(), sub.begin(), sub.end());
        }
        if (const auto* t = std::get_if<render::TextItem>(&item.item)) out.push_back(t->origin);
        if (const auto* i = std::get_if<render::ImageItem>(&item.item)) out.push_back(i->corner0);
        return out;
    }
    static render::Rgba color(const render::DisplayItem& item) {
        if (const auto* m = std::get_if<render::MarkerItem>(&item.item)) {
            if (!m->faces.empty()) return m->faces[0];
            if (!m->edges.empty()) return m->edges[0];
        }
        if (const auto* p = std::get_if<render::PathItem>(&item.item)) {
            if (p->stroke) return p->stroke->color;
            if (p->fill) return *p->fill;
        }
        if (const auto* t = std::get_if<render::TextItem>(&item.item)) return t->color;
        return {0.5, 0.25, 0.75, 1.0};
    }
    void set(std::vector<std::uint8_t>& px, render::Point p, render::Rgba c) const {
        if (!std::isfinite(p.x) || !std::isfinite(p.y)) return;
        set(px, static_cast<int>(std::floor(p.x)), static_cast<int>(std::floor(p.y)), c);
    }
    void set(std::vector<std::uint8_t>& px, int x, int y, render::Rgba c) const {
        if (x < 0 || y < 0 || x >= w || y >= h) return;
        std::uint8_t* q = &px[(static_cast<std::size_t>(y) * static_cast<std::size_t>(w) + static_cast<std::size_t>(x)) * 4];
        q[0] = byte(c.r);
        q[1] = byte(c.g);
        q[2] = byte(c.b);
        q[3] = 255;
    }
};

render::panels::BuiltFigure demoFigure(const core::Project& project) {
    return render::panels::buildFigure(demoDataset(), project, figureFonts());
}

struct LoopRun {
    std::vector<std::size_t> order;
    std::vector<std::vector<std::uint8_t>> frames;
    std::vector<std::pair<double, std::string>> progress;
};

LoopRun runLoop(const core::Project& project, std::span<const std::size_t> positions, unsigned threads,
            TestCompositor* compositor = nullptr) {
    render::panels::BuiltFigure fig = demoFigure(project);
    TestCompositor own;
    TestCompositor& comp = compositor ? *compositor : own;
    LoopRun run;
    const core::Reporter reporter({}, [&](double f, const std::string& m) { run.progress.emplace_back(f, m); });
    render::video::runFrameLoop(
        fig, figureFonts(), positions, comp,
        [&](std::size_t k, std::vector<std::uint8_t>& rgba) {
            run.order.push_back(k);
            run.frames.push_back(std::move(rgba));
        },
        reporter, {threads});
    return run;
}

}  // namespace

TEST(VideoPositions, StrideAsMdatk) {
    using render::video::videoPositions;
    EXPECT_EQ(videoPositions(24, 1).size(), 24u);
    EXPECT_EQ(videoPositions(24, 5), (std::vector<std::size_t>{0, 5, 10, 15, 20}));  // RN-17: 23 left out
    EXPECT_EQ(videoPositions(24, 23), (std::vector<std::size_t>{0, 23}));
    EXPECT_EQ(videoPositions(24, 100), (std::vector<std::size_t>{0}));
    EXPECT_EQ(videoPositions(3, 0), (std::vector<std::size_t>{0, 1, 2}));  // max(1, stride)
    EXPECT_EQ(videoPositions(3, -4), (std::vector<std::size_t>{0, 1, 2}));
    EXPECT_TRUE(videoPositions(0, 1).empty());
}

TEST(FrameLoop, ThreadCount) {
    EXPECT_EQ(render::video::frameLoopThreads({3}), 3u);
    const unsigned automatic = render::video::frameLoopThreads({});
    EXPECT_GE(automatic, 1u);
    EXPECT_LE(automatic, 8u);
}

TEST(FrameLoop, FramesInOrderAndIdenticalForAnyThreadCount) {
    const core::Project project = renderCliProject(renderCliCase("render"));
    const std::vector<std::size_t> positions = render::video::videoPositions(demoDataset().frameCount(), 1);
    const LoopRun one = runLoop(project, positions, 1);
    ASSERT_EQ(one.frames.size(), 24u);
    for (std::size_t k = 0; k < one.order.size(); ++k) EXPECT_EQ(one.order[k], k);
    EXPECT_EQ(one.frames[0].size(), 320u * 180u * 4u);
    EXPECT_NE(one.frames[0], one.frames[12]);  // the frames move
    for (const unsigned threads : {2u, 5u, 8u}) {
        const LoopRun many = runLoop(project, positions, threads);
        ASSERT_EQ(many.frames.size(), one.frames.size()) << threads;
        for (std::size_t k = 0; k < many.order.size(); ++k) EXPECT_EQ(many.order[k], k);
        for (std::size_t k = 0; k < one.frames.size(); ++k) {
            EXPECT_TRUE(many.frames[k] == one.frames[k]) << "frame " << k << " with " << threads << " threads";
        }
        EXPECT_EQ(many.progress, one.progress);
    }
}

TEST(FrameLoop, ProgressEveryNOver200Frames) {
    const core::Project project = renderCliProject(renderCliCase("render"));
    // 24 frames: every frame reported, (k + 1) / n
    const std::vector<std::size_t> positions = render::video::videoPositions(24, 1);
    const LoopRun run = runLoop(project, positions, 3);
    ASSERT_EQ(run.progress.size(), 24u);
    EXPECT_EQ(run.progress.front().second, "Video frame 1/24");
    EXPECT_DOUBLE_EQ(run.progress.front().first, 1.0 / 24.0);
    EXPECT_EQ(run.progress.back().second, "Video frame 24/24");
    EXPECT_DOUBLE_EQ(run.progress.back().first, 1.0);
    // 450 frames: every 2nd (450 // 200) and the last
    std::vector<std::size_t> many(450);
    for (std::size_t k = 0; k < many.size(); ++k) many[k] = (k * 5) % 24;
    const LoopRun big = runLoop(project, many, 4);
    ASSERT_EQ(big.frames.size(), 450u);
    ASSERT_EQ(big.progress.size(), 226u);
    EXPECT_EQ(big.progress[1].second, "Video frame 3/450");
    EXPECT_EQ(big.progress[224].second, "Video frame 449/450");
    EXPECT_EQ(big.progress[225].second, "Video frame 450/450");
    // a position shown twice is the same frame
    EXPECT_TRUE(big.frames[0] == big.frames[24]);
}

TEST(FrameLoop, CancelStopsBeforeTheNextFrame) {
    const core::Project project = renderCliProject(renderCliCase("render"));
    const std::vector<std::size_t> positions = render::video::videoPositions(24, 1);
    for (const unsigned threads : {1u, 4u}) {
        render::panels::BuiltFigure fig = demoFigure(project);
        TestCompositor comp;
        const core::Reporter reporter;
        std::vector<std::size_t> written;
        EXPECT_THROW(render::video::runFrameLoop(
                         fig, figureFonts(), positions, comp,
                         [&](std::size_t k, std::vector<std::uint8_t>&) {
                             written.push_back(k);
                             if (k == 5) reporter.cancel();
                         },
                         reporter, {threads}),
                     core::Cancelled);
        EXPECT_EQ(written, (std::vector<std::size_t>{0, 1, 2, 3, 4, 5})) << threads;
    }
}

TEST(FrameLoop, ErrorsPropagateInOrder) {
    const core::Project project = renderCliProject(renderCliCase("render"));
    const std::vector<std::size_t> positions = render::video::videoPositions(24, 1);
    for (const unsigned threads : {1u, 4u}) {
        TestCompositor comp;
        comp.failAt = 3;  // the 4th compose call (frame 3 with one thread)
        std::vector<std::size_t> written;
        render::panels::BuiltFigure fig = demoFigure(project);
        try {
            render::video::runFrameLoop(
                fig, figureFonts(), positions, comp,
                [&](std::size_t k, std::vector<std::uint8_t>&) { written.push_back(k); }, core::Reporter(), {threads});
            ADD_FAILURE() << "no error";
        } catch (const std::runtime_error& e) {
            EXPECT_STREQ(e.what(), "compositor failed");
        }
        EXPECT_LE(written.size(), 3u + threads) << threads;
        for (std::size_t k = 0; k < written.size(); ++k) EXPECT_EQ(written[k], k);
        // a failing sink stops the loop too
        render::panels::BuiltFigure fig2 = demoFigure(project);
        TestCompositor ok;
        EXPECT_THROW(render::video::runFrameLoop(
                         fig2, figureFonts(), positions, ok,
                         [&](std::size_t k, std::vector<std::uint8_t>&) {
                             if (k == 2) throw render::video::VideoError("sink failed");
                         },
                         core::Reporter(), {threads}),
                     render::video::VideoError);
    }
}

// ---- whole runs through the encoder ---------------------------------------------------------------

#define REQUIRE_FFMPEG()                                                               \
    if (!render::video::videoEncodingAvailable()) {                                    \
        GTEST_SKIP() << "built without the FFmpeg libraries (packaging/ffmpeg)";       \
    }

namespace {

bool leftovers(const fs::path& folder) {
    if (!fs::exists(folder)) return false;
    for (const auto& e : fs::directory_iterator(folder)) {
        if (e.path().filename().string().rfind(".MDATK_tmp_", 0) == 0) return true;
    }
    return false;
}

}  // namespace

class VideoRun : public ::testing::TestWithParam<const char*> {};

TEST_P(VideoRun, AsMdatksVideo) {
    REQUIRE_FFMPEG();
    const core::Value& record = renderCliCase(GetParam());
    const core::Project project = renderCliProject(record);
    const ExpectedVideo want = expectedVideo(record);
    TempDir dir;
    const fs::path path = dir.path() / stamped(want.file);  // results/ does not exist yet: created
    TestCompositor comp;
    LogCapture logs;
    std::vector<std::string> logged;
    const core::Reporter reporter([&](const std::string& m) { logged.push_back(m); });
    EXPECT_EQ(render::video::renderVideo(demoDataset(), project, path, figureFonts(), comp, reporter), path);
    ASSERT_TRUE(fs::exists(path));
    EXPECT_FALSE(leftovers(path.parent_path()));
    expectVideoLike(render::video::probeVideo(path), want);
    // the log lines as MDATK's, file names as given
    std::vector<std::string> expected = caseLines(record, "Rendering ", "{work}");
    for (const std::string& l : caseLines(record, "Wrote {work}", "{work}")) expected.push_back(l);
    const std::string mdatkPath = "{work}/" + want.file;
    for (std::string& l : expected) {
        const std::size_t at = l.find(mdatkPath);
        ASSERT_NE(at, std::string::npos) << l;
        l.replace(at, mdatkPath.size(), core::pathToUtf8(path));
    }
    EXPECT_EQ(logged, expected);
}

INSTANTIATE_TEST_SUITE_P(RenderCliCases, VideoRun,
                         ::testing::Values("render", "render_stride", "render_mpeg4", "render_odd_size",
                                           "render_codec_fallback"));

TEST(VideoRender, CodecFallbackWarns) {
    REQUIRE_FFMPEG();
    const core::Project project = renderCliProject(renderCliCase("render_codec_fallback"));
    TempDir dir;
    TestCompositor comp;
    LogCapture logs;
    render::video::renderVideo(demoDataset(), project, dir / "v.mp4", figureFonts(), comp, core::Reporter());
    std::vector<std::string> warnings;
    for (const std::string& l : logs.lines()) {
        if (l.rfind("WARNING: ", 0) == 0) warnings.push_back(l);
    }
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_EQ(warnings[0].rfind("WARNING: video: encoder nosuchcodec is not available; using ", 0), 0u) << warnings[0];
}

TEST(VideoRender, CancelLeavesNoFile) {
    REQUIRE_FFMPEG();
    const core::Project project = renderCliProject(renderCliCase("render"));
    for (const unsigned threads : {1u, 3u}) {
        TempDir dir;
        const fs::path path = dir / "MDATK_video.mp4";
        TestCompositor comp;
        const auto flag = std::make_shared<std::atomic<bool>>(false);
        const core::Reporter reporter({}, [flag](double fraction, const std::string&) {
            if (fraction > 0.2) flag->store(true);
        }, flag);
        render::video::VideoRenderOptions options;
        options.loop.threads = threads;
        EXPECT_THROW(render::video::renderVideo(demoDataset(), project, path, figureFonts(), comp, reporter, options),
                     core::Cancelled);
        EXPECT_FALSE(fs::exists(path));
        EXPECT_FALSE(leftovers(dir.path()));
    }
}

TEST(VideoRender, ExistingFileKeptOnFailure) {
    REQUIRE_FFMPEG();
    const core::Project project = renderCliProject(renderCliCase("render"));
    TempDir dir;
    const fs::path path = dir / "MDATK_video.mp4";
    writeFile(path, "old");
    TestCompositor comp;
    comp.failAt = 4;
    EXPECT_THROW(render::video::renderVideo(demoDataset(), project, path, figureFonts(), comp, core::Reporter()),
                 std::runtime_error);
    EXPECT_EQ(readFile(path), "old");
    EXPECT_FALSE(leftovers(dir.path()));
}

TEST(VideoRender, NoFramesIsAnError) {  // RN-02
    const core::Project project = renderCliProject(renderCliCase("render"));
    const core::Dataset empty({}, {}, {}, "ns");
    TempDir dir;
    TestCompositor comp;
    try {
        render::video::renderVideo(empty, project, dir / "v.mp4", figureFonts(), comp, core::Reporter());
        ADD_FAILURE() << "no error";
    } catch (const render::video::VideoError& e) {
        EXPECT_STREQ(e.what(), "no analysed frames to render");
        EXPECT_EQ(e.pythonName(), "VideoError");
    }
    EXPECT_TRUE(fs::is_empty(dir.path()));
}
