// The GUI's Vulkan 3D view (src/render/vulkan_view.cpp) drawing a molecule with the GPU scene
// renderer on Qt's device: the window shows what the offscreen renderer draws (same pixels over
// the background), a fitted viewport keeps atoms round, and setting no molecule brings the
// placeholder back. Needs a display whose Qt platform supports Vulkan (X11 / Wayland; Xvfb with
// lavapipe works: xvfb-run -a ctest -R VulkanWindow); skipped elsewhere (e.g. the "offscreen"
// platform of CI's test job), also with MDBRIDGE_REQUIRE_VULKAN=1, which concerns the driver;
// MDBRIDGE_REQUIRE_VULKAN_WINDOW=1 (CI's Xvfb step) makes a missing Vulkan window a failure.

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <memory>

#include <QElapsedTimer>
#include <QEventLoop>
#include <QGuiApplication>
#include <QTimer>

#include "mdbridge/render/view.hpp"
#include "mdbridge/render/vulkan/scene_renderer.hpp"
#include "vulkan_harness.hpp"

namespace mdbridge::render {
namespace {

std::shared_ptr<MoleculeFrame> spheres() {
    auto frame = std::make_shared<MoleculeFrame>();
    molecule::Spheres s;
    s.xyz = {{0, 0, 0}, {6, 3, -2}, {-5, -4, 2}};
    s.radius = {4.0, 2.5, 3.0};
    s.owner = {0, 1, 2};
    s.atom = {-1, -1, -1};
    molecule::Colored c;
    c.prim = s;
    c.rgb = {{0.9, 0.2, 0.2}, {0.2, 0.7, 0.3}, {0.2, 0.3, 0.9}};
    frame->items = {c};
    frame->xlim = {-10, 10};
    frame->ylim = {-10, 10};
    frame->zrange = {-8, 8};
    frame->shading.fog_rgb = {0.102, 0.118, 0.153};
    frame->shading.fog_strength = 0.3;
    frame->shading.glow = std::vector<double>{0.0, 1.0, 0.0};
    return frame;
}

// waits until `view` has rendered `frames` more frames (false after 20 s)
bool waitFrames(View& view, int frames) {
    const int target = view.frames() + frames;
    QElapsedTimer clock;
    clock.start();
    while (view.frames() < target && clock.elapsed() < 20000) {
        QEventLoop loop;
        QTimer::singleShot(20, &loop, &QEventLoop::quit);
        loop.exec();
    }
    return view.frames() >= target;
}

TEST(VulkanWindow, DrawsTheMoleculeAsTheOffscreenRenderer) {
    std::unique_ptr<View> view = createView(Renderer::Vulkan).view;
    if (!view || view->renderer() != Renderer::Vulkan) {
        if (std::getenv("MDBRIDGE_REQUIRE_VULKAN_WINDOW")) {  // CI's Xvfb step
            FAIL() << "no Vulkan window on Qt platform " << QGuiApplication::platformName().toStdString();
        }
        GTEST_SKIP() << "no Vulkan window on Qt platform " << QGuiApplication::platformName().toStdString();
    }
    QWindow* window = view->window();
    window->resize(240, 180);
    window->show();
    ASSERT_TRUE(waitFrames(*view, 2)) << "the Vulkan window rendered no frame";
    const std::shared_ptr<MoleculeFrame> frame = spheres();
    ASSERT_TRUE(view->setMoleculeFrame(frame));
    ASSERT_TRUE(waitFrames(*view, 1));
    window->requestUpdate();  // a molecule is drawn on demand; one more frame after any in flight
    ASSERT_TRUE(waitFrames(*view, 1));
    const QImage shot = view->grabFrame().convertToFormat(QImage::Format_RGBA8888);
    ASSERT_FALSE(shot.isNull());
    const int w = shot.width(), h = shot.height();

    // the same scene offscreen on a device of our own, composed over the background
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    auto renderer = vulkan::SceneRenderer::create(*context, &problem);
    ASSERT_TRUE(renderer) << problem;
    const molecule::Viewport vp = fittedViewport(*frame, w, h);
    EXPECT_NEAR((vp.xlim[1] - vp.xlim[0]) / w, (vp.ylim[1] - vp.ylim[0]) / h, 1e-12);  // square pixels
    const molecule::RgbaImage img = renderer->render(frame->items, vp, frame->shading, frame->quality);
    ASSERT_EQ(img.width, static_cast<std::size_t>(w));
    ASSERT_EQ(img.height, static_cast<std::size_t>(h));
    int worst = 0;
    std::size_t covered = 0;
    for (int y = 0; y < h; ++y) {
        const uchar* row = shot.constScanLine(y);
        for (int x = 0; x < w; ++x) {
            const std::uint8_t* p = img.pixel(static_cast<std::size_t>(x), static_cast<std::size_t>(y));
            const double a = p[3] / 255.0;
            covered += p[3] > 0 ? 1 : 0;
            for (int k = 0; k < 3; ++k) {
                const double want = (p[k] / 255.0 * a + frame->background[static_cast<std::size_t>(k)] * (1 - a)) * 255.0;
                worst = std::max(worst, static_cast<int>(std::lround(std::fabs(want - row[4 * x + k]))));
            }
        }
    }
    EXPECT_GT(covered, static_cast<std::size_t>(w * h / 10));
    EXPECT_LE(worst, 2);  // 8-bit rounding of the composition

    // fast quality, then the placeholder again
    auto fast = std::make_shared<MoleculeFrame>(*frame);
    fast->quality = molecule::RenderQuality::Fast;
    view->setMoleculeFrame(fast);
    EXPECT_TRUE(waitFrames(*view, 1));
    view->setMoleculeFrame(nullptr);
    EXPECT_TRUE(waitFrames(*view, 2));  // the placeholder animates again
}

TEST(VulkanWindow, FittedViewportKeepsPixelsSquare) {
    MoleculeFrame f;
    f.xlim = {-10, 10};
    f.ylim = {0, 10};
    for (const auto& [w, h] : {std::pair{100, 100}, std::pair{300, 100}, std::pair{100, 300}}) {
        const molecule::Viewport vp = fittedViewport(f, w, h);
        EXPECT_NEAR((vp.xlim[1] - vp.xlim[0]) / w, (vp.ylim[1] - vp.ylim[0]) / h, 1e-12);
        EXPECT_LE(vp.xlim[0], -10.0);
        EXPECT_GE(vp.xlim[1], 10.0);
        EXPECT_LE(vp.ylim[0], 0.0);
        EXPECT_GE(vp.ylim[1], 10.0);
        EXPECT_DOUBLE_EQ(vp.xlim[0] + vp.xlim[1], 0.0);
        EXPECT_DOUBLE_EQ(vp.ylim[0] + vp.ylim[1], 10.0);
    }
}

}  // namespace
}  // namespace mdbridge::render

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    if (!std::getenv("DISPLAY") && !std::getenv("WAYLAND_DISPLAY") && !std::getenv("QT_QPA_PLATFORM")) {
        qputenv("QT_QPA_PLATFORM", "offscreen");  // no display: the window test skips
    }
    QGuiApplication app(argc, argv);
    return RUN_ALL_TESTS();
}
