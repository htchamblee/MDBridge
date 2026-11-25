// The Vulkan scene renderer (src/render/vulkan/scene_renderer.cpp) against the software renderer
// (src/render/molecule/software.cpp), the reference of its look: on every molecule case of
// data/render/molecule (the renders of make_molecule_cases.py and every scene, also larger, glowing
// and fast) the coverage (alpha > 0) IoU is >= 0.99 and the SSIM >= 0.97 ("same look, not same
// pixels", docs/specs/m4-plan.md section 3); the resolved G-buffer against SoftwareRenderer::resolve;
// errors as the software renderer's; one renderer from several threads; the GPU stage renderer and
// its renderer choice (RN-16). Runs on the machine's Vulkan 1.3 driver (lavapipe in CI); skipped
// without one unless MDBRIDGE_REQUIRE_VULKAN=1 (vulkan_harness.hpp).
//
// MDBRIDGE_VULKAN_BENCHMARK=1 also prints frames/s of both renderers for the demo at several sizes
// (VulkanScene.Benchmark).

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <string>
#include <thread>
#include <vector>

#include "mdbridge/render/molecule/software.hpp"
#include "mdbridge/render/orientation.hpp"
#include "mdbridge/render/panels/vulkan_stage.hpp"
#include "mdbridge/render/vulkan/scene_renderer.hpp"
#include "molecule_cases.hpp"
#include "vulkan_harness.hpp"

namespace mdbridge::render::molecule {
namespace {

using namespace molecule_cases;  // NOLINT(google-build-using-namespace): the case reader

std::unique_ptr<vulkan::SceneRenderer>& rendererSlot() {
    static std::unique_ptr<vulkan::SceneRenderer> renderer;
    return renderer;
}

vulkan::SceneRenderer& sharedRenderer(vulkan::Context& context) {
    std::unique_ptr<vulkan::SceneRenderer>& renderer = rendererSlot();
    if (!renderer) {
        std::string problem;
        renderer = vulkan::SceneRenderer::create(context, &problem);
        if (!renderer) {
            throw std::runtime_error("SceneRenderer::create: " + problem);
        }
    }
    return *renderer;
}

// The shared renderer goes before the shared context (registered after the harness's environment,
// so torn down before it).
class RendererEnvironment : public ::testing::Environment {
public:
    void TearDown() override { rendererSlot().reset(); }
};
const bool kRendererEnvironmentRegistered = (::testing::AddGlobalTestEnvironment(new RendererEnvironment), true);

// ---- image comparison -------------------------------------------------------------------------

struct Similarity {
    double iou = 1.0;   // coverage (alpha > 0)
    double ssim = 1.0;  // mean over the premultiplied R, G, B, A channels
    double ssimMin = 1.0;
    double meanAbs = 0.0;  // mean absolute channel difference (premultiplied, 0..255)
    std::size_t covered = 0;  // pixels with alpha > 0 in the reference
};

// SSIM (Wang et al. 2004): Gaussian window 11 x 11, sigma 1.5, K1 0.01, K2 0.03, L 255, over the
// positions where the window lies inside the image (the whole image as one window when smaller).
double ssimChannel(const std::vector<double>& a, const std::vector<double>& b, std::size_t w, std::size_t h) {
    const int r = 5;
    std::array<double, 11> g{};
    double gs = 0;
    for (int i = -r; i <= r; ++i) {
        g[static_cast<std::size_t>(i + r)] = std::exp(-0.5 * i * i / (1.5 * 1.5));
        gs += g[static_cast<std::size_t>(i + r)];
    }
    for (double& v : g) {
        v /= gs;
    }
    const double C1 = std::pow(0.01 * 255, 2), C2 = std::pow(0.03 * 255, 2);
    auto window = [&](std::size_t x0, std::size_t y0, bool gauss, std::size_t ww, std::size_t wh) {
        double sw = 0, ma = 0, mb = 0, saa = 0, sbb = 0, sab = 0;
        for (std::size_t y = 0; y < wh; ++y) {
            for (std::size_t x = 0; x < ww; ++x) {
                const double wt = gauss ? g[y] * g[x] : 1.0;
                const double va = a[(y0 + y) * w + x0 + x], vb = b[(y0 + y) * w + x0 + x];
                sw += wt;
                ma += wt * va;
                mb += wt * vb;
                saa += wt * va * va;
                sbb += wt * vb * vb;
                sab += wt * va * vb;
            }
        }
        ma /= sw;
        mb /= sw;
        const double va = saa / sw - ma * ma, vb = sbb / sw - mb * mb, cab = sab / sw - ma * mb;
        return ((2 * ma * mb + C1) * (2 * cab + C2)) / ((ma * ma + mb * mb + C1) * (va + vb + C2));
    };
    if (w < 11 || h < 11) {
        return window(0, 0, false, w, h);
    }
    double sum = 0;
    std::size_t n = 0;
    for (std::size_t y = 0; y + 11 <= h; ++y) {
        for (std::size_t x = 0; x + 11 <= w; ++x) {
            sum += window(x, y, true, 11, 11);
            ++n;
        }
    }
    return sum / static_cast<double>(n);
}

Similarity similarity(const RgbaImage& got, const RgbaImage& want) {
    Similarity s;
    const std::size_t n = want.width * want.height;
    std::size_t inter = 0, uni = 0;
    std::array<std::vector<double>, 4> ga, wa;
    double absSum = 0;
    for (auto* v : {&ga, &wa}) {
        for (std::vector<double>& c : *v) {
            c.resize(n);
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint8_t* p = got.rgba.data() + 4 * i;
        const std::uint8_t* q = want.rgba.data() + 4 * i;
        const bool a = p[3] > 0, b = q[3] > 0;
        inter += (a && b) ? 1 : 0;
        uni += (a || b) ? 1 : 0;
        s.covered += b ? 1 : 0;
        for (std::size_t k = 0; k < 4; ++k) {
            ga[k][i] = k == 3 ? p[3] : p[k] * (p[3] / 255.0);
            wa[k][i] = k == 3 ? q[3] : q[k] * (q[3] / 255.0);
            absSum += std::fabs(ga[k][i] - wa[k][i]);
        }
    }
    s.iou = uni ? static_cast<double>(inter) / static_cast<double>(uni) : 1.0;
    s.ssim = 0;
    for (std::size_t k = 0; k < 4; ++k) {
        const double v = ssimChannel(ga[k], wa[k], want.width, want.height);
        s.ssim += v / 4;
        s.ssimMin = std::min(s.ssimMin, v);
    }
    s.meanAbs = n ? absSum / static_cast<double>(4 * n) : 0.0;
    return s;
}

void expectSameLook(const RgbaImage& gpu, const RgbaImage& cpu, const std::string& what) {
    ASSERT_EQ(gpu.width, cpu.width) << what;
    ASSERT_EQ(gpu.height, cpu.height) << what;
    const Similarity s = similarity(gpu, cpu);
    std::printf("%-34s %7zu px  IoU %.4f  SSIM %.4f (min %.4f)  mean |d| %.3f\n", what.c_str(), s.covered, s.iou,
                s.ssim, s.ssimMin, s.meanAbs);
    EXPECT_GT(s.covered, 0u) << what;
    EXPECT_GE(s.iou, 0.99) << what;
    EXPECT_GE(s.ssim, 0.97) << what;
}

// a viewport fitting the items (rotation `R`, square limits around their extent)
Viewport fitted(const std::vector<Colored>& items, const std::array<Vec3, 3>& R, int width, int height) {
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    auto add = [&](const Vec3& p, double rad) {
        for (int k = 0; k < 2; ++k) {
            const double v = p[0] * R[static_cast<std::size_t>(k)][0] + p[1] * R[static_cast<std::size_t>(k)][1] +
                             p[2] * R[static_cast<std::size_t>(k)][2];
            if (std::isfinite(v)) {
                lo[k] = std::min(lo[k], v - rad);
                hi[k] = std::max(hi[k], v + rad);
            }
        }
    };
    for (const Colored& it : items) {
        if (const auto* m = std::get_if<Mesh>(&it.prim)) {
            for (const Vec3& p : m->verts) add(p, 0);
        } else if (const auto* s = std::get_if<Spheres>(&it.prim)) {
            for (std::size_t i = 0; i < s->xyz.size(); ++i) add(s->xyz[i], s->radius[i]);
        } else if (const auto* k = std::get_if<Sticks>(&it.prim)) {
            for (std::size_t i = 0; i < k->a.size(); ++i) {
                add(k->a[i], k->radius[i]);
                add(k->b[i], k->radius[i]);
            }
        }
    }
    const double cx = 0.5 * (lo[0] + hi[0]), cy = 0.5 * (lo[1] + hi[1]);
    const double half = 0.55 * std::max(hi[0] - lo[0], hi[1] - lo[1]);
    const double aspect = static_cast<double>(width) / height;
    Viewport vp;
    vp.rotation = R;
    vp.xlim = {cx - half * std::max(aspect, 1.0), cx + half * std::max(aspect, 1.0)};
    vp.ylim = {cy - half / std::min(aspect, 1.0), cy + half / std::min(aspect, 1.0)};
    vp.width = width;
    vp.height = height;
    const auto [z0, z1] = depthRange(items, R[2]);
    vp.zrange = {z0, z1};
    return vp;
}

std::array<Vec3, 3> rotation(double ax, double ay, double az) {
    const Mat3 m = eulerToMatrix(Vec3d{ax, ay, az});
    return {Vec3{m[0][0], m[0][1], m[0][2]}, Vec3{m[1][0], m[1][1], m[1][2]}, Vec3{m[2][0], m[2][1], m[2][2]}};
}

std::size_t residues(const std::vector<Colored>& items) {
    std::int32_t n = 0;
    for (const Colored& it : items) {
        const std::vector<std::int32_t>* owner = std::visit([](const auto& p) { return &p.owner; }, it.prim);
        for (std::int32_t o : *owner) n = std::max(n, o + 1);
    }
    return static_cast<std::size_t>(n);
}

// ---- tests ------------------------------------------------------------------------------------

TEST(VulkanScene, MatchesTheSoftwareRendererOnEveryRenderCase) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    vulkan::SceneRenderer& gpu = sharedRenderer(*context);
    const SoftwareRenderer cpu;
    const std::vector<std::string> names = renderNames();
    ASSERT_GE(names.size(), 37u);
    for (const std::string& name : names) {
        SCOPED_TRACE(name);
        const RenderCase c = renderCase(name);
        expectSameLook(gpu.render(c.items, c.vp, c.shading, c.quality), cpu.render(c.items, c.vp, c.shading, c.quality),
                       name);
    }
}

TEST(VulkanScene, MatchesTheSoftwareRendererOnEveryScene) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    vulkan::SceneRenderer& gpu = sharedRenderer(*context);
    const SoftwareRenderer cpu;
    const Value& scenes = at(cases(), "scenes");
    std::size_t n = 0;
    for (const auto& [name, rec] : *scenes.get<core::Table>()) {
        const std::vector<Colored> items = storedItems("scenes.npz", name, at(rec, "items"));
        const std::size_t nres = residues(items);
        Shading plain{{0.988, 0.988, 0.984}, 0.45, std::nullopt, {1.0, 0.69, 0.0}};
        Shading glow = plain;
        glow.glow = std::vector<double>(nres, 0.0);
        for (std::size_t r = 0; r < nres; r += 3) {
            (*glow.glow)[r] = r % 2 ? 1.0 : 0.6;
        }
        struct Variant {
            const char* label;
            std::array<Vec3, 3> R;
            int width, height;
            const Shading* shading;
            RenderQuality quality;
        };
        const Variant variants[] = {
            {"front", rotation(-90, 0, 0), 256, 200, &plain, RenderQuality::High},
            {"side glow", rotation(10, 60, 30), 300, 300, &glow, RenderQuality::High},
            {"fast", rotation(-90, 0, 0), 257, 193, &plain, RenderQuality::Fast},
        };
        for (const Variant& v : variants) {
            const std::string what = name + " " + v.label;
            SCOPED_TRACE(what);
            const Viewport vp = fitted(items, v.R, v.width, v.height);
            expectSameLook(gpu.render(items, vp, *v.shading, v.quality), cpu.render(items, vp, *v.shading, v.quality),
                           what);
        }
        ++n;
    }
    EXPECT_GE(n, 17u);
}

TEST(VulkanScene, ResolvedFragmentsMatchTheSoftwareRenderer) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    vulkan::SceneRenderer& gpu = sharedRenderer(*context);
    const SoftwareRenderer cpu;
    for (const std::string& name : renderNames()) {
        SCOPED_TRACE(name);
        const RenderCase c = renderCase(name);
        const ResolvedPixels g = gpu.resolve(c.items, c.vp);
        const ResolvedPixels s = cpu.resolve(c.items, c.vp);
        ASSERT_EQ(g.width, s.width);
        ASSERT_EQ(g.height, s.height);
        EXPECT_DOUBLE_EQ(g.scale, s.scale);
        std::size_t si = 0, common = 0, owners = 0, depths = 0, colours = 0;
        for (std::size_t gi = 0; gi < g.pixel.size(); ++gi) {
            while (si < s.pixel.size() && s.pixel[si] < g.pixel[gi]) ++si;
            if (si == s.pixel.size() || s.pixel[si] != g.pixel[gi]) continue;
            ++common;
            owners += g.owner[gi] == s.owner[si] ? 1 : 0;
            depths += std::fabs(g.depth[gi] - s.depth[si]) < 1e-2 ? 1 : 0;
            bool near = true;
            for (std::size_t k = 0; k < 3; ++k) near = near && std::fabs(g.rgb[gi][k] - s.rgb[si][k]) < 2.0f / 255;
            colours += near ? 1 : 0;
            EXPECT_NEAR(g.alpha[gi], s.alpha[si], 1e-4);
        }
        const double total = static_cast<double>(std::max(g.pixel.size(), s.pixel.size()));
        const double frac = total > 0 ? static_cast<double>(common) / total : 1.0;
        std::printf("%-26s pixels %zu / %zu common %.5f owner %.5f depth %.5f colour %.5f\n", name.c_str(),
                    g.pixel.size(), s.pixel.size(), frac, common ? double(owners) / double(common) : 1.0,
                    common ? double(depths) / double(common) : 1.0, common ? double(colours) / double(common) : 1.0);
        EXPECT_GE(frac, 0.995);
        if (common) {
            EXPECT_GE(static_cast<double>(owners) / static_cast<double>(common), 0.995);
            EXPECT_GE(static_cast<double>(depths) / static_cast<double>(common), 0.995);
            EXPECT_GE(static_cast<double>(colours) / static_cast<double>(common), 0.99);
        }
    }
}

TEST(VulkanScene, UnusableViewsAndItemsThrowAsTheSoftwareRenderer) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    vulkan::SceneRenderer& gpu = sharedRenderer(*context);
    const RenderCase c = renderCase("sphere");
    Viewport bad = c.vp;
    bad.xlim = {1.0, 1.0};
    EXPECT_THROW(gpu.render(c.items, bad, c.shading), std::invalid_argument);
    bad = c.vp;
    bad.rotation[0][0] = std::nan("");
    EXPECT_THROW(gpu.render(c.items, bad, c.shading), std::invalid_argument);
    bad = c.vp;
    bad.width = 0;
    EXPECT_THROW(gpu.render(c.items, bad, c.shading), std::invalid_argument);
    std::vector<Colored> items = c.items;
    items[0].rgb.push_back({1, 1, 1});
    EXPECT_THROW(gpu.render(items, c.vp, c.shading), std::invalid_argument);
    // still usable afterwards
    EXPECT_EQ(gpu.render(c.items, c.vp, c.shading).width, static_cast<std::size_t>(c.vp.width));
}

TEST(VulkanScene, EmptySceneSizeChangesAndOnePixelViews) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    vulkan::SceneRenderer& gpu = sharedRenderer(*context);
    const SoftwareRenderer cpu;
    const RenderCase c = renderCase("forty_spheres");
    for (const auto& [w, h] : {std::pair{1, 1}, std::pair{3, 2}, std::pair{640, 7}, std::pair{100, 100}}) {
        Viewport vp = c.vp;
        vp.width = w;
        vp.height = h;
        for (RenderQuality q : {RenderQuality::High, RenderQuality::Fast}) {
            const RgbaImage g = gpu.render(c.items, vp, c.shading, q);
            const RgbaImage s = cpu.render(c.items, vp, c.shading, q);
            EXPECT_EQ(g.width, s.width);
            EXPECT_EQ(g.height, s.height);
            if (w >= 100) {
                expectSameLook(g, s, "forty_spheres " + std::to_string(w) + "x" + std::to_string(h));
            }
        }
        const RgbaImage empty = gpu.render({}, vp, c.shading);
        EXPECT_TRUE(std::all_of(empty.rgba.begin(), empty.rgba.end(), [](std::uint8_t v) { return v == 0; }));
    }
}

TEST(VulkanScene, OneRendererServesSeveralThreads) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    vulkan::SceneRenderer& gpu = sharedRenderer(*context);
    const RenderCase a = renderCase("demo_mixed");
    const RenderCase b = renderCase("glow_on");
    const RgbaImage ra = gpu.render(a.items, a.vp, a.shading);
    const RgbaImage rb = gpu.render(b.items, b.vp, b.shading);
    std::vector<std::future<bool>> runs;
    for (int t = 0; t < 4; ++t) {
        runs.push_back(std::async(std::launch::async, [&, t] {
            bool same = true;
            for (int i = 0; i < 3; ++i) {
                const RenderCase& c = (t + i) % 2 ? a : b;
                same = same && gpu.render(c.items, c.vp, c.shading).rgba == ((t + i) % 2 ? ra : rb).rgba;
            }
            return same;
        }));
    }
    for (auto& r : runs) {
        EXPECT_TRUE(r.get());
    }
}

// ---- the structure panel's GPU renderer choice -------------------------------------------------

TEST(VulkanStage, ChoosesTheGpuRendererWhenUsable) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    panels::GpuStageOptions options;
    options.context = vulkan::vulkan_test::testOptions();
    panels::GpuStageRenderers gpu(options);
    ASSERT_TRUE(gpu.usable()) << gpu.problem();
    EXPECT_FALSE(gpu.deviceName().empty());
    const auto factory = gpu.factory();
    for (const char* name : {"auto", "vtk", "vulkan"}) {
        const panels::StageRendererChoice c = factory->choose(name);
        const bool software = std::string(name) == "auto" && gpu.cpuDevice();  // lavapipe: software is faster
        EXPECT_EQ(c.renderer->name(), software ? "SoftwareRenderer" : "VulkanRenderer") << name;
        EXPECT_TRUE(c.note.empty()) << name;
    }
    EXPECT_EQ(factory->choose("software").renderer->name(), "SoftwareRenderer");
    const RenderCase c = renderCase("demo_cartoon_glow");
    const auto renderer = factory->choose("vulkan").renderer;
    EXPECT_TRUE(renderer->drawsGlow());
    expectSameLook(renderer->render(c.items, c.vp, c.shading, c.quality),
                   SoftwareRenderer().render(c.items, c.vp, c.shading, c.quality), "stage demo_cartoon_glow");
    Viewport reversed = c.vp;
    reversed.xlim = {1.0, 0.0};
    EXPECT_THROW(renderer->render(c.items, reversed, c.shading, c.quality), std::invalid_argument);
}

TEST(VulkanStage, FallsBackToTheSoftwareRendererWithTheNote) {
    panels::GpuStageOptions options;
    options.context.deviceName = "no device is called like this";
    panels::GpuStageRenderers gpu(options);
    EXPECT_FALSE(gpu.usable());
    EXPECT_FALSE(gpu.problem().empty());
    const auto factory = gpu.factory();
    const panels::StageRendererChoice autoChoice = factory->choose("auto");
    EXPECT_EQ(autoChoice.renderer->name(), "SoftwareRenderer");
    EXPECT_TRUE(autoChoice.note.empty());
    for (const char* name : {"vtk", "vulkan"}) {
        const panels::StageRendererChoice c = factory->choose(name);
        EXPECT_EQ(c.renderer->name(), "SoftwareRenderer") << name;
        EXPECT_EQ(c.note, panels::kNoGpuRendererNote) << name;
    }
}

TEST(VulkanStage, InstallsAndUninstallsTheFactory) {
    ASSERT_EQ(panels::stageRendererFactory(), nullptr);
    {
        panels::GpuStageOptions options;
        options.context = vulkan::vulkan_test::testOptions();
        auto gpu = panels::installGpuStageRenderers(options);
        EXPECT_EQ(panels::stageRendererFactory(), gpu->factory());
        // a renderer that outlives the installation draws with the software renderer
        const auto renderer = gpu->factory()->choose("vulkan").renderer;
        const RenderCase c = renderCase("sphere");
        const RgbaImage before = renderer->render(c.items, c.vp, c.shading, c.quality);
        EXPECT_EQ(before.width, static_cast<std::size_t>(c.vp.width));
        std::shared_ptr<const panels::StageRenderer> kept = renderer;
        gpu.reset();
        EXPECT_EQ(panels::stageRendererFactory(), nullptr);
        EXPECT_EQ(kept->render(c.items, c.vp, c.shading, c.quality).rgba,
                  SoftwareRenderer().render(c.items, c.vp, c.shading, c.quality).rgba);
    }
}

// ---- frames per second (MDBRIDGE_VULKAN_BENCHMARK=1) --------------------------------------------

TEST(VulkanScene, Benchmark) {
    const char* flag = std::getenv("MDBRIDGE_VULKAN_BENCHMARK");
    if (!flag || std::string(flag) != "1") {
        GTEST_SKIP() << "set MDBRIDGE_VULKAN_BENCHMARK=1 to measure frames/s";
    }
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    vulkan::SceneRenderer& gpu = sharedRenderer(*context);
    const SoftwareRenderer cpu;
    std::printf("device: %s; software renderer: %u threads\n", context->deviceInfo().name.c_str(),
                std::max(1u, std::min(16u, std::thread::hardware_concurrency())));
    for (const char* scene : {"demo_cartoon", "demo_surf", "demo_mixed", "demo_vdw"}) {
        const std::vector<Colored> items = storedItems("scenes.npz", scene, at(at(at(cases(), "scenes"), scene), "items"));
        Shading glow{{0.988, 0.988, 0.984}, 0.45, std::vector<double>(residues(items), 0.0), {1.0, 0.69, 0.0}};
        for (std::size_t r = 0; r < glow.glow->size(); r += 4) (*glow.glow)[r] = 1.0;
        for (const auto& [w, h] : {std::pair{320, 180}, std::pair{640, 360}, std::pair{1280, 720}, std::pair{1920, 1080}}) {
            const Viewport vp = fitted(items, rotation(-90, 0, 0), w, h);
            auto fps = [&](auto&& draw) {
                draw();  // warm-up (pipelines, buffers of this size)
                const auto t0 = std::chrono::steady_clock::now();
                int n = 0;
                double s = 0;
                do {
                    draw();
                    ++n;
                    s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                } while (s < 1.5 && n < 200);
                return n / s;
            };
            const double g = fps([&] { (void)gpu.render(items, vp, glow); });
            const double c = fps([&] { (void)cpu.render(items, vp, glow); });
            std::printf("%-14s %4d x %-4d  Vulkan %7.1f frames/s  software %7.1f frames/s  (x %.2f)\n", scene, w, h, g, c,
                        g / c);
        }
    }
}

}  // namespace
}  // namespace mdbridge::render::molecule
