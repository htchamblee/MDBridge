#pragma once
// The GPU renderer of the molecule scene (Vulkan 1.3): the same picture as the built-in software
// renderer (render/molecule/software.hpp) for the same items, viewport and shading — "same look,
// not same pixels" (docs/specs/m4-plan.md section 3: coverage IoU >= 0.99, SSIM >= 0.97).
//
//   std::string problem;
//   auto context = vulkan::Context::create({}, &problem);
//   auto renderer = vulkan::SceneRenderer::create(*context, &problem);
//   molecule::RgbaImage image = renderer->render(items, viewport, shading);   // offscreen, read back
//
// How it draws (src/render/shaders/molecule_*):
//   * geometry passes into a G-buffer (front-most fragment per pixel: lit colour, view depth, owner
//     residue, opacity; depth test on the software renderer's depth key, ties to the later
//     primitive): spheres and sticks as impostors on screen-aligned quads (ray / sphere and
//     ray / cylinder per pixel centre, the spheres' quarter-pixel phase), mesh triangles on their
//     pixel boxes with the software renderer's widened edge test; the primitives are set up on the
//     CPU exactly as the software renderer sets them up, and the shaders repeat its per-pixel
//     arithmetic in float32. Opaque items and translucent items (opacity < 0.999) get a G-buffer
//     each; the front-most translucent layer in front of the opaque surface is blended once.
//   * compute passes for the compositing: depth cue (fog), glow tint and the glow halo (the
//     glowing pixels reduced to a quarter-resolution grid, scipy's Gaussian, the disc
//     normalisation, bilinear upsampling), ink outlines at depth jumps and silhouettes, alpha
//     anti-aliasing of the edges, the translucent layer.
//   * "fast" quality: half resolution without outlines and anti-aliasing, scaled back (as the
//     software renderer, RN-27 included).
// Offscreen (render, resolve) for batch work and tests; recordFrame / recordPresent draw into a
// window's command buffer (the GUI's QVulkanWindow, src/render/vulkan_view.cpp) without a read-back.
//
// Thread-safe: render() and resolve() may be called from several threads (the GPU work is
// serialised; the primitives are set up on the calling thread). Errors: std::invalid_argument for
// an unusable viewport or inconsistent items (as the software renderer, RN-28), std::runtime_error
// when the GPU fails or the scene does not fit the device's limits (callers fall back to the
// software renderer).

#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/software.hpp"
#include "mdbridge/render/vulkan/context.hpp"

namespace mdbridge::render::vulkan {

struct SceneRendererOptions {
    // upload buffers for frames in flight (a window's concurrent frame count); 1 for offscreen work
    std::uint32_t frameSlots = 1;
};

class SceneRenderer {
public:
    // nullptr with `problem` set when a pipeline cannot be created on the device
    static std::unique_ptr<SceneRenderer> create(Context& context, std::string* problem,
                                                 const SceneRendererOptions& options = {});
    ~SceneRenderer();
    SceneRenderer(const SceneRenderer&) = delete;
    SceneRenderer& operator=(const SceneRenderer&) = delete;

    // The items as an image of the viewport's size (straight alpha, row 0 at the top), as
    // molecule::SoftwareRenderer::render draws them. Waits for the GPU.
    molecule::RgbaImage render(const std::vector<molecule::Colored>& items, const molecule::Viewport& viewport,
                               const molecule::Shading& shading,
                               molecule::RenderQuality quality = molecule::RenderQuality::High);

    // The G-buffer of all items (whatever their opacity) at full resolution, as
    // molecule::SoftwareRenderer::resolve returns it (tests).
    molecule::ResolvedPixels resolve(const std::vector<molecule::Colored>& items, const molecule::Viewport& viewport);

    // ---- drawing into a window (the caller owns the command buffer and submits it) ----

    // Sets the frame up in upload slot `slot` (< frameSlots, a slot whose previous frame the GPU
    // has finished) and records the geometry and compositing passes into `cmd` (outside a render
    // pass). The image is then in the renderer's output buffer for recordPresent(). Throws as
    // render(). Waits for the device to be idle when the image size changes.
    void recordFrame(VkCommandBuffer cmd, std::uint32_t slot, const std::vector<molecule::Colored>& items,
                     const molecule::Viewport& viewport, const molecule::Shading& shading,
                     molecule::RenderQuality quality);
    // The pipeline recordPresent() draws with, for `renderPass` (subpass 0, `samples` samples).
    // False with `problem` set on failure.
    bool createPresentPipeline(VkRenderPass renderPass, VkSampleCountFlagBits samples, std::string* problem);
    // Inside the render pass: the last recorded frame over `background` (RGB 0..1, the colours
    // as stored: sRGB-encoded like everything QPainter draws), filling a target of the given size
    // from its top-left corner (fast quality: 2 x 2 target pixels per image pixel).
    void recordPresent(VkCommandBuffer cmd, std::uint32_t targetWidth, std::uint32_t targetHeight,
                       const std::array<float, 3>& background) const;

    Context& context() const { return m_context; }

    struct Impl;

private:
    explicit SceneRenderer(Context& context);

    Context& m_context;
    std::unique_ptr<Impl> m_impl;
    std::mutex m_mutex;
};

}  // namespace mdbridge::render::vulkan
