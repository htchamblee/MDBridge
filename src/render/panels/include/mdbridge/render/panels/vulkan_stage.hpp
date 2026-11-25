#pragma once
// The structure view's 3D layer on the GPU: a StageRenderer drawing with the Vulkan scene renderer
// (render/vulkan/scene_renderer.hpp) and the renderer choice that structure panels use once it is
// installed (setStageRendererFactory). No Qt.
//
//   int main() {
//       const auto gpu = mdbridge::render::panels::installGpuStageRenderers();  // until main returns
//       ...
//   }
//
// visual.renderer (MDATK's get_renderer, RN-16):
//   "software"         the built-in software renderer
//   "auto"             the Vulkan renderer when this machine can run it on a GPU, else the software
//                      renderer (silently: both draw the same picture). A Vulkan driver that runs on
//                      the CPU (Mesa's lavapipe) is left to the software renderer, which is about
//                      twice as fast there (VulkanScene.Benchmark)
//   "vtk", "vulkan"    the Vulkan renderer when it can run, else the software renderer with the note
//                      kNoGpuRendererNote (logged once per panel as a warning)
//   anything else      the software renderer (validate() reports unknown names)
// The Vulkan device is opened on the first GPU choice, once per process. A frame the GPU fails to
// draw (device lost, out of memory, a scene beyond the device's limits) is drawn by the software
// renderer instead, with one warning; later frames then use the software renderer too, so a video
// does not change its look between frames.

#include <memory>
#include <string>

#include "mdbridge/render/panels/structure.hpp"
#include "mdbridge/render/vulkan/context.hpp"

namespace mdbridge::render::panels {

struct GpuStageOptions {
    vulkan::ContextOptions context;  // device choice (deviceName), validation, pipeline cache
};

class GpuStageBackend;  // the shared device and renderer (vulkan_stage.cpp)

// The installed GPU renderer choice. Destroying it uninstalls the factory and releases the Vulkan
// device (renderers that panels still hold then draw with the software renderer): keep it alive
// for as long as figures are drawn, and destroy it before static destruction (some drivers, e.g.
// lavapipe, hang when a device is destroyed after main() has returned).
class GpuStageRenderers {
public:
    explicit GpuStageRenderers(const GpuStageOptions& options);
    ~GpuStageRenderers();
    GpuStageRenderers(const GpuStageRenderers&) = delete;
    GpuStageRenderers& operator=(const GpuStageRenderers&) = delete;

    // The factory (installed by installGpuStageRenderers; usable without installing it).
    std::shared_ptr<const StageRendererFactory> factory() const { return m_factory; }
    // Opens the device now if it is not open yet: true when the Vulkan renderer can draw here.
    bool usable() const;
    // Why it cannot ("" when usable() or not tried yet); the device's name when it can.
    std::string problem() const;
    // usable() on a Vulkan driver that runs on the CPU ("auto" then uses the software renderer)
    bool cpuDevice() const;
    std::string deviceName() const;

private:
    std::shared_ptr<GpuStageBackend> m_backend;
    std::shared_ptr<const StageRendererFactory> m_factory;
    bool m_installed = false;
    friend std::unique_ptr<GpuStageRenderers> installGpuStageRenderers(const GpuStageOptions& options);
};

// Installs the GPU renderer choice for every structure panel built from now on.
std::unique_ptr<GpuStageRenderers> installGpuStageRenderers(const GpuStageOptions& options = {});

}  // namespace mdbridge::render::panels
