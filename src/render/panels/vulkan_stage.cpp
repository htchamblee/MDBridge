// The structure view's 3D layer on the GPU (see vulkan_stage.hpp).

#include "mdbridge/render/panels/vulkan_stage.hpp"

#include <mutex>
#include <optional>
#include <shared_mutex>
#include <stdexcept>
#include <utility>

#include "mdbridge/core/log.hpp"
#include "mdbridge/render/vulkan/scene_renderer.hpp"

namespace mdbridge::render::panels {

// The process's Vulkan device and scene renderer, opened on first use and released by shutdown().
class GpuStageBackend {
public:
    explicit GpuStageBackend(GpuStageOptions options) : m_options(std::move(options)) {}
    ~GpuStageBackend() { shutdown(); }

    bool open() {
        const std::lock_guard<std::mutex> lock(m_mutex);
        if (!m_tried && !m_closed) {
            m_tried = true;
            m_context = vulkan::Context::create(m_options.context, &m_problem);
            if (m_context) {
                m_renderer = vulkan::SceneRenderer::create(*m_context, &m_problem);
                if (!m_renderer) {
                    m_context.reset();
                }
            }
        }
        return m_renderer != nullptr && !m_broken;
    }

    // nullopt when the GPU cannot draw this frame (then the caller uses the software renderer)
    std::optional<molecule::RgbaImage> render(const std::vector<molecule::Colored>& items,
                                              const molecule::Viewport& viewport, const molecule::Shading& shading,
                                              molecule::RenderQuality quality) {
        std::shared_lock<std::shared_mutex> use(m_use);
        vulkan::SceneRenderer* renderer = nullptr;
        {
            const std::lock_guard<std::mutex> lock(m_mutex);
            if (m_broken || !m_renderer) {
                return std::nullopt;
            }
            renderer = m_renderer.get();
        }
        try {
            return renderer->render(items, viewport, shading, quality);
        } catch (const std::invalid_argument&) {
            throw;  // the view or items are unusable: as the software renderer
        } catch (const std::exception& e) {
            const std::lock_guard<std::mutex> lock(m_mutex);
            if (!m_broken) {
                m_broken = true;
                core::logWarning(std::string("3D representations: the GPU renderer failed (") + e.what() +
                                 "); using the built-in renderer");
            }
            return std::nullopt;
        }
    }

    void shutdown() {
        const std::unique_lock<std::shared_mutex> exclusive(m_use);  // no render in flight
        const std::lock_guard<std::mutex> lock(m_mutex);
        m_closed = true;
        m_renderer.reset();
        m_context.reset();
    }

    std::string problem() const {
        const std::lock_guard<std::mutex> lock(m_mutex);
        return m_renderer ? std::string() : m_problem;
    }
    std::string deviceName() const {
        const std::lock_guard<std::mutex> lock(m_mutex);
        return m_context ? m_context->deviceInfo().name : std::string();
    }

private:
    GpuStageOptions m_options;
    mutable std::mutex m_mutex;
    std::shared_mutex m_use;
    std::unique_ptr<vulkan::Context> m_context;
    std::unique_ptr<vulkan::SceneRenderer> m_renderer;
    std::string m_problem;
    bool m_tried = false;
    bool m_closed = false;
    bool m_broken = false;
};

namespace {

class VulkanStageRenderer final : public StageRenderer {
public:
    explicit VulkanStageRenderer(std::shared_ptr<GpuStageBackend> backend) : m_backend(std::move(backend)) {}
    molecule::RgbaImage render(const std::vector<molecule::Colored>& items, const molecule::Viewport& viewport,
                               const molecule::Shading& shading, molecule::RenderQuality quality) const override {
        if (std::optional<molecule::RgbaImage> image = m_backend->render(items, viewport, shading, quality)) {
            return std::move(*image);
        }
        return m_software.render(items, viewport, shading, quality);
    }
    std::string name() const override { return "VulkanRenderer"; }

private:
    std::shared_ptr<GpuStageBackend> m_backend;
    molecule::SoftwareRenderer m_software;
};

class GpuFactory final : public StageRendererFactory {
public:
    explicit GpuFactory(std::shared_ptr<GpuStageBackend> backend) : m_backend(std::move(backend)) {}
    StageRendererChoice choose(std::string_view name) const override {
        const bool gpu = name == "auto" || name == "vtk" || name == "vulkan";
        if (gpu && m_backend->open()) {
            const std::lock_guard<std::mutex> lock(m_mutex);
            if (!m_vulkan) {
                m_vulkan = std::make_shared<VulkanStageRenderer>(m_backend);
            }
            return {m_vulkan, {}};
        }
        if (name == "auto") {
            return defaultStageRenderer("software");
        }
        return defaultStageRenderer(name);  // RN-16: "vtk" / "vulkan" with the note
    }

private:
    std::shared_ptr<GpuStageBackend> m_backend;
    mutable std::mutex m_mutex;
    mutable std::shared_ptr<const StageRenderer> m_vulkan;
};

}  // namespace

GpuStageRenderers::GpuStageRenderers(const GpuStageOptions& options)
    : m_backend(std::make_shared<GpuStageBackend>(options)), m_factory(std::make_shared<GpuFactory>(m_backend)) {}

GpuStageRenderers::~GpuStageRenderers() {
    if (m_installed && stageRendererFactory() == m_factory) {
        setStageRendererFactory(nullptr);
    }
    m_backend->shutdown();
}

bool GpuStageRenderers::usable() const { return m_backend->open(); }
std::string GpuStageRenderers::problem() const { return m_backend->problem(); }
std::string GpuStageRenderers::deviceName() const { return m_backend->deviceName(); }

std::unique_ptr<GpuStageRenderers> installGpuStageRenderers(const GpuStageOptions& options) {
    auto renderers = std::make_unique<GpuStageRenderers>(options);
    setStageRendererFactory(renderers->factory());
    renderers->m_installed = true;
    return renderers;
}

}  // namespace mdbridge::render::panels
