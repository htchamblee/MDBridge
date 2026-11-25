#include "mdbridge/render/view.hpp"

#include <algorithm>

#include "software_view.hpp"
#include "vulkan_view.hpp"

namespace mdbridge::render {

molecule::Viewport fittedViewport(const MoleculeFrame& frame, int width, int height) {
    molecule::Viewport vp;
    vp.rotation = frame.rotation;
    vp.zrange = frame.zrange;
    vp.width = std::max(width, 1);
    vp.height = std::max(height, 1);
    const double dx = frame.xlim[1] - frame.xlim[0], dy = frame.ylim[1] - frame.ylim[0];
    const double cx = 0.5 * (frame.xlim[0] + frame.xlim[1]), cy = 0.5 * (frame.ylim[0] + frame.ylim[1]);
    const double target = static_cast<double>(vp.width) / static_cast<double>(vp.height);
    double hx = 0.5 * dx, hy = 0.5 * dy;
    if (dy > 0 && dx / dy < target) {
        hx = hy * target;  // window wider than the data: widen x
    } else if (dx > 0) {
        hy = hx / target;
    }
    vp.xlim = {cx - hx, cx + hx};
    vp.ylim = {cy - hy, cy + hy};
    return vp;
}

std::optional<Renderer> rendererFromName(std::string_view name) {
    if (name == "auto") {
        return Renderer::Auto;
    }
    if (name == "vulkan") {
        return Renderer::Vulkan;
    }
    if (name == "software") {
        return Renderer::Software;
    }
    return std::nullopt;
}

const char* rendererName(Renderer renderer) {
    switch (renderer) {
    case Renderer::Auto:
        return "auto";
    case Renderer::Vulkan:
        return "vulkan";
    case Renderer::Software:
        return "software";
    }
    return "auto";
}

View::View(QWindow* window, Renderer renderer, QString description)
    : m_window(window), m_renderer(renderer), m_description(std::move(description)) {}

View::~View() { delete m_window; }  // QPointer: null if a window container deleted it already

void View::countFrame() { emit frameRendered(++m_frames); }

ViewChoice createView(Renderer requested) {
    QString reason = QStringLiteral("requested");
    if (requested != Renderer::Software) {
        if (std::unique_ptr<View> view = createVulkanView(&reason)) {
            return {std::move(view), QString()};
        }
    }
    return {createSoftwareView(reason), reason};
}

std::unique_ptr<View> createSoftwareView(const QString& reason) {
    return std::make_unique<SoftwareView>(reason);
}

}  // namespace mdbridge::render
