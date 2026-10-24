#include "mdbridge/render/view.hpp"

#include "software_view.hpp"
#include "vulkan_view.hpp"

namespace mdbridge::render {

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
