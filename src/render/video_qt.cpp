// The frame loop's QPainter compositor and render_video (video_qt.hpp).
#include "mdbridge/render/video_qt.hpp"

#include <cstring>

#include "mdbridge/render/figure_qt.hpp"
#include "mdbridge/render/text_qt.hpp"

namespace mdbridge::render {

void QtFrameCompositor::setBackground(const DisplayList& background) {
    // composeLayers paints on a premultiplied copy: convert once, so each frame only copies
    m_background = renderImage(background).convertToFormat(QImage::Format_ARGB32_Premultiplied);
}

std::vector<std::uint8_t> QtFrameCompositor::compose(const DisplayList& layer) const {
    const QImage frame = composeLayers(m_background, layer);  // RGBA8888
    const auto w = static_cast<std::size_t>(frame.width()), h = static_cast<std::size_t>(frame.height());
    std::vector<std::uint8_t> rgba(w * h * 4);
    for (std::size_t y = 0; y < h; ++y) {
        std::memcpy(rgba.data() + y * w * 4, frame.constScanLine(static_cast<int>(y)), w * 4);
    }
    return rgba;
}

std::filesystem::path renderVideo(const core::Dataset& dataset, const core::Project& project,
                                  const std::filesystem::path& path, const core::Reporter& reporter,
                                  const video::VideoRenderOptions& options) {
    QtFrameCompositor compositor;
    return video::renderVideo(dataset, project, path, bundledFontSet(), compositor, reporter, options);
}

}  // namespace mdbridge::render
