#pragma once
// The video frame loop (video/frame_loop.hpp) with the QPainter backend: frames composed as
// figure_qt.hpp's composeLayers composes them (the static background painted once into a raster
// image, each frame's animated layer painted over a copy of it, text filled from the bundled fonts'
// outlines), and render_video for the `render` command and the GUI.

#include <cstdint>
#include <filesystem>
#include <vector>

#include <QImage>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "video/frame_loop.hpp"

namespace mdbridge::render {

// FrameCompositor on QImage: setBackground renders the static layer once (renderImage); compose
// paints a frame's layer over a copy of it (composeLayers) and returns its RGBA8888 bytes.
// Thread-safe after setBackground (raster painting only, no QFont).
class QtFrameCompositor final : public video::FrameCompositor {
public:
    void setBackground(const DisplayList& background) override;
    std::vector<std::uint8_t> compose(const DisplayList& layer) const override;
    // the background as composeLayers paints on it (premultiplied ARGB32; null before setBackground)
    const QImage& background() const { return m_background; }

private:
    QImage m_background;
};

// MDATK's render_video(dataset, project, path, reporter): video::renderVideo with the bundled
// fonts and QtFrameCompositor. Raster painting only; the batch commands run it in their
// application (startBatchApplication).
std::filesystem::path renderVideo(const core::Dataset& dataset, const core::Project& project,
                                  const std::filesystem::path& path, const core::Reporter& reporter = core::Reporter(),
                                  const video::VideoRenderOptions& options = {});

}  // namespace mdbridge::render
