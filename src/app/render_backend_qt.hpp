#pragma once
// The `mdbridge` program's RenderBackend (commands/render_pipeline.hpp): the batch Qt application
// on the offscreen platform (render::startBatchApplication), the bundled fonts, the QPainter
// figure writers (render::writeFigure: PNG, JPEG, SVG, PDF; EPS, PS, TIFF from render core) and
// the video frame loop (render::renderVideo, the bundled LGPL FFmpeg libraries).

#include "render_pipeline.hpp"

namespace mdbridge::app {

// `program`: argv[0], for the Qt application.
RenderBackend qtRenderBackend(const char* program);

}  // namespace mdbridge::app
