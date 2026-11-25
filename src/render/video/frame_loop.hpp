#pragma once
// MDATK's video frame loop (mdatk/render/video.py: video_positions, render_video; MDATK's own code,
// ported) on MDBridge's figure (render/panels build.hpp) and encoder (encoder.hpp). Qt-free: the
// frames are composed by a FrameCompositor, whose Qt implementation (QPainter raster, the bundled
// fonts) is render::QtFrameCompositor in mdbridge::render (mdbridge/render/video_qt.hpp, which also
// has the render::renderVideo the `render` command calls).
//
// The loop, as MDATK's blitting:
//   1. build_figure once (the frozen layout); the artists of panels.animatedArtists() are marked
//      animated; the static layer (everything else: axes, ticks, full traces, heatmaps) is drawn
//      once and handed to the compositor as the background
//   2. for each position (0, s, 2s, ... below the frame count; output.video.frame_stride):
//      check for cancellation, update(pos) of every panel and the animated layer's display list (on
//      the calling thread, in order: the panels and the 3D stage renderer they call are not shared
//      between threads), then the frame composed (background + layer) on a pool of threads, and
//      written in order (on the calling thread) to the sink: the encoder, or a test
//   3. progress "Video frame k/n" every n/200 frames and after the last, as MDATK reports it
//
// Composition is a pure function of the background and the frame's display list, so the frames
// are the same bytes whatever the number of threads (tested). The 3D layer is drawn by the
// structure panel's stage renderer (panels::stageRendererFactory(): the built-in software renderer,
// or the GPU renderer once the program installs a factory that knows it).
//
// Differences from MDATK v0.0.4.1 (docs/DIFFERENCES.md, "Rendering: the video frame loop"):
// - RN-02: no analysed frames -> VideoError("no analysed frames to render") before anything is
//   written (MDATK failed in the panels or in FFmpeg).
// - RN-06: the canvas is exactly figure_size(w, h) pixels, so MDATK's "canvas size not even" /
//   "choose a DPI that divides the resolution" errors cannot happen.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <span>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/panels/build.hpp"

namespace mdbridge::render::video {

// video_positions(n_frames, stride): 0, s, 2s, ... below n_frames (s = max(1, stride); the last
// frame is left out unless it falls on the stride, RN-17: kept as MDATK has it).
std::vector<std::size_t> videoPositions(std::size_t frames, std::int64_t stride);

// Background + animated layer -> a frame.
class FrameCompositor {
public:
    virtual ~FrameCompositor();
    // Once, before the frames: the static layer (DisplayLayer::Static) of the figure.
    virtual void setBackground(const DisplayList& background) = 0;
    // A frame: the background with `layer` painted over it, as layer.width x layer.height straight
    // RGBA pixels (8 bits per channel), rows top first, no padding. Called from several threads at
    // once (after setBackground).
    virtual std::vector<std::uint8_t> compose(const DisplayList& layer) const = 0;
};

struct FrameLoopOptions {
    // threads composing frames: 0 = the hardware's (at most 8); 1 = on the calling thread
    unsigned threads = 0;
};

// The number of composing threads `options` asks for on this machine (>= 1).
unsigned frameLoopThreads(const FrameLoopOptions& options);

// Receives frame k (the k-th of the positions) on the calling thread, in order.
using FrameSink = std::function<void(std::size_t k, std::vector<std::uint8_t>& rgba)>;

// Steps 1-3 above on a built figure: marks its animated artists, sets the compositor's background,
// then the frames of `positions`. reporter.checkCancel() before each frame (core::Cancelled
// propagates; frames still being composed are dropped); errors of the sink, the panels or the
// compositor propagate the same way (a compositor's error when its frame's turn comes).
void runFrameLoop(panels::BuiltFigure& figure, const FontSet& fonts, std::span<const std::size_t> positions,
                  FrameCompositor& compositor, const FrameSink& sink, const core::Reporter& reporter,
                  const FrameLoopOptions& options = {});

struct VideoRenderOptions {
    FrameLoopOptions loop;
    int encoderThreads = 0;  // EncoderOptions::threads (0: the encoder's choice)
};

// render_video(dataset, project, path, reporter): build_figure at output.video's size and dpi,
// the folder of `path` created, "Rendering N video frame(s) at WxH, F fps -> PATH" logged, the
// encoder opened (output.video's codec, crf, fps; its warnings go to core::logWarning), the frame
// loop (the frames are encoded in order on a thread of their own while the next ones are drawn:
// the same file as encoding them on the calling thread), the file finished (temporary name renamed) and "Wrote PATH (N frames)" logged. Returns
// `path`. On failure or cancellation the temporary file is removed and nothing is left at `path`
// (an existing file there is kept). Throws VideoError for encoding errors and RN-02.
std::filesystem::path renderVideo(const core::Dataset& dataset, const core::Project& project,
                                  const std::filesystem::path& path, const FontSet& fonts, FrameCompositor& compositor,
                                  const core::Reporter& reporter, const VideoRenderOptions& options = {});

}  // namespace mdbridge::render::video
