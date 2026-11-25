#pragma once
// What the interface shows of an opened project, without Qt: the `mdbridge` interface
// (app/gui/project_loader.cpp) runs these on a worker thread, the tests run them directly
// (tests/cli/test_project_view.cpp).
//
//   loadReferenceFrame  MDATK's Pipeline.structure_preview() (the GUI's "Load real frame"): the
//                       system of the project's input files (with cleaning on, the cleaned files
//                       MDATK left for the same inputs and settings when they exist, else the raw
//                       files, as MDATK; cpptraj never runs) and the structure analyses (trace,
//                       atoms; proximity when the colouring or glow needs it) on the first
//                       analysed frame (inputs.start / stop / step)
//   structurePreview    that frame as the example dataset's StructurePreview (structure_preview's
//                       dictionary)
//   structureView       the 3D view's scene of that frame: the structure panel's 3D layer
//                       (render/panels/structure.hpp StructurePanel::stageFrame) with the
//                       project's representations, colours, rotation, depth cue and glow
//   previewFigure       the read-only figure preview: MDATK's GUI canvas (gui/canvas.py render,
//                       app.py refresh_canvas): the project's figure built with the example
//                       dataset at the output.video size scaled to the space available, the
//                       watermark "Example data – layout preview" (with " (structure: reference
//                       frame)" once the reference frame is loaded), panels updated to 35 % of the
//                       example time axis
//
// MDATK's GUI had no separate 3D view: its structure showed only in the figure canvas. MDBridge's
// 3D view draws a project without representations with a cartoon (kDefaultViewStyle), as the
// figure's 2D trace is no 3D picture (docs/DIFFERENCES.md "Interface").

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/render/example.hpp"
#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/software.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/text.hpp"

namespace mdbridge::app {

struct ReferenceFrame {
    core::Dataset dataset{{}, {}};  // one frame: the trace and atoms results (and proximity)
    std::int64_t frame = 0;         // its frame number in the trajectory
    std::string topology;           // the topology file read
    bool cleaned = false;           // the cleaned files of an earlier MDATK run
    std::vector<std::string> notes;  // reader warnings, the raw files used instead of cleaned ones
};

// Throws what io::loadSystem throws (missing or unreadable files), core::ValueError for a system
// without coordinates or an empty frame selection, core::Cancelled when the reporter is cancelled,
// and what the analyses throw (e.g. a selection error). The reporter receives the runner's log and
// progress.
ReferenceFrame loadReferenceFrame(const core::Project& project, const core::Reporter& reporter = core::Reporter());

// structure_preview()'s dictionary of a loaded reference frame (render::exampleDataset's input);
// without atoms when the dataset has no atoms result. Throws core::ValueError without a trace.
render::StructurePreview structurePreview(const core::Dataset& reference);

// The representation the 3D view draws for a project without any.
inline constexpr std::string_view kDefaultViewStyle = "cartoon";

// The project the 3D view's structure panel is built from: the layout is one structure panel, the
// representations are the project's (else one kDefaultViewStyle representation), the 3D layer
// renderer is "software" (the view draws the scene itself; no device is opened for it).
core::Project structureViewProject(const core::Project& project);

// The analyses loadReferenceFrame runs (panelRequiredAnalyses of structureViewProject).
std::vector<std::string> referenceAnalyses(const core::Project& project);

struct StructureView {
    std::vector<render::molecule::Colored> items;  // the scene's primitives (ions included)
    render::molecule::Viewport viewport;           // rotation, square data window, depth range
    render::molecule::Shading shading;             // depth cue toward the panel surface, glow
    std::array<float, 3> background{};             // the panel surface (RGB 0..1)
    std::size_t residues = 0;
    std::size_t atoms = 0;
};

// The 3D view of `reference` (frame 0 of the dataset). Throws core::ValueError when the frame has
// no structure to draw ("the structure view has no data: <reason>", "no 3D representation could be
// drawn").
StructureView structureView(const core::Project& project, const core::Dataset& reference);

// ---- the figure preview

inline constexpr double kPreviewFrame = 0.35;  // the part of the example time axis the preview shows
inline constexpr int kPreviewMinimumSize = 50;  // pixels; smaller spaces draw nothing (MDATK)

// "Example data – layout preview", with " (structure: reference frame)" when `referenceFrame`.
std::string previewWatermark(bool referenceFrame);

// The preview's figure size: output.video's size and dpi scaled to fit `width` x `height`
// (min of the two ratios), at least 2 x 2 pixels.
struct PreviewSize {
    int width = 2;
    int height = 2;
    double dpi = 72.0;
};
PreviewSize previewSize(const core::Project& project, int width, int height);

// The frame position the preview shows: int(kPreviewFrame * max(0, frames - 1)).
std::size_t previewPosition(std::size_t frames);

// The preview figure of `example` (render::exampleDataset) drawn (Figure::draw) at previewSize, the
// panels updated to previewPosition. Throws what render::panels::buildFigure throws (MDATK's GUI:
// "Cannot draw this layout: <message>").
render::panels::BuiltFigure previewFigure(const core::Project& project, const core::Dataset& example,
                                          const render::FontSet& fonts, int width, int height,
                                          const std::string& watermark);

}  // namespace mdbridge::app
