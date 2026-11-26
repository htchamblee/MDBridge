#pragma once
// The density part of MDATK's Pipeline (mdatk/pipeline.py density_frame_atoms, fit_density,
// ensure_density_fit) on MDBridge's analysis pipeline and chem's rigid-body fit (chem/density_fit.hpp):
// for `density fit` / `simulate`, `render` and `preview` (auto_fit), and the analysis window's
// Global search / Refine. Qt-free.
//
// Differences from MDATK v0.0.4.1 (docs/DIFFERENCES.md "Density fit"): non-finite map voxels carry
// no density in the fit (CH-11; the fit's message counts them), and the density centre of a local
// start and its level use the finite voxels (CH-11 / RN-33).

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/chem/density_fit.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/rigid_transform.hpp"

namespace mdbridge::app {

struct DensityFrameAtoms {
    std::vector<core::Vec3> xyz;        // the frame's atoms, superposed as the structure view shows them
    std::vector<std::string> elements;  // per atom
    std::int64_t frame = 0;             // the trajectory frame
    std::vector<core::Vec3> reference;  // the backbone trace of the reference frame (the anchor)
};

// density_frame_atoms(frame): the atoms analysis on the first analysed frame ("ref", "" or that
// frame's number) and on frame `frame`, superposed onto it as the video frames are, and the trace
// of the reference frame. Throws PipelineError("no frames selected"), PipelineError("density fit
// frame 999 is outside the trajectory (0-23)"), core::ValueError("invalid literal for int() with
// base 10: 'x'") for a frame that is no number, and what the analyses and readers throw. The
// analyses' warnings go to `reporter` (the pipeline's when null).
DensityFrameAtoms densityFrameAtoms(analysis::Pipeline& pipeline, std::string_view frame = "ref",
                                    const core::Reporter* reporter = nullptr);

// fit_density(search, start): the rigid-body fit of the fit frame into density.map_path (search
// "global" or "manual"; nullopt: density.search). A manual fit starts at `start` (12 numbers), else
// at the stored placement (anchored to the current reference), else at the map's density centre
// above density.level. Stores transform, anchor, fit_cc and fit_info in the pipeline's project
// (Pipeline::setDensity), logs "Density fit (...) of frame ref into MAP at 6 Å" and "Density fit:
// <fit_info>", and returns the fit. Throws PipelineError("no density map selected"),
// PipelineError("density map: <reason>"), PipelineError("density map: stored placement: ..."),
// PipelineError("density fit: <FitError>"), core::Cancelled, and what densityFrameAtoms throws.
chem::FitResult fitDensity(analysis::Pipeline& pipeline, std::optional<std::string> search = std::nullopt,
                           const std::optional<std::vector<double>>& start = std::nullopt,
                           const core::Reporter* reporter = nullptr);

// ensure_density_fit: before a video or preview, a shown map without a placement is fitted
// (density.auto_fit; off: "WARNING: the density map has no placement and automatic fitting is off;
// it is not drawn"); a failed fit is the warning "WARNING: density map not drawn: <reason>" and the
// figure is drawn without the map. Cancelled propagates.
void ensureDensityFit(analysis::Pipeline& pipeline, const core::Reporter* reporter = nullptr);

}  // namespace mdbridge::app
