#pragma once
// The plot panels' arithmetic, shared with the standalone charts (mdatk/render/charts.py uses the
// same running mean and matrix cells): MDATK's own code, ported, with MDBridge's improvements
// (docs/DIFFERENCES.md "Rendering: plot panels"):
//
//   runningMean   panels._running_mean (edge padding, window clipped to the series) - NaN-aware
//                 (RN-14): a window's mean of its finite values, NaN only where the whole window
//                 is not finite (MDATK: one NaN blanked every window it fell in)
//   columnEdges   the time-axis cell edges of a frames x rows matrix: midpoints between frame
//                 times, half a step outward at both ends (RN-01; MDATK: t[1] - t[0] for every
//                 column, cells misplaced for non-uniform times); one frame: t0 -/+ 0.5
//   heatmapNorm   the sequential heatmap's vmin / vmax: options vmin (0), vmax (0 or missing:
//                 the largest finite value, 1 without one), vmax at least 1e-9 as in MDATK; vmin
//                 above vmax swapped with a warning (RN-04; MDATK's colour bar then rewrote the
//                 norm to a callback-order artefact, e.g. 0.9 .. 5 for vmin 5, vmax 1)

#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace mdbridge::render::panels {

std::vector<double> runningMean(std::span<const double> y, long long window);

std::vector<double> columnEdges(std::span<const double> time);

struct Norm {
    double vmin = 0.0;
    double vmax = 1.0;
    bool swapped = false;  // vmin was above vmax (the caller warns)
};
// vminOption / vmaxOption: the options as given (nullopt: not set; vmax 0 counts as not set)
Norm heatmapNorm(std::span<const double> values, std::optional<double> vminOption, std::optional<double> vmaxOption);

}  // namespace mdbridge::render::panels
