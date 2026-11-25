#pragma once
// MDATK's standalone charts (mdatk/render/charts.py; MDATK's own code, ported): one figure per
// series (a time chart) or matrix (a residue x frame heatmap), written in every format of
// [output.charts]. Qt-free: the figures are drawn here and handed to a FigureWriter, the output
// backend (the Qt one is render::writeFigure, figure_qt.hpp); the `charts` and `render` commands
// call renderCharts.
//
// What MDATK does, kept: the chosen series (output.charts.series, or every series but trace/ and
// proximity/), the heatmaps of every matrix but those with include_heatmaps, the chart sizes and
// chrome, smoothing ("Per frame" + "Running mean (n frames)"), the target line, labels and the
// title template, file names <out_dir>/<prefix>_<safe key>[_<stamp>].<format>, the log
// "Wrote N chart file(s) to D", the progress "Chart <key>" and the "chart series K is not
// available" warning. MDBridge's improvements (docs/DIFFERENCES.md "Rendering: charts"):
//
//   RN-01  matrix chart cells at the midpoints between frame times (plots.hpp columnEdges)
//   RN-05  matrix chart height capped at 2^15 pixels, with a warning
//   RN-10  a matrix key in output.charts.series is charted as a matrix chart (pair maps: warning)
//   RN-11  formats checked before anything is drawn (chartFormatProblem)
//   RN-14  the NaN-aware running mean (plots.hpp runningMean)
//   RN-50  the title and the matrix charts' colour bar in the theme's colours
//   D5     styles: "mdatk" (or "") and Matplotlib's "default" look; other names warn, MDATK chrome

#include <cstdint>
#include <filesystem>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/render/figure.hpp"
#include "mdbridge/render/text.hpp"

namespace mdbridge::render::charts {

// The formats a chart can be written in (lower case): png, svg, pdf, jpg, jpeg, tif, tiff, eps, ps.
std::span<const std::string_view> chartFormats();
// A format as MDATK names the file: lower case, leading dots removed (".PNG" -> "png").
std::string normalizeFormat(std::string_view format);
// The error for output.charts.formats (empty: png), or nullopt when every one can be written:
// "chart format 'bogus' is not supported (choose from png, svg, ...)" (RN-11; check() and
// renderCharts report it).
std::optional<std::string> chartFormatProblem(const std::vector<std::string>& formats);

// MDATK's _safe(key): runs of characters other than A-Z a-z 0-9 _ . + - become "_", then leading and
// trailing "_" are removed ("rmsd@Nope/x" -> "rmsd_Nope_x").
std::string safeKey(std::string_view key);
// _format_title(template, label, key): the label when the template is empty, else {label} and
// {key} replaced.
std::string formatTitle(std::string_view tmpl, std::string_view label, std::string_view key);

// What renderCharts draws, in order: ("series", key) for each chosen series, then ("matrix", key)
// for each requested matrix (RN-10) and, with include_heatmaps, the other matrices. Logs the
// warnings for missing keys and pair maps.
struct ChartJob {
    std::string kind;  // "series" or "matrix"
    std::string key;
};
std::vector<ChartJob> chartJobs(const core::Dataset& dataset, const core::Project& project);

// The chart figures, drawn (Figure::draw(fonts)) and ready to be written. Throw std::out_of_range
// (Python's KeyError) for a missing key. The style warnings (D5) are logged by each call, as MDATK
// warned for each chart.
Figure seriesChart(const core::Dataset& dataset, const core::Project& project, const std::string& key,
                   const FontSet& fonts);
Figure matrixChart(const core::Dataset& dataset, const core::Project& project, const std::string& key,
                   const FontSet& fonts);

// The pixel height MDATK's matrix chart has (max(height, 0.18 * rows + 1.5) inches at dpi) and
// MDBridge's (capped at kMaxChartPixels, RN-05).
inline constexpr int kMaxChartPixels = 32768;
int matrixChartHeightPx(double heightInches, std::size_t rows, double dpi);

// Writes a drawn figure to `path` in `format` (one of chartFormats()); throws std::runtime_error
// (naming the file) when it cannot.
using FigureWriter = std::function<void(const Figure& figure, const std::filesystem::path& path, const std::string& format)>;

// render_charts(dataset, project, out_dir, reporter, stamp): out_dir is created, each job's chart is
// drawn and written in every format, the paths in order (one per job and format) are returned.
// Throws core::ValueError (chartFormatProblem) before anything is written; core::Cancelled when
// the reporter is cancelled (checked before each chart).
std::vector<std::string> renderCharts(const core::Dataset& dataset, const core::Project& project,
                                      const std::filesystem::path& outDir, const core::Reporter& reporter,
                                      const std::string& stamp, const FontSet& fonts, const FigureWriter& writer);

}  // namespace mdbridge::render::charts
