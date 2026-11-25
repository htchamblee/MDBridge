#pragma once
// Axis limits, tick locations and tick labels as Matplotlib 3.11 computes them for MDATK's figures
// (linear axes, rcParams at Matplotlib's defaults, which MDATK does not change), and MDATK's
// thin_ticks. No Qt: the figure framework and the tests use the same code.
//
//   nonsingular        matplotlib.transforms._nonsingular
//   tickSpace          XAxis / YAxis.get_tick_space (the 'auto' number of bins)
//   MaxNLocator        matplotlib.ticker.MaxNLocator, AutoLocator, MDATK's integer_axis
//   fixedLocatorTicks  matplotlib.ticker.FixedLocator.tick_values (colour bars with ticks=range(n))
//   ScalarFormatter    matplotlib.ticker.ScalarFormatter (offset, scientific notation, U+2212 minus)
//   setLimits          Axis._set_lim on a linear axis (what set_xlim / set_ylim store)
//   autoscale          Axes.autoscale_view of one linear axis (margins, sticky edges, view_limits)
//   ticksToDraw        Axis._update_ticks: which located ticks are drawn
//   colorbarLimits     Colorbar._process_values / _draw_all: the long axis of a linear colour bar
//   thinTicks          mdatk.render.panels.thin_ticks
//
// Ported from Matplotlib 3.11.2 (lib/matplotlib/ticker.py, transforms.py, axis.py, axes/_base.py,
// colorbar.py) under Matplotlib's licence (packaging/licenses/Matplotlib-LICENSE.txt):
// Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved.
// Summary of the changes made to Matplotlib (its licence, clause 3): the algorithms above were
// translated from Python to C++ for linear axes only (no log / function scales, no 3D axes, no
// "round_numbers" autolimit mode, no locale), with rcParams passed as options; the arithmetic keeps
// Python's and numpy's operations and their order so that results are equal bit for bit.
//
// Results equal Matplotlib's bit for bit only when the translation unit is compiled without
// floating-point contraction (-ffp-contract=off, MSVC /fp:precise), as mdbridge_render_core is.

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace mdbridge::render {

struct Range {
    double lo = 0.0;
    double hi = 0.0;
    bool operator==(const Range&) const = default;
};

// transforms._nonsingular: expands an empty or tiny interval, gives (-expander, expander) for
// non-finite or near-zero input, and swaps vmin > vmax unless `increasing` is false.
Range nonsingular(double vmin, double vmax, double expander = 0.001, double tiny = 1e-15,
                  bool increasing = true);

enum class AxisDirection { X, Y };

// FontProperties.get_size_in_points of a size in points: Matplotlib uses at least 1 pt.
double fontSizeInPoints(double size);

// XAxis / YAxis.get_tick_space: how many ticks fit on an axis `lengthPoints` long (the axes'
// width or height in points, axes size in inches x 72) with tick labels of `labelSizePoints`
// (the size is clamped as fontSizeInPoints does). Labels are assumed 3 x their size wide on x and
// 2 x their size high on y; a non-positive size gives 2^31 - 1.
int tickSpace(AxisDirection direction, double lengthPoints, double labelSizePoints);

enum class Prune { None, Lower, Upper, Both };

struct MaxNLocatorParams {
    std::optional<int> nbins = 10;  // nullopt: 'auto' (from the tick space, clipped to 1..9)
    std::vector<double> steps;      // empty: Matplotlib's default 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10
    bool integer = false;
    bool symmetric = false;
    Prune prune = Prune::None;
    int minNTicks = 2;
};

class MaxNLocator {
public:
    // Throws std::invalid_argument with Matplotlib's message for steps that are not an increasing
    // sequence within [1, 10], and for nbins < 1 (Matplotlib then computes garbage from a division
    // by zero or a negative bin count).
    explicit MaxNLocator(MaxNLocatorParams params = {});

    // matplotlib.ticker.AutoLocator: nbins 'auto', steps 1, 2, 2.5, 5, 10 (the default locator of
    // every axis and colour bar).
    static MaxNLocator autoLocator();
    // MaxNLocator(integer=True): MDATK's integer_axis (nbins 10, default steps).
    static MaxNLocator integerLocator();

    // tick_values(vmin, vmax). `tickSpace` is the axis' get_tick_space(), used with nbins 'auto';
    // without an axis (nullopt) Matplotlib uses 9 bins. Ticks beyond the interval are included,
    // as Matplotlib returns them.
    std::vector<double> tickValues(double vmin, double vmax, std::optional<int> tickSpace = std::nullopt) const;

    // view_limits (autolimit_mode 'data'): the interval made non-singular.
    Range viewLimits(double dmin, double dmax) const;

    // Locator.nonsingular: what autoscaling and set_xlim apply (expander 0.05).
    static Range nonsingularLimits(double v0, double v1);

    // The number of bins tickValues uses.
    int bins(std::optional<int> tickSpace = std::nullopt) const;

    bool isAuto() const { return !m_nbins.has_value(); }
    const MaxNLocatorParams& params() const { return m_params; }

private:
    std::vector<double> rawTicks(double vmin, double vmax, std::optional<int> tickSpace) const;

    MaxNLocatorParams m_params;
    std::optional<int> m_nbins;
    std::vector<double> m_extendedSteps;
};

// FixedLocator(locs, nbins).tick_values: every step-th location, step = ceil(n / max(nbins, 2)), starting
// at the offset whose ticks come closest to zero; all of them without nbins.
std::vector<double> fixedLocatorTicks(std::span<const double> locs, std::optional<int> nbins = std::nullopt);

struct ScalarFormatterOptions {
    bool useOffset = true;               // axes.formatter.useoffset
    std::optional<double> fixedOffset;   // set_useOffset(number): this offset, never computed
    int offsetThreshold = 4;             // axes.formatter.offset_threshold
    int powerLimitLow = -5;              // axes.formatter.limits
    int powerLimitHigh = 6;
    bool scientific = true;              // set_scientific
    bool useMathText = false;            // axes.formatter.use_mathtext (or text.usetex)
    bool unicodeMinus = true;            // axes.unicode_minus: '-' written as U+2212
};

// matplotlib.ticker.ScalarFormatter. Like Matplotlib's, it keeps the offset, order of magnitude
// and format of the last setLocs; formatTicks = setLocs + label of each value.
class ScalarFormatter {
public:
    explicit ScalarFormatter(ScalarFormatterOptions options = {});

    // set_locs: `view` is the axis' view interval (either order).
    void setLocs(std::span<const double> locs, Range view);
    // __call__: the label of one value ('' when the last setLocs had no locations).
    std::string operator()(double value) const;
    // format_ticks: setLocs, then the label of every value.
    std::vector<std::string> formatTicks(std::span<const double> locs, Range view);
    // get_offset: the offset / scientific-notation text ('1e6', '+1.234e3', '1e−5+2e2', ...).
    std::string offsetText() const;
    // format_data: '1.5e3', '−2e−4', '7' (exponent 0); mathtext: '1.5 \times 10^{3}'.
    // `numpyValue` selects how the significand is rounded to 10 decimals: numpy's round (the
    // value is a numpy float, as a computed offset is) or Python's round (a Python float or int).
    std::string formatData(double value, bool numpyValue = true) const;

    double offset() const { return m_offset; }
    int orderOfMagnitude() const { return m_orderOfMagnitude; }
    int decimals() const { return m_decimals; }  // the precision of the '%1.Nf' format
    const ScalarFormatterOptions& options() const { return m_options; }

private:
    void computeOffset(Range view);
    void setOrderOfMagnitude(Range view);
    void setFormat(Range view);
    std::string fixMinus(std::string text) const;

    ScalarFormatterOptions m_options;
    std::vector<double> m_locs;
    double m_offset = 0.0;
    bool m_offsetFromNumpy = true;
    int m_orderOfMagnitude = 0;
    int m_decimals = 0;
    bool m_hasFormat = false;
};

// Axis._set_lim on a linear axis: (v0, v1) made non-singular (expander 0.05), orientation kept.
// Throws std::invalid_argument("Axis limits cannot be NaN or Inf") as set_xlim does.
Range setLimits(double v0, double v1);

struct AutoscaleOptions {
    double margin = 0.05;               // axes.xmargin / axes.ymargin or Axes.margins
    std::vector<double> stickyEdges;    // the artists' sticky edges on this axis (images: extent)
    bool tight = false;                 // autoscale_view(tight=True): no view_limits step
    bool inverted = false;              // the axis is inverted (set_xbound keeps it so)
};

// Axes.autoscale_view of one linear, unshared axis with `locator` as its major locator:
// [dataMin, dataMax] are the data limits (non-finite ends ignored; none finite: Matplotlib's
// default (-inf, inf) handling), then nonsingular, the margin, the sticky edges, view_limits and
// set_xbound / set_xlim.
Range autoscale(double dataMin, double dataMax, const MaxNLocator& locator, const AutoscaleOptions& options = {});

// Axis._update_ticks: indices of the locations drawn, those within the view interval or within
// 1e-10 of its length outside it.
std::vector<std::size_t> ticksToDraw(std::span<const double> locs, Range view);

// The long axis of a colour bar for Normalize(vmin, vmax) and a colormap of `colormapSize`
// colours: the norm's limits made non-singular (expander 0.1), the first and last boundary
// (vmin + b (vmax - vmin) for b = 0 and 1), and set_ylim. Inverted axes get the limits swapped.
// For vmin <= vmax (the heatmap clamps vmax first, RN-04): with vmin > vmax Matplotlib's result
// depends on its norm callbacks (setting vmin re-enters the colour bar with the old vmax).
Range colorbarLimits(double vmin, double vmax, bool inverted = false);

// mdatk.render.panels.thin_ticks: 0, step, 2 step, ... below n, step = max(1, ceil(n / maxTicks)).
// maxTicks < 1 is taken as 1 (MDATK: ZeroDivisionError for 0, every label for negative values;
// docs/DIFFERENCES.md).
std::vector<long long> thinTicks(long long n, long long maxTicks = 20);

}  // namespace mdbridge::render
