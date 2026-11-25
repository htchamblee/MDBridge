#pragma once
// The figure model of MDBridge's charts and video frames: a figure of axes on a grid (GridSpec with
// row height ratios), twin axes, colour bars, legends, texts, lines, scatter points, line segments,
// images and rectangles, laid out as Matplotlib 3.11's "constrained" layout engine lays out MDATK's
// figures (mdatk/render/video.py build_figure, charts.py), then frozen. Qt-free: what is drawn comes
// out as a display list (display_list.hpp) that the Qt backends paint (figure_qt.hpp); the record
// writer (figure_dump.hpp) describes a drawn figure as tests/parity/mpl_walker.py describes MDATK's.
//
// Coordinates follow Matplotlib: figure pixels with the origin at the bottom left ("display"
// coordinates); a figure is width x height pixels at `dpi` (points -> pixels = pt * dpi / 72).
//
// Algorithms ported from Matplotlib 3.11.2 (lib/matplotlib/_constrained_layout.py, _layoutgrid.py,
// text.py, axis.py, axes/_base.py, spines.py, legend.py, offsetbox.py, gridspec.py, colorbar.py,
// transforms.py) under Matplotlib's licence (packaging/licenses/Matplotlib-LICENSE.txt):
// Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved.
// Summary of the changes made to Matplotlib (its licence, clause 3): the geometry of the artists
// MDATK uses (text boxes, tick labels, axis labels, titles, spines, legends with Line2D handles,
// colour bars on the right of one parent) and the constrained layout of one GridSpec were translated
// from Python to C++; the kiwisolver constraint system is replaced by its closed-form solution
// (MDATK's grids never over-constrain), the figure has an exact pixel size (docs/DIFFERENCES.md).

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "mdbridge/render/color.hpp"
#include "mdbridge/render/colormap.hpp"
#include "mdbridge/render/text.hpp"
#include "mdbridge/render/theme.hpp"
#include "mdbridge/render/ticks.hpp"

namespace mdbridge::render {

struct Point {
    double x = 0.0;
    double y = 0.0;
};

// A box x0 <= x1, y0 <= y1 (Matplotlib's Bbox; display pixels or figure fractions by context).
struct Box {
    double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;
    double width() const { return x1 - x0; }
    double height() const { return y1 - y0; }
    static Box fromBounds(double x, double y, double w, double h) { return {x, y, x + w, y + h}; }
    Box translated(double dx, double dy) const { return {x0 + dx, y0 + dy, x1 + dx, y1 + dy}; }
    Box padded(double p) const { return {x0 - p, y0 - p, x1 + p, y1 + p}; }
    // Bbox._is_finite: finite and not a single point (a line counts)
    bool isFinite() const;
    bool operator==(const Box&) const = default;
};
Box unionOf(const std::vector<Box>& boxes);  // Bbox.union (boxes must not be empty)

// Which transform an artist's coordinates go through.
enum class Coords : std::uint8_t {
    Data,    // ax.transData
    Axes,    // ax.transAxes (0..1 across the axes)
    Figure,  // fig.transFigure (0..1 across the figure)
    XAxis,   // x data, y axes (axvline)
    YAxis,   // x axes, y data (axhline)
    Pixels,  // display pixels
};

enum class HAlign : std::uint8_t { Left, Center, Right };
enum class VAlign : std::uint8_t { Top, Center, Baseline, CenterBaseline, Bottom };

struct FontSpec {
    std::vector<std::string> family{"sans-serif"};
    std::string style = "normal";  // normal, italic, oblique
    int weight = 400;              // 400 normal, 700 bold
    double size = 10.0;            // points
};

// What every artist has (Matplotlib's Artist): draw order, visibility, alpha, the frame loop's
// "animated" flag (drawn per frame, not into the static background), the legend label.
struct ArtistBase {
    double zorder = 0.0;
    bool visible = true;
    std::optional<double> alpha;
    bool animated = false;
    std::string label;
    bool clipOn = true;   // clipped to the axes box
    bool inLayout = true;
};

struct Text : ArtistBase {
    std::string text;
    Point xy;
    Coords coords = Coords::Data;
    HAlign ha = HAlign::Left;
    VAlign va = VAlign::Baseline;
    Rgba color{0.0, 0.0, 0.0, 1.0};
    FontSpec font;
    double rotation = 0.0;        // degrees, counter-clockwise
    bool rotationAnchor = false;  // rotation_mode "anchor" (align, then rotate) vs "default"
    bool wrap = false;
    Text() {
        zorder = 3.0;
        clipOn = false;
    }
};

// Line2D's dash pattern: unscaled (points at line width 1; Matplotlib scales by the width).
struct Dashes {
    double offset = 0.0;
    std::vector<double> pattern;  // on, off, on, off ... points; empty: solid
};
enum class LineStyle : std::uint8_t { Solid, Dashed, DashDot, Dotted, None, Custom };
enum class Marker : std::uint8_t { None, Circle, Square };
enum class CapStyle : std::uint8_t { Butt, Round, Projecting };

struct Line : ArtistBase {
    std::vector<double> x, y;
    Coords coords = Coords::Data;
    Rgba color{0.0, 0.0, 0.0, 1.0};
    double width = 1.5;  // points
    LineStyle style = LineStyle::Solid;
    Dashes dashes;       // LineStyle::Custom (and the named styles' patterns)
    Marker marker = Marker::None;
    double markerSize = 6.0;  // points (diameter)
    std::optional<Rgba> markerFace;  // default: the line colour
    std::optional<Rgba> markerEdge;
    double markerEdgeWidth = 1.0;
    CapStyle capStyle = CapStyle::Projecting;
    Line() { zorder = 2.0; }
    // The dash pattern at width 1 (rcParams lines.*_pattern for the named styles), points; empty
    // for solid lines.
    Dashes unscaledDashes() const;
    // The dash pattern drawn (scaled by the width, lines.scale_dashes), points.
    Dashes scaledDashes() const;
    bool drawsLine() const { return style != LineStyle::None; }
};

// A PathCollection of markers (Axes.scatter): sizes are areas in points^2.
struct Scatter : ArtistBase {
    std::vector<Point> offsets;
    std::vector<double> sizes;
    std::vector<Rgba> faces;  // per point, cycled (one entry: all points)
    std::vector<Rgba> edges;
    std::vector<double> lineWidths{1.0};
    Marker marker = Marker::Circle;
    Scatter() { zorder = 1.0; }
};

// A LineCollection.
struct Segments : ArtistBase {
    std::vector<std::vector<Point>> segments;
    std::vector<Rgba> colors;      // cycled
    std::vector<double> widths{1.5};  // points, cycled
    CapStyle capStyle = CapStyle::Butt;
    Segments() { zorder = 1.0; }
};

// An AxesImage: a 2-D array shown through a colormap and Normalize(vmin, vmax), or RGBA pixels.
struct Image : ArtistBase {
    int rows = 0, cols = 0;
    std::vector<double> values;  // rows x cols, row 0 first; NaN: masked (the colormap's bad colour)
    std::optional<Colormap> cmap;
    double vmin = 0.0, vmax = 1.0;
    std::vector<std::uint8_t> rgba;  // rows x cols x 4 when there is no colormap
    std::array<double, 4> extent{-0.5, 0.5, 0.5, -0.5};  // left, right, bottom, top (data)
    bool originUpper = true;
    Image() { zorder = 0.0; }
    // The RGBA bytes shown (cmap(norm(values), bytes=True) or rgba), rows x cols x 4, row 0 first.
    std::vector<std::uint8_t> pixels() const;
};

struct Rectangle : ArtistBase {
    double x = 0, y = 0, w = 0, h = 0;
    Coords coords = Coords::Data;
    std::optional<Rgba> face;  // nullopt: not filled
    std::optional<Rgba> edge;
    double lineWidth = 1.0;
    Rectangle() { zorder = 1.0; }
};

using Artist = std::variant<Line, Scatter, Segments, Image, Text, Rectangle>;

// Matplotlib's legend location codes ("best" 0, "upper right" 1 ... "center" 10).
enum class LegendLoc : int {
    Best = 0, UpperRight, UpperLeft, LowerLeft, LowerRight, Right, CenterLeft, CenterRight, LowerCenter,
    UpperCenter, Center
};

struct LegendEntry {
    std::string text;
    Line handle;  // the proxy line drawn in the handle box (style of the original artist)
};

// Axes.legend with Matplotlib's rcParams (handlelength 2, handleheight 0.7, handletextpad 0.8,
// labelspacing 0.5, borderpad 0.4, borderaxespad 0.5, fancybox), styled by style_legend.
struct Legend : ArtistBase {
    LegendLoc loc = LegendLoc::Best;
    std::vector<LegendEntry> entries;
    bool frameOn = true;
    Rgba face{1.0, 1.0, 1.0, 0.8};
    Rgba edge{0.8, 0.8, 0.8, 0.8};
    double edgeWidth = 1.0;  // points
    Rgba textColor{0.0, 0.0, 0.0, 1.0};
    double textSize = 10.0;  // points (style_legend's size)
    double fontSize = 10.0;  // the size the legend was made with: its paddings scale with it
    Legend() { zorder = 5.0; }
};
// style_legend(legend, theme, font_size)
void styleLegend(Legend& legend, const LegendStyle& style);

struct TickStyle {
    bool tick1 = true, tick2 = false, label1 = true, label2 = false;
    Rgba color{0.0, 0.0, 0.0, 1.0};
    Rgba labelColor{0.0, 0.0, 0.0, 1.0};
    double labelSize = 10.0;  // points
    double length = 3.5;      // points
    double width = 0.8;
    double pad = 3.5;  // between tick and label (xtick.major.pad); the label sits length + pad away
    double rotation = 0.0;
    bool gridOn = false;
    Rgba gridColor{0.6901960784313725, 0.6901960784313725, 0.6901960784313725, 1.0};  // "#b0b0b0"
    double gridWidth = 0.8;
};

struct Locator {
    enum class Kind : std::uint8_t { MaxN, Fixed };
    Kind kind = Kind::MaxN;
    MaxNLocatorParams params{std::nullopt, {1.0, 2.0, 2.5, 5.0, 10.0}};  // AutoLocator
    std::vector<double> locs;  // Fixed
    std::optional<int> fixedBins;
    bool isAuto() const { return kind == Kind::MaxN && !params.nbins; }
};

struct Formatter {
    enum class Kind : std::uint8_t { Scalar, Fixed };
    Kind kind = Kind::Scalar;
    ScalarFormatterOptions options;
    std::vector<std::string> seq;  // Fixed: labels in tick order
};

// One axis' ticks as drawn: location, label and the label boxes (display pixels).
struct DrawnTick {
    double loc = 0.0;
    std::string label;
    Box box1, box2;  // label1 / label2 extents (empty when not shown)
};

struct Axis {
    AxisDirection direction = AxisDirection::X;
    bool visible = true;
    Text label;
    bool labelOpposite = false;  // label on the top / right
    double labelPad = 4.0;       // points
    Locator locator;
    Formatter formatter;
    TickStyle ticks;
    Range limits{0.0, 1.0};
    bool autoscale = true;
    double margin = 0.05;
    Text offsetText;

    // results of the last draw
    std::vector<DrawnTick> drawn;
    int tickSpace = 0;
    std::optional<ScalarFormatter> scalar;  // the formatter's state (offset, format)
};

struct Spine {
    bool visible = true;
    Rgba color{0.0, 0.0, 0.0, 1.0};
    double width = 0.8;
};
enum SpineSide : std::uint8_t { kLeft = 0, kBottom, kRight, kTop };

struct SubplotSpec {
    int row0 = 0, row1 = 1, col0 = 0, col1 = 1;
};

// A colour bar on the right of its parent axes (fig.colorbar(mappable, ax=parent, fraction, pad)):
// constrained layout places it in the parent cell's right margin, aspect 20.
struct ColorbarInfo {
    std::size_t parent = 0;
    std::optional<std::size_t> mappable;  // artist index in the parent axes
    double fraction = 0.15, pad = 0.05, aspect = 20.0, shrink = 1.0;
    Point anchor{0.0, 0.5};
    std::optional<Colormap> cmap;
    double vmin = 0.0, vmax = 1.0;
    Spine outline;
};

struct Axes {
    std::optional<SubplotSpec> spec;  // grid cells (none for a colour bar)
    std::optional<ColorbarInfo> colorbar;
    std::optional<std::size_t> twinOf;  // twinx of this axes index
    Rgba face{1.0, 1.0, 1.0, 1.0};
    bool patchVisible = true;
    std::array<Spine, 4> spines;  // left, bottom, right, top
    Axis xaxis, yaxis;
    std::array<Text, 3> titles;  // left, center, right
    std::vector<Artist> artists;
    std::optional<Legend> legend;
    double zorder = 0.0;
    bool visible = true;
    bool axison = true;
    AxisBelow axisBelow = AxisBelow::Line;
    std::optional<double> boxAspect;  // height / width of the drawn box
    Point anchor{0.5, 0.5};
    // tests and dumps: which panel and role ("main", "twin", "colorbar") the axes belongs to
    int panel = -1;
    std::string role;

    // positions in figure fractions: as laid out ("original") and as drawn ("active", after the
    // box aspect)
    Box position{0.125, 0.11, 0.9, 0.88};
    Box active{0.125, 0.11, 0.9, 0.88};

    Axes();
    void setXLim(double lo, double hi);
    void setYLim(double lo, double hi);
    void setXTicks(std::vector<double> locs, std::optional<std::vector<std::string>> labels = std::nullopt);
    void setYTicks(std::vector<double> locs, std::optional<std::vector<std::string>> labels = std::nullopt);
    // Axes.set_title(text, loc=..)
    Text& setTitle(std::string text, HAlign loc = HAlign::Center);
    // Adds an artist and, as Matplotlib's add_line / add_collection / imshow do, its data to the
    // data limits autoscaling uses (changing an artist's data later does not change them: relim()).
    std::size_t add(Artist artist);
    // Axes.relim(): the data limits from the artists as they are now
    void relim();
    // Axes.cla() for the "No data" panels: artists, legend, ticks and limits as new
    void clear();

    // the data limits (Axes.dataLim, per axis) and sticky edges (images' extents)
    struct DataLimits {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        std::vector<double> sticky;
        void add(double v);
        bool empty() const { return !(lo <= hi); }
    };
    DataLimits dataX, dataY;
    void updateDataLimits(const Artist& artist);
};

// style_axes(ax, theme, grid, font_size) on an axes
void styleAxes(Axes& ax, const AxesStyle& style);
// Axes.tick_params(axis, colors / labelcolor / labelsize / width / length) on one axis
struct TickParams {
    std::optional<Rgba> color, labelColor;
    std::optional<double> labelSize, width, length, rotation;
};
void tickParams(Axis& axis, const TickParams& params);

struct GridSpec {
    int rows = 1, cols = 1;
    std::vector<double> heightRatios;  // empty: equal rows
    std::vector<double> widthRatios;
};

// fig.set_layout_engine("constrained", h_pad, w_pad, hspace, wspace): pads in inches, spaces in
// fractions; MDATK's figures use 0.08, 0.08, 0.06, 0.04.
struct LayoutParams {
    bool constrained = true;
    double hPad = 0.04167, wPad = 0.04167;
    double hspace = 0.02, wspace = 0.02;
};

class Figure {
public:
    // An exactly width x height pixel figure at dpi (Matplotlib: figsize = w / dpi inches).
    Figure(int widthPx, int heightPx, double dpi);

    int widthPx() const { return m_width; }
    int heightPx() const { return m_height; }
    double dpi() const { return m_dpi; }
    double widthInches() const { return m_width / m_dpi; }
    double heightInches() const { return m_height / m_dpi; }
    double pointsToPixels(double pt) const { return pt * m_dpi / 72.0; }

    Rgba facecolor{1.0, 1.0, 1.0, 1.0};
    GridSpec grid;
    LayoutParams layout;
    std::optional<Text> suptitle;
    std::vector<Text> texts;  // fig.text (figure coordinates by default)
    std::vector<Axes> axes;   // fig.axes order: drawn and dumped in this order

    std::size_t addSubplot(SubplotSpec spec);
    // ax.twinx(): the same cell, y axis on the right, x axis hidden, no background patch
    std::size_t twinx(std::size_t ax);
    // fig.colorbar(mappable, ax=parent, fraction, pad): a vertical colour bar right of the parent
    std::size_t colorbar(std::size_t parent, std::optional<std::size_t> mappable, Colormap cmap, double vmin,
                         double vmax, double fraction = 0.15, double pad = 0.05);
    // fig.suptitle(text, color, fontsize)
    Text& setSuptitle(std::string text, Rgba color, double size);
    // fig.text(x, y, text, ...)
    Text& addText(double x, double y, std::string text);

    // Lays the figure out (unless frozen) and computes what a draw computes: limits, ticks and tick
    // labels, label, title and legend positions, text extents. Call after every change.
    void draw(const FontSet& fonts);
    // fig.set_layout_engine("none") after the first draw: positions stay as they are
    void freeze() { m_frozen = true; }
    bool frozen() const { return m_frozen; }
    // warnings of the last draw (e.g. the layout collapsed: Matplotlib's warn_collapsed)
    const std::vector<std::string>& warnings() const { return m_warnings; }

    // Geometry after draw(): the drawn box of an axes in pixels, a text's extent, a legend's box
    // and handle / text offsets, axis label positions.
    Box axesBox(std::size_t ax) const;
    Box figureBox() const { return {0.0, 0.0, static_cast<double>(m_width), static_cast<double>(m_height)}; }

    // ---- transforms (display pixels)
    Point toPixels(const Axes& ax, Coords coords, Point p) const;

    // ---- the internals the layout, dump and display list share
    struct TextGeometry {
        Box box;                 // the rotated text's axis-aligned extent
        std::vector<std::pair<std::string, Point>> lines;  // each line and its baseline origin
        std::vector<TextLayout> layouts;
        FontFace face = FontFace::Sans;
        double pixelSize = 0.0;
        double descent = 0.0;    // of the last line
    };
    // The text's geometry at display position `at` (Text._get_layout + get_window_extent).
    TextGeometry textGeometry(const Text& text, Point at, const FontSet& fonts) const;
    // Text._get_wrapped_text: a wrap=True text broken into lines that fit between its position and
    // the figure's edge (the text itself otherwise)
    std::string wrappedText(const Text& text, Point at, const FontSet& fonts) const;
    // get_window_extent of a text placed by its own coordinates in `ax` (or the figure)
    Box textExtent(const Text& text, const Axes* ax, const FontSet& fonts) const;

    struct LegendGeometry {
        Box box;  // the legend box (frame bounds)
        struct Item {
            Point handleOrigin;  // the handle box' origin (lower left + descent), pixels
            double handleWidth = 0, handleHeight = 0, handleDescent = 0;
            Point textOrigin;    // the text's anchor (left, baseline)
        };
        std::vector<Item> items;
    };
    LegendGeometry legendGeometry(std::size_t ax) const;  // after draw()
    // Legend._find_best_position: the location code ("upper right" 1 ... "center" 10) where a legend
    // box of width x height pixels covers the fewest data points, line segments, rectangles and
    // texts of the axes (ties: the lowest code). After draw().
    LegendLoc bestLegendLocation(std::size_t ax, double width, double height, const FontSet& fonts) const;
    // Where a legend box of width x height pixels goes for a location code (_get_anchored_bbox):
    // its lower left corner.
    Point legendAnchor(std::size_t ax, LegendLoc loc, double width, double height, double fontSize) const;

    // the axis label / title / offset text positions after draw (display; Text::xy is not changed)
    Point axisLabelPosition(std::size_t ax, AxisDirection d) const;
    Point offsetTextPosition(std::size_t ax, AxisDirection d) const;
    Point titlePosition(std::size_t ax, int which) const;
    Point suptitlePosition() const;

    // the tight box (Axes.get_tightbbox) of an axes after draw, display pixels
    std::optional<Box> tightBox(std::size_t ax) const;

private:
    void autoscaleAxes(Axes& ax);
    void computeTicks(std::size_t a, const FontSet& fonts);
    void updateAxisLabels(std::size_t a);
    void updateTitles(std::size_t a, const FontSet& fonts);
    void applyAspect(Axes& ax);
    std::optional<Box> computeTightBox(std::size_t a, const FontSet& fonts, bool forLayout);
    std::optional<Box> axisTightBox(std::size_t a, AxisDirection d, const FontSet& fonts, bool forLayout);
    void placeLegend(std::size_t a, const FontSet& fonts);
    void constrainedLayout(const FontSet& fonts);
    void initialPositions();
    Box spineExtent(std::size_t a, SpineSide side) const;

    int m_width;
    int m_height;
    double m_dpi;
    bool m_frozen = false;
    bool m_positioned = false;
    std::vector<std::string> m_warnings;

    struct Derived {
        Point xLabel, yLabel, xOffset, yOffset;
        std::array<Point, 3> titles;
        std::optional<Box> tight;
        LegendGeometry legend;
    };
    std::vector<Derived> m_derived;
};

}  // namespace mdbridge::render
