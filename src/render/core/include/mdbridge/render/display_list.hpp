#pragma once
// A drawn figure (figure.hpp) as a list of painting primitives in Matplotlib's draw order, for the
// Qt backends (figure_qt.hpp: QImage, SVG, PDF) and the frame loop. Qt-free.
//
// Coordinates are device pixels of the figure with the origin at the TOP left and y down (paint
// coordinates); widths and sizes are pixels (points already converted at the figure's dpi). Colours
// are straight (not premultiplied) RGBA in [0, 1], artist alphas applied.
//
// The frame loop splits a figure in two layers (MDATK's blitting, mdatk/render/video.py): the static
// background (everything not animated, drawn once) and the animated artists (drawn over the
// background for every frame, in draw order).

#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "mdbridge/render/figure.hpp"

namespace mdbridge::render {

struct DisplayStroke {
    Rgba color;
    double width = 1.0;
    CapStyle cap = CapStyle::Butt;
    bool roundJoin = true;
    std::vector<double> dashes;  // on, off, ... pixels; empty: solid
    double dashOffset = 0.0;
};

struct FillRectItem {
    Box rect;  // x0 <= x1, y0 <= y1 (top-left origin)
    Rgba color;
};

// A polyline (open) or polygon (closed); `fill` for closed shapes.
struct PathItem {
    std::vector<std::vector<Point>> subpaths;  // NaN breaks give several subpaths
    bool closed = false;
    std::optional<Rgba> fill;
    std::optional<DisplayStroke> stroke;
};

struct MarkerItem {
    Marker shape = Marker::Circle;
    std::vector<Point> centers;
    std::vector<double> sizes{6.0};  // diameter / side per marker (cycled), pixels
    std::vector<Rgba> faces;  // per marker (cycled); empty: unfilled
    std::vector<Rgba> edges;  // cycled; empty: no edge
    double edgeWidth = 1.0;
};

struct RoundRectItem {
    Box rect;
    double radius = 0.0;
    std::optional<Rgba> fill;
    std::optional<DisplayStroke> stroke;
};

struct TextItem {
    std::string text;  // one line
    FontFace face = FontFace::Regular;
    double pixelSize = 10.0;  // the em in pixels
    Point origin;             // baseline, left end of the line
    double rotation = 0.0;    // degrees, counter-clockwise on screen
    Rgba color;
    TextLayout layout;  // the glyph run (font.hpp) the raster backend fills
};

struct ImageItem {
    int width = 0, height = 0;
    std::vector<std::uint8_t> rgba;  // straight RGBA, row 0 first
    // where pixel (0, 0)'s outer corner and pixel (w - 1, h - 1)'s outer corner land: the image may be
    // mirrored in x or y (origin / inverted limits)
    Point corner0, corner1;
    bool smooth = false;  // nearest-neighbour unless set
};

struct DisplayItem {
    std::variant<FillRectItem, PathItem, MarkerItem, RoundRectItem, TextItem, ImageItem> item;
    std::optional<Box> clip;  // device pixels; nullopt: the whole figure
    bool animated = false;
};

struct DisplayList {
    int width = 0, height = 0;  // device pixels
    double dpi = 72.0;
    std::vector<DisplayItem> items;
};

enum class DisplayLayer : std::uint8_t {
    All,       // a still image: everything
    Static,    // the frame loop's background: everything not animated
    Animated,  // the frame loop's per-frame layer: the animated artists only
};

// The figure as drawn after Figure::draw(fonts).
DisplayList buildDisplayList(const Figure& figure, const FontSet& fonts, DisplayLayer layer = DisplayLayer::All);

}  // namespace mdbridge::render
