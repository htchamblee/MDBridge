#pragma once
// MDATK's figure themes and chart chrome (mdatk/render/theme.py; MDATK's own code, ported): two
// themes (light / dark) of seven colour tokens, the 8-colour categorical palette assigned in fixed
// order, the sequential blue ramp and the ordinal ramp for class matrices, and the styles that
// style_axes / style_legend apply. Qt-free; the figure code (W1-A) applies the styles.

#include <cstdint>
#include <span>
#include <string_view>

#include "mdbridge/render/color.hpp"
#include "mdbridge/render/colormap.hpp"

namespace mdbridge::render {

struct Theme {
    std::string_view name;     // "light" or "dark"
    std::string_view surface;  // axes face, legend face, the ramps' zero
    std::string_view page;     // figure face
    std::string_view ink;      // titles
    std::string_view ink2;     // tick labels, axis labels, legend texts
    std::string_view muted;    // ticks, "No data" texts, watermark
    std::string_view grid;     // grid lines, legend edge, heatmap NaN (bad) colour
    std::string_view axis;     // spines

    // The categorical palette ("#2a78d6", ... 8 colours of the theme).
    std::span<const std::string_view> series() const;
    // series()[k % 8] (Python's modulo: never cycled to new colours, negative k wraps).
    std::string_view seriesColor(std::int64_t k) const;
    // LinearSegmentedColormap.from_list("mdatk_seq_<name>", [surface] + 13 blues (reversed for dark)).
    Colormap sequential() const;
    // ListedColormap for n ordered classes; class 0 recedes ("#e1e0d9" / "#2c2c2a"); n <= 1: the
    // surface alone. Up to 5 classes as MDATK; more classes get n - 1 distinct steps along the ramp
    // (MDATK repeated its 4 ramp colours, so classes shared colours; docs/DIFFERENCES.md).
    Colormap ordinal(int n) const;
};

const Theme& lightTheme();
const Theme& darkTheme();
// get_theme(name): the dark theme for "dark", the light one for anything else (as MDATK).
const Theme& theme(std::string_view name);
// The names visual.theme means as written ("light", "dark"); check() warns about others.
bool isThemeName(std::string_view name);

// The 13 colours of MDATK's SEQUENTIAL_BLUE ramp.
std::span<const std::string_view> sequentialBlue();

enum class AxisBelow { Line, True };  // Axes.set_axisbelow: rcParams' "line" or True

struct SpineStyle {
    bool visible = true;
    Rgba color;
    double width = 0.8;
};

// The axes after style_axes(ax, theme, grid, font_size) (Matplotlib defaults where it sets nothing).
struct AxesStyle {
    Rgba face;
    SpineStyle left, bottom, top, right;
    Rgba tickColor;        // tick marks
    Rgba tickLabelColor;
    double tickLabelSize = 0.0;  // points
    double tickWidth = 0.0;
    double tickLength = 0.0;
    Rgba axisLabelColor;
    double axisLabelSize = 0.0;
    Rgba titleColor;
    double titleSize = 0.0;
    bool grid = false;  // grid lines on both axes
    Rgba gridColor;
    double gridWidth = 0.0;
    std::string_view gridLineStyle;  // "-"
    AxisBelow axisBelow = AxisBelow::Line;
};
AxesStyle axesStyle(const Theme& theme, bool grid = true, double fontSize = 10.0);

// The legend after style_legend(legend, theme, font_size).
struct LegendStyle {
    Rgba face;  // surface, alpha 0.8 (legend.framealpha)
    Rgba edge;  // grid, alpha 0.8
    double edgeWidth = 0.0;
    Rgba textColor;
    double textSize = 0.0;
};
LegendStyle legendStyle(const Theme& theme, double fontSize = 10.0);

}  // namespace mdbridge::render
