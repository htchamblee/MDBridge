// MDATK's themes and chart chrome (mdatk/render/theme.py, MDATK v0.0.4.1; same owner, ported).

#include "mdbridge/render/theme.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace mdbridge::render {

namespace {

constexpr std::string_view kLightSeries[] = {"#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                                             "#e87ba4", "#008300", "#6250d6", "#e34948"};
constexpr std::string_view kDarkSeries[] = {"#3987e5", "#d95926", "#199e70", "#c98500",
                                            "#d55181", "#008300", "#9085e9", "#e66767"};
constexpr std::string_view kSequentialBlue[] = {"#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
                                                "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
                                                "#184f95", "#104281", "#0d366b"};
constexpr std::string_view kLightOrdinal[] = {"#86b6ef", "#3987e5", "#1c5cab", "#0d366b"};
constexpr std::string_view kDarkOrdinal[] = {"#184f95", "#3987e5", "#86b6ef", "#cde2fb"};

constexpr Theme kLight{"light", "#fcfcfb", "#f9f9f7", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"};
constexpr Theme kDark{"dark", "#1a1a19", "#0d0d0d", "#f0efec", "#c3c2b7", "#898781", "#2c2c2a", "#383835"};

bool isDark(const Theme& t) { return t.name == "dark"; }

}  // namespace

std::span<const std::string_view> Theme::series() const {
    using Colors = std::span<const std::string_view>;
    return isDark(*this) ? Colors(kDarkSeries) : Colors(kLightSeries);
}

std::string_view Theme::seriesColor(std::int64_t k) const {
    const auto colors = series();
    const auto n = static_cast<std::int64_t>(colors.size());
    return colors[static_cast<std::size_t>(((k % n) + n) % n)];
}

Colormap Theme::sequential() const {
    std::vector<Rgba> steps{toRgba(surface)};
    if (isDark(*this)) {
        for (auto it = std::rbegin(kSequentialBlue); it != std::rend(kSequentialBlue); ++it) {
            steps.push_back(toRgba(*it));
        }
    } else {
        for (std::string_view c : kSequentialBlue) {
            steps.push_back(toRgba(c));
        }
    }
    return Colormap::fromList("mdatk_seq_" + std::string(name), steps);
}

Colormap Theme::ordinal(int n) const {
    if (n <= 1) {
        return Colormap::listed("unnamed", {toRgba(surface)});
    }
    using Colors = std::span<const std::string_view>;
    const Colors ramp = isDark(*this) ? Colors(kDarkOrdinal) : Colors(kLightOrdinal);
    std::vector<Rgba> colors{toRgba(isDark(*this) ? "#2c2c2a" : "#e1e0d9")};
    const auto steps = static_cast<std::size_t>(n - 1);
    if (steps <= ramp.size()) {
        if (n > 2) {
            for (int k = 0; k < n - 1; ++k) {
                // round(k * (len(ramp) - 1) / max(1, n - 2)): Python's round, ties to even
                const double at = static_cast<double>(k) * static_cast<double>(ramp.size() - 1) / std::max(1, n - 2);
                colors.push_back(toRgba(ramp[static_cast<std::size_t>(std::nearbyint(at))]));
            }
        } else {
            colors.push_back(toRgba(ramp[1]));
        }
    } else {
        // more classes than ramp colours: n - 1 evenly spaced steps along the ramp (MDATK gave the
        // four ramp colours to n - 1 classes, so several classes shared one; DIFFERENCES.md RN-18)
        std::vector<Rgba> stops;
        for (std::string_view c : ramp) {
            stops.push_back(toRgba(c));
        }
        const Colormap along = Colormap::fromList("ordinal", stops, n - 1);
        for (std::size_t i = 0; i < steps; ++i) {
            colors.push_back(along.lut()[i]);
        }
    }
    return Colormap::listed("unnamed", std::move(colors));
}

const Theme& lightTheme() { return kLight; }
const Theme& darkTheme() { return kDark; }
const Theme& theme(std::string_view name) { return name == "dark" ? kDark : kLight; }
bool isThemeName(std::string_view name) { return name == "light" || name == "dark"; }

std::span<const std::string_view> sequentialBlue() { return kSequentialBlue; }

AxesStyle axesStyle(const Theme& t, bool grid, double fontSize) {
    AxesStyle s;
    s.face = toRgba(t.surface);
    s.left = s.bottom = {true, toRgba(t.axis), 0.8};
    s.top = s.right = {false, toRgba("black"), 0.8};  // rcParams axes.edgecolor / axes.linewidth
    s.tickColor = toRgba(t.muted);
    s.tickLabelColor = toRgba(t.ink2);
    s.tickLabelSize = fontSize * 0.85;
    s.tickWidth = 0.6;
    s.tickLength = 3.0;
    s.axisLabelColor = toRgba(t.ink2);
    s.axisLabelSize = fontSize * 0.9;
    s.titleColor = toRgba(t.ink);
    s.titleSize = fontSize;
    s.grid = grid;
    if (grid) {
        s.gridColor = toRgba(t.grid);
        s.gridWidth = 0.6;
        s.gridLineStyle = "-";
        s.axisBelow = AxisBelow::True;
    } else {  // rcParams grid.color, grid.linewidth, grid.linestyle, axes.axisbelow
        s.gridColor = toRgba("#b0b0b0");
        s.gridWidth = 0.8;
        s.gridLineStyle = "-";
        s.axisBelow = AxisBelow::Line;
    }
    return s;
}

LegendStyle legendStyle(const Theme& t, double fontSize) {
    // the frame keeps rcParams' legend.framealpha (0.8) over the colours style_legend sets
    return {toRgba(t.surface, 0.8), toRgba(t.grid, 0.8), 0.6, toRgba(t.ink2), fontSize * 0.8};
}

}  // namespace mdbridge::render
