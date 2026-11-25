#pragma once
// Colour specifications as Matplotlib 3.11 reads them (matplotlib.colors.to_rgba / to_hex), for the
// colours of MDATK's project files (visual.base_color, glow_color, group colours, chart colour ...)
// and of the figures. Qt-free.
//
// Port of matplotlib.colors.to_rgba / _to_rgba_no_colorcycle / to_hex (Matplotlib 3.11.2,
// Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved; Matplotlib licence,
// packaging/licenses/Matplotlib-LICENSE.txt). Changes from Matplotlib: C++ API for strings and
// RGB(A) tuples only; a tuple with a NaN component is rejected (docs/DIFFERENCES.md).

#include <optional>
#include <span>
#include <string>
#include <string_view>

namespace mdbridge::render {

// RGBA, each channel in [0, 1] (what to_rgba returns).
struct Rgba {
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
    double a = 1.0;
    friend bool operator==(const Rgba&, const Rgba&) = default;
};

// to_rgba(color, alpha) of a str, Matplotlib's grammar:
//   "none" in any case: (0, 0, 0, 0), whatever alpha;
//   "C0", "C1", ... (C and ASCII digits, case-sensitive): colour n % 10 of the default cycle (tab10);
//   a name (base "r", "g" ..., CSS4 "red", "tab:blue", "xkcd:sky blue", "grey" / "gray" spellings),
//     exact or, for names longer than one character, after Python's str.lower();
//   "#rgb", "#rgba", "#rrggbb", "#rrggbbaa" (hex digits in any case);
//   a grey level: Python float() of the string, 0 <= v <= 1 ("0.5", " .25", "1e-1").
// alpha, when given, replaces the colour's alpha (must lie in [0, 1]). Throws core::ValueError with
// Matplotlib's message: "Invalid RGBA argument: 'x'", "Invalid hex color specifier: '#12'",
// "Invalid string grayscale value '2'. Value must be within 0-1 range",
// "'alpha' must be between 0 and 1, inclusive".
Rgba toRgba(std::string_view color, std::optional<double> alpha = std::nullopt);

// to_rgba of an RGB or RGBA tuple: length 3 or 4 ("RGBA sequence should have length 3 or 4"),
// components in [0, 1] ("RGBA values should be within 0-1 range"; NaN too, unlike Matplotlib).
Rgba toRgba(std::span<const double> components, std::optional<double> alpha = std::nullopt);

// The message toRgba would throw, or nothing for a valid colour (check() of project files, RN-12).
std::optional<std::string> colorProblem(std::string_view color);
// matplotlib.colors.is_color_like(color) for a str.
bool isColorLike(std::string_view color);

// to_hex(color, keep_alpha): "#rrggbb" (or "#rrggbbaa"), each channel round(v * 255) (ties to even).
std::string toHex(const Rgba& color, bool keepAlpha = false);

// The colours of "C0" ... "C9" (rcParamsDefault["axes.prop_cycle"], Matplotlib's tab10).
std::span<const Rgba> defaultColorCycle();

// Every name to_rgba knows, sorted by bytes (1163 in Matplotlib 3.11.2).
std::span<const std::string_view> namedColorNames();

}  // namespace mdbridge::render
