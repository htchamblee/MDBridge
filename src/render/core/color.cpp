// Port of matplotlib.colors.to_rgba / _to_rgba_no_colorcycle / to_hex (Matplotlib 3.11.2).
// Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved. Used under the Matplotlib
// licence (packaging/licenses/Matplotlib-LICENSE.txt). Changes: rewritten in C++ for str and
// RGB(A) tuple inputs; tuples with NaN components are rejected.

#include "mdbridge/render/color.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mpl_tables.hpp"

namespace mdbridge::render {

namespace {

std::string repr(std::string_view text) { return core::pythonStringRepr(std::string(text)); }

[[noreturn]] void fail(const std::string& message) { throw core::ValueError(message); }

void checkAlpha(const std::optional<double>& alpha) {
    if (alpha && !(0.0 <= *alpha && *alpha <= 1.0)) {
        fail("'alpha' must be between 0 and 1, inclusive");
    }
}

std::size_t codePoints(std::string_view text) {
    return static_cast<std::size_t>(std::count_if(text.begin(), text.end(), [](char c) {
        return (static_cast<unsigned char>(c) & 0xC0) != 0x80;
    }));
}

const detail::NamedColor* findNamed(std::string_view name) {
    const auto* begin = std::begin(detail::kNamedColors);
    const auto* end = std::end(detail::kNamedColors);
    const auto* it = std::lower_bound(begin, end, name, [](const detail::NamedColor& entry, std::string_view key) {
        return entry.name < key;  // byte order, as the generator sorts
    });
    return it != end && it->name == name ? it : nullptr;
}

bool isHexDigit(char c) { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'); }

int hexValue(char c) {
    if (c >= '0' && c <= '9') {
        return c - '0';
    }
    if (c >= 'a' && c <= 'f') {
        return c - 'a' + 10;
    }
    return c - 'A' + 10;
}

// \AC[0-9]+\Z: an entry of the colour cycle
bool isNthColor(std::string_view c) {
    return c.size() >= 2 && c[0] == 'C' &&
           std::all_of(c.begin() + 1, c.end(), [](char d) { return d >= '0' && d <= '9'; });
}

Rgba fromRgb(const detail::Rgb& rgb, const std::optional<double>& alpha) {
    return {rgb.r, rgb.g, rgb.b, alpha ? *alpha : 1.0};
}

// _to_rgba_no_colorcycle for a str
Rgba stringToRgba(std::string_view c, const std::optional<double>& alpha) {
    checkAlpha(alpha);
    const std::string lower = core::pythonLower(c);
    if (lower == "none") {
        return {0.0, 0.0, 0.0, 0.0};
    }
    const detail::NamedColor* named = findNamed(c);
    if (named == nullptr && codePoints(c) != 1) {
        named = findNamed(lower);
    }
    if (named != nullptr) {
        return fromRgb(named->rgb, alpha);
    }
    if (c.size() >= 2 && c[0] == '#' && std::all_of(c.begin() + 1, c.end(), isHexDigit)) {
        std::vector<double> channels;
        if (c.size() == 7 || c.size() == 9) {  // #rrggbb, #rrggbbaa
            for (std::size_t i = 1; i < c.size(); i += 2) {
                channels.push_back((hexValue(c[i]) * 16 + hexValue(c[i + 1])) / 255.0);
            }
        } else if (c.size() == 4 || c.size() == 5) {  // #rgb, #rgba
            for (std::size_t i = 1; i < c.size(); ++i) {
                channels.push_back(hexValue(c[i]) / 15.0);
            }
        } else {
            fail("Invalid hex color specifier: " + repr(c));
        }
        const double a = channels.size() == 4 ? channels[3] : 1.0;
        return {channels[0], channels[1], channels[2], alpha ? *alpha : a};
    }
    if (const auto grey = core::tryPythonFloat(c)) {
        const double v = *grey;
        if (!(0.0 <= v && v <= 1.0)) {
            fail("Invalid string grayscale value " + repr(c) + ". Value must be within 0-1 range");
        }
        return {v, v, v, alpha ? *alpha : 1.0};
    }
    fail("Invalid RGBA argument: " + repr(c));
}

}  // namespace

Rgba toRgba(std::string_view color, std::optional<double> alpha) {
    if (isNthColor(color)) {
        // int(color[1:]) % len(cycle), digit by digit (any length)
        const auto n = static_cast<std::uint64_t>(std::size(detail::kDefaultCycle));
        std::uint64_t index = 0;
        for (char d : color.substr(1)) {
            index = (index * 10 + static_cast<std::uint64_t>(d - '0')) % n;
        }
        checkAlpha(alpha);
        return fromRgb(detail::kDefaultCycle[index], alpha);
    }
    return stringToRgba(color, alpha);
}

Rgba toRgba(std::span<const double> components, std::optional<double> alpha) {
    checkAlpha(alpha);
    if (components.size() != 3 && components.size() != 4) {
        fail("RGBA sequence should have length 3 or 4");
    }
    Rgba out{components[0], components[1], components[2], components.size() == 4 ? components[3] : 1.0};
    if (alpha) {
        out.a = *alpha;
    }
    for (double v : {out.r, out.g, out.b, out.a}) {
        if (!(v >= 0.0 && v <= 1.0)) {  // NaN too (Matplotlib lets it through)
            fail("RGBA values should be within 0-1 range");
        }
    }
    return out;
}

std::optional<std::string> colorProblem(std::string_view color) {
    try {
        toRgba(color);
    } catch (const core::ValueError& error) {
        return std::string(error.what());
    }
    return std::nullopt;
}

bool isColorLike(std::string_view color) { return !colorProblem(color); }

std::string toHex(const Rgba& color, bool keepAlpha) {
    static constexpr char kDigits[] = "0123456789abcdef";
    std::string out = "#";
    const double channels[] = {color.r, color.g, color.b, color.a};
    for (std::size_t i = 0; i < (keepAlpha ? 4u : 3u); ++i) {
        // round(v * 255): ties to even, as std::nearbyint in the default rounding mode
        const auto v = static_cast<int>(std::nearbyint(channels[i] * 255.0));
        const int clamped = std::clamp(v, 0, 255);
        out += kDigits[clamped / 16];
        out += kDigits[clamped % 16];
    }
    return out;
}

std::span<const Rgba> defaultColorCycle() {
    static const std::vector<Rgba> cycle = [] {
        std::vector<Rgba> out;
        for (const auto& rgb : detail::kDefaultCycle) {
            out.push_back({rgb.r, rgb.g, rgb.b, 1.0});
        }
        return out;
    }();
    return cycle;
}

std::span<const std::string_view> namedColorNames() {
    static const std::vector<std::string_view> names = [] {
        std::vector<std::string_view> out;
        for (const auto& entry : detail::kNamedColors) {
            out.push_back(entry.name);
        }
        return out;
    }();
    return names;
}

}  // namespace mdbridge::render
