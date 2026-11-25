// Port of matplotlib.colors.Colormap, ListedColormap, LinearSegmentedColormap and
// _create_lookup_table (Matplotlib 3.11.2). Copyright (c) 2012- Matplotlib Development Team; All
// Rights Reserved. Used under the Matplotlib licence (packaging/licenses/Matplotlib-LICENSE.txt).
// numpy.linspace and numpy.searchsorted's left search (NumPy 2.5.5, Copyright (c) 2005-2025, NumPy
// Developers; BSD-3-Clause, packaging/licenses/numpy-LICENSE.txt). Changes: rewritten in C++ for
// scalar values; segment data is validated when the map is made; at least one colour.
//
// Compiled without fused multiply-adds (src/render/core/CMakeLists.txt): the lookup tables are
// rounded as numpy rounds them.

#include "mdbridge/render/colormap.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mpl_tables.hpp"

namespace mdbridge::render {

namespace {

[[noreturn]] void fail(const std::string& message) { throw core::ValueError(message); }

void validate(const std::vector<Colormap::Segment>& rows) {
    if (rows.empty() || rows.front().x != 0.0 || rows.back().x != 1.0) {
        fail("data mapping points must start with x=0 and end with x=1");
    }
    for (std::size_t i = 1; i < rows.size(); ++i) {
        if (rows[i].x - rows[i - 1].x < 0.0) {  // np.diff(x) < 0
            fail("data mapping points must have x in increasing order");
        }
    }
}

double clip01(double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }  // np.clip (NaN stays)

// _create_lookup_table(N, data, gamma) for (x, y0, y1) rows
std::vector<double> lookupTable(int n, const std::vector<Colormap::Segment>& rows, double gamma) {
    std::vector<double> lut(static_cast<std::size_t>(n));
    if (n == 1) {
        lut[0] = clip01(rows.back().y0);  // convention: y = f(x=1)
        return lut;
    }
    const double scale = static_cast<double>(n - 1);
    std::vector<double> x(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        x[i] = rows[i].x * scale;
    }
    std::vector<double> xind = linspace(0.0, 1.0, n);
    for (double& v : xind) {
        v = scale * (gamma == 1.0 ? v : std::pow(v, gamma));
    }
    lut.front() = rows.front().y1;
    lut.back() = rows.back().y0;
    for (std::size_t j = 1; j + 1 < xind.size(); ++j) {
        // np.searchsorted(x, xind) (side "left"): the first x >= xind
        const auto ind = static_cast<std::size_t>(std::lower_bound(x.begin(), x.end(), xind[j]) - x.begin());
        const double distance = (xind[j] - x[ind - 1]) / (x[ind] - x[ind - 1]);
        lut[j] = distance * (rows[ind].y0 - rows[ind - 1].y1) + rows[ind - 1].y1;
    }
    for (double& v : lut) {
        v = clip01(v);
    }
    return lut;
}

std::uint8_t toByte(double v) { return static_cast<std::uint8_t>(v * 255.0); }  // astype(uint8): truncation

std::string joinedNames() {
    std::string out;
    for (std::string_view name : bundledColormapNames()) {
        out += (out.empty() ? "" : ", ") + std::string(name);
    }
    return out;
}

std::string lowerAscii(std::string_view text) {
    std::string out(text);
    std::transform(out.begin(), out.end(), out.begin(),
                   [](char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; });
    return out;
}

}  // namespace

std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> y;
    if (num <= 0) {
        return y;
    }
    y.resize(static_cast<std::size_t>(num));
    const int div = num - 1;
    const double delta = stop - start;
    if (div > 0) {
        const double step = delta / div;
        for (int i = 0; i < num; ++i) {
            if (step == 0.0) {  // numpy's handling of denormal steps
                y[static_cast<std::size_t>(i)] = static_cast<double>(i) / div * delta;
            } else {
                y[static_cast<std::size_t>(i)] = static_cast<double>(i) * step;
            }
        }
    } else {
        y[0] = 0.0 * delta;
    }
    for (double& v : y) {
        v += start;
    }
    if (num > 1) {
        y.back() = stop;
    }
    return y;
}

Colormap Colormap::listed(std::string name, std::vector<Rgba> colors) {
    if (colors.empty()) {
        fail("a colormap needs at least one colour");
    }
    Colormap map;
    map.m_name = std::move(name);
    map.m_n = static_cast<int>(colors.size());
    map.m_listed = true;
    map.m_colors = std::move(colors);
    map.build();
    return map;
}

Colormap Colormap::linearSegmented(std::string name, SegmentData data, int n, double gamma) {
    if (n < 1) {
        fail("a colormap needs at least one entry (N >= 1)");
    }
    validate(data.red);
    validate(data.green);
    validate(data.blue);
    if (!data.alpha.empty()) {
        validate(data.alpha);
    }
    Colormap map;
    map.m_name = std::move(name);
    map.m_n = n;
    map.m_listed = false;
    map.m_gamma = gamma;
    map.m_segments = std::move(data);
    map.build();
    return map;
}

Colormap Colormap::fromList(std::string name, const std::vector<Rgba>& colors, int n, double gamma) {
    std::vector<std::pair<double, Rgba>> stops;
    const std::vector<double> values = linspace(0.0, 1.0, static_cast<int>(colors.size()));
    for (std::size_t i = 0; i < colors.size(); ++i) {
        stops.emplace_back(values[i], colors[i]);
    }
    SegmentData data;
    for (const auto& [x, c] : stops) {
        data.red.push_back({x, c.r, c.r});
        data.green.push_back({x, c.g, c.g});
        data.blue.push_back({x, c.b, c.b});
        data.alpha.push_back({x, c.a, c.a});
    }
    return linearSegmented(std::move(name), std::move(data), n, gamma);
}

Colormap Colormap::fromList(std::string name, const std::vector<std::pair<double, Rgba>>& stops, int n,
                            double gamma) {
    bool ok = !stops.empty();
    for (std::size_t i = 0; ok && i < stops.size(); ++i) {
        const double v = stops[i].first;
        ok = !(v < 0.0 || v > 1.0) && (i == 0 || !(v - stops[i - 1].first < 0.0));
    }
    if (!ok) {
        fail("the values passed in the (value, color) pairs must increase monotonically from 0 to 1.");
    }
    SegmentData data;
    for (const auto& [x, c] : stops) {
        data.red.push_back({x, c.r, c.r});
        data.green.push_back({x, c.g, c.g});
        data.blue.push_back({x, c.b, c.b});
        data.alpha.push_back({x, c.a, c.a});
    }
    return linearSegmented(std::move(name), std::move(data), n, gamma);
}

void Colormap::build() {
    const auto n = static_cast<std::size_t>(m_n);
    if (m_listed) {
        m_lut.assign(n + 3, Rgba{0.0, 0.0, 0.0, 0.0});
        std::copy(m_colors.begin(), m_colors.end(), m_lut.begin());
    } else {
        m_lut.assign(n + 3, Rgba{1.0, 1.0, 1.0, 1.0});
        const auto red = lookupTable(m_n, m_segments.red, m_gamma);
        const auto green = lookupTable(m_n, m_segments.green, m_gamma);
        const auto blue = lookupTable(m_n, m_segments.blue, m_gamma);
        for (std::size_t i = 0; i < n; ++i) {
            m_lut[i].r = red[i];
            m_lut[i].g = green[i];
            m_lut[i].b = blue[i];
        }
        if (!m_segments.alpha.empty()) {
            const auto alpha = lookupTable(m_n, m_segments.alpha, 1.0);
            for (std::size_t i = 0; i < n; ++i) {
                m_lut[i].a = alpha[i];
            }
        }
    }
    // _update_lut_extremes
    m_lut[n] = m_underColor ? *m_underColor : m_lut[0];
    m_lut[n + 1] = m_overColor ? *m_overColor : m_lut[n - 1];
    m_lut[n + 2] = m_badColor;
}

namespace {

// the entry of _get_rgba_and_mask for a float (index N + 2 for bad)
std::size_t floatIndex(double x, int n) {
    double xa = x * n;
    if (xa == n) {
        xa = n - 1;  // x == 1 is the last entry, not "over"
    }
    if (std::isnan(xa)) {
        return static_cast<std::size_t>(n) + 2;
    }
    if (xa < 0.0) {
        return static_cast<std::size_t>(n);
    }
    if (xa >= n) {
        return static_cast<std::size_t>(n) + 1;
    }
    return static_cast<std::size_t>(xa);  // truncation, as astype(int)
}

std::size_t intIndex(std::int64_t i, int n) {
    if (i < 0) {
        return static_cast<std::size_t>(n);
    }
    if (i >= n) {
        return static_cast<std::size_t>(n) + 1;
    }
    return static_cast<std::size_t>(i);
}

}  // namespace

Rgba Colormap::operator()(double x, std::optional<double> alpha) const {
    const std::size_t i = floatIndex(x, m_n);
    Rgba out = m_lut[i];
    if (alpha) {
        out.a = clip01(*alpha);
        const Rgba& bad = m_lut.back();
        if (i == static_cast<std::size_t>(m_n) + 2 && bad.r == 0.0 && bad.g == 0.0 && bad.b == 0.0 && bad.a == 0.0) {
            out = {0.0, 0.0, 0.0, 0.0};
        }
    }
    return out;
}

Rgba Colormap::at(std::int64_t index, std::optional<double> alpha) const {
    Rgba out = m_lut[intIndex(index, m_n)];
    if (alpha) {
        out.a = clip01(*alpha);
    }
    return out;
}

namespace {

std::array<std::uint8_t, 4> bytesOf(const std::vector<Rgba>& lut, std::size_t i, bool isBad,
                                    const std::optional<double>& alpha) {
    const Rgba& c = lut[i];
    std::array<std::uint8_t, 4> out{toByte(c.r), toByte(c.g), toByte(c.b), toByte(c.a)};
    if (alpha) {
        const double a = clip01(*alpha);
        out[3] = std::isnan(a) ? std::uint8_t{0} : toByte(a);
        const Rgba& bad = lut.back();
        if (isBad && toByte(bad.r) == 0 && toByte(bad.g) == 0 && toByte(bad.b) == 0 && toByte(bad.a) == 0) {
            out = {0, 0, 0, 0};
        }
    }
    return out;
}

}  // namespace

std::array<std::uint8_t, 4> Colormap::bytes(double x, std::optional<double> alpha) const {
    const std::size_t i = floatIndex(x, m_n);
    return bytesOf(m_lut, i, i == static_cast<std::size_t>(m_n) + 2, alpha);
}

std::array<std::uint8_t, 4> Colormap::bytesAt(std::int64_t index, std::optional<double> alpha) const {
    return bytesOf(m_lut, intIndex(index, m_n), false, alpha);
}

Colormap Colormap::withExtremes(std::optional<Rgba> bad, std::optional<Rgba> under, std::optional<Rgba> over) const {
    Colormap map = *this;
    if (bad) {
        map.m_badColor = *bad;
    }
    if (under) {
        map.m_underColor = under;
    }
    if (over) {
        map.m_overColor = over;
    }
    map.build();
    return map;
}

Colormap Colormap::reversed(std::optional<std::string> name) const {
    Colormap map = *this;
    map.m_name = name ? *name : m_name + "_r";
    if (m_listed) {
        std::reverse(map.m_colors.begin(), map.m_colors.end());
    } else {
        SegmentData& d = map.m_segments;
        for (auto* channel : {&d.red, &d.green, &d.blue, &d.alpha}) {
            std::vector<Segment> rows;
            for (auto it = channel->rbegin(); it != channel->rend(); ++it) {
                rows.push_back({1.0 - it->x, it->y1, it->y0});
            }
            *channel = std::move(rows);
        }
    }
    std::swap(map.m_underColor, map.m_overColor);
    map.build();
    return map;
}

std::span<const std::string_view> bundledColormapNames() {
    static const std::vector<std::string_view> names = [] {
        std::vector<std::string_view> out;
        for (const auto& row : detail::kListedMaps) {
            out.push_back(row.name);
        }
        return out;
    }();
    return names;
}

std::optional<Colormap> findColormap(std::string_view name) {
    const bool reverse = name.size() > 2 && name.substr(name.size() - 2) == "_r";
    const std::string_view base = reverse ? name.substr(0, name.size() - 2) : name;
    for (const auto& row : detail::kListedMaps) {
        if (row.name == base) {
            std::vector<Rgba> colors;
            for (int i = 0; i < row.count; ++i) {
                const auto& c = detail::kListedMapColors[row.first + i];
                colors.push_back({c.r, c.g, c.b, 1.0});
            }
            Colormap map = Colormap::listed(std::string(row.name), std::move(colors));
            return reverse ? map.reversed() : map;
        }
    }
    return std::nullopt;
}

std::optional<std::string> colormapProblem(std::string_view name) {
    if (findColormap(name)) {
        return std::nullopt;
    }
    std::string message = core::pythonStringRepr(std::string(name)) + " is not a bundled colormap; choose from " +
                          joinedNames() + " (or one of them reversed, with \"_r\")";
    const std::string lower = lowerAscii(name);
    for (std::string_view candidate : bundledColormapNames()) {
        for (const std::string& spelled : {std::string(candidate), std::string(candidate) + "_r"}) {
            if (lowerAscii(spelled) == lower) {
                return message + "; did you mean " + core::pythonStringRepr(spelled) + "?";
            }
        }
    }
    return message;
}

Colormap colormap(std::string_view name) {
    if (auto map = findColormap(name)) {
        return *map;
    }
    fail(*colormapProblem(name));
}

std::vector<Rgba> palette(std::string_view name, int n) {
    const Colormap map = colormap(name);
    std::vector<Rgba> out;
    if (map.size() <= 20) {
        for (int k = 0; k < n; ++k) {
            out.push_back(map.at(k % map.size()));
        }
    } else {
        for (double v : linspace(0.05, 0.95, std::max(n, 1))) {
            out.push_back(map(v));
        }
    }
    return out;
}

}  // namespace mdbridge::render
