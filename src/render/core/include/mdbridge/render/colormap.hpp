#pragma once
// Colormaps as Matplotlib 3.11 builds and evaluates them: ListedColormap,
// LinearSegmentedColormap (segment data and from_list), with_extremes (bad / under / over),
// reversed, the lookup by value (float in [0, 1]) or index (integer), the bundled named maps of
// owner decision D6 and MDATK's chain palettes (mdatk/render/coloring.py palette). Qt-free.
//
// Port of matplotlib.colors.Colormap / ListedColormap / LinearSegmentedColormap /
// _create_lookup_table (Matplotlib 3.11.2, Copyright (c) 2012- Matplotlib Development Team; All
// Rights Reserved; Matplotlib licence, packaging/licenses/Matplotlib-LICENSE.txt) and of
// numpy.linspace (NumPy 2.5.5, BSD-3-Clause, packaging/licenses/numpy-LICENSE.txt). Changes from
// Matplotlib: scalar C++ API; invalid segment data is reported when the map is made (Matplotlib:
// at its first use); a map needs at least one colour.

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/render/color.hpp"

namespace mdbridge::render {

class Colormap {
public:
    // One row (x, y0, y1) of a LinearSegmentedColormap channel: between x[i] and x[i+1] the
    // channel goes linearly from y1[i] to y0[i+1].
    struct Segment {
        double x = 0.0;
        double y0 = 0.0;
        double y1 = 0.0;
    };
    // red, green, blue and (optionally: empty = opaque) alpha
    struct SegmentData {
        std::vector<Segment> red, green, blue, alpha;
    };

    // ListedColormap(colors, name): one entry per colour. Throws core::ValueError for no colours.
    static Colormap listed(std::string name, std::vector<Rgba> colors);
    // LinearSegmentedColormap(name, segmentdata, N, gamma). Throws core::ValueError with
    // Matplotlib's "data mapping points must start with x=0 and end with x=1" / "data mapping
    // points must have x in increasing order", and for N < 1.
    static Colormap linearSegmented(std::string name, SegmentData data, int n = 256, double gamma = 1.0);
    // LinearSegmentedColormap.from_list(name, colors, N, gamma): colours equally spaced on [0, 1].
    static Colormap fromList(std::string name, const std::vector<Rgba>& colors, int n = 256, double gamma = 1.0);
    // from_list with (value, colour) pairs; values must increase from 0 to 1 ("the values passed
    // in the (value, color) pairs must increase monotonically from 0 to 1.").
    static Colormap fromList(std::string name, const std::vector<std::pair<double, Rgba>>& stops, int n = 256,
                             double gamma = 1.0);

    const std::string& name() const { return m_name; }
    int size() const { return m_n; }  // Matplotlib's N
    bool isListed() const { return m_listed; }

    // cmap(x) for a float: x * N truncated to an entry; x == 1 is the last entry; x < 0 gives the
    // under colour, x > 1 the over colour, NaN the bad colour. alpha (clipped to [0, 1]) replaces
    // the alpha, except that a bad value stays (0, 0, 0, 0) while the bad colour is all zeros.
    Rgba operator()(double x, std::optional<double> alpha = std::nullopt) const;
    // cmap(i) for an integer: entry i; i < 0 under, i >= N over.
    Rgba at(std::int64_t index, std::optional<double> alpha = std::nullopt) const;
    // cmap(x, bytes=True): the entry's channels as (v * 255) truncated (alpha: clip(alpha) * 255).
    std::array<std::uint8_t, 4> bytes(double x, std::optional<double> alpha = std::nullopt) const;
    std::array<std::uint8_t, 4> bytesAt(std::int64_t index, std::optional<double> alpha = std::nullopt) const;

    Rgba bad() const { return m_lut[static_cast<std::size_t>(m_n) + 2]; }    // get_bad()
    Rgba under() const { return m_lut[static_cast<std::size_t>(m_n)]; }      // get_under()
    Rgba over() const { return m_lut[static_cast<std::size_t>(m_n) + 1]; }   // get_over()
    // The N entries followed by under, over and bad (Matplotlib's _lut).
    const std::vector<Rgba>& lut() const { return m_lut; }

    // with_extremes(bad=, under=, over=): a copy with the given colours set.
    Colormap withExtremes(std::optional<Rgba> bad, std::optional<Rgba> under = std::nullopt,
                          std::optional<Rgba> over = std::nullopt) const;
    // reversed(name): name + "_r" by default; under and over swap, bad stays.
    Colormap reversed(std::optional<std::string> name = std::nullopt) const;

private:
    Colormap() = default;
    void build();

    std::string m_name;
    int m_n = 0;
    bool m_listed = true;
    double m_gamma = 1.0;
    std::vector<Rgba> m_colors;  // listed
    SegmentData m_segments;      // linear segmented
    Rgba m_badColor{0.0, 0.0, 0.0, 0.0};
    std::optional<Rgba> m_underColor;
    std::optional<Rgba> m_overColor;
    std::vector<Rgba> m_lut;
};

// The colormaps MDBridge bundles (owner decision D6), in this order: tab10, tab20, tab20b, tab20c,
// Set1, Set2, Set3, Pastel1, Pastel2, Dark2, Accent, Paired, viridis, plasma, inferno, magma,
// cividis. Each is also known with "_r" (reversed), as in Matplotlib's registry.
std::span<const std::string_view> bundledColormapNames();
// The bundled map of that (case-sensitive) name, or nothing.
std::optional<Colormap> findColormap(std::string_view name);
// The same; throws core::ValueError naming the bundled maps for any other name (RN-13).
Colormap colormap(std::string_view name);
// The message colormap(name) would throw, or nothing (check() of visual.chain_palette).
std::optional<std::string> colormapProblem(std::string_view name);

// MDATK's chain palette (mdatk/render/coloring.py palette(name, n)): a map with at most 20 entries
// gives entries k % N for k < n; a larger one samples numpy.linspace(0.05, 0.95, max(n, 1)).
std::vector<Rgba> palette(std::string_view name, int n);

// numpy.linspace(start, stop, num) (endpoint included), bit for bit.
std::vector<double> linspace(double start, double stop, int num);

}  // namespace mdbridge::render
