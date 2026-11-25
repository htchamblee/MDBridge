// Colours, colormaps, chain palettes and MDATK's themes (src/render/core) against what Matplotlib
// 3.11.2 and MDATK v0.0.4.1 give (tests/data/render/colors/color_cases.json, written by
// tests/parity/make_color_cases.py; MDBridge's improvements from tests/parity/improvements_colors.py).
// Every number must agree bit for bit.

#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/color.hpp"
#include "mdbridge/render/colormap.hpp"
#include "mdbridge/render/theme.hpp"

namespace {

using mdbridge::core::Array;
using mdbridge::core::Value;
using namespace mdbridge::render;

const Value& cases() {
    static const Value doc = [] {
        const std::filesystem::path path = std::filesystem::path(MDBRIDGE_TEST_DATA) / "render/colors/color_cases.json";
        std::ifstream in(path, std::ios::binary);
        std::ostringstream text;
        text << in.rdbuf();
        return mdbridge::core::parseJson(text.str(), path.string());
    }();
    return doc;
}

const Value& at(const Value& v, std::string_view key) {
    const Value* found = v.find(key);
    if (found == nullptr) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *found;
}
const Array& arr(const Value& v) { return *v.get<Array>(); }
const std::string& str(const Value& v) { return *v.get<std::string>(); }
double num(const Value& v) {
    if (const auto* i = v.get<std::int64_t>()) {
        return static_cast<double>(*i);
    }
    return *v.get<double>();
}
std::optional<double> optNum(const Value& v) { return v.isNull() ? std::nullopt : std::optional<double>(num(v)); }
bool boolean(const Value& v) { return *v.get<bool>(); }

Rgba rgba(const Value& v) {
    const Array& a = arr(v);
    return {num(a[0]), num(a[1]), num(a[2]), num(a[3])};
}

bool same(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}
bool same(const Rgba& a, const Rgba& b) { return same(a.r, b.r) && same(a.g, b.g) && same(a.b, b.b) && same(a.a, b.a); }

std::string show(const Rgba& c) {
    std::ostringstream out;
    out.precision(17);
    out << "(" << c.r << ", " << c.g << ", " << c.b << ", " << c.a << ")";
    return out.str();
}

std::string showValue(const Value& v) { return mdbridge::core::toJson(v, -1); }

// Expect `expected` (an RGBA list or {"error": message}) from `compute`.
template <class F>
void expectRgba(const Value& expected, F compute, const std::string& what) {
    if (const Value* error = expected.find("error")) {
        try {
            const Rgba got = compute();
            ADD_FAILURE() << what << ": expected error " << str(*error) << ", got " << show(got);
        } catch (const mdbridge::core::ValueError& e) {
            EXPECT_EQ(std::string(e.what()), str(*error)) << what;
        }
        return;
    }
    try {
        const Rgba got = compute();
        EXPECT_TRUE(same(got, rgba(expected)))
            << what << ": got " << show(got) << ", expected " << show(rgba(expected));
    } catch (const mdbridge::core::ValueError& e) {
        ADD_FAILURE() << what << ": unexpected error " << e.what();
    }
}

std::optional<Rgba> parseExtreme(const Value& spec, std::string_view key) {
    const Value* v = spec.find(key);
    return v == nullptr ? std::nullopt : std::optional<Rgba>(toRgba(str(*v)));
}

Colormap build(const Value& spec) {
    const std::string kind = str(at(spec, "kind"));
    const std::string name = str(at(spec, "name"));
    std::optional<Colormap> map;
    if (kind == "bundled") {
        map = colormap(name);
    } else if (kind == "from_list") {
        std::vector<Rgba> colors;
        for (const Value& c : arr(at(spec, "colors"))) {
            colors.push_back(toRgba(str(c)));
        }
        const int n = static_cast<int>(num(at(spec, "N")));
        const double gamma = num(at(spec, "gamma"));
        const Value* values = spec.find("values");
        if (values != nullptr && !values->isNull()) {
            std::vector<std::pair<double, Rgba>> stops;
            for (std::size_t i = 0; i < colors.size(); ++i) {
                stops.emplace_back(num(arr(*values)[i]), colors[i]);
            }
            map = Colormap::fromList(name, stops, n, gamma);
        } else {
            map = Colormap::fromList(name, colors, n, gamma);
        }
    } else if (kind == "segments") {
        Colormap::SegmentData data;
        const Value& d = at(spec, "data");
        auto rows = [&](std::string_view key) {
            std::vector<Colormap::Segment> out;
            if (const Value* channel = d.find(key)) {
                for (const Value& r : arr(*channel)) {
                    out.push_back({num(arr(r)[0]), num(arr(r)[1]), num(arr(r)[2])});
                }
            }
            return out;
        };
        data.red = rows("red");
        data.green = rows("green");
        data.blue = rows("blue");
        data.alpha = rows("alpha");
        map = Colormap::linearSegmented(name, data, static_cast<int>(num(at(spec, "N"))), num(at(spec, "gamma")));
    } else if (kind == "listed") {
        std::vector<Rgba> colors;
        for (const Value& c : arr(at(spec, "colors"))) {
            colors.push_back(toRgba(str(c)));
        }
        map = Colormap::listed(name, colors);
    } else {
        throw std::runtime_error("unknown spec kind " + kind);
    }
    if (const Value* r = spec.find("reversed"); r != nullptr && boolean(*r)) {
        map = map->reversed();
    }
    if (const Value* ext = spec.find("extremes")) {
        map = map->withExtremes(parseExtreme(*ext, "bad"), parseExtreme(*ext, "under"), parseExtreme(*ext, "over"));
    }
    return *map;
}

TEST(Color, ToRgbaStrings) {
    const Array& list = arr(at(cases(), "to_rgba"));
    ASSERT_GT(list.size(), 2000u);
    for (const Value& c : list) {
        const std::string color = str(arr(c)[0]);
        const auto alpha = optNum(arr(c)[1]);
        expectRgba(arr(c)[2], [&] { return toRgba(color, alpha); }, "to_rgba" + showValue(c));
        EXPECT_EQ(colorProblem(color).has_value(), isColorLike(color) == false);
    }
}

TEST(Color, ToRgbaTuples) {
    for (const Value& c : arr(at(cases(), "to_rgba_tuple"))) {
        std::vector<double> comps;
        for (const Value& v : arr(arr(c)[0])) {
            comps.push_back(num(v));
        }
        const auto alpha = optNum(arr(c)[1]);
        expectRgba(arr(c)[2], [&] { return toRgba(std::span<const double>(comps), alpha); }, showValue(c));
    }
}

TEST(Color, ToHex) {
    for (const Value& c : arr(at(cases(), "to_hex"))) {
        EXPECT_EQ(toHex(rgba(arr(c)[0]), boolean(arr(c)[1])), str(arr(c)[2])) << showValue(c);
    }
}

TEST(Color, NamesAndCycle) {
    EXPECT_EQ(namedColorNames().size(), static_cast<std::size_t>(num(at(cases(), "names"))));
    for (std::string_view name : namedColorNames()) {
        EXPECT_TRUE(isColorLike(name)) << name;
    }
    ASSERT_EQ(defaultColorCycle().size(), 10u);
    EXPECT_EQ(toHex(defaultColorCycle()[0]), "#1f77b4");
    EXPECT_TRUE(same(defaultColorCycle()[3], toRgba("tab:red")));
}

TEST(Colormap, Samples) {
    const Array& maps = arr(at(cases(), "colormaps"));
    std::size_t count = 0;
    for (const Value& m : maps) {
        const Value& spec = at(m, "spec");
        const std::string what = showValue(spec);
        const Colormap map = build(spec);
        EXPECT_EQ(map.name(), str(at(m, "name"))) << what;
        ASSERT_EQ(map.size(), static_cast<int>(num(at(m, "N")))) << what;
        EXPECT_TRUE(same(map.under(), rgba(at(m, "under")))) << what;
        EXPECT_TRUE(same(map.over(), rgba(at(m, "over")))) << what;
        EXPECT_TRUE(same(map.bad(), rgba(at(m, "bad")))) << what;
        if (const Value* lut = m.find("lut")) {
            ASSERT_EQ(map.lut().size(), arr(*lut).size()) << what;
            for (std::size_t i = 0; i < map.lut().size(); ++i) {
                EXPECT_TRUE(same(map.lut()[i], rgba(arr(*lut)[i])))
                    << what << " row " << i << ": " << show(map.lut()[i]);
            }
        }
        for (const Value& s : arr(at(m, "samples"))) {
            const double x = num(arr(s)[0]);
            const Rgba got = map(x, optNum(arr(s)[1]));
            EXPECT_TRUE(same(got, rgba(arr(s)[2]))) << what << " x " << showValue(s) << ": " << show(got);
            ++count;
        }
        for (const Value& s : arr(at(m, "int_samples"))) {
            const auto i = *arr(s)[0].get<std::int64_t>();
            const Rgba got = map.at(i, optNum(arr(s)[1]));
            EXPECT_TRUE(same(got, rgba(arr(s)[2]))) << what << " i " << showValue(s) << ": " << show(got);
            ++count;
        }
        for (const Value& s : arr(at(m, "bytes"))) {
            const auto alpha = optNum(arr(s)[1]);
            const auto got = arr(s).size() == 4 ? map.bytesAt(*arr(s)[0].get<std::int64_t>(), alpha)
                                                : map.bytes(num(arr(s)[0]), alpha);
            const Array& want = arr(arr(s)[2]);
            for (std::size_t k = 0; k < 4; ++k) {
                EXPECT_EQ(static_cast<int>(got[k]), static_cast<int>(num(want[k])))
                    << what << " bytes " << showValue(s);
            }
            ++count;
        }
    }
    EXPECT_GT(count, 2000u);
}

TEST(Colormap, Errors) {
    for (const Value& e : arr(at(cases(), "colormap_errors"))) {
        try {
            build(arr(e)[0]);
            ADD_FAILURE() << "no error for " << showValue(arr(e)[0]);
        } catch (const mdbridge::core::ValueError& error) {
            EXPECT_EQ(std::string(error.what()), str(arr(e)[1])) << showValue(arr(e)[0]);
        }
    }
    EXPECT_THROW(Colormap::listed("empty", {}), mdbridge::core::ValueError);
    EXPECT_THROW(Colormap::fromList("none", std::vector<Rgba>{}), mdbridge::core::ValueError);
    EXPECT_THROW(Colormap::fromList("n0", {Rgba{}, Rgba{}}, 0), mdbridge::core::ValueError);
}

TEST(Colormap, ChainPalettes) {
    for (const Value& p : arr(at(cases(), "palettes"))) {
        const auto got = palette(str(arr(p)[0]), static_cast<int>(num(arr(p)[1])));
        const Array& want = arr(arr(p)[2]);
        ASSERT_EQ(got.size(), want.size()) << showValue(arr(p)[0]) << " " << num(arr(p)[1]);
        for (std::size_t i = 0; i < got.size(); ++i) {
            EXPECT_TRUE(same(got[i], rgba(want[i]))) << showValue(arr(p)[0]) << " n " << num(arr(p)[1]) << " #" << i;
        }
    }
}

TEST(Colormap, Names) {
    for (const Value& c : arr(at(cases(), "colormap_names"))) {
        const std::string name = str(arr(c)[0]);
        const Value& problem = arr(c)[2];
        if (problem.isNull()) {
            EXPECT_FALSE(colormapProblem(name).has_value()) << name;
            EXPECT_EQ(colormap(name).name(), name);
        } else {
            EXPECT_EQ(colormapProblem(name).value_or("<none>"), str(problem));
            EXPECT_THROW(colormap(name), mdbridge::core::ValueError);
        }
    }
    EXPECT_EQ(bundledColormapNames().size(), 17u);
}

TEST(Colormap, Linspace) {
    // numpy.linspace's rounding: i * step + start, the end point exact
    const auto y = linspace(0.05, 0.95, 7);
    ASSERT_EQ(y.size(), 7u);
    const double step = (0.95 - 0.05) / 6;
    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(same(y[static_cast<std::size_t>(i)], i * step + 0.05));
    }
    EXPECT_EQ(y.back(), 0.95);
    EXPECT_EQ(linspace(0.05, 0.95, 1), std::vector<double>{0.05});
    EXPECT_TRUE(linspace(0.0, 1.0, 0).empty());
}

TEST(Theme, Tokens) {
    for (const Theme* t : {&lightTheme(), &darkTheme()}) {
        const Value& rec = at(at(cases(), "themes"), t->name);
        EXPECT_EQ(t->name, str(at(rec, "name")));
        EXPECT_EQ(t->surface, str(at(rec, "surface")));
        EXPECT_EQ(t->page, str(at(rec, "page")));
        EXPECT_EQ(t->ink, str(at(rec, "ink")));
        EXPECT_EQ(t->ink2, str(at(rec, "ink2")));
        EXPECT_EQ(t->muted, str(at(rec, "muted")));
        EXPECT_EQ(t->grid, str(at(rec, "grid")));
        EXPECT_EQ(t->axis, str(at(rec, "axis")));
        const Array& series = arr(at(rec, "series"));
        ASSERT_EQ(t->series().size(), series.size());
        for (std::size_t i = 0; i < series.size(); ++i) {
            EXPECT_EQ(t->series()[i], str(series[i]));
        }
        for (const Value& s : arr(at(rec, "series_color"))) {
            EXPECT_EQ(t->seriesColor(*arr(s)[0].get<std::int64_t>()), str(arr(s)[1])) << showValue(s);
        }
    }
    for (const Value& g : arr(at(cases(), "get_theme"))) {
        EXPECT_EQ(theme(str(arr(g)[0])).name, str(arr(g)[1])) << showValue(g);
    }
    EXPECT_TRUE(isThemeName("light"));
    EXPECT_TRUE(isThemeName("dark"));
    EXPECT_FALSE(isThemeName("Dark"));
    EXPECT_EQ(sequentialBlue().size(), 13u);
}

TEST(Theme, Ramps) {
    for (const Theme* t : {&lightTheme(), &darkTheme()}) {
        const Value& rec = at(at(cases(), "themes"), t->name);
        const Colormap seq = t->sequential();
        const Value& want = at(rec, "sequential");
        EXPECT_EQ(seq.name(), str(at(want, "name")));
        ASSERT_EQ(seq.lut().size(), arr(at(want, "lut")).size());
        for (std::size_t i = 0; i < seq.lut().size(); ++i) {
            EXPECT_TRUE(same(seq.lut()[i], rgba(arr(at(want, "lut"))[i]))) << t->name << " row " << i;
        }
        const Colormap bad = seq.withExtremes(toRgba(t->grid));
        const Value& wantBad = at(rec, "sequential_bad");
        EXPECT_TRUE(same(bad.bad(), rgba(at(wantBad, "bad"))));
        for (const Value& s : arr(at(wantBad, "samples"))) {
            EXPECT_TRUE(same(bad(num(arr(s)[0])), rgba(arr(s)[1]))) << showValue(s);
        }
        for (const Value& o : arr(at(rec, "ordinal"))) {
            const int n = static_cast<int>(num(at(o, "n")));
            const Colormap map = t->ordinal(n);
            EXPECT_EQ(map.name(), str(at(o, "name")));
            const Array& colors = arr(at(o, "colors"));
            ASSERT_EQ(static_cast<std::size_t>(map.size()), colors.size()) << t->name << " n " << n;
            for (std::size_t i = 0; i < colors.size(); ++i) {
                EXPECT_TRUE(same(map.lut()[i], rgba(colors[i]))) << t->name << " n " << n << " #" << i;
            }
            if (n > 5) {  // RN-18: every class its own colour
                for (std::size_t i = 1; i < colors.size(); ++i) {
                    EXPECT_FALSE(same(map.lut()[i], map.lut()[i - 1])) << t->name << " n " << n;
                }
            }
        }
    }
}

void expectTick(const Value& want, const AxesStyle& s, const std::string& what) {
    EXPECT_TRUE(same(s.tickColor, rgba(at(want, "tick_color")))) << what;
    EXPECT_EQ(s.tickWidth, num(at(want, "tick_width"))) << what;
    EXPECT_EQ(s.tickLength, num(at(want, "tick_length"))) << what;
    EXPECT_TRUE(same(s.tickLabelColor, rgba(at(want, "label_color")))) << what;
    EXPECT_EQ(s.tickLabelSize, num(at(want, "label_size"))) << what;
    EXPECT_EQ(s.grid, boolean(at(want, "grid_visible"))) << what;
    EXPECT_TRUE(same(s.gridColor, rgba(at(want, "grid_color")))) << what;
    EXPECT_EQ(s.gridWidth, num(at(want, "grid_width"))) << what;
    EXPECT_EQ(s.gridLineStyle, str(at(want, "grid_linestyle"))) << what;
    EXPECT_TRUE(same(s.axisLabelColor, rgba(at(want, "axis_label_color")))) << what;
    EXPECT_EQ(s.axisLabelSize, num(at(want, "axis_label_size"))) << what;
}

TEST(Theme, Chrome) {
    for (const Value& c : arr(at(cases(), "style_axes"))) {
        const std::string what = showValue(Value(Array{arr(c)[0], arr(c)[1], arr(c)[2]}));
        const AxesStyle s = axesStyle(theme(str(arr(c)[0])), boolean(arr(c)[1]), num(arr(c)[2]));
        const Value& want = arr(c)[3];
        EXPECT_TRUE(same(s.face, rgba(at(want, "face")))) << what;
        const Value& spines = at(want, "spines");
        using Side = std::pair<const char*, const SpineStyle*>;
        for (const auto& [side, spine] : {Side{"left", &s.left}, Side{"right", &s.right}, Side{"bottom", &s.bottom},
                                          Side{"top", &s.top}}) {
            const Array& w = arr(at(spines, side));
            EXPECT_EQ(spine->visible, boolean(w[0])) << what << side;
            EXPECT_TRUE(same(spine->color, rgba(w[1]))) << what << side;
            EXPECT_EQ(spine->width, num(w[2])) << what << side;
        }
        expectTick(at(want, "x"), s, what + " x");
        expectTick(at(want, "y"), s, what + " y");
        EXPECT_TRUE(same(s.titleColor, rgba(at(want, "title_color")))) << what;
        EXPECT_EQ(s.titleSize, num(at(want, "title_size"))) << what;
        EXPECT_EQ(s.axisBelow == AxisBelow::True ? "True" : "line", str(at(want, "axisbelow"))) << what;
    }
    for (const Value& c : arr(at(cases(), "style_legend"))) {
        const LegendStyle s = legendStyle(theme(str(arr(c)[0])), num(arr(c)[1]));
        const Value& want = arr(c)[2];
        EXPECT_TRUE(same(s.face, rgba(at(want, "face"))));
        EXPECT_TRUE(same(s.edge, rgba(at(want, "edge"))));
        EXPECT_EQ(s.edgeWidth, num(at(want, "edge_width")));
        EXPECT_TRUE(same(s.textColor, rgba(at(want, "text_color"))));
        EXPECT_EQ(s.textSize, num(at(want, "text_size")));
    }
}

}  // namespace
