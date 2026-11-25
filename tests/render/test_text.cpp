// The bundled fonts, Matplotlib's font choice and the text metrics (src/render/core/font.cpp,
// text.cpp) against Matplotlib 3.11.2 (data/render/text/text_cases.json,
// tests/parity/make_text_cases.py): the committed data/fonts files are Matplotlib's; the face its
// font manager picks; width / height / descent and glyph positions exactly as Matplotlib measures
// with text.hinting "no_hinting" (26.6 fixed point), and within a bound of its default hinting.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/render/text.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using core::Array;
using core::Value;
using render::FontFace;

namespace {

const std::filesystem::path kFonts = MDBRIDGE_FONT_DIR;

const Value& data() {
    static const Value value = test::readJson(test::kData / "render" / "text" / "text_cases.json");
    return value;
}
const Value& at(const Value& v, std::string_view key) { return *v.find(key); }
const Array& list(const Value& v) { return *v.get<Array>(); }
std::string str(const Value& v) { return *v.get<std::string>(); }
std::int64_t integer(const Value& v) { return *v.get<std::int64_t>(); }
double num(const Value& v) {
    if (const auto* d = v.get<double>()) {
        return *d;
    }
    return static_cast<double>(*v.get<std::int64_t>());
}

const render::FontSet& fonts() {
    static const render::FontSet set =
        render::FontSet::load([](std::string_view file) { return test::readFile(kFonts / std::string(file)); });
    return set;
}

FontFace faceOfFile(const std::string& file) {
    for (const render::BundledFont& f : render::bundledFonts()) {
        if (f.file == file) {
            return f.face;
        }
    }
    throw std::runtime_error("not a bundled face: " + file);
}

}  // namespace

TEST(BundledFonts, AreMatplotlibsFiles) {
    EXPECT_EQ(str(at(data(), "matplotlib")), "3.11.2");
    const Array& expected = list(at(data(), "fonts"));
    ASSERT_EQ(expected.size(), render::kFontFaceCount);
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const Value& e = expected[i];
        const render::BundledFont& f = render::bundledFonts()[i];
        SCOPED_TRACE(f.file);
        EXPECT_EQ(str(at(e, "file")), f.file);
        EXPECT_EQ(str(at(e, "sha256")), f.sha256);
        EXPECT_EQ(core::sha256Hex(test::readFile(kFonts / std::string(f.file))), f.sha256);  // data/fonts
        const render::Font& font = fonts().font(f.face);
        EXPECT_EQ(font.family(), str(at(e, "family")));
        EXPECT_EQ(font.family(), f.family);
        EXPECT_EQ(font.subfamily(), f.style);
        EXPECT_EQ(font.version(), str(at(e, "version")));
        EXPECT_EQ(font.version(), "Version " + std::string(render::kDejaVuVersion));
        EXPECT_EQ(font.unitsPerEm(), integer(at(e, "units_per_em")));
        EXPECT_EQ(font.ascender(), integer(at(e, "ascender")));
        EXPECT_EQ(font.descender(), integer(at(e, "descender")));
        EXPECT_EQ(font.glyphCount(), integer(at(e, "glyphs")));
        EXPECT_EQ(render::bundledFont(f.face).file, f.file);
    }
}

TEST(BundledFonts, RejectsOtherData) {
    EXPECT_THROW(render::Font::fromData("not a font"), render::FontError);
    EXPECT_THROW(render::Font::fromData(std::string(64, '\0')), render::FontError);
    std::string truncated = test::readFile(kFonts / "DejaVuSans.ttf").substr(0, 4000);
    EXPECT_THROW(render::Font::fromData(truncated), render::FontError);
    EXPECT_THROW(render::FontSet::load([](std::string_view file) {
                     // the Bold face where Book belongs: rejected by its checksum
                     return test::readFile(kFonts / (file == "DejaVuSans.ttf" ? std::string("DejaVuSans-Bold.ttf")
                                                                               : std::string(file)));
                 }),
                 render::FontError);
}

TEST(FontChoice, AsMatplotlibsFontManager) {
    int n = 0;
    for (const Value& c : list(at(data(), "choice"))) {
        std::vector<std::string> families;
        for (const Value& f : list(at(c, "family"))) {
            families.push_back(str(f));
        }
        const Value& w = at(c, "weight");
        const int weight = w.get<std::string>() ? *render::fontWeight(str(w)) : static_cast<int>(integer(w));
        const std::string style = str(at(c, "style"));
        SCOPED_TRACE(core::toJson(c, -1));
        EXPECT_EQ(render::bundledFont(render::chooseFace(families, style, weight)).file, str(at(c, "file")));
        ++n;
    }
    EXPECT_GT(n, 700);
    EXPECT_EQ(render::chooseFace("monospace", "normal", 700), FontFace::MonoBold);  // MDATK's sequence letters
    EXPECT_EQ(render::chooseFace("sans-serif", "italic", 400), FontFace::SansOblique);
    EXPECT_EQ(render::fontWeight("bold"), 700);
    EXPECT_EQ(render::fontWeight("650"), 650);
    EXPECT_EQ(render::fontWeight("bolder"), std::nullopt);
    EXPECT_EQ(render::fontWeight("0"), std::nullopt);
}

TEST(TextUnits, AsMatplotlib) {
    EXPECT_DOUBLE_EQ(render::pointsToPixels(10.0, 100.0), 10.0 * 100.0 / 72.0);
    EXPECT_DOUBLE_EQ(render::pointsToPixels(72.0, 150.0), 150.0);
    EXPECT_DOUBLE_EQ(render::pixelsToPoints(150.0, 150.0), 72.0);
    EXPECT_DOUBLE_EQ(render::qtPointSize(16.0, 96.0), 12.0);
    // FreeType's size for FT_Set_Char_Size(10 pt * 64, 100 dpi): (640 * 100 + 36) / 72 in 26.6
    const render::FontScale s = render::fontScale(2048, 10.0, 100.0);
    EXPECT_EQ(s.pixels26, (640 * 100 + 36) / 72);
    EXPECT_EQ(render::fontScale(2048, 10.0, 72.9).pixels26, 640);  // the dpi's integer part
    EXPECT_EQ(render::fontScale(2048, 10.0, 0.0).pixels26, 640);   // FreeType's default 72 dpi
}

TEST(TextMetrics, AsMatplotlibUnhinted) {
    int n = 0;
    for (const Value& c : list(at(data(), "metrics"))) {
        const std::string text = str(at(c, "text"));
        const double points = num(at(c, "points")), dpi = num(at(c, "dpi"));
        SCOPED_TRACE(str(at(c, "file")) + " '" + text + "' " + std::to_string(points) + " pt " + std::to_string(dpi) + " dpi");
        const render::TextLayout t = fonts().layout(faceOfFile(str(at(c, "file"))), text, points, dpi);
        EXPECT_EQ(t.advance, integer(at(c, "width")));
        EXPECT_EQ(t.hasInk ? t.inkYMax - t.inkYMin : 0, integer(at(c, "height")));
        EXPECT_EQ(t.hasInk ? -t.inkYMin : 0, integer(at(c, "descent")));
        const Array& xs = list(at(c, "x"));
        ASSERT_EQ(t.glyphs.size(), xs.size());
        for (std::size_t i = 0; i < xs.size(); ++i) {
            EXPECT_EQ(t.glyphs[i].x, integer(xs[i])) << "glyph " << i;
        }
        // Matplotlib's default hinting (MDATK's figures) rounds each glyph to the pixel grid
        const Array& hinted = list(at(c, "hinted"));
        const double glyphs = static_cast<double>(std::max<std::size_t>(xs.size(), 1));
        EXPECT_LE(std::abs(t.width() - num(hinted[0])), 0.7 * glyphs);
        EXPECT_LE(std::abs(t.height() - num(hinted[1])), 2.0);
        EXPECT_LE(std::abs(t.descent() - num(hinted[2])), 1.5);
        ++n;
    }
    EXPECT_GT(n, 1000);
}

TEST(TextMetrics, LayoutDetails) {
    const render::TextLayout blank = fonts().layout(FontFace::Sans, "  ", 10, 100);
    EXPECT_FALSE(blank.hasInk);
    EXPECT_GT(blank.width(), 0);
    EXPECT_EQ(blank.height(), 0);
    EXPECT_EQ(blank.descent(), 0);
    const render::TextLayout empty = fonts().layout(FontFace::Sans, "", 10, 100);
    EXPECT_TRUE(empty.glyphs.empty());
    EXPECT_EQ(empty.width(), 0);
    // ligature: "fi" is one glyph whose cluster is the 'f'
    const render::TextLayout fi = fonts().layout(FontFace::Sans, "afi", 10, 100);
    ASSERT_EQ(fi.glyphs.size(), 2u);
    EXPECT_EQ(fi.glyphs[1].cluster, 1u);
    // kerning: "AV" is narrower than A and V apart
    const render::TextLayout av = fonts().layout(FontFace::Sans, "AV", 10, 100);
    EXPECT_LT(av.advance, fonts().layout(FontFace::Sans, "A", 10, 100).advance +
                              fonts().layout(FontFace::Sans, "V", 10, 100).advance);
    // monospace: every glyph the same advance, no kerning
    const render::TextLayout mono = fonts().layout(FontFace::MonoBold, "AVWi", 10, 100);
    EXPECT_EQ(mono.advance, 4 * fonts().layout(FontFace::MonoBold, "A", 10, 100).advance);
    // a missing glyph is .notdef; invalid UTF-8 becomes U+FFFD
    EXPECT_EQ(fonts().font(FontFace::Sans).glyphIndex(0x10FFFD), 0);
    EXPECT_EQ(render::decodeUtf8("a\xff" "b\xe2\x82"), std::u32string({U'a', 0xFFFD, U'b', 0xFFFD}));
    EXPECT_EQ(render::decodeUtf8("\xc3\xa5\xe2\x88\x92"), std::u32string({0xE5, 0x2212}));
    EXPECT_EQ(render::decodeUtf8("\xed\xa0\x80"), std::u32string({0xFFFD, 0xFFFD, 0xFFFD}));
}

TEST(BatchPlatform, OffscreenWithoutDisplayServer) {
#if defined(_WIN32) || defined(__APPLE__)
    EXPECT_EQ(render::batchQtPlatform(), std::nullopt);
#else
    EXPECT_EQ(render::batchQtPlatform(), "offscreen");
#endif
    for (const char* c : {"render", "charts", "preview", "fonts"}) {
        EXPECT_TRUE(render::drawsWithoutWindow(c)) << c;
    }
    for (const char* c : {"gui", "analyze", "info", ""}) {
        EXPECT_FALSE(render::drawsWithoutWindow(c)) << c;
    }
}

TEST(TextPath, OutlinesAtThePenPositions) {
    using K = render::PathCommand::Kind;
    const render::Font& sans = fonts().font(FontFace::Sans);
    // 'l' is a rectangle: four lines
    const render::TextLayout l = sans.layout("l", 10, 100);
    const std::vector<render::PathCommand> lp = render::textPath(sans, l);
    ASSERT_FALSE(lp.empty());
    EXPECT_EQ(lp.front().kind, K::Move);
    EXPECT_EQ(lp.back().kind, K::Close);
    for (const auto& c : lp) {
        EXPECT_TRUE(c.kind != K::Quad);
    }
    // every point of a text's path lies in its ink box (y down), and curves are closed per contour
    const render::TextLayout t = sans.layout("Gyroscope 0.5 Å", 12, 150);
    const std::vector<render::PathCommand> path = render::textPath(sans, t);
    int moves = 0, closes = 0, quads = 0;
    for (const auto& c : path) {
        moves += c.kind == K::Move;
        closes += c.kind == K::Close;
        quads += c.kind == K::Quad;
        if (c.kind != K::Close) {
            EXPECT_GE(c.x, t.inkXMin / 64.0 - 1e-9);
            EXPECT_LE(c.x, t.inkXMax / 64.0 + 1e-9);
            EXPECT_GE(-c.y, t.inkYMin / 64.0 - 1e-9);
            EXPECT_LE(-c.y, t.inkYMax / 64.0 + 1e-9);
        }
    }
    EXPECT_EQ(moves, closes);
    EXPECT_GT(moves, 14);
    EXPECT_GT(quads, 40);
    // the outline's box is the glyph box the metrics use
    const std::vector<render::OutlinePoint> o = sans.outline(sans.glyphIndex(U'g'), l.pixelSize);
    std::int64_t ymin = 0;
    for (const auto& q : o) {
        ymin = std::min(ymin, q.y);
    }
    EXPECT_EQ(-ymin, sans.layout("g", 10, 100).inkYMin * -1);
    EXPECT_TRUE(o.back().contourEnd);
}
