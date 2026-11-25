#pragma once
// TrueType fonts read and laid out without Qt or FreeType: the metrics of the bundled DejaVu faces
// (text.hpp) as Matplotlib 3.11 measures text with `text.hinting: no_hinting`, i.e. FreeType's
// unhinted outlines and HarfBuzz' shaping (cmap, the GSUB ligatures, the GPOS kerning), in 26.6 fixed
// point. Device independent: the same numbers on every system and at every dpi (Matplotlib's default
// hinting snaps glyphs to the pixel grid; docs/DIFFERENCES.md). The figure layout measures text with
// it; the Qt side draws the glyph runs it produces (text_qt.hpp).
//
// Implemented from the OpenType specification (tables head, hhea, maxp, hmtx, loca, glyf, cmap
// formats 4 / 12, name, GDEF glyph classes, GSUB lookup type 4, GPOS lookup type 2, extension
// lookups) and from the fixed-point rules FreeType and HarfBuzz document for scaling; checked against
// Matplotlib (tests/parity/make_text_cases.py).

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::render {

class FontError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// A point of a glyph outline (TrueType: quadratic curves; consecutive off-curve points have an
// implied on-curve point midway), 26.6 pixels, y up; `contourEnd` marks the last point of a contour.
struct OutlinePoint {
    std::int64_t x = 0;
    std::int64_t y = 0;
    bool onCurve = true;
    bool contourEnd = false;
};

// A glyph's control box (the box of its outline points) in font units (empty for glyphs without an
// outline, e.g. space).
struct GlyphBox {
    int xMin = 0;
    int yMin = 0;
    int xMax = 0;
    int yMax = 0;
    bool empty = true;
};

// The size of a layout as FreeType sets it for FT_Set_Char_Size(points * 64, dpi): the 26.6 pixel
// size (points * 64 truncated, the dpi's integer part, rounded as FreeType rounds) and the 16.16
// scale from font units to 26.6 pixels.
struct FontScale {
    std::int64_t pixels26 = 0;  // the em in 26.6 pixels
    std::int64_t scale = 0;     // FreeType's x_scale (= y_scale): units -> 26.6, 16.16 fixed point
    double pixelSize() const { return static_cast<double>(pixels26) / 64.0; }  // the em in pixels
};
FontScale fontScale(int unitsPerEm, double points, double dpi);

// One glyph of a laid-out line: its pen position in 26.6 pixels from the origin (x offsets included)
// and the index of its first code point in the text (ligatures cover several).
struct PositionedGlyph {
    std::uint16_t glyph = 0;
    std::int64_t x = 0;
    std::int64_t advance = 0;  // 26.6, kerning included
    std::size_t cluster = 0;
};

// A single line of text laid out (newlines are not special: the figure layout splits lines first,
// as Matplotlib's Text does).
struct TextLayout {
    std::vector<PositionedGlyph> glyphs;
    std::int64_t advance = 0;  // pen position after the last glyph, 26.6
    // ink box (glyph control boxes at their pen positions), 26.6; hasInk false for blank text
    std::int64_t inkXMin = 0, inkYMin = 0, inkXMax = 0, inkYMax = 0;
    bool hasInk = false;
    double pixelSize = 0;  // the em in pixels (FontScale::pixelSize)

    // Matplotlib's RendererAgg.get_text_width_height_descent: width = the advance, height = the
    // ink's height, descent = how far the ink goes below the baseline (0 for blank text), in pixels.
    double width() const { return static_cast<double>(advance) / 64.0; }
    double height() const { return hasInk ? static_cast<double>(inkYMax - inkYMin) / 64.0 : 0.0; }
    double descent() const { return hasInk ? static_cast<double>(-inkYMin) / 64.0 : 0.0; }
};

class Font {
public:
    // Parses a TrueType font file (its bytes are kept); throws FontError when it is not one or is
    // malformed.
    static Font fromData(std::string data);

    const std::string& data() const;
    const std::string& family() const;     // name ID 1, e.g. "DejaVu Sans"
    const std::string& subfamily() const;  // name ID 2, e.g. "Book", "Bold", "Oblique"
    const std::string& version() const;    // name ID 5, e.g. "Version 2.35"
    int unitsPerEm() const;
    int ascender() const;   // hhea, font units
    int descender() const;  // hhea, font units (negative)
    int lineGap() const;
    int glyphCount() const;

    std::uint16_t glyphIndex(char32_t codepoint) const;  // 0 (.notdef) when the font lacks it
    int advanceUnits(std::uint16_t glyph) const;
    GlyphBox box(std::uint16_t glyph) const;
    // The glyph's outline at `pixelSize` pixels per em (TextLayout::pixelSize), scaled as FreeType's
    // unhinted loader scales it; drawn by the Qt side as a filled path (the same pixels on every
    // system: no font engine, hinting or system font settings involved).
    std::vector<OutlinePoint> outline(std::uint16_t glyph, double pixelSize) const;

    // Glyphs of `codepoints` after the default ligatures (GSUB liga / clig / rlig / ccmp, type 4),
    // with the index of each glyph's first code point.
    std::vector<std::pair<std::uint16_t, std::size_t>> substitute(const std::u32string& codepoints) const;
    // GPOS kern (type 2) adjustments of a glyph sequence in font units: x advance per glyph and x
    // placement per glyph, applied in lookup order as HarfBuzz does.
    void kern(const std::vector<std::uint16_t>& glyphs, std::vector<int>& xAdvance, std::vector<int>& xPlacement) const;

    // Lays out one line of UTF-8 text (invalid bytes become U+FFFD) at `points` and `dpi`.
    TextLayout layout(std::string_view utf8, double points, double dpi) const;

    struct Impl;

private:
    explicit Font(std::shared_ptr<const Impl> impl);
    std::shared_ptr<const Impl> impl_;
};

// A laid-out line as path commands in pixels, origin on the baseline at the left, y down (paint
// coordinates): each glyph's outline at its pen position, TrueType's implied on-curve points made
// explicit. Fill with the non-zero winding rule.
struct PathCommand {
    enum class Kind : std::uint8_t { Move, Line, Quad, Close };
    Kind kind = Kind::Move;
    double x = 0, y = 0;    // end point
    double cx = 0, cy = 0;  // control point (Quad)
};
std::vector<PathCommand> textPath(const Font& font, const TextLayout& layout);

// UTF-8 -> code points; invalid or truncated sequences become U+FFFD (as Python's "replace").
std::u32string decodeUtf8(std::string_view text);

}  // namespace mdbridge::render
