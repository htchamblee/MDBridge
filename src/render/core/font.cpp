// TrueType reading and single-line layout (font.hpp). The table layouts follow the OpenType
// specification; the fixed-point scaling is FreeType's (FT_Set_Char_Size, FT_MulFix / FT_MulDiv /
// FT_DivFix rounding, unhinted advances) and HarfBuzz' (font scale and em_mult for GPOS values), as
// both libraries document them. Matplotlib 3.11 lays text out with them (FT2Font + libraqm).
#include "mdbridge/render/font.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace mdbridge::render {

namespace {

constexpr char32_t kReplacement = 0xFFFD;

// --- fixed point, as FreeType rounds (sign and magnitude) ---------------------------------------

std::int64_t mulDiv(std::int64_t a, std::int64_t b, std::int64_t c) {
    const bool negative = (a < 0) != (b < 0);
    const std::int64_t ua = a < 0 ? -a : a;
    const std::int64_t ub = b < 0 ? -b : b;
    const std::int64_t d = c > 0 ? (ua * ub + c / 2) / c : std::numeric_limits<std::int32_t>::max();
    return negative ? -d : d;
}
std::int64_t mulFix(std::int64_t a, std::int64_t b) {
    const bool negative = (a < 0) != (b < 0);
    const std::int64_t ua = a < 0 ? -a : a;
    const std::int64_t ub = b < 0 ? -b : b;
    const std::int64_t c = (ua * ub + 0x8000) >> 16;
    return negative ? -c : c;
}

// --- big-endian reading with bounds checks ---------------------------------------------------

class Reader {
public:
    explicit Reader(const std::string& data) : data_(data) {}
    std::size_t size() const { return data_.size(); }
    void need(std::size_t offset, std::size_t n) const {
        if (offset > data_.size() || n > data_.size() - offset) {
            throw FontError("malformed font: a table reaches past the end of the file");
        }
    }
    std::uint8_t u8(std::size_t o) const {
        need(o, 1);
        return static_cast<std::uint8_t>(data_[o]);
    }
    std::uint16_t u16(std::size_t o) const {
        need(o, 2);
        return static_cast<std::uint16_t>((static_cast<unsigned>(static_cast<std::uint8_t>(data_[o])) << 8) |
                                          static_cast<std::uint8_t>(data_[o + 1]));
    }
    std::int16_t s16(std::size_t o) const { return static_cast<std::int16_t>(u16(o)); }
    std::uint32_t u32(std::size_t o) const {
        return (static_cast<std::uint32_t>(u16(o)) << 16) | u16(o + 2);
    }

private:
    const std::string& data_;
};

std::uint32_t tag(const char (&t)[5]) {
    return (static_cast<std::uint32_t>(static_cast<unsigned char>(t[0])) << 24) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(t[1])) << 16) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(t[2])) << 8) |
           static_cast<std::uint32_t>(static_cast<unsigned char>(t[3]));
}

void appendUtf8(std::string& out, char32_t c) {
    if (c < 0x80) {
        out += static_cast<char>(c);
    } else if (c < 0x800) {
        out += static_cast<char>(0xC0 | (c >> 6));
        out += static_cast<char>(0x80 | (c & 0x3F));
    } else if (c < 0x10000) {
        out += static_cast<char>(0xE0 | (c >> 12));
        out += static_cast<char>(0x80 | ((c >> 6) & 0x3F));
        out += static_cast<char>(0x80 | (c & 0x3F));
    } else {
        out += static_cast<char>(0xF0 | (c >> 18));
        out += static_cast<char>(0x80 | ((c >> 12) & 0x3F));
        out += static_cast<char>(0x80 | ((c >> 6) & 0x3F));
        out += static_cast<char>(0x80 | (c & 0x3F));
    }
}

// The scripts that select OpenType features (the others count as common: they take the script of
// the text around them, as libraqm itemises runs).
enum class Script { Common, Latin, Greek, Cyrillic };

Script scriptOf(char32_t c) {
    if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == 0xAA || c == 0xBA ||
        (c >= 0xC0 && c <= 0x24F && c != 0xD7 && c != 0xF7) || (c >= 0x1E00 && c <= 0x1EFF)) {
        return Script::Latin;
    }
    if ((c >= 0x370 && c <= 0x3FF && c != 0x37E && c != 0x387) || (c >= 0x1F00 && c <= 0x1FFF)) {
        return Script::Greek;
    }
    if (c >= 0x400 && c <= 0x52F) {
        return Script::Cyrillic;
    }
    return Script::Common;
}

}  // namespace

std::u32string decodeUtf8(std::string_view text) {
    std::u32string out;
    out.reserve(text.size());
    std::size_t i = 0;
    while (i < text.size()) {
        const auto b0 = static_cast<unsigned char>(text[i]);
        if (b0 < 0x80) {
            out += static_cast<char32_t>(b0);
            ++i;
            continue;
        }
        int n = 0;
        char32_t c = 0;
        char32_t minimum = 0;
        if (b0 >= 0xC2 && b0 <= 0xDF) {
            n = 1, c = b0 & 0x1F, minimum = 0x80;
        } else if (b0 >= 0xE0 && b0 <= 0xEF) {
            n = 2, c = b0 & 0x0F, minimum = 0x800;
        } else if (b0 >= 0xF0 && b0 <= 0xF4) {
            n = 3, c = b0 & 0x07, minimum = 0x10000;
        } else {
            out += kReplacement;
            ++i;
            continue;
        }
        int k = 1;
        for (; k <= n && i + static_cast<std::size_t>(k) < text.size(); ++k) {
            const auto b = static_cast<unsigned char>(text[i + static_cast<std::size_t>(k)]);
            if ((b & 0xC0) != 0x80) {
                break;
            }
            c = (c << 6) | (b & 0x3F);
        }
        if (k <= n) {  // truncated: the lead byte and the continuation bytes seen are one error
            out += kReplacement;
            i += static_cast<std::size_t>(k);
            continue;
        }
        if (c < minimum || c > 0x10FFFF || (c >= 0xD800 && c <= 0xDFFF)) {
            out += kReplacement;
            ++i;
            continue;
        }
        out += c;
        i += static_cast<std::size_t>(n) + 1;
    }
    return out;
}

FontScale fontScale(int unitsPerEm, double points, double dpi) {
    // FT_Set_Char_Size(face, (FT_F26Dot6)(points * 64), 0, (FT_UInt)dpi, (FT_UInt)dpi), as
    // Matplotlib's FT2Font::set_size calls it: both values truncated
    const double p = points * 64.0;
    const auto charWidth = static_cast<std::int64_t>(std::isfinite(p) ? std::clamp(p, 0.0, 1e12) : 0.0);
    std::int64_t resolution = static_cast<std::int64_t>(std::isfinite(dpi) ? std::clamp(dpi, 0.0, 1e9) : 0.0);
    if (resolution == 0) {
        resolution = 72;  // FreeType's default
    }
    FontScale s;
    s.pixels26 = (charWidth * resolution + 36) / 72;  // FT_REQUEST_WIDTH
    s.scale = unitsPerEm > 0 ? mulDiv(s.pixels26, 0x10000, unitsPerEm) : 0;  // FT_DivFix
    return s;
}

// --- the parsed font ---------------------------------------------------------------------------

struct Font::Impl {
    std::string data;
    std::string family, subfamily, version;
    int unitsPerEm = 0, ascender = 0, descender = 0, lineGap = 0, numGlyphs = 0;
    int numHMetrics = 0;
    bool longLoca = false;
    std::size_t hmtx = 0, loca = 0, glyf = 0, glyfSize = 0;
    std::size_t cmap = 0;  // the chosen subtable
    int cmapFormat = 0;
    std::size_t gdefClasses = 0, gdefMarkAttach = 0;  // 0: none
    std::size_t gsub = 0, gpos = 0;

    Reader r() const { return Reader(data); }

    std::uint16_t glyphIndex(char32_t c) const;
    int advance(std::uint16_t g) const;
    using Point = OutlinePoint;
    void points(std::uint16_t g, std::int64_t scale, int depth, std::vector<Point>& out) const;
    GlyphBox box(std::uint16_t g, std::int64_t scale) const;  // the points' box, 26.6 (font units at 0x10000)

    // OpenType layout
    int coverage(std::size_t table, std::uint16_t g) const;  // -1 when not covered
    int classOf(std::size_t classDef, std::uint16_t g) const;
    int glyphClass(std::uint16_t g) const { return gdefClasses ? classOf(gdefClasses, g) : 0; }
    bool ignored(std::uint16_t g, std::uint16_t flags) const;
    std::vector<std::uint16_t> lookups(std::size_t table, Script script, std::initializer_list<const char*> features) const;
    // the subtable offsets of a lookup with its type, extension lookups resolved
    std::vector<std::size_t> subtables(std::size_t table, std::uint16_t lookup, int wantedType, int extensionType,
                                       std::uint16_t& flags) const;
};

namespace {

std::string nameString(const Reader& r, std::size_t name, std::uint16_t id) {
    const std::uint16_t count = r.u16(name + 2);
    const std::size_t strings = name + r.u16(name + 4);
    std::optional<std::string> mac;
    for (std::uint16_t i = 0; i < count; ++i) {
        const std::size_t rec = name + 6 + 12u * i;
        const std::uint16_t platform = r.u16(rec), encoding = r.u16(rec + 2), language = r.u16(rec + 4);
        if (r.u16(rec + 6) != id) {
            continue;
        }
        const std::size_t length = r.u16(rec + 8), offset = strings + r.u16(rec + 10);
        r.need(offset, length);
        if ((platform == 3 && (encoding == 1 || encoding == 10) && language == 0x409) || platform == 0) {
            std::string out;
            for (std::size_t k = 0; k + 1 < length; k += 2) {
                char32_t c = r.u16(offset + k);
                if (c >= 0xD800 && c <= 0xDBFF && k + 3 < length) {
                    const char32_t low = r.u16(offset + k + 2);
                    c = 0x10000 + ((c - 0xD800) << 10) + (low - 0xDC00);
                    k += 2;
                }
                appendUtf8(out, c);
            }
            return out;
        }
        if (platform == 1 && encoding == 0 && !mac) {
            std::string out;
            for (std::size_t k = 0; k < length; ++k) {
                const auto c = static_cast<unsigned char>(r.u8(offset + k));
                appendUtf8(out, c < 0x80 ? c : kReplacement);
            }
            mac = out;
        }
    }
    return mac.value_or("");
}

}  // namespace

std::uint16_t Font::Impl::glyphIndex(char32_t c) const {
    const Reader rd = r();
    if (cmapFormat == 12) {
        const std::uint32_t groups = rd.u32(cmap + 12);
        std::uint32_t lo = 0, hi = groups;
        while (lo < hi) {
            const std::uint32_t mid = lo + (hi - lo) / 2;
            const std::size_t g = cmap + 16 + 12u * mid;
            if (c < rd.u32(g)) {
                hi = mid;
            } else if (c > rd.u32(g + 4)) {
                lo = mid + 1;
            } else {
                const std::uint32_t glyph = rd.u32(g + 8) + (static_cast<std::uint32_t>(c) - rd.u32(g));
                return glyph < static_cast<std::uint32_t>(numGlyphs) ? static_cast<std::uint16_t>(glyph) : 0;
            }
        }
        return 0;
    }
    if (cmapFormat == 4 && c <= 0xFFFF) {
        const std::size_t segments = rd.u16(cmap + 6) / 2u;
        const std::size_t ends = cmap + 14, starts = ends + 2 * segments + 2, deltas = starts + 2 * segments,
                          ranges = deltas + 2 * segments;
        std::size_t lo = 0, hi = segments;
        while (lo < hi) {  // the first segment whose end >= c
            const std::size_t mid = lo + (hi - lo) / 2;
            if (rd.u16(ends + 2 * mid) < c) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if (lo == segments || rd.u16(starts + 2 * lo) > c) {
            return 0;
        }
        const std::uint16_t delta = rd.u16(deltas + 2 * lo), range = rd.u16(ranges + 2 * lo);
        std::uint32_t glyph = 0;
        if (range == 0) {
            glyph = (static_cast<std::uint32_t>(c) + delta) & 0xFFFF;
        } else {
            glyph = rd.u16(ranges + 2 * lo + range + 2 * (static_cast<std::uint32_t>(c) - rd.u16(starts + 2 * lo)));
            if (glyph != 0) {
                glyph = (glyph + delta) & 0xFFFF;
            }
        }
        return glyph < static_cast<std::uint32_t>(numGlyphs) ? static_cast<std::uint16_t>(glyph) : 0;
    }
    return 0;
}

int Font::Impl::advance(std::uint16_t g) const {
    const int k = std::min<int>(g, numHMetrics - 1);
    return r().u16(hmtx + 4u * static_cast<unsigned>(k));
}

// The outline points of a glyph (on- and off-curve; FreeType's control box is their box), scaled
// to 26.6 with `scale` as FreeType's unhinted loader scales them: each simple glyph's points with
// FT_MulFix, composite components moved by their scaled offsets (or matched points) after their
// 2 x 2 transform. scale 0x10000 gives font units.
void Font::Impl::points(std::uint16_t g, std::int64_t scale, int depth, std::vector<Point>& out) const {
    if (g >= numGlyphs || depth > 8) {
        return;
    }
    const Reader rd = r();
    const std::size_t start = longLoca ? rd.u32(loca + 4u * g) : 2u * rd.u16(loca + 2u * g);
    const std::size_t end = longLoca ? rd.u32(loca + 4u * (g + 1u)) : 2u * rd.u16(loca + 2u * (g + 1u));
    if (end <= start || end > glyfSize) {
        return;
    }
    const std::size_t o = glyf + start;
    const std::int16_t contours = rd.s16(o);
    if (contours > 0) {
        const std::size_t n = static_cast<std::size_t>(rd.u16(o + 10 + 2u * static_cast<unsigned>(contours - 1))) + 1;
        std::size_t p = o + 10 + 2u * static_cast<unsigned>(contours);
        p += 2 + rd.u16(p);  // instructions
        std::vector<std::uint8_t> flags;
        flags.reserve(n);
        while (flags.size() < n) {
            const std::uint8_t f = rd.u8(p++);
            flags.push_back(f);
            if (f & 0x08) {
                for (int k = rd.u8(p++); k > 0 && flags.size() < n; --k) {
                    flags.push_back(f);
                }
            }
        }
        std::vector<std::int64_t> xs(n), ys(n);
        std::int64_t v = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (flags[i] & 0x02) {
                const int d = rd.u8(p++);
                v += (flags[i] & 0x10) ? d : -d;
            } else if (!(flags[i] & 0x10)) {
                v += rd.s16(p);
                p += 2;
            }
            xs[i] = v;
        }
        v = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (flags[i] & 0x04) {
                const int d = rd.u8(p++);
                v += (flags[i] & 0x20) ? d : -d;
            } else if (!(flags[i] & 0x20)) {
                v += rd.s16(p);
                p += 2;
            }
            ys[i] = v;
        }
        std::vector<bool> ends(n, false);
        for (int c = 0; c < contours; ++c) {
            const std::size_t last = rd.u16(o + 10 + 2u * static_cast<unsigned>(c));
            if (last < n) {
                ends[last] = true;
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            out.push_back({mulFix(xs[i], scale), mulFix(ys[i], scale), (flags[i] & 0x01) != 0, ends[i]});
        }
        return;
    }
    if (contours == 0) {
        return;
    }
    std::size_t p = o + 10;  // composite glyph
    for (;;) {
        const std::uint16_t flags = rd.u16(p), component = rd.u16(p + 2);
        p += 4;
        std::int64_t a1 = 0, a2 = 0;
        if (flags & 0x0001) {  // ARG_1_AND_2_ARE_WORDS
            a1 = (flags & 0x0002) ? rd.s16(p) : rd.u16(p);
            a2 = (flags & 0x0002) ? rd.s16(p + 2) : rd.u16(p + 2);
            p += 4;
        } else {
            a1 = (flags & 0x0002) ? static_cast<std::int8_t>(rd.u8(p)) : rd.u8(p);
            a2 = (flags & 0x0002) ? static_cast<std::int8_t>(rd.u8(p + 1)) : rd.u8(p + 1);
            p += 2;
        }
        std::int64_t xx = 0x10000, xy = 0, yx = 0, yy = 0x10000;  // 16.16 (F2Dot14 << 2)
        bool transformed = false;
        if (flags & 0x0008) {  // WE_HAVE_A_SCALE
            xx = yy = static_cast<std::int64_t>(rd.s16(p)) * 4;
            p += 2;
            transformed = true;
        } else if (flags & 0x0040) {  // WE_HAVE_AN_X_AND_Y_SCALE
            xx = static_cast<std::int64_t>(rd.s16(p)) * 4;
            yy = static_cast<std::int64_t>(rd.s16(p + 2)) * 4;
            p += 4;
            transformed = true;
        } else if (flags & 0x0080) {  // WE_HAVE_A_TWO_BY_TWO
            xx = static_cast<std::int64_t>(rd.s16(p)) * 4;
            yx = static_cast<std::int64_t>(rd.s16(p + 2)) * 4;
            xy = static_cast<std::int64_t>(rd.s16(p + 4)) * 4;
            yy = static_cast<std::int64_t>(rd.s16(p + 6)) * 4;
            p += 8;
            transformed = true;
        }
        const std::size_t first = out.size();
        points(component, scale, depth + 1, out);
        if (transformed) {
            for (std::size_t i = first; i < out.size(); ++i) {
                const Point q = out[i];
                out[i].x = mulFix(q.x, xx) + mulFix(q.y, xy);
                out[i].y = mulFix(q.x, yx) + mulFix(q.y, yy);
            }
        }
        std::int64_t dx = 0, dy = 0;
        if (flags & 0x0002) {  // ARGS_ARE_XY_VALUES
            dx = mulFix(a1, scale);
            dy = mulFix(a2, scale);
        } else {  // matched points: parent point a1 meets component point a2
            const std::size_t k1 = static_cast<std::size_t>(a1), k2 = first + static_cast<std::size_t>(a2);
            if (k1 < first && k2 < out.size()) {
                dx = out[k1].x - out[k2].x;
                dy = out[k1].y - out[k2].y;
            }
        }
        for (std::size_t i = first; i < out.size(); ++i) {
            out[i].x += dx;
            out[i].y += dy;
        }
        if (!(flags & 0x0020)) {  // MORE_COMPONENTS
            break;
        }
    }
}

GlyphBox Font::Impl::box(std::uint16_t g, std::int64_t scale) const {
    std::vector<Point> pts;
    points(g, scale, 0, pts);
    if (pts.empty()) {
        return {};
    }
    GlyphBox b{std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::min(),
               std::numeric_limits<int>::min(), false};
    for (const Point& q : pts) {
        b.xMin = std::min(b.xMin, static_cast<int>(q.x));
        b.yMin = std::min(b.yMin, static_cast<int>(q.y));
        b.xMax = std::max(b.xMax, static_cast<int>(q.x));
        b.yMax = std::max(b.yMax, static_cast<int>(q.y));
    }
    return b;
}

int Font::Impl::coverage(std::size_t table, std::uint16_t g) const {
    const Reader rd = r();
    const std::uint16_t format = rd.u16(table);
    if (format == 1) {
        int lo = 0, hi = rd.u16(table + 2);
        while (lo < hi) {
            const int mid = lo + (hi - lo) / 2;
            const std::uint16_t v = rd.u16(table + 4 + 2u * static_cast<unsigned>(mid));
            if (v == g) {
                return mid;
            }
            (v < g ? lo = mid + 1 : hi = mid);
        }
        return -1;
    }
    if (format == 2) {
        int lo = 0, hi = rd.u16(table + 2);
        while (lo < hi) {
            const int mid = lo + (hi - lo) / 2;
            const std::size_t rec = table + 4 + 6u * static_cast<unsigned>(mid);
            if (g < rd.u16(rec)) {
                hi = mid;
            } else if (g > rd.u16(rec + 2)) {
                lo = mid + 1;
            } else {
                return rd.u16(rec + 4) + (g - rd.u16(rec));
            }
        }
    }
    return -1;
}

int Font::Impl::classOf(std::size_t table, std::uint16_t g) const {
    const Reader rd = r();
    const std::uint16_t format = rd.u16(table);
    if (format == 1) {
        const std::uint16_t first = rd.u16(table + 2), count = rd.u16(table + 4);
        return g >= first && g - first < count ? rd.u16(table + 6 + 2u * (g - first)) : 0;
    }
    if (format == 2) {
        int lo = 0, hi = rd.u16(table + 2);
        while (lo < hi) {
            const int mid = lo + (hi - lo) / 2;
            const std::size_t rec = table + 4 + 6u * static_cast<unsigned>(mid);
            if (g < rd.u16(rec)) {
                hi = mid;
            } else if (g > rd.u16(rec + 2)) {
                lo = mid + 1;
            } else {
                return rd.u16(rec + 4);
            }
        }
    }
    return 0;
}

bool Font::Impl::ignored(std::uint16_t g, std::uint16_t flags) const {
    if ((flags & 0xFF1E) == 0) {
        return false;
    }
    const int cls = glyphClass(g);
    if ((cls == 1 && (flags & 0x2)) || (cls == 2 && (flags & 0x4))) {
        return true;
    }
    if (cls == 3) {
        if (flags & 0x8) {
            return true;
        }
        if (flags & 0x10) {  // mark filtering set: not supported, marks ignored
            return true;
        }
        const int attach = flags >> 8;
        if (attach != 0 && gdefMarkAttach != 0 && classOf(gdefMarkAttach, g) != attach) {
            return true;
        }
    }
    return false;
}

std::vector<std::uint16_t> Font::Impl::lookups(std::size_t table, Script script,
                                               std::initializer_list<const char*> features) const {
    std::vector<std::uint16_t> out;
    if (table == 0) {
        return out;
    }
    const Reader rd = r();
    const std::size_t scripts = table + rd.u16(table + 4), featureList = table + rd.u16(table + 6);
    const char* wanted = script == Script::Latin ? "latn" : script == Script::Greek ? "grek"
                         : script == Script::Cyrillic ? "cyrl" : "DFLT";
    std::size_t scriptTable = 0;
    for (const char* t : {wanted, "DFLT", "dflt", "latn"}) {
        const std::uint32_t want = (static_cast<std::uint32_t>(static_cast<unsigned char>(t[0])) << 24) |
                                   (static_cast<std::uint32_t>(static_cast<unsigned char>(t[1])) << 16) |
                                   (static_cast<std::uint32_t>(static_cast<unsigned char>(t[2])) << 8) |
                                   static_cast<unsigned char>(t[3]);
        const std::uint16_t count = rd.u16(scripts);
        for (std::uint16_t i = 0; i < count && !scriptTable; ++i) {
            if (rd.u32(scripts + 2 + 6u * i) == want) {
                scriptTable = scripts + rd.u16(scripts + 2 + 6u * i + 4);
            }
        }
        if (scriptTable) {
            break;
        }
    }
    if (!scriptTable || rd.u16(scriptTable) == 0) {
        return out;
    }
    const std::size_t langSys = scriptTable + rd.u16(scriptTable);
    std::vector<std::uint16_t> featureIndices;
    if (rd.u16(langSys + 2) != 0xFFFF) {
        featureIndices.push_back(rd.u16(langSys + 2));  // required feature
    }
    for (std::uint16_t i = 0, n = rd.u16(langSys + 4); i < n; ++i) {
        featureIndices.push_back(rd.u16(langSys + 6 + 2u * i));
    }
    const std::uint16_t featureCount = rd.u16(featureList);
    for (const std::uint16_t f : featureIndices) {
        if (f >= featureCount) {
            continue;
        }
        const std::size_t rec = featureList + 2 + 6u * f;
        const std::uint32_t t = rd.u32(rec);
        bool use = false;
        for (const char* name : features) {
            use = use || t == ((static_cast<std::uint32_t>(static_cast<unsigned char>(name[0])) << 24) |
                               (static_cast<std::uint32_t>(static_cast<unsigned char>(name[1])) << 16) |
                               (static_cast<std::uint32_t>(static_cast<unsigned char>(name[2])) << 8) |
                               static_cast<unsigned char>(name[3]));
        }
        if (!use) {
            continue;
        }
        const std::size_t feature = featureList + rd.u16(rec + 4);
        for (std::uint16_t i = 0, n = rd.u16(feature + 2); i < n; ++i) {
            out.push_back(rd.u16(feature + 4 + 2u * i));
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<std::size_t> Font::Impl::subtables(std::size_t table, std::uint16_t lookup, int wantedType,
                                               int extensionType, std::uint16_t& flags) const {
    std::vector<std::size_t> out;
    const Reader rd = r();
    const std::size_t list = table + rd.u16(table + 8);
    if (lookup >= rd.u16(list)) {
        return out;
    }
    const std::size_t l = list + rd.u16(list + 2 + 2u * lookup);
    const std::uint16_t type = rd.u16(l);
    flags = rd.u16(l + 2);
    for (std::uint16_t i = 0, n = rd.u16(l + 4); i < n; ++i) {
        std::size_t sub = l + rd.u16(l + 6 + 2u * i);
        int subType = type;
        if (type == extensionType) {
            subType = rd.u16(sub + 2);
            sub += rd.u32(sub + 4);
        }
        if (subType == wantedType) {
            out.push_back(sub);
        }
    }
    return out;
}

// --- Font ---------------------------------------------------------------------------------------

Font::Font(std::shared_ptr<const Impl> impl) : impl_(std::move(impl)) {}

Font Font::fromData(std::string data) {
    auto impl = std::make_shared<Impl>();
    impl->data = std::move(data);
    const Reader rd = impl->r();
    if (rd.size() < 12) {
        throw FontError("not a TrueType font: the file is too short");
    }
    const std::uint32_t version = rd.u32(0);
    if (version != 0x00010000 && version != tag("true")) {
        throw FontError("not a TrueType font (sfnt version " + std::to_string(version) + ")");
    }
    std::size_t head = 0, hhea = 0, maxp = 0, name = 0, cmapTable = 0, gdef = 0;
    for (std::uint16_t i = 0, n = rd.u16(4); i < n; ++i) {
        const std::size_t rec = 12 + 16u * i;
        const std::uint32_t t = rd.u32(rec);
        const std::size_t offset = rd.u32(rec + 8), length = rd.u32(rec + 12);
        rd.need(offset, length);
        if (t == tag("head")) head = offset;
        else if (t == tag("hhea")) hhea = offset;
        else if (t == tag("maxp")) maxp = offset;
        else if (t == tag("hmtx")) impl->hmtx = offset;
        else if (t == tag("loca")) impl->loca = offset;
        else if (t == tag("glyf")) impl->glyf = offset, impl->glyfSize = length;
        else if (t == tag("cmap")) cmapTable = offset;
        else if (t == tag("name")) name = offset;
        else if (t == tag("GDEF")) gdef = offset;
        else if (t == tag("GSUB")) impl->gsub = offset;
        else if (t == tag("GPOS")) impl->gpos = offset;
    }
    if (!head || !hhea || !maxp || !impl->hmtx || !impl->loca || !impl->glyf || !cmapTable) {
        throw FontError("not a TrueType font with outlines: a required table is missing");
    }
    impl->unitsPerEm = rd.u16(head + 18);
    impl->longLoca = rd.s16(head + 50) == 1;
    impl->ascender = rd.s16(hhea + 4);
    impl->descender = rd.s16(hhea + 6);
    impl->lineGap = rd.s16(hhea + 8);
    impl->numHMetrics = rd.u16(hhea + 34);
    impl->numGlyphs = rd.u16(maxp + 4);
    if (impl->unitsPerEm < 16 || impl->numHMetrics < 1 || impl->numGlyphs < 1) {
        throw FontError("malformed font: bad head / hhea / maxp values");
    }
    rd.need(impl->hmtx, 4u * static_cast<unsigned>(impl->numHMetrics));
    rd.need(impl->loca, (impl->longLoca ? 4u : 2u) * (static_cast<unsigned>(impl->numGlyphs) + 1));
    // the Unicode cmap FreeType selects: UCS-4 (format 12) before UCS-2 (format 4)
    int best = 0;
    for (std::uint16_t i = 0, n = rd.u16(cmapTable + 2); i < n; ++i) {
        const std::size_t rec = cmapTable + 4 + 8u * i;
        const std::uint16_t platform = rd.u16(rec), encoding = rd.u16(rec + 2);
        const std::size_t sub = cmapTable + rd.u32(rec + 4);
        const std::uint16_t format = rd.u16(sub);
        int rank = 0;
        if (format == 12 && ((platform == 3 && encoding == 10) || (platform == 0 && (encoding == 4 || encoding == 6)))) {
            rank = 2;
        } else if (format == 4 && ((platform == 3 && encoding == 1) || platform == 0)) {
            rank = 1;
        }
        if (rank > best) {
            best = rank;
            impl->cmap = sub;
            impl->cmapFormat = format;
        }
    }
    if (best == 0) {
        throw FontError("font has no Unicode character map (cmap format 4 or 12)");
    }
    if (impl->cmapFormat == 12) {
        rd.need(impl->cmap, 16 + 12u * rd.u32(impl->cmap + 12));
    } else {
        rd.need(impl->cmap, rd.u16(impl->cmap + 2));
    }
    if (gdef) {
        if (const std::uint16_t o = rd.u16(gdef + 4)) impl->gdefClasses = gdef + o;
        if (const std::uint16_t o = rd.u16(gdef + 10)) impl->gdefMarkAttach = gdef + o;
    }
    if (name) {
        impl->family = nameString(rd, name, 1);
        impl->subfamily = nameString(rd, name, 2);
        impl->version = nameString(rd, name, 5);
    }
    return Font(std::move(impl));
}

const std::string& Font::data() const { return impl_->data; }
const std::string& Font::family() const { return impl_->family; }
const std::string& Font::subfamily() const { return impl_->subfamily; }
const std::string& Font::version() const { return impl_->version; }
int Font::unitsPerEm() const { return impl_->unitsPerEm; }
int Font::ascender() const { return impl_->ascender; }
int Font::descender() const { return impl_->descender; }
int Font::lineGap() const { return impl_->lineGap; }
int Font::glyphCount() const { return impl_->numGlyphs; }
std::uint16_t Font::glyphIndex(char32_t c) const { return impl_->glyphIndex(c); }
int Font::advanceUnits(std::uint16_t g) const { return impl_->advance(g); }
GlyphBox Font::box(std::uint16_t g) const { return impl_->box(g, 0x10000); }

std::vector<OutlinePoint> Font::outline(std::uint16_t g, double pixelSize) const {
    std::vector<OutlinePoint> out;
    const std::int64_t scale = mulDiv(static_cast<std::int64_t>(std::llround(pixelSize * 64.0)), 0x10000, impl_->unitsPerEm);
    impl_->points(g, scale, 0, out);
    return out;
}

namespace {

// substitute / kern for one script run
std::vector<std::pair<std::uint16_t, std::size_t>> substituteRun(const Font::Impl& f, const std::u32string& text,
                                                                 std::size_t begin, std::size_t end, Script script) {
    std::vector<std::pair<std::uint16_t, std::size_t>> glyphs;
    for (std::size_t i = begin; i < end; ++i) {
        glyphs.emplace_back(f.glyphIndex(text[i]), i);
    }
    if (!f.gsub) {
        return glyphs;
    }
    const Reader rd = f.r();
    for (const std::uint16_t lookup : f.lookups(f.gsub, script, {"ccmp", "locl", "rlig", "liga", "clig", "calt"})) {
        std::uint16_t flags = 0;
        const std::vector<std::size_t> subs = f.subtables(f.gsub, lookup, 4, 7, flags);
        if (subs.empty()) {
            continue;
        }
        for (std::size_t i = 0; i < glyphs.size(); ++i) {
            if (f.ignored(glyphs[i].first, flags)) {
                continue;
            }
            bool done = false;
            for (const std::size_t sub : subs) {
                if (done || rd.u16(sub) != 1) {
                    continue;
                }
                const int cov = f.coverage(sub + rd.u16(sub + 2), glyphs[i].first);
                if (cov < 0 || cov >= rd.u16(sub + 4)) {
                    continue;
                }
                const std::size_t set = sub + rd.u16(sub + 6 + 2u * static_cast<unsigned>(cov));
                for (std::uint16_t k = 0, n = rd.u16(set); k < n && !done; ++k) {
                    const std::size_t lig = set + rd.u16(set + 2 + 2u * k);
                    const std::uint16_t components = rd.u16(lig + 2);
                    std::vector<std::size_t> matched{i};
                    std::size_t j = i + 1;
                    for (std::uint16_t c = 1; c < components; ++c) {
                        while (j < glyphs.size() && f.ignored(glyphs[j].first, flags)) {
                            ++j;
                        }
                        if (j >= glyphs.size() || glyphs[j].first != rd.u16(lig + 4 + 2u * (c - 1u))) {
                            break;
                        }
                        matched.push_back(j++);
                    }
                    if (matched.size() != components) {
                        continue;
                    }
                    glyphs[i].first = rd.u16(lig);
                    for (std::size_t m = matched.size(); m-- > 1;) {
                        glyphs.erase(glyphs.begin() + static_cast<std::ptrdiff_t>(matched[m]));
                    }
                    done = true;
                }
            }
        }
    }
    return glyphs;
}

struct ValueRecord {
    int xPlacement = 0;
    int xAdvance = 0;
};
std::size_t valueSize(std::uint16_t format) {
    std::size_t n = 0;
    for (int b = 0; b < 8; ++b) {
        n += (format >> b) & 1u;
    }
    return 2 * n;
}
ValueRecord readValue(const Reader& rd, std::size_t o, std::uint16_t format) {
    ValueRecord v;
    if (format & 0x1) {
        v.xPlacement = rd.s16(o);
        o += 2;
    }
    if (format & 0x2) {
        o += 2;
    }
    if (format & 0x4) {
        v.xAdvance = rd.s16(o);
    }
    return v;
}

void kernRun(const Font::Impl& f, const std::vector<std::uint16_t>& glyphs, std::size_t begin, std::size_t end,
             Script script, std::vector<int>& xAdvance, std::vector<int>& xPlacement) {
    if (!f.gpos) {
        return;
    }
    const Reader rd = f.r();
    for (const std::uint16_t lookup : f.lookups(f.gpos, script, {"kern"})) {
        std::uint16_t flags = 0;
        const std::vector<std::size_t> subs = f.subtables(f.gpos, lookup, 2, 9, flags);
        if (subs.empty()) {
            continue;
        }
        std::size_t i = begin;
        while (i < end) {
            if (f.ignored(glyphs[i], flags)) {
                ++i;
                continue;
            }
            std::size_t j = i + 1;
            while (j < end && f.ignored(glyphs[j], flags)) {
                ++j;
            }
            if (j >= end) {
                break;
            }
            bool applied = false;
            std::uint16_t format2 = 0;
            for (const std::size_t sub : subs) {
                const std::uint16_t format = rd.u16(sub);
                const int cov = f.coverage(sub + rd.u16(sub + 2), glyphs[i]);
                if (cov < 0) {
                    continue;
                }
                const std::uint16_t vf1 = rd.u16(sub + 4), vf2 = rd.u16(sub + 6);
                const std::size_t s1 = valueSize(vf1), s2 = valueSize(vf2);
                std::size_t record = 0;
                if (format == 1) {
                    if (cov >= rd.u16(sub + 8)) {
                        continue;
                    }
                    const std::size_t set = sub + rd.u16(sub + 10 + 2u * static_cast<unsigned>(cov));
                    const std::size_t stride = 2 + s1 + s2;
                    int lo = 0, hi = rd.u16(set);
                    while (lo < hi) {
                        const int mid = lo + (hi - lo) / 2;
                        const std::size_t rec = set + 2 + stride * static_cast<unsigned>(mid);
                        const std::uint16_t second = rd.u16(rec);
                        if (second == glyphs[j]) {
                            record = rec + 2;
                            break;
                        }
                        (second < glyphs[j] ? lo = mid + 1 : hi = mid);
                    }
                } else if (format == 2) {
                    const int c1 = f.classOf(sub + rd.u16(sub + 8), glyphs[i]);
                    const int c2 = f.classOf(sub + rd.u16(sub + 10), glyphs[j]);
                    const std::uint16_t n1 = rd.u16(sub + 12), n2 = rd.u16(sub + 14);
                    if (c1 < n1 && c2 < n2) {
                        record = sub + 16 + (static_cast<std::size_t>(c1) * n2 + static_cast<std::size_t>(c2)) * (s1 + s2);
                    }
                }
                if (record == 0) {
                    continue;
                }
                const ValueRecord v1 = readValue(rd, record, vf1), v2 = readValue(rd, record + s1, vf2);
                xAdvance[i] += v1.xAdvance;
                xPlacement[i] += v1.xPlacement;
                xAdvance[j] += v2.xAdvance;
                xPlacement[j] += v2.xPlacement;
                applied = true;
                format2 = vf2;
                break;
            }
            i = applied ? (format2 != 0 ? j + 1 : j) : i + 1;
        }
    }
}

// script runs: common characters take the script of the run before them (leading ones that of the
// first strong character after them)
std::vector<std::pair<std::size_t, Script>> scriptRuns(const std::u32string& text) {
    std::vector<std::pair<std::size_t, Script>> runs;  // (start, script)
    Script current = Script::Common;
    for (const char32_t c : text) {
        if (const Script s = scriptOf(c); s != Script::Common) {
            current = s;
            break;
        }
    }
    runs.emplace_back(0, current);
    for (std::size_t i = 0; i < text.size(); ++i) {
        const Script s = scriptOf(text[i]);
        if (s != Script::Common && s != current) {
            current = s;
            runs.emplace_back(i, s);
        }
    }
    return runs;
}

}  // namespace

std::vector<std::pair<std::uint16_t, std::size_t>> Font::substitute(const std::u32string& text) const {
    std::vector<std::pair<std::uint16_t, std::size_t>> out;
    const auto runs = scriptRuns(text);
    for (std::size_t r = 0; r < runs.size(); ++r) {
        const std::size_t end = r + 1 < runs.size() ? runs[r + 1].first : text.size();
        auto part = substituteRun(*impl_, text, runs[r].first, end, runs[r].second);
        out.insert(out.end(), part.begin(), part.end());
    }
    return out;
}

void Font::kern(const std::vector<std::uint16_t>& glyphs, std::vector<int>& xAdvance, std::vector<int>& xPlacement) const {
    xAdvance.assign(glyphs.size(), 0);
    xPlacement.assign(glyphs.size(), 0);
    kernRun(*impl_, glyphs, 0, glyphs.size(), Script::Latin, xAdvance, xPlacement);
}

TextLayout Font::layout(std::string_view utf8, double points, double dpi) const {
    const Impl& f = *impl_;
    const std::u32string text = decodeUtf8(utf8);
    const FontScale size = fontScale(f.unitsPerEm, points, dpi);
    // HarfBuzz' font scale from FreeType's size (hb-ft) and its em_mult for GPOS values
    const std::int64_t hbScale = (size.scale * f.unitsPerEm + (1 << 15)) >> 16;
    const std::int64_t emMult = (hbScale << 16) / f.unitsPerEm;
    const auto em = [emMult](int v) { return (static_cast<std::int64_t>(v) * emMult + 32768) >> 16; };

    std::vector<std::uint16_t> glyphs;
    std::vector<std::size_t> clusters;
    std::vector<int> kernAdvance, kernPlacement;
    const auto runs = scriptRuns(text);
    for (std::size_t r = 0; r < runs.size(); ++r) {
        const std::size_t end = r + 1 < runs.size() ? runs[r + 1].first : text.size();
        const auto part = substituteRun(f, text, runs[r].first, end, runs[r].second);
        const std::size_t first = glyphs.size();
        for (const auto& [g, cluster] : part) {
            glyphs.push_back(g);
            clusters.push_back(cluster);
        }
        kernAdvance.resize(glyphs.size(), 0);
        kernPlacement.resize(glyphs.size(), 0);
        kernRun(f, glyphs, first, glyphs.size(), runs[r].second, kernAdvance, kernPlacement);
    }

    TextLayout out;
    out.pixelSize = size.pixelSize();
    std::int64_t pen = 0;
    for (std::size_t i = 0; i < glyphs.size(); ++i) {
        // FT_Get_Advance (unhinted: units x x_scale / 64 in 16.16) -> 26.6, as hb-ft rounds it
        const std::int64_t advance = ((mulDiv(f.advance(glyphs[i]), size.scale, 64) + (1 << 9)) >> 10) + em(kernAdvance[i]);
        const std::int64_t x = pen + em(kernPlacement[i]);
        const GlyphBox b = f.box(glyphs[i], size.scale);
        if (!b.empty) {
            const std::int64_t x0 = x + b.xMin, x1 = x + b.xMax, y0 = b.yMin, y1 = b.yMax;
            if (!out.hasInk) {
                out.inkXMin = x0, out.inkXMax = x1, out.inkYMin = y0, out.inkYMax = y1;
                out.hasInk = true;
            } else {
                out.inkXMin = std::min(out.inkXMin, x0), out.inkXMax = std::max(out.inkXMax, x1);
                out.inkYMin = std::min(out.inkYMin, y0), out.inkYMax = std::max(out.inkYMax, y1);
            }
        }
        out.glyphs.push_back({glyphs[i], x, advance, clusters[i]});
        pen += advance;
    }
    out.advance = pen;
    return out;
}

std::vector<PathCommand> textPath(const Font& font, const TextLayout& layout) {
    using K = PathCommand::Kind;
    struct P {
        double x, y;
        bool on;
    };
    std::vector<PathCommand> out;
    for (const PositionedGlyph& g : layout.glyphs) {
        const std::vector<OutlinePoint> pts = font.outline(g.glyph, layout.pixelSize);
        std::vector<P> contour;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            contour.push_back({static_cast<double>(g.x + pts[i].x) / 64.0, -static_cast<double>(pts[i].y) / 64.0,
                               pts[i].onCurve});
            if (!pts[i].contourEnd && i + 1 < pts.size()) {
                continue;
            }
            const std::size_t n = contour.size();
            if (n >= 2) {
                std::size_t first = 0;
                while (first < n && !contour[first].on) {
                    ++first;
                }
                P start{};
                std::vector<P> rest;  // the points after the start, in order
                if (first < n) {
                    start = contour[first];
                    for (std::size_t k = 1; k < n; ++k) {
                        rest.push_back(contour[(first + k) % n]);
                    }
                } else {  // only off-curve points: start midway between the last and the first
                    start = {(contour[n - 1].x + contour[0].x) / 2, (contour[n - 1].y + contour[0].y) / 2, true};
                    rest = contour;
                }
                out.push_back({K::Move, start.x, start.y});
                std::optional<P> control;
                for (const P& q : rest) {
                    if (q.on) {
                        out.push_back(control ? PathCommand{K::Quad, q.x, q.y, control->x, control->y}
                                              : PathCommand{K::Line, q.x, q.y});
                        control.reset();
                    } else {
                        if (control) {
                            out.push_back({K::Quad, (control->x + q.x) / 2, (control->y + q.y) / 2, control->x, control->y});
                        }
                        control = q;
                    }
                }
                if (control) {
                    out.push_back({K::Quad, start.x, start.y, control->x, control->y});
                }
                out.push_back({K::Close});
            }
            contour.clear();
        }
    }
    return out;
}

}  // namespace mdbridge::render
