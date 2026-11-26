// PostScript / EPS and TIFF output of figures (figure_formats.hpp). Written from the PostScript
// Language Reference (3rd edition), the EPSF 3.0 and DSC 3.0 specifications and the TIFF 6.0
// specification.

#include "mdbridge/render/figure_formats.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "mdbridge/render/font.hpp"

namespace mdbridge::render {

namespace {

// A number for the program: fixed point, trailing zeros dropped ("12.5", "-3", "0")
std::string num(double v) {
    if (!std::isfinite(v)) v = 0.0;
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.4f", v);
    std::string s = buf;
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    if (s == "-0") s = "0";
    return s;
}

// a DSC text value: printable ASCII, parentheses and backslashes escaped
std::string dscText(std::string_view text) {
    std::string out;
    for (unsigned char c : text) {
        if (c < 0x20 || c > 0x7e) continue;
        if (c == '(' || c == ')' || c == '\\') out.push_back('\\');
        out.push_back(static_cast<char>(c));
    }
    return out;
}

class PostScript {
public:
    explicit PostScript(const FontSet& fonts) : m_fonts(fonts) {}

    std::string out;

    // PostScript has no transparency: alpha 0 is not drawn, anything else opaque (Matplotlib's PS
    // backend does the same)
    bool color(const Rgba& c) {
        if (!(c.a > 0.0)) return false;
        out += num(std::clamp(c.r, 0.0, 1.0)) + " " + num(std::clamp(c.g, 0.0, 1.0)) + " " +
               num(std::clamp(c.b, 0.0, 1.0)) + " rg\n";
        return true;
    }

    void stroke(const DisplayStroke& s) {
        out += num(s.width) + " setlinewidth " + std::to_string(s.cap == CapStyle::Butt ? 0 : s.cap == CapStyle::Round ? 1 : 2) +
               " setlinecap " + (s.roundJoin ? "1" : "0") + " setlinejoin\n";
        out += "[";
        if (!s.dashes.empty()) {
            for (std::size_t i = 0; i < s.dashes.size(); ++i) out += (i ? " " : "") + num(std::max(s.dashes[i], 1e-6));
        }
        out += "] " + num(s.dashOffset) + " setdash\n";
    }

    void operator()(const FillRectItem& r) {
        if (!color(r.color)) return;
        out += num(r.rect.x0) + " " + num(r.rect.y0) + " " + num(r.rect.width()) + " " + num(r.rect.height()) + " rectfill\n";
    }

    void operator()(const PathItem& item) {
        auto path = [&] {
            out += "newpath\n";
            for (const auto& sub : item.subpaths) {
                if (sub.empty()) continue;
                out += num(sub[0].x) + " " + num(sub[0].y) + " m\n";
                for (std::size_t i = 1; i < sub.size(); ++i) out += num(sub[i].x) + " " + num(sub[i].y) + " l\n";
                if (item.closed) out += "cp\n";
            }
        };
        if (item.fill && item.fill->a > 0.0) {
            path();
            color(*item.fill);
            out += "fill\n";
        }
        if (item.stroke && item.stroke->width > 0 && item.stroke->color.a > 0.0) {
            path();
            color(item.stroke->color);
            stroke(*item.stroke);
            out += "stroke\n";
        }
    }

    void operator()(const MarkerItem& m) {
        for (std::size_t i = 0; i < m.centers.size(); ++i) {
            const double size = m.sizes.empty() ? 0.0 : m.sizes[i % m.sizes.size()];
            const Point c = m.centers[i];
            auto shape = [&] {
                out += "newpath\n";
                if (m.shape == Marker::Square) {
                    out += num(c.x - size / 2) + " " + num(c.y - size / 2) + " " + num(size) + " " + num(size) + " re\n";
                } else {
                    out += num(c.x + size / 2) + " " + num(c.y) + " m " + num(c.x) + " " + num(c.y) + " " + num(size / 2) +
                           " 0 360 arc cp\n";
                }
            };
            if (!m.faces.empty() && m.faces[i % m.faces.size()].a > 0.0) {
                shape();
                color(m.faces[i % m.faces.size()]);
                out += "fill\n";
            }
            if (!m.edges.empty() && m.edgeWidth > 0 && m.edges[i % m.edges.size()].a > 0.0) {
                shape();
                color(m.edges[i % m.edges.size()]);
                out += num(m.edgeWidth) + " setlinewidth 0 setlinejoin [] 0 setdash stroke\n";
            }
        }
    }

    void operator()(const RoundRectItem& r) {
        auto shape = [&] {
            const double rad = std::min({r.radius, r.rect.width() / 2, r.rect.height() / 2});
            const Box& b = r.rect;
            out += "newpath " + num(b.x0 + rad) + " " + num(b.y0) + " m\n";
            out += num(b.x1) + " " + num(b.y0) + " " + num(b.x1) + " " + num(b.y1) + " " + num(rad) + " arct\n";
            out += num(b.x1) + " " + num(b.y1) + " " + num(b.x0) + " " + num(b.y1) + " " + num(rad) + " arct\n";
            out += num(b.x0) + " " + num(b.y1) + " " + num(b.x0) + " " + num(b.y0) + " " + num(rad) + " arct\n";
            out += num(b.x0) + " " + num(b.y0) + " " + num(b.x1) + " " + num(b.y0) + " " + num(rad) + " arct cp\n";
        };
        if (r.fill && r.fill->a > 0.0) {
            shape();
            color(*r.fill);
            out += "fill\n";
        }
        if (r.stroke && r.stroke->width > 0 && r.stroke->color.a > 0.0) {
            shape();
            color(r.stroke->color);
            stroke(*r.stroke);
            out += "stroke\n";
        }
    }

    void operator()(const TextItem& t) {
        if (!(t.color.a > 0.0)) return;
        const std::vector<PathCommand> cmds = m_fonts.path(t.face, t.layout);
        if (cmds.empty()) return;
        out += "gsave " + num(t.origin.x) + " " + num(t.origin.y) + " translate";
        if (t.rotation != 0.0) out += " " + num(-t.rotation) + " rotate";
        out += "\nnewpath\n";
        double x = 0, y = 0;  // the current point (quadratic -> cubic Bézier)
        for (const PathCommand& c : cmds) {
            switch (c.kind) {
            case PathCommand::Kind::Move:
                out += num(c.x) + " " + num(c.y) + " m\n";
                break;
            case PathCommand::Kind::Line:
                out += num(c.x) + " " + num(c.y) + " l\n";
                break;
            case PathCommand::Kind::Quad:
                out += num(x + 2.0 / 3.0 * (c.cx - x)) + " " + num(y + 2.0 / 3.0 * (c.cy - y)) + " " +
                       num(c.x + 2.0 / 3.0 * (c.cx - c.x)) + " " + num(c.y + 2.0 / 3.0 * (c.cy - c.y)) + " " + num(c.x) +
                       " " + num(c.y) + " c\n";
                break;
            case PathCommand::Kind::Close:
                out += "cp\n";
                break;
            }
            if (c.kind != PathCommand::Kind::Close) {
                x = c.x;
                y = c.y;
            }
        }
        color(t.color);
        out += "fill grestore\n";
    }

    void operator()(const ImageItem& im) {
        if (im.width <= 0 || im.height <= 0) return;
        // the unit square onto the image's corners (user space is y down, as the image's rows)
        out += "gsave " + num(im.corner0.x) + " " + num(im.corner0.y) + " translate " + num(im.corner1.x - im.corner0.x) +
               " " + num(im.corner1.y - im.corner0.y) + " scale\n";
        out += "/DeviceRGB setcolorspace\n<< /ImageType 1 /Width " + std::to_string(im.width) + " /Height " +
               std::to_string(im.height) + " /BitsPerComponent 8 /Decode [0 1 0 1 0 1] /ImageMatrix [" +
               std::to_string(im.width) + " 0 0 " + std::to_string(im.height) + " 0 0] /Interpolate " +
               (im.smooth ? "true" : "false") + " /DataSource currentfile /ASCIIHexDecode filter >> image\n";
        static const char hex[] = "0123456789abcdef";
        std::size_t column = 0;
        const std::size_t pixels = std::min(im.rgba.size() / 4, static_cast<std::size_t>(im.width) * static_cast<std::size_t>(im.height));
        for (std::size_t p = 0; p < 4 * pixels; p += 4) {
            const double a = im.rgba[p + 3] / 255.0;
            for (int ch = 0; ch < 3; ++ch) {
                // over white
                const double v = im.rgba[p + static_cast<std::size_t>(ch)] * a + 255.0 * (1.0 - a);
                const int b = std::clamp(static_cast<int>(std::lround(v)), 0, 255);
                out.push_back(hex[b >> 4]);
                out.push_back(hex[b & 15]);
            }
            if (++column == 40) {
                out.push_back('\n');
                column = 0;
            }
        }
        out += ">\ngrestore\n";
    }

private:
    const FontSet& m_fonts;
};

}  // namespace

std::string renderPostScript(const DisplayList& list, const FontSet& fonts, bool eps, std::string_view title) {
    const double k = 72.0 / list.dpi;  // pixels -> points
    const double wpt = list.width * k, hpt = list.height * k;
    PostScript ps(fonts);
    std::string& o = ps.out;
    o += eps ? "%!PS-Adobe-3.0 EPSF-3.0\n" : "%!PS-Adobe-3.0\n";
    if (!title.empty()) o += "%%Title: " + dscText(title) + "\n";
    o += "%%Creator: MDBridge\n";
    o += "%%BoundingBox: 0 0 " + std::to_string(static_cast<long long>(std::ceil(wpt - 1e-9))) + " " +
         std::to_string(static_cast<long long>(std::ceil(hpt - 1e-9))) + "\n";
    o += "%%HiResBoundingBox: 0 0 " + num(wpt) + " " + num(hpt) + "\n";
    o += "%%LanguageLevel: 2\n%%Pages: 1\n";
    if (!eps) o += "%%DocumentMedia: figure " + num(wpt) + " " + num(hpt) + " 0 () ()\n";
    o += "%%EndComments\n%%BeginProlog\n";
    o += "/m {moveto} bind def /l {lineto} bind def /c {curveto} bind def /cp {closepath} bind def\n";
    o += "/re {4 2 roll moveto 1 index 0 rlineto 0 exch rlineto neg 0 rlineto closepath} bind def\n";
    o += "/rg {setrgbcolor} bind def\n%%EndProlog\n";
    if (!eps) o += "%%BeginSetup\n<< /PageSize [" + num(wpt) + " " + num(hpt) + "] >> setpagedevice\n%%EndSetup\n";
    o += "%%Page: 1 1\ngsave\n";
    // figure pixels, origin top left, y down
    o += "0 " + num(hpt) + " translate " + num(k) + " " + num(-k) + " scale\n";
    for (const DisplayItem& item : list.items) {
        if (item.clip) {
            const Box& c = *item.clip;
            o += "gsave " + num(c.x0) + " " + num(c.y0) + " " + num(c.width()) + " " + num(c.height()) + " rectclip\n";
        }
        std::visit(ps, item.item);
        if (item.clip) o += "grestore\n";
    }
    o += "grestore\nshowpage\n%%EOF\n";
    return std::move(ps.out);
}

std::string encodeTiff(int width, int height, std::span<const std::uint8_t> rgba, double dpi) {
    if (width <= 0 || height <= 0 || rgba.size() < static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4) {
        throw std::invalid_argument("TIFF: the pixels do not match the size");
    }
    const std::uint32_t bytes = static_cast<std::uint32_t>(static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4);
    std::string out;
    auto u16 = [&](std::uint32_t v) {
        out.push_back(static_cast<char>(v & 0xff));
        out.push_back(static_cast<char>((v >> 8) & 0xff));
    };
    auto u32 = [&](std::uint32_t v) {
        u16(v & 0xffff);
        u16(v >> 16);
    };
    struct Entry {
        std::uint16_t tag, type;
        std::uint32_t count, value;
    };
    constexpr std::uint16_t kShort = 3, kLong = 4, kRational = 5;
    constexpr std::uint32_t kEntries = 14;
    const std::uint32_t ifd = 8;
    const std::uint32_t extra = ifd + 2 + kEntries * 12 + 4;  // after the IFD: bits per sample, resolutions
    const std::uint32_t bitsAt = extra, xResAt = extra + 8, yResAt = extra + 16, dataAt = extra + 24;
    const std::uint32_t denominator = 10000;
    const std::uint32_t numerator = static_cast<std::uint32_t>(std::clamp(std::lround(dpi * denominator), 1L, 2000000000L));
    const Entry entries[kEntries] = {
        {256, kLong, 1, static_cast<std::uint32_t>(width)},   // ImageWidth
        {257, kLong, 1, static_cast<std::uint32_t>(height)},  // ImageLength
        {258, kShort, 4, bitsAt},                             // BitsPerSample 8, 8, 8, 8
        {259, kShort, 1, 1},                                  // Compression: none
        {262, kShort, 1, 2},                                  // PhotometricInterpretation: RGB
        {273, kLong, 1, dataAt},                              // StripOffsets
        {274, kShort, 1, 1},                                  // Orientation: top left
        {277, kShort, 1, 4},                                  // SamplesPerPixel
        {278, kLong, 1, static_cast<std::uint32_t>(height)},  // RowsPerStrip
        {279, kLong, 1, bytes},                               // StripByteCounts
        {282, kRational, 1, xResAt},                          // XResolution
        {283, kRational, 1, yResAt},                          // YResolution
        {296, kShort, 1, 2},                                  // ResolutionUnit: inch
        {338, kShort, 1, 2},                                  // ExtraSamples: unassociated alpha
    };
    out += "II";
    u16(42);
    u32(ifd);
    u16(kEntries);
    for (const Entry& e : entries) {
        u16(e.tag);
        u16(e.type);
        u32(e.count);
        if (e.type == kShort && e.count == 1) {
            u16(e.value);
            u16(0);
        } else {
            u32(e.value);
        }
    }
    u32(0);  // no next IFD
    for (int i = 0; i < 4; ++i) u16(8);
    u32(numerator);
    u32(denominator);
    u32(numerator);
    u32(denominator);
    out.append(reinterpret_cast<const char*>(rgba.data()), bytes);
    return out;
}

}  // namespace mdbridge::render
