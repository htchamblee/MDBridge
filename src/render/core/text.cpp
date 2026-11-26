// The bundled figure fonts and Matplotlib's font choice (text.hpp).
//
// chooseFace follows Matplotlib's FontManager.findfont scoring (score_family, score_style,
// score_weight; variant, stretch and size score 0 for the bundled scalable faces) and the default
// font lists of its rcParams (font.sans-serif, font.monospace, ...). Matplotlib 3.11.2, Copyright (c)
// 2012- Matplotlib Development Team; Copyright (c) 2002-2011 John D. Hunter; Matplotlib licence
// (PSF-style), packaging/licenses/Matplotlib-LICENSE.txt. Changes: rewritten in C++ for the bundled
// faces.
#include "mdbridge/render/text.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <utility>

#include "mdbridge/core/hash.hpp"

namespace mdbridge::render {

namespace {

constexpr std::array<BundledFont, kFontFaceCount> kFonts{{
    {{"IBMPlexSerif-Regular.ttf", "IBM Plex Serif", "Regular", "Version 3.006",
      "6c490ccd38e856c674e0d163f3f40742a75058b7b935998dbfbc1745ba54e474"},
     FontFace::Regular, 400, false},
    {{"IBMPlexSerif-Bold.ttf", "IBM Plex Serif", "Bold", "Version 3.006",
      "5fa66833dc80adb6252c793adb555716ec10833d2b1e5b0235b68c72cd69adf0"},
     FontFace::Bold, 700, false},
    {{"IBMPlexSerif-Italic.ttf", "IBM Plex Serif", "Italic", "Version 3.006",
      "d28aa63cf618439712d6b459c5649af1d6d789cb911994f2357f1ce44128f163"},
     FontFace::Italic, 400, true},
    {{"IBMPlexSerif-BoldItalic.ttf", "IBM Plex Serif", "Bold Italic", "Version 3.006",
      "e33c9d0a692bb26c9c1ece4a35fd2a351317f5a60009d0e9945ae30d0297bf42"},
     FontFace::BoldItalic, 700, true},
    {{"IBMPlexMono-Regular.ttf", "IBM Plex Mono", "Regular", "Version 2.005",
      "7c6fbddca4b700be918f5f6183d9bd4464fa427fe435f0b480d77fe2bb8c5a43"},
     FontFace::Mono, 400, false, true},
    {{"IBMPlexMono-Bold.ttf", "IBM Plex Mono", "Bold", "Version 2.005",
      "74e5eedcfa4596497d34e19023cabdabd3a8c852b903007a5654a59591a72ffb"},
     FontFace::MonoBold, 700, false, true},
}};

constexpr BundledFile kFallback{"IBMPlexMath-Regular.ttf", "IBM Plex Math", "Regular", "Version 1.000",
                                "10d07a23cb322725f83ca51748a4c7a61cf14caf2c4d880d793c65d3a1402533"};

std::string lower(std::string_view s) {
    std::string out(s);
    for (char& c : out) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return out;
}

// Matplotlib's rcParams defaults for the generic families (font.sans-serif, ...), with IBM Plex
// Serif first in sans-serif and serif and IBM Plex Mono first in monospace (MDBridge's figures;
// tests/parity/improvements_fonts.py gives Matplotlib the same lists)
const std::vector<std::string_view>* aliasList(std::string_view generic) {
    static const std::vector<std::string_view> sans{
        kFigureFamily, "DejaVu Sans", "Bitstream Vera Sans", "Computer Modern Sans Serif", "Lucida Grande", "Verdana",
        "Geneva", "Lucid", "Arial", "Helvetica", "Avant Garde", "sans-serif"};
    static const std::vector<std::string_view> mono{
        kSequenceFamily, "DejaVu Sans Mono", "Bitstream Vera Sans Mono", "Computer Modern Typewriter", "Andale Mono",
        "Nimbus Mono L", "Courier New", "Courier", "Fixed", "Terminal", "monospace"};
    static const std::vector<std::string_view> serif{
        kFigureFamily, "DejaVu Serif", "Bitstream Vera Serif", "Computer Modern Roman", "New Century Schoolbook",
        "Century Schoolbook L", "Utopia", "ITC Bookman", "Bookman", "Nimbus Roman No9 L", "Times New Roman",
        "Times", "Palatino", "Charter", "serif"};
    static const std::vector<std::string_view> cursive{
        "Apple Chancery", "Textile", "Zapf Chancery", "Sand", "Script MT", "Felipa", "Comic Neue",
        "Comic Sans MS", "cursive"};
    static const std::vector<std::string_view> fantasy{
        "Chicago", "Charcoal", "Impact", "Western", "xkcd script", "fantasy"};
    if (generic == "sans-serif" || generic == "sans" || generic == "sans serif") return &sans;
    if (generic == "monospace") return &mono;
    if (generic == "serif") return &serif;
    if (generic == "cursive") return &cursive;
    if (generic == "fantasy") return &fantasy;
    return nullptr;
}

double scoreFamily(std::span<const std::string> families, std::string_view name) {
    if (families.empty()) {
        return 1.0;
    }
    const std::string family2 = lower(name);
    const double step = 1.0 / static_cast<double>(families.size());
    for (std::size_t i = 0; i < families.size(); ++i) {
        const std::string family1 = lower(families[i]);
        if (const auto* options = aliasList(family1)) {
            for (std::size_t k = 0; k < options->size(); ++k) {
                if (lower((*options)[k]) == family2) {
                    return (static_cast<double>(i) + static_cast<double>(k) / static_cast<double>(options->size())) * step;
                }
            }
        } else if (family1 == family2) {
            return static_cast<double>(i) * step;
        }
    }
    return 1.0;
}

double scoreStyle(std::string_view a, std::string_view b) {
    if (a == b) {
        return 0.0;
    }
    const auto slanted = [](std::string_view s) { return s == "italic" || s == "oblique"; };
    return slanted(a) && slanted(b) ? 0.1 : 1.0;
}

FontFace bestFace(std::span<const std::string> families, std::string_view style, int weight, double& best) {
    best = 1e64;
    FontFace face = FontFace::Regular;
    for (const BundledFont& f : kFonts) {
        const double score = scoreFamily(families, f.family) * 10 + scoreStyle(style, f.italic ? "italic" : "normal") +
                             0.95 * (std::abs(weight - f.weight) / 1000.0) + 0.05;
        if (score < best) {
            best = score;
            face = f.face;
        }
    }
    return face;
}

}  // namespace

std::span<const BundledFont> bundledFonts() { return kFonts; }

const BundledFont& bundledFont(FontFace face) { return kFonts[static_cast<std::size_t>(face)]; }

const BundledFile& fallbackFont() { return kFallback; }


std::optional<int> fontWeight(std::string_view name) {
    static constexpr std::pair<std::string_view, int> kWeights[] = {
        {"ultralight", 100}, {"light", 200}, {"normal", 400}, {"regular", 400}, {"book", 400},
        {"medium", 500}, {"roman", 500}, {"semibold", 600}, {"demibold", 600}, {"demi", 600},
        {"bold", 700}, {"heavy", 800}, {"extra bold", 800}, {"black", 900}};
    for (const auto& [n, w] : kWeights) {
        if (n == name) {
            return w;
        }
    }
    int value = 0;
    const auto [end, error] = std::from_chars(name.data(), name.data() + name.size(), value);
    if (error == std::errc() && end == name.data() + name.size() && value >= 1 && value <= 1000) {
        return value;
    }
    return std::nullopt;
}

FontFace chooseFace(std::span<const std::string> families, std::string_view style, int weight) {
    double best = 0;
    const FontFace face = bestFace(families, style, weight, best);
    if (best >= 10.0) {  // no family matched: Matplotlib retries with its default family
        const std::string fallback[] = {std::string(kFigureFamily)};
        return bestFace(fallback, style, weight, best);
    }
    return face;
}

FontFace chooseFace(std::string_view family, std::string_view style, int weight) {
    const std::string families[] = {std::string(family)};
    return chooseFace(families, style, weight);
}

double pointsToPixels(double points, double dpi) { return points * dpi / 72.0; }
double pixelsToPoints(double pixels, double dpi) { return pixels * 72.0 / dpi; }
double qtPointSize(double pixels, double deviceDpi) { return pixels * 72.0 / deviceDpi; }

Font checkedFont(const BundledFile& f, std::string bytes) {
    if (bytes.empty()) {
        throw FontError("bundled font " + std::string(f.file) + " is missing");
    }
    if (core::sha256Hex(bytes) != f.sha256) {
        throw FontError("bundled font " + std::string(f.file) + " is not " + std::string(f.family) + " " +
                        std::string(f.style) + " as IBM's release has it (sha256 differs)");
    }
    Font font = Font::fromData(std::move(bytes));
    if (font.family() != f.family || font.subfamily() != f.style) {
        throw FontError("bundled font " + std::string(f.file) + " is " + font.family() + " " + font.subfamily());
    }
    return font;
}

FontSet::FontSet(std::array<std::string, kFontFaceCount> data, std::string fallbackData) {
    for (const BundledFont& f : kFonts) {
        fonts_.push_back(checkedFont(f, std::move(data[static_cast<std::size_t>(f.face)])));
    }
    fonts_.push_back(checkedFont(kFallback, std::move(fallbackData)));
}

std::optional<std::string_view> batchQtPlatform() {
#if defined(_WIN32) || defined(__APPLE__)
    return std::nullopt;
#else
    return std::string_view("offscreen");
#endif
}

bool drawsWithoutWindow(std::string_view command) {
    return command == "render" || command == "charts" || command == "preview" || command == "fonts" ||
           command == "demo";  // demo --run renders
}

}  // namespace mdbridge::render
