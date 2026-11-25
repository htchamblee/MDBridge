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
    {FontFace::Sans, "DejaVuSans.ttf", "DejaVu Sans", "Book", 400, false, false,
     "3fdf69cabf06049ea70a00b5919340e2ce1e6d02b0cc3c4b44fb6801bd1e0d22"},
    {FontFace::SansBold, "DejaVuSans-Bold.ttf", "DejaVu Sans", "Bold", 700, false, false,
     "b184b89e3c1075f22f6b71575b6fc20d4972b3cfd3b23322ca6fd596dcaef167"},
    {FontFace::SansOblique, "DejaVuSans-Oblique.ttf", "DejaVu Sans", "Oblique", 400, true, false,
     "ccdf74b350f11fd3dd5774de50e5e6346a1a5da1f5b7d5fb83590665e97a5213"},
    {FontFace::Mono, "DejaVuSansMono.ttf", "DejaVu Sans Mono", "Book", 400, false, true,
     "602ec86b8948cfcd956482fe64f94c36c867770149ef2f791d4613f443bcecb3"},
    {FontFace::MonoBold, "DejaVuSansMono-Bold.ttf", "DejaVu Sans Mono", "Bold", 700, false, true,
     "baada9a5172fe20886251aff0433fc38461912d5daf07287e7bee56620a8da96"},
}};

std::string lower(std::string_view s) {
    std::string out(s);
    for (char& c : out) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return out;
}

// Matplotlib's rcParams defaults for the generic families (font.sans-serif, ...)
const std::vector<std::string_view>* aliasList(std::string_view generic) {
    static const std::vector<std::string_view> sans{
        "DejaVu Sans", "Bitstream Vera Sans", "Computer Modern Sans Serif", "Lucida Grande", "Verdana", "Geneva",
        "Lucid", "Arial", "Helvetica", "Avant Garde", "sans-serif"};
    static const std::vector<std::string_view> mono{
        "DejaVu Sans Mono", "Bitstream Vera Sans Mono", "Computer Modern Typewriter", "Andale Mono",
        "Nimbus Mono L", "Courier New", "Courier", "Fixed", "Terminal", "monospace"};
    static const std::vector<std::string_view> serif{
        "DejaVu Serif", "Bitstream Vera Serif", "Computer Modern Roman", "New Century Schoolbook",
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
    FontFace face = FontFace::Sans;
    for (const BundledFont& f : kFonts) {
        const double score = scoreFamily(families, f.family) * 10 + scoreStyle(style, f.oblique ? "oblique" : "normal") +
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
        const std::string fallback[] = {"DejaVu Sans"};
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

FontSet::FontSet(std::array<std::string, kFontFaceCount> data) {
    for (const BundledFont& f : kFonts) {
        std::string& bytes = data[static_cast<std::size_t>(f.face)];
        if (bytes.empty()) {
            throw FontError("bundled font " + std::string(f.file) + " is missing");
        }
        if (core::sha256Hex(bytes) != f.sha256) {
            throw FontError("bundled font " + std::string(f.file) + " is not DejaVu " + std::string(kDejaVuVersion) +
                            " as Matplotlib 3.11.2 ships it (sha256 differs)");
        }
        Font font = Font::fromData(std::move(bytes));
        if (font.family() != f.family || font.subfamily() != f.style) {
            throw FontError("bundled font " + std::string(f.file) + " is " + font.family() + " " + font.subfamily());
        }
        fonts_.push_back(std::move(font));
    }
}

std::optional<std::string_view> batchQtPlatform() {
#if defined(_WIN32) || defined(__APPLE__)
    return std::nullopt;
#else
    return std::string_view("offscreen");
#endif
}

bool drawsWithoutWindow(std::string_view command) {
    return command == "render" || command == "charts" || command == "preview" || command == "fonts";
}

}  // namespace mdbridge::render
