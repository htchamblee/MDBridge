#pragma once
// The figure fonts: the DejaVu faces MDBridge bundles (data/fonts, the files Matplotlib 3.11.2 ships,
// DejaVu 2.35), how Matplotlib's font properties choose among them, point / pixel conversion, and
// the Qt platform the batch commands use. No Qt: the Qt side (text_qt.hpp) registers the same files
// from the program's resources and draws the glyph runs of font.hpp.
//
// Figures use the bundled faces only, never the system's fonts: the same metrics, layout and look on
// every machine, also on systems without any font installed (a bare cluster node).

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/render/font.hpp"

namespace mdbridge::render {

// The bundled faces. MDATK's figures use DejaVu Sans (labels, ticks, titles, legends), its Oblique
// (italic notes of the sequence panels), DejaVu Sans Mono (sequence letters) and its Bold (the
// sequence panel's residue letters); DejaVu Sans Bold for bold text.
enum class FontFace : std::uint8_t { Sans, SansBold, SansOblique, Mono, MonoBold };
inline constexpr std::size_t kFontFaceCount = 5;

struct BundledFont {
    FontFace face;
    std::string_view file;     // data/fonts/<file>; Qt resource kFontResourcePrefix + file
    std::string_view family;   // the font's family name ("DejaVu Sans")
    std::string_view style;    // its style name ("Book", "Bold", "Oblique")
    int weight;                // Matplotlib / CSS weight: 400 or 700
    bool oblique;
    bool monospace;
    std::string_view sha256;   // of the file: Matplotlib 3.11.2's copy (tests/parity/make_text_cases.py)
};

// The faces in registry order (data/fonts' order, also the order the font choice breaks ties in).
std::span<const BundledFont> bundledFonts();
const BundledFont& bundledFont(FontFace face);
inline constexpr std::string_view kFontResourcePrefix = ":/mdbridge/fonts/";
inline constexpr std::string_view kDejaVuVersion = "2.35";

// Matplotlib's weight names (font_manager.weight_dict: "normal" 400, "bold" 700, "light" 200 ...)
// or a number 1..1000 as text; nullopt for anything else.
std::optional<int> fontWeight(std::string_view name);

// The bundled face Matplotlib's font manager (FontManager.findfont) picks for font properties when
// only the bundled faces are installed: its score over family (generic names sans-serif / sans /
// monospace / serif / cursive / fantasy with Matplotlib's default lists, or names; case-insensitive),
// style ("normal", "italic", "oblique") and weight, the first face with the lowest score, and its
// fallback to DejaVu Sans when no family matches. Bold oblique and mono oblique are not bundled:
// they give the closest bundled face, as Matplotlib's scoring would.
FontFace chooseFace(std::span<const std::string> families, std::string_view style, int weight);
FontFace chooseFace(std::string_view family, std::string_view style, int weight);

// Points <-> pixels at a figure's dpi, as Matplotlib converts (RendererBase.points_to_pixels:
// points * dpi / 72).
double pointsToPixels(double points, double dpi);
double pixelsToPoints(double pixels, double dpi);
// The point size to give a QFont for `pixels` on a paint device of `deviceDpi` logical dpi (Qt turns
// a point size into pixels with the device's dpi).
double qtPointSize(double pixels, double deviceDpi);

// The bundled faces parsed (font.hpp), from bytes the caller provides: the Qt resources in the
// programs, data/fonts in tests. Throws FontError when a face is missing, malformed or not the
// bundled file (family / style / sha256 differ).
class FontSet {
public:
    // `load(file)` returns the bytes of a bundled file (BundledFont::file).
    template <typename Loader>
    static FontSet load(Loader&& loader) {
        std::array<std::string, kFontFaceCount> data;
        for (const BundledFont& f : bundledFonts()) {
            data[static_cast<std::size_t>(f.face)] = loader(f.file);
        }
        return FontSet(std::move(data));
    }
    const Font& font(FontFace face) const { return fonts_[static_cast<std::size_t>(face)]; }

    // Matplotlib's get_text_width_height_descent for one line of text (font.hpp, unhinted).
    TextLayout layout(FontFace face, std::string_view utf8, double points, double dpi) const {
        return font(face).layout(utf8, points, dpi);
    }

private:
    explicit FontSet(std::array<std::string, kFontFaceCount> data);
    std::vector<Font> fonts_;
};

// The Qt platform plugin (QT_QPA_PLATFORM) for commands that draw without a window (`render`,
// `charts`, `preview`, `fonts`): "offscreen" on Linux and other X11 / Wayland systems, so they run
// without a display (ssh, batch jobs, containers) whatever DISPLAY / WAYLAND_DISPLAY or a desktop's
// QT_QPA_PLATFORM say; nullopt on Windows and macOS, whose native platforms need no display server.
std::optional<std::string_view> batchQtPlatform();
// The commands that draw without a window.
bool drawsWithoutWindow(std::string_view command);

}  // namespace mdbridge::render
