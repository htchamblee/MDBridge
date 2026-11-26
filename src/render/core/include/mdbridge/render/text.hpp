#pragma once
// The figure fonts: the IBM Plex faces MDBridge bundles (data/fonts, unmodified files of IBM's
// releases, data/sources/fonts/README.md): IBM Plex Serif for all text, IBM Plex Math as the one
// fallback for the characters Plex Serif lacks (superscript + / -, ...), IBM Plex Mono for the
// interface's fixed-width text only; how Matplotlib's font properties choose among the Serif faces,
// point / pixel conversion, and the Qt platform the batch commands use. No Qt: the Qt side
// (text_qt.hpp) registers the same files from the program's resources and draws the glyph runs of
// font.hpp.
//
// Figures use the bundled faces only, never the system's fonts: the same metrics, layout and look on
// every machine, also on systems without any font installed (a bare cluster node). MDATK draws in
// Matplotlib's DejaVu Sans; MDBridge in IBM Plex Serif, whose I, l and 1 cannot be confused
// (docs/DIFFERENCES.md, "Rendering: fonts and text"); the parity generators give Matplotlib the same
// files and font lists (tests/parity/improvements_fonts.py).

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/render/font.hpp"

namespace mdbridge::render {

// The figure faces: IBM Plex Serif Regular (labels, ticks, titles, legends, MDATK's "monospace"
// sequence letters), Bold, Italic (the sequence panels' chain notes) and Bold Italic.
enum class FontFace : std::uint8_t { Regular, Bold, Italic, BoldItalic };
inline constexpr std::size_t kFontFaceCount = 4;

// A bundled font file: data/fonts/<file>, Qt resource kFontResourcePrefix + file.
struct BundledFile {
    std::string_view file;
    std::string_view family;   // the font's family name (name ID 1, "IBM Plex Serif")
    std::string_view style;    // its style name (name ID 2: "Regular", "Bold", "Italic", "Bold Italic")
    std::string_view version;  // name ID 5 ("Version 3.006")
    std::string_view sha256;   // of the file as IBM's release has it (data/sources/fonts/README.md)
};

struct BundledFont : BundledFile {
    FontFace face = FontFace::Regular;
    int weight = 400;  // Matplotlib / CSS weight: 400 or 700
    bool italic = false;
};

// The figure faces in registry order (also the order the font choice breaks ties in).
std::span<const BundledFont> bundledFonts();
const BundledFont& bundledFont(FontFace face);
// IBM Plex Math: the fallback of every figure face (and of the interface's fonts) for the
// characters they have no glyph for.
const BundledFile& fallbackFont();
// IBM Plex Mono Regular and Bold: the interface's fixed-width text (logs, multi-line fields); not
// used in figures.
std::span<const BundledFile> interfaceMonoFonts();
inline constexpr std::string_view kFontResourcePrefix = ":/mdbridge/fonts/";
// IBM's release tags of the bundled files (github.com/IBM/plex)
inline constexpr std::string_view kPlexSerifRelease = "@ibm/plex-serif@2.0.0";
inline constexpr std::string_view kPlexMathRelease = "@ibm/plex-math@1.1.0";
inline constexpr std::string_view kPlexMonoRelease = "@ibm/plex-mono@2.5.0";
// The family the figures use whatever family a text asks for (Matplotlib's generic families
// sans-serif, serif and monospace list it first; unknown families fall back to it).
inline constexpr std::string_view kFigureFamily = "IBM Plex Serif";

// Matplotlib's weight names (font_manager.weight_dict: "normal" 400, "bold" 700, "light" 200 ...)
// or a number 1..1000 as text; nullopt for anything else.
std::optional<int> fontWeight(std::string_view name);

// The bundled face Matplotlib's font manager (FontManager.findfont) picks for font properties when
// only the bundled faces are installed and its font lists are MDBridge's (IBM Plex Serif first in
// font.sans-serif, font.serif and font.monospace; tests/parity/improvements_fonts.py): its score over
// family (generic names sans-serif / sans / monospace / serif / cursive / fantasy, or names;
// case-insensitive), style ("normal", "italic", "oblique") and weight, the first face with the lowest
// score, and its fallback to IBM Plex Serif when no family matches.
FontFace chooseFace(std::span<const std::string> families, std::string_view style, int weight);
FontFace chooseFace(std::string_view family, std::string_view style, int weight);

// Points <-> pixels at a figure's dpi, as Matplotlib converts (RendererBase.points_to_pixels:
// points * dpi / 72).
double pointsToPixels(double points, double dpi);
double pixelsToPoints(double pixels, double dpi);
// The point size to give a QFont for `pixels` on a paint device of `deviceDpi` logical dpi (Qt turns
// a point size into pixels with the device's dpi).
double qtPointSize(double pixels, double deviceDpi);

// The figure faces and their fallback parsed (font.hpp), from bytes the caller provides: the Qt
// resources in the programs, data/fonts in tests. Throws FontError when a file is missing, malformed
// or not the bundled file (family / style / sha256 differ).
class FontSet {
public:
    // `load(file)` returns the bytes of a bundled file (BundledFile::file).
    template <typename Loader>
    static FontSet load(Loader&& loader) {
        std::array<std::string, kFontFaceCount> data;
        for (const BundledFont& f : bundledFonts()) {
            data[static_cast<std::size_t>(f.face)] = loader(f.file);
        }
        return FontSet(std::move(data), loader(fallbackFont().file));
    }
    const Font& font(FontFace face) const { return fonts_[static_cast<std::size_t>(face)]; }
    const Font& fallback() const { return fonts_.back(); }

    // Matplotlib's get_text_width_height_descent for one line of text (font.hpp, unhinted), the
    // characters the face lacks from the fallback font.
    TextLayout layout(FontFace face, std::string_view utf8, double points, double dpi) const {
        return font(face).layout(utf8, points, dpi, &fallback());
    }
    // The outlines of a layout of `face` (fallback glyphs from the fallback font).
    std::vector<PathCommand> path(FontFace face, const TextLayout& line) const {
        return textPath(font(face), line, &fallback());
    }

private:
    FontSet(std::array<std::string, kFontFaceCount> data, std::string fallbackData);
    std::vector<Font> fonts_;  // the faces in FontFace order, then the fallback
};

// Throws FontError when `bytes` are not the bundled file `f` (sha256, family, style); returns the
// parsed font.
Font checkedFont(const BundledFile& f, std::string bytes);

// The Qt platform plugin (QT_QPA_PLATFORM) for commands that draw without a window (`render`,
// `charts`, `preview`, `fonts`): "offscreen" on Linux and other X11 / Wayland systems, so they run
// without a display (ssh, batch jobs, containers) whatever DISPLAY / WAYLAND_DISPLAY or a desktop's
// QT_QPA_PLATFORM say; nullopt on Windows and macOS, whose native platforms need no display server.
std::optional<std::string_view> batchQtPlatform();
// The commands that draw without a window (also `demo`, for `demo --run`).
bool drawsWithoutWindow(std::string_view command);

}  // namespace mdbridge::render
