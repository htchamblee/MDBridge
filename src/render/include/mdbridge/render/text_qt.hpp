#pragma once
// The bundled fonts in the Qt programs: the IBM Plex faces compiled in as Qt resources
// (:/mdbridge/fonts/, from data/fonts), figure text drawn from their outlines, QFonts of the same
// faces (figures' vector output and the interface), and the Qt application the batch commands draw
// with.
//
// Raster figures draw text as filled paths of the glyph outlines that render core lays out
// (font.hpp): positions and shapes are the same numbers on every system, and the pixels come from
// Qt's own path rasteriser, not from a font engine (FreeType / DirectWrite), its hinting or the
// system's font settings, so a text image is identical on a bare server and on a desktop. QFonts of
// the bundled faces (registered with QFontDatabase) serve the interface and vector output that keeps
// text as text; IBM Plex Math is their fallback for missing glyphs (in figures per glyph as
// Matplotlib falls back, in the interface as the next family of Qt's font family list).

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <QColor>
#include <QFont>
#include <QImage>
#include <QPainterPath>
#include <QPointF>
#include <QString>
#include <QStringList>

#include "mdbridge/render/text.hpp"

class QGuiApplication;
class QPainter;
class QPaintDevice;

namespace mdbridge::render {

// The bundled faces read from the program's resources (loaded and checked once; thread-safe). Throws
// FontError when the resources are missing or not the bundled files.
const FontSet& bundledFontSet();

// Registers the bundled fonts with QFontDatabase (once; needs a QGuiApplication): the figure faces,
// the fallback and the interface's fixed-width faces; returns the family names Qt gives them ("IBM
// Plex Serif", "IBM Plex Math", "IBM Plex Mono"). Throws FontError when one is missing, not the
// bundled file or rejected by Qt.
QStringList registerBundledFonts();

// A QFont of a figure face at `pixels` pixels per em on `device` (point size from the device's
// logical dpi), unhinted, without fallback to other fonts; with `fallbackGlyphs`, the fallback font
// (IBM Plex Math Regular) for the glyphs a layout takes from it. registerBundledFonts first.
QFont figureFont(FontFace face, double pixels, const QPaintDevice& device, bool fallbackGlyphs = false);

// The interface's fonts (registerBundledFonts first): IBM Plex Serif, and IBM Plex Mono for
// fixed-width text (logs, multi-line fields), each with IBM Plex Math as the next family for the
// characters it lacks; at `pointSize` points (Qt's default size when <= 0).
QFont interfaceFont(double pointSize = 0);
QFont interfaceFixedFont(double pointSize = 0);

// The runs of a laid-out line drawn with one QFont each (vector output that keeps text as text): the
// pen position where the run starts (pixels), whether its glyphs are the fallback font's, and its
// characters [begin, end) as code point indices of the line.
struct TextRun {
    double x = 0;
    bool fallback = false;
    std::size_t begin = 0;
    std::size_t end = 0;
};
std::vector<TextRun> fontRuns(std::string_view utf8, const TextLayout& layout);
QString runText(std::string_view utf8, const TextRun& run);

// A laid-out line as a path, origin on the baseline at the left (textPath, non-zero winding).
QPainterPath textPainterPath(FontFace face, const TextLayout& layout);
// Fills that path at `origin` with `color`, antialiased (the painter's transform applies: rotated
// labels rotate the painter).
void drawText(QPainter& painter, QPointF origin, FontFace face, const TextLayout& layout, const QColor& color);

// The text check of `mdbridge fonts --render` and the tests: every bundled face at several sizes,
// rotated and coloured text, on white, 720 x 400 pixels (RGBA8888, premultiplied). Its pixels must
// be the same wherever the same build runs (CI compares the Linux build machine and a bare system).
QImage textSampleImage();
// SHA-256 of an image's pixels as RGBA8888 bytes, row by row (no padding, no file format).
std::string imagePixelsSha256(const QImage& image);

// The Qt application of a command that draws without a window (render, charts, preview, fonts): the
// "offscreen" platform where batchQtPlatform() says so (QT_QPA_PLATFORM is set for this process,
// whatever it was), so it runs without a display. Qt reads no options from the command line.
std::unique_ptr<QGuiApplication> startBatchApplication(const char* program);

}  // namespace mdbridge::render
