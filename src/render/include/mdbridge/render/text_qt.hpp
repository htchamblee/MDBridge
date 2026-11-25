#pragma once
// The figure fonts in the Qt programs: the bundled DejaVu faces compiled in as Qt resources
// (:/mdbridge/fonts/, from data/fonts), text drawn from their outlines, QFonts of the same faces, and
// the Qt application the batch commands draw with.
//
// Raster figures draw text as filled paths of the glyph outlines that render core lays out
// (font.hpp): positions and shapes are the same numbers on every system, and the pixels come from
// Qt's own path rasteriser, not from a font engine (FreeType / DirectWrite), its hinting or the
// system's font settings, so a text image is identical on a bare server and on a desktop. QFonts of
// the bundled faces (registered with QFontDatabase) serve the interface and vector output that keeps
// text as text.

#include <memory>
#include <string>

#include <QColor>
#include <QFont>
#include <QImage>
#include <QPainterPath>
#include <QPointF>
#include <QStringList>

#include "mdbridge/render/text.hpp"

class QGuiApplication;
class QPainter;
class QPaintDevice;

namespace mdbridge::render {

// The bundled faces read from the program's resources (loaded and checked once; thread-safe). Throws
// FontError when the resources are missing or not the bundled files.
const FontSet& bundledFontSet();

// Registers the bundled faces with QFontDatabase (once; needs a QGuiApplication) and returns the
// family names Qt gives them ("DejaVu Sans", "DejaVu Sans Mono"). Throws FontError when Qt rejects
// one.
QStringList registerBundledFonts();

// A QFont of a bundled face at `pixels` pixels per em on `device` (point size from the device's
// logical dpi), unhinted, without fallback to other fonts. registerBundledFonts first.
QFont figureFont(FontFace face, double pixels, const QPaintDevice& device);

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
