#pragma once
// The QPainter backends of the figure framework: a display list (display_list.hpp) painted into a
// raster image, an SVG document or a PDF page of exactly the figure's size.
//
//   raster  QImage, RGBA8888 (straight alpha), width x height pixels; text filled from the bundled
//           fonts' outlines (text_qt.hpp): the same pixels on every system
//   SVG     QSvgGenerator, width / height = the figure's size in inches at its dpi, view box in
//           pixels; text kept as text (<text> elements in the bundled IBM Plex families: one per
//           run of the figure face or of its fallback, IBM Plex Math)
//   PDF     QPdfWriter, one page of the figure's size in points, no margins; text kept as text
//           (the bundled faces embedded by Qt)
//   EPS, PS, TIFF  render core's writers (figure_formats.hpp)
//
// SVG and PDF need a QGuiApplication (fonts); the batch commands start one on the offscreen
// platform (startBatchApplication). No Qt Charts / Qt Graphs (GPL): plain QPainter.

#include <filesystem>
#include <string>

#include <QByteArray>
#include <QImage>

#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/figure.hpp"

class QPainter;

namespace mdbridge::render {

enum class TextOutput {
    Outlines,  // glyph outlines filled as paths (raster images)
    Fonts,     // QFont text (vector output keeps text as text)
};

// Paints `list` with `painter`, whose coordinates are the list's pixels.
void paintDisplayList(QPainter& painter, const DisplayList& list, TextOutput text);

// The list as an image of exactly list.width x list.height pixels (QImage::Format_RGBA8888).
QImage renderImage(const DisplayList& list);
// Composes the frame loop's layers: a copy of `background` with `layer` painted over it.
QImage composeLayers(const QImage& background, const DisplayList& layer);

// SVG / PDF documents of the list.
QByteArray renderSvg(const DisplayList& list, const QString& title = {});
QByteArray renderPdf(const DisplayList& list, const QString& title = {});

// The image formats a figure can be written as: png, svg, pdf, eps, ps, tif / tiff (render core's
// writers, figure_formats.hpp), and the raster formats Qt's image plugins write here (jpg / jpeg
// with the qjpeg plugin, bmp ...), lower-case.
bool canWriteFigureFormat(const std::string& format);
// Writes the list to `path` in the format of its extension (png, svg, pdf, jpg ...), or in
// `format` when one is given; throws std::runtime_error naming the file when the format is unknown
// or the file cannot be written. Raster files record the figure's dpi; JPEG has no alpha (the
// figure over white, as Matplotlib's savefig).
void writeFigure(const DisplayList& list, const std::filesystem::path& path, const std::string& format = {});
// A drawn figure (Figure::draw) written as writeFigure writes its display list, text from the
// bundled fonts: the charts' writer (render/panels charts.hpp FigureWriter).
void writeFigure(const Figure& figure, const std::filesystem::path& path, const std::string& format = {});

}  // namespace mdbridge::render
