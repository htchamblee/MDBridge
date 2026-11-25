#pragma once
// Figure file formats written without Qt (the QPainter backends, figure_qt.hpp, write PNG, JPEG,
// SVG and PDF; these complete the formats MDATK's charts write through Matplotlib and pillow):
//
//   PostScript  the display list (display_list.hpp) as a PostScript Level 2 program: an EPS file
//               (EPSF-3.0, %%BoundingBox of the figure in points) or a one-page PS document whose
//               page is the figure. Paths, markers and images as vectors and image data, text as
//               the filled outlines of the bundled fonts (font.hpp textPath: no fonts embedded or
//               needed). PostScript has no transparency: partially transparent colours are drawn
//               opaque and fully transparent ones not at all, as Matplotlib's PS backend does;
//               images are composited over white.
//   TIFF        RGBA pixels (straight alpha) as a baseline TIFF: little endian, one uncompressed
//               strip, 8 bits per sample, unassociated alpha, the resolution in dpi.

#include <cstdint>
#include <span>
#include <string>
#include <string_view>

#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/text.hpp"

namespace mdbridge::render {

std::string renderPostScript(const DisplayList& list, const FontSet& fonts, bool eps, std::string_view title = {});

// `rgba`: width x height x 4 bytes, row 0 (the top) first.
std::string encodeTiff(int width, int height, std::span<const std::uint8_t> rgba, double dpi);

}  // namespace mdbridge::render
