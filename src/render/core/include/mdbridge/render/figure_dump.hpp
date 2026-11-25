#pragma once
// A drawn figure (figure.hpp) as data, in the record structure tests/parity/mpl_walker.py writes for
// MDATK's Matplotlib figures, so that tests compare MDBridge's figures with MDATK's like with like
// (owner decision D10: contents only, no pixels):
//
//   figure   size_px, dpi, facecolor, suptitle, texts, axes
//   axes     id "ax<k>", position_px, tight_px, facecolor, visible, axison, xlim, ylim, aspect,
//            zorder, titles {left, center, right} (non-empty ones), spines {side: visible, color,
//            lw}, xaxis / yaxis, artists (each with id "ax<k>/<n>"; the legend as "ax<k>/legend"),
//            panel / role when set, colorbar {orientation, vmin, vmax, cmap {name, N}, outline,
//            label}
//   axis     label {text, color, size, extent_px}, ticks [[loc, label, label1 shown, label2
//            shown]], tick_style, offset_text, visible, inverted, locator, formatter
//   artists  line, scatter, segments, image, text, patch (Rectangle), legend as mpl_walker.py
//            describes them
//
// Colours are RGBA floats in [0, 1]; positions in figure pixels, origin bottom left.

#include "mdbridge/core/value.hpp"
#include "mdbridge/render/figure.hpp"

namespace mdbridge::render {

// The figure after Figure::draw(fonts).
core::Value dumpFigure(const Figure& figure, const FontSet& fonts);

}  // namespace mdbridge::render
