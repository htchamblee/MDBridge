#pragma once
// The tables of mpl_colors.inc (tools/gen_mpl_tables.py: Matplotlib 3.11.2's named colours,
// default colour cycle and bundled colormaps; the notices are in the .inc).

#include <iterator>
#include <string_view>

namespace mdbridge::render::detail {

struct Rgb {
    double r, g, b;
};
struct NamedColor {
    std::string_view name;
    Rgb rgb;
};
struct ListedMapRow {
    std::string_view name;
    int first;  // row in kListedMapColors
    int count;  // Matplotlib's N
};

#include "mpl_colors.inc"

}  // namespace mdbridge::render::detail
