#pragma once
// The rendering part of check() (docs/specs/m4-plan.md RN-11 ... RN-15): the settings that MDATK
// v0.0.4.1 only found wrong in the middle of a video, a chart run or a preview, after the analyses
// had run (a Matplotlib ValueError / TypeError traceback), are reported up front with the
// problems of core::validate, each naming its field. Qt-free.
//
//   RN-11  output.charts.formats: every format one of charts::chartFormats()
//   RN-12  colours (Matplotlib's to_rgba grammar, render/color.hpp): groups[k].color,
//          visual.base_color, visual.proximity_color, visual.glow_color, visual.ion_color,
//          density.color (when a map is set), output.charts.color (when set: "" is the palette's)
//   RN-13  visual.chain_palette: a bundled colormap (render/colormap.hpp, owner decision D6)
//   RN-15  panel options the panels cannot read (render/panels panel.cpp optionNumber / optionInt):
//          timeseries smooth (a whole number), target (a number), ylim (two numbers); heatmap vmin,
//          vmax (numbers, RN-04), max_ticks; contact_map mode (instant, frame or frequency), frame,
//          max_ticks (whole numbers of at least 1 for max_ticks). Options of other types, and the
//          sequence / secondary panels' own options (their "No data" messages, RN-40, RN-41), are
//          left to the panels.
//
// Messages: "<field>: <problem>", e.g. "visual.glow_color: Invalid RGBA argument: 'gold-ish'",
// "layout.panels[2].options.ylim: expected two numbers, got 5" (k counts from 0, as --set's
// layout.panels.2.options.ylim).

#include <optional>
#include <string>
#include <vector>

#include "mdbridge/core/project.hpp"

namespace mdbridge::render::panels {

// Every problem above, in this order: chart formats, colours, chain palette, panel options.
std::vector<std::string> renderProblems(const core::Project& project);

// The error for a preview (or other single figure) file name: its format, from the extension
// (lower case, as Matplotlib's savefig reads it; no extension means png), must be one of
// charts::chartFormats(): "preview format 'gif' is not supported (choose from png, svg, ...)".
// nullopt when it can be written.
std::optional<std::string> figureFileProblem(const std::string& path, const std::string& what = "preview");
// The format figureFileProblem reads from the name ("png" when it has no extension).
std::string figureFileFormat(const std::string& path);

}  // namespace mdbridge::render::panels
