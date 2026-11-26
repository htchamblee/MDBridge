#pragma once
// The panel types and layout presets of MDATK's figure (mdatk/render/panels.py PANEL_TYPES,
// mdatk/render/layout.py PRESETS, PRESET_LABELS, HEIGHT_RATIOS), generated from MDATK by
// tools/gen_panel_tables.py into src/core/panel_tables.inc: what `list` prints, the panel types
// the pipeline knows, and the grid the renderer lays panels out on.

#include <span>
#include <string_view>

namespace mdbridge::core {

struct PanelType {
    std::string_view name;         // type_name ("structure", "timeseries", ...)
    std::string_view label;        // e.g. "Time series"
    std::string_view description;  // one paragraph
    std::string_view dataKind;     // what spec.series refers to: "series", "matrices" or "pairmaps"
    int minSeries;
    int maxSeries;  // 0: the panel takes no series
};

// A range of grid cells, ends exclusive: rows row0 .. row1 - 1, columns col0 .. col1 - 1.
struct LayoutSlot {
    int row0;
    int row1;
    int col0;
    int col1;
};

struct LayoutPreset {
    std::string_view name;   // e.g. "large_plus_3"
    std::string_view label;  // e.g. "1 large + 3 small (right)"
    int rows;
    int cols;
    std::span<const LayoutSlot> slots;     // filled in order by the layout's panels
    std::span<const double> heightRatios;  // relative row heights; empty: equal rows
};

// The slots of a preset (Qt code, where `slots` is a macro, includes this header with the macro
// pushed away and reads the slots through this).
inline std::span<const LayoutSlot> presetSlots(const LayoutPreset& preset) { return preset.slots; }

// PANEL_TYPES in registration order: structure, timeseries, dual_axis, heatmap, contact_map,
// sequence, secondary.
std::span<const PanelType> panelTypes();
const PanelType* findPanelType(std::string_view name);

// PRESETS in MDATK's order (single, large_plus_2, ...).
std::span<const LayoutPreset> layoutPresets();
const LayoutPreset* findLayoutPreset(std::string_view name);

}  // namespace mdbridge::core
