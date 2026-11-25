// The panel types and layout presets of MDATK's figure (data: panel_tables.inc).
#include "mdbridge/core/panels.hpp"

#include <array>
#include <cstddef>
#include <iterator>

namespace mdbridge::core {

namespace {

struct LayoutPresetRow {
    std::string_view name;
    std::string_view label;
    int rows;
    int cols;
    std::size_t firstSlot;
    std::size_t slotCount;
    std::size_t firstRatio;
    std::size_t ratioCount;
};

#include "panel_tables.inc"

constexpr std::size_t kPresetCount = std::size(kLayoutPresets);

const std::array<LayoutPreset, kPresetCount>& presets() {
    static const std::array<LayoutPreset, kPresetCount> table = [] {
        std::array<LayoutPreset, kPresetCount> out{};
        for (std::size_t i = 0; i < kPresetCount; ++i) {
            const LayoutPresetRow& r = kLayoutPresets[i];
            out[i] = LayoutPreset{r.name,
                                  r.label,
                                  r.rows,
                                  r.cols,
                                  std::span<const LayoutSlot>(kLayoutSlots).subspan(r.firstSlot, r.slotCount),
                                  std::span<const double>(kLayoutHeightRatios).subspan(r.firstRatio, r.ratioCount)};
        }
        return out;
    }();
    return table;
}

}  // namespace

std::span<const PanelType> panelTypes() { return kPanelTypes; }

const PanelType* findPanelType(std::string_view name) {
    for (const PanelType& type : kPanelTypes) {
        if (type.name == name) {
            return &type;
        }
    }
    return nullptr;
}

std::span<const LayoutPreset> layoutPresets() { return presets(); }

const LayoutPreset* findLayoutPreset(std::string_view name) {
    for (const LayoutPreset& preset : presets()) {
        if (preset.name == name) {
            return &preset;
        }
    }
    return nullptr;
}

}  // namespace mdbridge::core
