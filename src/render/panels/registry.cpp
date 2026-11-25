// The panel registry (panel.hpp): mdatk/render/panels.py PANEL_TYPES / panel_class.
//
// Each panel family adds its table below, one line per family (wave 2: W2-A plots, W2-B structure,
// W2-C sequence and secondary); the registry keeps MDATK's order (core::panelTypes()).

#include <map>
#include <string>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/render/panels/panel.hpp"
#include "mdbridge/render/panels/seqpanels.hpp"
#include "mdbridge/render/panels/structure.hpp"

namespace mdbridge::render::panels {

namespace {

const std::map<std::string, PanelFactory, std::less<>>& factories() {
    static const std::map<std::string, PanelFactory, std::less<>> table = [] {
        std::map<std::string, PanelFactory, std::less<>> m;
        auto add = [&m](std::span<const PanelRegistration> family) {
            for (const PanelRegistration& r : family) m.emplace(std::string(r.type), r.make);
        };

        add(structurePanels());  // W2-B: structure
        add(plotPanels());  // W2-A: timeseries, dual_axis, heatmap, contact_map
        add(residuePanels());  // W2-C: sequence, secondary

        return m;
    }();
    return table;
}

}  // namespace

std::vector<std::string_view> registeredPanelTypes() {
    std::vector<std::string_view> out;
    for (const core::PanelType& t : core::panelTypes()) {
        if (factories().count(t.name)) out.push_back(t.name);
    }
    return out;
}

std::unique_ptr<Panel> createPanel(Figure& fig, SubplotSpec slot, const core::PanelSpec& spec, Context& ctx) {
    const auto it = factories().find(spec.type);
    if (it == factories().end()) {
        std::string names;
        for (const core::PanelType& t : core::panelTypes()) names += (names.empty() ? "" : ", ") + std::string(t.name);
        throw core::ValueError("unknown panel type " + core::pythonRepr(core::Value(spec.type)) + " (choose from " + names +
                               ")");
    }
    return it->second(fig, slot, spec, ctx);
}

}  // namespace mdbridge::render::panels
