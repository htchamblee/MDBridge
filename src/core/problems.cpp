// Located problems: the presentation of mdatk/problems.py of MDATK 0.1.1 (owner's code, ported).

#include "mdbridge/core/problems.hpp"

#include <array>
#include <utility>

#include "mdbridge/core/project.hpp"
#include "mdbridge/core/value.hpp"
#include "project_internal.hpp"

namespace mdbridge::core {

namespace {

using Pair = std::pair<std::string_view, std::string_view>;

// SECTION_TABS: first component of a field path -> the page that shows it
constexpr std::array<Pair, 9> kSectionPages{{
    {"notes", "input"}, {"inputs", "input"}, {"cleaning", "input"}, {"groups", "groups"},
    {"analyses", "analyses"}, {"visual", "structure"}, {"density", "structure"}, {"layout", "layout"},
    {"output", "outputs"},
}};

// _LABELS: readable names of field paths ('#' stands for a list index)
constexpr std::array<Pair, 43> kLabels{{
    {"notes", "Project notes"},
    {"inputs.topology", "Input files › Topology"},
    {"inputs.trajectories", "Input files › Trajectories"},
    {"inputs.start", "Frames › Start"},
    {"inputs.stop", "Frames › Stop"},
    {"inputs.step", "Frames › Step"},
    {"inputs.dt_ps", "Frames › Time step override"},
    {"inputs.external", "External data"},
    {"inputs.external.#", "External file {n}"},
    {"cleaning", "Trajectory cleaning"},
    {"groups.#", "Group {n}{name}"},
    {"groups.#.selection", "Selection (overrides)"},
    {"groups.#.color", "Colour"},
    {"analyses.rmsd", "RMSD"},
    {"analyses.rmsd.targets.#", "Target {n}{name}"},
    {"analyses.contacts", "Pairwise contacts"},
    {"analyses.hbonds", "Hydrogen bonds"},
    {"analyses.mg", "Mg²⁺ coordination and occupancy"},
    {"analyses.mg.occupancy_selection", "Occupancy residues"},
    {"analyses.stacking", "π–π stacking"},
    {"analyses.basepairs", "Base pairs"},
    {"analyses.energy", "Intermolecular energy"},
    {"visual.representations.#", "Representation {n}"},
    {"visual.structure_groups", "Appearance › Groups shown"},
    {"visual.proximity_set", "Appearance › Proximity set"},
    {"visual.color_mode", "Appearance › Colouring"},
    {"visual.glow_mode", "Appearance › Glow"},
    {"visual.center_on", "Orientation › Centre on"},
    {"visual.zoom", "Orientation › Zoom"},
    {"visual.label_residues", "Residue number tags › Label residues"},
    {"visual.rep_rest", "Representations › Other residues"},
    {"visual.backbone_atom", "Representations › Backbone atom"},
    {"visual.renderer", "Representations › Renderer"},
    {"density", "Density map"},
    {"density.map_path", "Map file"},
    {"layout.preset", "Grid preset"},
    {"layout.panels.#", "Panel {n}{name}"},
    {"layout.panels.#.series", "Data"},
    {"output.video", "Video"},
    {"output.charts", "Time charts"},
    {"output.pdb", "Frame PDBs"},
    {"output.pdb.selection", "Atom subset"},
    {"output.video.fps", "FPS"},
}};

// _LEAF_LABELS: names of the last path component wherever it occurs
constexpr std::array<Pair, 12> kLeafLabels{{
    {"group_a", "Group A"}, {"group_b", "Group B"}, {"groups", "Groups"}, {"selection", "Selection"},
    {"fit_selection", "Fit selection"}, {"mg_selection", "Mg selection"},
    {"ligand_selection", "Ligand selection"}, {"residues", "Residues"}, {"name", "Name"},
    {"options", "Options"}, {"style", "Style"}, {"color", "Colouring"},
}};

// PANEL_NAMES: short names of the panel types for "Panel 3 (Contact map)"
constexpr std::array<Pair, 7> kPanelNames{{
    {"structure", "Structure"}, {"timeseries", "Time series"}, {"dual_axis", "Dual-axis time series"},
    {"heatmap", "Residue heatmap"}, {"contact_map", "Contact map"}, {"sequence", "Sequence"},
    {"secondary", "Secondary structure"},
}};

constexpr std::string_view kSeparator = " › ";

template <std::size_t N>
std::optional<std::string_view> lookup(const std::array<Pair, N>& table, std::string_view key) {
    for (const auto& [name, text] : table) {
        if (name == key) {
            return text;
        }
    }
    return std::nullopt;
}

bool isIndex(std::string_view part) {
    if (part.empty()) {
        return false;
    }
    for (const char c : part) {
        if (c < '0' || c > '9') {
            return false;
        }
    }
    return true;
}

std::vector<std::string_view> split(std::string_view path) {
    std::vector<std::string_view> parts;
    std::size_t start = 0;
    for (;;) {
        const std::size_t end = path.find('.', start);
        parts.push_back(path.substr(start, end == std::string_view::npos ? std::string_view::npos : end - start));
        if (end == std::string_view::npos) {
            return parts;
        }
        start = end + 1;
    }
}

// str.capitalize() of an ASCII key
std::string capitalize(std::string_view text) {
    std::string out(text);
    for (std::size_t i = 0; i < out.size(); ++i) {
        char& c = out[i];
        if (i == 0 && c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        } else if (i > 0 && c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return out;
}

// _item_name: " (Aptamer)", what the list item parts[0..k] is called
std::string itemName(const std::vector<std::string_view>& parts, std::size_t k, const Value* plain) {
    if (!plain) {
        return {};
    }
    const Value* node = plain;
    for (std::size_t i = 0; i <= k && node; ++i) {
        if (const Array* array = node->get<Array>()) {
            if (!isIndex(parts[i])) {
                return {};
            }
            const std::size_t index = std::stoul(std::string(parts[i]));
            node = index < array->size() ? &(*array)[index] : nullptr;
        } else if (node->isTable() && !isIndex(parts[i])) {
            node = node->find(parts[i]);
        } else {
            node = nullptr;
        }
    }
    if (!node || !node->isTable()) {
        return {};
    }
    std::string name;
    if (parts[0] == "layout") {  // a panel is known by its type
        const Value* type = node->find("type");
        const std::string kind = type && type->isString() ? type->stringOr() : std::string();
        name = std::string(lookup(kPanelNames, kind).value_or(kind));
    } else {
        const Value* own = node->find("name");
        name = own && own->isString() ? own->stringOr() : std::string();
        if (name.empty()) {
            const Value* label = node->find("label");
            name = label && label->isString() ? label->stringOr() : std::string();
        }
    }
    return name.empty() ? std::string() : " (" + name + ")";
}

std::string format(std::string_view label, std::size_t n, const std::string& name) {
    std::string out;
    for (std::size_t i = 0; i < label.size(); ++i) {
        if (label.compare(i, 3, "{n}") == 0) {
            out += std::to_string(n);
            i += 2;
        } else if (label.compare(i, 6, "{name}") == 0) {
            out += name;
            i += 5;
        } else {
            out += label[i];
        }
    }
    return out;
}

}  // namespace

namespace detail {

std::string describeField(std::string_view field, const Value* plain) {
    if (field.empty()) {
        return {};
    }
    const std::vector<std::string_view> parts = split(field);
    const std::string_view page = problemPage(field);
    std::string out(page.empty() ? capitalize(parts[0]) : std::string(problemPageTitle(page)));
    std::string pattern;
    for (std::size_t k = 0; k < parts.size(); ++k) {
        if (k) {
            pattern += '.';
        }
        pattern += isIndex(parts[k]) ? std::string("#") : std::string(parts[k]);
        std::optional<std::string> label;
        if (const auto known = lookup(kLabels, pattern)) {
            label = std::string(*known);
        } else if (k == parts.size() - 1 && k > 0 && !isIndex(parts[k])) {
            if (const auto leaf = lookup(kLeafLabels, parts[k])) {
                label = std::string(*leaf);
            } else {
                std::string spaced(parts[k]);
                for (char& c : spaced) {
                    if (c == '_') {
                        c = ' ';
                    }
                }
                label = capitalize(spaced);
            }
        }
        if (label && !label->empty()) {
            const std::size_t n = isIndex(parts[k]) ? std::stoul(std::string(parts[k])) + 1 : 0;
            const std::string name =
                label->find("{name}") != std::string::npos ? itemName(parts, k, plain) : std::string();
            out += kSeparator;
            out += format(*label, n, name);
        }
    }
    return out;
}

}  // namespace detail

std::string_view Problem::page() const { return problemPage(field); }

std::string Problem::text() const { return location.empty() ? message : location + ": " + message; }

std::string_view problemPage(std::string_view field) {
    if (field.empty()) {
        return {};
    }
    return lookup(kSectionPages, field.substr(0, field.find('.'))).value_or(std::string_view());
}

std::string_view problemPageTitle(std::string_view page) {
    for (const ProblemPage& p : kProblemPages) {
        if (p.id == page) {
            return p.title;
        }
    }
    return {};
}

std::string describeField(std::string_view field, const Project* project) {
    if (!project) {
        return detail::describeField(field, nullptr);
    }
    const Value plain = toPlain(*project);
    return detail::describeField(field, &plain);
}

Problem locatedProblem(std::string field, std::string message, const Project* project, std::string level) {
    Problem out;
    out.location = describeField(field, project);
    out.field = std::move(field);
    out.message = std::move(message);
    out.level = std::move(level);
    return out;
}

std::vector<Problem> projectProblems(const Project& project) {
    std::vector<Problem> out = detail::checkProblems(project);
    for (std::string& text : layoutWarnings(project)) {
        Problem warning;
        warning.field = "layout.panels";
        warning.message = std::move(text);
        warning.level = "warning";
        out.push_back(std::move(warning));
    }
    const Value plain = toPlain(project);
    for (Problem& p : out) {
        p.location = detail::describeField(p.field, &plain);
    }
    return out;
}

}  // namespace mdbridge::core
