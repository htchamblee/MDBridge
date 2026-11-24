// Internal analyses feeding the structural view (mdatk/analysis/structure.py): backbone trace,
// proximity to a residue set, atomic coordinates.
// Descriptors ported; the computations land with chunk W2.5.

#include <set>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

std::vector<std::string> structureGroups(const core::Project& project) {
    const std::vector<std::string> all = core::groupNames(project);
    const std::vector<std::string>& shown = project.visual.structure_groups;
    const std::vector<std::string>& names = shown.empty() ? all : shown;
    const std::set<std::string> known(all.begin(), all.end());
    std::vector<std::string> out;
    for (const std::string& n : names) {
        if (known.contains(n)) {
            out.push_back(n);
        }
    }
    return out;
}

namespace {

core::Array stringArray(const std::vector<std::string>& items) {
    return core::Array(items.begin(), items.end());
}

class TraceDescriptor final : public AnalysisDescriptor {
public:
    // version 4: v2 residue numbers for tags; v3 group membership; v4 chain and residue names
    TraceDescriptor() : AnalysisDescriptor({"trace", "Backbone trace", 4, true, 0, {}}) {}

    core::Table paramsFor(const core::Project& project) const override {
        const core::VisualParams& v = project.visual;
        return {{"groups", stringArray(structureGroups(project))},
                {"trace_atoms", v.trace_atoms},
                {"align", v.align_frames},
                {"show_ions", v.show_ions},
                {"ions", v.show_ions ? project.analyses.mg.mg_selection : std::string()},
                {"group_order", stringArray(core::groupNames(project))}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

class ProximityDescriptor final : public AnalysisDescriptor {
public:
    ProximityDescriptor() : AnalysisDescriptor({"proximity", "Proximity to residue set", 1, true, 0, {}}) {}

    ExpectedOutputs expectedOutputs(const core::Project&) const override { return {{}, {"min_distance"}, {}}; }

    core::Table paramsFor(const core::Project& project) const override {
        return {{"groups", stringArray(structureGroups(project))},
                {"set", project.visual.proximity_set},
                {"radius", project.analyses.proximity_search_radius},
                {"ions", project.analyses.mg.mg_selection}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.visual.proximity_set);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

class AtomsDescriptor final : public AnalysisDescriptor {
public:
    // version 2: elements guessed per residue (protein CA is carbon, not calcium)
    AtomsDescriptor() : AnalysisDescriptor({"atoms", "Atomic coordinates (structure view)", 2, true, 0, {}}) {}

    core::Table paramsFor(const core::Project& project) const override {
        const core::VisualParams& v = project.visual;
        return {{"groups", stringArray(structureGroups(project))},
                {"trace_atoms", v.trace_atoms},
                {"align", v.align_frames},
                {"hydrogens", v.rep_hydrogens},
                {"ions", project.analyses.mg.mg_selection},
                {"group_order", stringArray(core::groupNames(project))}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& traceAnalysis() {
    static const TraceDescriptor descriptor;
    return descriptor;
}

const AnalysisDescriptor& proximityAnalysis() {
    static const ProximityDescriptor descriptor;
    return descriptor;
}

const AnalysisDescriptor& atomsAnalysis() {
    static const AtomsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
