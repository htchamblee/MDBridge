// Mg2+ coordination, bridges and occupancy (mdatk/analysis/mg.py).
// Descriptor ported; the computation lands with chunk W2.2.

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

namespace {

class MgDescriptor final : public AnalysisDescriptor {
public:
    MgDescriptor()
        : AnalysisDescriptor({"mg",
                              "Mg coordination",
                              1,
                              false,
                              2,
                              {{"contacts", "Mg–O inner-sphere contacts"},
                               {"outer_contacts", "Mg–O outer-sphere contacts"},
                               {"bridges", "O–Mg–O bridges"},
                               {"inter_bridges", "{a}–Mg–{b} bridges"},
                               {"bridging_ions", "Bridging Mg ions ({a}–{b})"},
                               {"inner_residues", "Residues with inner-sphere Mg"},
                               {"occupancy", "Mg occupancy"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.mg.enabled; }

    core::Table paramsFor(const core::Project& project) const override { return detail::sectionParams(project, "mg"); }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"contacts", "outer_contacts", "bridges", "inter_bridges", "bridging_ions", "inner_residues"},
                {"occupancy"},
                {}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::MgParams& m = project.analyses.mg;
        return {m.mg_selection, m.ligand_selection, m.occupancy_selection, m.group_a, m.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& mgAnalysis() {
    static const MgDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
