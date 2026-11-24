// Non-bonded interaction energy between two groups (mdatk/analysis/energy.py).
// Descriptor ported; the computation lands with chunk W2.3.

#include <algorithm>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

namespace {

class EnergyDescriptor final : public AnalysisDescriptor {
public:
    EnergyDescriptor()
        : AnalysisDescriptor({"energy",
                              "Intermolecular energy",
                              1,
                              false,
                              2,
                              {{"total", "{a}–{b} interaction energy"},
                               {"elec", "{a}–{b} electrostatic energy"},
                               {"vdw", "{a}–{b} van der Waals energy"},
                               {"net_force", "Net force of {b} on {a}"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.energy.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        core::Table params = detail::sectionParams(project, "energy");
        // target_value is display only
        std::erase_if(params, [](const auto& entry) { return entry.first == "target_value"; });
        return params;
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"total", "elec", "vdw", "net_force"}, {}, {}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::EnergyParams& e = project.analyses.energy;
        return {e.group_a, e.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& energyAnalysis() {
    static const EnergyDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
