// Hydrogen bonds between two groups (mdatk/analysis/hbonds.py).
// Descriptor ported; the computation lands with chunk W2.1.

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

namespace {

class HBondsDescriptor final : public AnalysisDescriptor {
public:
    HBondsDescriptor()
        : AnalysisDescriptor({"hbonds",
                              "Hydrogen bonds",
                              1,
                              false,
                              2,
                              {{"count", "{a}–{b} H-bonds"},
                               {"per_residue_a", "{a} residue H-bonds"},
                               {"per_residue_b", "{b} residue H-bonds"},
                               {"map", "{a}–{b} H-bond map"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.hbonds.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        return detail::sectionParams(project, "hbonds");
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"count"}, {"per_residue_a", "per_residue_b"}, {"map"}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::HBondParams& h = project.analyses.hbonds;
        return {h.group_a, h.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& hbondsAnalysis() {
    static const HBondsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
