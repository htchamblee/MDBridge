// Base pairs of the structural-view residues (mdatk/analysis/basepairs.py).
// Descriptor ported; the computation lands with chunk W2.6.

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

namespace {

class BasePairsDescriptor final : public AnalysisDescriptor {
public:
    BasePairsDescriptor()
        : AnalysisDescriptor(
              {"basepairs", "Base pairs", 1, true, 0, {{"map", "Base pairs"}, {"count", "Base pairs"}}}) {}

    core::Table paramsFor(const core::Project& project) const override {
        const core::BasePairParams& bp = project.analyses.basepairs;
        const std::vector<std::string> names = structureGroups(project);
        return {{"groups", core::Array(names.begin(), names.end())},
                {"ions", project.analyses.mg.mg_selection},
                {"engine", bp.engine},
                {"stride", bp.stride},
                {"cutoff", bp.cutoff},
                {"wobble", bp.wobble},
                {"noncanonical", bp.noncanonical}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override { return {{"count"}, {}, {"map"}}; }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& basepairsAnalysis() {
    static const BasePairsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
