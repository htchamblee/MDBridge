// Base stacking from ring geometry (mdatk/analysis/stacking.py).
// Descriptor ported; the computation lands with chunk W2.2.

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

namespace {

class StackingDescriptor final : public AnalysisDescriptor {
public:
    StackingDescriptor()
        : AnalysisDescriptor({"stacking",
                              "π–π stacking",
                              1,
                              false,
                              2,
                              {{"count", "π–π stacked pairs"},
                               {"inter", "{a}–{b} stacked pairs"},
                               {"intra", "Intra-group stacked pairs"},
                               {"map", "Stacking map"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.stacking.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        return detail::sectionParams(project, "stacking");
    }

    ExpectedOutputs expectedOutputs(const core::Project& project) const override {
        ExpectedOutputs out{{"count", "inter"}, {}, {"map"}};
        if (project.analyses.stacking.scope != "inter") {
            out.series.emplace_back("intra");
        }
        return out;
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::StackingParams& s = project.analyses.stacking;
        return {s.selection, s.group_a, s.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& stackingAnalysis() {
    static const StackingDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
