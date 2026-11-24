#include "mdbridge/analysis/registry.hpp"

#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

const std::vector<const AnalysisDescriptor*>& registeredAnalyses() {
    // MDATK's import order (mdatk/analysis/__init__.py): contacts, energy, hbonds, mg, rmsd,
    // stacking, structure (trace, proximity, atoms), then basepairs
    static const std::vector<const AnalysisDescriptor*> registry{
        &contactsAnalysis(), &energyAnalysis(), &hbondsAnalysis(),    &mgAnalysis(),    &rmsdAnalysis(),
        &stackingAnalysis(), &traceAnalysis(),  &proximityAnalysis(), &atomsAnalysis(), &basepairsAnalysis()};
    return registry;
}

const AnalysisDescriptor* findAnalysis(std::string_view name) {
    for (const AnalysisDescriptor* d : registeredAnalyses()) {
        if (d->name() == name) {
            return d;
        }
    }
    return nullptr;
}

int registryIndex(std::string_view name) {
    const auto& registry = registeredAnalyses();
    for (std::size_t i = 0; i < registry.size(); ++i) {
        if (registry[i]->name() == name) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

}  // namespace mdbridge::analysis
