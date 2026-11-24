// Root mean square deviation against a reference frame (mdatk/analysis/rmsd.py).
// Descriptor ported; the computation lands with chunk W2.4.

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

std::string rmsdOutputKey(std::string_view label, std::size_t k) {
    // re.sub(r"[^A-Za-z0-9_.+\-]+", "_", label).strip("_")
    const auto keep = [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
               c == '+' || c == '-';
    };
    std::string key;
    bool inRun = false;
    for (const char c : label) {
        if (keep(c)) {
            key += c;
            inRun = false;
        } else if (!inRun) {
            key += '_';
            inRun = true;
        }
    }
    const std::size_t first = key.find_first_not_of('_');
    if (first == std::string::npos) {
        return "target" + std::to_string(k);
    }
    return key.substr(first, key.find_last_not_of('_') - first + 1);
}

namespace {

class RmsdDescriptor final : public AnalysisDescriptor {
public:
    RmsdDescriptor() : AnalysisDescriptor({"rmsd", "RMSD", 1, false, 1, {}}) {}

    void withGroups(core::Project& project, const std::vector<std::string>& groups) const override {
        const std::string& g = groups.at(0);
        project.analyses.rmsd.targets = {{"Heavy", "group:" + g + " and not name H*", ""},
                                         {"Backbone", "group:" + g + " and name C4'", ""}};
    }

    std::string outputLabel(const core::Project& project, std::string_view out) const override {
        const auto& targets = project.analyses.rmsd.targets;
        for (std::size_t k = 0; k < targets.size(); ++k) {
            if (rmsdOutputKey(targets[k].label, k) == out) {
                return "RMSD " + targets[k].label;
            }
        }
        return "RMSD " + std::string(out);
    }

    bool enabled(const core::Project& project) const override {
        return project.analyses.rmsd.enabled && !project.analyses.rmsd.targets.empty();
    }

    core::Table paramsFor(const core::Project& project) const override {
        const core::Table section = detail::sectionParams(project, "rmsd");
        core::Table out;
        for (const char* key : {"targets", "reference_frame", "superposition", "mass_weighted"}) {
            for (const auto& [k, v] : section) {
                if (k == key) {
                    out.emplace_back(k, v);
                }
            }
        }
        return out;
    }

    ExpectedOutputs expectedOutputs(const core::Project& project) const override {
        ExpectedOutputs out;
        const auto& targets = project.analyses.rmsd.targets;
        for (std::size_t k = 0; k < targets.size(); ++k) {
            out.series.push_back(rmsdOutputKey(targets[k].label, k));
        }
        return out;
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out;
        for (const core::RMSDTarget& t : project.analyses.rmsd.targets) {
            out.push_back(t.selection);
            out.push_back(t.fit_selection);
        }
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& rmsdAnalysis() {
    static const RmsdDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
