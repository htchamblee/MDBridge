// Heavy-atom and residue contacts between two groups (mdatk/analysis/contacts.py).
// Descriptor ported; the computation lands with chunk W2.1.

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"

namespace mdbridge::analysis {

namespace {

class ContactsDescriptor final : public AnalysisDescriptor {
public:
    ContactsDescriptor()
        : AnalysisDescriptor({"contacts",
                              "Pairwise contacts",
                              1,
                              false,
                              2,
                              {{"atom_pairs", "{a}–{b} heavy-atom contacts"},
                               {"residue_pairs", "{a}–{b} residue contacts"},
                               {"q_native", "{a}–{b} fraction of native contacts (Q)"},
                               {"per_residue_a", "{a} residue contacts with {b}"},
                               {"per_residue_b", "{b} residue contacts with {a}"},
                               {"map", "{a}–{b} contact map"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.contacts.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        return detail::sectionParams(project, "contacts");
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"atom_pairs", "residue_pairs", "q_native"}, {"per_residue_a", "per_residue_b"}, {"map"}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::ContactParams& c = project.analyses.contacts;
        return {c.group_a, c.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { detail::notImplemented(*this); }
};

}  // namespace

const AnalysisDescriptor& contactsAnalysis() {
    static const ContactsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
