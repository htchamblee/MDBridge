#pragma once
// The registered analyses and the data keys built from them (mdatk/analysis/base.py and
// __init__.py): instance ids ("contacts@Aptamer+Mg" is the contacts analysis run on the groups
// Aptamer and Mg), data keys ("<instance>/<output>"), labels, expected keys and the residue groups
// a cache key covers.
//
// Errors carry MDATK's messages: an unknown analysis is std::out_of_range (MDATK's KeyError,
// without the quotes Python adds when printing it), a bad instance core::ValueError.

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::analysis {

inline constexpr char kInstanceSeparator = '@';
inline constexpr char kGroupSeparator = '+';

// Every registered analysis in MDATK's order: contacts, energy, hbonds, mg, rmsd, stacking,
// trace, proximity, atoms, basepairs.
const std::vector<const AnalysisDescriptor*>& registeredAnalyses();
// The analysis named `name` (no instance part), or nullptr.
const AnalysisDescriptor* findAnalysis(std::string_view name);
// Position of `name` in the registry, or -1.
int registryIndex(std::string_view name);

struct Instance {
    std::string base;
    std::vector<std::string> groups;
    bool operator==(const Instance&) const = default;
};

// parse_instance: "contacts@A+B" -> {"contacts", {"A", "B"}}; "rmsd" -> {"rmsd", {}}; empty group
// names are dropped ("contacts@A++B" -> {"A", "B"}, "contacts@" -> {}).
Instance parseInstance(std::string_view instance);
// instance_id(base, groups): empty names dropped; the base alone without groups.
std::string instanceId(std::string_view base, const std::vector<std::string>& groups = {});
// split_key: "contacts@A+B/atom_pairs" -> {"contacts@A+B", "atom_pairs"} (at the first '/'; the
// output is "" without one).
std::pair<std::string, std::string> splitKey(std::string_view key);

// analysis_class(instance): throws std::out_of_range("unknown analysis 'x' (known: contacts,
// energy, ...)").
const AnalysisDescriptor& analysisFor(std::string_view instance);

// instance_project(project, instance): the project an instance runs with, a copy with its groups
// substituted (the project unchanged without groups). Throws std::out_of_range (unknown
// analysis) or core::ValueError: "<name> does not take residue groups ('<instance>')",
// "'<instance>': <name> takes N group(s)", "'<instance>': unknown residue group(s) A, B".
core::Project instanceProject(const core::Project& project, std::string_view instance);

// key_label(key, project): the label of a data key without computing anything; the key itself
// when the analysis or the instance is invalid.
std::string keyLabel(std::string_view key, const core::Project& project);

// expected_keys(project, include_disabled): "<name>/<output>" per kind, in registry order;
// without includeDisabled, user analyses that are not enabled are left out (internal ones never).
ExpectedOutputs expectedKeys(const core::Project& project, bool includeDisabled = true);

// referenced_groups(texts, project): the definitions {"name", "residues", "selection"} of every
// residue group named by group:<Name> tokens or as a whole (stripped) text, followed through the
// groups' own selections, sorted by name (a repeated name: its last definition).
core::Array referencedGroups(const std::vector<std::string>& texts, const core::Project& project);

}  // namespace mdbridge::analysis
