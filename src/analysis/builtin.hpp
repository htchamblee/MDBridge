#pragma once
// Private helpers of the built-in analyses' descriptors.

#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::analysis::detail {

// dict(vars(project.analyses.<section>), enabled=None): the section's fields in declaration
// order (toPlain), with "enabled" Null in its place.
core::Table sectionParams(const core::Project& project, std::string_view section);

// The group_a / group_b fields of analyses.<section> (contacts, hbonds, mg, stacking, energy), or
// nulls for sections without them.
struct GroupFields {
    const std::string* a = nullptr;
    const std::string* b = nullptr;
};
GroupFields groupFields(const core::Project& project, std::string_view section);
struct MutableGroupFields {
    std::string* a = nullptr;
    std::string* b = nullptr;
};
MutableGroupFields groupFields(core::Project& project, std::string_view section);

// The computation of an analysis whose port has not landed yet: throws
// core::NamedError("NotImplementedError", "<name>: not implemented yet").
[[noreturn]] void notImplemented(const AnalysisDescriptor& descriptor);

}  // namespace mdbridge::analysis::detail
