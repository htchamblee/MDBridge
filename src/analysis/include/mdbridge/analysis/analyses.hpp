#pragma once
// The built-in analyses (mdatk/analysis/*.py), in registry order, and the class-level helpers of
// their modules that other code needs.

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/core/project.hpp"

namespace mdbridge::analysis {

const AnalysisDescriptor& contactsAnalysis();   // contacts.cpp
const AnalysisDescriptor& energyAnalysis();     // energy.cpp
const AnalysisDescriptor& hbondsAnalysis();     // hbonds.cpp
const AnalysisDescriptor& mgAnalysis();         // mg.cpp
const AnalysisDescriptor& rmsdAnalysis();       // rmsd.cpp
const AnalysisDescriptor& stackingAnalysis();   // stacking.cpp
const AnalysisDescriptor& traceAnalysis();      // structure.cpp
const AnalysisDescriptor& proximityAnalysis();  // structure.cpp
const AnalysisDescriptor& atomsAnalysis();      // structure.cpp
const AnalysisDescriptor& basepairsAnalysis();  // basepairs.cpp

// rmsd._out_key(label, k): runs of characters other than [A-Za-z0-9_.+-] become "_", leading and
// trailing "_" are stripped; "target<k>" when nothing is left.
std::string rmsdOutputKey(std::string_view label, std::size_t k);

// energy.COULOMB (kcal mol^-1 Å e^-2)
inline constexpr double kEnergyCoulomb = 332.0522;

// energy._bonded_exclusions(u, a, b, depth): the sorted keys a_local * b.size() + b_local of the
// A-B pairs separated by <= depth bonds (none for depth <= 0 or a topology without bonds). The
// search starts from the A atoms within depth - 1 bonds (through A atoms) of a bond between A and
// B, and reaches B atoms through any atoms.
std::vector<std::int64_t> energyBondedExclusions(const core::AtomGroup& a, const core::AtomGroup& b,
                                                 std::int64_t depth);

// structure.structure_groups(project): visual.structure_groups (all group names when empty),
// unknown names left out.
std::vector<std::string> structureGroups(const core::Project& project);

}  // namespace mdbridge::analysis
