#pragma once
// The built-in analyses (mdatk/analysis/*.py), in registry order, and the class-level helpers of
// their modules that other code needs.

#include <cstddef>
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

// structure.structure_groups(project): visual.structure_groups (all group names when empty),
// unknown names left out.
std::vector<std::string> structureGroups(const core::Project& project);

}  // namespace mdbridge::analysis
