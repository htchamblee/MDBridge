#pragma once
// The built-in analyses (mdatk/analysis/*.py), in registry order, and the class-level helpers of
// their modules that other code needs.

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"

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

// structure.structure_selection_text(project): "group:A or group:B" of structureGroups(); throws
// core::SelectionError("no residue groups available for the structural view") when there are none.
std::string structureSelectionText(const core::Project& project);

// structure.structure_residues(ctx): the residues of the structural view (ions are drawn
// separately) and their atoms.
struct StructureResidues {
    std::vector<std::int32_t> residues;  // sorted unique residue indices of `atoms`
    core::AtomGroup atoms;  // select(structureSelectionText) minus select(mg_selection, allow_empty)
};
// Throws core::SelectionError from the selections (no structure groups, no atoms selected, invalid
// selections); the residues may be empty (every atom an ion).
StructureResidues structureResidues(const AnalysisContext& context);

// structure.guess_element(name): the element the atoms analysis tells hydrogens by: the name
// stripped of surrounding whitespace and leading digits; "Mg", "Na", "Cl", "Zn", "Ca" for exactly
// those two letters (any case), else its first character upper-cased ("X" for an empty name).
// Upper-casing is ASCII only (Python's str.upper also maps other letters).
std::string structureGuessElement(std::string_view name);

}  // namespace mdbridge::analysis
