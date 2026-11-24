#pragma once
// The built-in analyses (mdatk/analysis/*.py), in registry order, and the class-level helpers of
// their modules that other code needs.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/rigid_transform.hpp"  // Vec3
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

// stacking.ring_atoms(residue): the atoms of the residue's aromatic ring in ring order (purine
// N9 ... C4, else pyrimidine N1 ... C6, else the side-chain ring of PHE / TYR / TRP / HIS and its
// protonation names, by the upper-cased residue name), or std::nullopt; a repeated atom name
// refers to the residue's last atom of that name.
std::optional<std::vector<std::int32_t>> stackingRingAtoms(const core::Topology& topology, std::int32_t residue);

// stacking.ring_geometry(positions): per ring of `ringSize` consecutive points the centroid (mean
// over the ring atoms), the unit normal of the least-squares plane (the last right singular
// vector, up to sign) and the planarity s_min / sqrt(ringSize). Throws std::invalid_argument
// unless positions.size() is a multiple of ringSize > 0.
struct RingGeometry {
    std::vector<core::Vec3> centroids;
    std::vector<core::Vec3> normals;
    std::vector<double> planarity;
};
RingGeometry ringGeometry(std::span<const core::Vec3> positions, std::size_t ringSize);

}  // namespace mdbridge::analysis
