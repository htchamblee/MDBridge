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
#include "mdbridge/analysis/context.hpp"
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

// energy.COULOMB (kcal mol^-1 Å e^-2)
inline constexpr double kEnergyCoulomb = 332.0522;

// energy._bonded_exclusions(u, a, b, depth): the sorted keys a_local * b.size() + b_local of the
// A-B pairs joined by a path of <= depth bonds through any atoms (AN-08; MDATK only follows paths
// that start at a direct A-B bond), along `bonds` (the topology's when null). None for depth <= 0
// or without bonds. A repeated A atom keys every position; a repeated B atom its last one.
std::vector<std::int64_t> energyBondedExclusions(const core::AtomGroup& a, const core::AtomGroup& b,
                                                 std::int64_t depth, const std::vector<core::Bond>* bonds = nullptr);

// The energy analysis' Lennard-Jones parameters and partial charges without parm7 parameters
// (src/analysis/energy_parameters.cpp; docs/DIFFERENCES.md AN-25 - AN-27, CIT-04), per atom of
// `atoms`, decided on `frame` (the first analysed frame) where bonds must be guessed:
// * Template: the residue matches a residue template of Amber's ff14SB, OL15, OL3 or TIP3P water /
//   ion libraries (OpenMM 8.1.2's conversion): R* / epsilon of the template atom; its charge when
//   the residue has every atom of the template (templateCharge), else none.
// * Typed: one Amber atom type by element (parm10.dat; Joung-Cheatham / Li-Merz TIP3P ions for
//   monatomic ions and metals), hydrogens by the element they are bonded to, sp3 carbon and
//   hydroxyl / ether / water oxygen by their bonded neighbours; elements without a type take
//   carbon's (defaultElements, defaultAtoms).
struct EnergyParameters {
    enum class Layer : std::uint8_t { Template, Typed };
    std::vector<double> rstar;    // R* = Rmin / 2, Å
    std::vector<double> epsilon;  // kcal/mol
    std::vector<double> charge;   // e: the template's where templateCharge, else 0
    std::vector<Layer> layer;
    std::vector<bool> templateCharge;
    std::vector<std::string> forceFields;      // template force fields used (ff14SB, OL15, OL3, TIP3P), first use
    std::vector<std::string> references;       // citation keys of what was used, first use
    std::vector<std::string> defaultElements;  // elements that took carbon's parameters, first use
    std::size_t defaultAtoms = 0;
    std::size_t chargelessResidues = 0;  // residues of `atoms` with atoms without a template charge
    bool guessedBonds = false;  // bonds were guessed (AN-14, Cordero et al.'s radii) to decide a type or a chain end
};
EnergyParameters energyParameters(const core::AtomGroup& atoms, const core::Frame& frame);

// "Cornell et al. 1995, doi:10.1021/ja00124a002" for a key of EnergyParameters::references.
std::string energyShortCitation(std::string_view key);

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

// AN-14 (docs/DIFFERENCES.md): bonds among `atoms` guessed from their positions in `frame`, as
// local index pairs (i < j) [m, 2] row-major, sorted: atoms no farther apart than the sum of their
// covalent radii (Cordero et al. 2008; carbon's for an unknown element) + 0.4 Å, never two
// hydrogens, never an atom of a one-atom residue (ions), and within one residue, between
// neighbouring residues (residue indices differing by 1, whatever their chains or segments: O3'-P
// across a chain break) or between two sulfurs of any residues (a disulfide); MDATK 0.1.2's
// residue rule. Elements from core::atomElements. (MDATK 0.0.4.1: 1.95 Å, H within 1.25 Å,
// neighbouring residues only; 0.1.2: covalent radii + 0.4 Å but at least 1.95 Å, H within 1.25 Å.)
std::vector<std::int64_t> guessBonds(const core::AtomGroup& atoms, const core::Frame& frame);

// structure.guess_element(name): the element the atoms analysis tells hydrogens by: the name
// stripped of surrounding whitespace and leading digits; "Mg", "Na", "Cl", "Zn", "Ca" for exactly
// those two letters (any case), else its first character upper-cased ("X" for an empty name).
// Upper-casing is ASCII only (Python's str.upper also maps other letters).
std::string structureGuessElement(std::string_view name);

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
