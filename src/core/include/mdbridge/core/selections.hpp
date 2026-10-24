#pragma once
// MDATK's selection layer (mdatk/selections.py): residue groups (group:<Name> in selections),
// residue ranges, element guessing, residue labels and the frame / time lists of the PDB export.

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"

namespace mdbridge::core {

class SelectionError : public ValueError {
public:
    using ValueError::ValueError;
    const char* pythonName() const noexcept override { return "SelectionError"; }
};

struct ResidueRange {
    std::optional<std::string> chain;
    std::int64_t first = 0;
    std::int64_t last = 0;
    bool operator==(const ResidueRange&) const = default;
};

// "1-40, 45, A:50-60" -> {(none, 1, 40), (none, 45, 45), (A, 50, 60)}; reversed ranges are
// swapped. Throws SelectionError with MDATK's message.
std::vector<ResidueRange> parseResidueRanges(std::string_view text);

// The ranges as a selection: "(resid 1:40) or (chainID A and resid 50:60)"; chains are chainID
// when the topology has chain IDs, else segid.
std::string rangesToSelection(std::string_view text, const Topology* topology = nullptr);

// "(<selection>)" of a group (its selection wins over its residues).
std::string groupSelection(const ResidueGroup& group, const Topology* topology = nullptr);

// Every group:<Name> replaced by the group's selection, recursively (at most 8 levels).
std::string expandGroups(std::string_view text, const std::vector<ResidueGroup>& groups,
                         const Topology* topology = nullptr);

// A bare group name, "sel:<selection>" or a selection, with the groups expanded.
std::string resolveSelection(std::string_view text, const std::vector<ResidueGroup>& groups,
                             const Topology* topology = nullptr);

// select(universe, text, groups, allow_empty): the atoms, sorted and unique. Throws SelectionError
// ("invalid selection 'x' -> 'y': <MDAnalysis error>", "selection 'x' matched no atoms", ...).
AtomGroup select(const std::shared_ptr<const Topology>& topology, const Frame* current, std::string_view text,
                 const std::vector<ResidueGroup>& groups = {}, bool allowEmpty = false);
AtomGroup select(const System& system, std::string_view text, const std::vector<ResidueGroup>& groups = {},
                 bool allowEmpty = false);

// Upper-case element of an atom name ("C4'" -> "C"; ions of single-atom residues -> "MG", "NA" ...).
std::string guessElement(std::string_view name, std::string_view resname = {}, bool singleAtomResidue = false);

// Upper-case element per atom: the topology's when every atom of the group has one, else guessed.
std::vector<std::string> atomElements(const AtomGroup& atoms);
std::vector<bool> heavyMask(const AtomGroup& atoms);
AtomGroup heavy(const AtomGroup& atoms);
AtomGroup elementFilter(const AtomGroup& atoms, const std::vector<std::string>& elements);

// "a, ,b" -> {"a", "b"}
std::vector<std::string> splitCsv(std::string_view text);

// "G12" per residue, "A:G12" when the residues span several segids.
std::vector<std::string> residueLabels(const Topology& topology, std::span<const std::int32_t> residues);

// "0, 10, 20-50:10, last, -1" -> sorted unique 0-based frames of a trajectory of frameCount frames.
std::vector<std::int64_t> parseFrameSpec(std::string_view text, std::int64_t frameCount);

// "1.5ns, 2000ps, 3" -> times in ps (numbers without unit in defaultUnit).
std::vector<double> parseTimeSpec(std::string_view text, std::string_view defaultUnit = "ps");

// ps per unit (ps, ns, us, fs). Throws ValueError for other units.
double timeUnitFactor(std::string_view unit);

}  // namespace mdbridge::core
