#pragma once
// MDAnalysis' atom selection language (select_atoms), which MDATK passes its selections to and its
// users write: the same tokens, precedence (and / or equal and left-associative, not tighter),
// case-sensitive glob matching, resid ranges with insertion codes, and the same error texts.
// Not chemfiles' language, which differs in keywords, precedence and matching.

#include <cstdint>
#include <string_view>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/topology.hpp"

namespace mdbridge::core {

// A selection MDAnalysis rejects; what() is MDAnalysis' str(exc), e.g. "Unknown selection token:
// 'foo'" or "'Topology' object has no attribute 'chainIDs'".
class AtomSelectionError : public ValueError {
public:
    using ValueError::ValueError;
    const char* pythonName() const noexcept override { return "SelectionError"; }
};

// select_atoms(selection) on all atoms: sorted unique atom indices. Distance keywords (around,
// sphzone, ...) use `current` (positions and, periodic, its box); without a frame they throw.
// An empty string selects nothing (MDAnalysis warns).
std::vector<std::int32_t> selectAtoms(const Topology& topology, const Frame* current, std::string_view selection);

}  // namespace mdbridge::core
