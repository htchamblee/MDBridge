#pragma once
// MDAnalysis' atom selection language (select_atoms), which MDATK passes its selections to and its
// users write: the same tokens, precedence (and / or equal and left-associative, not tighter),
// case-sensitive glob matching, resid ranges with insertion codes, and the same error texts,
// except the deliberate differences of docs/DIFFERENCES.md "Selections": reversed prop operands
// mirror the operator, float range bounds take the single-value syntax, reversed resid ranges
// select nothing, and errors name what is missing instead of Python's internal texts.
// Not chemfiles' language, which differs in keywords, precedence and matching.
// Periodic distances are computed independently of MDAnalysis, exactly for every cutoff; pairs at
// exactly the cutoff may differ.

#include <cstdint>
#include <string_view>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/topology.hpp"

namespace mdbridge::core {

// A rejected selection; what() is MDAnalysis' str(exc) where MDBridge keeps it, e.g. "Unknown
// selection token: 'foo'" or "'Topology' object has no attribute 'chainIDs'", or MDBridge's own
// text ("Unexpected end of selection after 'and'", "around: distance must not be negative").
class AtomSelectionError : public ValueError {
public:
    using ValueError::ValueError;
    const char* pythonName() const noexcept override { return "SelectionError"; }
};

// select_atoms(selection) on all atoms: sorted unique atom indices. Distance keywords (around,
// sphzone, ...) use `current` (positions and, periodic, its box); without a frame they throw.
// An empty string selects nothing (MDAnalysis warns); one of blanks only is "empty selection".
std::vector<std::int32_t> selectAtoms(const Topology& topology, const Frame* current, std::string_view selection);

}  // namespace mdbridge::core
