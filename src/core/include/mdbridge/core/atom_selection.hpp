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
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/topology.hpp"

namespace mdbridge::core {

// A rejected selection; what() is MDAnalysis' str(exc) where MDBridge keeps it, e.g. "Unknown
// selection token: 'foo'" or "'Topology' object has no attribute 'chainIDs'", or MDBridge's own
// text ("Unexpected end of selection after 'and'", "around: distance must not be negative").
// Errors the selection layer adds a hint to (core::select, SEL-19) say what kind they are and
// carry the offending token.
class AtomSelectionError : public ValueError {
public:
    enum class Problem {
        Other,
        UnknownToken,  // "Unknown selection token: 'x'" (token: x)
        UnknownGroup,  // "Selection failed: 'Failed to find group: x'" (token: x)
        UnclosedParenthesis,  // "Unexpected token: 'x' Expected: ')'" (token: x)
    };
    using ValueError::ValueError;
    // `plain`: the message without the parser's own hint (what() has it)
    AtomSelectionError(const std::string& message, Problem problem, std::string token, std::string plain)
        : ValueError(message), m_problem(problem), m_token(std::move(token)), m_plain(std::move(plain)) {}
    const char* pythonName() const noexcept override { return "SelectionError"; }
    Problem problem() const noexcept { return m_problem; }
    const std::string& token() const noexcept { return m_token; }
    std::string plainMessage() const { return m_plain.empty() ? std::string(what()) : m_plain; }

private:
    Problem m_problem = Problem::Other;
    std::string m_token;
    std::string m_plain;
};

// The keywords of the language and "and" / "or", in table order (for typo hints).
const std::vector<std::string>& selectionKeywords();

// The hint for an unknown token: "keywords are lower case: did you mean 'and'?" when its lower-case
// spelling is a keyword, else "did you mean 'resname'?" for a close keyword (suggestName), else ""
// (also for "and" / "or" themselves, unknown where a selection should start).
std::string keywordHint(const std::string& token);

// select_atoms(selection) on all atoms: sorted unique atom indices. Distance keywords (around,
// sphzone, ...) use `current` (positions and, periodic, its box); without a frame they throw.
// An empty string selects nothing (MDAnalysis warns); one of blanks only is "empty selection".
std::vector<std::int32_t> selectAtoms(const Topology& topology, const Frame* current, std::string_view selection);

}  // namespace mdbridge::core
