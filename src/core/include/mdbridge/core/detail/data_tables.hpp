#pragma once
// Residue, atom-name, element and mass tables rebuilt from public sources (wwPDB CCD, PDB format,
// CIAAW, Amber / GROMACS manuals, OpenMM force-field files; see data/sources/ATTRIBUTION.md),
// shared by the core (selection macros) and io (topology guessing, PDB writer) libraries. Internal
// to MDBridge: not part of the public API.

#include <algorithm>
#include <span>
#include <string_view>

namespace mdbridge::core::detail::tables {

#include "data_tables_independent.inc"

// The generator writes every name list sorted (byte order), so membership is a binary search.
constexpr bool contains(std::span<const std::string_view> sorted, std::string_view value) {
    return std::binary_search(sorted.begin(), sorted.end(), value);
}

static_assert(std::is_sorted(std::begin(kIndProteinResidues), std::end(kIndProteinResidues)));
static_assert(std::is_sorted(std::begin(kIndNucleicResidues), std::end(kIndNucleicResidues)));
static_assert(std::is_sorted(std::begin(kIndWaterResidues), std::end(kIndWaterResidues)));
static_assert(std::is_sorted(std::begin(kIndProteinBackboneAtoms), std::end(kIndProteinBackboneAtoms)));
static_assert(std::is_sorted(std::begin(kIndNucleicBackboneAtoms), std::end(kIndNucleicBackboneAtoms)));
static_assert(std::is_sorted(std::begin(kIndNucleicBaseAtoms), std::end(kIndNucleicBaseAtoms)));
static_assert(std::is_sorted(std::begin(kIndNucleicSugarAtoms), std::end(kIndNucleicSugarAtoms)));
static_assert(std::is_sorted(std::begin(kIndGuessOneLetterElements), std::end(kIndGuessOneLetterElements)));
static_assert(std::is_sorted(std::begin(kIndGuessTwoLetterElements), std::end(kIndGuessTwoLetterElements)));
static_assert(std::is_sorted(std::begin(kIndMonatomicIonNames), std::end(kIndMonatomicIonNames)));
static_assert(std::is_sorted(std::begin(kIndMercuryHgResidues), std::end(kIndMercuryHgResidues)));
static_assert(std::is_sorted(std::begin(kIndChlorineClResidues), std::end(kIndChlorineClResidues)));
static_assert(std::size(kIndElements) == 118 && kIndElements[0].z == 1 && kIndElements[117].z == 118);

// ASCII case-insensitive equality (element symbols, residue and atom names).
constexpr bool equalsIgnoringCase(std::string_view a, std::string_view b) {
    if (a.size() != b.size()) {
        return false;
    }
    const auto up = [](char c) { return c >= 'a' && c <= 'z' ? static_cast<char>(c - 'a' + 'A') : c; };
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (up(a[i]) != up(b[i])) {
            return false;
        }
    }
    return true;
}

// The element with this symbol (exact spelling: "Ca", not "CA"); nullptr when there is none.
constexpr const IndElement* elementBySymbol(std::string_view symbol) {
    for (const IndElement& e : kIndElements) {
        if (e.symbol == symbol) {
            return &e;
        }
    }
    return nullptr;
}

// The element whose symbol matches in any letter case ("CA", "Ca", "ca"); nullptr when none does.
constexpr const IndElement* elementBySymbolIgnoringCase(std::string_view symbol) {
    for (const IndElement& e : kIndElements) {
        if (equalsIgnoringCase(e.symbol, symbol)) {
            return &e;
        }
    }
    return nullptr;
}

}  // namespace mdbridge::core::detail::tables
