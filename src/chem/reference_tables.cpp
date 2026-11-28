// The generated reference tables (tools/gen_chem_tables.py) behind mdbridge/chem/reference_tables.hpp.

#include "mdbridge/chem/reference_tables.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>

namespace mdbridge::chem::reference {

namespace {

#include "component_tables.inc"
#include "helix_tables.inc"
#include "triplex_tables.inc"

// The contiguous run of rows whose key equals `key` (the tables are grouped by their first field).
template <typename Row, std::size_t N, typename Key>
std::span<const Row> run(const Row (&rows)[N], std::string_view key, Key keyOf) {
    const Row* first = std::find_if(rows, rows + N, [&](const Row& row) { return keyOf(row) == key; });
    const Row* last = std::find_if(first, rows + N, [&](const Row& row) { return keyOf(row) != key; });
    return {first, static_cast<std::size_t>(last - first)};
}

}  // namespace

std::span<const CylindricalAtom> arnottBackbone() { return kArnottBackbone; }

std::span<const CylindricalBaseAtom> arnottBase(char base) {
    const char key[1] = {base};
    return run(kArnottBases, std::string_view(key, 1),
               [](const CylindricalBaseAtom& row) { return std::string_view(&row.base, 1); });
}

std::span<const CylindricalAtom> jrdNucleotide() { return kJrdNucleotide; }
double jrdTwist() { return kJrdTwist; }
double jrdRise() { return kJrdRise; }

std::span<const TriplexAtom> triplexLevel(std::string_view level, int strand) {
    const auto inLevel = [&](const TriplexAtom& a) { return a.level == level && a.strand == strand; };
    const TriplexAtom* first = std::find_if(std::begin(kTriplexLevels), std::end(kTriplexLevels), inLevel);
    const TriplexAtom* last =
        std::find_if(first, std::end(kTriplexLevels), [&](const TriplexAtom& a) { return !inLevel(a); });
    return {first, static_cast<std::size_t>(last - first)};
}

const TriplexScrew* triplexScrew(std::string_view kind) {
    for (const TriplexScrew& s : kTriplexScrews) {
        if (s.kind == kind) {
            return &s;
        }
    }
    return nullptr;
}

std::span<const TriplexHbond> triplexHbonds() { return kTriplexHbonds; }
std::span<const TriplexPoseSpread> triplexPoseSpreads() { return kTriplexPoseSpreads; }

std::span<const ComponentInfo> components() { return kComponents; }

std::span<const ComponentAtom> componentAtoms(std::string_view code) {
    return run(kComponentAtoms, code, [](const ComponentAtom& row) { return row.code; });
}

std::span<const ComponentBond> componentBonds(std::string_view code) {
    return run(kComponentBonds, code, [](const ComponentBond& row) { return row.code; });
}

std::span<const ComponentHydrogen> componentHydrogens(std::string_view code) {
    return run(kComponentHydrogens, code, [](const ComponentHydrogen& row) { return row.code; });
}

std::span<const HydrogenPosition> hydrogenPositions(std::string_view residue) {
    return run(kHydrogenPositions, residue, [](const HydrogenPosition& row) { return row.residue; });
}

std::span<const ChiralCentre> lRiboseChirality() { return kLChirality; }
const ProtonatedBase& protonatedCytidine() { return kProtonatedCytidine; }

}  // namespace mdbridge::chem::reference
