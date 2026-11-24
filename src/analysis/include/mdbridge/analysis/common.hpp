#pragma once
// Helpers shared by several analyses (mdatk/analysis/common.py). Residue labels are
// core::residueLabels (common.labels).

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/topology.hpp"

namespace mdbridge::analysis {

struct ResidueRows {
    std::vector<std::int32_t> residues;  // sorted unique residue indices of the atoms
    std::vector<std::int64_t> rows;  // per atom: its residue's position in `residues`
};

// residue_rows(atoms)
ResidueRows residueRows(const core::AtomGroup& atoms);

// without(atoms, other): the atoms of `atoms` that are not in `other`, order and duplicates kept
// (`atoms` itself when `other` is empty).
core::AtomGroup without(const core::AtomGroup& atoms, const core::AtomGroup& other);

// per_row_counts(offsets, rows, values, n_rows): dense float32 rows x frames matrix (frames =
// offsets.size() - 1) of the CSR triplets, values summed in order in float32 (np.add.at).
// Throws std::invalid_argument when the sizes disagree, std::out_of_range for a row >= nRows.
core::NumArray perRowCounts(std::span<const std::int64_t> offsets, std::span<const std::int32_t> rows,
                            std::span<const float> values, std::size_t nRows);

}  // namespace mdbridge::analysis
