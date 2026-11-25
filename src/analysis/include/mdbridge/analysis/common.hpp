#pragma once
// Helpers shared by several analyses (mdatk/analysis/common.py). Residue labels are
// core::residueLabels (common.labels).

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
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

// One sparse frame of entries per analysed frame, kept by frame position so that process() can
// run for several frames at once (analysis.hpp), joined in frame order by csr().
struct FrameEntries {
    std::vector<std::int64_t> i, j;
    std::vector<double> v;
};
class FrameCsr {
public:
    // `frames` empty frames
    void reset(std::size_t frames) { m_frames.assign(frames, {}); }
    // frame k's entries (frames are independent: concurrent calls for different k are safe)
    void set(std::size_t k, FrameEntries entries) { m_frames.at(k) = std::move(entries); }
    // CsrBuilder::add of every frame in order
    core::CsrBuilder csr() const;

private:
    std::vector<FrameEntries> m_frames;
};

}  // namespace mdbridge::analysis
