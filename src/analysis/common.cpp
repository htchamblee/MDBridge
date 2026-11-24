#include "mdbridge/analysis/common.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

namespace mdbridge::analysis {

ResidueRows residueRows(const core::AtomGroup& atoms) {
    ResidueRows out;
    out.residues = atoms.residues();  // np.unique(atoms.resindices)
    out.rows.reserve(atoms.size());
    const core::Topology& topology = atoms.topology();
    for (const std::int32_t atom : atoms.indices()) {
        const std::int32_t residue = topology.atomResidue[static_cast<std::size_t>(atom)];
        const auto it = std::lower_bound(out.residues.begin(), out.residues.end(), residue);
        out.rows.push_back(static_cast<std::int64_t>(it - out.residues.begin()));
    }
    return out;
}

core::AtomGroup without(const core::AtomGroup& atoms, const core::AtomGroup& other) {
    if (other.empty()) {
        return atoms;
    }
    const std::unordered_set<std::int32_t> drop(other.indices().begin(), other.indices().end());
    std::vector<bool> keep;
    keep.reserve(atoms.size());
    for (const std::int32_t atom : atoms.indices()) {
        keep.push_back(!drop.contains(atom));
    }
    return atoms.masked(keep);
}

core::NumArray perRowCounts(std::span<const std::int64_t> offsets, std::span<const std::int32_t> rows,
                            std::span<const float> values, std::size_t nRows) {
    if (offsets.empty()) {
        throw std::invalid_argument("perRowCounts: offsets must hold frames + 1 entries");
    }
    const std::size_t nFrames = offsets.size() - 1;
    if (rows.size() != values.size() || offsets.back() != static_cast<std::int64_t>(rows.size())) {
        throw std::invalid_argument("perRowCounts: the CSR arrays disagree");
    }
    std::vector<float> out(nRows * nFrames, 0.0f);
    for (std::size_t f = 0; f < nFrames; ++f) {
        const std::int64_t begin = offsets[f];
        const std::int64_t end = offsets[f + 1];
        if (begin < 0 || end < begin) {
            throw std::invalid_argument("perRowCounts: offsets must not decrease");
        }
        for (auto e = static_cast<std::size_t>(begin); e < static_cast<std::size_t>(end); ++e) {
            const std::int32_t row = rows[e];
            if (row < 0 || static_cast<std::size_t>(row) >= nRows) {
                throw std::out_of_range("index " + std::to_string(row) + " is out of bounds for axis 0 with size " +
                                        std::to_string(nRows));
            }
            float& cell = out[static_cast<std::size_t>(row) * nFrames + f];
            cell = cell + values[e];
        }
    }
    return core::NumArray(std::move(out), {nRows, nFrames});
}

}  // namespace mdbridge::analysis
