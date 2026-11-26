#pragma once
// scipy.spatial.cKDTree(points).query_pairs(r) for 3-D points (MDATK's helix builder: contacts of a
// closure and of the contact report, mdatk/builder.py). A cell grid (core/detail/neighbors.hpp), not a
// k-d tree: the same set of pairs, returned sorted rather than in the tree's order (cKDTree returns a
// Python set, or with output_type="ndarray" an array in an order that depends on the tree; MDBridge
// documents a fixed order instead, so that callers' results do not depend on it).

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace mdbridge::numerics {

// Every pair (i, j), i < j, of the points (x, y, z interleaved) whose squared distance, summed as
// ((dx^2 + dy^2) + dz^2) as cKDTree sums it, is at most r * r; sorted by i, then j. A non-finite
// coordinate throws std::invalid_argument with cKDTree's message ("data must be finite, check for
// nan or inf values"); r < 0 or NaN gives no pairs (cKDTree squares r, so a negative radius found the
// pairs within |r|: docs/DIFFERENCES.md NU-01). Single-threaded and deterministic.
std::vector<std::pair<std::size_t, std::size_t>> queryPairs(std::span<const double> xyz, double r);

}  // namespace mdbridge::numerics
