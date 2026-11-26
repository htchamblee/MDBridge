// The pair search of cKDTree.query_pairs (scipy/spatial/ckdtree/src/query_pairs.cxx: the pairs whose
// squared distance, sqeuclidean_distance_double, is within r^2) on MDBridge's own cell grid
// (core/detail/neighbors.hpp); no SciPy code. Built without floating-point contraction, so that the
// squared distances round as cKDTree's do.

#include "mdbridge/numerics/pairs.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "mdbridge/core/detail/neighbors.hpp"

namespace mdbridge::numerics {

namespace {

double squaredDistance(const double* a, const double* b) {
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    double s = 0.0;
    s += dx * dx;
    s += dy * dy;
    s += dz * dz;
    return s;
}

}  // namespace

std::vector<std::pair<std::size_t, std::size_t>> queryPairs(std::span<const double> xyz, double r) {
    if (xyz.size() % 3 != 0) {
        throw std::invalid_argument("queryPairs: coordinates come in threes (x, y, z)");
    }
    for (double v : xyz) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("data must be finite, check for nan or inf values");
        }
    }
    std::vector<std::pair<std::size_t, std::size_t>> out;
    if (!(r >= 0)) {
        return out;
    }
    const double r2 = r * r;
    const std::size_t n = xyz.size() / 3;
    const auto grid = core::detail::CellGrid<double>::build(xyz, r, nullptr);
    if (grid) {
        for (std::size_t a = 0; a < n; ++a) {
            const double* p = &xyz[3 * a];
            const std::size_t first = out.size();
            grid->forEach(p, [&](std::size_t b) {
                if (b > a && squaredDistance(p, &xyz[3 * b]) <= r2) {
                    out.emplace_back(a, b);
                }
            });
            std::sort(out.begin() + static_cast<std::ptrdiff_t>(first), out.end());
        }
    } else {  // no grid (a spread with far more cells than points): every pair
        for (std::size_t a = 0; a < n; ++a) {
            for (std::size_t b = a + 1; b < n; ++b) {
                if (squaredDistance(&xyz[3 * a], &xyz[3 * b]) <= r2) {
                    out.emplace_back(a, b);
                }
            }
        }
    }
    return out;
}

}  // namespace mdbridge::numerics
