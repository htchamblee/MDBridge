#include "mdbridge/analysis/spatial.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "mdbridge/core/detail/neighbors.hpp"
#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/linalg.hpp"

// Compiled without floating-point contraction (see CMakeLists.txt): distances, centres and
// products are rounded as scipy and numpy round them.

namespace mdbridge::analysis {

namespace {

std::size_t pointCount(std::size_t values, const char* what) {
    if (values % 3 != 0) {
        throw std::invalid_argument(std::string(what) + ": coordinates must come in x, y, z triples");
    }
    return values / 3;
}

std::vector<double> widened(std::span<const float> xyz) {
    return {xyz.begin(), xyz.end()};
}

std::vector<float> narrowed(std::span<const double> xyz) {
    std::vector<float> out(xyz.size());
    std::transform(xyz.begin(), xyz.end(), out.begin(),
                   [](double v) { return static_cast<float>(v); });
    return out;
}

void requireFinite(std::span<const double> xyz) {
    for (const double v : xyz) {
        if (!std::isfinite(v)) {
            throw core::ValueError("data must be finite, check for nan or inf values");
        }
    }
}

// np.sign, with MDATK's 0 -> 1
double signOrOne(double x) {
    if (x > 0) {
        return 1.0;
    }
    if (x < 0) {
        return -1.0;
    }
    return std::isnan(x) ? x : 1.0;
}

}  // namespace

// ---- boxes -------------------------------------------------------------------------------------

std::optional<Box> usableBox(std::span<const double> box) {
    if (box.size() != 6) {
        return std::nullopt;
    }
    Box out{};
    std::copy(box.begin(), box.end(), out.begin());
    return usableBox(std::optional<Box>(out));
}

std::optional<Box> usableBox(const std::optional<Box>& box) {
    if (!box) {
        return std::nullopt;
    }
    for (std::size_t k = 0; k < 3; ++k) {
        if ((*box)[k] <= 0) {  // NaN passes, as np.any(box[:3] <= 0) lets it
            return std::nullopt;
        }
    }
    return box;
}

std::optional<Box> usableBox(const std::optional<std::array<float, 6>>& dims) {
    if (!dims) {
        return std::nullopt;
    }
    Box out{};
    std::copy(dims->begin(), dims->end(), out.begin());
    return usableBox(std::optional<Box>(out));
}

// ---- pairs_within ------------------------------------------------------------------------------

PairList pairsWithinOpen(std::span<const double> a, std::span<const double> b, double cutoff) {
    const std::size_t na = pointCount(a.size(), "pairsWithin");
    const std::size_t nb = pointCount(b.size(), "pairsWithin");
    PairList out;
    if (na == 0 || nb == 0) {
        return out;
    }
    requireFinite(a);  // cKDTree(a) is built first, then cKDTree(b): the same message either way
    requireFinite(b);
    // scipy keeps d2 <= r ** 2 (RectRectDistanceTracker: upper_bound = r * r for p = 2)
    const double bound = cutoff * cutoff;
    if (std::isnan(bound)) {
        return out;
    }
    const auto keep = [&](std::size_t i, std::size_t j) {
        // sqeuclidean_distance_double for 3 dimensions: the 4-way unrolled loop does not run,
        // s = (0 + 0 + 0 + 0), then s += d * d per coordinate
        double s = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            const double diff = a[3 * i + k] - b[3 * j + k];
            s += diff * diff;
        }
        if (s <= bound) {
            out.i.push_back(static_cast<std::int64_t>(i));
            out.j.push_back(static_cast<std::int64_t>(j));
            out.d.push_back(std::sqrt(s));
        }
    };
    const std::optional<core::detail::CellGrid<double>> grid =
        core::detail::CellGrid<double>::build(b, std::abs(cutoff), nullptr);
    if (!grid) {
        for (std::size_t i = 0; i < na; ++i) {
            for (std::size_t j = 0; j < nb; ++j) {
                keep(i, j);
            }
        }
        return out;
    }
    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < na; ++i) {
        candidates.clear();
        grid->forEach(&a[3 * i], [&candidates](std::size_t j) { candidates.push_back(j); });
        std::sort(candidates.begin(), candidates.end());
        for (const std::size_t j : candidates) {
            keep(i, j);
        }
    }
    return out;
}

PairList pairsWithinPeriodic(std::span<const float> /*a*/, std::span<const float> /*b*/,
                             double /*cutoff*/, const Box& /*box*/) {
    // W1.1 implements MDAnalysis' capped_distance from docs/specs/m3-capped-distance.md (clean room)
    throw std::logic_error("periodic neighbour search not implemented yet");
}

PairList pairsWithin(std::span<const double> a, std::span<const double> b, double cutoff,
                     const std::optional<Box>& box) {
    pointCount(a.size(), "pairsWithin");
    pointCount(b.size(), "pairsWithin");
    if (a.empty() || b.empty()) {
        return {};
    }
    if (const std::optional<Box> usable = usableBox(box)) {
        // MDATK passes np.asarray(a, np.float32) to capped_distance
        return pairsWithinPeriodic(narrowed(a), narrowed(b), cutoff, *usable);
    }
    return pairsWithinOpen(a, b, cutoff);
}

PairList pairsWithin(std::span<const float> a, std::span<const float> b, double cutoff,
                     const std::optional<Box>& box) {
    pointCount(a.size(), "pairsWithin");
    pointCount(b.size(), "pairsWithin");
    if (a.empty() || b.empty()) {
        return {};
    }
    if (const std::optional<Box> usable = usableBox(box)) {
        return pairsWithinPeriodic(a, b, cutoff, *usable);
    }
    return pairsWithinOpen(widened(a), widened(b), cutoff);
}

// ---- minimum_image_vectors ---------------------------------------------------------------------

namespace {

// usable and np.allclose(box[3:], 90.0): |x - 90| <= atol + rtol * |90| (NaN and inf fail)
bool orthorhombic(const std::optional<Box>& box) {
    if (!usableBox(box)) {
        return false;
    }
    constexpr double tolerance = 1e-8 + 1e-5 * 90.0;
    for (std::size_t k = 3; k < 6; ++k) {
        if (!(std::abs((*box)[k] - 90.0) <= tolerance)) {
            return false;
        }
    }
    return true;
}

}  // namespace

std::vector<double> minimumImageVectors(std::span<const double> vectors,
                                        const std::optional<Box>& box) {
    const std::size_t n = pointCount(vectors.size(), "minimumImageVectors");
    std::vector<double> out(vectors.begin(), vectors.end());
    if (!orthorhombic(box)) {
        return out;
    }
    for (std::size_t p = 0; p < n; ++p) {
        for (std::size_t k = 0; k < 3; ++k) {
            const double length = (*box)[k];
            double& v = out[3 * p + k];
            v = v - length * std::nearbyint(v / length);  // np.round: half to even
        }
    }
    return out;
}

Vec3 minimumImageVector(const Vec3& vector, const std::optional<Box>& box) {
    const std::vector<double> out = minimumImageVectors(vector, box);
    return {out[0], out[1], out[2]};
}

// ---- kabsch ------------------------------------------------------------------------------------

KabschResult kabsch(std::span<const Vec3> mobile, std::span<const Vec3> reference,
                    std::optional<std::span<const double>> weights) {
    const std::size_t n = mobile.size();
    if (reference.size() != n || (weights && weights->size() != n)) {
        throw std::invalid_argument("kabsch: mobile, reference and weights differ in length");
    }
    std::vector<double> w(n);
    if (!weights) {
        std::fill(w.begin(), w.end(), 1.0 / static_cast<double>(std::max<std::size_t>(1, n)));
    } else {
        const double total = core::detail::numpySum(*weights);
        for (std::size_t k = 0; k < n; ++k) {
            w[k] = (*weights)[k] / total;
        }
    }
    KabschResult out;
    out.rotation = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    // w @ mobile, w @ reference
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t c = 0; c < 3; ++c) {
            out.mobileCenter[c] += w[k] * mobile[k][c];
            out.referenceCenter[c] += w[k] * reference[k][c];
        }
    }
    if (n < 3) {
        return out;
    }
    // h = ((mobile - mc) * w[:, None]).T @ (reference - rc)
    Mat3 h{};
    for (std::size_t k = 0; k < n; ++k) {
        Vec3 x{}, y{};
        for (std::size_t c = 0; c < 3; ++c) {
            x[c] = (mobile[k][c] - out.mobileCenter[c]) * w[k];
            y[c] = reference[k][c] - out.referenceCenter[c];
        }
        for (std::size_t p = 0; p < 3; ++p) {
            for (std::size_t q = 0; q < 3; ++q) {
                h[p][q] += x[p] * y[q];
            }
        }
    }
    for (const auto& row : h) {
        for (const double v : row) {
            if (!std::isfinite(v)) {
                throw LinAlgError("SVD did not converge");
            }
        }
    }
    const core::Svd3 svd = core::svd3x3(h);
    // d = sign(det(vt.T @ u.T))
    Mat3 vu{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            vu[i][j] = svd.vt[0][i] * svd.u[j][0] + svd.vt[1][i] * svd.u[j][1] +
                       svd.vt[2][i] * svd.u[j][2];
        }
    }
    const double d = signOrOne(core::det3(vu));
    // r = (vt.T @ diag(1, 1, d)) @ u.T
    for (std::size_t i = 0; i < 3; ++i) {
        const Vec3 v{svd.vt[0][i], svd.vt[1][i], svd.vt[2][i] * d};
        for (std::size_t j = 0; j < 3; ++j) {
            out.rotation[i][j] = v[0] * svd.u[j][0] + v[1] * svd.u[j][1] + v[2] * svd.u[j][2];
        }
    }
    return out;
}

}  // namespace mdbridge::analysis
