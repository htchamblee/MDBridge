#include "mdbridge/analysis/spatial.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "mdbridge/core/detail/neighbors.hpp"
#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/json.hpp"
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

// AN-20: the points whose three coordinates are finite (their indices), or std::nullopt when all are
std::optional<std::vector<std::size_t>> finitePoints(std::span<const double> xyz) {
    const std::size_t n = xyz.size() / 3;
    bool all = true;
    for (const double v : xyz) {
        all = all && std::isfinite(v);
    }
    if (all) {
        return std::nullopt;
    }
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < n; ++i) {
        if (std::isfinite(xyz[3 * i]) && std::isfinite(xyz[3 * i + 1]) && std::isfinite(xyz[3 * i + 2])) {
            out.push_back(i);
        }
    }
    return out;
}

std::vector<double> pointsAt(std::span<const double> xyz, const std::vector<std::size_t>& rows) {
    std::vector<double> out;
    out.reserve(3 * rows.size());
    for (const std::size_t i : rows) {
        out.insert(out.end(), xyz.begin() + static_cast<std::ptrdiff_t>(3 * i),
                   xyz.begin() + static_cast<std::ptrdiff_t>(3 * i + 3));
    }
    return out;
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
    // AN-20: points with a non-finite coordinate are in no pair (MDATK: scipy's "data must be
    // finite" error); the search runs on the finite points and the indices are mapped back
    const std::optional<std::vector<std::size_t>> finiteA = finitePoints(a);
    const std::optional<std::vector<std::size_t>> finiteB = finitePoints(b);
    if (finiteA || finiteB) {
        const std::vector<double> fa = finiteA ? pointsAt(a, *finiteA) : std::vector<double>(a.begin(), a.end());
        const std::vector<double> fb = finiteB ? pointsAt(b, *finiteB) : std::vector<double>(b.begin(), b.end());
        out = pairsWithinOpen(fa, fb, cutoff);
        for (std::size_t p = 0; p < out.size(); ++p) {
            if (finiteA) {
                out.i[p] = static_cast<std::int64_t>((*finiteA)[static_cast<std::size_t>(out.i[p])]);
            }
            if (finiteB) {
                out.j[p] = static_cast<std::int64_t>((*finiteB)[static_cast<std::size_t>(out.j[p])]);
            }
        }
        return out;
    }
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

// ---- periodic cells and minimum images (docs/DIFFERENCES.md, AN-01 to AN-04) -------------------
//
// MDBridge's own periodic geometry; it replaces the reproduction of MDAnalysis' capped_distance
// (docs/specs/m3-capped-distance.md documents that reference behaviour). Every step is a float64
// operation in the order written here (no contraction, see above), so that
// tests/parity/improvements_periodic.py reproduces the results bit for bit with numpy.

namespace {

constexpr double kRadiansPerDegree = 3.14159265358979323846 / 180.0;

// At most this many reductions by whole cells along one vector: each one shrinks the component
// by a factor of about 2^-52, so float64 values from float32 inputs never need more than ~25.
constexpr int kMaxReductions = 64;

double cosineOf(double degrees) {
    return degrees == 90.0 ? 0.0 : std::cos(degrees * kRadiansPerDegree);
}

double sineOf(double degrees) {
    return degrees == 90.0 ? 1.0 : std::sin(degrees * kRadiansPerDegree);
}

double squaredLength(const Vec3& v) {
    return (v[0] * v[0] + v[1] * v[1]) + v[2] * v[2];
}

}  // namespace

std::optional<PeriodicCell> PeriodicCell::from(const Box& box) {
    PeriodicCell cell;
    std::array<double, 3> length{}, angle{};
    for (std::size_t k = 0; k < 3; ++k) {
        length[k] = static_cast<float>(box[k]);
        angle[k] = static_cast<float>(box[3 + k]);
    }
    cell.m_orthogonal = angle[0] == 90.0 && angle[1] == 90.0 && angle[2] == 90.0;
    if (cell.m_orthogonal) {
        for (std::size_t k = 0; k < 3; ++k) {
            cell.m_periodic[k] = std::isfinite(length[k]) && length[k] > 0;
            cell.m_vectors[k][k] = cell.m_periodic[k] ? length[k] : 0.0;
        }
        cell.flatten();
        return cell;
    }
    const double a = length[0], b = length[1], c = length[2];
    for (std::size_t k = 0; k < 3; ++k) {
        if (!(std::isfinite(length[k]) && length[k] > 0 && angle[k] > 0 && angle[k] < 180)) {
            return std::nullopt;
        }
    }
    const double cosAlpha = cosineOf(angle[0]), cosBeta = cosineOf(angle[1]), cosGamma = cosineOf(angle[2]);
    const double sinGamma = sineOf(angle[2]);
    const double bx = b * cosGamma;
    const double by = b * sinGamma;
    const double cx = c * cosBeta;
    const double cy = c * (cosAlpha - cosBeta * cosGamma) / sinGamma;
    const double cz = std::sqrt(c * c - cx * cx - cy * cy);  // NaN for a negative radicand
    if (!(by > 0 && cz > 0 && std::isfinite(cz))) {
        return std::nullopt;
    }
    cell.m_vectors = {Vec3{a, 0.0, 0.0}, Vec3{bx, by, 0.0}, Vec3{cx, cy, cz}};
    cell.m_periodic = {true, true, true};
    cell.flatten();
    return cell;
}

std::optional<PeriodicCell> PeriodicCell::from(const std::optional<Box>& box) {
    const std::optional<Box> usable = usableBox(box);
    return usable ? from(*usable) : std::nullopt;
}

std::optional<PeriodicCell> PeriodicCell::from(const std::optional<std::array<float, 6>>& dims) {
    const std::optional<Box> usable = usableBox(dims);
    return usable ? from(*usable) : std::nullopt;
}

void PeriodicCell::flatten() {
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            m_flat[3 * r + c] = m_vectors[r][c];
        }
    }
}

// Written with plain scalars and m_flat (no std::array accessors, no calls): the pair search calls
// it for every candidate pair (1e8 in the tests), also in unoptimised Debug builds. |v| > h is
// written v > h || v < -h (the same for every v, NaN included, as h > 0).
Vec3 PeriodicCell::minimumImage(const Vec3& d) const {
    const double* v = d.data();
    double x = v[0], y = v[1], z = v[2];
    const double* u = m_flat;  // a = u[0..2], b = u[3..5], c = u[6..8]
    if (m_orthogonal) {
        // per periodic axis (length > 0; 0 when not periodic): v - L * nearbyint(v / L), repeated
        // while |v| > L / 2
        const double lx = u[0], ly = u[4], lz = u[8];
        const double hx = lx / 2, hy = ly / 2, hz = lz / 2;
        for (int step = 0; lx > 0 && step < kMaxReductions && (x > hx || x < -hx); ++step) {
            x = x - lx * std::nearbyint(x / lx);
        }
        for (int step = 0; ly > 0 && step < kMaxReductions && (y > hy || y < -hy); ++step) {
            y = y - ly * std::nearbyint(y / ly);
        }
        for (int step = 0; lz > 0 && step < kMaxReductions && (z > hz || z < -hz); ++step) {
            z = z - lz * std::nearbyint(z / lz);
        }
        return {x, y, z};
    }
    const double ax = u[0], bx = u[3], by = u[4], cx = u[6], cy = u[7], cz = u[8];
    // along c, then b, then a: v - n * u for n = nearbyint(v_k / u_k) (half to even), repeated
    // while |v_k| > u_k / 2 (the components of u after k are 0)
    const double hz = cz / 2, hy = by / 2, hx = ax / 2;
    for (int step = 0; step < kMaxReductions && (z > hz || z < -hz); ++step) {
        const double n = std::nearbyint(z / cz);
        x = x - n * cx;
        y = y - n * cy;
        z = z - n * cz;
    }
    for (int step = 0; step < kMaxReductions && (y > hy || y < -hy); ++step) {
        const double n = std::nearbyint(y / by);
        x = x - n * bx;
        y = y - n * by;
    }
    for (int step = 0; step < kMaxReductions && (x > hx || x < -hx); ++step) {
        const double n = std::nearbyint(x / ax);
        x = x - n * ax;
    }
    // the shortest of the 27 neighbouring images; the image itself unless another is strictly
    // shorter (a first, b, c innermost)
    double bestX = x, bestY = y, bestZ = z;
    double bestSquare = (x * x + y * y) + z * z;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                if (i == 0 && j == 0 && k == 0) {
                    continue;
                }
                const double fi = i, fj = j, fk = k;
                const double wx = ((x + fi * ax) + fj * bx) + fk * cx;
                const double wy = (y + fj * by) + fk * cy;
                const double wz = z + fk * cz;
                const double s = (wx * wx + wy * wy) + wz * wz;
                if (s < bestSquare) {
                    bestX = wx;
                    bestY = wy;
                    bestZ = wz;
                    bestSquare = s;
                }
            }
        }
    }
    return {bestX, bestY, bestZ};
}

std::array<double, 3> PeriodicCell::widths() const {
    if (m_orthogonal) {
        return {m_vectors[0][0], m_vectors[1][1], m_vectors[2][2]};
    }
    const Vec3& a = m_vectors[0];
    const Vec3& b = m_vectors[1];
    const Vec3& c = m_vectors[2];
    const auto cross = [](const Vec3& x, const Vec3& y) {
        return Vec3{x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]};
    };
    const double volume = a[0] * b[1] * c[2];
    return {volume / std::sqrt(squaredLength(cross(b, c))), volume / std::sqrt(squaredLength(cross(c, a))),
            volume / std::sqrt(squaredLength(cross(a, b)))};
}

std::array<double, 3> PeriodicCell::fractions(const Vec3& x) const {
    if (m_orthogonal) {
        std::array<double, 3> s{};
        for (std::size_t k = 0; k < 3; ++k) {
            s[k] = m_periodic[k] ? x[k] / m_vectors[k][k] : 0.0;
        }
        return s;
    }
    const Vec3& a = m_vectors[0];
    const Vec3& b = m_vectors[1];
    const Vec3& c = m_vectors[2];
    const double sc = x[2] / c[2];
    const double sb = (x[1] - sc * c[1]) / b[1];
    const double sa = (x[0] - sb * b[0] - sc * c[0]) / a[0];
    return {sa, sb, sc};
}

// ---- periodic pairs ----------------------------------------------------------------------------

namespace {

bool finitePoint(std::span<const double> xyz, std::size_t p) {
    return std::isfinite(xyz[3 * p]) && std::isfinite(xyz[3 * p + 1]) && std::isfinite(xyz[3 * p + 2]);
}

Vec3 pointAt(std::span<const double> xyz, std::size_t p) {
    return {xyz[3 * p], xyz[3 * p + 1], xyz[3 * p + 2]};
}

// Plain distances for a cell without periodic geometry: the open path on the finite points.
PairList plainPairs(std::span<const double> a, std::span<const double> b, double cutoff) {
    std::vector<double> fa, fb;
    std::vector<std::int64_t> ia, ib;
    for (std::size_t p = 0; p < a.size() / 3; ++p) {
        if (finitePoint(a, p)) {
            fa.insert(fa.end(), a.begin() + static_cast<std::ptrdiff_t>(3 * p),
                      a.begin() + static_cast<std::ptrdiff_t>(3 * p + 3));
            ia.push_back(static_cast<std::int64_t>(p));
        }
    }
    for (std::size_t p = 0; p < b.size() / 3; ++p) {
        if (finitePoint(b, p)) {
            fb.insert(fb.end(), b.begin() + static_cast<std::ptrdiff_t>(3 * p),
                      b.begin() + static_cast<std::ptrdiff_t>(3 * p + 3));
            ib.push_back(static_cast<std::int64_t>(p));
        }
    }
    PairList out = pairsWithinOpen(fa, fb, cutoff);
    for (std::size_t k = 0; k < out.size(); ++k) {
        out.i[k] = ia[static_cast<std::size_t>(out.i[k])];
        out.j[k] = ib[static_cast<std::size_t>(out.j[k])];
    }
    return out;
}

}  // namespace

void visitPairsWithinPeriodic(std::span<const double> a, std::span<const double> b, double cutoff,
                              const Box& box, const PairRowVisitor& visit) {
    const std::size_t n = pointCount(a.size(), "pairsWithinPeriodic");
    const std::size_t m = pointCount(b.size(), "pairsWithinPeriodic");
    const double bound = cutoff * cutoff;
    if (n == 0 || m == 0 || std::isnan(bound)) {
        return;
    }
    const std::optional<PeriodicCell> cell = PeriodicCell::from(box);
    if (!cell) {
        const PairList plain = plainPairs(a, b, cutoff);  // sorted by (i, j): one row at a time
        for (std::size_t s = 0; s < plain.size();) {
            std::size_t e = s;
            while (e < plain.size() && plain.i[e] == plain.i[s]) {
                ++e;
            }
            visit(plain.i[s], std::span(plain.j).subspan(s, e - s), std::span(plain.d).subspan(s, e - s));
            s = e;
        }
        return;
    }

    // Candidate pairs from a cell list over the fractional coordinates wrapped into [0, 1): along
    // each periodic direction, bins at least |cutoff| (plus a margin for rounding) wide, so that
    // every image within the cutoff lies in the same or a neighbouring bin. A direction narrower
    // than three bins is one bin (every candidate along it: cutoffs beyond a third of the cell).
    // Points very far outside the cell (fractions beyond 1e4) are compared with every point.
    const std::array<double, 3> width = cell->widths();
    const double reach = std::abs(cutoff);
    std::array<std::int64_t, 3> bins{};
    for (std::size_t k = 0; k < 3; ++k) {
        bins[k] = 1;
        if (cell->periodic()[k] && std::isfinite(reach)) {
            const double size = reach + 1e-7 * (reach + width[k]);
            const double q = std::floor(width[k] / size);
            if (q >= 3) {
                bins[k] = static_cast<std::int64_t>(std::min(q, 256.0));
            }
        }
    }
    const auto binOf = [&](std::span<const double> xyz, std::size_t p, std::int64_t& key) {
        const std::array<double, 3> s = cell->fractions(pointAt(xyz, p));
        std::array<std::int64_t, 3> idx{};
        for (std::size_t k = 0; k < 3; ++k) {
            if (bins[k] == 1) {
                continue;
            }
            if (!(std::abs(s[k]) <= 1e4)) {
                return false;  // far outside: compared with everything
            }
            const double f = s[k] - std::floor(s[k]);
            idx[k] = std::min<std::int64_t>(static_cast<std::int64_t>(f * static_cast<double>(bins[k])), bins[k] - 1);
        }
        key = (idx[0] * bins[1] + idx[1]) * bins[2] + idx[2];
        return true;
    };

    std::vector<std::pair<std::int64_t, std::size_t>> binned;  // (bin, j)
    std::vector<std::size_t> farB;
    for (std::size_t j = 0; j < m; ++j) {
        if (!finitePoint(b, j)) {
            continue;
        }
        std::int64_t key = 0;
        if (binOf(b, j, key)) {
            binned.emplace_back(key, j);
        } else {
            farB.push_back(j);
        }
    }
    std::sort(binned.begin(), binned.end());
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> members;
    for (std::size_t s = 0; s < binned.size();) {
        std::size_t e = s;
        while (e < binned.size() && binned[e].first == binned[s].first) {
            ++e;
        }
        members.emplace(binned[s].first, std::pair(s, e));
        s = e;
    }
    std::vector<std::size_t> binnedJ;  // the j of binned: each bin's members, increasing
    for (const auto& entry : binned) {
        binnedJ.push_back(entry.second);
    }
    std::vector<std::size_t> everyB = binnedJ;
    everyB.insert(everyB.end(), farB.begin(), farB.end());
    std::sort(everyB.begin(), everyB.end());

    std::vector<std::size_t> candidates;
    std::vector<std::int64_t> rowJ;
    std::vector<double> rowD;
    for (std::size_t i = 0; i < n; ++i) {
        if (!finitePoint(a, i)) {
            continue;
        }
        std::int64_t key = 0;
        std::span<const std::size_t> list;  // the candidates j, increasing
        if (!binOf(a, i, key)) {
            list = everyB;
        } else {
            candidates = farB;
            int ranges = farB.empty() ? 0 : 1;  // increasing runs appended
            const std::int64_t home[3] = {key / (bins[1] * bins[2]), (key / bins[2]) % bins[1], key % bins[2]};
            const auto offsets = [&](std::size_t k) {
                return bins[k] == 1 ? std::vector<std::int64_t>{0} : std::vector<std::int64_t>{-1, 0, 1};
            };
            for (const std::int64_t dx : offsets(0)) {
                for (const std::int64_t dy : offsets(1)) {
                    for (const std::int64_t dz : offsets(2)) {
                        const auto wrap = [&](std::size_t k, std::int64_t v) { return (v + bins[k]) % bins[k]; };
                        const std::int64_t neighbour =
                            (wrap(0, home[0] + dx) * bins[1] + wrap(1, home[1] + dy)) * bins[2] + wrap(2, home[2] + dz);
                        const auto found = members.find(neighbour);
                        if (found == members.end()) {
                            continue;
                        }
                        candidates.insert(candidates.end(),
                                          binnedJ.begin() + static_cast<std::ptrdiff_t>(found->second.first),
                                          binnedJ.begin() + static_cast<std::ptrdiff_t>(found->second.second));
                        ++ranges;
                    }
                }
            }
            if (ranges > 1) {
                std::sort(candidates.begin(), candidates.end());
            }
            list = candidates;
        }
        // plain pointers in the loop (it runs for every candidate pair, also in Debug builds)
        const double* p = a.data() + 3 * i;
        const double px = p[0], py = p[1], pz = p[2];
        const double* bxyz = b.data();
        const std::size_t* candidate = list.data();
        const std::size_t count = list.size();
        if (rowJ.size() < count) {
            rowJ.resize(count);
            rowD.resize(count);
        }
        std::int64_t* keptJ = rowJ.data();
        double* keptD = rowD.data();
        std::size_t kept = 0;
        for (std::size_t c = 0; c < count; ++c) {
            const std::size_t j = candidate[c];
            const double* q = bxyz + 3 * j;
            const Vec3 image = cell->minimumImage({px - q[0], py - q[1], pz - q[2]});
            const double* w = image.data();
            const double s = (w[0] * w[0] + w[1] * w[1]) + w[2] * w[2];
            if (s <= bound) {
                keptJ[kept] = static_cast<std::int64_t>(j);
                keptD[kept] = std::sqrt(s);
                ++kept;
            }
        }
        if (kept > 0) {
            visit(static_cast<std::int64_t>(i), std::span(rowJ).first(kept), std::span(rowD).first(kept));
        }
    }
}

PairList pairsWithinPeriodic(std::span<const double> a, std::span<const double> b, double cutoff,
                             const Box& box) {
    PairList out;
    visitPairsWithinPeriodic(a, b, cutoff, box,
                             [&](std::int64_t i, std::span<const std::int64_t> j, std::span<const double> d) {
                                 out.i.insert(out.i.end(), j.size(), i);
                                 out.j.insert(out.j.end(), j.begin(), j.end());
                                 out.d.insert(out.d.end(), d.begin(), d.end());
                             });
    return out;
}

void visitPairsWithinPeriodic(std::span<const float> a, std::span<const float> b, double cutoff,
                              const Box& box, const PairRowVisitor& visit) {
    visitPairsWithinPeriodic(widened(a), widened(b), cutoff, box, visit);
}

PairList pairsWithinPeriodic(std::span<const float> a, std::span<const float> b, double cutoff,
                             const Box& box) {
    return pairsWithinPeriodic(widened(a), widened(b), cutoff, box);
}

PairList pairsWithin(std::span<const double> a, std::span<const double> b, double cutoff,
                     const std::optional<Box>& box) {
    pointCount(a.size(), "pairsWithin");
    pointCount(b.size(), "pairsWithin");
    if (a.empty() || b.empty()) {
        return {};
    }
    if (const std::optional<Box> usable = usableBox(box)) {
        // MDATK passes np.asarray(a, np.float32) to the periodic search
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

std::vector<double> minimumImageVectors(std::span<const double> vectors, const std::optional<Box>& box) {
    const std::size_t n = pointCount(vectors.size(), "minimumImageVectors");
    std::vector<double> out(vectors.begin(), vectors.end());
    const std::optional<PeriodicCell> cell = PeriodicCell::from(box);
    if (!cell) {
        return out;
    }
    for (std::size_t p = 0; p < n; ++p) {
        const Vec3 v = cell->minimumImage({out[3 * p], out[3 * p + 1], out[3 * p + 2]});
        std::copy(v.begin(), v.end(), out.begin() + static_cast<std::ptrdiff_t>(3 * p));
    }
    return out;
}

Vec3 minimumImageVector(const Vec3& vector, const std::optional<Box>& box) {
    const std::optional<PeriodicCell> cell = PeriodicCell::from(box);
    return cell ? cell->minimumImage(vector) : vector;
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
