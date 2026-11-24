#include "mdbridge/analysis/spatial.hpp"

#include <algorithm>
#include <cmath>
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

// ---- periodic pairs (MDAnalysis' capped_distance) ----------------------------------------------
//
// Written from the behaviour specification docs/specs/m3-capped-distance.md (clean room); the
// section numbers below refer to it. Single-precision steps use float variables: each float
// operation rounds to binary32, each double one to binary64 (no contraction, see above).

namespace {

constexpr double kRadiansPerDegree = 3.14159265358979323846 / 180.0;

using Point = std::array<float, 3>;

// The whole number of cells a far point is moved by. MDAnalysis holds it in a 32-bit integer and
// its result is undefined beyond that (§4.1: |x| beyond about 2^31 L); MDBridge then leaves the
// point where it is along that axis (0 cells), which is what MDAnalysis returned in the one probed
// case (a 1e-10 Å triclinic length, points 1-9 Å from the origin).
double wholeCells(double fraction) {
    const double n = std::floor(fraction);
    return std::abs(n) < 0x1p31 ? n : 0.0;
}

// §2: the cell rounded to float32; for a cell that is not orthogonal, its lower-triangular box
// vectors a = (ax, 0, 0), b = (bx, by, 0), c = (cx, cy, cz), all zero when the cell is invalid
struct PeriodicCell {
    std::array<float, 3> length{};
    std::array<float, 3> angle{};
    bool orthogonal = false;
    bool valid = false;  // orthogonal, or a triclinic cell with box vectors
    float ax = 0, bx = 0, by = 0, cx = 0, cy = 0, cz = 0;
};

PeriodicCell periodicCell(const Box& box) {
    PeriodicCell cell;
    for (std::size_t k = 0; k < 3; ++k) {
        cell.length[k] = static_cast<float>(box[k]);
        cell.angle[k] = static_cast<float>(box[3 + k]);
    }
    cell.orthogonal = cell.angle[0] == 90.0f && cell.angle[1] == 90.0f && cell.angle[2] == 90.0f;
    if (cell.orthogonal) {
        cell.valid = true;
        cell.ax = cell.length[0];
        cell.by = cell.length[1];
        cell.cz = cell.length[2];
        return cell;
    }
    const double a = cell.length[0], b = cell.length[1], c = cell.length[2];
    const double alpha = cell.angle[0], beta = cell.angle[1], gamma = cell.angle[2];
    const auto cosine = [](double degrees) {
        return degrees == 90.0 ? 0.0 : std::cos(degrees * kRadiansPerDegree);
    };
    const double cosAlpha = cosine(alpha), cosBeta = cosine(beta), cosGamma = cosine(gamma);
    const double sinGamma = gamma == 90.0 ? 1.0 : std::sin(gamma * kRadiansPerDegree);
    const double bx = b * cosGamma;
    const double by = b * sinGamma;
    const double cx = c * cosBeta;
    const double cy = c * (cosAlpha - cosBeta * cosGamma) / sinGamma;
    const double cz = std::sqrt(c * c - cx * cx - cy * cy);  // NaN for a negative radicand
    const bool positive = a > 0 && b > 0 && c > 0 && alpha > 0 && beta > 0 && gamma > 0;
    const bool below180 = alpha < 180 && beta < 180 && gamma < 180;
    cell.valid = positive && below180 && cz > 0;
    if (cell.valid) {
        cell.ax = static_cast<float>(a);
        cell.bx = static_cast<float>(bx);
        cell.by = static_cast<float>(by);
        cell.cx = static_cast<float>(cx);
        cell.cy = static_cast<float>(cy);
        cell.cz = static_cast<float>(cz);
    }
    return cell;
}

// §3: the automatic method choice
bool useBruteForce(std::size_t n, std::size_t m, double cutoff, const PeriodicCell& cell) {
    if (n < 10 || m < 10) {
        return true;
    }
    if (static_cast<std::uint64_t>(n) * static_cast<std::uint64_t>(m) >= 100000000u) {
        return false;
    }
    std::array<float, 3> size = cell.length;
    if (!cell.orthogonal) {
        // per axis, the largest minus the smallest component of the three box vectors
        const auto extent = [](float u, float v, float w) {
            return std::max({u, v, w}) - std::min({u, v, w});
        };
        size = {extent(cell.ax, cell.bx, cell.cx), extent(0.0f, cell.by, cell.cy),
                extent(0.0f, 0.0f, cell.cz)};
    }
    const float c = static_cast<float>(cutoff);
    return std::any_of(size.begin(), size.end(), [c](float s) { return c > 0.3f * s; });
}

bool finitePoint(std::span<const float> xyz, std::size_t p) {
    return std::isfinite(xyz[3 * p]) && std::isfinite(xyz[3 * p + 1]) &&
           std::isfinite(xyz[3 * p + 2]);
}

void addPair(PairList& out, std::size_t i, std::size_t j, double d) {
    out.i.push_back(static_cast<std::int64_t>(i));
    out.j.push_back(static_cast<std::int64_t>(j));
    out.d.push_back(d);
}

// §4.1, brute force on an orthogonal cell: the difference scaled by the float32 reciprocal of
// the length, reduced to the nearest image (round half away from zero), scaled back
PairList bruteForceOrthogonal(std::span<const float> a, std::span<const float> b, double cutoff,
                              const PeriodicCell& cell) {
    const std::size_t n = a.size() / 3, m = b.size() / 3;
    std::array<double, 3> length{}, reciprocal{};
    std::array<bool, 3> periodic{};
    for (std::size_t k = 0; k < 3; ++k) {
        length[k] = cell.length[k];
        reciprocal[k] = static_cast<float>(1.0 / length[k]);
        periodic[k] = length[k] > 0x1p-23;
    }
    PairList out;
    for (std::size_t i = 0; i < n; ++i) {
        if (!finitePoint(a, i)) {
            continue;
        }
        for (std::size_t j = 0; j < m; ++j) {
            if (!finitePoint(b, j)) {
                continue;
            }
            std::array<double, 3> square{};
            for (std::size_t k = 0; k < 3; ++k) {
                const float difference = b[3 * j + k] - a[3 * i + k];
                double d = difference;
                if (periodic[k]) {
                    const double s = d * reciprocal[k];
                    d = length[k] * (s - std::round(s));
                }
                square[k] = d * d;
            }
            const double distance = std::sqrt((square[0] + square[1]) + square[2]);
            if (distance <= cutoff) {
                addPair(out, i, j, distance);
            }
        }
    }
    return out;
}

// §4.1, grid: x wrapped into [0, L) in single precision; an infinite length is not periodic
float wrapOrthogonal(float x, float length) {
    if (!std::isfinite(length)) {
        return x;
    }
    const double inverse = 1.0 / static_cast<double>(length);
    const auto reduced = [&](float v) {
        const double n = wholeCells(static_cast<double>(v) * inverse);
        const float multiple = static_cast<float>(n * static_cast<double>(length));
        return v - multiple;
    };
    float v = x;
    if (x < 0.0f) {
        const float y = x + length;
        if (!(y >= 0.0f)) {
            float r = reduced(x);
            if (r < 0.0f) {
                r = r + length;
            }
            return r;
        }
        v = y;
    }
    if (v >= length) {
        const float w = v - length;
        if (w < length) {
            return w;
        }
        float r = reduced(v);
        if (r >= length) {
            r = r - length;
        }
        return r;
    }
    return v;
}

// §4.3 step 1: a point moved into the primary cell by whole box vectors, c first (0 <= z < cz),
// then b (bounds tilted by c), then a (bounds tilted by b and c). Along each vector: a point less
// than one cell outside moves by that vector once; one further out first moves by floor of its
// fractional coordinate (float32 multiple) and then by one more vector if still outside.
Point wrapTriclinic(Point p, const PeriodicCell& cell) {
    const Point va{cell.ax, 0.0f, 0.0f}, vb{cell.bx, cell.by, 0.0f}, vc{cell.cx, cell.cy, cell.cz};
    const double cyz = static_cast<double>(cell.cy) / cell.cz;
    const double cxz = static_cast<double>(cell.cx) / cell.cz;
    const double bxy = static_cast<double>(cell.bx) / cell.by;
    const auto plus = [](Point q, const Point& v) {
        for (std::size_t k = 0; k < 3; ++k) {
            q[k] = q[k] + v[k];
        }
        return q;
    };
    const auto minus = [](Point q, const Point& v) {
        for (std::size_t k = 0; k < 3; ++k) {
            q[k] = q[k] - v[k];
        }
        return q;
    };
    const auto minusMultiple = [](Point q, const Point& v, double n) {
        for (std::size_t k = 0; k < 3; ++k) {
            const float multiple = static_cast<float>(n * static_cast<double>(v[k]));
            q[k] = q[k] - multiple;
        }
        return q;
    };
    // k: the axis of the vector's own component; lower(q): the lower bound of q[k]
    const auto along = [&](Point q, const Point& v, std::size_t k, const auto& lower) {
        const double extent = v[k];
        if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2])) {
            return q;  // an infinite cell length: not periodic along this vector
        }
        if (static_cast<double>(q[k]) < lower(q)) {
            const Point once = plus(q, v);
            if (!(static_cast<double>(once[k]) >= lower(once))) {
                q = minusMultiple(q, v, wholeCells((q[k] - lower(q)) / extent));
                if (static_cast<double>(q[k]) < lower(q)) {
                    q = plus(q, v);
                }
                return q;
            }
            q = once;
        }
        if (static_cast<double>(q[k]) >= lower(q) + extent) {
            const Point once = minus(q, v);
            if (static_cast<double>(once[k]) < lower(once) + extent) {
                return once;
            }
            q = minusMultiple(q, v, wholeCells((q[k] - lower(q)) / extent));
            if (static_cast<double>(q[k]) >= lower(q) + extent) {
                q = minus(q, v);
            }
        }
        return q;
    };
    p = along(p, vc, 2, [](const Point&) { return 0.0; });
    p = along(p, vb, 1, [cyz](const Point& q) { return q[2] * cyz; });
    p = along(p, va, 0, [cyz, cxz, bxy](const Point& q) {
        return (q[1] - q[2] * cyz) * bxy + q[2] * cxz;
    });
    return p;
}

// §4.3 steps 2-3: the squared length of the shortest of the 27 images of p - q (the first
// strictly shorter one, i outermost, k innermost); an invalid cell's zero vectors give |p - q|^2.
// A box vector with an infinite component (an infinite cell length) is not periodic: only its
// 0 multiple is tried, and it adds nothing (MDAnalysis' result is meaningless there: distinct
// points at distance 0).
double squaredTriclinic(const Point& p, const Point& q, const PeriodicCell& cell) {
    const float fx = p[0] - q[0], fy = p[1] - q[1], fz = p[2] - q[2];
    const double dx = fx, dy = fy, dz = fz;
    const double ax = cell.ax, bx = cell.bx, by = cell.by, cx = cell.cx, cy = cell.cy,
                 cz = cell.cz;
    const double spanA = std::isfinite(ax) ? 1.0 : 0.0;
    const double spanB = std::isfinite(bx) && std::isfinite(by) ? 1.0 : 0.0;
    const double spanC = std::isfinite(cx) && std::isfinite(cy) && std::isfinite(cz) ? 1.0 : 0.0;
    const auto times = [](double n, double v) { return n == 0.0 ? 0.0 : n * v; };
    double best = 0.0;
    bool first = true;
    for (double i = -spanA; i <= spanA; ++i) {
        for (double j = -spanB; j <= spanB; ++j) {
            for (double k = -spanC; k <= spanC; ++k) {
                const double x = dx + times(i, ax) + times(j, bx) + times(k, cx);
                const double y = dy + times(j, by) + times(k, cy);
                const double z = dz + times(k, cz);
                const double s = (x * x + y * y) + z * z;
                if (first || s < best) {
                    best = s;
                    first = false;
                }
            }
        }
    }
    return best;
}

// §4.1, grid: the squared minimum-image distance of wrapped points (reference minus configuration)
double squaredOrthogonal(const Point& p, const Point& q, const std::array<float, 3>& length) {
    std::array<double, 3> square{};
    for (std::size_t k = 0; k < 3; ++k) {
        const float difference = p[k] - q[k];
        double d = difference;
        const double box = length[k];
        if (d > box / 2) {
            d -= box;
        } else if (d <= -box / 2) {
            d += box;
        }
        square[k] = d * d;
    }
    return (square[0] + square[1]) + square[2];
}

// both point sets wrapped as the method needs (orthogonal per axis, triclinic by box vectors,
// nothing for an invalid cell's brute force); non-finite points are left as they are (skipped)
std::vector<Point> wrappedPoints(std::span<const float> xyz, const PeriodicCell& cell,
                                 bool orthogonal) {
    std::vector<Point> out(xyz.size() / 3);
    for (std::size_t p = 0; p < out.size(); ++p) {
        Point q{xyz[3 * p], xyz[3 * p + 1], xyz[3 * p + 2]};
        if (finitePoint(xyz, p)) {
            if (orthogonal) {
                for (std::size_t k = 0; k < 3; ++k) {
                    q[k] = wrapOrthogonal(q[k], cell.length[k]);
                }
            } else if (cell.valid) {
                q = wrapTriclinic(q, cell);
            }
        }
        out[p] = q;
    }
    return out;
}

// §4.3, brute force on a triclinic or invalid cell
PairList bruteForceTriclinic(std::span<const float> a, std::span<const float> b, double cutoff,
                             const PeriodicCell& cell) {
    const std::vector<Point> pa = wrappedPoints(a, cell, false);
    const std::vector<Point> pb = wrappedPoints(b, cell, false);
    PairList out;
    for (std::size_t i = 0; i < pa.size(); ++i) {
        if (!finitePoint(a, i)) {
            continue;
        }
        for (std::size_t j = 0; j < pb.size(); ++j) {
            if (!finitePoint(b, j)) {
                continue;
            }
            const double distance = std::sqrt(squaredTriclinic(pa[i], pb[j], cell));
            if (distance <= cutoff) {
                addPair(out, i, j, distance);
            }
        }
    }
    return out;
}

// §3.1: the largest cutoff the grid accepts (float32)
float gridLimit(const PeriodicCell& cell) {
    if (cell.orthogonal) {
        return std::min({cell.length[0], cell.length[1], cell.length[2]}) / 2.0f;
    }
    const auto sine = [](float degrees) {
        const float radians = static_cast<float>(static_cast<double>(degrees) * kRadiansPerDegree);
        return static_cast<float>(std::sin(static_cast<double>(radians)));
    };
    const float sa = sine(cell.angle[0]), sb = sine(cell.angle[1]), sg = sine(cell.angle[2]);
    const std::array<float, 6> products = {cell.by * sa, cell.cz * sa, cell.ax * sb,
                                           cell.cz * sb, cell.ax * sg, cell.by * sg};
    float smallest = products[0];
    for (const float v : products) {
        smallest = std::min(smallest, v);
    }
    return smallest / 2.0f;
}

// §3.1, §4.1, §4.3, §5: the grid method
PairList gridSearch(std::span<const float> a, std::span<const float> b, double cutoff,
                    const PeriodicCell& cell) {
    if (cutoff < 0) {
        throw core::ValueError("Cutoff must be positive");
    }
    if (cell.length[0] == 0.0f || cell.length[1] == 0.0f || cell.length[2] == 0.0f) {
        throw core::ValueError("Any of the box dimensions cannot be 0");
    }
    const float limit = gridLimit(cell);
    if (cutoff > static_cast<double>(limit)) {
        throw core::ValueError("Cutoff " + core::pythonFloatRepr(cutoff) +
                               " too large for box (max " +
                               core::pythonFloatRepr(static_cast<double>(limit)) + ")");
    }
    PairList out;
    if (std::isnan(cutoff)) {
        return out;
    }
    const double bound = cutoff * cutoff;
    // an invalid triclinic cell is wrapped and searched as an orthogonal one with its lengths
    const bool orthogonal = cell.orthogonal || !cell.valid;
    const std::vector<Point> pa = wrappedPoints(a, cell, orthogonal);
    const std::vector<Point> pb = wrappedPoints(b, cell, orthogonal);

    // cells per axis: the cell's perpendicular widths (from the box vectors; 0 for an invalid
    // cell) over max(cutoff, 1) + 0.001; up to 3 cells is one cell, not wrapped
    const double ax = cell.ax, bx = cell.bx, by = cell.by, cx = cell.cx, cy = cell.cy,
                 cz = cell.cz;
    std::array<double, 3> width{};
    if (cell.orthogonal) {
        width = {ax, by, cz};
    } else if (cell.valid) {
        const double volume = ax * by * cz;
        const double t = bx * cy - by * cx;
        width = {volume / std::sqrt((by * cz) * (by * cz) + (bx * cz) * (bx * cz) + t * t),
                 by * cz / std::sqrt(cy * cy + cz * cz), cz};
    }
    const float reach = static_cast<float>(std::max(cutoff, 1.0));
    const float h = static_cast<float>(static_cast<double>(reach) + 0.001);
    std::array<std::int64_t, 3> count{};
    std::array<bool, 3> wrapped{};
    for (std::size_t k = 0; k < 3; ++k) {
        const double q = std::floor(width[k] / h);
        wrapped[k] = q >= 4;
        count[k] = wrapped[k] ? static_cast<std::int64_t>(std::min(q, 1290.0)) : 1;
    }
    // the grid cell of a wrapped point: integer parts truncated, then the remainder modulo the
    // count with the dividend's sign; for an orthogonal search the tilt terms vanish
    const std::array<double, 3> extent =
        orthogonal ? std::array<double, 3>{cell.length[0], cell.length[1], cell.length[2]}
                   : std::array<double, 3>{ax, by, cz};
    const double tiltYz = orthogonal ? 0.0 : cy / cz;
    const double tiltXy = orthogonal ? 0.0 : bx / by;
    const double tiltXz = orthogonal ? 0.0 : cx / cz;
    const auto index = [&](double v, std::size_t k) -> std::int64_t {
        const double size = extent[k] / static_cast<double>(count[k]);
        const double r = std::fmod(std::trunc(v / size), static_cast<double>(count[k]));
        return std::isfinite(r) ? static_cast<std::int64_t>(r) : 0;
    };
    using Cell = std::array<std::int64_t, 3>;
    const auto cellOf = [&](const Point& p) -> Cell {
        const double x = p[0], y = p[1], z = p[2];
        const double y1 = y - z * tiltYz;
        return {index(x - y1 * tiltXy - z * tiltXz, 0), index(y1, 1), index(z, 2)};
    };
    const auto key = [](const Cell& c) {
        return ((c[0] + 2048) * 4096 + (c[1] + 2048)) * 4096 + (c[2] + 2048);
    };

    // configuration points by cell, each cell in descending index order
    std::vector<std::pair<std::int64_t, std::size_t>> sorted;
    sorted.reserve(pb.size());
    for (std::size_t j = 0; j < pb.size(); ++j) {
        if (finitePoint(b, j)) {
            sorted.emplace_back(key(cellOf(pb[j])), j);
        }
    }
    std::sort(sorted.begin(), sorted.end(), [](const auto& l, const auto& r) {
        return l.first != r.first ? l.first < r.first : l.second > r.second;
    });
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> members;
    for (std::size_t s = 0; s < sorted.size();) {
        std::size_t e = s;
        while (e < sorted.size() && sorted[e].first == sorted[s].first) {
            ++e;
        }
        members.emplace(sorted[s].first, std::pair(s, e));
        s = e;
    }

    // each reference point in turn; its 27 neighbour cells x outermost, z innermost (offsets
    // -1, 0, +1; only 0 along an axis that is not wrapped)
    for (std::size_t i = 0; i < pa.size(); ++i) {
        if (!finitePoint(a, i)) {
            continue;
        }
        const Cell home = cellOf(pa[i]);
        const auto neighbour = [&](std::size_t k, std::int64_t offset) {
            std::int64_t v = home[k] + offset;
            if (v == -1) {
                v = count[k] - 1;
            } else if (v == count[k]) {
                v = 0;
            }
            return v;
        };
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            if (!wrapped[0] && dx != 0) {
                continue;
            }
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                if (!wrapped[1] && dy != 0) {
                    continue;
                }
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    if (!wrapped[2] && dz != 0) {
                        continue;
                    }
                    const auto found =
                        members.find(key({neighbour(0, dx), neighbour(1, dy), neighbour(2, dz)}));
                    if (found == members.end()) {
                        continue;
                    }
                    for (std::size_t s = found->second.first; s < found->second.second; ++s) {
                        const std::size_t j = sorted[s].second;
                        const double d2 = orthogonal ? squaredOrthogonal(pa[i], pb[j], cell.length)
                                                     : squaredTriclinic(pa[i], pb[j], cell);
                        if (d2 <= bound) {
                            addPair(out, i, j, std::sqrt(d2));
                        }
                    }
                }
            }
        }
    }
    return out;
}

}  // namespace

PairList pairsWithinPeriodic(std::span<const float> a, std::span<const float> b, double cutoff,
                             const Box& box) {
    const std::size_t n = pointCount(a.size(), "pairsWithinPeriodic");
    const std::size_t m = pointCount(b.size(), "pairsWithinPeriodic");
    if (n == 0 || m == 0) {
        return {};
    }
    const PeriodicCell cell = periodicCell(box);
    if (useBruteForce(n, m, cutoff, cell)) {
        return cell.orthogonal ? bruteForceOrthogonal(a, b, cutoff, cell)
                               : bruteForceTriclinic(a, b, cutoff, cell);
    }
    return gridSearch(a, b, cutoff, cell);
}

bool periodicSearchUsesGrid(std::size_t n, std::size_t m, double cutoff, const Box& box) {
    return !useBruteForce(n, m, cutoff, periodicCell(box));
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
