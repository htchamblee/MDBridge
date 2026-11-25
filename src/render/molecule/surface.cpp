// Port of MDATK's mdatk/render/molecule/surface.py (MDATK v0.0.4.1, the owner's code). numpy's and
// scipy's rounding is followed where it decides the mesh: the splat is summed in float64 in numpy's
// order and stored as float32, the Gaussian blur is scipy's (ndimage.cpp), the level and the edge
// crossings are float32 (numpy 2 computes `float32 array op Python float` in float32), and the
// gradient is np.gradient's float32 central / one-sided differences, interpolated by scipy's
// order-1 map_coordinates. scipy.spatial.cKDTree's nearest-atom query is replaced by a cell grid
// (ties: the lowest atom index).

#include "mdbridge/render/molecule/surface.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <unordered_map>


namespace mdbridge::render::molecule {
namespace {

// cube corners as (dz, dy, dx) bits (corner c = 4 dz + 2 dy + dx) and the 12 cube edges between them
constexpr std::array<std::array<int, 3>, 8> kCorners = {
    {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}};
constexpr std::array<std::array<int, 2>, 12> kEdges = {
    {{0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7}}};

std::string wholeNumber(double x) {
    char text[64];
    std::snprintf(text, sizeof text, "%.0f", x);
    return text;
}

}  // namespace

DensityGrid densityGrid(std::span<const Vec3> xyz, double spacing, double sigma, std::optional<double> pad_opt) {
    if (xyz.empty()) {
        throw std::invalid_argument("zero-size array to reduction operation minimum which has no identity");
    }
    const double h = spacing;
    const double pad = pad_opt ? *pad_opt : 3.2 * sigma + 1.5;
    Vec3 lo{}, hi{};
    for (std::size_t c = 0; c < 3; ++c) {
        double mn = xyz[0][c];
        double mx = xyz[0][c];
        for (const Vec3& p : xyz) {
            mn = std::min(mn, p[c]);
            mx = std::max(mx, p[c]);
        }
        lo[c] = mn - pad;
        hi[c] = mx;
    }
    double dims[3];
    for (std::size_t c = 0; c < 3; ++c) {
        dims[c] = std::ceil((hi[c] + pad - lo[c]) / h) + 1.0;
        if (!std::isfinite(dims[c]) || dims[c] < 1.0) {
            throw std::runtime_error("surface grid: coordinates or spacing are not usable");
        }
    }
    const double total = dims[0] * dims[1] * dims[2];
    if (total > kMaxSurfaceGridPoints) {
        throw std::runtime_error("surface grid of " + wholeNumber(dims[0]) + " x " + wholeNumber(dims[1]) + " x " +
                                 wholeNumber(dims[2]) + " points is too large (at most 67108864 points)");
    }
    const auto nx = static_cast<std::size_t>(dims[0]);
    const auto ny = static_cast<std::size_t>(dims[1]);
    const auto nz = static_cast<std::size_t>(dims[2]);
    std::vector<double> acc(nx * ny * nz, 0.0);
    const std::size_t n = xyz.size();
    std::vector<std::array<std::int64_t, 3>> i0(n);
    std::vector<Vec3> f(n);
    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t c = 0; c < 3; ++c) {
            const double g = (xyz[a][c] - lo[c]) / h;
            const double fl = std::floor(g);
            i0[a][c] = static_cast<std::int64_t>(fl);
            f[a][c] = g - static_cast<double>(i0[a][c]);
        }
    }
    const auto sny = static_cast<std::int64_t>(ny);
    const auto snx = static_cast<std::int64_t>(nx);
    for (int dz = 0; dz < 2; ++dz) {
        for (int dy = 0; dy < 2; ++dy) {
            for (int dx = 0; dx < 2; ++dx) {
                for (std::size_t a = 0; a < n; ++a) {
                    const double wz = dz ? f[a][2] : 1 - f[a][2];
                    const double wy = dy ? f[a][1] : 1 - f[a][1];
                    const double wx = dx ? f[a][0] : 1 - f[a][0];
                    const std::int64_t idx = ((i0[a][2] + dz) * sny + (i0[a][1] + dy)) * snx + (i0[a][0] + dx);
                    acc[static_cast<std::size_t>(idx)] += wx * wy * wz;
                }
            }
        }
    }
    DensityGrid grid;
    grid.shape = {nz, ny, nx};
    grid.origin = lo;
    grid.spacing = h;
    grid.values.resize(acc.size());
    for (std::size_t i = 0; i < acc.size(); ++i) {
        grid.values[i] = static_cast<float>(acc[i]);
    }
    acc.clear();
    acc.shrink_to_fit();
    const double s = sigma / h;
    gaussianFilter(grid.values, grid.shape, s, 3.0, 0.0);
    const auto norm = static_cast<float>(std::pow(2 * std::numbers::pi * s * s, 1.5));
    for (float& v : grid.values) {
        v *= norm;
    }
    return grid;
}

SurfaceNets surfaceNets(const DensityGrid& grid, double level) {
    const float lvl = static_cast<float>(level);
    const auto [nz, ny, nx] = grid.shape;
    SurfaceNets out;
    if (nz < 2 || ny < 2 || nx < 2) {
        return out;
    }
    const std::vector<float>& vol = grid.values;
    auto vidx = [&](std::size_t z, std::size_t y, std::size_t x) { return (z * ny + y) * nx + x; };
    auto inside = [&](std::size_t z, std::size_t y, std::size_t x) { return vol[vidx(z, y, x)] >= lvl; };
    const std::size_t cz = nz - 1, cy = ny - 1, cx = nx - 1;
    std::vector<std::int32_t> vid(cz * cy * cx, -1);
    std::vector<std::size_t> cells;
    for (std::size_t z = 0; z < cz; ++z) {
        for (std::size_t y = 0; y < cy; ++y) {
            for (std::size_t x = 0; x < cx; ++x) {
                int n_in = 0;
                for (const auto& c : kCorners) {
                    n_in += inside(z + static_cast<std::size_t>(c[0]), y + static_cast<std::size_t>(c[1]),
                                   x + static_cast<std::size_t>(c[2]))
                                ? 1
                                : 0;
                }
                if (n_in > 0 && n_in < 8) {
                    vid[(z * cy + y) * cx + x] = static_cast<std::int32_t>(cells.size());
                    cells.push_back((z * cy + y) * cx + x);
                }
            }
        }
    }
    if (cells.empty()) {
        return out;
    }
    out.verts.reserve(cells.size());
    out.coords.reserve(cells.size());
    for (std::size_t cell : cells) {
        const std::size_t kz = cell / (cy * cx);
        const std::size_t rem = cell % (cy * cx);
        const std::size_t ky = rem / cx;
        const std::size_t kx = rem % cx;
        float cv[8];
        bool cm[8];
        for (std::size_t c = 0; c < 8; ++c) {
            cv[c] = vol[vidx(kz + static_cast<std::size_t>(kCorners[c][0]), ky + static_cast<std::size_t>(kCorners[c][1]),
                             kx + static_cast<std::size_t>(kCorners[c][2]))];
            cm[c] = cv[c] >= lvl;
        }
        double acc[3] = {0.0, 0.0, 0.0};
        double cnt = 0.0;
        for (const auto& [a, b] : kEdges) {
            if (cm[a] == cm[b]) {
                continue;
            }
            const float t = (lvl - cv[a]) / (cv[b] - cv[a]);
            for (std::size_t k = 0; k < 3; ++k) {
                const auto ca = static_cast<double>(kCorners[static_cast<std::size_t>(a)][k]);
                const auto cb = static_cast<double>(kCorners[static_cast<std::size_t>(b)][k]);
                acc[k] += ca + static_cast<double>(t) * (cb - ca);
            }
            cnt += 1.0;
        }
        const double gz = static_cast<double>(kz) + acc[0] / cnt;
        const double gy = static_cast<double>(ky) + acc[1] / cnt;
        const double gx = static_cast<double>(kx) + acc[2] / cnt;
        out.verts.push_back({grid.origin[0] + gx * grid.spacing, grid.origin[1] + gy * grid.spacing,
                             grid.origin[2] + gz * grid.spacing});
        out.coords.push_back({gz, gy, gx});
    }
    auto cellId = [&](std::size_t z, std::size_t y, std::size_t x) { return vid[(z * cy + y) * cx + x]; };
    // every grid edge that crosses the level is shared by four cells: one quad
    std::vector<std::array<std::int32_t, 4>> quads;
    auto addQuad = [&](bool in, std::array<std::int32_t, 4> q) {
        quads.push_back(in ? q : std::array<std::int32_t, 4>{q[3], q[2], q[1], q[0]});
    };
    for (std::size_t z = 1; z + 1 < nz; ++z) {  // x edges
        for (std::size_t y = 1; y + 1 < ny; ++y) {
            for (std::size_t x = 0; x + 1 < nx; ++x) {
                const bool in = inside(z, y, x);
                if (in != inside(z, y, x + 1)) {
                    addQuad(in, {cellId(z - 1, y - 1, x), cellId(z - 1, y, x), cellId(z, y, x), cellId(z, y - 1, x)});
                }
            }
        }
    }
    for (std::size_t z = 1; z + 1 < nz; ++z) {  // y edges
        for (std::size_t y = 0; y + 1 < ny; ++y) {
            for (std::size_t x = 1; x + 1 < nx; ++x) {
                const bool in = inside(z, y, x);
                if (in != inside(z, y + 1, x)) {
                    addQuad(in, {cellId(z - 1, y, x - 1), cellId(z, y, x - 1), cellId(z, y, x), cellId(z - 1, y, x)});
                }
            }
        }
    }
    for (std::size_t z = 0; z + 1 < nz; ++z) {  // z edges
        for (std::size_t y = 1; y + 1 < ny; ++y) {
            for (std::size_t x = 1; x + 1 < nx; ++x) {
                const bool in = inside(z, y, x);
                if (in != inside(z + 1, y, x)) {
                    addQuad(in, {cellId(z, y - 1, x - 1), cellId(z, y - 1, x), cellId(z, y, x), cellId(z, y, x - 1)});
                }
            }
        }
    }
    out.tris.reserve(2 * quads.size());
    for (const auto& q : quads) {
        out.tris.push_back({q[0], q[1], q[2]});
    }
    for (const auto& q : quads) {
        out.tris.push_back({q[0], q[2], q[3]});
    }
    return out;
}

GaussianSurface gaussianSurface(std::span<const Vec3> xyz, double spacing, double sigma, double level) {
    const DensityGrid grid = densityGrid(xyz, spacing, sigma);
    SurfaceNets nets = surfaceNets(grid, level);
    GaussianSurface out;
    out.verts = std::move(nets.verts);
    if (nets.tris.empty()) {
        return out;
    }
    out.tris = std::move(nets.tris);
    const std::vector<float>& vol = grid.values;
    const std::size_t strides[3] = {grid.shape[1] * grid.shape[2], grid.shape[2], 1};
    // np.gradient(vol)[axis] at a flat index: float32 central differences, one-sided at the borders
    auto gradient = [&](std::size_t axis, std::size_t i) -> float {
        const std::size_t len = grid.shape[axis];
        const std::size_t p = (i / strides[axis]) % len;
        const std::size_t st = strides[axis];
        if (len < 2) {
            return 0.0F;
        }
        if (p == 0) {
            return (vol[i + st] - vol[i]) / 1.0F;
        }
        if (p == len - 1) {
            return (vol[i] - vol[i - st]) / 1.0F;
        }
        return (vol[i + st] - vol[i - st]) / 2.0F;
    };
    out.normals.reserve(out.verts.size());
    const auto floor_norm = static_cast<float>(1e-9);
    for (const Vec3& c : nets.coords) {
        float g[3];
        for (std::size_t k = 0; k < 3; ++k) {  // components x, y, z = gradients along axes 2, 1, 0
            const std::size_t axis = 2 - k;
            g[k] = -mapCoordinatesLinearNearest([&](std::size_t i) { return gradient(axis, i); }, grid.shape, c[0],
                                                c[1], c[2]);
        }
        const float len = std::max(std::sqrt((g[0] * g[0] + g[1] * g[1]) + g[2] * g[2]), floor_norm);
        out.normals.push_back({static_cast<double>(g[0] / len), static_cast<double>(g[1] / len),
                               static_cast<double>(g[2] / len)});
    }
    out.nearest = nearestPoints(xyz, out.verts);
    // wind each triangle so its face normal agrees with the gradient normals
    for (Tri& t : out.tris) {
        const Vec3& v0 = out.verts[static_cast<std::size_t>(t[0])];
        const Vec3& v1 = out.verts[static_cast<std::size_t>(t[1])];
        const Vec3& v2 = out.verts[static_cast<std::size_t>(t[2])];
        const Vec3 e1{v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        const Vec3 e2{v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
        const Vec3 fn{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]};
        double d[3];
        for (std::size_t k = 0; k < 3; ++k) {
            const float sum = (static_cast<float>(out.normals[static_cast<std::size_t>(t[0])][k]) +
                               static_cast<float>(out.normals[static_cast<std::size_t>(t[1])][k])) +
                              static_cast<float>(out.normals[static_cast<std::size_t>(t[2])][k]);
            d[k] = fn[k] * static_cast<double>(sum);
        }
        if ((d[0] + d[1]) + d[2] < 0) {
            std::swap(t[0], t[2]);
        }
    }
    return out;
}

std::vector<std::int32_t> nearestPoints(std::span<const Vec3> points, std::span<const Vec3> queries) {
    std::vector<std::int32_t> out(queries.size(), -1);
    if (points.empty()) {
        return out;
    }
    auto d2 = [](const Vec3& a, const Vec3& b) {
        const double dx = a[0] - b[0];
        const double dy = a[1] - b[1];
        const double dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    };
    auto better = [](double d, std::int32_t i, double best, std::int32_t best_i) {
        return d < best || (d == best && i < best_i);
    };
    if (points.size() < 32) {
        for (std::size_t q = 0; q < queries.size(); ++q) {
            double best = std::numeric_limits<double>::infinity();
            std::int32_t best_i = -1;
            for (std::size_t i = 0; i < points.size(); ++i) {
                const double d = d2(points[i], queries[q]);
                if (best_i < 0 || better(d, static_cast<std::int32_t>(i), best, best_i)) {
                    best = d;
                    best_i = static_cast<std::int32_t>(i);
                }
            }
            out[q] = best_i;
        }
        return out;
    }
    // a hashed cell grid (cells of 4 Å), searched in growing shells around the query's cell
    constexpr double cell = 4.0;
    Vec3 lo = points[0];
    for (const Vec3& p : points) {
        for (std::size_t c = 0; c < 3; ++c) {
            lo[c] = std::min(lo[c], p[c]);
        }
    }
    auto cellOf = [&](const Vec3& p) {
        std::array<std::int64_t, 3> k{};
        for (std::size_t c = 0; c < 3; ++c) {
            k[c] = static_cast<std::int64_t>(std::floor((p[c] - lo[c]) / cell));
        }
        return k;
    };
    auto key = [](std::int64_t x, std::int64_t y, std::int64_t z) {
        return (static_cast<std::uint64_t>(x) * 73856093ULL) ^ (static_cast<std::uint64_t>(y) * 19349663ULL) ^
               (static_cast<std::uint64_t>(z) * 83492791ULL);
    };
    struct Bucket {
        std::array<std::int64_t, 3> cell;
        std::vector<std::int32_t> members;
    };
    std::unordered_multimap<std::uint64_t, std::size_t> index;
    std::vector<Bucket> buckets;
    std::array<std::int64_t, 3> kmin{0, 0, 0}, kmax{0, 0, 0};
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto k = cellOf(points[i]);
        for (std::size_t c = 0; c < 3; ++c) {
            kmax[c] = std::max(kmax[c], k[c]);
        }
        const std::uint64_t h = key(k[0], k[1], k[2]);
        std::size_t found = buckets.size();
        for (auto [it, end] = index.equal_range(h); it != end; ++it) {
            if (buckets[it->second].cell == k) {
                found = it->second;
                break;
            }
        }
        if (found == buckets.size()) {
            index.emplace(h, buckets.size());
            buckets.push_back({k, {}});
        }
        buckets[found].members.push_back(static_cast<std::int32_t>(i));
    }
    for (std::size_t q = 0; q < queries.size(); ++q) {
        const Vec3& p = queries[q];
        if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
            continue;
        }
        const auto k = cellOf(p);
        // shells beyond this one hold no cells
        std::int64_t r_max = 0;
        for (std::size_t c = 0; c < 3; ++c) {
            r_max = std::max({r_max, std::abs(k[c] - kmin[c]), std::abs(k[c] - kmax[c])});
        }
        double best = std::numeric_limits<double>::infinity();
        std::int32_t best_i = -1;
        for (std::int64_t r = 0; r <= r_max; ++r) {
            for (std::int64_t x = k[0] - r; x <= k[0] + r; ++x) {
                for (std::int64_t y = k[1] - r; y <= k[1] + r; ++y) {
                    const bool edge_xy = std::abs(x - k[0]) == r || std::abs(y - k[1]) == r;
                    const std::int64_t step = edge_xy ? 1 : 2 * r;
                    for (std::int64_t z = k[2] - r; z <= k[2] + r; z += std::max<std::int64_t>(step, 1)) {
                        const std::array<std::int64_t, 3> kc{x, y, z};
                        for (auto [it, end] = index.equal_range(key(x, y, z)); it != end; ++it) {
                            const Bucket& b = buckets[it->second];
                            if (b.cell != kc) {
                                continue;
                            }
                            for (std::int32_t i : b.members) {
                                const double d = d2(points[static_cast<std::size_t>(i)], p);
                                if (best_i < 0 || better(d, i, best, best_i)) {
                                    best = d;
                                    best_i = i;
                                }
                            }
                        }
                    }
                }
            }
            // every point not yet seen is at least r cells away
            const double reach = static_cast<double>(r) * cell;
            if (best_i >= 0 && best < reach * reach) {
                break;
            }
        }
        out[q] = best_i;
    }
    return out;
}

}  // namespace mdbridge::render::molecule
