#pragma once
// A uniform cell grid for neighbour searches: MDBridge's own code (written for the atom selections'
// around / sphzone / ... keywords, computed independently of MDAnalysis), shared by the selections
// and the analyses' pair searches.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace mdbridge::core::detail {

// Cells of a uniform grid holding reference points (x, y, z interleaved, float or double): every
// point within `reach` of a query lies in the 27 cells around the query's cell. With a periodic
// cell the grid covers the cell in fractional coordinates; otherwise the points' bounding box.
template <class T>
class CellGrid {
public:
    // box: the cell's row vectors a, b, c (row-major; lower triangular: a along x, b in the xy
    // plane) for a periodic grid, nullptr for open boundaries. No grid (std::nullopt) for a
    // non-positive or non-finite reach, no points, a periodic cell less than 3 cells thick along
    // some axis, non-finite points (open), or an open grid that would have far more cells than
    // points; callers then compare every pair.
    static std::optional<CellGrid> build(std::span<const T> ref, double reach,
                                         const std::array<double, 9>* box) {
        CellGrid grid;
        const std::size_t count = ref.size() / 3;
        if (!(reach > 0.0) || !std::isfinite(reach) || count == 0) {
            return std::nullopt;
        }
        const double cell = reach * (1.0 + 1e-4) + 1e-4;
        if (box == nullptr) {
            return buildOpen(ref, cell, grid);
        }
        const std::array<double, 9>& m = *box;
        grid.m_periodic = true;
        grid.m_box = m;
        const double volume = m[0] * m[4] * m[8];
        const double bc = std::sqrt(std::pow(m[4] * m[8], 2) + std::pow(m[3] * m[8], 2) +
                                    std::pow(m[3] * m[7] - m[4] * m[6], 2));
        const double ca = m[0] * std::sqrt(m[8] * m[8] + m[7] * m[7]);
        const double heights[3] = {volume / bc, volume / ca, m[8]};
        for (std::size_t i = 0; i < 3; ++i) {
            const double cells = std::floor(heights[i] / cell);
            if (!(cells >= 3)) {
                return std::nullopt;
            }
            grid.m_dims[i] = static_cast<std::int64_t>(std::min(cells, 256.0));
        }
        grid.fill(ref);
        return grid;
    }

    // Calls visit(reference index) for the candidates near p until it returns true; whether one
    // did. Candidates come cell by cell (x, then y, then z of the 27 cells), in each cell in
    // increasing reference index.
    template <class F>
    bool any(const T* p, F&& visit) const {
        std::array<std::int64_t, 3> c{};
        if (!cellOf(p, c, true)) {
            return false;
        }
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    std::array<std::int64_t, 3> q = {c[0] + dx, c[1] + dy, c[2] + dz};
                    bool inside = true;
                    for (std::size_t i = 0; i < 3; ++i) {
                        if (m_periodic) {
                            q[i] = (q[i] % m_dims[i] + m_dims[i]) % m_dims[i];
                        } else {
                            inside = inside && q[i] >= 0 && q[i] < m_dims[i];
                        }
                    }
                    if (!inside) {
                        continue;
                    }
                    const auto k =
                        static_cast<std::size_t>((q[0] * m_dims[1] + q[1]) * m_dims[2] + q[2]);
                    for (std::size_t s = m_start[k]; s < m_start[k + 1]; ++s) {
                        if (visit(m_items[s])) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    // Calls visit(reference index) for every candidate near p. A periodic grid with fewer than 3
    // cells along an axis is never built, so no candidate repeats.
    template <class F>
    void forEach(const T* p, F&& visit) const {
        any(p, [&visit](std::size_t r) {
            visit(r);
            return false;
        });
    }

private:
    static std::optional<CellGrid> buildOpen(std::span<const T> ref, double cell, CellGrid& grid) {
        const std::size_t count = ref.size() / 3;
        std::array<double, 3> lo{}, hi{};
        for (std::size_t i = 0; i < 3; ++i) {
            lo[i] = std::numeric_limits<double>::infinity();
            hi[i] = -std::numeric_limits<double>::infinity();
        }
        for (std::size_t a = 0; a < count; ++a) {
            for (std::size_t i = 0; i < 3; ++i) {
                const double v = ref[3 * a + i];
                if (!std::isfinite(v)) {
                    return std::nullopt;
                }
                lo[i] = std::min(lo[i], v);
                hi[i] = std::max(hi[i], v);
            }
        }
        // at most 1024 cells per axis: wider cells along an axis that would need more (every cell
        // at least `cell` wide, so the 27 cells around a query still hold every neighbour)
        double total = 1;
        for (std::size_t i = 0; i < 3; ++i) {
            const double extent = hi[i] - lo[i];
            grid.m_cell[i] = std::max(cell, extent / 1000.0);
            grid.m_dims[i] =
                static_cast<std::int64_t>(std::min(std::floor(extent / grid.m_cell[i]) + 1, 1024.0));
            total *= static_cast<double>(grid.m_dims[i]);
        }
        if (total > 4.0 * static_cast<double>(count) + 1e6) {
            return std::nullopt;
        }
        grid.m_low = lo;
        grid.fill(ref);
        return grid;
    }

    bool cellOf(const T* p, std::array<std::int64_t, 3>& c, bool query) const {
        if (m_periodic) {
            const double z = p[2] / m_box[8];
            const double y = (p[1] - z * m_box[7]) / m_box[4];
            const double x = (p[0] - y * m_box[3] - z * m_box[6]) / m_box[0];
            const double s[3] = {x, y, z};
            for (std::size_t i = 0; i < 3; ++i) {
                if (!std::isfinite(s[i])) {
                    return false;
                }
                const double f = s[i] - std::floor(s[i]);
                c[i] = std::min(m_dims[i] - 1,
                                static_cast<std::int64_t>(f * static_cast<double>(m_dims[i])));
            }
            return true;
        }
        for (std::size_t i = 0; i < 3; ++i) {
            const double v = p[i];
            if (!std::isfinite(v)) {
                return false;
            }
            const double k = std::floor((v - m_low[i]) / m_cell[i]);
            if (query && (k < -1 || k > static_cast<double>(m_dims[i]))) {
                return false;
            }
            c[i] = static_cast<std::int64_t>(k);
        }
        return true;
    }

    void fill(std::span<const T> ref) {
        const std::size_t count = ref.size() / 3;
        const auto cells = static_cast<std::size_t>(m_dims[0] * m_dims[1] * m_dims[2]);
        std::vector<std::size_t> cellOfPoint(count);
        m_start.assign(cells + 1, 0);
        for (std::size_t a = 0; a < count; ++a) {
            std::array<std::int64_t, 3> c{};
            cellOf(&ref[3 * a], c, false);
            const auto k = static_cast<std::size_t>((c[0] * m_dims[1] + c[1]) * m_dims[2] + c[2]);
            cellOfPoint[a] = k;
            ++m_start[k + 1];
        }
        for (std::size_t k = 0; k < cells; ++k) {
            m_start[k + 1] += m_start[k];
        }
        m_items.resize(count);
        std::vector<std::size_t> fillAt(m_start.begin(), m_start.end() - 1);
        for (std::size_t a = 0; a < count; ++a) {
            m_items[fillAt[cellOfPoint[a]]++] = a;
        }
    }

    bool m_periodic = false;
    std::array<double, 9> m_box{};
    std::array<double, 3> m_low{};
    std::array<double, 3> m_cell{1, 1, 1};
    std::array<std::int64_t, 3> m_dims{1, 1, 1};
    std::vector<std::size_t> m_start;
    std::vector<std::size_t> m_items;
};

}  // namespace mdbridge::core::detail
