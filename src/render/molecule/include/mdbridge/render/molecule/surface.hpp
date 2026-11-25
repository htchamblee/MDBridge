#pragma once
// Smooth molecular surface (like VMD QuickSurf): MDATK's mdatk/render/molecule/surface.py.
//
// Every atom is splatted trilinearly onto a grid and blurred with a Gaussian, normalised so an
// isolated atom peaks at 1. The iso-surface is extracted with naive surface nets: one vertex per
// grid cell that straddles the level (the mean of its edge crossings) and one quad per grid edge
// that crosses it. Normals are the trilinearly sampled density gradient, so the mesh is smooth with
// about a third of the triangles of marching cubes. The mesh does not depend on the view, so it is
// built once per frame. The arithmetic follows numpy / scipy's (float32 grid, ndimage.hpp).
//
// Deliberate difference (docs/DIFFERENCES.md RN-18): a grid of more than kMaxSurfaceGridPoints
// points (a blown-up frame, or an enormous system) is refused with std::runtime_error before
// anything is allocated, where MDATK tried to allocate it.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "mdbridge/render/molecule/geometry.hpp"
#include "mdbridge/render/molecule/ndimage.hpp"

namespace mdbridge::render::molecule {

inline constexpr double kMaxSurfaceGridPoints = 67108864.0;  // 2^26 (about 0.25 GB of float32)

// Summed unit Gaussians on a grid [z, y, x] with `spacing` (Å).
struct DensityGrid {
    std::vector<float> values;  // C order, shape {nz, ny, nx}
    Shape3 shape{};
    Vec3 origin{};  // Å, (x, y, z) of grid point [0, 0, 0]
    double spacing = 1.0;
};
// pad defaults to 3.2 sigma + 1.5 Å. Throws std::invalid_argument without points, and
// std::runtime_error for grids above kMaxSurfaceGridPoints.
DensityGrid densityGrid(std::span<const Vec3> xyz, double spacing = 1.0, double sigma = 1.4,
                        std::optional<double> pad = std::nullopt);

// Iso-surface vol == level: verts in Å (x, y, z), triangles and the vertices' fractional grid
// coordinates (z, y, x). The grid must be below `level` on its border (densityGrid's padding).
struct SurfaceNets {
    std::vector<Vec3> verts;
    std::vector<Tri> tris;
    std::vector<Vec3> coords;  // (z, y, x) per vertex
};
SurfaceNets surfaceNets(const DensityGrid& grid, double level);

// verts, outward normals, triangles (wound outward) and the nearest atom per vertex.
struct GaussianSurface {
    std::vector<Vec3> verts;
    std::vector<Vec3> normals;
    std::vector<Tri> tris;
    std::vector<std::int32_t> nearest;
};
GaussianSurface gaussianSurface(std::span<const Vec3> xyz, double spacing = 1.0, double sigma = 1.4,
                                double level = 0.5);

// Index of the point of `points` nearest to `query` (Euclidean; the lowest index among equally
// near ones), for many queries.
std::vector<std::int32_t> nearestPoints(std::span<const Vec3> points, std::span<const Vec3> queries);

}  // namespace mdbridge::render::molecule
