// density_centroid of mdatk/density.py (the owner's code; moved from render/panels/structure.cpp).

#include "mdbridge/numerics/density.hpp"

#include <cmath>
#include <cstddef>

namespace mdbridge::numerics {

core::Vec3 densityCentroid(const io::DensityMap& map, double level) {
    const auto nx = static_cast<std::size_t>(map.shapeXyz[0]), ny = static_cast<std::size_t>(map.shapeXyz[1]),
               nz = static_cast<std::size_t>(map.shapeXyz[2]);
    double total = 0.0, cx = 0.0, cy = 0.0, cz = 0.0;
    for (std::size_t z = 0; z < nz; ++z) {
        for (std::size_t y = 0; y < ny; ++y) {
            for (std::size_t x = 0; x < nx; ++x) {
                const double v = static_cast<double>(map.data[(z * ny + y) * nx + x]) - level;
                if (!(v > 0.0) || !std::isfinite(v)) continue;
                total += v;
                cx += v * static_cast<double>(x);
                cy += v * static_cast<double>(y);
                cz += v * static_cast<double>(z);
            }
        }
    }
    if (!(total > 0.0)) {
        return {map.origin[0] + 0.5 * static_cast<double>(map.shapeXyz[0] - 1) * map.voxel[0],
                map.origin[1] + 0.5 * static_cast<double>(map.shapeXyz[1] - 1) * map.voxel[1],
                map.origin[2] + 0.5 * static_cast<double>(map.shapeXyz[2] - 1) * map.voxel[2]};
    }
    return map.indexToXyz({cx / total, cy / total, cz / total});
}

}  // namespace mdbridge::numerics
