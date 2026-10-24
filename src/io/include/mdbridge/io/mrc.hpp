#pragma once
// Density maps: MRC / CCP4 files as MDATK reads and writes them (mdatk/io/mrc.py), plain or gzip,
// and the map functions of mdatk/density.py that need no fitting (statistics, contour level,
// simulated maps of a model).

#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::io {

// mrc.MapError (a ValueError): unreadable or unsupported maps, with MDATK's messages.
class MapError : public core::ValueError {
public:
    using core::ValueError::ValueError;
    const char* pythonName() const noexcept override { return "MapError"; }
};

// numpy-exact float32 statistics of a map, widened to double.
struct MapStats {
    double min = 0, max = 0, mean = 0, rms = 0;
};

struct DensityMap {
    std::array<std::int64_t, 3> shapeXyz{};  // nx, ny, nz
    std::vector<float> data;  // [z][y][x], x fastest
    core::Vec3 origin{};
    core::Vec3 voxel{};
    std::string path;  // as given to readMap

    float at(std::int64_t x, std::int64_t y, std::int64_t z) const;
    core::Vec3 extent() const;
    MapStats stats() const;
    core::Vec3 indexToXyz(const core::Vec3& index) const;
    core::Vec3 xyzToIndex(const core::Vec3& xyz) const;
};

// read_map: both byte orders, modes 0 / 1 / 2 / 6 / 12, every axis order, the NSTART / ORIGIN /
// cell rules, gzip. Throws MapError (MDATK's messages; truncated or corrupt gzip and absurd sizes
// also as MapError, a documented deviation), std::filesystem::filesystem_error.
DensityMap readMap(const std::filesystem::path& path);

// write_map: float32 MRC2014, little-endian, byte-exact with MDATK.
void writeMap(const std::filesystem::path& path, const DensityMap& map);

// load_map_cached: one map object per path (and modification time), shared.
std::shared_ptr<const DensityMap> loadMapCached(const std::filesystem::path& path);

// map_level(map, sigmas): mean + sigmas * rms.
double mapLevel(const DensityMap& map, double sigmas);

// The ordered table `mdatk density info MAP` prints (as JSON), for the CLI.
core::Value densityInfo(const std::filesystem::path& path);

// simulate_map / atom_weights: a model's density at a resolution (Gaussian splatting).
inline constexpr double kSigmaFactor = 0.22507907903927651;  // 1 / (pi sqrt 2), ChimeraX molmap
std::vector<double> atomWeights(std::span<const std::string> elements);
DensityMap simulateMap(std::span<const core::Vec3> xyz, std::span<const double> weights, double resolution,
                       std::optional<double> spacing = std::nullopt, double sigmaFactor = kSigmaFactor,
                       double padSigmas = 5.0);

}  // namespace mdbridge::io
