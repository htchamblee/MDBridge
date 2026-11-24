#pragma once
// Density maps: MRC / CCP4 files as MDATK reads and writes them (mdatk/io/mrc.py), plain or gzip,
// and the map functions of mdatk/density.py that need no fitting (statistics, contour level,
// simulated maps of a model).
//
// Deliberate deviations from MDATK (where Python raises other exceptions or allocates first):
// - Truncated gzip data ends the data instead of raising EOFError: MapError "{path}: file too
//   short for an MRC header" or "{path}: file is shorter than its header says". Corrupt deflate
//   data raises MapError "{path}: Error -3 while decompressing data: <zlib's message>" (Python:
//   zlib.error without the path); a bad member trailer or a damaged following member raises
//   MapError "{path}: CRC check failed 0x.. != 0x..", "{path}: Incorrect length of data
//   produced", "{path}: Not a gzipped file (b'..')" or "{path}: Unknown compression method"
//   (Python: gzip.BadGzipFile). Members are decompressed only as far as the map needs, as in
//   Python, but a member's trailer is checked only when the map needs data beyond it; Python's
//   128 KiB read-ahead may check it (and raise) earlier.
// - Header sizes beyond the file (up to 2^93 bytes) raise "{path}: file is shorter than its header
//   says" without allocating them (Python: MemoryError or that message).
// - Statistics: min / max of maps holding both +0.0 and -0.0 as the extreme may differ from numpy
//   in the sign of the zero (numpy's SIMD reduction order depends on the CPU); a NaN result keeps
//   the first NaN's bits.
// - simulateMap uses the C library's exp(), as numpy does without SIMD exp kernels; numpy on
//   AVX-512 CPUs uses Intel's SVML exp, which differs in the last bit for a few percent of
//   arguments, so such maps could differ in the last float32 bit of a voxel (none of the parity
//   cases does after the float32 rounding).

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

// density.FitError (a ValueError): simulateMap's "no atoms to simulate a map from" and "the
// resolution must be positive".
class FitError : public core::ValueError {
public:
    using core::ValueError::ValueError;
    const char* pythonName() const noexcept override { return "FitError"; }
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

// load_map_cached: one map object per path (and size and modification time), shared; the two
// most recently read maps are kept.
std::shared_ptr<const DensityMap> loadMapCached(const std::filesystem::path& path);

// map_level(map, sigmas): mean + sigmas * rms.
double mapLevel(const DensityMap& map, double sigmas);

// The ordered table `mdatk density info MAP` prints (as JSON), for the CLI.
core::Value densityInfo(const std::filesystem::path& path);

// simulate_map / atom_weights: a model's density at a resolution (Gaussian splatting). Elements
// are looked up as str(e).strip().upper() in MDATK's table (unknown ones count as carbon).
inline constexpr double kSigmaFactor = 0.22507907903927651;  // 1 / (pi sqrt 2), ChimeraX molmap
std::vector<double> atomWeights(std::span<const std::string> elements);
DensityMap simulateMap(std::span<const core::Vec3> xyz, std::span<const double> weights,
                       double resolution, std::optional<double> spacing = std::nullopt,
                       double sigmaFactor = kSigmaFactor, double padSigmas = 5.0);

}  // namespace mdbridge::io
