#pragma once
// Amber NetCDF trajectories (netCDF classic, 64-bit offset), read as MDAnalysis' NCDFReader reads
// them and written as its NCDFWriter (scipy backend) writes them. A self-contained reader of the
// classic format: no netCDF-C dependency.

#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "mdbridge/core/topology.hpp"

namespace mdbridge::io {

// The Amber convention checks and messages of MDAnalysis' reader ("NCDF trajectory <path> is
// missing Conventions", "... does not conform to AMBER specifications ...", atom-count mismatch,
// units other than angstrom / picosecond) are core::ValueError. Times: the file's time variable
// (in its own type, or the scale_factor product's, widened) - `dt` never changes them; without a
// time variable frame * dt (dt default 1.0). dt(): the given dt, else time[1] - time[0] of the
// unscaled values in the variable's type, 1.0 without a time variable or for one frame. Boxes from cell_lengths / cell_angles (rounded to float32), none without them.
// `expectedAtoms`: the topology's atom count (0: take the file's).
std::unique_ptr<core::Trajectory> openNetcdfTrajectory(const std::filesystem::path& path,
                                                       std::optional<double> dt = std::nullopt,
                                                       std::size_t expectedAtoms = 0);

// Writes frames as MDAnalysis' NCDFWriter does without the netCDF4 package (scipy's classic
// writer), byte for byte: the Amber attributes, float32 coordinates, float32 time, float64 cell
// lengths/angles when the first frame has a box (then every frame needs one). Nothing is written
// for no frames. `program`, `programVersion` and `application` are the global attributes of
// MDATK's files; `title` is not stored (MDAnalysis 2.11's NCDFWriter drops its remarks).
struct NetcdfWriteOptions {
    std::string title;
    std::string program = "MDAnalysis.coordinates.TRJ.NCDFWriter";
    std::string programVersion = "2.11.0";
    std::string application = "MDAnalysis";
};
void writeNetcdfTrajectory(const std::filesystem::path& path, const std::vector<core::Frame>& frames,
                           std::size_t atoms, const NetcdfWriteOptions& options = {});

}  // namespace mdbridge::io
