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
// (float32 widened, scale_factor applied) - `dt` never changes them; without a time variable
// frame * dt (dt default 1.0). dt(): the given dt, else float32(time[1]) - float32(time[0]), 1.0
// for one frame. Boxes from cell_lengths / cell_angles (rounded to float32), none without them.
// `expectedAtoms`: the topology's atom count (0: take the file's).
std::unique_ptr<core::Trajectory> openNetcdfTrajectory(const std::filesystem::path& path,
                                                       std::optional<double> dt = std::nullopt,
                                                       std::size_t expectedAtoms = 0);

// Writes frames as MDAnalysis' NCDFWriter does without the netCDF4 package (scipy's classic
// writer): the Amber attributes, float32 coordinates, float32 time, float64 cell lengths/angles
// when the frames have boxes. `programVersion` is MDAnalysis' version string in MDATK's files.
struct NetcdfWriteOptions {
    std::string title;
    std::string program = "MDAnalysis";
    std::string programVersion = "2.11.0";
};
void writeNetcdfTrajectory(const std::filesystem::path& path, const std::vector<core::Frame>& frames,
                           std::size_t atoms, const NetcdfWriteOptions& options = {});

}  // namespace mdbridge::io
