#pragma once
// Coordinate files as MDAnalysis reads them for MDATK: the reader per extension, times, dt, boxes
// and float32 precision of each format (io-trajectories spec 2.2), and several files concatenated
// with ChainReader's rules (per-file times for formats that store them).

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/topology.hpp"

namespace mdbridge::io {

// MDAnalysis' format names: NCDF, DCD, XTC, TRR, TRJ (Amber ASCII mdcrd), INPCRD (rst7), CRD
// (CHARMM card), PDB, GRO, ... "" when the extension is unknown.
std::string coordinateFormatForExtension(std::string_view path);

// One coordinate file in `format` (MDAnalysis name; "" guesses from the extension). `dt`: the
// project's inputs.dt_ps when set. `atoms`: the topology's atom count, checked against the file.
// Throws core::ValueError (unknown format, malformed file, atom-count mismatch) or
// std::filesystem::filesystem_error.
std::unique_ptr<core::Trajectory> openTrajectory(const std::filesystem::path& path, std::string_view format,
                                                 std::optional<double> dt, std::size_t atoms);

// Several trajectories one after the other (MDAnalysis ChainReader, continuous=False): global
// frame numbers, each file keeping its own times where it stores them; dt() and info() of the
// first file (no info at all: MDAnalysis' ChainReader has no header / remarks).
std::unique_ptr<core::Trajectory> chainTrajectories(std::vector<std::unique_ptr<core::Trajectory>> parts);

}  // namespace mdbridge::io
