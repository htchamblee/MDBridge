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
// (CHARMM card), PDB, GRO, ... "" when the extension is unknown. MDATK's table of
// mdatk/io/loader.py first (.nc .ncdf .netcdf NCDF, .mdcrd .crd .trj TRJ, .rst7 .inpcrd INPCRD,
// .ncrst NCRST), else MDAnalysis' guess (the upper-cased extension).
std::string coordinateFormatForExtension(std::string_view path);

// One coordinate file in `format` (MDAnalysis name, any case; "" guesses from the extension as
// MDAnalysis does per file: .rst7 / .netcdf / .ncrst are then unknown, .crd is CHARMM CRD). `dt`:
// the project's inputs.dt_ps when set. `atoms`: the topology's atom count, checked against the
// file (0: not checked). Throws core::ValueError (unknown or unsupported format, malformed file,
// atom-count mismatch; pythonName() names MDAnalysis' exception) or
// std::filesystem::filesystem_error.
std::unique_ptr<core::Trajectory> openTrajectory(const std::filesystem::path& path, std::string_view format,
                                                 std::optional<double> dt, std::size_t atoms);

// The coordinates of Universe(topology, *paths, format=..., dt=...): one file is openTrajectory;
// several are opened as MDAnalysis' ChainReader opens them (a file's ValueError becomes TypeError
// "Unable to read <path> with <class '...'>.", the files must agree on the atom count) and chained,
// then checked against `atoms` with the Universe's message.
std::unique_ptr<core::Trajectory> openTrajectories(const std::vector<std::filesystem::path>& paths,
                                                   std::string_view format, std::optional<double> dt,
                                                   std::size_t atoms);

// Several trajectories one after the other (MDAnalysis ChainReader, continuous=False): global
// frame numbers, each file keeping its own times where it stores them (NCDF with a time variable,
// XTC / TRR, DCD, INPCRD with a time), else global frame * that file's dt. dt() is the dt of the
// file holding the current frame (the first file's after opening); no info() (MDAnalysis'
// ChainReader has no header / remarks).
std::unique_ptr<core::Trajectory> chainTrajectories(std::vector<std::unique_ptr<core::Trajectory>> parts);

}  // namespace mdbridge::io
