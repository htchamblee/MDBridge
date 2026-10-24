#pragma once
// The per-format readers behind openTrajectory / readTopologyFile (internal to io). Each follows
// MDAnalysis 2.11.0's reader of that format for what MDATK sees: times, dt, boxes, float32 values,
// attributes, and the reader metadata the PDB export writes.

#include <filesystem>
#include <memory>
#include <optional>

#include "mdbridge/core/topology.hpp"
#include "mdbridge/io/loader.hpp"

namespace mdbridge::io::detail {

// amber_ascii.cpp: Amber ASCII trajectory (TRJ: .mdcrd .crd .trj; title line as info().header,
// box-line detection) and restart (INPCRD: .rst7 .inpcrd; time from line 2, never a box)
std::unique_ptr<core::Trajectory> openAmberTrajectory(const std::filesystem::path& path, std::optional<double> dt,
                                                      std::size_t atoms);
std::unique_ptr<core::Trajectory> openAmberRestart(const std::filesystem::path& path, std::optional<double> dt,
                                                   std::size_t atoms);

// dcd.cpp: CHARMM / NAMD DCD
std::unique_ptr<core::Trajectory> openDcd(const std::filesystem::path& path, std::optional<double> dt,
                                          std::size_t atoms);

// gromacs.cpp: XTC and TRR
std::unique_ptr<core::Trajectory> openXtc(const std::filesystem::path& path, std::optional<double> dt,
                                          std::size_t atoms);
std::unique_ptr<core::Trajectory> openTrr(const std::filesystem::path& path, std::optional<double> dt,
                                          std::size_t atoms);

// pdb_topology.cpp: PDB and GRO as topologies (MDAnalysis PDBParser / GROParser) and as
// coordinates (PDBReader with header / compounds / remarks, multi-model; GROReader)
TopologyFile readPdbTopology(const std::filesystem::path& path, std::optional<double> dt);
TopologyFile readGroTopology(const std::filesystem::path& path, std::optional<double> dt);
std::unique_ptr<core::Trajectory> openPdbCoordinates(const std::filesystem::path& path, std::optional<double> dt,
                                                     std::size_t atoms);
std::unique_ptr<core::Trajectory> openGroCoordinates(const std::filesystem::path& path, std::optional<double> dt,
                                                     std::size_t atoms);

}  // namespace mdbridge::io::detail
