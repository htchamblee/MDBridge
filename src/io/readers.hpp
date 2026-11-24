#pragma once
// The per-format readers behind openTrajectory / readTopologyFile (internal to io). Each follows
// MDAnalysis 2.11.0's reader of that format for what MDATK sees: times, dt, boxes, float32 values,
// attributes, and the reader metadata the PDB export writes.

#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <utility>

#include "mdbridge/core/topology.hpp"
#include "mdbridge/io/loader.hpp"

namespace mdbridge::io::detail {

// A reader error whose MDAnalysis counterpart is a Python exception other than ValueError
// (TypeError, OSError, NotImplementedError ...): MDAnalysis' ChainReader rewraps only ValueErrors.
class ReaderError : public core::ValueError {
public:
    ReaderError(std::string pythonName, const std::string& message)
        : core::ValueError(message), m_name(std::move(pythonName)) {}
    const char* pythonName() const noexcept override { return m_name.c_str(); }

private:
    std::string m_name;
};

// A trajectory whose frames are decoded from its file on demand. Times are either stored in the
// file (kept per file when chained) or frame * dt (renumbered by chainTrajectories).
class FileTrajectory : public core::Trajectory {
public:
    std::int64_t frameCount() const override { return m_frameCount; }
    double dt() const override { return m_dt; }
    const core::Frame& seek(std::int64_t index) override;
    const core::Frame& current() const override { return m_frame; }
    // Whether frame `index`'s time comes from the file (MDAnalysis sets ts.time) rather than
    // index * dt.
    virtual bool storedTime(std::int64_t index) const = 0;

protected:
    // Fills positions, box and time of frame `index` (0 <= index < frameCount()).
    virtual void readFrame(std::int64_t index, core::Frame& frame) = 0;

    std::int64_t m_frameCount = 0;
    double m_dt = 1.0;
    core::Frame m_frame;
};

// MDAnalysis' triclinic_box of three box vectors (float64): lengths and angles rounded to float32,
// all zero for an invalid cell (trajectory.cpp).
std::array<float, 6> triclinicBox(const std::array<double, 3>& x, const std::array<double, 3>& y,
                                  const std::array<double, 3>& z);

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
