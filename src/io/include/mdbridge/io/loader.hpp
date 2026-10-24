#pragma once
// Loading an MD system as MDATK does (mdatk/io/loader.py): one topology file, zero or more
// trajectory files (concatenated), the frame selection and the summary `mdatk info` prints.

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "mdbridge/core/topology.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::io {

// The topology MDAnalysis builds from a file of a topology format (parm7, PDB, GRO, ...), and the
// coordinates the file itself holds (PDB, GRO), as MDAnalysis' parsers do. `format`: MDAnalysis'
// topology format name ("PARM7", "PDB", "GRO", ...), "" guesses from the extension.
struct TopologyFile {
    std::shared_ptr<core::Topology> topology;
    std::unique_ptr<core::Trajectory> coordinates;  // nullptr for formats without coordinates (parm7)
};
TopologyFile readTopologyFile(const std::filesystem::path& path, std::string_view format = {},
                              std::optional<double> dt = std::nullopt);

// load_universe(topology, trajectories, dt_ps). Throws core::ValueError("no topology file given"),
// std::filesystem::filesystem_error for "topology not found: <path>" / "trajectory not found:
// <path>" (what() starts with that text), core::ValueError for files the readers reject, and for
// a parm7 topology without trajectories ("No trajectory loaded into Universe"). dtPs: 0 or
// std::nullopt means unset. Format choice: parm7 extensions -> PARM7; one shared trajectory
// extension of the MDATK table -> that reader for every file; otherwise per file.
core::System loadSystem(const std::filesystem::path& topology, const std::vector<std::filesystem::path>& trajectories,
                        std::optional<double> dtPs = std::nullopt);

// frame_indices(n_frames, start, stop, step): Python slice semantics. Throws core::ValueError
// ("frame step must be >= 1", "frame selection start=5 stop=None step=1 selects no frames
// (trajectory has 24)").
std::vector<std::int64_t> frameIndices(std::int64_t frameCount, std::int64_t start = 0,
                                       std::optional<std::int64_t> stop = std::nullopt, std::int64_t step = 1);

// system_summary(u): {n_atoms, n_residues, n_segments, n_frames, dt_ps (null for one frame),
// has_box, top_resnames (12 most frequent residue names; ties by name)}.
core::Value systemSummary(const core::System& system);

}  // namespace mdbridge::io
