#pragma once
// The PDB export of trajectory frames (mdatk/io/pdb_export.py), whose files are MDAnalysis'
// PDBWriter's (single frame, reindexed, CONECT for non-guessed bonds): HEADER / TITLE / COMPND /
// REMARK / CRYST1 (unitary cell and REMARK 285 without a box), 80-column ATOM / HETATM records
// with MDAnalysis' atom-name alignment, chain-ID validation and wrapping serials and resids,
// CONECT, END. Byte-exact with MDATK v0.0.4.1 (MDAnalysis 2.11.0).

#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "mdbridge/core/topology.hpp"

namespace mdbridge::io {

// The text of one PDB file of `atoms` at `frame`, with the reader metadata `info`. Throws
// core::ValueError("PDB files must have coordinate values between -1000.000 and 9999.999
// Angstroem: file writing was aborted."), for record types other than ATOM / HETATM and formal
// charges beyond 9.
std::string pdbFrameText(const core::AtomGroup& atoms, const core::Frame& frame, const core::TrajectoryInfo& info);

// export_frames: one file <out_dir>/<prefix>_<frame, zero-padded to max(5, digits of the frame
// count)>[_<stamp>].pdb per frame, in the given order; progress((k + 1) / n, "PDB frame <f>").
// Throws SelectionError("no frames selected for PDB export"); a coordinate error deletes that file
// and stops (earlier files stay).
std::vector<std::filesystem::path> exportFrames(core::System& system, const std::vector<std::int64_t>& frames,
                                                const core::AtomGroup& atoms, const std::filesystem::path& outDir,
                                                const std::string& prefix,
                                                const std::function<void(double, const std::string&)>& progress = {},
                                                const std::string& stamp = {});

}  // namespace mdbridge::io
