// PDB and GRO files as MDAnalysis 2.11.0 reads them for MDATK: the topologies of its PDBParser and
// GROParser (with the types and masses the Universe then guesses) and the coordinates of its
// PDBReader (models as frames, CRYST1, HEADER / COMPND / REMARK) and GROReader.

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

#include "pdb_tables.inc"

// Readers: CLEAN ROOM: reimplement from scratchpad/cleanroom/pdb-topology.md and gro.md

}  // namespace

TopologyFile readPdbTopology(const std::filesystem::path& path, std::optional<double> dt) {
    (void)path; (void)dt;
    throw std::logic_error("clean room: PDB topology not implemented");
}

TopologyFile readGroTopology(const std::filesystem::path& path, std::optional<double> dt) {
    (void)path; (void)dt;
    throw std::logic_error("clean room: GRO topology not implemented");
}

std::unique_ptr<core::Trajectory> openPdbCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t atoms) {
    (void)path; (void)dt; (void)atoms;
    throw std::logic_error("clean room: PDB coordinates not implemented");
}

std::unique_ptr<core::Trajectory> openGroCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t atoms) {
    (void)path; (void)dt; (void)atoms;
    throw std::logic_error("clean room: GRO coordinates not implemented");
}

}  // namespace mdbridge::io::detail
