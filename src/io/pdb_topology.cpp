// TEMPORARY STUB (replaced by the PDB / GRO readers)
#include "mdbridge/core/errors.hpp"
#include "readers.hpp"

namespace mdbridge::io::detail {

TopologyFile readPdbTopology(const std::filesystem::path&, std::optional<double>) {
    throw core::ValueError("PDB topologies are not implemented yet");
}

TopologyFile readGroTopology(const std::filesystem::path&, std::optional<double>) {
    throw core::ValueError("GRO topologies are not implemented yet");
}

std::unique_ptr<core::Trajectory> openPdbCoordinates(const std::filesystem::path&,
                                                     std::optional<double>, std::size_t) {
    throw core::ValueError("PDB coordinates are not implemented yet");
}

std::unique_ptr<core::Trajectory> openGroCoordinates(const std::filesystem::path&,
                                                     std::optional<double>, std::size_t) {
    throw core::ValueError("GRO coordinates are not implemented yet");
}

}  // namespace mdbridge::io::detail
