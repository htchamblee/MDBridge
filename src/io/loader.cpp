#include "mdbridge/io/loader.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <system_error>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/parm7.hpp"
#include "mdbridge/io/trajectory.hpp"
#include "mdbridge/io/warnings.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace {

using core::ValueError;
using detail::ReaderError;

// FileNotFoundError: what() is the message alone
class NotFoundError : public std::filesystem::filesystem_error {
public:
    NotFoundError(const std::string& message, const std::filesystem::path& path)
        : std::filesystem::filesystem_error(
              message, path, std::make_error_code(std::errc::no_such_file_or_directory)),
          m_message(message) {}
    const char* what() const noexcept override { return m_message.c_str(); }

private:
    std::string m_message;
};

std::string upper(std::string_view text) {
    std::string out(text);
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return out;
}

std::string lowerExtension(std::string_view path) {
    std::string ext = core::posixpath::splitext(core::pathBasename(path)).second;
    for (char& c : ext) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return ext;
}

// MDAnalysis' util.guess_format for a topology: the extension (after .gz / .bz2), upper-cased
std::string guessTopologyFormat(std::string_view path) {
    const std::string base = core::pathBasename(path);
    auto [root, ext] = core::posixpath::splitext(base);
    const std::string compressed = upper(ext);
    if (compressed == ".GZ" || compressed == ".BZ2") {
        ext = core::posixpath::splitext(root).second;
    }
    return upper(ext.empty() ? ext : ext.substr(1));
}

// MDAnalysis' topology parser registry in its order (its error lists them)
constexpr std::array<std::string_view, 37> kTopologyFormats = {
    "CRD",      "CONFIG", "HISTORY",  "DMS",        "PDB",     "ENT",
    "XPDB",     "IN",     "FHIAIMS",  "GMS",        "GRO",     "GSD",
    "XML",      "ITP",    "DATA",     "LAMMPSDUMP", "MINIMAL", "CIF",
    "CIF.GZ",   "MMCIF",  "MMCIF.GZ", "MMTF",       "MOL2",    "PDBQT",
    "PQR",      "PSF",    "TOP",      "PRMTOP",     "PARM7",   "TPR",
    "TXYZ",     "ARC",    "XYZ",      "PARMED",     "RDKIT",   "OPENMMTOPOLOGY",
    "OPENMMAPP"};

std::string unknownTopologyText(const std::string& format) {
    std::string keys;
    for (const std::string_view key : kTopologyFormats) {
        keys += (keys.empty() ? "'" : ", '") + std::string(key) + "'";
    }
    return "'" + format +
           "' isn't a valid topology format, nor a coordinate format\n"
           "   from which a topology can be minimally inferred.\n"
           "   You can use 'Universe(topology, ..., topology_format=FORMAT)'\n"
           "   to explicitly specify the format and\n"
           "   override automatic detection. Known FORMATs are:\n"
           "   dict_keys([" +
           keys +
           "])\n"
           "   See "
           "https://docs.mdanalysis.org/documentation_pages/topology/"
           "init.html#supported-topology-formats\n"
           "   For missing formats, raise an issue at \n"
           "   https://github.com/MDAnalysis/mdanalysis/issues";
}

// MDATK's extension table (mdatk/io/loader.py), applied per file; "" lets the reader guess
std::string mdatkFormat(std::string_view path) {
    static const std::array<std::pair<std::string_view, std::string_view>, 9> kFormats = {{
        {".nc", "NCDF"},
        {".ncdf", "NCDF"},
        {".netcdf", "NCDF"},
        {".mdcrd", "TRJ"},
        {".crd", "TRJ"},
        {".trj", "TRJ"},
        {".rst7", "INPCRD"},
        {".inpcrd", "INPCRD"},
        {".ncrst", "NCRST"},
    }};
    const std::string ext = lowerExtension(path);
    for (const auto& [suffix, name] : kFormats) {
        if (suffix == ext) {
            return std::string(name);
        }
    }
    return {};
}

// A .top file read as parm7 (MDATK's mapping) that is not one: GROMACS topologies share the
// extension.
void checkTopExtension(const std::filesystem::path& path, const std::string& text) {
    if (lowerExtension(text) != ".top") {
        return;
    }
    std::ifstream in(path, std::ios::binary);
    std::string first;
    std::getline(in, first);
    if (!first.starts_with("%VE")) {
        throw detail::FormatError(text + ": not an Amber parm7 topology (GROMACS .top files are not supported)");
    }
}

std::string noCoordinates(const std::string& topology) {
    return topology + " has no coordinates; give a trajectory with --traj";
}

}  // namespace

TopologyFile readTopologyFile(const std::filesystem::path& path, std::string_view format,
                              std::optional<double> dt) {
    const std::string text = detail::pathText(path);
    const std::string name = format.empty() ? guessTopologyFormat(text) : upper(format);
    ReaderWarningScope scope;
    TopologyFile out;
    if (name == "PARM7" || name == "PRMTOP" || name == "TOP") {
        checkTopExtension(path, text);
        out.topology = readParm7Topology(path);
    } else if (name == "PDB" || name == "ENT") {
        out = detail::readPdbTopology(path, dt);
    } else if (name == "GRO") {
        out = detail::readGroTopology(path, dt);
    }
    if (out.topology) {
        out.warnings = scope.take();
        return out;
    }
    if (std::find(kTopologyFormats.begin(), kTopologyFormats.end(), name) !=
        kTopologyFormats.end()) {
        throw ValueError("MDBridge cannot read " + name + " topology files ('" + text +
                         "'); it reads parm7 (PARM7 / PRMTOP / TOP), PDB and GRO topologies");
    }
    throw ValueError(unknownTopologyText(name));
}

core::System loadSystem(const std::filesystem::path& topology,
                        const std::vector<std::filesystem::path>& trajectories,
                        std::optional<double> dtPs) {
    if (topology.empty()) {
        throw ValueError("no topology file given");
    }
    if (!std::filesystem::exists(topology)) {
        throw NotFoundError("topology not found: " + detail::pathText(topology), topology);
    }
    for (const auto& path : trajectories) {
        if (!std::filesystem::exists(path)) {
            throw NotFoundError("trajectory not found: " + detail::pathText(path), path);
        }
    }
    const std::string topologyText = detail::pathText(topology);
    const bool parm7 = isParm7(topologyText);
    std::vector<std::string> formats;  // MDATK's table per file (a chain may mix formats)
    for (const auto& path : trajectories) {
        formats.push_back(mdatkFormat(detail::pathText(path)));
    }
    std::optional<double> dt;
    if (dtPs && !(*dtPs >= 0.0 && std::isfinite(*dtPs))) {
        throw ValueError("inputs.dt_ps must be a positive number of picoseconds (0: the files' time step), got " +
                         core::pythonFloatRepr(*dtPs));
    }
    if (dtPs && *dtPs != 0.0) {
        dt = *dtPs;
    }
    TopologyFile file = readTopologyFile(topology, parm7 ? "PARM7" : "", dt);
    core::System system;
    system.topology = file.topology;
    system.warnings = std::move(file.warnings);
    ReaderWarningScope scope;
    if (!trajectories.empty()) {
        system.trajectory = openTrajectories(trajectories, formats, dt, file.topology->atomCount());
    } else if (file.coordinates) {
        system.trajectory = std::move(file.coordinates);
    } else {
        throw ValueError(noCoordinates(topologyText));
    }
    for (std::string& warning : scope.take()) {
        if (std::find(system.warnings.begin(), system.warnings.end(), warning) == system.warnings.end()) {
            system.warnings.push_back(std::move(warning));
        }
    }
    return system;
}

std::vector<std::int64_t> frameIndices(std::int64_t frameCount, std::int64_t start,
                                       std::optional<std::int64_t> stop, std::int64_t step) {
    if (step < 1) {
        throw ValueError("frame step must be >= 1");
    }
    const std::int64_t n = std::max<std::int64_t>(frameCount, 0);
    const auto clamp = [n](std::int64_t index) {
        if (index < 0) {
            index += n;
            return std::max<std::int64_t>(index, 0);
        }
        return std::min(index, n);
    };
    const std::int64_t first = clamp(start);
    const std::int64_t last = stop ? clamp(*stop) : n;
    std::vector<std::int64_t> out;
    for (std::int64_t i = first; i < last; i += step) {
        out.push_back(i);
    }
    if (out.empty()) {
        throw ValueError("frame selection start=" + std::to_string(start) +
                         " stop=" + (stop ? std::to_string(*stop) : std::string("end")) +
                         " step=" + std::to_string(step) + " selects no frames (trajectory has " +
                         std::to_string(frameCount) + ")");
    }
    return out;
}

core::Value systemSummary(const core::System& system) {
    if (!system.trajectory) {
        throw ValueError("the system has no coordinates; give a trajectory with --traj");
    }
    const core::Topology& top = *system.topology;
    std::map<std::string, std::int64_t> counts;
    for (std::size_t r = 0; r < top.residueCount(); ++r) {
        ++counts[top.resname(static_cast<std::int32_t>(r))];
    }
    std::vector<std::pair<std::string, std::int64_t>> order(counts.begin(), counts.end());
    std::stable_sort(order.begin(), order.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });
    core::Table names;
    for (std::size_t k = 0; k < order.size() && k < 12; ++k) {
        names.emplace_back(order[k].first, core::Value(order[k].second));
    }
    const core::Trajectory& trajectory = *system.trajectory;
    const std::int64_t frames = trajectory.frameCount();
    core::Table out;
    out.emplace_back("n_atoms", core::Value(static_cast<std::int64_t>(top.atomCount())));
    out.emplace_back("n_residues", core::Value(static_cast<std::int64_t>(top.residueCount())));
    out.emplace_back("n_segments", core::Value(static_cast<std::int64_t>(top.segmentCount())));
    out.emplace_back("n_frames", core::Value(frames));
    out.emplace_back("dt_ps", frames > 1 ? core::Value(trajectory.dt()) : core::Value());
    out.emplace_back("has_box", core::Value(frames > 0 && trajectory.current().box.has_value()));
    out.emplace_back("top_resnames", core::Value(std::move(names)));
    return core::Value(std::move(out));
}

}  // namespace mdbridge::io
