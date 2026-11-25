#include "mdbridge/io/trajectory.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/netcdf.hpp"
#include "mdbridge/io/warnings.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace {

using core::ValueError;
using detail::ReaderError;

enum class Kind { Ncdf, Trj, Inpcrd, Dcd, Xtc, Trr, Pdb, Gro };

struct ReaderClass {
    Kind kind;
    const char* className;  // MDAnalysis' reader class, as ChainReader's messages name it
    const char* form;       // str(Reader.format), as the Universe's atom-count message shows it
};

std::optional<ReaderClass> readerFor(std::string_view format) {
    static const std::array<std::pair<std::string_view, ReaderClass>, 13> readers = {{
        {"NCDF", {Kind::Ncdf, "MDAnalysis.coordinates.TRJ.NCDFReader", "['NCDF', 'NC']"}},
        {"NC", {Kind::Ncdf, "MDAnalysis.coordinates.TRJ.NCDFReader", "['NCDF', 'NC']"}},
        {"TRJ", {Kind::Trj, "MDAnalysis.coordinates.TRJ.TRJReader", "['TRJ', 'MDCRD', 'CRDBOX']"}},
        {"MDCRD",
         {Kind::Trj, "MDAnalysis.coordinates.TRJ.TRJReader", "['TRJ', 'MDCRD', 'CRDBOX']"}},
        {"CRDBOX",
         {Kind::Trj, "MDAnalysis.coordinates.TRJ.TRJReader", "['TRJ', 'MDCRD', 'CRDBOX']"}},
        {"INPCRD",
         {Kind::Inpcrd, "MDAnalysis.coordinates.INPCRD.INPReader", "['INPCRD', 'RESTRT']"}},
        {"RESTRT",
         {Kind::Inpcrd, "MDAnalysis.coordinates.INPCRD.INPReader", "['INPCRD', 'RESTRT']"}},
        {"DCD", {Kind::Dcd, "MDAnalysis.coordinates.DCD.DCDReader", "DCD"}},
        {"XTC", {Kind::Xtc, "MDAnalysis.coordinates.XTC.XTCReader", "XTC"}},
        {"TRR", {Kind::Trr, "MDAnalysis.coordinates.TRR.TRRReader", "TRR"}},
        {"PDB", {Kind::Pdb, "MDAnalysis.coordinates.PDB.PDBReader", "['PDB', 'ENT']"}},
        {"ENT", {Kind::Pdb, "MDAnalysis.coordinates.PDB.PDBReader", "['PDB', 'ENT']"}},
        {"GRO", {Kind::Gro, "MDAnalysis.coordinates.GRO.GROReader", "GRO"}},
    }};
    for (const auto& [name, reader] : readers) {
        if (name == format) {
            return reader;
        }
    }
    return std::nullopt;
}

// MDAnalysis' reader registry in its order (its "Unknown format" message lists them).
constexpr std::array<std::string_view, 51> kRegistry = {
    "CHAIN",    "CHEMFILES", "CRD",    "DCD",        "CONFIG",
    "HISTORY",  "DMS",       "GMS",    "GRO",        "IMD",
    "INPCRD",   "RESTRT",    "CIF",    "CIF.GZ",     "MMCIF",
    "MMCIF.GZ", "LAMMPS",    "DATA",   "LAMMPSDUMP", "MOL2",
    "PDB",      "ENT",       "XPDB",   "PDBQT",      "PQR",
    "TPR",      "TRC",       "TRJ",    "MDCRD",      "CRDBOX",
    "NCDF",     "NC",        "TRR",    "H5MD",       "TRZ",
    "XTC",      "XYZ",       "TXYZ",   "ARC",        "MEMORY",
    "MMTF",     "GSD",       "COOR",   "NAMDBIN",    "IN",
    "FHIAIMS",  "TNG",       "PARMED", "RDKIT",      "OPENMMSIMULATION",
    "OPENMMAPP"};

std::string upper(std::string_view text) {
    std::string out(text);
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return out;
}

std::string lower(std::string_view text) {
    std::string out(text);
    for (char& c : out) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return out;
}

std::string unknownFormatText(const std::string& format, const std::string& path) {
    std::string keys;
    for (const std::string_view key : kRegistry) {
        keys += (keys.empty() ? "'" : ", '") + std::string(key) + "'";
    }
    return "Unknown coordinate trajectory format '" + format + "' for '" + path +
           "'. The FORMATs \n           dict_keys([" + keys +
           "])\n           are implemented in MDAnalysis.\n"
           "           See "
           "https://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1\n"
           "           Use the format keyword to explicitly set the format: "
           "'Universe(...,format=FORMAT)'\n"
           "           For missing formats, raise an issue at "
           "https://github.com/MDAnalysis/mdanalysis/issues";
}

// MDAnalysis' util.guess_format: the extension after a .gz / .bz2 suffix, upper-cased.
std::string guessFormat(std::string_view path) {
    const std::string base = core::pathBasename(path);
    auto [root, ext] = core::posixpath::splitext(base);
    const std::string compressed = lower(ext);
    if (compressed == ".gz" || compressed == ".bz2") {
        ext = core::posixpath::splitext(root).second;
    }
    return upper(ext.empty() ? ext : ext.substr(1));
}

// The reader of `format` ("" guesses), or the error MDAnalysis gives; `chained`: the form of the
// message inside a ChainReader (no "Cannot find an appropriate coordinate reader" preamble).
ReaderClass resolve(const std::filesystem::path& path, std::string_view format, bool chained) {
    const std::string text = detail::pathText(path);
    const std::string name = format.empty() ? guessFormat(text) : upper(format);
    if (const auto reader = readerFor(name)) {
        return *reader;
    }
    if (std::find(kRegistry.begin(), kRegistry.end(), name) != kRegistry.end()) {
        throw ValueError("MDBridge cannot read " + name + " coordinate files ('" + text +
                         "'); it reads NCDF, TRJ (mdcrd), INPCRD (rst7), DCD, XTC, TRR, PDB "
                         "and GRO coordinates");
    }
    if (chained) {
        throw ValueError(unknownFormatText(name, text));
    }
    throw ReaderError("TypeError", "Cannot find an appropriate coordinate reader for file '" +
                                       text + "'.\n           " + unknownFormatText(name, text));
}

std::size_t atomCountOf(const core::Trajectory& trajectory) {
    return trajectory.frameCount() > 0 ? trajectory.current().positions.size() / 3 : 0;
}

// The reader alone: `atoms` reaches the readers that take n_atoms (NCDF checks it, TRJ needs it).
std::unique_ptr<core::Trajectory> openWith(const ReaderClass& reader,
                                           const std::filesystem::path& path,
                                           std::optional<double> dt, std::size_t atoms) {
    switch (reader.kind) {
    case Kind::Ncdf:
        return openNetcdfTrajectory(path, dt, atoms);
    case Kind::Trj:
        return detail::openAmberTrajectory(path, dt, atoms);
    case Kind::Inpcrd:
        return detail::openAmberRestart(path, dt, atoms);
    case Kind::Dcd:
        return detail::openDcd(path, dt, atoms);
    case Kind::Xtc:
        return detail::openXtc(path, dt, atoms);
    case Kind::Trr:
        return detail::openTrr(path, dt, atoms);
    case Kind::Pdb:
        return detail::openPdbCoordinates(path, dt, atoms);
    case Kind::Gro:
        return detail::openGroCoordinates(path, dt, atoms);
    }
    throw std::logic_error("unknown reader");
}

void checkAtoms(const core::Trajectory& trajectory, std::size_t atoms, const std::string& form,
                const std::string& fileText) {
    const std::size_t found = atomCountOf(trajectory);
    if (atoms != 0 && found != atoms) {
        throw ValueError("The topology and " + form +
                         " trajectory files don't have the same number of atoms!\n"
                         "Topology number of atoms " +
                         std::to_string(atoms) + "\nTrajectory: " + fileText + " Number of atoms " +
                         std::to_string(found));
    }
}

class ChainTrajectory : public core::Trajectory {
public:
    explicit ChainTrajectory(std::vector<std::unique_ptr<core::Trajectory>> parts)
        : m_parts(std::move(parts)) {
        m_starts.push_back(0);
        for (const auto& part : m_parts) {
            m_starts.push_back(m_starts.back() + part->frameCount());
        }
        if (frameCount() > 0) {
            seek(0);
        }
    }

    std::int64_t frameCount() const override { return m_starts.back(); }
    double dt() const override { return m_parts[m_active]->dt(); }

    const core::Frame& seek(std::int64_t index) override {
        const std::int64_t n = frameCount();
        if (index < 0) {
            index += n;
        }
        if (index < 0 || index >= n) {
            throw std::out_of_range("Frame index " + std::to_string(index) +
                                    " is out of range for a trajectory of " + std::to_string(n) +
                                    " frames");
        }
        const auto next = std::upper_bound(m_starts.begin(), m_starts.end(), index);
        const auto part = static_cast<std::size_t>(next - m_starts.begin() - 1);
        const std::int64_t local = index - m_starts[part];
        core::Trajectory& reader = *m_parts[part];
        m_frame = reader.seek(local);
        m_active = part;
        m_frame.index = index;
        const auto* file = dynamic_cast<const detail::FileTrajectory*>(&reader);
        if (file == nullptr || !file->storedTime(local)) {
            m_frame.time = static_cast<double>(index) * reader.dt();
        }
        return m_frame;
    }

    const core::Frame& current() const override {
        if (frameCount() == 0) {
            throw std::out_of_range("the trajectory has no frames");
        }
        return m_frame;
    }

private:
    std::vector<std::unique_ptr<core::Trajectory>> m_parts;
    std::vector<std::int64_t> m_starts;
    std::size_t m_active = 0;
    core::Frame m_frame;
};

// str() of a numpy array of these strings: repr()s separated by spaces, wrapped at 75 columns.
std::string numpyStrList(const std::vector<std::string>& items) {
    std::string out;
    std::string line = "[";
    for (std::size_t i = 0; i < items.size(); ++i) {
        const bool last = i + 1 == items.size();
        const std::string word = core::pythonStringRepr(items[i]);
        if (line.size() + word.size() > (last ? 74u : 75u) && line.size() > 1) {
            while (!line.empty() && line.back() == ' ') {
                line.pop_back();
            }
            out += line + "\n";
            line = " ";
        }
        line += word;
        if (!last) {
            line += ' ';
        }
    }
    return out + line + "]";
}

std::string pythonListRepr(const std::vector<std::string>& items) {
    std::string out = "[";
    for (std::size_t i = 0; i < items.size(); ++i) {
        out += (i ? ", " : "") + core::pythonStringRepr(items[i]);
    }
    return out + "]";
}

}  // namespace

namespace detail {

namespace {

// u . v as a fused chain: fma(u2, v2, fma(u1, v1, u0 * v0))
double dot3(const std::array<double, 3>& u, const std::array<double, 3>& v) {
    return std::fma(u[2], v[2], std::fma(u[1], v[1], u[0] * v[0]));
}

// The angle (degrees) between two vectors of the given lengths; NaN when undefined.
double angleBetween(const std::array<double, 3>& u, const std::array<double, 3>& v, double lu,
                    double lv) {
    constexpr double kDegreesPerRadian = 57.29577951308232;  // 180 / pi as one double
    return std::acos(dot3(u, v) / (lu * lv)) * kDegreesPerRadian;
}

}  // namespace

std::array<float, 6> triclinicBox(const std::array<double, 3>& x, const std::array<double, 3>& y,
                                  const std::array<double, 3>& z) {
    const double a = std::sqrt(dot3(x, x));
    const double b = std::sqrt(dot3(y, y));
    const double c = std::sqrt(dot3(z, z));
    const std::array<double, 3> angles = {angleBetween(y, z, b, c), angleBetween(x, z, a, c),
                                          angleBetween(x, y, a, b)};
    const std::array<double, 6> cell = {a, b, c, angles[0], angles[1], angles[2]};
    std::array<float, 6> out{};
    for (std::size_t i = 0; i < 6; ++i) {
        out[i] = static_cast<float>(cell[i]);
        // NaN fails the comparison too
        if (!(out[i] > 0.0f)) {
            return {};
        }
    }
    for (const double angle : angles) {
        if (!(angle < 180.0)) {
            return {};
        }
    }
    return out;
}

const core::Frame& FileTrajectory::seek(std::int64_t index) {
    const std::int64_t n = frameCount();
    if (index < 0) {
        index += n;
    }
    if (index < 0 || index >= n) {
        throw std::out_of_range("Frame index " + std::to_string(index) +
                                " is out of range for a trajectory of " + std::to_string(n) +
                                " frames");
    }
    core::Frame frame;
    readFrame(index, frame);
    frame.index = index;
    m_frame = std::move(frame);
    return m_frame;
}

}  // namespace detail

std::string coordinateFormatForExtension(std::string_view path) {
    static const std::array<std::pair<std::string_view, std::string_view>, 9> table = {{
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
    const std::string ext = lower(core::posixpath::splitext(core::pathBasename(path)).second);
    for (const auto& [suffix, format] : table) {
        if (suffix == ext) {
            return std::string(format);
        }
    }
    return guessFormat(path);
}

namespace {

// inputs.dt_ps given for a file that stores its own times (NetCDF with a time variable): the
// override replaces them (frame x dt_ps, MDATK 0.1.0), so say so when the file's step is another one.
void warnTimeStep(core::Trajectory& trajectory, const std::string& fileText, double dt) {
    if (trajectory.frameCount() < 2) {
        return;
    }
    const double t0 = trajectory.seek(0).time;
    const double t1 = trajectory.seek(1).time;
    trajectory.seek(0);
    const double step = t1 - t0;
    if (std::abs(step - dt) > 1e-4 * std::max(std::abs(dt), std::abs(step))) {
        readerWarning("inputs.dt_ps (" + core::pythonFloatRepr(dt) + " ps) replaces the time step in " + fileText +
                      " (" + core::pythonFloatRepr(step) + " ps)");
    }
}

// The time-step override (MDATK 0.1.0, io/loader.py _FrameTimes): the time of frame k is
// k x dt, whatever the file's time stamps. MDAnalysis' own dt only fills in for trajectories
// without time stamps; a NetCDF file keeps the times stored in it, and a second file of a chain
// starts again at its own first time.
class FrameTimesTrajectory final : public core::Trajectory {
public:
    FrameTimesTrajectory(std::unique_ptr<core::Trajectory> inner, double dt) : m_inner(std::move(inner)), m_dt(dt) {}

    core::TrajectoryInfo info() const override { return m_inner->info(); }
    std::int64_t frameCount() const override { return m_inner->frameCount(); }
    double dt() const override { return m_inner->dt(); }
    const core::Frame& seek(std::int64_t index) override {
        m_valid = false;
        retime(m_inner->seek(index));
        return m_frame;
    }
    const core::Frame& current() const override {
        if (!m_valid) {  // the frame the reader opened on (only seek() moves the inner reader)
            retime(m_inner->current());
        }
        return m_frame;
    }

private:
    void retime(const core::Frame& frame) const {
        m_frame = frame;
        m_frame.time = static_cast<double>(frame.index) * m_dt;
        m_valid = true;
    }

    std::unique_ptr<core::Trajectory> m_inner;
    double m_dt;
    mutable core::Frame m_frame;
    mutable bool m_valid = false;
};

}  // namespace

std::unique_ptr<core::Trajectory> withFrameTimes(std::unique_ptr<core::Trajectory> trajectory,
                                                 std::optional<double> dt) {
    if (!dt || *dt == 0.0) {
        return trajectory;
    }
    return std::make_unique<FrameTimesTrajectory>(std::move(trajectory), *dt);
}

std::unique_ptr<core::Trajectory> openTrajectory(const std::filesystem::path& path,
                                                 std::string_view format, std::optional<double> dt,
                                                 std::size_t atoms) {
    const ReaderClass reader = resolve(path, format, false);
    auto trajectory = openWith(reader, path, dt, atoms);
    if (reader.kind != Kind::Ncdf && reader.kind != Kind::Trj) {
        checkAtoms(*trajectory, atoms, reader.form, detail::pathText(path));
    }
    if (dt) {
        warnTimeStep(*trajectory, detail::pathText(path), *dt);
    }
    return withFrameTimes(std::move(trajectory), dt);
}

std::unique_ptr<core::Trajectory> openTrajectories(const std::vector<std::filesystem::path>& paths,
                                                   std::string_view format,
                                                   std::optional<double> dt, std::size_t atoms) {
    return openTrajectories(paths, std::vector<std::string>(paths.size(), std::string(format)), dt, atoms);
}

std::unique_ptr<core::Trajectory> openTrajectories(const std::vector<std::filesystem::path>& paths,
                                                   const std::vector<std::string>& formats,
                                                   std::optional<double> dt, std::size_t atoms) {
    if (formats.size() != paths.size()) {
        throw std::invalid_argument("openTrajectories: one format per path");
    }
    if (paths.size() == 1) {
        return openTrajectory(paths.front(), formats.front(), dt, atoms);
    }
    std::vector<std::unique_ptr<core::Trajectory>> parts;
    std::vector<std::string> texts;
    for (std::size_t p = 0; p < paths.size(); ++p) {
        const std::filesystem::path& path = paths[p];
        const ReaderClass reader = resolve(path, formats[p], true);
        texts.push_back(detail::pathText(path));
        try {
            // NCDF checks the atom count itself; in a chain the chain's checks name the files
            parts.push_back(openWith(reader, path, dt, reader.kind == Kind::Ncdf ? 0 : atoms));
        } catch (const detail::FormatError&) {
            throw;  // names the file and place already
        } catch (const ValueError& error) {
            if (std::string_view(error.pythonName()) != "ValueError") {
                throw;
            }
            throw ReaderError("TypeError", "Unable to read " + texts.back() + " with <class '" +
                                               reader.className + "'>.");
        }
        if (dt) {
            warnTimeStep(*parts.back(), texts.back(), *dt);
        }
    }
    const std::size_t first = atomCountOf(*parts.front());
    std::vector<std::string> wrong;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (atomCountOf(*parts[i]) != first) {
            wrong.push_back(texts[i]);
        }
    }
    if (!wrong.empty()) {
        throw ValueError("The following trajectories do not have the correct n_atoms  (" +
                         std::to_string(first) + "):\n" + numpyStrList(wrong));
    }
    auto chain = chainTrajectories(std::move(parts));
    checkAtoms(*chain, atoms, "CHAIN", pythonListRepr(texts));
    return withFrameTimes(std::move(chain), dt);
}

std::unique_ptr<core::Trajectory> chainTrajectories(
    std::vector<std::unique_ptr<core::Trajectory>> parts) {
    if (parts.empty()) {
        throw ValueError("no trajectories to chain");
    }
    return std::make_unique<ChainTrajectory>(std::move(parts));
}

}  // namespace mdbridge::io
