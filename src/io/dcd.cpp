// CHARMM / NAMD / X-PLOR DCD trajectories as MDAnalysis 2.11.0's DCDReader reads them: Fortran
// records of either byte order, CHARMM unit cells interpreted as cosines, degrees or box vectors,
// fixed atoms, frame count from the file size, times (frame + istart / nsavc) * dt.

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "binary_file.hpp"
#include "mdbridge/core/errors.hpp"
#include "readers.hpp"

namespace mdbridge::io::detail {

namespace {

using core::ValueError;

constexpr double kPi = 3.141592653589793;

[[noreturn]] void headerError(const char* reason) {
    throw ReaderError("OSError", std::string("Reading DCD header failed: ") + reason);
}

constexpr const char* kPrematureEof = "premature EOF found in DCD file";
constexpr const char* kBadFormat = "format of DCD file is wrong";

class DcdTrajectory : public FileTrajectory {
public:
    DcdTrajectory(const std::filesystem::path& path, std::optional<double> dt) : m_file(path) {
        readHeader();
        const std::uint64_t dims = m_fourDims ? 4 : 3;
        const std::uint64_t extra = m_unitCell ? 56 : 0;
        m_firstFrameBytes = (m_atoms + 2) * dims * 4 + extra;
        m_frameBytes = (m_atoms - fixedCount() + 2) * dims * 4 + extra;
        const double frames =
            (static_cast<double>(m_file.size()) - static_cast<double>(m_headerBytes) -
             static_cast<double>(m_firstFrameBytes)) /
                static_cast<double>(m_frameBytes) +
            1.0;
        m_frameCount = static_cast<std::int64_t>(frames);
        if (m_frameCount < 0) {
            m_frameCount = 0;
        }
        if (m_frameCount == 0) {
            throw ReaderError("OSError", "opened empty file. No frames are saved");
        }
        // the first frame holds the fixed atoms' positions
        try {
            m_firstFrame = readRaw(0);
        } catch (const ReaderError&) {
            throw ReaderError("OSError", "DCD is corrupted");
        }
        // delta (AKMA time units) in ps, as MDAnalysis.units.convert computes it
        const double akma = 1.0 / 4.888821e-2;
        const double delta = m_delta * (1.0 / akma);
        m_dt = dt ? *dt : delta * static_cast<double>(m_nsavc);
        if (m_nsavc == 0) {
            throw ReaderError("ZeroDivisionError", "division by zero");
        }
        if (m_fourDims) {
            throw ValueError("could not broadcast input array from shape (" +
                             std::to_string(m_atoms) + ",4) into shape (" +
                             std::to_string(m_atoms) + ",3)");
        }
        m_timeOffset = static_cast<double>(m_istart) / static_cast<double>(m_nsavc);
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return true; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        const Raw raw = index == 0 ? m_firstFrame : readRaw(index);
        frame.positions = raw.positions;
        frame.time = (static_cast<double>(index) + m_timeOffset) * m_dt;
        frame.box = unitCell(raw.cell);
    }

private:
    struct Raw {
        std::vector<float> positions;
        std::array<double, 6> cell{0.0, 90.0, 0.0, 90.0, 90.0, 0.0};  // A, gamma, B, beta, alpha, C
    };

    std::size_t fixedCount() const { return m_free.empty() ? 0 : m_atoms - m_free.size(); }

    // The payload of the Fortran record at m_pos (and moves past it).
    std::string record(std::uint64_t& pos, std::optional<std::uint32_t> expected = std::nullopt) {
        const std::string head = m_file.read(pos, 4);
        if (head.size() != 4) {
            headerError(kPrematureEof);
        }
        const std::uint32_t size = loadU32(head.data(), m_bigEndian);
        if (expected && size != *expected) {
            headerError(kBadFormat);
        }
        const std::string payload = m_file.read(pos + 4, size);
        const std::string tail = m_file.read(pos + 4 + size, 4);
        if (payload.size() != size || tail.size() != 4) {
            headerError(kPrematureEof);
        }
        if (loadU32(tail.data(), m_bigEndian) != size) {
            headerError(kBadFormat);
        }
        pos += 8 + static_cast<std::uint64_t>(size);
        return payload;
    }

    void readHeader() {
        const std::string first = m_file.read(0, 4);
        if (first.size() != 4) {
            headerError(kPrematureEof);
        }
        if (loadU32(first.data(), false) == 84) {
            m_bigEndian = false;
        } else if (loadU32(first.data(), true) == 84) {
            m_bigEndian = true;
        } else {
            headerError(kBadFormat);
        }
        std::uint64_t pos = 0;
        const std::string head = record(pos, 84);
        if (head.compare(0, 4, "CORD") != 0) {
            headerError(kBadFormat);
        }
        const auto control = [&](int i) { return loadI32(head.data() + 4 + 4 * i, m_bigEndian); };
        m_istart = control(1);
        m_nsavc = control(2);
        const std::int32_t fixed = control(8);
        const bool charmm = control(19) != 0;
        if (charmm) {
            // MDAnalysis' libdcd takes the CHARMM delta in the host's (little-endian) byte order
            m_delta = static_cast<double>(loadF32(head.data() + 4 + 36, false));
            m_unitCell = control(10) != 0;
            m_fourDims = control(11) == 1;
        } else {
            m_delta = loadF64(head.data() + 4 + 36, m_bigEndian);
        }
        const std::string title = record(pos);
        if (title.size() < 4) {
            headerError(kBadFormat);
        }
        const std::string atoms = record(pos, 4);
        const std::int32_t n = loadI32(atoms.data(), m_bigEndian);
        if (n <= 0 || fixed < 0 || fixed > n) {
            headerError(kBadFormat);
        }
        m_atoms = static_cast<std::size_t>(n);
        if (fixed > 0) {
            const std::size_t freeCount = m_atoms - static_cast<std::size_t>(fixed);
            const std::string indices = record(pos, static_cast<std::uint32_t>(4 * freeCount));
            for (std::size_t i = 0; i < freeCount; ++i) {
                const std::int32_t atom = loadI32(indices.data() + 4 * i, m_bigEndian);
                if (atom < 1 || static_cast<std::size_t>(atom) > m_atoms) {
                    headerError(kBadFormat);
                }
                m_free.push_back(static_cast<std::size_t>(atom - 1));
            }
        }
        m_headerBytes = pos;
    }

    Raw readRaw(std::int64_t index) {
        std::uint64_t pos = m_headerBytes;
        if (index > 0) {
            pos += m_firstFrameBytes + m_frameBytes * static_cast<std::uint64_t>(index - 1);
        }
        Raw raw;
        if (m_unitCell) {
            const std::string cell = record(pos, 48);
            for (int i = 0; i < 6; ++i) {
                raw.cell[static_cast<std::size_t>(i)] = loadF64(cell.data() + 8 * i, m_bigEndian);
            }
        }
        const bool partial = index > 0 && !m_free.empty();
        const std::size_t count = partial ? m_free.size() : m_atoms;
        raw.positions = partial ? m_firstFrame.positions : std::vector<float>(3 * m_atoms);
        for (std::size_t d = 0; d < (m_fourDims ? 4u : 3u); ++d) {
            const std::string values = record(pos, static_cast<std::uint32_t>(4 * count));
            if (d == 3) {
                break;
            }
            for (std::size_t i = 0; i < count; ++i) {
                const std::size_t atom = partial ? m_free[i] : i;
                raw.positions[3 * atom + d] = loadF32(values.data() + 4 * i, m_bigEndian);
            }
        }
        return raw;
    }

    // DCDReader._frame_to_ts: [A, gamma, B, beta, alpha, C] as lengths and angles.
    static std::optional<std::array<float, 6>> unitCell(const std::array<double, 6>& cell) {
        std::array<double, 6> uc = {cell[0], cell[2], cell[5], cell[4], cell[3], cell[1]};
        const auto inUnit = [](double x) { return -1.0 <= x && x <= 1.0; };
        std::array<float, 6> box{};
        if (inUnit(uc[3]) && inUnit(uc[4]) && inUnit(uc[5])) {
            for (std::size_t i = 3; i < 6; ++i) {
                uc[i] = 90.0 - std::asin(uc[i]) * 90.0 / (kPi / 2);
            }
            for (std::size_t i = 0; i < 6; ++i) {
                box[i] = static_cast<float>(uc[i]);
            }
        } else if (std::any_of(uc.begin(), uc.end(), [](double x) { return x < 0.0; }) ||
                   uc[3] > 180.0 || uc[4] > 180.0 || uc[5] > 180.0) {
            const std::array<double, 3> e1 = {cell[0], cell[1], cell[3]};
            const std::array<double, 3> e2 = {cell[1], cell[2], cell[4]};
            const std::array<double, 3> e3 = {cell[3], cell[4], cell[5]};
            box = triclinicBox(e1, e2, e3);
        } else {
            for (std::size_t i = 0; i < 6; ++i) {
                box[i] = static_cast<float>(uc[i]);
            }
        }
        if (box[0] == 0.0f && box[1] == 0.0f && box[2] == 0.0f) {
            return std::nullopt;
        }
        return box;
    }

    BinaryFile m_file;
    bool m_bigEndian = false;
    std::int32_t m_istart = 0;
    std::int32_t m_nsavc = 0;
    double m_delta = 0.0;  // CHARMM: float32; X-PLOR: float64
    bool m_unitCell = false;
    bool m_fourDims = false;
    std::size_t m_atoms = 0;
    std::vector<std::size_t> m_free;
    std::uint64_t m_headerBytes = 0;
    std::uint64_t m_firstFrameBytes = 0;
    std::uint64_t m_frameBytes = 0;
    double m_timeOffset = 0.0;
    Raw m_firstFrame;
};

}  // namespace

std::unique_ptr<core::Trajectory> openDcd(const std::filesystem::path& path,
                                          std::optional<double> dt, std::size_t) {
    return std::make_unique<DcdTrajectory>(path, dt);
}

}  // namespace mdbridge::io::detail
