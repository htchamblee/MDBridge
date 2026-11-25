// Amber ASCII coordinates: mdcrd trajectories (TRJReader of MDAnalysis 2.11.0: a title line, 10F8.3
// lines, a box line after each frame when the line after the first frame holds three F8.3 fields)
// and rst7 / inpcrd restarts (INPReader: one frame, time from the second line, never a box).
// MDBridge checks the layout against the atom count at open and names the file, the line and the
// columns in its errors (docs/DIFFERENCES.md, "Amber ASCII coordinates").

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/warnings.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

using core::ValueError;

// A line without its line end.
U32 body(const U32& line) {
    U32 out = line;
    while (!out.empty() && (out.back() == U'\n' || out.back() == U'\r')) {
        out.pop_back();
    }
    return out;
}

bool isBlank(const U32& line) { return rstrip(line).empty(); }

// The value of the field at columns [first, first + width) of `line`, or the RM error naming it.
double fieldValue(const std::string& file, const U32& line, std::size_t lineNumber,
                  std::size_t first, std::size_t width) {
    const U32 text = slice(body(line), first, first + width);
    if (const auto value = toFloat(text)) {
        return *value;
    }
    std::string message = file + " line " + std::to_string(lineNumber) + ", columns " +
                          std::to_string(first + 1) + "-" + std::to_string(first + width) + ": " +
                          core::pythonStringRepr(toUtf8(text)) + " is not a number";
    const U32 stripped = strip(text);
    if (!stripped.empty() && stripped.find_first_not_of(U'*') == U32::npos) {
        message += " (Fortran overflow: the value does not fit the field)";
    }
    throw ReaderError("OSError", message);
}

// ---------------------------------------------------------------- mdcrd (TRJReader)

constexpr std::size_t kFieldWidth = 8;
constexpr std::size_t kFieldsPerLine = 10;

// The 8-column fields a line holds: its length without trailing blanks, in fields, rounded up.
std::size_t fieldCount(const U32& line) {
    return (rstrip(line).size() + kFieldWidth - 1) / kFieldWidth;
}

class TrjTrajectory : public FileTrajectory {
public:
    TrjTrajectory(const std::filesystem::path& path, std::optional<double> dt, std::size_t atoms)
        : m_file(pathText(path)), m_atoms(atoms) {
        if (atoms == 0) {
            throw ValueError("AMBER TRJ reader REQUIRES the n_atoms keyword");
        }
        const std::string bytes = readBytes(path);
        if (bytes.empty()) {
            fail("file is empty");
        }
        m_lines = splitPythonLines(decodeUtf8Replace(bytes));
        if (rstrip(m_lines.front()).size() > 80) {
            throw ReaderError("OSError",
                              "Header of AMBER formatted trajectory has more than 80 chars. This "
                              "is probably not a AMBER trajectory.");
        }
        while (m_lines.size() > 1 && isBlank(m_lines.back())) {  // RD-A04
            m_lines.pop_back();
        }
        if (m_lines.size() < 2) {
            fail("no frames after the title");
        }
        m_lineCount = (3 * atoms + kFieldsPerLine - 1) / kFieldsPerLine;
        m_lastFields = (3 * atoms - 1) % kFieldsPerLine + 1;
        checkFrame(1);
        detectBox();
        const std::size_t perFrame = m_lineCount + (m_periodic ? 1 : 0);
        std::size_t first = 1;
        while (first + perFrame <= m_lines.size()) {
            checkFrame(first);
            if (m_periodic) {
                checkBoxLine(first + m_lineCount);
            }
            first += perFrame;
        }
        if (first < m_lines.size()) {  // RD-A04
            readerWarning(m_file + ": the last frame is incomplete (lines " + std::to_string(first + 1) + "-" +
                          std::to_string(m_lines.size()) + ") and is ignored");
        }
        m_frameCount = static_cast<std::int64_t>((m_lines.size() - 1) / perFrame);
        if (m_frameCount == 0) {
            fail("the file ends inside the first frame (" + std::to_string(atoms) + " atoms need " +
                 std::to_string(m_lineCount) + " coordinate lines)");
        }
        m_dt = dt.value_or(1.0);
        seek(0);
    }

    core::TrajectoryInfo info() const override {
        core::TrajectoryInfo out;
        out.header = toUtf8(m_lines.front());
        return out;
    }

    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        const std::size_t perFrame = m_lineCount + (m_periodic ? 1 : 0);
        const std::size_t first = 1 + static_cast<std::size_t>(index) * perFrame;
        frame.positions.clear();
        for (std::size_t i = 0; i < m_lineCount; ++i) {
            const std::size_t fields = i + 1 == m_lineCount ? m_lastFields : kFieldsPerLine;
            for (std::size_t k = 0; k < fields; ++k) {
                frame.positions.push_back(static_cast<float>(
                    fieldValue(m_file, m_lines[first + i], first + i + 1, k * kFieldWidth,
                               kFieldWidth)));
            }
        }
        if (m_periodic) {
            const std::size_t at = first + m_lineCount;
            std::array<float, 6> box = {0.0f, 0.0f, 0.0f, 90.0f, 90.0f, 90.0f};
            for (std::size_t k = 0; k < 3; ++k) {
                box[k] = static_cast<float>(
                    fieldValue(m_file, m_lines[at], at + 1, k * kFieldWidth, kFieldWidth));
            }
            if (box[0] != 0.0f || box[1] != 0.0f || box[2] != 0.0f) {
                frame.box = box;
            }
        }
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    [[noreturn]] void fail(const std::string& message) const {
        throw ReaderError("OSError", m_file + ": " + message);
    }

    [[noreturn]] void failAt(std::size_t index, const std::string& message) const {
        throw ReaderError("OSError", m_file + " line " + std::to_string(index + 1) + ": " + message);
    }

    // The coordinate lines of the frame starting at line index `first` (RD-A02): ten fields each,
    // 3N mod 10 (or ten) on the last. Lines past the end of the file are left for the caller.
    void checkFrame(std::size_t first) const {
        for (std::size_t i = 0; i < m_lineCount && first + i < m_lines.size(); ++i) {
            const std::size_t expected = i + 1 == m_lineCount ? m_lastFields : kFieldsPerLine;
            const std::size_t found = fieldCount(m_lines[first + i]);
            if (found != expected) {
                failAt(first + i, "expected " + std::to_string(expected) +
                                      (expected == 1 ? " coordinate for " : " coordinates for ") +
                                      std::to_string(m_atoms) + " atoms, found " +
                                      std::to_string(found) + " (does the topology match?)");
            }
        }
    }

    void checkBoxLine(std::size_t index) const {
        const std::size_t found = fieldCount(m_lines[index]);
        if (found != 3) {
            failAt(index, "expected a box line (3 values), found " + std::to_string(found) +
                              " values (does the topology match?)");
        }
    }

    // Periodic when the line after the first frame is a box line (RD-A01, RD-A03): three fields
    // where the next frame would start with another count. With one atom (a frame is one line of
    // three) every second line is the box when there is an even number of lines and all of those
    // hold three positive values.
    void detectBox() {
        const std::size_t next = 1 + m_lineCount;
        if (next >= m_lines.size()) {
            m_periodic = false;  // one frame, no box line
            return;
        }
        const std::size_t startFields = m_lineCount == 1 ? m_lastFields : kFieldsPerLine;
        const std::size_t found = fieldCount(m_lines[next]);
        if (startFields == 3) {
            const std::size_t data = m_lines.size() - 1;
            m_periodic = data % 2 == 0;
            for (std::size_t i = 2; m_periodic && i < m_lines.size(); i += 2) {
                for (std::size_t k = 0; k < 3; ++k) {
                    const auto value = toFloat(slice(body(m_lines[i]), k * kFieldWidth,
                                                     (k + 1) * kFieldWidth));
                    m_periodic = m_periodic && fieldCount(m_lines[i]) == 3 && value && *value > 0.0;
                }
            }
            return;
        }
        if (found == 3) {
            m_periodic = true;
            return;
        }
        if (found != startFields) {
            failAt(next, "expected a box line (3 values) or the first " +
                             std::to_string(startFields) + " coordinates of the next frame, found " +
                             std::to_string(found) + " values (does the topology match?)");
        }
        m_periodic = false;
    }

    std::string m_file;
    std::size_t m_atoms;
    Lines m_lines;
    std::size_t m_lineCount = 0;   // coordinate lines per frame
    std::size_t m_lastFields = 0;  // fields on a frame's last coordinate line
    bool m_periodic = false;
};

// ---------------------------------------------------------------- inpcrd / rst7 (INPReader)

constexpr std::size_t kRestartFieldWidth = 12;

class RestartTrajectory : public FileTrajectory {
public:
    RestartTrajectory(const std::filesystem::path& path, std::optional<double> dt) {
        const std::string file = pathText(path);
        const auto fail = [&](const std::string& message) {
            throw ReaderError("OSError", file + ": " + message);
        };
        const Lines lines = readTextLines(path);
        if (lines.empty()) {
            fail("file is empty");
        }
        if (lines.size() < 2) {
            fail("no atom count line after the title");
        }
        const std::vector<U32> tokens = split(lines[1]);
        if (tokens.empty()) {
            fail("line 2 has no atom count");
        }
        std::optional<std::int64_t> count;
        try {
            count = core::tryPythonInt(toUtf8(tokens[0]));
        } catch (const std::overflow_error&) {
            fail("line 2: atom count " + core::pythonStringRepr(toUtf8(tokens[0])) +
                 " is too large");
        }
        if (!count || *count < 0) {
            throw ReaderError("OSError", file + " line 2: atom count " +
                                             core::pythonStringRepr(toUtf8(tokens[0])) +
                                             " is not a non-negative integer");
        }
        if (tokens.size() > 1) {
            m_time = core::tryPythonFloat(toUtf8(tokens[1]));
            if (!m_time) {
                throw ReaderError("OSError", file + " line 2: time " +
                                                 core::pythonStringRepr(toUtf8(tokens[1])) +
                                                 " is not a number");
            }
        }
        const auto atoms = static_cast<std::size_t>(*count);
        for (std::size_t atom = 0; atom < atoms; atom += 2) {
            const std::size_t lineIndex = 2 + atom / 2;
            const std::size_t fields = atom + 1 < atoms ? 6 : 3;
            if (lineIndex >= lines.size()) {
                fail("the file ends at line " + std::to_string(lines.size()) + "; " +
                     std::to_string(atoms) + " atoms need " +
                     std::to_string(2 + (atoms + 1) / 2) + " lines");
            }
            const U32& line = lines[lineIndex];
            const std::size_t found = (rstrip(line).size() + kRestartFieldWidth - 1) / kRestartFieldWidth;
            if (found < fields) {
                throw ReaderError("OSError", file + " line " + std::to_string(lineIndex + 1) +
                                                 ": expected " + std::to_string(fields) +
                                                 " coordinates, found " + std::to_string(found) +
                                                 " (the file has " + std::to_string(atoms) +
                                                 " atoms)");
            }
            for (std::size_t k = 0; k < fields; ++k) {
                m_positions.push_back(static_cast<float>(fieldValue(
                    file, line, lineIndex + 1, k * kRestartFieldWidth, kRestartFieldWidth)));
            }
        }
        m_frameCount = 1;
        m_dt = dt.value_or(1.0);
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return m_time.has_value(); }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        frame.positions = m_positions;
        frame.time = m_time.value_or(static_cast<double>(index) * m_dt);
    }

private:
    std::vector<float> m_positions;
    std::optional<double> m_time;
};

}  // namespace

std::unique_ptr<core::Trajectory> openAmberTrajectory(const std::filesystem::path& path,
                                                      std::optional<double> dt, std::size_t atoms) {
    return std::make_unique<TrjTrajectory>(path, dt, atoms);
}

std::unique_ptr<core::Trajectory> openAmberRestart(const std::filesystem::path& path,
                                                   std::optional<double> dt, std::size_t) {
    return std::make_unique<RestartTrajectory>(path, dt);
}

}  // namespace mdbridge::io::detail
