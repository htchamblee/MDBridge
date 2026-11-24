// Amber ASCII coordinates as MDAnalysis 2.11.0 reads them: TRJReader (mdcrd: 10F8.3 lines, a
// box line after each frame when the line after the first frame has exactly three F8.3 fields)
// and INPReader (rst7 / inpcrd: one frame, time from the second line, never a box).

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

using core::ValueError;

// ---------------------------------------------------------------- mdcrd (TRJReader)

constexpr std::size_t kFieldWidth = 8;
constexpr std::size_t kFieldsPerLine = 10;

std::string reshapeMessage(std::size_t count, std::size_t atoms) {
    return "cannot reshape array of size " + std::to_string(count) + " into shape (" +
           std::to_string(atoms) + ",3)";
}

// Field `k` (8 characters) of a line: its text and value (nullopt: not a number).
struct Field {
    U32 text;
    std::optional<double> value;
};

Field field(const U32& line, std::size_t k) {
    Field f;
    f.text = slice(line, k * kFieldWidth, (k + 1) * kFieldWidth);
    f.value = toFloat(f.text);
    return f;
}

[[noreturn]] void fieldError(std::size_t k, const U32& text) {
    throw ValueError("FixedcolumnEntry(" + std::to_string(k * kFieldWidth) + "," +
                     std::to_string((k + 1) * kFieldWidth) + ",'F'): Failed to read&convert " +
                     toUtf8(text));
}

// The values of the first `count` fields; the first that does not convert is an error.
std::vector<double> requireFields(const U32& line, std::size_t count) {
    std::vector<double> out;
    for (std::size_t k = 0; k < count; ++k) {
        Field f = field(line, k);
        if (!f.value) {
            fieldError(k, f.text);
        }
        out.push_back(*f.value);
    }
    return out;
}

class TrjTrajectory : public FileTrajectory {
public:
    TrjTrajectory(const std::filesystem::path& path, std::optional<double> dt, std::size_t atoms)
        : m_atoms(atoms) {
        if (atoms == 0) {
            throw ValueError("AMBER TRJ reader REQUIRES the n_atoms keyword");
        }
        const std::string bytes = readBytes(path);
        if (bytes.empty()) {
            throw ReaderError("EOFError",
                              "Compressed file ended before the end-of-stream marker was reached");
        }
        m_lines = splitPythonLines(decodeUtf8Replace(bytes));
        if (rstrip(m_lines.front()).size() > 80) {
            throw ReaderError("OSError",
                              "Header of AMBER formatted trajectory has more than 80 chars. This "
                              "is probably not a AMBER trajectory.");
        }
        m_lineCount = (3 * atoms + kFieldsPerLine - 1) / kFieldsPerLine;
        m_remainder = 3 * atoms % kFieldsPerLine;
        if (m_lines.size() < 2) {
            throw ReaderError("EOFError", "");
        }
        if (atoms > 1) {
            detectBox();
        }
        const std::size_t perFrame = m_lineCount + (m_periodic ? 1 : 0);
        m_frameCount = static_cast<std::int64_t>((m_lines.size() - 1) / perFrame);
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
        const std::vector<double> values = coordinates(first);
        frame.positions.assign(values.size(), 0.0f);
        for (std::size_t i = 0; i < values.size(); ++i) {
            frame.positions[i] = static_cast<float>(values[i]);
        }
        if (m_periodic) {
            const std::vector<double> abc = requireFields(m_lines[first + m_lineCount], 3);
            const std::array<float, 6> box = {static_cast<float>(abc[0]),
                                              static_cast<float>(abc[1]),
                                              static_cast<float>(abc[2]),
                                              90.0f,
                                              90.0f,
                                              90.0f};
            if (box[0] != 0.0f || box[1] != 0.0f || box[2] != 0.0f) {
                frame.box = box;
            }
        }
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    // The coordinates of the frame whose lines start at `first` (lines past the end of the file
    // are not read): a line contributes all ten fields when they all convert, else its first
    // 3N mod 10 fields.
    std::vector<double> coordinates(std::size_t first) const {
        std::vector<double> values;
        for (std::size_t i = first; i < first + m_lineCount && i < m_lines.size(); ++i) {
            const U32& line = m_lines[i];
            std::vector<double> full;
            for (std::size_t k = 0; k < kFieldsPerLine; ++k) {
                const auto value = field(line, k).value;
                if (!value) {
                    break;
                }
                full.push_back(*value);
            }
            if (full.size() != kFieldsPerLine) {
                full = requireFields(line, m_remainder);
            }
            values.insert(values.end(), full.begin(), full.end());
        }
        if (values.size() != 3 * m_atoms) {
            throw ValueError(reshapeMessage(values.size(), m_atoms));
        }
        return values;
    }

    // Periodic when the line after the first frame has exactly three numeric fields.
    void detectBox() {
        coordinates(1);
        const std::size_t next = 1 + m_lineCount;
        if (next >= m_lines.size()) {
            throw ReaderError("StopIteration", "");
        }
        std::size_t numbers = 0;
        for (std::size_t k = 0; k < kFieldsPerLine; ++k) {
            if (field(m_lines[next], k).value) {
                ++numbers;
            }
        }
        m_periodic = numbers == 3;
        if (m_periodic) {
            requireFields(m_lines[next], 3);
        }
    }

    std::size_t m_atoms;
    Lines m_lines;
    std::size_t m_lineCount = 0;  // coordinate lines per frame
    std::size_t m_remainder = 0;  // 3N mod 10
    bool m_periodic = false;
};

// ---------------------------------------------------------------- inpcrd / rst7 (INPReader)

constexpr std::size_t kRestartFieldWidth = 12;

class RestartTrajectory : public FileTrajectory {
public:
    RestartTrajectory(const std::filesystem::path& path, std::optional<double> dt) {
        const Lines lines = readTextLines(path);
        if (lines.size() < 2) {
            throw ReaderError("IndexError", "list index out of range");
        }
        const std::vector<U32> tokens = split(lines[1]);
        if (tokens.empty()) {
            throw ReaderError("IndexError", "list index out of range");
        }
        std::int64_t count = 0;
        try {
            count = core::pythonInt(toUtf8(tokens[0]));
        } catch (const std::overflow_error& error) {
            throw ReaderError("OverflowError", error.what());
        }
        if (count < 0) {
            throw ReaderError("OverflowError", "can't convert negative value to uint64_t");
        }
        if (tokens.size() > 1) {
            m_time = core::pythonFloat(toUtf8(tokens[1]));
        }
        const auto atoms = static_cast<std::size_t>(count);
        for (std::size_t atom = 0; atom < atoms; atom += 2) {
            const std::size_t lineIndex = 2 + atom / 2;
            const U32 line = lineIndex < lines.size() ? lines[lineIndex] : U32();
            const std::size_t fields = atom + 1 < atoms ? 6 : 3;
            for (std::size_t k = 0; k < fields; ++k) {
                const U32 text = slice(line, k * kRestartFieldWidth, (k + 1) * kRestartFieldWidth);
                m_positions.push_back(static_cast<float>(core::pythonFloat(toUtf8(text))));
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
