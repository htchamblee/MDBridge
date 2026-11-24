// Amber ASCII coordinates as MDAnalysis 2.11.0 reads them: TRJReader (mdcrd: 10F8.3 lines, a
// box line after each frame when the line after the first frame has exactly three F8.3 fields)
// and INPReader (rst7 / inpcrd: one frame, time from the second line, never a box).

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "binary_file.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

using core::ValueError;

std::string utf8(std::u32string_view text) {
    return core::encodeUtf8(text);
}

std::u32string rstrip(std::u32string_view text) {
    std::size_t last = text.size();
    while (last > 0 && core::isPythonSpace(text[last - 1])) {
        --last;
    }
    return std::u32string(text.substr(0, last));
}

// Lines of a text file as Python's universal-newline readline() yields them ("\n"-terminated but
// for a last line without one), read in chunks, with the byte offset where each line starts.
class LineReader {
public:
    LineReader(BinaryFile& file, std::uint64_t offset)
        : m_file(file), m_offset(offset), m_bufferStart(offset) {}

    // false at the end of the file
    bool next(std::string& line, std::uint64_t& start) {
        line.clear();
        start = m_offset;
        for (;;) {
            const int c = get();
            if (c < 0) {
                return !line.empty();
            }
            if (c == '\n' || c == '\r') {
                if (c == '\r' && peek() == '\n') {
                    get();
                }
                line += '\n';
                return true;
            }
            line += static_cast<char>(c);
        }
    }

    std::uint64_t offset() const { return m_offset; }

private:
    int peek() {
        if (m_at == m_buffer.size()) {
            m_bufferStart += m_buffer.size();
            m_buffer = m_file.read(m_bufferStart, 1 << 16);
            m_at = 0;
        }
        return m_at < m_buffer.size() ? static_cast<unsigned char>(m_buffer[m_at]) : -1;
    }
    int get() {
        const int c = peek();
        if (c >= 0) {
            ++m_at;
            ++m_offset;
        }
        return c;
    }

    BinaryFile& m_file;
    std::uint64_t m_offset;
    std::uint64_t m_bufferStart;
    std::string m_buffer;
    std::size_t m_at = 0;
};

// MDAnalysis' FORTRANReader("<count>F8.3") entries: float(line[8 * i:8 * (i + 1)]).
std::optional<double> fieldValue(const std::u32string& line, std::size_t i, std::string* failed) {
    const std::size_t start = std::min(line.size(), 8 * i);
    const std::u32string text = line.substr(start, 8);
    if (auto value = core::tryPythonFloat(utf8(text))) {
        return value;
    }
    if (failed != nullptr) {
        *failed = "FixedcolumnEntry(" + std::to_string(8 * i) + "," + std::to_string(8 * i + 8) +
                  ",'F'): Failed to read&convert " + utf8(text);
    }
    return std::nullopt;
}

// FORTRANReader.read: every entry or the first entry's ValueError.
std::vector<double> readFields(const std::u32string& line, std::size_t count) {
    std::vector<double> out;
    std::string failed;
    for (std::size_t i = 0; i < count; ++i) {
        const auto value = fieldValue(line, i, &failed);
        if (!value) {
            throw ValueError(failed);
        }
        out.push_back(*value);
    }
    return out;
}

std::size_t matchingFields(const std::u32string& line, std::size_t count) {
    std::size_t matches = 0;
    for (std::size_t i = 0; i < count; ++i) {
        matches += fieldValue(line, i, nullptr) ? 1 : 0;
    }
    return matches;
}

std::u32string decodeLine(const std::string& bytes) {
    return decodeUtf8Replace(bytes);
}

class TrjTrajectory : public FileTrajectory {
public:
    TrjTrajectory(const std::filesystem::path& path, std::optional<double> dt, std::size_t atoms)
        : m_file(path), m_atoms(atoms) {
        if (atoms == 0) {
            throw ValueError("AMBER TRJ reader REQUIRES the n_atoms keyword");
        }
        if (m_file.size() == 0) {
            throw ReaderError("EOFError",
                              "Compressed file ended before the end-of-stream marker was reached");
        }
        m_dt = dt.value_or(1.0);
        m_linesPerFrame = (3 * atoms + 9) / 10;
        m_lastPerLine = 3 * atoms % 10;

        LineReader reader(m_file, 0);
        std::string line;
        std::uint64_t start = 0;
        reader.next(line, start);
        m_header = decodeLine(line);
        if (rstrip(m_header).size() > 80) {
            throw ReaderError(
                "OSError",
                "Header of AMBER formatted trajectory has more than 80 chars. This is "
                "probably not a AMBER trajectory.");
        }
        m_dataStart = reader.offset();
        detectBox();

        // frame offsets: every linesPerFrame (+1 with a box) lines after the header
        const std::size_t perFrame = m_linesPerFrame + (m_periodic ? 1 : 0);
        LineReader scan(m_file, m_dataStart);
        std::size_t count = 0;
        while (scan.next(line, start)) {
            if (count % perFrame == 0) {
                m_offsets.push_back(start);
            }
            ++count;
        }
        m_frameCount = static_cast<std::int64_t>(count / perFrame);
        m_offsets.resize(static_cast<std::size_t>(m_frameCount));
        if (m_frameCount == 0) {
            core::Frame frame;
            readAt(m_dataStart, frame);  // the error MDAnalysis' first read gives
            throw ReaderError("EOFError", "");
        }
        seek(0);
    }

    core::TrajectoryInfo info() const override {
        core::TrajectoryInfo info;
        info.header = utf8(m_header);
        return info;
    }

    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        readAt(m_offsets[static_cast<std::size_t>(index)], frame);
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    // One frame's coordinates from `offset` (TRJReader._read_next_timestep); the box when periodic.
    void readAt(std::uint64_t offset, core::Frame& frame) {
        LineReader reader(m_file, offset);
        std::string bytes;
        std::uint64_t start = 0;
        std::vector<double> values;
        for (std::size_t number = 0; number < m_linesPerFrame && reader.next(bytes, start);
             ++number) {
            const std::u32string line = decodeLine(bytes);
            if (matchingFields(line, 10) == 10) {
                const auto fields = readFields(line, 10);
                values.insert(values.end(), fields.begin(), fields.end());
            } else {
                const auto fields = readFields(line, m_lastPerLine);
                values.insert(values.end(), fields.begin(), fields.end());
            }
        }
        if (values.empty()) {
            throw ReaderError("EOFError", "");
        }
        if (m_periodic) {
            if (!reader.next(bytes, start)) {
                throw ReaderError("StopIteration", "");
            }
            const auto box = readFields(decodeLine(bytes), 3);
            frame.box = std::array<float, 6>{static_cast<float>(box[0]),
                                             static_cast<float>(box[1]),
                                             static_cast<float>(box[2]),
                                             90.0f,
                                             90.0f,
                                             90.0f};
            if ((*frame.box)[0] == 0.0f && (*frame.box)[1] == 0.0f && (*frame.box)[2] == 0.0f) {
                frame.box.reset();
            }
        }
        if (values.size() != 3 * m_atoms) {
            throw ValueError("cannot reshape array of size " + std::to_string(values.size()) +
                             " into shape (" + std::to_string(m_atoms) + ",3)");
        }
        frame.positions.resize(values.size());
        std::transform(values.begin(), values.end(), frame.positions.begin(),
                       [](double v) { return static_cast<float>(v); });
    }

    void detectBox() {
        if (m_atoms == 1) {
            return;  // a box line cannot be told from a one-atom frame
        }
        core::Frame frame;
        LineReader reader(m_file, m_dataStart);
        std::string bytes;
        std::uint64_t start = 0;
        std::size_t lines = 0;
        while (lines < m_linesPerFrame && reader.next(bytes, start)) {
            ++lines;
        }
        readAt(m_dataStart, frame);
        if (!reader.next(bytes, start)) {
            throw ReaderError("StopIteration", "");
        }
        const std::u32string line = decodeLine(bytes);
        if (matchingFields(line, 10) == 3) {
            m_periodic = true;
            readFields(line, 3);
        }
    }

    BinaryFile m_file;
    std::size_t m_atoms;
    std::size_t m_linesPerFrame = 0;
    std::size_t m_lastPerLine = 0;
    std::u32string m_header;
    std::uint64_t m_dataStart = 0;
    bool m_periodic = false;
    std::vector<std::uint64_t> m_offsets;
};

std::vector<std::string> splitWhitespace(std::u32string_view text) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && core::isPythonSpace(text[i])) {
            ++i;
        }
        const std::size_t start = i;
        while (i < text.size() && !core::isPythonSpace(text[i])) {
            ++i;
        }
        if (i > start) {
            out.push_back(utf8(text.substr(start, i - start)));
        }
    }
    return out;
}

class RestartTrajectory : public FileTrajectory {
public:
    RestartTrajectory(const std::filesystem::path& path, std::optional<double> dt) {
        const Lines lines = readTextLines(path);
        m_dt = dt.value_or(1.0);
        m_frameCount = 1;
        const std::vector<std::string> second =
            lines.size() > 1 ? splitWhitespace(lines[1]) : std::vector<std::string>{};
        if (second.empty()) {
            throw ReaderError("IndexError", "list index out of range");
        }
        const std::int64_t atoms = core::pythonInt(second[0]);
        if (atoms < 0) {
            throw ValueError("negative dimensions are not allowed");
        }
        core::Frame frame;
        if (second.size() > 1) {
            m_time = core::pythonFloat(second[1]);
        }
        const auto n = static_cast<std::size_t>(atoms);
        frame.positions.assign(3 * n, 0.0f);
        std::size_t next = 2;
        const auto field = [&](const std::u32string& line, std::size_t i) {
            const std::size_t start = std::min(line.size(), 12 * i);
            return static_cast<float>(core::pythonFloat(utf8(line.substr(start, 12))));
        };
        for (std::size_t p = 0; p < n / 2; ++p) {
            const std::u32string line = next < lines.size() ? lines[next++] : std::u32string();
            for (std::size_t i = 0; i < 6; ++i) {
                frame.positions[6 * p + i] = field(line, i);
            }
        }
        if (n % 2 == 1) {
            const std::u32string line = next < lines.size() ? lines[next++] : std::u32string();
            for (std::size_t i = 0; i < 3; ++i) {
                frame.positions[3 * (n - 1) + i] = field(line, i);
            }
        }
        frame.time = m_time.value_or(0.0);
        m_positions = std::move(frame.positions);
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return m_time.has_value(); }

protected:
    void readFrame(std::int64_t, core::Frame& frame) override {
        frame.positions = m_positions;
        frame.time = m_time.value_or(0.0);
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
