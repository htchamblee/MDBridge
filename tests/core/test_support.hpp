#pragma once
// Helpers shared by the core tests: files, temporary folders, a fake clock, a log capture and the
// tagged values of the vector files that tests/parity/make_writer_vectors.py writes.

#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::test {

namespace fs = std::filesystem;

inline const fs::path kData = MDBRIDGE_TEST_DATA;

inline std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("cannot read " + path.string());
    }
    std::ostringstream text;
    text << in.rdbuf();
    return std::move(text).str();
}

inline void writeFile(const fs::path& path, std::string_view text) {
    std::ofstream out(path, std::ios::binary);
    out.write(text.data(), static_cast<std::streamsize>(text.size()));
}

// The text MDATK wrote on Linux, with the line ends this platform's text files get.
inline std::string nativeLines(std::string_view text) {
#ifdef _WIN32
    std::string out;
    for (const char c : text) {
        if (c == '\n') {
            out += '\r';
        }
        out += c;
    }
    return out;
#else
    return std::string(text);
#endif
}

// A fresh folder, removed again at the end of the test.
class TempDir {
public:
    TempDir() {
        static int counter = 0;
        const fs::path base = fs::temp_directory_path();
        for (;;) {
            m_path = base / ("mdbridge-test-" + std::to_string(std::rand()) + "-" + std::to_string(++counter));
            if (fs::create_directory(m_path)) {
                break;
            }
        }
    }
    ~TempDir() {
        std::error_code ec;
        fs::remove_all(m_path, ec);
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;

    const fs::path& path() const { return m_path; }
    fs::path operator/(const std::string& name) const { return m_path / name; }

private:
    fs::path m_path;
};

// The clocks of make_project_file_golden.py: wall time from `start`, one second later on every
// call; monotonic readings from a list.
class FakeClock {
public:
    FakeClock(core::DateTime start, std::vector<double> monotonic = {})
        : m_state(std::make_shared<State>(State{start, std::move(monotonic), 0})) {
        auto state = m_state;
        m_previous = core::setClock(core::Clock{
            [state] {
                core::DateTime now = state->next;
                if (++state->next.second == 60) {
                    state->next.second = 0;
                    if (++state->next.minute == 60) {
                        state->next.minute = 0;
                        ++state->next.hour;
                    }
                }
                return now;
            },
            [state] {
                if (state->index >= state->monotonic.size()) {
                    throw std::logic_error("FakeClock: no monotonic reading left");
                }
                return state->monotonic[state->index++];
            }});
    }
    ~FakeClock() { core::setClock(m_previous); }
    FakeClock(const FakeClock&) = delete;
    FakeClock& operator=(const FakeClock&) = delete;

    std::size_t monotonicReadings() const { return m_state->index; }

private:
    struct State {
        core::DateTime next;
        std::vector<double> monotonic;
        std::size_t index;
    };
    std::shared_ptr<State> m_state;
    core::Clock m_previous;
};

inline core::DateTime localDateTime(int year, int month, int day, int hour, int minute, int second) {
    core::DateTime t;
    t.kind = core::DateTime::Kind::LocalDateTime;
    t.year = year;
    t.month = month;
    t.day = day;
    t.hour = hour;
    t.minute = minute;
    t.second = second;
    return t;
}

// Collects the messages sent to the core log while it lives.
class LogCapture {
public:
    LogCapture() {
        m_previous = core::setLogSink([this](core::LogLevel level, std::string_view message) {
            m_lines.push_back(std::string(core::logLevelName(level)) + ": " + std::string(message));
        });
    }
    ~LogCapture() { core::setLogSink(m_previous); }
    LogCapture(const LogCapture&) = delete;
    LogCapture& operator=(const LogCapture&) = delete;

    const std::vector<std::string>& lines() const { return m_lines; }

private:
    std::vector<std::string> m_lines;
    core::LogSink m_previous;
};

inline core::Value readJson(const fs::path& path) { return core::parseJson(readFile(path), path.string()); }

inline const std::string& text(const core::Value& value, std::string_view key) {
    static const std::string empty;
    const core::Value* entry = value.find(key);
    const std::string* s = entry ? entry->get<std::string>() : nullptr;
    return s ? *s : empty;
}

inline std::string fromHex(std::string_view hex) {
    std::string out;
    for (std::size_t i = 0; i + 1 < hex.size(); i += 2) {
        out += static_cast<char>(std::stoi(std::string(hex.substr(i, 2)), nullptr, 16));
    }
    return out;
}

inline core::DateTime parseIso(std::string_view iso, core::DateTime::Kind kind) {
    core::DateTime t;
    t.kind = kind;
    const std::string s(iso);
    if (kind == core::DateTime::Kind::LocalTime) {
        std::sscanf(s.c_str(), "%d:%d:%d", &t.hour, &t.minute, &t.second);
    } else if (kind == core::DateTime::Kind::LocalDate) {
        std::sscanf(s.c_str(), "%d-%d-%d", &t.year, &t.month, &t.day);
    } else {
        std::sscanf(s.c_str(), "%d-%d-%dT%d:%d:%d", &t.year, &t.month, &t.day, &t.hour, &t.minute, &t.second);
    }
    if (const auto dot = s.find('.'); dot != std::string::npos) {
        std::string digits;
        for (std::size_t i = dot + 1; i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])); ++i) {
            digits += s[i];
        }
        digits.resize(6, '0');
        t.microsecond = std::stoi(digits);
    }
    const auto zone = s.find_first_of("+-Z", s.find('T') == std::string::npos ? s.size() : s.find('T'));
    if (kind == core::DateTime::Kind::LocalDateTime && zone != std::string::npos) {
        t.kind = core::DateTime::Kind::OffsetDateTime;
        if (s[zone] != 'Z') {
            int h = 0, m = 0;
            std::sscanf(s.c_str() + zone + 1, "%d:%d", &h, &m);
            t.offsetMinutes = (s[zone] == '-' ? -1 : 1) * (h * 60 + m);
        }
    }
    return t;
}

// A value of the vector files: {"t": "none" | "bool" | "int" | "float" | "str" | "datetime" | "date" |
// "time" | "array" | "table", ...} (see make_writer_vectors.py).
inline core::Value tagged(const core::Value& v) {
    const std::string& t = text(v, "t");
    if (t == "none") {
        return {};
    }
    if (t == "bool") {
        return *v.find("v");
    }
    if (t == "int") {
        return core::Value(std::int64_t{std::stoll(text(v, "v"))});
    }
    if (t == "float") {
        const std::string& hex = text(v, "hex");
        if (hex == "nan") {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (hex == "inf" || hex == "-inf") {
            return (hex[0] == '-' ? -1 : 1) * std::numeric_limits<double>::infinity();
        }
        return std::strtod(hex.c_str(), nullptr);
    }
    if (t == "str") {
        if (const core::Value* repeat = v.find("repeat")) {
            const auto count = static_cast<std::size_t>(*v.find("count")->get<std::int64_t>());
            std::string out;
            for (std::size_t i = 0; i < count; ++i) {
                out += *repeat->get<std::string>();
            }
            return out;
        }
        return fromHex(text(v, "utf8"));
    }
    if (t == "datetime") {
        return parseIso(text(v, "iso"), core::DateTime::Kind::LocalDateTime);
    }
    if (t == "date") {
        return parseIso(text(v, "iso"), core::DateTime::Kind::LocalDate);
    }
    if (t == "time") {
        return parseIso(text(v, "iso"), core::DateTime::Kind::LocalTime);
    }
    if (t == "array") {
        core::Array out;
        for (const core::Value& item : *v.find("v")->get<core::Array>()) {
            out.push_back(tagged(item));
        }
        return out;
    }
    if (t == "table") {
        core::Table out;
        for (const core::Value& pair : *v.find("v")->get<core::Array>()) {
            const core::Array& kv = *pair.get<core::Array>();
            out.emplace_back(*kv[0].get<std::string>(), tagged(kv[1]));
        }
        return out;
    }
    throw std::runtime_error("unknown tag " + t);
}

}  // namespace mdbridge::test
