// File names of generated workflow files (MDATK's mdatk/workflows/naming.py) and the replaceable
// clocks (expected results made with MDATK v0.0.4.1 by tests/parity/python_cases.py; see
// data/README.md).

#include <atomic>
#include <charconv>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <thread>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/naming.hpp"
#include "mdbridge/core/project_file.hpp"

namespace fs = std::filesystem;
using namespace mdbridge::core;

namespace {

const fs::path kData = fs::path(MDBRIDGE_TEST_DATA) / "python";

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

Value readVectors(const std::string& name) { return parseJson(readFile(kData / name), name); }

const Array& items(const Value& value, std::string_view key) {
    return *value.find(key)->get<Array>();
}
const Value& at(const Value& array, std::size_t i) { return (*array.get<Array>())[i]; }
std::string text(const Value& value) { return value.stringOr(); }  // None -> ""
std::int64_t integer(const Value& value) { return *value.get<std::int64_t>(); }

std::string fromHex(std::string_view hex) {
    std::string out;
    for (std::size_t i = 0; i + 1 < hex.size(); i += 2) {
        unsigned byte = 0;
        std::from_chars(hex.data() + i, hex.data() + i + 2, byte, 16);
        out += static_cast<char>(byte);
    }
    return out;
}

std::optional<ParsedName> expectedParts(const Value& parts) {
    if (parts.isNull()) {
        return std::nullopt;
    }
    return ParsedName{text(*parts.find("project")), text(*parts.find("step")),
                      text(*parts.find("tag")), text(*parts.find("stamp")),
                      text(*parts.find("ext"))};
}

DateTime localDateTime(const std::string& iso) {
    return *parseToml("v = " + iso, "datetime").find("v")->get<DateTime>();
}

// Restores the clocks when the test ends.
class ClockGuard {
public:
    explicit ClockGuard(Clock clock) : m_previous(setClock(std::move(clock))) {}
    ~ClockGuard() { setClock(m_previous); }
    ClockGuard(const ClockGuard&) = delete;
    ClockGuard& operator=(const ClockGuard&) = delete;

private:
    Clock m_previous;
};

auto fields(const DateTime& t) {
    return std::make_tuple(t.year, t.month, t.day, t.hour, t.minute, t.second);
}

DateTime fromTm(const std::tm& tm) {
    return {DateTime::Kind::LocalDateTime,
            tm.tm_year + 1900,
            tm.tm_mon + 1,
            tm.tm_mday,
            tm.tm_hour,
            tm.tm_min,
            tm.tm_sec};
}

std::tm localTm(std::time_t t) {
    std::tm out{};
#if defined(_WIN32)
    localtime_s(&out, &t);
#else
    localtime_r(&t, &out);
#endif
    return out;
}

}  // namespace

TEST(Naming, LikeMdatk) {
    const Value vectors = readVectors("naming_cases.json");
    for (const Value& c : items(vectors, "sanitize_project")) {
        EXPECT_EQ(sanitizeProject(text(at(c, 0))), text(at(c, 1))) << text(at(c, 0));
    }
    for (const Value& c : items(vectors, "sanitize_tag")) {
        EXPECT_EQ(sanitizeTag(text(at(c, 0))), text(at(c, 1))) << text(at(c, 0));
    }
    for (const Value& c : items(vectors, "step_label")) {
        EXPECT_EQ(stepLabel(integer(at(c, 0))), text(at(c, 1)));
    }
    for (const Value& c : items(vectors, "step_folder")) {
        EXPECT_EQ(stepFolder(integer(at(c, 0)), text(at(c, 1))), text(at(c, 2)));
    }
    for (const Value& c : items(vectors, "file_name")) {
        EXPECT_EQ(fileName(text(at(c, 0)), integer(at(c, 1)), text(at(c, 2)), text(at(c, 3)),
                           text(at(c, 4))),
                  text(at(c, 5)));
    }
    for (const Value& c : items(vectors, "parse_name")) {
        EXPECT_EQ(parseName(text(at(c, 0))), expectedParts(at(c, 1))) << text(at(c, 0));
    }
    EXPECT_EQ(newStamp(localDateTime("2026-10-15T14:32:59")),
              text(*vectors.find("new_stamp_example")));
    EXPECT_EQ(kStampFormat, "%Y%m%d-%H%M");
}

TEST(Naming, EdgeCasesLikeMdatk) {
    // Unicode digits and whitespace, "$" before a final newline, greedy matching, undecodable bytes
    const Value vectors = readVectors("naming_extra_cases.json");
    for (const Value& c : items(vectors, "parse_name")) {
        const std::string name = fromHex(text(at(c, 0)));
        EXPECT_EQ(parseName(name), expectedParts(at(c, 1))) << name;
    }
    for (const Value& c : items(vectors, "sanitize_project")) {
        EXPECT_EQ(sanitizeProject(fromHex(text(at(c, 0)))), text(at(c, 1))) << text(at(c, 0));
    }
    for (const Value& c : items(vectors, "sanitize_tag")) {
        EXPECT_EQ(sanitizeTag(fromHex(text(at(c, 0)))), text(at(c, 1))) << text(at(c, 0));
    }
    for (const Value& c : items(vectors, "new_stamp")) {
        EXPECT_EQ(newStamp(localDateTime(text(at(c, 0)))), text(at(c, 1))) << text(at(c, 0));
    }
    EXPECT_EQ(stepLabel(-5), "-5");
    EXPECT_EQ(stepLabel(-10), "-10");
    EXPECT_EQ(fileName("p", 1, "t", "", "..."), "p_01_t");
}

TEST(Naming, SanitizeProjectOfTheModelVectors) {
    const Value vectors = readVectors("dirs_cases.json");
    for (const Value& c : items(vectors, "sanitize_project")) {
        EXPECT_EQ(sanitizeProject(text(*c.find("name"))), text(*c.find("sanitized")));
    }
}

TEST(Clock, FixedClocksForTests) {
    const DateTime fixed{DateTime::Kind::LocalDateTime, 2026, 10, 15, 14, 32, 5};
    {
        const ClockGuard guard({[&] { return fixed; }, [] { return 12.5; }});
        EXPECT_EQ(clockNow(), fixed);
        EXPECT_EQ(clockMonotonic(), 12.5);
        EXPECT_EQ(newStamp(), "20261015-1432");
        // setClock hands back what was installed
        Clock mine = setClock({});
        ASSERT_TRUE(mine.now && mine.monotonic);
        EXPECT_EQ(mine.now(), fixed);
        setClock(std::move(mine));
        EXPECT_EQ(clockNow(), fixed);
    }
    EXPECT_NE(clockNow(), fixed);
}

TEST(Clock, EmptyMembersAreTheSystemClocks) {
    const ClockGuard guard({{}, [] { return 1.0; }});
    const std::time_t before = std::time(nullptr);
    const DateTime now = clockNow();
    const std::time_t after = std::time(nullptr);
    EXPECT_EQ(now.kind, DateTime::Kind::LocalDateTime);
    EXPECT_EQ(now.microsecond, 0);
    EXPECT_LE(fields(fromTm(localTm(before))), fields(now));
    EXPECT_LE(fields(now), fields(fromTm(localTm(after))));
    EXPECT_EQ(clockMonotonic(), 1.0);
    setClock({});
    const double a = clockMonotonic();
    const double b = clockMonotonic();
    EXPECT_LE(a, b);
}

TEST(Clock, ReplaceableWhileOtherThreadsReadIt) {
    std::atomic<bool> stop = false;
    std::vector<std::thread> readers;
    std::atomic<int> reads = 0;
    for (int t = 0; t < 4; ++t) {
        readers.emplace_back([&] {
            while (!stop) {
                const DateTime now = clockNow();
                EXPECT_TRUE(now.year == 2000 || now.year == 2001 || now.year > 2020);
                ++reads;
            }
        });
    }
    Clock previous;
    for (int k = 0; k < 2000; ++k) {
        const int year = 2000 + k % 2;
        Clock replaced =
            setClock({[year] { return DateTime{DateTime::Kind::LocalDateTime, year, 1, 1}; }, {}});
        if (k == 0) {
            previous = std::move(replaced);
        }
    }
    stop = true;
    for (std::thread& reader : readers) {
        reader.join();
    }
    setClock(std::move(previous));
    EXPECT_GT(reads.load(), 0);
}
