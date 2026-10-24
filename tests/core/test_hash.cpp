// SHA-256, stable_hash and file_fingerprint as MDATK's config.py computes them (expected results
// made with MDATK v0.0.4.1 by tests/parity/python_cases.py; see data/README.md).

#include <charconv>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

#if !defined(_WIN32)
#include <fcntl.h>
#include <sys/stat.h>
#endif

#include <gtest/gtest.h>

#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/python.hpp"

namespace fs = std::filesystem;
using namespace mdbridge::core;

namespace {

const fs::path kData = fs::path(MDBRIDGE_TEST_DATA) / "python";
constexpr std::size_t kLengths[] = {0, 8, 16, 64, 100};  // the hash lengths of the vectors

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

Value readVectors(const std::string& name) { return parseJson(readFile(kData / name), name); }

std::string text(const Value* value) { return value ? value->stringOr() : std::string(); }

std::string fromHex(std::string_view hex) {
    std::string out;
    for (std::size_t i = 0; i + 1 < hex.size(); i += 2) {
        unsigned byte = 0;
        std::from_chars(hex.data() + i, hex.data() + i + 2, byte, 16);
        out += static_cast<char>(byte);
    }
    return out;
}

// The Value of a tagged input of stable_hash_inputs.json; false for an integer beyond 64 bits.
bool untag(const Value& in, Value& out) {
    if (const Table* table = in.get<Table>()) {
        if (table->size() == 1) {
            const auto& [tag, value] = table->front();
            if (tag == "$float") {
                const std::string name = value.stringOr();
                const double inf = std::numeric_limits<double>::infinity();
                out = name == "nan"   ? std::numeric_limits<double>::quiet_NaN()
                      : name == "inf" ? inf
                                      : -inf;
                return true;
            }
            if (tag == "$datetime" || tag == "$date") {
                out = *parseToml("v = " + value.stringOr(), "tag").find("v");
                return true;
            }
            if (tag == "$int") {
                return false;
            }
        }
        Table result;
        for (const auto& [key, value] : *table) {
            Value item;
            if (!untag(value, item)) {
                return false;
            }
            result.emplace_back(key, std::move(item));
        }
        out = std::move(result);
        return true;
    }
    if (const Array* array = in.get<Array>()) {
        Array result;
        for (const Value& value : *array) {
            Value item;
            if (!untag(value, item)) {
                return false;
            }
            result.push_back(std::move(item));
        }
        out = std::move(result);
        return true;
    }
    out = in;
    return true;
}

// what _to_plain does to the None values of dicts (the test's own copy of the rule)
Value withoutNulls(const Value& value) {
    if (const Table* table = value.get<Table>()) {
        Table out;
        for (const auto& [key, item] : *table) {
            if (!item.isNull()) {
                out.emplace_back(key, withoutNulls(item));
            }
        }
        return out;
    }
    if (const Array* array = value.get<Array>()) {
        Array out;
        for (const Value& item : *array) {
            out.push_back(withoutNulls(item));
        }
        return out;
    }
    return value;
}

// Sets the modification time of a file in nanoseconds since the Unix epoch (Windows: in 100 ns).
void setModified(const fs::path& path, std::int64_t ns) {
#if defined(_WIN32)
    const auto ticks = ns / 100 + 116'444'736'000'000'000;
    fs::last_write_time(path, fs::file_time_type(fs::file_time_type::duration(ticks)));
#else
    timespec times[2];
    times[0].tv_sec = times[1].tv_sec = static_cast<time_t>(ns / 1'000'000'000);
    times[0].tv_nsec = times[1].tv_nsec = static_cast<long>(ns % 1'000'000'000);
    ASSERT_EQ(::utimensat(AT_FDCWD, path.c_str(), times, 0), 0);
#endif
}

class TempDir {
public:
    TempDir() {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        m_path = fs::temp_directory_path() / ("mdbridge-hash-" + std::to_string(stamp));
        fs::create_directories(m_path);
    }
    ~TempDir() {
        std::error_code ec;
        fs::remove_all(m_path, ec);
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    const fs::path& path() const { return m_path; }

private:
    fs::path m_path;
};

}  // namespace

TEST(Sha256, Fips180Examples) {
    EXPECT_EQ(sha256Hex(""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    EXPECT_EQ(sha256Hex("abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    EXPECT_EQ(sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
              "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    EXPECT_EQ(
        sha256Hex(
            "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmnoijklmnopjklmnopqklmnop"
            "qrlmnopqrsmnopqrstnopqrstu"),
        "cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1");
    EXPECT_EQ(sha256Hex(std::string(1'000'000, 'a')),  // 15625 blocks
              "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

TEST(Sha256, LikeHashlibAtEveryBlockBoundary) {
    const Value vectors = readVectors("sha256_cases.json");
    for (const Value& c : *vectors.find("cases")->get<Array>()) {
        const Array& pair = *c.get<Array>();
        EXPECT_EQ(sha256Hex(fromHex(pair[0].stringOr())), pair[1].stringOr()) << pair[0].stringOr();
    }
    EXPECT_GT(vectors.find("cases")->size(), 200u);
}

TEST(StableHash, LikeMdatkForEveryCase) {
    const Value hashes = readVectors("stable_hash_cases.json");
    const Value inputs = readVectors("stable_hash_inputs.json");
    const Array& cases = *hashes.find("cases")->get<Array>();
    const Array& data = *inputs.find("cases")->get<Array>();
    ASSERT_EQ(cases.size(), data.size());
    std::size_t checked = 0;
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const std::string name = text(cases[i].find("name"));
        ASSERT_EQ(text(data[i].find("name")), name);
        const std::string blob = text(cases[i].find("blob"));
        EXPECT_EQ(sha256Hex(blob), text(cases[i].find("sha256"))) << name;
        Value input;
        if (!untag(*data[i].find("input"), input)) {
            EXPECT_EQ(name, "big_int");  // a Value holds 64-bit integers only
            continue;
        }
        if (*data[i].find("dataclass")->get<bool>()) {
            // _to_plain of a dataclass keeps its None fields: hashed as given, nothing dropped
            EXPECT_EQ(pythonJsonCompact(input), blob) << name;
            for (const std::size_t n : kLengths) {
                EXPECT_EQ(sha256Hex(pythonJsonCompact(input)).substr(0, n),
                          text(cases[i].find("stable_hash_" + std::to_string(n))))
                    << name;
            }
        } else {
            EXPECT_EQ(pythonJsonCompact(withoutNulls(input)), blob) << name;
            EXPECT_EQ(stableHash(input), text(cases[i].find("stable_hash_16"))) << name;
            for (const std::size_t n : kLengths) {
                EXPECT_EQ(stableHash(input, n),
                          text(cases[i].find("stable_hash_" + std::to_string(n))))
                    << name;
            }
        }
        ++checked;
    }
    EXPECT_EQ(checked, cases.size() - 1);
}

TEST(StableHash, DropsNullsOfTablesAtEveryDepthButNotOfArrays) {
    const Value withNulls =
        Table{{"a", Value()}, {"b", Array{Value(), Table{{"c", Value()}, {"d", 1}}}}};
    const Value without = Table{{"b", Array{Value(), Table{{"d", 1}}}}};
    EXPECT_EQ(stableHash(withNulls), stableHash(without));
    EXPECT_EQ(stableHash(without, 64), sha256Hex(R"({"b":[null,{"d":1}]})"));
    EXPECT_EQ(
        pythonJsonCompact(Table{{"t", DateTime{DateTime::Kind::LocalTime, 0, 0, 0, 7, 8, 0}}}),
        R"({"t":"07:08:00"})");
}

TEST(FileFingerprint, LikeMdatk) {
    const Value vectors = readVectors("stable_hash_cases.json");
    const Value& expected = *vectors.find("file_fingerprint");
    const TempDir dir;
    const fs::path file = dir.path() / "fp-file";
    std::ofstream(file, std::ios::binary) << "hello";
    const std::int64_t ns = *expected.find("mtime_set_ns")->get<std::int64_t>();
    setModified(file, ns);

    const Value fingerprint = fileFingerprint(file);
    std::vector<std::string> keys;
    for (const auto& [key, value] : *fingerprint.get<Table>()) {
        keys.push_back(key);
    }
    EXPECT_EQ(keys, (std::vector<std::string>{"path", "size", "mtime_ns"}));
    EXPECT_EQ(*fingerprint.find("size"), *expected.path({"result", "size"}));
#if defined(_WIN32)
    EXPECT_EQ(*fingerprint.find("mtime_ns"), Value(ns / 100 * 100));  // FILETIME has 100 ns steps
    EXPECT_EQ(text(fingerprint.find("path")), pathToUtf8(fs::absolute(file)));
#else
    EXPECT_EQ(*fingerprint.find("mtime_ns"), *expected.path({"result", "mtime_ns"}));
    EXPECT_EQ(text(fingerprint.find("path")), pathToUtf8(file));
    // os.path.abspath: normalised as text, symbolic links not resolved
    fs::create_directory(dir.path() / "sub");
    const fs::path roundabout = dir.path() / "." / "sub" / ".." / "fp-file";
    EXPECT_EQ(text(fileFingerprint(roundabout).find("path")), pathToUtf8(file));
#endif
}

TEST(FileFingerprint, MissingFileThrows) {
    const TempDir dir;
    EXPECT_THROW(fileFingerprint(dir.path() / "missing"), fs::filesystem_error);
}
