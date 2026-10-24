// Python's repr() / str() / round(), json.loads and posixpath, and the log sink, as MDATK's Python
// side has them (expected results made with Python 3.13 and MDATK v0.0.4.1 by
// tests/parity/python_cases.py; see data/README.md).

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#if !defined(_WIN32)
#include <pwd.h>
#include <unistd.h>
#endif

#include <gtest/gtest.h>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/naming.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

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
std::string text(const Value& value) { return value.stringOr(); }
std::string at(const Value& array, std::size_t i) { return (*array.get<Array>())[i].stringOr(); }

std::string fromHex(std::string_view hex) {
    std::string out;
    for (std::size_t i = 0; i + 1 < hex.size(); i += 2) {
        unsigned byte = 0;
        std::from_chars(hex.data() + i, hex.data() + i + 2, byte, 16);
        out += static_cast<char>(byte);
    }
    return out;
}

// float.fromhex of float.hex() output, e.g. "-0x1.8p+1"
double fromPythonHex(std::string_view hex) {
    const bool negative = hex.front() == '-';
    hex.remove_prefix(negative ? 3 : 2);
    double value = 0;
    std::from_chars(hex.data(), hex.data() + hex.size(), value, std::chars_format::hex);
    return negative ? -value : value;
}

// Sets (or, with std::nullopt, removes) an environment variable for the life of the object.
class ScopedEnv {
public:
    ScopedEnv(std::string name, const std::optional<std::string>& value) : m_name(std::move(name)) {
        if (const char* old = std::getenv(m_name.c_str())) {
            m_old = old;
        }
        set(value);
    }
    ~ScopedEnv() { set(m_old); }
    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;

private:
    void set(const std::optional<std::string>& value) {
#if defined(_WIN32)
        const std::wstring name(m_name.begin(), m_name.end());
        const std::wstring wide =
            value ? std::wstring(value->begin(), value->end()) : std::wstring();
        _wputenv_s(name.c_str(), wide.c_str());  // an empty value removes the variable
#else
        if (value) {
            ::setenv(m_name.c_str(), value->c_str(), 1);
        } else {
            ::unsetenv(m_name.c_str());
        }
#endif
    }

    std::string m_name;
    std::optional<std::string> m_old;
};

// The variable expanduser("~") reads.
#if defined(_WIN32)
constexpr const char* kHomeVariable = "USERPROFILE";
#else
constexpr const char* kHomeVariable = "HOME";
#endif

// What "~root" expands to here (the vectors were made where it is /root).
std::optional<std::string> rootHome() {
#if defined(_WIN32)
    return std::nullopt;
#else
    const passwd* entry = ::getpwnam("root");
    return entry && entry->pw_dir ? std::optional<std::string>(entry->pw_dir) : std::nullopt;
#endif
}

#if defined(_WIN32)
bool startsWithUserName(std::string_view path) {
    return path.size() > 1 && path[0] == '~' && path[1] != '/';
}
#endif

// The expected message of a ProjectError, or "" when `f` does not throw one.
template <class F>
std::string projectError(F&& f) {
    try {
        f();
    } catch (const ProjectError& error) {
        return error.what();
    }
    return {};
}

}  // namespace

// ------------------------------------------------------------------------------- repr / str / type

TEST(PythonRepr, TomlValuesLikePython) {
    const Value vectors = readVectors("python_repr_cases.json");
    // a Value holds 64-bit integers; i3 = 12345678901234567890 is beyond them
    std::string toml;
    std::istringstream lines(text(*vectors.find("input_toml")));
    for (std::string line; std::getline(lines, line);) {
        if (line.rfind("i3 = ", 0) != 0) {
            toml += line + "\n";
        }
    }
    const Value values = parseToml(toml, "python_repr_cases");
    std::size_t checked = 0;
    for (const Value& c : items(vectors, "cases")) {
        const Value* value = values.find(text(*c.find("key")));
        if (!value) {
            EXPECT_EQ(text(*c.find("key")), "i3");
            continue;
        }
        EXPECT_EQ(pythonRepr(*value), text(*c.find("repr"))) << text(*c.find("key"));
        EXPECT_EQ(pythonStr(*value), text(*c.find("str"))) << text(*c.find("key"));
        EXPECT_EQ(pythonTypeName(*value), text(*c.find("type_name"))) << text(*c.find("key"));
        ++checked;
    }
    EXPECT_EQ(checked, items(vectors, "cases").size() - 1);
}

TEST(PythonRepr, StringsOfEveryPrintableBoundary) {
    // repr() and json.dumps() (ensure_ascii) of each text; undecodable bytes as surrogateescape
    const Value vectors = readVectors("string_repr_cases.json");
    for (const Value& c : items(vectors, "cases")) {
        const std::string s = fromHex(at(c, 0));
        EXPECT_EQ(pythonStringRepr(s), at(c, 1)) << at(c, 0);
        EXPECT_EQ(pythonJsonCompact(s), at(c, 2)) << at(c, 0);
    }
    EXPECT_GT(items(vectors, "cases").size(), 1500u);
}

TEST(PythonRepr, EdgeCases) {
    EXPECT_EQ(pythonRepr(Value()), "None");
    EXPECT_EQ(pythonStr(Value()), "None");
    EXPECT_EQ(pythonTypeName(Value()), "NoneType");
    EXPECT_EQ(pythonRepr(Array{Value(), Table{{"k", Array{}}}, -0.0}), "[None, {'k': []}, -0.0]");
    EXPECT_EQ(pythonStr(Table{{"it's", "x"}}), "{\"it's\": 'x'}");
    EXPECT_EQ(pythonStr(std::int64_t{-9223372036854775807 - 1}), "-9223372036854775808");
    // microseconds keep a zero second; offsets normalise as timedelta does
    DateTime dt{DateTime::Kind::LocalDateTime, 2026, 10, 15, 14, 32, 0, 5, 0};
    EXPECT_EQ(pythonRepr(dt), "datetime.datetime(2026, 10, 15, 14, 32, 0, 5)");
    dt = {DateTime::Kind::OffsetDateTime, 2026, 1, 2, 0, 0, 0, 0, 23 * 60 + 59};
    EXPECT_EQ(pythonRepr(dt),
              "datetime.datetime(2026, 1, 2, 0, 0, "
              "tzinfo=datetime.timezone(datetime.timedelta(seconds=86340)))");
    dt.offsetMinutes = -1;
    EXPECT_EQ(pythonRepr(dt),
              "datetime.datetime(2026, 1, 2, 0, 0, "
              "tzinfo=datetime.timezone(datetime.timedelta(days=-1, seconds=86340)))");
    EXPECT_EQ(pythonStr(dt), "2026-01-02 00:00:00-00:01");
    EXPECT_EQ(pythonRepr(DateTime{DateTime::Kind::LocalTime, 0, 0, 0, 0, 0, 0, 7, 0}),
              "datetime.time(0, 0, 0, 7)");
    // milestone-1 messages use the same quoting
    EXPECT_EQ(pythonStringRepr("4.5 \xc3\x85 \xe2\x80\x93 \x07"),
              "'4.5 \xc3\x85 \xe2\x80\x93 \\x07'");
    EXPECT_EQ(pythonStringRepr("\xc2\xa0\xe2\x80\x8b\xf3\xa0\x80\x81"),
              "'\\xa0\\u200b\\U000e0001'");
}

// ------------------------------------------------------------------------------------------ round

TEST(PythonRound, TwoDigitsLikePython) {
    const Value vectors = readVectors("round2_vectors.json");
    for (const Value& c : *vectors.get<Array>()) {
        const double x = fromPythonHex(text(*c.find("hex")));
        EXPECT_EQ(pythonFloatRepr(x), text(*c.find("repr_in")));
        EXPECT_EQ(pythonFloatRepr(pythonRound(x, 2)), text(*c.find("round2")))
            << text(*c.find("repr_in"));
    }
    EXPECT_EQ(vectors.size(), 732u);
}

TEST(PythonRound, AnyDigitsLikePython) {
    const Value vectors = readVectors("round_cases.json");
    for (const Value& c : items(vectors, "cases")) {
        const Array& v = *c.get<Array>();
        const double x = fromPythonHex(text(v[0]));
        const auto digits = static_cast<int>(*v[1].get<std::int64_t>());
        const std::string expected = text(v[2]);
        if (expected.rfind("OverflowError: ", 0) == 0) {
            EXPECT_THROW(pythonRound(x, digits), std::overflow_error)
                << text(v[0]) << " " << digits;
            continue;
        }
        EXPECT_EQ(pythonFloatRepr(pythonRound(x, digits)), expected) << text(v[0]) << " " << digits;
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(std::isnan(pythonRound(nan, 2)));
    EXPECT_EQ(pythonRound(-std::numeric_limits<double>::infinity(), -400),
              -std::numeric_limits<double>::infinity());
    try {
        pythonRound(std::numeric_limits<double>::max(), -308);
        ADD_FAILURE() << "no overflow";
    } catch (const std::overflow_error& error) {
        EXPECT_STREQ(error.what(), "rounded value too large to represent");
    }
}

// -------------------------------------------------------------------------------------- posixpath

TEST(Posixpath, StringFunctionsLikePython) {
    const Value vectors = readVectors("posixpath_cases.json");
    for (const Value& c : items(vectors, "normpath")) {
        EXPECT_EQ(posixpath::normpath(at(c, 0)), at(c, 1)) << at(c, 0);
    }
    for (const Value& c : items(vectors, "isabs")) {
        EXPECT_EQ(posixpath::isabs(at(c, 0)), *(*c.get<Array>())[1].get<bool>()) << at(c, 0);
    }
    for (const Value& c : items(vectors, "dirname")) {
        EXPECT_EQ(posixpath::dirname(at(c, 0)), at(c, 1)) << at(c, 0);
    }
    for (const Value& c : items(vectors, "basename")) {
        EXPECT_EQ(posixpath::basename(at(c, 0)), at(c, 1)) << at(c, 0);
    }
    for (const Value& c : items(vectors, "splitext")) {
        const Value& parts = (*c.get<Array>())[1];
        EXPECT_EQ(posixpath::splitext(at(c, 0)), std::make_pair(at(parts, 0), at(parts, 1)))
            << at(c, 0);
    }
    for (const Value& c : items(vectors, "join")) {
        EXPECT_EQ(posixpath::join(at(c, 0), at(c, 1)), at(c, 2)) << at(c, 0) << " + " << at(c, 1);
    }
    for (const Value& c : items(vectors, "abspath")) {
        EXPECT_EQ(posixpath::abspath(at(c, 0), at(c, 1)), at(c, 2))
            << at(c, 1) << " + " << at(c, 0);
    }
}

TEST(Posixpath, ExpanduserLikePython) {
    const Value vectors = readVectors("posixpath_cases.json");
    const std::optional<std::string> root = rootHome();
    std::size_t checked = 0;
    for (const Value& c : items(vectors, "expanduser")) {
        const std::string home = at(c, 0);
        const std::string path = at(c, 1);
        if (path.rfind("~root", 0) == 0 && root != "/root") {
            continue;  // made where root's home is /root
        }
#if defined(_WIN32)
        if (home.empty() || startsWithUserName(path)) {
            continue;  // an empty variable cannot be set there; there is no passwd database
        }
#endif
        const ScopedEnv env(kHomeVariable, home);
        EXPECT_EQ(posixpath::expanduser(path), at(c, 2)) << "HOME=" << home << " " << path;
        ++checked;
    }
    EXPECT_GT(checked, 40u);
}

#if !defined(_WIN32)
TEST(Posixpath, ExpanduserWithoutHomeUsesThePasswdEntry) {
    const passwd* entry = ::getpwuid(::getuid());
    ASSERT_NE(entry, nullptr);
    std::string home = entry->pw_dir;
    while (!home.empty() && home.back() == '/') {
        home.pop_back();
    }
    const ScopedEnv env("HOME", std::nullopt);
    EXPECT_EQ(posixpath::expanduser("~/x"), home + "/x");
    EXPECT_EQ(posixpath::expanduser("~"), home.empty() ? "/" : home);
    EXPECT_EQ(posixpath::expanduser("~nosuchuser_xyz/x"), "~nosuchuser_xyz/x");
    EXPECT_EQ(posixpath::expanduser(std::string("~a\0b/x", 6)), std::string("~a\0b/x", 6));
}
#endif

TEST(Posixpath, ResolvesProjectPathsLikeMdatk) {
    // Project.resolve_paths' fix(): "" stays, ~ expanded, relative paths joined to the base and
    // normalised, absolute ones kept as they are
    const Value vectors = readVectors("resolve_paths_cases.json");
    const ScopedEnv env(kHomeVariable, text(*vectors.find("home")));
    const std::optional<std::string> root = rootHome();
    const auto fix = [](const std::string& base, const std::string& p) {
        if (p.empty()) {
            return p;
        }
        const std::string expanded = posixpath::expanduser(p);
        return posixpath::isabs(expanded) ? expanded
                                          : posixpath::normpath(posixpath::join(base, expanded));
    };
    const std::string base = text(*vectors.find("base_dir"));
    const auto skip = [&](const std::string& p) {
#if defined(_WIN32)
        return startsWithUserName(p) && p.rfind("~nosuchuser", 0) != 0;
#else
        return p.rfind("~root", 0) == 0 && root != "/root";
#endif
    };
    for (const Value& c : items(vectors, "fix")) {
        if (!skip(text(*c.find("before")))) {
            EXPECT_EQ(fix(base, text(*c.find("before"))), text(*c.find("after")))
                << text(*c.find("before"));
        }
    }
    for (const Value& c : items(vectors, "fix_relative_base")) {
        EXPECT_EQ(fix(text(*c.find("base_dir")), text(*c.find("before"))), text(*c.find("after")));
    }
    for (const Value& c : items(vectors, "cpptraj_path")) {  // only paths holding a "/" are fixed
        const std::string before = text(*c.find("before"));
        const std::string after =
            before.find('/') != std::string::npos ? fix(base, before) : before;
        EXPECT_EQ(after, text(*c.find("after"))) << before;
    }
    std::size_t changes = 0;
    for (const Value& c : items(vectors, "changes")) {
        if (text(*c.find("path")).find(".files.") == std::string::npos) {  // files: split at ";"
            EXPECT_EQ(fix(base, text(*c.find("before"))), text(*c.find("after")))
                << text(*c.find("path"));
            ++changes;
        }
    }
    EXPECT_GT(changes, 15u);
}

TEST(Posixpath, ProjectDirectoriesLikeMdatk) {
    // Project.out_dir / cache_dir / workflow_dir with the working directory {cwd}
    const Value vectors = readVectors("dirs_cases.json");
    const std::string cwd = "/work/dir";
    const auto expected = [&](const Value& c, std::string_view key) {
        std::string s = text(*c.find(key));
        if (s.rfind("{cwd}", 0) == 0) {
            s.replace(0, 5, cwd);
        }
        return s;
    };
    for (const Value& c : items(vectors, "dirs")) {
        const std::string out = text(*c.find("output.out_dir"));
        const std::string cache = text(*c.find("output.cache_dir"));
        const std::string workflow = text(*c.find("workflow.out_dir"));
        const std::string outDir = posixpath::abspath(out.empty() ? "MDATK_output" : out, cwd);
        EXPECT_EQ(outDir, expected(c, "out_dir"));
        EXPECT_EQ(
            cache.empty() ? posixpath::join(outDir, "MDATK_cache") : posixpath::abspath(cache, cwd),
            expected(c, "cache_dir"));
        const std::string name =
            sanitizeProject(text(*c.find("workflow.project_name"))) + "_workflow";
        EXPECT_EQ(
            workflow.empty() ? posixpath::join(outDir, name) : posixpath::abspath(workflow, cwd),
            expected(c, "workflow_dir"));
    }
}

TEST(Posixpath, CurrentDirectory) {
    EXPECT_EQ(pathFromUtf8(currentDirectory()), fs::current_path());
#if !defined(_WIN32)
    EXPECT_TRUE(posixpath::isabs(currentDirectory()));
#endif
}

// ------------------------------------------------------------------------------------- json.loads

TEST(ParseJson, LikePythonJsonLoads) {
    const Value vectors = readVectors("json_loads_cases.json");
    std::size_t errors = 0;
    for (const Value& c : items(vectors, "cases")) {
        const std::string input = text(*c.find("text"));
        if (const Value* error = c.find("error")) {
            EXPECT_EQ(projectError([&] { parseJson(input, "case"); }),
                      "case: not valid JSON: " + text(*error))
                << input;
            ++errors;
            continue;
        }
        const std::string lossy = c.find("lossy") ? text(*c.find("lossy")) : "";
        if (lossy == "int64") {
            const std::string message = projectError([&] { parseJson(input, "case"); });
            const std::string prefix = "case: not valid JSON: Integer outside the 64-bit range: ";
            EXPECT_EQ(message.rfind(prefix, 0), 0u) << input;
            continue;
        }
        const std::string expected = text(*c.find(lossy == "surrogate" ? "repr_mdbridge" : "repr"));
        try {
            EXPECT_EQ(pythonRepr(parseJson(input, "case")), expected) << input;
        } catch (const ProjectError& e) {
            ADD_FAILURE() << input << ": " << e.what();
        }
    }
    EXPECT_GT(errors, 300u);
}

TEST(ParseJson, KeepsTypesOrderAndExactNumbers) {
    const Value v =
        parseJson(R"({"b": 1, "a": [1.0, -0.0, 2.5e-5, 1e400, NaN], "b": {"c": null}})", "t");
    ASSERT_TRUE(v.isTable());
    // a repeated key keeps its first place and takes its last value
    EXPECT_EQ(v.get<Table>()->front().first, "b");
    EXPECT_EQ(*v.find("b"), (Value(Table{{"c", Value()}})));
    const Array& a = *v.find("a")->get<Array>();
    EXPECT_EQ(a[0].type(), Value::Type::Float);
    EXPECT_TRUE(std::signbit(*a[1].get<double>()));
    EXPECT_EQ(*a[2].get<double>(), 2.5e-5);
    EXPECT_EQ(*a[3].get<double>(), std::numeric_limits<double>::infinity());
    EXPECT_TRUE(std::isnan(*a[4].get<double>()));
    EXPECT_EQ(parseJson("-9223372036854775808", "t"),
              Value(std::int64_t{-9223372036854775807 - 1}));
    EXPECT_EQ(parseJson("\"caf\xe9\"", "t"), Value("caf\xe9"));  // bytes kept as they are
    EXPECT_EQ(parseJson(R"("\u0000")", "t"), Value(std::string(1, '\0')));
}

TEST(ParseJson, LimitsNesting) {
    const auto nested = [](int depth, char open, char close) {
        return std::string(static_cast<std::size_t>(depth), open) + "1" +
               std::string(static_cast<std::size_t>(depth), close);
    };
    EXPECT_NO_THROW(parseJson(nested(1000, '[', ']'), "deep"));
    EXPECT_EQ(projectError([&] { parseJson(nested(1001, '[', ']'), "deep"); }),
              "deep: not valid JSON: maximum recursion depth exceeded while decoding a JSON array "
              "from a unicode string");
    std::string objects;
    for (int k = 0; k < 1001; ++k) {
        objects += "{\"a\":";
    }
    EXPECT_EQ(projectError([&] { parseJson(objects, "deep"); }),
              "deep: not valid JSON: maximum recursion depth exceeded while decoding a JSON object "
              "from a unicode string");
}

TEST(ParseJson, ReadsEveryVectorFile) {
    for (const auto& entry : fs::directory_iterator(kData)) {
        if (entry.path().extension() == ".json") {
            EXPECT_NO_THROW(parseJson(readFile(entry.path()), entry.path().filename().string()))
                << entry.path();
        }
    }
}

// -------------------------------------------------------------------------------------------- log

TEST(Log, DefaultSinkPrintsWarningsAndErrors) {
    const LogSink previous = setLogSink({});
    testing::internal::CaptureStderr();
    logWarning("project: ignoring unknown key 'bogus'");
    logMessage(LogLevel::Error, "broken");
    logInfo("project log moved to x");
    logMessage(LogLevel::Debug, "detail");
    EXPECT_EQ(testing::internal::GetCapturedStderr(),
              "WARNING: project: ignoring unknown key 'bogus'\nERROR: broken\n");
    setLogSink(previous);
    EXPECT_EQ(logLevelName(LogLevel::Debug), "DEBUG");
    EXPECT_EQ(logLevelName(LogLevel::Info), "INFO");
}

TEST(Log, SinksReceiveEveryMessageInOrder) {
    std::vector<std::pair<LogLevel, std::string>> seen;
    const LogSink previous = setLogSink([&](LogLevel level, std::string_view message) {
        seen.emplace_back(level, std::string(message));
    });
    logInfo("one");
    logWarning("two");
    const LogSink mine = setLogSink(previous);
    testing::internal::CaptureStderr();
    logInfo("not captured");  // the default sink again (which drops INFO)
    EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
    ASSERT_EQ(seen.size(), 2u);
    EXPECT_EQ(seen[0], std::make_pair(LogLevel::Info, std::string("one")));
    EXPECT_EQ(seen[1], std::make_pair(LogLevel::Warning, std::string("two")));
    ASSERT_TRUE(mine);
    mine(LogLevel::Error, "direct");
    EXPECT_EQ(seen.size(), 3u);
}

TEST(Log, SinksMayLogAndReplaceThemselves) {
    std::vector<std::string> seen;
    LogSink previous;
    previous = setLogSink([&](LogLevel, std::string_view message) {
        seen.emplace_back(message);
        if (message == "first") {
            logInfo("from the sink");
            setLogSink([&](LogLevel, std::string_view m) {
                seen.push_back("second sink: " + std::string(m));
            });
        }
    });
    logInfo("first");
    logInfo("after");
    setLogSink(previous);
    EXPECT_EQ(seen, (std::vector<std::string>{"first", "from the sink", "second sink: after"}));
}

TEST(Log, MessagesFromManyThreadsReachTheSinkOneAtATime) {
    std::vector<std::string> seen;  // not synchronised itself
    const LogSink previous =
        setLogSink([&](LogLevel, std::string_view m) { seen.emplace_back(m); });
    std::vector<std::thread> threads;
    for (int t = 0; t < 8; ++t) {
        threads.emplace_back([t] {
            for (int k = 0; k < 500; ++k) {
                logWarning("thread " + std::to_string(t) + " message " + std::to_string(k));
            }
        });
    }
    for (std::thread& thread : threads) {
        thread.join();
    }
    setLogSink(previous);
    EXPECT_EQ(seen.size(), 4000u);
}
