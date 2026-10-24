// The step log of a project, as MDATK v0.0.4.1 writes it: entry texts (vectors of
// tests/parity/make_writer_vectors.py) and C++ ports of tests/test_projectfile.py.

#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/core/text.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

std::vector<std::string> strings(const Value& array) {
    std::vector<std::string> out;
    for (const Value& item : *array.get<Array>()) {
        out.push_back(*item.get<std::string>());
    }
    return out;
}

std::vector<std::string> steps(const fs::path& project) {
    std::vector<std::string> out;
    for (const Value& entry : readLog(project)) {
        out.push_back(entry.find("step") ? entry.find("step")->stringOr() : std::string());
    }
    return out;
}

// Restores the default split size at the end of a test.
struct SplitBytes {
    explicit SplitBytes(std::uintmax_t bytes) { setLogSplitBytes(bytes); }
    ~SplitBytes() { setLogSplitBytes(kDefaultLogSplitBytes); }
};

}  // namespace

TEST(StepLogVectors, EntryText) {
    const Value cases = readJson(kData / "project-file" / "vectors" / "entry_vectors.json");
    ASSERT_EQ(cases.size(), 14u);
    for (const Value& c : *cases.get<Array>()) {
        const Value& fields = *c.find("entry");
        LogEntry entry;
        entry.step = *tagged(*fields.find("step")).get<std::string>();
        entry.title = *tagged(*fields.find("title")).get<std::string>();
        entry.status = *tagged(*fields.find("status")).get<std::string>();
        entry.time = *tagged(*fields.find("time")).get<DateTime>();
        entry.duration_s = *tagged(*fields.find("duration_s")).get<double>();
        entry.command = *tagged(*fields.find("command")).get<std::string>();
        entry.folder = *tagged(*fields.find("folder")).get<std::string>();
        entry.files = strings(tagged(*fields.find("files")));
        entry.summary = *tagged(*fields.find("summary")).get<std::string>();
        entry.citations = strings(tagged(*fields.find("citations")));
        entry.output = *tagged(*fields.find("output")).get<std::string>();
        const bool withDoc = *c.find("with_doc")->get<bool>();

        const std::string out = logEntryText(entry, withDoc);
        if (c.find("output")) {
            EXPECT_EQ(out, text(c, "output")) << "step " << entry.step << " with_doc " << withDoc;
        } else {
            EXPECT_EQ(characterCount(out), static_cast<std::size_t>(*c.find("output_len_codepoints")->get<std::int64_t>()));
            const std::u32string codes = decodeUtf8(out);
            const std::u32string head = decodeUtf8(text(c, "output_head"));
            const std::u32string tail = decodeUtf8(text(c, "output_tail"));
            EXPECT_EQ(codes.substr(0, head.size()), head) << "step " << entry.step;
            EXPECT_EQ(codes.substr(codes.size() - tail.size()), tail) << "step " << entry.step;
        }
    }
}

TEST(StepLog, PlainRoundsAndDropsEmptyFields) {
    LogEntry entry("02", "tleap build", {"key1"});
    entry.duration_s = 2.675;  // 2.67499999... in binary: round() gives 2.67
    entry.status = "";
    const Value plain = entry.plain();
    std::vector<std::string> keys;
    for (const auto& [key, value] : *plain.get<Table>()) {
        keys.push_back(key);
    }
    EXPECT_EQ(keys, (std::vector<std::string>{"step", "title", "time", "duration_s", "citations"}));
    EXPECT_EQ(*plain.find("duration_s")->get<double>(), 2.67);
}

// tests/test_projectfile.py::test_step_log_append_keep_and_split
TEST(StepLog, AppendKeepAndSplit) {
    TempDir dir;
    const fs::path p = dir / "p.mdatk.toml";
    saveProject(p, Project{});
    logged(p, "02", "tleap build", {"key1"}, [](LogEntry& entry) {
        entry.command = "bash build.sh";
        entry.output = std::string(20050, 'x');
    });
    EXPECT_THROW(logged(p, "03", "md", {}, [](LogEntry&) { throw std::runtime_error("pmemd crashed"); }),
                 std::runtime_error);

    Array log = readLog(p);
    ASSERT_EQ(steps(p), (std::vector<std::string>{"02", "03"}));
    EXPECT_EQ(log[0].find("status")->stringOr(), "ok");
    EXPECT_EQ(*log[0].find("citations"), Value(Array{Value("key1")}));
    const std::string output = log[0].find("output")->stringOr();
    EXPECT_EQ(output.rfind("[...]", 0), 0u);
    EXPECT_EQ(characterCount(output), 20006u);
    EXPECT_EQ(log[1].find("status")->stringOr(), "failed");
    EXPECT_NE(log[1].find("summary")->stringOr().find("pmemd crashed"), std::string::npos);

    Project q = loadProject(p);
    q.notes = "edited";
    saveProject(p, q);
    EXPECT_EQ(steps(p).size(), 2u);
    EXPECT_EQ(loadProject(p).notes, "edited");

    {
        SplitBytes small(1000);
        logged(p, "04", "analysis", {}, [](LogEntry&) {});
    }
    EXPECT_TRUE(fs::exists(dir / "p.log.toml"));
    const std::string projectText = readFile(p);
    EXPECT_NE(projectText.find("log_file = \"p.log.toml\""), std::string::npos);
    EXPECT_EQ(projectText.find(kLogMarker), std::string::npos);
    EXPECT_EQ(steps(p), (std::vector<std::string>{"02", "03", "04"}));

    saveProject(p, loadProject(p));
    logged(p, "05", "video", {}, [](LogEntry&) {});
    EXPECT_EQ(steps(p), (std::vector<std::string>{"02", "03", "04", "05"}));
    EXPECT_EQ(logFileName(p), "p.log.toml");
}

// tests/test_projectfile.py::test_params_log_is_never_logged_into
TEST(StepLog, ParamsLogIsNeverLoggedInto) {
    TempDir dir;
    const fs::path path = dir / "MDATK_params_log.toml";
    writeParamsLog(path, Project{}, Table{{"timestamp", "t"}});
    EXPECT_TRUE(isParamsLog(path));
    const std::string before = readFile(path);
    logged(path, "render", "", {}, [](LogEntry&) {});
    EXPECT_EQ(readFile(path), before);
}

TEST(StepLog, NoProjectNoEntry) {
    TempDir dir;
    bool ran = false;
    logged({}, "render", "", {}, [&](LogEntry&) { ran = true; });
    EXPECT_TRUE(ran);
    logged(dir / "missing.mdatk.toml", "render", "", {}, [](LogEntry&) {});
    EXPECT_FALSE(fs::exists(dir / "missing.mdatk.toml"));
    EXPECT_THROW(appendLog(dir / "missing.mdatk.toml", LogEntry("x")), ProjectError);
}

TEST(StepLog, StatusAndSummaryOfFailures) {
    TempDir dir;
    const fs::path p = dir / "s.mdatk.toml";
    saveProject(p, Project{});
    EXPECT_THROW(logged(p, "a", "", {}, [](LogEntry&) { throw Cancelled("stopped by the user"); }), Cancelled);
    EXPECT_THROW(logged(p, "b", "", {}, [](LogEntry&) { throw std::invalid_argument(""); }), std::invalid_argument);
    EXPECT_THROW(logged(p, "c", "", {}, [](LogEntry& e) {
                     e.summary = "exit status 3";
                     throw std::runtime_error("ignored");
                 }),
                 std::runtime_error);
    EXPECT_THROW(logged(p, "d", "", {}, [](LogEntry&) { throw NamedError("SystemExit", "3"); }), NamedError);
    EXPECT_THROW(logged(p, "e", "", {}, [&](LogEntry&) { (void)fs::file_size(dir / "a.pdb"); }), fs::filesystem_error);
    const Array log = readLog(p);
    ASSERT_EQ(log.size(), 5u);
    EXPECT_EQ(log[0].find("status")->stringOr(), "cancelled");
    EXPECT_EQ(log[0].find("summary"), nullptr);
    EXPECT_EQ(log[1].find("summary")->stringOr(), "ValueError: ");
    EXPECT_EQ(log[2].find("summary")->stringOr(), "exit status 3");
    EXPECT_EQ(log[3].find("summary")->stringOr(), "SystemExit: 3");
#ifndef _WIN32
    EXPECT_EQ(log[4].find("summary")->stringOr(),
              "FileNotFoundError: [Errno 2] No such file or directory: '" + (dir / "a.pdb").string() + "'");
#endif
}

// tests/test_projectfile.py::test_concurrent_saves_and_log_appends_keep_every_entry
TEST(StepLog, ConcurrentSavesAndAppendsKeepEveryEntry) {
    TempDir dir;
    const fs::path p = dir / "c.mdatk.toml";
    saveProject(p, Project{});
    std::thread saver([&] {
        for (int i = 0; i < 40; ++i) {
            saveProject(p, Project{});
        }
    });
    for (int i = 0; i < 40; ++i) {
        appendLog(p, LogEntry(std::to_string(i)));
    }
    saver.join();
    std::vector<std::string> expected;
    for (int i = 0; i < 40; ++i) {
        expected.push_back(std::to_string(i));
    }
    EXPECT_EQ(steps(p), expected);
}

// tests/test_projectfile.py::test_undecodable_file_names_do_not_break_the_log
TEST(StepLog, UndecodableFileNamesDoNotBreakTheLog) {
    TempDir dir;
    const fs::path p = dir / "u.mdatk.toml";
    saveProject(p, Project{});
    LogEntry entry("pdb clean");
    entry.files = {"caf\xE9.pdb"};
    appendLog(p, entry);
    const Array log = readLog(p);
    ASSERT_EQ(log.size(), 1u);
    EXPECT_EQ(*log[0].find("files"), Value(Array{Value("caf\xEF\xBF\xBD.pdb")}));
}
