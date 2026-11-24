// The MDATK commands of the command line (init, info, export-pdb, density info), run in-process on
// the inputs of tests/data/cli/cli_cases.json and compared with what MDATK v0.0.4.1 printed and
// wrote there (tests/parity/cli_parity.py --write; normalised the same way: {work}, {time},
// {stamp}).
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include "commands.hpp"
#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/text.hpp"
#include "test_support.hpp"

namespace mdbridge {
namespace {

namespace fs = std::filesystem;
using core::Value;
using test::kData;

const Value& field(const Value& v, std::string_view key) {
    const Value* f = v.find(key);
    if (f == nullptr) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *f;
}

std::string text(const Value& v, std::string_view key) { return field(v, key).stringOr(); }

std::string replaceAll(std::string s, const std::string& from, const std::string& to) {
    for (std::size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size()) {
        s.replace(pos, from.size(), to);
    }
    return s;
}

struct CommandRun {
    int exit = 0;
    std::string out;
    std::string err;
};

// The command line in `folder` (the working directory meanwhile).
CommandRun runIn(const fs::path& folder, const std::vector<std::string>& args) {
    std::vector<const char*> argv{"mdbridge"};
    for (const std::string& a : args) {
        argv.push_back(a.c_str());
    }
    std::ostringstream out;
    std::ostringstream err;
    app::Console console{out, err, false};
    const fs::path previous = fs::current_path();
    fs::current_path(folder);
    CommandRun run;
    try {
        run.exit = app::runCoreCommandLine(static_cast<int>(argv.size()), argv.data(), console);
    } catch (...) {
        fs::current_path(previous);
        throw;
    }
    fs::current_path(previous);
    run.out = out.str();
    run.err = err.str();
    return run;
}

class Normaliser {
public:
    explicit Normaliser(const fs::path& work) : m_work(core::pathToUtf8(fs::canonical(work))) {}

    std::string operator()(std::string s) const {
        s = replaceAll(std::move(s), m_work, "{work}");
        return std::regex_replace(s, m_stamp, "{stamp}");
    }
    std::string file(std::string s) const {
#ifdef _WIN32
        s = replaceAll(std::move(s), "\r\n", "\n");
#endif
        return std::regex_replace((*this)(std::move(s)), m_saved, "$1saved = {time}");
    }

private:
    std::string m_work;
    std::regex m_stamp{R"(\d{4}-\d{2}-\d{2}_\d{6})"};
    std::regex m_saved{R"((^|\n)saved = \d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?(?=\n|$))"};
};

bool validUtf8(const std::string& s) {
    return core::encodeUtf8(core::decodeUtf8(s)) == s;
}

TEST(Commands, MatchMdatk) {
    const Value data = test::readJson(kData / "cli" / "cli_cases.json");
    std::size_t run = 0;
    for (const Value& c : *field(data, "cases").get<core::Array>()) {
        const std::string name = text(c, "name");
        SCOPED_TRACE(name);
#ifdef _WIN32
        // absolute paths and OS error texts as on Linux, where parity is defined (posix_only)
        if (*field(c, "posix_only").get<bool>()) {
            continue;
        }
#endif
        test::TempDir work;
        for (const auto& [dest, source] : *field(c, "inputs").get<core::Table>()) {
            const fs::path target = work.path() / core::pathFromUtf8(dest);
            fs::create_directories(target.parent_path());
            fs::copy_file(kData / core::pathFromUtf8(source.stringOr()), target);
        }
        for (const auto& [dest, content] : *field(c, "setup").get<core::Table>()) {
            test::writeFile(work.path() / core::pathFromUtf8(dest), content.stringOr());
        }
        std::vector<std::string> args;
        for (const Value& a : *field(c, "argv").get<core::Array>()) {
            args.push_back(a.stringOr());
        }
        const CommandRun result = runIn(work.path(), args);
        ++run;
        const Normaliser normalise(work.path());
        if (field(c, "exit_only").get<bool>() != nullptr && *field(c, "exit_only").get<bool>()) {
            EXPECT_EQ(result.exit, *field(c, "exit").get<std::int64_t>());
            continue;
        }
        EXPECT_EQ(result.exit, *field(c, "exit").get<std::int64_t>()) << result.err;
        EXPECT_EQ(normalise(result.out), text(c, "stdout"));
        if (field(c, "stderr").isNull()) {
            // MDATK crashed with a traceback (mdatk_crash): MDBridge reports the error instead
            EXPECT_TRUE(result.err.starts_with("error: ")) << text(c, "mdatk_crash") << " / " << result.err;
        } else {
            EXPECT_EQ(normalise(result.err), text(c, "stderr"));
        }

        // the files left in the folder: everything but the unchanged inputs
        core::Table expectedFiles = *field(c, "written").get<core::Table>();
#ifdef _WIN32
        // MDATK's project lock file (.NAME.lock) needs fcntl, which Windows lacks: there neither MDATK
        // nor MDBridge creates one (projectfile._locked, detail::ProjectLock)
        static const std::regex projectLock{R"((^|/)\.[^/]+\.lock$)"};
        std::erase_if(expectedFiles,
                      [](const auto& entry) { return std::regex_search(entry.first, projectLock); });
#endif
        std::vector<std::string> seen;
        for (const auto& entry : fs::recursive_directory_iterator(work.path())) {
            if (!entry.is_regular_file()) {
                continue;
            }
            std::string rel = core::pathToUtf8(entry.path().lexically_relative(work.path()));
            rel = replaceAll(rel, "\\", "/");
            const std::string content = test::readFile(entry.path());
            if (field(c, "inputs").find(rel) != nullptr) {
                continue;
            }
            if (const Value* given = field(c, "setup").find(rel); given != nullptr && given->stringOr() == content) {
                continue;
            }
            rel = normalise(rel);
            seen.push_back(rel);
            const Value* expected = nullptr;
            for (const auto& [key, value] : expectedFiles) {
                if (key == rel) {
                    expected = &value;
                }
            }
            if (expected == nullptr) {
                ADD_FAILURE() << "unexpected file " << rel;
                continue;
            }
            if (validUtf8(content)) {
                const std::string got = normalise.file(content);
                if (const Value* t = expected->find("text")) {
                    EXPECT_EQ(got, t->stringOr()) << rel;
                } else {
                    EXPECT_EQ(core::sha256Hex(got), text(*expected, "sha256")) << rel;
                    EXPECT_EQ(static_cast<std::int64_t>(got.size()), *field(*expected, "size").get<std::int64_t>())
                        << rel;
                }
            } else {
                EXPECT_EQ(core::sha256Hex(content), text(*expected, "sha256")) << rel;
            }
        }
        for (const auto& [key, value] : expectedFiles) {
            if (std::find(seen.begin(), seen.end(), key) == seen.end()) {
                ADD_FAILURE() << "missing file " << key;
            }
        }
    }
    EXPECT_GT(run, 0u);
}

// ---- pieces ------------------------------------------------------------------------------------

TEST(Commands, OverridesFollowMdatk) {
    core::Project p;
    app::applyOverride(p, "analyses.contacts.cutoff=5");
    EXPECT_EQ(p.analyses.contacts.cutoff, 5.0);
    app::applyOverride(p, " groups.-1.name =Last");  // the key is stripped, indices may be negative
    EXPECT_EQ(p.groups.back().name, "Last");
    app::applyOverride(p, "notes=[1, 2]");  // a text setting takes the text
    EXPECT_EQ(p.notes, "[1, 2]");
    app::applyOverride(p, "notes=\"quoted\"");
    EXPECT_EQ(p.notes, "quoted");
    app::applyOverride(p, "inputs.stop=7");
    EXPECT_EQ(p.inputs.stop, 7);
    app::applyOverride(p, "inputs.stop=null");
    EXPECT_EQ(p.inputs.stop, std::nullopt);
    app::applyOverride(p, "groups=[{\"name\": \"X\"}]");
    ASSERT_EQ(p.groups.size(), 1u);
    EXPECT_EQ(p.groups[0].residues, "");

    const auto message = [&p](const std::string& assignment) {
        try {
            app::applyOverride(p, assignment);
        } catch (const core::ValueError& e) {
            return std::string(e.what());
        }
        return std::string("no error");
    };
    EXPECT_EQ(message("x"), "--set expects key=value, got 'x'");
    EXPECT_EQ(message(" . =1"), "empty key in ' . =1'");
    EXPECT_EQ(message("inputs.bogus=1"), "unknown parameter 'inputs.bogus'");
    EXPECT_EQ(message("layout.panels.3.options.mode=x"), "unknown parameter 'layout.panels.3.options.mode'");
    EXPECT_EQ(message("validate.x=1"), "unknown parameter 'validate.x'");
    EXPECT_EQ(message("cleaning.enabled=1"), "cleaning.enabled expects true/false");
    EXPECT_EQ(message("groups.1.name=x"), "groups.1.name: no list element '1'");
    EXPECT_EQ(message("groups.99999999999999999999999.name=x"),
              "groups.99999999999999999999999.name: no list element '99999999999999999999999'");
    EXPECT_EQ(message("groups.--1.name=x"), "invalid literal for int() with base 10: '--1'");
    EXPECT_EQ(message("analyses.contacts.cutoff=far"),
              "project.analyses.contacts.cutoff: expected a number, got 'far'");
}

TEST(Commands, JsonIsAscii) {
    EXPECT_EQ(app::jsonForConsole(Value(core::Table{{"k\xc3\xa9", "\xf0\x9f\x98\x80\x7f"}})),
              "{\n  \"k\\u00e9\": \"\\ud83d\\ude00\\u007f\"\n}");
    EXPECT_EQ(app::jsonForConsole(Value(std::string("a\xff"))), "\"a\\udcff\"");
}

TEST(Commands, OnlyTheseCommandsAreOffered) {
    CLI::App root("test", "mdbridge");
    app::CoreCommands commands(root);
    std::vector<std::string> names;
    for (const CLI::App* sub : root.get_subcommands({})) {
        names.push_back(sub->get_name());
    }
    EXPECT_EQ(names, (std::vector<std::string>{"init", "info", "export-pdb", "density"}));
    const std::string help = root.help();
    EXPECT_EQ(help.find("analyze"), std::string::npos);
    EXPECT_EQ(help.find("render"), std::string::npos);
}

TEST(Commands, UsageErrorsExitWith2) {
    test::TempDir work;
    EXPECT_EQ(runIn(work.path(), {"init"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"density"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"analyze"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"--help"}).exit, 0);
}

TEST(Commands, InitWritesTheProjectFileWithTheClock) {
    test::TempDir work;
    const core::Clock previous = core::setClock({[] {
                                                     core::DateTime t;
                                                     t.year = 2026, t.month = 10, t.day = 15, t.hour = 14, t.minute = 32;
                                                     return t;
                                                 },
                                                 {}});
    const CommandRun run = runIn(work.path(), {"init-params", "a.mdatk.toml", "-s", "notes=hello"});
    core::setClock(previous);
    EXPECT_EQ(run.exit, 0);
    EXPECT_EQ(run.out, "Wrote a project file with default settings to a.mdatk.toml\n");
    const std::string written = test::readFile(work.path() / "a.mdatk.toml");
    EXPECT_NE(written.find("saved = 2026-10-15T14:32:00"), std::string::npos);
    EXPECT_NE(written.find("notes = \"hello\""), std::string::npos);
    EXPECT_EQ(runIn(work.path(), {"init", "a.mdatk.toml"}).err, "a.mdatk.toml exists (use --overwrite)\n");
}

}  // namespace
}  // namespace mdbridge
