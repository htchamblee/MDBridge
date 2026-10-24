// Project files are read exactly as MDATK v0.0.4.1 reads them. The fixtures in data/projects were
// written by MDATK itself or by hand, and data/expected holds what MDATK's projectfile.read() and
// read_log() return for each (tests/parity/reference.py; see data/README.md).

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/project_file.hpp"

namespace fs = std::filesystem;
using namespace mdbridge::core;

namespace {

const fs::path kData = MDBRIDGE_TEST_DATA;
const fs::path kProjects = kData / "projects";
const fs::path kExpected = kData / "expected";

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

std::vector<std::string> fixtureNames() {
    std::vector<std::string> names;
    for (const auto& entry : fs::directory_iterator(kProjects)) {
        names.push_back(entry.path().filename().string());
    }
    std::sort(names.begin(), names.end());
    return names;
}

class MdatkFixture : public testing::TestWithParam<std::string> {};

}  // namespace

TEST_P(MdatkFixture, ReadsLikeMdatk) {
    const std::string name = GetParam();
    const fs::path path = kProjects / name;
    const fs::path error = kExpected / (name + ".error.txt");
    if (fs::exists(error)) {
        try {
            readProject(path);
            FAIL() << name << " was read, MDATK rejects it with: " << readFile(error);
        } catch (const ProjectError& e) {
            const std::string message = e.what();
            EXPECT_EQ(message.rfind(path.string() + ": ", 0), 0u) << message;  // names the file
            const std::string expected = readFile(error);
            if (expected.find("not an MDATK project") != std::string::npos) {  // the very same text
                EXPECT_EQ(message.substr(path.string().size()), expected.substr(name.size(), expected.size() - name.size() - 1));
            } else if (expected.find("not valid TOML") != std::string::npos) {  // the parser's words differ
                EXPECT_NE(message.find(": not valid TOML: "), std::string::npos) << message;
            } else {
                EXPECT_NE(message.find("JSON presets are no longer supported"), std::string::npos) << message;
            }
        }
        return;
    }
    const ProjectFile file = readProject(path);
    EXPECT_EQ(toJson(file.settings) + "\n", readFile(kExpected / (name + ".settings.json")));
    EXPECT_EQ(toJson(readLog(path)) + "\n", readFile(kExpected / (name + ".log.json")));
}

INSTANTIATE_TEST_SUITE_P(Fixtures, MdatkFixture, testing::ValuesIn(fixtureNames()),
                         [](const testing::TestParamInfo<std::string>& info) {
                             std::string id = info.param;
                             std::replace_if(id.begin(), id.end(), [](char c) { return !std::isalnum(static_cast<unsigned char>(c)); }, '_');
                             return id;
                         });

TEST(ProjectFile, EveryFixtureHasItsReference) {
    for (const std::string& name : fixtureNames()) {
        const bool error = fs::exists(kExpected / (name + ".error.txt"));
        const bool data = fs::exists(kExpected / (name + ".settings.json")) && fs::exists(kExpected / (name + ".log.json"));
        EXPECT_TRUE(error != data) << name << ": run tests/parity/reference.py (data/README.md)";
    }
}

TEST(ProjectFile, HeaderOfAProject) {
    const ProjectFile file = readProject(kProjects / "default.mdatk.toml");
    EXPECT_EQ(file.kind, FileKind::Project);
    std::vector<std::string> keys;
    for (const auto& entry : *file.header.get<Table>()) {
        keys.push_back(entry.first);
    }
    EXPECT_EQ(keys, (std::vector<std::string>{"format", "format_version", "mdatk_version", "saved", "log_file"}));
    EXPECT_EQ(file.header.find("mdatk_version")->stringOr(), "0.0.4.1");
    EXPECT_EQ(*file.header.find("format_version")->get<std::int64_t>(), kProjectFormatVersion);
    EXPECT_NE(file.header.find("saved")->get<DateTime>(), nullptr);
    EXPECT_EQ(file.settings.find("format"), nullptr);
    EXPECT_EQ(file.settings.path({"inputs", "step"})->get<std::int64_t>() != nullptr, true);
}

TEST(ProjectFile, ParameterLogOpensAsItsParameters) {
    const fs::path path = kProjects / "MDATK_params_log.toml";
    EXPECT_TRUE(isParamsLog(path));
    EXPECT_FALSE(isParamsLog(kProjects / "default.mdatk.toml"));
    EXPECT_FALSE(isParamsLog(kProjects / "no-such-file.toml"));
    const ProjectFile file = readProject(path);
    EXPECT_EQ(file.kind, FileKind::ParamsLog);
    EXPECT_EQ(file.header.find("format")->stringOr(), kParamsLogFormat);
    EXPECT_NE(file.header.find("run"), nullptr);
    EXPECT_EQ(file.settings.path({"inputs", "topology"})->stringOr(), "/data/run1/system.parm7");

    const ProjectFile empty = readProject(kProjects / "empty_params_log.toml");
    EXPECT_EQ(empty.kind, FileKind::ParamsLog);
    EXPECT_EQ(empty.settings.size(), 0u);
}

TEST(ProjectFile, SeparateStepLog) {
    const fs::path project = kProjects / "split.mdatk.toml";
    EXPECT_EQ(logFileName(project), "split.log.toml");
    EXPECT_EQ(fs::canonical(logPath(project)), fs::canonical(kProjects / "split.log.toml"));
    EXPECT_EQ(logFileName(kProjects / "full.mdatk.toml"), "");
    EXPECT_EQ(logPath(kProjects / "full.mdatk.toml"), kProjects / "full.mdatk.toml");
    const Array log = readLog(project);
    std::vector<std::string> steps;
    for (const Value& entry : log) {
        steps.push_back(entry.find("step")->stringOr());
    }
    EXPECT_EQ(steps, (std::vector<std::string>{"pdb clean", "04", "05"}));
}

TEST(ProjectFile, MissingLogGivesAWarningAndNoEntries) {
    std::string warning;
    EXPECT_TRUE(readLog(kProjects / "missing_log.mdatk.toml", &warning).empty());
    EXPECT_NE(warning.find("could not read the step log"), std::string::npos) << warning;
    EXPECT_NE(warning.find("missing.log.toml"), std::string::npos) << warning;
}

TEST(ProjectFile, MissingFile) {
    EXPECT_THROW(readProject(kProjects / "no-such-file.mdatk.toml"), ProjectError);
    EXPECT_THROW(readProject(kProjects), ProjectError);  // a folder
}

TEST(ProjectFile, ParseTomlKeepsTheOrderOfTheText) {
    const Value data = parseToml("z = 1\na = 2\n[m]\ny = 3\nb = 4\n[c]\nx = { q = 1, p = 2 }\n", "text");
    EXPECT_EQ(toJson(data, -1), "{\"z\": 1, \"a\": 2, \"m\": {\"y\": 3, \"b\": 4}, \"c\": {\"x\": {\"q\": 1, \"p\": 2}}}");
    try {
        parseToml("a = \n", "notes.toml");
        FAIL() << "invalid TOML was accepted";
    } catch (const ProjectError& e) {
        EXPECT_EQ(std::string(e.what()).rfind("notes.toml: not valid TOML: ", 0), 0u) << e.what();
        EXPECT_NE(std::string(e.what()).find("(at line 1, column"), std::string::npos) << e.what();
    }
}
