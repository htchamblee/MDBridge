// Saving and loading projects, parameter logs and the recent-projects list, byte for byte as
// MDATK v0.0.4.1 writes them: the golden files of tests/parity/make_project_file_golden.py are
// replayed step by step with the same clocks, and the tests of MDATK's test_projectfile.py and
// test_config_selections.py that save or load files are ported.

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/settings.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/core/step_types.hpp"
#include "mdbridge/core/toml_writer.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

const fs::path kGolden = kData / "project-file" / "golden";
const fs::path kWorkflowGolden = kData / "project-file" / "workflow";

// start / end readings of time.monotonic() for each logged() block of make_project_file_golden.py
const std::vector<double> kMonotonic{100.0, 101.23456, 102.0, 102.004, 103.0, 105.555, 106.0, 106.125, 107.0,
                                     107.135, 108.0,   115.0,     116.0, 116.5,   117.0, 117.995, 118.0, 118.005,
                                     119.0,   119.015, 120.0,     120.25, 121.0,  121.5, 122.0,   122.75, 123.0,
                                     123.1};

void setEnvironment(const char* name, const std::string& value) {
#ifdef _WIN32
    _putenv_s(name, value.c_str());
#else
    if (value.empty()) {
        unsetenv(name);
    } else {
        setenv(name, value.c_str(), 1);
    }
#endif
}

std::string environment(const char* name) {
    const char* value = std::getenv(name);
    return value ? value : "";
}

class ScopedEnvironment {
public:
    ScopedEnvironment(const char* name, const std::string& value) : m_name(name), m_old(environment(name)) {
        setEnvironment(name, value);
    }
    ~ScopedEnvironment() { setEnvironment(m_name, m_old); }

private:
    const char* m_name;
    std::string m_old;
};

// A path as UTF-8 text in generic form ('/' separators): how loaded projects write absolute paths
// on Windows (MDBridge's Windows rules, ntpath.hpp); the plain text on POSIX.
std::string genericUtf8(const fs::path& path) {
    const std::u8string text = path.generic_u8string();
    return std::string(text.begin(), text.end());
}

Project fromFile(const fs::path& path) { return projectFromPlain(readProject(path).settings); }

std::string replaceAll(std::string text, const std::string& from, const std::string& to) {
    for (std::size_t at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size())) {
        text.replace(at, from.size(), to);
    }
    return text;
}

void expectSameFile(const fs::path& written, const std::string& golden) {
    EXPECT_EQ(readFile(written), nativeLines(readFile(kGolden / golden))) << golden;
}

}  // namespace

// The whole sequence of make_project_file_golden.py: the default project, a rich project and its
// step log (every failure kind), a truncated output, the log moving to a side file, a parameter
// log and the recent-projects list. Every file must match MDATK's byte for byte.
TEST(ProjectFileGolden, ReplaysMdatkSequence) {
    FakeClock clock(localDateTime(2026, 10, 15, 14, 32, 5), kMonotonic);
    TempDir dir;

    // 1. default project
    saveProject(dir / "default/default.mdatk.toml", Project{});
    expectSameFile(dir / "default/default.mdatk.toml", "default.mdatk.toml");
#ifndef _WIN32
    EXPECT_TRUE(fs::exists(dir / "default/.default.mdatk.toml.lock"));
#endif
    EXPECT_FALSE(fs::exists(dir / "default/default.mdatk.toml.tmp"));

    // 2. rich project and its step log
    const Project rich = fromFile(kGolden / "rich.mdatk.toml");
    const fs::path path = dir / "rich/rich.mdatk.toml";
    saveProject(path, rich);
    expectSameFile(path, "rich.mdatk.toml");

    logged(path, "02", "tleap build", {"case-2023-amber", "key1"}, [](LogEntry& entry) {
        entry.command = "bash apt1_02_tleap.sh";
        entry.folder = "/home/user/projects/apt1_workflow/02_amber-build";
        entry.files = {"apt1_02_tleap_20261015-1432.parm7", "apt1_02_tleap_20261015-1432.rst7"};
        entry.summary = "exit status 0";
        entry.output = "Welcome to LEaP!\nLoading parameters\n\tdone\n";
    });
    EXPECT_THROW(logged(path, "03", "md", {}, [](LogEntry&) { throw std::runtime_error("pmemd crashed"); }),
                 std::runtime_error);
    EXPECT_THROW(logged(path, "04", "slow run", {}, [](LogEntry&) { throw Cancelled("stopped by the user"); }),
                 Cancelled);
    EXPECT_THROW(logged(path, "render", "Analyses and outputs", {}, [](LogEntry&) { throw Cancelled(); }),
                 Cancelled);  // KeyboardInterrupt in MDATK
    EXPECT_THROW(logged(path, "pdb clean", "", {}, [](LogEntry&) { throw std::invalid_argument(""); }),
                 std::invalid_argument);
    EXPECT_THROW(logged(path, "05", "summary set before the failure", {},
                        [](LogEntry& entry) {
                            entry.summary = "exit status 3";
                            throw std::runtime_error("ignored: the summary is already set");
                        }),
                 std::runtime_error);
    logged(path, "pdb clean", "PDB clean-up of caf\xE9.pdb", {}, [](LogEntry& entry) {
        entry.files = {"caf\xE9.pdb", "clean/caf\xE9_clean.pdb"};  // undecodable file names
        entry.summary = "renamed 3 residues\nrenumbered chain B\n";
    });
    expectSameFile(path, "rich.with-log.mdatk.toml");

    Project q = fromFile(path);  // no path resolution
    q.notes = "edited notes";
    saveProject(path, q);  // keeps the log
    expectSameFile(path, "rich.resaved.mdatk.toml");

    // 3. output truncated to its last 20000 characters
    const fs::path t = dir / "truncation/t.mdatk.toml";
    saveProject(t, Project{});
    logged(t, "06", "long output", {}, [](LogEntry& entry) {
        for (int k = 0; k < 2000; ++k) {
            char line[16];
            std::snprintf(line, sizeof line, "line %05d\n", k);
            entry.output += line;
        }
    });
    expectSameFile(t, "truncation.mdatk.toml");

    // 4. the log moves to p.log.toml once the project file is too large
    const fs::path p = dir / "split/p.mdatk.toml";
    saveProject(p, Project{});
    logged(p, "02", "tleap build", {}, [](LogEntry&) {});
    logged(p, "03", "md", {}, [](LogEntry&) {});
    expectSameFile(p, "split.before.mdatk.toml");
    setLogSplitBytes(1000);
    logged(p, "04", "analysis", {}, [](LogEntry&) {});
    setLogSplitBytes(kDefaultLogSplitBytes);
    expectSameFile(p, "split.after.mdatk.toml");
    expectSameFile(dir / "split/p.log.toml", "split.after.log.toml");
    saveProject(p, fromFile(p));
    logged(p, "05", "video", {}, [](LogEntry&) {});
    expectSameFile(p, "split.resaved.mdatk.toml");
    expectSameFile(dir / "split/p.log.toml", "split.resaved.log.toml");

    const fs::path p2 = dir / "split2/q.mdatk.toml";  // over the limit before it has a log
    saveProject(p2, Project{});
    setLogSplitBytes(1000);
    logged(p2, "01", "first entry of a big project", {}, [](LogEntry&) {});
    setLogSplitBytes(kDefaultLogSplitBytes);
    expectSameFile(p2, "split-nolog.after.mdatk.toml");
    expectSameFile(dir / "split2/q.log.toml", "split-nolog.after.log.toml");

    // 5. parameter log of a run: never logged into, never saved over
    fs::create_directories(dir / "params");
    const fs::path params = dir / "params/MDATK_params_log.toml";
    const Value run = Table{
        {"project", "Molecular Dynamic Analysis Tool Kit (MDATK)"},
        {"environment", Table{{"mdatk", "0.0.4.1"},
                              {"python", "3.13.17"},
                              {"MDAnalysis", "2.9.0"},
                              {"numpy", "2.3.1"},
                              {"platform", "Linux-6.18-x86_64"},
                              {"cpptraj", Value()}}},
        {"inputs", Table{{"topology", Value()}, {"trajectories", Value()}, {"cleaning", Value()}}},
        {"cache_dir", "/home/user/projects/results/MDATK_cache"},
        {"analysis_keys", Table{{"rmsd", "0123456789abcdef"}, {"contacts", "fedcba9876543210"}}},
        {"timestamp", "2026-10-15_143205"},
        {"outputs", Table{{"video", ""}, {"charts", Array{}}, {"pdbs", Array{}}, {"params_log", Value()}}},
        {"citation_keys", Array{"mdanalysis-2016", "mdanalysis-2011", "numpy-2020", "matplotlib-2007", "mdatk"}},
        {"citations",
         Array{"Gowers, R. J. et al. MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics "
               "Simulations. Proc. 15th Python in Science Conf. (2016).",
               "Harris, C. R. et al. Array programming with NumPy. Nature 585, 357-362 (2020)."}},
    };
    writeParamsLog(params, rich, run);
    expectSameFile(params, "params_log.toml");
    const std::string before = readFile(params);
    logged(params, "render", "", {}, [](LogEntry&) {});
    EXPECT_EQ(readFile(params), before);
    try {
        saveProject(params, Project{});
        ADD_FAILURE() << "saved over a parameter log";
    } catch (const ProjectError& error) {
        EXPECT_EQ(std::string(error.what()), params.string() +
                                                 " is the parameter log of a run; save the project under another name");
    }
    EXPECT_EQ(clock.monotonicReadings(), kMonotonic.size());

    // 6. recent projects
    const ScopedEnvironment config("XDG_CONFIG_HOME", (dir / "config").string());
    const ScopedEnvironment appData("APPDATA", (dir / "config").string());  // Windows (D16)
    std::vector<fs::path> paths;
    for (std::size_t k = 0; k < kRecentMax + 2; ++k) {
        const fs::path project = dir / ("recent/p" + std::to_string(k) + ".mdatk.toml");
        saveProject(project, Project{});
        addRecent(project);
        paths.push_back(project);
    }
#ifndef _WIN32
    EXPECT_EQ(replaceAll(readFile(dir / "config/mdatk/settings.toml"), dir.path().string(), "/home/user/projects"),
              readFile(kGolden / "settings.toml"));
#endif
    ASSERT_FALSE(recentProjects().empty());
    EXPECT_EQ(recentProjects()[0], paths.back());
}

// The planner goldens of the step-type spec: saving what was loaded gives the same text, so every
// settings comment, the record quirk and the unnormalised values come out as MDATK writes them.
TEST(ProjectFileGolden, PlannerProjectsRewriteIdentically) {
    for (const char* name : {"golden_all_tools.mdatk.toml", "golden_edge_cases.mdatk.toml",
                             "golden_test_every_entry.mdatk.toml", "golden_test_saved_in_project.mdatk.toml",
                             "golden_nondefault.mdatk.toml", "golden_written_stamp.mdatk.toml"}) {
        FakeClock clock(localDateTime(2026, 10, 15, 14, 32, 0));
        EXPECT_EQ(projectText(fromFile(kWorkflowGolden / name)), readFile(kWorkflowGolden / name)) << name;
    }
}

// tests/test_projectfile.py::test_every_entry_is_commented, building the project as MDATK does
TEST(ProjectIo, EveryEntryIsCommented) {
    Project p;
    p.notes = "Aim: compare\nthree \"quoted\" runs";
    PlanNode node;
    node.id = "n1";
    node.tool = "amber-build";
    node.step = 1;
    node.settings = stepType("amber-build").normalise(Table{{"mg_count", 4}});
    node.files = Table{{"structure", "model.pdb"}};
    p.workflow.nodes.push_back(node);
    p.layout.panels[1].options = Table{{"smooth", 5}, {"progressive", true}, {"custom", 1}};

    FakeClock clock(localDateTime(2026, 10, 15, 14, 32, 0));
    const std::string text = projectText(p);
    EXPECT_EQ(text, readFile(kWorkflowGolden / "golden_test_every_entry.mdatk.toml"));

    // every "key = value" line outside multi-line strings and arrays has a comment line above it
    std::vector<std::string> lines;
    for (std::size_t start = 0; start < text.size();) {
        const std::size_t end = text.find('\n', start);
        lines.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    bool inString = false;
    int depth = 0;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string& line = lines[i];
        if (!inString && depth == 0 && !line.empty() && line[0] != '#' && line[0] != '[') {
            const std::size_t eq = line.find(" = ");
            ASSERT_NE(eq, std::string::npos) << line;
            EXPECT_TRUE(i > 0 && !lines[i - 1].empty() && lines[i - 1][0] == '#') << line;
        }
        if (!inString && depth == 0 && line.rfind("[", 0) == 0 && line.rfind("[[", 0) != 0) {
            EXPECT_TRUE(i > 0 && !lines[i - 1].empty() && lines[i - 1][0] == '#') << line;
        }
        for (std::size_t at = line.find("\"\"\""); at != std::string::npos; at = line.find("\"\"\"", at + 3)) {
            inString = !inString;
        }
        if (!inString && line.size() > 2 && line.compare(line.size() - 2, 2, " [") == 0) {
            ++depth;
        } else if (!inString && depth > 0 && line == "]") {
            --depth;
        }
    }
    EXPECT_EQ(text.find("\n# key ="), std::string::npos);

    TempDir dir;
    saveProject(dir / "e.mdatk.toml", p);
    EXPECT_EQ(loadProject(dir / "e.mdatk.toml").notes, p.notes);
}

// tests/test_projectfile.py::test_optional_values_and_round_trip
TEST(ProjectIo, OptionalValuesAndRoundTrip) {
    TempDir dir;
    Project p;
    p.inputs.stop = std::nullopt;
    p.inputs.dt_ps = 2.5;
    p.analyses.energy.target_value = std::nullopt;
    PanelSpec panel;
    panel.type = "contact_map";
    panel.series = {"contacts/map"};
    panel.title = "Map";
    panel.options = Table{{"transpose", true}};
    p.layout.panels.push_back(panel);
    const fs::path path = dir / "o.mdatk.toml";
    saveProject(path, p);
    EXPECT_NE(readFile(path).find(nativeLines("\n# stop =   (not set)\n")), std::string::npos);
    const Project q = loadProject(path);
    EXPECT_FALSE(q.inputs.stop.has_value());
    EXPECT_EQ(q.inputs.dt_ps, 2.5);
    EXPECT_EQ(toPlain(loadProject(path)), toPlain(q));
    saveProject(path, q);
    EXPECT_EQ(toPlain(loadProject(path)), toPlain(q));
}

// tests/test_projectfile.py::test_examples_load_and_rewrite_identically (the examples are empty)
TEST(ProjectIo, EmptyFileLoadsAsTheDefaults) {
    TempDir dir;
    writeFile(dir / "empty.mdatk.toml", "");
    const Project p = loadProject(dir / "empty.mdatk.toml");
    saveProject(dir / "again.mdatk.toml", p);
    EXPECT_EQ(toPlain(loadProject(dir / "again.mdatk.toml")), toPlain(p));
    EXPECT_EQ(fromFile(dir / "empty.mdatk.toml"), Project{});
}

// tests/test_projectfile.py::test_params_log_is_never_overwritten_by_save
TEST(ProjectIo, ParamsLogIsNeverOverwrittenBySave) {
    TempDir dir;
    const fs::path path = dir / "MDATK_params_log.toml";
    writeParamsLog(path, Project{}, Table{{"timestamp", "t"}});
    EXPECT_THROW(saveProject(path, Project{}), ProjectError);
    EXPECT_TRUE(isParamsLog(path));
}

// tests/test_config_selections.py::test_project_roundtrip
TEST(ProjectIo, ProjectRoundTrip) {
    TempDir dir;
    Project p;
    p.inputs.topology = genericUtf8(dir / "sys.parm7");
    p.inputs.trajectories = {genericUtf8(dir / "a.nc"), genericUtf8(dir / "b.nc")};
    p.output.out_dir = genericUtf8(dir / "out");
    p.analyses.contacts.cutoff = 5.25;
    PanelSpec panel;
    panel.type = "contact_map";
    panel.series = {"contacts/map"};
    panel.title = "Map";
    panel.options = Table{{"mode", "frequency"}};
    p.layout.panels.push_back(panel);
    saveProject(dir / "p.mdatk.toml", p);
    const Project q = loadProject(dir / "p.mdatk.toml");
    EXPECT_EQ(toPlain(q), toPlain(p));
    EXPECT_EQ(stableHash(toPlain(q)), stableHash(toPlain(p)));
}

// tests/test_config_selections.py::test_relative_paths_and_aliases
TEST(ProjectIo, RelativePathsAndAliases) {
    TempDir dir;
    writeFile(dir / "p.mdatk.toml",
              "not_a_key = 1\n[input]\ntopology = \"top.parm7\"\ntrajectory = \"traj.nc\"\n[output]\nout_dir = \"out\"\n");
    std::vector<std::string> warnings;
    const Project p = loadProject(dir / "p.mdatk.toml", &warnings);
    const std::string base = genericUtf8(fs::absolute(dir.path()).lexically_normal());
    EXPECT_EQ(p.inputs.topology, base + "/top.parm7");
    EXPECT_EQ(p.inputs.trajectories, std::vector<std::string>{base + "/traj.nc"});
    EXPECT_EQ(p.output.out_dir, base + "/out");
    ASSERT_EQ(warnings.size(), 1u);
    // PRJ-09 (D9): kept and saved back, one warning (MDATK: "project: ignoring unknown key 'not_a_key'")
    EXPECT_EQ(warnings[0], "project: 1 unknown key kept and saved back unchanged (from a newer MDATK, or "
                           "misspelt?): project.not_a_key");
    EXPECT_EQ(p.unknownKeys, (Table{{"not_a_key", 1}}));
}

// tests/test_config_selections.py::test_params_log_can_be_replayed
TEST(ProjectIo, ParamsLogCanBeReplayed) {
    TempDir dir;
    Project p;
    p.analyses.hbonds.angle_cutoff = 150.0;
    const fs::path path = dir / "MDATK_params_log.toml";
    writeParamsLog(path, p, Table{{"timestamp", "now"}, {"citations", Array{"x"}}});
    EXPECT_EQ(loadProject(path).analyses.hbonds.angle_cutoff, 150.0);
}

// tests/test_config_selections.py::test_json_presets_are_rejected
TEST(ProjectIo, JsonPresetsAreRejected) {
    TempDir dir;
    for (const char* name : {"old_preset.json", "old_preset.JSON"}) {
        try {
            loadProject(dir / name);
            ADD_FAILURE() << "loaded " << name;
        } catch (const ProjectError& error) {
            EXPECT_NE(std::string(error.what()).find("TOML"), std::string::npos);
        }
        try {
            saveProject(dir / name, Project{});
            ADD_FAILURE() << "saved " << name;
        } catch (const ProjectError& error) {
            EXPECT_NE(std::string(error.what()).find("TOML"), std::string::npos);
        }
    }
}

// tests/test_workflows.py::test_workflow_saved_in_project
TEST(ProjectIo, WorkflowSavedInProject) {
    TempDir dir;
    const fs::path path = dir / "w.mdatk.toml";
    saveProject(path, fromFile(kWorkflowGolden / "golden_test_saved_in_project.mdatk.toml"));
    const Value raw = parseToml(readFile(path), "w.mdatk.toml");
    const Value* nodes = raw.path({"workflow", "nodes"});
    ASSERT_TRUE(nodes && nodes->isArray());
    EXPECT_EQ(*(*nodes->get<Array>())[0].path({"settings", "mg_count"}), Value(4));

    const Project q = loadProject(path);
    const std::string base = genericUtf8(fs::absolute(dir.path()).lexically_normal());
    ASSERT_GE(q.workflow.nodes.size(), 3u);
    EXPECT_EQ(*q.workflow.nodes[0].files.front().second.get<std::string>(), base + "/inputs/model.pdb");
    EXPECT_EQ(q.workflow.nodes[1].settings.front().second, Value(base + "/maps/map.mrc"));
    EXPECT_EQ(q.workflow.edges.size(), 1u);
    ASSERT_EQ(q.workflow.scratch.size(), 1u);
    const Project original = fromFile(kWorkflowGolden / "golden_test_saved_in_project.mdatk.toml");
    for (std::size_t i = 0; i < q.workflow.nodes.size(); ++i) {
        EXPECT_EQ(q.workflow.nodes[i].id, original.workflow.nodes[i].id);
        EXPECT_EQ(q.workflow.nodes[i].tool, original.workflow.nodes[i].tool);
    }
    bool tube = false;
    for (const auto& [key, value] : q.workflow.scratch[0].settings) {
        tube = tube || (key == "rna_style" && value == Value("Tube"));
    }
    EXPECT_TRUE(tube);
}

TEST(ProjectIo, InvalidUtf8InAnExistingProjectIsAnError) {
    TempDir dir;
    const fs::path path = dir / "bad.mdatk.toml";
    writeFile(path, "notes = \"caf\xE9\"\n");
    EXPECT_THROW(saveProject(path, Project{}), ProjectError);  // as MDATK: the old log cannot be kept
}

// tests/test_projectfile.py::test_recent_projects
TEST(ProjectIo, RecentProjects) {
    TempDir dir;
    const ScopedEnvironment config("XDG_CONFIG_HOME", (dir / "config").string());
    const ScopedEnvironment appData("APPDATA", (dir / "config").string());
#ifdef _WIN32
    EXPECT_EQ(userSettingsPath(), dir / "config" / "MDBridge" / "settings.toml");  // D16: not shared
#else
    EXPECT_EQ(userSettingsPath(), dir / "config" / "mdatk" / "settings.toml");  // shared with MDATK
#endif
    EXPECT_TRUE(recentProjects().empty());
    std::vector<fs::path> paths;
    for (int k = 0; k < 12; ++k) {
        const fs::path project = dir / ("p" + std::to_string(k) + ".mdatk.toml");
        saveProject(project, Project{});
        addRecent(project);
        paths.push_back(project);
    }
    std::vector<fs::path> recent = recentProjects();
    ASSERT_EQ(recent.size(), 10u);
    EXPECT_EQ(recent[0], paths[11]);
    fs::remove(paths[11]);
    recent = recentProjects();
    EXPECT_EQ(recent[0], paths[10]);
}

// PRJ-12 (AP-17, MDATK 0.0.14): the recent list keeps every other setting of the shared file
// (RNAdo's cat_proof and whatever else MDATK keeps there) and is written atomically. The text is
// MDATK 0.1.1's add_recent for the same file.
TEST(ProjectIo, RecentListKeepsOtherSettings) {
    TempDir dir;
    const ScopedEnvironment config("XDG_CONFIG_HOME", (dir / "config").string());
    const ScopedEnvironment appData("APPDATA", (dir / "config").string());
    const fs::path settings = userSettingsPath();
    fs::create_directories(settings.parent_path());
    writeFile(settings, "# x\ncat_proof = true\nrecent = [\"/gone.mdatk.toml\"]\n[window]\nwidth = 3\n");
    const fs::path project = dir / "p.mdatk.toml";
    saveProject(project, Project{});
    addRecent(project);
    const std::string written = replaceAll(readFile(settings), "\r\n", "\n");
    EXPECT_EQ(written, "# Molecular Dynamic Analysis Tool Kit (MDATK) user settings.\n\n"
                       "# RNAdo's Cat-proof mode: the keyboard is ignored (untick it with the mouse).\n"
                       "cat_proof = true\n"
                       "# Recently used project files, newest first.\n"
                       "recent = " + tomlValue(Array{Value(pathToUtf8(recentProjects().at(0)))}) + "\n\n"
                       "[window]\nwidth = 3\n");
    EXPECT_EQ(recentProjects(), std::vector<fs::path>{project});
    EXPECT_FALSE(fs::exists(settings.string() + ".tmp"));

}

// AP-17 (MDATK 0.1.2, _settings_to_change): a settings file that is no TOML is kept as
// settings.toml.unreadable, with a warning, and a new one starts; one that cannot be read at all is
// left alone and nothing is written.
TEST(ProjectIo, UnreadableSettingsFileIsKept) {
    TempDir dir;
    const ScopedEnvironment config("XDG_CONFIG_HOME", (dir / "config").string());
    const ScopedEnvironment appData("APPDATA", (dir / "config").string());
    const fs::path settings = userSettingsPath();
    fs::create_directories(settings.parent_path());
    const fs::path project = dir / "p.mdatk.toml";
    saveProject(project, Project{});
    const fs::path kept = settings.string() + ".unreadable";
    // cut off; Latin-1 text
    const std::vector<std::string> files{"cat_proof = true\nrecent = [\n", "notes = \"caf\xE9\"\n"};
    for (const std::string& broken : files) {
        writeFile(settings, broken);
        const LogCapture log;
        addRecent(project);
        EXPECT_EQ(recentProjects(), std::vector<fs::path>{project});
        EXPECT_EQ(readFile(kept), broken);  // every byte kept
        ASSERT_EQ(log.lines().size(), 1u);
        const std::string& warning = log.lines()[0];
        EXPECT_TRUE(warning.starts_with("WARNING: " + pathToUtf8(settings) + " cannot be read (")) << warning;
        EXPECT_TRUE(warning.ends_with("): it is kept as " + pathToUtf8(kept) + " and a new settings file starts"))
            << warning;
        EXPECT_EQ(warning.find(pathToUtf8(settings) + ": "), std::string::npos) << warning;  // named once
    }

    // cannot be read (here: a folder): nothing is written, nothing is moved
    fs::remove(settings);
    fs::remove(kept);
    fs::create_directories(settings);
    const LogCapture log;
    addRecent(project);
    EXPECT_TRUE(fs::is_directory(settings));
    EXPECT_FALSE(fs::exists(kept));
    EXPECT_FALSE(fs::exists(settings.string() + ".tmp"));
    ASSERT_EQ(log.lines().size(), 1u);
    EXPECT_TRUE(log.lines()[0].starts_with("WARNING: could not update the recent projects list: "))
        << log.lines()[0];
}

// D16: Linux keeps the recent list in MDATK's shared file and MDBridge's own settings in its own
// folder; Windows keeps both in %APPDATA%\MDBridge.
TEST(ProjectIo, ConfigFolders) {
    TempDir dir;
    const ScopedEnvironment config("XDG_CONFIG_HOME", (dir / "xdg").string());
    const ScopedEnvironment appData("APPDATA", (dir / "roaming").string());
#ifdef _WIN32
    EXPECT_EQ(mdbridgeConfigDir(), dir / "roaming" / "MDBridge");
    EXPECT_EQ(userSettingsPath(), dir / "roaming" / "MDBridge" / "settings.toml");
#else
    EXPECT_EQ(mdbridgeConfigDir(), dir / "xdg" / "mdbridge");
    EXPECT_EQ(userSettingsPath(), dir / "xdg" / "mdatk" / "settings.toml");
    const ScopedEnvironment noXdg("XDG_CONFIG_HOME", "");
    const ScopedEnvironment home("HOME", (dir / "home").string());
    EXPECT_EQ(mdbridgeConfigDir(), dir / "home" / ".config" / "mdbridge");
    EXPECT_EQ(userSettingsPath(), dir / "home" / ".config" / "mdatk" / "settings.toml");
#endif
}

// ---- project files of a later MDATK (D9; tests/parity/make_mdatk_latest_projects.py) ---------

namespace {

const fs::path kLatest = kData / "project-file" / "mdatk-0.1.1";

std::string unixLines(const std::string& text) { return replaceAll(text, "\r\n", "\n"); }

}  // namespace

// PRJ-09: MDBridge reads what MDATK 0.1.1 reads from its own files: the settings, the step log
// (statuses waiting and stopped, the step "workflow supersede") and the separate log file; every key
// 0.1.x added to the plan is modelled, so there is no warning.
TEST(ProjectCompat, ReadsWhatMdatkLatestReads) {
    const Value expected = readJson(kLatest / "read.json");
    for (const auto& [name, c] : *expected.find("files")->get<Table>()) {
        const fs::path path = kLatest / name;
        const Value settings = readProject(path).settings;
        EXPECT_EQ(settings, tagged(*c.find("settings"))) << name;
        std::string warning;
        EXPECT_EQ(Value(readLog(path, &warning)), tagged(*c.find("log"))) << name;
        EXPECT_EQ(warning, "") << name;
        EXPECT_EQ(logFileName(path), c.find("log_file")->stringOr()) << name;
        std::vector<std::string> warnings;
        const Project p = projectFromPlain(settings, {}, &warnings);
        EXPECT_EQ(warnings, std::vector<std::string>{}) << name;
        EXPECT_EQ(p.workflow.supersede_runs, false);
        EXPECT_EQ(p.workflow.parallel, 3);
        EXPECT_EQ(p.workflow.cpu_cores, 8);
        EXPECT_EQ(p.workflow.gpus, 1);
        EXPECT_EQ(p.workflow.last_node, 7);
        ASSERT_EQ(p.workflow.nodes.size(), 4u);
        EXPECT_EQ(p.workflow.nodes[1].adopted, "20261101-0900");
        EXPECT_EQ(p.workflow.nodes[1].folder_tag, "model12");
        EXPECT_EQ(p.workflow.nodes[0].adopted, "");  // 0.1.1 writes adopted always, folder_tag when set
        EXPECT_EQ(p.workflow.nodes[0].folder_tag, std::nullopt);
    }
}

// PRJ-09: saving a project MDATK 0.1.1 wrote keeps every key and value in the order of the file (the
// plan's budgets, adopted runs and folder tags, the options of step types 0.0.4.1 does not know) and
// the step log byte for byte, so MDATK 0.1.1 reads the saved file as it read its own.
TEST(ProjectCompat, MdatkLatestProjectsSurviveASave) {
    TempDir dir;
    for (const char* name : {"plan.mdatk.toml", "sidecar.mdatk.toml", "sidecar.log.toml"}) {
        fs::copy_file(kLatest / name, dir / name);
    }
    for (const char* name : {"plan.mdatk.toml", "sidecar.mdatk.toml"}) {
        const fs::path path = dir / name;
        const std::string before = unixLines(readFile(path));
        const Value settings = readProject(path).settings;
        const Array log = readLog(path);
        std::vector<std::string> warnings;
        Project project = projectFromPlain(settings, {}, &warnings);  // as `init -p` and a GUI save
        project.mdatkVersion = "0.1.1";  // what loadProject sets from the file
        saveProject(path, project);
        EXPECT_EQ(warnings, std::vector<std::string>{}) << name;
        EXPECT_EQ(*readProject(path).header.find("mdatk_version"), Value("0.1.1")) << name;
        EXPECT_EQ(readProject(path).settings, settings) << name;
        EXPECT_EQ(readLog(path), log) << name;
        EXPECT_EQ(logFileName(path), logFileName(kLatest / name)) << name;
        const std::string after = unixLines(readFile(path));
        const std::string marker = "\n" + std::string(kLogMarker);
        if (std::string(name) == "plan.mdatk.toml") {  // the notes and a script hold the marker too
            const std::string logPart = before.substr(before.rfind(marker));
            EXPECT_EQ(after.substr(after.size() - logPart.size()), logPart);
        }
    }
    EXPECT_EQ(readFile(dir / "sidecar.log.toml"), readFile(kLatest / "sidecar.log.toml"));
}

// An unknown key of any table, at any depth, is kept with its value and written back; one warning
// names them all (D9).
TEST(ProjectCompat, UnknownKeysAreKeptAtEveryDepth) {
    TempDir dir;
    const fs::path path = dir / "p.mdatk.toml";
    writeFile(path, "format = \"mdatk-project\"\nfuture = { a = 1, b = [1.5, \"x\"] }\n"
                    "[inputs]\ntopology = \"t.parm7\"\nnew_option = 2026-11-25\n"
                    "[[groups]]\nname = \"A\"\nresidues = \"1-3\"\nshade = \"light\"\n"
                    "[analyses.contacts]\ncutoff = 5.0\nmode = \"atoms\"\n"
                    "[[workflow.nodes]]\nid = \"n1\"\ntool = \"note\"\nnew_flag = true\n"
                    "[[workflow.nodes.lanes]]\nname = \"x\"\n");
    std::vector<std::string> warnings;
    const Project p = projectFromPlain(readProject(path).settings, {}, &warnings);
    EXPECT_EQ(warnings, (std::vector<std::string>{
                            "project: 6 unknown keys kept and saved back unchanged (from a newer MDATK, or "
                            "misspelt?): project.future, project.inputs.new_option, project.groups[0].shade, "
                            "project.analyses.contacts.mode, project.workflow.nodes[0].new_flag, "
                            "project.workflow.nodes[0].lanes"}));
    EXPECT_EQ(unknownKeyPaths(p).size(), 6u);
    saveProject(path, p);
    const Value saved = readProject(path).settings;
    EXPECT_EQ(*saved.path({"future"}), (Table{{"a", 1}, {"b", Array{1.5, "x"}}}));
    EXPECT_EQ(saved.path({"inputs", "new_option"})->type(), Value::Type::DateTime);
    EXPECT_EQ(*saved.path({"analyses", "contacts", "mode"}), Value("atoms"));
    EXPECT_EQ(*(*saved.find("groups")->get<Array>())[0].find("shade"), Value("light"));
    const Value& node = (*saved.path({"workflow", "nodes"})->get<Array>())[0];
    EXPECT_EQ(*node.find("new_flag"), Value(true));
    EXPECT_EQ(*node.find("lanes"), (Array{Table{{"name", "x"}}}));
    std::vector<std::string> again;
    EXPECT_EQ(projectFromPlain(saved, {}, &again), p);  // nothing lost or added on the second round
}

// PRJ-10 (MDATK 0.1.0): the step-log marker pasted into the notes or a script does not start the
// log. MDATK 0.0.4.1 took it for the log: the next entry was appended without a marker of its own and
// the next save cut the settings at the notes, leaving a file that was no TOML.
TEST(ProjectCompat, LogMarkerInsideTheSettingsIsNoLog) {
    TempDir dir;
    const fs::path path = dir / "p.mdatk.toml";
    Project p;
    p.notes = "pasted:\n" + std::string(kLogMarker) + "\n[[log]]\nstep = \"x\"\n";
    saveProject(path, p);
    EXPECT_TRUE(readLog(path).empty());
    appendLog(path, LogEntry("render", "first"));
    ASSERT_EQ(readLog(path).size(), 1u);
    p.analyses.contacts.cutoff = 5.0;
    saveProject(path, p);
    EXPECT_EQ(fromFile(path).notes, p.notes);
    EXPECT_EQ(fromFile(path).analyses.contacts.cutoff, 5.0);
    ASSERT_EQ(readLog(path).size(), 1u);
    appendLog(path, LogEntry("render", "second"));
    EXPECT_EQ(readLog(path).size(), 2u);
    const fs::path sidecar = splitLog(path);
    EXPECT_EQ(sidecar, dir / "p.log.toml");
    EXPECT_EQ(fromFile(path).notes, p.notes);
    EXPECT_EQ(readLog(path).size(), 2u);
}

// PRJ-10 (MDATK 0.1.0): a separate log file that is gone (deleted, or the project copied without it)
// starts again with the next entry instead of failing the step's log.
TEST(ProjectCompat, MissingSeparateLogStartsAgain) {
    TempDir dir;
    const fs::path path = dir / "sidecar.mdatk.toml";
    fs::copy_file(kLatest / "sidecar.mdatk.toml", path);
    EXPECT_EQ(appendLog(path, LogEntry("render", "after the copy")), dir / "sidecar.log.toml");
    const Array log = readLog(path);
    ASSERT_EQ(log.size(), 1u);
    EXPECT_EQ(*log[0].find("title"), Value("after the copy"));
    EXPECT_EQ(unixLines(readFile(dir / "sidecar.log.toml")).rfind("\n" + std::string(kLogMarker) + "\n\n# One entry", 0),
              0u);
}

// PRJ-09: the mdatk_version of a later MDATK's file is kept on save; files without one, or of
// 0.0.4.1 and older, are written as 0.0.4.1.
TEST(ProjectCompat, MdatkVersionOfALaterMdatkIsKept) {
    TempDir dir;
    const fs::path path = dir / "plan.mdatk.toml";
    fs::copy_file(kLatest / "plan.mdatk.toml", path);
    Project p = loadProject(path);
    EXPECT_EQ(p.mdatkVersion, "0.1.1");
    saveProject(dir / "copy.mdatk.toml", p);
    EXPECT_EQ(*readProject(dir / "copy.mdatk.toml").header.find("mdatk_version"), Value("0.1.1"));
    const auto written = [](const std::string& version) {
        Project q;
        q.mdatkVersion = version;
        return writtenMdatkVersion(q);
    };
    for (const char* older : {"", "0.0.4.1", "0.0.4", "0.0.3d", "0.0.1", "x", ".1", "abc.0.2"}) {
        EXPECT_EQ(written(older), "0.0.4.1") << older;
    }
    for (const char* newer : {"0.1.1", "0.0.10", "0.0.4.2", "0.0.4.1a", "1", "0.1.1.dev3"}) {
        EXPECT_EQ(written(newer), newer) << newer;
    }
    p.mdatkVersion.clear();
    saveProject(path, p);
    EXPECT_EQ(*readProject(path).header.find("mdatk_version"), Value("0.0.4.1"));
}

// PRJ-13: a save that cannot replace the file leaves the old file as it was and no <name>.tmp.
TEST(ProjectSave, AFailedReplaceLeavesNoTemporaryFile) {
    TempDir dir;
    const fs::path path = dir / "p.mdatk.toml";
    fs::create_directories(path / "inside");  // a folder of that name: the rename must fail
    EXPECT_THROW(saveProject(path, Project{}), fs::filesystem_error);
    EXPECT_TRUE(fs::is_directory(path / "inside"));
    EXPECT_FALSE(fs::exists(dir / "p.mdatk.toml.tmp"));
}

#ifdef _WIN32
// PRJ-13: another program reading the project while it is saved, as virus scanners and the search
// indexer open files (sharing read, write and delete). Windows refuses to replace an open file;
// the save waits for it (GuiTests ap13_unsaved_changes failed so in CI run 45790979168).
TEST(ProjectSave, WaitsForAProgramReadingTheFile) {
    TempDir dir;
    const fs::path path = dir / "p.mdatk.toml";
    saveProject(path, Project{});
    const HANDLE held = CreateFileW(path.c_str(), GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                                    nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    ASSERT_NE(held, INVALID_HANDLE_VALUE);
    std::thread release([held] {
        std::this_thread::sleep_for(std::chrono::milliseconds(300));
        CloseHandle(held);
    });
    Project q;
    q.notes = "saved while another program read the file";
    EXPECT_NO_THROW(saveProject(path, q));
    release.join();
    EXPECT_EQ(fromFile(path).notes, q.notes);
    EXPECT_FALSE(fs::exists(dir / "p.mdatk.toml.tmp"));
}

// PRJ-13: a file held open (without delete sharing) for longer than the 2 s the save waits: the
// save fails with the system's message and the file stays as it was.
TEST(ProjectSave, GivesUpOnAFileHeldOpen) {
    TempDir dir;
    const fs::path path = dir / "p.mdatk.toml";
    Project first;
    first.notes = "first";
    saveProject(path, first);
    const HANDLE held = CreateFileW(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                                    FILE_ATTRIBUTE_NORMAL, nullptr);
    ASSERT_NE(held, INVALID_HANDLE_VALUE);
    const auto start = std::chrono::steady_clock::now();
    EXPECT_THROW(saveProject(path, Project{}), fs::filesystem_error);
    EXPECT_GE(std::chrono::steady_clock::now() - start, std::chrono::milliseconds(1900));
    CloseHandle(held);
    EXPECT_EQ(fromFile(path).notes, "first");
    EXPECT_FALSE(fs::exists(dir / "p.mdatk.toml.tmp"));
}
#endif
