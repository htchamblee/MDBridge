// The demo system (src/app/commands/demo.cpp), the process runner (src/core/process.cpp), the
// cpptraj cleaning (src/analysis/cleaning.cpp) and the citations of a project
// (src/app/commands/citations.cpp): MDBridge's own behaviour and improvements (docs/DIFFERENCES.md
// CMD-01 ... CMD-06). What MDATK prints and writes for demo and clean is compared by
// tests/cli/test_commands.cpp (cli_cases.json, tests/parity/cli_parity.py).
#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "citations.hpp"
#include "commands.hpp"
#include "demo.hpp"
#include "mdbridge/analysis/cache.hpp"
#include "mdbridge/analysis/cleaning.hpp"
#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/process.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/version.hpp"
#include "test_support.hpp"

namespace mdbridge {
namespace {

namespace fs = std::filesystem;
using test::kData;

struct CommandRun {
    int exit = 0;
    std::string out;
    std::string err;
};

CommandRun runIn(const fs::path& folder, const std::vector<std::string>& args) {
    std::vector<const char*> argv{"mdbridge"};
    for (const std::string& a : args) {
        argv.push_back(a.c_str());
    }
    std::ostringstream out, err;
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

std::string utf8(const fs::path& p) { return core::pathToUtf8(p); }

// the separator os.path.join puts between the demo's folder and its files ('\\' on Windows, ntpath)
#ifdef _WIN32
constexpr const char* kJoin = "\\";
#else
constexpr const char* kJoin = "/";
#endif

// ---- demo ---------------------------------------------------------------------------------------

TEST(Demo, DrySystemIsMdatksByteForByte) {
    // tests/data's demo_dry files are MDATK's build_demo(n_frames=24, solvent=False)
    // (tests/parity/make_selection_fixtures.py); the solvated systems are compared by cli_cases.json
    test::TempDir work;
    app::DemoOptions options;
    options.frames = 24;
    options.solvent = false;
    const app::DemoFiles files = app::buildDemo(utf8(work / "d"), options);
    EXPECT_EQ(files.topology, utf8(work / "d") + kJoin + "demo.parm7");
    EXPECT_EQ(files.entries().size(), 4u);
    EXPECT_EQ(files.entries()[3].second, files.project);
    EXPECT_TRUE(test::readFile(work / "d" / "demo.parm7") == test::readFile(kData / "topology" / "demo_dry.parm7"));
    EXPECT_TRUE(test::readFile(work / "d" / "demo.nc") == test::readFile(kData / "io" / "trajectories" / "demo_dry.nc"));
    const core::Project p = core::loadProject(work / "d" / "demo.mdatk.toml");
    EXPECT_FALSE(p.cleaning.enabled);  // no solvent: nothing to clean
    EXPECT_EQ(p.layout.preset, "large_plus_4");
    EXPECT_EQ(p.layout.panels.back().title, "Mg occupancy");
    EXPECT_EQ(p.output.pdb.times, "120ps");
    EXPECT_EQ(p.groups[1].residues, "13-24");
}

TEST(Demo, ProjectSettings) {
    const core::Project p = app::demoProject("t.parm7", "t.nc", "out", 4, true, 12.5);
    EXPECT_EQ(p.groups[0].residues, "1-4");
    EXPECT_EQ(p.groups[1].residues, "5-8");
    EXPECT_TRUE(p.cleaning.enabled);
    EXPECT_EQ(p.output.pdb.times, "12.5ps");
    EXPECT_EQ(p.analyses.energy.target_value, -60.0);
    EXPECT_EQ(p.visual.rotation, (std::vector<double>{-90.0, 0.0, 0.0}));
}

// CMD-01: impossible arguments are refused before anything is written
TEST(Demo, RefusesImpossibleArguments) {
    test::TempDir work;
    const auto refused = [&](app::DemoOptions options, const std::string& message) {
        try {
            app::buildDemo(utf8(work / "x"), options);
            ADD_FAILURE() << message;
        } catch (const core::ValueError& e) {
            EXPECT_EQ(std::string(e.what()), message);
        }
        EXPECT_FALSE(fs::exists(work / "x"));
    };
    app::DemoOptions o;
    o.frames = 0;
    refused(o, "the demo needs at least one frame (--n-frames 0)");
    o = {};
    o.dtPs = -2.5;
    refused(o, "the time between frames must be a positive number of ps (--dt -2.5)");
    o.dtPs = std::nan("");
    refused(o, "the time between frames must be a positive number of ps (--dt nan)");
    o = {};
    o.sequence = "GCT";
    refused(o, "demo sequence 'GCT': at least four of the letters A, C, G and U");
}

TEST(Demo, ProgressAndCancel) {
    test::TempDir work;
    std::vector<double> fractions;
    core::Reporter reporter({}, [&](double f, const std::string&) { fractions.push_back(f); });
    app::DemoOptions options;
    options.frames = 4;
    options.solvent = false;
    app::buildDemo(utf8(work / "a"), options, reporter);
    EXPECT_EQ(fractions, (std::vector<double>{0.25, 0.5, 0.75, 1.0}));
    core::Reporter cancelled;
    cancelled.cancel();
    EXPECT_THROW(app::buildDemo(utf8(work / "b"), options, cancelled), core::Cancelled);
    EXPECT_FALSE(fs::exists(work / "b" / "demo.parm7"));
}

TEST(Demo, Command) {
    test::TempDir work;
    const CommandRun run = runIn(work.path(), {"demo", "sub/d", "--n-frames", "3", "--no-solvent"});
    EXPECT_EQ(run.exit, 0) << run.err;
    // os.path.abspath of the working folder as the process sees it: on Windows the folder as given
    // (8.3 short names such as RUNNER~1 are not expanded, unlike fs::canonical), with '\\'
#ifdef _WIN32
    const std::string outDir = utf8(fs::absolute(work.path() / "sub" / "d").lexically_normal());
#else
    const std::string outDir = utf8(fs::canonical(work.path()) / "sub" / "d");
#endif
    const std::string d = std::string("sub/d") + kJoin;
    EXPECT_EQ(run.out, "Synthetic system written to " + outDir + ":\n" +
                           "  topology   " + d + "demo.parm7\n" +
                           "  trajectory " + d + "demo.nc\n" +
                           "  project    " + d + "demo.mdatk.toml\n" +
                           "  preset     " + d + "demo.mdatk.toml\n" +
                           "Next: mdbridge render -p " + d + "demo.mdatk.toml   (or: mdbridge gui " + d +
                           "demo.mdatk.toml)\n");
    // --run without a renderer (this test's command line has none): the demo is written, render fails
    // (the counterpart of MDATK's test_cli_demo_render_and_tools is in test_demo_commands.cpp)
    const CommandRun rendered = runIn(work.path(), {"demo", "r", "--n-frames", "2", "--no-solvent", "--run", "-q"});
    EXPECT_EQ(rendered.exit, 1);
    EXPECT_TRUE(fs::exists(work / "r" / "demo.nc"));
    EXPECT_TRUE(rendered.err.starts_with("error: ")) << rendered.err;
}

// ---- the process runner -------------------------------------------------------------------------

TEST(Process, LineSplitter) {
    core::process::LineSplitter s;
    EXPECT_TRUE(s.add("a b  \r\nc").size() == 1);
    std::vector<std::string> lines = s.add("\xC3");  // half of "ä"
    EXPECT_TRUE(lines.empty());
    lines = s.add("\xA4 \t\nbad \xFF\xFE x\n");
    ASSERT_EQ(lines.size(), 2u);
    EXPECT_EQ(lines[0], "c\xC3\xA4");
    EXPECT_EQ(lines[1], "bad \xEF\xBF\xBD\xEF\xBF\xBD x");
    s.add("rest  ");
    EXPECT_EQ(s.finish(), std::optional<std::string>("rest"));
    EXPECT_EQ(s.finish(), std::nullopt);
    EXPECT_EQ(core::process::windowsCommandLine({"a b", "c\"d", "e\\", "f\\ g\\", ""}),
              "\"a b\" c\\\"d e\\ \"f\\ g\\\\\" \"\"");
    EXPECT_EQ(core::process::commandLineText({"conda", "run", "-n", "x y"}), "conda run -n x y");
}

#ifndef _WIN32

void writeScript(const fs::path& path, const std::string& text) {
    fs::create_directories(path.parent_path());
    test::writeFile(path, text);
    fs::permissions(path, fs::perms::owner_all, fs::perm_options::add);
}

TEST(Process, RunsAndStreams) {
    test::TempDir work;
    writeScript(work / "p.sh", "#!/bin/sh\necho out; echo err >&2; pwd; printf '%s\\n' \"$MDB_X\"; exit 7\n");
    std::string output;
    core::process::RunOptions options;
    options.workingDirectory = work.path() / "";
    options.environment = std::vector<std::pair<std::string, std::string>>{{"MDB_X", "given"}, {"PATH", "/usr/bin:/bin"}};
    options.onOutput = [&](std::string_view chunk) { output += chunk; };
    const core::process::Result r = core::process::run({utf8(work / "p.sh")}, options);
    EXPECT_EQ(r.exitCode, 7);
    EXPECT_FALSE(r.cancelled);
    EXPECT_EQ(output, "out\nerr\n" + utf8(fs::canonical(work.path())) + "\ngiven\n");
    EXPECT_THROW(core::process::run({"mdbridge-no-such-program"}), core::process::StartError);
    try {
        core::process::run({utf8(work / "missing.sh")});
        ADD_FAILURE();
    } catch (const core::process::StartError& e) {
        EXPECT_EQ(std::string(e.what()), "[Errno 2] No such file or directory: '" + utf8(work / "missing.sh") + "'");
    }
    writeScript(work / "kill.sh", "#!/bin/sh\nkill -9 $$\n");
    EXPECT_EQ(core::process::run({utf8(work / "kill.sh")}).exitCode, -9);
}

// CMD-03: a cancel ends a program that prints nothing, and its children (MDATK: only between reads)
TEST(Process, CancelEndsTheProcessGroup) {
    test::TempDir work;
    writeScript(work / "slow.sh", "#!/bin/sh\necho started\nsleep 30 &\necho $! > child.pid\nwait\n");
    std::atomic<bool> cancel{false};
    core::process::RunOptions options;
    options.workingDirectory = work.path();
    options.terminateGrace = std::chrono::milliseconds(2000);
    std::thread canceller([&cancel, &work] {
        for (int k = 0; k < 200 && !fs::exists(work / "child.pid"); ++k) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        cancel = true;
    });
    options.cancelled = [&] { return cancel.load(); };
    const auto start = std::chrono::steady_clock::now();
    const core::process::Result r = core::process::run({utf8(work / "slow.sh")}, options);
    canceller.join();
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    EXPECT_TRUE(r.cancelled);
    EXPECT_LT(seconds, 5.0);
    const std::string pid = test::readFile(work / "child.pid");
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    EXPECT_FALSE(fs::exists("/proc/" + pid.substr(0, pid.find('\n')) + "/cmdline") &&
                 test::readFile("/proc/" + pid.substr(0, pid.find('\n')) + "/cmdline").starts_with("sleep"))
        << "the background child survived";
}

// CMD-07: started from the AppImage, programs do not inherit its library folders
TEST(Process, ChildEnvironmentLeavesTheAppImageOut) {
    const auto value = [](const std::vector<std::pair<std::string, std::string>>& env, const std::string& name) {
        for (const auto& [n, v] : env) {
            if (n == name) {
                return std::optional<std::string>(v);
            }
        }
        return std::optional<std::string>();
    };
    std::vector<std::pair<std::string, std::optional<std::string>>> saved;
    for (const char* name : {"APPDIR", "LD_LIBRARY_PATH", "FONTCONFIG_FILE", "MDB_KEPT"}) {
        const char* v = std::getenv(name);
        saved.emplace_back(name, v ? std::optional<std::string>(v) : std::nullopt);
    }
    setenv("APPDIR", "/tmp/.mount_MDBxyz", 1);
    setenv("LD_LIBRARY_PATH", "/tmp/.mount_MDBxyz/usr/lib/fallback:/tmp/.mount_MDBxyz/usr/lib:/opt/amber/lib", 1);
    setenv("FONTCONFIG_FILE", "/tmp/.mount_MDBxyz/usr/share/mdbridge/fontconfig/fonts.conf", 1);
    setenv("MDB_KEPT", "x", 1);
    auto env = core::process::childEnvironment();
    EXPECT_EQ(value(env, "LD_LIBRARY_PATH"), std::optional<std::string>("/opt/amber/lib"));
    EXPECT_EQ(value(env, "FONTCONFIG_FILE"), std::nullopt);
    EXPECT_EQ(value(env, "MDB_KEPT"), std::optional<std::string>("x"));
    setenv("LD_LIBRARY_PATH", "/tmp/.mount_MDBxyz/usr/lib", 1);
    EXPECT_EQ(value(core::process::childEnvironment(), "LD_LIBRARY_PATH"), std::nullopt);
    unsetenv("APPDIR");  // not from the AppImage: everything as it is
    env = core::process::childEnvironment();
    EXPECT_EQ(value(env, "LD_LIBRARY_PATH"), std::optional<std::string>("/tmp/.mount_MDBxyz/usr/lib"));
    for (const auto& [name, v] : saved) {
        if (v) {
            setenv(name.c_str(), v->c_str(), 1);
        } else {
            unsetenv(name.c_str());
        }
    }
}

TEST(Process, Which) {
    test::TempDir work;
    writeScript(work / "bin" / "tool", "#!/bin/sh\n");
    test::writeFile(work / "bin" / "plain", "not executable");
    const char* old = std::getenv("PATH");
    const std::string saved = old ? old : "";
    setenv("PATH", (utf8(work / "nothing") + ":" + utf8(work / "bin")).c_str(), 1);
    EXPECT_EQ(core::process::which("tool"), utf8(work / "bin" / "tool"));
    EXPECT_EQ(core::process::which("plain"), std::nullopt);
    EXPECT_EQ(core::process::which(utf8(work / "bin" / "tool")), utf8(work / "bin" / "tool"));
    EXPECT_EQ(core::process::which("bin/tool"), std::nullopt);  // relative to the current folder
    setenv("PATH", saved.c_str(), 1);
}

// ---- cleaning -----------------------------------------------------------------------------------

const char* kFakeCpptraj = R"SH(#!/bin/sh
script="$2"
parm=$(sed -n 's/^parm "\(.*\)"$/\1/p' "$script")
traj=$(sed -n 's/^trajin "\([^"]*\)".*/\1/p' "$script" | head -n 1)
parmout=$(sed -n 's/.* parmout "\([^"]*\)"$/\1/p' "$script")
nc=$(sed -n 's/^trajout "\([^"]*\)" netcdf$/\1/p' "$script")
[ -n "$parmout" ] && cp "$parm" "$parmout"
[ -n "$MDB_SLOW" ] && { echo waiting; sleep 30; }
cp "$traj" "$nc"
printf '  0%% 50%%' ; printf ' 100%% Complete.\n'
)SH";

struct CleaningFixture {
    test::TempDir work;
    core::Project project;
    CleaningFixture() {
        fs::copy_file(kData / "topology" / "demo_dry.parm7", work / "demo_dry.parm7");
        fs::copy_file(kData / "io" / "trajectories" / "demo_dry.nc", work / "demo_dry.nc");
        writeScript(work / "bin" / "cpptraj", kFakeCpptraj);
        project.inputs.topology = utf8(work / "demo_dry.parm7");
        project.inputs.trajectories = {utf8(work / "demo_dry.nc")};
        project.output.out_dir = utf8(work / "out");
        project.cleaning.enabled = true;
        project.cleaning.cpptraj_path = utf8(work / "bin" / "cpptraj");
    }
};

TEST(Cleaning, RunsReusesAndForces) {
    CleaningFixture f;
    std::vector<std::string> log;
    std::vector<std::pair<double, std::string>> progress;
    const core::Reporter reporter([&](const std::string& m) { log.push_back(m); },
                                  [&](double x, const std::string& m) { progress.emplace_back(x, m); });
    const analysis::CleanResult first = analysis::clean(f.project.cleaning, f.project.inputs.topology,
                                                        f.project.inputs.trajectories, core::cacheDir(f.project), reporter);
    EXPECT_FALSE(first.reused);
    EXPECT_EQ(first.outDir, utf8(f.work / "out" / "MDATK_cache" / "cleaned" / first.key));
    EXPECT_TRUE(fs::exists(first.topology));
    EXPECT_TRUE(test::readFile(first.trajectory) == test::readFile(f.work / "demo_dry.nc"));
    EXPECT_EQ(test::readFile(first.log), "  0% 50% 100% Complete.\n");
    ASSERT_FALSE(log.empty());
    EXPECT_EQ(log.front(), "cpptraj script written to " + first.script + ":");
    EXPECT_EQ(log.back(), "  0% 50% 100% Complete.");
    EXPECT_EQ(progress.back(), (std::pair<double, std::string>(1.0, "cpptraj")));
    const core::Value marker = test::readJson(fs::path(first.outDir) / "MDATK_clean.json");
    EXPECT_EQ(marker.find("key")->stringOr(), first.key);
    EXPECT_EQ(marker.find("cpptraj")->stringOr(), f.project.cleaning.cpptraj_path);
    EXPECT_TRUE(marker.find("params")->find("stop")->isNull());

    log.clear();
    const analysis::CleanResult again = analysis::clean(f.project.cleaning, f.project.inputs.topology,
                                                        f.project.inputs.trajectories, core::cacheDir(f.project), reporter);
    EXPECT_TRUE(again.reused);
    EXPECT_EQ(log, (std::vector<std::string>{"Reusing cleaned trajectory from " + first.outDir}));
    const analysis::CleanResult forced = analysis::clean(f.project.cleaning, f.project.inputs.topology,
                                                         f.project.inputs.trajectories, core::cacheDir(f.project),
                                                         reporter, true);
    EXPECT_FALSE(forced.reused);

    // the pipeline uses the cleaned files; previews find them without running anything
    analysis::Pipeline pipeline(f.project);
    const analysis::ResolvedInputs& in = pipeline.resolveInputs();
    EXPECT_EQ(in.topology, first.topology);
    ASSERT_TRUE(in.cleaned.has_value());
    EXPECT_TRUE(in.cleaned->reused);
    EXPECT_EQ(analysis::Pipeline(f.project).cleanedIfAvailable().trajectories,
              std::vector<std::string>{first.trajectory});
}

// MDATK 0.1.0: cleaned files in a folder of cleaning.output_dir count only with this key
TEST(Cleaning, OutputDirNeedsTheKey) {
    CleaningFixture f;
    f.project.cleaning.output_dir = utf8(f.work / "cleaned");
    const core::Reporter quiet;
    analysis::clean(f.project.cleaning, f.project.inputs.topology, f.project.inputs.trajectories,
                    core::cacheDir(f.project), quiet);
    EXPECT_TRUE(analysis::cleanedIfAvailable(f.project.cleaning, f.project.inputs.topology,
                                             f.project.inputs.trajectories, core::cacheDir(f.project)));
    core::CleaningParams other = f.project.cleaning;
    other.offset = 2;
    EXPECT_FALSE(analysis::cleanedIfAvailable(other, f.project.inputs.topology, f.project.inputs.trajectories,
                                              core::cacheDir(f.project)));
    EXPECT_FALSE(analysis::clean(other, f.project.inputs.topology, f.project.inputs.trajectories,
                                 core::cacheDir(f.project), quiet)
                     .reused);
}

// CMD-03: a cancel stops cpptraj at once although it prints nothing; no marker is left
TEST(Cleaning, Cancel) {
    CleaningFixture f;
    setenv("MDB_SLOW", "1", 1);
    core::Reporter reporter;
    std::thread canceller([&reporter] {
        std::this_thread::sleep_for(std::chrono::milliseconds(300));
        reporter.cancel();
    });
    const auto start = std::chrono::steady_clock::now();
    try {
        analysis::clean(f.project.cleaning, f.project.inputs.topology, f.project.inputs.trajectories,
                        core::cacheDir(f.project), reporter);
        ADD_FAILURE();
    } catch (const core::Cancelled& e) {
        EXPECT_EQ(std::string(e.what()), "cpptraj cancelled");
    }
    canceller.join();
    unsetenv("MDB_SLOW");
    EXPECT_LT(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(), 10.0);
    const std::string key = analysis::cleaningKey(f.project.cleaning, f.project.inputs.topology,
                                                  f.project.inputs.trajectories);
    EXPECT_FALSE(fs::exists(f.work / "out" / "MDATK_cache" / "cleaned" / key / "MDATK_clean.json"));
}

TEST(Cleaning, Script) {
    core::CleaningParams p;
    p.fit_reference = "file";
    EXPECT_THROW(analysis::buildCpptrajScript(p, "a.parm7", {"a.nc"}, "/o"), analysis::CpptrajError);
    p.reference_file = "ref.pdb";
    p.solvent_mask = " :WAT, HOH ";
    p.ion_mask = "";
    p.extra_strip_mask = "@H\t=";
    p.autoimage_anchor = ":1 - 12";
    p.fit_mask = "  ";
    p.offset = 0;
    p.stop = 50;
    p.extra_commands = "\n  center :1  \r\nimage\n\n";
    const std::string script = analysis::buildCpptrajScript(p, "/w/a.PRMTOP", {"/w/a.nc", "/w/b.nc"}, "/o");
    const std::string body = script.substr(script.find('\n') + 1);
    EXPECT_TRUE(script.starts_with(std::string("# Generated by MDBridge ") + kVersion + " on ")) << script;  // CMD-08
    EXPECT_EQ(body,
              "parm \"/w/a.PRMTOP\"\n"
              "trajin \"/w/a.nc\" 1 50 1\n"
              "trajin \"/w/b.nc\" 1 50 1\n"
              "reference \"" + core::posixpath::abspath("ref.pdb", core::currentDirectory()) + "\" [mdatkref]\n"
              "strip (:WAT,HOH)|(@H=) parmout \"/o/MDATK_cleaned.parm7\"\n"
              "autoimage anchor :1-12\n"
              "box nobox\n"
              "rms MDATKfit !@H= ref [mdatkref] mass\n"
              "center :1\n"
              "image\n"
              "trajout \"/o/MDATK_cleaned.nc\" netcdf\n"
              "trajout \"/o/MDATK_cleaned_ref.pdb\" pdb onlyframes 1\n"
              "run\n"
              "quit\n");
    EXPECT_THROW(analysis::buildCpptrajScript({}, "a.parm7", {}, "/o"), analysis::CpptrajError);
    EXPECT_THROW(analysis::buildCpptrajScript({}, "a\".parm7", {"x.nc"}, "/o"), analysis::CpptrajError);
}

#endif  // !_WIN32

// ---- citations (D14) ----------------------------------------------------------------------------

TEST(Citations, SettingsAndRunKeys) {
    test::TempDir work;
    core::Project project;
    project.output.out_dir = utf8(work / "out");
    project.cleaning.enabled = true;
    app::ProjectCitations c = app::projectCitations(project);
    EXPECT_TRUE(c.cacheKeys.empty());
    EXPECT_EQ(c.keys, c.settingsKeys);
    ASSERT_FALSE(c.keys.empty());
    const std::string text = app::citationText(c);
    EXPECT_TRUE(text.starts_with("Resources this project uses (" + std::to_string(c.keys.size()) +
                                 " references; each generated script lists only its own):\n\n"))
        << text;
    EXPECT_NE(text.find("\nCPPTRAJ\n    Roe D. R."), std::string::npos) << text;

    // a cached result recording a run-time citation (the covalent radii of a bond guess)
    const analysis::CacheStore cache(core::cacheDir(project));
    cache.ensure();
    core::AnalysisResult result;
    result.name = "atoms";
    result.key = "0123456789abcdef";
    result.meta.emplace_back("citation_keys", core::Array{"cordero2008", "not-a-key"});
    cache.saveResult(result, std::vector<std::int64_t>{0}, std::vector<double>{0.0});
    const core::Dataset dataset({0}, {0.0});
    cache.writeManifest(core::Table{}, core::Table{}, dataset, core::Array{}, {{"atoms", "0123456789abcdef"}});
    c = app::projectCitations(project);
    EXPECT_EQ(c.cacheKeys, (std::vector<std::string>{"cordero2008", "not-a-key"}));
    EXPECT_EQ(c.runOnly, std::vector<std::string>{"cordero2008"});
    EXPECT_EQ(c.keys.back(), "cordero2008");
    EXPECT_NE(app::citationText(c).find("\nRecorded by the cached results of the last run:\n    "), std::string::npos);
    const core::Value json = app::citationJson(c);
    EXPECT_EQ(json.find("run_keys")->get<core::Array>()->size(), 2u);
    EXPECT_EQ(json.find("keys")->get<core::Array>()->size(), c.keys.size());

    // an entry the manifest lists but the cache lacks: nothing is recorded
    cache.writeManifest(core::Table{}, core::Table{}, dataset, core::Array{},
                        {{"atoms", "0123456789abcdef"}, {"rmsd", "fedcba9876543210"}});
    EXPECT_TRUE(app::projectCitations(project).cacheKeys.empty());

    core::Project none;
    none.output.out_dir = utf8(work / "none");
    none.analyses.rmsd.enabled = false;
    none.output.video.enabled = false;
    const app::ProjectCitations empty = app::projectCitations(none);
    EXPECT_TRUE(empty.keys.empty());
    EXPECT_EQ(app::citationText(empty), "This project uses no external resource with a reference yet.\n");
}

TEST(Citations, Command) {
    test::TempDir work;
    test::writeFile(work / "p.mdatk.toml", "format = \"mdatk-project\"\n[cleaning]\nenabled = true\n");
    const CommandRun text = runIn(work.path(), {"citations", "-p", "p.mdatk.toml"});
    EXPECT_EQ(text.exit, 0) << text.err;
    EXPECT_TRUE(text.out.starts_with("Resources this project uses (")) << text.out;
    const CommandRun json = runIn(work.path(), {"citations", "-p", "p.mdatk.toml", "--json", "--bibtex", "refs.bib"});
    EXPECT_EQ(json.exit, 0) << json.err;
    const core::Value value = core::parseJson(json.out, "stdout");
    std::vector<std::string> keys;
    for (const core::Value& k : *value.find("keys")->get<core::Array>()) {
        keys.push_back(k.stringOr());
    }
    EXPECT_EQ(test::readFile(work / "refs.bib"), core::citations::bibtex(keys));
    const CommandRun bib = runIn(work.path(), {"citations", "-p", "p.mdatk.toml", "--bibtex", "b.bib"});
    EXPECT_EQ(bib.out, "Wrote " + std::to_string(keys.size()) + " references to b.bib\n");
    const CommandRun missing = runIn(work.path(), {"citations", "-p", "nope.mdatk.toml"});
    EXPECT_EQ(missing.exit, 1);
}

}  // namespace
}  // namespace mdbridge
