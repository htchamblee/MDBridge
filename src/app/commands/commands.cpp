// MDATK's milestone-2 commands (mdatk/cli.py cmd_init, cmd_info, cmd_export_pdb, cmd_density_info
// and the parts of load_project, make_reporter and Pipeline they use).
#include "commands.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <CLI/CLI.hpp>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/loader.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/io/pdb_export.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;
using core::Value;

namespace {

// The options MDATK's _common() adds to init, info and export-pdb.
struct CommonArguments {
    std::string params;
    std::string top;
    std::vector<std::string> traj;
    std::string frames;
    bool clean = false;
    bool noClean = false;
    std::string cpptraj;
    std::string condaEnv;
    std::string outDir;
    std::string cacheDir;
    std::vector<std::string> set;
    bool quiet = false;
    bool verbose = false;
    CLI::Option* cleanOption = nullptr;
    CLI::Option* noCleanOption = nullptr;
};

}  // namespace

struct CoreCommands::Arguments {
    CLI::App* init = nullptr;
    CLI::App* info = nullptr;
    CLI::App* exportPdb = nullptr;
    CLI::App* density = nullptr;
    CLI::App* densityInfo = nullptr;

    CommonArguments initCommon;
    std::string initOutput;
    bool overwrite = false;

    CommonArguments infoCommon;

    CommonArguments exportCommon;
    // given (even as "") or not: CLI11 leaves an optional unset for an empty value
    std::string frameList, times, selection;
    CLI::Option* frameListOption = nullptr;
    CLI::Option* timesOption = nullptr;
    CLI::Option* selectionOption = nullptr;

    std::string map;
};

namespace {

void addCommon(CLI::App& command, CommonArguments& a) {
    command.add_option("--project,-p,--params", a.params,
                       "project file (.mdatk.toml) or a run's MDATK_params_log.toml");
    command.add_option("--top", a.top, "topology file (.parm7, .pdb, ...)");
    command.add_option("--traj", a.traj, "trajectory file(s) (.nc, .dcd, ...)");
    command.add_option("--frames", a.frames, "frame selection start:stop:step (0-based, stop exclusive)");
    a.cleanOption = command.add_flag("--clean", a.clean, "enable cpptraj cleaning");
    a.noCleanOption = command.add_flag("--no-clean", a.noClean, "disable cleaning");
    command.add_option("--cpptraj", a.cpptraj, "path to the cpptraj executable");
    command.add_option("--conda-env", a.condaEnv,
                       "conda environment providing cpptraj (e.g. AmberTools26); cpptraj then runs via `conda run`");
    command.add_option("--out-dir,-o", a.outDir, "output directory");
    command.add_option("--cache-dir", a.cacheDir, "cache directory (default <out-dir>/MDATK_cache)");
    command.add_option("--set,-s", a.set, "override a parameter, e.g. analyses.contacts.cutoff=5 (repeatable)")
        ->type_name("KEY=VALUE")
        ->allow_extra_args(false);
    command.add_flag("--quiet,-q", a.quiet, "only print the result");
    command.add_flag("--verbose,-v", a.verbose, "also print full cpptraj output");
}

// os.path.abspath
std::string absolutePath(const std::string& path) {
#ifdef _WIN32
    return core::pathToUtf8(fs::absolute(core::pathFromUtf8(path)).lexically_normal());
#else
    return core::posixpath::abspath(path, core::currentDirectory());
#endif
}

// --clean / --no-clean (store_true / store_false into one destination): the last one given wins
std::optional<bool> cleanChoice(const CLI::App& command, const CommonArguments& a) {
    std::optional<bool> choice;
    for (const CLI::Option* option : command.parse_order()) {
        if (option == a.cleanOption) {
            choice = true;
        } else if (option == a.noCleanOption) {
            choice = false;
        }
    }
    return choice;
}

// load_project(args, resolve)
core::Project loadProject(const CLI::App& command, const CommonArguments& a, bool resolve) {
    core::Project project = a.params.empty() ? core::Project{} : core::loadProject(core::pathFromUtf8(a.params));
    if (!a.top.empty()) {
        project.inputs.topology = absolutePath(a.top);
    }
    if (!a.traj.empty()) {
        project.inputs.trajectories.clear();
        for (const std::string& t : a.traj) {
            project.inputs.trajectories.push_back(absolutePath(t));
        }
    }
    if (!a.outDir.empty()) {
        project.output.out_dir = absolutePath(a.outDir);
    }
    if (!a.cacheDir.empty()) {
        project.output.cache_dir = absolutePath(a.cacheDir);
    }
    if (!a.cpptraj.empty()) {
        project.cleaning.cpptraj_path = a.cpptraj;
    }
    if (!a.condaEnv.empty()) {
        project.cleaning.conda_env = a.condaEnv;
    }
    if (const std::optional<bool> clean = cleanChoice(command, a)) {
        project.cleaning.enabled = *clean;
    }
    if (!a.frames.empty()) {
        std::vector<std::string> parts;  // (frames.split(":") + ["", "", ""])[:3]
        std::size_t start = 0;
        while (true) {
            const std::size_t colon = a.frames.find(':', start);
            parts.push_back(a.frames.substr(start, colon == std::string::npos ? std::string::npos : colon - start));
            if (colon == std::string::npos) {
                break;
            }
            start = colon + 1;
        }
        parts.resize(std::max<std::size_t>(parts.size(), 3));
        project.inputs.start = parts[0].empty() ? 0 : core::pythonInt(parts[0]);
        project.inputs.stop =
            parts[1].empty() ? std::nullopt : std::optional<std::int64_t>(core::pythonInt(parts[1]));
        project.inputs.step = parts[2].empty() ? 1 : core::pythonInt(parts[2]);
    }
    for (const std::string& assignment : a.set) {
        applyOverride(project, assignment);
    }
    if (resolve) {
        core::resolvePaths(project, core::currentDirectory());
    }
    return project;
}

// The console reporter of MDATK's command line (make_reporter, ConsoleProgress).
class ConsoleReporter {
public:
    ConsoleReporter(Console& console, bool quiet, bool verbose)
        : m_console(console), m_quiet(quiet), m_verbose(verbose) {}

    core::Reporter reporter() {
        return core::Reporter([this](const std::string& text) { logLine(text); },
                              [this](double fraction, const std::string& message) { progress(fraction, message); });
    }

private:
    void logLine(const std::string& text) {
        const bool detail = text.starts_with("    ") || text.starts_with("\t") || text.starts_with("|");
        if ((m_verbose || !detail) && !m_quiet) {
            m_console.err << (m_console.errIsTerminal ? "\r\033[K" : "") << text << '\n';
            m_console.err.flush();
        }
    }

    void progress(double fraction, const std::string& message) {
        if (m_quiet) {
            return;
        }
        const double now =
            std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
        if (fraction < 1.0 && fraction - m_last < 0.01 && now - m_lastTime < 2.0) {
            return;
        }
        m_last = fraction;
        m_lastTime = now;
        char percent[32];
        std::snprintf(percent, sizeof percent, "%5.1f", fraction * 100);
        const std::string line = "[" + std::string(percent) + "%] " + message;
        if (m_console.errIsTerminal) {
            m_console.err << "\r\033[K" << line << (fraction >= 1.0 ? "\n" : "");
        } else if (fraction >= 1.0) {
            // MDATK also prints when int(fraction * 10) != int(self.last * 10), compared after
            // self.last = fraction: never, so only final lines reach a file or a pipe
            m_console.err << line << '\n';
        }
        m_console.err.flush();
    }

    Console& m_console;
    bool m_quiet;
    bool m_verbose;
    double m_last = -1.0;
    double m_lastTime = 0.0;
};

// pipeline.PipelineError (a RuntimeError)
class PipelineError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

std::vector<fs::path> toPaths(const std::vector<std::string>& texts) {
    std::vector<fs::path> out;
    for (const std::string& t : texts) {
        out.push_back(core::pathFromUtf8(t));
    }
    return out;
}

// Pipeline.stamp(): the time stamp of a run's file names ("" when disabled)
std::string runStamp(const core::Project& project) {
    if (!project.output.timestamp_files) {
        return {};
    }
    const core::DateTime now = core::clockNow();
    char text[64];
    std::snprintf(text, sizeof text, "%04d-%02d-%02d_%02d%02d%02d", now.year, now.month, now.day, now.hour,
                  now.minute, now.second);
    return text;
}

int commandInit(CLI::App& command, const CommonArguments& common, const std::string& outputText, bool overwrite,
                Console& console) {
    const core::Project project = loadProject(command, common, false);  // relative paths stay relative to the file
    const fs::path output = core::pathFromUtf8(outputText);
    std::error_code ec;
    if (fs::exists(output, ec) && !overwrite) {
        console.err << outputText << " exists (use --overwrite)\n";
        return 1;
    }
    core::saveProject(output, project);
    console.out << "Wrote a project file with default settings to " << outputText << '\n';
    return 0;
}

int commandInfo(CLI::App& command, const CommonArguments& common, Console& console) {
    const core::Project project = loadProject(command, common, true);
    if (project.inputs.topology.empty()) {
        console.err << "info needs --project or --top\n";
        return 2;
    }
    const core::System system = io::loadSystem(core::pathFromUtf8(project.inputs.topology),
                                               toPaths(project.inputs.trajectories), project.inputs.dt_ps);
    Value info = io::systemSummary(system);
    core::Table groups;
    for (const core::ResidueGroup& g : project.groups) {
        Value entry;
        try {
            const core::AtomGroup atoms = core::select(system, "group:" + g.name, project.groups);
            const std::vector<std::int32_t> residues = atoms.residues();
            const core::Topology& t = atoms.topology();
            const auto label = [&t](std::int32_t r) {
                return t.resname(r) + std::to_string(t.resids[static_cast<std::size_t>(r)]);
            };
            entry = core::Table{{"atoms", static_cast<std::int64_t>(atoms.size())},
                                {"residues", static_cast<std::int64_t>(residues.size())},
                                {"first", label(residues.front())},
                                {"last", label(residues.back())}};
        } catch (const std::exception& error) {  // every problem is reported
            entry = core::Table{{"error", pythonErrorText(error)}};
        }
        bool replaced = false;  // a dict: a repeated name keeps its place and takes the new value
        for (auto& [name, value] : groups) {
            if (name == g.name) {
                value = entry;
                replaced = true;
            }
        }
        if (!replaced) {
            groups.emplace_back(g.name, std::move(entry));
        }
    }
    core::Array problems;
    for (const std::string& p : core::validate(project)) {
        problems.emplace_back(p);
    }
    core::Table& table = *info.get<core::Table>();
    table.emplace_back("groups", std::move(groups));
    table.emplace_back("problems", std::move(problems));
    console.out << jsonForConsole(info) << '\n';
    return 0;
}

int commandExportPdb(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    const CommonArguments& common = a.exportCommon;
    core::Project project = loadProject(command, common, true);
    ConsoleReporter consoleReporter(console, common.quiet, common.verbose);
    const core::Reporter reporter = consoleReporter.reporter();
    core::PDBParams& pdb = project.output.pdb;
    if (a.frameListOption->count() > 0) {
        pdb.frames = a.frameList;
    }
    if (a.timesOption->count() > 0) {
        pdb.times = a.times;
    }
    if (a.selectionOption->count() > 0) {
        pdb.selection = a.selection;
    }

    // Pipeline.check()
    const std::vector<std::string> problems = core::validate(project);
    if (!problems.empty()) {
        std::string text = "Invalid parameters:";
        for (const std::string& p : problems) {
            text += "\n  - " + p;
        }
        throw PipelineError(text);
    }
    for (const std::string& warning : core::layoutWarnings(project)) {
        reporter.log("WARNING: " + warning);
    }
    const std::string outDir = core::outDir(project);
    fs::create_directories(core::pathFromUtf8(outDir));

    // Pipeline.export_pdbs(): universe() (resolve_inputs, load_universe), frames, atoms, files
    if (project.cleaning.enabled) {
        throw PipelineError(
            "cleaning with cpptraj is not available in MDBridge yet: turn it off (--no-clean or "
            "cleaning.enabled = false) and give the cleaned files as --top / --traj");
    }
    reporter.progress(1.0, "inputs ready");
    core::System system = io::loadSystem(core::pathFromUtf8(project.inputs.topology),
                                         toPaths(project.inputs.trajectories), project.inputs.dt_ps);
    if (!system.trajectory) {
        throw core::ValueError("No trajectory loaded into Universe");
    }
    const std::vector<std::int64_t> frames = core::framesFromSpecs(
        *system.trajectory, pdb.frames, pdb.times, project.inputs.time_unit, project.inputs.dt_ps);
    if (frames.empty()) {
        throw PipelineError("PDB export enabled but no frames or times were given");
    }
    const core::AtomGroup atoms = core::resolveExportSelection(system, pdb.selection, project.groups);
    reporter.log("Exporting " + std::to_string(frames.size()) + " frame(s) (" + std::to_string(atoms.size()) +
                  " atoms) as PDB");
    const std::vector<fs::path> written = io::exportFrames(
        system, frames, atoms, core::pathFromUtf8(outDir), pdb.prefix,
        [&reporter](double fraction, const std::string& message) { reporter.progress(fraction, message); },
        runStamp(project));
    core::Array files;
    for (const fs::path& path : written) {
        files.emplace_back(core::pathToUtf8(path));
    }
    console.out << jsonForConsole(core::Table{{"pdbs", std::move(files)}}) << '\n';
    return 0;
}

int commandDensityInfo(const std::string& map, Console& console) {
    console.out << jsonForConsole(io::densityInfo(core::pathFromUtf8(map))) << '\n';
    return 0;
}

}  // namespace

// ---- the command line ------------------------------------------------------------------------

Console standardConsole() {
#ifdef _WIN32
    const bool tty = _isatty(_fileno(stderr)) != 0;
#else
    const bool tty = isatty(fileno(stderr)) != 0;
#endif
    return Console{std::cout, std::cerr, tty};
}

CoreCommands::CoreCommands(CLI::App& app) : m_args(std::make_unique<Arguments>()) {
    Arguments& a = *m_args;

    a.init = app.add_subcommand("init", "write a project file with default settings");
    a.init->alias("init-params");
    a.init->add_option("output", a.initOutput, "project file to write (.mdatk.toml)")->required();
    a.init->add_flag("--overwrite", a.overwrite, "replace an existing file");
    addCommon(*a.init, a.initCommon);

    a.info = app.add_subcommand("info", "summarise a system and check residue groups");
    addCommon(*a.info, a.infoCommon);

    a.exportPdb = app.add_subcommand("export-pdb", "write selected frames as PDB files");
    addCommon(*a.exportPdb, a.exportCommon);
    a.frameListOption = a.exportPdb->add_option("--frame-list", a.frameList, "frames, e.g. \"0, 10, 20-50:10, last\"");
    a.timesOption = a.exportPdb->add_option("--times", a.times, "times, e.g. \"1.5ns, 2000ps\"");
    a.selectionOption = a.exportPdb->add_option("--selection", a.selection, "group names (comma separated) or a selection");

    a.density = app.add_subcommand("density", "density maps (MRC / CCP4): header info");
    a.density->require_subcommand(0, 1);
    a.densityInfo = a.density->add_subcommand("info", "print a map's grid, voxel size, origin and statistics");
    a.densityInfo->add_option("map", a.map, "MRC / CCP4 map (.mrc, .map, .ccp4; gzip-compressed too)")->required();
}

CoreCommands::~CoreCommands() = default;

std::optional<int> CoreCommands::run(Console& console) {
    Arguments& a = *m_args;
    const bool any = a.init->parsed() || a.info->parsed() || a.exportPdb->parsed() || a.density->parsed();
    if (!any) {
        return std::nullopt;
    }
    if (a.density->parsed() && !a.densityInfo->parsed()) {
        console.out << a.density->help();
        return 2;
    }
    // logging.basicConfig(level=WARNING, format="%(levelname)s: %(message)s")
    const core::LogSink previous = core::setLogSink([&console](core::LogLevel level, std::string_view message) {
        if (level == core::LogLevel::Warning || level == core::LogLevel::Error) {
            console.err << core::logLevelName(level) << ": " << message << '\n';
            console.err.flush();
        }
    });
    struct RestoreSink {
        const core::LogSink& sink;
        ~RestoreSink() { core::setLogSink(sink); }
    } restore{previous};
    try {
        if (a.init->parsed()) {
            return commandInit(*a.init, a.initCommon, a.initOutput, a.overwrite, console);
        }
        if (a.info->parsed()) {
            return commandInfo(*a.info, a.infoCommon, console);
        }
        if (a.exportPdb->parsed()) {
            return commandExportPdb(*a.exportPdb, a, console);
        }
        return commandDensityInfo(a.map, console);
    } catch (const core::Cancelled&) {
        console.err << "cancelled\n";
        return 130;
    } catch (const std::exception& error) {
        console.out.flush();
        console.err << "error: " << pythonErrorText(error) << '\n';
        return 1;
    }
}

int runCoreCommandLine(int argc, const char* const* argv, Console& console) {
    CLI::App app("MDBridge: the MDATK commands of the core and io libraries (no interface)",
                 argc > 0 ? core::pathBasename(argv[0]) : "mdbridge");
    app.require_subcommand(0, 1);
    CoreCommands commands(app);
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& error) {
        const int code = app.exit(error, console.out, console.err);
        return code == 0 ? 0 : 2;
    }
    if (std::optional<int> code = commands.run(console)) {
        return *code;
    }
    console.out << app.help();
    return 2;
}

// ---- Python's texts --------------------------------------------------------------------------

std::string jsonForConsole(const Value& value) {
    const std::string text = core::toJson(value, 2);
    std::string out;
    out.reserve(text.size());
    const auto escape = [&out](std::uint32_t code) {
        char buffer[16];
        std::snprintf(buffer, sizeof buffer, "\\u%04x", static_cast<unsigned>(code));
        out += buffer;
    };
    for (char32_t c : core::decodeUtf8Surrogateescape(text)) {
        const auto code = static_cast<std::uint32_t>(c);
        if (code < 0x7f) {
            out += static_cast<char>(code);
        } else if (code < 0x10000) {
            escape(code);
        } else {
            escape(0xd800 + ((code - 0x10000) >> 10));
            escape(0xdc00 + ((code - 0x10000) & 0x3ff));
        }
    }
    return out;
}

std::string pythonErrorText(const std::exception& error) {
    if (dynamic_cast<const fs::filesystem_error*>(&error) != nullptr) {
        // io's FileNotFoundError("topology not found: <path>") carries its message alone
        const std::string what = error.what();
        if (what.starts_with("topology not found: ") || what.starts_with("trajectory not found: ")) {
            return what;
        }
    }
    if (dynamic_cast<const PipelineError*>(&error) != nullptr) {
        return error.what();
    }
    const std::string summary = core::exceptionSummary(error);  // "<Python class>: <str(exc)>"
    const std::size_t colon = summary.find(": ");
    return colon == std::string::npos ? std::string() : summary.substr(colon + 2);
}

}  // namespace mdbridge::app
