#include "command_support.hpp"

#include <algorithm>
#include <filesystem>
#include <system_error>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::app::detail {

namespace fs = std::filesystem;

namespace {

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

}  // namespace

CLI::Option* lastWins(CLI::Option* option) { return option->multi_option_policy(CLI::MultiOptionPolicy::TakeLast); }

void addCommon(CLI::App& command, CommonArguments& a) {
    lastWins(command.add_option("--project,-p,--params", a.params,
                                "project file (.mdatk.toml) or a run's MDATK_params_log.toml"));
    lastWins(command.add_option("--top", a.top, "topology file (.parm7, .pdb, ...)"));
    // nargs="+": a repeated --traj replaces the earlier list (each occurrence starts afresh)
    command.add_option("--traj", a.traj, "trajectory file(s) (.nc, .dcd, ...)")->trigger_on_parse();
    lastWins(command.add_option("--frames", a.frames, "frame selection start:stop:step (0-based, stop exclusive)"));
    a.cleanOption = command.add_flag("--clean", a.clean, "enable cpptraj cleaning");
    a.noCleanOption = command.add_flag("--no-clean", a.noClean, "disable cleaning");
    lastWins(command.add_option("--cpptraj", a.cpptraj, "path to the cpptraj executable"));
    lastWins(command.add_option(
        "--conda-env", a.condaEnv,
        "conda environment providing cpptraj (e.g. AmberTools26); cpptraj then runs via `conda run`"));
    lastWins(command.add_option("--out-dir,-o", a.outDir, "output directory"));
    lastWins(command.add_option("--cache-dir", a.cacheDir, "cache directory (default <out-dir>/MDATK_cache)"));
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

// load_project(args, resolve)
core::Project loadProject(const CLI::App& command, const CommonArguments& a, bool resolve) {
    if (!a.params.empty()) {  // Project.load's open(): Python's FileNotFoundError text
        std::error_code ec;
        if (!fs::exists(core::pathFromUtf8(a.params), ec)) {
            throw fs::filesystem_error("open", core::pathFromUtf8(a.params),
                                       std::make_error_code(std::errc::no_such_file_or_directory));
        }
    }
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
        // MDATK reports int()'s "invalid literal for int() with base 10: 'x'"; MDBridge names the option
        const auto number = [&a](const std::string& text) {
            try {
                return core::pythonInt(text);
            } catch (const core::ValueError&) {
                throw core::ValueError("--frames: " + core::pythonStringRepr(a.frames) +
                                       " is not a frame number or range");
            }
        };
        project.inputs.start = parts[0].empty() ? 0 : number(parts[0]);
        project.inputs.stop = parts[1].empty() ? std::nullopt : std::optional<std::int64_t>(number(parts[1]));
        project.inputs.step = parts[2].empty() ? 1 : number(parts[2]);
    }
    for (const std::string& assignment : a.set) {
        applyOverride(project, assignment);
    }
    if (resolve) {
        core::resolvePaths(project, core::currentDirectory());
    }
    return project;
}

std::string commandText(const std::vector<std::string>& arguments, const std::string& group) {
    if (arguments.empty() || arguments.front() != group) {
        return {};
    }
    std::string out;
    for (const std::string& a : arguments) {
        out += (out.empty() ? "" : " ") + a;
    }
    return out;
}

}  // namespace mdbridge::app::detail
