#pragma once
// MDATK project files (.mdatk.toml) and per-run parameter logs (MDATK_params_log.toml).
//
// The specification is mdatk/projectfile.py of MDATK v0.0.4.1 (read(), log_path(), read_log()):
//
//     format = "mdatk-project"
//     notes = "..."     # free-text project notes
//     [inputs]          # analysis settings ...
//     [workflow]        # planner steps ...
//     # ---- step log ----
//     [[log]]           # one entry per step run (appended, newest last)
//
// Once a project file would grow past 50 MB its log moves to <name>.log.toml next to it and
// log_file names that file. Milestone 1 reads these files; the typed parameter model and the
// commented writer follow with the `core` milestone.

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

inline constexpr std::string_view kProjectFormat = "mdatk-project";
inline constexpr std::int64_t kProjectFormatVersion = 1;
inline constexpr std::string_view kParamsLogFormat = "mdatk-params-log";
inline constexpr std::string_view kProjectSuffix = ".mdatk.toml";
inline constexpr std::string_view kLogMarker = "# ---- step log (entries are appended; newest last) ----";

// A file that is not a readable project; the message names the file, as MDATK's do.
class ProjectError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class FileKind { Project, ParamsLog };

struct ProjectFile {
    std::filesystem::path path;
    FileKind kind = FileKind::Project;
    // Entries about the file rather than the project: format, format_version, mdatk_version,
    // saved, log_file (project file); format, mdatk_version, created, run (parameter log).
    Value header = Table{};
    // The project settings exactly as mdatk.projectfile.read() returns them: meta entries and
    // the step log removed; for a parameter log, its "parameters" table.
    Value settings = Table{};
};

// Reads a project file, or a run's parameter log (which opens like a project, to replay the run).
// Throws ProjectError for JSON presets, unreadable files, invalid TOML and other formats.
ProjectFile readProject(const std::filesystem::path& path);

// True for a run's parameter log: format = "mdatk-params-log" in its first 4096 characters.
bool isParamsLog(const std::filesystem::path& path);

// The log_file entry of a project ("" when the log is at the end of the project file itself).
std::string logFileName(const std::filesystem::path& project);

// Where the step log of `project` is written: the project file or the file log_file names.
std::filesystem::path logPath(const std::filesystem::path& project);

// The step log entries of `project`, oldest first. A missing or unreadable log gives no entries
// (MDATK warns and carries on); `warning`, if given, receives the reason.
Array readLog(const std::filesystem::path& project, std::string* warning = nullptr);

// TOML text as plain data, keys in the order of the text.
// Throws ProjectError("<source>: not valid TOML: <reason> (at line L, column C)").
Value parseToml(std::string_view text, std::string_view source);

}  // namespace mdbridge::core
