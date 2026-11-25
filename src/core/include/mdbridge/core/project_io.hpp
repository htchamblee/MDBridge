#pragma once
// Loading and saving projects (Project.load / Project.save) and writing parameter logs: MDATK's
// mdatk/projectfile.py.
//
// A saved project file is the commented settings text (projectText) followed, when the old file
// had one, by its step log copied byte for byte. The log starts at the marker line that ends complete
// TOML settings, so the marker pasted into the notes is no log (MDATK 0.1.0; DIFFERENCES PRJ-10). Saves are atomic (<path>.tmp, then rename) and
// exclude each other, and appends to the step log, across threads and processes: a global mutex
// and, on Linux, flock() on the lock file .<name>.lock next to the project (the lock MDATK takes,
// so MDATK and MDBridge can work on the same project).
//
// Text files are UTF-8 with the platform's line ends (LF; CR LF on Windows), as MDATK writes them.

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// Project.load: readProject(path) (projects and parameter logs), then projectFromPlain with the
// paths resolved against the file's folder; Project::mdatkVersion is the file's mdatk_version.
// Throws ProjectError.
Project loadProject(const std::filesystem::path& path, std::vector<std::string>* warnings = nullptr);

// Project.save. Throws ProjectError for a .json path ("<path>: JSON presets are no longer
// supported; MDATK 0.0.4.1 uses TOML project files (.mdatk.toml)") and for a parameter log
// ("<path> is the parameter log of a run; save the project under another name"), and
// std::filesystem::filesystem_error when the file cannot be written. The step log and log_file of
// an existing file are kept.
void saveProject(const std::filesystem::path& path, const Project& project);

// The settings part of a project file: the 3-line header, format, format_version, mdatk_version
// (writtenMdatkVersion), saved (clockNow()), log_file, then toPlain(project), every entry with its
// comment.
std::string projectText(const Project& project, std::string_view logFile = {});

// The mdatk_version a save writes: project.mdatkVersion when it is a later version than
// kMdatkReference ("0.1.1" is kept), else kMdatkReference ("0.0.4.1"; also for "" and texts that
// are no version). Versions compare component by component (split at '.'), each by its leading number and
// then the rest as text ("0.0.3d" > "0.0.3"), a missing component lowest ("0.0.4" < "0.0.4.1").
std::string writtenMdatkVersion(const Project& project);

// config_doc(path, stack): the comment above an entry (or table header) of a project file.
// `path` holds the keys from the root without list indices, `stack` the enclosing tables, root
// first (see TomlDocFunction). Model fields, meta entries, panel options, node files and planner
// settings (from the step types) are commented; anything else is not.
std::optional<std::string> projectDoc(const std::vector<std::string>& path,
                                      const std::vector<const Value*>& stack);

// write_params_log: the per-run parameter log (MDATK_params_log.toml) with `run` (what was run,
// written as given) and the complete settings. Atomic, not locked.
void writeParamsLog(const std::filesystem::path& path, const Project& project, const Value& run);

}  // namespace mdbridge::core
