#pragma once
// Internals shared by the project-file sources (project_file.cpp, project_io.cpp, step_log.cpp,
// settings.cpp); not installed.

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

namespace mdbridge::core::detail {

// A path as UTF-8 text, for messages ("{path}" in MDATK's messages is the path as given).
std::string displayPath(const std::filesystem::path& path);

// os.path.abspath of a path, as UTF-8 text (POSIX string rules; std::filesystem on Windows).
std::string absolutePath(const std::filesystem::path& path);

// projectfile._reject_json: throws ProjectError for a .json path (any case).
void rejectJson(const std::filesystem::path& path);

// The bytes of a file. Throws std::filesystem::filesystem_error.
std::string readFileBytes(const std::filesystem::path& path);

// The whole file read as Python's text mode reads it: strict UTF-8 (ProjectError naming the
// first undecodable byte otherwise) and universal newlines (CR LF and CR become LF).
std::string readTextFile(const std::filesystem::path& path);

// Writes (or appends) text as Python's text mode does: LF becomes the platform's line end.
// Throws std::filesystem::filesystem_error.
void writeTextFile(const std::filesystem::path& path, std::string_view text, bool append = false);

// Writes <path>.tmp and renames it over `path` (os.replace).
void replaceFile(const std::filesystem::path& path, std::string_view text);

// projectfile._locked(path): the process-wide lock, and on POSIX flock() on .<name>.lock next to
// the project (only the process-wide lock when that file cannot be opened). Not re-entrant across
// projects in one thread: never nest two of them.
class ProjectLock {
public:
    explicit ProjectLock(const std::filesystem::path& path);
    ~ProjectLock();
    ProjectLock(const ProjectLock&) = delete;
    ProjectLock& operator=(const ProjectLock&) = delete;

private:
    struct State;
    std::unique_ptr<State> m_state;
};

}  // namespace mdbridge::core::detail
