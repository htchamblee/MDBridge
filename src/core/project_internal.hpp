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

// os.path.abspath of a path against the current directory, as UTF-8 text: posixpath's rules; on
// Windows ntpath::abspath (MDBridge's Windows rules, ntpath.hpp), in generic form ("C:/p/x.toml").
std::string absolutePath(const std::filesystem::path& path);

// os.path.dirname(os.path.abspath(path)): the folder of a project, as UTF-8 text (posixpath; on
// Windows ntpath::dirname of the generic absolute path, "C:/projects/p"), the base against which
// resolvePaths makes the project's relative paths absolute.
std::string absoluteDirectory(const std::filesystem::path& path);

// projectfile._reject_json: throws ProjectError for a .json path (any case).
void rejectJson(const std::filesystem::path& path);

// The bytes of a file. Throws std::filesystem::filesystem_error.
std::string readFileBytes(const std::filesystem::path& path);

// The whole file read as Python's text mode reads it: strict UTF-8 (ProjectError naming the
// first undecodable byte otherwise) and universal newlines (CR LF and CR become LF).
std::string readTextFile(const std::filesystem::path& path);

// The whole file read as Python's text mode with errors="replace" reads it: bytes that are not
// UTF-8 become U+FFFD, universal newlines. Throws std::filesystem::filesystem_error.
std::string readTextFileReplacing(const std::filesystem::path& path);

// projectfile._log_start of MDATK 0.1.0: where the step log starts in the text of a project file,
// the position of the line break before the marker line that ends the settings (npos: no log).
// The marker's text may also be part of the settings (pasted into the notes or a custom script):
// there it stands inside a string, so the text before it is no complete TOML and it does not
// start the log. When no marker ends complete settings and the whole text is valid TOML, there is
// no log; when the whole text is not valid TOML either, the first marker (as MDATK 0.0.4.1).
std::size_t logStart(std::string_view text);

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
