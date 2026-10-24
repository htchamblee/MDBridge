#include "mdbridge/core/project_io.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdio>
#include <mutex>
#include <system_error>
#include <utility>

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/step_types.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/core/toml_writer.hpp"
#include "mdbridge/version.hpp"
#include "project_internal.hpp"

#ifdef _WIN32
#include <cwchar>
#else
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#endif

namespace mdbridge::core {

namespace fs = std::filesystem;

namespace detail {

namespace {

std::recursive_mutex& projectMutex() {
    static std::recursive_mutex mutex;
    return mutex;
}

bool endsWithIgnoringCase(std::string_view text, std::string_view suffix) {
    if (text.size() < suffix.size()) {
        return false;
    }
    return std::equal(suffix.begin(), suffix.end(), text.end() - static_cast<std::ptrdiff_t>(suffix.size()),
                      [](char a, char b) {
                          auto lower = [](char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; };
                          return lower(a) == lower(b);
                      });
}

struct FileCloser {
    void operator()(std::FILE* file) const { std::fclose(file); }
};
using File = std::unique_ptr<std::FILE, FileCloser>;

File openFile(const fs::path& path, const char* mode) {
#ifdef _WIN32
    std::wstring wmode(mode, mode + std::char_traits<char>::length(mode));
    return File(_wfopen(path.c_str(), wmode.c_str()));
#else
    return File(std::fopen(path.c_str(), mode));
#endif
}

[[noreturn]] void throwFileError(const char* what, const fs::path& path, int error) {
    throw fs::filesystem_error(what, path, std::error_code(error ? error : EIO, std::generic_category()));
}

}  // namespace

std::string readFileBytes(const fs::path& path) {
    errno = 0;
    File file = openFile(path, "rb");
    if (!file) {
        throwFileError("cannot open the file", path, errno);
    }
    std::string bytes;
    char buffer[65536];
    std::size_t n = 0;
    while ((n = std::fread(buffer, 1, sizeof buffer, file.get())) > 0) {
        bytes.append(buffer, n);
    }
    if (std::ferror(file.get())) {
        throwFileError("cannot read the file", path, errno);
    }
    return bytes;
}

std::string displayPath(const fs::path& path) { return pathToUtf8(path); }

std::string absolutePath(const fs::path& path) {
#ifdef _WIN32
    std::error_code ec;
    const fs::path absolute = fs::absolute(path, ec);
    return pathToUtf8((ec ? path : absolute).lexically_normal());
#else
    return posixpath::abspath(pathToUtf8(path), currentDirectory());
#endif
}

void rejectJson(const fs::path& path) {
    if (endsWithIgnoringCase(displayPath(path), ".json")) {
        throw ProjectError(displayPath(path) + ": JSON presets are no longer supported; MDBridge uses TOML "
                                               "project files (" + std::string(kProjectSuffix) + ")");
    }
}

std::string readTextFile(const fs::path& path) {
    const std::string bytes = readFileBytes(path);
    const std::u32string codes = decodeUtf8(bytes);
    if (encodeUtf8(codes) != bytes) {  // a byte was replaced: Python's strict decoding fails
        std::size_t offset = 0;
        for (const char32_t c : codes) {
            std::string one;
            appendUtf8(one, c);
            if (bytes.compare(offset, one.size(), one) != 0) {
                break;
            }
            offset += one.size();
        }
        char byte[8];
        std::snprintf(byte, sizeof byte, "0x%02x", static_cast<unsigned char>(bytes[offset]));
        throw ProjectError(displayPath(path) + ": 'utf-8' codec can't decode byte " + byte + " in position " +
                           std::to_string(offset));
    }
    std::string text;
    text.reserve(bytes.size());
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        if (bytes[i] == '\r') {
            text += '\n';
            if (i + 1 < bytes.size() && bytes[i + 1] == '\n') {
                ++i;
            }
        } else {
            text += bytes[i];
        }
    }
    return text;
}

void writeTextFile(const fs::path& path, std::string_view text, bool append) {
#ifdef _WIN32
    std::string data;
    data.reserve(text.size() + text.size() / 32);
    for (const char c : text) {
        if (c == '\n') {
            data += '\r';
        }
        data += c;
    }
#else
    const std::string_view data = text;
#endif
    errno = 0;
    File file = openFile(path, append ? "ab" : "wb");
    if (!file) {
        throwFileError("cannot open the file for writing", path, errno);
    }
    if (std::fwrite(data.data(), 1, data.size(), file.get()) != data.size() || std::fflush(file.get()) != 0) {
        throwFileError("cannot write the file", path, errno);
    }
}

void replaceFile(const fs::path& path, std::string_view text) {
    fs::path tmp = path;
    tmp += ".tmp";
    writeTextFile(tmp, text);
    fs::rename(tmp, path);
}

struct ProjectLock::State {
    std::unique_lock<std::recursive_mutex> thread;
    int fd = -1;
};

ProjectLock::ProjectLock(const fs::path& path) : m_state(std::make_unique<State>()) {
    m_state->thread = std::unique_lock<std::recursive_mutex>(projectMutex());
#ifndef _WIN32
    // MDATK has no lock file on Windows (no fcntl there) and neither has MDBridge.
    const fs::path lock = pathFromUtf8(posixpath::dirname(absolutePath(path))) /
                          ("." + pathToUtf8(path.filename()) + ".lock");
    const int fd = ::open(lock.c_str(), O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC, 0666);
    if (fd < 0) {
        return;  // a read-only folder: the thread lock only, as MDATK
    }
    int rc = 0;
    while ((rc = ::flock(fd, LOCK_EX)) != 0 && errno == EINTR) {
    }
    if (rc != 0) {
        ::close(fd);
        return;
    }
    m_state->fd = fd;
#endif
}

ProjectLock::~ProjectLock() {
#ifndef _WIN32
    if (m_state->fd >= 0) {
        ::flock(m_state->fd, LOCK_UN);
        ::close(m_state->fd);
    }
#endif
}

}  // namespace detail

namespace {

constexpr std::array<std::pair<std::string_view, std::string_view>, 5> kMetaDocs{{
    {"format", "File type marker of an MDATK project; do not change."},
    {"format_version", "Version of the project file layout."},
    {"mdatk_version", "Molecular Dynamic Analysis Tool Kit (MDATK) version that wrote this file."},
    {"saved", "When the settings were last saved."},
    {"log_file",
     "File next to the project that holds its step log (\"\": the log is at the end of this file; it moves "
     "out once the project file would grow too large)."},
}};

constexpr std::array<std::pair<std::string_view, std::string_view>, 13> kPanelOptionDocs{{
    {"mode",
     "Contact map: 'instant' (evolves with the movie), 'frame' (one fixed frame) or 'frequency' (fraction of "
     "frames in contact)."},
    {"frame", "Contact map in 'frame' mode: the frame shown (0-based)."},
    {"transpose", "Contact map: true puts the first group (e.g. the aptamer) along x and the second along y."},
    {"progressive", "Time series: draw the curves up to the current movie frame only."},
    {"smooth", "Time series: running mean over this many frames (0 or 1 = raw data)."},
    {"target", "Time series: value drawn as a horizontal reference line."},
    {"vmin", "Heatmap: lower end of the colour scale."},
    {"vmax", "Heatmap: upper end of the colour scale (0 = the data maximum)."},
    {"max_ticks", "Heatmap: at most this many residue labels on the axis."},
    {"residues",
     "Sequence / secondary structure: residues shown, e.g. 'A:1-20, B:5-30' (empty: all residues of the "
     "structure view)."},
    {"per_row", "Sequence: residues per row (0: as many as fit)."},
    {"source",
     "Secondary structure: 'dynamic' (pairs re-detected every n frames), 'first' (the first frame's pairs) or "
     "'dotbracket' (the 'dotbracket' option)."},
    {"dotbracket", "Secondary structure: dot-bracket of the shown residues ('&' or a space between chains)."},
}};

constexpr std::array<std::string_view, 3> kProjectHeader{
    "Molecular Dynamic Analysis Tool Kit (MDATK) project file.",
    "Plain TOML: every entry is explained by the comment above it. Edit with care; MDATK rewrites this",
    "file (keeping the step log) when the project is saved."};

template <std::size_t N>
std::optional<std::string> lookup(const std::array<std::pair<std::string_view, std::string_view>, N>& table,
                                  std::string_view key) {
    for (const auto& [name, text] : table) {
        if (name == key) {
            return std::string(text);
        }
    }
    return std::nullopt;
}

bool startsWith(const std::vector<std::string>& path, std::string_view a, std::string_view b,
                std::string_view c) {
    return path.size() >= 3 && path[0] == a && path[1] == b && path[2] == c;
}

// projectfile._free_dict_doc / _settings_doc: entries below the free-form tables
std::optional<std::string> freeDictDoc(const std::vector<std::string>& path, const std::vector<const Value*>& stack) {
    if (startsWith(path, "layout", "panels", "options") && path.size() == 4) {
        if (auto text = lookup(kPanelOptionDocs, path[3])) {
            return text;
        }
        return "Panel option '" + path[3] + "' (see the panel type's description).";
    }
    if (path.size() == 4 && (path[1] == "nodes" || path[1] == "scratch") && path[2] == "files") {
        return "File chosen by hand for the step's '" + path[3] +
               "' input (paths relative to this file; ';' separates several).";
    }
    if (!startsWith(path, "workflow", "nodes", "settings") && !startsWith(path, "workflow", "scratch", "settings")) {
        return std::nullopt;
    }
    // The step is the table two levels up. Below the settings themselves (records of a repeated
    // section, nested tables) that is the settings table, which names no step type: MDATK leaves
    // those entries uncommented.
    const Value* node = stack.size() >= 2 ? stack[stack.size() - 2] : nullptr;
    const Value* tool = node ? node->find("tool") : nullptr;
    const std::string* name = tool ? tool->get<std::string>() : nullptr;
    return settingComment(name ? *name : std::string_view(), path[3]);
}

// The table of the dataclass model below `path`, as config_doc walks it: true for a free-form
// table (PanelSpec.options, PlanNode.settings, PlanNode.files) at path[0..2].
bool belowFreeDict(const std::vector<std::string>& path) {
    return path.size() >= 4 && (startsWith(path, "layout", "panels", "options") ||
                                startsWith(path, "workflow", "nodes", "settings") ||
                                startsWith(path, "workflow", "nodes", "files") ||
                                startsWith(path, "workflow", "scratch", "settings") ||
                                startsWith(path, "workflow", "scratch", "files"));
}

std::string joined(const std::vector<std::string>& path, std::size_t first = 0) {
    std::string out;
    for (std::size_t i = first; i < path.size(); ++i) {
        if (i > first) {
            out += '.';
        }
        out += path[i];
    }
    return out;
}

Table metaEntries(std::string_view logFile) {
    return Table{{"format", std::string(kProjectFormat)},
                 {"format_version", kProjectFormatVersion},
                 {"mdatk_version", std::string(kMdatkReference)},
                 {"saved", clockNow()},
                 {"log_file", std::string(logFile)}};
}

// _existing_log: (the log section from its marker on, the project's log_file)
std::pair<std::string, std::string> existingLog(const fs::path& path) {
    std::error_code ec;
    if (!fs::is_regular_file(path, ec)) {
        return {};
    }
    const std::string text = detail::readTextFile(path);
    std::string logFile;
    const std::vector<std::string> lines = pythonSplitlines(text);
    for (std::size_t i = 0; i < lines.size() && i < 60; ++i) {
        if (lines[i].rfind("log_file = ", 0) == 0) {
            try {
                const Value entry = parseToml(lines[i], detail::displayPath(path));
                // MDATK would write any other type back verbatim and fail later; "" keeps the log here
                logFile = entry.find("log_file") ? entry.find("log_file")->stringOr() : std::string();
            } catch (const ProjectError&) {
                logFile.clear();
            }
            break;
        }
    }
    const std::size_t at = text.find("\n" + std::string(kLogMarker));
    return {at == std::string::npos ? std::string() : text.substr(at + 1), logFile};
}

}  // namespace

std::optional<std::string> projectDoc(const std::vector<std::string>& path, const std::vector<const Value*>& stack) {
    if (path.empty()) {
        return std::nullopt;
    }
    if (path.size() == 1) {
        if (auto text = lookup(kMetaDocs, path[0])) {
            return text;
        }
    }
    if (belowFreeDict(path)) {
        return freeDictDoc(path, stack);
    }
    if (const auto text = modelComment(joined(path))) {
        return std::string(*text);
    }
    return std::nullopt;
}

std::string projectText(const Project& project, std::string_view logFile) {
    Table data = metaEntries(logFile);
    Value settings = toPlain(project);
    for (auto& entry : *settings.get<Table>()) {
        data.push_back(std::move(entry));
    }
    const std::vector<std::string> header(kProjectHeader.begin(), kProjectHeader.end());
    return dumpToml(Value(std::move(data)), projectDoc, header);
}

Project loadProject(const fs::path& path, std::vector<std::string>* warnings) {
    const ProjectFile file = readProject(path);
    return projectFromPlain(file.settings, posixpath::dirname(detail::absolutePath(path)), warnings);
}

void saveProject(const fs::path& path, const Project& project) {
    detail::rejectJson(path);
    if (isParamsLog(path)) {
        throw ProjectError(detail::displayPath(path) +
                           " is the parameter log of a run; save the project under another name");
    }
    fs::create_directories(pathFromUtf8(posixpath::dirname(detail::absolutePath(path))));
    detail::ProjectLock lock(path);
    const auto [logPart, logFile] = existingLog(path);
    std::string text = projectText(project, logFile);
    if (!logPart.empty()) {
        text += "\n";
        text += logPart;
    }
    detail::replaceFile(path, text);
}

void writeParamsLog(const fs::path& path, const Project& project, const Value& run) {
    Table data{{"format", std::string(kParamsLogFormat)},
               {"mdatk_version", std::string(kMdatkReference)},
               {"created", clockNow()},
               {"run", run},
               {"parameters", toPlain(project)}};
    auto doc = [](const std::vector<std::string>& p, const std::vector<const Value*>& stack) -> std::optional<std::string> {
        if (p.size() == 1) {
            if (p[0] == "format") {
                return "Per-run parameter log of MDATK; open it as a project to replay the run.";
            }
            if (p[0] == "created") {
                return "When the run finished.";
            }
            if (p[0] == "run") {
                return "What was run: environment, input files, cleaning, cache keys, outputs, citations.";
            }
            if (p[0] == "parameters") {
                return "The complete settings of the run.";
            }
            return lookup(kMetaDocs, p[0]);
        }
        if (!p.empty() && p[0] == "parameters") {
            const std::vector<std::string> rest(p.begin() + 1, p.end());
            const std::vector<const Value*> inner(stack.empty() ? stack.begin() : stack.begin() + 1, stack.end());
            return projectDoc(rest, inner);
        }
        return std::nullopt;
    };
    detail::replaceFile(path, dumpToml(Value(std::move(data)), doc,
                                       {"Molecular Dynamic Analysis Tool Kit (MDATK) parameter log of one run."}));
}

}  // namespace mdbridge::core
