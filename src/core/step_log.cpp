#include "mdbridge/core/step_log.hpp"

#include <array>
#include <atomic>
#include <new>
#include <system_error>

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/core/toml_writer.hpp"
#include "project_internal.hpp"

namespace mdbridge::core {

namespace fs = std::filesystem;

namespace {

std::atomic<std::uintmax_t> g_splitBytes{kDefaultLogSplitBytes};

constexpr std::array<std::pair<std::string_view, std::string_view>, 11> kLogDocs{{
    {"time", "When the step started."},
    {"step",
     "What ran: 'render', 'workflow write', 'pdb clean', or the planner step number of a generated script "
     "that MDATK ran (e.g. '03')."},
    {"title", "Short human-readable description."},
    {"status", "ok, failed or cancelled."},
    {"duration_s", "Wall time in seconds."},
    {"command", "Command line that was run (for generated scripts)."},
    {"folder", "Working folder of the command."},
    {"files", "Files written by the step."},
    {"summary", "Result in a few lines."},
    {"citations", "Citation keys of the resources the step used (see the BibTeX export)."},
    {"output", "Program output (the last 20000 characters)."},
}};

constexpr std::string_view kEntryDoc =
    "One entry per step run in this project: analyses, renders, generated scripts that MDATK ran, PDB "
    "clean-ups, density fits.";

Value stringList(const std::vector<std::string>& items) {
    Array out;
    out.reserve(items.size());
    for (const std::string& item : items) {
        out.emplace_back(item);
    }
    return out;
}

std::string lastCharacters(std::string_view text, std::size_t count) {
    const std::u32string codes = decodeUtf8(text);
    return encodeUtf8(std::u32string_view(codes).substr(codes.size() - count));
}

// _split_log without the lock (append_log holds it already)
fs::path splitLogUnlocked(const fs::path& project) {
    const std::string base = pathToUtf8(project.filename());
    std::string stem;
    if (base.size() >= kProjectSuffix.size() &&
        base.compare(base.size() - kProjectSuffix.size(), kProjectSuffix.size(), kProjectSuffix) == 0) {
        stem = base.substr(0, base.size() - kProjectSuffix.size());
    } else {
        stem = posixpath::splitext(base).first;
    }
    const std::string name = stem + ".log.toml";
    const fs::path sidecar = pathFromUtf8(posixpath::join(posixpath::dirname(detail::absolutePath(project)), name));

    const std::string text = detail::readTextFile(project);
    const std::size_t at = text.find("\n" + std::string(kLogMarker));
    std::string head = at == std::string::npos ? text : text.substr(0, at + 1);
    const std::string logPart = at == std::string::npos ? std::string(kLogMarker) + "\n" : text.substr(at + 1);

    std::error_code ec;
    const bool empty = !fs::exists(sidecar, ec) || fs::file_size(sidecar, ec) == 0;
    std::string append;
    if (empty) {
        append = "# Step log of the MDATK project " + base + " (moved here when the project file grew too large).\n\n";
    }
    append += logPart;
    detail::writeTextFile(sidecar, append, true);

    const std::string from = "\nlog_file = \"\"\n";
    if (const std::size_t pos = head.find(from); pos != std::string::npos) {
        head.replace(pos, from.size(), "\nlog_file = " + tomlQuote(name) + "\n");
    }
    detail::replaceFile(project, head);
    logInfo("project log moved to " + detail::displayPath(sidecar));
    return sidecar;
}

bool containsMarker(const fs::path& target) {
    // MDATK decodes with errors="replace"; the marker is ASCII, which that decoding never alters
    return detail::readFileBytes(target).find(kLogMarker) != std::string::npos;
}

}  // namespace

std::uintmax_t logSplitBytes() { return g_splitBytes.load(); }
void setLogSplitBytes(std::uintmax_t bytes) { g_splitBytes.store(bytes); }

LogEntry::LogEntry(std::string step_, std::string title_, std::vector<std::string> citations_)
    : step(std::move(step_)), title(std::move(title_)), time(clockNow()), citations(std::move(citations_)) {}

Value LogEntry::plain() const {
    Table data;
    auto text = [&](std::string_view key, const std::string& value) {
        if (!value.empty()) {
            data.emplace_back(std::string(key), value);
        }
    };
    auto list = [&](std::string_view key, const std::vector<std::string>& value) {
        if (!value.empty()) {
            data.emplace_back(std::string(key), stringList(value));
        }
    };
    text("step", step);
    text("title", title);
    text("status", status);
    data.emplace_back("time", time);
    data.emplace_back("duration_s", pythonRound(duration_s, 2));
    text("command", command);
    text("folder", folder);
    list("files", files);
    text("summary", summary);
    list("citations", citations);
    if (characterCount(output) > kLogOutputChars) {
        data.emplace_back("output", "[...]\n" + lastCharacters(output, kLogOutputChars));
    } else {
        text("output", output);
    }
    return data;
}

std::string logEntryText(const LogEntry& entry, bool withDoc) {
    TomlDocFunction doc;
    if (withDoc) {
        doc = [](const std::vector<std::string>& path, const std::vector<const Value*>&) -> std::optional<std::string> {
            if (path.size() == 1) {
                return std::string(kEntryDoc);
            }
            for (const auto& [key, text] : kLogDocs) {
                if (path.size() == 2 && path[1] == key) {
                    return std::string(text);
                }
            }
            return std::nullopt;
        };
    }
    TomlWriter writer(doc);
    writer.table({"log"}, entry.plain(), {}, true, withDoc);
    std::string out;
    for (const std::string& line : writer.lines()) {
        out += line;
        out += '\n';
    }
    // "\n".join(lines).lstrip("\n") + "\n"
    out.erase(0, out.find_first_not_of('\n'));
    return out.empty() ? "\n" : out;
}

fs::path splitLog(const fs::path& project) {
    detail::ProjectLock lock(project);
    return splitLogUnlocked(project);
}

fs::path appendLog(const fs::path& project, const LogEntry& entry) {
    std::error_code ec;
    if (!fs::is_regular_file(project, ec)) {
        throw ProjectError(detail::displayPath(project));
    }
    detail::ProjectLock lock(project);
    fs::path target = logPath(project);
    if (target == project && fs::file_size(project) > logSplitBytes()) {
        target = splitLogUnlocked(project);
    }
    const bool hasMarker = containsMarker(target);
    std::string text;
    if (!hasMarker) {
        text += "\n" + std::string(kLogMarker) + "\n";
    }
    text += "\n" + logEntryText(entry, !hasMarker);
    detail::writeTextFile(target, text, true);
    return target;
}

std::string exceptionSummary(const std::exception& error) {
    if (const auto* named = dynamic_cast<const NamedError*>(&error)) {
        return named->pythonName() + ": " + error.what();
    }
    if (dynamic_cast<const Cancelled*>(&error)) {
        return std::string("Cancelled: ") + error.what();
    }
    if (const auto* fsError = dynamic_cast<const fs::filesystem_error*>(&error)) {
        // Python's OSError text: [Errno 2] No such file or directory: 'a.pdb'
        const std::error_code code = fsError->code();
        std::string name = "OSError";
        if (code == std::errc::no_such_file_or_directory) {
            name = "FileNotFoundError";
        } else if (code == std::errc::permission_denied) {
            name = "PermissionError";
        } else if (code == std::errc::file_exists) {
            name = "FileExistsError";
        } else if (code == std::errc::is_a_directory) {
            name = "IsADirectoryError";
        } else if (code == std::errc::not_a_directory) {
            name = "NotADirectoryError";
        }
        std::string text = name + ": [Errno " + std::to_string(code.value()) + "] " + code.message();
        if (!fsError->path1().empty()) {
            text += ": " + pythonRepr(Value(pathToUtf8(fsError->path1())));
            if (!fsError->path2().empty()) {
                text += " -> " + pythonRepr(Value(pathToUtf8(fsError->path2())));
            }
        }
        return text;
    }
    if (dynamic_cast<const ProjectError*>(&error) || dynamic_cast<const std::invalid_argument*>(&error)) {
        return std::string("ValueError: ") + error.what();
    }
    if (dynamic_cast<const std::out_of_range*>(&error)) {
        return "KeyError: " + pythonRepr(Value(std::string(error.what())));  // str(KeyError) is its repr
    }
    if (dynamic_cast<const std::bad_alloc*>(&error)) {
        return "MemoryError: ";
    }
    if (dynamic_cast<const std::runtime_error*>(&error)) {
        return std::string("RuntimeError: ") + error.what();
    }
    return std::string("Exception: ") + error.what();
}

void logged(const fs::path& project, const std::string& step, const std::string& title,
            const std::vector<std::string>& citations, const std::function<void(LogEntry&)>& body) {
    LogEntry entry(step, title, citations);
    const double t0 = clockMonotonic();
    auto finish = [&] {
        entry.duration_s = clockMonotonic() - t0;
        std::error_code ec;
        if (project.empty() || !fs::is_regular_file(project, ec) || isParamsLog(project)) {
            return;
        }
        try {
            appendLog(project, entry);
        } catch (const std::exception& error) {  // a failed log write never fails the step
            logWarning("could not write the step log of " + detail::displayPath(project) + ": " + error.what());
        }
    };
    try {
        body(entry);
    } catch (const Cancelled&) {
        entry.status = "cancelled";
        finish();
        throw;
    } catch (const std::exception& error) {
        entry.status = "failed";
        if (entry.summary.empty()) {
            entry.summary = exceptionSummary(error);
        }
        finish();
        throw;
    } catch (...) {
        entry.status = "failed";
        if (entry.summary.empty()) {
            entry.summary = "Exception: unknown error";
        }
        finish();
        throw;
    }
    finish();
}

}  // namespace mdbridge::core
