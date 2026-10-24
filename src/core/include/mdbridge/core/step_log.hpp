#pragma once
// The step log of a project: one [[log]] entry per step run (analyses, renders, generated scripts
// that ran, PDB clean-ups, density fits), appended after the settings of the project file or, once
// that would grow past 50 MB, to <name>.log.toml next to it. MDATK's projectfile.py (LogEntry,
// append_log, _split_log, logged).

#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// Program output kept in an entry (the last characters, as code points).
inline constexpr std::size_t kLogOutputChars = 20000;
inline constexpr std::uintmax_t kDefaultLogSplitBytes = 50ull * 1024 * 1024;

// The size above which the next append moves the log out of the project file. Settable for tests.
std::uintmax_t logSplitBytes();
void setLogSplitBytes(std::uintmax_t bytes);

struct LogEntry {
    std::string step;  // "render", "workflow write", "pdb clean", or a planner step number ("03")
    std::string title;
    std::string status = "ok";  // ok, failed or cancelled
    DateTime time;  // when the step started
    double duration_s = 0.0;
    std::string command;
    std::string folder;
    std::vector<std::string> files;
    std::string summary;
    std::vector<std::string> citations;
    std::string output;

    // An entry that starts now (clockNow()).
    explicit LogEntry(std::string step = {}, std::string title = {},
                      std::vector<std::string> citations = {});

    // The entry as written: duration rounded to 2 decimals (Python round()), output cut to its last
    // kLogOutputChars code points behind "[...]\n", empty strings and lists left out.
    Value plain() const;
};

// The text of one entry, with comments (the first entry of a log) or without.
std::string logEntryText(const LogEntry& entry, bool withDoc);

// append_log: appends the entry to the project's log (moving the log out first when the project
// file is too large) and returns the file written to. Throws ProjectError when `project` is not a
// file, std::filesystem::filesystem_error when the log cannot be written.
std::filesystem::path appendLog(const std::filesystem::path& project, const LogEntry& entry);

// _split_log: moves the log of `project` to <stem>.log.toml, sets log_file, returns the side file.
std::filesystem::path splitLog(const std::filesystem::path& project);

// "<class>: <message>" of an exception, the summary of a failed step: NamedError's name; the
// Python name of a ValueError (ValueError for ProjectError, SelectionError, NoDataError ...);
// ValueError for std::invalid_argument; KeyError (message as its repr) for std::out_of_range;
// FileNotFoundError and the other OSError names, with Python's "[Errno N] text: 'path'", for
// filesystem errors; RuntimeError for other std::runtime_error; Exception otherwise.
std::string exceptionSummary(const std::exception& error);

// logged(path, step, title, citations): runs `body` with an entry that started now. When `body`
// throws, the entry's status becomes "cancelled" (Cancelled) or "failed" (with the summary
// exceptionSummary(), unless the body set one) and the exception propagates after the entry is
// written. The duration is always measured. Nothing is written when `project` is empty, not a
// file, or a parameter log; a failed write is only logged as a warning ("could not write the step
// log of <path>: <reason>").
void logged(const std::filesystem::path& project, const std::string& step, const std::string& title,
            const std::vector<std::string>& citations, const std::function<void(LogEntry&)>& body);

}  // namespace mdbridge::core
