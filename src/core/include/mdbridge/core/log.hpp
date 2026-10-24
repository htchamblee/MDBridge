#pragma once
// Messages that MDATK sends through Python's logging (warnings about unknown keys, unreadable step
// logs, ...). The CLI and GUI print them as "WARNING: <message>", as MDATK's front ends do.

#include <functional>
#include <string_view>

namespace mdbridge::core {

enum class LogLevel { Debug, Info, Warning, Error };

// "DEBUG", "INFO", "WARNING", "ERROR" (Python's level names).
std::string_view logLevelName(LogLevel level);

using LogSink = std::function<void(LogLevel level, std::string_view message)>;

// Installs the process-wide sink and returns the previous one. The default sink writes warnings
// and errors to stderr as "<LEVEL>: <message>" and drops the rest. Thread-safe.
LogSink setLogSink(LogSink sink);

void logMessage(LogLevel level, std::string_view message);
inline void logWarning(std::string_view message) { logMessage(LogLevel::Warning, message); }
inline void logInfo(std::string_view message) { logMessage(LogLevel::Info, message); }

}  // namespace mdbridge::core
