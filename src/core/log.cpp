#include "mdbridge/core/log.hpp"

#include <cstdio>
#include <mutex>
#include <string>

namespace mdbridge::core {

namespace {

// logging.basicConfig(level=WARNING, format="%(levelname)s: %(message)s") of MDATK's front ends
void defaultSink(LogLevel level, std::string_view message) {
    if (level != LogLevel::Warning && level != LogLevel::Error) {
        return;
    }
    std::string line(logLevelName(level));
    line += ": ";
    line += message;
    line += '\n';
    std::fwrite(line.data(), 1, line.size(), stderr);
    std::fflush(stderr);
}

struct SinkState {
    // recursive: a sink may log, or install another sink, while it runs
    std::recursive_mutex mutex;
    LogSink sink = defaultSink;
};

SinkState& sinkState() {
    static SinkState state;
    return state;
}

}  // namespace

std::string_view logLevelName(LogLevel level) {
    switch (level) {
    case LogLevel::Debug:
        return "DEBUG";
    case LogLevel::Info:
        return "INFO";
    case LogLevel::Warning:
        return "WARNING";
    case LogLevel::Error:
        return "ERROR";
    }
    return "UNKNOWN";
}

LogSink setLogSink(LogSink sink) {
    SinkState& state = sinkState();
    const std::lock_guard lock(state.mutex);
    LogSink previous = std::move(state.sink);
    state.sink = sink ? std::move(sink) : LogSink(defaultSink);
    return previous;
}

void logMessage(LogLevel level, std::string_view message) {
    SinkState& state = sinkState();
    // messages reach the sink one at a time, as Python's handlers take their lock in emit()
    const std::lock_guard lock(state.mutex);
    const LogSink sink = state.sink;
    sink(level, message);
}

}  // namespace mdbridge::core
