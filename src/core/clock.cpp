#include "mdbridge/core/clock.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <mutex>

namespace mdbridge::core {

namespace {

// datetime.now().replace(microsecond=0): local time, naive, whole seconds
DateTime systemNow() {
    const auto now = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());
    const std::time_t seconds = std::chrono::system_clock::to_time_t(now);
    std::tm local{};
#if defined(_WIN32)
    localtime_s(&local, &seconds);
#else
    localtime_r(&seconds, &local);
#endif
    DateTime out;
    out.kind = DateTime::Kind::LocalDateTime;
    out.year = local.tm_year + 1900;
    out.month = local.tm_mon + 1;
    out.day = local.tm_mday;
    out.hour = local.tm_hour;
    out.minute = local.tm_min;
    out.second = std::min(local.tm_sec, 59);  // a leap second, as CPython clamps it
    return out;
}

// time.monotonic()
double systemMonotonic() {
    const auto since = std::chrono::steady_clock::now().time_since_epoch();
    return std::chrono::duration<double>(since).count();
}

struct ClockState {
    std::mutex mutex;
    Clock clock;
};

ClockState& clockState() {
    static ClockState state;
    return state;
}

}  // namespace

Clock setClock(Clock clock) {
    ClockState& state = clockState();
    const std::lock_guard lock(state.mutex);
    std::swap(state.clock, clock);
    return clock;
}

DateTime clockNow() {
    std::function<DateTime()> now;
    {
        ClockState& state = clockState();
        const std::lock_guard lock(state.mutex);
        now = state.clock.now;
    }
    return now ? now() : systemNow();
}

double clockMonotonic() {
    std::function<double()> monotonic;
    {
        ClockState& state = clockState();
        const std::lock_guard lock(state.mutex);
        monotonic = state.clock.monotonic;
    }
    return monotonic ? monotonic() : systemMonotonic();
}

}  // namespace mdbridge::core
