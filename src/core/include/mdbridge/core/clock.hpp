#pragma once
// The clocks of project files and step logs, replaceable so that tests can write files with
// fixed time stamps and durations.

#include <functional>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

struct Clock {
    // Local wall-clock time in whole seconds (MDATK's projectfile._now()).
    std::function<DateTime()> now;
    // Seconds of a monotonic clock (time.monotonic()), for durations.
    std::function<double()> monotonic;
};

// Installs `clock` (empty members: the system clocks) and returns the previous one.
Clock setClock(Clock clock);

DateTime clockNow();
double clockMonotonic();

}  // namespace mdbridge::core
