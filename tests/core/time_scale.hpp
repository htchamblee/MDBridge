#pragma once
// The factor by which a test's time limits grow: MDBRIDGE_TIME_SCALE (a number >= 1; default 1).
// Builds that run slower by design set it, e.g. the sanitizer CI job (AddressSanitizer and
// UndefinedBehaviorSanitizer run the code several times slower), so a limit that checks that
// something is fast still checks it there, against a proportionally larger budget, and a limit
// that checks something happens at once (a cancel, a stall of the event loop) is not skipped.

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

namespace mdbridge::test {

inline double timeScale() {
    static const double scale = [] {
        const char* text = std::getenv("MDBRIDGE_TIME_SCALE");
        if (text == nullptr || *text == '\0') {
            return 1.0;
        }
        try {
            std::size_t used = 0;
            const double value = std::stod(text, &used);
            if (used == std::string(text).size() && value >= 1.0) {
                return value;
            }
        } catch (const std::exception&) {
        }
        // a mistyped scale must not quietly test against the wrong budget
        std::fprintf(stderr, "MDBRIDGE_TIME_SCALE=%s is not a number of at least 1\n", text);
        std::abort();
    }();
    return scale;
}

}  // namespace mdbridge::test
