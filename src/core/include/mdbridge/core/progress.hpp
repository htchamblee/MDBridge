#pragma once
// Progress and cancellation shared by the CLI and the GUI (mdatk/progress.py Reporter): messages,
// a fraction of the work done, sub-ranges for parts of the work, and one cancel flag for the whole
// family of reporters.

#include <atomic>
#include <functional>
#include <memory>
#include <string>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::core {

class Reporter {
public:
    using LogCallback = std::function<void(const std::string& message)>;
    using ProgressCallback = std::function<void(double fraction, const std::string& message)>;

    explicit Reporter(LogCallback log = {}, ProgressCallback progress = {},
                      std::shared_ptr<std::atomic<bool>> cancelFlag = {});

    // A message: logged at info level and passed to the log callback.
    void log(const std::string& message) const;
    // Forwarded program output: logged at debug level and passed to the log callback.
    void rawLog(const std::string& message) const;
    // `fraction` of this reporter's range, clamped to [0, 1] (NaN as 0), mapped onto the whole.
    void progress(double fraction, const std::string& message = {}) const;
    // A reporter for the part [start, end] of this one's range (not clamped), same callbacks and
    // cancel flag.
    Reporter sub(double start, double end) const;

    // Throws Cancelled("operation cancelled") once any reporter of the family was cancelled.
    void checkCancel() const;
    void cancel() const;
    bool cancelled() const;

private:
    LogCallback m_log;
    ProgressCallback m_progress;
    std::shared_ptr<std::atomic<bool>> m_cancel;
    double m_start = 0.0;
    double m_end = 1.0;
};

}  // namespace mdbridge::core
