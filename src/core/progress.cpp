#include "mdbridge/core/progress.hpp"

#include <algorithm>
#include <utility>

#include "mdbridge/core/log.hpp"

namespace mdbridge::core {

Reporter::Reporter(LogCallback log, ProgressCallback progress, std::shared_ptr<std::atomic<bool>> cancelFlag)
    : m_log(std::move(log)), m_progress(std::move(progress)), m_cancel(std::move(cancelFlag)) {
    if (!m_cancel) {
        m_cancel = std::make_shared<std::atomic<bool>>(false);
    }
}

void Reporter::log(const std::string& message) const {
    logMessage(LogLevel::Info, message);
    if (m_log) {
        m_log(message);
    }
}

void Reporter::rawLog(const std::string& message) const {
    logMessage(LogLevel::Debug, message);
    if (m_log) {
        m_log(message);
    }
}

void Reporter::progress(double fraction, const std::string& message) const {
    // min(1.0, max(0.0, fraction)) as Python evaluates it: NaN becomes 0.0
    fraction = std::min(1.0, std::max(0.0, fraction));
    if (m_progress) {
        m_progress(m_start + (m_end - m_start) * fraction, message);
    }
}

Reporter Reporter::sub(double start, double end) const {
    Reporter out(m_log, m_progress, m_cancel);
    out.m_start = m_start + (m_end - m_start) * start;
    out.m_end = m_start + (m_end - m_start) * end;
    return out;
}

void Reporter::checkCancel() const {
    if (m_cancel->load()) {
        throw Cancelled("operation cancelled");
    }
}

void Reporter::cancel() const { m_cancel->store(true); }

bool Reporter::cancelled() const { return m_cancel->load(); }

}  // namespace mdbridge::core
