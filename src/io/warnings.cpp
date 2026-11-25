#include "mdbridge/io/warnings.hpp"

#include <algorithm>
#include <utility>

#include "mdbridge/core/log.hpp"

namespace mdbridge::io {

namespace {
thread_local ReaderWarningScope* t_scope = nullptr;
}  // namespace

ReaderWarningScope::ReaderWarningScope() : m_outer(t_scope) {
    t_scope = this;
}

ReaderWarningScope::~ReaderWarningScope() {
    t_scope = m_outer;
}

std::vector<std::string> ReaderWarningScope::take() {
    return std::exchange(m_warnings, {});
}

void readerWarning(std::string message) {
    if (t_scope == nullptr) {
        core::logWarning(message);
        return;
    }
    std::vector<std::string>& list = t_scope->m_warnings;
    if (std::find(list.begin(), list.end(), message) == list.end()) {
        list.push_back(std::move(message));
    }
}

}  // namespace mdbridge::io
