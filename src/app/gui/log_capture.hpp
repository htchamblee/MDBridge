#pragma once
// The messages (core::logMessage) of one thread kept apart while a ThreadLogCapture lives on it: the
// window's log sink hands them to the capture instead of the log. The figure preview draws with
// one, so the window can log a drawing's warnings only when they change (as MDATK's canvas logs
// what a redraw after an edit warns about, not every redraw).

#include <string>
#include <vector>

namespace mdbridge::app {

class ThreadLogCapture {
public:
    explicit ThreadLogCapture(std::vector<std::string>* lines) : m_previous(slot()) { slot() = lines; }
    ~ThreadLogCapture() { slot() = m_previous; }
    ThreadLogCapture(const ThreadLogCapture&) = delete;
    ThreadLogCapture& operator=(const ThreadLogCapture&) = delete;

    // The lines of the capture of this thread (null without one).
    static std::vector<std::string>* current() { return slot(); }

private:
    static std::vector<std::string>*& slot() {
        thread_local std::vector<std::string>* lines = nullptr;
        return lines;
    }
    std::vector<std::string>* m_previous;
};

}  // namespace mdbridge::app
