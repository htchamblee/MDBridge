#pragma once
// Error types whose MDATK counterparts are Python exceptions. Callers of MDATK catch ValueError
// around most operations; the C++ code catches mdbridge::core::ValueError at the same places.

#include <stdexcept>
#include <string>
#include <utility>

namespace mdbridge::core {

// Python's ValueError: bad input (a project entry, a selection, a residue range ...).
class ValueError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
    // The Python class MDATK raises (step-log summaries name it).
    virtual const char* pythonName() const noexcept { return "ValueError"; }
};

// MDAnalysis' NoDataError: the topology lacks an attribute (bonds, elements, ...).
class NoDataError : public ValueError {
public:
    using ValueError::ValueError;
    const char* pythonName() const noexcept override { return "NoDataError"; }
};

// The user stopped the work (mdatk.progress.Cancelled); a logged step ends "cancelled".
class Cancelled : public std::runtime_error {
public:
    Cancelled() : std::runtime_error("operation cancelled") {}
    using std::runtime_error::runtime_error;
};

// An error that names the Python exception class MDATK reports in the same situation; the summary
// of a failed step is "<pythonName>: <message>", e.g. "RuntimeError: pmemd crashed".
class NamedError : public std::runtime_error {
public:
    NamedError(std::string pythonName, const std::string& message)
        : std::runtime_error(message), m_name(std::move(pythonName)) {}
    const std::string& pythonName() const { return m_name; }

private:
    std::string m_name;
};

}  // namespace mdbridge::core
