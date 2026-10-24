#pragma once
// A project file in a few lines, for `mdbridge project show` and the interface's sidebar.

#include <string>
#include <vector>

#include "mdbridge/core/project_file.hpp"

namespace mdbridge::app {

struct OpenedProject {
    core::ProjectFile file;
    std::filesystem::path logPath;
    core::Array log;         // step log entries, oldest first
    std::string logWarning;  // why the log could not be read ("" if it could)
};

// Reads the project and its step log (throws core::ProjectError like core::readProject).
OpenedProject openProject(const std::filesystem::path& path);

struct SummaryRow {
    std::string label;
    std::string value;
};

std::vector<SummaryRow> summarize(const OpenedProject& project);

// One line for a scalar or a short array (e.g. a log entry's field in a table cell).
std::string valueText(const core::Value& value);

// Everything MDBridge read: header, settings, log (JSON, see core::toJson).
core::Value projectData(const OpenedProject& project);

}  // namespace mdbridge::app
