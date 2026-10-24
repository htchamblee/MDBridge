#include "project_summary.hpp"

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"

namespace mdbridge::app {

using core::Array;
using core::Table;
using core::Value;

namespace {

std::string text(const Value* value, const std::string& fallback = "(not set)") {
    if (!value || value->isNull()) {
        return fallback;
    }
    const std::string out = valueText(*value);
    return out.empty() ? fallback : out;
}

std::string plural(std::size_t n, const char* one, const char* many) {
    return std::to_string(n) + " " + (n == 1 ? one : many);
}

std::string notesText(const Value* notes) {
    const std::string all = notes ? notes->stringOr() : "";
    if (all.empty()) {
        return "(none)";
    }
    const std::size_t end = all.find('\n');
    if (end == std::string::npos) {
        return all;
    }
    std::size_t lines = 1;
    for (std::size_t at = end; at != std::string::npos; at = all.find('\n', at + 1)) {
        ++lines;
    }
    return all.substr(0, end) + " (" + std::to_string(lines) + " lines)";
}

std::string groupsText(const Value* groups) {
    const Array* list = groups ? groups->get<Array>() : nullptr;
    if (!list) {
        return "(MDATK's defaults)";
    }
    if (list->empty()) {
        return "(none)";
    }
    std::string out;
    for (const Value& group : *list) {
        const std::string residues = text(group.find("residues"), "");
        const std::string selection = text(group.find("selection"), "");
        out += (out.empty() ? "" : ", ") + text(group.find("name"), "Group") + " (" +
               (selection.empty() ? residues : selection) + ")";
    }
    return out;
}

std::string plannerText(const Value* workflow) {
    const std::size_t steps = workflow ? (workflow->find("nodes") ? workflow->find("nodes")->size() : 0) : 0;
    const std::size_t edges = workflow ? (workflow->find("edges") ? workflow->find("edges")->size() : 0) : 0;
    if (steps == 0) {
        return "(empty)";
    }
    return plural(steps, "step", "steps") + ", " + plural(edges, "connection", "connections");
}

std::string layoutText(const Value* layout) {
    const Array* panels = layout && layout->find("panels") ? layout->find("panels")->get<Array>() : nullptr;
    if (!panels) {
        return "(MDATK's defaults)";
    }
    std::string types;
    for (const Value& panel : *panels) {
        types += (types.empty() ? "" : ", ") + text(panel.find("type"), "structure");
    }
    return plural(panels->size(), "panel", "panels") + (types.empty() ? "" : ": " + types);
}

std::string logText(const OpenedProject& project) {
    if (!project.logWarning.empty()) {
        return project.logWarning;
    }
    if (project.log.empty()) {
        return "(empty)";
    }
    const Value& last = project.log.back();
    std::string out = plural(project.log.size(), "entry", "entries") + "; last: " + text(last.find("time"), "") +
                      " " + text(last.find("step"), "") + " " + text(last.find("title"), "") + " (" +
                      text(last.find("status"), "ok") + ")";
    if (project.logPath != project.file.path) {
        out += "; in " + core::pathToUtf8(project.logPath.filename());
    }
    return out;
}

}  // namespace

OpenedProject openProject(const std::filesystem::path& path) {
    OpenedProject project;
    project.file = core::readProject(path);
    project.logPath = core::logPath(path);
    if (project.file.kind == core::FileKind::Project) {  // a run's parameter log has no step log
        project.log = core::readLog(path, &project.logWarning);
    }
    return project;
}

std::string valueText(const Value& value) {
    switch (value.type()) {
    case Value::Type::Null:
        return "";
    case Value::Type::Bool:
        return *value.get<bool>() ? "true" : "false";
    case Value::Type::Int:
        return std::to_string(*value.get<std::int64_t>());
    case Value::Type::Float:
        return core::pythonFloatRepr(*value.get<double>());
    case Value::Type::String:
        return *value.get<std::string>();
    case Value::Type::DateTime:
        return value.get<core::DateTime>()->isoformat();
    case Value::Type::Array:
    case Value::Type::Table:
        break;
    }
    return core::toJson(value, -1);
}

std::vector<SummaryRow> summarize(const OpenedProject& project) {
    const core::ProjectFile& file = project.file;
    const Value& s = file.settings;
    const bool paramsLog = file.kind == core::FileKind::ParamsLog;
    std::vector<SummaryRow> rows;
    rows.push_back({"File", core::pathToUtf8(file.path)});
    if (paramsLog) {
        rows.push_back({"Type", "MDATK parameter log of one run (opens as a project, to repeat the run)"});
        rows.push_back({"Written by", "MDATK " + text(file.header.find("mdatk_version"), "?") + ", run finished " +
                                          text(file.header.find("created"), "?")});
    } else {
        const Value* version = file.header.find("format_version");
        rows.push_back({"Type", "MDATK project" + (version ? " (file layout " + text(version) + ")" : std::string())});
        rows.push_back({"Written by", file.header.find("mdatk_version")
                                          ? "MDATK " + text(file.header.find("mdatk_version")) + ", saved " +
                                                text(file.header.find("saved"), "?")
                                          : "(no version recorded: written by hand?)"});
    }
    rows.push_back({"Notes", notesText(s.find("notes"))});
    const Value* inputs = s.find("inputs") ? s.find("inputs") : s.find("input");  // MDATK accepts both
    rows.push_back({"Topology", text(inputs ? inputs->find("topology") : nullptr)});
    const Value* trajectories = inputs ? inputs->find("trajectories") : nullptr;
    if (!trajectories && inputs) {
        trajectories = inputs->find("trajectory");
    }
    std::string trajectoryText = "(not set)";
    if (const Array* list = trajectories ? trajectories->get<Array>() : nullptr; list && !list->empty()) {
        trajectoryText = plural(list->size(), "file", "files") + ": " + text(&list->front()) +
                         (list->size() > 1 ? ", ..." : "");
    } else if (trajectories && trajectories->isString() && !trajectories->stringOr().empty()) {
        trajectoryText = "1 file: " + trajectories->stringOr();
    }
    rows.push_back({"Trajectories", trajectoryText});
    rows.push_back({"Residue groups", groupsText(s.find("groups"))});
    rows.push_back({"Figure layout", layoutText(s.find("layout"))});
    rows.push_back({"Project Planner", plannerText(s.find("workflow"))});
    if (!paramsLog) {
        rows.push_back({"Step log", logText(project)});
    }
    return rows;
}

Value projectData(const OpenedProject& project) {
    Table data;
    data.emplace_back("path", core::pathToUtf8(project.file.path));
    data.emplace_back("kind", project.file.kind == core::FileKind::ParamsLog ? "params-log" : "project");
    data.emplace_back("header", project.file.header);
    data.emplace_back("settings", project.file.settings);
    data.emplace_back("log_path", core::pathToUtf8(project.logPath));
    data.emplace_back("log", project.log);
    if (!project.logWarning.empty()) {
        data.emplace_back("log_warning", project.logWarning);
    }
    return data;
}

}  // namespace mdbridge::app
