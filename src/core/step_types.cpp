// TEMPORARY STUB (replaced by the step-type registry)
// Knows only "input-file" and its path setting, so that the parameter model links and resolves
// planner paths; the rest exists only for the core tests to link.

#include "mdbridge/core/step_types.hpp"

#include <stdexcept>
#include <string>

namespace mdbridge::core {

namespace {

const StepType& inputFile() {
    static const StepType type = [] {
        StepType t;
        t.name = "input-file";
        t.generates = false;
        t.pathSettings = {"path"};
        return t;
    }();
    return type;
}

}  // namespace

const StepType* findStepType(std::string_view name) {
    return name == "input-file" ? &inputFile() : nullptr;
}

const StepType& stepType(std::string_view name) {
    if (const StepType* type = findStepType(name)) {
        return *type;
    }
    throw std::out_of_range("unknown step type '" + std::string(name) + "' (known: input-file)");
}

Table StepType::normalise(const Table& settings, std::vector<Issue>*,
                          const std::function<void(const std::string&)>&) const {
    return settings;
}

std::optional<std::string> settingComment(std::string_view, std::string_view) {
    return std::nullopt;
}

}  // namespace mdbridge::core
