// --set KEY=VALUE: MDATK's apply_override (mdatk/cli.py).
#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "commands.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::app {

namespace {

using core::Value;

// Attributes of MDATK's Project object besides its fields (methods and properties): hasattr()
// finds them, so `--set out_dir=x` passes the check and reaches from_dict as an unknown key.
constexpr std::array<std::string_view, 14> kProjectAttributes{
    "cache_dir", "copy",          "density_active", "density_problems", "from_dict", "group", "group_names",
    "load",      "out_dir",       "resolve_paths",  "save",             "to_dict",   "validate", "workflow_dir"};
// The str-valued ones (properties): a non-string JSON value is replaced by the text, as for fields.
constexpr std::array<std::string_view, 3> kTextAttributes{"cache_dir", "out_dir", "workflow_dir"};

// The free-form tables of the model (dicts, not dataclasses): no parameters inside them.
constexpr std::array<std::string_view, 3> kFreeFormFields{"options", "settings", "files"};

template <std::size_t N>
bool contains(const std::array<std::string_view, N>& names, std::string_view name) {
    for (std::string_view n : names) {
        if (n == name) {
            return true;
        }
    }
    return false;
}

std::string repr(std::string_view text) { return core::pythonStringRepr(std::string(text)); }

// _parse_value: JSON when the text is JSON, else the text.
Value parseValue(std::string_view raw) {
    try {
        return core::parseJson(raw, "--set");
    } catch (const core::ProjectError&) {
        return Value(std::string(raw));
    }
}

// part.lstrip("-").isdigit()
bool looksLikeIndex(std::string_view part) {
    std::size_t start = 0;
    while (start < part.size() && part[start] == '-') {
        ++start;
    }
    const std::u32string digits = core::decodeUtf8(part.substr(start));
    if (digits.empty()) {
        return false;
    }
    for (char32_t c : digits) {
        if (!core::isPythonDigit(c)) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> splitPath(std::string_view path) {
    std::vector<std::string> parts;
    const std::string stripped = core::pythonStrip(path);
    std::size_t start = 0;
    while (start <= stripped.size()) {
        std::size_t end = stripped.find('.', start);
        if (end == std::string::npos) {
            end = stripped.size();
        }
        if (end > start) {
            parts.push_back(stripped.substr(start, end - start));
        }
        start = end + 1;
    }
    return parts;
}

}  // namespace

void applyOverride(core::Project& project, std::string_view assignment) {
    const std::size_t eq = assignment.find('=');
    if (eq == std::string_view::npos) {
        throw core::ValueError("--set expects key=value, got " + repr(assignment));
    }
    const std::string path(assignment.substr(0, eq));
    const std::string raw(assignment.substr(eq + 1));
    const std::vector<std::string> parts = splitPath(path);
    if (parts.empty()) {
        throw core::ValueError("empty key in " + repr(assignment));
    }

    // Walk the plain form of the project as MDATK walks the objects: tables of the model
    // (dataclasses) have attributes, lists have elements, anything else has neither.
    Value data = core::toPlain(project);
    Value* parent = &data;
    bool parentIsModel = true;
    bool parentIsRoot = true;
    Value* current = nullptr;
    bool currentIsModel = false;
    bool currentIsText = false;  // a str attribute that is no field
    bool currentIsAttribute = false;  // a method or property of Project: not in the plain data
    std::size_t index = 0;  // the list element of `current` when `parent` is a list
    std::string modelPath;  // the model path of `parent` (keys, list indices left out)
    for (std::size_t n = 0; n < parts.size(); ++n) {
        const std::string& part = parts[n];
        if (n > 0 && currentIsAttribute) {
            throw core::ValueError("unknown parameter " + repr(path));
        }
        if (core::Array* list = parent->get<core::Array>()) {
            const auto size = static_cast<std::int64_t>(list->size());
            bool inRange = false;
            std::int64_t k = 0;
            if (looksLikeIndex(part)) {
                try {
                    k = core::pythonInt(part);
                    inRange = -size <= k && k < size;
                } catch (const std::overflow_error&) {
                    inRange = false;  // beyond 64 bits: outside any list
                } catch (const core::ValueError&) {
                    inRange = false;  // "--1": MDATK reports int()'s "invalid literal" (PRJ-04)
                }
            }
            if (!inRange) {
                throw core::ValueError(path + ": no list element " + repr(part));
            }
            index = static_cast<std::size_t>(k < 0 ? k + size : k);
            current = &(*list)[index];
            currentIsModel = current->isTable();  // the model's lists of tables hold dataclasses
            currentIsText = false;
            currentIsAttribute = false;
        } else {
            // a field of the model, not a key the project keeps without knowing it (D9): hasattr()
            const std::string fieldPath = modelPath.empty() ? part : modelPath + "." + part;
            Value* field = parentIsModel && core::modelComment(fieldPath) ? parent->find(part) : nullptr;
            if (field != nullptr) {
                modelPath = fieldPath;
                current = field;
                currentIsModel = field->isTable() && !contains(kFreeFormFields, part);
                currentIsText = false;
                currentIsAttribute = false;
            } else if (parentIsRoot && contains(kProjectAttributes, part)) {
                current = nullptr;
                currentIsModel = false;
                currentIsText = contains(kTextAttributes, part);
                currentIsAttribute = true;
            } else {
                throw core::ValueError("unknown parameter " + repr(path));
            }
        }
        if (n + 1 < parts.size()) {
            if (currentIsAttribute) {
                continue;  // the next part fails: a method or a str has no parameters
            }
            parent = current;
            parentIsModel = currentIsModel;
            parentIsRoot = false;
        }
    }

    Value value = parseValue(raw);
    const bool currentIsBool = current != nullptr && current->get<bool>() != nullptr;
    if (currentIsBool && value.get<bool>() == nullptr) {
        throw core::ValueError(path + " expects true/false");
    }
    const bool currentIsString = currentIsText || (current != nullptr && current->isString());
    if (currentIsString && !value.isString()) {
        value = Value(raw);
    }
    if (core::Array* list = parent->get<core::Array>()) {
        (*list)[index] = std::move(value);
    } else if (current != nullptr) {
        *current = std::move(value);
    } else {
        parent->get<core::Table>()->emplace_back(parts.back(), std::move(value));  // a new key, reported by from_dict
    }
    // the unknown keys of the project were reported when it was read: warn about new ones only
    const std::vector<std::string> before = core::unknownKeyPaths(project);
    const std::string version = project.mdatkVersion;  // not a setting: not in the plain data
    std::vector<std::string> ignored;
    project = core::projectFromPlain(data, {}, &ignored);
    project.mdatkVersion = version;
    std::vector<std::string> added;
    for (const std::string& key : core::unknownKeyPaths(project)) {
        if (std::find(before.begin(), before.end(), key) == before.end()) {
            added.push_back(key);
        }
    }
    if (!added.empty()) {
        core::logWarning(core::unknownKeysWarning(added));
    }
}

}  // namespace mdbridge::app
