#pragma once
// Step types of the Project Planner and script wizard: their settings schemas, as MDATK declares
// them in mdatk/workflows (registry.py, schema.py, base.py).
//
// The schemas are data of the project file: every entry of [workflow.nodes.settings] is commented
// with the option text of its step type, and settings that hold file paths are resolved against
// the project's folder. They therefore live in core, generated from MDATK v0.0.4.1
// (tools/gen_step_types.py -> step_types.inc); the workflow module adds checks and script
// generators for the same step types and may register further ones.

#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// KINDS of schema.py, in that order.
enum class OptionKind { Bool, Choice, Radio, Int, Float, Text, Multiline, File };
std::string_view optionKindName(OptionKind kind);

// A condition on other settings (schema.py requirements_met), checked in order:
// a value of true: the setting must be truthy; false: it must not be; an array: the setting must
// equal one of its items (Python ==, so 1 == true == 1.0); anything else: it must equal the value.
using Requirements = std::vector<std::pair<std::string, Value>>;
bool requirementsMet(const Requirements& requirements, const Table& values);

struct Choice {
    std::string value;
    std::string label;
};

struct Option {
    std::string key;
    std::string label;
    OptionKind kind = OptionKind::Text;
    Value defaultValue;  // with MDATK's Python type: Null, Bool, Int, Float or String
    std::string tooltip;
    std::vector<Choice> choices;  // choice and radio options only
    std::optional<double> minimum, maximum, step;
    int digits = 0;
    bool optional = false;  // int / float: an empty value means "not set"
    std::string placeholder;
    std::string unit;
    Requirements requirements;  // "requires" in MDATK

    // Option.coerce: the value converted to the option's type. Throws std::invalid_argument with
    // MDATK's message, e.g. "engine: 'x' is not one of pmemd.cuda, pmemd, pmemd.MPI".
    Value coerce(const Value& value) const;
};

struct RecordTemplate {
    std::string id;
    std::string label;
    Table values;  // overrides of the record defaults
};

struct Section {
    std::string key;
    std::string title;
    std::string tooltip;
    bool repeat = false;  // a list of records stored under `key`
    std::string itemLabel = "name";
    std::vector<Option> options;
    std::vector<RecordTemplate> templates;
    std::vector<std::string> defaultItems;  // template ids of the default records
    Requirements requirements;  // "requires" in MDATK

    const Option* option(std::string_view key) const;
    Table recordDefaults() const;
    // Throws std::out_of_range("<key>: unknown template '<id>'").
    Table recordTemplate(std::string_view id) const;
    Array defaultRecords() const;
};

struct Port {
    std::string name;
    std::string label;
    std::vector<std::string> kinds;  // file kinds (FILE_KINDS ids); the first is the port's kind
    bool optional = false;
    bool multiple = false;
    std::string tooltip;
};

// A warning or error about settings (schema.Issue). str(issue) = "WARNING: <message>".
struct Issue {
    std::string level;  // error, warning, info
    std::string message;
    std::vector<std::string> keys;
    std::string text() const;
};

struct StepType {
    std::string name;
    std::string label;
    std::string category;
    std::string tag;  // step folder name and default file tag
    std::string description;
    bool generates = true;  // false for planner helpers (note, input-file, mdatk)
    std::vector<std::string> pathSettings;  // settings holding paths relative to the project file
    std::vector<Port> inputs;
    std::vector<Port> outputs;
    std::vector<Section> sections;

    const Section* section(std::string_view key) const;
    // The option of a non-repeated section with this key (Tool.flat_options).
    const Option* flatOption(std::string_view key) const;
    // Tool.defaults(): options of plain sections, default records of repeated ones, in order.
    Table defaults() const;
    // Tool.normalise: defaults overlaid with the coerced given settings. A value that cannot be
    // used becomes "<message> - using the default": an Issue in `issues`, or, without a list, a
    // warning "<key>: <message> - using the default" to `warn`. Unknown keys always go to `warn`
    // as "<name>: ignoring unknown settings a, b" (sorted).
    Table normalise(const Table& settings, std::vector<Issue>* issues = nullptr,
                    const std::function<void(const std::string&)>& warn = {}) const;
};

// All step types, in registration order (MDATK's TOOLS).
const std::vector<StepType>& stepTypes();

// The step type `name`, or nullptr.
const StepType* findStepType(std::string_view name);

// get_tool(name): throws std::out_of_range("unknown step type 'x' (known: a, b, ...)").
const StepType& stepType(std::string_view name);

// Adds a step type (registry.register). Throws std::invalid_argument for an empty or taken name
// ("bad or duplicate tool name 'x'") and for option keys used twice ("x: duplicate option keys
// ['a', 'b']").
void registerStepType(StepType type);

// CATEGORIES of registry.py, and the step types grouped by them (tools_by_category; step types of
// other categories last, under "Other").
const std::vector<std::string>& stepCategories();
std::vector<std::pair<std::string, std::vector<const StepType*>>> stepTypesByCategory();

// The comment MDATK writes above setting `key` of a step of type `tool` (projectfile._tool_option
// and _option_text): the option's "label: tooltip[ Choices: a, b.]", or a repeated section's
// tooltip. std::nullopt for an unknown step type or setting.
std::optional<std::string> settingComment(std::string_view tool, std::string_view key);

// FILE_KINDS of schema.py, in order.
struct FileKindInfo {
    std::string id;
    std::string label;
    std::vector<std::string> extensions;
};
const std::vector<FileKindInfo>& fileKinds();
std::string kindLabel(std::string_view kind);    // the label, or `kind` itself when unknown
std::string guessKind(std::string_view path);    // by extension (case-insensitive), else "text"

}  // namespace mdbridge::core
