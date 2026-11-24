#include "mdbridge/core/step_types.hpp"

#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

#include "mdbridge/core/log.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

#include "step_types.inc"

struct Registry {
    std::mutex mutex;
    std::vector<StepType> types = builtinStepTypes();
    std::vector<std::string> categories = builtinCategories();
    std::vector<FileKindInfo> kinds = builtinFileKinds();
};

Registry& registry() {
    static Registry instance;  // built on first use; thread-safe initialisation
    return instance;
}

std::string repr(std::string_view text) { return pythonRepr(Value(std::string(text))); }

std::string joined(const std::vector<std::string>& items, std::string_view separator = ", ") {
    std::string out;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i > 0) {
            out += separator;
        }
        out += items[i];
    }
    return out;
}

// dict lookup in a table built from pairs: a repeated key holds its last value
const Value* lookup(const Table& table, std::string_view key) {
    const Value* found = nullptr;
    for (const auto& [k, v] : table) {
        if (k == key) {
            found = &v;
        }
    }
    return found;
}

// d[key] = value: a new key goes last, an existing one keeps its place
void assign(Table& table, const std::string& key, Value value) {
    for (auto& [k, v] : table) {
        if (k == key) {
            v = std::move(value);
            return;
        }
    }
    table.emplace_back(key, std::move(value));
}

// ---- Python semantics ------------------------------------------------------------------------

bool truthy(const Value& value) {
    switch (value.type()) {
    case Value::Type::Null:
        return false;
    case Value::Type::Bool:
        return *value.get<bool>();
    case Value::Type::Int:
        return *value.get<std::int64_t>() != 0;
    case Value::Type::Float:
        return *value.get<double>() != 0.0;  // NaN is true
    case Value::Type::String:
        return !value.get<std::string>()->empty();
    case Value::Type::DateTime:
        return true;  // datetime, date and time objects are always true
    case Value::Type::Array:
    case Value::Type::Table:
        return value.size() != 0;
    }
    return true;
}

bool isNumber(const Value& v) {
    const Value::Type t = v.type();
    return t == Value::Type::Bool || t == Value::Type::Int || t == Value::Type::Float;
}

std::int64_t integerOf(const Value& v) {
    if (const bool* b = v.get<bool>()) {
        return *b ? 1 : 0;
    }
    return *v.get<std::int64_t>();
}

std::int64_t utcMinutes(const DateTime& t) {
    // days from civil (proleptic Gregorian), enough to order instants
    const int y = t.year - (t.month <= 2 ? 1 : 0);
    const int era = (y >= 0 ? y : y - 399) / 400;
    const int yoe = y - era * 400;
    const int doy = (153 * (t.month + (t.month > 2 ? -3 : 9)) + 2) / 5 + t.day - 1;
    const int doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    const std::int64_t days = std::int64_t{era} * 146097 + doe;
    return (days * 24 + t.hour) * 60 + t.minute - t.offsetMinutes;
}

bool dateTimesEqual(const DateTime& a, const DateTime& b) {
    if (a.kind != b.kind) {
        return false;  // date != datetime; naive != aware
    }
    if (a.kind == DateTime::Kind::OffsetDateTime) {
        return utcMinutes(a) == utcMinutes(b) && a.second == b.second &&
               a.microsecond == b.microsecond;
    }
    return a == b;
}

// Python's ==: bool, int and float compare as numbers; lists element-wise; dicts by keys
bool pythonEqual(const Value& a, const Value& b) {
    if (isNumber(a) && isNumber(b)) {
        const double* fa = a.get<double>();
        const double* fb = b.get<double>();
        if (!fa && !fb) {
            return integerOf(a) == integerOf(b);
        }
        if (fa && fb) {
            return *fa == *fb;
        }
        const double f = fa ? *fa : *fb;
        const std::int64_t i = integerOf(fa ? b : a);
        // exact comparison: only an integral double within the 64-bit range can equal an int
        if (!(f >= -0x1p63 && f < 0x1p63) || f != std::trunc(f)) {
            return false;
        }
        return static_cast<std::int64_t>(f) == i;
    }
    if (a.type() != b.type()) {
        return false;
    }
    switch (a.type()) {
    case Value::Type::Null:
        return true;
    case Value::Type::String:
        return *a.get<std::string>() == *b.get<std::string>();
    case Value::Type::DateTime:
        return dateTimesEqual(*a.get<DateTime>(), *b.get<DateTime>());
    case Value::Type::Array: {
        const Array& x = *a.get<Array>();
        const Array& y = *b.get<Array>();
        if (x.size() != y.size()) {
            return false;
        }
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (!pythonEqual(x[i], y[i])) {
                return false;
            }
        }
        return true;
    }
    case Value::Type::Table: {
        const Table& x = *a.get<Table>();
        const Table& y = *b.get<Table>();
        if (x.size() != y.size()) {
            return false;
        }
        for (const auto& [key, value] : x) {
            const Value* other = lookup(y, key);
            if (!other || !pythonEqual(value, *other)) {
                return false;
            }
        }
        return true;
    }
    default:
        return false;
    }
}

// format(x, "g")
std::string formatG(double x) {
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%g", x);
    return buffer;
}

std::string issueText(const std::string& message) { return message + " - using the default"; }

// projectfile._option_text: f"{label}: {text}".strip(": ").strip()
std::string optionText(const std::string& label, std::string text,
                       const std::vector<Choice>& choices) {
    if (!choices.empty()) {
        std::vector<std::string> values;
        for (const Choice& c : choices) {
            values.push_back(c.value);
        }
        text += " Choices: " + joined(values) + ".";
    }
    std::string out = label + ": " + text;
    const auto first = out.find_first_not_of(": ");
    if (first == std::string::npos) {
        return {};
    }
    out = out.substr(first, out.find_last_not_of(": ") - first + 1);
    return pythonStrip(out);
}

}  // namespace

std::string_view optionKindName(OptionKind kind) {
    switch (kind) {
    case OptionKind::Bool:
        return "bool";
    case OptionKind::Choice:
        return "choice";
    case OptionKind::Radio:
        return "radio";
    case OptionKind::Int:
        return "int";
    case OptionKind::Float:
        return "float";
    case OptionKind::Text:
        return "text";
    case OptionKind::Multiline:
        return "multiline";
    case OptionKind::File:
        return "file";
    }
    return "text";
}

bool requirementsMet(const Requirements& requirements, const Table& values) {
    static const Value none;
    for (const auto& [key, wanted] : requirements) {
        const Value* found = lookup(values, key);
        const Value& have = found ? *found : none;
        if (const bool* flag = wanted.get<bool>()) {
            if (truthy(have) != *flag) {
                return false;
            }
        } else if (const Array* items = wanted.get<Array>()) {
            bool any = false;
            for (const Value& item : *items) {
                any = any || pythonEqual(have, item);
            }
            if (!any) {
                return false;
            }
        } else if (!pythonEqual(have, wanted)) {
            return false;
        }
    }
    return true;
}

Value Option::coerce(const Value& value) const {
    switch (kind) {
    case OptionKind::Bool: {
        if (value.type() == Value::Type::Bool) {
            return value;
        }
        if (const std::string* text = value.get<std::string>()) {
            // str.lower(): no other character lowers to these ASCII words
            std::string low = *text;
            for (char& c : low) {
                if (c >= 'A' && c <= 'Z') {
                    c = static_cast<char>(c - 'A' + 'a');
                }
            }
            if (low == "true" || low == "yes" || low == "on" || low == "1") {
                return true;
            }
            if (low == "false" || low == "no" || low == "off" || low == "0" || low.empty()) {
                return false;
            }
        }
        if (const std::int64_t* i = value.get<std::int64_t>(); i && (*i == 0 || *i == 1)) {
            return *i == 1;
        }
        if (const double* f = value.get<double>(); f && (*f == 0.0 || *f == 1.0)) {
            return *f == 1.0;
        }
        throw std::invalid_argument(key + ": expected true/false, got " + pythonRepr(value));
    }
    case OptionKind::Int:
    case OptionKind::Float: {
        const std::string* text = value.get<std::string>();
        if (value.isNull() || (text && pythonStrip(*text).empty())) {
            if (optional) {
                return Value();
            }
            throw std::invalid_argument(key + ": a value is required");
        }
        const std::string notNumber = key + ": expected a number, got " + pythonRepr(value);
        double number = 0;
        if (const std::int64_t* i = value.get<std::int64_t>()) {
            number = static_cast<double>(*i);
        } else if (const double* f = value.get<double>()) {
            number = *f;
        } else if (text) {
            const std::optional<double> parsed = tryPythonFloat(*text);
            if (!parsed) {
                throw std::invalid_argument(notNumber);
            }
            number = *parsed;
        } else {
            throw std::invalid_argument(notNumber);  // bool, datetime, list, table
        }
        if (kind == OptionKind::Int) {
            // int(number): Python raises ValueError for NaN and OverflowError (which normalise
            // does not catch) for infinity; both are reported like the NaN case here
            if (std::isnan(number)) {
                throw std::invalid_argument("cannot convert float NaN to integer");
            }
            if (std::isinf(number)) {
                throw std::invalid_argument("cannot convert float infinity to integer");
            }
            if (number != std::trunc(number)) {
                throw std::invalid_argument(key + ": expected an integer, got " +
                                            pythonRepr(value));
            }
        }
        if (minimum && number < *minimum) {
            throw std::invalid_argument(key + ": must be >= " + formatG(*minimum));
        }
        if (maximum && number > *maximum) {
            throw std::invalid_argument(key + ": must be <= " + formatG(*maximum));
        }
        if (kind == OptionKind::Float) {
            return number;
        }
        // Python keeps a big integer; a Value cannot hold one
        if (!(number >= -0x1p63 && number < 0x1p63)) {
            throw std::invalid_argument(key + ": expected an integer, got " + pythonRepr(value));
        }
        return static_cast<std::int64_t>(number);
    }
    case OptionKind::Choice:
    case OptionKind::Radio: {
        const std::string text = value.isNull() ? std::string() : pythonStr(value);
        for (const Choice& c : choices) {
            if (c.value == text) {
                return text;
            }
        }
        std::vector<std::string> values;
        for (const Choice& c : choices) {
            values.push_back(c.value);
        }
        throw std::invalid_argument(key + ": " + repr(text) + " is not one of " + joined(values));
    }
    case OptionKind::Text:
    case OptionKind::Multiline:
    case OptionKind::File:
        break;
    }
    return value.isNull() ? std::string() : pythonStr(value);
}

const Option* Section::option(std::string_view optionKey) const {
    for (const Option& o : options) {
        if (o.key == optionKey) {
            return &o;
        }
    }
    return nullptr;
}

Table Section::recordDefaults() const {
    Table out;
    for (const Option& o : options) {
        assign(out, o.key, o.defaultValue);
    }
    return out;
}

Table Section::recordTemplate(std::string_view id) const {
    for (const RecordTemplate& t : templates) {
        if (t.id == id) {
            Table record = recordDefaults();
            for (const auto& [k, v] : t.values) {
                assign(record, k, v);
            }
            return record;
        }
    }
    throw std::out_of_range(key + ": unknown template " + repr(id));
}

Array Section::defaultRecords() const {
    Array out;
    for (const std::string& id : defaultItems) {
        out.emplace_back(recordTemplate(id));
    }
    return out;
}

std::string Issue::text() const {
    std::string upper = level;
    for (char& c : upper) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return upper + ": " + message;
}

const Section* StepType::section(std::string_view key) const {
    for (const Section& s : sections) {
        if (s.key == key) {
            return &s;
        }
    }
    return nullptr;
}

const Option* StepType::flatOption(std::string_view key) const {
    const Option* found = nullptr;  // a dict: the last option of a key wins
    for (const Section& s : sections) {
        if (!s.repeat) {
            for (const Option& o : s.options) {
                if (o.key == key) {
                    found = &o;
                }
            }
        }
    }
    return found;
}

Table StepType::defaults() const {
    Table out;
    for (const Section& s : sections) {
        if (s.repeat) {
            assign(out, s.key, s.defaultRecords());
        } else {
            for (const Option& o : s.options) {
                assign(out, o.key, o.defaultValue);
            }
        }
    }
    return out;
}

Table StepType::normalise(const Table& settings, std::vector<Issue>* issues,
                          const std::function<void(const std::string&)>& warn) const {
    const auto warning = [&warn](const std::string& message) {
        if (warn) {
            warn(message);
        } else {
            logWarning(message);
        }
    };
    const auto bad = [&](const std::string& key, const std::string& message) {
        if (issues) {
            issues->push_back(Issue{"warning", issueText(message), {key}});
        } else {
            warning(key + ": " + issueText(message));
        }
    };
    const auto coerceInto = [&](Table& target, const Option& o, const Value& raw,
                                const std::string& key) {
        try {
            assign(target, o.key, o.coerce(raw));
        } catch (const std::invalid_argument& e) {
            bad(key, e.what());
        }
    };

    Table out = defaults();
    std::set<std::string> known;
    for (const auto& [k, v] : out) {
        known.insert(k);
    }
    for (const Section& s : sections) {
        if (!s.repeat) {
            for (const Option& o : s.options) {
                if (const Value* raw = lookup(settings, o.key)) {
                    coerceInto(out, o, *raw, o.key);
                }
            }
            continue;
        }
        const Value* raw = lookup(settings, s.key);
        if (!raw || raw->isNull()) {
            continue;
        }
        const Array* list = raw->get<Array>();
        if (!list) {
            bad(s.key, "expected a list");
            continue;
        }
        Array records;
        for (std::size_t k = 0; k < list->size(); ++k) {
            const std::string prefix = s.key + "." + std::to_string(k);
            const Table* record = (*list)[k].get<Table>();
            if (!record) {
                bad(prefix, "expected an object");
                continue;
            }
            Table clean = s.recordDefaults();
            for (const Option& o : s.options) {
                if (const Value* value = lookup(*record, o.key)) {
                    coerceInto(clean, o, *value, prefix + "." + o.key);
                }
            }
            records.emplace_back(std::move(clean));
        }
        assign(out, s.key, std::move(records));
    }
    std::set<std::string> unknown;  // sorted by code point, as Python sorts str
    for (const auto& [k, v] : settings) {
        if (!known.count(k)) {
            unknown.insert(k);
        }
    }
    if (!unknown.empty()) {
        warning(name + ": ignoring unknown settings " +
                joined(std::vector<std::string>(unknown.begin(), unknown.end())));
    }
    return out;
}

const std::vector<StepType>& stepTypes() { return registry().types; }

const StepType* findStepType(std::string_view name) {
    for (const StepType& t : stepTypes()) {
        if (t.name == name) {
            return &t;
        }
    }
    return nullptr;
}

const StepType& stepType(std::string_view name) {
    if (const StepType* t = findStepType(name)) {
        return *t;
    }
    std::set<std::string> names;
    for (const StepType& t : stepTypes()) {
        names.insert(t.name);
    }
    throw std::out_of_range("unknown step type " + repr(name) + " (known: " +
                            joined(std::vector<std::string>(names.begin(), names.end())) + ")");
}

void registerStepType(StepType type) {
    Registry& r = registry();
    const std::lock_guard<std::mutex> lock(r.mutex);
    bool taken = type.name.empty();
    for (const StepType& t : r.types) {
        taken = taken || t.name == type.name;
    }
    if (taken) {
        throw std::invalid_argument("bad or duplicate tool name " + repr(type.name));
    }
    std::vector<std::string> keys;
    for (const Section& s : type.sections) {
        if (!s.repeat) {
            for (const Option& o : s.options) {
                keys.push_back(o.key);
            }
        }
    }
    for (const Section& s : type.sections) {
        if (s.repeat) {
            keys.push_back(s.key);
        }
    }
    std::set<std::string> seen, dupes;
    for (const std::string& k : keys) {
        if (!seen.insert(k).second) {
            dupes.insert(k);
        }
    }
    if (!dupes.empty()) {
        Array list;
        for (const std::string& k : dupes) {
            list.emplace_back(k);
        }
        throw std::invalid_argument(type.name + ": duplicate option keys " + pythonRepr(list));
    }
    r.types.push_back(std::move(type));
}

const std::vector<std::string>& stepCategories() { return registry().categories; }

std::vector<std::pair<std::string, std::vector<const StepType*>>> stepTypesByCategory() {
    std::vector<std::pair<std::string, std::vector<const StepType*>>> out;
    const std::vector<std::string>& categories = stepCategories();
    for (const std::string& category : categories) {
        std::vector<const StepType*> types;
        for (const StepType& t : stepTypes()) {
            if (t.category == category) {
                types.push_back(&t);
            }
        }
        if (!types.empty()) {
            out.emplace_back(category, std::move(types));
        }
    }
    std::vector<const StepType*> rest;
    for (const StepType& t : stepTypes()) {
        bool listed = false;
        for (const std::string& category : categories) {
            listed = listed || t.category == category;
        }
        if (!listed) {
            rest.push_back(&t);
        }
    }
    if (!rest.empty()) {
        out.emplace_back("Other", std::move(rest));
    }
    return out;
}

std::optional<std::string> settingComment(std::string_view tool, std::string_view key) {
    const StepType* type = findStepType(tool);
    if (!type) {
        return std::nullopt;
    }
    for (const Section& s : type->sections) {
        if (s.key == key && s.repeat) {
            return optionText({}, s.tooltip, {});  // a Section has no label or choices
        }
        for (const Option& o : s.options) {
            if (o.key == key) {
                return optionText(o.label, o.tooltip, o.choices);
            }
        }
    }
    return std::nullopt;
}

const std::vector<FileKindInfo>& fileKinds() { return registry().kinds; }

std::string kindLabel(std::string_view kind) {
    for (const FileKindInfo& k : fileKinds()) {
        if (k.id == kind) {
            return k.label;
        }
    }
    return std::string(kind);
}

std::string guessKind(std::string_view path) {
    // str.lower(): of the characters that lower to ASCII only the Kelvin sign (to "k") is not
    // ASCII itself, and no extension holds a "k"
    std::string low(path);
    for (char& c : low) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    for (const FileKindInfo& k : fileKinds()) {
        for (const std::string& ext : k.extensions) {
            std::string e = ext;
            for (char& c : e) {
                if (c >= 'A' && c <= 'Z') {
                    c = static_cast<char>(c - 'A' + 'a');
                }
            }
            if (low.size() >= e.size() && low.compare(low.size() - e.size(), e.size(), e) == 0) {
                return k.id;
            }
        }
    }
    return "text";
}

}  // namespace mdbridge::core
