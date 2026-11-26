// The edit model of the GUI (project_document.hpp). The group rename cascade ports MDATK 0.1.1's
// Project.rename_group / group_references (config.py, selections.py rename_group_tokens; owner's
// code).

#include "mdbridge/core/project_document.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <utility>

#include <toml++/toml.hpp>

#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/ntpath.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/settings.hpp"
#include "mdbridge/core/step_types.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/core/toml_writer.hpp"
#include "mdbridge/version.hpp"
#include "project_internal.hpp"

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

namespace mdbridge::core {

namespace fs = std::filesystem;

namespace {

#include "project_model.inc"

template <class T>
concept Model = requires(const T* p) { fieldsOf(p); };

template <class T>
struct IsVector : std::false_type {};
template <class T>
struct IsVector<std::vector<T>> : std::true_type {
    using Item = T;
};

template <class T>
struct IsOptional : std::false_type {};
template <class T>
struct IsOptional<std::optional<T>> : std::true_type {
    using Item = T;
};

// Calls f(name, member) for every field, in declaration order.
template <class S, class F>
void forEachField(S& s, F&& f) {
    using T = std::remove_const_t<S>;
    std::apply([&](const auto&... fd) { (f(fd.name, s.*(fd.member)), ...); },
               fieldsOf(static_cast<const T*>(nullptr)));
}

// Calls f(member) for the field named `key`; false when there is none.
template <class S, class F>
bool withField(S& s, std::string_view key, F&& f) {
    using T = std::remove_const_t<S>;
    return std::apply([&](const auto&... fd) { return ((fd.name == key ? (f(s.*(fd.member)), true) : false) || ...); },
                      fieldsOf(static_cast<const T*>(nullptr)));
}

std::string repr(std::string_view text) { return pythonRepr(Value(std::string(text))); }

std::string joinPath(std::string_view base, std::string_view key) {
    if (base.empty()) {
        return std::string(key);
    }
    std::string out(base);
    out += '.';
    out += key;
    return out;
}

std::vector<std::string> splitPath(std::string_view path) {
    std::vector<std::string> parts;
    if (path.empty()) {
        return parts;
    }
    std::size_t start = 0;
    for (;;) {
        const std::size_t end = path.find('.', start);
        parts.emplace_back(path.substr(start, end == std::string_view::npos ? std::string_view::npos : end - start));
        if (end == std::string_view::npos) {
            return parts;
        }
        start = end + 1;
    }
}

std::optional<std::size_t> indexOf(std::string_view part) {
    if (part.empty() || part.size() > 9) {
        return std::nullopt;
    }
    std::size_t n = 0;
    for (const char c : part) {
        if (c < '0' || c > '9') {
            return std::nullopt;
        }
        n = n * 10 + static_cast<std::size_t>(c - '0');
    }
    return n;
}

// ---- kinds -----------------------------------------------------------------------------------

template <class T>
constexpr FieldKind kindOf() {
    if constexpr (std::is_same_v<T, std::string>) {
        return FieldKind::Text;
    } else if constexpr (std::is_same_v<T, std::int64_t>) {
        return FieldKind::Integer;
    } else if constexpr (std::is_same_v<T, double>) {
        return FieldKind::Number;
    } else if constexpr (std::is_same_v<T, bool>) {
        return FieldKind::Boolean;
    } else if constexpr (std::is_same_v<T, std::optional<std::int64_t>>) {
        return FieldKind::OptionalInteger;
    } else if constexpr (std::is_same_v<T, std::optional<double>>) {
        return FieldKind::OptionalNumber;
    } else if constexpr (std::is_same_v<T, std::optional<std::string>>) {
        return FieldKind::OptionalText;
    } else if constexpr (std::is_same_v<T, std::optional<bool>>) {
        return FieldKind::OptionalBoolean;
    } else if constexpr (std::is_same_v<T, std::vector<std::string>>) {
        return FieldKind::TextList;
    } else if constexpr (std::is_same_v<T, std::vector<double>>) {
        return FieldKind::NumberList;
    } else if constexpr (std::is_same_v<T, std::vector<std::vector<double>>>) {
        return FieldKind::NumberRows;
    } else if constexpr (std::is_same_v<T, Table>) {
        return FieldKind::FreeTable;
    } else if constexpr (Model<T>) {
        return FieldKind::Record;
    } else {
        static_assert(IsVector<T>::value, "a field type of the model without a kind");
        return FieldKind::Records;
    }
}

FieldInfo makeInfo(FieldKind kind, std::string modelPath, bool free) {
    FieldInfo info;
    info.kind = kind;
    info.comment = free ? std::nullopt : modelComment(modelPath);
    info.modelPath = std::move(modelPath);
    return info;
}

std::optional<FieldInfo> freeInfo(const Value& value, const std::vector<std::string>& parts, std::size_t i,
                                  std::string modelPath) {
    if (i == parts.size()) {
        return makeInfo(value.isTable() ? FieldKind::FreeTable : FieldKind::FreeValue, std::move(modelPath), true);
    }
    if (const Array* array = value.get<Array>()) {
        const auto index = indexOf(parts[i]);
        if (!index || *index >= array->size()) {
            return std::nullopt;
        }
        return freeInfo((*array)[*index], parts, i + 1, std::move(modelPath));
    }
    if (const Value* entry = value.find(parts[i])) {
        return freeInfo(*entry, parts, i + 1, joinPath(modelPath, parts[i]));
    }
    return std::nullopt;
}

template <class T>
std::optional<FieldInfo> infoOf(const T& member, const std::vector<std::string>& parts, std::size_t i,
                                std::string modelPath) {
    if (i == parts.size()) {
        return makeInfo(kindOf<T>(), std::move(modelPath), false);
    }
    if constexpr (Model<T>) {
        std::optional<FieldInfo> out;
        withField(member, parts[i], [&](const auto& m) { out = infoOf(m, parts, i + 1, joinPath(modelPath, parts[i])); });
        return out;
    } else if constexpr (std::is_same_v<T, Table>) {
        return freeInfo(Value(member), parts, i, std::move(modelPath));
    } else if constexpr (IsVector<T>::value) {
        const auto index = indexOf(parts[i]);
        if (!index || *index >= member.size()) {
            return std::nullopt;
        }
        return infoOf(member[*index], parts, i + 1, std::move(modelPath));
    } else {
        return std::nullopt;
    }
}

template <class T>
void collectPaths(const T& member, const std::string& path, std::vector<std::string>& out) {
    if constexpr (Model<T>) {
        forEachField(member, [&](std::string_view name, const auto& m) { collectPaths(m, joinPath(path, name), out); });
    } else if constexpr (IsVector<T>::value) {
        if constexpr (Model<typename IsVector<T>::Item>) {
            for (std::size_t i = 0; i < member.size(); ++i) {
                collectPaths(member[i], joinPath(path, std::to_string(i)), out);
            }
        } else {
            out.push_back(path);
        }
    } else {
        out.push_back(path);
    }
}

// ---- plain values ----------------------------------------------------------------------------

bool sameValue(const Value& a, const Value& b) {
    if (a.type() != b.type()) {
        return false;
    }
    if (const double* x = a.get<double>()) {
        const double y = *b.get<double>();
        if (std::isnan(*x) || std::isnan(y)) {
            return std::isnan(*x) && std::isnan(y);
        }
        return *x == y && std::signbit(*x) == std::signbit(y);
    }
    if (const Array* x = a.get<Array>()) {
        const Array& y = *b.get<Array>();
        return x->size() == y.size() && std::equal(x->begin(), x->end(), y.begin(), sameValue);
    }
    if (const Table* x = a.get<Table>()) {
        const Table& y = *b.get<Table>();
        return x->size() == y.size() && std::equal(x->begin(), x->end(), y.begin(), [](const auto& p, const auto& q) {
                   return p.first == q.first && sameValue(p.second, q.second);
               });
    }
    return a == b;
}

// The fields whose values differ, as paths below `prefix`.
void diffValues(const Value& a, const Value& b, const std::string& prefix, std::vector<std::string>& out) {
    if (sameValue(a, b)) {
        return;
    }
    const Table* x = a.get<Table>();
    const Table* y = b.get<Table>();
    if (x && y) {
        bool sameKeys = x->size() == y->size();
        for (std::size_t i = 0; sameKeys && i < x->size(); ++i) {
            sameKeys = (*x)[i].first == (*y)[i].first;
        }
        if (sameKeys) {
            for (std::size_t i = 0; i < x->size(); ++i) {
                diffValues((*x)[i].second, (*y)[i].second, joinPath(prefix, (*x)[i].first), out);
            }
            return;
        }
    }
    const Array* u = a.get<Array>();
    const Array* v = b.get<Array>();
    if (u && v && u->size() == v->size()) {
        for (std::size_t i = 0; i < u->size(); ++i) {
            diffValues((*u)[i], (*v)[i], joinPath(prefix, std::to_string(i)), out);
        }
        return;
    }
    out.push_back(prefix);
}

std::vector<std::string> changedFields(const Project& before, const Project& after) {
    std::vector<std::string> out;
    diffValues(toPlain(before), toPlain(after), "", out);
    return out;
}

// "project.analyses.contacts.cutoff: expected a number, got 'far'" -> "expected a number, got 'far'"
std::string withoutWhere(const std::string& message) {
    if (message.rfind("project", 0) == 0) {
        const std::size_t colon = message.find(": ");
        if (colon != std::string::npos && message.find(' ') > colon) {
            return message.substr(colon + 2);
        }
    }
    return message;
}

[[noreturn]] void noSuchField(std::string_view path) {
    throw FieldError(std::string(path), "no setting " + repr(path));
}

FieldInfo requireInfo(const Project& project, std::string_view path) {
    auto info = fieldInfo(project, path);
    if (!info) {
        noSuchField(path);
    }
    return *info;
}

// The plain value at `parts` (nullptr when absent).
Value* findPlain(Value& root, const std::vector<std::string>& parts, std::size_t count) {
    Value* node = &root;
    for (std::size_t i = 0; i < count && node; ++i) {
        if (Array* array = node->get<Array>()) {
            const auto index = indexOf(parts[i]);
            node = index && *index < array->size() ? &(*array)[*index] : nullptr;
        } else {
            node = node->find(parts[i]);
        }
    }
    return node;
}

std::string pythonFloatList(const Value& value) {
    std::string out;
    if (const Array* array = value.get<Array>()) {
        for (std::size_t i = 0; i < array->size(); ++i) {
            if (i) {
                out += ", ";
            }
            const Value& item = (*array)[i];
            if (const double* d = item.get<double>()) {
                out += pythonFloatRepr(*d);
            } else {
                out += pythonStr(item);
            }
        }
    }
    return out;
}

std::string asciiLower(std::string_view text) {
    std::string out(text);
    for (char& c : out) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return out;
}

std::vector<std::string> splitCsv(std::string_view text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (;;) {
        const std::size_t end = text.find(',', start);
        std::string item = pythonStrip(text.substr(start, end == std::string_view::npos ? std::string_view::npos : end - start));
        if (!item.empty()) {
            out.push_back(std::move(item));
        }
        if (end == std::string_view::npos) {
            return out;
        }
        start = end + 1;
    }
}

Value parseJsonText(std::string_view path, std::string_view text) {
    try {
        return parseJson(text, "");
    } catch (const ProjectError& error) {
        const std::string message = error.what();
        const std::size_t at = message.find("not valid JSON");
        throw FieldError(std::string(path), at == std::string::npos ? message : message.substr(at));
    }
}

}  // namespace

// ---- fields ----------------------------------------------------------------------------------

Problem FieldError::problem(const Project* project) const { return locatedProblem(m_field, what(), project); }

std::optional<FieldInfo> fieldInfo(const Project& project, std::string_view path) {
    return infoOf(project, splitPath(path), 0, std::string());
}

std::vector<std::string> fieldPaths(const Project& project) {
    std::vector<std::string> out;
    collectPaths(project, std::string(), out);
    return out;
}

std::string fieldText(const Project& project, std::string_view path) {
    const FieldInfo info = requireInfo(project, path);
    Value plain = toPlain(project);
    const Value* value = findPlain(plain, splitPath(path), splitPath(path).size());
    const Value none;
    const Value& v = value ? *value : none;
    switch (info.kind) {
    case FieldKind::Text:
    case FieldKind::OptionalText:
        return v.stringOr();
    case FieldKind::Integer:
    case FieldKind::OptionalInteger:
        return v.isNull() ? std::string() : pythonStr(v);
    case FieldKind::Number:
    case FieldKind::OptionalNumber:
        if (const double* d = v.get<double>()) {
            return pythonFloatRepr(*d);
        }
        return v.isNull() ? std::string() : pythonStr(v);
    case FieldKind::Boolean:
    case FieldKind::OptionalBoolean:
        if (const bool* b = v.get<bool>()) {
            return *b ? "true" : "false";
        }
        return {};
    case FieldKind::TextList: {
        std::string out;
        if (const Array* array = v.get<Array>()) {
            for (std::size_t i = 0; i < array->size(); ++i) {
                out += (i ? ", " : "") + (*array)[i].stringOr();
            }
        }
        return out;
    }
    case FieldKind::NumberList:
        return pythonFloatList(v);
    case FieldKind::FreeTable:
        return v.size() == 0 ? std::string() : toJson(v, -1);
    case FieldKind::FreeValue:
        return v.isString() ? v.stringOr() : toJson(v, -1);
    case FieldKind::NumberRows:
    case FieldKind::Records:
    case FieldKind::Record:
        break;
    }
    throw FieldError(std::string(path), repr(path) + " has no text form");
}

Value parseFieldText(const Project& project, std::string_view path, std::string_view text) {
    const FieldInfo info = requireInfo(project, path);
    const std::string stripped = pythonStrip(text);
    const auto fail = [&](const std::string& message) -> FieldError { return FieldError(std::string(path), message); };
    switch (info.kind) {
    case FieldKind::Text:
        return std::string(text);
    case FieldKind::OptionalText: {  // "" is unset for a field that is not set (a set "" stays "")
        const std::vector<std::string> parts = splitPath(path);
        Value plain = toPlain(project);
        const Value* current = findPlain(plain, parts, parts.size());
        return text.empty() && (!current || current->isNull()) ? Value() : Value(std::string(text));
    }
    case FieldKind::Integer:
    case FieldKind::OptionalInteger:
        if (stripped.empty() && info.kind == FieldKind::OptionalInteger) {
            return {};
        }
        if (const auto n = tryPythonInt(stripped)) {
            return *n;
        }
        throw fail(repr(stripped) + " is not a whole number");
    case FieldKind::Number:
    case FieldKind::OptionalNumber:
        if (stripped.empty() && info.kind == FieldKind::OptionalNumber) {
            return {};
        }
        if (const auto d = tryPythonFloat(stripped)) {
            return *d;
        }
        throw fail(repr(stripped) + " is not a number");
    case FieldKind::Boolean:
    case FieldKind::OptionalBoolean: {
        if (stripped.empty() && info.kind == FieldKind::OptionalBoolean) {
            return {};
        }
        const std::string lower = asciiLower(stripped);
        if (lower == "true" || lower == "false") {
            return lower == "true";
        }
        throw fail(repr(stripped) + " is not true or false");
    }
    case FieldKind::TextList: {
        Array out;
        for (std::string& item : splitCsv(text)) {
            out.emplace_back(std::move(item));
        }
        return out;
    }
    case FieldKind::NumberList: {
        Array out;
        for (const std::string& item : splitCsv(text)) {
            const auto d = tryPythonFloat(item);
            if (!d) {
                throw fail(repr(item) + " is not a number");
            }
            out.emplace_back(*d);
        }
        return out;
    }
    case FieldKind::FreeTable: {
        if (stripped.empty()) {
            return Table{};
        }
        Value value = parseJsonText(path, text);
        if (!value.isTable()) {
            const std::vector<std::string> parts = splitPath(path);
            throw fail((parts.empty() ? std::string("value") : parts.back()) + " must be a JSON object");
        }
        return value;
    }
    case FieldKind::FreeValue: {
        Value plain = toPlain(project);
        const std::vector<std::string> parts = splitPath(path);
        const Value* current = findPlain(plain, parts, parts.size());
        if (current && current->isString()) {
            return std::string(text);
        }
        return parseJsonText(path, text);
    }
    case FieldKind::NumberRows:
    case FieldKind::Records:
    case FieldKind::Record:
        break;
    }
    throw fail(repr(path) + " has no text form");
}

bool DocumentChange::affects(std::string_view path) const {
    for (const std::string& f : fields) {
        if (f.empty() || path.empty() || f == path) {
            return true;
        }
        const auto inside = [](std::string_view a, std::string_view b) {  // a inside b
            return a.size() > b.size() && a.compare(0, b.size(), b) == 0 && a[b.size()] == '.';
        };
        if (inside(f, path) || inside(path, f)) {
            return true;
        }
    }
    return false;
}

// ---- the document ----------------------------------------------------------------------------

namespace {

std::atomic<std::uint64_t> g_documents{0};

Project projectFromPlainQuietly(const Value& plain) {
    std::vector<std::string> ignored;
    return projectFromPlain(plain, {}, &ignored);
}

std::string settingsText(const fs::path& path) {
    const std::string text = detail::readTextFile(path);
    const std::size_t at = detail::logStart(text);
    return at == std::string::npos ? text : text.substr(0, at);
}

// ---- the comments of a file ----

// The lines (1-based) that continue a value begun on an earlier line (multi-line strings and arrays).
void markContinuations(const toml::node& node, std::vector<bool>& inside) {
    const auto mark = [&inside](const toml::node& value) {
        const auto& where = value.source();
        for (auto line = where.begin.line + 1; line <= where.end.line && line < inside.size(); ++line) {
            inside[line] = true;
        }
    };
    if (const toml::table* table = node.as_table()) {
        if (table->is_inline()) {
            mark(node);
            return;
        }
        for (auto&& [key, value] : *table) {
            markContinuations(value, inside);
        }
    } else if (const toml::array* array = node.as_array()) {
        if (array->is_array_of_tables() && !array->empty() && !array->front().as_table()->is_inline()) {
            for (const toml::node& item : *array) {
                markContinuations(item, inside);
            }
        } else {
            mark(node);
        }
    } else {
        mark(node);
    }
}

// The keys of a TOML key path ("a.\"b c\".d"), unquoted (basic strings with simple escapes).
std::vector<std::string> keyPath(std::string_view text) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && (text[i] == ' ' || text[i] == '\t' || text[i] == '.')) {
            ++i;
        }
        if (i >= text.size()) {
            break;
        }
        std::string key;
        if (text[i] == '"' || text[i] == '\'') {
            const char quote = text[i++];
            while (i < text.size() && text[i] != quote) {
                if (quote == '"' && text[i] == '\\' && i + 1 < text.size()) {
                    ++i;
                }
                key += text[i++];
            }
            ++i;
        } else {
            while (i < text.size() && text[i] != '.' && text[i] != ' ' && text[i] != '\t') {
                key += text[i++];
            }
        }
        out.push_back(std::move(key));
    }
    return out;
}

std::string joinKeys(const std::vector<std::string>& keys) {
    std::string out;
    for (std::size_t i = 0; i < keys.size(); ++i) {
        out += (i ? "." : "") + keys[i];
    }
    return out;
}

bool isNodeTable(const std::vector<std::string>& path) {
    return path.size() == 2 && path[0] == "workflow" && (path[1] == "nodes" || path[1] == "scratch");
}

bool isNodeSetting(const std::vector<std::string>& path) {
    return path.size() == 4 && path[0] == "workflow" && (path[1] == "nodes" || path[1] == "scratch") &&
           path[2] == "settings";
}

// The comment block above every entry and table header of a project file written by tomlio, keyed by
// the path without list indices (planner settings with "@<tool>"); the first one of each key.
std::vector<std::pair<std::string, std::string>> fileComments(const std::string& settings) {
    std::vector<std::pair<std::string, std::string>> out;
    toml::table root;
    try {
        root = toml::parse(settings);
    } catch (const toml::parse_error&) {
        return out;
    }
    std::vector<std::string> lines;
    for (std::size_t start = 0; start <= settings.size();) {
        const std::size_t end = settings.find('\n', start);
        lines.push_back(settings.substr(start, end == std::string::npos ? std::string::npos : end - start));
        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    std::vector<bool> inside(lines.size() + 2, false);
    markContinuations(root, inside);
    std::vector<std::string> table;
    std::string tool;
    std::vector<std::string> pending;
    const auto record = [&](std::vector<std::string> path) {
        std::string key = joinKeys(path);
        if (isNodeSetting(path)) {
            key += "@" + tool;
        }
        if (std::none_of(out.begin(), out.end(), [&](const auto& e) { return e.first == key; })) {
            std::string text;
            for (std::size_t i = 0; i < pending.size(); ++i) {
                text += (i ? "\n" : "") + pending[i];
            }
            out.emplace_back(std::move(key), std::move(text));
        }
        pending.clear();
    };
    static constexpr std::string_view kNotSet = " =   (not set)";
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string& line = lines[i];
        if (inside[i + 1]) {
            pending.clear();
            continue;
        }
        if (line.empty()) {
            pending.clear();
        } else if (line[0] == '#') {
            if (line.size() > 2 + kNotSet.size() && line.compare(0, 2, "# ") == 0 &&
                line.compare(line.size() - kNotSet.size(), kNotSet.size(), kNotSet) == 0) {
                std::vector<std::string> path = table;
                for (std::string& k : keyPath(std::string_view(line).substr(2, line.size() - 2 - kNotSet.size()))) {
                    path.push_back(std::move(k));
                }
                record(std::move(path));
            } else {
                pending.push_back(line.size() >= 2 && line[1] == ' ' ? line.substr(2) : line.substr(1));
            }
        } else if (line[0] == '[') {
            const bool array = line.size() > 1 && line[1] == '[';
            const std::size_t close = line.rfind(array ? "]]" : "]");
            table = keyPath(std::string_view(line).substr(array ? 2 : 1, close == std::string::npos ? 0 : close - (array ? 2 : 1)));
            if (isNodeTable(table)) {
                tool.clear();
            }
            record(table);
        } else {
            const std::size_t equals = line.find(" = ");
            if (equals == std::string::npos) {
                pending.clear();
                continue;
            }
            std::vector<std::string> path = table;
            for (std::string& k : keyPath(std::string_view(line).substr(0, equals))) {
                path.push_back(std::move(k));
            }
            if (isNodeTable(table) && path.size() == 3 && path[2] == "tool") {
                try {
                    const toml::table entry = toml::parse(line);
                    tool = entry["tool"].value_or(std::string());
                } catch (const toml::parse_error&) {
                    tool.clear();
                }
            }
            record(std::move(path));
        }
    }
    return out;
}

}  // namespace

TomlDocFunction ProjectDocument::comments() const {
    if (m_fileComments.empty() || writtenMdatkVersion(m_project) == std::string_view(kMdatkReference)) {
        return {};
    }
    auto entries = std::make_shared<const std::vector<std::pair<std::string, std::string>>>(m_fileComments);
    return [entries](const std::vector<std::string>& path, const std::vector<const Value*>& stack) -> std::optional<std::string> {
        std::string key = joinKeys(path);
        if (isNodeSetting(path)) {  // the step is the table two levels up (as projectDoc reads it)
            const Value* node = stack.size() >= 2 ? stack[stack.size() - 2] : nullptr;
            const Value* tool = node ? node->find("tool") : nullptr;
            key += "@" + (tool ? tool->stringOr() : std::string());
        }
        for (const auto& [name, text] : *entries) {
            if (name == key) {
                return text;
            }
        }
        return projectDoc(path, stack);
    };
}

ProjectDocument::ProjectDocument() : m_serial(++g_documents) {
    m_currentId = m_savedId = m_nextId++;
}

Project ProjectDocument::resolvedProject() const {
    Project out = m_project;
    resolvePaths(out, folder());
    return out;
}

std::string ProjectDocument::folder() const {
    return hasPath() ? detail::absoluteDirectory(m_path) : detail::absolutePath(".");
}

void ProjectDocument::rememberDisk() {
    m_diskSize.reset();
    m_diskBytes.clear();
    std::error_code ec;
    if (!hasPath() || !fs::is_regular_file(m_path, ec)) {
        return;
    }
    const auto size = fs::file_size(m_path, ec);
    const auto time = fs::last_write_time(m_path, ec);
    if (ec) {
        return;
    }
    try {
        m_diskBytes = settingsText(m_path);
    } catch (const std::exception&) {
        return;
    }
    m_diskSize = size;
    m_diskTime = time;
}

bool ProjectDocument::changedOnDisk() const {
    if (!hasPath()) {
        return false;
    }
    std::error_code ec;
    if (!fs::is_regular_file(m_path, ec)) {
        return true;
    }
    if (!m_diskSize) {
        return true;
    }
    const auto size = fs::file_size(m_path, ec);
    const auto time = fs::last_write_time(m_path, ec);
    if (!ec && size == *m_diskSize && time == m_diskTime) {
        return false;
    }
    try {
        return settingsText(m_path) != m_diskBytes;  // a step log appended by a run is no change
    } catch (const std::exception&) {
        return true;
    }
}

void ProjectDocument::newProject() { replace(Project{}, {}, false); }

void ProjectDocument::open(const fs::path& path, std::vector<std::string>* warnings) {
    const ProjectFile file = readProject(path);
    std::vector<std::string> found;
    Project project = projectFromPlain(file.settings, {}, &found);
    const Value* version = file.header.find("mdatk_version");
    project.mdatkVersion = version && version->isString() ? version->stringOr() : std::string();
    Project check = project;  // what loadProject would refuse
    resolvePaths(check, detail::absoluteDirectory(path));
    std::vector<std::pair<std::string, std::string>> commentsOfFile = fileComments(settingsText(path));
    assign(std::move(project), path, false, file.kind);
    m_fileComments = std::move(commentsOfFile);
    if (warnings) {
        warnings->insert(warnings->end(), found.begin(), found.end());
    }
}

void ProjectDocument::replace(Project project, fs::path path, bool modified) {
    assign(std::move(project), std::move(path), modified, FileKind::Project);
}

void ProjectDocument::assign(Project project, fs::path path, bool modified, FileKind kind) {
    const bool wasModified = isModified();
    m_project = std::move(project);
    m_path = path.empty() ? fs::path() : pathFromUtf8(detail::absolutePath(path));
    m_kind = kind;
    m_fileComments.clear();
    m_undo.clear();
    m_redo.clear();
    m_groupDepth = 0;
    m_groupStart.reset();
    m_groupChanged = false;
    m_currentId = m_nextId++;
    m_savedId = modified ? 0 : m_currentId;
    rememberDisk();
    DocumentChange change;
    change.kind = DocumentChange::Kind::Replaced;
    change.fields = {""};
    change.modifiedChanged = wasModified != isModified();
    notify(std::move(change));
}

void ProjectDocument::save() {
    if (!hasPath()) {
        throw ProjectError("the project has no file yet; use Save As");
    }
    saveProject(m_path, m_project, comments());
    const bool wasModified = isModified();
    m_savedId = m_currentId;
    rememberDisk();
    DocumentChange change;
    change.kind = DocumentChange::Kind::Saved;
    change.modifiedChanged = wasModified != isModified();
    notify(std::move(change));
}

void ProjectDocument::saveAs(const fs::path& path) {
    const fs::path target = pathFromUtf8(detail::absolutePath(path));
    const std::string from = folder();
    const std::string to = detail::absoluteDirectory(target);
    Project next = m_project;
    const std::vector<std::string> changed = rebasePaths(next, from, to);
    saveProject(target, next, comments());
    for (Step& step : m_undo) {
        rebasePaths(step.project, from, to);
    }
    for (Step& step : m_redo) {
        rebasePaths(step.project, from, to);
    }
    if (m_groupStart) {
        rebasePaths(m_groupStart->project, from, to);
    }
    const bool wasModified = isModified();
    m_project = std::move(next);
    m_path = target;
    m_kind = FileKind::Project;
    m_savedId = m_currentId;
    rememberDisk();
    DocumentChange change;
    change.kind = DocumentChange::Kind::Saved;
    change.fields = changed;
    change.modifiedChanged = wasModified != isModified();
    notify(std::move(change));
}

Value ProjectDocument::get(std::string_view path) const {
    requireInfo(m_project, path);
    Value plain = toPlain(m_project);
    const std::vector<std::string> parts = splitPath(path);
    const Value* value = findPlain(plain, parts, parts.size());
    return value ? *value : Value();
}

bool ProjectDocument::set(std::string_view path, const Value& value, std::string_view label) {
    const std::vector<std::string> parts = splitPath(path);
    const std::optional<FieldInfo> info = fieldInfo(m_project, path);
    if (parts.empty()) {
        if (!value.isTable()) {
            throw FieldError(std::string(path), "expected a table, got " + pythonRepr(value));
        }
    } else if (!info) {  // only a new entry of a free-form table
        const std::size_t dot = path.rfind('.');
        const auto parentInfo =
            dot == std::string_view::npos ? std::nullopt : fieldInfo(m_project, path.substr(0, dot));
        if (!parentInfo || parentInfo->kind != FieldKind::FreeTable || indexOf(parts.back())) {
            noSuchField(path);
        }
    } else if ((info->kind == FieldKind::Text && !value.isString()) ||
               (info->kind == FieldKind::OptionalText && !value.isString() && !value.isNull())) {
        throw FieldError(std::string(path), "expected a string, got " + pythonRepr(value));
    } else if (info->kind == FieldKind::TextList && value.isArray()) {
        for (const Value& item : *value.get<Array>()) {
            if (!item.isString()) {
                throw FieldError(std::string(path), "expected a list of strings, got " + pythonRepr(value));
            }
        }
    }
    const Value before = toPlain(m_project);
    Value plain = before;
    if (parts.empty()) {
        plain = value;
    } else {
        Value* parent = findPlain(plain, parts, parts.size() - 1);
        if (!parent) {
            noSuchField(path);
        }
        if (Array* array = parent->get<Array>()) {
            (*array)[*indexOf(parts.back())] = value;
        } else if (Table* table = parent->get<Table>()) {
            auto it = std::find_if(table->begin(), table->end(), [&](const auto& e) { return e.first == parts.back(); });
            // an entry of a free-form table (not the table a model field holds): Null removes it
            const bool freeEntry = !info || info->kind == FieldKind::FreeValue ||
                                   (info->kind == FieldKind::FreeTable && !info->comment);
            if (value.isNull() && freeEntry) {
                if (it != table->end()) {
                    table->erase(it);
                }
            } else if (it != table->end()) {
                it->second = value;
            } else {
                table->emplace_back(parts.back(), value);
            }
        } else {
            noSuchField(path);
        }
    }
    Project next;
    try {
        next = projectFromPlainQuietly(plain);
        if (parts.empty() || parts[0] == "workflow") {
            Project check = next;
            resolvePaths(check, "/");  // planner files must be paths, as loadProject requires
        }
    } catch (const ProjectError& error) {
        throw FieldError(std::string(path), withoutWhere(error.what()));
    }
    next.mdatkVersion = m_project.mdatkVersion;
    if (sameValue(toPlain(next), before)) {
        return false;
    }
    commit(std::move(next), label.empty() ? "Set " + describeField(path, &m_project) : std::string(label));
    return true;
}

bool ProjectDocument::setText(std::string_view path, std::string_view text, std::string_view label) {
    if (fieldText(m_project, path) == text) {
        return false;
    }
    return set(path, parseFieldText(m_project, path, text), label);
}

bool ProjectDocument::edit(std::string_view label, const std::function<void(Project&)>& change) {
    Project next = m_project;
    change(next);
    if (sameValue(toPlain(next), toPlain(m_project))) {
        return false;
    }
    commit(std::move(next), label);
    return true;
}

namespace {

Array listAt(const ProjectDocument& doc, std::string_view path) {
    const auto info = doc.info(path);
    if (!info || (info->kind != FieldKind::Records && info->kind != FieldKind::TextList &&
                  info->kind != FieldKind::NumberList && info->kind != FieldKind::NumberRows)) {
        throw FieldError(std::string(path), repr(path) + " is not a list");
    }
    const Value value = doc.get(path);
    return value.get<Array>() ? *value.get<Array>() : Array{};
}

[[noreturn]] void badIndex(std::string_view path, std::size_t index) {
    throw FieldError(std::string(path), "no item " + std::to_string(index) + " in " + repr(path));
}

}  // namespace

void ProjectDocument::insertItem(std::string_view listPath, std::size_t index, const Value& item, std::string_view label) {
    Array list = listAt(*this, listPath);
    if (index > list.size()) {
        badIndex(listPath, index);
    }
    list.insert(list.begin() + static_cast<std::ptrdiff_t>(index), item);
    set(listPath, list, label.empty() ? "Add to " + describeField(listPath, &m_project) : std::string(label));
}

void ProjectDocument::removeItem(std::string_view listPath, std::size_t index, std::string_view label) {
    Array list = listAt(*this, listPath);
    if (index >= list.size()) {
        badIndex(listPath, index);
    }
    list.erase(list.begin() + static_cast<std::ptrdiff_t>(index));
    set(listPath, list, label.empty() ? "Remove from " + describeField(listPath, &m_project) : std::string(label));
}

void ProjectDocument::moveItem(std::string_view listPath, std::size_t from, std::size_t to, std::string_view label) {
    Array list = listAt(*this, listPath);
    if (from >= list.size()) {
        badIndex(listPath, from);
    }
    if (to >= list.size()) {
        badIndex(listPath, to);
    }
    Value item = list[from];
    list.erase(list.begin() + static_cast<std::ptrdiff_t>(from));
    list.insert(list.begin() + static_cast<std::ptrdiff_t>(to), std::move(item));
    set(listPath, list, label.empty() ? "Move in " + describeField(listPath, &m_project) : std::string(label));
}

std::vector<std::string> ProjectDocument::renameGroup(std::string_view oldName, std::string_view newName) {
    Project next = m_project;
    std::vector<std::string> changed = core::renameGroup(next, oldName, newName);
    if (!sameValue(toPlain(next), toPlain(m_project))) {
        commit(std::move(next), "Rename group " + repr(oldName) + " to " + repr(pythonStrip(newName)));
    }
    return changed;
}

std::vector<std::string> ProjectDocument::groupReferences(std::string_view name) const {
    return core::groupReferences(m_project, name);
}

void ProjectDocument::commit(Project next, std::string_view label) {
    const bool wasModified = isModified();
    DocumentChange change;
    change.kind = DocumentChange::Kind::Edit;
    change.fields = changedFields(m_project, next);
    change.label = std::string(label);
    if (m_groupDepth > 0) {
        if (!m_groupChanged) {
            m_redo.clear();
        }
        m_groupChanged = true;
    } else {
        m_undo.push_back(Step{std::move(m_project), m_currentId, std::string(label)});
        if (m_undo.size() > kUndoLimit) {
            m_undo.erase(m_undo.begin());
        }
        m_redo.clear();
    }
    m_project = std::move(next);
    m_currentId = m_nextId++;
    change.modifiedChanged = wasModified != isModified();
    notify(std::move(change));
}

std::string ProjectDocument::undoLabel() const { return m_undo.empty() ? std::string() : m_undo.back().label; }

std::string ProjectDocument::redoLabel() const { return m_redo.empty() ? std::string() : m_redo.back().label; }

void ProjectDocument::undo() {
    while (m_groupDepth > 0) {
        endGroup();
    }
    if (m_undo.empty()) {
        return;
    }
    const bool wasModified = isModified();
    Step step = std::move(m_undo.back());
    m_undo.pop_back();
    DocumentChange change;
    change.kind = DocumentChange::Kind::Undo;
    change.fields = changedFields(m_project, step.project);
    change.label = step.label;
    m_redo.push_back(Step{std::move(m_project), m_currentId, step.label});
    m_project = std::move(step.project);
    m_currentId = step.id;
    change.modifiedChanged = wasModified != isModified();
    notify(std::move(change));
}

void ProjectDocument::redo() {
    while (m_groupDepth > 0) {
        endGroup();
    }
    if (m_redo.empty()) {
        return;
    }
    const bool wasModified = isModified();
    Step step = std::move(m_redo.back());
    m_redo.pop_back();
    DocumentChange change;
    change.kind = DocumentChange::Kind::Redo;
    change.fields = changedFields(m_project, step.project);
    change.label = step.label;
    m_undo.push_back(Step{std::move(m_project), m_currentId, step.label});
    m_project = std::move(step.project);
    m_currentId = step.id;
    change.modifiedChanged = wasModified != isModified();
    notify(std::move(change));
}

void ProjectDocument::clearUndo() {
    m_undo.clear();
    m_redo.clear();
}

void ProjectDocument::beginGroup(std::string_view label) {
    if (m_groupDepth++ == 0) {
        m_groupStart = Step{m_project, m_currentId, std::string(label)};
        m_groupLabel = std::string(label);
        m_groupChanged = false;
    }
}

void ProjectDocument::endGroup() {
    if (m_groupDepth == 0 || --m_groupDepth > 0) {
        return;
    }
    if (m_groupChanged && m_groupStart) {
        m_undo.push_back(std::move(*m_groupStart));
        m_undo.back().label = m_groupLabel;
        if (m_undo.size() > kUndoLimit) {
            m_undo.erase(m_undo.begin());
        }
    }
    m_groupStart.reset();
    m_groupChanged = false;
}

void ProjectDocument::cancelGroup() {
    if (m_groupDepth == 0) {
        return;
    }
    m_groupDepth = 0;
    std::optional<Step> start = std::move(m_groupStart);
    m_groupStart.reset();
    const bool changed = m_groupChanged;
    m_groupChanged = false;
    if (!changed || !start) {
        return;
    }
    const bool wasModified = isModified();
    DocumentChange change;
    change.kind = DocumentChange::Kind::Undo;
    change.fields = changedFields(m_project, start->project);
    change.label = m_groupLabel;
    m_project = std::move(start->project);
    m_currentId = start->id;
    change.modifiedChanged = wasModified != isModified();
    notify(std::move(change));
}

std::vector<Problem> ProjectDocument::problems() const {
    try {
        return projectProblems(resolvedProject());
    } catch (const ProjectError&) {
        return projectProblems(m_project);
    }
}

ProjectDocument::ListenerId ProjectDocument::addListener(Listener listener) {
    const ListenerId id = m_nextListener++;
    m_listeners.emplace_back(id, std::move(listener));
    return id;
}

void ProjectDocument::removeListener(ListenerId id) {
    std::erase_if(m_listeners, [id](const auto& entry) { return entry.first == id; });
}

void ProjectDocument::notify(DocumentChange change) {
    const auto listeners = m_listeners;  // a listener may add or remove listeners
    for (const auto& [id, listener] : listeners) {
        if (listener) {
            listener(change);
        }
    }
}

// ---- residue groups --------------------------------------------------------------------------

namespace {

// GroupUse of MDATK 0.1.1's config.py: how a string field uses group names, by model path
constexpr std::array<std::pair<std::string_view, std::string_view>, 23> kGroupUses{{
    {"notes", "text"},
    {"analyses.rmsd.targets.selection", "selection"},
    {"analyses.rmsd.targets.fit_selection", "selection"},
    {"analyses.contacts.group_a", "name"},
    {"analyses.contacts.group_b", "name"},
    {"analyses.hbonds.group_a", "name"},
    {"analyses.hbonds.group_b", "name"},
    {"analyses.mg.mg_selection", "selection"},
    {"analyses.mg.ligand_selection", "selection"},
    {"analyses.mg.group_a", "name"},
    {"analyses.mg.group_b", "name"},
    {"analyses.mg.occupancy_selection", "selection"},
    {"analyses.stacking.selection", "selection"},
    {"analyses.stacking.group_a", "name"},
    {"analyses.stacking.group_b", "name"},
    {"analyses.energy.group_a", "name"},
    {"analyses.energy.group_b", "name"},
    {"visual.representations.groups", "name"},
    {"visual.structure_groups", "name"},
    {"visual.proximity_set", "selection"},
    {"layout.panels.series", "key"},
    {"output.charts.series", "key"},
    {"output.pdb.selection", "selection"},
}};

std::string_view groupUse(std::string_view modelPath) {
    for (const auto& [path, kind] : kGroupUses) {
        if (path == modelPath) {
            return kind;
        }
    }
    return {};
}

using StringFunction = std::function<std::string(const std::string& field, std::string_view kind, const std::string& text)>;

void walkValue(Value& value, const std::string& path, const StringFunction& fn) {
    if (std::string* text = value.get<std::string>()) {
        *text = fn(path, "", *text);
    } else if (Table* table = value.get<Table>()) {
        for (auto& [key, item] : *table) {
            walkValue(item, joinPath(path, key), fn);
        }
    } else if (Array* array = value.get<Array>()) {
        for (std::size_t i = 0; i < array->size(); ++i) {
            walkValue((*array)[i], joinPath(path, std::to_string(i)), fn);
        }
    }
}

// _walk_strings: every string of the model below `member` through fn(field, kind, text)
template <class T>
void walkStrings(T& member, const std::string& path, const std::string& modelPath, std::string_view kind,
                 const StringFunction& fn) {
    if constexpr (std::is_same_v<T, std::string>) {
        member = fn(path, kind, member);
    } else if constexpr (std::is_same_v<T, std::optional<std::string>>) {
        if (member) {
            *member = fn(path, kind, *member);
        }
    } else if constexpr (std::is_same_v<T, Table>) {
        for (auto& [key, item] : member) {
            walkValue(item, joinPath(path, key), fn);
        }
    } else if constexpr (Model<T>) {
        forEachField(member, [&](std::string_view name, auto& m) {
            const std::string inner = joinPath(modelPath, name);
            walkStrings(m, joinPath(path, name), inner, groupUse(inner), fn);
        });
    } else if constexpr (IsVector<T>::value) {
        for (std::size_t i = 0; i < member.size(); ++i) {
            walkStrings(member[i], joinPath(path, std::to_string(i)), modelPath, kind, fn);
        }
    }
}

bool isGroupTokenCharacter(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
           c == '+' || c == '-' || c == '\'';
}

// rename_group_tokens of MDATK 0.1.1: group:<old> tokens renamed, other group:<name> tokens read whole
std::string renameGroupTokens(const std::string& text, std::string_view oldName, std::string_view newName) {
    static constexpr std::string_view kToken = "group:";
    std::string out;
    std::size_t done = 0;
    std::size_t at = text.find(kToken);
    while (at != std::string::npos) {
        std::size_t end = at + kToken.size();
        while (end < text.size() && isGroupTokenCharacter(text[end])) {
            ++end;
        }
        if (end == at + kToken.size()) {  // no name: the regex tries again one character later
            at = text.find(kToken, at + 1);
            continue;
        }
        if (std::string_view(text).substr(at + kToken.size(), end - at - kToken.size()) == oldName) {
            out.append(text, done, at - done);
            out += kToken;
            out += newName;
            done = end;
        }
        at = text.find(kToken, end);
    }
    out.append(text, done, std::string::npos);
    return out;
}

bool contains(const std::vector<std::string>& names, std::string_view name) {
    return std::find(names.begin(), names.end(), name) != names.end();
}

// _renamed_text of MDATK 0.1.1
std::string renamedText(std::string_view kind, const std::string& text, std::string_view oldName,
                        std::string_view newName, const std::vector<std::string>& names) {
    if (kind == "text") {
        return text;
    }
    if (kind == "name") {
        return pythonStrip(text) == oldName ? std::string(newName) : text;
    }
    if (kind == "key") {  // base@GroupA+GroupB/output
        const std::size_t at = text.find('@');
        if (at == std::string::npos || at == 0 || text.find('/') < at) {
            return text;
        }
        const std::size_t slash = text.find('/', at + 1);
        const std::string groups = text.substr(at + 1, slash == std::string::npos ? std::string::npos : slash - at - 1);
        std::vector<std::string> items;
        std::size_t start = 0;
        for (;;) {
            const std::size_t plus = groups.find('+', start);
            items.push_back(groups.substr(start, plus == std::string::npos ? std::string::npos : plus - start));
            if (plus == std::string::npos) {
                break;
            }
            start = plus + 1;
        }
        if (!contains(items, oldName)) {
            return text;
        }
        std::string out = text.substr(0, at + 1);
        for (std::size_t i = 0; i < items.size(); ++i) {
            out += (i ? "+" : "") + (items[i] == oldName ? std::string(newName) : items[i]);
        }
        if (slash != std::string::npos) {
            out += text.substr(slash);
        }
        return out;
    }
    if (kind == "selection") {  // a group name, or comma-separated group names
        std::vector<std::string> items;
        std::size_t start = 0;
        for (;;) {
            const std::size_t comma = text.find(',', start);
            items.push_back(text.substr(start, comma == std::string::npos ? std::string::npos : comma - start));
            if (comma == std::string::npos) {
                break;
            }
            start = comma + 1;
        }
        bool named = false;
        bool allNames = true;
        for (const std::string& item : items) {
            const std::string name = pythonStrip(item);
            if (name.empty()) {
                continue;  // empty items are skipped, as split_csv does
            }
            named = named || name == oldName;
            allNames = allNames && contains(names, name);
        }
        if (named && allNames) {
            std::string out;
            for (std::size_t i = 0; i < items.size(); ++i) {
                const std::string& item = items[i];
                if (i) {
                    out += ',';
                }
                if (pythonStrip(item) == oldName) {  // keep the spaces around the name
                    const std::size_t first = item.find(oldName);
                    out += item.substr(0, first) + std::string(newName) + item.substr(first + oldName.size());
                } else {
                    out += item;
                }
            }
            return out;
        }
    }
    return renameGroupTokens(text, oldName, newName);
}

bool isPlainName(std::string_view name) {
    return !name.empty() && std::all_of(name.begin(), name.end(), [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
               c == '-';
    });
}

}  // namespace

std::string groupNameProblem(std::string_view name, const std::vector<std::string>& taken) {
    if (pythonStrip(name).empty()) {
        return "A residue group needs a name.";
    }
    if (!isPlainName(name)) {
        return "Group name " + repr(name) + " may only contain letters, digits, '_', '.' and '-' (no spaces).";
    }
    if (contains(taken, name)) {
        return "There is already a residue group named " + repr(name) + ".";
    }
    return {};
}

std::vector<std::string> renameGroup(Project& project, std::string_view oldName, std::string_view newName) {
    group(project, oldName);  // std::out_of_range for an unknown group
    const std::string name = pythonStrip(newName);
    if (name == oldName) {
        return {};
    }
    const std::vector<std::string> names = groupNames(project);
    if (std::count(names.begin(), names.end(), oldName) > 1) {
        throw ValueError("several residue groups are named " + repr(oldName) + "; give them distinct names first");
    }
    std::vector<std::string> others;
    for (const std::string& n : names) {
        if (n != oldName) {
            others.push_back(n);
        }
    }
    if (const std::string problem = groupNameProblem(name, others); !problem.empty()) {
        throw ValueError(problem);
    }
    std::vector<std::string> changed;
    walkStrings(project, std::string(), std::string(), {},
                [&](const std::string& field, std::string_view kind, const std::string& text) {
                    std::string out = renamedText(kind, text, oldName, name, names);
                    if (out != text) {
                        changed.push_back(field);
                    }
                    return out;
                });
    for (ResidueGroup& g : project.groups) {
        if (g.name == oldName) {
            g.name = name;
            break;
        }
    }
    return changed;
}

std::vector<std::string> groupReferences(const Project& project, std::string_view name) {
    std::vector<std::string> names = groupNames(project);
    if (!contains(names, name)) {
        names.emplace_back(name);
    }
    std::vector<std::string> found;
    Project copy = project;
    walkStrings(copy, std::string(), std::string(), {},
                [&](const std::string& field, std::string_view kind, const std::string& text) {
                    if (renamedText(kind, text, name, std::string_view("\0", 1), names) != text) {
                        found.push_back(field);
                    }
                    return text;
                });
    return found;
}

// ---- paths -----------------------------------------------------------------------------------

namespace {

#ifdef _WIN32
// the drive of a generic absolute path: "C:" or "//server/share" ("" for none)
std::string driveOf(std::string_view path) {
    if (path.size() >= 2 && path[1] == ':') {
        return std::string(path.substr(0, 2));
    }
    if (path.size() > 2 && path[0] == '/' && path[1] == '/') {
        const std::size_t server = path.find('/', 2);
        if (server == std::string_view::npos) {
            return std::string(path);
        }
        const std::size_t share = path.find('/', server + 1);
        return std::string(path.substr(0, share));
    }
    return {};
}
#endif

std::vector<std::string_view> components(std::string_view path) {
    std::vector<std::string_view> out;
    std::size_t start = 0;
    while (start <= path.size()) {
        std::size_t end = path.find('/', start);
        if (end == std::string_view::npos) {
            end = path.size();
        }
        if (end > start) {
            out.push_back(path.substr(start, end - start));
        }
        start = end + 1;
    }
    return out;
}

bool isRelative(std::string_view path) {
    if (path.empty() || path[0] == '~') {
        return false;
    }
#ifdef _WIN32
    return !ntpath::isabs(path) && !ntpath::isDevicePath(path);
#else
    return !posixpath::isabs(path);
#endif
}

std::string absoluteAgainst(std::string_view path, std::string_view folder) {
#ifdef _WIN32
    return ntpath::resolve(path, folder);
#else
    return posixpath::normpath(posixpath::join(folder, path));
#endif
}

std::string normalFolder(std::string_view folder) {
#ifdef _WIN32
    return ntpath::normpath(folder);
#else
    return posixpath::normpath(folder);
#endif
}

}  // namespace

std::optional<std::string> relativePath(std::string_view path, std::string_view start) {
#ifdef _WIN32
    const std::string a = ntpath::generic(ntpath::normpath(path));
    const std::string b = ntpath::generic(ntpath::normpath(start));
    const std::string da = driveOf(a);
    const std::string db = driveOf(b);
    if (asciiLower(da) != asciiLower(db)) {
        return std::nullopt;
    }
    const std::vector<std::string_view> pa = components(std::string_view(a).substr(da.size()));
    const std::vector<std::string_view> pb = components(std::string_view(b).substr(db.size()));
#else
    const std::string a = posixpath::normpath(path);
    const std::string b = posixpath::normpath(start);
    const std::vector<std::string_view> pa = components(a);
    const std::vector<std::string_view> pb = components(b);
#endif
    // ntpath.relpath compares the components with normcase (Windows names ignore case; ASCII here)
#ifdef _WIN32
    const auto same = [](std::string_view x, std::string_view y) { return asciiLower(x) == asciiLower(y); };
#else
    const auto same = [](std::string_view x, std::string_view y) { return x == y; };
#endif
    std::size_t common = 0;
    while (common < pa.size() && common < pb.size() && same(pa[common], pb[common])) {
        ++common;
    }
    std::string out;
    for (std::size_t i = common; i < pb.size(); ++i) {
        out += out.empty() ? ".." : "/..";
    }
    for (std::size_t i = common; i < pa.size(); ++i) {
        if (!out.empty()) {
            out += '/';
        }
        out += pa[i];
    }
    return out.empty() ? std::string(".") : out;
}

std::vector<std::string> rebasePaths(Project& project, std::string_view fromFolder, std::string_view toFolder) {
    std::vector<std::string> changed;
    const std::string from = normalFolder(fromFolder);
    const std::string to = normalFolder(toFolder);
    if (from == to) {
        return changed;
    }
    const auto rebased = [&](const std::string& path) {
        if (!isRelative(path)) {
            return path;
        }
        const std::string absolute = absoluteAgainst(path, from);
        return relativePath(absolute, to).value_or(absolute);
    };
    const auto fix = [&](std::string& path, const std::string& field) {
        std::string next = rebased(path);
        if (next != path) {
            path = std::move(next);
            changed.push_back(field);
        }
    };
    Project& p = project;
    fix(p.inputs.topology, "inputs.topology");
    for (std::size_t i = 0; i < p.inputs.trajectories.size(); ++i) {
        fix(p.inputs.trajectories[i], "inputs.trajectories." + std::to_string(i));
    }
    for (std::size_t i = 0; i < p.inputs.external.size(); ++i) {
        fix(p.inputs.external[i].path, "inputs.external." + std::to_string(i) + ".path");
    }
    fix(p.cleaning.reference_file, "cleaning.reference_file");
    fix(p.cleaning.output_dir, "cleaning.output_dir");
#ifdef _WIN32
    const bool cpptrajIsAPath = p.cleaning.cpptraj_path.find_first_of("/\\") != std::string::npos;
#else
    const bool cpptrajIsAPath = p.cleaning.cpptraj_path.find('/') != std::string::npos;
#endif
    if (cpptrajIsAPath) {
        fix(p.cleaning.cpptraj_path, "cleaning.cpptraj_path");
    }
    fix(p.output.out_dir, "output.out_dir");
    fix(p.output.cache_dir, "output.cache_dir");
    fix(p.density.map_path, "density.map_path");
    fix(p.workflow.out_dir, "workflow.out_dir");
    const auto nodes = [&](std::vector<PlanNode>& list, const std::string& where) {
        for (std::size_t i = 0; i < list.size(); ++i) {
            PlanNode& node = list[i];
            const std::string at = where + "." + std::to_string(i);
            for (auto& [port, value] : node.files) {
                std::string* text = value.get<std::string>();
                if (!text) {
                    continue;
                }
                std::vector<std::string> pieces;
                bool any = false;
                std::size_t start = 0;
                for (;;) {
                    const std::size_t end = text->find(';', start);
                    const std::string piece = pythonStrip(std::string_view(*text).substr(
                        start, end == std::string::npos ? std::string::npos : end - start));
                    if (!piece.empty()) {
                        std::string next = rebased(piece);
                        any = any || next != piece;
                        pieces.push_back(std::move(next));
                    }
                    if (end == std::string::npos) {
                        break;
                    }
                    start = end + 1;
                }
                if (any) {
                    std::string joined;
                    for (std::size_t k = 0; k < pieces.size(); ++k) {
                        joined += (k ? ";" : "") + pieces[k];
                    }
                    *text = joined;
                    changed.push_back(at + ".files." + port);
                }
            }
            if (const StepType* type = findStepType(node.tool)) {
                for (const std::string& key : type->pathSettings) {
                    for (auto& [name, value] : node.settings) {
                        if (name == key) {
                            std::string* text = value.get<std::string>();
                            if (text && !text->empty()) {
                                fix(*text, at + ".settings." + name);
                            }
                            break;
                        }
                    }
                }
            }
        }
    };
    nodes(p.workflow.nodes, "workflow.nodes");
    nodes(p.workflow.scratch, "workflow.scratch");
    return changed;
}

// ---- recovery copies -------------------------------------------------------------------------

namespace {

constexpr std::string_view kRecoveryOf = "recovery_of";
constexpr std::string_view kRecoveryFolder = "recovery_folder";

int processId() {
#ifdef _WIN32
    return _getpid();
#else
    return static_cast<int>(::getpid());
#endif
}

}  // namespace

fs::path recoveryDir() { return mdbridgeConfigDir() / "recovery"; }

fs::path recoveryPath(const ProjectDocument& document, const fs::path& dir) {
    if (document.hasPath()) {
        return dir / ("recovery-" + sha256Hex(pathToUtf8(document.path())).substr(0, 16) + std::string(kProjectSuffix));
    }
    return dir / ("recovery-untitled-" + std::to_string(processId()) + "-" + std::to_string(document.serial()) +
                  std::string(kProjectSuffix));
}

fs::path writeRecovery(const ProjectDocument& document, const fs::path& dir) {
    fs::create_directories(dir);
    const fs::path file = recoveryPath(document, dir);
    std::string text = tomlKey(kRecoveryOf) + " = " +
                       tomlQuote(document.hasPath() ? pathToUtf8(document.path()) : std::string()) + "\n" +
                       tomlKey(kRecoveryFolder) + " = " + tomlQuote(document.folder()) + "\n";
    text += projectText(document.project(), {}, document.comments());
    detail::replaceFile(file, text);
    return file;
}

namespace {

std::optional<RecoveryFile> readRecoveryHeader(const fs::path& file, Value* data) {
    try {
        Value parsed = parseToml(detail::readTextFile(file), detail::displayPath(file));
        const Value* of = parsed.find(kRecoveryOf);
        const Value* folder = parsed.find(kRecoveryFolder);
        const Value* format = parsed.find("format");
        if (!of || !of->isString() || !folder || !folder->isString() || !format ||
            *format != Value(std::string(kProjectFormat))) {
            return std::nullopt;
        }
        RecoveryFile out;
        out.file = file;
        out.original = of->stringOr().empty() ? fs::path() : pathFromUtf8(of->stringOr());
        out.folder = folder->stringOr();
        if (const Value* saved = parsed.find("saved"); saved && saved->get<DateTime>()) {
            out.saved = *saved->get<DateTime>();
        }
        if (data) {
            *data = std::move(parsed);
        }
        return out;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

}  // namespace

std::vector<RecoveryFile> findRecoveries(const fs::path& dir) {
    std::vector<RecoveryFile> out;
    std::error_code ec;
    for (fs::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec)) {
        const std::string name = pathToUtf8(it->path().filename());
        if (name.rfind("recovery-", 0) != 0 || name.size() < kProjectSuffix.size() ||
            name.compare(name.size() - kProjectSuffix.size(), kProjectSuffix.size(), kProjectSuffix) != 0) {
            continue;
        }
        if (auto found = readRecoveryHeader(it->path(), nullptr)) {
            out.push_back(std::move(*found));
        }
    }
    const auto key = [](const DateTime& t) {
        return std::make_tuple(t.year, t.month, t.day, t.hour, t.minute, t.second, t.microsecond);
    };
    std::sort(out.begin(), out.end(), [&](const RecoveryFile& a, const RecoveryFile& b) {
        if (key(a.saved) != key(b.saved)) {
            return key(a.saved) > key(b.saved);
        }
        return a.file < b.file;
    });
    return out;
}

void openRecovery(ProjectDocument& document, const RecoveryFile& recovery, std::vector<std::string>* warnings) {
    Value data;
    const auto header = readRecoveryHeader(recovery.file, &data);
    if (!header) {
        throw ProjectError(detail::displayPath(recovery.file) + ": not a recovery copy of a project");
    }
    Table settings;
    std::string version;
    for (auto& [key, value] : *data.get<Table>()) {
        if (key == kRecoveryOf || key == kRecoveryFolder) {
            continue;
        }
        if (key == "mdatk_version") {
            version = value.stringOr();
        }
        settings.emplace_back(key, std::move(value));
    }
    std::vector<std::string> found;
    Project project = projectFromPlain(Value(std::move(settings)), {}, &found);
    project.mdatkVersion = version;
    const std::string folder =
        header->original.empty() ? detail::absolutePath(".") : detail::absoluteDirectory(header->original);
    rebasePaths(project, header->folder, folder);
    document.assign(std::move(project), header->original, true, FileKind::Project);
    document.m_fileComments = fileComments(detail::readTextFile(recovery.file));
    if (warnings) {
        warnings->insert(warnings->end(), found.begin(), found.end());
    }
}

bool removeRecovery(const fs::path& file) {
    std::error_code ec;
    return fs::remove(file, ec) && !ec;
}

}  // namespace mdbridge::core
