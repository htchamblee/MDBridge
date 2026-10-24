#include "mdbridge/core/toml_writer.hpp"

#include <cmath>
#include <cstdio>

#include "mdbridge/core/json.hpp"

namespace mdbridge::core {

namespace {

constexpr int kWidth = 100;               // comment wrap width (tomlio._WIDTH)
constexpr std::size_t kInlineMax = 96;    // longer arrays are written one element per line

bool isBareKey(std::string_view name) {
    if (name.empty()) {
        return false;
    }
    for (const char c : name) {
        const bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
        if (!ok) {
            return false;
        }
    }
    return true;
}

void appendUnicodeEscape(std::string& out, char32_t code) {
    char escape[8];
    std::snprintf(escape, sizeof escape, "\\u%04X", static_cast<unsigned>(code));
    out += escape;
}

bool isControl(char32_t c) { return c < 0x20 || c == 0x7F; }

// tomlio._multiline: a multi-line basic string; the newline after the opening quotes is not part of it
std::string multiline(std::string_view text) {
    std::string out = "\"\"\"\n";
    for (const char32_t c : decodeUtf8(text)) {
        if (c == U'\\') {
            out += "\\\\";
        } else if (c == U'"') {
            out += "\\\"";  // never three quotes in a row
        } else if (c == U'\n' || c == U'\t') {
            out += static_cast<char>(c);
        } else if (c == U'\r') {
            out += "\\r";
        } else if (isControl(c)) {
            appendUnicodeEscape(out, c);
        } else {
            appendUtf8(out, c);
        }
    }
    return out + "\"\"\"";
}

std::string floatText(double v) {
    if (std::isnan(v)) {
        return "nan";
    }
    if (std::isinf(v)) {
        return v > 0 ? "inf" : "-inf";
    }
    return pythonFloatRepr(v);  // has ".0", "." or an exponent, as tomlio requires
}

bool isTableArray(const Value& v) {
    const Array* array = v.get<Array>();
    if (!array || array->empty()) {
        return false;
    }
    for (const Value& item : *array) {
        if (!item.isTable()) {
            return false;
        }
    }
    return true;
}

std::string joinPath(const std::vector<std::string>& path) {
    std::string out;
    for (std::size_t i = 0; i < path.size(); ++i) {
        out += (i ? "." : "") + tomlKey(path[i]);
    }
    return out;
}

}  // namespace

std::string tomlKey(std::string_view name) { return isBareKey(name) ? std::string(name) : tomlQuote(name); }

std::string tomlQuote(std::string_view text) {
    std::string out = "\"";
    for (const char32_t c : decodeUtf8(text)) {
        switch (c) {
        case U'\\':
            out += "\\\\";
            break;
        case U'"':
            out += "\\\"";
            break;
        case U'\n':
            out += "\\n";
            break;
        case U'\t':
            out += "\\t";
            break;
        case U'\r':
            out += "\\r";
            break;
        default:
            if (isControl(c)) {
                appendUnicodeEscape(out, c);
            } else {
                appendUtf8(out, c);
            }
        }
    }
    return out + "\"";
}

std::string tomlValue(const Value& value) {
    switch (value.type()) {
    case Value::Type::Bool:
        return *value.get<bool>() ? "true" : "false";
    case Value::Type::Int:
        return std::to_string(*value.get<std::int64_t>());
    case Value::Type::Float:
        return floatText(*value.get<double>());
    case Value::Type::String: {
        const std::string& s = *value.get<std::string>();
        return s.find('\n') != std::string::npos ? multiline(s) : tomlQuote(s);
    }
    case Value::Type::DateTime:
        return value.get<DateTime>()->isoformat();
    case Value::Type::Table: {
        const Table& table = *value.get<Table>();
        if (table.empty()) {
            return "{}";
        }
        std::string items;
        for (const auto& [k, x] : table) {
            if (!x.isNull()) {
                items += (items.empty() ? "" : ", ") + tomlKey(k) + " = " + tomlValue(x);
            }
        }
        return "{ " + items + " }";
    }
    case Value::Type::Array: {
        std::vector<std::string> items;
        for (const Value& x : *value.get<Array>()) {
            if (!x.isNull()) {
                items.push_back(tomlValue(x));
            }
        }
        std::string flat = "[";
        for (std::size_t i = 0; i < items.size(); ++i) {
            flat += (i ? ", " : "") + items[i];
        }
        flat += "]";
        if (characterCount(flat) <= kInlineMax && flat.find('\n') == std::string::npos) {
            return flat;
        }
        std::string out = "[\n";
        for (const std::string& item : items) {
            out += "    " + item + ",\n";
        }
        return out + "]";
    }
    case Value::Type::Null:
        break;
    }
    return "";  // None has no TOML form: the writer comments such entries out
}

TomlWriter::TomlWriter(TomlDocFunction doc) : m_doc(std::move(doc)) {}

std::optional<std::string> TomlWriter::doc(const std::vector<std::string>& path,
                                           const std::vector<const Value*>& stack) const {
    return m_doc ? m_doc(path, stack) : std::nullopt;
}

void TomlWriter::comment(const std::optional<std::string>& text, std::string_view indent) {
    if (!text || text->empty()) {
        return;
    }
    const int width = kWidth - static_cast<int>(characterCount(indent)) - 2;
    std::vector<std::string> paragraphs = pythonSplitlines(*text);
    if (paragraphs.empty()) {
        paragraphs.emplace_back();
    }
    for (const std::string& paragraph : paragraphs) {
        std::vector<std::string> wrapped = pythonTextwrap(paragraph, width);
        if (wrapped.empty()) {
            wrapped.emplace_back();
        }
        for (const std::string& line : wrapped) {
            m_lines.push_back(pythonRstrip(std::string(indent) + "# " + line));
        }
    }
}

void TomlWriter::blank() {
    if (!m_lines.empty() && !m_lines.back().empty()) {
        m_lines.emplace_back();
    }
}

std::string TomlWriter::document(const Value& data, const std::vector<std::string>& header) {
    for (const std::string& line : header) {
        m_lines.push_back(line.empty() ? "#" : pythonRstrip("# " + line));
    }
    if (!header.empty()) {
        m_lines.emplace_back();
    }
    body({}, data, {&data});
    std::string out;
    for (std::size_t i = 0; i < m_lines.size(); ++i) {
        out += (i ? "\n" : "") + m_lines[i];
    }
    while (!out.empty() && out.back() == '\n') {
        out.pop_back();
    }
    return out + "\n";
}

void TomlWriter::table(const std::vector<std::string>& path, const Value& data, std::vector<const Value*> stack,
                       bool array, bool first) {
    blank();
    if (first) {
        comment(doc(path, stack));
    }
    const std::string name = joinPath(path);
    m_lines.push_back(array ? "[[" + name + "]]" : "[" + name + "]");
    stack.push_back(&data);
    body(path, data, stack);
}

void TomlWriter::body(const std::vector<std::string>& path, const Value& data, const std::vector<const Value*>& stack) {
    const Table* table = data.get<Table>();
    if (!table) {
        return;
    }
    // scalars (and arrays of scalars) and empty tables first, then tables, then arrays of tables
    const auto here = [&path](const std::string& key) {
        std::vector<std::string> p = path;
        p.push_back(key);
        return p;
    };
    const auto writeEntry = [&](const std::string& key, const Value& v) {
        const std::vector<std::string> p = here(key);
        comment(doc(p, stack));
        if (v.isNull()) {
            m_lines.push_back("# " + tomlKey(key) + " =   (not set)");
        } else {
            m_lines.push_back(tomlKey(key) + " = " + tomlValue(v));
        }
    };
    for (const auto& [key, v] : *table) {
        if (!v.isTable() && !isTableArray(v)) {
            writeEntry(key, v);
        }
    }
    for (const auto& [key, v] : *table) {
        if (v.isTable() && v.size() == 0) {
            writeEntry(key, v);
        }
    }
    for (const auto& [key, v] : *table) {
        if (v.isTable() && v.size() > 0) {
            this->table(here(key), v, stack);
        }
    }
    for (const auto& [key, v] : *table) {
        if (isTableArray(v)) {
            const Array& items = *v.get<Array>();
            for (std::size_t i = 0; i < items.size(); ++i) {
                this->table(here(key), items[i], stack, true, i == 0);
            }
        }
    }
}

std::string dumpToml(const Value& data, TomlDocFunction doc, const std::vector<std::string>& header) {
    return TomlWriter(std::move(doc)).document(data, header);
}

}  // namespace mdbridge::core
