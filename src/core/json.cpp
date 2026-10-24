#include "mdbridge/core/json.hpp"

#include <charconv>
#include <cmath>
#include <cstdio>
#include <string_view>

namespace mdbridge::core {

namespace {

// json's encode_basestring (ensure_ascii=False): only quotes, backslashes and control characters
void appendJsonString(std::string& out, std::string_view text) {
    out += '"';
    for (const char ch : text) {
        switch (ch) {
        case '"':
            out += "\\\"";
            break;
        case '\\':
            out += "\\\\";
            break;
        case '\b':
            out += "\\b";
            break;
        case '\f':
            out += "\\f";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (static_cast<unsigned char>(ch) < 0x20) {
                char escape[8];
                std::snprintf(escape, sizeof escape, "\\u%04x", static_cast<unsigned>(ch));
                out += escape;
            } else {
                out += ch;
            }
        }
    }
    out += '"';
}

std::string jsonFloat(double value) {
    if (std::isnan(value)) {
        return "NaN";
    }
    if (std::isinf(value)) {
        return value > 0 ? "Infinity" : "-Infinity";
    }
    return pythonFloatRepr(value);
}

void newline(std::string& out, int indent, int depth) {
    if (indent >= 0) {
        out += '\n';
        out.append(static_cast<std::size_t>(indent) * static_cast<std::size_t>(depth), ' ');
    }
}

void appendJson(std::string& out, const Value& value, int indent, int depth) {
    const char* separator = indent >= 0 ? "," : ", ";
    switch (value.type()) {
    case Value::Type::Null:
        out += "null";
        return;
    case Value::Type::Bool:
        out += *value.get<bool>() ? "true" : "false";
        return;
    case Value::Type::Int:
        out += std::to_string(*value.get<std::int64_t>());
        return;
    case Value::Type::Float:
        out += jsonFloat(*value.get<double>());
        return;
    case Value::Type::String:
        appendJsonString(out, *value.get<std::string>());
        return;
    case Value::Type::DateTime:
        appendJsonString(out, value.get<DateTime>()->isoformat());
        return;
    case Value::Type::Array: {
        const Array& array = *value.get<Array>();
        if (array.empty()) {
            out += "[]";
            return;
        }
        out += '[';
        for (std::size_t i = 0; i < array.size(); ++i) {
            if (i) {
                out += separator;
            }
            newline(out, indent, depth + 1);
            appendJson(out, array[i], indent, depth + 1);
        }
        newline(out, indent, depth);
        out += ']';
        return;
    }
    case Value::Type::Table: {
        const Table& table = *value.get<Table>();
        if (table.empty()) {
            out += "{}";
            return;
        }
        out += '{';
        for (std::size_t i = 0; i < table.size(); ++i) {
            if (i) {
                out += separator;
            }
            newline(out, indent, depth + 1);
            appendJsonString(out, table[i].first);
            out += ": ";
            appendJson(out, table[i].second, indent, depth + 1);
        }
        newline(out, indent, depth);
        out += '}';
        return;
    }
    }
}

}  // namespace

std::string toJson(const Value& value, int indent) {
    std::string out;
    appendJson(out, value, indent, 0);
    return out;
}

std::string pythonFloatRepr(double value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    // shortest round-trip digits, read off the scientific form "[-]d[.ddd]e±XX"
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof buffer, value, std::chars_format::scientific);
    std::string_view text(buffer, static_cast<std::size_t>(result.ptr - buffer));
    std::string out;
    if (text.front() == '-') {
        out += '-';
        text.remove_prefix(1);
    }
    const std::size_t e = text.find('e');
    std::string digits(1, text[0]);
    if (e > 1) {
        digits += text.substr(2, e - 2);
    }
    int exponent = 0;
    std::string_view expText = text.substr(e + 1);
    if (expText.front() == '+') {
        expText.remove_prefix(1);
    }
    std::from_chars(expText.data(), expText.data() + expText.size(), exponent);

    // CPython's float_repr_style "short": exponent notation when the decimal point would sit
    // more than 4 places left of the digits or more than 16 places right of the first digit
    const int point = exponent + 1;  // position of the decimal point after the first digit
    const int n = static_cast<int>(digits.size());
    if (point <= -4 || point > 16) {
        out += digits[0];
        if (n > 1) {
            out += '.';
            out += digits.substr(1);
        }
        char exp[16];
        std::snprintf(exp, sizeof exp, "e%c%02d", exponent < 0 ? '-' : '+', std::abs(exponent));
        out += exp;
    } else if (point <= 0) {
        out += "0.";
        out.append(static_cast<std::size_t>(-point), '0');
        out += digits;
    } else if (point >= n) {
        out += digits;
        out.append(static_cast<std::size_t>(point - n), '0');
        out += ".0";
    } else {
        out += digits.substr(0, static_cast<std::size_t>(point));
        out += '.';
        out += digits.substr(static_cast<std::size_t>(point));
    }
    return out;
}

std::string pythonStringRepr(const std::string& text) {
    const bool hasSingle = text.find('\'') != std::string::npos;
    const bool hasDouble = text.find('"') != std::string::npos;
    const char quote = hasSingle && !hasDouble ? '"' : '\'';
    std::string out(1, quote);
    for (const char ch : text) {
        const auto code = static_cast<unsigned char>(ch);
        if (ch == quote || ch == '\\') {
            out += '\\';
            out += ch;
        } else if (ch == '\n') {
            out += "\\n";
        } else if (ch == '\r') {
            out += "\\r";
        } else if (ch == '\t') {
            out += "\\t";
        } else if (code < 0x20 || code == 0x7f) {
            char escape[8];
            std::snprintf(escape, sizeof escape, "\\x%02x", code);
            out += escape;
        } else {
            out += ch;
        }
    }
    out += quote;
    return out;
}

}  // namespace mdbridge::core
