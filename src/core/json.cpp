#include "mdbridge/core/json.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "mdbridge/core/text.hpp"

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


void appendHex(std::string& out, char32_t code, int digits) {
    static constexpr char kHex[] = "0123456789abcdef";
    for (int shift = 4 * (digits - 1); shift >= 0; shift -= 4) {
        out += kHex[(code >> shift) & 0xF];
    }
}

// Deeper nesting is refused before the stack (or the destructors of the Value tree) can run out.
constexpr int kMaxJsonDepth = 1000;

bool isJsonSpace(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }
bool isAsciiDigit(char c) { return c >= '0' && c <= '9'; }

// json.loads as CPython's C scanner (Modules/_json.c) does it, on the UTF-8 bytes: every
// character the grammar looks at is ASCII, so byte offsets work wherever Python counts characters,
// and only the positions in messages are converted to characters.
class JsonParser {
public:
    JsonParser(std::string_view text, std::string_view source) : m_text(text), m_source(source) {}

    Value parse() {
        if (m_text.substr(0, 3) == "\xEF\xBB\xBF") {
            fail("Unexpected UTF-8 BOM (decode using utf-8-sig)", 0);
        }
        std::size_t i = skipSpace(0);
        BigInt big;
        Value value = scan(i, 0, big);
        i = skipSpace(i);
        if (i != m_text.size()) {
            fail("Extra data", i);
        }
        if (big) {
            fail("Integer outside the 64-bit range", *big);
        }
        return value;
    }

private:
    // Where the first integer beyond 64 bits starts that the result would hold (Python holds a
    // big int). It is reported only once the whole text has been read, so that json's own errors
    // come first and a value that a repeated key replaces does not count, as in Python.
    using BigInt = std::optional<std::size_t>;

    static void keepFirst(BigInt& into, const BigInt& from) {
        if (from && (!into || *from < *into)) {
            into = from;
        }
    }

    // json.JSONDecodeError: "<message>: line L column C (char P)", P counted in characters
    [[noreturn]] void fail(std::string_view message, std::size_t at) const {
        const std::string_view before = m_text.substr(0, at);
        const std::size_t newline = before.rfind('\n');
        const std::size_t pos = characterCount(before);
        const auto lines = std::count(before.begin(), before.end(), '\n');
        const std::size_t line = 1 + static_cast<std::size_t>(lines);
        const std::size_t column = newline == std::string_view::npos
                                       ? pos + 1
                                       : pos - characterCount(before.substr(0, newline));
        throw ProjectError(std::string(m_source) + ": not valid JSON: " + std::string(message) +
                           ": line " + std::to_string(line) + " column " + std::to_string(column) +
                           " (char " + std::to_string(pos) + ")");
    }

    std::size_t skipSpace(std::size_t i) const {
        while (i < m_text.size() && isJsonSpace(m_text[i])) {
            ++i;
        }
        return i;
    }

    bool follows(std::size_t i, std::string_view word) const {
        return m_text.substr(i, word.size()) == word;
    }

    // scan_once: the value at i; i moves past it
    Value scan(std::size_t& i, int depth, BigInt& big) {
        const std::size_t n = m_text.size();
        if (i >= n) {
            fail("Expecting value", i);
        }
        switch (m_text[i]) {
        case '"':
            ++i;
            return scanString(i);
        case '{':
        case '[': {
            const bool object = m_text[i] == '{';
            if (depth >= kMaxJsonDepth) {  // Python's RecursionError, which has no position
                throw ProjectError(std::string(m_source) + ": not valid JSON: maximum recursion " +
                                   "depth exceeded while decoding a JSON " +
                                   (object ? "object" : "array") + " from a unicode string");
            }
            ++i;
            return object ? scanObject(i, depth + 1, big) : scanArray(i, depth + 1, big);
        }
        case 'n':
            if (follows(i, "null")) {
                i += 4;
                return {};
            }
            break;
        case 't':
            if (follows(i, "true")) {
                i += 4;
                return true;
            }
            break;
        case 'f':
            if (follows(i, "false")) {
                i += 5;
                return false;
            }
            break;
        case 'N':
            if (follows(i, "NaN")) {
                i += 3;
                return std::numeric_limits<double>::quiet_NaN();
            }
            break;
        case 'I':
            if (follows(i, "Infinity")) {
                i += 8;
                return std::numeric_limits<double>::infinity();
            }
            break;
        case '-':
            if (follows(i, "-Infinity")) {
                i += 9;
                return -std::numeric_limits<double>::infinity();
            }
            break;
        default:
            break;
        }
        return scanNumber(i, big);
    }

    // i: the character after '{'
    Value scanObject(std::size_t& i, int depth, BigInt& big) {
        const std::size_t n = m_text.size();
        Table table;
        std::vector<BigInt> bigs;                             // of each entry
        std::unordered_map<std::string, std::size_t> index;  // a repeated key keeps its first place
        i = skipSpace(i);
        if (i >= n || m_text[i] != '}') {
            while (true) {
                if (i >= n || m_text[i] != '"') {
                    fail("Expecting property name enclosed in double quotes", i);
                }
                ++i;
                std::string key = scanString(i);
                i = skipSpace(i);
                if (i >= n || m_text[i] != ':') {
                    fail("Expecting ':' delimiter", i);
                }
                i = skipSpace(i + 1);
                BigInt valueBig;
                Value value = scan(i, depth, valueBig);
                if (const auto found = index.find(key); found != index.end()) {
                    table[found->second].second = std::move(value);
                    bigs[found->second] = valueBig;
                } else {
                    index.emplace(key, table.size());
                    table.emplace_back(std::move(key), std::move(value));
                    bigs.push_back(valueBig);
                }
                i = skipSpace(i);
                if (i < n && m_text[i] == '}') {
                    break;
                }
                if (i >= n || m_text[i] != ',') {
                    fail("Expecting ',' delimiter", i);
                }
                const std::size_t comma = i;
                i = skipSpace(i + 1);
                if (i < n && m_text[i] == '}') {
                    fail("Illegal trailing comma before end of object", comma);
                }
            }
        }
        ++i;
        for (const BigInt& entry : bigs) {
            keepFirst(big, entry);
        }
        return table;
    }

    // i: the character after '['
    Value scanArray(std::size_t& i, int depth, BigInt& big) {
        const std::size_t n = m_text.size();
        Array array;
        i = skipSpace(i);
        if (i >= n || m_text[i] != ']') {
            while (true) {
                BigInt itemBig;
                array.push_back(scan(i, depth, itemBig));
                keepFirst(big, itemBig);
                i = skipSpace(i);
                if (i < n && m_text[i] == ']') {
                    break;
                }
                if (i >= n || m_text[i] != ',') {
                    fail("Expecting ',' delimiter", i);
                }
                const std::size_t comma = i;
                i = skipSpace(i + 1);
                if (i < n && m_text[i] == ']') {
                    fail("Illegal trailing comma before end of array", comma);
                }
            }
        }
        ++i;
        return array;
    }

    // Four hex digits at i, or -1.
    long hex4(std::size_t i) const {
        long code = 0;
        for (std::size_t k = i; k < i + 4; ++k) {
            const char c = m_text[k];
            int digit = 0;
            if (c >= '0' && c <= '9') {
                digit = c - '0';
            } else if (c >= 'a' && c <= 'f') {
                digit = c - 'a' + 10;
            } else if (c >= 'A' && c <= 'F') {
                digit = c - 'A' + 10;
            } else {
                return -1;
            }
            code = code * 16 + digit;
        }
        return code;
    }

    // scanstring (strict): `end` is the character after the opening quote; it moves past the
    // closing one
    std::string scanString(std::size_t& end) {
        const std::size_t n = m_text.size();
        const std::size_t begin = end - 1;
        std::string out;
        while (true) {
            std::size_t next = end;
            while (next < n && m_text[next] != '"' && m_text[next] != '\\') {
                if (static_cast<unsigned char>(m_text[next]) <= 0x1F) {
                    fail("Invalid control character at", next);
                }
                ++next;
            }
            if (next >= n) {
                fail("Unterminated string starting at", begin);
            }
            out.append(m_text.substr(end, next - end));  // bytes as they are, valid UTF-8 or not
            if (m_text[next] == '"') {
                end = next + 1;
                return out;
            }
            ++next;  // past the backslash
            if (next == n) {
                fail("Unterminated string starting at", begin);
            }
            if (m_text[next] != 'u') {
                end = next + 1;
                switch (m_text[next]) {
                case '"':
                case '\\':
                case '/':
                    out += m_text[next];
                    break;
                case 'b':
                    out += '\b';
                    break;
                case 'f':
                    out += '\f';
                    break;
                case 'n':
                    out += '\n';
                    break;
                case 'r':
                    out += '\r';
                    break;
                case 't':
                    out += '\t';
                    break;
                default:
                    fail("Invalid \\escape", end - 2);
                }
                continue;
            }
            ++next;
            end = next + 4;
            const long first = end > n ? -1 : hex4(next);
            if (first < 0) {
                fail("Invalid \\uXXXX escape", next - 1);
            }
            auto code = static_cast<char32_t>(first);
            // a high surrogate joins a following \u low surrogate (only if more text follows it)
            if (code >= 0xD800 && code <= 0xDBFF && end + 6 < n && m_text[end] == '\\' &&
                m_text[end + 1] == 'u') {
                const long second = hex4(end + 2);
                if (second < 0) {
                    fail("Invalid \\uXXXX escape", end + 1);
                }
                if (second >= 0xDC00 && second <= 0xDFFF) {
                    code = 0x10000 + ((code - 0xD800) << 10) +
                           static_cast<char32_t>(second - 0xDC00);
                    end += 6;
                }
            }
            // Python keeps a lone surrogate; UTF-8 text cannot hold one
            appendUtf8(out, code >= 0xD800 && code <= 0xDFFF ? char32_t{0xFFFD} : code);
        }
    }

    // _match_number: -?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?, an int unless it has a
    // fraction or an exponent
    Value scanNumber(std::size_t& i, BigInt& big) {
        const std::size_t n = m_text.size();
        const std::size_t start = i;
        std::size_t idx = i;
        if (m_text[idx] == '-' && ++idx >= n) {
            fail("Expecting value", start);
        }
        if (m_text[idx] >= '1' && m_text[idx] <= '9') {
            while (++idx < n && isAsciiDigit(m_text[idx])) {
            }
        } else if (m_text[idx] == '0') {
            ++idx;
        } else {
            fail("Expecting value", start);
        }
        bool isFloat = false;
        if (idx + 1 < n && m_text[idx] == '.' && isAsciiDigit(m_text[idx + 1])) {
            isFloat = true;
            idx += 2;
            while (idx < n && isAsciiDigit(m_text[idx])) {
                ++idx;
            }
        }
        if (idx + 1 < n && (m_text[idx] == 'e' || m_text[idx] == 'E')) {
            const std::size_t exponentStart = idx++;
            if (idx + 1 < n && (m_text[idx] == '-' || m_text[idx] == '+')) {
                ++idx;
            }
            while (idx < n && isAsciiDigit(m_text[idx])) {
                ++idx;
            }
            if (isAsciiDigit(m_text[idx - 1])) {
                isFloat = true;
            } else {
                idx = exponentStart;  // "1e" or "1e+": the number ends before the e
            }
        }
        const std::string_view number = m_text.substr(start, idx - start);
        i = idx;
        if (isFloat) {
            return parseFloat(number);
        }
        std::int64_t value = 0;
        const auto result = std::from_chars(number.data(), number.data() + number.size(), value);
        if (result.ec != std::errc()) {
            big = start;
            return {};
        }
        return value;
    }

    // float(text) for a number of the JSON grammar: correctly rounded; too large -> inf, too
    // small -> 0, keeping the sign
    static double parseFloat(std::string_view number) {
        double value = 0;
        const auto result = std::from_chars(number.data(), number.data() + number.size(), value);
        if (result.ec != std::errc::result_out_of_range) {
            return value;
        }
        const bool negative = number.front() == '-';
        // decimal exponent of the first significant digit
        const std::size_t e = number.find_first_of("eE");
        const std::string_view mantissa = number.substr(negative ? 1 : 0, e - (negative ? 1 : 0));
        const std::size_t point = std::min(mantissa.find('.'), mantissa.size());
        const std::size_t first = mantissa.find_first_not_of("0.");
        if (first == std::string_view::npos) {
            return negative ? -0.0 : 0.0;
        }
        long long magnitude = first < point ? static_cast<long long>(point - first) - 1
                                            : -static_cast<long long>(first - point);
        if (e != std::string_view::npos) {
            std::string_view digits = number.substr(e + 1);
            const bool down = digits.front() == '-';
            if (digits.front() == '-' || digits.front() == '+') {
                digits.remove_prefix(1);
            }
            long long exponent = 0;
            for (const char c : digits) {
                exponent = std::min(exponent * 10 + (c - '0'), 1'000'000'000LL);
            }
            magnitude += down ? -exponent : exponent;
        }
        const double out = magnitude > 0 ? std::numeric_limits<double>::infinity() : 0.0;
        return negative ? -out : out;
    }

    std::string_view m_text;
    std::string_view m_source;
};

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
    // unicode_repr of CPython
    const std::u32string s = decodeUtf8Surrogateescape(text);
    const bool hasSingle = s.find(U'\'') != std::u32string::npos;
    const bool hasDouble = s.find(U'"') != std::u32string::npos;
    const char quote = hasSingle && !hasDouble ? '"' : '\'';
    std::string out(1, quote);
    for (const char32_t c : s) {
        if (c == static_cast<char32_t>(quote) || c == U'\\') {
            out += '\\';
            out += static_cast<char>(c);
        } else if (c == U'\t') {
            out += "\\t";
        } else if (c == U'\n') {
            out += "\\n";
        } else if (c == U'\r') {
            out += "\\r";
        } else if (c < 0x20 || c == 0x7F) {
            out += "\\x";
            appendHex(out, c, 2);
        } else if (c < 0x7F || isPythonPrintable(c)) {
            appendUtf8(out, c);
        } else if (c <= 0xFF) {
            out += "\\x";
            appendHex(out, c, 2);
        } else if (c <= 0xFFFF) {
            out += "\\u";
            appendHex(out, c, 4);
        } else {
            out += "\\U";
            appendHex(out, c, 8);
        }
    }
    out += quote;
    return out;
}

Value parseJson(std::string_view text, std::string_view source) {
    return JsonParser(text, source).parse();
}

}  // namespace mdbridge::core
