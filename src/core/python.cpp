#include "mdbridge/core/python.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#if !defined(_WIN32)
#include <cerrno>
#include <pwd.h>
#include <unistd.h>
#endif

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

void appendInts(std::string& out, const std::vector<int>& fields) {
    for (std::size_t i = 0; i < fields.size(); ++i) {
        out += (i ? ", " : "") + std::to_string(fields[i]);
    }
}

// datetime.timedelta.__repr__ of a whole number of seconds: days=, seconds= (normalised as
// timedelta does: -05:30 is days=-1, seconds=66600)
std::string timedeltaRepr(long long totalSeconds) {
    long long days = totalSeconds / 86400;
    long long seconds = totalSeconds % 86400;
    if (seconds < 0) {
        seconds += 86400;
        --days;
    }
    std::string args;
    if (days) {
        args = "days=" + std::to_string(days);
    }
    if (seconds) {
        args += (args.empty() ? "" : ", ") + std::string("seconds=") + std::to_string(seconds);
    }
    return "datetime.timedelta(" + (args.empty() ? std::string("0") : args) + ")";
}

// datetime.date / time / datetime __repr__: trailing zero microseconds and then trailing zero
// seconds are left out (hour and minute stay)
std::string dateTimeRepr(const DateTime& dt) {
    if (dt.kind == DateTime::Kind::LocalDate) {
        std::string out = "datetime.date(";
        appendInts(out, {dt.year, dt.month, dt.day});
        return out + ")";
    }
    std::vector<int> fields;
    if (dt.kind != DateTime::Kind::LocalTime) {
        fields = {dt.year, dt.month, dt.day};
    }
    fields.insert(fields.end(), {dt.hour, dt.minute, dt.second, dt.microsecond});
    for (int k = 0; k < 2 && fields.back() == 0; ++k) {
        fields.pop_back();
    }
    const bool time = dt.kind == DateTime::Kind::LocalTime;
    std::string out = time ? "datetime.time(" : "datetime.datetime(";
    appendInts(out, fields);
    if (dt.kind == DateTime::Kind::OffsetDateTime) {
        // tomllib gives timezone.utc for a zero offset (Z, +00:00 and -00:00)
        out += ", tzinfo=";
        out += dt.offsetMinutes == 0 ? std::string("datetime.timezone.utc")
                                     : "datetime.timezone(" +
                                           timedeltaRepr(60LL * dt.offsetMinutes) + ")";
    }
    return out + ")";
}

void appendRepr(std::string& out, const Value& value) {
    switch (value.type()) {
    case Value::Type::Null:
        out += "None";
        return;
    case Value::Type::Bool:
        out += *value.get<bool>() ? "True" : "False";
        return;
    case Value::Type::Int:
        out += std::to_string(*value.get<std::int64_t>());
        return;
    case Value::Type::Float:
        out += pythonFloatRepr(*value.get<double>());
        return;
    case Value::Type::String:
        out += pythonStringRepr(*value.get<std::string>());
        return;
    case Value::Type::DateTime:
        out += dateTimeRepr(*value.get<DateTime>());
        return;
    case Value::Type::Array: {
        out += '[';
        const Array& array = *value.get<Array>();
        for (std::size_t i = 0; i < array.size(); ++i) {
            if (i) {
                out += ", ";
            }
            appendRepr(out, array[i]);
        }
        out += ']';
        return;
    }
    case Value::Type::Table: {
        out += '{';
        const Table& table = *value.get<Table>();
        for (std::size_t i = 0; i < table.size(); ++i) {
            if (i) {
                out += ", ";
            }
            out += pythonStringRepr(table[i].first);
            out += ": ";
            appendRepr(out, table[i].second);
        }
        out += '}';
        return;
    }
    }
}

// \uXXXX of a UTF-16 code unit, lower-case hex digits
void appendUnicodeEscape(std::string& out, char32_t unit) {
    static constexpr char kHex[] = "0123456789abcdef";
    out += "\\u";
    for (int shift = 12; shift >= 0; shift -= 4) {
        out += kHex[(unit >> shift) & 0xF];
    }
}

// json's encode_basestring_ascii: short escapes, everything outside ' '..'~' as \uXXXX
void appendAsciiJsonString(std::string& out, std::string_view text) {
    out += '"';
    for (const char32_t c : decodeUtf8Surrogateescape(text)) {
        switch (c) {
        case U'"':
            out += "\\\"";
            break;
        case U'\\':
            out += "\\\\";
            break;
        case U'\n':
            out += "\\n";
            break;
        case U'\r':
            out += "\\r";
            break;
        case U'\t':
            out += "\\t";
            break;
        case U'\b':
            out += "\\b";
            break;
        case U'\f':
            out += "\\f";
            break;
        default:
            if (c >= 0x20 && c <= 0x7E) {
                out += static_cast<char>(c);
            } else if (c >= 0x10000) {  // a UTF-16 surrogate pair
                appendUnicodeEscape(out, 0xD800 | ((c - 0x10000) >> 10));
                appendUnicodeEscape(out, 0xDC00 | ((c - 0x10000) & 0x3FF));
            } else {
                appendUnicodeEscape(out, c);
            }
        }
    }
    out += '"';
}

void appendJsonCompact(std::string& out, const Value& value) {
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
    case Value::Type::Float: {
        const double x = *value.get<double>();
        if (std::isnan(x)) {
            out += "NaN";
        } else if (std::isinf(x)) {
            out += x > 0 ? "Infinity" : "-Infinity";
        } else {
            out += pythonFloatRepr(x);
        }
        return;
    }
    case Value::Type::String:
        appendAsciiJsonString(out, *value.get<std::string>());
        return;
    case Value::Type::DateTime:  // default=str
        appendAsciiJsonString(out, pythonStr(value));
        return;
    case Value::Type::Array: {
        out += '[';
        const Array& array = *value.get<Array>();
        for (std::size_t i = 0; i < array.size(); ++i) {
            if (i) {
                out += ',';
            }
            appendJsonCompact(out, array[i]);
        }
        out += ']';
        return;
    }
    case Value::Type::Table: {
        // sort_keys: Python compares str by code point (not by UTF-8 byte once lone surrogates
        // of undecodable bytes take part)
        const Table& table = *value.get<Table>();
        std::vector<std::pair<std::u32string, const std::pair<std::string, Value>*>> entries;
        entries.reserve(table.size());
        for (const auto& entry : table) {
            entries.emplace_back(decodeUtf8Surrogateescape(entry.first), &entry);
        }
        std::stable_sort(entries.begin(), entries.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        out += '{';
        for (std::size_t i = 0; i < entries.size(); ++i) {
            if (i) {
                out += ',';
            }
            appendAsciiJsonString(out, entries[i].second->first);
            out += ':';
            appendJsonCompact(out, entries[i].second->second);
        }
        out += '}';
        return;
    }
    }
}

// round(x, digits) for digits < 0: half-even at 10^-digits on the exact value, read off the exact
// decimal digits of the integer part (the fraction only decides ties)
std::string roundToTens(double x, int digits) {
    const double whole = std::trunc(std::fabs(x));
    char buffer[400];  // the digits of an integral double: at most 309
    const auto result =
        std::to_chars(buffer, buffer + sizeof buffer, whole, std::chars_format::fixed, 0);
    const std::string integer(buffer, result.ptr);
    const bool fraction = whole != std::fabs(x);
    const auto drop = static_cast<std::size_t>(-digits);
    std::string kept = integer.size() > drop ? integer.substr(0, integer.size() - drop) : "0";
    const std::string rest =
        integer.size() > drop ? integer.substr(integer.size() - drop)
                              : std::string(drop - integer.size(), '0') + integer;
    bool up = rest[0] > '5';
    if (rest[0] == '5') {
        const bool above = fraction || rest.find_first_not_of('0', 1) != std::string::npos;
        up = above || (kept.back() - '0') % 2 == 1;
    }
    if (up) {
        std::size_t k = kept.size();
        while (k > 0 && kept[k - 1] == '9') {
            kept[--k] = '0';
        }
        if (k == 0) {
            kept.insert(kept.begin(), '1');
        } else {
            ++kept[k - 1];
        }
    }
    return (std::signbit(x) ? "-" : "") + kept + "e" + std::to_string(drop);
}

#if defined(_WIN32)
// MDBridge's Windows rule (ntpath.hpp; MDATK's Windows build uses Python's ntpath, which CI cannot
// run): "~" is the profile folder, as in Python's ntpath; "~user" stays as it is (ntpath guesses a
// sibling folder of the profile).
std::optional<std::string> currentUserHome() {
    if (const wchar_t* profile = _wgetenv(L"USERPROFILE")) {
        return pathToUtf8(std::filesystem::path(profile));
    }
    return std::nullopt;
}

std::optional<std::string> userHome(const std::string&) { return std::nullopt; }
#else
// pw_dir of a passwd entry (by name, or of the process's user)
std::optional<std::string> passwdHome(const std::string* name) {
    const long hint = ::sysconf(_SC_GETPW_R_SIZE_MAX);
    std::vector<char> buffer(hint > 0 ? static_cast<std::size_t>(hint) : 16384);
    passwd entry{};
    passwd* found = nullptr;
    while (true) {
        const int rc =
            name ? ::getpwnam_r(name->c_str(), &entry, buffer.data(), buffer.size(), &found)
                 : ::getpwuid_r(::getuid(), &entry, buffer.data(), buffer.size(), &found);
        if (rc != ERANGE || buffer.size() > (std::size_t{1} << 24)) {
            break;
        }
        buffer.resize(buffer.size() * 2);
    }
    if (!found || !found->pw_dir) {
        return std::nullopt;
    }
    return std::string(found->pw_dir);
}

std::optional<std::string> currentUserHome() {
    if (const char* home = std::getenv("HOME")) {  // set counts, even when empty
        return std::string(home);
    }
    return passwdHome(nullptr);
}

std::optional<std::string> userHome(const std::string& name) {
    if (name.find('\0') != std::string::npos) {  // no such user (Python raises ValueError)
        return std::nullopt;
    }
    return passwdHome(&name);
}
#endif

bool isAsciiSpace(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
}
bool isAsciiDigit(char c) { return c >= '0' && c <= '9'; }

int decimalValue(char32_t code) {
    // Unicode decimal digits come in runs of ten that start at a zero
    char32_t first = code;
    while (first > 0 && isPythonDigit(first - 1)) {
        --first;
    }
    return static_cast<int>((code - first) % 10);
}

// The string _PyUnicode_TransformDecimalAndSpaceToASCII and the underscore rule leave: ASCII with
// Unicode digits as 0-9 and Unicode whitespace as spaces, single underscores between digits
// removed. std::nullopt for any other character, an embedded NUL or a misplaced underscore.
std::optional<std::string> asciiNumberText(std::string_view text) {
    std::string out;
    char previous = '\0';
    for (const char32_t code : decodeUtf8(text)) {
        char c;
        if (code < 127) {
            c = static_cast<char>(code);
        } else if (isPythonSpace(code)) {
            c = ' ';
        } else if (isPythonDigit(code)) {
            c = static_cast<char>('0' + decimalValue(code));
        } else {
            return std::nullopt;
        }
        if (c == '\0') {
            return std::nullopt;
        }
        if (c == '_') {
            if (!isAsciiDigit(previous)) {
                return std::nullopt;
            }
        } else {
            if (previous == '_' && !isAsciiDigit(c)) {
                return std::nullopt;
            }
            out += c;
        }
        previous = c;
    }
    if (previous == '_') {
        return std::nullopt;
    }
    return out;
}

bool startsWithNoCase(std::string_view text, std::size_t at, std::string_view word) {
    if (text.size() - at < word.size()) {
        return false;
    }
    for (std::size_t i = 0; i < word.size(); ++i) {
        char c = text[at + i];
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
        if (c != word[i]) {
            return false;
        }
    }
    return true;
}

// round-half-even of digits * 10^exponent / 2^-1074, times 2^-1074: the exact value of a decimal
// below the normal range (some from_chars implementations report those as out of range)
double subnormalValue(const std::string& digits, std::int64_t exponent) {
    std::vector<unsigned char> big;  // decimal digits, least significant first
    for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
        big.push_back(static_cast<unsigned char>(*it - '0'));
    }
    for (int k = 0; k < 1074; ++k) {
        int carry = 0;
        for (unsigned char& d : big) {
            const int v = d * 2 + carry;
            d = static_cast<unsigned char>(v % 10);
            carry = v / 10;
        }
        if (carry) {
            big.push_back(static_cast<unsigned char>(carry));
        }
    }
    const auto fraction = static_cast<std::size_t>(-exponent);
    std::uint64_t n = 0;
    for (std::size_t i = big.size(); i-- > fraction;) {
        n = n * 10 + big[i];
    }
    if (fraction > 0 && fraction <= big.size()) {
        const int top = big[fraction - 1];
        bool rest = false;
        for (std::size_t i = 0; i + 1 < fraction; ++i) {
            rest = rest || big[i] != 0;
        }
        if (top > 5 || (top == 5 && (rest || (n & 1)))) {
            ++n;
        }
    }
    return std::ldexp(static_cast<double>(n), -1074);
}

// digits (no leading zeros, at least one non-zero) * 10^exponent, correctly rounded
double decimalToDouble(std::string digits, std::int64_t exponent) {
    while (digits.size() > 1 && digits.back() == '0') {
        digits.pop_back();
        ++exponent;
    }
    constexpr std::size_t kMaxDigits = 800;  // more than any halfway case needs
    if (digits.size() > kMaxDigits) {
        const bool sticky = digits.find_first_not_of('0', kMaxDigits) != std::string::npos;
        exponent += static_cast<std::int64_t>(digits.size() - kMaxDigits);
        digits.resize(kMaxDigits);
        if (sticky) {
            digits += '1';
            --exponent;
        }
    }
    const std::int64_t lead = exponent + static_cast<std::int64_t>(digits.size());
    if (lead > 310) {
        return std::numeric_limits<double>::infinity();
    }
    if (lead < -330) {
        return 0.0;
    }
    const std::string text = digits + "e" + std::to_string(exponent);
    double value = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    if (result.ec == std::errc::result_out_of_range) {
        return lead > 0 ? std::numeric_limits<double>::infinity()
                        : subnormalValue(digits, exponent);
    }
    return value;
}

// float(text) of a Python str; std::nullopt where Python raises ValueError
std::optional<double> parseFloat(std::string_view text) {
    const std::optional<std::string> ascii = asciiNumberText(text);
    if (!ascii) {
        return std::nullopt;
    }
    std::string_view s = *ascii;
    while (!s.empty() && isAsciiSpace(s.front())) {
        s.remove_prefix(1);
    }
    while (!s.empty() && isAsciiSpace(s.back())) {
        s.remove_suffix(1);
    }
    if (s.empty()) {
        return std::nullopt;
    }
    std::size_t i = 0;
    const bool negative = s[0] == '-';
    if (s[0] == '-' || s[0] == '+') {
        ++i;
    }
    std::string digits;
    std::size_t intDigits = 0, fracDigits = 0;
    while (i < s.size() && isAsciiDigit(s[i])) {
        digits += s[i++];
        ++intDigits;
    }
    if (i < s.size() && s[i] == '.') {
        ++i;
        while (i < s.size() && isAsciiDigit(s[i])) {
            digits += s[i++];
            ++fracDigits;
        }
    }
    const double sign = negative ? -1.0 : 1.0;
    if (intDigits + fracDigits == 0) {  // _Py_parse_inf_or_nan
        if (startsWithNoCase(s, i, "inf")) {
            i += 3;
            if (startsWithNoCase(s, i, "inity")) {
                i += 5;
            }
            if (i != s.size()) {
                return std::nullopt;
            }
            return sign * std::numeric_limits<double>::infinity();
        }
        if (startsWithNoCase(s, i, "nan") && i + 3 == s.size()) {
            return std::copysign(std::numeric_limits<double>::quiet_NaN(), sign);
        }
        return std::nullopt;
    }
    std::int64_t exponent = 0;
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        std::size_t j = i + 1;
        bool expNegative = false;
        if (j < s.size() && (s[j] == '+' || s[j] == '-')) {
            expNegative = s[j] == '-';
            ++j;
        }
        if (j < s.size() && isAsciiDigit(s[j])) {
            while (j < s.size() && isAsciiDigit(s[j])) {
                if (exponent < 1'000'000'000) {
                    exponent = exponent * 10 + (s[j] - '0');
                }
                ++j;
            }
            if (expNegative) {
                exponent = -exponent;
            }
            i = j;
        }
    }
    if (i != s.size()) {
        return std::nullopt;
    }
    exponent -= static_cast<std::int64_t>(fracDigits);
    const std::size_t nonZero = digits.find_first_not_of('0');
    if (nonZero == std::string::npos) {
        return sign * 0.0;
    }
    return sign * decimalToDouble(digits.substr(nonZero), exponent);
}


// int(text) of a Python str: the digits as a decimal string with its sign; std::nullopt where
// Python raises ValueError
std::optional<std::pair<bool, std::string>> parseInt(std::string_view text) {
    const std::optional<std::string> ascii = asciiNumberText(text);
    if (!ascii) {
        return std::nullopt;
    }
    std::string_view s = *ascii;
    while (!s.empty() && isAsciiSpace(s.front())) {
        s.remove_prefix(1);
    }
    while (!s.empty() && isAsciiSpace(s.back())) {
        s.remove_suffix(1);
    }
    bool negative = false;
    if (!s.empty() && (s[0] == '-' || s[0] == '+')) {
        negative = s[0] == '-';
        s.remove_prefix(1);
    }
    if (s.empty()) {
        return std::nullopt;
    }
    for (const char c : s) {
        if (!isAsciiDigit(c)) {
            return std::nullopt;
        }
    }
    const std::size_t nonZero = s.find_first_not_of('0');
    const std::string digits =
        nonZero == std::string_view::npos ? "0" : std::string(s.substr(nonZero));
    return std::pair<bool, std::string>(negative, digits);
}

// "%.200R": the repr cut to 200 characters
std::string truncatedRepr(std::string_view text) {
    const std::u32string repr = decodeUtf8(pythonRepr(Value(std::string(text))));
    return encodeUtf8(std::u32string_view(repr).substr(0, 200));
}

std::optional<std::int64_t> toInt64(bool negative, const std::string& digits) {
    if (digits.size() > 19) {
        return std::nullopt;
    }
    std::uint64_t magnitude = 0;
    for (const char c : digits) {
        magnitude = magnitude * 10 + static_cast<std::uint64_t>(c - '0');
    }
    constexpr auto kMax = static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    if (digits.size() == 19 && digits > "9223372036854775808") {
        return std::nullopt;
    }
    if (negative) {
        if (magnitude > kMax + 1) {
            return std::nullopt;
        }
        return magnitude == kMax + 1 ? std::numeric_limits<std::int64_t>::min()
                                     : -static_cast<std::int64_t>(magnitude);
    }
    if (magnitude > kMax) {
        return std::nullopt;
    }
    return static_cast<std::int64_t>(magnitude);
}

}  // namespace

std::string pythonRepr(const Value& value) {
    std::string out;
    appendRepr(out, value);
    return out;
}

std::string pythonStr(const Value& value) {
    if (const std::string* text = value.get<std::string>()) {
        return *text;
    }
    if (const DateTime* dt = value.get<DateTime>()) {
        std::string text = dt->isoformat();
        const bool date = dt->kind == DateTime::Kind::LocalDate;
        if (!date && dt->kind != DateTime::Kind::LocalTime) {
            text[10] = ' ';  // datetime.__str__ is isoformat(sep=" ")
        }
        return text;
    }
    return pythonRepr(value);
}

std::string_view pythonTypeName(const Value& value) {
    switch (value.type()) {
    case Value::Type::Null:
        return "NoneType";
    case Value::Type::Bool:
        return "bool";
    case Value::Type::Int:
        return "int";
    case Value::Type::Float:
        return "float";
    case Value::Type::String:
        return "str";
    case Value::Type::DateTime:
        switch (value.get<DateTime>()->kind) {
        case DateTime::Kind::LocalDate:
            return "date";
        case DateTime::Kind::LocalTime:
            return "time";
        case DateTime::Kind::LocalDateTime:
        case DateTime::Kind::OffsetDateTime:
            break;
        }
        return "datetime";
    case Value::Type::Array:
        return "list";
    case Value::Type::Table:
        return "dict";
    }
    return "object";
}

double pythonRound(double x, int digits) {
    // float.__round__: non-finite values and digits beyond what a double holds round to
    // themselves, digits below -308 give a zero of x's sign
    if (!std::isfinite(x) || digits > 323) {
        return x;
    }
    if (digits < -308) {
        return 0.0 * x;
    }
    std::string text;
    if (digits >= 0) {
        // to_chars with a precision rounds the exact binary value half to even, as Python's dtoa
        char buffer[700];  // sign, 309 integer digits, the point and up to 323 decimals
        const auto result =
            std::to_chars(buffer, buffer + sizeof buffer, x, std::chars_format::fixed, digits);
        text.assign(buffer, result.ptr);
    } else {
        text = roundToTens(x, digits);
    }
    double rounded = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), rounded);
    if (result.ec == std::errc::result_out_of_range) {
        throw std::overflow_error("rounded value too large to represent");
    }
    return rounded;
}

std::string pythonJsonCompact(const Value& value) {
    std::string out;
    appendJsonCompact(out, value);
    return out;
}

std::int64_t pythonInt(std::string_view text) {
    const auto parsed = parseInt(text);
    if (!parsed) {
        throw ValueError("invalid literal for int() with base 10: " + truncatedRepr(text));
    }
    const std::optional<std::int64_t> value = toInt64(parsed->first, parsed->second);
    if (!value) {
        throw std::overflow_error("Python int too large to convert to C long");
    }
    return *value;
}

double pythonFloat(std::string_view text) {
    const std::optional<double> value = parseFloat(text);
    if (!value) {
        throw ValueError("could not convert string to float: " +
                         pythonRepr(Value(std::string(text))));
    }
    return *value;
}

std::optional<std::int64_t> tryPythonInt(std::string_view text) {
    const auto parsed = parseInt(text);
    if (!parsed) {
        return std::nullopt;
    }
    return toInt64(parsed->first, parsed->second);
}

std::optional<double> tryPythonFloat(std::string_view text) { return parseFloat(text); }

namespace posixpath {

std::string expanduser(std::string_view path) {
    if (path.empty() || path.front() != '~') {
        return std::string(path);
    }
    std::size_t i = path.find('/', 1);
    if (i == std::string_view::npos) {
        i = path.size();
    }
    const std::optional<std::string> home =
        i == 1 ? currentUserHome() : userHome(std::string(path.substr(1, i - 1)));
    if (!home) {
        return std::string(path);
    }
    std::string out = *home;
    while (!out.empty() && out.back() == '/') {
        out.pop_back();
    }
    out += path.substr(i);
    return out.empty() ? std::string("/") : out;
}

bool isabs(std::string_view path) { return !path.empty() && path.front() == '/'; }

std::string join(std::string_view base, std::string_view path) {
    if (isabs(path)) {
        return std::string(path);
    }
    std::string out(base);
    if (!out.empty() && out.back() != '/') {
        out += '/';
    }
    out += path;
    return out;
}

std::string normpath(std::string_view path) {
    if (path.empty()) {
        return ".";
    }
    // POSIX keeps exactly two leading slashes ("//" may mean something else); more become one
    std::size_t slashes = 0;
    if (path.front() == '/') {
        slashes = path.size() >= 2 && path[1] == '/' && (path.size() < 3 || path[2] != '/') ? 2 : 1;
    }
    std::vector<std::string_view> parts;
    std::size_t start = 0;
    while (start <= path.size()) {
        std::size_t end = path.find('/', start);
        if (end == std::string_view::npos) {
            end = path.size();
        }
        const std::string_view part = path.substr(start, end - start);
        start = end + 1;
        if (part.empty() || part == ".") {
            continue;
        }
        const bool keep = part != ".." || (!slashes && parts.empty()) ||
                          (!parts.empty() && parts.back() == "..");
        if (keep) {
            parts.push_back(part);
        } else if (!parts.empty()) {
            parts.pop_back();
        }
    }
    std::string out(slashes, '/');
    for (std::size_t k = 0; k < parts.size(); ++k) {
        if (k) {
            out += '/';
        }
        out += parts[k];
    }
    return out.empty() ? std::string(".") : out;
}

std::string abspath(std::string_view path, std::string_view cwd) {
    return normpath(isabs(path) ? std::string(path) : join(cwd, path));
}

std::string dirname(std::string_view path) {
    const std::size_t slash = path.rfind('/');
    std::string head(slash == std::string_view::npos ? std::string_view()
                                                     : path.substr(0, slash + 1));
    if (head.find_first_not_of('/') != std::string::npos) {  // not only slashes: strip the last
        while (head.back() == '/') {
            head.pop_back();
        }
    }
    return head;
}

std::string basename(std::string_view path) {
    const std::size_t slash = path.rfind('/');
    return std::string(slash == std::string_view::npos ? path : path.substr(slash + 1));
}

std::pair<std::string, std::string> splitext(std::string_view path) {
    // genericpath._splitext: the last dot of the last component, unless only dots precede it
    const std::size_t slash = path.rfind('/');
    const std::size_t dot = path.rfind('.');
    const std::size_t name = slash == std::string_view::npos ? 0 : slash + 1;
    if (dot != std::string_view::npos && (slash == std::string_view::npos || dot > slash)) {
        for (std::size_t k = name; k < dot; ++k) {
            if (path[k] != '.') {
                return {std::string(path.substr(0, dot)), std::string(path.substr(dot))};
            }
        }
    }
    return {std::string(path), std::string()};
}

}  // namespace posixpath

namespace {

struct OsErrorEntry {
    std::errc condition;
    int number;        // errno on Linux
    const char* text;  // glibc's strerror()
};

constexpr OsErrorEntry kOsErrors[] = {
    {std::errc::operation_not_permitted, 1, "Operation not permitted"},
    {std::errc::no_such_file_or_directory, 2, "No such file or directory"},
    {std::errc::io_error, 5, "Input/output error"},
    {std::errc::bad_file_descriptor, 9, "Bad file descriptor"},
    {std::errc::not_enough_memory, 12, "Cannot allocate memory"},
    {std::errc::permission_denied, 13, "Permission denied"},
    {std::errc::device_or_resource_busy, 16, "Device or resource busy"},
    {std::errc::file_exists, 17, "File exists"},
    {std::errc::cross_device_link, 18, "Invalid cross-device link"},
    {std::errc::not_a_directory, 20, "Not a directory"},
    {std::errc::is_a_directory, 21, "Is a directory"},
    {std::errc::invalid_argument, 22, "Invalid argument"},
    {std::errc::too_many_files_open, 24, "Too many open files"},
    {std::errc::file_too_large, 27, "File too large"},
    {std::errc::no_space_on_device, 28, "No space left on device"},
    {std::errc::read_only_file_system, 30, "Read-only file system"},
    {std::errc::filename_too_long, 36, "File name too long"},
    {std::errc::directory_not_empty, 39, "Directory not empty"},
    {std::errc::too_many_symbolic_link_levels, 40, "Too many levels of symbolic links"},
};

}  // namespace

std::string pythonOsErrorText(const std::error_code& code, const std::filesystem::path& path1,
                              const std::filesystem::path& path2) {
    int number = code.value();
    std::string text = code.message();
    const std::error_condition condition = code.default_error_condition();
    if (condition.category() == std::generic_category()) {
        for (const OsErrorEntry& entry : kOsErrors) {
            if (condition.value() == static_cast<int>(entry.condition)) {
                number = entry.number;
                text = entry.text;
                break;
            }
        }
    }
    std::string out = "[Errno " + std::to_string(number) + "] " + text;
    if (!path1.empty()) {
        out += ": " + pythonRepr(Value(pathToUtf8(path1)));
        if (!path2.empty()) {
            out += " -> " + pythonRepr(Value(pathToUtf8(path2)));
        }
    }
    return out;
}

std::string pythonOsErrorText(const std::filesystem::filesystem_error& error) {
    return pythonOsErrorText(error.code(), error.path1(), error.path2());
}

std::string currentDirectory() { return pathToUtf8(std::filesystem::current_path()); }

}  // namespace mdbridge::core
