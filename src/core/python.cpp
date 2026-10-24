#include "mdbridge/core/python.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <vector>

#if !defined(_WIN32)
#include <cerrno>
#include <pwd.h>
#include <unistd.h>
#endif

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
// No reference behaviour: MDATK never ran on Windows. "~" is the profile folder, as in Python's
// ntpath; "~user" stays as it is.
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

std::string currentDirectory() { return pathToUtf8(std::filesystem::current_path()); }

}  // namespace mdbridge::core
