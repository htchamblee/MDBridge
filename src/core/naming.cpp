#include "mdbridge/core/naming.hpp"

#include <cstdio>
#include <vector>

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

bool isAsciiAlnum(char32_t c) {
    return (c >= U'A' && c <= U'Z') || (c >= U'a' && c <= U'z') || (c >= U'0' && c <= U'9');
}

// re.sub(r"[^<allowed>]+", "-", text.strip()).strip(trim) or fallback. The classes are ASCII, so
// working on UTF-8 bytes gives Python's result: every byte of another character is outside them
// and a run of such bytes becomes one "-".
std::string sanitize(std::string_view text, bool keepDotAndUnderscore, std::string_view trim,
                     std::string_view fallback) {
    std::string out;
    bool inRun = false;
    for (const char c : pythonStrip(text)) {
        const auto u = static_cast<unsigned char>(c);
        const bool allowed =
            isAsciiAlnum(u) || c == '-' || (keepDotAndUnderscore && (c == '.' || c == '_'));
        if (allowed) {
            out += c;
        } else if (!inRun) {
            out += '-';
        }
        inRun = !allowed;
    }
    const std::size_t first = out.find_first_not_of(trim);
    if (first == std::string::npos) {
        return std::string(fallback);
    }
    return out.substr(first, out.find_last_not_of(trim) - first + 1);
}

// UTF-8 length of a code point of decodeUtf8Surrogateescape (an escaped byte is one byte)
std::size_t encodedLength(char32_t c) {
    if (c >= 0xDC80 && c <= 0xDCFF) {
        return 1;
    }
    return c < 0x80 ? 1 : c < 0x800 ? 2 : c < 0x10000 ? 3 : 4;
}

}  // namespace

std::string newStamp(const std::optional<DateTime>& now) {
    const DateTime t = now ? *now : clockNow();
    char stamp[64];  // %Y is not padded (glibc): year 999 gives "999..."
    std::snprintf(stamp, sizeof stamp, "%d%02d%02d-%02d%02d", t.year, t.month, t.day, t.hour,
                  t.minute);
    return stamp;
}

std::string sanitizeProject(std::string_view name) {
    return sanitize(name, true, "-._", "project");
}

std::string sanitizeTag(std::string_view tag) { return sanitize(tag, false, "-", "step"); }

std::string stepLabel(std::int64_t step) {
    const std::string digits = std::to_string(step);
    return digits.size() < 2 ? "0" + digits : digits;
}

std::string fileName(std::string_view project, std::int64_t step, std::string_view tag,
                     std::string_view stamp, std::string_view ext) {
    std::string name = sanitizeProject(project) + "_" + stepLabel(step) + "_" + sanitizeTag(tag);
    if (!stamp.empty()) {
        name += "_";
        name += stamp;
    }
    const std::size_t dots = ext.find_first_not_of('.');
    if (dots != std::string_view::npos) {
        name += ".";
        name += ext.substr(dots);
    }
    return name;
}

std::string stepFolder(std::int64_t step, std::string_view tag) {
    return stepLabel(step) + "_" + sanitizeTag(tag);
}

std::optional<ParsedName> parseName(std::string_view name) {
    // re.match(r"^(?P<project>.+)_(?P<step>\d{2,})_(?P<tag>[A-Za-z0-9-]+)"
    //          r"(?:_(?P<stamp>\d{8}-\d{4}))?(?:\.(?P<ext>[A-Za-z0-9]+))?$", name)
    // with Python's backtracking order, \d = any Unicode decimal digit, "." = anything but a
    // newline, and "$" also before a final newline
    const std::u32string s = decodeUtf8Surrogateescape(name);
    const std::size_t n = s.size();
    std::vector<std::size_t> offset(n + 1, 0);  // byte offset of each character
    for (std::size_t k = 0; k < n; ++k) {
        offset[k + 1] = offset[k] + encodedLength(s[k]);
    }
    const auto part = [&](std::size_t from, std::size_t to) {
        return std::string(name.substr(offset[from], offset[to] - offset[from]));
    };
    const auto run = [&](std::size_t from, bool (*accept)(char32_t)) {
        std::size_t to = from;
        while (to < n && accept(s[to])) {
            ++to;
        }
        return to - from;
    };
    const auto isTag = [](char32_t c) { return isAsciiAlnum(c) || c == U'-'; };
    const auto isExt = [](char32_t c) { return isAsciiAlnum(c); };
    const auto isEnd = [&](std::size_t at) { return at == n || (at + 1 == n && s[at] == U'\n'); };
    const auto isStamp = [&](std::size_t at) {  // _\d{8}-\d{4}
        if (at + 14 > n || s[at] != U'_' || s[at + 9] != U'-') {
            return false;
        }
        for (std::size_t k = at + 1; k < at + 14; ++k) {
            if (k != at + 9 && !isPythonDigit(s[k])) {
                return false;
            }
        }
        return true;
    };

    std::size_t longest = 0;  // .+ stops at the first newline
    while (longest < n && s[longest] != U'\n') {
        ++longest;
    }
    for (std::size_t p = longest; p >= 1; --p) {
        if (p >= n || s[p] != U'_') {
            continue;
        }
        const std::size_t step0 = p + 1;
        for (std::size_t digits = run(step0, isPythonDigit); digits >= 2; --digits) {
            const std::size_t step1 = step0 + digits;
            if (step1 >= n || s[step1] != U'_') {
                continue;
            }
            const std::size_t tag0 = step1 + 1;
            for (std::size_t tagLength = run(tag0, isTag); tagLength >= 1; --tagLength) {
                const std::size_t tag1 = tag0 + tagLength;
                for (const bool withStamp : {true, false}) {
                    if (withStamp && !isStamp(tag1)) {
                        continue;
                    }
                    const std::size_t at = withStamp ? tag1 + 14 : tag1;
                    ParsedName parsed{part(0, p), part(step0, step1), part(tag0, tag1),
                                      withStamp ? part(tag1 + 1, at) : std::string(), {}};
                    if (at < n && s[at] == U'.') {
                        for (std::size_t length = run(at + 1, isExt); length >= 1; --length) {
                            if (isEnd(at + 1 + length)) {
                                parsed.ext = part(at + 1, at + 1 + length);
                                return parsed;
                            }
                        }
                    }
                    if (isEnd(at)) {
                        return parsed;
                    }
                }
            }
        }
    }
    return std::nullopt;
}

}  // namespace mdbridge::core
