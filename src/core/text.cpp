#include "mdbridge/core/text.hpp"

#include <algorithm>
#include <iterator>

namespace mdbridge::core {

namespace {

#include "unicode_tables.inc"

template <std::size_t N>
bool inRanges(const CodeRange (&table)[N], char32_t code) {
    const auto* end = std::end(table);
    const auto* it = std::upper_bound(std::begin(table), end, code,
                                      [](char32_t c, const CodeRange& r) { return c < r.first; });
    return it != std::begin(table) && code <= std::prev(it)->last;
}

constexpr char32_t kReplacement = 0xFFFD;

bool isContinuation(unsigned char byte) { return (byte & 0xC0) == 0x80; }

// Decodes the sequence at text[i]; returns its length, or 0 when it is not valid UTF-8.
std::size_t decodeOne(std::string_view text, std::size_t i, char32_t& code) {
    const auto b0 = static_cast<unsigned char>(text[i]);
    const auto at = [&](std::size_t k) { return static_cast<unsigned char>(text[i + k]); };
    const std::size_t left = text.size() - i;
    if (b0 < 0x80) {
        code = b0;
        return 1;
    }
    if (b0 >= 0xC2 && b0 <= 0xDF && left >= 2 && isContinuation(at(1))) {
        code = (static_cast<char32_t>(b0 & 0x1F) << 6) | (at(1) & 0x3F);
        return 2;
    }
    if (b0 >= 0xE0 && b0 <= 0xEF && left >= 3 && isContinuation(at(1)) && isContinuation(at(2))) {
        const bool ok = (b0 == 0xE0 && at(1) >= 0xA0) || (b0 == 0xED && at(1) <= 0x9F) || (b0 != 0xE0 && b0 != 0xED);
        if (ok) {
            code = (static_cast<char32_t>(b0 & 0x0F) << 12) | (static_cast<char32_t>(at(1) & 0x3F) << 6) | (at(2) & 0x3F);
            return 3;
        }
    }
    if (b0 >= 0xF0 && b0 <= 0xF4 && left >= 4 && isContinuation(at(1)) && isContinuation(at(2)) &&
        isContinuation(at(3))) {
        const bool ok = (b0 == 0xF0 && at(1) >= 0x90) || (b0 == 0xF4 && at(1) <= 0x8F) || (b0 > 0xF0 && b0 < 0xF4);
        if (ok) {
            code = (static_cast<char32_t>(b0 & 0x07) << 18) | (static_cast<char32_t>(at(1) & 0x3F) << 12) |
                   (static_cast<char32_t>(at(2) & 0x3F) << 6) | (at(3) & 0x3F);
            return 4;
        }
    }
    return 0;
}

bool isAsciiWhitespace(char32_t c) { return c == U' ' || c == U'\t' || c == U'\n' || c == U'\v' || c == U'\f' || c == U'\r'; }

// textwrap's character classes: \w, letter = [^\d\W], word_punct = [\w!"'&.,?]
bool isLetter(char32_t c) { return isPythonWordCharacter(c) && !isPythonDigit(c); }
bool isWordPunct(char32_t c) {
    return isPythonWordCharacter(c) || c == U'!' || c == U'"' || c == U'\'' || c == U'&' || c == U'.' || c == U',' ||
           c == U'?';
}

bool isBlank(std::u32string_view chunk) {
    return std::all_of(chunk.begin(), chunk.end(), [](char32_t c) { return isPythonSpace(c); });
}

// textwrap.TextWrapper._munge_whitespace: expandtabs(8), then each whitespace character -> ' '
std::u32string mungeWhitespace(std::u32string_view text) {
    std::u32string out;
    std::size_t column = 0;
    for (const char32_t c : text) {
        if (c == U'\t') {
            const std::size_t spaces = 8 - column % 8;
            out.append(spaces, U' ');
            column += spaces;
            continue;
        }
        column = (c == U'\n' || c == U'\r') ? 0 : column + 1;
        out.push_back(isAsciiWhitespace(c) ? U' ' : c);
    }
    return out;
}

// textwrap.TextWrapper.wordsep_re (break_on_hyphens): whitespace runs, em-dashes between words, and
// words, split after a hyphen between letters ("self-|contained")
std::vector<std::u32string> splitChunks(const std::u32string& s) {
    std::vector<std::u32string> chunks;
    const std::size_t n = s.size();
    const auto ws = [&](std::size_t k) { return k < n && isAsciiWhitespace(s[k]); };
    const auto letter = [&](std::size_t k) { return k < n && isLetter(s[k]); };
    const auto hyphen = [&](std::size_t k) { return k < n && s[k] == U'-'; };
    // -{2,} followed by \w, starting at k
    const auto emDashAhead = [&](std::size_t k) {
        std::size_t j = k;
        while (hyphen(j)) {
            ++j;
        }
        return j - k >= 2 && j < n && isPythonWordCharacter(s[j]);
    };
    std::size_t i = 0;
    while (i < n) {
        std::size_t end = 0;
        if (ws(i)) {  // whitespace+
            end = i;
            while (ws(end)) {
                ++end;
            }
        } else if (i > 0 && isWordPunct(s[i - 1]) && emDashAhead(i)) {  // (?<=wp) -{2,} (?=\w)
            end = i;
            while (hyphen(end)) {
                ++end;
            }
        } else {  // nws+? followed by a hyphen break, the end of the word, or an em-dash
            for (std::size_t k = i + 1; k <= n; ++k) {
                const bool hyphenBreak =
                    hyphen(k) && ((k >= 2 && letter(k - 2) && letter(k - 1)) ||
                                  (k >= 3 && letter(k - 3) && hyphen(k - 2) && letter(k - 1))) &&
                    letter(k + 1) && (letter(k + 2) || (hyphen(k + 2) && letter(k + 3)));
                if (hyphenBreak) {
                    end = k + 1;
                    break;
                }
                if (k == n || ws(k)) {
                    end = k;
                    break;
                }
                if (isWordPunct(s[k - 1]) && emDashAhead(k)) {
                    end = k;
                    break;
                }
            }
        }
        chunks.emplace_back(s.substr(i, end - i));
        i = end;
    }
    return chunks;
}

}  // namespace

bool isPythonSpace(char32_t code) { return inRanges(kSpaceRanges, code); }
bool isPythonWordCharacter(char32_t code) { return inRanges(kWordRanges, code); }
bool isPythonDigit(char32_t code) { return inRanges(kDigitRanges, code); }
bool isPythonPrintable(char32_t code) { return inRanges(kPrintableRanges, code); }

std::u32string decodeUtf8(std::string_view text) {
    std::u32string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size();) {
        char32_t code = 0;
        const std::size_t length = decodeOne(text, i, code);
        out.push_back(length ? code : kReplacement);
        i += length ? length : 1;
    }
    return out;
}

std::u32string decodeUtf8Surrogateescape(std::string_view text) {
    std::u32string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size();) {
        char32_t code = 0;
        const std::size_t length = decodeOne(text, i, code);
        out.push_back(length ? code : char32_t{0xDC00} + static_cast<unsigned char>(text[i]));
        i += length ? length : 1;
    }
    return out;
}

void appendUtf8(std::string& out, char32_t code) {
    if (code < 0x80) {
        out += static_cast<char>(code);
    } else if (code < 0x800) {
        out += static_cast<char>(0xC0 | (code >> 6));
        out += static_cast<char>(0x80 | (code & 0x3F));
    } else if (code < 0x10000) {
        out += static_cast<char>(0xE0 | (code >> 12));
        out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
        out += static_cast<char>(0x80 | (code & 0x3F));
    } else {
        out += static_cast<char>(0xF0 | (code >> 18));
        out += static_cast<char>(0x80 | ((code >> 12) & 0x3F));
        out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
        out += static_cast<char>(0x80 | (code & 0x3F));
    }
}

std::string encodeUtf8(std::u32string_view text) {
    std::string out;
    out.reserve(text.size());
    for (const char32_t c : text) {
        appendUtf8(out, c);
    }
    return out;
}

std::size_t characterCount(std::string_view text) {
    std::size_t count = 0;
    for (std::size_t i = 0; i < text.size(); ++count) {
        char32_t code = 0;
        const std::size_t length = decodeOne(text, i, code);
        i += length ? length : 1;
    }
    return count;
}

std::vector<std::string> pythonSplitlines(std::string_view text) {
    const std::u32string s = decodeUtf8(text);
    std::vector<std::string> lines;
    std::size_t start = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char32_t c = s[i];
        const bool boundary = c == U'\n' || c == U'\r' || c == U'\v' || c == U'\f' || c == 0x1C || c == 0x1D ||
                              c == 0x1E || c == 0x85 || c == 0x2028 || c == 0x2029;
        if (!boundary) {
            continue;
        }
        lines.push_back(encodeUtf8(std::u32string_view(s).substr(start, i - start)));
        if (c == U'\r' && i + 1 < s.size() && s[i + 1] == U'\n') {
            ++i;
        }
        start = i + 1;
    }
    if (start < s.size()) {
        lines.push_back(encodeUtf8(std::u32string_view(s).substr(start)));
    }
    return lines;
}

std::string pythonRstrip(std::string_view text) {
    std::u32string s = decodeUtf8(text);
    while (!s.empty() && isPythonSpace(s.back())) {
        s.pop_back();
    }
    return encodeUtf8(s);
}

std::string pythonStrip(std::string_view text) {
    const std::u32string s = decodeUtf8(text);
    std::size_t a = 0;
    std::size_t b = s.size();
    while (a < b && isPythonSpace(s[a])) {
        ++a;
    }
    while (b > a && isPythonSpace(s[b - 1])) {
        --b;
    }
    return encodeUtf8(std::u32string_view(s).substr(a, b - a));
}

std::vector<std::string> pythonTextwrap(std::string_view text, int width) {
    // TextWrapper._wrap_chunks with the default options (no indents, drop_whitespace)
    std::vector<std::u32string> chunks = splitChunks(mungeWhitespace(decodeUtf8(text)));
    std::reverse(chunks.begin(), chunks.end());
    std::vector<std::string> lines;
    const auto w = static_cast<std::size_t>(std::max(width, 1));
    while (!chunks.empty()) {
        std::vector<std::u32string> line;
        std::size_t length = 0;
        if (isBlank(chunks.back()) && !lines.empty()) {  // no whitespace at the start of later lines
            chunks.pop_back();
        }
        while (!chunks.empty() && length + chunks.back().size() <= w) {
            length += chunks.back().size();
            line.push_back(std::move(chunks.back()));
            chunks.pop_back();
        }
        if (!chunks.empty() && chunks.back().size() > w) {  // _handle_long_word (break_long_words)
            const std::size_t spaceLeft = width < 1 ? 1 : w - length;
            std::u32string& chunk = chunks.back();
            std::size_t end = spaceLeft;
            if (chunk.size() > spaceLeft) {  // break after the last hyphen, if non-hyphens precede it
                const std::size_t hyphenAt = spaceLeft == 0 ? std::u32string::npos : chunk.rfind(U'-', spaceLeft - 1);
                if (hyphenAt != std::u32string::npos && hyphenAt > 0 &&
                    std::any_of(chunk.begin(), chunk.begin() + static_cast<std::ptrdiff_t>(hyphenAt),
                                [](char32_t c) { return c != U'-'; })) {
                    end = hyphenAt + 1;
                }
            }
            if (end > 0) {  // Python appends nothing when the line has no room left
                line.push_back(chunk.substr(0, end));
                chunk.erase(0, end);
            }
        }
        if (!line.empty() && isBlank(line.back())) {
            line.pop_back();
        }
        if (!line.empty()) {
            std::u32string joined;
            for (const std::u32string& part : line) {
                joined += part;
            }
            lines.push_back(encodeUtf8(joined));
        }
    }
    return lines;
}

}  // namespace mdbridge::core
