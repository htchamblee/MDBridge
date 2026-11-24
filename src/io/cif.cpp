// The glob and fnmatch pattern handling of find_dictionary (Lib/glob.py, Lib/fnmatch.py)
// is ported from CPython (Python Software Foundation License Version 2; full text in
// packaging/licenses/Python-LICENSE.txt, see THIRD_PARTY_NOTICES.md):
// Copyright (c) 2001-2024 Python Software Foundation; All Rights Reserved

#include "mdbridge/io/cif.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

#if !defined(_WIN32)
#include <sys/stat.h>
#endif

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace fs = std::filesystem;
using core::decodeUtf8;
using core::encodeUtf8;
using core::isPythonSpace;
using core::pythonLower;

namespace {

// Text as Python holds it after decoding with errors="replace": every maximal invalid UTF-8
// subpart becomes one U+FFFD.
std::string normalised(std::string_view text) {
    return encodeUtf8(detail::decodeUtf8Replace(text));
}

std::string_view beforeDot(std::string_view tag) { return tag.substr(0, tag.find('.')); }

std::string_view afterDot(std::string_view tag) {
    const std::size_t dot = tag.find('.');
    return dot == std::string_view::npos ? tag : tag.substr(dot + 1);
}

bool startsWith(std::u32string_view text, std::u32string_view prefix) {
    return text.substr(0, prefix.size()) == prefix;
}

// str.lower() of a word, as far as the keyword tests need it: of the characters that lower to
// ASCII only the Kelvin sign (to "k") and U+0130 (to "i" + U+0307) are not ASCII, and no keyword
// holds k or i.
std::u32string keywordLower(std::u32string_view word) {
    std::u32string out(word.substr(0, 8));
    for (char32_t& c : out) {
        if (c >= U'A' && c <= U'Z') {
            c = static_cast<char32_t>(c + 32);
        }
    }
    return out;
}

// _cif_tokens, one token at a time (parse_cif meets its errors in the order of the text).
class Tokenizer {
public:
    explicit Tokenizer(std::string_view text) : m_lines(core::pythonSplitlines(normalised(text))) {}

    bool next(CifToken& token) {
        for (;;) {
            if (!m_inLine) {
                if (m_index >= m_lines.size()) {
                    return false;
                }
                if (!m_lines[m_index].empty() && m_lines[m_index][0] == ';') {
                    token = textField();
                    return true;
                }
                m_line = decodeUtf8(m_lines[m_index]);
                m_pos = 0;
                m_inLine = true;
            }
            if (word(token)) {
                return true;
            }
            m_inLine = false;
            ++m_index;
        }
    }

private:
    CifToken textField() {
        std::vector<std::string> buf{m_lines[m_index].substr(1)};
        ++m_index;
        while (m_index < m_lines.size() &&
               (m_lines[m_index].empty() || m_lines[m_index][0] != ';')) {
            buf.push_back(m_lines[m_index]);
            ++m_index;
        }
        if (m_index >= m_lines.size()) {
            throw DictionaryError("mmCIF: unterminated text field");
        }
        ++m_index;
        std::size_t first = core::pythonStrip(buf[0]).empty() ? 1 : 0;
        std::string value;
        for (std::size_t k = first; k < buf.size(); ++k) {
            value += (k > first ? "\n" : "") + buf[k];
        }
        return {CifToken::Kind::Value, std::move(value)};
    }

    // The next token of the current line, false at its end.
    bool word(CifToken& token) {
        const std::u32string& line = m_line;
        const std::size_t n = line.size();
        while (m_pos < n) {
            const char32_t c = line[m_pos];
            if (isPythonSpace(c)) {
                ++m_pos;
                continue;
            }
            if (c == U'#') {
                token = {CifToken::Kind::Comment, encodeUtf8(line.substr(m_pos))};
                m_pos = n;
                return true;
            }
            if (c == U'\'' || c == U'"') {
                std::size_t end = m_pos + 1;
                for (;;) {
                    end = line.find(c, end);
                    if (end == std::u32string::npos) {
                        throw DictionaryError("mmCIF line " + std::to_string(m_index + 1) +
                                              ": unterminated quote");
                    }
                    if (end + 1 >= n || isPythonSpace(line[end + 1])) {
                        break;
                    }
                    ++end;
                }
                token = {CifToken::Kind::Value,
                         encodeUtf8(std::u32string_view(line).substr(m_pos + 1, end - m_pos - 1))};
                m_pos = end + 1;
                return true;
            }
            std::size_t end = m_pos;
            while (end < n && !isPythonSpace(line[end])) {
                ++end;
            }
            const std::u32string_view w = std::u32string_view(line).substr(m_pos, end - m_pos);
            m_pos = end;
            const std::u32string low = keywordLower(w);
            if (startsWith(low, U"data_")) {
                token = {CifToken::Kind::Data, encodeUtf8(w.substr(5))};
            } else if (low == U"loop_" && w.size() == 5) {
                token = {CifToken::Kind::Loop, ""};
            } else if (startsWith(low, U"save_") || (low == U"global_" && w.size() == 7) ||
                       (low == U"stop_" && w.size() == 5)) {
                token = {CifToken::Kind::Other, encodeUtf8(w)};
            } else if (w[0] == U'_') {
                token = {CifToken::Kind::Tag, encodeUtf8(w)};
            } else {
                token = {CifToken::Kind::Value, encodeUtf8(w)};
            }
            return true;
        }
        return false;
    }

    std::vector<std::string> m_lines;
    std::size_t m_index = 0;
    std::u32string m_line;
    std::size_t m_pos = 0;
    bool m_inLine = false;
};

std::string ljust(const std::string& text, std::size_t width) {
    const std::size_t length = core::characterCount(text);
    return length >= width ? text : text + std::string(width - length, ' ');
}

}  // namespace

// ---- documents -------------------------------------------------------------------------------

std::string CifLoop::category() const {
    return tags.empty() ? std::string() : pythonLower(beforeDot(tags[0]));
}

std::ptrdiff_t CifLoop::col(std::string_view item) const {
    const std::string want = pythonLower(item);
    for (std::size_t k = 0; k < tags.size(); ++k) {
        if (pythonLower(afterDot(tags[k])) == want) {
            return static_cast<std::ptrdiff_t>(k);
        }
    }
    return -1;
}

std::vector<std::string> CifLoop::column(std::string_view item) const {
    const std::ptrdiff_t k = col(item);
    std::vector<std::string> out;
    if (k >= 0) {
        for (const std::vector<std::string>& row : rows) {
            out.push_back(row.at(static_cast<std::size_t>(k)));
        }
    }
    return out;
}

const CifLoop* CifBlock::get(std::string_view category) const {
    const std::string want = pythonLower(category);
    for (const CifLoop& item : items) {
        if (item.category() == want) {
            return &item;
        }
    }
    return nullptr;
}

CifLoop* CifBlock::get(std::string_view category) {
    return const_cast<CifLoop*>(std::as_const(*this).get(category));
}

const CifBlock* CifDocument::block(std::string_view name) const {
    const std::string want = pythonLower(name);
    for (const CifBlock& b : blocks) {
        if (pythonLower(b.name) == want) {
            return &b;
        }
    }
    return nullptr;
}

CifBlock* CifDocument::block(std::string_view name) {
    return const_cast<CifBlock*>(std::as_const(*this).block(name));
}

// ---- reading ---------------------------------------------------------------------------------

std::vector<CifToken> cifTokens(std::string_view text) {
    Tokenizer tokenizer(text);
    std::vector<CifToken> out;
    CifToken token;
    while (tokenizer.next(token)) {
        out.push_back(std::move(token));
    }
    return out;
}

CifDocument parseCif(std::string_view text) {
    CifDocument doc;
    CifBlock* block = nullptr;
    CifLoop* loop = nullptr;
    std::vector<std::string> values;
    std::optional<std::string> tag;

    const auto closeLoop = [&] {
        if (!loop) {
            return;
        }
        if (loop->tags.empty() || values.size() % loop->tags.size() != 0) {
            const std::string category = loop->category();
            throw DictionaryError("mmCIF: loop " + (category.empty() ? "?" : category) + " has " +
                                  std::to_string(values.size()) + " values for " +
                                  std::to_string(loop->tags.size()) + " columns");
        }
        const std::size_t w = loop->tags.size();
        loop->rows.clear();
        for (std::size_t k = 0; k < values.size(); k += w) {
            loop->rows.emplace_back(values.begin() + static_cast<std::ptrdiff_t>(k),
                                    values.begin() + static_cast<std::ptrdiff_t>(k + w));
        }
        loop = nullptr;
        values.clear();
    };
    const auto noValue = [&] { return DictionaryError("mmCIF: " + *tag + " has no value"); };

    Tokenizer tokenizer(text);
    CifToken token;
    while (tokenizer.next(token)) {
        std::string& value = token.text;
        switch (token.kind) {
        case CifToken::Kind::Comment:
            if (!block) {
                doc.header.push_back(std::move(value));
            }
            continue;
        case CifToken::Kind::Data:
            closeLoop();
            doc.blocks.push_back(CifBlock{std::move(value), {}});
            block = &doc.blocks.back();
            continue;
        default:
            break;
        }
        if (!block) {
            // a monomer-library global_ section before the first data_ block: dropped
            continue;
        }
        if (token.kind == CifToken::Kind::Loop) {
            closeLoop();
            if (tag) {
                throw noValue();
            }
            block->items.emplace_back();
            loop = &block->items.back();
        } else if (token.kind == CifToken::Kind::Tag) {
            if (loop && values.empty()) {
                loop->tags.push_back(std::move(value));
                continue;
            }
            closeLoop();
            if (tag) {
                throw noValue();
            }
            tag = std::move(value);
        } else if (token.kind == CifToken::Kind::Value) {
            if (loop) {
                values.push_back(std::move(value));
            } else if (tag) {
                const std::string category = pythonLower(beforeDot(*tag));
                CifLoop* last = block->items.empty() ? nullptr : &block->items.back();
                if (last && last->pairs && last->category() == category &&
                    last->col(afterDot(*tag)) < 0) {
                    last->tags.push_back(*tag);
                    last->rows[0].push_back(std::move(value));
                } else {
                    CifLoop item;
                    item.tags = {*tag};
                    item.rows = {{std::move(value)}};
                    item.pairs = true;
                    block->items.push_back(std::move(item));
                }
                tag.reset();
            } else {
                throw DictionaryError("mmCIF: value " + core::pythonStringRepr(value) +
                                      " without a tag");
            }
        }
    }
    closeLoop();
    if (tag) {
        throw noValue();
    }
    return doc;
}

// ---- writing ---------------------------------------------------------------------------------

std::string cifValue(std::string_view value) {
    if (value.find('\n') != std::string_view::npos) {
        return "\n;" + normalised(value) + "\n;\n";
    }
    if (value.empty()) {
        return "''";
    }
    const std::u32string s = detail::decodeUtf8Replace(value);
    const bool hasSingle = value.find('\'') != std::string_view::npos;
    const bool hasDouble = value.find('"') != std::string_view::npos;
    const std::string low = pythonLower(value);
    bool needs = std::any_of(s.begin(), s.end(), isPythonSpace) ||
                 std::u32string_view(U"_#$'\"[];").find(s[0]) != std::u32string_view::npos ||
                 hasSingle || hasDouble;
    for (const std::string_view reserved : {"data_", "save_", "loop_", "global_", "stop_"}) {
        needs = needs || low.starts_with(reserved);
    }
    const std::string text = encodeUtf8(s);
    if (!needs) {
        return text;
    }
    if (!hasDouble) {
        return "\"" + text + "\"";
    }
    if (!hasSingle) {
        return "'" + text + "'";
    }
    return "\n;" + text + "\n;\n";
}

std::string writeCif(const CifDocument& doc) {
    std::vector<std::string> out = doc.header;
    for (const CifBlock& block : doc.blocks) {
        if (out.empty() || out.back() != "#") {
            out.emplace_back("#");
        }
        out.push_back("data_" + block.name);
        out.emplace_back("#");
        for (const CifLoop& item : block.items) {
            if (item.tags.empty()) {
                continue;
            }
            if (item.pairs && item.rows.size() == 1) {
                std::size_t width = 0;
                for (const std::string& t : item.tags) {
                    width = std::max(width, core::characterCount(t));
                }
                const std::size_t count = std::min(item.tags.size(), item.rows[0].size());
                for (std::size_t k = 0; k < count; ++k) {
                    out.push_back(core::pythonRstrip(ljust(item.tags[k], width) + "   " +
                                                     cifValue(item.rows[0][k])));
                }
                out.emplace_back("#");
                continue;
            }
            out.emplace_back("loop_");
            out.insert(out.end(), item.tags.begin(), item.tags.end());
            std::vector<std::vector<std::string>> cells;
            std::size_t columns = item.rows.empty() ? 0 : std::string::npos;
            for (const std::vector<std::string>& row : item.rows) {
                std::vector<std::string> c;
                for (const std::string& v : row) {
                    c.push_back(cifValue(v));
                }
                columns = std::min(columns, c.size());
                cells.push_back(std::move(c));
            }
            // zip(*cells): as many columns as the shortest row
            std::vector<std::size_t> widths(columns, 0);
            for (std::size_t k = 0; k < columns; ++k) {
                bool any = false;
                for (const std::vector<std::string>& row : cells) {
                    if (row[k].find('\n') == std::string::npos) {
                        widths[k] = std::max(widths[k], core::characterCount(row[k]));
                        any = true;
                    }
                }
                if (!any) {
                    widths[k] = 1;
                }
            }
            for (const std::vector<std::string>& row : cells) {
                std::string line;
                if (std::any_of(row.begin(), row.end(), [](const std::string& c) {
                        return c.find('\n') != std::string::npos;
                    })) {
                    for (std::size_t k = 0; k < row.size(); ++k) {
                        line += (k ? " " : "") + row[k];
                    }
                    out.push_back(std::move(line));
                } else {
                    for (std::size_t k = 0; k < std::min(row.size(), columns); ++k) {
                        line += (k ? "  " : "") + ljust(row[k], widths[k]);
                    }
                    out.push_back(core::pythonRstrip(line));
                }
            }
            out.emplace_back("#");
        }
    }
    std::string text;
    for (std::size_t k = 0; k < out.size(); ++k) {
        text += (k ? "\n" : "") + out[k];
    }
    return text + "\n";
}

// ---- dictionary files ------------------------------------------------------------------------

namespace {

constexpr std::size_t kDictionaryCharacters = 2'000'000;

// ---- fnmatch.fnmatchcase and glob.glob (POSIX) -----------------------------------------------

struct GlobItem {
    enum class Kind { Star, Any, Literal, Set, Never } kind = Kind::Literal;
    char32_t ch = 0;
    bool negated = false;
    std::vector<std::pair<char32_t, char32_t>> ranges;

    bool matches(char32_t c) const {
        switch (kind) {
        case Kind::Any:
            return true;
        case Kind::Literal:
            return c == ch;
        case Kind::Set: {
            bool in = false;
            for (const auto& [lo, hi] : ranges) {
                in = in || (c >= lo && c <= hi);
            }
            return in != negated;
        }
        case Kind::Never:
        case Kind::Star:
            break;
        }
        return false;
    }
};

std::u32string escapeSetText(const std::u32string& text, bool hyphens) {
    std::u32string out;
    for (const char32_t c : text) {
        if (c == U'\\' || (hyphens && c == U'-')) {
            out += U'\\';
        }
        out += c;
    }
    return out;
}

// The character class Python's re builds from fnmatch's translation of "[...]".
GlobItem setItem(const std::u32string& pat, std::size_t i, std::size_t j) {
    std::u32string stuff;
    const std::u32string inner = pat.substr(i, j - i);
    if (inner.find(U'-') == std::u32string::npos) {
        stuff = escapeSetText(inner, false);
    } else {
        std::vector<std::u32string> chunks;
        std::size_t k = pat[i] == U'!' ? i + 2 : i + 1;
        for (;;) {
            const std::size_t dash = k < j ? pat.find(U'-', k) : std::u32string::npos;
            if (dash == std::u32string::npos || dash >= j) {
                break;
            }
            chunks.push_back(pat.substr(i, dash - i));
            i = dash + 1;
            k = dash + 3;
        }
        const std::u32string chunk = pat.substr(i, j - i);
        if (!chunk.empty()) {
            chunks.push_back(chunk);
        } else {
            chunks.back() += U'-';
        }
        for (std::size_t m = chunks.size() - 1; m > 0; --m) {
            if (chunks[m - 1].empty() || chunks[m].empty()) {  // Python's IndexError (defensive)
                throw core::NamedError("IndexError", "string index out of range");
            }
            if (chunks[m - 1].back() > chunks[m].front()) {
                chunks[m - 1] =
                    chunks[m - 1].substr(0, chunks[m - 1].size() - 1) + chunks[m].substr(1);
                chunks.erase(chunks.begin() + static_cast<std::ptrdiff_t>(m));
            }
        }
        for (std::size_t m = 0; m < chunks.size(); ++m) {
            stuff += (m ? U"-" : U"") + escapeSetText(chunks[m], true);
        }
    }
    GlobItem item;
    if (stuff.empty()) {
        item.kind = GlobItem::Kind::Never;
        return item;
    }
    if (stuff == U"!") {
        item.kind = GlobItem::Kind::Any;
        return item;
    }
    item.kind = GlobItem::Kind::Set;
    std::size_t pos = 0;
    if (stuff[0] == U'!') {
        item.negated = true;
        pos = 1;
    }
    const auto atom = [&](std::size_t& at) {
        char32_t c = stuff[at++];
        if (c == U'\\' && at < stuff.size()) {
            c = stuff[at++];
        }
        return c;
    };
    while (pos < stuff.size()) {
        const char32_t lo = atom(pos);
        if (pos + 1 < stuff.size() && stuff[pos] == U'-') {
            ++pos;
            item.ranges.emplace_back(lo, atom(pos));
        } else {
            item.ranges.emplace_back(lo, lo);
        }
    }
    return item;
}

std::vector<GlobItem> compileGlob(const std::u32string& pat) {
    std::vector<GlobItem> out;
    std::size_t i = 0;
    const std::size_t n = pat.size();
    while (i < n) {
        const char32_t c = pat[i++];
        GlobItem item;
        if (c == U'*') {
            if (out.empty() || out.back().kind != GlobItem::Kind::Star) {
                item.kind = GlobItem::Kind::Star;
                out.push_back(item);
            }
            continue;
        }
        if (c == U'?') {
            item.kind = GlobItem::Kind::Any;
        } else if (c == U'[') {
            std::size_t j = i;
            if (j < n && pat[j] == U'!') {
                ++j;
            }
            if (j < n && pat[j] == U']') {
                ++j;
            }
            while (j < n && pat[j] != U']') {
                ++j;
            }
            if (j >= n) {
                item.ch = U'[';
            } else {
                item = setItem(pat, i, j);
                i = j + 1;
            }
        } else {
            item.ch = c;
        }
        out.push_back(item);
    }
    return out;
}

bool globMatch(const std::vector<GlobItem>& pattern, const std::u32string& text) {
    std::size_t p = 0, t = 0;
    std::size_t starP = std::u32string::npos, starT = 0;
    while (t < text.size()) {
        if (p < pattern.size() && pattern[p].kind == GlobItem::Kind::Star) {
            starP = p++;
            starT = t;
        } else if (p < pattern.size() && pattern[p].matches(text[t])) {
            ++p;
            ++t;
        } else if (starP != std::u32string::npos) {
            p = starP + 1;
            t = ++starT;
        } else {
            return false;
        }
    }
    while (p < pattern.size() && pattern[p].kind == GlobItem::Kind::Star) {
        ++p;
    }
    return p == pattern.size();
}

bool hasMagic(std::string_view text) { return text.find_first_of("*?[") != std::string_view::npos; }

fs::file_status fileStatus(const std::string& path) {
    std::error_code ec;
    return path.empty() ? fs::file_status() : fs::status(core::pathFromUtf8(path), ec);
}

bool isDir(const std::string& path) { return fs::is_directory(fileStatus(path)); }

bool isFile(const std::string& path) { return fs::is_regular_file(fileStatus(path)); }

bool lexists(const std::string& path) {
    std::error_code ec;
    return !path.empty() && fs::exists(fs::symlink_status(core::pathFromUtf8(path), ec));
}

// os.path.split
std::pair<std::string, std::string> splitPath(const std::string& path) {
    return {core::posixpath::dirname(path), core::posixpath::basename(path)};
}

// The names in a folder ("" is the current one), only folders with dirOnly; none on errors.
std::vector<std::string> listDir(const std::string& dirname, bool dirOnly) {
    std::vector<std::string> names;
    std::error_code ec;
    fs::directory_iterator it(core::pathFromUtf8(dirname.empty() ? "." : dirname), ec);
    for (; !ec && it != fs::directory_iterator(); it.increment(ec)) {
        std::error_code dirError;
        if (!dirOnly || it->is_directory(dirError)) {
            names.push_back(core::pathToUtf8(it->path().filename()));
        }
    }
    return names;
}

std::vector<std::string> glob1(const std::string& dirname, const std::string& pattern,
                               bool dirOnly) {
    const std::vector<GlobItem> compiled = compileGlob(core::decodeUtf8Surrogateescape(pattern));
    std::vector<std::string> out;
    for (std::string& name : listDir(dirname, dirOnly)) {
        if (pattern[0] != '.' && !name.empty() && name[0] == '.') {
            continue;
        }
        if (globMatch(compiled, core::decodeUtf8Surrogateescape(name))) {
            out.push_back(std::move(name));
        }
    }
    return out;
}

std::vector<std::string> glob0(const std::string& dirname, const std::string& basename) {
    if (!basename.empty()) {
        if (lexists(dirname.empty() ? basename : core::posixpath::join(dirname, basename))) {
            return {basename};
        }
    } else if (isDir(dirname)) {
        return {basename};
    }
    return {};
}

// glob._iglob(pathname, root_dir="", dir_fd=None, recursive=False, dironly)
void iglob(const std::string& pathname, bool dirOnly, std::vector<std::string>& out) {
    const auto [dirname, basename] = splitPath(pathname);
    if (!hasMagic(pathname)) {
        if (!basename.empty() ? lexists(pathname) : isDir(dirname)) {
            out.push_back(pathname);
        }
        return;
    }
    if (dirname.empty()) {
        for (std::string& name : glob1("", basename, dirOnly)) {
            out.push_back(std::move(name));
        }
        return;
    }
    std::vector<std::string> dirs;
    if (dirname != pathname && hasMagic(dirname)) {
        iglob(dirname, true, dirs);
    } else {
        dirs.push_back(dirname);
    }
    for (const std::string& dir : dirs) {
        const std::vector<std::string> names =
            hasMagic(basename) ? glob1(dir, basename, dirOnly) : glob0(dir, basename);
        for (const std::string& name : names) {
            out.push_back(core::posixpath::join(dir, name));
        }
    }
}

std::vector<std::string> glob(const std::string& pathname) {
    std::vector<std::string> out;
    iglob(pathname, false, out);
    return out;
}

// glob.escape (no drive on POSIX)
std::string globEscape(std::string_view path) {
    std::string out;
    for (const char c : path) {
        if (c == '*' || c == '?' || c == '[') {
            out += '[';
            out += c;
            out += ']';
        } else {
            out += c;
        }
    }
    return out;
}

std::string absolutePath(const std::string& path) {
#if defined(_WIN32)
    return core::pathToUtf8(fs::absolute(core::pathFromUtf8(path)).lexically_normal());
#else
    return core::posixpath::abspath(path, core::currentDirectory());
#endif
}

// os.path.getmtime: st_mtime as Python computes it (sec + nsec * 1e-9).
double modificationTime(const std::string& path) {
#if defined(_WIN32)
    std::error_code ec;
    const auto time = fs::last_write_time(core::pathFromUtf8(path), ec);
    return ec ? 0.0 : static_cast<double>(time.time_since_epoch().count());
#else
    struct stat st {};
    if (::stat(path.c_str(), &st) != 0) {
        return 0.0;
    }
    return static_cast<double>(st.st_mtim.tv_sec) + static_cast<double>(st.st_mtim.tv_nsec) * 1e-9;
#endif
}

bool contains(const std::string& text, std::string_view part) {
    return text.find(part) != std::string::npos;
}

}  // namespace

bool isDictionaryFile(const std::string& path) {
    const fs::path file = core::pathFromUtf8(path);
    std::error_code ec;
    if (path.empty() || fs::is_directory(file, ec)) {
        return false;
    }
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        return false;
    }
    // enough bytes for 2,000,000 characters (at most 4 bytes each, or "\r\n") and one more
    const std::size_t limit = 4 * kDictionaryCharacters + 1;
    std::string bytes;
    std::vector<char> buffer(std::size_t{1} << 16);
    while (bytes.size() < limit) {
        const std::size_t want = std::min(buffer.size(), limit - bytes.size());
        in.read(buffer.data(), static_cast<std::streamsize>(want));
        const auto got = static_cast<std::size_t>(in.gcount());
        bytes.append(buffer.data(), got);
        if (got < want) {
            break;
        }
    }
    if (in.bad()) {
        return false;
    }
    const std::u32string decoded = detail::decodeUtf8Replace(bytes);
    std::u32string text;
    for (std::size_t i = 0; i < decoded.size() && text.size() < kDictionaryCharacters; ++i) {
        if (decoded[i] == U'\r') {
            if (i + 1 < decoded.size() && decoded[i + 1] == U'\n') {
                ++i;
            }
            text += U'\n';
        } else {
            text += decoded[i];
        }
    }
    return text.find(U"_chem_comp_atom.") != std::u32string::npos &&
           text.find(U"_chem_comp_bond.") != std::u32string::npos;
}

std::string findDictionary(std::span<const std::string> sources, std::string_view code,
                           std::span<const std::string> exclude) {
    std::set<std::string> skip;
    for (const std::string& p : exclude) {
        skip.insert(absolutePath(p));
    }
    std::vector<std::string> candidates;
    for (const std::string& src : sources) {
        std::vector<std::string> found;
        if (isDir(src)) {
            found = glob(core::posixpath::join(globEscape(src), "*.cif"));
        } else if (hasMagic(src)) {
            found = glob(src);
        } else if (isFile(src)) {
            found.push_back(src);
        }
        for (std::string& p : found) {
            if (isFile(p) && !skip.contains(absolutePath(p))) {
                candidates.push_back(std::move(p));
            }
        }
    }
    // sorted(set(...)): str order is code point order
    const auto codePoints = [](const std::string& s) { return core::decodeUtf8Surrogateescape(s); };
    std::sort(candidates.begin(), candidates.end(),
              [&](const std::string& a, const std::string& b) {
                  return codePoints(a) < codePoints(b);
              });
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    struct Candidate {
        bool otherName;
        bool tmp;
        double mtime;
        std::u32string path;
        std::string text;
    };
    std::vector<Candidate> good;
    const std::string lowCode = pythonLower(code);
    for (const std::string& p : candidates) {
        if (isDictionaryFile(p)) {
            const std::string base = pythonLower(core::posixpath::basename(p));
            good.push_back({!code.empty() && !contains(base, lowCode),
                            contains(pythonLower(p), "tmp"), modificationTime(p), codePoints(p),
                            p});
        }
    }
    if (good.empty()) {
        std::string joined;
        for (std::size_t k = 0; k < sources.size(); ++k) {
            joined += (k ? ", " : "") + sources[k];
        }
        throw DictionaryError(
            "no restraint dictionary (.cif with _chem_comp_atom and _chem_comp_bond loops) in " +
            joined);
    }
    // key: (code not in the name, "tmp" in the path, -mtime, path)
    const auto best =
        std::min_element(good.begin(), good.end(), [](const Candidate& a, const Candidate& b) {
            return std::tie(a.otherName, a.tmp, b.mtime, a.path) <
                   std::tie(b.otherName, b.tmp, a.mtime, b.path);
        });
    return best->text;
}

}  // namespace mdbridge::io
