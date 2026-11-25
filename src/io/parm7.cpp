#include "mdbridge/io/parm7.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <system_error>

#include "mdbridge/core/detail/data_tables.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/warnings.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace {

using core::ValueError;
using detail::Lines;

std::string utf8(std::u32string_view text) { return core::encodeUtf8(text); }

bool startsWith(std::u32string_view text, std::u32string_view prefix) {
    return text.substr(0, prefix.size()) == prefix;
}

std::u32string strip(std::u32string_view text) {
    std::size_t first = 0;
    std::size_t last = text.size();
    while (first < last && core::isPythonSpace(text[first])) {
        ++first;
    }
    while (last > first && core::isPythonSpace(text[last - 1])) {
        --last;
    }
    return std::u32string(text.substr(first, last - first));
}

// str.split() without arguments
std::vector<std::u32string> splitWhitespace(std::u32string_view text) {
    std::vector<std::u32string> out;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && core::isPythonSpace(text[i])) {
            ++i;
        }
        const std::size_t start = i;
        while (i < text.size() && !core::isPythonSpace(text[i])) {
            ++i;
        }
        if (i > start) {
            out.emplace_back(text.substr(start, i - start));
        }
    }
    return out;
}

std::string repr(const std::string& text) { return core::pythonStringRepr(text); }

// ---------------------------------------------------------------- MDATK's section reader

struct Format {
    char kind = 'a';
    std::size_t width = 80;
};

// _FMT.search(line): \(\s*(\d+)\s*([aAiIeEfF])\s*(\d+)(?:\.(\d+))?\s*\)
std::optional<Format> searchFormat(std::u32string_view line) {
    const auto spaces = [&](std::size_t& i) {
        while (i < line.size() && core::isPythonSpace(line[i])) {
            ++i;
        }
    };
    const auto digits = [&](std::size_t& i) {
        const std::size_t start = i;
        while (i < line.size() && core::isPythonDigit(line[i])) {
            ++i;
        }
        return line.substr(start, i - start);
    };
    for (std::size_t start = 0; start < line.size(); ++start) {
        if (line[start] != U'(') {
            continue;
        }
        std::size_t i = start + 1;
        spaces(i);
        if (digits(i).empty()) {
            continue;
        }
        spaces(i);
        if (i >= line.size() ||
            std::u32string_view(U"aAiIeEfF").find(line[i]) == std::u32string_view::npos) {
            continue;
        }
        const char kind = static_cast<char>(line[i] | 0x20);
        ++i;
        spaces(i);
        const std::u32string_view width = digits(i);
        if (width.empty()) {
            continue;
        }
        if (i + 1 < line.size() && line[i] == U'.' && core::isPythonDigit(line[i + 1])) {
            ++i;
            digits(i);
        }
        spaces(i);
        if (i >= line.size() || line[i] != U')') {
            continue;
        }
        Format format;
        format.kind = kind;
        try {
            format.width = static_cast<std::size_t>(core::pythonInt(utf8(width)));
        } catch (const std::overflow_error&) {
            format.width = std::numeric_limits<std::size_t>::max();
        }
        return format;
    }
    return std::nullopt;
}

struct RawSection {
    std::u32string name;
    Format format;
    std::vector<std::u32string> lines;
};

Parm7Values sectionValues(const RawSection& section) {
    std::vector<std::u32string> raw;
    for (std::u32string_view line : section.lines) {
        while (!line.empty() && line.back() == U'\n') {
            line.remove_suffix(1);
        }
        if (section.format.width == 0) {
            throw ValueError("range() arg 3 must not be zero");
        }
        for (std::size_t k = 0; k < line.size(); k += section.format.width) {
            const std::u32string_view chunk = line.substr(k, section.format.width);
            if (section.format.kind == 'a' || !strip(chunk).empty()) {
                raw.emplace_back(chunk);
            }
            if (section.format.width > line.size()) {
                break;
            }
        }
    }
    if (section.format.kind == 'a') {
        std::vector<std::string> out;
        for (const std::u32string& s : raw) {
            out.push_back(utf8(strip(s)));
        }
        return out;
    }
    if (section.format.kind == 'i') {
        std::vector<std::int64_t> out;
        for (const std::u32string& s : raw) {
            out.push_back(core::pythonInt(utf8(s)));
        }
        return out;
    }
    std::vector<double> out;
    for (std::u32string s : raw) {
        for (char32_t& c : s) {
            c = c == U'D' ? U'E' : c == U'd' ? U'e' : c;
        }
        out.push_back(core::pythonFloat(utf8(s)));
    }
    return out;
}

// ---------------------------------------------------------------- conversions of written values

std::string pyStr(const Parm7Values& values, std::size_t k) {
    if (const auto* s = std::get_if<std::vector<std::string>>(&values)) {
        return (*s)[k];
    }
    if (const auto* i = std::get_if<std::vector<std::int64_t>>(&values)) {
        return std::to_string((*i)[k]);
    }
    return core::pythonFloatRepr(std::get<std::vector<double>>(values)[k]);
}

// str(int(v))
std::string pyIntText(const Parm7Values& values, std::size_t k) {
    if (const auto* s = std::get_if<std::vector<std::string>>(&values)) {
        return std::to_string(core::pythonInt((*s)[k]));
    }
    if (const auto* i = std::get_if<std::vector<std::int64_t>>(&values)) {
        return std::to_string((*i)[k]);
    }
    const double x = std::get<std::vector<double>>(values)[k];
    if (std::isnan(x)) {
        throw ValueError("cannot convert float NaN to integer");
    }
    if (std::isinf(x)) {
        throw std::overflow_error("cannot convert float infinity to integer");
    }
    const double t = std::trunc(x);
    if (t == 0.0) {
        return "0";
    }
    char buffer[400];
    const auto result =
        std::to_chars(buffer, buffer + sizeof buffer, t, std::chars_format::fixed, 0);
    return std::string(buffer, result.ptr);
}

double pyFloat(const Parm7Values& values, std::size_t k) {
    if (const auto* s = std::get_if<std::vector<std::string>>(&values)) {
        return core::pythonFloat((*s)[k]);
    }
    if (const auto* i = std::get_if<std::vector<std::int64_t>>(&values)) {
        return static_cast<double>((*i)[k]);
    }
    return std::get<std::vector<double>>(values)[k];
}

std::size_t valueCount(const Parm7Values& values) {
    return std::visit([](const auto& v) { return v.size(); }, values);
}

// f"{text[:n]:<n}" (code points)
std::string leftField(const std::string& text, std::size_t n) {
    std::u32string s = core::decodeUtf8(text);
    if (s.size() > n) {
        s.resize(n);
    }
    std::string out = utf8(s);
    out.append(n - s.size(), ' ');
    return out;
}

std::string rightField(const std::string& text, std::size_t n) {
    return text.size() >= n ? text : std::string(n - text.size(), ' ') + text;
}

// f"{x:16.8E}"
std::string formatE(double x) {
    std::string text;
    if (std::isnan(x)) {
        text = "NAN";  // Python drops the sign of a NaN
    } else if (std::isinf(x)) {
        text = x < 0 ? "-INF" : "INF";
    } else {
        char buffer[64];
        const auto result =
            std::to_chars(buffer, buffer + sizeof buffer, x, std::chars_format::scientific, 8);
        text.assign(buffer, result.ptr);
        for (char& c : text) {
            c = c == 'e' ? 'E' : c;
        }
    }
    return rightField(text, 16);
}

// What Python's strict ASCII encoder reports for the first run of non-ASCII characters.
void checkAscii(const std::string& text) {
    const std::u32string s = core::decodeUtf8(text);
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] < 0x80) {
            continue;
        }
        std::size_t end = i;
        while (end < s.size() && s[end] >= 0x80) {
            ++end;
        }
        std::string message = "'ascii' codec can't encode ";
        if (end == i + 1) {
            char escape[16];
            const auto c = static_cast<unsigned long>(s[i]);
            std::snprintf(escape, sizeof escape,
                          c <= 0xFF     ? "\\x%02lx"
                          : c <= 0xFFFF ? "\\u%04lx"
                                        : "\\U%08lx",
                          c);
            message += "character '" + std::string(escape) + "' in position " + std::to_string(i);
        } else {
            message +=
                "characters in position " + std::to_string(i) + "-" + std::to_string(end - 1);
        }
        throw ValueError(message + ": ordinal not in range(128)");
    }
}

// ---------------------------------------------------------------- topology (MDAnalysis TOPParser)
//
// One pass over the lines: the two header lines, the POINTERS counts, then the %FLAG sections.
// The sections the topology uses hold the number of values their POINTERS count says; their data
// lines are counted from the fields per line of their %FORMAT (RD-P06). Unknown sections are
// skipped wherever they are (RD-P01). Every error names the file and, where one is known, the
// line, the section and the field (docs/DIFFERENCES.md, "parm7 topologies").

// Element symbols by atomic number: kIndElements[Z - 1] (mdbridge/core/detail/data_tables.hpp).
namespace tables = core::detail::tables;

// How the field texts of a section become values.
enum class FieldKind {
    Text,     // the stripped text
    Charge,   // x / 18.2223, in float64 (RD-P12)
    Real,     // double
    Element,  // atomic number -> element symbol
    Integer,  // int
    Pointer,  // int - 1
    Index,    // int, checked and divided by 3 when the topology is built
};

struct SectionRule {
    std::u32string_view flag;
    std::string_view key;
    std::size_t pointer;  // the POINTERS entry with the record count
    std::int64_t perRecord;
    FieldKind kind;
};

constexpr std::array<SectionRule, 15> kSectionRules = {{
    {U"ATOM_NAME", "name", 0, 1, FieldKind::Text},
    {U"CHARGE", "charge", 0, 1, FieldKind::Charge},
    {U"ATOMIC_NUMBER", "elements", 0, 1, FieldKind::Element},
    {U"MASS", "mass", 0, 1, FieldKind::Real},
    {U"ATOM_TYPE_INDEX", "type_indices", 0, 1, FieldKind::Integer},
    {U"AMBER_ATOM_TYPE", "types", 0, 1, FieldKind::Text},
    {U"RESIDUE_LABEL", "resname", 11, 1, FieldKind::Text},
    {U"RESIDUE_POINTER", "respoint", 11, 1, FieldKind::Pointer},
    {U"BONDS_INC_HYDROGEN", "bondh", 2, 3, FieldKind::Index},
    {U"BONDS_WITHOUT_HYDROGEN", "bonda", 3, 3, FieldKind::Index},
    {U"ANGLES_INC_HYDROGEN", "angh", 4, 4, FieldKind::Index},
    {U"ANGLES_WITHOUT_HYDROGEN", "anga", 5, 4, FieldKind::Index},
    {U"DIHEDRALS_INC_HYDROGEN", "dihh", 6, 5, FieldKind::Index},
    {U"DIHEDRALS_WITHOUT_HYDROGEN", "diha", 7, 5, FieldKind::Index},
    {U"RESIDUE_CHAINID", "chainid", 11, 1, FieldKind::Text},
}};

const SectionRule* findRule(std::u32string_view flag) {
    for (const SectionRule& rule : kSectionRules) {
        if (rule.flag == flag) {
            return &rule;
        }
    }
    return nullptr;
}

const SectionRule& ruleOfKey(std::string_view key) {
    for (const SectionRule& rule : kSectionRules) {
        if (rule.key == key) {
            return rule;
        }
    }
    throw std::logic_error("parm7: unknown section key");
}

std::string flagName(const SectionRule& rule) { return "%FLAG " + utf8(rule.flag); }

// The values one section gave: texts (Text, Element), integers or reals.
struct Collected {
    std::vector<std::string> texts;
    std::vector<std::int64_t> integers;
    std::vector<double> reals;

    std::size_t size() const { return texts.size() + integers.size() + reals.size(); }
};

// One item of a Fortran format: `repeat` fields of `width` columns, or `repeat` skipped columns.
struct FormatItem {
    std::size_t repeat = 1;
    std::size_t width = 0;
    bool skip = false;
};

bool isAsciiDigit(char32_t c) { return c >= U'0' && c <= U'9'; }

std::size_t parseCount(std::u32string_view digits) {
    try {
        return static_cast<std::size_t>(core::pythonInt(utf8(digits)));
    } catch (const std::overflow_error&) {
        return std::numeric_limits<std::size_t>::max();
    }
}

// The items of a %FORMAT text: comma-separated [repeat]letter[width][.decimals] with a letter of
// I F E A X in either case; only the start of an item has to have this shape. Empty: not a format.
std::optional<std::vector<FormatItem>> parseFortranFormat(std::u32string_view text) {
    std::vector<FormatItem> items;
    std::size_t begin = 0;
    while (true) {
        const std::size_t comma = text.find(U',', begin);
        const std::u32string_view item = comma == std::u32string_view::npos
                                             ? text.substr(begin)
                                             : text.substr(begin, comma - begin);
        std::size_t i = 0;
        while (i < item.size() && isAsciiDigit(item[i])) {
            ++i;
        }
        const std::u32string_view repeat = item.substr(0, i);
        char32_t letter = i < item.size() ? item[i] : U' ';
        if (letter >= U'A' && letter <= U'Z') {
            letter = static_cast<char32_t>(letter - U'A' + U'a');
        }
        if (std::u32string_view(U"ifeax").find(letter) == std::u32string_view::npos) {
            return std::nullopt;
        }
        ++i;
        const std::size_t widthStart = i;
        while (i < item.size() && isAsciiDigit(item[i])) {
            ++i;
        }
        FormatItem parsed;
        parsed.repeat = repeat.empty() ? 1 : parseCount(repeat);
        if (letter == U'x') {
            parsed.skip = true;
        } else {
            const std::u32string_view width = item.substr(widthStart, i - widthStart);
            parsed.width = width.empty() ? 0 : parseCount(width);
        }
        items.push_back(parsed);
        if (comma == std::u32string_view::npos) {
            break;
        }
        begin = comma + 1;
    }
    return items;
}

// The data fields of one line of a format (saturating at SIZE_MAX).
std::size_t fieldsPerLine(const std::vector<FormatItem>& items) {
    std::size_t total = 0;
    for (const FormatItem& item : items) {
        if (!item.skip && item.width > 0) {
            total = item.repeat > std::numeric_limits<std::size_t>::max() - total
                        ? std::numeric_limits<std::size_t>::max()
                        : total + item.repeat;
        }
    }
    return total;
}

// Moves `column` on by `step`, at most to `end`.
std::size_t advance(std::size_t column, std::size_t step, std::size_t end) {
    return step > end - column ? end : column + step;
}

// The stripped texts of the fields of one line laid out by `items` (blank fields as "") up to the
// end of the line.
std::vector<std::u32string> splitFields(std::u32string_view line,
                                        const std::vector<FormatItem>& items) {
    std::vector<std::u32string> out;
    std::size_t column = 0;
    for (const FormatItem& item : items) {
        if (column >= line.size()) {
            break;
        }
        if (item.skip) {
            column = advance(column, item.repeat, line.size());
            continue;
        }
        if (item.width == 0) {
            continue;
        }
        for (std::size_t r = 0; r < item.repeat && column < line.size(); ++r) {
            out.push_back(strip(line.substr(column, item.width)));
            column = advance(column, item.width, line.size());
        }
    }
    return out;
}

// ceil(a / b) for a >= 0, b > 0
std::int64_t ceilDiv(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b;
    if (a % b != 0) {
        ++q;
    }
    return q;
}

constexpr double kAmberChargeFactor = 18.2223;

class Parm7TopologyReader {
public:
    explicit Parm7TopologyReader(const std::filesystem::path& path)
        : m_file(detail::pathText(path)), m_lines(detail::readTextLines(path)) {}

    std::shared_ptr<core::Topology> read() {
        readHeader();
        const std::u32string* line = readPointers();
        walkSections(line);
        return build();
    }

private:
    std::string m_file;
    Lines m_lines;
    std::size_t m_next = 0;
    std::vector<std::int64_t> m_pointers;
    std::map<std::string, Collected, std::less<>> m_sections;

    [[noreturn]] void fail(const std::string& message) const {
        throw ValueError(m_file + ": " + message);
    }
    [[noreturn]] void failAt(std::size_t line, const std::string& message) const {
        throw ValueError(m_file + " line " + std::to_string(line) + ": " + message);
    }

    // The 1-based number of the line nextLine() returned last.
    std::size_t lineNumber() const { return m_next; }

    const std::u32string* nextLine() {
        return m_next < m_lines.size() ? &m_lines[m_next++] : nullptr;
    }

    void readHeader() {
        if (m_lines.empty()) {
            fail("file is empty");
        }
        if (!startsWith(*nextLine(), U"%VE")) {
            throw ValueError(m_file + " is not a valid TOP file. %VE Missing in header");
        }
        const std::u32string* second = nextLine();
        const std::vector<std::u32string> tokens =
            second ? splitWhitespace(*second) : std::vector<std::u32string>();
        if (tokens.size() >= 2 && tokens[1] == U"CTITLE") {
            throw ValueError(m_file +
                             " is detected as a Chamber-style TOP file. At this time MDAnalysis "
                             "does not support such topologies");
        }
        if (tokens.size() < 2 || tokens[1] != U"TITLE") {
            throw ValueError(m_file + " is not a valid TOP file. 'TITLE' missing in header");
        }
    }

    // Reads %FLAG POINTERS (anywhere after the header) up to the next %FLAG line, which it
    // returns (nullptr at the end of the file).
    const std::u32string* readPointers() {
        while (true) {
            const std::u32string* line = nextLine();
            if (!line) {
                fail("required section %FLAG POINTERS is missing");
            }
            if (startsWith(*line, U"%FLAG POINTERS")) {
                break;
            }
        }
        readSectionFormat(U"POINTERS");
        while (const std::u32string* line = nextLine()) {
            if (startsWith(*line, U"%FLAG")) {
                checkPointerCount();
                return line;
            }
            for (const std::u32string& token : splitWhitespace(*line)) {
                const std::optional<std::int64_t> value = core::tryPythonInt(utf8(token));
                if (!value) {
                    failAt(lineNumber(), "%FLAG POINTERS value " + repr(utf8(token)) +
                                             " is not an integer");
                }
                m_pointers.push_back(*value);
            }
        }
        checkPointerCount();
        return nullptr;
    }

    void checkPointerCount() const {
        if (m_pointers.size() < 12) {
            fail("%FLAG POINTERS needs at least 12 integers (found " +
                 std::to_string(m_pointers.size()) + ")");
        }
    }

    // The text after the first %FLAG (at `at`) up to any later %FLAG, stripped.
    static std::u32string sectionName(std::u32string_view line, std::size_t at) {
        std::u32string_view rest = line.substr(at + 5);
        const std::size_t again = rest.find(U"%FLAG");
        if (again != std::u32string_view::npos) {
            rest = rest.substr(0, again);
        }
        return strip(rest);
    }

    void walkSections(const std::u32string* line) {
        // the known section just read: its flag and data line count (a blank line may follow an
        // empty section)
        std::optional<std::pair<const SectionRule*, std::int64_t>> previous;
        while (line) {
            if (previous && previous->second == 0 && strip(*line).empty()) {
                line = nextLine();
                if (!line) {
                    return;
                }
            }
            const std::size_t at = line->find(U"%FLAG");
            if (at == std::u32string::npos) {
                std::string after = "the data lines";
                if (previous) {
                    after = "the " + std::to_string(previous->second) + " data lines of " +
                            flagName(*previous->first);
                }
                failAt(lineNumber(), "expected a %FLAG line after " + after + ", found " +
                                         repr(utf8(strip(*line))));
            }
            const SectionRule* rule = findRule(sectionName(*line, at));
            if (!rule) {
                line = nextLine();
                while (line && !startsWith(*line, U"%FLAG")) {
                    line = nextLine();
                }
                previous.reset();
                continue;
            }
            previous = std::make_pair(rule, readKnownSection(*rule));
            line = nextLine();
        }
    }

    std::vector<FormatItem> readSectionFormat(std::u32string_view flag) {
        while (true) {
            const std::u32string* line = nextLine();
            if (!line) {
                fail("unexpected end of file in %FLAG " + utf8(flag) +
                     " (expected a %FORMAT line)");
            }
            if (startsWith(*line, U"%COMMENT")) {
                continue;
            }
            if (startsWith(*line, U"%FORMAT")) {
                std::u32string_view text = std::u32string_view(*line).substr(
                    std::min<std::size_t>(8, line->size()));
                const std::size_t close = text.find(U')');
                if (close != std::u32string_view::npos) {
                    text = text.substr(0, close);
                }
                const auto items = parseFortranFormat(text);
                if (!items) {
                    failAt(lineNumber(), "cannot parse the format of %FLAG " + utf8(flag) + ": " +
                                             repr(utf8(strip(*line))));
                }
                return *items;
            }
            failAt(lineNumber(), "expected %FORMAT or %COMMENT after %FLAG " + utf8(flag) +
                                     ", found " + repr(utf8(strip(*line))));
        }
    }

    std::int64_t pointer(std::size_t index) const { return m_pointers[index]; }

    void convertField(const SectionRule& rule, const std::u32string& field, std::size_t line,
                      std::size_t column, Collected& into) const {
        const std::string text = utf8(field);
        const auto where = [&] {
            return flagName(rule) + " field " + std::to_string(column) + ": " + repr(text);
        };
        const auto bad = [&](const std::string& what) {
            throw ValueError(m_file + " line " + std::to_string(line) + ", " + where() + what);
        };
        const auto integer = [&] {
            std::optional<std::int64_t> value;
            try {
                value = core::tryPythonInt(text);
            } catch (const std::overflow_error&) {
                bad(" is too large");
            }
            if (!value) {
                bad(" is not an integer");
            }
            return *value;
        };
        const auto real = [&] {
            std::u32string s = field;
            for (char32_t& c : s) {  // Fortran D exponents (RD-P07)
                c = c == U'D' ? U'E' : c == U'd' ? U'e' : c;
            }
            const std::optional<double> value = core::tryPythonFloat(utf8(s));
            if (!value) {
                bad(" is not a number");
            }
            return *value;
        };
        switch (rule.kind) {
        case FieldKind::Text: into.texts.push_back(text); break;
        case FieldKind::Charge: into.reals.push_back(real() / kAmberChargeFactor); break;
        case FieldKind::Real: into.reals.push_back(real()); break;
        case FieldKind::Element: {
            const std::int64_t z = integer();
            if (z > static_cast<std::int64_t>(std::size(tables::kIndElements))) {
                fail("ATOMIC_NUMBER " + std::to_string(z) + " of atom " +
                     std::to_string(into.texts.size() + 1) + " is not an element");
            }
            into.texts.emplace_back(z <= 0 ? std::string_view()
                                           : tables::kIndElements[static_cast<std::size_t>(z - 1)]
                                                 .symbol);
            break;
        }
        case FieldKind::Integer:
        case FieldKind::Index: into.integers.push_back(integer()); break;
        case FieldKind::Pointer: into.integers.push_back(integer() - 1); break;
        }
    }

    // Reads the header and data lines of a known section; returns its data line count.
    std::int64_t readKnownSection(const SectionRule& rule) {
        const std::size_t flagLine = lineNumber();
        const std::vector<FormatItem> format = readSectionFormat(rule.flag);
        const std::size_t perLine = fieldsPerLine(format);
        const std::int64_t records = pointer(rule.pointer);
        if (records < 0) {
            fail("%FLAG POINTERS gives " + flagName(rule) + " a negative count (" +
                 std::to_string(records) + ")");
        }
        if (records > std::numeric_limits<std::int64_t>::max() / rule.perRecord) {
            fail("%FLAG POINTERS gives " + flagName(rule) + " too large a count (" +
                 std::to_string(records) + ")");
        }
        const std::int64_t count = records * rule.perRecord;
        if (perLine == 0 && count > 0) {
            failAt(flagLine, "the %FORMAT of " + flagName(rule) + " has no data fields");
        }
        constexpr auto kMaxPerLine = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max());
        const auto perLine64 = static_cast<std::int64_t>(std::min(perLine, kMaxPerLine));
        const std::int64_t lineCount = count == 0 ? 0 : ceilDiv(count, perLine64);
        Collected collected;
        const bool text = rule.kind == FieldKind::Text;
        for (std::int64_t k = 0; k < lineCount; ++k) {
            const std::u32string* line = nextLine();
            if (!line) {
                fail("unexpected end of file in " + flagName(rule) + " (expected " +
                     std::to_string(lineCount) + " data lines, found " + std::to_string(k) + ")");
            }
            const std::vector<std::u32string> fields = splitFields(*line, format);
            std::size_t used = fields.size();  // up to the last non-blank field
            while (used > 0 && fields[used - 1].empty()) {
                --used;
            }
            for (std::size_t f = 0; f < used; ++f) {
                if (fields[f].empty() && !text) {  // RD-P05
                    failAt(lineNumber(), "blank field " + std::to_string(f + 1) + " in " +
                                             flagName(rule));
                }
                convertField(rule, fields[f], lineNumber(), f + 1, collected);
            }
        }
        m_sections[std::string(rule.key)] = std::move(collected);
        return lineCount;
    }

    const Collected* section(std::string_view key) const {
        const auto it = m_sections.find(key);
        return it == m_sections.end() ? nullptr : &it->second;
    }

    void expectCount(std::string_view key, std::size_t expected, const char* what) const {
        const std::size_t have = section(key)->size();
        if (have != expected) {
            fail(flagName(ruleOfKey(key)) + " has " + std::to_string(have) + " values for " +
                 std::to_string(expected) + " " + what);
        }
    }

    // The atom indices of a bond / angle / dihedral section, checked (RD-P10): each entry has
    // `n` values, the first n - 1 of them atom indices 3 * (atom - 1) (dihedrals: the third and
    // fourth may be negative, as Amber flags them).
    std::vector<std::vector<std::int32_t>> recordTuples(std::string_view key,
                                                        std::size_t atomCount) const {
        const SectionRule& rule = ruleOfKey(key);
        const std::vector<std::int64_t>& v = section(key)->integers;
        const auto n = static_cast<std::size_t>(rule.perRecord);
        if (v.size() % n != 0) {
            fail(flagName(rule) + " entry " + std::to_string(v.size() / n + 1) + " has " +
                 std::to_string(v.size() % n) + " of " + std::to_string(n) + " values");
        }
        std::vector<std::vector<std::int32_t>> out;
        for (std::size_t k = 0; k < v.size(); k += n) {
            std::vector<std::int32_t> atoms;
            for (std::size_t i = 0; i + 1 < n; ++i) {
                std::int64_t value = v[k + i];
                const std::string entry =
                    flagName(rule) + " entry " + std::to_string(k / n + 1) + ": atom index " +
                    std::to_string(value);
                if (n == 5 && i >= 2 && value < 0 && value != std::numeric_limits<std::int64_t>::min()) {
                    value = -value;
                }
                if (value < 0) {
                    fail(entry + " is negative");
                }
                if (value % 3 != 0) {
                    fail(entry + " is not a multiple of 3");
                }
                if (static_cast<std::uint64_t>(value / 3) >= atomCount) {
                    fail(entry + " refers to atom " + std::to_string(value / 3 + 1) + " of " +
                         std::to_string(atomCount));
                }
                atoms.push_back(static_cast<std::int32_t>(value / 3));
            }
            out.push_back(std::move(atoms));
        }
        return out;
    }

    std::shared_ptr<core::Topology> build() const {
        static constexpr std::array<std::string_view, 9> kRequired = {
            "name", "respoint", "resname", "bonda", "bondh", "anga", "angh", "diha", "dihh"};
        for (const std::string_view required : kRequired) {
            if (!section(required)) {
                fail("required section " + flagName(ruleOfKey(required)) + " is missing");
            }
        }
        const auto atomCount = static_cast<std::size_t>(pointer(0));
        const auto residueCount = static_cast<std::size_t>(pointer(11));
        // per-atom and per-residue sections against the POINTERS counts (RD-P14)
        for (const std::string_view key : {"name", "charge", "elements", "mass", "types",
                                           "type_indices"}) {
            if (section(key)) {
                expectCount(key, atomCount, "atoms");
            }
        }
        expectCount("resname", residueCount, "residues");
        expectCount("respoint", residueCount, "residues");

        // residues: names from RESIDUE_LABEL, membership from RESIDUE_POINTER (RD-P08)
        const std::vector<std::int64_t>& starts = section("respoint")->integers;
        for (std::size_t r = 0; r < starts.size(); ++r) {
            const std::string at = "%FLAG RESIDUE_POINTER residue " + std::to_string(r + 1) +
                                   " starts at atom " + std::to_string(starts[r] + 1);
            if (r == 0 && starts[0] != 0) {
                fail(at + "; the first residue must start at atom 1");
            }
            if (r > 0 && starts[r] <= starts[r - 1]) {
                fail("%FLAG RESIDUE_POINTER is not increasing at residue " +
                     std::to_string(r + 1) + " (" + std::to_string(starts[r - 1] + 1) + ", then " +
                     std::to_string(starts[r] + 1) + ")");
            }
            if (starts[r] >= static_cast<std::int64_t>(atomCount)) {
                fail(at + " but there are " + std::to_string(atomCount) + " atoms");
            }
        }

        // bonds, angles and dihedrals: the topology keeps the bonds; the others are only checked
        std::vector<std::vector<std::int32_t>> bondTuples = recordTuples("bonda", atomCount);
        for (std::vector<std::int32_t>& t : recordTuples("bondh", atomCount)) {
            bondTuples.push_back(std::move(t));
        }
        for (const std::string_view key : {"anga", "angh", "diha", "dihh"}) {
            (void)recordTuples(key, atomCount);
        }

        auto topology = std::make_shared<core::Topology>();
        core::Topology& top = *topology;
        for (const std::string& name : section("name")->texts) {
            top.names.push_back(top.strings.intern(name));
        }
        if (const Collected* charges = section("charge")) {
            top.charges = charges->reals;
        }
        if (const Collected* elements = section("elements")) {
            std::vector<std::uint32_t> symbols;
            for (const std::string& symbol : elements->texts) {
                symbols.push_back(top.strings.intern(symbol));
            }
            top.elements = std::move(symbols);
        }
        if (const Collected* masses = section("mass")) {
            top.masses = masses->reals;
        }
        if (const Collected* types = section("types")) {
            for (const std::string& type : types->texts) {
                top.types.push_back(top.strings.intern(type));
            }
        }
        std::vector<std::int64_t> serials(atomCount);
        for (std::size_t i = 0; i < atomCount; ++i) {
            serials[i] = static_cast<std::int64_t>(i) + 1;
        }
        top.ids = std::move(serials);

        std::vector<core::Bond> bonds;
        for (const std::vector<std::int32_t>& t : bondTuples) {
            if (t[0] != t[1]) {  // MDBridge's bonds join two different atoms
                bonds.push_back({std::min(t[0], t[1]), std::max(t[0], t[1])});
            }
        }
        top.bonds = std::move(bonds);

        const std::vector<std::string>& resnames = section("resname")->texts;
        for (std::size_t r = 0; r < residueCount; ++r) {
            top.resnames.push_back(top.strings.intern(resnames[r]));
            top.resids.push_back(static_cast<std::int64_t>(r) + 1);
        }
        top.resnums = top.resids;
        top.atomResidue.assign(atomCount, 0);
        for (std::size_t r = 0; r < residueCount; ++r) {
            const auto last = r + 1 < residueCount ? static_cast<std::size_t>(starts[r + 1]) : atomCount;
            for (auto a = static_cast<std::size_t>(starts[r]); a < last; ++a) {
                top.atomResidue[a] = static_cast<std::int32_t>(r);
            }
        }

        // segments: runs of equal RESIDUE_CHAINID values, or a single segment SYSTEM
        const Collected* chains = section("chainid");
        if (chains && chains->texts.size() == residueCount) {
            std::vector<std::uint32_t> perAtom;
            for (const std::int32_t r : top.atomResidue) {
                perAtom.push_back(top.strings.intern(chains->texts[static_cast<std::size_t>(r)]));
            }
            top.chainIDs = std::move(perAtom);
            for (std::size_t r = 0; r < residueCount; ++r) {
                if (r == 0 || chains->texts[r] != chains->texts[r - 1]) {
                    top.segids.push_back(top.strings.intern(chains->texts[r]));
                }
                top.residueSegment.push_back(static_cast<std::int32_t>(top.segids.size() - 1));
            }
        } else {
            if (chains) {  // LD-03 (MDAnalysis: "... Skipping section", a Python warning MDATK hides)
                readerWarning(m_file + ": %FLAG RESIDUE_CHAINID has " + std::to_string(chains->texts.size()) +
                              " values for " + std::to_string(residueCount) + " residues; chain IDs are not read");
            }
            top.segids.push_back(top.strings.intern("SYSTEM"));
            top.residueSegment.assign(residueCount, 0);
        }
        top.finish();
        return topology;
    }
};

}  // namespace

Parm7Sections readParm7Sections(const std::filesystem::path& path,
                                const std::optional<std::set<std::string>>& wanted) {
    const std::string pathText = detail::pathText(path);
    const Lines lines = detail::readTextLines(path);
    Parm7Sections out;
    std::optional<RawSection> current;
    const auto isWanted = [&](const std::u32string& name) {
        return !wanted || wanted->count(utf8(name)) > 0;
    };
    const auto flush = [&] {
        if (current && isWanted(current->name)) {
            Parm7Values values = sectionValues(*current);
            const std::string name = utf8(current->name);
            for (auto& [flag, v] : out) {
                if (flag == name) {
                    v = std::move(values);
                    return;
                }
            }
            out.emplace_back(name, std::move(values));
        }
    };
    for (const std::u32string& line : lines) {
        if (startsWith(line, U"%FLAG")) {
            flush();
            current = RawSection{strip(std::u32string_view(line).substr(5)), Format{}, {}};
            continue;
        }
        if (startsWith(line, U"%FORMAT")) {
            if (!current) {
                continue;
            }
            const std::optional<Format> format = searchFormat(line);
            if (!format) {
                throw ValueError(pathText + ": cannot parse format line " +
                                 repr(utf8(strip(line))));
            }
            current->format = *format;
            continue;
        }
        if (startsWith(line, U"%VERSION") || startsWith(line, U"%COMMENT")) {
            continue;
        }
        if (current && isWanted(current->name)) {
            current->lines.push_back(line);
        }
    }
    flush();
    return out;
}

const Parm7Values* findSection(const Parm7Sections& sections, std::string_view flag) {
    for (const auto& [name, values] : sections) {
        if (name == flag) {
            return &values;
        }
    }
    return nullptr;
}

namespace {

std::vector<std::int64_t> asIntegers(const Parm7Values& values) {
    std::vector<std::int64_t> out;
    for (std::size_t k = 0; k < valueCount(values); ++k) {
        out.push_back(std::stoll(pyIntText(values, k)));
    }
    return out;
}

std::vector<double> asReals(const Parm7Values& values) {
    std::vector<double> out;
    for (std::size_t k = 0; k < valueCount(values); ++k) {
        out.push_back(pyFloat(values, k));
    }
    return out;
}

// numpy indexing: negative positions count from the end
std::size_t numpyIndex(std::int64_t index, std::size_t size) {
    const auto n = static_cast<std::int64_t>(size);
    if (index < -n || index >= n) {
        throw std::out_of_range("index " + std::to_string(index) +
                                " is out of bounds for axis 0 with size " + std::to_string(size));
    }
    return static_cast<std::size_t>(index < 0 ? index + n : index);
}

}  // namespace

LjTable LjTable::fromParm7(const std::filesystem::path& path) {
    static const char* const kFlags[] = {"POINTERS", "ATOM_TYPE_INDEX", "NONBONDED_PARM_INDEX",
                                         "LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF"};
    const Parm7Sections sections =
        readParm7Sections(path, std::set<std::string>(std::begin(kFlags), std::end(kFlags)));
    std::vector<std::string> missing;
    for (const char* flag : kFlags) {
        if (!findSection(sections, flag)) {
            missing.push_back(flag);
        }
    }
    if (!missing.empty()) {
        std::string list;
        for (const std::string& m : missing) {
            list += (list.empty() ? "'" : ", '") + m + "'";
        }
        throw ValueError(detail::pathText(path) + ": missing parm7 sections [" + list + "]");
    }
    LjTable table;
    const Parm7Values& pointers = *findSection(sections, "POINTERS");
    if (valueCount(pointers) < 2) {
        throw std::out_of_range("index 1 is out of bounds for axis 0 with size " +
                                std::to_string(valueCount(pointers)));
    }
    table.ntypes = std::stoll(pyIntText(pointers, 1));
    table.typeIndex = asIntegers(*findSection(sections, "ATOM_TYPE_INDEX"));
    table.nbIndex = asIntegers(*findSection(sections, "NONBONDED_PARM_INDEX"));
    table.acoef = asReals(*findSection(sections, "LENNARD_JONES_ACOEF"));
    table.bcoef = asReals(*findSection(sections, "LENNARD_JONES_BCOEF"));
    return table;
}

std::pair<std::vector<double>, std::vector<double>> LjTable::pairCoefficients(
    std::span<const std::int64_t> i, std::span<const std::int64_t> j) const {
    if (i.size() != j.size()) {
        throw std::invalid_argument("pairCoefficients: i and j differ in length");
    }
    std::vector<double> a(i.size(), 0.0);
    std::vector<double> b(i.size(), 0.0);
    for (std::size_t k = 0; k < i.size(); ++k) {
        const std::int64_t ti = typeIndex[numpyIndex(i[k], typeIndex.size())] - 1;
        const std::int64_t tj = typeIndex[numpyIndex(j[k], typeIndex.size())] - 1;
        const std::int64_t idx = nbIndex[numpyIndex(ti * ntypes + tj, nbIndex.size())];
        if (idx > 0) {
            a[k] = acoef[numpyIndex(idx - 1, acoef.size())];
            b[k] = bcoef[numpyIndex(idx - 1, bcoef.size())];
        }
    }
    return {std::move(a), std::move(b)};
}

std::vector<std::string> formatParm7Values(char kind, const Parm7Values& values) {
    const std::size_t n = valueCount(values);
    std::vector<std::string> lines;
    if (kind == 's') {
        for (std::size_t k = 0; k < n; ++k) {
            lines.push_back(leftField(pyStr(values, k), 80));
        }
        if (lines.empty()) {
            lines.emplace_back();
        }
        return lines;
    }
    std::size_t perLine = 0;
    switch (kind) {
    case 'a': perLine = 20; break;
    case 'i': perLine = 10; break;
    case 'e': perLine = 5; break;
    default: throw std::invalid_argument(std::string("parm7: unknown section kind '") + kind + "'");
    }
    std::vector<std::string> cells;
    for (std::size_t k = 0; k < n; ++k) {
        if (kind == 'a') {
            cells.push_back(leftField(pyStr(values, k), 4));
        } else if (kind == 'i') {
            cells.push_back(rightField(pyIntText(values, k), 8));
        } else {
            cells.push_back(formatE(pyFloat(values, k)));
        }
    }
    for (std::size_t k = 0; k < cells.size(); k += perLine) {
        std::string line;
        for (std::size_t c = k; c < std::min(cells.size(), k + perLine); ++c) {
            line += cells[c];
        }
        lines.push_back(std::move(line));
    }
    if (lines.empty()) {
        lines.emplace_back();
    }
    return lines;
}

void writeParm7(const std::filesystem::path& path, const std::vector<Parm7Section>& sections,
                std::string_view title) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        throw std::filesystem::filesystem_error("cannot write parm7 file", path,
                                                std::make_error_code(std::errc::io_error));
    }
    std::string out;
    // each piece is encoded as one write() of MDATK; what came before an error stays in the file
    const auto write = [&](const std::string& text) {
        checkAscii(text);
        out += text;
    };
    try {
        write("%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/26  00:00:00\n");
        write("%FLAG TITLE\n%FORMAT(20a4)\n");
        write(leftField(std::string(title), 80) + "\n");
        for (const Parm7Section& section : sections) {
            const char* fmt = section.kind == 'a'   ? "20a4"
                              : section.kind == 'i' ? "10I8"
                              : section.kind == 'e' ? "5E16.8"
                              : section.kind == 's' ? "1a80"
                                                    : nullptr;
            if (!fmt) {
                throw std::invalid_argument(std::string("parm7: unknown section kind '") +
                                            section.kind + "'");
            }
            write("%FLAG " + section.flag + "\n%FORMAT(" + fmt + ")\n");
            for (const std::string& line : formatParm7Values(section.kind, section.values)) {
                write(line + "\n");
            }
        }
    } catch (...) {
        file.write(out.data(), static_cast<std::streamsize>(out.size()));
        throw;
    }
    file.write(out.data(), static_cast<std::streamsize>(out.size()));
    if (!file.flush()) {
        throw std::filesystem::filesystem_error("cannot write parm7 file", path,
                                                std::make_error_code(std::errc::io_error));
    }
}

std::shared_ptr<core::Topology> readParm7Topology(const std::filesystem::path& path) {
    return Parm7TopologyReader(path).read();
}

bool isParm7(std::string_view path) {
    std::string ext = core::posixpath::splitext(core::pathBasename(path)).second;
    for (char& c : ext) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return ext == ".parm7" || ext == ".prmtop" || ext == ".top";
}

}  // namespace mdbridge::io
