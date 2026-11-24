#include "mdbridge/io/parm7.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <system_error>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
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

// Topology reader: CLEAN ROOM: reimplement from scratchpad/cleanroom/parm7-topology.md

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
    (void)path;
    throw std::logic_error("clean room: parm7 topology not implemented");
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
