#include "mdbridge/io/pdb.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <utility>

#include "mdbridge/core/detail/data_tables.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "python_text.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace detail {

namespace {
#include "pdb_tables.inc"
}  // namespace

U32 toU32(std::string_view utf8) {
    return core::decodeUtf8(utf8);
}

std::string toUtf8(U32View text) {
    return core::encodeUtf8(text);
}

U32 slice(U32View text, std::size_t first, std::size_t last) {
    last = std::min(last, text.size());
    if (first >= last) {
        return {};
    }
    return U32(text.substr(first, last - first));
}

bool isSpace(char32_t c) {
    return core::isPythonSpace(c);
}

U32 lstrip(U32View text) {
    std::size_t first = 0;
    while (first < text.size() && isSpace(text[first])) {
        ++first;
    }
    return U32(text.substr(first));
}

U32 rstrip(U32View text) {
    std::size_t last = text.size();
    while (last > 0 && isSpace(text[last - 1])) {
        --last;
    }
    return U32(text.substr(0, last));
}

U32 strip(U32View text) {
    return lstrip(rstrip(text));
}

std::vector<U32> split(U32View text) {
    std::vector<U32> words;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && isSpace(text[i])) {
            ++i;
        }
        const std::size_t first = i;
        while (i < text.size() && !isSpace(text[i])) {
            ++i;
        }
        if (i > first) {
            words.emplace_back(text.substr(first, i - first));
        }
    }
    return words;
}

U32 ljust(U32View text, std::size_t width) {
    U32 out(text);
    if (out.size() < width) {
        out.append(width - out.size(), U' ');
    }
    return out;
}

U32 rjust(U32View text, std::size_t width) {
    if (text.size() >= width) {
        return U32(text);
    }
    return U32(width - text.size(), U' ') + U32(text);
}

U32 upper(U32View text) {
    U32 out;
    out.reserve(text.size());
    for (const char32_t c : text) {
        if (c < 0x80) {
            out += c >= U'a' && c <= U'z' ? static_cast<char32_t>(c - 32) : c;
            continue;
        }
        const auto special = std::find_if(std::begin(kUpperSpecial), std::end(kUpperSpecial),
                                          [c](const UpperSpecial& s) { return s.code == c; });
        if (special != std::end(kUpperSpecial)) {
            out += toU32(special->upper);
            continue;
        }
        const auto run =
            std::upper_bound(std::begin(kUpperRuns), std::end(kUpperRuns), c,
                             [](char32_t code, const UpperRun& r) { return code < r.first; });
        if (run != std::begin(kUpperRuns)) {
            const UpperRun& r = *(run - 1);
            if (c <= r.last && (r.step == 0 ? c == r.first : (c - r.first) % r.step == 0)) {
                out += static_cast<char32_t>(static_cast<std::int64_t>(c) + r.delta);
                continue;
            }
        }
        out += c;
    }
    return out;
}

bool isDigit(char32_t c) {
    return core::isPythonDigit(c);
}

bool isAlnum(char32_t c) {
    return c != U'_' && core::isPythonWordCharacter(c);
}

bool isAlpha(char32_t c) {
    if (!isAlnum(c) || isDigit(c)) {
        return false;
    }
    return std::none_of(std::begin(kWordNotAlpha), std::end(kWordNotAlpha),
                        [c](const std::array<char32_t, 2>& r) { return c >= r[0] && c <= r[1]; });
}

bool isAllAlpha(U32View text) {
    return !text.empty() &&
           std::all_of(text.begin(), text.end(), [](char32_t c) { return isAlpha(c); });
}

std::optional<std::int64_t> toInt(U32View text) {
    try {
        return core::tryPythonInt(toUtf8(text));
    } catch (const std::overflow_error&) {
        return std::nullopt;
    }
}

std::optional<double> toFloat(U32View text) {
    return core::tryPythonFloat(toUtf8(text));
}

std::string fixed(double x, int width, int precision) {
    std::string body;
    if (std::isnan(x)) {
        body = "nan";
    } else if (std::isinf(x)) {
        body = x > 0 ? "inf" : "-inf";
    } else {
        const int n = std::snprintf(nullptr, 0, "%.*f", precision, x);
        body.resize(static_cast<std::size_t>(n) + 1);
        std::snprintf(body.data(), body.size(), "%.*f", precision, x);
        body.resize(static_cast<std::size_t>(n));
    }
    if (body.size() < static_cast<std::size_t>(width)) {
        body.insert(0, static_cast<std::size_t>(width) - body.size(), ' ');
    }
    return body;
}

std::string padInt(std::int64_t n, int width) {
    std::string body = std::to_string(n);
    if (body.size() < static_cast<std::size_t>(width)) {
        body.insert(0, static_cast<std::size_t>(width) - body.size(), ' ');
    }
    return body;
}

}  // namespace detail

namespace {

using namespace detail;

// Python's str.splitlines(keepends=True)
std::vector<U32> splitlinesKeepends(U32View text) {
    std::vector<U32> lines;
    std::size_t first = 0;
    std::size_t i = 0;
    while (i < text.size()) {
        const char32_t c = text[i];
        const bool lineBreak = c == U'\n' || c == U'\r' || c == 0x0B || c == 0x0C || c == 0x1C ||
                               c == 0x1D || c == 0x1E || c == 0x85 || c == 0x2028 || c == 0x2029;
        ++i;
        if (!lineBreak) {
            continue;
        }
        if (c == U'\r' && i < text.size() && text[i] == U'\n') {
            ++i;
        }
        lines.emplace_back(text.substr(first, i - first));
        first = i;
    }
    if (first < text.size()) {
        lines.emplace_back(text.substr(first));
    }
    return lines;
}

std::vector<U32> splitlines(U32View text) {
    std::vector<U32> lines = splitlinesKeepends(text);
    for (U32& line : lines) {
        std::size_t end = line.size();
        if (end > 0) {
            const char32_t c = line[end - 1];
            if (c == U'\n' && end > 1 && line[end - 2] == U'\r') {
                end -= 2;
            } else if (c == U'\n' || c == U'\r' || c == 0x0B || c == 0x0C || c == 0x1C ||
                       c == 0x1D || c == 0x1E || c == 0x85 || c == 0x2028 || c == 0x2029) {
                end -= 1;
            }
        }
        line.resize(end);
    }
    return lines;
}

bool startsWith(U32View text, std::u32string_view prefix) {
    return text.substr(0, prefix.size()) == prefix;
}

std::string repr(U32View text) {
    return core::pythonStringRepr(toUtf8(text));
}

U32 recordOf(U32View line) {
    return upper(strip(slice(line, 0, 6)));
}

bool coordsOk(U32View tail) {
    return toFloat(slice(tail, 3, 11)) && toFloat(slice(tail, 11, 19)) &&
           toFloat(slice(tail, 19, 27));
}

struct Fused {
    U32 alt;
    U32 resname;
    U32 chain;
    std::size_t colon = 0;
};

// _fused_field: ChimeraX 'resname+chain:' field
std::optional<Fused> fusedField(U32View line) {
    std::size_t c = U32::npos;
    for (std::size_t i = 18; i < 30 && i < line.size(); ++i) {
        if (line[i] == U':') {
            c = i;
            break;
        }
    }
    if (c == U32::npos || isSpace(line[c - 1])) {
        return std::nullopt;
    }
    const std::vector<U32> words = split(slice(line, 16, c - 1));
    if (words.size() == 1) {
        return Fused{U" ", words[0], U32(1, line[c - 1]), c};
    }
    if (words.size() == 2 && words[0] == U32(1, line[16])) {
        return Fused{words[0], words[1], U32(1, line[c - 1]), c};
    }
    return std::nullopt;
}

// regex helpers of _tokens_tail (str patterns: \s and \d are Unicode)
bool matchCoord(U32View text, std::size_t at, std::size_t& end, bool exactlyThree) {
    std::size_t i = at;
    if (i < text.size() && text[i] == U'-') {
        ++i;
    }
    const std::size_t digits = i;
    while (i < text.size() && isDigit(text[i])) {
        ++i;
    }
    if (i == digits || i >= text.size() || text[i] != U'.') {
        return false;
    }
    ++i;
    const std::size_t fraction = i;
    if (exactlyThree) {
        while (i < text.size() && i - fraction < 3 && isDigit(text[i])) {
            ++i;
        }
        if (i - fraction != 3) {
            return false;
        }
    } else {
        while (i < text.size() && isDigit(text[i])) {
            ++i;
        }
        if (i == fraction) {
            return false;
        }
    }
    end = i;
    return true;
}

// re.finditer / findall of a pattern (non-overlapping, left to right) from `from`
std::vector<std::pair<std::size_t, std::size_t>> findAll(U32View text, std::size_t from,
                                                         bool exactlyThree) {
    std::vector<std::pair<std::size_t, std::size_t>> out;
    std::size_t i = from;
    while (i <= text.size()) {
        std::size_t end = 0;
        if (i < text.size() && matchCoord(text, i, end, exactlyThree)) {
            out.emplace_back(i, end);
            i = end;
        } else {
            ++i;
        }
    }
    return out;
}

struct TokensTail {
    std::int64_t resseq = 0;
    U32 icode;
    U32 tail;
    U32 element;
};

TokensTail tokensTail(U32View rest, std::size_t lineno) {
    const auto damaged = [&] {
        return PdbError("line " + std::to_string(lineno) + ": cannot read the damaged ATOM record");
    };
    // r"\s*(-?\d+)([A-Za-z]?)\s"
    std::size_t i = 0;
    while (i < rest.size() && isSpace(rest[i])) {
        ++i;
    }
    // \s* is greedy but may backtrack; a digit or '-' cannot be whitespace, so no backtracking
    // helps
    const std::size_t numberStart = i;
    if (i < rest.size() && rest[i] == U'-') {
        ++i;
    }
    const std::size_t digitStart = i;
    while (i < rest.size() && isDigit(rest[i])) {
        ++i;
    }
    std::optional<std::size_t> matchEnd;
    std::size_t numberEnd = i;
    U32 letter;
    if (i > digitStart) {
        // (-?\d+) may backtrack over digits; ([A-Za-z]?)\s must follow
        for (std::size_t end = i; end > digitStart && !matchEnd; --end) {
            const char32_t next = end < rest.size() ? rest[end] : 0;
            const bool ascii = (next >= U'A' && next <= U'Z') || (next >= U'a' && next <= U'z');
            if (ascii && end + 1 < rest.size() && isSpace(rest[end + 1])) {
                numberEnd = end;
                letter = U32(1, next);
                matchEnd = end + 2;
            } else if (end < rest.size() && isSpace(next)) {
                numberEnd = end;
                matchEnd = end + 1;
            }
        }
    }
    if (!matchEnd) {
        throw damaged();
    }
    const auto coords = findAll(rest, *matchEnd, true);
    if (coords.size() < 3) {
        throw damaged();
    }
    std::array<double, 3> xyz{};
    for (int k = 0; k < 3; ++k) {
        xyz[k] = *toFloat(slice(rest, coords[k].first, coords[k].second));
    }
    const U32 after = slice(rest, coords[2].second);
    std::vector<double> nums;
    for (const auto& [first, last] : findAll(after, 0, false)) {
        if (nums.size() == 2) {
            break;
        }
        nums.push_back(*toFloat(slice(after, first, last)));
    }
    const double occupancy = nums.size() > 0 ? nums[0] : 1.0;
    const double bfactor = nums.size() > 1 ? nums[1] : 0.0;
    const std::vector<U32> words = split(after);
    TokensTail out;
    if (!words.empty() && isAllAlpha(words.back()) && words.back().size() <= 2) {
        out.element = words.back();
    }
    const std::string tail = "   " + fixed(xyz[0], 8, 3) + fixed(xyz[1], 8, 3) +
                             fixed(xyz[2], 8, 3) + fixed(occupancy, 6, 2) + fixed(bfactor, 6, 2);
    out.tail = ljust(toU32(tail), 49);
    const auto resseq = toInt(slice(rest, numberStart, numberEnd));
    if (!resseq) {
        throw core::ValueError("line " + std::to_string(lineno) +
                               ": residue number beyond the 64-bit range");
    }
    out.resseq = *resseq;
    out.icode = letter;
    return out;
}

struct ParsedAtom {
    PdbAtom atom;
    U32 resname;
    U32 chain;
    std::int64_t resseq = 0;
    U32 icode;
    bool fused = false;
};

ParsedAtom parseAtom(U32View raw, std::size_t lineno) {
    const U32 line = ljust(raw, 80);
    ParsedAtom out;
    U32 alt;
    U32 tail;
    U32 element;
    U32 charge;
    if (const auto fused = fusedField(line)) {
        alt = fused->alt;
        out.resname = fused->resname;
        out.chain = strip(fused->chain);
        out.fused = true;
        const U32 rest = ljust(slice(line, fused->colon + 1), 58);
        bool ok = false;
        if (const auto resseq = toInt(slice(rest, 0, 4))) {
            out.resseq = *resseq;
            out.icode = strip(slice(rest, 4, 5));
            tail = slice(rest, 5, 54);
            element = strip(slice(rest, 54, 56));
            charge = strip(slice(rest, 56, 58));
            ok = coordsOk(tail);
        }
        if (!ok) {
            TokensTail t = tokensTail(rest, lineno);
            out.resseq = t.resseq;
            out.icode = t.icode;
            tail = t.tail;
            element = t.element;
            charge.clear();
        }
    } else {
        alt = slice(line, 16, 17);
        out.resname = strip(slice(line, 17, 21));
        out.chain = strip(slice(line, 21, 22));
        const auto resseq = toInt(slice(line, 22, 26));
        if (!resseq) {
            throw PdbError("line " + std::to_string(lineno) + ": unreadable residue number " +
                           repr(slice(line, 22, 26)));
        }
        out.resseq = *resseq;
        out.icode = strip(slice(line, 26, 27));
        tail = slice(line, 27, 76);
        element = strip(slice(line, 76, 78));
        charge = strip(slice(line, 78, 80));
        if (!coordsOk(tail)) {
            // CH-01: a number below -999 as MDATK wrote it ("A-1234"), every later column shifted
            const U32 number = slice(line, 22, 27);
            const bool wide = number.size() == 5 && number[0] == U'-' &&
                              std::all_of(number.begin() + 1, number.end(),
                                          [](char32_t c) { return c >= U'0' && c <= U'9'; });
            if (wide && coordsOk(slice(line, 28, 77))) {
                throw PdbError("line " + std::to_string(lineno) + ": residue number " +
                               std::to_string(*toInt(number)) +
                               " does not fit the PDB columns (-999 to 9999); the later columns are shifted");
            }
            throw PdbError("line " + std::to_string(lineno) + ": unreadable coordinates");
        }
    }
    out.atom.serial = toInt(slice(line, 6, 11));
    out.atom.record = toUtf8(strip(slice(line, 0, 6)));
    out.atom.name = toUtf8(strip(slice(line, 12, 16)));
    out.atom.alt = toUtf8(strip(alt));
    out.atom.tail = toUtf8(tail);
    out.atom.element = toUtf8(element);
    out.atom.charge = toUtf8(charge);
    return out;
}

// --- writing

bool twoLetterElement(U32View bare) {
    static constexpr std::array<std::u32string_view, 23> kTwoLetter = {
        U"MG", U"NA", U"CL", U"ZN", U"MN", U"CA", U"LI", U"RB", U"CS", U"SR", U"BA", U"CD",
        U"FE", U"CO", U"NI", U"CU", U"BR", U"AG", U"HG", U"PB", U"TL", U"AU", U"PT"};
    return std::find(kTwoLetter.begin(), kTwoLetter.end(), bare) != kTwoLetter.end();
}

// re.sub(r"[\d+\-]", "", text)
U32 withoutDigitsAndSigns(U32View text) {
    U32 out;
    for (const char32_t c : text) {
        if (!isDigit(c) && c != U'+' && c != U'-') {
            out += c;
        }
    }
    return out;
}

U32 elementOf(const PdbAtom& atom, const std::string& resname) {
    if (!atom.element.empty()) {
        return toU32(atom.element);
    }
    const U32 name = strip(toU32(atom.name));
    const U32 bare = upper(withoutDigitsAndSigns(name));
    if (twoLetterElement(bare) && bare == upper(withoutDigitsAndSigns(toU32(resname)))) {
        return bare;
    }
    std::size_t first = 0;
    while (first < name.size() && isDigit(name[first])) {
        ++first;
    }
    return first < name.size() ? upper(name.substr(first, 1)) : U32();
}

U32 atomNameField(U32View name, U32View element) {
    if (name.size() >= 4 || element.size() == 2) {
        return ljust(slice(name, 0, 4), 4);
    }
    return U" " + ljust(name, 3);
}

U32 residueField(U32View resname) {
    return resname.size() <= 3 ? rjust(resname, 3) + U" " : slice(resname, 0, 4);
}

U32 firstOr(const std::string& text, char32_t fallback) {
    const U32 value = toU32(text);
    return value.empty() ? U32(1, fallback) : value.substr(0, 1);
}

U32 idColumns(const PdbResidue& res) {
    if (res.resseq < -999) {  // CH-01: "-1002" would shift every later column
        throw PdbError("residue number " + std::to_string(res.resseq) + " (" +
                       (res.chain.empty() ? std::string("_") : res.chain) +
                       ") does not fit the PDB columns (-999 to 9999)");
    }
    const std::int64_t resseq = res.resseq > 9999 ? res.resseq % 10000 : res.resseq;
    return residueField(toU32(res.name)) + firstOr(res.chain, U' ') + toU32(padInt(resseq, 4)) +
           firstOr(res.icode, U' ');
}

U32 atomLine(const PdbAtom& atom, const PdbResidue& res, std::int64_t serial, bool anisou) {
    const U32 element = elementOf(atom, res.name);
    const U32 body = anisou ? slice(ljust(toU32(atom.anisou), 53), 0, 49)
                            : slice(ljust(toU32(atom.tail), 49), 0, 49);
    const U32 record = anisou ? U32(U"ANISOU") : toU32(atom.record);
    return ljust(record, 6) + toU32(padInt(serial % 100000, 5)) + U" " +
           atomNameField(toU32(atom.name), element) + firstOr(atom.alt, U' ') + idColumns(res) +
           body + rjust(element, 2) + ljust(toU32(atom.charge), 2);
}

std::vector<U32> remapConect(const std::vector<std::string>& lines,
                             const std::map<std::int64_t, std::vector<std::int64_t>>& serialMap) {
    for (const auto& [serial, values] : serialMap) {
        if (values.size() > 1) {
            return {};
        }
    }
    std::vector<U32> out;
    for (const std::string& text : lines) {
        const U32 line = ljust(toU32(text), 31);
        std::vector<std::int64_t> ids;
        bool bad = false;
        for (std::size_t i = 6; i < 31; i += 5) {
            const U32 field = slice(line, i, i + 5);
            if (strip(field).empty()) {
                continue;
            }
            const auto value = toInt(field);
            if (!value) {
                bad = true;
                break;
            }
            ids.push_back(*value);
        }
        if (bad || ids.empty() || !serialMap.contains(ids[0])) {
            continue;
        }
        std::vector<std::int64_t> mapped = {serialMap.at(ids[0])[0]};
        for (std::size_t k = 1; k < ids.size(); ++k) {
            if (const auto it = serialMap.find(ids[k]); it != serialMap.end()) {
                mapped.push_back(it->second[0]);
            }
        }
        if (mapped.size() > 1) {
            U32 record = U"CONECT";
            for (const std::int64_t serial : mapped) {
                record += toU32(padInt(serial % 100000, 5));
            }
            out.push_back(std::move(record));
        }
    }
    return out;
}

std::size_t seqresInsertIndex(const std::vector<std::string>& header) {
    static constexpr std::array<std::u32string_view, 15> kOrder = {
        U"MODRES", U"HET",    U"HETNAM", U"HETSYN", U"FORMUL", U"HELIX", U"SHEET", U"SSBOND",
        U"LINK",   U"CISPEP", U"SITE",   U"CRYST1", U"ORIGX",  U"SCALE", U"MTRIX"};
    for (std::size_t k = 0; k < header.size(); ++k) {
        const U32 record = recordOf(toU32(header[k]));
        for (const auto prefix : kOrder) {
            if (startsWith(record, prefix)) {
                return k;
            }
        }
    }
    return header.size();
}

// numpy's round-half-even to an integer
std::int64_t roundHalfEven(double v) {
    return static_cast<std::int64_t>(std::nearbyint(v));
}

U32 transformAnisou(const U32& line, const std::array<std::array<double, 3>, 3>& r) {
    std::array<std::int64_t, 6> u{};
    for (std::size_t k = 0; k < 6; ++k) {
        const auto value = toInt(slice(line, 28 + 7 * k, 35 + 7 * k));
        if (!value) {
            return line;
        }
        u[k] = *value;
    }
    const double m[3][3] = {{double(u[0]), double(u[3]), double(u[4])},
                            {double(u[3]), double(u[1]), double(u[5])},
                            {double(u[4]), double(u[5]), double(u[2])}};
    double rm[3][3];
    double out[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rm[i][j] = r[i][0] * m[0][j] + r[i][1] * m[1][j] + r[i][2] * m[2][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out[i][j] = rm[i][0] * r[j][0] + rm[i][1] * r[j][1] + rm[i][2] * r[j][2];
        }
    }
    const double values[6] = {out[0][0], out[1][1], out[2][2], out[0][1], out[0][2], out[1][2]};
    U32 result = slice(line, 0, 28);
    for (const double v : values) {
        result += toU32(padInt(roundHalfEven(v), 7));
    }
    return result + slice(line, 70);
}

U32 rstripLineEnd(U32View line) {
    std::size_t end = line.size();
    while (end > 0 && (line[end - 1] == U'\r' || line[end - 1] == U'\n')) {
        --end;
    }
    return U32(line.substr(0, end));
}

}  // namespace

// --- model

std::array<double, 3> PdbAtom::xyz() const {
    const U32 t = toU32(tail);
    std::array<double, 3> out{};
    for (std::size_t k = 0; k < 3; ++k) {
        out[k] = core::pythonFloat(toUtf8(slice(t, 3 + 8 * k, 11 + 8 * k)));
    }
    return out;
}

std::string PdbResidue::label() const {
    return (chain.empty() ? std::string("_") : chain) + ":" + std::to_string(resseq) + icode + " " +
           name;
}

std::vector<std::pair<std::string, std::array<double, 3>>> PdbResidue::coords() const {
    std::vector<std::pair<std::string, std::array<double, 3>>> out;
    for (const PdbAtom& atom : atoms) {
        std::string key = atom.name;
        std::replace(key.begin(), key.end(), '*', '\'');
        const bool seen = std::any_of(out.begin(), out.end(),
                                      [&](const auto& entry) { return entry.first == key; });
        if (!seen) {
            out.emplace_back(std::move(key), atom.xyz());
        }
    }
    return out;
}

std::vector<const PdbResidue*> PdbStructure::residues() const {
    std::vector<const PdbResidue*> out;
    for (const auto& segment : segments) {
        for (const PdbResidue& residue : segment) {
            out.push_back(&residue);
        }
    }
    return out;
}

std::vector<std::string> PdbStructure::chainIds() const {
    std::vector<std::string> out;
    for (const PdbResidue* residue : residues()) {
        if (std::find(out.begin(), out.end(), residue->chain) == out.end()) {
            out.push_back(residue->chain);
        }
    }
    return out;
}

// --- reading

PdbStructure readPdbText(std::string_view text) {
    PdbStructure s;
    std::vector<PdbResidue> seg;
    PdbResidue* res = nullptr;
    bool started = false;
    bool ended = false;
    bool skipping = false;
    std::int64_t models = 0;
    std::int64_t repaired = 0;
    std::int64_t dropped = 0;
    std::int64_t afterEnd = 0;

    const auto closeSegment = [&] {
        if (!seg.empty()) {
            s.segments.push_back(std::move(seg));
        }
        seg.clear();
        res = nullptr;
    };

    const std::vector<U32> lines = splitlines(toU32(text));
    for (std::size_t k = 0; k < lines.size(); ++k) {
        const std::size_t lineno = k + 1;
        const U32 line = rstrip(lines[k]);
        const U32 rec = recordOf(line);
        if (line.empty()) {
            continue;
        }
        const bool atomRecord = rec == U"ATOM" || rec == U"HETATM";
        if (ended) {
            afterEnd += atomRecord ? 1 : 0;
            continue;
        }
        if (rec == U"MODEL") {
            ++models;
            if (models == 1) {
                s.model = toUtf8(line);
            } else {
                skipping = true;
            }
            continue;
        }
        if (rec == U"ENDMDL") {
            closeSegment();
            skipping = true;
            continue;
        }
        const bool sig = rec == U"SIGATM" || rec == U"SIGUIJ" || rec == U"ANISOU";
        if ((atomRecord || rec == U"TER" || sig) && skipping) {
            continue;
        }
        if (atomRecord) {
            started = true;
            ParsedAtom parsed = parseAtom(line, lineno);
            repaired += parsed.fused ? 1 : 0;
            const std::string resname = toUtf8(parsed.resname);
            const std::string chain = toUtf8(parsed.chain);
            const std::string icode = toUtf8(parsed.icode);
            if (!seg.empty() && seg.back().chain != chain) {
                closeSegment();
            }
            if (res == nullptr || res->name != resname || res->chain != chain ||
                res->resseq != parsed.resseq || res->icode != icode) {
                PdbResidue residue;
                residue.name = resname;
                residue.chain = chain;
                residue.resseq = parsed.resseq;
                residue.icode = icode;
                seg.push_back(std::move(residue));
                res = &seg.back();
            }
            res->atoms.push_back(std::move(parsed.atom));
        } else if (rec == U"TER") {
            closeSegment();
        } else if (rec == U"ANISOU" && res != nullptr) {
            res->atoms.back().anisou = toUtf8(slice(ljust(line, 80), 27, 80));
        } else if (sig) {
            ++dropped;
        } else if (rec == U"CONECT") {
            s.conect.push_back(toUtf8(line));
        } else if (rec == U"MASTER") {
            s.notes.emplace_back("Dropped the MASTER record (its counts no longer apply)");
        } else if (rec == U"END") {
            closeSegment();
            ended = true;
        } else if (started) {
            s.trailer.push_back(toUtf8(line));
        } else {
            s.header.push_back(toUtf8(line));
        }
    }
    closeSegment();
    if (s.segments.empty()) {
        throw PdbError("no ATOM/HETATM records found (not a PDB file, or the path does not exist)");
    }
    if (repaired) {
        s.notes.push_back(
            "Repaired " + std::to_string(repaired) +
            " ChimeraX-damaged ATOM/HETATM lines (residue name fused with the chain id and "
            "a colon, e.g. 'GC:' -> residue G, chain C)");
    }
    if (models > 1) {
        s.notes.push_back("Kept model 1 of " + std::to_string(models));
    }
    if (afterEnd) {
        s.notes.push_back("Ignored " + std::to_string(afterEnd) +
                          " ATOM/HETATM lines after END (only the first model is kept)");
    }
    if (dropped) {
        s.notes.push_back("Dropped " + std::to_string(dropped) +
                          " SIGATM/SIGUIJ/orphan ANISOU records");
    }
    return s;
}

namespace detail {

std::optional<std::string> utf8Error(std::string_view bytes) {
    const auto byte = [&](std::size_t i) { return static_cast<unsigned char>(bytes[i]); };
    std::size_t i = 0;
    while (i < bytes.size()) {
        const unsigned char b = byte(i);
        if (b < 0x80) {
            ++i;
            continue;
        }
        std::size_t length = 0;
        unsigned char low = 0x80;
        unsigned char high = 0xBF;
        if (b >= 0xC2 && b <= 0xDF) {
            length = 2;
        } else if (b >= 0xE0 && b <= 0xEF) {
            length = 3;
            low = b == 0xE0 ? 0xA0 : 0x80;
            high = b == 0xED ? 0x9F : 0xBF;
        } else if (b >= 0xF0 && b <= 0xF4) {
            length = 4;
            low = b == 0xF0 ? 0x90 : 0x80;
            high = b == 0xF4 ? 0x8F : 0xBF;
        }
        std::size_t k = 1;
        std::string reason = "invalid start byte";
        if (length > 0) {
            for (; k < length && i + k < bytes.size(); ++k) {
                const unsigned char c = byte(i + k);
                if (c < (k == 1 ? low : 0x80) || c > (k == 1 ? high : 0xBF)) {
                    break;
                }
            }
            if (k == length) {
                i += length;
                continue;
            }
            reason = i + k >= bytes.size() ? "unexpected end of data" : "invalid continuation byte";
        }
        char hex[8];
        std::snprintf(hex, sizeof hex, "0x%02x", b);
        if (k == 1) {
            return "'utf-8' codec can't decode byte " + std::string(hex) + " in position " +
                   std::to_string(i) + ": " + reason;
        }
        return "'utf-8' codec can't decode bytes in position " + std::to_string(i) + "-" +
               std::to_string(i + k - 1) + ": " + reason;
    }
    return std::nullopt;
}

U32 decodeStrict(std::string_view bytes) {
    if (const auto error = utf8Error(bytes)) {
        throw UnicodeDecodeError(*error);
    }
    return toU32(bytes);
}

U32 universalNewlines(U32View text) {
    U32 out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == U'\r') {
            out += U'\n';
            if (i + 1 < text.size() && text[i + 1] == U'\n') {
                ++i;
            }
        } else {
            out += text[i];
        }
    }
    return out;
}

std::vector<U32> fileLines(U32View text) {
    std::vector<U32> lines;
    std::size_t first = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == U'\n') {
            lines.emplace_back(text.substr(first, i + 1 - first));
            first = i + 1;
        }
    }
    if (first < text.size()) {
        lines.emplace_back(text.substr(first));
    }
    return lines;
}

}  // namespace detail

PdbStructure readPdbFile(const std::filesystem::path& path) {
    return readPdbText(readPdbFileText(path));
}

std::string readPdbFileText(const std::filesystem::path& path) {
    // universal newlines turn \r\n and \r into \n, which splitlines treats alike
    return toUtf8(decodeStrict(readBytes(path)));
}

// --- writing

std::vector<std::string> seqresLines(const std::vector<SeqresChain>& chains) {
    std::vector<std::string> out;
    for (const SeqresChain& chain : chains) {
        const std::size_t n = chain.names.size();
        for (std::size_t k = 0; k < n; k += 13) {
            U32 body;
            for (std::size_t i = k; i < std::min(n, k + 13); ++i) {
                if (core::characterCount(chain.names[i]) > 3) {  // CH-07
                    throw PdbError("SEQRES residue names have at most three characters: " +
                                   core::pythonStringRepr(chain.names[i]));
                }
                body += (i > k ? U" " : U"") + rjust(toU32(chain.names[i]), 3);
            }
            const U32 line = U"SEQRES " + toU32(padInt(static_cast<std::int64_t>(k / 13 + 1), 3)) +
                             U" " + firstOr(chain.chain, U' ') + U" " +
                             toU32(padInt(static_cast<std::int64_t>(n), 4)) + U"  " + body;
            out.push_back(toUtf8(line));
        }
    }
    return out;
}

std::string writePdb(const PdbStructure& structure, bool stripConect,
                     const std::optional<std::vector<SeqresChain>>& seqres) {
    std::vector<std::string> header = structure.header;
    if (seqres) {
        std::erase_if(header, [](const std::string& line) { return line.starts_with("SEQRES"); });
        const std::size_t k = seqresInsertIndex(header);
        const std::vector<std::string> records = seqresLines(*seqres);
        header.insert(header.begin() + static_cast<std::ptrdiff_t>(k), records.begin(),
                      records.end());
    }
    std::vector<U32> out;
    for (const std::string& line : header) {
        out.push_back(toU32(line));
    }
    if (!structure.model.empty()) {
        out.push_back(toU32(structure.model));
    }
    std::int64_t serial = 0;
    std::map<std::int64_t, std::vector<std::int64_t>> serialMap;
    for (const auto& seg : structure.segments) {
        if (seg.empty()) {
            continue;
        }
        for (const PdbResidue& res : seg) {
            for (const PdbAtom& atom : res.atoms) {
                ++serial;
                if (atom.serial) {
                    serialMap[*atom.serial].push_back(serial);
                }
                out.push_back(atomLine(atom, res, serial, false));
                if (!atom.anisou.empty()) {
                    out.push_back(atomLine(atom, res, serial, true));
                }
            }
        }
        ++serial;
        out.push_back(U"TER   " + toU32(padInt(serial % 100000, 5)) + U"      " +
                      idColumns(seg.back()));
    }
    if (!structure.model.empty()) {
        out.push_back(U"ENDMDL");
    }
    for (const std::string& line : structure.trailer) {
        out.push_back(toU32(line));
    }
    if (!stripConect) {
        for (U32& line : remapConect(structure.conect, serialMap)) {
            out.push_back(std::move(line));
        }
    }
    out.push_back(U"END");
    std::string text;
    for (const U32& line : out) {
        text += toUtf8(rstrip(line));
        text += '\n';
    }
    return text;
}

// --- residue names

namespace {

// L-DNA (MDATK 0.1.x _LDNA_RE): the CCD's L-deoxyribonucleotides and LD[ACGT][53N]?
std::optional<ResidueName> parseLDnaName(const U32& n) {
    const std::string name = toUtf8(n);
    for (const auto& component : core::detail::tables::kIndLDeoxyribonucleotides) {
        if (name == component.residue) {
            return ResidueName{component.base, "dna", 'L', ""};
        }
    }
    if (n.size() >= 3 && n.size() <= 4 && n[0] == U'L' && n[1] == U'D' &&
        std::u32string_view(U"ACGT").find(n[2]) != std::u32string_view::npos &&
        (n.size() == 3 || std::u32string_view(U"53N").find(n[3]) != std::u32string_view::npos)) {
        return ResidueName{static_cast<char>(n[2]), "dna", 'L', n.size() == 4 ? std::string(1, static_cast<char>(n[3])) : ""};
    }
    return std::nullopt;
}

}  // namespace

std::optional<ResidueName> parseResidueName(std::string_view name) {
    const U32 n = upper(strip(toU32(name)));
    static const std::array<std::pair<std::u32string_view, std::pair<char, const char*>>, 6>
        kAliases = {{
            {U"ADE", {'A', "rna"}},
            {U"CYT", {'C', "rna"}},
            {U"GUA", {'G', "rna"}},
            {U"URA", {'U', "rna"}},
            {U"URI", {'U', "rna"}},
            {U"THY", {'T', "dna"}},
        }};
    for (const auto& [alias, value] : kAliases) {
        if (n == alias) {
            return ResidueName{value.first, value.second, 'D', ""};
        }
    }
    // CH-22: the modified ribonucleotides MDBridge writes (the CCD's CH: C+), read as their parent
    const std::string utf8 = toUtf8(n);
    for (const auto& component : core::detail::tables::kIndModifiedRibonucleotides) {
        if (utf8 == component.residue) {
            return ResidueName{component.base, "rna", 'D', "", std::string(component.residue)};
        }
    }
    // ^([RL0D]?)([ACGTU])([53N]?)$
    std::size_t i = 0;
    char pre = 0;
    const auto in = [](char32_t c, std::u32string_view set) {
        return set.find(c) != std::u32string_view::npos;
    };
    if (i < n.size() && in(n[i], U"RL0D") && i + 1 < n.size() && in(n[i + 1], U"ACGTU")) {
        pre = static_cast<char>(n[i]);
        ++i;
    }
    if (i >= n.size() || !in(n[i], U"ACGTU")) {
        return parseLDnaName(n);
    }
    const char base = static_cast<char>(n[i]);
    ++i;
    std::string term;
    if (i < n.size() && in(n[i], U"53N")) {
        term = std::string(1, static_cast<char>(n[i]));
        ++i;
    }
    if (i != n.size()) {
        return std::nullopt;
    }
    if (pre == 'D' || (base == 'T' && pre == 0)) {
        return ResidueName{base, "dna", 'D', term};
    }
    if (base == 'T') {
        return std::nullopt;
    }
    return ResidueName{base, "rna", (pre == 'L' || pre == '0') ? 'L' : 'D', term};
}

std::string pythonUpper(std::string_view text) {
    return toUtf8(upper(toU32(text)));
}

bool isBackboneAtomName(std::string_view atomName) {
    return atomName.find('\'') != std::string_view::npos || atomName.find('*') != std::string_view::npos ||
           core::pythonStrip(atomName) == "P";
}

std::string standardAtomName(std::string_view resname, std::string_view name) {
    const std::string residue = toUtf8(upper(strip(toU32(resname))));
    for (const auto& atom : core::detail::tables::kIndComponentAtomNames) {
        if (atom.residue == residue && atom.name == name) {
            return std::string(atom.standard);
        }
    }
    return std::string(name);
}

// --- annotation and builder text utilities

std::string insertHeaderRecords(std::string_view text, const std::vector<std::string>& records) {
    std::vector<U32> lines;
    for (U32& line : splitlines(toU32(text))) {
        if (!startsWith(line, U"REMARK 950") && !startsWith(line, U"MASTER")) {
            lines.push_back(std::move(line));
        }
    }
    std::size_t at = lines.size();
    for (std::size_t k = 0; k < lines.size(); ++k) {
        const U32 rec = recordOf(lines[k]);
        if (rec == U"MODEL" || rec == U"ATOM" || rec == U"HETATM") {
            at = k;
            break;
        }
    }
    std::vector<U32> all(lines.begin(), lines.begin() + static_cast<std::ptrdiff_t>(at));
    for (const std::string& record : records) {
        all.push_back(toU32(record));
    }
    all.insert(all.end(), lines.begin() + static_cast<std::ptrdiff_t>(at), lines.end());
    U32 out;
    for (std::size_t k = 0; k < all.size(); ++k) {
        out += (k ? U"\n" : U"") + all[k];
    }
    return toUtf8(out + U"\n");
}

std::optional<std::string> transformPdbText(std::string_view text,
                                            const std::array<std::array<double, 3>, 3>& rotation,
                                            const std::array<double, 3>& translation) {
    U32 out;
    for (U32 line : splitlinesKeepends(toU32(text))) {
        const U32 rec = slice(line, 0, 6);
        if (rec == U"ATOM  " || rec == U"HETATM") {
            const U32 body = rstripLineEnd(line);
            std::array<double, 3> xyz{};
            for (std::size_t k = 0; k < 3; ++k) {
                const auto value = toFloat(slice(body, 30 + 8 * k, 38 + 8 * k));
                if (!value) {
                    return std::nullopt;
                }
                xyz[k] = *value;
            }
            std::array<double, 3> moved{};
            for (std::size_t i = 0; i < 3; ++i) {
                moved[i] = rotation[i][0] * xyz[0] + rotation[i][1] * xyz[1] +
                           rotation[i][2] * xyz[2] + translation[i];
            }
            // numpy's min / max propagate NaN, which then fails both comparisons
            const bool nan = std::isnan(moved[0]) || std::isnan(moved[1]) || std::isnan(moved[2]);
            const double low = std::min({moved[0], moved[1], moved[2]});
            const double high = std::max({moved[0], moved[1], moved[2]});
            if (!nan && (low <= -999.9995 || high >= 9999.9995)) {
                throw BuildError("the fitted coordinates do not fit the PDB columns");
            }
            const U32 end = line.substr(body.size());
            line = slice(body, 0, 30) +
                   toU32(fixed(moved[0], 8, 3) + fixed(moved[1], 8, 3) + fixed(moved[2], 8, 3)) +
                   slice(body, 54) + end;
        } else if (rec == U"ANISOU") {
            const U32 body = rstripLineEnd(line);
            const U32 end = line.substr(body.size());
            line = rstrip(transformAnisou(ljust(body, 70), rotation)) + end;
        }
        out += line;
    }
    return toUtf8(out);
}

std::string newAtomTail(const std::array<double, 3>& xyz, double occupancy, double bfactor) {
    // Python's max() / min() keep the first of incomparable (NaN) values
    double high = std::abs(xyz[0]);
    double low = xyz[0];
    for (std::size_t k = 1; k < 3; ++k) {
        if (std::abs(xyz[k]) > high) {
            high = std::abs(xyz[k]);
        }
        if (xyz[k] < low) {
            low = xyz[k];
        }
    }
    if (high >= 9999.9995 || low <= -999.9995) {
        throw BuildError("coordinates do not fit the PDB columns");
    }
    const std::string tail = "   " + fixed(xyz[0], 8, 3) + fixed(xyz[1], 8, 3) +
                             fixed(xyz[2], 8, 3) + fixed(occupancy, 6, 2) + fixed(bfactor, 6, 2);
    return toUtf8(ljust(toU32(tail), 49));
}

}  // namespace mdbridge::io
