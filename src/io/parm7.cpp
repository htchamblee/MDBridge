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

// ---------------------------------------------------------------- MDAnalysis' TOPParser

constexpr std::array<const char*, 119> kElements = {
    "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
    "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

enum class Attr {
    Names,
    Charges,
    Elements,
    Masses,
    TypeIndices,
    Types,
    Resnames,
    ResPointers,
    BondH,
    BondA,
    AngH,
    AngA,
    DihH,
    DihA,
    ChainIds
};

struct KnownSection {
    std::u32string_view flag;
    int perRecord;
    int perLine;
    Attr attr;
    std::size_t pointer;  // index into POINTERS
};

constexpr KnownSection kKnown[] = {
    {U"ATOM_NAME", 1, 20, Attr::Names, 0},
    {U"CHARGE", 1, 5, Attr::Charges, 0},
    {U"ATOMIC_NUMBER", 1, 10, Attr::Elements, 0},
    {U"MASS", 1, 5, Attr::Masses, 0},
    {U"ATOM_TYPE_INDEX", 1, 10, Attr::TypeIndices, 0},
    {U"AMBER_ATOM_TYPE", 1, 20, Attr::Types, 0},
    {U"RESIDUE_LABEL", 1, 20, Attr::Resnames, 11},
    {U"RESIDUE_POINTER", 1, 10, Attr::ResPointers, 11},
    {U"BONDS_INC_HYDROGEN", 3, 10, Attr::BondH, 2},
    {U"BONDS_WITHOUT_HYDROGEN", 3, 10, Attr::BondA, 3},
    {U"ANGLES_INC_HYDROGEN", 4, 10, Attr::AngH, 4},
    {U"ANGLES_WITHOUT_HYDROGEN", 4, 10, Attr::AngA, 5},
    {U"DIHEDRALS_INC_HYDROGEN", 5, 10, Attr::DihH, 6},
    {U"DIHEDRALS_WITHOUT_HYDROGEN", 5, 10, Attr::DihA, 7},
    {U"RESIDUE_CHAINID", 1, 20, Attr::ChainIds, 11},
};

// The names MDAnalysis reports in its messages (KeyError of a missing section, length checks).
const char* attrName(Attr attr) {
    switch (attr) {
    case Attr::Names: return "name";
    case Attr::Charges: return "charge";
    case Attr::Elements: return "elements";
    case Attr::Masses: return "mass";
    case Attr::TypeIndices: return "type_indices";
    case Attr::Types: return "types";
    case Attr::Resnames: return "resname";
    case Attr::ResPointers: return "respoint";
    case Attr::BondH: return "bondh";
    case Attr::BondA: return "bonda";
    case Attr::AngH: return "angh";
    case Attr::AngA: return "anga";
    case Attr::DihH: return "dihh";
    case Attr::DihA: return "diha";
    case Attr::ChainIds: return "segids";
    }
    return "";
}

// The attribute names of the Universe's length check.
const char* topologyAttrName(Attr attr) {
    switch (attr) {
    case Attr::Names: return "names";
    case Attr::Charges: return "charges";
    case Attr::Elements: return "elements";
    case Attr::Masses: return "masses";
    case Attr::TypeIndices: return "type_indices";
    case Attr::Types: return "types";
    default: return "";
    }
}

std::int64_t floorDiv(std::int64_t a, std::int64_t b) {
    const std::int64_t q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

struct Field {
    std::size_t start;
    std::size_t stop;
};

// FORTRANReader(fmt).entries: comma-separated [count]<letter><width>[.d] items and [count]X skips.
std::vector<Field> fortranFields(std::u32string_view fmt) {
    std::vector<Field> fields;
    std::size_t position = 0;
    std::size_t begin = 0;
    const auto bad = [&] {
        return ValueError("unrecognized FORTRAN format " + repr(utf8(fmt)));
    };
    while (begin <= fmt.size()) {
        std::size_t end = fmt.find(U',', begin);
        if (end == std::u32string_view::npos) {
            end = fmt.size();
        }
        std::u32string_view item = fmt.substr(begin, end - begin);
        while (!item.empty() && core::isPythonSpace(item.back())) {
            item.remove_suffix(1);
        }
        std::size_t i = 0;
        std::size_t count = 0;
        bool hasCount = false;
        while (i < item.size() && item[i] >= U'0' && item[i] <= U'9') {
            count = count * 10 + static_cast<std::size_t>(item[i] - U'0');
            hasCount = true;
            ++i;
        }
        if (!hasCount) {
            count = 1;
        }
        if (i >= item.size()) {
            throw bad();
        }
        const char32_t letter = item[i] | 0x20;
        ++i;
        if (letter == U'x' && i == item.size()) {
            position += count;
        } else if (letter == U'a' || letter == U'i' || letter == U'e' || letter == U'f' ||
                   letter == U'd') {
            std::size_t width = 0;
            const std::size_t widthStart = i;
            while (i < item.size() && item[i] >= U'0' && item[i] <= U'9') {
                width = width * 10 + static_cast<std::size_t>(item[i] - U'0');
                ++i;
            }
            if (i == widthStart) {
                throw bad();
            }
            if (i < item.size() && item[i] == U'.') {
                ++i;
                while (i < item.size() && item[i] >= U'0' && item[i] <= U'9') {
                    ++i;
                }
            }
            if (i != item.size()) {
                throw bad();
            }
            for (std::size_t k = 0; k < count; ++k) {
                fields.push_back({position, position + width});
                position += width;
            }
        } else {
            throw bad();
        }
        begin = end + 1;
    }
    return fields;
}

// The attributes TOPParser collects, by section.
struct TopAttrs {
    std::vector<Attr> order;  // first insertion of each attribute
    std::map<Attr, std::vector<std::string>> strings;
    std::map<Attr, std::vector<double>> reals;
    std::map<Attr, std::vector<std::int64_t>> ints;

    void touch(Attr attr) {
        if (std::find(order.begin(), order.end(), attr) == order.end()) {
            order.push_back(attr);
        }
    }
    bool has(Attr attr) const { return std::find(order.begin(), order.end(), attr) != order.end(); }
    std::size_t size(Attr attr) const {
        if (auto it = strings.find(attr); it != strings.end()) {
            return it->second.size();
        }
        if (auto it = reals.find(attr); it != reals.end()) {
            return it->second.size();
        }
        return ints.at(attr).size();
    }
};

class TopReader {
public:
    explicit TopReader(const std::filesystem::path& path)
        : m_path(detail::pathText(path)), m_lines(detail::readTextLines(path)) {}

    std::shared_ptr<core::Topology> parse();

private:
    std::string m_path;
    Lines m_lines;
    std::size_t m_next = 0;

    bool atEnd() const { return m_next >= m_lines.size(); }
    const std::u32string& next() {
        if (atEnd()) {
            throw ValueError(m_path + ": unexpected end of file");
        }
        return m_lines[m_next++];
    }
    void readSection(const KnownSection& known, std::int64_t numlines, TopAttrs& attrs);
};

void TopReader::readSection(const KnownSection& known, std::int64_t numlines, TopAttrs& attrs) {
    std::u32string fmt;
    for (;;) {
        const std::u32string& line = next();
        if (startsWith(line, U"%FORMAT")) {
            const std::u32string_view rest =
                std::u32string_view(line).substr(std::min<std::size_t>(8, line.size()));
            fmt = rest.substr(0, rest.find(U')'));
            break;
        }
        if (!startsWith(line, U"%COMMENT")) {
            throw ValueError(
                "Invalid header line. Does not begin with either %FLAG, %FORMAT nor %COMMENT:\n" +
                utf8(line));
        }
    }
    const std::vector<Field> fields = fortranFields(fmt);
    std::vector<std::string> values;
    for (std::int64_t n = 0; n < numlines; ++n) {
        const std::u32string& line = next();
        for (const Field& field : fields) {
            if (field.start >= line.size()) {
                continue;
            }
            const std::u32string value =
                strip(std::u32string_view(line).substr(field.start, field.stop - field.start));
            if (!value.empty()) {
                values.push_back(utf8(value));
            }
        }
    }
    const Attr attr = known.attr;
    attrs.touch(attr);
    attrs.strings.erase(attr);
    attrs.reals.erase(attr);
    attrs.ints.erase(attr);
    switch (attr) {
    case Attr::Names:
    case Attr::Types:
    case Attr::Resnames:
    case Attr::ChainIds: attrs.strings[attr] = std::move(values); break;
    case Attr::Charges: {
        std::vector<double> out;
        for (const std::string& v : values) {
            const float charge = static_cast<float>(core::pythonFloat(v)) / 18.2223f;
            out.push_back(static_cast<double>(charge));
        }
        attrs.reals[attr] = std::move(out);
        break;
    }
    case Attr::Masses: {
        std::vector<double> out;
        for (const std::string& v : values) {
            out.push_back(core::pythonFloat(v));
        }
        attrs.reals[attr] = std::move(out);
        break;
    }
    case Attr::Elements: {
        std::vector<std::string> out;
        for (const std::string& v : values) {
            const std::int64_t z = core::pythonInt(v);
            if (z > static_cast<std::int64_t>(kElements.size()) - 1) {
                throw ValueError(std::to_string(z));  // KeyError
            }
            out.push_back(z > 0 ? kElements[static_cast<std::size_t>(z)] : "");
        }
        attrs.strings[attr] = std::move(out);
        break;
    }
    default: {
        std::vector<std::int64_t> out;
        for (const std::string& v : values) {
            const std::int64_t x = core::pythonInt(v);
            out.push_back(attr == Attr::ResPointers   ? x - 1
                          : attr == Attr::TypeIndices ? x
                                                      : floorDiv(x, 3));
        }
        if (known.perRecord >
            1) {  // parse_chunks: the atoms of each record, without the parameter index
            std::vector<std::int64_t> atoms;
            const auto per = static_cast<std::size_t>(known.perRecord);
            for (std::size_t k = 0; k + per <= out.size(); k += per) {
                atoms.insert(atoms.end(), out.begin() + static_cast<std::ptrdiff_t>(k),
                             out.begin() + static_cast<std::ptrdiff_t>(k + per - 1));
            }
            out = std::move(atoms);
        }
        attrs.ints[attr] = std::move(out);
        break;
    }
    }
}

// header.split("%FLAG")[1].strip()
std::optional<std::u32string> flagName(std::u32string_view line) {
    const std::size_t at = line.find(U"%FLAG");
    if (at == std::u32string_view::npos) {
        return std::nullopt;
    }
    std::u32string_view rest = line.substr(at + 5);
    return strip(rest.substr(0, rest.find(U"%FLAG")));
}

std::shared_ptr<core::Topology> TopReader::parse() {
    if (!startsWith(next(), U"%VE")) {
        throw ValueError(m_path + " is not a valid TOP file. %VE Missing in header");
    }
    const std::vector<std::u32string> title = splitWhitespace(next());
    if (title.size() < 2) {
        throw ValueError("list index out of range");
    }
    if (title[1] != U"TITLE") {
        if (title[1] == U"CTITLE") {
            throw ValueError(
                m_path +
                " is detected as a Chamber-style TOP file. At this time MDAnalysis does "
                "not support such topologies");
        }
        throw ValueError(m_path + " is not a valid TOP file. 'TITLE' missing in header");
    }
    while (!startsWith(next(), U"%FLAG POINTERS")) {
    }
    next();
    std::vector<std::int64_t> pointers;
    for (int k = 0; k < 4; ++k) {
        for (const std::u32string& token : splitWhitespace(next())) {
            pointers.push_back(core::pythonInt(utf8(token)));
        }
    }

    TopAttrs attrs;
    std::optional<std::int64_t> numlines;
    std::optional<std::u32string> section = flagName(next());
    if (!section) {
        throw ValueError("list index out of range");
    }
    while (section) {
        const KnownSection* known = nullptr;
        for (const KnownSection& k : kKnown) {
            if (k.flag == *section) {
                known = &k;
            }
        }
        const std::u32string* line = nullptr;
        if (known) {
            if (known->pointer >= pointers.size()) {
                throw ValueError("list index out of range");
            }
            const std::int64_t num = pointers[known->pointer] * known->perRecord;
            numlines = floorDiv(num, known->perLine) + (num % known->perLine != 0 ? 1 : 0);
            readSection(*known, *numlines, attrs);
            if (!atEnd()) {
                line = &next();
                if (*numlines == 0 && strip(*line).empty()) {
                    line = atEnd() ? nullptr : &next();
                }
            }
        } else {
            while (!atEnd()) {
                line = &next();
                if (startsWith(*line, U"%FLAG")) {
                    break;
                }
                line = nullptr;
            }
            if (line) {
                if (!numlines) {
                    throw ValueError(
                        "cannot access local variable 'numlines' where it is not associated with a "
                        "value");
                }
                if (*numlines == 0 && strip(*line).empty()) {
                    line = atEnd() ? nullptr : &next();
                }
            }
        }
        if (!line) {
            break;
        }
        section = flagName(*line);
        if (!section) {
            throw ValueError("%FLAG section not found, formatting error for PARM7 file " + m_path +
                             " ");
        }
    }

    const auto require = [&](Attr attr) {
        if (!attrs.has(attr)) {
            throw ValueError(std::string("'") + attrName(attr) + "'");  // KeyError
        }
    };
    require(Attr::Names);
    const std::size_t nAtoms = attrs.size(Attr::Names);
    require(Attr::ResPointers);
    std::vector<std::int64_t> pointersOfResidues = attrs.ints[Attr::ResPointers];
    pointersOfResidues.push_back(static_cast<std::int64_t>(nAtoms));
    std::vector<std::int32_t> residx(nAtoms, 0);
    const auto clampIndex = [&](std::int64_t x) {
        const auto n = static_cast<std::int64_t>(nAtoms);
        if (x < 0) {
            x = std::max<std::int64_t>(x + n, 0);
        }
        return std::min(x, n);
    };
    for (std::size_t r = 0; r + 1 < pointersOfResidues.size(); ++r) {
        const std::int64_t from = clampIndex(pointersOfResidues[r]);
        const std::int64_t to = clampIndex(pointersOfResidues[r + 1]);
        for (std::int64_t a = from; a < to; ++a) {
            residx[static_cast<std::size_t>(a)] = static_cast<std::int32_t>(r);
        }
    }
    require(Attr::Resnames);
    const std::size_t nResidues = attrs.size(Attr::Resnames);
    for (const Attr attr :
         {Attr::BondA, Attr::BondH, Attr::AngA, Attr::AngH, Attr::DihA, Attr::DihH}) {
        require(attr);
    }

    // the Universe's length checks, in the order of the attributes
    for (const Attr attr : attrs.order) {
        const std::string name = topologyAttrName(attr);
        if (!name.empty() && attrs.size(attr) != nAtoms) {
            throw ValueError("Length of " + name + " does not match number of atoms.\nExpect: " +
                             std::to_string(nAtoms) + " Have: " + std::to_string(attrs.size(attr)));
        }
    }

    auto topology = std::make_shared<core::Topology>();
    core::StringPool& pool = topology->strings;
    const auto intern = [&](const std::vector<std::string>& values) {
        std::vector<std::uint32_t> out;
        out.reserve(values.size());
        for (const std::string& v : values) {
            out.push_back(pool.intern(v));
        }
        return out;
    };
    topology->names = intern(attrs.strings[Attr::Names]);
    if (attrs.has(Attr::Types)) {
        topology->types = intern(attrs.strings[Attr::Types]);
    }
    if (attrs.has(Attr::Elements)) {
        topology->elements = intern(attrs.strings[Attr::Elements]);
    }
    if (attrs.has(Attr::Charges)) {
        topology->charges = attrs.reals[Attr::Charges];
    }
    if (attrs.has(Attr::Masses)) {
        topology->masses = attrs.reals[Attr::Masses];
    }
    topology->atomResidue = residx;
    std::vector<std::int64_t> ids(nAtoms);
    for (std::size_t a = 0; a < nAtoms; ++a) {
        ids[a] = static_cast<std::int64_t>(a) + 1;
    }
    topology->ids = std::move(ids);
    topology->resnames = intern(attrs.strings[Attr::Resnames]);
    for (std::size_t r = 0; r < nResidues; ++r) {
        topology->resids.push_back(static_cast<std::int64_t>(r) + 1);
    }
    topology->resnums = topology->resids;

    const std::vector<std::string>* chains =
        attrs.has(Attr::ChainIds) ? &attrs.strings[Attr::ChainIds] : nullptr;
    if (chains && chains->size() == nResidues) {
        topology->residueSegment.resize(nResidues);
        for (std::size_t r = 0; r < nResidues; ++r) {
            if (r == 0 || (*chains)[r] != (*chains)[r - 1]) {
                topology->segids.push_back(pool.intern((*chains)[r]));
            }
            topology->residueSegment[r] = static_cast<std::int32_t>(topology->segids.size()) - 1;
        }
        std::vector<std::uint32_t> chainIDs;
        for (const std::int32_t r : residx) {
            chainIDs.push_back(pool.intern((*chains)[static_cast<std::size_t>(r)]));
        }
        topology->chainIDs = std::move(chainIDs);
    } else {
        topology->segids.push_back(pool.intern("SYSTEM"));
        topology->residueSegment.assign(nResidues, 0);
    }

    std::vector<core::Bond> bonds;
    for (const Attr attr : {Attr::BondA, Attr::BondH}) {
        const std::vector<std::int64_t>& atoms = attrs.ints[attr];
        for (std::size_t k = 0; k + 1 < atoms.size(); k += 2) {
            if (atoms[k] != atoms[k + 1]) {  // core::Topology has no bonds of an atom to itself
                bonds.push_back(
                    {static_cast<std::int32_t>(atoms[k]), static_cast<std::int32_t>(atoms[k + 1])});
            }
        }
    }
    topology->bonds = std::move(bonds);
    try {
        topology->finish();
    } catch (const std::invalid_argument& e) {
        throw ValueError(m_path + ": " + e.what());
    }
    return topology;
}

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
    return TopReader(path).parse();
}

bool isParm7(std::string_view path) {
    std::string ext = core::posixpath::splitext(path).second;
    for (char& c : ext) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return ext == ".parm7" || ext == ".prmtop" || ext == ".top";
}

}  // namespace mdbridge::io
