// PDB and GRO files as MDAnalysis 2.11.0 reads them for MDATK: the topologies of its PDBParser and
// GROParser (with the types and masses the Universe then guesses) and the coordinates of its
// PDBReader (models as frames, CRYST1, HEADER / COMPND / REMARK) and GROReader.

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

#include "pdb_tables.inc"

using core::ValueError;

// --- MDAnalysis' guessers (DefaultGuesser.guess_atom_element / get_atom_mass)

std::optional<std::string_view> atomElementEntry(const std::string& key) {
    for (const auto& [name, element] : kAtomElements) {
        if (name == key) {
            return element;
        }
    }
    return std::nullopt;
}

bool guessElement(const std::string& name) {
    return std::find(std::begin(kGuessElements), std::end(kGuessElements), name) !=
           std::end(kGuessElements);
}

std::string guessAtomElement(const std::string& atomName) {
    if (atomName.empty()) {
        return "";
    }
    const U32 name32 = toU32(atomName);
    if (const auto element = atomElementEntry(toUtf8(upper(name32)))) {
        return std::string(*element);
    }
    U32 noSymbols;
    for (const char32_t c : name32) {
        if (c != U'*' && c != U'+' && c != U'-') {
            noSymbols += c;
        }
    }
    U32 bare;
    for (const char32_t c : noSymbols) {
        if (c < U'0' || c > U'9') {
            bare += c;
        }
    }
    U32 name = upper(bare);
    if (const auto element = atomElementEntry(toUtf8(name))) {
        return std::string(*element);
    }
    while (!name.empty()) {
        if (guessElement(toUtf8(name))) {
            return toUtf8(name);
        }
        const U32 head = name.substr(0, name.size() - 1);
        if (guessElement(toUtf8(head))) {
            return toUtf8(head);
        }
        const U32 tail = name.substr(1);
        if (guessElement(toUtf8(tail))) {
            return toUtf8(tail);
        }
        if (name.size() <= 2) {
            return toUtf8(name.substr(0, 1));
        }
        name.pop_back();
    }
    return toUtf8(noSymbols);
}

double atomMass(const std::string& element) {
    const auto find = [](const std::string& key) -> std::optional<double> {
        for (const auto& [name, mass] : kMasses) {
            if (name == key) {
                return mass;
            }
        }
        return std::nullopt;
    };
    if (const auto mass = find(element)) {
        return *mass;
    }
    if (const auto mass = find(toUtf8(upper(toU32(element))))) {
        return *mass;
    }
    return 0.0;
}

// elem.capitalize() when the result is ASCII (element symbols are), else std::nullopt
std::optional<std::string> capitalizeAscii(U32View text) {
    std::string out;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char32_t c = text[i];
        if (i == 0) {
            if (c < 0x80) {
                out += static_cast<char>(c >= U'a' && c <= U'z' ? c - 32 : c);
                continue;
            }
            static const std::array<std::pair<char32_t, const char*>, 10> kTitle = {{
                {0xDF, "Ss"},
                {0x131, "I"},
                {0x17F, "S"},
                {0xFB00, "Ff"},
                {0xFB01, "Fi"},
                {0xFB02, "Fl"},
                {0xFB03, "Ffi"},
                {0xFB04, "Ffl"},
                {0xFB05, "St"},
                {0xFB06, "St"},
            }};
            const auto it = std::find_if(kTitle.begin(), kTitle.end(),
                                         [c](const auto& e) { return e.first == c; });
            if (it == kTitle.end()) {
                return std::nullopt;
            }
            out += it->second;
            continue;
        }
        if (c < 0x80) {
            out += static_cast<char>(c >= U'A' && c <= U'Z' ? c + 32 : c);
        } else if (c == 0x212A) {  // KELVIN SIGN
            out += 'k';
        } else {
            return std::nullopt;
        }
    }
    return out;
}

bool knownElement(const std::string& symbol) {
    return std::find(std::begin(kElementSymbols), std::end(kElementSymbols), symbol) !=
           std::end(kElementSymbols);
}

// int(text) raising Python's ValueError
std::int64_t intOrThrow(U32View text) {
    return core::pythonInt(toUtf8(text));
}

// MDAnalysis' anyopen takes an empty file for a truncated bzip2 stream
std::string readNonEmpty(const std::filesystem::path& path) {
    std::string bytes = readBytes(path);
    if (bytes.empty()) {
        throw ReaderError("EOFError",
                          "Compressed file ended before the end-of-stream marker was reached");
    }
    return bytes;
}

// The atom count of a GRO file's second line (numpy refuses negative array sizes)
std::size_t atomCount(U32View line) {
    const std::int64_t count = intOrThrow(line);
    if (count < 0) {
        throw ValueError("negative dimensions are not allowed");
    }
    return static_cast<std::size_t>(count);
}

U32 readStrictLines(const std::filesystem::path& path) {
    return universalNewlines(decodeStrict(readNonEmpty(path)));
}

// --- PDB topology

constexpr std::u32string_view kUpperDigits = U"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
constexpr std::u32string_view kLowerDigits = U"0123456789abcdefghijklmnopqrstuvwxyz";

// hy36decode(5, s): std::nullopt for MDAnalysis' ValueError; a character outside the digits of
// the chosen base raises KeyError
std::optional<std::int64_t> hy36Decode(U32View s) {
    if (s.size() != 5) {
        return std::nullopt;
    }
    const char32_t f = s[0];
    if (f == U'-' || f == U' ' || isDigit(f)) {
        return toInt(s);
    }
    const auto decode = [&](std::u32string_view digits) {
        std::int64_t result = 0;
        for (const char32_t c : s) {
            const std::size_t value = digits.find(c);
            if (value == std::u32string_view::npos) {
                throw ReaderError("KeyError", core::pythonStringRepr(toUtf8(U32(1, c))));
            }
            result = result * 36 + static_cast<std::int64_t>(value);
        }
        return result;
    };
    if (kUpperDigits.find(f) != std::u32string_view::npos) {
        return decode(kUpperDigits) - 10 * 36 * 36 * 36 * 36 + 100000;
    }
    if (kLowerDigits.find(f) != std::u32string_view::npos) {
        return decode(kLowerDigits) + 16 * 36 * 36 * 36 * 36 + 100000;
    }
    return std::nullopt;
}

float float32(double value) {
    return static_cast<float>(value);
}

struct PdbRows {
    std::vector<std::string> records, names, altLocs, resnames, chainIDs, elements, charges, icodes,
        segids;
    std::vector<std::int64_t> serials, resids;
    std::vector<double> occupancies, tempFactors;
    bool wrappedSerials = false;
};

PdbRows parsePdbAtoms(const std::vector<U32>& lines) {
    PdbRows rows;
    std::int64_t residPrevious = 0;
    std::int64_t lastWrapped = 100000;
    for (const U32& raw : lines) {
        const U32 line = strip(raw);
        if (line.empty()) {
            continue;
        }
        if (line.starts_with(U"END")) {
            break;
        }
        if (!line.starts_with(U"ATOM") && !line.starts_with(U"HETATM")) {
            continue;
        }
        rows.records.push_back(toUtf8(strip(slice(line, 0, 6))));
        const U32 serialText = slice(line, 6, 11);
        std::optional<std::int64_t> serial = toInt(serialText);
        if (!serial) {
            try {
                serial = hy36Decode(serialText);
            } catch (const ReaderError&) {
                // the KeyError leaves `serial` unbound for the finally clause that appends it
                throw ReaderError("UnboundLocalError",
                                  "cannot access local variable 'serial' where it is not "
                                  "associated with a value");
            }
        }
        if (!serial) {
            rows.wrappedSerials = true;
            serial = lastWrapped++;
        }
        rows.serials.push_back(static_cast<std::int32_t>(*serial));
        rows.names.push_back(toUtf8(strip(slice(line, 12, 16))));
        rows.altLocs.push_back(toUtf8(strip(slice(line, 16, 17))));
        rows.resnames.push_back(toUtf8(strip(slice(line, 17, 21))));
        rows.chainIDs.push_back(toUtf8(strip(slice(line, 21, 22))));
        rows.elements.push_back(toUtf8(strip(slice(line, 76, 78))));
        rows.charges.push_back(toUtf8(strip(slice(line, 78, 80))));
        if (auto resid = toInt(slice(line, 22, 26))) {
            while (*resid - residPrevious < -5000) {
                *resid += 10000;
            }
            residPrevious = *resid;
            rows.resids.push_back(*resid);
        } else {
            rows.resids.push_back(1);
        }
        rows.icodes.push_back(toUtf8(strip(slice(line, 26, 27))));
        rows.occupancies.push_back(float32(toFloat(slice(line, 54, 60)).value_or(0.0)));
        rows.tempFactors.push_back(float32(toFloat(slice(line, 60, 66)).value_or(1.0)));
        rows.segids.push_back(toUtf8(strip(slice(line, 72, 76))));
    }
    return rows;
}

// formal charges "2+" -> 2, "1-" -> -1, "" / "0" -> 0; std::nullopt (attribute dropped) otherwise
std::optional<std::vector<std::int64_t>> formalCharges(const std::vector<std::string>& entries) {
    std::vector<std::int64_t> out;
    for (const std::string& entry : entries) {
        if (entry.empty() || entry == "0") {
            out.push_back(0);
            continue;
        }
        if (entry.find_first_of("+-") == std::string::npos) {
            return std::nullopt;
        }
        U32 reversed = toU32(entry);
        std::reverse(reversed.begin(), reversed.end());
        const auto value = toInt(reversed);
        if (!value) {
            return std::nullopt;
        }
        out.push_back(*value);
    }
    return out;
}

// _parsebonds: CONECT records anywhere in the file; std::nullopt for no bonds
std::optional<std::vector<core::Bond>> pdbBonds(const std::vector<U32>& lines,
                                                const std::vector<std::int64_t>& serials) {
    std::map<std::int64_t, std::int32_t> mapping;
    for (std::size_t i = 0; i < serials.size(); ++i) {
        mapping[serials[i]] = static_cast<std::int32_t>(i);
    }
    std::set<std::pair<std::int32_t, std::int32_t>> bonds;
    for (const U32& raw : lines) {
        if (slice(raw, 0, 6) != U"CONECT") {
            continue;
        }
        const U32 conect = strip(raw);
        const std::int64_t atom = intOrThrow(slice(conect, 6, 11));
        const U32 rest = slice(conect, 11);
        const std::size_t count = rest.size() / 5;
        if (count == 0) {
            continue;
        }
        if (rest.size() % count != 0) {
            return std::nullopt;  // corrupt CONECT records: no bonds at all
        }
        for (std::size_t i = 0; i < count; ++i) {
            const std::int64_t other = intOrThrow(slice(rest, 5 * i, 5 * i + 5));
            const auto a = mapping.find(atom);
            const auto b = mapping.find(other);
            if (a != mapping.end() && b != mapping.end()) {
                bonds.emplace(a->second, b->second);
            }
        }
    }
    if (bonds.empty()) {
        return std::nullopt;
    }
    std::vector<core::Bond> out;
    for (const auto& [a, b] : bonds) {
        if (a != b) {
            out.push_back({std::min(a, b), std::max(a, b)});
        }
    }
    return out;
}

std::shared_ptr<core::Topology> pdbTopology(const std::vector<U32>& lines) {
    PdbRows rows = parsePdbAtoms(lines);
    auto top = std::make_shared<core::Topology>();
    core::StringPool& pool = top->strings;
    const std::size_t n = rows.names.size();
    const auto intern = [&](const std::vector<std::string>& values) {
        std::vector<std::uint32_t> out;
        out.reserve(values.size());
        for (const std::string& v : values) {
            out.push_back(pool.intern(v));
        }
        return out;
    };
    const bool anySegid = std::any_of(rows.segids.begin(), rows.segids.end(),
                                      [](const auto& s) { return !s.empty(); });
    std::vector<std::string> segids = anySegid ? rows.segids : rows.chainIDs;

    top->names = intern(rows.names);
    top->altLocs = intern(rows.altLocs);
    top->chainIDs = intern(rows.chainIDs);
    top->recordTypes = intern(rows.records);
    top->ids = rows.serials;
    top->tempFactors = rows.tempFactors;
    top->occupancies = rows.occupancies;

    const bool anyElement = std::any_of(rows.elements.begin(), rows.elements.end(),
                                        [](const auto& e) { return !e.empty(); });
    std::vector<std::string> types(n);
    std::vector<double> masses(n);
    if (anyElement) {
        std::vector<std::string> validated(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto symbol = capitalizeAscii(toU32(rows.elements[i]));
            validated[i] = symbol && knownElement(*symbol) ? *symbol : "";
            types[i] = rows.elements[i];  // the raw column; blank entries are not guessed
            masses[i] = atomMass(validated[i]);
        }
        top->elements = intern(validated);
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            types[i] = guessAtomElement(rows.names[i]);
            masses[i] = atomMass(types[i]);
        }
    }
    top->types = intern(types);
    top->masses = std::move(masses);
    if (std::any_of(rows.charges.begin(), rows.charges.end(),
                    [](const auto& c) { return !c.empty(); })) {
        top->formalCharges = formalCharges(rows.charges);
    }

    // residues: consecutive runs of equal (resid, resname, icode, segid)
    std::vector<std::string> residueSegids;
    for (std::size_t i = 0; i < n; ++i) {
        const bool same = i > 0 && rows.resids[i] == rows.resids[i - 1] &&
                          rows.resnames[i] == rows.resnames[i - 1] &&
                          rows.icodes[i] == rows.icodes[i - 1] && segids[i] == segids[i - 1];
        if (!same) {
            top->resids.push_back(static_cast<std::int32_t>(rows.resids[i]));
            top->resnames.push_back(pool.intern(rows.resnames[i]));
            if (!top->icodes) {
                top->icodes.emplace();
            }
            top->icodes->push_back(pool.intern(rows.icodes[i]));
            residueSegids.push_back(segids[i]);
        }
        top->atomResidue.push_back(static_cast<std::int32_t>(top->resids.size() - 1));
    }
    if (!top->icodes) {
        top->icodes.emplace();
    }
    top->resnums = top->resids;
    if (std::any_of(residueSegids.begin(), residueSegids.end(),
                    [](const auto& s) { return !s.empty(); })) {
        for (std::size_t r = 0; r < residueSegids.size(); ++r) {
            if (r == 0 || residueSegids[r] != residueSegids[r - 1]) {
                top->segids.push_back(pool.intern(residueSegids[r]));
            }
            top->residueSegment.push_back(static_cast<std::int32_t>(top->segids.size() - 1));
        }
    } else {
        top->segids.push_back(pool.intern("SYSTEM"));
        top->residueSegment.assign(residueSegids.size(), 0);
    }
    if (!rows.wrappedSerials) {
        top->bonds = pdbBonds(lines, rows.serials);
    }
    top->finish();
    return top;
}

// --- PDB coordinates

std::string bytesStrip(std::string_view bytes) {
    const auto space = [](char c) {
        return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
    };
    std::size_t first = 0;
    std::size_t last = bytes.size();
    while (first < last && space(bytes[first])) {
        ++first;
    }
    while (last > first && space(bytes[last - 1])) {
        --last;
    }
    return std::string(bytes.substr(first, last - first));
}

std::string byteSlice(std::string_view bytes, std::size_t first,
                      std::size_t last = std::string_view::npos) {
    if (first >= bytes.size()) {
        return {};
    }
    return std::string(bytes.substr(first, last == std::string_view::npos ? last : last - first));
}

std::string decodedText(std::string_view bytes) {
    return toUtf8(decodeStrict(bytes));
}

std::vector<U32> splitlines(U32View text) {
    std::vector<U32> lines;
    std::size_t first = 0;
    std::size_t i = 0;
    while (i < text.size()) {
        const char32_t c = text[i];
        const bool lineBreak = c == U'\n' || c == U'\r' || c == 0x0B || c == 0x0C || c == 0x1C ||
                               c == 0x1D || c == 0x1E || c == 0x85 || c == 0x2028 || c == 0x2029;
        if (!lineBreak) {
            ++i;
            continue;
        }
        lines.emplace_back(text.substr(first, i - first));
        ++i;
        if (c == U'\r' && i < text.size() && text[i] == U'\n') {
            ++i;
        }
        first = i;
    }
    if (first < text.size()) {
        lines.emplace_back(text.substr(first));
    }
    return lines;
}

float coordinate(U32View text) {
    const auto value = toFloat(text);
    if (!value) {
        throw ValueError("could not convert string to float: " +
                         core::pythonStringRepr(toUtf8(text)));
    }
    return static_cast<float>(*value);
}

class PdbTrajectory : public FileTrajectory {
public:
    PdbTrajectory(const std::filesystem::path& path, std::optional<double> dt, std::size_t atoms)
        : m_path(pathText(path)), m_atoms(atoms) {
        m_bytes = readNonEmpty(path);
        if (m_atoms == 0) {
            m_atoms =
                parsePdbAtoms(fileLines(universalNewlines(decodeStrict(m_bytes)))).names.size();
        }
        m_dt = dt.value_or(1.0);
        std::vector<std::size_t> models;
        std::vector<std::size_t> crysts;
        std::size_t at = 0;
        while (at < m_bytes.size()) {
            const std::size_t newline = m_bytes.find('\n', at);
            const std::size_t end = newline == std::string::npos ? m_bytes.size() : newline + 1;
            const std::string_view line(m_bytes.data() + at, end - at);
            if (line.starts_with("MODEL")) {
                models.push_back(end);
            } else if (line.starts_with("CRYST")) {
                crysts.push_back(at);
            } else if (line.starts_with("HEADER")) {
                m_info.header = decodedText(bytesStrip(byteSlice(line, 10, 66)));
            } else if (line.starts_with("TITLE")) {
                decodedText(bytesStrip(byteSlice(line, 8, 80)));
            } else if (line.starts_with("COMPND")) {
                m_info.compounds->push_back(decodedText(bytesStrip(byteSlice(line, 7, 80))));
            } else if (line.starts_with("REMARK")) {
                m_info.remarks->push_back(decodedText(bytesStrip(byteSlice(line, 6))));
            }
            at = end;
        }
        if (models.empty()) {
            models.push_back(0);
        }
        m_starts = models;
        if (crysts.size() == models.size()) {
            for (std::size_t k = 0; k < models.size(); ++k) {
                m_starts[k] = std::min(models[k], crysts[k]);
            }
        }
        m_stops.assign(m_starts.begin() + 1, m_starts.end());
        m_stops.push_back(m_bytes.size());
        m_frameCount = static_cast<std::int64_t>(m_starts.size());
        seek(0);
    }

    core::TrajectoryInfo info() const override { return m_info; }
    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        const auto k = static_cast<std::size_t>(index);
        const std::size_t start = std::min(m_starts[k], m_bytes.size());
        const std::size_t stop = std::max(start, std::min(m_stops[k], m_bytes.size()));
        const U32 chunk = decodeStrict(std::string_view(m_bytes).substr(start, stop - start));
        std::vector<float> positions;
        std::size_t count = 0;
        for (const U32& line : splitlines(chunk)) {
            const U32 record = slice(line, 0, 6);
            if (record == U"ATOM  " || record == U"HETATM") {
                if (count >= m_atoms) {
                    throw ReaderError("IndexError", "index " + std::to_string(count) +
                                                        " is out of bounds for axis 0 with size " +
                                                        std::to_string(m_atoms));
                }
                positions.push_back(coordinate(slice(line, 30, 38)));
                positions.push_back(coordinate(slice(line, 38, 46)));
                positions.push_back(coordinate(slice(line, 46, 54)));
                ++count;
            } else if (record == U"CRYST1") {
                std::array<float, 6> cell{};
                bool ok = true;
                for (std::size_t i = 0; i < 6 && ok; ++i) {
                    static constexpr std::array<std::size_t, 7> kColumns = {6,  15, 24, 33,
                                                                            40, 47, 54};
                    const auto value = toFloat(slice(line, kColumns[i], kColumns[i + 1]));
                    ok = value.has_value();
                    cell[i] = ok ? static_cast<float>(*value) : 0.0f;
                }
                if (!ok || placeholder(cell)) {
                    m_box.reset();
                } else {
                    m_box = cell;
                }
            }
        }
        if (count != m_atoms) {
            throw ValueError(
                "Inconsistency in file '" + m_path + "': The number of atoms (" +
                std::to_string(count) + ") in trajectory frame " + std::to_string(index) +
                " differs from the number of atoms (" + std::to_string(m_atoms) +
                ") in the corresponding topology.\nTrajectories with varying numbers of atoms are "
                "currently not supported.");
        }
        frame.positions = std::move(positions);
        frame.box = m_box;
        if (m_box && (*m_box)[0] == 0.0f && (*m_box)[1] == 0.0f && (*m_box)[2] == 0.0f) {
            frame.box.reset();
        }
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    // np.allclose(cell, [1, 1, 1, 90, 90, 90])
    static bool placeholder(const std::array<float, 6>& cell) {
        static constexpr std::array<double, 6> kUnit = {1.0, 1.0, 1.0, 90.0, 90.0, 90.0};
        for (std::size_t i = 0; i < 6; ++i) {
            const double v = cell[i];
            if (!(std::abs(v - kUnit[i]) <= 1e-8 + 1e-5 * kUnit[i])) {
                return false;
            }
        }
        return true;
    }

    std::string m_path;
    std::size_t m_atoms;
    std::string m_bytes;
    std::vector<std::size_t> m_starts;
    std::vector<std::size_t> m_stops;
    core::TrajectoryInfo m_info{std::string(), std::vector<std::string>(),
                                std::vector<std::string>()};
    std::optional<std::array<float, 6>> m_box;  // kept from frame to frame, as MDAnalysis' Timestep
};

// --- GRO

class GroTrajectory : public FileTrajectory {
public:
    GroTrajectory(const std::filesystem::path& path, std::optional<double> dt) {
        m_dt = dt.value_or(1.0);
        m_frameCount = 1;
        const std::vector<U32> lines = fileLines(readStrictLines(path));
        if (lines.size() < 2) {
            throw ValueError("invalid literal for int() with base 10: ''");
        }
        const std::size_t atoms = atomCount(lines[1]);
        m_positions.assign(3 * atoms, 0.0f);
        const U32 first = lines.size() > 2 ? lines[2] : U32();
        const std::size_t dot = slice(first, 25).find(U'.');
        const std::size_t cs = dot == U32::npos ? 0 : dot + 1;
        const auto setPosition = [&](std::size_t atom, const U32& line) {
            if (atom >= atoms) {
                throw ReaderError("IndexError", "index " + std::to_string(atom) +
                                                    " is out of bounds for axis 0 with size " +
                                                    std::to_string(atoms));
            }
            for (std::size_t i = 0; i < 3; ++i) {
                m_positions[3 * atom + i] = coordinate(slice(line, 20 + cs * i, 20 + cs * (i + 1)));
            }
        };
        setPosition(0, first);
        std::optional<std::vector<float>> cell;
        for (std::size_t pos = 1; 2 + pos < lines.size(); ++pos) {
            const U32& line = lines[2 + pos];
            if (pos == atoms) {
                cell = unitCell(line);
                break;
            }
            setPosition(pos, line);
        }
        if (!cell) {
            throw ReaderError(
                "UnboundLocalError",
                "cannot access local variable 'unitcell' where it is not associated with a value");
        }
        for (float& x : m_positions) {
            x *= 10.0f;
        }
        if (cell->size() == 3) {
            const bool zero = std::all_of(cell->begin(), cell->end(),
                                          [](float v) { return std::abs(double(v)) <= 1e-8; });
            if (!zero) {
                m_box =
                    std::array<float, 6>{(*cell)[0], (*cell)[1], (*cell)[2], 90.0f, 90.0f, 90.0f};
            }
        } else if (cell->size() == 9) {
            const auto& c = *cell;
            m_box = triclinicBox({c[0], c[3], c[4]}, {c[5], c[1], c[6]}, {c[7], c[8], c[2]});
        } else {
            throw ValueError("GRO unitcell has neither 3 nor 9 entries.");
        }
        if (m_box) {
            for (std::size_t i = 0; i < 3; ++i) {
                (*m_box)[i] *= 10.0f;
            }
            if ((*m_box)[0] == 0.0f && (*m_box)[1] == 0.0f && (*m_box)[2] == 0.0f) {
                m_box.reset();
            }
        }
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t, core::Frame& frame) override {
        frame.positions = m_positions;
        frame.box = m_box;
        frame.time = 0.0;
    }

private:
    static std::vector<float> unitCell(const U32& line) {
        std::vector<float> out;
        bool ok = true;
        for (const U32& word : split(line)) {
            const auto value = toFloat(word);
            if (!value) {
                ok = false;
                break;
            }
            out.push_back(static_cast<float>(*value));
        }
        if (ok) {
            return out;
        }
        // re.findall(r"(\d+\.\d{5})", line)
        out.clear();
        std::size_t i = 0;
        while (i < line.size()) {
            std::size_t j = i;
            while (j < line.size() && isDigit(line[j])) {
                ++j;
            }
            if (j > i && j < line.size() && line[j] == U'.') {
                std::size_t k = j + 1;
                while (k < line.size() && k - j - 1 < 5 && isDigit(line[k])) {
                    ++k;
                }
                if (k - j - 1 == 5) {
                    out.push_back(static_cast<float>(*toFloat(slice(line, i, k))));
                    i = k;
                    continue;
                }
            }
            ++i;
        }
        return out;
    }

    std::vector<float> m_positions;
    std::optional<std::array<float, 6>> m_box;
};

std::shared_ptr<core::Topology> groTopology(const std::filesystem::path& path) {
    const std::vector<U32> lines = fileLines(readStrictLines(path));
    if (lines.size() < 2) {
        throw ReaderError("StopIteration", "");
    }
    const std::size_t atoms = atomCount(lines[1]);
    std::vector<std::int64_t> resids(atoms, 0);
    std::vector<std::string> resnames(atoms);
    std::vector<std::string> names(atoms);
    std::vector<std::int64_t> indices(atoms, 0);
    std::size_t read = 0;
    for (std::size_t i = 0; 2 + i < lines.size() && i < atoms; ++i) {
        const U32& line = lines[2 + i];
        const auto resid = toInt(slice(line, 0, 5));
        const auto index = toInt(slice(line, 15, 20));
        if (!resid || !index) {
            throw ReaderError(
                "OSError", "Couldn't read the following line of the .gro file:\n" + toUtf8(line));
        }
        resids[i] = static_cast<std::int32_t>(*resid);
        resnames[i] = toUtf8(strip(slice(line, 5, 10)));
        names[i] = toUtf8(strip(slice(line, 10, 15)));
        indices[i] = static_cast<std::int32_t>(*index);
        ++read;
    }
    for (std::size_t i = 0; i < read; ++i) {
        if (names[i].empty()) {
            throw ReaderError("OSError", "Missing atom name on line: " + std::to_string(i + 3));
        }
    }
    if (read < atoms) {
        throw ReaderError("IndexError", "index 0 is out of bounds for axis 0 with size 0");
    }
    // resid wraps: each block of consecutive zero resids, but one at the first atom, adds 100000
    // from there on
    std::vector<std::size_t> starts;
    std::optional<std::size_t> previousZero;
    for (std::size_t i = 0; i < atoms; ++i) {
        if (resids[i] != 0) {
            continue;
        }
        if (!previousZero || *previousZero + 1 != i) {
            starts.push_back(i);
        }
        previousZero = i;
    }
    if (!starts.empty() && starts.front() == 0) {
        starts.erase(starts.begin());
    }
    for (const std::size_t s : starts) {
        for (std::size_t i = s; i < atoms; ++i) {
            resids[i] = static_cast<std::int32_t>(resids[i] + 100000);
        }
    }
    auto top = std::make_shared<core::Topology>();
    core::StringPool& pool = top->strings;
    std::vector<std::string> types(atoms);
    std::vector<double> masses(atoms);
    for (std::size_t i = 0; i < atoms; ++i) {
        top->names.push_back(pool.intern(names[i]));
        types[i] = guessAtomElement(names[i]);
        top->types.push_back(pool.intern(types[i]));
        masses[i] = atomMass(types[i]);
        if (i == 0 || resids[i] != resids[i - 1] || resnames[i] != resnames[i - 1]) {
            top->resids.push_back(resids[i]);
            top->resnames.push_back(pool.intern(resnames[i]));
        }
        top->atomResidue.push_back(static_cast<std::int32_t>(top->resids.size() - 1));
    }
    top->masses = std::move(masses);
    top->ids = indices;
    top->resnums = top->resids;
    top->segids.push_back(pool.intern("SYSTEM"));
    top->residueSegment.assign(top->resids.size(), 0);
    top->finish();
    return top;
}

// The Universe's wrapping of a parser's ValueError / OSError
template <class Fn>
std::shared_ptr<core::Topology> parseTopology(const std::filesystem::path& path, const char* parser,
                                              Fn&& fn) {
    const std::string name = std::string("<class '") + parser + "'>";
    try {
        return fn();
    } catch (const ValueError& error) {
        const std::string python = error.pythonName();
        if (python == "ValueError" || python == "UnicodeDecodeError") {
            throw ValueError("Failed to construct topology from file " + pathText(path) +
                             " with parser " + name + ".\nError: " + error.what());
        }
        if (python == "OSError") {
            throw ReaderError("OSError", "Failed to load from the topology file " + pathText(path) +
                                             " with parser " + name + ".\nError: " + error.what());
        }
        throw;
    }
}

}  // namespace

TopologyFile readPdbTopology(const std::filesystem::path& path, std::optional<double> dt) {
    TopologyFile out;
    out.topology = parseTopology(path, "MDAnalysis.topology.PDBParser.PDBParser", [&] {
        auto top = pdbTopology(fileLines(readStrictLines(path)));
        if (top->atomCount() == 0) {
            throw ReaderError("IndexError", "index 0 is out of bounds for axis 0 with size 0");
        }
        return top;
    });
    out.coordinates = openPdbCoordinates(path, dt, out.topology->atomCount());
    return out;
}

TopologyFile readGroTopology(const std::filesystem::path& path, std::optional<double> dt) {
    TopologyFile out;
    out.topology = parseTopology(path, "MDAnalysis.topology.GROParser.GROParser",
                                 [&] { return groTopology(path); });
    out.coordinates = openGroCoordinates(path, dt, out.topology->atomCount());
    return out;
}

std::unique_ptr<core::Trajectory> openPdbCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t atoms) {
    return std::make_unique<PdbTrajectory>(path, dt, atoms);
}

std::unique_ptr<core::Trajectory> openGroCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t) {
    return std::make_unique<GroTrajectory>(path, dt);
}

}  // namespace mdbridge::io::detail
