// PDB and GRO files as MDAnalysis 2.11.0 reads them for MDATK: the topologies of its PDBParser and
// GROParser (with the types and masses the Universe then guesses) and the coordinates of its
// PDBReader (models as frames, CRYST1, HEADER / COMPND / REMARK) and GROReader. The guessing
// follows MDAnalysis' procedure on MDBridge's own element / name tables (data_tables.hpp).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>

#include "mdbridge/core/detail/data_tables.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

namespace tables = core::detail::tables;

constexpr const char* kPdbParser = "<class 'MDAnalysis.topology.PDBParser.PDBParser'>";
constexpr const char* kGroParser = "<class 'MDAnalysis.topology.GROParser.GROParser'>";

// ---------------------------------------------------------------------------------------------
// Errors

ReaderError emptyFileError() {
    return ReaderError("EOFError", "Compressed file ended before the end-of-stream marker was reached");
}

ReaderError outOfBounds(std::size_t index, std::size_t size) {
    return ReaderError("IndexError", "index " + std::to_string(index) +
                                         " is out of bounds for axis 0 with size " +
                                         std::to_string(size));
}

ReaderError unboundLocal(const std::string& variable) {
    return ReaderError("UnboundLocalError", "cannot access local variable '" + variable +
                                                "' where it is not associated with a value");
}

// The Universe's rewrapping of parser errors: ValueError / TypeError (UnicodeDecodeError is a
// ValueError) become "Failed to construct topology ..." ValueErrors; other errors pass through.
template <class Fn>
auto wrapParserErrors(const std::string& fileText, const char* parser, Fn&& parse) {
    try {
        return parse();
    } catch (const core::ValueError& error) {
        const std::string_view name = error.pythonName();
        if (name != "ValueError" && name != "UnicodeDecodeError" && name != "TypeError") {
            throw;
        }
        throw core::ValueError("Failed to construct topology from file " + fileText +
                               " with parser " + parser + ".\nError: " + error.what());
    }
}

// An OSError of a parser, as the Universe rewraps it
ReaderError topologyLoadError(const std::string& fileText, const char* parser,
                              const std::string& message) {
    return ReaderError("OSError", "Failed to load from the topology file " + fileText +
                                      " with parser " + parser + ".\nError: " + message);
}

// ---------------------------------------------------------------------------------------------
// Types and masses the Universe guesses
//
// The procedure is MDAnalysis' (guess_atom_type / guess_masses); the data are the independent
// tables of mdbridge/core/detail/data_tables.hpp, so names those tables treat differently get
// other types and masses than MDAnalysis gives them.

// What the force-field name table says about an atom of a residue: its element, or "" for a
// virtual site. Names and residues compare in any letter case; "*water*" entries apply to every
// water residue.
std::optional<std::string_view> forceFieldElement(std::string_view name, std::string_view resname) {
    const bool water = tables::contains(tables::kIndWaterResidues, toUtf8(upper(toU32(resname))));
    for (const tables::IndNameElement& entry : tables::kIndForceFieldNameElements) {
        const bool residue = entry.residue == "*water*"
                                 ? water
                                 : tables::equalsIgnoringCase(entry.residue, resname);
        if (residue && tables::equalsIgnoringCase(entry.name, name)) {
            return entry.element;
        }
    }
    return std::nullopt;
}

// A guessed type is the element symbol in upper case ("CL", "NA"), as MDAnalysis spells guessed
// types; a virtual site gets the empty type (no element, mass 0).
std::optional<std::string> specialElement(std::string_view name, std::string_view resname) {
    if (const auto element = forceFieldElement(name, resname)) {
        return toUtf8(upper(toU32(*element)));
    }
    return std::nullopt;
}

bool isOneLetterElement(const U32& upperOption) {
    return upperOption.size() == 1 &&
           tables::contains(tables::kIndGuessOneLetterElements, toUtf8(upperOption));
}

// Whether a prefix option of the name may be guessed as an element: a one-letter element in any
// case; a two-letter element as the name spells it in the symbol's own capitalisation ("Cl",
// "Ca"), or in upper case when its first letter is no one-letter element (MG, ZN, LI). An
// upper-case CA, CD, CL, NA, SE reads as a one-letter element and a locant (carbons alpha and
// delta, a heme nitrogen); upper-case ion names of known ion residues are the force-field table's.
bool isShortElement(const U32& upperOption, const U32& spelled) {
    if (upperOption.size() == 1) {
        return isOneLetterElement(upperOption);
    }
    if (upperOption.size() != 2) {
        return false;
    }
    if (spelled.size() == 2 && tables::contains(tables::kIndGuessTwoLetterElements, toUtf8(spelled))) {
        return true;
    }
    const std::string key = toUtf8(upperOption);
    return !isOneLetterElement(upperOption.substr(0, 1)) &&
           std::any_of(std::begin(tables::kIndGuessTwoLetterElements),
                       std::end(tables::kIndGuessTwoLetterElements),
                       [&](std::string_view e) { return tables::equalsIgnoringCase(e, key); });
}

bool isElementSymbol(std::string_view key) { return tables::elementBySymbol(key) != nullptr; }

std::string guessTypeFromName(const std::string& atomName, std::string_view resname) {
    const U32 text = toU32(atomName);
    if (text.empty()) {
        return {};
    }
    if (const auto special = specialElement(atomName, resname)) {
        return *special;
    }
    U32 bare;
    for (const char32_t c : text) {
        if (c != U'*' && c != U'+' && c != U'-') {
            bare += c;
        }
    }
    U32 letters;
    for (const char32_t c : bare) {
        if (!isDigit(c)) {
            letters += c;
        }
    }
    const U32 stem = upper(letters);
    if (stem.empty()) {
        return toUtf8(bare);
    }
    if (const auto special = specialElement(toUtf8(letters), resname)) {
        return *special;
    }
    // upper() can change the length (U+00DF -> SS); the name then has no spelling to compare
    // with a two-letter symbol and only the one-letter test applies.
    const bool sameLength = letters.size() == stem.size();
    for (std::size_t length = stem.size(); length > 0; --length) {
        const U32 candidate = stem.substr(0, length);
        const U32 written = sameLength ? letters.substr(0, length) : U32();
        const std::pair<U32, U32> options[] = {
            {candidate, written},
            {candidate.substr(0, length - 1), sameLength ? written.substr(0, length - 1) : U32()},
            {candidate.substr(1), sameLength ? written.substr(1) : U32()}};
        for (const auto& [option, spelled] : options) {
            if (isShortElement(option, spelled)) {
                return toUtf8(option);
            }
        }
        if (length <= 2) {
            return toUtf8(candidate.substr(0, 1));
        }
    }
    return {};
}

// The atomic weight of an element symbol or a guessed type, the symbol in any letter case; 0 when
// the key is no element (an unknown or empty type, a virtual site).
double massOf(const std::string& key) {
    if (const tables::IndElement* element = tables::elementBySymbol(key)) {
        return element->mass;
    }
    if (const tables::IndElement* element = tables::elementBySymbolIgnoringCase(key)) {
        return element->mass;
    }
    return 0.0;
}

// str.capitalize() as far as element symbols are concerned
std::string capitalised(const U32& text) {
    if (text.empty()) {
        return {};
    }
    U32 out = upper(text.substr(0, 1));
    for (std::size_t i = 1; i < text.size(); ++i) {
        const char32_t c = text[i];
        out += (c >= U'A' && c <= U'Z') ? static_cast<char32_t>(c - U'A' + U'a') : c;
    }
    return toUtf8(out);
}

// Fills types (guessed from names unless given) and masses (from elements when present, else
// from types).
void guessTypesAndMasses(core::Topology& top, const std::optional<std::vector<std::string>>& givenTypes,
                         const std::optional<std::vector<std::string>>& elements) {
    const std::size_t n = top.names.size();
    std::vector<double> masses(n);
    top.types.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::string type =
            givenTypes ? (*givenTypes)[i]
                       : guessTypeFromName(top.strings[top.names[i]],
                                           top.resname(top.atomResidue[i]));
        top.types[i] = top.strings.intern(type);
        masses[i] = massOf(elements ? (*elements)[i] : type);
    }
    top.masses = std::move(masses);
}

// Splits atoms into residues at every change of `key` and fills atomResidue; returns the first
// atom of each residue.
template <class Key>
std::vector<std::size_t> residueRuns(core::Topology& top, std::size_t n, Key&& key) {
    std::vector<std::size_t> firsts;
    top.atomResidue.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || key(i) != key(i - 1)) {
            firsts.push_back(i);
        }
        top.atomResidue[i] = static_cast<std::int32_t>(firsts.size() - 1);
    }
    return firsts;
}

// ---------------------------------------------------------------------------------------------
// PDB topology

struct PdbAtoms {
    std::vector<std::string> recordTypes, names, altLocs, resnames, chainIDs, icodes, segids,
        elementColumn, chargeColumn;
    std::vector<std::int64_t> serials, resids;
    std::vector<double> occupancies, tempFactors;
    bool guessedSerials = false;
};

std::optional<int> hybrid36Digit(char32_t c, bool lowerCase) {
    if (c >= U'0' && c <= U'9') {
        return static_cast<int>(c - U'0');
    }
    const char32_t first = lowerCase ? U'a' : U'A';
    if (c >= first && c <= first + 25) {
        return static_cast<int>(c - first) + 10;
    }
    return std::nullopt;
}

// Columns 7-11 as a serial: a Python int, else hybrid-36 of 5 characters; std::nullopt when
// neither applies (the caller makes one up).
std::optional<std::int64_t> decodeSerial(const U32& field, bool firstAtom) {
    if (const auto value = core::tryPythonInt(toUtf8(field))) {
        return value;
    }
    if (field.size() != 5) {
        return std::nullopt;
    }
    const char32_t lead = field[0];
    const bool upperCase = lead >= U'A' && lead <= U'Z';
    const bool lowerCase = lead >= U'a' && lead <= U'z';
    if (!upperCase && !lowerCase) {
        return std::nullopt;
    }
    std::int64_t value = 0;
    for (const char32_t c : field) {
        const auto digit = hybrid36Digit(c, lowerCase);
        if (!digit) {
            if (firstAtom) {
                throw unboundLocal("serial");
            }
            throw ReaderError("KeyError", core::pythonStringRepr(toUtf8(U32(1, c))));
        }
        value = value * 36 + *digit;
    }
    constexpr std::int64_t kPow4 = 36LL * 36 * 36 * 36;
    return upperCase ? value - 10 * kPow4 + 100000 : value + 16 * kPow4 + 100000;
}

float toFloat32(double value) {
    return static_cast<float>(value);
}

double numberOr(const U32& field, double fallback) {
    const auto value = core::tryPythonFloat(toUtf8(field));
    return value ? static_cast<double>(toFloat32(*value)) : fallback;
}

PdbAtoms readPdbAtoms(const std::vector<U32>& lines) {
    PdbAtoms atoms;
    std::int64_t madeUpSerial = 100000;
    std::int64_t previousResid = 0;
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
        const auto field = [&](std::size_t first, std::size_t last) {
            return toUtf8(strip(slice(line, first, last)));
        };
        const auto serial = decodeSerial(slice(line, 6, 11), atoms.serials.empty());
        if (serial) {
            atoms.serials.push_back(*serial);
        } else {
            atoms.serials.push_back(madeUpSerial++);
            atoms.guessedSerials = true;
        }
        atoms.recordTypes.push_back(field(0, 6));
        atoms.names.push_back(field(12, 16));
        atoms.altLocs.push_back(field(16, 17));
        atoms.resnames.push_back(field(17, 21));
        atoms.chainIDs.push_back(field(21, 22));
        if (auto resid = core::tryPythonInt(toUtf8(slice(line, 22, 26)))) {
            while (*resid - previousResid < -5000) {
                *resid += 10000;
            }
            previousResid = *resid;
            atoms.resids.push_back(*resid);
        } else {
            atoms.resids.push_back(1);
        }
        atoms.icodes.push_back(field(26, 27));
        atoms.occupancies.push_back(numberOr(slice(line, 54, 60), 0.0));
        atoms.tempFactors.push_back(numberOr(slice(line, 60, 66), 1.0));
        atoms.segids.push_back(field(72, 76));
        atoms.elementColumn.push_back(field(76, 78));
        atoms.chargeColumn.push_back(field(78, 80));
    }
    return atoms;
}

// Columns 79-80: "2-" is -2; std::nullopt for the first entry that is not a charge
std::optional<std::vector<std::int64_t>> formalCharges(const std::vector<std::string>& column) {
    if (std::all_of(column.begin(), column.end(), [](const std::string& q) { return q.empty(); })) {
        return std::nullopt;
    }
    std::vector<std::int64_t> out;
    for (const std::string& q : column) {
        if (q.empty() || q == "0") {
            out.push_back(0);
            continue;
        }
        if (q.find_first_of("+-") == std::string::npos) {
            return std::nullopt;
        }
        U32 reversed = toU32(q);
        std::reverse(reversed.begin(), reversed.end());
        const auto value = core::tryPythonInt(toUtf8(reversed));
        if (!value) {
            return std::nullopt;
        }
        out.push_back(*value);
    }
    return out;
}

// CONECT records of the whole file; std::nullopt when they are corrupt or give no bond.
std::optional<std::vector<core::Bond>> conectBonds(const std::vector<U32>& lines,
                                                   const std::vector<std::int64_t>& serials) {
    std::map<std::int64_t, std::int32_t> atomOf;
    for (std::size_t i = 0; i < serials.size(); ++i) {
        atomOf[serials[i]] = static_cast<std::int32_t>(i);
    }
    std::set<core::Bond> bonds;
    bool any = false;
    for (const U32& raw : lines) {
        if (slice(raw, 0, 6) != U"CONECT") {
            continue;
        }
        const U32 line = strip(raw);
        const std::int64_t base = core::pythonInt(toUtf8(slice(line, 6, 11)));
        const std::size_t length = line.size() > 11 ? line.size() - 11 : 0;
        const std::size_t count = length / 5;
        if (count == 0) {
            continue;
        }
        if (length % count != 0) {
            return std::nullopt;
        }
        for (std::size_t k = 0; k < count; ++k) {
            const std::int64_t partner = core::pythonInt(toUtf8(slice(line, 11 + 5 * k, 16 + 5 * k)));
            const auto a = atomOf.find(base);
            const auto b = atomOf.find(partner);
            if (a == atomOf.end() || b == atomOf.end()) {
                continue;
            }
            any = true;
            // the topology model has no self-bonds
            if (a->second != b->second) {
                bonds.insert({std::min(a->second, b->second), std::max(a->second, b->second)});
            }
        }
    }
    if (!any) {
        return std::nullopt;
    }
    return std::vector<core::Bond>(bonds.begin(), bonds.end());
}

std::shared_ptr<core::Topology> parsePdbTopology(const std::string& bytes) {
    const std::vector<U32> lines = fileLines(universalNewlines(decodeStrict(bytes)));
    PdbAtoms atoms = readPdbAtoms(lines);
    const std::size_t n = atoms.names.size();
    if (n == 0) {
        throw outOfBounds(0, 0);
    }
    auto top = std::make_shared<core::Topology>();
    core::StringPool& pool = top->strings;
    const auto intern = [&](const std::vector<std::string>& values) {
        std::vector<std::uint32_t> ids;
        ids.reserve(values.size());
        for (const std::string& v : values) {
            ids.push_back(pool.intern(v));
        }
        return ids;
    };

    // per-atom segid: columns 73-76, or the chain ID when they are blank throughout
    const bool noSegids = std::all_of(atoms.segids.begin(), atoms.segids.end(),
                                      [](const std::string& s) { return s.empty(); });
    const std::vector<std::string>& atomSegids = noSegids ? atoms.chainIDs : atoms.segids;

    top->names = intern(atoms.names);
    top->chainIDs = intern(atoms.chainIDs);
    top->altLocs = intern(atoms.altLocs);
    top->recordTypes = intern(atoms.recordTypes);
    top->ids = atoms.serials;
    top->occupancies = atoms.occupancies;
    top->tempFactors = atoms.tempFactors;
    top->formalCharges = formalCharges(atoms.chargeColumn);

    const auto firsts = residueRuns(*top, n, [&](std::size_t i) {
        return std::tie(atoms.resids[i], atoms.resnames[i], atoms.icodes[i], atomSegids[i]);
    });
    std::vector<std::uint32_t> icodes;
    std::vector<std::string> residueSegids;
    for (const std::size_t first : firsts) {
        top->resids.push_back(atoms.resids[first]);
        top->resnames.push_back(pool.intern(atoms.resnames[first]));
        icodes.push_back(pool.intern(atoms.icodes[first]));
        residueSegids.push_back(atomSegids[first]);
    }
    top->resnums = top->resids;
    top->icodes = std::move(icodes);

    const bool segmented = std::any_of(residueSegids.begin(), residueSegids.end(),
                                       [](const std::string& s) { return !s.empty(); });
    for (std::size_t r = 0; r < residueSegids.size(); ++r) {
        if (!segmented) {
            if (r == 0) {
                top->segids.push_back(pool.intern("SYSTEM"));
            }
        } else if (r == 0 || residueSegids[r] != residueSegids[r - 1]) {
            top->segids.push_back(pool.intern(residueSegids[r]));
        }
        top->residueSegment.push_back(static_cast<std::int32_t>(top->segids.size() - 1));
    }

    std::optional<std::vector<std::string>> types;
    std::optional<std::vector<std::string>> elements;
    if (std::any_of(atoms.elementColumn.begin(), atoms.elementColumn.end(),
                    [](const std::string& e) { return !e.empty(); })) {
        types = atoms.elementColumn;
        elements.emplace();
        for (const std::string& e : atoms.elementColumn) {
            const std::string symbol = capitalised(toU32(e));
            elements->push_back(isElementSymbol(symbol) ? symbol : std::string());
        }
        top->elements = intern(*elements);
    }
    guessTypesAndMasses(*top, types, elements);

    if (!atoms.guessedSerials) {
        top->bonds = conectBonds(lines, atoms.serials);
    }
    top->finish();
    return top;
}

// ---------------------------------------------------------------------------------------------
// PDB coordinates

std::string readByteRange(const std::filesystem::path& path, std::uint64_t first, std::uint64_t last) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::filesystem::filesystem_error(
            "cannot read file", path, std::make_error_code(std::errc::no_such_file_or_directory));
    }
    std::string out(static_cast<std::size_t>(last - first), '\0');
    in.seekg(static_cast<std::streamoff>(first));
    in.read(out.data(), static_cast<std::streamsize>(out.size()));
    out.resize(static_cast<std::size_t>(in.gcount()));
    return out;
}

std::vector<U32> splitLines(const U32& text) {
    std::vector<U32> lines;
    std::size_t first = 0;
    for (std::size_t i = 0; i <= text.size(); ++i) {
        if (i == text.size() || text[i] == U'\n') {
            if (i < text.size() || first < text.size()) {
                lines.push_back(text.substr(first, i - first));
            }
            first = i + 1;
        }
    }
    return lines;
}

// A CRYST1 record's cell; std::nullopt for an unreadable record, the 1 A^3 placeholder and a
// cell without lengths.
std::optional<std::array<float, 6>> crystCell(const U32& line) {
    static constexpr std::array<std::pair<std::size_t, std::size_t>, 6> kFields = {
        {{6, 15}, {15, 24}, {24, 33}, {33, 40}, {40, 47}, {47, 54}}};
    static constexpr std::array<double, 6> kPlaceholder = {1, 1, 1, 90, 90, 90};
    std::array<float, 6> cell{};
    bool placeholder = true;
    for (std::size_t i = 0; i < 6; ++i) {
        const auto value = core::tryPythonFloat(toUtf8(slice(line, kFields[i].first, kFields[i].second)));
        if (!value) {
            return std::nullopt;
        }
        cell[i] = toFloat32(*value);
        const double r = kPlaceholder[i];
        placeholder = placeholder && std::abs(static_cast<double>(cell[i]) - r) <= 1e-8 + 1e-5 * std::abs(r);
    }
    if (placeholder || (cell[0] == 0.0f && cell[1] == 0.0f && cell[2] == 0.0f)) {
        return std::nullopt;
    }
    return cell;
}

class PdbTrajectory : public FileTrajectory {
public:
    PdbTrajectory(const std::filesystem::path& path, std::optional<double> dt, std::size_t atoms)
        : m_path(path), m_pathText(pathText(path)), m_atoms(atoms) {
        const std::string bytes = readBytes(path);
        if (bytes.empty()) {
            throw emptyFileError();
        }
        m_dt = dt.value_or(1.0);
        m_fileSize = bytes.size();
        std::vector<std::uint64_t> models;
        std::vector<std::uint64_t> crysts;
        std::string header;
        std::vector<std::string> compounds;
        std::vector<std::string> remarks;
        const auto text = [](std::string_view line, std::size_t first, std::size_t last) {
            if (first >= line.size()) {
                return std::string();
            }
            return toUtf8(strip(decodeStrict(line.substr(first, last - first))));
        };
        std::size_t start = 0;
        while (start < bytes.size()) {
            std::size_t end = bytes.find('\n', start);
            end = end == std::string::npos ? bytes.size() : end + 1;
            const std::string_view line(bytes.data() + start, end - start);
            if (line.starts_with("MODEL")) {
                models.push_back(end);
            }
            if (line.starts_with("CRYST")) {
                crysts.push_back(start);
            }
            if (line.starts_with("HEADER")) {
                header = text(line, 10, 66);
            } else if (line.starts_with("TITLE")) {
                (void)text(line, 8, 80);  // decoded (and so checked) but not kept: no export writes it
            } else if (line.starts_with("COMPND")) {
                compounds.push_back(text(line, 7, 80));
            } else if (line.starts_with("REMARK")) {
                remarks.push_back(text(line, 6, std::string_view::npos));
            }
            start = end;
        }
        m_info.header = header;
        m_info.compounds = std::move(compounds);
        m_info.remarks = std::move(remarks);
        if (models.empty()) {
            m_starts = {0};
        } else {
            m_starts = models;
            if (crysts.size() == models.size()) {
                for (std::size_t k = 0; k < models.size(); ++k) {
                    m_starts[k] = std::min(models[k], crysts[k]);
                }
            }
        }
        m_frameCount = static_cast<std::int64_t>(m_starts.size());
        if (m_atoms == 0) {
            m_atoms = countAtomLines(0);
        }
        seek(0);
    }

    core::TrajectoryInfo info() const override { return m_info; }
    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        frame.positions.assign(3 * m_atoms, 0.0f);
        std::size_t count = 0;
        for (const U32& line : frameLines(index)) {
            const U32 record = slice(line, 0, 6);
            if (record == U"CRYST1") {
                m_box = crystCell(line);
                continue;
            }
            if (record != U"ATOM  " && record != U"HETATM") {
                continue;
            }
            std::array<float, 3> xyz{};
            for (std::size_t k = 0; k < 3; ++k) {
                xyz[k] = toFloat32(core::pythonFloat(toUtf8(slice(line, 30 + 8 * k, 38 + 8 * k))));
            }
            if (count >= m_atoms) {
                throw outOfBounds(m_atoms, m_atoms);
            }
            std::copy(xyz.begin(), xyz.end(), frame.positions.begin() + static_cast<std::ptrdiff_t>(3 * count));
            ++count;
        }
        if (count < m_atoms) {
            throw core::ValueError("Inconsistency in file '" + m_pathText + "': The number of atoms (" +
                                   std::to_string(count) + ") in trajectory frame " +
                                   std::to_string(index) + " differs from the number of atoms (" +
                                   std::to_string(m_atoms) +
                                   ") in the corresponding topology.\nTrajectories with varying "
                                   "numbers of atoms are currently not supported.");
        }
        frame.box = m_box;
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    std::vector<U32> frameLines(std::int64_t index) const {
        const auto k = static_cast<std::size_t>(index);
        const std::uint64_t last = k + 1 < m_starts.size() ? m_starts[k + 1] : m_fileSize;
        const std::uint64_t first = std::min<std::uint64_t>(m_starts[k], last);
        return splitLines(universalNewlines(decodeStrict(readByteRange(m_path, first, last))));
    }

    std::size_t countAtomLines(std::int64_t index) const {
        std::size_t count = 0;
        for (const U32& line : frameLines(index)) {
            const U32 record = slice(line, 0, 6);
            count += (record == U"ATOM  " || record == U"HETATM") ? 1 : 0;
        }
        return count;
    }

    std::filesystem::path m_path;
    std::string m_pathText;
    std::size_t m_atoms;
    std::uint64_t m_fileSize = 0;
    std::vector<std::uint64_t> m_starts;
    core::TrajectoryInfo m_info;
    std::optional<std::array<float, 6>> m_box;  // kept from frame to frame without CRYST1
};

// ---------------------------------------------------------------------------------------------
// GRO

std::vector<U32> groLines(const std::filesystem::path& path) {
    const std::string bytes = readBytes(path);
    if (bytes.empty()) {
        throw emptyFileError();
    }
    std::vector<U32> lines = splitPythonLines(decodeUtf8Replace(bytes));
    if (lines.size() < 2) {
        throw ReaderError("StopIteration", "");
    }
    return lines;
}

std::shared_ptr<core::Topology> parseGroTopology(const std::vector<U32>& lines,
                                                 const std::string& fileText) {
    const std::int64_t count = wrapParserErrors(fileText, kGroParser, [&] {
        return core::pythonInt(toUtf8(lines[1]));
    });
    if (count <= 0) {
        throw outOfBounds(0, 0);
    }
    const auto n = static_cast<std::size_t>(count);
    std::vector<std::int64_t> rawResids;
    std::vector<std::int64_t> ids;
    std::vector<std::string> resnames;
    std::vector<std::string> names;
    for (std::size_t i = 0; i < n; ++i) {
        if (2 + i >= lines.size()) {
            throw outOfBounds(0, 0);
        }
        const U32& line = lines[2 + i];
        const auto resid = core::tryPythonInt(toUtf8(slice(line, 0, 5)));
        const auto id = core::tryPythonInt(toUtf8(slice(line, 15, 20)));
        if (!resid || !id) {
            throw topologyLoadError(fileText, kGroParser,
                                    "Couldn't read the following line of the .gro file:\n" +
                                        toUtf8(line));
        }
        rawResids.push_back(*resid);
        ids.push_back(*id);
        resnames.push_back(toUtf8(strip(slice(line, 5, 10))));
        names.push_back(toUtf8(strip(slice(line, 10, 15))));
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (names[i].empty()) {
            throw topologyLoadError(fileText, kGroParser,
                                    "Missing atom name on line: " + std::to_string(i + 3));
        }
    }
    // a run of resid 0 after the first atom starts the next block of 100000
    std::vector<std::int64_t> resids(n);
    std::int64_t offset = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0 && rawResids[i] == 0 && rawResids[i - 1] != 0) {
            offset += 100000;
        }
        resids[i] = rawResids[i] + offset;
    }

    auto top = std::make_shared<core::Topology>();
    for (const std::string& name : names) {
        top->names.push_back(top->strings.intern(name));
    }
    top->ids = std::move(ids);
    const auto firsts = residueRuns(*top, n, [&](std::size_t i) { return std::tie(resids[i], resnames[i]); });
    for (const std::size_t first : firsts) {
        top->resids.push_back(resids[first]);
        top->resnames.push_back(top->strings.intern(resnames[first]));
        top->residueSegment.push_back(0);
    }
    top->resnums = top->resids;
    top->segids.push_back(top->strings.intern("SYSTEM"));
    guessTypesAndMasses(*top, std::nullopt, std::nullopt);
    top->finish();
    return top;
}

// nm (float32) to Å in float32 arithmetic
float nmToAngstrom(float value) {
    return value * 10.0f;
}

std::vector<float> groBoxValues(const U32& line) {
    std::vector<float> values;
    bool valid = true;
    for (const U32& token : split(line)) {
        const auto value = core::tryPythonFloat(toUtf8(token));
        if (!value) {
            valid = false;
            break;
        }
        values.push_back(toFloat32(*value));
    }
    if (!valid) {
        // every "digits . five digits" in the line
        values.clear();
        std::size_t i = 0;
        while (i < line.size()) {
            std::size_t j = i;
            while (j < line.size() && isDigit(line[j])) {
                ++j;
            }
            if (j > i && j + 6 <= line.size() && line[j] == U'.' &&
                std::all_of(line.begin() + static_cast<std::ptrdiff_t>(j + 1),
                            line.begin() + static_cast<std::ptrdiff_t>(j + 6), isDigit)) {
                values.push_back(toFloat32(core::pythonFloat(toUtf8(line.substr(i, j + 6 - i)))));
                i = j + 6;
            } else {
                i = j > i ? j : i + 1;
            }
        }
    }
    if (values.size() != 3 && values.size() != 9) {
        throw core::ValueError("GRO unitcell has neither 3 nor 9 entries.");
    }
    return values;
}

std::optional<std::array<float, 6>> groBox(const U32& line) {
    const std::vector<float> v = groBoxValues(line);
    std::array<float, 6> cell{};
    if (v.size() == 3) {
        if (std::all_of(v.begin(), v.end(), [](float x) { return std::abs(x) <= 1e-8; })) {
            return std::nullopt;
        }
        cell = {v[0], v[1], v[2], 90.0f, 90.0f, 90.0f};
    } else {
        const auto d = [&](std::size_t i) { return static_cast<double>(v[i]); };
        cell = triclinicBox({d(0), d(3), d(4)}, {d(5), d(1), d(6)}, {d(7), d(8), d(2)});
    }
    for (std::size_t i = 0; i < 3; ++i) {
        cell[i] = nmToAngstrom(cell[i]);
    }
    if (cell[0] == 0.0f && cell[1] == 0.0f && cell[2] == 0.0f) {
        return std::nullopt;
    }
    return cell;
}

class GroTrajectory : public FileTrajectory {
public:
    GroTrajectory(const std::filesystem::path& path, std::optional<double> dt) {
        const std::vector<U32> lines = groLines(path);
        const auto n = static_cast<std::size_t>(std::max<std::int64_t>(0, core::pythonInt(toUtf8(lines[1]))));
        m_dt = dt.value_or(1.0);
        m_frameCount = 1;
        // field width: from the first '.' at or after column 26 of the first atom line
        std::size_t width = 0;
        if (lines.size() > 2) {
            const std::size_t dot = lines[2].find(U'.', 25);
            width = dot == U32::npos ? 0 : dot - 24;
        }
        const std::size_t available = std::min(n, lines.size() - 2);
        m_first.positions.reserve(3 * available);
        for (std::size_t i = 0; i < available; ++i) {
            const U32& line = lines[2 + i];
            for (std::size_t k = 0; k < 3; ++k) {
                const U32 field = slice(line, 20 + k * width, 20 + (k + 1) * width);
                m_first.positions.push_back(nmToAngstrom(toFloat32(core::pythonFloat(toUtf8(field)))));
            }
        }
        if (2 + n >= lines.size()) {
            throw unboundLocal("unitcell");
        }
        m_first.box = groBox(lines[2 + n]);
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        frame = m_first;
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    core::Frame m_first;
};

}  // namespace

TopologyFile readPdbTopology(const std::filesystem::path& path, std::optional<double> dt) {
    const std::string bytes = readBytes(path);
    if (bytes.empty()) {
        throw emptyFileError();
    }
    TopologyFile out;
    out.topology = wrapParserErrors(pathText(path), kPdbParser, [&] { return parsePdbTopology(bytes); });
    out.coordinates = openPdbCoordinates(path, dt, out.topology->atomCount());
    return out;
}

TopologyFile readGroTopology(const std::filesystem::path& path, std::optional<double> dt) {
    TopologyFile out;
    out.topology = parseGroTopology(groLines(path), pathText(path));
    out.coordinates = openGroCoordinates(path, dt, out.topology->atomCount());
    return out;
}

std::unique_ptr<core::Trajectory> openPdbCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t atoms) {
    return std::make_unique<PdbTrajectory>(path, dt, atoms);
}

std::unique_ptr<core::Trajectory> openGroCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t atoms) {
    (void)atoms;  // the GRO reader takes its atom count from the file
    return std::make_unique<GroTrajectory>(path, dt);
}

}  // namespace mdbridge::io::detail
