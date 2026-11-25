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
#include "mdbridge/io/warnings.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

namespace tables = core::detail::tables;


// ---------------------------------------------------------------------------------------------
// Errors

// RD-P11: a parser error that is not one of MDBridge's own (FormatError: it names the file and
// line already) is worded "<file>: cannot read PDB topology: <message>" (MDAnalysis: "Failed to
// construct topology from file ... with parser <class 'MDAnalysis...PDBParser'>.\nError: ...").
// Errors other than ValueError / TypeError pass through.
template <class Fn>
auto wrapParserErrors(const std::string& fileText, const char* format, Fn&& parse) {
    try {
        return parse();
    } catch (const FormatError&) {
        throw;
    } catch (const core::ValueError& error) {
        const std::string_view name = error.pythonName();
        if (name != "ValueError" && name != "UnicodeDecodeError" && name != "TypeError") {
            throw;
        }
        throw FormatError(fileText + ": cannot read " + format + " topology: " + error.what());
    }
}

// ---------------------------------------------------------------------------------------------
// Types and masses the Universe guesses
//
// The procedure is MDAnalysis' (guess_atom_type / guess_masses); the data are the independent
// tables of mdbridge/core/detail/data_tables.hpp, so names those tables treat differently get
// other types and masses than MDAnalysis gives them.

// The element of an exact (upper-case) atom name of the force-field table, upper case as guessed
// types are spelled ("CL", "NA"); "" for a virtual site (no element, mass 0).
std::optional<std::string> specialElement(std::string_view key) {
    if (const tables::IndNameElement* entry = tables::forceFieldNameElement(key)) {
        return toUtf8(upper(toU32(entry->element)));
    }
    return std::nullopt;
}

// An upper-case prefix that is a guessable element ("C", "NA").
bool isShortElement(std::string_view key) {
    const auto matches = [&](std::string_view element) {
        return tables::equalsIgnoringCase(element, key);
    };
    return std::any_of(std::begin(tables::kIndGuessOneLetterElements),
                       std::end(tables::kIndGuessOneLetterElements), matches) ||
           std::any_of(std::begin(tables::kIndGuessTwoLetterElements),
                       std::end(tables::kIndGuessTwoLetterElements), matches);
}

bool isElementSymbol(std::string_view key) { return tables::elementBySymbol(key) != nullptr; }

// MDBridge's residue rule first (a deliberate improvement over MDAnalysis, which reads the name
// alone): an atom whose (residue, name), in upper case, is in kIndResidueAtomElements is that
// element ("CA" of residue CA is calcium); every other atom gets MDAnalysis' name-only guess.
std::string guessTypeFromName(const std::string& atomName, const std::string& resname) {
    const U32 text = toU32(atomName);
    if (text.empty()) {
        return {};
    }
    if (const tables::IndResidueAtomElement* ion =
            tables::residueAtomElement(toUtf8(upper(toU32(resname))), toUtf8(upper(text)))) {
        return toUtf8(upper(toU32(ion->element)));
    }
    if (const auto special = specialElement(toUtf8(upper(text)))) {
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
    if (const auto special = specialElement(toUtf8(stem))) {
        return *special;
    }
    for (std::size_t length = stem.size(); length > 0; --length) {
        const U32 candidate = stem.substr(0, length);
        for (const U32& option : {candidate, candidate.substr(0, length - 1), candidate.substr(1)}) {
            if (isShortElement(toUtf8(option))) {
                return toUtf8(option);
            }
        }
        if (length <= 2) {
            return toUtf8(candidate.substr(0, 1));
        }
    }
    return {};
}

// The atomic weight of an element symbol or a guessed type: the symbol as spelled, else in upper
// case ("NA" for Na); 0 when the key is no element (an unknown or empty type, a virtual site).
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

// ---------------------------------------------------------------------------------------------
// Messages: errors that follow the reader-message rule (file, line, field, what was expected) and
// reader warnings, one per kind with a count and the first line (docs/DIFFERENCES.md)

std::string lineText(std::size_t line) {
    return " line " + std::to_string(line);
}

// "1 atom" / "3 atoms"
std::string counted(std::size_t n, std::string_view one, std::string_view many) {
    return std::to_string(n) + " " + std::string(n == 1 ? one : many);
}

// One kind of warning: how often, and where first.
struct Tally {
    std::size_t count = 0;
    std::size_t line = 0;
    std::string detail;  // what was found first ("'XX'")
    void add(std::size_t at, std::string what = {}) {
        if (count++ == 0) {
            line = at;
            detail = std::move(what);
        }
    }
    // " (line 4)" / " (first: line 4)"; with the detail: " (first: 'XX' on line 4)"
    std::string where() const {
        const std::string place = (detail.empty() ? "line " : detail + " on line ") + std::to_string(line);
        return count == 1 ? " (" + place + ")" : " (first: " + place + ")";
    }
};

// One UTF-8 character at bytes[i]: its length and code point; 0 for a byte that does not start a
// valid UTF-8 character there.
std::size_t utf8Character(std::string_view bytes, std::size_t i, char32_t& code) {
    const auto byte = [&](std::size_t k) { return static_cast<unsigned char>(bytes[k]); };
    const unsigned char b = byte(i);
    if (b < 0x80) {
        code = b;
        return 1;
    }
    std::size_t length = 0;
    unsigned char low = 0x80;
    unsigned char high = 0xBF;
    if (b >= 0xC2 && b <= 0xDF) {
        length = 2;
        code = b & 0x1Fu;
    } else if (b >= 0xE0 && b <= 0xEF) {
        length = 3;
        code = b & 0x0Fu;
        low = b == 0xE0 ? 0xA0 : 0x80;
        high = b == 0xED ? 0x9F : 0xBF;
    } else if (b >= 0xF0 && b <= 0xF4) {
        length = 4;
        code = b & 0x07u;
        low = b == 0xF0 ? 0x90 : 0x80;
        high = b == 0xF4 ? 0x8F : 0xBF;
    } else {
        return 0;
    }
    if (i + length > bytes.size()) {
        return 0;
    }
    for (std::size_t k = 1; k < length; ++k) {
        const unsigned char c = byte(i + k);
        if (c < (k == 1 ? low : 0x80) || c > (k == 1 ? high : 0xBF)) {
            return 0;
        }
        code = (code << 6) | (c & 0x3Fu);
    }
    return length;
}

// The text of a file: UTF-8, where every byte that is not part of valid UTF-8 is read as Latin-1
// (old PDB files with Latin-1 remarks or names). `bad`: the offset of the first such byte (npos:
// none).
U32 decodeLenient(std::string_view bytes, std::size_t* bad = nullptr) {
    U32 out;
    out.reserve(bytes.size());
    if (bad) {
        *bad = std::string_view::npos;
    }
    std::size_t i = 0;
    while (i < bytes.size()) {
        char32_t code = 0;
        const std::size_t length = utf8Character(bytes, i, code);
        if (length == 0) {
            if (bad && *bad == std::string_view::npos) {
                *bad = i;
            }
            out += static_cast<char32_t>(static_cast<unsigned char>(bytes[i]));
            ++i;
        } else {
            out += code;
            i += length;
        }
    }
    return out;
}

std::size_t firstNonUtf8(std::string_view bytes) {
    std::size_t i = 0;
    while (i < bytes.size()) {
        char32_t code = 0;
        const std::size_t length = utf8Character(bytes, i, code);
        if (length == 0) {
            return i;
        }
        i += length;
    }
    return std::string_view::npos;
}

// The 1-based line (\n, \r\n and \r end lines) that holds bytes[offset]
std::size_t lineAt(std::string_view bytes, std::size_t offset) {
    std::size_t line = 1;
    for (std::size_t i = 0; i < offset && i < bytes.size(); ++i) {
        if (bytes[i] == '\n' || (bytes[i] == '\r' && (i + 1 >= bytes.size() || bytes[i + 1] != '\n'))) {
            ++line;
        }
    }
    return line;
}

void warnIfNotUtf8(const std::string& fileText, std::string_view bytes) {
    const std::size_t bad = firstNonUtf8(bytes);
    if (bad != std::string_view::npos) {
        readerWarning(fileText + lineText(lineAt(bytes, bad)) +
                      ": the text is not valid UTF-8; bytes that are not UTF-8 are read as Latin-1");
    }
}

FormatError emptyFile(const std::string& fileText) {
    return FormatError(fileText + ": file is empty");
}

// ---------------------------------------------------------------------------------------------
// Types and masses

// The type of an atom (given, or guessed from its name and residue) and its mass; counts the atoms
// whose mass is unknown (mass 0 for a type that is not empty; virtual sites and nameless atoms
// have an empty type).
void warnUnknownMasses(const core::Topology& top, const std::string& fileText,
                       const std::vector<std::size_t>& lines) {
    Tally unknown;
    for (std::size_t i = 0; i < top.names.size(); ++i) {
        if (top.masses && (*top.masses)[i] == 0.0 && !top.strings[top.types[i]].empty()) {
            unknown.add(lines[i], core::pythonStringRepr(top.strings[top.names[i]]));
        }
    }
    if (unknown.count > 0) {
        readerWarning(fileText + ": " + counted(unknown.count, "atom has", "atoms have") +
                      " no known element mass and get mass 0" + unknown.where());
    }
}

// Types guessed from the names, masses from the types (GRO, PDB without element column).
void guessTypesAndMasses(core::Topology& top) {
    const std::size_t n = top.names.size();
    std::vector<double> masses(n);
    top.types.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::string type = guessTypeFromName(top.strings[top.names[i]], top.resname(top.atomResidue[i]));
        top.types[i] = top.strings.intern(type);
        masses[i] = massOf(type);
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
    std::vector<std::size_t> lines;  // 1-based file line of each atom
    bool guessedSerials = false;
};

// The one atom-line rule of both passes: the record name in columns 1-6, not indented.
bool isAtomRecord(U32View line) {
    const U32 record = slice(line, 0, 6);
    return record == U"ATOM  " || record == U"HETATM";
}

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

enum class SerialKind { Value, Undecodable, BadHybrid36 };
struct Serial {
    SerialKind kind = SerialKind::Undecodable;
    std::int64_t value = 0;
};

// Columns 7-11 (an ATOM serial, a CONECT field) as a serial: a Python int, else hybrid-36 of 5
// characters (upper-case A0000 = 100000, lower-case a0000 = 43770016).
Serial decodeSerial(const U32& field) {
    if (const auto value = core::tryPythonInt(toUtf8(field))) {
        return {SerialKind::Value, *value};
    }
    if (field.size() != 5) {
        return {};
    }
    const char32_t lead = field[0];
    const bool upperCase = lead >= U'A' && lead <= U'Z';
    const bool lowerCase = lead >= U'a' && lead <= U'z';
    if (!upperCase && !lowerCase) {
        return {};
    }
    std::int64_t value = 0;
    for (const char32_t c : field) {
        const auto digit = hybrid36Digit(c, lowerCase);
        if (!digit) {
            return {SerialKind::BadHybrid36, 0};
        }
        value = value * 36 + *digit;
    }
    constexpr std::int64_t kPow4 = 36LL * 36 * 36 * 36;
    return {SerialKind::Value, upperCase ? value - 10 * kPow4 + 100000 : value + 16 * kPow4 + 100000};
}

float toFloat32(double value) {
    return static_cast<float>(value);
}

double numberOr(const U32& field, double fallback) {
    const auto value = core::tryPythonFloat(toUtf8(field));
    return value ? static_cast<double>(toFloat32(*value)) : fallback;
}

PdbAtoms readPdbAtoms(const std::vector<U32>& lines, const std::string& fileText) {
    PdbAtoms atoms;
    std::int64_t madeUpSerial = 100000;
    std::int64_t previousResid = 0;
    Tally ignored;
    Tally madeUp;
    Tally noResid;
    for (std::size_t k = 0; k < lines.size(); ++k) {
        const U32& raw = lines[k];
        const std::size_t lineNumber = k + 1;
        const U32 line = strip(raw);
        if (line.empty()) {
            continue;
        }
        if (line.starts_with(U"END")) {
            break;
        }
        if (!isAtomRecord(raw)) {
            if (line.starts_with(U"ATOM") || line.starts_with(U"HETATM")) {
                ignored.add(lineNumber);
            }
            continue;
        }
        const auto field = [&](std::size_t first, std::size_t last) {
            return toUtf8(strip(slice(line, first, last)));
        };
        const U32 serialField = slice(line, 6, 11);
        const Serial serial = decodeSerial(serialField);
        if (serial.kind == SerialKind::BadHybrid36) {
            throw FormatError(fileText + lineText(lineNumber) + ": invalid hybrid-36 atom serial " +
                              core::pythonStringRepr(toUtf8(serialField)));
        }
        if (serial.kind == SerialKind::Value) {
            atoms.serials.push_back(serial.value);
        } else {
            atoms.serials.push_back(madeUpSerial++);
            atoms.guessedSerials = true;
            madeUp.add(lineNumber);
        }
        atoms.lines.push_back(lineNumber);
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
            noResid.add(lineNumber);
        }
        atoms.icodes.push_back(field(26, 27));
        atoms.occupancies.push_back(numberOr(slice(line, 54, 60), 0.0));
        atoms.tempFactors.push_back(numberOr(slice(line, 60, 66), 1.0));
        atoms.segids.push_back(field(72, 76));
        atoms.elementColumn.push_back(field(76, 78));
        atoms.chargeColumn.push_back(field(78, 80));
    }
    if (ignored.count > 0) {
        readerWarning(fileText + ": " + counted(ignored.count, "line starts", "lines start") +
                      " with ATOM or HETATM but not in columns 1-6 ('ATOM  ' or 'HETATM'); ignored" +
                      ignored.where());
    }
    if (madeUp.count > 0) {
        readerWarning(fileText + ": " + (madeUp.count == 1 ? std::string("1 atom serial is not a number")
                                         : std::to_string(madeUp.count) + " atom serials are not numbers") +
                      "; serials from 100000 are used and CONECT records are not read" +
                      madeUp.where());
    }
    if (noResid.count > 0) {
        readerWarning(fileText + ": " + counted(noResid.count, "atom has", "atoms have") +
                      " no residue number in columns 23-26; resid 1 is used" + noResid.where());
    }
    return atoms;
}

// Columns 79-80: "2-" is -2; std::nullopt (and the index of the first invalid entry) for an entry
// that is not a charge
std::optional<std::vector<std::int64_t>> formalCharges(const std::vector<std::string>& column, std::size_t& invalid) {
    if (std::all_of(column.begin(), column.end(), [](const std::string& q) { return q.empty(); })) {
        invalid = column.size();
        return std::nullopt;
    }
    std::vector<std::int64_t> out;
    for (std::size_t i = 0; i < column.size(); ++i) {
        const std::string& q = column[i];
        if (q.empty() || q == "0") {
            out.push_back(0);
            continue;
        }
        invalid = i;
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
    invalid = column.size();
    return out;
}

// CONECT records of the whole file (raw columns 1-6 "CONECT"); std::nullopt when no record names
// two atoms. Serials as in ATOM records (hybrid-36 too); a blank partner field is no partner. A
// malformed record (a base or partner serial that is no number, a length MDAnalysis rejects) is
// skipped alone.
std::optional<std::vector<core::Bond>> conectBonds(const std::vector<U32>& lines,
                                                   const std::vector<std::int64_t>& serials,
                                                   const std::string& fileText) {
    std::map<std::int64_t, std::int32_t> atomOf;
    for (std::size_t i = 0; i < serials.size(); ++i) {
        atomOf[serials[i]] = static_cast<std::int32_t>(i);
    }
    std::set<core::Bond> bonds;
    bool any = false;
    Tally malformed;
    Tally single;
    Tally notAtoms;
    for (std::size_t k = 0; k < lines.size(); ++k) {
        const U32& raw = lines[k];
        if (slice(raw, 0, 6) != U"CONECT") {
            continue;
        }
        const std::size_t lineNumber = k + 1;
        const U32 line = strip(raw);
        const Serial base = decodeSerial(slice(line, 6, 11));
        const std::size_t length = line.size() > 11 ? line.size() - 11 : 0;
        const std::size_t count = length / 5;
        if (base.kind != SerialKind::Value) {
            malformed.add(lineNumber);
            continue;
        }
        if (count == 0) {
            single.add(lineNumber);
            continue;
        }
        if (length % count != 0) {
            malformed.add(lineNumber);
            continue;
        }
        std::vector<std::int64_t> partners;
        bool good = true;
        for (std::size_t f = 0; f < count; ++f) {
            const U32 field = slice(line, 11 + 5 * f, 16 + 5 * f);
            if (strip(field).empty()) {
                continue;
            }
            const Serial partner = decodeSerial(field);
            if (partner.kind != SerialKind::Value) {
                good = false;
                break;
            }
            partners.push_back(partner.value);
        }
        if (!good) {
            malformed.add(lineNumber);
            continue;
        }
        for (const std::int64_t partner : partners) {
            const auto a = atomOf.find(base.value);
            const auto b = atomOf.find(partner);
            if (a == atomOf.end() || b == atomOf.end()) {
                notAtoms.add(lineNumber);
                continue;
            }
            any = true;
            // the topology model has no self-bonds
            if (a->second != b->second) {
                bonds.insert({std::min(a->second, b->second), std::max(a->second, b->second)});
            }
        }
    }
    if (malformed.count > 0) {
        readerWarning(fileText + ": " + counted(malformed.count, "CONECT record is", "CONECT records are") +
                      " malformed and skipped" + malformed.where());
    }
    if (single.count > 0) {
        readerWarning(fileText + ": " + counted(single.count, "CONECT record names", "CONECT records name") +
                      " no bonded atom and " + (single.count == 1 ? "is" : "are") + " ignored" + single.where());
    }
    if (notAtoms.count > 0) {
        readerWarning(fileText + ": " +
                      counted(notAtoms.count, "CONECT bond names a serial", "CONECT bonds name serials") +
                      " of no atom (TER records, missing atoms) and " + (notAtoms.count == 1 ? "is" : "are") +
                      " left out" + notAtoms.where());
    }
    if (!any) {
        return std::nullopt;
    }
    return std::vector<core::Bond>(bonds.begin(), bonds.end());
}

std::shared_ptr<core::Topology> parsePdbTopology(const std::string& bytes, const std::string& fileText) {
    const std::vector<U32> lines = fileLines(universalNewlines(decodeLenient(bytes)));
    PdbAtoms atoms = readPdbAtoms(lines, fileText);
    const std::size_t n = atoms.names.size();
    if (n == 0) {
        throw FormatError(fileText + ": no ATOM or HETATM records");
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
    std::size_t invalidCharge = n;
    top->formalCharges = formalCharges(atoms.chargeColumn, invalidCharge);
    if (invalidCharge < n) {
        readerWarning(fileText + lineText(atoms.lines[invalidCharge]) + ": formal charge " +
                      core::pythonStringRepr(atoms.chargeColumn[invalidCharge]) +
                      " in columns 79-80 is not a charge (such as '2-'); no formal charges are read");
    }

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

    // Types, elements and masses: an element column with at least one entry gives types (as
    // written) and elements; atoms whose column is blank get their type, element and mass guessed
    // from the name, as in a file without the column.
    if (std::any_of(atoms.elementColumn.begin(), atoms.elementColumn.end(),
                    [](const std::string& e) { return !e.empty(); })) {
        std::vector<std::string> elements(n);
        std::vector<double> masses(n);
        top->types.resize(n);
        Tally blank;
        Tally unknown;
        for (std::size_t i = 0; i < n; ++i) {
            const std::string& column = atoms.elementColumn[i];
            std::string type = column;
            if (column.empty()) {
                type = guessTypeFromName(atoms.names[i], top->resname(top->atomResidue[i]));
                blank.add(atoms.lines[i]);
            }
            const std::string symbol = capitalised(toU32(type));
            elements[i] = isElementSymbol(symbol) ? symbol : std::string();
            if (!column.empty() && elements[i].empty()) {
                unknown.add(atoms.lines[i], core::pythonStringRepr(column));
            }
            top->types[i] = pool.intern(type);
            masses[i] = column.empty() ? massOf(type) : massOf(elements[i]);
        }
        top->elements = intern(elements);
        top->masses = std::move(masses);
        if (blank.count > 0) {
            readerWarning(fileText + ": " + counted(blank.count, "atom has", "atoms have") +
                          " no element symbol in columns 77-78; guessed from the atom name" + blank.where());
        }
        if (unknown.count > 0) {
            readerWarning(fileText + ": " + counted(unknown.count, "atom has", "atoms have") +
                          " an unknown element symbol in columns 77-78; the element is left blank" +
                          unknown.where());
        }
    } else {
        readerWarning(fileText + ": no element symbols in columns 77-78; atom types and masses are guessed "
                                 "from the atom names");
        guessTypesAndMasses(*top);
    }
    warnUnknownMasses(*top, fileText, atoms.lines);

    if (!atoms.guessedSerials) {
        top->bonds = conectBonds(lines, atoms.serials, fileText);
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

enum class CellRecord { Cell, Unreadable, Placeholder, Zero };

// A CRYST1 record's cell; std::nullopt for an unreadable record, the 1 A^3 placeholder and a
// cell without lengths (`kind` says which).
std::optional<std::array<float, 6>> crystCell(const U32& line, CellRecord& kind) {
    static constexpr std::array<std::pair<std::size_t, std::size_t>, 6> kFields = {
        {{6, 15}, {15, 24}, {24, 33}, {33, 40}, {40, 47}, {47, 54}}};
    static constexpr std::array<double, 6> kPlaceholder = {1, 1, 1, 90, 90, 90};
    std::array<float, 6> cell{};
    bool placeholder = true;
    for (std::size_t i = 0; i < 6; ++i) {
        const auto value = core::tryPythonFloat(toUtf8(slice(line, kFields[i].first, kFields[i].second)));
        if (!value) {
            kind = CellRecord::Unreadable;
            return std::nullopt;
        }
        cell[i] = toFloat32(*value);
        const double r = kPlaceholder[i];
        placeholder = placeholder && std::abs(static_cast<double>(cell[i]) - r) <= 1e-8 + 1e-5 * std::abs(r);
    }
    if (placeholder) {
        kind = CellRecord::Placeholder;
        return std::nullopt;
    }
    if (cell[0] == 0.0f && cell[1] == 0.0f && cell[2] == 0.0f) {
        kind = CellRecord::Zero;
        return std::nullopt;
    }
    kind = CellRecord::Cell;
    return cell;
}

std::string noTimesWarning(const std::string& fileText) {
    return fileText + ": the file has no time information; frames are taken to be 1 ps apart (set "
                      "inputs.dt_ps for the real time step)";
}

std::string inconsistentAtoms(const std::string& fileText, std::size_t count, std::int64_t frame, std::size_t atoms) {
    return "Inconsistency in file '" + fileText + "': The number of atoms (" + std::to_string(count) +
           ") in trajectory frame " + std::to_string(frame) + " differs from the number of atoms (" +
           std::to_string(atoms) +
           ") in the corresponding topology.\nTrajectories with varying numbers of atoms are currently not "
           "supported.";
}

class PdbTrajectory : public FileTrajectory {
public:
    PdbTrajectory(const std::filesystem::path& path, std::optional<double> dt, std::size_t atoms)
        : m_path(path), m_pathText(pathText(path)), m_atoms(atoms) {
        const std::string bytes = readBytes(path);
        if (bytes.empty()) {
            throw emptyFile(m_pathText);
        }
        warnIfNotUtf8(m_pathText, bytes);
        m_dt = dt.value_or(1.0);
        m_fileSize = bytes.size();
        std::vector<std::pair<std::uint64_t, std::size_t>> models;  // byte offset, line of the frame start
        std::vector<std::pair<std::uint64_t, std::size_t>> crysts;
        std::string header;
        std::vector<std::string> compounds;
        std::vector<std::string> remarks;
        // a record's text: the line decoded, then columns first..last (characters), stripped
        const auto text = [](std::string_view line, std::size_t first, std::size_t last) {
            return toUtf8(strip(slice(decodeLenient(line), first, last)));
        };
        std::size_t start = 0;
        std::size_t lineNumber = 1;
        while (start < bytes.size()) {
            std::size_t end = bytes.find('\n', start);
            end = end == std::string::npos ? bytes.size() : end + 1;
            const std::string_view line(bytes.data() + start, end - start);
            std::size_t breaks = 0;
            for (std::size_t i = 0; i < line.size(); ++i) {
                breaks += (line[i] == '\n' || (line[i] == '\r' && (i + 1 >= line.size() || line[i + 1] != '\n'))) ? 1 : 0;
            }
            if (line.starts_with("MODEL")) {
                models.emplace_back(end, lineNumber + breaks);
            }
            if (line.starts_with("CRYST")) {
                crysts.emplace_back(start, lineNumber);
            }
            if (line.starts_with("HEADER")) {
                header = text(line, 10, 66);
            } else if (line.starts_with("COMPND")) {
                compounds.push_back(text(line, 7, 80));
            } else if (line.starts_with("REMARK")) {
                remarks.push_back(text(line, 6, U32::npos));
            }
            lineNumber += breaks;
            start = end;
        }
        m_info.header = header;
        m_info.compounds = std::move(compounds);
        m_info.remarks = std::move(remarks);
        if (models.empty()) {
            m_frames.push_back({0, 1});
        } else {
            for (std::size_t k = 0; k < models.size(); ++k) {
                Frame f{models[k].first, models[k].second};
                if (crysts.size() == models.size() && crysts[k].first < f.start) {
                    f = {crysts[k].first, crysts[k].second};
                }
                m_frames.push_back(f);
            }
        }
        m_frameCount = static_cast<std::int64_t>(m_frames.size());
        indexFrames();
        if (m_frameCount > 1 && !dt) {
            readerWarning(noTimesWarning(m_pathText));
        }
        seek(0);
    }

    core::TrajectoryInfo info() const override { return m_info; }
    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        const Frame& f = m_frames[static_cast<std::size_t>(index)];
        frame.positions.assign(3 * m_atoms, 0.0f);
        std::size_t count = 0;
        const std::vector<U32> lines = frameLines(index);
        for (std::size_t l = 0; l < lines.size() && count < m_atoms; ++l) {
            const U32& line = lines[l];
            if (!isAtomRecord(line)) {
                continue;
            }
            for (std::size_t k = 0; k < 3; ++k) {
                const U32 field = slice(line, 30 + 8 * k, 38 + 8 * k);
                const auto value = core::tryPythonFloat(toUtf8(field));
                if (!value) {
                    throw FormatError(m_pathText + lineText(f.line + l) +
                                      (m_frameCount > 1 ? " (model " + std::to_string(index + 1) + ")" : "") +
                                      ": " + std::string(1, static_cast<char>('x' + k)) + " coordinate " +
                                      core::pythonStringRepr(toUtf8(field)) + " is not a number");
                }
                frame.positions[3 * count + k] = toFloat32(*value);
            }
            ++count;
        }
        frame.box = f.cell;
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    struct Frame {
        std::uint64_t start = 0;
        std::size_t line = 1;  // the file line where the frame starts
        std::optional<std::array<float, 6>> cell;
    };

    std::vector<U32> frameLines(std::int64_t index) const {
        const auto k = static_cast<std::size_t>(index);
        const std::uint64_t last = k + 1 < m_frames.size() ? m_frames[k + 1].start : m_fileSize;
        const std::uint64_t first = std::min<std::uint64_t>(m_frames[k].start, last);
        return splitLines(universalNewlines(decodeLenient(readByteRange(m_path, first, last))));
    }

    // Every frame's atom count (checked) and cell: its last CRYST1 record, else the cell of the
    // nearest earlier frame in the file that has one (CRYST1 records before the first frame count
    // as earlier), whatever order frames are read in.
    void indexFrames() {
        std::optional<std::array<float, 6>> cell;
        Tally unreadable;
        Tally placeholder;
        if (m_frames.front().start > 0) {  // the text before the first frame
            const std::vector<U32> lines = splitLines(
                universalNewlines(decodeLenient(readByteRange(m_path, 0, m_frames.front().start))));
            for (std::size_t l = 0; l < lines.size(); ++l) {
                if (slice(lines[l], 0, 6) == U"CRYST1") {
                    CellRecord kind = CellRecord::Cell;
                    cell = crystCell(lines[l], kind);
                    if (kind == CellRecord::Unreadable) {
                        unreadable.add(1 + l);
                    } else if (kind == CellRecord::Placeholder) {
                        placeholder.add(1 + l);
                    }
                }
            }
        }
        for (std::size_t k = 0; k < m_frames.size(); ++k) {
            Frame& f = m_frames[k];
            const std::vector<U32> lines = frameLines(static_cast<std::int64_t>(k));
            std::size_t count = 0;
            std::size_t endLine = 0;    // the first END / ENDMDL line of the frame
            std::size_t afterEnd = 0;   // the first atom record after it
            for (std::size_t l = 0; l < lines.size(); ++l) {
                const U32& line = lines[l];
                if (slice(line, 0, 6) == U"CRYST1") {
                    CellRecord kind = CellRecord::Cell;
                    cell = crystCell(line, kind);
                    if (kind == CellRecord::Unreadable) {
                        unreadable.add(f.line + l);
                    } else if (kind == CellRecord::Placeholder) {
                        placeholder.add(f.line + l);
                    }
                } else if (isAtomRecord(line)) {
                    ++count;
                    if (endLine != 0 && afterEnd == 0) {
                        afterEnd = f.line + l;
                    }
                } else if (endLine == 0 && strip(line).starts_with(U"END")) {
                    endLine = f.line + l;
                }
            }
            f.cell = cell;
            if (k == 0 && m_atoms == 0) {
                m_atoms = count;
            }
            if (count == m_atoms) {
                continue;
            }
            if (afterEnd != 0) {
                throw FormatError(m_pathText + lineText(afterEnd) + ": ATOM or HETATM record after END (line " +
                                  std::to_string(endLine) + ")");
            }
            if (k == 0) {
                throw core::ValueError(inconsistentAtoms(m_pathText, count, 0, m_atoms));
            }
            throw FormatError(m_pathText + ": model " + std::to_string(k + 1) + " has " +
                              counted(count, "atom", "atoms") + ", model 1 has " + std::to_string(m_atoms) +
                              " (line " + std::to_string(f.line) + ")");
        }
        if (unreadable.count > 0) {
            readerWarning(m_pathText + ": " + counted(unreadable.count, "CRYST1 record", "CRYST1 records") +
                          " could not be read (six numbers in columns 7-54 expected); no unit cell is used" +
                          unreadable.where());
        }
        if (placeholder.count > 0) {
            readerWarning(m_pathText + ": " + counted(placeholder.count, "CRYST1 record is", "CRYST1 records are") +
                          " the 1 A^3 placeholder cell; no unit cell is used" + placeholder.where());
        }
    }

    std::filesystem::path m_path;
    std::string m_pathText;
    std::size_t m_atoms;
    std::uint64_t m_fileSize = 0;
    std::vector<Frame> m_frames;
    core::TrajectoryInfo m_info;
};

// ---------------------------------------------------------------------------------------------
// GRO

// The lines of a text in bytes: offsets of each line's start and end (without its terminator);
// \n, \r\n and \r end lines.
struct ByteLines {
    explicit ByteLines(std::string_view bytes) : m_bytes(bytes) {}
    // the next line from `pos`: [first, last) of its text; false at the end of the text
    bool next(std::size_t& pos, std::size_t& first, std::size_t& last) const {
        if (pos >= m_bytes.size()) {
            return false;
        }
        first = pos;
        std::size_t i = pos;
        while (i < m_bytes.size() && m_bytes[i] != '\n' && m_bytes[i] != '\r') {
            ++i;
        }
        last = i;
        if (i < m_bytes.size()) {
            i += (m_bytes[i] == '\r' && i + 1 < m_bytes.size() && m_bytes[i + 1] == '\n') ? 2 : 1;
        }
        pos = i;
        return true;
    }

private:
    std::string_view m_bytes;
};

// nm (float32) to Å in float32 arithmetic
float nmToAngstrom(float value) {
    return value * 10.0f;
}

std::optional<std::vector<float>> numbersOf(const std::vector<U32>& tokens) {
    std::vector<float> values;
    for (const U32& token : tokens) {
        const auto value = core::tryPythonFloat(toUtf8(token));
        if (!value) {
            return std::nullopt;
        }
        values.push_back(toFloat32(*value));
    }
    return values;
}

// The values of a box line: split on whitespace; else fixed fields of 10 characters (GROMACS'
// %10.5f, whose values run together beyond 9999 nm); else numbers with five decimals written
// without separators (signs kept). std::nullopt when none of these gives 3 or 9 values.
std::optional<std::vector<float>> groBoxValues(const U32& raw) {
    const U32 line = rstrip(raw);
    const auto usable = [](const std::optional<std::vector<float>>& v) {
        return v && (v->size() == 3 || v->size() == 9);
    };
    if (auto values = numbersOf(split(line)); usable(values)) {
        return values;
    }
    if (!line.empty() && line.size() % 10 == 0) {
        std::vector<U32> fields;
        for (std::size_t i = 0; i < line.size(); i += 10) {
            fields.push_back(line.substr(i, 10));
        }
        if (auto values = numbersOf(fields); usable(values)) {
            return values;
        }
    }
    std::vector<U32> tokens;
    std::size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && isSpace(line[i])) {
            ++i;
        }
        if (i == line.size()) {
            break;
        }
        std::size_t j = i + (line[i] == U'-' ? 1 : 0);
        const std::size_t digits = j;
        while (j < line.size() && isDigit(line[j])) {
            ++j;
        }
        if (j == digits || j + 6 > line.size() || line[j] != U'.' ||
            !std::all_of(line.begin() + static_cast<std::ptrdiff_t>(j + 1),
                         line.begin() + static_cast<std::ptrdiff_t>(j + 6), isDigit)) {
            return std::nullopt;
        }
        tokens.push_back(line.substr(i, j + 6 - i));
        i = j + 6;
    }
    if (auto values = numbersOf(tokens); usable(values)) {
        return values;
    }
    return std::nullopt;
}

enum class BoxKind { Cell, Zero, Invalid };

std::optional<std::array<float, 6>> groBox(const std::vector<float>& v, BoxKind& kind) {
    std::array<float, 6> cell{};
    kind = BoxKind::Cell;
    if (v.size() == 3) {
        if (std::all_of(v.begin(), v.end(), [](float x) { return std::abs(x) <= 1e-8; })) {
            kind = BoxKind::Zero;
            return std::nullopt;
        }
        cell = {v[0], v[1], v[2], 90.0f, 90.0f, 90.0f};
    } else {
        const auto d = [&](std::size_t i) { return static_cast<double>(v[i]); };
        cell = triclinicBox({d(0), d(3), d(4)}, {d(5), d(1), d(6)}, {d(7), d(8), d(2)});
        if (cell[0] == 0.0f && cell[1] == 0.0f && cell[2] == 0.0f) {
            kind = std::all_of(v.begin(), v.end(), [](float x) { return x == 0.0f; }) ? BoxKind::Zero
                                                                                     : BoxKind::Invalid;
            return std::nullopt;
        }
    }
    for (std::size_t i = 0; i < 3; ++i) {
        cell[i] = nmToAngstrom(cell[i]);
    }
    if (cell[0] == 0.0f && cell[1] == 0.0f && cell[2] == 0.0f) {
        return std::nullopt;
    }
    return cell;
}

// A GRO file's frames (title, atom count, atoms, box line each; GROMACS writes trajectories as
// such frames one after the other): where each frame's atom lines are, and its box.
struct GroIndex {
    std::size_t atoms = 0;
    struct Frame {
        std::uint64_t first = 0;  // byte offset of the first atom line
        std::uint64_t last = 0;   // byte offset just after the box line
        std::size_t line = 3;     // file line of the first atom line
        std::optional<std::array<float, 6>> box;
    };
    std::vector<Frame> frames;
    std::vector<U32> firstLines;  // the first frame: title, count, atom lines, box line
};

GroIndex indexGro(const std::string& bytes, const std::string& fileText) {
    if (bytes.empty()) {
        throw emptyFile(fileText);
    }
    warnIfNotUtf8(fileText, bytes);
    const ByteLines reader(bytes);
    const auto decoded = [&](std::size_t first, std::size_t last) {
        return decodeLenient(std::string_view(bytes).substr(first, last - first));
    };
    GroIndex index;
    std::size_t pos = 0;
    std::size_t first = 0;
    std::size_t last = 0;
    std::size_t lineNumber = 0;  // of the line read last
    Tally zero;
    Tally invalid;
    for (std::size_t frame = 0;; ++frame) {
        const std::size_t frameStart = pos;
        const std::size_t titleLine = lineNumber + 1;
        // a later frame: only when the rest of the file is more than white space
        if (frame > 0) {
            bool blank = true;
            for (std::size_t i = pos; i < bytes.size() && blank; ++i) {
                blank = isSpace(static_cast<unsigned char>(bytes[i]));
            }
            if (blank) {
                break;
            }
        }
        const auto incomplete = [&] {
            readerWarning(fileText + lineText(titleLine) + ": the text after frame " + std::to_string(frame) +
                          " is not a complete GRO frame and is ignored");
        };
        if (!reader.next(pos, first, last)) {  // title
            throw ReaderError("StopIteration", "");
        }
        ++lineNumber;
        U32 title = decoded(first, last);
        if (!reader.next(pos, first, last)) {
            if (frame > 0) {
                incomplete();
                break;
            }
            throw FormatError(fileText + ": the file ends after its title line; line 2 should hold the atom count");
        }
        ++lineNumber;
        const U32 countText = decoded(first, last);
        const auto count = core::tryPythonInt(toUtf8(countText));
        if (frame > 0 && (!count || *count != static_cast<std::int64_t>(index.atoms))) {
            if (count) {
                throw FormatError(fileText + lineText(lineNumber) + ": frame " + std::to_string(frame + 1) + " has " +
                                  counted(static_cast<std::size_t>(std::max<std::int64_t>(*count, 0)), "atom",
                                          "atoms") +
                                  ", frame 1 has " + std::to_string(index.atoms));
            }
            incomplete();
            break;
        }
        if (!count) {
            throw FormatError(fileText + lineText(lineNumber) + ": atom count " +
                              core::pythonStringRepr(toUtf8(countText)) + " is not an integer");
        }
        if (*count <= 0) {
            throw FormatError(fileText + ": no atoms (the atom count on line 2 is " + std::to_string(*count) + ")");
        }
        const auto n = static_cast<std::size_t>(*count);
        index.atoms = n;
        if (frame == 0) {
            index.firstLines.push_back(std::move(title));
            index.firstLines.push_back(countText + U"\n");
        }
        GroIndex::Frame f;
        f.first = pos;
        f.line = lineNumber + 1;
        std::size_t read = 0;
        while (read <= n && reader.next(pos, first, last)) {
            ++lineNumber;
            ++read;
            if (frame == 0) {
                index.firstLines.push_back(decoded(first, last) + U"\n");
            }
        }
        if (read <= n) {  // the file ended early
            if (frame > 0) {
                pos = frameStart;
                incomplete();
                break;
            }
            if (read == n) {
                throw FormatError(fileText + ": missing box line after atom " + std::to_string(n) + " (line " +
                                  std::to_string(lineNumber) + ")");
            }
            const bool lastIsBox = read > 0 && groBoxValues(index.firstLines.back()).has_value();
            throw FormatError(fileText + ": the header (line 2) says " + counted(n, "atom", "atoms") +
                              " but the file has " + counted(read - (lastIsBox ? 1 : 0), "atom line", "atom lines"));
        }
        const U32 boxLine = decoded(first, last);
        f.last = pos;
        const auto values = groBoxValues(boxLine);
        if (!values) {
            throw FormatError(fileText + lineText(lineNumber) + ": GRO box line " +
                              core::pythonStringRepr(toUtf8(boxLine)) + " has neither 3 nor 9 values");
        }
        BoxKind kind = BoxKind::Cell;
        f.box = groBox(*values, kind);
        if (kind == BoxKind::Zero) {
            zero.add(lineNumber);
        } else if (kind == BoxKind::Invalid) {
            invalid.add(lineNumber);
        }
        index.frames.push_back(f);
    }
    if (zero.count > 0) {
        readerWarning(fileText + ": " + counted(zero.count, "box is", "boxes are") + " zero; no periodic box" +
                      zero.where());
    }
    if (invalid.count > 0) {
        readerWarning(fileText + ": " + counted(invalid.count, "box has", "boxes have") +
                      " vectors that do not form a cell; no periodic box" + invalid.where());
    }
    return index;
}

std::shared_ptr<core::Topology> parseGroTopology(const std::vector<U32>& lines, const std::string& fileText) {
    const auto n = static_cast<std::size_t>(core::pythonInt(toUtf8(lines[1])));
    std::vector<std::int64_t> rawResids;
    std::vector<std::int64_t> ids;
    std::vector<std::string> resnames;
    std::vector<std::string> names;
    for (std::size_t i = 0; i < n; ++i) {
        const U32& line = lines[2 + i];
        const auto resid = core::tryPythonInt(toUtf8(slice(line, 0, 5)));
        const auto id = core::tryPythonInt(toUtf8(slice(line, 15, 20)));
        if (!resid || !id) {
            // RD-P11, RM (MDAnalysis: "Failed to load from the topology file ... with parser
            // <class '...GROParser'>.\nError: Couldn't read the following line of the .gro file: ...")
            const bool residue = !resid;
            throw FormatError(fileText + " line " + std::to_string(i + 3) + ": " +
                              (residue ? "residue number " : "atom number ") +
                              core::pythonStringRepr(toUtf8(residue ? slice(line, 0, 5) : slice(line, 15, 20))) +
                              (residue ? " (columns 1-5)" : " (columns 16-20)") + " is not an integer");
        }
        rawResids.push_back(*resid);
        ids.push_back(*id);
        resnames.push_back(toUtf8(strip(slice(line, 5, 10))));
        names.push_back(toUtf8(strip(slice(line, 10, 15))));
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (names[i].empty()) {
            throw FormatError(fileText + " line " + std::to_string(i + 3) +
                              ": no atom name in columns 11-15");  // RD-P11, RM
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
    guessTypesAndMasses(*top);
    std::vector<std::size_t> lineNumbers(n);
    for (std::size_t i = 0; i < n; ++i) {
        lineNumbers[i] = i + 3;
    }
    warnUnknownMasses(*top, fileText, lineNumbers);
    top->finish();
    return top;
}

class GroTrajectory : public FileTrajectory {
public:
    GroTrajectory(const std::filesystem::path& path, std::optional<double> dt, GroIndex index)
        : m_path(path), m_pathText(pathText(path)), m_index(std::move(index)) {
        m_dt = dt.value_or(1.0);
        m_frameCount = static_cast<std::int64_t>(m_index.frames.size());
        if (m_frameCount > 1 && !dt) {
            readerWarning(noTimesWarning(m_pathText));
        }
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return false; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        const GroIndex::Frame& f = m_index.frames[static_cast<std::size_t>(index)];
        const std::vector<U32> lines =
            splitLines(universalNewlines(decodeLenient(readByteRange(m_path, f.first, f.last))));
        const std::size_t n = m_index.atoms;
        const std::string frameText = m_frameCount > 1 ? " (frame " + std::to_string(index + 1) + ")" : "";
        // the field width: the distance between the decimal points of y and z (GROMACS writes
        // fixed %w.df fields), the x field's point at the same place in its field
        std::array<std::size_t, 3> dots{};
        std::size_t from = 20;
        bool found = true;
        for (std::size_t k = 0; k < 3 && found; ++k) {
            dots[k] = lines[0].find(U'.', from);
            found = dots[k] != U32::npos;
            from = dots[k] + 1;
        }
        const std::size_t width = found ? dots[2] - dots[1] : 0;
        if (!found || width == 0 || dots[1] - dots[0] != width || dots[0] >= 20 + width) {
            throw FormatError(m_pathText + lineText(f.line) + frameText +
                              ": cannot determine the coordinate field width (decimal points of x, y and z)");
        }
        frame.positions.resize(3 * n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < 3; ++k) {
                const U32 field = slice(lines[i], 20 + k * width, 20 + (k + 1) * width);
                const auto value = core::tryPythonFloat(toUtf8(field));
                if (!value) {
                    throw FormatError(m_pathText + lineText(f.line + i) + frameText + ": " +
                                      std::string(1, static_cast<char>('x' + k)) + " coordinate " +
                                      core::pythonStringRepr(toUtf8(field)) + " is not a number");
                }
                frame.positions[3 * i + k] = nmToAngstrom(toFloat32(*value));
            }
        }
        frame.box = f.box;
        frame.time = static_cast<double>(index) * m_dt;
    }

private:
    std::filesystem::path m_path;
    std::string m_pathText;
    GroIndex m_index;
};

}  // namespace

TopologyFile readPdbTopology(const std::filesystem::path& path, std::optional<double> dt) {
    const std::string bytes = readBytes(path);
    const std::string fileText = pathText(path);
    if (bytes.empty()) {
        throw emptyFile(fileText);
    }
    warnIfNotUtf8(fileText, bytes);
    TopologyFile out;
    out.topology = wrapParserErrors(fileText, "PDB", [&] { return parsePdbTopology(bytes, fileText); });
    out.coordinates = openPdbCoordinates(path, dt, out.topology->atomCount());
    return out;
}

TopologyFile readGroTopology(const std::filesystem::path& path, std::optional<double> dt) {
    const std::string fileText = pathText(path);
    GroIndex index = indexGro(readBytes(path), fileText);
    TopologyFile out;
    out.topology = parseGroTopology(index.firstLines, fileText);
    out.coordinates = std::make_unique<GroTrajectory>(path, dt, std::move(index));
    return out;
}

std::unique_ptr<core::Trajectory> openPdbCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t atoms) {
    return std::make_unique<PdbTrajectory>(path, dt, atoms);
}

std::unique_ptr<core::Trajectory> openGroCoordinates(const std::filesystem::path& path,
                                                     std::optional<double> dt, std::size_t atoms) {
    (void)atoms;  // the GRO reader takes its atom count from the file
    return std::make_unique<GroTrajectory>(path, dt, indexGro(readBytes(path), pathText(path)));
}

}  // namespace mdbridge::io::detail
