// MDATK v0.0.4.1's PDB tools (mdatk/pdbtools.py, MDATK's own code) on io's PDB model, with the
// fixes of docs/DIFFERENCES.md "Chemistry" (CH-01, CH-03 ... CH-07, CH-13, CH-15). Every loop
// follows MDATK's order, so notes and outputs are MDATK's byte for byte (tests/chem/
// test_pdbtools.cpp).

#include "mdbridge/chem/pdbtools.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::chem {

using io::PdbAtom;
using io::PdbError;
using io::PdbResidue;
using io::PdbStructure;

namespace {

// --- tables (pdbtools.py:65-79; MDATK's own short lists) ------------------------------------------

constexpr std::array<std::pair<std::string_view, char>, 20> kAmino = {{
    {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'},
    {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'}, {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'},
    {"PRO", 'P'}, {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'},
}};
// Amber protonation / disulfide variants: one-letter codes only, not standard PDB names
constexpr std::array<std::pair<std::string_view, char>, 8> kAminoVariants = {{
    {"HID", 'H'}, {"HIE", 'H'}, {"HIP", 'H'}, {"CYX", 'C'}, {"CYM", 'C'}, {"ASH", 'D'}, {"GLH", 'E'},
    {"LYN", 'K'},
}};
constexpr std::array<std::string_view, 13> kStandardNucleic = {"A",  "C",  "G",  "U",  "I",  "N",  "DA",
                                                               "DC", "DG", "DT", "DI", "DN", "UNK"};
constexpr std::array<std::string_view, 14> kWater = {"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3",
                                                     "TIP4", "TIP5", "T3P", "T4P", "T4E", "SPC", "OPC"};

// L_CHIRAL_SIGNS: (centre, (atom 1, atom 2, atom 3), sign of an L-ribose); "N" is N9 / N1
struct ChiralCentre {
    std::string_view centre;
    std::array<std::string_view, 3> partners;
    double lSign;
};
constexpr std::array<ChiralCentre, 3> kLChiralSigns = {{
    {"C4'", {"O4'", "C3'", "C5'"}, 1.0},
    {"C3'", {"O3'", "C4'", "C2'"}, -1.0},
    {"C1'", {"O4'", "N", "C2'"}, 1.0},
}};

template <std::size_t N>
std::optional<char> lookup(const std::array<std::pair<std::string_view, char>, N>& table, std::string_view key) {
    for (const auto& [name, code] : table) {
        if (name == key) {
            return code;
        }
    }
    return std::nullopt;
}

template <std::size_t N>
bool contains(const std::array<std::string_view, N>& table, std::string_view key) {
    return std::find(table.begin(), table.end(), key) != table.end();
}

std::string withPrimes(std::string name) {
    std::replace(name.begin(), name.end(), '*', '\'');
    return name;
}

bool isPhosphateAtom(const PdbAtom& atom) {
    return contains(kPhosphate5p, withPrimes(atom.name));
}

bool isAmino(const PdbResidue& residue) {
    const std::string name = io::pythonUpper(residue.name);
    return lookup(kAmino, name).has_value() || lookup(kAminoVariants, name).has_value();
}

std::string chainOrBlank(const std::string& chain) {
    return chain.empty() ? std::string("_") : chain;
}

std::string join(const std::vector<std::string>& items, std::string_view sep) {
    std::string out;
    for (std::size_t k = 0; k < items.size(); ++k) {
        if (k) {
            out += sep;
        }
        out += items[k];
    }
    return out;
}

// Python's f"{x:+.2f}" / f"{x:8.3f}" (nan / inf as Python spells them)
std::string formatFixed(double x, int width, int precision, bool plus) {
    std::string body;
    if (std::isnan(x)) {
        body = plus ? "+nan" : "nan";
    } else if (std::isinf(x)) {
        body = x > 0 ? (plus ? "+inf" : "inf") : "-inf";
    } else {
        char buf[64];
        std::snprintf(buf, sizeof buf, plus ? "%+.*f" : "%.*f", precision, x);
        body = buf;
    }
    if (body.size() < static_cast<std::size_t>(width)) {
        body.insert(0, static_cast<std::size_t>(width) - body.size(), ' ');
    }
    return body;
}

// the code points of `text` from `first` on (Python's text[first:])
std::string codePointsFrom(std::string_view text, std::size_t first) {
    const std::u32string u = core::decodeUtf8(text);
    return first >= u.size() ? std::string() : core::encodeUtf8(std::u32string_view(u).substr(first));
}

const Vec3* find(const AtomCoords& coords, std::string_view name) {
    for (const auto& [key, xyz] : coords) {
        if (key == name) {
            return &xyz;
        }
    }
    return nullptr;
}

template <class Fn>
void forEachResidue(PdbStructure& structure, Fn&& fn) {
    for (auto& segment : structure.segments) {
        for (PdbResidue& residue : segment) {
            fn(residue);
        }
    }
}

std::vector<PdbResidue*> allResidues(PdbStructure& structure) {
    std::vector<PdbResidue*> out;
    forEachResidue(structure, [&](PdbResidue& r) { out.push_back(&r); });
    return out;
}

bool allSolvent(const std::vector<PdbResidue>& segment) {
    return std::all_of(segment.begin(), segment.end(), [](const PdbResidue& r) { return isSolvent(r); });
}

// _polymer_segments (indices)
std::vector<std::size_t> polymerSegments(const PdbStructure& structure) {
    std::vector<std::size_t> out;
    for (std::size_t k = 0; k < structure.segments.size(); ++k) {
        if (!allSolvent(structure.segments[k])) {
            out.push_back(k);
        }
    }
    return out;
}

// _resegment: segments split wherever the chain id changes
void resegment(PdbStructure& structure) {
    std::vector<std::vector<PdbResidue>> segments;
    for (auto& segment : structure.segments) {
        std::vector<PdbResidue> current;
        for (PdbResidue& residue : segment) {
            if (!current.empty() && current.back().chain != residue.chain) {
                segments.push_back(std::move(current));
                current.clear();
            }
            current.push_back(std::move(residue));
        }
        if (!current.empty()) {
            segments.push_back(std::move(current));
        }
    }
    structure.segments = std::move(segments);
}

// a Counter keyed by T in order of first appearance
template <class T>
struct Counter {
    std::vector<std::pair<T, std::int64_t>> items;
    void add(const T& key, std::int64_t n = 1) {
        for (auto& [k, v] : items) {
            if (k == key) {
                v += n;
                return;
            }
        }
        items.emplace_back(key, n);
    }
    std::int64_t total() const {
        std::int64_t n = 0;
        for (const auto& item : items) {
            n += item.second;
        }
        return n;
    }
};

std::string volumesText(const std::vector<std::pair<std::string, double>>& volumes) {
    std::vector<std::string> parts;
    for (const auto& [centre, v] : volumes) {
        parts.push_back(centre + " " + formatFixed(v, 0, 2, true));
    }
    return join(parts, ", ");
}

// _residue_ranges: "C:45-57, C:64-73" (consecutive numbers merged per chain)
std::string residueRanges(const std::vector<const PdbResidue*>& residues) {
    std::vector<std::tuple<std::string, std::int64_t, std::int64_t>> parts;
    for (const PdbResidue* r : residues) {
        if (!parts.empty() && std::get<0>(parts.back()) == r->chain && r->resseq == std::get<2>(parts.back()) + 1) {
            std::get<2>(parts.back()) = r->resseq;
        } else {
            parts.emplace_back(r->chain, r->resseq, r->resseq);
        }
    }
    std::vector<std::string> out;
    for (const auto& [chain, a, b] : parts) {
        out.push_back(chainOrBlank(chain) + ":" + std::to_string(a) + (b != a ? "-" + std::to_string(b) : ""));
    }
    return join(out, ", ");
}

// _compress: "3-4, 7"
std::string compress(const std::vector<std::int64_t>& numbers) {
    std::vector<std::pair<std::int64_t, std::int64_t>> parts;
    for (const std::int64_t n : numbers) {
        if (!parts.empty() && n == parts.back().second + 1) {
            parts.back().second = n;
        } else {
            parts.emplace_back(n, n);
        }
    }
    std::vector<std::string> out;
    for (const auto& [a, b] : parts) {
        out.push_back(a == b ? std::to_string(a) : std::to_string(a) + "-" + std::to_string(b));
    }
    return join(out, ", ");
}

std::string schemeList() {
    return "amber, ccd, rosetta, pdb";
}

void checkScheme(std::string_view scheme) {
    if (!contains(kSchemes, scheme)) {
        throw PdbError("unknown residue name scheme " + core::pythonStringRepr(std::string(scheme)) +
                       " (choose from " + schemeList() + ")");
    }
}

// CH-06 (MDATK 0.1.0 _bonded_on): O3' of `before` within kO3PBonded of P of `after`
bool bondedOn(const PdbResidue& before, const PdbResidue& after) {
    const AtomCoords a = before.coords();
    const AtomCoords b = after.coords();
    const Vec3* o3 = find(a, "O3'");
    const Vec3* p = find(b, "P");
    if (o3 == nullptr || p == nullptr) {
        return false;
    }
    const double dx = (*o3)[0] - (*p)[0];
    const double dy = (*o3)[1] - (*p)[1];
    const double dz = (*o3)[2] - (*p)[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz) <= kO3PBonded;
}

bool phosphateOnly(const PdbResidue& residue) {
    return !residue.atoms.empty() && std::all_of(residue.atoms.begin(), residue.atoms.end(), isPhosphateAtom);
}

std::string endsState(bool five, bool three) {
    return five && three ? "N" : five ? "5" : three ? "3" : "";
}

// the CCD component of a nucleotide name (no terminal variant)
std::string ccdComponent(const io::ResidueName& info) {
    if (info.sugar == "dna") {
        return std::string(info.chirality == 'L' ? "0D" : "D") + info.base;
    }
    return schemeName(info.base, info.chirality, "", "ccd");
}

// the free single-character chain ids, lazily (pdbtools._free_ids)
class FreeIds {
public:
    explicit FreeIds(const PdbStructure& structure) : m_structure(structure) {}
    std::string next() {
        if (!m_started) {
            m_started = true;
            for (const std::string& c : m_structure.chainIds()) {
                m_used.insert(c);
            }
        }
        while (m_pos < io::kChainIds.size()) {
            const std::string c(1, io::kChainIds[m_pos++]);
            if (!m_used.count(c)) {
                m_used.insert(c);
                return c;
            }
        }
        throw PdbError("no unused single-character chain id left");
    }

private:
    const PdbStructure& m_structure;
    bool m_started = false;
    std::set<std::string> m_used;
    std::size_t m_pos = 0;
};

// one item of parse_ranges' text form: ^([^:\s]?)(?::(-?\d+)(?:-(-?\d+))?)?$
std::optional<ResidueRange> matchRangeItem(const std::string& item) {
    std::u32string u = core::decodeUtf8(item);
    if (!u.empty() && u.back() == U'\n') {  // `$` matches before a final newline
        u.pop_back();
    }
    ResidueRange out;
    std::size_t i = 0;
    if (i < u.size() && u[i] != U':' && !core::isPythonSpace(u[i])) {
        out.chain = core::encodeUtf8(u.substr(0, 1));
        ++i;
    }
    if (i == u.size()) {
        return out;
    }
    if (u[i] != U':') {
        return std::nullopt;
    }
    ++i;
    const auto number = [&](std::size_t& at) -> std::optional<std::string> {
        const std::size_t first = at;
        if (at < u.size() && u[at] == U'-') {
            ++at;
        }
        const std::size_t digits = at;
        while (at < u.size() && core::isPythonDigit(u[at])) {
            ++at;
        }
        if (at == digits) {
            return std::nullopt;
        }
        return core::encodeUtf8(u.substr(first, at - first));
    };
    const std::optional<std::string> a = number(i);
    if (!a) {
        return std::nullopt;
    }
    std::optional<std::string> b;
    if (i < u.size() && u[i] == U'-') {
        ++i;
        b = number(i);
        if (!b) {
            return std::nullopt;
        }
    }
    if (i != u.size()) {
        return std::nullopt;
    }
    try {
        out.start = core::pythonInt(*a);
        out.end = b ? core::pythonInt(*b) : *out.start;
    } catch (const std::overflow_error&) {
        return std::nullopt;  // beyond 64 bits: no residue number
    }
    return out;
}

ResidueRange readRangeText(const std::string& item) {
    std::optional<ResidueRange> r = matchRangeItem(item);
    if (!r) {
        throw PdbError("cannot read residue range " + core::pythonStringRepr(item) +
                       " (expected e.g. C:45-73, C:45 or C)");
    }
    return *r;
}

// re.split(r"[,;\s]+", text.strip())
std::vector<std::string> splitRangeText(const std::string& text) {
    const std::u32string u = core::decodeUtf8(core::pythonStrip(text));
    std::vector<std::string> out;
    std::u32string current;
    bool inSeparator = false;
    for (const char32_t c : u) {
        const bool sep = c == U',' || c == U';' || core::isPythonSpace(c);
        if (sep) {
            if (!inSeparator) {
                out.push_back(core::encodeUtf8(current));
                current.clear();
            }
            inSeparator = true;
        } else {
            current += c;
            inSeparator = false;
        }
    }
    out.push_back(core::encodeUtf8(current));
    return out;
}

}  // namespace

// --- residues --------------------------------------------------------------------------------------

std::optional<io::ResidueName> nucleotide(const PdbResidue& residue) {
    std::optional<io::ResidueName> info = io::parseResidueName(residue.name);
    if (!info) {
        return std::nullopt;
    }
    const bool backbone = std::any_of(residue.atoms.begin(), residue.atoms.end(), [](const PdbAtom& a) {
        return a.name.find('\'') != std::string::npos || a.name.find('*') != std::string::npos || a.name == "P";
    });
    if (!backbone) {
        return std::nullopt;
    }
    if (info->sugar == "rna") {  // CH-05 (MDATK 0.1.0): a 2'-deoxy sugar is DNA
        bool c2 = false;
        bool o2 = false;
        for (const PdbAtom& a : residue.atoms) {
            const std::string name = withPrimes(a.name);
            c2 = c2 || name == "C2'";
            o2 = o2 || name == "O2'";
        }
        if (c2 && !o2) {
            info->sugar = "dna";
        }
    }
    return info;
}

bool isPolymer(const PdbResidue& residue) {
    return nucleotide(residue).has_value() || isAmino(residue);
}

bool isSolvent(const PdbResidue& residue) {
    return !isPolymer(residue) && (contains(kWater, io::pythonUpper(residue.name)) || residue.atoms.size() == 1);
}

char residueCode(const PdbResidue& residue) {
    if (const auto nt = nucleotide(residue)) {
        return nt->base;
    }
    const std::string name = io::pythonUpper(residue.name);
    if (const auto code = lookup(kAmino, name)) {
        return *code;
    }
    return lookup(kAminoVariants, name).value_or('X');
}

std::optional<std::string> residueChirality(const PdbResidue& residue) {
    if (residue.chirality && !residue.chirality->empty()) {
        return residue.chirality;
    }
    if (const auto nt = nucleotide(residue)) {
        return std::string(1, nt->chirality);
    }
    return std::nullopt;
}

std::string schemeName(char base, char chirality, std::string_view terminal, std::string_view scheme) {
    const std::string b(1, base);
    if (scheme == "amber") {
        if (chirality == 'L') {  // nucleic12_LRNA.lib defines LX, LX5 and LX3 only
            return "L" + b + (terminal == "5" || terminal == "N" ? "5" : terminal == "3" ? "3" : "");
        }
        return b + std::string(terminal);
    }
    if (scheme == "ccd") {
        return chirality == 'L' ? "0" + b + (terminal == "5" ? "5" : terminal == "3" ? "3" : "") : b;
    }
    if (scheme == "rosetta" || scheme == "pdb") {
        return chirality == 'L' ? "0" + b : b;
    }
    checkScheme(scheme);
    return b;
}

std::vector<std::string> segmentTerminals(const std::vector<PdbResidue>& segment) {
    std::vector<std::string> out(segment.size());
    std::vector<std::size_t> nts;
    for (std::size_t k = 0; k < segment.size(); ++k) {
        if (nucleotide(segment[k])) {
            nts.push_back(k);
        }
    }
    if (nts.empty()) {
        return out;
    }
    const std::size_t k5 = nts.front();
    const std::size_t k3 = nts.back();
    const bool five = !(k5 > 0 && bondedOn(segment[k5 - 1], segment[k5]));
    const bool three = !(k3 + 1 < segment.size() && !phosphateOnly(segment[k3 + 1]) &&
                         bondedOn(segment[k3], segment[k3 + 1]));
    if (nts.size() == 1) {
        out[k5] = endsState(five, three);
    } else {
        out[k5] = endsState(five, false);
        out[k3] = endsState(false, three);
    }
    return out;
}

// --- chirality -------------------------------------------------------------------------------------

double chiralVolume(const Vec3& centre, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 va{};
    Vec3 vb{};
    Vec3 vc{};
    for (std::size_t i = 0; i < 3; ++i) {
        va[i] = a[i] - centre[i];
        vb[i] = b[i] - centre[i];
        vc[i] = c[i] - centre[i];
    }
    const Vec3 cross = {vb[1] * vc[2] - vb[2] * vc[1], vb[2] * vc[0] - vb[0] * vc[2],
                        vb[0] * vc[1] - vb[1] * vc[0]};
    return va[0] * cross[0] + va[1] * cross[1] + va[2] * cross[2];
}

SugarChirality sugarChirality(const AtomCoords& coords, std::optional<char> base) {
    const bool purine = base == 'A' || base == 'G';
    const bool pyrimidine = base == 'C' || base == 'U' || base == 'T';
    const std::string_view glyco = (purine || (!pyrimidine && find(coords, "N9") != nullptr)) ? "N9" : "N1";
    SugarChirality out;
    int lVotes = 0;
    int dVotes = 0;
    int flat = 0;
    for (const ChiralCentre& centre : kLChiralSigns) {
        const Vec3* c = find(coords, centre.centre);
        std::array<const Vec3*, 3> partners{};
        bool complete = c != nullptr;
        for (std::size_t k = 0; k < 3; ++k) {
            partners[k] = find(coords, centre.partners[k] == "N" ? glyco : centre.partners[k]);
            complete = complete && partners[k] != nullptr;
        }
        if (!complete) {
            continue;
        }
        const double v = chiralVolume(*c, *partners[0], *partners[1], *partners[2]);
        out.volumes.emplace_back(std::string(centre.centre), v);
        if (std::abs(v) >= kMinChiralVolume) {
            (v * centre.lSign > 0 ? lVotes : dVotes) += 1;
        } else {
            ++flat;  // a nearly planar centre: distorted geometry
        }
    }
    if (!(lVotes || dVotes) || (lVotes && dVotes) || flat) {
        return out;
    }
    out.call = lVotes > dVotes ? 'L' : 'D';
    return out;
}

// --- ranges ----------------------------------------------------------------------------------------

bool RangeSpec::truthy() const {
    if (const std::string* t = text()) {
        return !t->empty();
    }
    if (const auto* list = items()) {
        return !list->empty();
    }
    return false;
}

std::vector<ResidueRange> parseRanges(const RangeSpec& spec) {
    std::vector<ResidueRange> out;
    if (const std::string* text = spec.text()) {
        for (const std::string& item : splitRangeText(*text)) {
            if (!item.empty()) {
                out.push_back(readRangeText(item));
            }
        }
    } else if (const auto* items = spec.items()) {
        for (const RangeSpec::Item& item : *items) {
            if (const auto* tuple = std::get_if<ResidueRange>(&item)) {
                ResidueRange r{core::pythonStrip(tuple->chain), tuple->start, tuple->end ? tuple->end : tuple->start};
                out.push_back(std::move(r));
            } else if (const std::string& t = std::get<std::string>(item); !t.empty()) {
                out.push_back(readRangeText(t));
            }
        }
    }
    for (const ResidueRange& r : out) {  // CH-15
        if (r.start && r.end && *r.end < *r.start) {
            throw PdbError("residue range " + r.chain + ":" + std::to_string(*r.start) + "-" + std::to_string(*r.end) +
                           " ends before it starts");
        }
    }
    return out;
}

bool inRanges(const PdbResidue& residue, const std::vector<ResidueRange>& ranges) {
    for (const ResidueRange& r : ranges) {
        if (residue.chain == r.chain &&
            (!r.start || (*r.start <= residue.resseq && (!r.end || residue.resseq <= *r.end)))) {
            return true;
        }
    }
    return false;
}

std::string listing(const std::vector<std::string>& items, std::size_t limit) {
    const std::vector<std::string> shown(items.begin(), items.begin() + static_cast<std::ptrdiff_t>(std::min(limit, items.size())));
    std::string out = join(shown, ", ");
    if (items.size() > limit) {
        out += " and " + std::to_string(items.size() - limit) + " more";
    }
    return out;
}

// --- names and records -----------------------------------------------------------------------------

Notes assignChirality(PdbStructure& structure, bool detect, const RangeSpec& lRanges) {
    const std::vector<ResidueRange> ranges = parseRanges(lRanges);
    Notes notes;
    std::vector<const PdbResidue*> lInDNames;
    std::vector<std::string> dInLNames;
    std::vector<std::string> undetermined;
    std::vector<std::string> forcedD;
    std::vector<std::string> lDna;
    std::int64_t countL = 0;
    std::int64_t countD = 0;
    std::int64_t inRange = 0;
    for (PdbResidue* res : allResidues(structure)) {
        const auto nt = nucleotide(*res);
        if (!nt) {
            continue;
        }
        const SugarChirality sugar = detect ? sugarChirality(res->coords(), nt->base) : SugarChirality{};
        const bool forced = !ranges.empty() && inRanges(*res, ranges);
        inRange += forced;
        if (forced) {
            res->chirality = "L";
            if (sugar.call == 'D') {
                forcedD.push_back(res->label() + " (" + volumesText(sugar.volumes) + ")");
            }
        } else {
            res->chirality = std::string(1, sugar.call.value_or(nt->chirality));
            if (detect && !sugar.call) {
                undetermined.push_back(res->label());
            } else if (sugar.call == 'L' && nt->chirality == 'D') {
                lInDNames.push_back(res);
            } else if (sugar.call == 'D' && nt->chirality == 'L') {
                dInLNames.push_back(res->label() + " (" + volumesText(sugar.volumes) + ")");
            }
        }
        if (nt->sugar == "dna" && res->chirality == "L") {
            lDna.push_back(res->label());
        }
        (*res->chirality == "L" ? countL : countD) += 1;
    }
    notes.push_back(std::string("Chirality (") + (detect ? "sugar geometry" : "residue names") + "): " +
                    std::to_string(countL) + " L and " + std::to_string(countD) + " D nucleotides" +
                    (!ranges.empty() ? ", " + std::to_string(inRange) + " forced L by ranges" : ""));
    if (!lInDNames.empty()) {
        notes.push_back("L sugars in " + std::to_string(lInDNames.size()) +
                        " residues with D names (e.g. cpptraj output), L names when a scheme is applied: " +
                        residueRanges(lInDNames));
    }
    if (!dInLNames.empty()) {
        notes.push_back("WARNING: L-named residues with D sugar geometry (D names when a scheme is applied; "
                        "force L with ranges if intended): " + listing(dInLNames));
    }
    if (!forcedD.empty()) {
        notes.push_back("WARNING: forced L although the sugar geometry is D: " + listing(forcedD));
    }
    if (!undetermined.empty()) {
        notes.push_back("Chirality undetermined (missing sugar atoms or centres that disagree), using the residue "
                        "name for " + listing(undetermined));
    }
    if (!lDna.empty()) {
        notes.push_back("L-DNA residues have no name in the supported schemes and keep their names: " + listing(lDna));
    }
    if (!ranges.empty() && !inRange) {
        notes.push_back("WARNING: the L ranges match no nucleotide");
    }
    return notes;
}

Notes convertNames(PdbStructure& structure, std::string_view scheme, const RangeSpec& only, bool terminals) {
    checkScheme(scheme);
    const std::vector<ResidueRange> ranges = parseRanges(only);
    Counter<std::pair<std::string, std::string>> changes;
    std::vector<std::string> termMismatch;
    std::vector<std::string> phosphate;
    std::vector<std::string> singleL;
    std::vector<std::string> keptPhosphate;
    for (auto& segment : structure.segments) {
        const std::vector<std::string> terms = segmentTerminals(segment);
        for (std::size_t k = 0; k < segment.size(); ++k) {
            PdbResidue& res = segment[k];
            const auto nt = nucleotide(res);
            if (!nt || nt->sugar != "rna" || (!ranges.empty() && !inRanges(res, ranges))) {
                continue;
            }
            const std::string segTerm = terminals ? terms[k] : "";
            std::string term = segTerm;
            const char chir = res.chirality && !res.chirality->empty() ? (*res.chirality)[0] : nt->chirality;
            const bool hasP = std::any_of(res.atoms.begin(), res.atoms.end(), [](const PdbAtom& a) { return a.name == "P"; });
            if (scheme == "ccd" && chir == 'L' && (term == "5" || term == "N") && hasP) {
                // 0A5 ... stand for a 5'-OH end; a 5' end with its phosphate is the CCD component itself
                term.clear();
                keptPhosphate.push_back(res.label());
            }
            const std::string name = schemeName(nt->base, chir, term, scheme);
            const bool singleLUnit = segTerm == "N" && chir == 'L' && nt->terminal == "5";  // LA5 for LAN
            if (terminals && !nt->terminal.empty() && nt->terminal != segTerm && !singleLUnit) {
                termMismatch.push_back(res.label());  // the old name shows the misplaced terminal
            }
            if (name != res.name) {
                changes.add({res.name, name});
                res.name = name;
            }
            if (scheme == "amber" && (terms[k] == "5" || terms[k] == "N") && hasP) {
                phosphate.push_back(res.label());
            }
            if (scheme == "amber" && term == "N" && chir == 'L') {
                singleL.push_back(res.label());
            }
        }
    }
    Notes notes;
    if (!changes.items.empty()) {
        auto sorted = changes.items;
        std::stable_sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
        std::vector<std::string> parts;
        for (const auto& [key, n] : sorted) {
            parts.push_back(key.first + "->" + key.second + " x" + std::to_string(n));
        }
        notes.push_back("Renamed " + std::to_string(changes.total()) + " residues to the " + std::string(scheme) +
                        " scheme: " + join(parts, ", "));
    } else {
        notes.push_back("Residue names already follow the " + std::string(scheme) + " scheme");
    }
    if (!termMismatch.empty()) {
        notes.push_back("Terminal names that are not at a segment end (missing TER?): " + listing(termMismatch));
    }
    if (!phosphate.empty()) {
        notes.push_back("WARNING: 5'-terminal residues carry P/OP1/OP2, which the Amber 5' units do not have "
                        "(tleap cannot type them; remove them with --strip-5p or cap the end): " + listing(phosphate));
    }
    if (!keptPhosphate.empty()) {
        notes.push_back("5'-terminal L residues with a 5'-phosphate keep the CCD component names (0A ..., which "
                        "include P/OP1/OP2/OP3; the 5' variants 0A5 ... have a 5'-OH): " + listing(keptPhosphate));
    }
    if (!singleL.empty()) {
        notes.push_back("nucleic12_LRNA.lib has no unit for a single-residue L segment; wrote the 5' unit "
                        "(its O3' lacks a hydrogen): " + listing(singleL));
    }
    return notes;
}

Notes setRecords(PdbStructure& structure, std::string_view policy) {
    if (!contains(kRecordPolicies, policy)) {
        throw PdbError("unknown record policy " + core::pythonStringRepr(std::string(policy)) +
                       " (choose from keep, rosetta, pdb)");
    }
    if (policy == "keep") {
        return {};
    }
    std::map<std::pair<std::string, std::string>, std::int64_t> changed;
    forEachResidue(structure, [&](PdbResidue& res) {
        std::string record;
        if (policy == "rosetta") {
            record = isPolymer(res) ? "ATOM" : "";
        } else {
            const bool standard = lookup(kAmino, res.name).has_value() || contains(kStandardNucleic, res.name);
            record = standard ? "ATOM" : "HETATM";
        }
        if (!record.empty() &&
            std::any_of(res.atoms.begin(), res.atoms.end(), [&](const PdbAtom& a) { return a.record != record; })) {
            changed[{res.name, record}] += 1;
            for (PdbAtom& a : res.atoms) {
                a.record = record;
            }
        }
    });
    if (changed.empty()) {
        return {};
    }
    std::vector<std::string> parts;
    for (const auto& [key, n] : changed) {
        parts.push_back(key.first + "->" + key.second + " x" + std::to_string(n));
    }
    return {"Record policy " + std::string(policy) + ": " + join(parts, ", ")};
}

// --- chains ----------------------------------------------------------------------------------------

Notes restoreChainsFromReference(PdbStructure& structure, const PdbStructure& reference, bool copyNumbers) {
    std::vector<PdbResidue*> mine;
    for (PdbResidue* r : allResidues(structure)) {
        if (!isSolvent(*r)) {
            mine.push_back(r);
        }
    }
    std::vector<const PdbResidue*> theirs;
    for (const PdbResidue* r : reference.residues()) {
        if (!isSolvent(*r)) {
            theirs.push_back(r);
        }
    }
    if (mine.size() != theirs.size()) {
        throw PdbError("chain restore: the structure has " + std::to_string(mine.size()) +
                       " residues but the reference has " + std::to_string(theirs.size()) +
                       " (water and ions excluded)");
    }
    std::vector<std::string> bad;
    for (std::size_t k = 0; k < mine.size(); ++k) {
        const char a = residueCode(*mine[k]);
        const char b = residueCode(*theirs[k]);
        if (a != b || (a == 'X' && io::pythonUpper(mine[k]->name) != io::pythonUpper(theirs[k]->name))) {
            bad.push_back("#" + std::to_string(k + 1) + ": " + mine[k]->label() + " vs reference " + theirs[k]->label());
        }
    }
    if (!bad.empty()) {
        throw PdbError("chain restore: sequence differs from the reference at " + std::to_string(bad.size()) +
                       " residue(s): " + listing(bad, 5));
    }
    for (std::size_t k = 0; k < mine.size(); ++k) {
        mine[k]->chain = theirs[k]->chain;
        if (copyNumbers) {
            mine[k]->resseq = theirs[k]->resseq;
            mine[k]->icode = theirs[k]->icode;
        }
    }
    const std::size_t n = mine.size();
    resegment(structure);
    std::vector<std::string> chains;
    for (const std::string& c : reference.chainIds()) {
        chains.push_back(chainOrBlank(c));
    }
    return {std::string("Restored chain ids") + (copyNumbers ? " and residue numbers" : "") + " of " +
            std::to_string(n) + " residues from the reference (chains " + join(chains, ", ") + ")"};
}

Notes restoreChains(PdbStructure& structure, const std::optional<std::vector<std::string>>& chainIds) {
    const std::vector<std::size_t> segments = polymerSegments(structure);
    std::vector<std::string> ids;
    if (chainIds) {
        for (const std::string& c : *chainIds) {
            ids.push_back(core::pythonStrip(c));
        }
        if (ids.size() != segments.size()) {
            throw PdbError("chain restore: " + std::to_string(ids.size()) + " chain ids given for " +
                           std::to_string(segments.size()) + " TER segments");
        }
        if (std::any_of(ids.begin(), ids.end(), [](const std::string& c) { return core::characterCount(c) > 1; })) {
            throw PdbError("chain ids must be single characters");
        }
    } else {
        if (segments.size() > io::kChainIds.size()) {
            throw PdbError(std::to_string(segments.size()) + " segments but only " +
                           std::to_string(io::kChainIds.size()) + " single-character chain ids");
        }
        for (std::size_t k = 0; k < segments.size(); ++k) {
            ids.emplace_back(1, io::kChainIds[k]);
        }
    }
    for (std::size_t k = 0; k < segments.size(); ++k) {
        for (PdbResidue& res : structure.segments[segments[k]]) {
            res.chain = ids[k];
        }
    }
    resegment(structure);
    std::vector<std::string> shown;
    for (const std::string& c : ids) {
        shown.push_back(chainOrBlank(c));
    }
    return {"Assigned chain ids " + join(shown, ", ") + " to " + std::to_string(segments.size()) + " TER segments"};
}

Notes splitChains(PdbStructure& structure, const RangeSpec& at) {
    FreeIds freeIds(structure);
    Notes notes;
    if (at.isNone()) {
        std::set<std::string> seen;
        for (const std::size_t k : polymerSegments(structure)) {
            auto& segment = structure.segments[k];
            const std::string chain = segment.front().chain;
            if (seen.count(chain)) {
                const std::string id = freeIds.next();
                notes.push_back("Segment starting at " + segment.front().label() + " -> chain " + id);
                for (PdbResidue& res : segment) {
                    res.chain = id;
                }
            }
            seen.insert(chain);
        }
        if (notes.empty()) {
            notes.push_back("No chain spans several TER segments; nothing split");
        }
        return notes;
    }
    // residues to split before, by their input chain and number, as positions in file order (which
    // splitting segments keeps)
    const std::vector<const PdbResidue*> residues = structure.residues();
    std::vector<std::size_t> targets;
    for (const ResidueRange& r : parseRanges(at)) {
        if (!r.start) {
            throw PdbError("split position needs a residue number, e.g. " + r.chain + ":45");
        }
        const auto it = std::find_if(residues.begin(), residues.end(), [&](const PdbResidue* res) {
            return res->chain == r.chain && res->resseq == *r.start;
        });
        if (it == residues.end()) {
            throw PdbError("split position " + r.chain + ":" + std::to_string(*r.start) + " not found");
        }
        targets.push_back(static_cast<std::size_t>(it - residues.begin()));
    }
    std::stable_sort(targets.begin(), targets.end());
    for (const std::size_t target : targets) {
        std::size_t i = 0;
        std::size_t k = target;
        while (k >= structure.segments[i].size()) {
            k -= structure.segments[i].size();
            ++i;
        }
        const std::string chain = structure.segments[i][k].chain;
        const std::int64_t start = structure.segments[i][k].resseq;
        if (k > 0) {
            auto& segment = structure.segments[i];
            std::vector<PdbResidue> rest(std::make_move_iterator(segment.begin() + static_cast<std::ptrdiff_t>(k)),
                                         std::make_move_iterator(segment.end()));
            segment.erase(segment.begin() + static_cast<std::ptrdiff_t>(k), segment.end());
            structure.segments.insert(structure.segments.begin() + static_cast<std::ptrdiff_t>(i) + 1, std::move(rest));
            ++i;
        }
        const std::string id = freeIds.next();
        for (std::size_t s = i; s < structure.segments.size(); ++s) {
            for (PdbResidue& res : structure.segments[s]) {
                if (res.chain == chain) {
                    res.chain = id;
                }
            }
        }
        notes.push_back("Split chain " + chainOrBlank(chain) + " at residue " + std::to_string(start) +
                        ": the rest is chain " + id);
    }
    return notes;
}

Notes combineChains(PdbStructure& structure, const std::vector<std::string>& chainList,
                    const std::optional<std::string>& into) {
    if (into && core::characterCount(core::pythonStrip(*into)) > 1) {  // CH-03
        throw PdbError("chain ids must be single characters");
    }
    std::vector<std::string> chains;
    for (const std::string& c : chainList) {
        chains.push_back(core::pythonStrip(c));
    }
    if (chains.empty()) {
        return {};
    }
    const std::string target = into ? core::pythonStrip(*into) : chains.front();
    const std::vector<std::string> present = structure.chainIds();
    std::vector<std::string> missing;
    for (const std::string& c : chains) {
        if (std::find(present.begin(), present.end(), c) == present.end()) {
            missing.push_back(c);
        }
    }
    std::int64_t n = 0;
    forEachResidue(structure, [&](PdbResidue& res) {
        if (std::find(chains.begin(), chains.end(), res.chain) != chains.end()) {
            res.chain = target;
            ++n;
        }
    });
    std::vector<std::string> listed;
    for (const std::string& c : chains) {
        listed.push_back(chainOrBlank(c));
    }
    Notes notes = {"Combined chains " + join(listed, ", ") + " into chain " + chainOrBlank(target) + " (" +
                   std::to_string(n) + " residues)"};
    if (!missing.empty()) {
        notes.push_back("WARNING: chains not present: " + join(missing, ", "));
    }
    Counter<std::tuple<std::string, std::int64_t, std::string>> keys;
    for (const PdbResidue* r : structure.residues()) {
        if (r->chain == target) {
            keys.add({r->chain, r->resseq, r->icode});
        }
    }
    std::vector<std::int64_t> dup;
    for (const auto& [key, count] : keys.items) {
        if (count > 1) {
            dup.push_back(std::get<1>(key));
        }
    }
    std::sort(dup.begin(), dup.end());
    if (!dup.empty()) {
        std::vector<std::string> numbers;
        for (const std::int64_t d : dup) {
            numbers.push_back(std::to_string(d));
        }
        notes.push_back("WARNING: chain " + chainOrBlank(target) + " now has duplicate residue numbers (" +
                        listing(numbers) + "); renumber the residues");
    }
    return notes;
}

Notes orderChains(PdbStructure& structure, const std::optional<std::vector<std::string>>& order) {
    const std::vector<std::string> before = structure.chainIds();
    std::map<std::string, std::int64_t> rank;
    if (order) {
        for (std::size_t k = 0; k < order->size(); ++k) {
            rank[(*order)[k]] = static_cast<std::int64_t>(k);  // a repeated id: its last place
        }
    }
    using Key = std::tuple<int, int, std::string, std::int64_t>;
    std::vector<std::pair<Key, std::size_t>> keys;
    for (std::size_t k = 0; k < structure.segments.size(); ++k) {
        const auto& seg = structure.segments[k];
        const auto pos = static_cast<std::int64_t>(k);
        Key key;
        if (allSolvent(seg)) {
            key = {2, 0, "", pos};
        } else if (!order) {
            const std::string& c = seg.front().chain;
            key = {0, c.empty() ? 1 : 0, c, pos};
        } else if (const auto it = rank.find(seg.front().chain); it != rank.end()) {
            key = {0, 0, "", it->second};
        } else {
            key = {1, 0, "", pos};
        }
        keys.emplace_back(std::move(key), k);
    }
    std::stable_sort(keys.begin(), keys.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<std::vector<PdbResidue>> segments;
    for (const auto& [key, k] : keys) {
        segments.push_back(std::move(structure.segments[k]));
    }
    structure.segments = std::move(segments);
    const std::vector<std::string> after = structure.chainIds();
    if (after == before) {
        return {};
    }
    std::vector<std::string> shown;
    for (const std::string& c : after) {
        shown.push_back(chainOrBlank(c));
    }
    return {"Chain order: " + join(shown, " ")};
}

// --- numbering, phosphates, sequence ---------------------------------------------------------------

Notes strip5primePhosphates(PdbStructure& structure) {
    std::vector<std::string> changed;
    std::vector<std::string> skipped;
    std::vector<std::string> emptied;
    for (auto& seg : structure.segments) {
        const auto firstPolymer = [&seg]() {
            return std::find_if(seg.begin(), seg.end(), [](const PdbResidue& r) { return !isSolvent(r); });
        };
        auto res = firstPolymer();
        while (res != seg.end() && phosphateOnly(*res)) {
            emptied.push_back(res->label());  // a lone 5'-phosphate modelled as its own residue
            seg.erase(res);
            res = firstPolymer();
        }
        if (res == seg.end()) {
            continue;
        }
        if (!nucleotide(*res)) {
            skipped.push_back(res->label());  // the segment starts with something else: its phosphates stay
            continue;
        }
        if (std::none_of(res->atoms.begin(), res->atoms.end(), isPhosphateAtom)) {
            continue;
        }
        const auto p = std::find_if(res->atoms.begin(), res->atoms.end(), [](const PdbAtom& a) {
            return isPhosphateAtom(a) && a.name == "P";
        });
        const std::optional<PdbAtom> pAtom = p != res->atoms.end() ? std::optional<PdbAtom>(*p) : std::nullopt;
        const bool hasH = std::any_of(res->atoms.begin(), res->atoms.end(), [](const PdbAtom& a) {
            const std::size_t first = a.name.find_first_not_of("0123456789");
            return first != std::string::npos && a.name[first] == 'H';
        });
        std::erase_if(res->atoms, isPhosphateAtom);
        const auto o5 = std::find_if(res->atoms.begin(), res->atoms.end(), [](const PdbAtom& a) {
            return withPrimes(a.name) == "O5'";
        });
        const bool hasHo5 = std::any_of(res->atoms.begin(), res->atoms.end(), [](const PdbAtom& a) {
            return a.name == "HO5'";
        });
        if (hasH && pAtom && o5 != res->atoms.end() && !hasHo5) {
            const Vec3 o = o5->xyz();
            const Vec3 pXyz = pAtom->xyz();
            Vec3 d{};
            for (std::size_t k = 0; k < 3; ++k) {
                d[k] = pXyz[k] - o[k];
            }
            double n = std::sqrt(0.0 + d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            if (n == 0.0) {
                n = 1.0;
            }
            std::string tail = "   ";
            for (std::size_t k = 0; k < 3; ++k) {
                tail += formatFixed(o[k] + 0.96 * d[k] / n, 8, 3, false);
            }
            tail += codePointsFrom(o5->tail, 27);
            PdbAtom h;
            h.record = o5->record;
            h.name = "HO5'";
            h.tail = std::move(tail);
            h.element = "H";
            res->atoms.insert(o5 + 1, std::move(h));
        }
        changed.push_back(res->label());
    }
    Notes notes;
    if (!changed.empty()) {
        notes.push_back("Removed the 5'-terminal phosphate of " + std::to_string(changed.size()) +
                        " residue(s): " + listing(changed));
    }
    if (!emptied.empty()) {
        notes.push_back("Removed " + std::to_string(emptied.size()) +
                        " 5'-terminal residue(s) that held only a phosphate: " + listing(emptied));
    }
    if (!skipped.empty()) {
        notes.push_back("5' phosphates kept where a segment does not start with a recognised nucleotide: " +
                        listing(skipped));
    }
    return notes;
}

Notes renumberResidues(PdbStructure& structure, std::int64_t start, bool perChain, bool keepGaps) {
    std::map<std::optional<std::string>, std::int64_t> last;
    const PdbResidue* prev = nullptr;
    std::string prevChain;
    std::pair<std::int64_t, std::string> prevOld{0, ""};
    std::int64_t cleared = 0;
    std::int64_t disorder = 0;
    std::int64_t highest = 0;
    bool any = false;
    forEachResidue(structure, [&](PdbResidue& res) {
        const std::optional<std::string> key = perChain ? std::optional<std::string>(res.chain) : std::nullopt;
        const std::pair<std::int64_t, std::string> old{res.resseq, res.icode};
        const bool sameChain = prev != nullptr && prevChain == res.chain;
        std::int64_t number = 0;
        const auto it = last.find(key);
        if (it == last.end()) {
            number = start;
        } else if (keepGaps && sameChain && old > prevOld) {
            number = it->second + std::max<std::int64_t>(0, old.first - prevOld.first);
        } else {
            number = it->second + 1;
            disorder += keepGaps && sameChain;
        }
        if (!keepGaps && !res.icode.empty()) {
            ++cleared;
            res.icode.clear();
        }
        prev = &res;
        prevChain = res.chain;
        prevOld = old;
        res.resseq = number;
        last[key] = number;
        highest = any ? std::max(highest, number) : number;
        any = true;
    });
    Notes notes = {std::string("Renumbered residues ") + (perChain ? "per chain" : "continuously") + " from " +
                   std::to_string(start) + (keepGaps ? " (gaps kept)" : "")};
    if (cleared) {
        notes.push_back("Cleared " + std::to_string(cleared) + " insertion codes");
    }
    if (disorder) {
        notes.push_back("WARNING: " + std::to_string(disorder) + " residues out of order; numbered sequentially there");
    }
    if ((any ? highest : 0) > 9999) {
        notes.push_back("WARNING: residue numbers above 9999 wrap in the PDB columns");
    }
    return notes;
}

std::vector<std::pair<std::string, std::vector<const PdbResidue*>>> sequenceChains(const PdbStructure& structure) {
    std::vector<std::pair<std::string, std::vector<const PdbResidue*>>> chains;
    for (const std::size_t k : polymerSegments(structure)) {
        const auto& seg = structure.segments[k];
        if (std::none_of(seg.begin(), seg.end(), [](const PdbResidue& r) { return isPolymer(r); })) {
            continue;  // ligand-only segment
        }
        for (const PdbResidue& res : seg) {
            if (isSolvent(res)) {
                continue;
            }
            auto it = std::find_if(chains.begin(), chains.end(), [&](const auto& c) { return c.first == res.chain; });
            if (it == chains.end()) {
                chains.emplace_back(res.chain, std::vector<const PdbResidue*>{});
                it = chains.end() - 1;
            }
            it->second.push_back(&res);
        }
    }
    return chains;
}

SeqresNames seqresNames(const PdbStructure& structure) {
    SeqresNames out;
    std::map<std::pair<std::string, std::string>, std::int64_t> mapped;
    std::vector<std::string> unknown;
    for (const auto& [chain, residues] : sequenceChains(structure)) {
        io::SeqresChain entry{chain, {}};
        for (const PdbResidue* res : residues) {
            std::string name = core::pythonStrip(res->name);
            const auto nt = nucleotide(*res);
            const std::u32string up = core::decodeUtf8(io::pythonUpper(name));
            const bool amberLDna = nt && nt->sugar == "dna" && nt->chirality == 'L' && up.substr(0, 2) == U"LD";
            if (core::characterCount(name) > 3 || amberLDna) {
                std::string replacement;
                if (const auto info = io::parseResidueName(name)) {
                    replacement = ccdComponent(*info);
                    mapped[{name, replacement}] += 1;
                } else {
                    replacement = "UNK";
                    unknown.push_back(res->label());
                }
                name = std::move(replacement);
            }
            entry.names.push_back(std::move(name));
        }
        out.chains.push_back(std::move(entry));
    }
    if (!mapped.empty()) {
        std::vector<std::string> parts;
        for (const auto& [key, n] : mapped) {
            parts.push_back(key.first + "->" + key.second + " x" + std::to_string(n));
        }
        out.notes.push_back("SEQRES names are the CCD's three-character ones: wrote the CCD components of MDATK's "
                            "Amber-style L-DNA names and of other four-character names " + join(parts, ", "));
    }
    if (!unknown.empty()) {
        out.notes.push_back("WARNING: SEQRES names have three characters: wrote UNK for the residue name(s) of " +
                            listing(unknown));
    }
    return out;
}

std::vector<std::string> seqresRecords(const PdbStructure& structure) {
    return io::seqresLines(seqresNames(structure).chains);
}

std::string fasta(const PdbStructure& structure, std::string_view title) {
    std::vector<std::string> out;
    for (const auto& [chain, residues] : sequenceChains(structure)) {
        std::vector<const PdbResidue*> nts;
        std::set<std::string> sugars;
        std::int64_t nAmino = 0;
        for (const PdbResidue* r : residues) {
            if (const auto nt = nucleotide(*r)) {
                nts.push_back(r);
                sugars.insert(nt->sugar);
            }
            nAmino += isAmino(*r);
        }
        std::string kind;
        if (!nts.empty() && !nAmino) {
            kind = sugars == std::set<std::string>{"rna"}   ? "RNA"
                   : sugars == std::set<std::string>{"dna"} ? "DNA"
                                                            : "nucleic acid";
        } else {
            kind = nAmino && nts.empty() ? "protein" : "polymer";
        }
        std::vector<std::int64_t> lNumbers;
        for (const PdbResidue* r : nts) {
            if (residueChirality(*r) == "L") {
                lNumbers.push_back(r->resseq);
            }
        }
        if (!lNumbers.empty() && lNumbers.size() == nts.size()) {
            kind = "L-" + kind;
        } else if (!lNumbers.empty()) {
            kind += ", L-nucleotides at " + compress(lNumbers);
        }
        out.push_back(">" + std::string(title) + "_" + chainOrBlank(chain) + " chain " + chainOrBlank(chain) + " | " +
                      kind + " | " + std::to_string(residues.size()) + " residues");
        std::string sequence;
        for (const PdbResidue* r : residues) {
            sequence += residueCode(*r);
        }
        out.push_back(std::move(sequence));
    }
    return join(out, "\n") + "\n";
}

std::string writePdb(const PdbStructure& structure, bool stripConect, bool seqres) {
    if (!seqres) {
        return io::writePdb(structure, stripConect);
    }
    return io::writePdb(structure, stripConect, seqresNames(structure).chains);
}

core::Value summary(const PdbStructure& structure) {
    core::Table chains;
    for (const auto& [chain, residues] : sequenceChains(structure)) {
        Counter<std::string> calls;
        std::set<std::string> names;
        std::string sequence;
        for (const PdbResidue* r : residues) {
            if (nucleotide(*r)) {
                calls.add(residueChirality(*r).value_or("?"));
            }
            names.insert(r->name);
            sequence += residueCode(*r);
        }
        core::Table chirality;
        for (const auto& [call, n] : calls.items) {
            chirality.emplace_back(call, core::Value(n));
        }
        core::Array nameList(names.begin(), names.end());
        core::Table entry = {
            {"residues", core::Value(static_cast<std::int64_t>(residues.size()))},
            {"numbers", core::Value(residues.empty() ? std::string()
                                                     : std::to_string(residues.front()->resseq) + "-" +
                                                           std::to_string(residues.back()->resseq))},
            {"sequence", core::Value(sequence)},
            {"names", core::Value(std::move(nameList))},
            {"chirality", core::Value(std::move(chirality))},
        };
        const std::string key = chainOrBlank(chain);
        const auto it = std::find_if(chains.begin(), chains.end(), [&](const auto& kv) { return kv.first == key; });
        if (it != chains.end()) {
            it->second = core::Value(std::move(entry));  // a blank chain and a chain "_" share the key
        } else {
            chains.emplace_back(key, core::Value(std::move(entry)));
        }
    }
    Counter<std::string> other;
    for (const PdbResidue* r : structure.residues()) {
        if (!isPolymer(*r)) {
            other.add(r->name);
        }
    }
    core::Table otherTable;
    for (const auto& [name, n] : other.items) {
        otherTable.emplace_back(name, core::Value(n));
    }
    core::Array notes(structure.notes.begin(), structure.notes.end());
    return core::Value(core::Table{
        {"chains", core::Value(std::move(chains))},
        {"segments", core::Value(static_cast<std::int64_t>(structure.segments.size()))},
        {"other_residues", core::Value(std::move(otherTable))},
        {"conect_records", core::Value(static_cast<std::int64_t>(structure.conect.size()))},
        {"notes", core::Value(std::move(notes))},
    });
}

// --- whole files -----------------------------------------------------------------------------------

PdbSource PdbSource::text(std::string text) {
    PdbSource s;
    s.m_text = std::move(text);
    return s;
}

PdbSource PdbSource::file(std::filesystem::path path) {
    PdbSource s;
    s.m_isFile = true;
    s.m_path = std::move(path);
    return s;
}

io::PdbStructure PdbSource::read() const {
    return m_isFile ? io::readPdbFile(m_path) : io::readPdbText(m_text);
}

std::string PdbSource::content() const {
    return m_isFile ? io::readPdbFileText(m_path) : m_text;
}

std::string PdbSource::defaultTitle() const {
    if (!m_isFile) {
        return "structure";
    }
    return core::posixpath::splitext(core::pathToUtf8(m_path.filename())).first;
}

CleanResult clean(const PdbSource& source, const CleanOptions& o) {
    PdbStructure s = source.read();
    Notes notes = s.notes;
    const auto stale = std::count_if(s.header.begin(), s.header.end(),
                                     [](const std::string& line) { return line.starts_with("REMARK 950"); });
    if (stale) {  // base-pair annotation of the old names / numbers
        std::erase_if(s.header, [](const std::string& line) { return line.starts_with("REMARK 950"); });
        notes.push_back("Dropped " + std::to_string(stale) + " REMARK 950 base-pair annotation lines (they name "
                        "the input residues); annotate the cleaned file again if needed");
    }
    const auto append = [&notes](Notes more) { notes.insert(notes.end(), more.begin(), more.end()); };
    if (o.chainReference) {
        append(restoreChainsFromReference(s, o.chainReference->read(), o.copyNumbers));
    } else if (o.chainIds || o.assignChains) {
        append(restoreChains(s, o.chainIds));
    }
    if (o.detectL || o.lRanges.truthy()) {
        append(assignChirality(s, o.detectL, o.lRanges));
    }
    if (o.splitSegments) {
        append(splitChains(s));
    }
    if (o.splitAt.truthy()) {
        append(splitChains(s, o.splitAt));
    }
    if (!o.combine.empty()) {
        append(combineChains(s, o.combine, o.combineInto));
    }
    if (!o.chainOrder.empty()) {
        append(orderChains(s, o.chainOrder));
    } else if (o.sortChains) {
        append(orderChains(s));
    }
    if (o.strip5p) {
        append(strip5primePhosphates(s));
    }
    const std::string scheme = o.scheme.value_or("");
    if (!scheme.empty()) {
        append(convertNames(s, scheme, {}, o.terminalNames));
    }
    if (o.renumber) {
        append(renumberResidues(s, o.start, o.perChain, o.keepGaps));
    }
    std::string policy = o.records.value_or("");
    if (policy.empty()) {
        policy = scheme == "rosetta" ? "rosetta" : scheme == "pdb" ? "pdb" : "keep";
    }
    append(setRecords(s, policy));
    if (!s.conect.empty()) {
        if (o.stripConect) {
            notes.push_back("Stripped " + std::to_string(s.conect.size()) + " CONECT records");
        } else {
            std::set<std::int64_t> serials;
            bool repeated = false;
            for (const PdbResidue* r : s.residues()) {
                for (const PdbAtom& a : r->atoms) {
                    if (a.serial) {
                        repeated = repeated || !serials.insert(*a.serial).second;
                    }
                }
            }
            notes.push_back(repeated ? "WARNING: dropped " + std::to_string(s.conect.size()) +
                                           " CONECT records: the input's atom serials repeat, so they cannot be remapped"
                                     : "Kept CONECT records (remapped to the new atom serials)");  // CH-04
        }
    }
    std::size_t nAtoms = 0;
    for (const PdbResidue* r : s.residues()) {
        nAtoms += r->atoms.size();
    }
    if (nAtoms + s.segments.size() > 99999) {
        notes.push_back("WARNING: more than 99999 atoms; serial numbers wrap in the PDB columns");
    }
    if (o.seqres) {
        notes.push_back("Wrote SEQRES records for " + std::to_string(sequenceChains(s).size()) + " chain(s)");
        append(seqresNames(s).notes);  // CH-07
    }
    std::string title = o.title ? (o.title->empty() ? std::string("structure") : *o.title) : source.defaultTitle();
    CleanResult result;
    result.text = chem::writePdb(s, o.stripConect, o.seqres);
    if (o.makeFasta) {
        result.fasta = fasta(s, title);
    }
    result.notes = std::move(notes);
    result.structure = std::move(s);
    return result;
}

CleanResult toLrna(const PdbSource& source, const RangeSpec& ranges, std::string_view scheme, bool stripConect,
                   bool terminals) {
    PdbStructure s = source.read();
    Notes notes = s.notes;
    Notes more = assignChirality(s, true, ranges);
    notes.insert(notes.end(), more.begin(), more.end());
    more = convertNames(s, scheme, ranges, terminals);
    notes.insert(notes.end(), more.begin(), more.end());
    CleanResult result;
    result.text = chem::writePdb(s, stripConect);
    result.notes = std::move(notes);
    result.structure = std::move(s);
    return result;
}

core::Value check(const PdbSource& source, const RangeSpec& lRanges) {
    PdbStructure s = source.read();
    const Notes notes = assignChirality(s, true, lRanges);
    core::Value info = summary(s);
    core::Array* all = info.find("notes")->get<core::Array>();
    all->insert(all->end(), notes.begin(), notes.end());
    return info;
}

}  // namespace mdbridge::chem
