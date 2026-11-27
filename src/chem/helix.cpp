// MDATK's structure builder, part 1 (mdatk/builder.py v0.0.4.1, MDATK's own code, ported): the
// Arnott A-RNA templates and their hydrogens, the sequence functions, check_helix and the
// Watson-Crick build_helix. Arithmetic in NumPy's order (no fused multiply-adds: CMakeLists.txt), so
// that the written coordinates are MDATK's to the last printed digit. The fibre monomer comes from
// the generated reference tables (data/sources/helix/arnott_arna.tsv, reference_tables.hpp).

#include "mdbridge/chem/helix.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/chem/reference_tables.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "helix_internal.hpp"

namespace mdbridge::chem {

namespace helix_detail {

// --- templates -------------------------------------------------------------------------------

constexpr std::string_view kBackboneNames[] = {"P",   "OP1", "OP2", "O5'", "C5'", "C4'",
                                               "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"};
constexpr std::string_view kBaseA[] = {"N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"};
constexpr std::string_view kBaseC[] = {"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"};
constexpr std::string_view kBaseG[] = {"N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"};
constexpr std::string_view kBaseU[] = {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"};

enum class HRule { Ch, Ch2, Sp2, Nh, Nh2, Oh };

double bondH(HRule rule) {
    switch (rule) {
        case HRule::Ch:
        case HRule::Ch2:
            return 1.09;
        case HRule::Sp2:
            return 1.08;
        case HRule::Nh:
        case HRule::Nh2:
            return 1.01;
        case HRule::Oh:
            return 0.96;
    }
    return 0.0;
}

struct HydrogenRule {
    std::string_view name;
    std::string_view parent;
    HRule rule;
    std::string_view refs[3];  // "N": the glycosidic nitrogen
    std::size_t refCount;
};

// _H_RULES; guanine's amino hydrogens follow the CCD (CHD-01, D15): H21 parallel to N1 -> C2 (trans
// to N1, on the N3 side), H22 parallel to N3 -> C2 (cis to N1, the Watson-Crick edge)
constexpr HydrogenRule kSugarH[] = {
    {"H5'", "C5'", HRule::Ch2, {"O5'", "C4'", ""}, 2},  {"H5''", "C5'", HRule::Ch2, {"O5'", "C4'", ""}, 2},
    {"H4'", "C4'", HRule::Ch, {"C5'", "O4'", "C3'"}, 3}, {"H3'", "C3'", HRule::Ch, {"C4'", "O3'", "C2'"}, 3},
    {"H2'", "C2'", HRule::Ch, {"C3'", "O2'", "C1'"}, 3}, {"HO2'", "O2'", HRule::Oh, {"C2'", "C3'", ""}, 2},
    {"H1'", "C1'", HRule::Ch, {"O4'", "C2'", "N"}, 3}};
constexpr HydrogenRule kHydA[] = {{"H8", "C8", HRule::Sp2, {"N7", "N9", ""}, 2},
                                  {"H61", "N6", HRule::Nh2, {"C6", "C5", ""}, 2},
                                  {"H62", "N6", HRule::Nh2, {"C6", "N1", ""}, 2},
                                  {"H2", "C2", HRule::Sp2, {"N1", "N3", ""}, 2}};
constexpr HydrogenRule kHydG[] = {{"H8", "C8", HRule::Sp2, {"N7", "N9", ""}, 2},
                                  {"H1", "N1", HRule::Nh, {"C2", "C6", ""}, 2},
                                  {"H21", "N2", HRule::Nh2, {"C2", "N1", ""}, 2},
                                  {"H22", "N2", HRule::Nh2, {"C2", "N3", ""}, 2}};
constexpr HydrogenRule kHydC[] = {{"H41", "N4", HRule::Nh2, {"C4", "C5", ""}, 2},
                                  {"H42", "N4", HRule::Nh2, {"C4", "N3", ""}, 2},
                                  {"H5", "C5", HRule::Sp2, {"C4", "C6", ""}, 2},
                                  {"H6", "C6", HRule::Sp2, {"C5", "N1", ""}, 2}};
constexpr HydrogenRule kHydU[] = {{"H3", "N3", HRule::Nh, {"C2", "C4", ""}, 2},
                                  {"H5", "C5", HRule::Sp2, {"C4", "C6", ""}, 2},
                                  {"H6", "C6", HRule::Sp2, {"C5", "N1", ""}, 2}};

constexpr std::string_view kSugarHNames[] = {"H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'"};
constexpr std::string_view kHNamesA[] = {"H8", "H61", "H62", "H2"};
constexpr std::string_view kHNamesG[] = {"H8", "H1", "H21", "H22"};
constexpr std::string_view kHNamesC[] = {"H41", "H42", "H5", "H6"};
constexpr std::string_view kHNamesU[] = {"H3", "H5", "H6"};

constexpr double kHo2Torsion = 180.0;  // C3'-C2'-O2'-HO2'

double wobbleTurn(char base) {
    return base == 'G' ? -13.3 : base == 'U' ? 11.4 : 0.0;
}

std::span<const HydrogenRule> hydrogenRules(char base) {
    switch (base) {
        case 'A':
            return kHydA;
        case 'C':
            return kHydC;
        case 'G':
            return kHydG;
        case 'U':
            return kHydU;
        default:
            return {};
    }
}

Glycosidic glycosidic(char base) {
    if (base == 'A' || base == 'G') {
        return {"N9", "C4", "C8"};
    }
    return {"N1", "C2", "C6"};
}

// _hydrogens: ideal hydrogens of a D-nucleotide
Atoms hydrogens(const Atoms& heavy, char base) {
    Atoms out;
    const std::string_view glyco = glycosidic(base).n;
    std::vector<HydrogenRule> rules(std::begin(kSugarH), std::end(kSugarH));
    for (const HydrogenRule& rule : hydrogenRules(base)) {
        rules.push_back(rule);
    }
    for (const HydrogenRule& rule : rules) {
        const Vec3& c = heavy.at(rule.parent);
        Vec3 ref[3]{};
        for (std::size_t k = 0; k < rule.refCount; ++k) {
            ref[k] = heavy.at(rule.refs[k] == "N" ? glyco : rule.refs[k]);
        }
        const double d = bondH(rule.rule);
        switch (rule.rule) {
            case HRule::Ch:
            case HRule::Sp2:
            case HRule::Nh: {
                Vec3 sum = unit(ref[0] - c);
                for (std::size_t k = 1; k < rule.refCount; ++k) {
                    sum = sum + unit(ref[k] - c);
                }
                out.set(rule.name, c - d * unit(sum));
                break;
            }
            case HRule::Nh2:  // parallel to the ring bond refs[1] -> refs[0]
                out.set(rule.name, c + d * unit(ref[0] - ref[1]));
                break;
            case HRule::Oh:
                out.set(rule.name, place(ref[1], ref[0], c, d, 109.5, kHo2Torsion));
                break;
            case HRule::Ch2: {
                if (out.find(rule.name)) {
                    break;
                }
                const Vec3 u1 = unit(ref[0] - c);
                const Vec3 u2 = unit(ref[1] - c);
                const Vec3 bis = -unit(u1 + u2);
                const Vec3 perp = unit(cross(u1, u2));
                const double half = radians(109.5 / 2);
                Vec3 h1 = c + d * (std::cos(half) * bis + std::sin(half) * perp);
                Vec3 h2 = c + d * (std::cos(half) * bis - std::sin(half) * perp);
                // PDB v3 names of a D-ribose: H5' gives a positive volume C5'(O5', C4', H5')
                if (helixChiralVolume(c, ref[0], ref[1], h1) < 0) {
                    std::swap(h1, h2);
                }
                out.set("H5'", h1);
                out.set("H5''", h2);
                break;
            }
        }
    }
    return out;
}

char checkedBase(std::string_view base) {
    if (base.size() == 1) {
        const char b = base[0] >= 'a' && base[0] <= 'z' ? static_cast<char>(base[0] - 'a' + 'A') : base[0];
        if (b == 'A' || b == 'C' || b == 'G' || b == 'U') {
            return b;
        }
    }
    std::string upper(base);
    for (char& ch : upper) {
        if (ch >= 'a' && ch <= 'z') {
            ch = static_cast<char>(ch - 'a' + 'A');
        }
    }
    throw BuildError("no template for base " + core::pythonStringRepr(upper) + " (A, C, G, U)");
}

// _template: heavy atoms, then hydrogens
const std::vector<TemplateAtom>& fullTemplate(char base) {
    static const std::map<char, std::vector<TemplateAtom>> cache = [] {
        std::map<char, std::vector<TemplateAtom>> out;
        for (const char b : std::string_view("ACGU")) {
            Atoms heavy;
            for (const reference::CylindricalAtom& a : reference::arnottBackbone()) {
                heavy.set(a.name, cartesian(a.r, a.phi, a.z));
            }
            for (const reference::CylindricalBaseAtom& a : reference::arnottBase(b)) {
                heavy.set(a.name, cartesian(a.r, a.phi, a.z));
            }
            std::vector<TemplateAtom> all = heavy.items;
            for (const TemplateAtom& h : hydrogens(heavy, b).items) {
                all.push_back(h);
            }
            out[b] = std::move(all);
        }
        return out;
    }();
    return cache.at(base);
}

Atoms templateAtoms(char base, bool withHydrogens) {
    Atoms out;
    for (const TemplateAtom& a : fullTemplate(base)) {
        if (withHydrogens || a.name[0] != 'H') {
            out.items.push_back(a);
        }
    }
    return out;
}

// _turn_base: the base turned in its plane about C1' (in place)
void turnBase(Atoms& coords, char base, double degrees, bool withHydrogens) {
    if (degrees == 0.0) {
        return;
    }
    const Glycosidic g = glycosidic(base);
    const Vec3& n = coords.at(g.n);
    const Vec3 normal = unit(cross(coords.at(g.a) - n, coords.at(g.b) - n));
    const Mat3 r = turnMatrix(normal, degrees);
    const Vec3 c1 = coords.at("C1'");
    std::vector<std::string_view> names(baseAtomNames(base).begin(), baseAtomNames(base).end());
    if (withHydrogens) {
        for (const std::string_view h : baseHydrogenNames(base)) {
            names.push_back(h);
        }
    }
    for (const std::string_view name : names) {
        if (Vec3* v = coords.find(name)) {
            *v = matVec(r, *v - c1) + c1;
        }
    }
}

// --- text helpers ----------------------------------------------------------------------------

std::string lower(std::string_view text) {
    return core::pythonLower(text);
}
std::string strip(std::string_view text) {
    return core::pythonStrip(text);
}
std::string upperAscii(std::string_view text) {
    std::string out(text);
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return out;
}
std::string repr(std::string_view text) {
    return core::pythonStringRepr(std::string(text));
}
std::string fixed(double x, int precision) {
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%.*f", precision, x);
    return buffer;
}
template <class Items>
std::string join(const Items& items, std::string_view sep) {
    std::string out;
    bool first = true;
    for (const auto& item : items) {
        if (!first) {
            out += sep;
        }
        out += item;
        first = false;
    }
    return out;
}
// ", ".join(items[:8]) + (" ..." if more)
std::string firstEight(const std::vector<std::string>& items) {
    std::vector<std::string> shown(items.begin(), items.begin() + static_cast<std::ptrdiff_t>(std::min<std::size_t>(8, items.size())));
    return join(shown, ", ") + (items.size() > 8 ? " ..." : "");
}

// --- pairings --------------------------------------------------------------------------------

bool parallel(std::string_view pairing) {
    return pairing == "reverse-wc" || pairing == "hoogsteen" || pairing == "hh";
}

std::string_view pairingLw(std::string_view pairing) {
    if (pairing == "wc") return "cWW";
    if (pairing == "reverse-wc") return "tWW";
    if (pairing == "hoogsteen") return "cHW";
    if (pairing == "reverse-hoogsteen") return "tHW";
    return "tHH";
}

std::string_view pairingLabel(std::string_view pairing) {
    if (pairing == "wc") return "Watson-Crick (cis WC/WC), antiparallel A-form";
    if (pairing == "reverse-wc") return "reverse Watson-Crick (trans WC/WC), parallel";
    if (pairing == "hoogsteen") return "Hoogsteen (cis WC/Hoogsteen: A·U, G·C+), parallel";
    if (pairing == "reverse-hoogsteen")
        return "reverse Hoogsteen (trans WC/Hoogsteen: A·U, A·A, G·G), antiparallel";
    return "trans Hoogsteen/Hoogsteen A·A (as PDB 4JRD), parallel";
}

std::string_view triplexLabel(std::string_view motif) {
    if (motif == "hoogsteen")
        return "pyrimidine motif: Hoogsteen third strand parallel to the purine strand (U·A-U, C+·G-C)";
    return "purine motif: reverse-Hoogsteen third strand antiparallel to the purine strand (A·A-U, G·G-C, also U·A-U)";
}

std::string_view referencePair(std::string_view pairing) {
    return pairing == "hh" ? "AA" : "AU";
}

char complementBase(char b) {
    switch (b) {
        case 'A':
            return 'U';
        case 'U':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        default:
            return b;
    }
}

char autoPartnerBase(std::string_view pairing, char b) {
    if (pairing == "reverse-hoogsteen") {
        return b == 'A' ? 'U' : b == 'G' ? 'G' : b == 'U' ? 'A' : 'G';
    }
    if (pairing == "hh") {
        return b;
    }
    return complementBase(b);
}

std::string autoPartnerText(std::string_view pairing) {
    // ", ".join(f"{a}·{b}" for a, b in AUTO_PARTNER[pairing].items()) in the dict's order
    const std::string_view order = pairing == "wc" || pairing == "reverse-wc" || pairing == "hoogsteen" ? "AUGC" : "AGUC";
    std::vector<std::string> parts;
    for (const char b : (pairing == "hh" ? std::string_view("ACGU") : order)) {
        parts.push_back(std::string(1, b) + "·" + std::string(1, autoPartnerBase(pairing, b)));
    }
    return join(parts, ", ");
}

char autoThirdBase(std::string_view motif, char b) {
    if (motif == "hoogsteen") {
        return b == 'A' ? 'U' : b == 'G' ? 'C' : b;
    }
    return b == 'A' ? 'A' : b == 'G' ? 'G' : b == 'U' ? 'A' : 'G';
}

struct HbondEntry {
    std::string_view rule;
    char a;
    char b;
    std::vector<PairHbond> bonds;
};

const std::vector<HbondEntry>& hbondTable() {
    static const std::vector<HbondEntry> table = [] {
        std::vector<HbondEntry> out = {
            {"wc", 'A', 'U', {{"N1", "N3", 2, "H3"}, {"N6", "O4", 1, "H61"}}},
            {"wc", 'G', 'C', {{"N1", "N3", 1, "H1"}, {"N2", "O2", 1, "H22"}, {"O6", "N4", 2, "H41"}}},
            {"wc", 'G', 'U', {{"O6", "N3", 2, "H3"}, {"N1", "O2", 1, "H1"}}},
            {"reverse-wc", 'A', 'U', {{"N1", "N3", 2, "H3"}, {"N6", "O2", 1, "H61"}}},
            {"reverse-wc", 'G', 'C', {{"N1", "O2", 1, "H1"}, {"N2", "N3", 1, "H22"}}},
            {"hoogsteen", 'A', 'U', {{"N7", "N3", 2, "H3"}, {"N6", "O4", 1, "H62"}}},
            {"hoogsteen", 'G', 'C', {{"N7", "N3", 2, "H3+"}, {"O6", "N4", 2, "H41"}}},
            {"reverse-hoogsteen", 'A', 'U', {{"N7", "N3", 2, "H3"}, {"N6", "O2", 1, "H62"}}},
            {"reverse-hoogsteen", 'A', 'A', {{"N7", "N6", 2, "H61"}, {"N6", "N1", 1, "H62"}}},
            {"reverse-hoogsteen", 'G', 'G', {{"N7", "N1", 2, "H1"}, {"O6", "N2", 2, "H22"}}},
            {"hh", 'A', 'A', {{"N6", "N7", 1, "H62"}, {"N7", "N6", 2, "H62"}}},
        };
        const std::size_t base = out.size();
        for (std::size_t k = 0; k < base; ++k) {  // the same pairs with the strands swapped
            if (out[k].a != out[k].b) {
                HbondEntry swapped{out[k].rule, out[k].b, out[k].a, {}};
                for (const PairHbond& h : out[k].bonds) {
                    swapped.bonds.push_back({h.atom2, h.atom1, 3 - h.donor, h.hydrogen});
                }
                out.push_back(std::move(swapped));
            }
        }
        return out;
    }();
    return table;
}

// _canonical: a canonical pair of `rule`; the Hoogsteen rules need the purine on strand 1
bool canonical(std::string_view rule, char a, char b) {
    return !pairHbonds(rule, a, b).empty() &&
           ((rule != "hoogsteen" && rule != "reverse-hoogsteen") || a == 'A' || a == 'G');
}

std::string reversed(std::string_view text) {
    return std::string(text.rbegin(), text.rend());
}

// _parse_pair / _parse_ids on strings: the parts as text (strip()ped) or an error
std::vector<std::string> splitParts(std::string_view value) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = value.find(',', start);
        parts.push_back(strip(value.substr(start, comma == std::string_view::npos ? std::string_view::npos : comma - start)));
        if (comma == std::string_view::npos) {
            break;
        }
        start = comma + 1;
    }
    return parts;
}

std::string listRepr(const std::vector<std::string>& parts) {
    std::vector<std::string> items;
    for (const std::string& p : parts) {
        items.push_back(repr(p));
    }
    return "[" + join(items, ", ") + "]";
}

template <class Kind>
std::vector<Kind> convertParts(const std::vector<std::string>& parts, const std::string& shown);

template <>
std::vector<std::string> convertParts<std::string>(const std::vector<std::string>& parts, const std::string&) {
    return parts;
}

template <>
std::vector<std::int64_t> convertParts<std::int64_t>(const std::vector<std::string>& parts, const std::string& shown) {
    std::vector<std::int64_t> out;
    for (const std::string& p : parts) {
        std::optional<std::int64_t> v;
        try {
            v = core::tryPythonInt(p);
        } catch (const std::exception&) {
            v.reset();
        }
        if (!v) {
            throw BuildError("cannot read " + shown);
        }
        out.push_back(*v);
    }
    return out;
}

// _parse_pair(parts, kind); `shown` is repr(value)
template <class Kind>
std::vector<Kind> parsePair(std::vector<std::string> parts, const std::string& shown) {
    if (parts.size() == 1) {
        parts.push_back(parts[0]);
    }
    if (parts.size() != 2) {
        throw BuildError("expected two values (strand 1, strand 2), got " + shown);
    }
    return convertParts<Kind>(parts, shown);
}

template <class Kind>
std::vector<Kind> parseIds(std::string_view value, int count) {
    const std::vector<std::string> parts = splitParts(value);
    const std::string shown = repr(value);
    if (count == 2) {
        return parsePair<Kind>(parts, shown);
    }
    if (parts.size() == 3) {
        return convertParts<Kind>(parts, shown);
    }
    if (parts.size() > 3) {
        throw BuildError("expected three values (strand 1, 2, 3), got " + shown);
    }
    // two parts: _parse_pair(list) shows the list; one part: _parse_pair(str) splits it again
    std::vector<Kind> two = parts.size() == 2 ? parsePair<Kind>(parts, listRepr(parts))
                                              : parsePair<Kind>(splitParts(parts[0]), repr(parts[0]));
    if constexpr (std::is_same_v<Kind, std::string>) {
        for (const char c : std::string_view("CDEFGHIJKLMNOPQRSTUVWXYZ")) {
            const std::string id(1, c);
            if (std::find(two.begin(), two.end(), id) == two.end()) {
                two.push_back(id);
                break;
            }
        }
    } else {
        two.push_back(1);
    }
    return two;
}

}  // namespace helix_detail

using namespace helix_detail;

// --- public: templates ---------------------------------------------------------------------------

std::vector<TemplateAtom> helixTemplate(std::string_view base, bool withHydrogens) {
    return templateAtoms(checkedBase(base), withHydrogens).items;
}

std::span<const std::string_view> baseAtomNames(char base) {
    switch (base) {
        case 'A':
            return kBaseA;
        case 'C':
            return kBaseC;
        case 'G':
            return kBaseG;
        case 'U':
            return kBaseU;
        default:
            return {};
    }
}

std::span<const std::string_view> baseHydrogenNames(char base) {
    switch (base) {
        case 'A':
            return kHNamesA;
        case 'C':
            return kHNamesC;
        case 'G':
            return kHNamesG;
        case 'U':
            return kHNamesU;
        default:
            return {};
    }
}

std::span<const std::string_view> sugarHydrogenNames() {
    return kSugarHNames;
}

std::span<const std::string_view> helix_detail::backboneNames() {
    return kBackboneNames;
}

// --- public: sequences ---------------------------------------------------------------------------

ParsedSequence parseSequence(std::string_view text, std::string_view what) {
    const std::string raw = strip(text);
    std::u32string seq;
    if (!raw.empty()) {
        std::u32string s = core::decodeUtf8(upperAscii(raw));
        std::replace(s.begin(), s.end(), U'’', U'\'');
        // re.sub(r"^5'?-?|-?3'?$", "", s)
        std::u32string kept;
        std::size_t p = 0;
        while (p < s.size()) {
            if (p == 0 && s[0] == U'5') {
                std::size_t q = 1;
                if (q < s.size() && s[q] == U'\'') {
                    ++q;
                }
                if (q < s.size() && s[q] == U'-') {
                    ++q;
                }
                p = q;
                continue;
            }
            const std::u32string tail = s.substr(p);
            if (tail == U"-3'" || tail == U"-3" || tail == U"3'" || tail == U"3") {
                break;
            }
            kept += s[p];
            ++p;
        }
        for (const char32_t c : kept) {
            if (c != U'-' && !core::isPythonSpace(c)) {
                seq += c;
            }
        }
    }
    ParsedSequence out;
    if (seq.find(U'T') != std::u32string::npos) {
        std::replace(seq.begin(), seq.end(), U'T', U'U');
        out.notes.push_back(std::string(what) + ": T read as U (RNA)");
    }
    std::set<char32_t> bad;
    for (const char32_t c : seq) {
        if (c != U'A' && c != U'C' && c != U'G' && c != U'U') {
            bad.insert(c);
        }
    }
    if (!bad.empty()) {
        std::vector<std::string> letters;
        for (const char32_t c : bad) {
            letters.push_back(core::encodeUtf8(std::u32string(1, c)));
        }
        throw BuildError(std::string(what) + ": invalid letter(s) " + join(letters, ", ") + " (RNA: A, C, G, U)");
    }
    if (seq.empty()) {
        throw BuildError(std::string(what) + " is empty");
    }
    out.sequence = core::encodeUtf8(seq);
    return out;
}

std::string complement(std::string_view sequence) {
    std::string out;
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        out += complementBase(*it);
    }
    return out;
}

std::string checkPairing(std::optional<std::string_view> pairing) {
    const std::string_view given = pairing && !pairing->empty() ? *pairing : std::string_view("wc");
    const std::string value = lower(strip(given));
    if (std::find(kPairings.begin(), kPairings.end(), value) == kPairings.end()) {
        throw BuildError("unknown pairing " + repr(*pairing) + " (choose from " + join(kPairings, ", ") + ")");
    }
    return value;
}

std::optional<std::string> checkTriplex(std::optional<std::string_view> triplex) {
    const std::string value = lower(strip(triplex ? *triplex : std::string_view()));
    if (value.empty() || value == "none" || value == "no") {
        return std::nullopt;
    }
    if (std::find(kTriplexMotifs.begin(), kTriplexMotifs.end(), value) == kTriplexMotifs.end()) {
        throw BuildError("unknown triplex motif " + repr(*triplex) + " (choose from none, " +
                         join(kTriplexMotifs, ", ") + ")");
    }
    return value;
}

std::string helixKind(std::optional<std::string_view> pairing, std::optional<std::string_view> triplex) {
    const std::string p = checkPairing(pairing);
    const std::optional<std::string> motif = checkTriplex(triplex);
    if (motif) {
        return *motif + " triplex";
    }
    if (p == "wc") {
        return "A-form";
    }
    return std::string(parallel(p) ? "parallel" : "antiparallel") + " " + std::string(pairingLw(p));
}

std::string autoPartner(std::string_view sequence, std::string_view pairing) {
    std::string mapped;
    for (const char b : sequence) {
        mapped += autoPartnerBase(pairing, b);
    }
    return parallel(pairing) ? mapped : reversed(mapped);
}

ResolvedPartner resolvePartner(std::string_view sequence, std::optional<std::string_view> partner,
                               std::optional<std::string_view> pairingName) {
    const std::string pairing = checkPairing(pairingName);
    ParsedSequence parsed = parseSequence(sequence);
    const std::string mode = partner ? strip(*partner) : std::string("auto");
    if (mode.empty()) {
        throw BuildError("type the partner sequence (5'->3'), or choose the Watson-Crick complement or a "
                         "single strand");
    }
    ResolvedPartner out{parsed.sequence, std::nullopt, parsed.notes};
    const std::string m = lower(mode);
    if (m == "none") {
        return out;
    }
    if (m == "auto") {
        if (pairing == "wc") {
            out.partner = complement(out.sequence);
            out.notes.push_back("partner strand: Watson-Crick complement");
            return out;
        }
        out.partner = autoPartner(out.sequence, pairing);
        out.notes.push_back("partner strand: " + std::string(pairingLw(pairing)) + " partners of strand 1 (" +
                            autoPartnerText(pairing) + "), " + (parallel(pairing) ? "parallel" : "antiparallel"));
        return out;
    }
    ParsedSequence other = parseSequence(mode, "partner");
    for (std::string& note : other.notes) {
        out.notes.push_back(std::move(note));
    }
    if (other.sequence.size() != out.sequence.size()) {
        throw BuildError("the partner has " + std::to_string(other.sequence.size()) + " nucleotides but the sequence has " +
                         std::to_string(out.sequence.size()) +
                         "; an explicit partner must have the same length (5'->3')");
    }
    out.partner = other.sequence;
    return out;
}

ResolvedThird resolveThird(std::string_view sequence, std::optional<std::string_view> third, std::string_view motif) {
    const std::string mode = strip(third && !third->empty() ? *third : std::string_view("auto"));
    const bool par = motif == "hoogsteen";
    if (motif != "hoogsteen" && motif != "reverse-hoogsteen") {
        throw std::invalid_argument("resolveThird: unknown motif " + std::string(motif));
    }
    if (lower(mode) == "auto") {
        std::string mapped;
        for (const char b : sequence) {
            mapped += autoThirdBase(motif, b);
        }
        return {par ? mapped : reversed(mapped),
                {std::string("third strand: ") + (par ? "U·A, C+·G" : "A·A, G·G") + " on the strand-1 purines, " +
                 (par ? "parallel" : "antiparallel") + " to strand 1"}};
    }
    ParsedSequence other = parseSequence(mode, "third strand");
    if (other.sequence.size() != sequence.size()) {
        throw BuildError("the third strand has " + std::to_string(other.sequence.size()) +
                         " nucleotides but the sequence has " + std::to_string(sequence.size()) +
                         "; it must have the same length (5'->3')");
    }
    return {other.sequence, other.notes};
}

std::vector<std::string> pairTypes(std::string_view sequence, std::string_view partner) {
    std::vector<std::string> out;
    const std::size_t n = std::min(sequence.size(), partner.size());
    for (std::size_t k = 0; k < n; ++k) {
        const char a = sequence[k];
        const char b = partner[partner.size() - 1 - k];
        const std::string pair{a, b};
        if (pair == "AU" || pair == "UA" || pair == "GC" || pair == "CG") {
            out.emplace_back("wc");
        } else if (pair == "GU" || pair == "UG") {
            out.emplace_back("wobble");
        } else {
            out.emplace_back("mismatch");
        }
    }
    return out;
}

std::vector<PairKind> pairingKinds(std::string_view sequence, std::string_view other,
                                   std::optional<std::string_view> pairingName) {
    const std::string pairing = checkPairing(pairingName);
    const std::string partners = parallel(pairing) ? std::string(other) : reversed(other);
    std::vector<PairKind> out;
    const std::size_t n = std::min(sequence.size(), partners.size());
    if (pairing == "wc") {
        const std::vector<std::string> kinds = pairTypes(sequence, other);
        for (std::size_t k = 0; k < std::min(n, kinds.size()); ++k) {
            out.push_back({sequence[k], partners[k], kinds[k]});
        }
        return out;
    }
    for (std::size_t k = 0; k < n; ++k) {
        out.push_back({sequence[k], partners[k], canonical(pairing, sequence[k], partners[k]) ? "canonical" : "mismatch"});
    }
    return out;
}

std::vector<PairKind> triplexKinds(std::string_view sequence, std::string_view third, std::string_view motif) {
    const std::string partners = motif == "hoogsteen" ? std::string(third) : reversed(third);
    std::vector<PairKind> out;
    for (std::size_t k = 0; k < std::min(sequence.size(), partners.size()); ++k) {
        out.push_back({sequence[k], partners[k], canonical(motif, sequence[k], partners[k]) ? "canonical" : "mismatch"});
    }
    return out;
}

std::span<const PairHbond> pairHbonds(std::string_view rule, char a, char b) {
    for (const HbondEntry& entry : hbondTable()) {
        if (entry.rule == rule && entry.a == a && entry.b == b) {
            return entry.bonds;
        }
    }
    return {};
}

namespace {

std::string pairText(std::size_t k, char a, char b) {
    return std::to_string(k + 1) + ":" + std::string(1, a) + "·" + std::string(1, b);
}

// _pairing_messages
void pairingMessages(std::string_view seq, std::string_view other, std::string_view pairing,
                     std::vector<HelixMessage>& out) {
    const std::vector<PairKind> kinds = pairingKinds(seq, other, pairing);
    const std::string lw(pairingLw(pairing));
    std::size_t good = 0;
    std::vector<std::string> mism, pyr;
    bool gc = false;
    for (std::size_t k = 0; k < kinds.size(); ++k) {
        const PairKind& p = kinds[k];
        if (p.kind != "mismatch") {
            ++good;
            const std::string ab{p.first, p.partner};
            gc = gc || (pairing == "hoogsteen" && (ab == "GC" || ab == "CG"));
        } else {
            mism.push_back(pairText(k, p.first, p.partner));
        }
        if (p.first == 'C' || p.first == 'U') {
            pyr.push_back(pairText(k, p.first, p.partner));
        }
    }
    out.push_back({"info", std::string(pairingLabel(pairing)) + ": " + std::to_string(seq.size()) + " pairs, " +
                               std::to_string(good) + " canonical " + lw + ", " + std::to_string(kinds.size() - good) +
                               " mismatched"});
    if (!mism.empty()) {
        const std::string_view ref = referencePair(pairing);
        out.push_back({"warning", "mismatched pairs (no " + lw + " pair) start from the reference " +
                                      std::string(1, ref[0]) + "·" + std::string(1, ref[1]) +
                                      " pair, are moved off the bases they would hit and are reported: " +
                                      firstEight(mism)});
    }
    if ((pairing == "hoogsteen" || pairing == "reverse-hoogsteen") && !pyr.empty()) {
        out.push_back({"warning", "strand 1 must be the purine strand: a Hoogsteen pair needs the Hoogsteen edge of a "
                                  "strand-1 purine, so these pairs are built as mismatches: " +
                                      firstEight(pyr)});
    }
    if (gc) {
        out.push_back({"info", "G·C+ Hoogsteen pairs need N3-protonated cytosine (C+); written as plain C "
                               "(protonate or restrain for the simulation)"});
    }
}

}  // namespace

std::optional<std::string> rnaRefusal(std::string_view pairing, std::optional<std::string_view> motif) {
    std::string_view what;
    if (motif) {
        if (*motif == "reverse-hoogsteen") {
            what = "the purine-motif (reverse-Hoogsteen) triplex";
        }
    } else if (pairing == "reverse-wc") {
        what = "the parallel reverse Watson-Crick (tWW) duplex";
    } else if (pairing == "reverse-hoogsteen") {
        what = "the reverse-Hoogsteen (tHW) duplex";
    }
    if (what.empty()) {
        return std::nullopt;
    }
    return (motif ? "triplex " + std::string(*motif) : "pairing " + std::string(pairing)) +
           ": RNA is not known to form " + std::string(what) +
           "; DNA and DNA-RNA hybrid builds of it will come with DNA structure building";
}

std::vector<HelixMessage> checkHelix(std::string_view sequence, std::optional<std::string_view> partner,
                                     std::optional<std::string_view> pairingName,
                                     std::optional<std::string_view> triplex, std::optional<std::string_view> third) {
    std::string pairing;
    std::optional<std::string> motif;
    ResolvedPartner resolved;
    ResolvedThird thirdStrand;
    try {
        pairing = checkPairing(pairingName);
        motif = checkTriplex(triplex);
        resolved = resolvePartner(sequence, partner, pairing);
        if (!resolved.partner && (pairing != "wc" || motif)) {
            throw BuildError(!motif ? "a single strand has no pairing: choose a partner strand (or pairing wc)"
                                    : "a triplex needs the Watson-Crick duplex: choose a partner strand");
        }
        if (motif && pairing != "wc") {
            throw BuildError("a triplex is built on the Watson-Crick duplex: choose pairing wc");
        }
        if (const std::optional<std::string> refused = rnaRefusal(pairing, motif)) {  // CH-18
            throw BuildError(*refused);
        }
        if (motif) {
            thirdStrand = resolveThird(resolved.sequence, third, *motif);
        }
    } catch (const BuildError& error) {
        return {{"error", error.what()}};
    }
    const std::string& seq = resolved.sequence;
    std::vector<HelixMessage> out;
    for (const std::string& note : resolved.notes) {
        out.push_back({"info", note});
    }
    if (!resolved.partner) {
        out.push_back({"info", "single strand of " + std::to_string(seq.size()) + " nucleotides in the A-form helix geometry"});
        return out;
    }
    const std::string& other = *resolved.partner;
    const std::string mode = lower(strip(partner && !partner->empty() ? *partner : std::string_view("auto")));
    const bool explicitPartner = mode != "auto" && mode != "none";
    if (pairing != "wc") {
        pairingMessages(seq, other, pairing, out);
        if (explicitPartner && other == autoPartner(seq, pairing)) {
            out.push_back({"info", "the explicit partner is the " + std::string(pairingLw(pairing)) +
                                       " partner strand of 'auto'"});
        }
        return out;
    }
    if (motif) {
        for (const std::string& note : thirdStrand.notes) {
            out.push_back({"info", note});
        }
        const std::vector<PairKind> kinds = triplexKinds(seq, thirdStrand.sequence, *motif);
        std::size_t good = 0;
        std::vector<std::string> mism;
        bool cgc = false;
        for (std::size_t k = 0; k < kinds.size(); ++k) {
            if (kinds[k].kind == "canonical") {
                ++good;
                cgc = cgc || (*motif == "hoogsteen" && kinds[k].first == 'G' && kinds[k].partner == 'C');
            } else {
                mism.push_back(pairText(k, kinds[k].partner, kinds[k].first));
            }
        }
        out.push_back({"info", "triplex (" + std::string(triplexLabel(*motif)) + "): " + std::to_string(good) +
                                   " canonical base triple(s), " + std::to_string(kinds.size() - good) + " mismatched"});
        if (!mism.empty()) {
            out.push_back({"warning", "third-strand mismatches (no Hoogsteen-edge pair; strand-1 pyrimidines have "
                                      "no purine Hoogsteen edge) start from the reference triple, are moved off the "
                                      "bases they would hit and are reported: " +
                                          firstEight(mism)});
        }
        if (cgc) {
            out.push_back({"info", "C+·G-C triples need N3-protonated cytosine (C+) in the third strand; written "
                                   "as plain C (protonate or restrain for the simulation)"});
        }
    }
    const std::vector<std::string> kinds = pairTypes(seq, other);
    const auto nWc = std::count(kinds.begin(), kinds.end(), "wc");
    const auto nWob = std::count(kinds.begin(), kinds.end(), "wobble");
    std::vector<std::string> mism;
    for (std::size_t k = 0; k < kinds.size(); ++k) {
        if (kinds[k] == "mismatch") {
            mism.push_back(pairText(k, seq[k], other[seq.size() - 1 - k]));
        }
    }
    out.push_back({"info", std::to_string(seq.size()) + " base pairs: " + std::to_string(nWc) + " Watson-Crick, " +
                               std::to_string(nWob) + " G·U wobble, " + std::to_string(mism.size()) + " mismatched"});
    if (explicitPartner && other == complement(seq)) {
        out.push_back({"info", "the explicit partner is the Watson-Crick complement"});
    }
    if (!mism.empty()) {
        out.push_back({"warning", "mismatched pairs keep the Watson-Crick frame (bases may collide; minimise "
                                  "or remodel): " +
                                      firstEight(mism)});
    }
    return out;
}

std::vector<std::string> parseChainIds(std::string_view value, int count) {
    return parseIds<std::string>(value, count);
}

std::vector<std::int64_t> parseStartNumbers(std::string_view value, int count) {
    return parseIds<std::int64_t>(value, count);
}

// --- public: Watson-Crick helices ------------------------------------------------------------------

namespace helix_detail {

struct BuiltAtom {
    std::string name;
    Vec3 xyz;
};

// _strand_residue: one nucleotide in the helix frame
std::vector<BuiltAtom> strandResidue(char base, int level, bool partner, bool five, bool three, bool phosphate5,
                                     bool withHydrogens, bool wobble) {
    Atoms t = templateAtoms(base, withHydrogens);
    if (wobble) {
        turnBase(t, base, wobbleTurn(base), withHydrogens);
    }
    std::vector<BuiltAtom> first, last;
    if (five) {
        if (phosphate5) {  // OP3 towards the O3' the previous residue would have
            const Vec3 prevO3 = screw(t.at("O3'"), -1);
            first.push_back({"OP3", t.at("P") + kPOp3 * unit(prevO3 - t.at("P"))});
        } else {
            if (withHydrogens) {
                last.push_back({"HO5'", place(t.at("C4'"), t.at("C5'"), t.at("O5'"), kBondOH, 109.5, 180.0)});
            }
            for (const std::string_view name : {"P", "OP1", "OP2"}) {
                t.erase(name);
            }
        }
    }
    if (three && withHydrogens) {
        last.push_back({"HO3'", place(t.at("C4'"), t.at("C3'"), t.at("O3'"), kBondOH, 109.5, 180.0)});
    }
    std::vector<BuiltAtom> atoms = first;
    for (const std::string_view name : kBackboneNames) {
        if (const Vec3* v = t.find(name)) {
            atoms.push_back({std::string(name), *v});
        }
    }
    for (const std::string_view name : baseAtomNames(base)) {
        if (const Vec3* v = t.find(name)) {
            atoms.push_back({std::string(name), *v});
        }
    }
    if (withHydrogens) {
        for (const std::string_view name : kSugarHNames) {
            atoms.push_back({std::string(name), t.at(name)});
        }
        for (const std::string_view name : baseHydrogenNames(base)) {
            atoms.push_back({std::string(name), t.at(name)});
        }
    }
    for (BuiltAtom& a : last) {
        atoms.push_back(std::move(a));
    }
    for (BuiltAtom& a : atoms) {
        a.xyz = screw(partner ? dyad(a.xyz) : a.xyz, level);
    }
    return atoms;
}

std::string elementOf(std::string_view name) {
    const std::size_t k = name.find_first_not_of("0123456789");
    return k == std::string_view::npos ? std::string() : std::string(name.substr(k, 1));
}

// _mismatch_notes: mismatched pairs whose bases collide (on the written, rounded coordinates)
std::vector<std::string> mismatchNotes(const io::PdbStructure& s, const std::string& seq, const std::string& other,
                                       const std::vector<std::string>& kinds) {
    std::vector<std::string> out;
    const std::size_t n = seq.size();
    for (std::size_t k = 0; k < kinds.size(); ++k) {
        if (kinds[k] != "mismatch") {
            continue;
        }
        const io::PdbResidue& r1 = s.segments[0][k];
        const io::PdbResidue& r2 = s.segments[1][n - 1 - k];
        const auto inBase = [](const io::PdbAtom& a, char base) {
            const auto names = baseAtomNames(base);
            return std::find(names.begin(), names.end(), a.name) != names.end();
        };
        double d = std::numeric_limits<double>::infinity();
        for (const io::PdbAtom& a : r1.atoms) {
            if (!inBase(a, seq[k])) {
                continue;
            }
            for (const io::PdbAtom& b : r2.atoms) {
                if (!inBase(b, other[n - 1 - k])) {
                    continue;
                }
                const std::array<double, 3> pa = a.xyz();
                const std::array<double, 3> pb = b.xyz();
                d = std::min(d, norm(Vec3{pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]}));
            }
        }
        if (d < kClash) {
            out.push_back("WARNING: mismatch " + r1.chain + ":" + std::to_string(r1.resseq) + " " + std::string(1, seq[k]) +
                          " · " + r2.chain + ":" + std::to_string(r2.resseq) + " " + std::string(1, other[n - 1 - k]) +
                          ": the bases collide (closest heavy atoms " + fixed(d, 1) + " Å)");
        }
    }
    return out;
}

// CH-18: the base pairs (the strand-1 base with each partner at its level) whose bases have no N/O
// donor-acceptor pair closer than kHbondNote, as one warning (none: nullopt)
constexpr double kHbondNote = 3.5;

std::optional<std::string> hbondNote(const io::PdbStructure& s, const std::string& seq, const std::string& partner,
                                     const std::string& third, bool partnerParallel, bool thirdParallel) {
    static const std::map<char, std::vector<std::string_view>> donors = {
        {'A', {"N6"}}, {'G', {"N1", "N2"}}, {'C', {"N4"}}, {'U', {"N3"}}};
    static const std::map<char, std::vector<std::string_view>> acceptors = {
        {'A', {"N1", "N3", "N7"}}, {'G', {"O6", "N3", "N7"}}, {'C', {"O2", "N3"}}, {'U', {"O2", "O4"}}};
    const auto coords = [](const io::PdbResidue& r) {
        std::map<std::string, Vec3> out;
        for (const io::PdbAtom& a : r.atoms) {
            const std::array<double, 3> x = a.xyz();
            out.emplace(a.name, Vec3{x[0], x[1], x[2]});
        }
        return out;
    };
    const auto bonded = [&](const std::map<std::string, Vec3>& r1, char b1, const std::map<std::string, Vec3>& r2,
                            char b2) {
        for (int side = 0; side < 2; ++side) {
            const auto& dres = side == 0 ? r1 : r2;
            const auto& ares = side == 0 ? r2 : r1;
            for (const std::string_view x : donors.at(side == 0 ? b1 : b2)) {
                for (const std::string_view y : acceptors.at(side == 0 ? b2 : b1)) {
                    const auto dx = dres.find(std::string(x));
                    const auto ay = ares.find(std::string(y));
                    if (dx != dres.end() && ay != ares.end() && norm(dx->second - ay->second) < kHbondNote) {
                        return true;
                    }
                }
            }
        }
        return false;
    };
    const std::size_t n = seq.size();
    std::vector<std::string> rows;
    for (std::size_t k = 0; k < n; ++k) {
        const io::PdbResidue& r1 = s.segments[0][k];
        const auto c1 = coords(r1);
        for (std::size_t idx = 1; idx < s.segments.size(); ++idx) {
            const bool par = idx == 1 ? partnerParallel : thirdParallel;
            const std::size_t j = par ? k : n - 1 - k;
            const char b2 = (idx == 1 ? partner : third)[j];
            const io::PdbResidue& r2 = s.segments[idx][j];
            if (!bonded(c1, seq[k], coords(r2), b2)) {
                rows.push_back(r1.chain + ":" + std::to_string(r1.resseq) + " " + std::string(1, seq[k]) + " · " +
                               r2.chain + ":" + std::to_string(r2.resseq) + " " + std::string(1, b2));
            }
        }
    }
    if (rows.empty()) {
        return std::nullopt;
    }
    std::string text = "WARNING: " + std::to_string(rows.size()) +
                       " base pair(s) without a base-base H-bond (no N/O donor-acceptor pair within " +
                       fixed(kHbondNote, 1) + " Å): ";
    for (std::size_t k = 0; k < rows.size() && k < 8; ++k) {
        text += (k ? ", " : "") + rows[k];
    }
    if (rows.size() > 8) {
        text += " ...";
    }
    return text;
}

}  // namespace helix_detail

std::string buildHelix(std::string_view sequence, const HelixOptions& options, std::vector<std::string>* notes) {
    const std::string pairing = checkPairing(options.pairing ? std::optional<std::string_view>(*options.pairing)
                                                             : std::nullopt);
    const std::optional<std::string> motif =
        checkTriplex(options.triplex ? std::optional<std::string_view>(*options.triplex) : std::nullopt);
    ResolvedPartner resolved = resolvePartner(
        sequence, options.partner ? std::optional<std::string_view>(*options.partner) : std::nullopt, pairing);
    const std::string chirality = upperAscii(strip(options.chirality.empty() ? std::string("D") : options.chirality));
    if (chirality != "D" && chirality != "L") {
        throw BuildError("chirality must be D or L, got " + repr(chirality));
    }
    const bool duplex = resolved.partner.has_value();
    if (!duplex && pairing != "wc") {
        throw BuildError("a single strand has no pairing: choose a partner strand (or pairing wc)");
    }
    if (motif && (!duplex || pairing != "wc")) {
        throw BuildError("a triplex is built on the Watson-Crick duplex: choose a partner strand and pairing wc");
    }
    if (const std::optional<std::string> refused = rnaRefusal(pairing, motif)) {  // CH-18
        throw BuildError(*refused);
    }
    std::string thirdSeq;
    if (motif) {
        ResolvedThird thirdStrand = resolveThird(
            resolved.sequence, options.third ? std::optional<std::string_view>(*options.third) : std::nullopt, *motif);
        resolved.notes.insert(resolved.notes.end(), thirdStrand.notes.begin(), thirdStrand.notes.end());
        thirdSeq = std::move(thirdStrand.sequence);
    }
    const int count = motif ? 3 : 2;
    const std::vector<std::string> ids = parseChainIds(options.chains, count);
    const std::vector<std::int64_t> firsts = parseStartNumbers(options.start, count);
    for (const std::string& c : ids) {
        if (core::characterCount(c) > 1) {
            throw BuildError("chain ids must be single characters");
        }
    }
    if (duplex && ids[0] == ids[1]) {
        throw BuildError("the two strands need different chain ids");
    }
    if (motif && (ids[2] == ids[0] || ids[2] == ids[1])) {
        throw BuildError("the three strands need different chain ids");
    }
    const bool chosen = options.scheme.has_value() && !options.scheme->empty();
    const std::string scheme = chosen ? *options.scheme : std::string("rosetta");
    if (std::find(kSchemes.begin(), kSchemes.end(), scheme) == kSchemes.end()) {
        throw BuildError("unknown residue name scheme " + repr(scheme) + " (choose from " + join(kSchemes, ", ") + ")");
    }
    const bool rebuilt = pairing != "wc" || motif.has_value();
    const bool experimental = pairing == "hoogsteen" || (motif && *motif == "hoogsteen");  // CH-18

    const std::string& seq = resolved.sequence;
    const std::size_t n = seq.size();
    const std::vector<std::string> kinds =
        duplex && pairing == "wc" ? pairTypes(seq, *resolved.partner) : std::vector<std::string>(n, "none");
    const char hand = chirality[0];

    io::PdbStructure s;
    std::vector<Vec3> coords;
    std::string kind;
    FittedHelix fit;
    const auto addStrand = [&](const std::string& chain, std::int64_t first, const std::string& bases, bool partner) {
        std::vector<io::PdbResidue> seg;
        for (std::size_t j = 0; j < bases.size(); ++j) {
            const std::size_t level = partner ? n - 1 - j : j;
            io::PdbResidue res;
            res.name = std::string(1, bases[j]);
            res.chain = chain;
            res.resseq = first + static_cast<std::int64_t>(j);
            res.chirality = std::string(1, hand);
            for (BuiltAtom& a : strandResidue(bases[j], static_cast<int>(level), partner, j == 0, j == bases.size() - 1,
                                              options.phosphate5, options.hydrogens, kinds[level] == "wobble")) {
                io::PdbAtom atom;
                atom.record = "ATOM";
                atom.name = a.name;
                atom.element = elementOf(a.name);
                res.atoms.push_back(std::move(atom));
                coords.push_back(a.xyz);
            }
            seg.push_back(std::move(res));
        }
        s.segments.push_back(std::move(seg));
    };
    if (!rebuilt) {
        addStrand(ids[0], firsts[0], seq, false);
        if (duplex) {  // strand 2 pairs residue n-1-j of strand 1 at the same level
            addStrand(ids[1], firsts[1], *resolved.partner, true);
        }
    } else {
        kind = pairing != "wc" ? pairing : "triplex-" + *motif;
        std::vector<std::string> sequences = {seq, *resolved.partner};
        if (motif) {
            sequences.push_back(thirdSeq);
        }
        fit = fittedHelix(kind, sequences, options.phosphate5, options.hydrogens, options);
        for (std::size_t strand = 0; strand < fit.strands.size(); ++strand) {
            std::vector<io::PdbResidue> seg;
            for (std::size_t j = 0; j < fit.strands[strand].size(); ++j) {
                FittedResidue& built = fit.strands[strand][j];
                io::PdbResidue res;
                res.name = std::string(1, built.base);
                res.chain = ids[strand];
                res.resseq = firsts[strand] + static_cast<std::int64_t>(j);
                res.chirality = std::string(1, hand);
                for (FittedAtom& a : built.atoms) {
                    io::PdbAtom atom;
                    atom.record = "ATOM";
                    atom.name = a.name;
                    atom.element = elementOf(a.name);
                    res.atoms.push_back(std::move(atom));
                    coords.push_back(a.xyz);
                }
                seg.push_back(std::move(res));
            }
            s.segments.push_back(std::move(seg));
        }
        for (Vec3& v : coords) {  // helix axis along z through the origin
            v = matVec(fit.axisRotation, v) + fit.axisShift;
        }
    }
    if (hand == 'L') {
        for (Vec3& v : coords) {
            v[0] *= -1.0;  // the mirror image: left-handed helix of L-nucleotides
        }
    }
    if (options.center && !coords.empty()) {
        double low = coords[0][2];
        double high = coords[0][2];
        for (const Vec3& v : coords) {
            low = std::min(low, v[2]);
            high = std::max(high, v[2]);
        }
        const double mid = 0.5 * (low + high);
        for (Vec3& v : coords) {
            v = Vec3{v[0] - 0.0, v[1] - 0.0, v[2] - mid};
        }
    }
    std::size_t at = 0;
    for (std::vector<io::PdbResidue>& seg : s.segments) {
        for (io::PdbResidue& res : seg) {
            for (io::PdbAtom& atom : res.atoms) {
                atom.tail = io::newAtomTail(coords[at++]);
            }
        }
    }

    const std::vector<std::string> geometry = rebuilt ? modelNotes(s, fit.levels, kind) : std::vector<std::string>{};
    std::vector<std::string> out = resolved.notes;
    std::vector<std::string> renamed;
    for (std::string& note : convertNames(s, scheme, {}, options.terminals)) {
        if (note.rfind("WARNING: 5'-terminal", 0) != 0) {
            renamed.push_back(std::move(note));
        }
    }
    if (!chosen) {
        renamed = {"Residue names " +
                   std::string(hand == 'D' ? "A, C, G, U" : "0A, 0C, 0G, 0U (the CCD L-ribonucleotides)") +
                   ", ATOM records"};
    }
    out.insert(out.end(), renamed.begin(), renamed.end());
    if (scheme == "pdb") {
        for (std::string& note : setRecords(s, "pdb")) {
            out.push_back(std::move(note));
        }
    }
    if (options.phosphate5 && scheme == "amber") {
        out.emplace_back("WARNING: the Amber 5' units have no 5'-phosphate (tleap cannot type P/OP1/OP2/OP3 "
                         "of A5/LA5 ...); build without the 5'-phosphate for tleap");
    }
    const std::string handText = hand == 'D' ? "right-handed, D" : "left-handed, L: mirror image";
    std::string what = duplex ? "duplex of " + std::to_string(n) + " base pairs"
                              : "single strand of " + std::to_string(n) + " nucleotides";
    if (motif) {
        what = "triplex of " + std::to_string(n) + " base triples";
    }
    // D13: MDBridge names itself in the record
    if (!rebuilt) {
        s.header = {"REMARK 951 MDBridge IDEAL A-FORM RNA (" + handText + "), " + what,
                    "REMARK 951 ARNOTT A-RNA FIBER MONOMER, TWIST " + fixed(kHelixTwist, 1) + " DEG, RISE " +
                        fixed(kHelixRise, 2) + " A, HELIX AXIS Z",
                    "REMARK 951 STRAND " + ids[0] + " " + seq};
    } else {
        const std::string form =
            motif ? "WATSON-CRICK DUPLEX + " + upperAscii(*motif) + " THIRD STRAND"
                  : "PAIRING " + std::string(pairingLw(pairing)) + " (" +
                        (parallel(pairing) ? "PARALLEL" : "ANTIPARALLEL") + " STRANDS)";
        s.header = {"REMARK 951 MDBridge IDEAL RNA (" + handText + "), " + what + ", " + form,
                    std::string(experimental ? "REMARK 951 PDB 3P22 4PLX 6SVS 7JNH, TWIST "
                                             : "REMARK 951 NUCLEOTIDE OF PDB 4JRD, TWIST ") +
                        fixed(fit.twist, 1) + " DEG, RISE " +
                        fixed(fit.rise, 2) + " A, HELIX AXIS Z",
                    "REMARK 951 STRAND " + ids[0] + " " + seq};
    }
    if (duplex) {
        s.header.push_back("REMARK 951 STRAND " + ids[1] + " " + *resolved.partner);
        if (!rebuilt) {
            for (std::string& note : mismatchNotes(s, seq, *resolved.partner, kinds)) {
                out.push_back(std::move(note));
            }
            const auto wobbles = std::count(kinds.begin(), kinds.end(), "wobble");
            if (wobbles) {
                out.push_back(std::to_string(wobbles) +
                              " G·U pair(s) built with an approximate wobble geometry (bases turned in their plane "
                              "about C1'); minimise before production");
            }
        }
    }
    if (motif) {
        s.header.push_back("REMARK 951 STRAND " + ids[2] + " " + thirdSeq + " (THIRD STRAND, " +
                           (*motif == "hoogsteen" ? "PARALLEL" : "ANTIPARALLEL") + " TO STRAND " + ids[0] + ")");
    }
    if (rebuilt) {
        const std::string which =
            motif ? "Triplex, " + std::string(triplexLabel(*motif)) +
                        "; the Watson-Crick duplex and the third strand are fitted together (the Arnott A-form duplex "
                        "leaves no room for a third strand)"
                  : "Pairing " + pairing + ": " + std::string(pairingLabel(pairing));
        const std::string source =
            experimental ? "nucleotides (C3'-endo, anti), base pairs and stacking from the RNA triple helices of PDB "
                           "3P22, 4PLX, 6SVS and 7JNH (mean U·A-U triples), the helix, the pairs and the "
                           "backbone torsions fitted to them and to close the O3'-P links"
                         : "nucleotides (C3'-endo, anti) from PDB 4JRD with backbone torsions fitted to close the "
                           "O3'-P links";
        out.push_back(which + "; helix twist " + fixed(fit.twist, 1) + "°, rise " + fixed(fit.rise, 2) + " Å, " +
                      source + " (reference helix: fit residual " + fixed(fit.cost, 2) + ", closest contacts " +
                      fixed(fit.contact, 2) + " Å heavy atoms, " + fixed(fit.contactH, 2) + " Å H-H)");
    }
    out.insert(out.end(), geometry.begin(), geometry.end());
    std::size_t atomCount = 0;
    for (const auto& seg : s.segments) {
        for (const auto& res : seg) {
            atomCount += res.atoms.size();
        }
    }
    const std::string atomsNote =
        std::to_string(atomCount) + " atoms" + (options.hydrogens ? ", hydrogens included" : ", heavy atoms only");
    if (duplex) {  // CH-18: pairs without a base-base H-bond
        if (std::optional<std::string> note = hbondNote(s, seq, *resolved.partner, thirdSeq, parallel(pairing),
                                                        motif && *motif == "hoogsteen")) {
            out.push_back(std::move(*note));
        }
    }
    if (!rebuilt) {
        out.push_back("Built an ideal A-form RNA " + what + " (" + handText + "): " + atomsNote);
    } else {
        const std::string kindText = motif ? "with a " + *motif + " third strand"
                                           : "with " + std::string(pairingLw(pairing)) + " pairs";
        out.push_back("Built an RNA " + what + " " + kindText + " (" + handText + "): " + atomsNote +
                      "; a starting model: minimise before production");
    }
    // CH-01: a residue number below -999 is refused by io::writePdb (below), before anything is written
    // CH-02: numbers above 9999 wrap (a PDB convention), now with renumber_residues' warning
    std::int64_t highest = firsts[0];
    if (duplex) {
        highest = std::max(highest, firsts[1]);
    }
    if (motif) {
        highest = std::max(highest, firsts[2]);
    }
    if (highest + static_cast<std::int64_t>(n) - 1 > 9999) {
        out.emplace_back("WARNING: residue numbers above 9999 wrap in the PDB columns");
    }
    std::string text = io::writePdb(s);
    if (notes) {
        notes->insert(notes->end(), out.begin(), out.end());
    }
    return text;
}

}  // namespace mdbridge::chem
