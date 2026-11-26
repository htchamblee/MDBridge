// Counterparts of MDATK v0.0.4.1's tests/test_pdbtools.py (32 tests; MDATK's own code, ported),
// one TEST per MDATK test, named after it. The fixtures (AmberTools' RNA.lib nucleotides, make_pdb,
// atom_line, numbered) are ported from that file. Where MDBridge deliberately differs (CH items),
// the expectation is MDBridge's. test_embedded_copy_runs_standalone belongs to milestone 7 (the
// generated scripts' embedded copy); test_mdatk_pdb_cli_check_clean_and_log is checked here through
// the library calls the commands make (the command line and its step log: W2-C).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/value.hpp"
#include "test_support.hpp"

using namespace mdbridge;
namespace pt = mdbridge::chem;

namespace {

struct AtomXyz {
    const char* name;
    double x, y, z;
};

const std::vector<AtomXyz> kBackbone = {
    {"P", 2.314, 1.777, -0.478},  {"OP1", 1.933, 3.073, -1.082}, {"OP2", 3.157, 0.872, -1.292},
    {"O5'", 2.994, 2.142, 0.923}, {"C5'", 2.156, 2.552, 2.020},  {"H5''", 1.605, 3.450, 1.741},
    {"C4'", 3.008, 2.848, 3.243}, {"O4'", 3.228, 1.605, 3.977},  {"C1'", 4.559, 1.151, 3.778},
    {"C3'", 4.416, 3.378, 2.961}, {"C2'", 5.230, 2.096, 2.779},  {"O2'", 6.273, 2.050, 3.757},
    {"HO2'", 6.759, 1.235, 3.611}, {"O3'", 5.030, 4.129, 3.998},
};

const std::map<char, std::vector<AtomXyz>>& bases() {
    static const std::map<char, std::vector<AtomXyz>> b = {
        {'A', {{"N9", 4.539, -0.297, 3.316}, {"C8", 4.337, -0.846, 2.077}, {"N7", 4.453, -2.141, 2.049},
               {"C5", 4.754, -2.480, 3.363}, {"C6", 4.998, -3.708, 3.989}, {"N6", 4.977, -4.883, 3.346},
               {"N1", 5.265, -3.689, 5.302}, {"C2", 5.285, -2.520, 5.935}, {"N3", 5.073, -1.309, 5.456},
               {"C4", 4.806, -1.356, 4.133}}},
        {'C', {{"N1", 4.501, -0.268, 3.327}, {"C6", 4.258, -0.582, 2.026}, {"C5", 4.206, -1.878, 1.616},
               {"C4", 4.416, -2.881, 2.614}, {"N4", 4.379, -4.159, 2.286}, {"N3", 4.653, -2.564, 3.884},
               {"C2", 4.701, -1.262, 4.273}, {"O2", 4.919, -0.928, 5.447}}},
        {'G', {{"N9", 4.500, -0.268, 3.327}, {"C8", 4.268, -0.773, 2.063}, {"N7", 4.279, -2.081, 2.002},
               {"C5", 4.535, -2.475, 3.310}, {"C6", 4.665, -3.772, 3.873}, {"O6", 4.578, -4.865, 3.316},
               {"N1", 4.926, -3.709, 5.247}, {"C2", 5.048, -2.547, 5.982}, {"N2", 5.299, -2.702, 7.289},
               {"N3", 4.928, -1.333, 5.453}, {"C4", 4.674, -1.377, 4.118}}},
        {'U', {{"N1", 4.539, -0.307, 3.313}, {"C6", 4.317, -0.673, 2.011}, {"C5", 4.360, -1.958, 1.634},
               {"C4", 4.643, -2.997, 2.590}, {"O4", 4.703, -4.199, 2.337}, {"N3", 4.853, -2.538, 3.874},
               {"C2", 4.815, -1.224, 4.292}, {"O2", 5.015, -0.915, 5.455}}},
    };
    return b;
}

std::vector<AtomXyz> residueAtoms(char base, bool mirror = false, double dz = 0.0) {
    const double sign = mirror ? -1.0 : 1.0;
    std::vector<AtomXyz> out;
    const auto add = [&](const AtomXyz& a) { out.push_back({a.name, sign * a.x, a.y, a.z + dz}); };
    for (std::size_t k = 0; k < 9; ++k) {
        add(kBackbone[k]);
    }
    for (const AtomXyz& a : bases().at(base)) {
        add(a);
    }
    for (std::size_t k = 9; k < kBackbone.size(); ++k) {
        add(kBackbone[k]);
    }
    return out;
}

std::string fmt(const char* format, auto... args) {
    char buf[256];
    std::snprintf(buf, sizeof buf, format, args...);
    return buf;
}

std::string atomLine(int serial, const std::string& name, const std::string& resname, const std::string& chain,
                     int resseq, double x, double y, double z, const std::string& record = "ATOM",
                     std::string element = "") {
    if (element.empty()) {
        element = name.substr(name.find_first_not_of("0123456789"), 1);
    }
    const std::string field = name.size() == 4 ? name : " " + name + std::string(3 - std::min<std::size_t>(3, name.size()), ' ');
    return fmt("%-6s%5d %s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s", record.c_str(), serial,
               field.c_str(), resname.c_str(), chain.empty() ? " " : chain.c_str(), resseq, x, y, z, element.c_str());
}

struct Res {
    std::string name;
    int resseq;
    char base;
    bool mirror;
};
using Segments = std::vector<std::pair<std::string, std::vector<Res>>>;

std::string makePdb(const Segments& segments, const std::map<std::string, std::string>& records = {},
                    const std::vector<std::string>& extra = {}, const std::vector<std::string>& header = {},
                    int firstSerial = 1) {
    std::vector<std::string> lines = header;
    int serial = firstSerial - 1;
    int k = 0;
    for (const auto& [chain, residues] : segments) {
        const Res* last = nullptr;
        for (const Res& r : residues) {
            for (const AtomXyz& a : residueAtoms(r.base, r.mirror, 7.0 * k)) {
                ++serial;
                const auto it = records.find(r.name);
                lines.push_back(atomLine(serial, a.name, r.name, chain, r.resseq, a.x, a.y, a.z,
                                         it == records.end() ? "ATOM" : it->second));
            }
            ++k;
            last = &r;
        }
        ++serial;
        lines.push_back(fmt("TER   %5d      %3s %1s%4d", serial, last->name.c_str(), chain.empty() ? " " : chain.c_str(),
                            last->resseq));
    }
    lines.insert(lines.end(), extra.begin(), extra.end());
    lines.push_back("END");
    std::string text;
    for (const std::string& l : lines) {
        text += l + "\n";
    }
    return text;
}

pt::AtomCoords coords(char base, bool mirror = false) {
    pt::AtomCoords out;
    for (const AtomXyz& a : residueAtoms(base, mirror)) {
        out.emplace_back(a.name, pt::Vec3{a.x, a.y, a.z});
    }
    return out;
}

using Names = std::vector<std::tuple<std::string, std::int64_t, std::string>>;
Names names(const io::PdbStructure& s) {
    Names out;
    for (const io::PdbResidue* r : s.residues()) {
        out.emplace_back(r->chain, r->resseq, r->name);
    }
    return out;
}
std::vector<std::string> resnames(const io::PdbStructure& s) {
    std::vector<std::string> out;
    for (const io::PdbResidue* r : s.residues()) {
        out.push_back(r->name);
    }
    return out;
}
std::vector<std::string> chainsOf(const io::PdbStructure& s) {
    std::vector<std::string> out;
    for (const io::PdbResidue* r : s.residues()) {
        out.push_back(r->chain);
    }
    return out;
}

bool anyNote(const pt::Notes& notes, const std::function<bool(const std::string&)>& pred) {
    return std::any_of(notes.begin(), notes.end(), pred);
}
bool hasNote(const pt::Notes& notes, std::string_view part) {
    return anyNote(notes, [&](const std::string& n) { return n.find(part) != std::string::npos; });
}
bool hasNote(const pt::Notes& notes, std::string_view part1, std::string_view part2) {
    return anyNote(notes, [&](const std::string& n) {
        return n.find(part1) != std::string::npos && n.find(part2) != std::string::npos;
    });
}
bool startsNote(const pt::Notes& notes, std::string_view prefix) {
    return anyNote(notes, [&](const std::string& n) { return n.starts_with(prefix); });
}

pt::CleanResult clean(const std::string& text, pt::CleanOptions o = {}) {
    return pt::clean(pt::PdbSource::text(text), o);
}
pt::CleanOptions scheme(const char* s) {
    pt::CleanOptions o;
    o.scheme = s;
    return o;
}

std::vector<std::string> lines(const std::string& text) {
    std::vector<std::string> out;
    std::size_t first = 0;
    while (first < text.size()) {
        const std::size_t end = text.find('\n', first);
        out.push_back(text.substr(first, end - first));
        first = end == std::string::npos ? text.size() : end + 1;
    }
    return out;
}
std::string joinLines(const std::vector<std::string>& ls) {
    std::string out;
    for (const std::string& l : ls) {
        out += l + "\n";
    }
    return out;
}
std::string trimmed(const std::string& s) {
    const std::size_t a = s.find_first_not_of(' ');
    const std::size_t b = s.find_last_not_of(' ');
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}
std::size_t count(const std::string& text, std::string_view part) {
    std::size_t n = 0;
    for (std::size_t at = text.find(part); at != std::string::npos; at = text.find(part, at + 1)) {
        ++n;
    }
    return n;
}

const Segments kHybrid = {
    {"A", {{"G", 1, 'G', false}, {"C", 2, 'C', false}, {"A", 3, 'A', false}, {"U", 4, 'U', false}}},
    {"B", {{"0G", 5, 'G', true}, {"0C", 6, 'C', true}, {"0A", 7, 'A', true}, {"0U", 8, 'U', true}}},
    {"C", {{"0A", 9, 'A', true}}},
    {"D", {{"U", 10, 'U', false}}},
    {"E", {{"DG", 11, 'G', false}, {"PSU", 12, 'U', false}}}};
const Segments kCpptraj = {{"", {{"G5", 1, 'G', false}, {"C", 2, 'C', false}, {"U3", 3, 'U', false}}},
                           {"", {{"G5", 4, 'G', true}, {"A", 5, 'A', true}, {"C3", 6, 'C', true}}}};
const Segments kReference = {{"A", {{"G", 1, 'G', false}, {"C", 2, 'C', false}, {"U", 3, 'U', false}}},
                             {"B", {{"LG", 1, 'G', true}, {"LA", 2, 'A', true}, {"LC", 3, 'C', true}}}};
const std::vector<std::string>& ions() {
    static const std::vector<std::string> v = {atomLine(500, "Na+", "Na+", "", 7, 1, 1, 1, "ATOM", "NA"), "TER",
                                               atomLine(501, "O", "WAT", "", 8, 2, 2, 2),
                                               atomLine(502, "H1", "WAT", "", 8, 2, 2, 2), "TER"};
    return v;
}

template <class Fn>
std::string errorOf(Fn&& fn) {
    try {
        fn();
    } catch (const io::PdbError& e) {
        return e.what();
    }
    return "<no error>";
}

double volume(const pt::SugarChirality& s, const char* centre) {
    for (const auto& [k, v] : s.volumes) {
        if (k == centre) {
            return v;
        }
    }
    return std::nan("");
}

std::string chimeraxDamage(const std::string& line) {
    const std::string fused = trimmed(line.substr(17, 3)) + line[21] + ":";
    return line.substr(0, 22 - fused.size()) + fused + line.substr(22);
}

// tests/test_pdbtools.py numbered: G residues with (resseq, icode) per chain
std::string numbered(const std::vector<std::pair<std::string, std::vector<std::pair<int, std::string>>>>& spec) {
    std::vector<std::string> ls;
    int serial = 0;
    for (const auto& [chain, residues] : spec) {
        for (const auto& [resseq, icode] : residues) {
            for (const AtomXyz& a : residueAtoms('G')) {
                ++serial;
                std::string l = atomLine(serial, a.name, "G", chain, resseq, a.x, a.y, a.z);
                l[26] = icode.empty() ? ' ' : icode[0];
                ls.push_back(l);
            }
        }
        ls.push_back("TER");
    }
    return joinLines(ls) + "END\n";
}
const std::vector<std::pair<std::string, std::vector<std::pair<int, std::string>>>> kGapped = {
    {"A", {{5, ""}, {6, ""}, {10, ""}, {10, "A"}, {11, ""}}}, {"B", {{20, ""}, {21, ""}}}};

}  // namespace

// --- chirality ---------------------------------------------------------------------------------------

TEST(PdbToolsMdatk, test_chiral_volume_convention) {
    EXPECT_DOUBLE_EQ(pt::chiralVolume({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}), 1.0);
    EXPECT_DOUBLE_EQ(pt::chiralVolume({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}), -1.0);
}

TEST(PdbToolsMdatk, test_sugar_chirality_d_and_mirrored) {
    for (const char base : std::string("ACGU")) {
        const pt::SugarChirality d = pt::sugarChirality(coords(base), base);
        EXPECT_EQ(d.call, 'D') << base;
        ASSERT_EQ(d.volumes.size(), 3u) << base;
        EXPECT_LT(volume(d, "C4'"), 0) << base;
        EXPECT_GT(volume(d, "C3'"), 0) << base;
        EXPECT_LT(volume(d, "C1'"), 0) << base;
        const pt::SugarChirality l = pt::sugarChirality(coords(base, true), base);
        EXPECT_EQ(l.call, 'L') << base;
        for (const char* c : {"C4'", "C3'", "C1'"}) {
            EXPECT_NEAR(volume(l, c), -volume(d, c), 1e-9) << base << c;
        }
    }
}

TEST(PdbToolsMdatk, test_sugar_chirality_missing_atoms_and_tie) {
    const auto without = [](pt::AtomCoords c, std::initializer_list<std::string_view> drop) {
        std::erase_if(c, [&](const auto& e) { return std::find(drop.begin(), drop.end(), e.first) != drop.end(); });
        return c;
    };
    const pt::SugarChirality partial = pt::sugarChirality(without(coords('G'), {"C5'", "O3'"}), 'G');
    EXPECT_EQ(partial.call, 'D');
    ASSERT_EQ(partial.volumes.size(), 1u);
    EXPECT_EQ(partial.volumes[0].first, "C1'");
    const pt::AtomCoords noN9 = without(coords('G', true), {"N9"});
    const pt::SugarChirality l = pt::sugarChirality(noN9, 'G');
    EXPECT_EQ(l.call, 'L');  // purine: N1 is not used in place of N9
    ASSERT_EQ(l.volumes.size(), 2u);
    EXPECT_EQ(l.volumes[0].first, "C4'");
    EXPECT_EQ(l.volumes[1].first, "C3'");
    const pt::SugarChirality none = pt::sugarChirality({{"P", {0, 0, 0}}}, 'A');
    EXPECT_FALSE(none.call);
    EXPECT_TRUE(none.volumes.empty());
    const pt::AtomCoords tie = {{"C4'", {0, 0, 0}}, {"O4'", {1, 0, 0}}, {"C3'", {0, 1, 0}}, {"C5'", {0, 0, 1}},
                                {"C1'", {2, 0, 0}}, {"N9", {2, 1, 0}},  {"C2'", {2, 0, 1}}};
    const pt::SugarChirality t = pt::sugarChirality(tie, 'A');
    EXPECT_FALSE(t.call);
    EXPECT_GT(volume(t, "C4'"), 0);
    EXPECT_LT(volume(t, "C1'"), 0);
}

TEST(PdbToolsMdatk, test_assign_chirality_detection_override_and_fallback) {
    const std::string text = makePdb({{"A", {{"G", 1, 'G', false}, {"C", 2, 'C', false}}},
                                      {"B", {{"G", 1, 'G', true}, {"U", 2, 'U', true}}},
                                      {"C", {{"0A", 1, 'A', false}}}});
    const auto hands = [](const io::PdbStructure& s) {
        std::string out;
        for (const io::PdbResidue* r : s.residues()) {
            out += r->chirality.value_or("?");
        }
        return out;
    };
    io::PdbStructure s = io::readPdbText(text);
    pt::Notes notes = pt::assignChirality(s);
    EXPECT_EQ(hands(s), "DDLLD");
    EXPECT_TRUE(hasNote(notes, "L sugars in 2 residues with D names", "B:1-2"));
    EXPECT_TRUE(anyNote(notes, [](const std::string& n) {
        return n.starts_with("WARNING: L-named residues with D sugar") && n.find("C:1 0A") != std::string::npos;
    }));
    s = io::readPdbText(text);
    notes = pt::assignChirality(s, true, "A:2, C");
    EXPECT_EQ(hands(s), "DLLLL");
    EXPECT_TRUE(hasNote(notes, "forced L although the sugar geometry is D", "A:2 C"));
    s = io::readPdbText(text);
    pt::assignChirality(s, false);
    EXPECT_EQ(hands(s), "DDDDL");  // names only

    std::vector<std::string> kept;  // no sugar ring: undetermined -> residue name
    for (const std::string& l : lines(makePdb({{"A", {{"LG", 1, 'G', false}}}}))) {
        const std::string atom = l.size() > 16 ? trimmed(l.substr(12, 4)) : "";
        if (!l.starts_with("ATOM") || atom == "P" || atom == "OP1" || atom == "OP2" || atom == "O5'" || atom == "C5'") {
            kept.push_back(l);
        }
    }
    s = io::readPdbText(joinLines(kept));
    notes = pt::assignChirality(s);
    EXPECT_EQ(s.segments[0][0].chirality, "L");
    EXPECT_TRUE(anyNote(notes, [](const std::string& n) {
        return n.starts_with("Chirality undetermined") && n.find("A:1 LG") != std::string::npos;
    }));
}

TEST(PdbToolsMdatk, test_parse_ranges) {
    using R = pt::ResidueRange;
    EXPECT_EQ(pt::parseRanges("C:45-73, A:1-5"), (std::vector<R>{{"C", 45, 73}, {"A", 1, 5}}));
    EXPECT_EQ(pt::parseRanges(std::vector<pt::RangeSpec::Item>{std::string("C:45"), std::string("B"), std::string(":3-4")}),
              (std::vector<R>{{"C", 45, 45}, {"B", std::nullopt, std::nullopt}, {"", 3, 4}}));
    EXPECT_EQ(pt::parseRanges(std::vector<pt::RangeSpec::Item>{R{"A", 2, 9}, R{"B", 4, std::nullopt}}),
              (std::vector<R>{{"A", 2, 9}, {"B", 4, 4}}));
    EXPECT_EQ(pt::parseRanges("A:-3--1"), (std::vector<R>{{"A", -3, -1}}));
    EXPECT_TRUE(pt::parseRanges(pt::RangeSpec()).empty());
    EXPECT_THROW(pt::parseRanges("45-73"), io::PdbError);
    // CH-15: a reversed range is an error, not a range that matches nothing
    EXPECT_EQ(errorOf([] { pt::parseRanges("A:4-1"); }), "residue range A:4-1 ends before it starts");
}

// --- residue names -----------------------------------------------------------------------------------

TEST(PdbToolsMdatk, test_parse_residue_name) {
    struct Case {
        const char* name;
        std::optional<io::ResidueName> expected;
    };
    const auto n = [](char b, const char* sugar, char hand, const char* term) {
        return std::optional<io::ResidueName>(io::ResidueName{b, sugar, hand, term});
    };
    const std::vector<Case> cases = {
        {"A", n('A', "rna", 'D', "")},   {"G5", n('G', "rna", 'D', "5")},  {"U3", n('U', "rna", 'D', "3")},
        {"CN", n('C', "rna", 'D', "N")}, {"RA5", n('A', "rna", 'D', "5")}, {"rG", n('G', "rna", 'D', "")},
        {"RAN", n('A', "rna", 'D', "N")}, {"ADE", n('A', "rna", 'D', "")}, {"GUA", n('G', "rna", 'D', "")},
        {"CYT", n('C', "rna", 'D', "")}, {"URA", n('U', "rna", 'D', "")},  {"URI", n('U', "rna", 'D', "")},
        {"LA", n('A', "rna", 'L', "")},  {"LG3", n('G', "rna", 'L', "3")}, {"0G5", n('G', "rna", 'L', "5")},
        {"0C", n('C', "rna", 'L', "")},  {"DA", n('A', "dna", 'D', "")},   {"DT3", n('T', "dna", 'D', "3")},
        {"THY", n('T', "dna", 'D', "")}, {"HOH", std::nullopt},           {"MG", std::nullopt},
        {"PSU", std::nullopt},           {"LT", std::nullopt},            {"RT", std::nullopt},
    };
    for (const Case& c : cases) {
        EXPECT_EQ(io::parseResidueName(c.name), c.expected) << c.name;
    }
}

TEST(PdbToolsMdatk, test_scheme_name) {
    struct Case {
        char chir;
        const char* term;
        const char* amber;
        const char* ccd;
        const char* rosetta;
    };
    for (const Case& c : std::vector<Case>{{'D', "", "G", "G", "G"},       {'D', "5", "G5", "G", "G"},
                                           {'D', "3", "G3", "G", "G"},     {'D', "N", "GN", "G", "G"},
                                           {'L', "", "LG", "0G", "0G"},    {'L', "5", "LG5", "0G5", "0G"},
                                           {'L', "3", "LG3", "0G3", "0G"}, {'L', "N", "LG5", "0G", "0G"}}) {
        EXPECT_EQ(pt::schemeName('G', c.chir, c.term, "amber"), c.amber);
        EXPECT_EQ(pt::schemeName('G', c.chir, c.term, "ccd"), c.ccd);
        EXPECT_EQ(pt::schemeName('G', c.chir, c.term, "rosetta"), c.rosetta);
        EXPECT_EQ(pt::schemeName('G', c.chir, c.term, "pdb"), c.rosetta);
    }
}

TEST(PdbToolsMdatk, test_scheme_name_unknown) {
    EXPECT_THROW(pt::schemeName('G', 'D', "", "charmm"), io::PdbError);
}

TEST(PdbToolsMdatk, test_convert_names_all_schemes_and_round_trip) {
    const std::string text = makePdb(kHybrid);
    const pt::CleanResult amber = clean(text, scheme("amber"));
    EXPECT_EQ(resnames(amber.structure), (std::vector<std::string>{"G5", "C", "A", "U3", "LG5", "LC", "LA", "LU3",
                                                                    "LA5", "UN", "DG", "PSU"}));
    EXPECT_TRUE(hasNote(amber.notes, "no unit for a single-residue L segment", "C:9 LA5"));
    EXPECT_TRUE(hasNote(amber.notes, "5'-terminal residues carry P/OP1/OP2"));
    const pt::CleanResult ccd = clean(amber.text, scheme("ccd"));
    EXPECT_EQ(resnames(ccd.structure), (std::vector<std::string>{"G", "C", "A", "U", "0G", "0C", "0A", "0U3", "0A",
                                                                  "U", "DG", "PSU"}));
    EXPECT_TRUE(hasNote(ccd.notes, "keep the CCD component names", "B:5"));
    EXPECT_FALSE(hasNote(ccd.notes, "missing TER?"));
    std::vector<std::string> noP;  // a 5'-OH end gets the variant
    for (const std::string& l : lines(text)) {
        const bool phosphate0G = l.size() > 20 && trimmed(l.substr(17, 3)) == "0G" &&
                                 std::set<std::string>{"P", "OP1", "OP2", "OP3"}.count(trimmed(l.substr(12, 4)));
        if (!phosphate0G) {
            noP.push_back(l);
        }
    }
    EXPECT_EQ(resnames(clean(joinLines(noP), scheme("ccd")).structure)[4], "0G5");
    const pt::CleanResult rosetta = clean(ccd.text, scheme("rosetta"));
    EXPECT_EQ(resnames(rosetta.structure), (std::vector<std::string>{"G", "C", "A", "U", "0G", "0C", "0A", "0U", "0A",
                                                                      "U", "DG", "PSU"}));
    EXPECT_EQ(clean(rosetta.text, scheme("amber")).text, amber.text);
    EXPECT_EQ(clean(amber.text, scheme("pdb")).text, clean(rosetta.text, scheme("pdb")).text);
    std::vector<std::string> a;
    std::vector<std::string> r;
    for (const std::string& l : lines(amber.text)) {
        if (l.starts_with("ATOM")) {
            a.push_back(l.substr(0, 17) + l.substr(20));
        }
    }
    for (const std::string& l : lines(rosetta.text)) {
        if (l.starts_with("ATOM")) {
            r.push_back(l.substr(0, 17) + l.substr(20));
        }
    }
    EXPECT_EQ(a, r);
}

TEST(PdbToolsMdatk, test_convert_names_mirrored_plain_names_and_terminal_mismatch) {
    const std::string text = makePdb({{"", {{"G5", 1, 'G', true}, {"C", 2, 'C', true}, {"C3", 3, 'C', true}}},
                                      {"", {{"A5", 4, 'A', false}, {"U3", 5, 'U', false}, {"G", 6, 'G', false}}}});
    const pt::CleanResult r = clean(text, scheme("amber"));
    EXPECT_EQ(resnames(r.structure), (std::vector<std::string>{"LG5", "LC", "LC3", "A5", "U", "G3"}));
    EXPECT_TRUE(anyNote(r.notes, [](const std::string& n) {
        return n.starts_with("Terminal names that are not at a segment end") && n.find(":5 U3") != std::string::npos;
    }));
    pt::CleanOptions o = scheme("amber");
    o.detectL = false;
    EXPECT_EQ(clean(text, o).structure.segments[0][0].name, "G5");
}

TEST(PdbToolsMdatk, test_to_lrna_renames_only_ranges_and_keeps_terminals) {
    const std::string text =
        makePdb({{"A", {{"G5", 1, 'G', false}, {"C", 2, 'C', false}, {"A", 3, 'A', false}, {"U3", 4, 'U', false}}},
                 {"B", {{"G5", 1, 'G', false}, {"C3", 2, 'C', false}}}});
    pt::CleanResult r = pt::toLrna(pt::PdbSource::text(text), "A:3-4, B");
    EXPECT_EQ(names(r.structure), (Names{{"A", 1, "G5"}, {"A", 2, "C"}, {"A", 3, "LA"}, {"A", 4, "LU3"},
                                         {"B", 1, "LG5"}, {"B", 2, "LC3"}}));
    EXPECT_TRUE(hasNote(r.notes, "forced L although the sugar geometry is D"));
    r = pt::toLrna(pt::PdbSource::text(text), std::vector<pt::RangeSpec::Item>{pt::ResidueRange{"B", 1, 2}}, "ccd");
    const std::vector<std::string> n = resnames(r.structure);
    EXPECT_EQ(std::vector<std::string>(n.end() - 2, n.end()), (std::vector<std::string>{"0G", "0C3"}));
    EXPECT_EQ(std::vector<std::string>(n.begin(), n.begin() + 4), (std::vector<std::string>{"G5", "C", "A", "U3"}));
}

// --- parsing / writing ---------------------------------------------------------------------------------

TEST(PdbToolsMdatk, test_chimerax_repair) {
    const Segments segs = {{"C", {{"G", 1, 'G', true}, {"0A", 2, 'A', true}, {"LG5", 3, 'G', true}}},
                           {"B", {{"C", 4, 'C', false}, {"U", 5, 'U', false}}}};
    const std::string good = makePdb(segs);
    std::vector<std::string> damagedLines;
    for (const std::string& l : lines(good)) {
        damagedLines.push_back(l.starts_with("ATOM") || l.starts_with("TER") ? chimeraxDamage(l) : l);
    }
    const std::string damaged = joinLines(damagedLines);
    EXPECT_NE(damaged.find("  GC:   1"), std::string::npos);
    EXPECT_NE(damaged.find("LG5C:   3"), std::string::npos);
    const pt::CleanResult fixed = clean(damaged);
    EXPECT_EQ(fixed.text, clean(good).text);
    EXPECT_EQ(names(fixed.structure), (Names{{"C", 1, "G"}, {"C", 2, "0A"}, {"C", 3, "LG5"}, {"B", 4, "C"}, {"B", 5, "U"}}));
    const std::size_t nAtoms = count(good, "\nATOM") + (good.starts_with("ATOM") ? 1 : 0);
    EXPECT_TRUE(startsNote(fixed.notes, "Repaired " + std::to_string(nAtoms) + " ChimeraX-damaged"));
    const std::string line = lines(good)[0];
    const io::PdbResidue shifted = io::readPdbText(line.substr(0, 17) + "  LG5C:" + line.substr(22)).segments[0][0];
    EXPECT_EQ(std::tie(shifted.name, shifted.chain, shifted.resseq), std::make_tuple("LG5", "C", std::int64_t{1}));
    EXPECT_NEAR(shifted.atoms[0].xyz()[0], -2.314, 1e-9);
    const io::PdbResidue ragged =
        io::readPdbText("ATOM      1  P     0GC:  7   -2.314   1.777   -0.478 1.00 0.00 P").segments[0][0];
    EXPECT_EQ(std::tie(ragged.name, ragged.chain, ragged.resseq, ragged.atoms[0].element),
              std::make_tuple("0G", "C", std::int64_t{7}, "P"));
    EXPECT_NEAR(ragged.atoms[0].xyz()[2], -0.478, 1e-9);
    EXPECT_EQ(io::readPdbText(line.substr(0, 17) + "G   C" + line.substr(22)).segments[0][0].name, "G");
}

TEST(PdbToolsMdatk, test_writer_columns_alignment_and_elements) {
    const std::string raw = joinLines({
        "ATOM      7 C1'    G A   1       1.000   2.000   3.000  1.00  0.00",
        "ATOM      8 H5''   G A   1       1.000   2.000   3.000  1.00  0.00",
        "HETATM    9 MG    MG A   2       4.000   5.000   6.000  1.00  0.00",
        "HETATM   10 Na+  Na+ A   3       4.000   5.000   6.000  1.00  0.00",
        "ATOM     11  CA  ALA B   1       0.000   0.000   0.000  1.00  0.00           C",
    });
    const std::vector<std::string> out = lines(clean(raw).text);
    EXPECT_EQ(out[0], "ATOM      1  C1'   G A   1       1.000   2.000   3.000  1.00  0.00           C");
    EXPECT_EQ(out[1].substr(12, 4), "H5''");
    EXPECT_EQ(out[1].substr(76, 2), " H");
    EXPECT_EQ(out[2].substr(12, 4), "MG  ");
    EXPECT_EQ(out[2].substr(76, 2), "MG");
    EXPECT_EQ(out[2].substr(17, 3), " MG");
    EXPECT_EQ(out[3].substr(12, 4), "Na+ ");
    EXPECT_EQ(out[3].substr(76, 2), "NA");
    std::vector<std::string> records;
    for (const std::string& l : out) {
        records.push_back(l.substr(0, 6));
    }
    EXPECT_EQ(records, (std::vector<std::string>{"ATOM  ", "ATOM  ", "HETATM", "HETATM", "TER   ", "ATOM  ", "TER   ", "END"}));
    EXPECT_EQ(out[4], "TER       5      Na+ A   3");
    EXPECT_EQ(out[5].substr(12, 4), " CA ");
    EXPECT_EQ(out[5].substr(6, 5), "    6");
}

TEST(PdbToolsMdatk, test_models_header_trailer_and_anisou) {
    std::vector<std::string> body = lines(makePdb({{"A", {{"G", 1, 'G', false}}}}));
    body.pop_back();  // without END
    std::vector<std::string> ls = {"HEADER    RNA", "REMARK   1 kept", "SEQRES   1 A    1    G", "CRYST1   10.000",
                                   "MODEL        1", body[0],
                                   "ANISOU    1  P     G A   1     1000   2000   3000    100    200    300       P"};
    ls.insert(ls.end(), body.begin() + 1, body.end());
    ls.push_back("ENDMDL");
    ls.push_back("MODEL        2");
    for (std::string l : body) {
        const std::size_t at = l.find("  1.00  0.00");
        if (at != std::string::npos) {
            l.replace(at, 12, "  1.00  9.99");
        }
        ls.push_back(l);
    }
    ls.insert(ls.end(), {"ENDMDL", "MASTER        0", "END"});
    std::string text;
    for (std::size_t k = 0; k < ls.size(); ++k) {
        text += (k ? "\n" : "") + ls[k];
    }
    const pt::CleanResult r = clean(text);
    const std::vector<std::string> out = lines(r.text);
    EXPECT_EQ(std::vector<std::string>(out.begin(), out.begin() + 5),
              (std::vector<std::string>{"HEADER    RNA", "REMARK   1 kept", "SEQRES   1 A    1    G", "CRYST1   10.000",
                                        "MODEL        1"}));
    EXPECT_TRUE(out[6].starts_with("ANISOU    1  P     G A   1     1000"));
    EXPECT_TRUE(out[6].ends_with("P"));
    EXPECT_EQ(std::vector<std::string>(out.end() - 3, out.end()),
              (std::vector<std::string>{fmt("TER   %5d        G A   1", static_cast<int>(body.size())), "ENDMDL", "END"}));
    EXPECT_EQ(r.text.find("9.99"), std::string::npos);
    EXPECT_EQ(r.text.find("MASTER"), std::string::npos);
    EXPECT_TRUE(hasNote(r.notes, "Kept model 1 of 2"));
    EXPECT_TRUE(startsNote(r.notes, "Dropped the MASTER"));
    EXPECT_THROW(io::readPdbText("REMARK only\nEND\n"), io::PdbError);
    // CH-13: a path is no PDB text; a file that does not exist is a file error
    EXPECT_THROW(pt::PdbSource::file("/no/such/file.pdb").read(), std::filesystem::filesystem_error);
    EXPECT_THROW(pt::PdbSource::text("/no/such/file.pdb").read(), io::PdbError);
}

TEST(PdbToolsMdatk, test_record_policy) {
    const std::vector<std::string> water = {atomLine(90, "O", "HOH", "W", 1, 9, 9, 9, "HETATM"),
                                            atomLine(91, "H1", "HOH", "W", 1, 9, 9, 9, "HETATM"),
                                            atomLine(92, "MG", "MG", "W", 2, 8, 8, 8, "HETATM", "MG"), "TER"};
    const std::string text = makePdb({{"A", {{"G", 1, 'G', false}, {"0A", 2, 'A', true}, {"DG", 3, 'G', false}}}},
                                     {{"0A", "HETATM"}, {"G", "HETATM"}}, water);
    const auto records = [](const pt::CleanResult& r) {
        std::map<std::string, std::string> out;
        for (const std::string& l : lines(r.text)) {
            if (l.starts_with("ATOM") || l.starts_with("HETATM")) {
                out[trimmed(l.substr(17, 3))] = trimmed(l.substr(0, 6));
            }
        }
        return out;
    };
    using M = std::map<std::string, std::string>;
    EXPECT_EQ(records(clean(text)), (M{{"G", "HETATM"}, {"0A", "HETATM"}, {"DG", "ATOM"}, {"HOH", "HETATM"}, {"MG", "HETATM"}}));
    pt::CleanOptions o;
    o.records = "rosetta";
    const pt::CleanResult ros = clean(text, o);
    EXPECT_EQ(records(ros), (M{{"G", "ATOM"}, {"0A", "ATOM"}, {"DG", "ATOM"}, {"HOH", "HETATM"}, {"MG", "HETATM"}}));
    EXPECT_TRUE(startsNote(ros.notes, "Record policy rosetta"));
    o.records = "pdb";
    EXPECT_EQ(records(clean(text, o)), (M{{"G", "ATOM"}, {"0A", "HETATM"}, {"DG", "ATOM"}, {"HOH", "HETATM"}, {"MG", "HETATM"}}));
    EXPECT_EQ(records(clean(text, scheme("rosetta")))["0A"], "ATOM");
    EXPECT_EQ(records(clean(text, scheme("pdb")))["0A"], "HETATM");
    o.records = "charmm";
    EXPECT_THROW(clean(text, o), io::PdbError);
}

TEST(PdbToolsMdatk, test_conect_strip_and_remap) {
    const std::string text =
        makePdb({{"A", {{"G", 1, 'G', false}}}}, {}, {"CONECT  101  102  103", "CONECT  102  101"}, {}, 101);
    const pt::CleanResult r = clean(text);
    EXPECT_EQ(r.text.find("CONECT"), std::string::npos);
    EXPECT_TRUE(hasNote(r.notes, "Stripped 2 CONECT records"));
    pt::CleanOptions o;
    o.stripConect = false;
    const pt::CleanResult kept = clean(text, o);
    const std::vector<std::string> out = lines(kept.text);
    EXPECT_EQ(std::vector<std::string>(out.end() - 3, out.end()),
              (std::vector<std::string>{"CONECT    1    2    3", "CONECT    2    1", "END"}));
    // CH-04: repeated input serials: the CONECT records cannot be remapped, and the note says so
    const std::string twice = text.substr(0, text.size() - 4) + text.substr(0, text.size() - 4) + "END\n";
    const pt::CleanResult dropped = clean(twice, o);
    EXPECT_EQ(dropped.text.find("CONECT"), std::string::npos);
    EXPECT_TRUE(hasNote(dropped.notes, "WARNING: dropped 4 CONECT records: the input's atom serials repeat"));
}

// --- chains ----------------------------------------------------------------------------------------------

TEST(PdbToolsMdatk, test_restore_chains_from_reference) {
    const std::string text = makePdb(kCpptraj, {}, ions());
    const std::string ref = makePdb(kReference);
    pt::CleanOptions o = scheme("amber");
    o.chainReference = pt::PdbSource::text(ref);
    pt::CleanResult r = clean(text, o);
    EXPECT_EQ(names(r.structure), (Names{{"A", 1, "G5"}, {"A", 2, "C"}, {"A", 3, "U3"}, {"B", 4, "LG5"}, {"B", 5, "LA"},
                                         {"B", 6, "LC3"}, {"", 7, "Na+"}, {"", 8, "WAT"}}));
    EXPECT_TRUE(startsNote(r.notes, "Restored chain ids of 6 residues"));
    pt::CleanOptions numbers;
    numbers.chainReference = pt::PdbSource::text(ref);
    numbers.copyNumbers = true;
    r = clean(text, numbers);
    const Names n = names(r.structure);
    for (std::size_t k = 0; k < 6; ++k) {
        EXPECT_EQ(std::get<0>(n[k]), k < 3 ? "A" : "B");
        EXPECT_EQ(std::get<1>(n[k]), static_cast<std::int64_t>(k % 3 + 1));
    }
    Segments one = {{"", kCpptraj[0].second}};
    one[0].second.insert(one[0].second.end(), kCpptraj[1].second.begin(), kCpptraj[1].second.end());
    pt::CleanOptions plain;
    plain.chainReference = pt::PdbSource::text(ref);
    EXPECT_EQ(clean(makePdb(one), plain).structure.segments.size(), 2u);  // TER between chains
}

TEST(PdbToolsMdatk, test_restore_chains_mismatch_errors) {
    const std::string text = makePdb(kCpptraj);
    pt::CleanOptions o;
    o.chainReference = pt::PdbSource::text(
        makePdb({kReference[0], {"B", {kReference[1].second[0], kReference[1].second[1]}}}));
    EXPECT_NE(errorOf([&] { clean(text, o); }).find("has 6 residues but the reference has 5"), std::string::npos);
    o.chainReference = pt::PdbSource::text(
        makePdb({kReference[0], {"B", {{"LG", 1, 'G', true}, {"LU", 2, 'U', true}, {"LC", 3, 'C', true}}}}));
    EXPECT_NE(errorOf([&] { clean(text, o); }).find("sequence differs from the reference at 1 residue(s): #5"),
              std::string::npos);
}

TEST(PdbToolsMdatk, test_restore_chains_from_list_and_alphabetical) {
    const std::string text = makePdb(kCpptraj, {}, ions());
    pt::CleanOptions o;
    o.chainIds = std::vector<std::string>{"X", "Y"};
    EXPECT_EQ(chainsOf(clean(text, o).structure), (std::vector<std::string>{"X", "X", "X", "Y", "Y", "Y", "", ""}));
    pt::CleanOptions assign;
    assign.assignChains = true;
    const pt::CleanResult r = clean(text, assign);
    EXPECT_EQ(chainsOf(r.structure), (std::vector<std::string>{"A", "A", "A", "B", "B", "B", "", ""}));
    EXPECT_TRUE(hasNote(r.notes, "Assigned chain ids A, B to 2 TER segments"));
    o.chainIds = std::vector<std::string>{"A"};
    EXPECT_NE(errorOf([&] { clean(text, o); }).find("1 chain ids given for 2 TER segments"), std::string::npos);
}

TEST(PdbToolsMdatk, test_split_and_combine_chains) {
    const std::string text = makePdb({{"A", {{"G", 1, 'G', false}, {"C", 2, 'C', false}}},
                                      {"A", {{"A", 3, 'A', false}}},
                                      {"B", {{"U", 1, 'U', false}}}});
    pt::CleanOptions o;
    o.splitSegments = true;
    EXPECT_EQ(chainsOf(clean(text, o).structure), (std::vector<std::string>{"A", "A", "C", "B"}));
    std::vector<Res> gcaug;
    for (int k = 0; k < 5; ++k) {
        const char b = "GCAUG"[k];
        gcaug.push_back({std::string(1, b), k + 1, b, false});
    }
    const std::string five = makePdb({{"A", gcaug}});
    pt::CleanOptions split = scheme("amber");
    split.splitAt = "A:3";
    const pt::CleanResult r = clean(five, split);
    EXPECT_EQ(names(r.structure), (Names{{"A", 1, "G5"}, {"A", 2, "C3"}, {"B", 3, "A5"}, {"B", 4, "U"}, {"B", 5, "G3"}}));
    EXPECT_EQ(count(r.text, "TER"), 2u);
    split.splitAt = "A:9";
    EXPECT_THROW(clean(five, split), io::PdbError);
    pt::CleanOptions combine;
    combine.combine = {"A", "B"};
    const pt::CleanResult c = clean(text, combine);
    for (const std::string& chain : chainsOf(c.structure)) {
        EXPECT_EQ(chain, "A");
    }
    EXPECT_EQ(count(c.text, "TER"), 3u);
    EXPECT_TRUE(hasNote(c.notes, "duplicate residue numbers (1)"));
    combine.combineInto = "Z";
    combine.renumber = true;
    EXPECT_EQ(names(clean(text, combine).structure), (Names{{"Z", 1, "G"}, {"Z", 2, "C"}, {"Z", 3, "A"}, {"Z", 4, "U"}}));
    // CH-03: the merged chain's id has one character
    combine.combineInto = "XY";
    EXPECT_EQ(errorOf([&] { clean(text, combine); }), "chain ids must be single characters");
}

TEST(PdbToolsMdatk, test_order_chains) {
    const std::string text = makePdb({{"C", {{"G", 1, 'G', false}}}, {"A", {{"C", 2, 'C', false}}},
                                      {"B", {{"A", 3, 'A', false}}}}, {}, ions());
    pt::CleanOptions o;
    o.sortChains = true;
    EXPECT_EQ(clean(text, o).structure.chainIds(), (std::vector<std::string>{"A", "B", "C", ""}));
    pt::CleanOptions given;
    given.chainOrder = {"B"};
    EXPECT_EQ(clean(text, given).structure.chainIds(), (std::vector<std::string>{"B", "C", "A", ""}));
}

// --- numbering and sequence ------------------------------------------------------------------------------

TEST(PdbToolsMdatk, test_renumber_residues) {
    struct Case {
        bool perChain;
        bool keepGaps;
        std::int64_t start;
        std::vector<std::int64_t> expected;
    };
    for (const Case& c : std::vector<Case>{{false, false, 1, {1, 2, 3, 4, 5, 6, 7}},
                                           {true, false, 1, {1, 2, 3, 4, 5, 1, 2}},
                                           {false, true, 1, {1, 2, 6, 6, 7, 8, 9}},
                                           {true, true, 1, {1, 2, 6, 6, 7, 1, 2}},
                                           {false, false, 100, {100, 101, 102, 103, 104, 105, 106}}}) {
        pt::CleanOptions o;
        o.renumber = true;
        o.perChain = c.perChain;
        o.keepGaps = c.keepGaps;
        o.start = c.start;
        const pt::CleanResult r = clean(numbered(kGapped), o);
        std::vector<std::int64_t> got;
        for (const io::PdbResidue* res : r.structure.residues()) {
            got.push_back(res->resseq);
        }
        EXPECT_EQ(got, c.expected);
        EXPECT_EQ(r.structure.residues()[3]->icode, c.keepGaps ? "A" : "");
        if (!c.keepGaps) {
            EXPECT_TRUE(hasNote(r.notes, "Cleared 1 insertion codes"));
        }
    }
    // CH-01: numbers below -999 do not fit the PDB columns
    pt::CleanOptions o;
    o.renumber = true;
    o.start = -2000;
    EXPECT_EQ(errorOf([&] { clean(numbered(kGapped), o); }),
              "residue number -2000 (A) does not fit the PDB columns (-999 to 9999)");
}

TEST(PdbToolsMdatk, test_seqres_records_and_placement) {
    std::vector<Res> chainA;
    const std::string seq = "GCAUGCAUGCAUGCA";
    for (std::size_t k = 0; k < seq.size(); ++k) {
        chainA.push_back({std::string(1, seq[k]), static_cast<int>(k + 1), seq[k], false});
    }
    const std::string text = makePdb({{"A", chainA}, {"B", {{"0G", 1, 'G', true}, {"0C", 2, 'C', true}}}}, {}, ions(),
                                     {"REMARK   1 note", "SEQRES   1 Z    1    A", "CRYST1   10.000"});
    pt::CleanOptions o;
    o.seqres = true;
    const pt::CleanResult r = clean(text, o);
    const std::vector<std::string> out = lines(r.text);
    EXPECT_EQ(std::vector<std::string>(out.begin(), out.begin() + 5),
              (std::vector<std::string>{"REMARK   1 note",
                                        "SEQRES   1 A   15    G   C   A   U   G   C   A   U   G   C   A   U   G",
                                        "SEQRES   2 A   15    C   A", "SEQRES   1 B    2   0G  0C", "CRYST1   10.000"}));
    const std::string& line = out[1];
    EXPECT_EQ(line.substr(0, 6), "SEQRES");
    EXPECT_EQ(line.substr(7, 3), "  1");
    EXPECT_EQ(line[11], 'A');
    EXPECT_EQ(line.substr(13, 4), "  15");
    EXPECT_EQ(line.substr(19, 3), "  G");
    EXPECT_EQ(line.substr(67, 3), "  G");
    EXPECT_TRUE(hasNote(r.notes, "Wrote SEQRES records for 2 chain(s)"));
}

TEST(PdbToolsMdatk, test_fasta) {
    const std::string text =
        makePdb({{"A", {{"G", 1, 'G', false}, {"C", 2, 'C', false}, {"LA", 3, 'A', true}, {"LU", 4, 'U', true}}},
                 {"B", {{"G", 1, 'G', true}, {"C", 2, 'C', true}}},
                 {"C", {{"DG", 1, 'G', false}, {"PSU", 2, 'U', false}}}},
                {}, ions());
    pt::CleanOptions o;
    o.makeFasta = true;
    o.title = "demo";
    EXPECT_EQ(clean(text, o).fasta.value_or(""),
              ">demo_A chain A | RNA, L-nucleotides at 3-4 | 4 residues\nGCAU\n"
              ">demo_B chain B | L-RNA | 2 residues\nGC\n"
              ">demo_C chain C | DNA | 2 residues\nGX\n");
    EXPECT_FALSE(clean(text).fasta.has_value());
}

TEST(PdbToolsMdatk, test_clean_end_to_end_from_path) {
    test::TempDir dir;
    const auto path = dir / "solvated.pdb";
    std::vector<std::string> extra = ions();
    extra.push_back("CONECT  500  501");
    test::writeFile(path, makePdb(kCpptraj, {}, extra, {"CRYST1  105.626  105.626  105.626 109.47 109.47 109.47"}));
    const auto ref = dir / "ref.pdb";
    test::writeFile(ref, makePdb(kReference));
    pt::CleanOptions o = scheme("ccd");
    o.chainReference = pt::PdbSource::file(ref);
    o.renumber = true;
    o.perChain = true;
    o.seqres = true;
    o.makeFasta = true;
    const pt::CleanResult r = pt::clean(pt::PdbSource::file(path), o);
    const Names n = names(r.structure);
    EXPECT_EQ(Names(n.begin(), n.begin() + 6),
              (Names{{"A", 1, "G"}, {"A", 2, "C"}, {"A", 3, "U"}, {"B", 1, "0G"}, {"B", 2, "0A"}, {"B", 3, "0C3"}}));
    EXPECT_TRUE(r.fasta.value_or("").starts_with(">solvated_A chain A | RNA | 3 residues\nGCU\n>solvated_B chain B | L-RNA"));
    const std::vector<std::string> out = lines(r.text);
    EXPECT_EQ(out[0], "SEQRES   1 A    3    G   C   U");
    EXPECT_TRUE(out[2].starts_with("CRYST1"));
    EXPECT_EQ(out.back(), "END");
    EXPECT_EQ(r.text.find("CONECT"), std::string::npos);
    std::vector<int> serials;
    std::size_t ters = 0;
    for (const std::string& l : out) {
        ters += l.starts_with("TER");
        if (l.starts_with("ATOM") || l.starts_with("HETATM") || l.starts_with("TER")) {
            serials.push_back(std::stoi(l.substr(6, 5)));
        }
    }
    EXPECT_EQ(ters, 4u);
    for (std::size_t k = 0; k < serials.size(); ++k) {
        EXPECT_EQ(serials[k], static_cast<int>(k + 1));
    }
}

// --- command line (the library side of `pdb check / clean / to-lrna`; the commands: W2-C) ----------------

TEST(PdbToolsMdatk, test_mdatk_pdb_cli_check_clean_and_log) {
    const std::string text = makePdb({{"A", {{"G", 1, 'G', false}, {"C", 2, 'C', false}}},
                                      {"B", {{"G", 3, 'G', true}, {"U", 4, 'U', true}}}});
    const core::Value info = pt::check(pt::PdbSource::text(text));
    EXPECT_EQ(*info.path({"chains", "A", "chirality"}), core::Value(core::Table{{"D", core::Value(2)}}));
    EXPECT_EQ(*info.path({"chains", "B", "chirality"}), core::Value(core::Table{{"L", core::Value(2)}}));
    pt::CleanOptions o = scheme("amber");
    o.renumber = true;
    o.start = 10;
    o.seqres = true;
    o.makeFasta = true;
    const pt::CleanResult r = clean(text, o);
    EXPECT_EQ(names(io::readPdbText(r.text)), (Names{{"A", 10, "G5"}, {"A", 11, "C3"}, {"B", 12, "LG5"}, {"B", 13, "LU3"}}));
    EXPECT_NE(r.text.find("SEQRES"), std::string::npos);
    EXPECT_NE(r.fasta.value_or("").find('>'), std::string::npos);
    EXPECT_TRUE(startsNote(r.notes, "Chirality"));
    const pt::CleanResult l = pt::toLrna(pt::PdbSource::text(text), "A:2");
    const std::vector<std::string> n = resnames(io::readPdbText(l.text));
    EXPECT_EQ(std::vector<std::string>(n.begin(), n.begin() + 2), (std::vector<std::string>{"G", "LC3"}));
}

// test_embedded_copy_runs_standalone: milestone 7 (the embedded copy of the generated scripts)

TEST(PdbToolsMdatk, test_terminal_names_can_be_left_to_leap) {
    const std::string pdb = makePdb({{"C", {{"0G", 45, 'G', true}, {"0C", 46, 'C', true}, {"0A", 47, 'A', true}}}});
    pt::CleanOptions o = scheme("amber");
    o.terminalNames = false;
    EXPECT_EQ(resnames(clean(pdb, o).structure), (std::vector<std::string>{"LG", "LC", "LA"}));  // like to_lrna.sh
    EXPECT_EQ(resnames(clean(pdb, scheme("amber")).structure), (std::vector<std::string>{"LG5", "LC", "LA3"}));
    const pt::CleanResult out = pt::toLrna(pt::PdbSource::text(makePdb({{"A", {{"G", 1, 'G', true}, {"C", 2, 'C', true}}}})),
                                           "A:1-2", "amber", true, false);
    EXPECT_EQ(resnames(out.structure), (std::vector<std::string>{"LG", "LC"}));
}

TEST(PdbToolsMdatk, test_split_at_several_positions_in_file_order) {
    std::vector<Res> six;
    for (int i = 1; i < 7; ++i) {
        six.push_back({"G", i, 'G', false});
    }
    io::PdbStructure s = io::readPdbText(makePdb({{"A", six}}));
    pt::splitChains(s, "A:5, A:3");
    std::vector<std::pair<std::string, std::int64_t>> got;
    for (const io::PdbResidue* r : s.residues()) {
        got.emplace_back(r->chain, r->resseq);
    }
    EXPECT_EQ(got, (std::vector<std::pair<std::string, std::int64_t>>{{"A", 1}, {"A", 2}, {"B", 3}, {"B", 4}, {"C", 5}, {"C", 6}}));
}

TEST(PdbToolsMdatk, test_single_residue_l_segment_has_no_false_terminal_warning) {
    const std::string once = clean(makePdb({{"C", {{"0G", 45, 'G', true}}}}), scheme("amber")).text;
    const pt::CleanResult again = clean(once, scheme("amber"));
    EXPECT_EQ(resnames(again.structure), (std::vector<std::string>{"LG5"}));
    EXPECT_FALSE(hasNote(again.notes, "missing TER"));
}

TEST(PdbToolsMdatk, test_strip_5prime_phosphates_for_amber) {
    const std::string text = makePdb(kHybrid);
    pt::CleanOptions o = scheme("amber");
    o.strip5p = true;
    const pt::CleanResult r = clean(text, o);
    EXPECT_TRUE(hasNote(r.notes, "Removed the 5'-terminal phosphate of"));
    EXPECT_FALSE(hasNote(r.notes, "carry P/OP1/OP2"));  // the Amber warning is gone
    for (const auto& seg : r.structure.segments) {
        for (const io::PdbResidue& res : seg) {
            if (pt::nucleotide(res)) {
                for (const io::PdbAtom& a : res.atoms) {
                    EXPECT_TRUE(a.name != "P" && a.name != "OP1" && a.name != "OP2" && a.name != "OP3") << res.label();
                }
                break;
            }
        }
    }
    const auto& second = r.structure.segments[0][1].atoms;  // a later residue keeps its phosphate
    EXPECT_TRUE(std::any_of(second.begin(), second.end(), [](const io::PdbAtom& a) { return a.name == "P"; }));
    for (const std::string& n : clean(text, scheme("amber")).notes) {
        EXPECT_EQ(n.find("Removed"), std::string::npos);
    }
    std::vector<std::string> ls = lines(text);  // hydrogens: HO5' 0.96 A from O5', along the old O5'-P bond
    const auto p = std::find_if(ls.begin(), ls.end(), [](const std::string& l) { return trimmed(l.substr(12, 4)) == "P"; });
    const std::string h = p->substr(0, 12) + " H5'" + p->substr(16);
    ls.insert(p, h);
    pt::CleanOptions strip;
    strip.strip5p = true;
    const pt::CleanResult withH = clean(joinLines(ls), strip);
    const pt::AtomCoords c = withH.structure.segments[0][0].coords();
    const auto at = [&](const char* name) {
        return std::find_if(c.begin(), c.end(), [&](const auto& e) { return e.first == name; });
    };
    ASSERT_NE(at("HO5'"), c.end());
    const pt::Vec3 a = at("HO5'")->second;
    const pt::Vec3 b = at("O5'")->second;
    EXPECT_NEAR(std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])),
                0.96, 0.01);
}

TEST(PdbToolsMdatk, test_strip_5p_keeps_internal_phosphates_and_drops_lone_phosphates) {
    pt::CleanOptions o;
    o.strip5p = true;
    const pt::CleanResult r =
        clean(makePdb({{"A", {{"XYZ", 1, 'G', false}, {"G", 2, 'G', false}, {"C", 3, 'C', false}}}}), o);
    const auto& g2 = r.structure.segments[0][1].atoms;  // G2 is not a 5' end here
    EXPECT_TRUE(std::any_of(g2.begin(), g2.end(), [](const io::PdbAtom& a) { return a.name == "P"; }));
    EXPECT_TRUE(hasNote(r.notes, "does not start with a recognised nucleotide"));
    std::vector<std::string> ls;  // residue 1 keeps only its phosphate
    for (const std::string& l : lines(makePdb({{"A", {{"G", 1, 'G', false}, {"C", 2, 'C', false}}}}))) {
        const bool drop = l.starts_with("ATOM") && trimmed(l.substr(22, 4)) == "1" &&
                          !std::set<std::string>{"P", "OP1", "OP2"}.count(trimmed(l.substr(12, 4)));
        if (!drop) {
            ls.push_back(l);
        }
    }
    const pt::CleanResult lone = clean(joinLines(ls), o);
    const auto res = lone.structure.residues();
    ASSERT_EQ(res.size(), 1u);
    EXPECT_EQ(res[0]->resseq, 2);
    EXPECT_TRUE(std::none_of(res[0]->atoms.begin(), res[0]->atoms.end(), [](const io::PdbAtom& a) { return a.name == "P"; }));
    EXPECT_TRUE(hasNote(lone.notes, "held only a phosphate"));
}
