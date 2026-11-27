// The structure builder's base mutations and fit to a map (src/chem/mutate.cpp, fit_to_map.cpp)
// against what MDATK's mdatk/builder.py computes with MDBridge's improvements applied (CH-08, CH-09,
// CH-10, CH-11; tests/parity/make_builder_cases.py -> data/chem/builder_cases.json "mutations",
// builder_fitted_cases.json "fits" with the maps of data/chem/maps), and the counterparts of MDATK's
// tests/test_builder.py tests of mutations and fits (test_mutations_keep_the_frame_and_compensate,
// test_l_rna_mutation_keeps_l_names_and_chirality, test_mutation_options_and_errors,
// test_global_fit_recovers_the_pose, test_local_fit_and_anisou_rotation and the DNA-partner and
// altloc part of test_negative_start_empty_partner_dna_partner_and_altlocs).
//
// Mutated PDB text is compared byte for byte; the Kabsch rotation of a base swap comes from
// MDBridge's own SVD (NumPy: LAPACK), so a line whose rounded coordinate differs is accepted with
// every coordinate within 0.0011 A (docs/specs/m6-plan.md 4.1), and the number of such lines is
// printed. The fits are compared within D10's tolerance: coordinates RMSD <= 0.2 A and every atom
// within 0.5 A of MDATK's, CC within 1e-3, the notes equal with their numbers within the printed
// rounding of that tolerance.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/fit_to_map.hpp"
#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/mutate.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/io/pdb.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::chem;
using core::Array;
using core::Value;
using test::kData;
namespace fs = std::filesystem;

namespace {

const Value& builderCases() {
    static const Value data =
        core::parseJson(test::readFile(kData / "chem" / "builder_cases.json"), "builder_cases.json");
    return data;
}
const Value& fittedCases() {
    static const Value data =
        core::parseJson(test::readFile(kData / "chem" / "builder_fitted_cases.json"), "builder_fitted_cases.json");
    return data;
}

const Array& list(const Value& value) {
    return *value.get<Array>();
}
const Value& at(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
std::string str(const Value& value) {
    return *value.get<std::string>();
}
std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    for (const Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}
double hexDouble(const Value& value) {
    return std::strtod(str(value).c_str(), nullptr);
}

template <class Fn>
std::optional<std::pair<std::string, std::string>> errorOf(Fn&& fn) {
    try {
        fn();
    } catch (const core::ValueError& error) {
        return std::make_pair(std::string(error.pythonName()), std::string(error.what()));
    }
    return std::nullopt;
}

void expectError(const std::optional<std::pair<std::string, std::string>>& got, const Value& expected,
                 const std::string& where) {
    ASSERT_TRUE(got.has_value()) << where << ": expected " << str(at(expected, "message"));
    EXPECT_EQ(got->first, str(at(expected, "type"))) << where;
    EXPECT_EQ(got->second, str(at(expected, "message"))) << where;
}

std::vector<std::string> lines(const std::string& text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {
        const std::size_t end = text.find('\n', start);
        out.push_back(text.substr(start, end - start));
        start = end == std::string::npos ? text.size() : end + 1;
    }
    return out;
}

bool atomLine(const std::string& line) {
    return line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0;
}

// byte for byte, or the coordinates of a differing atom line within `tolerance`; returns the number
// of atom lines that differ
std::size_t expectPdbText(const std::string& got, const std::string& expected, const std::string& where,
                          double tolerance = 0.0011) {
    if (got == expected) {
        return 0;
    }
    const std::vector<std::string> a = lines(got);
    const std::vector<std::string> b = lines(expected);
    EXPECT_EQ(a.size(), b.size()) << where;
    std::size_t differ = 0;
    for (std::size_t k = 0; k < std::min(a.size(), b.size()); ++k) {
        if (a[k] == b[k]) {
            continue;
        }
        if (!atomLine(a[k]) || a[k].size() < 54 || b[k].size() < 54) {
            ADD_FAILURE() << where << ": line " << k << "\n  got      " << a[k] << "\n  expected " << b[k];
            continue;
        }
        ++differ;
        EXPECT_EQ(a[k].substr(0, 30), b[k].substr(0, 30)) << where;
        EXPECT_EQ(a[k].substr(54), b[k].substr(54)) << where;
        for (std::size_t c = 30; c < 54; c += 8) {
            EXPECT_NEAR(std::stod(a[k].substr(c, 8)), std::stod(b[k].substr(c, 8)), tolerance)
                << where << ": " << a[k];
        }
    }
    return differ;
}

MutateOptions mutateOptions(const Value& options) {
    MutateOptions out;
    for (const auto& [key, value] : *options.get<core::Table>()) {
        if (key == "compensatory") {
            out.compensatory = *value.get<bool>();
        } else if (key == "wobble") {
            out.wobble = *value.get<bool>();
        } else {
            throw std::runtime_error("unknown mutate option " + key);
        }
    }
    return out;
}

using Vec = std::array<double, 3>;

double dist(const Vec& a, const Vec& b) {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

std::map<std::string, Vec> coordsOf(const io::PdbResidue& res) {
    std::map<std::string, Vec> out;
    for (const auto& [name, xyz] : res.coords()) {
        out.emplace(name, xyz);
    }
    return out;
}

using Key = std::pair<std::string, std::int64_t>;
std::map<Key, io::PdbResidue> byKey(const io::PdbStructure& s) {
    std::map<Key, io::PdbResidue> out;
    for (const io::PdbResidue* r : s.residues()) {
        out[{r->chain, r->resseq}] = *r;
    }
    return out;
}

struct FoundPair {
    io::PdbResidue a;
    io::PdbResidue b;
    int kind = 0;
};
// tests/test_builder.py pairs_of: the builtin detector's pairs
std::vector<FoundPair> pairsOf(const std::string& text, bool wobble = true) {
    const io::PdbStructure s = io::readPdbText(text);
    const std::vector<const io::PdbResidue*> residues = s.residues();
    std::vector<std::string> names;
    std::vector<AtomIndex> index;
    std::vector<double> xyz;
    std::int64_t n = 0;
    for (const io::PdbResidue* r : residues) {
        AtomIndex rows;
        for (const auto& [name, p] : r->coords()) {
            rows.emplace(name, n++);
            xyz.insert(xyz.end(), p.begin(), p.end());
        }
        index.push_back(std::move(rows));
        names.push_back(r->name);
    }
    std::vector<FoundPair> out;
    for (const BasePair& p : PairDetector(names, index, 3.5, wobble).detect(xyz)) {
        out.push_back({*residues[static_cast<std::size_t>(p.i)], *residues[static_cast<std::size_t>(p.j)], p.kind});
    }
    return out;
}

bool anyContains(const Notes& notes, std::string_view text) {
    return std::any_of(notes.begin(), notes.end(), [&](const std::string& n) { return n.find(text) != std::string::npos; });
}

HelixOptions helix(std::optional<std::string> partner = "auto") {
    HelixOptions o;
    o.partner = std::move(partner);
    return o;
}

// The current directory, restored on leaving the scope
class WorkingDirectory {
public:
    explicit WorkingDirectory(const fs::path& dir) : m_old(fs::current_path()) { fs::current_path(dir); }
    ~WorkingDirectory() { fs::current_path(m_old); }
    WorkingDirectory(const WorkingDirectory&) = delete;
    WorkingDirectory& operator=(const WorkingDirectory&) = delete;

private:
    fs::path m_old;
};

}  // namespace

// --- recorded cases ----------------------------------------------------------------------------------

TEST(MutationCases, ParseMutations) {
    for (const Value& c : list(at(at(builderCases(), "mutations"), "parse_mutations"))) {
        const std::string spec = str(list(c)[0]);
        const Value& expected = list(c)[1];
        std::vector<Mutation> got;
        const auto error = errorOf([&] { got = parseMutations(spec); });
        if (const Value* e = expected.find("error")) {
            expectError(error, *e, spec);
            continue;
        }
        ASSERT_FALSE(error) << spec << ": " << error->second;
        const Array& want = list(at(expected, "value"));
        ASSERT_EQ(got.size(), want.size()) << spec;
        for (std::size_t k = 0; k < got.size(); ++k) {
            const Array& w = list(want[k]);
            EXPECT_EQ(got[k].chain, str(w[0])) << spec;
            EXPECT_EQ(got[k].resseq, *w[1].get<std::int64_t>()) << spec;
            EXPECT_EQ(got[k].icode, str(w[2])) << spec;
            EXPECT_EQ(std::string(1, got[k].base), str(w[3])) << spec;
            EXPECT_TRUE(got[k].isExplicit);
        }
    }
}

TEST(MutationCases, Mutate) {
    const Value& m = at(builderCases(), "mutations");
    const Value& inputs = at(m, "inputs");
    std::size_t cases = 0;
    std::size_t lineDiffs = 0;
    for (const Value& c : list(at(m, "mutate"))) {
        const std::string input = str(at(c, "input"));
        const std::string spec = str(at(c, "spec"));
        const std::string where = input + " " + spec;
        MutateResult got;
        const auto error = errorOf(
            [&] { got = mutate(PdbSource::text(str(at(inputs, input))), spec, mutateOptions(at(c, "options"))); });
        ++cases;
        if (const Value* e = c.find("error")) {
            expectError(error, *e, where);
            continue;
        }
        ASSERT_FALSE(error) << where << ": " << error->second;
        lineDiffs += expectPdbText(got.text, str(at(c, "text")), where);
        EXPECT_EQ(got.notes, strings(at(c, "notes"))) << where;
    }
    EXPECT_EQ(cases, 32u);
    std::printf("[mutations] %zu cases, %zu atom line(s) not byte-identical (within 0.0011 A)\n", cases, lineDiffs);
}

// The inputs the generator built are MDBridge's own Watson-Crick helices
TEST(MutationCases, InputsAreMdbridgeHelices) {
    const Value& inputs = at(at(builderCases(), "mutations"), "inputs");
    EXPECT_EQ(buildHelix("GGCAUCG"), str(at(inputs, "helix")));
    HelixOptions h;
    h.hydrogens = true;
    EXPECT_EQ(buildHelix("GGCAUCG", h), str(at(inputs, "helix_h")));
    HelixOptions l;
    l.chirality = "L";
    l.scheme = "ccd";
    l.hydrogens = true;
    EXPECT_EQ(buildHelix("GAUC", l), str(at(inputs, "l_ccd_h")));
}

TEST(MutationCases, MutatedNames) {
    EXPECT_EQ(mutatedName("LG5", 'C'), "LC5");
    EXPECT_EQ(mutatedName("GUA", 'U'), "URA");
    EXPECT_EQ(mutatedName("RA3", 'G'), "RG3");
    EXPECT_EQ(mutatedName("0A", 'U'), "0U");
    EXPECT_EQ(mutatedName(" g ", 'A'), "A");
    EXPECT_EQ(mutatedName("URI", 'C'), "CYT");
    EXPECT_EQ(mutatedName("C5", 'G'), "G5");
    EXPECT_EQ(mutatedName("UN", 'A'), "AN");
    const auto error = errorOf([] { (void)mutatedName("DA", 'G'); });
    ASSERT_TRUE(error);
    EXPECT_EQ(error->first, "BuildError");
    EXPECT_EQ(error->second, "cannot rename residue 'DA' (RNA nucleotides only)");
    EXPECT_TRUE(errorOf([] { (void)mutatedName("LGG", 'G'); }));
    EXPECT_TRUE(errorOf([] { (void)mutatedName("A55", 'G'); }));
}

TEST(MutationCases, ListsOfSpecs) {
    const std::vector<std::string> items = {"A:5G", "", "B:-2>C"};
    const std::vector<Mutation> got = parseMutations(items);
    ASSERT_EQ(got.size(), 2u);
    EXPECT_EQ(got[1].where(), "B:-2");
    EXPECT_EQ(got[1].base, 'C');
    EXPECT_TRUE(errorOf([] { (void)parseMutations(std::vector<std::string>{}); }));
    EXPECT_TRUE(errorOf([] { (void)parseMutations("A:99999999999999999999999G"); }));
}

// --- MDATK's tests -------------------------------------------------------------------------------

// test_mutations_keep_the_frame_and_compensate, with CH-09: the mutated bases get ideal hydrogens
TEST(Mutations, KeepTheFrameAndCompensate) {
    HelixOptions h;
    h.hydrogens = true;
    const std::string src = buildHelix("GGCAUCG", h);
    const std::map<Key, io::PdbResidue> before = byKey(io::readPdbText(src));
    const MutateResult r = mutate(PdbSource::text(src), "A:1C, A:4G");
    io::PdbStructure s = io::readPdbText(r.text);
    const std::map<Key, io::PdbResidue> after = byKey(s);
    EXPECT_EQ(after.at({"A", 1}).name, "C");
    EXPECT_EQ(after.at({"A", 4}).name, "G");
    EXPECT_EQ(after.at({"B", 7}).name, "G");  // the compensatory partners
    EXPECT_EQ(after.at({"B", 4}).name, "C");
    EXPECT_TRUE(anyContains(r.notes, "compensatory, compensates A:1 C"));
    EXPECT_TRUE(anyContains(r.notes, "A:4: replaced 4 hydrogen(s) of the old base (H8 H61 H62 H2) with 4 of the new "
                                     "base (H8 H1 H21 H22) at ideal positions"));
    // sugar / backbone atoms verbatim, the new glycosidic N on the old one
    const io::PdbResidue& old = before.at({"A", 4});
    const io::PdbResidue& now = after.at({"A", 4});
    const std::map<std::string, Vec> oldXyz = coordsOf(old);
    const std::map<std::string, Vec> newXyz = coordsOf(now);
    for (const io::PdbAtom& a : old.atoms) {
        if (a.name.find('\'') != std::string::npos || a.name == "P" || a.name == "OP1" || a.name == "OP2") {
            EXPECT_EQ(newXyz.at(a.name), oldXyz.at(a.name)) << a.name;
        }
    }
    EXPECT_LT(dist(newXyz.at("N9"), oldXyz.at("N9")), 0.05);
    EXPECT_LT(dist(coordsOf(after.at({"A", 1})).at("N1"), coordsOf(before.at({"A", 1})).at("N9")), 0.1);
    for (const io::PdbAtom& a : now.atoms) {  // the old base (adenine) and its hydrogens are gone
        EXPECT_TRUE(a.name != "H61" && a.name != "H62" && a.name != "H2" && a.name != "N6") << a.name;
    }
    // CH-09: the new base's hydrogens 1.0 A (C-H 1.08) from their parents, after the heavy atoms
    std::vector<std::string> names;
    for (const io::PdbAtom& a : now.atoms) {
        names.push_back(a.name);
    }
    const auto pos = [&](const std::string& n) { return std::find(names.begin(), names.end(), n) - names.begin(); };
    EXPECT_LT(pos("C4"), pos("H8"));
    EXPECT_LT(pos("H22"), static_cast<std::ptrdiff_t>(names.size()));
    EXPECT_NEAR(dist(newXyz.at("H8"), newXyz.at("C8")), 1.08, 0.01);
    EXPECT_NEAR(dist(newXyz.at("H1"), newXyz.at("N1")), 1.01, 0.02);
    EXPECT_NEAR(dist(newXyz.at("H21"), newXyz.at("N2")), 1.01, 0.02);
    EXPECT_NEAR(dist(newXyz.at("H22"), newXyz.at("N2")), 1.01, 0.02);
    for (const io::PdbAtom& a : now.atoms) {
        if (a.name == "H8") {
            EXPECT_EQ(a.element, "H");
        }
    }
    // the mutated helix has the same hydrogen set as a helix built with these bases
    HelixOptions mh;
    mh.hydrogens = true;
    const io::PdbStructure built = io::readPdbText(buildHelix("CGCGUCG", mh));
    const std::map<Key, io::PdbResidue> ideal = byKey(built);
    for (const auto& [key, res] : after) {
        std::vector<std::string> a;
        std::vector<std::string> b;
        for (const io::PdbAtom& atom : res.atoms) {
            a.push_back(atom.name);
        }
        for (const io::PdbAtom& atom : ideal.at(key).atoms) {
            b.push_back(atom.name);
        }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        EXPECT_EQ(a, b) << key.first << ":" << key.second;
    }
    const std::vector<FoundPair> found = pairsOf(r.text);
    EXPECT_EQ(found.size(), 7u);
    for (const FoundPair& p : found) {
        EXPECT_EQ(p.kind, kWatsonCrick);
    }
    const Notes chir = assignChirality(s);
    ASSERT_FALSE(chir.empty());
    const std::string_view tail = " 0 L and 14 D nucleotides";
    EXPECT_TRUE(chir[0].size() >= tail.size() && chir[0].substr(chir[0].size() - tail.size()) == tail) << chir[0];
}

// test_l_rna_mutation_keeps_l_names_and_chirality
TEST(Mutations, LRnaKeepsLNamesAndChirality) {
    const std::vector<std::pair<std::optional<std::string>, std::array<std::string, 4>>> runs = {
        {"amber", {"LG5", "LC5", "LC3", "LG3"}}, {"ccd", {"0G5", "0C5", "0C3", "0G3"}}, {std::nullopt, {"0G", "0C", "0C", "0G"}}};
    for (const auto& [scheme, names] : runs) {
        HelixOptions o;
        o.chirality = "L";
        o.scheme = scheme;
        const std::string src = buildHelix("GAUC", o);
        const MutateResult r = mutate(PdbSource::text(src), "A:1C");
        io::PdbStructure s = io::readPdbText(r.text);
        const std::map<Key, io::PdbResidue> res = byKey(s);
        const std::map<Key, io::PdbResidue> old = byKey(io::readPdbText(src));
        EXPECT_EQ(old.at({"A", 1}).name, names[0]);
        EXPECT_EQ(res.at({"A", 1}).name, names[1]);
        EXPECT_EQ(old.at({"B", 4}).name, names[2]);
        EXPECT_EQ(res.at({"B", 4}).name, names[3]);
        const Notes chir = assignChirality(s);
        EXPECT_EQ(chir.at(0).rfind("Chirality (sugar geometry): 8 L and 0 D", 0), 0u) << chir.at(0);
        for (const Key& key : {Key{"A", 1}, Key{"B", 4}}) {
            const io::PdbResidue& r2 = res.at(key);
            const SugarChirality c = sugarChirality(r2.coords(), io::parseResidueName(r2.name)->base);
            ASSERT_TRUE(c.call);
            EXPECT_EQ(*c.call, 'L');
            for (const auto& v : c.volumes) {
                EXPECT_GT(std::abs(v.second), 1.5);
            }
        }
        const std::vector<FoundPair> found = pairsOf(r.text);
        EXPECT_EQ(found.size(), 4u);
        for (const FoundPair& p : found) {
            EXPECT_EQ(p.kind, kWatsonCrick);
        }
    }
}

// test_mutation_options_and_errors (CH-08: a spec without a chain is the blank chain)
TEST(Mutations, OptionsAndErrors) {
    const std::string src = buildHelix("GGCAUCG");
    Notes notes = mutate(PdbSource::text(src), "A:2U, B:6A").notes;
    EXPECT_FALSE(anyContains(notes, "compensatory"));
    notes = mutate(PdbSource::text(src), "A:2A, B:6A").notes;
    EXPECT_TRUE(anyContains(notes, "do not pair canonically"));
    MutateOptions no;
    no.compensatory = false;
    MutateResult r = mutate(PdbSource::text(src), "A:2U", no);
    EXPECT_EQ(io::readPdbText(r.text).residues()[7 + 5]->name, "C");
    // wobble: a G facing U stays, a U facing G stays
    const std::string wob = buildHelix("GGCAUCG", helix("CGAUGUC"));  // pairs A:2 G · B:6 U
    MutateOptions w;
    w.wobble = true;
    notes = mutate(PdbSource::text(wob), "A:6U", w).notes;
    EXPECT_TRUE(anyContains(notes, "already pairs (G·U wobble)"));
    r = mutate(PdbSource::text(wob), "A:6U");
    EXPECT_TRUE(anyContains(r.notes, "compensates A:6 U"));
    std::vector<int> kinds;
    for (const FoundPair& p : pairsOf(r.text)) {
        if ((p.a.chain == "A" && p.a.resseq == 6) || (p.b.chain == "A" && p.b.resseq == 6)) {
            kinds.push_back(p.kind);
        }
    }
    EXPECT_EQ(kinds, std::vector<int>{kWatsonCrick});
    const auto message = [&](std::string_view spec) {
        const auto e = errorOf([&] { (void)mutate(PdbSource::text(src), spec); });
        return e ? e->second : std::string();
    };
    EXPECT_NE(message("A:99G").find("not found"), std::string::npos);
    EXPECT_NE(message("A5").find("cannot read"), std::string::npos);
    const auto twice = errorOf([] { (void)parseMutations("A:5G, A:5C"); });
    ASSERT_TRUE(twice);
    EXPECT_NE(twice->second.find("twice"), std::string::npos);
    const std::vector<Mutation> m = parseMutations("A:5g; B:-2>C :7U");
    ASSERT_EQ(m.size(), 3u);
    EXPECT_EQ(std::make_tuple(m[0].chain, m[0].resseq, m[0].base), std::make_tuple(std::string("A"), std::int64_t{5}, 'G'));
    EXPECT_EQ(std::make_tuple(m[1].chain, m[1].resseq, m[1].base), std::make_tuple(std::string("B"), std::int64_t{-2}, 'C'));
    EXPECT_EQ(std::make_tuple(m[2].chain, m[2].resseq, m[2].base), std::make_tuple(std::string(), std::int64_t{7}, 'U'));
    const std::vector<Mutation> blank = parseMutations("15G");  // CH-08
    EXPECT_EQ(std::make_tuple(blank[0].chain, blank[0].resseq), std::make_tuple(std::string(), std::int64_t{15}));
}

// the DNA-partner and altloc part of test_negative_start_empty_partner_dna_partner_and_altlocs
TEST(Mutations, DnaPartnerAndAltlocs) {
    const std::string text = buildHelix("GGCAU");
    std::string hybrid;
    for (std::string ln : lines(text)) {
        if (ln.rfind("ATOM", 0) == 0 && ln[21] == 'B') {
            if (core::pythonStrip(ln.substr(12, 4)) == "O2'") {
                continue;
            }
            std::string name = core::pythonStrip(ln.substr(17, 3));
            while (!name.empty() && (name.back() == '3' || name.back() == '5')) {
                name.pop_back();
            }
            std::string field = "D" + name;
            field.insert(0, 3 - std::min<std::size_t>(3, field.size()), ' ');
            ln = ln.substr(0, 17) + field + ln.substr(20);
        }
        hybrid += ln + "\n";
    }
    const MutateResult h = mutate(PdbSource::text(hybrid), "A:2C");
    EXPECT_TRUE(anyContains(h.notes, "is not an RNA nucleotide; no compensatory change"));
    std::string alt;
    for (const std::string& ln : lines(text)) {
        if (ln.rfind("ATOM", 0) == 0 && ln[21] == 'A' && core::pythonStrip(ln.substr(22, 4)) == "3") {
            alt += ln.substr(0, 16) + "A" + ln.substr(17, 37) + "  0.60" + ln.substr(60) + "\n";
            alt += ln.substr(0, 16) + "B" + ln.substr(17, 37) + "  0.40" + ln.substr(60) + "\n";
        } else {
            alt += ln + "\n";
        }
    }
    MutateOptions no;
    no.compensatory = false;
    const MutateResult r = mutate(PdbSource::text(alt), "A:3G", no);
    std::size_t base = 0;
    for (const std::string& ln : lines(r.text)) {
        const std::string name = ln.size() > 16 ? core::pythonStrip(ln.substr(12, 4)) : "";
        if (ln.rfind("ATOM", 0) == 0 && ln[21] == 'A' && core::pythonStrip(ln.substr(22, 4)) == "3" &&
            (name == "N9" || name == "C8" || name == "O6")) {
            ++base;
            EXPECT_EQ(ln[16], ' ');
            EXPECT_EQ(core::pythonStrip(ln.substr(54, 6)), "1.00");
        }
    }
    EXPECT_EQ(base, 3u);
    EXPECT_TRUE(anyContains(r.notes, "alternate locations A,B"));
}

// A mutation of a model without hydrogens adds none (CH-09 only replaces hydrogens that were there),
// and unknown atoms of the base are kept with a warning
TEST(Mutations, HydrogensOnlyWhereThereWereAndUnknownAtomsKept) {
    const std::string src = buildHelix("GGCAUCG");
    const MutateResult r = mutate(PdbSource::text(src), "A:4G");
    const io::PdbStructure mutated = io::readPdbText(r.text);  // residues() points into it
    for (const io::PdbResidue* res : mutated.residues()) {
        for (const io::PdbAtom& a : res->atoms) {
            EXPECT_NE(a.element, "H") << res->label() << " " << a.name;
        }
    }
    std::string extra;
    for (const std::string& ln : lines(src)) {
        extra += ln + "\n";
        if (ln.rfind("ATOM", 0) == 0 && ln.substr(12, 4) == " N9 " && ln[21] == 'A' &&
            core::pythonStrip(ln.substr(22, 4)) == "4") {
            extra += ln.substr(0, 12) + " XX " + ln.substr(16) + "\n";
        }
    }
    const MutateResult u = mutate(PdbSource::text(extra), "A:4G");
    EXPECT_TRUE(anyContains(u.notes, "WARNING: A:4 G: kept unrecognised atom(s) XX"));
}

// --- the fit to a map -----------------------------------------------------------------------------

namespace {

// numbers of `text` (as written: "0.973", "59.0", "4") and the text with each replaced by '#'
std::pair<std::string, std::vector<std::string>> numbersOf(const std::string& text) {
    static const std::regex number(R"(-?\d+(\.\d+)?)");
    std::vector<std::string> numbers;
    std::string masked;
    std::size_t last = 0;
    for (auto it = std::sregex_iterator(text.begin(), text.end(), number); it != std::sregex_iterator(); ++it) {
        masked += text.substr(last, static_cast<std::size_t>(it->position()) - last) + "#";
        numbers.push_back(it->str());
        last = static_cast<std::size_t>(it->position() + it->length());
    }
    masked += text.substr(last);
    return {masked, numbers};
}

// notes equal, their numbers within D10 (CC 1e-3, rotation 0.5 degrees, shift 0.2 A) plus rounding
void expectFitNotes(const Notes& got, const Notes& expected, const std::string& where) {
    ASSERT_EQ(got.size(), expected.size()) << where;
    for (std::size_t k = 0; k < got.size(); ++k) {
        if (got[k] == expected[k]) {
            continue;
        }
        const auto [a, na] = numbersOf(got[k]);
        const auto [b, nb] = numbersOf(expected[k]);
        EXPECT_EQ(a, b) << where << "\n  got      " << got[k] << "\n  expected " << expected[k];
        if (na.size() != nb.size()) {
            continue;
        }
        for (std::size_t j = 0; j < na.size(); ++j) {
            const std::size_t dot = nb[j].find('.');
            const std::size_t decimals = dot == std::string::npos ? 0 : nb[j].size() - dot - 1;
            const double tolerance = decimals >= 3 ? 0.0015 : decimals == 1 ? 0.55 : 0.0;
            EXPECT_NEAR(std::stod(na[j]), std::stod(nb[j]), tolerance) << where << ": " << got[k];
        }
    }
}

void expectFittedText(const std::string& got, const std::string& expected, const std::string& where) {
    const std::vector<std::string> a = lines(got);
    const std::vector<std::string> b = lines(expected);
    ASSERT_EQ(a.size(), b.size()) << where;
    double sum = 0.0;
    double worst = 0.0;
    std::size_t n = 0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        if (!atomLine(b[k])) {
            EXPECT_EQ(a[k], b[k]) << where;
            continue;
        }
        EXPECT_EQ(a[k].substr(0, 30), b[k].substr(0, 30)) << where;
        EXPECT_EQ(a[k].substr(54), b[k].substr(54)) << where;
        const Vec p{std::stod(a[k].substr(30, 8)), std::stod(a[k].substr(38, 8)), std::stod(a[k].substr(46, 8))};
        const Vec q{std::stod(b[k].substr(30, 8)), std::stod(b[k].substr(38, 8)), std::stod(b[k].substr(46, 8))};
        const double d = dist(p, q);
        sum += d * d;
        worst = std::max(worst, d);
        ++n;
    }
    const double rmsd = n ? std::sqrt(sum / static_cast<double>(n)) : 0.0;
    EXPECT_LE(rmsd, 0.2) << where;
    EXPECT_LE(worst, 0.5) << where;
    std::printf("[fit] %s: RMSD %.4f A, max %.4f A to MDATK's model\n", where.c_str(), rmsd, worst);
}

}  // namespace

TEST(FitCases, FitToMap) {
    const Value& fits = at(fittedCases(), "fits");
    const Value& models = at(fits, "models");
    const WorkingDirectory cwd(kData / "chem");  // map paths relative to it, as the generator ran
    std::size_t cases = 0;
    for (const Value& c : list(at(fits, "cases"))) {
        const std::string model = str(at(c, "model"));
        const std::string map = str(at(c, "map"));
        const std::string search = str(at(c, "search"));
        const double resolution = hexDouble(at(c, "resolution"));
        FitToMapOptions o;
        o.search = search;
        o.angularStep = hexDouble(at(c, "step"));
        const std::string where = model + " " + map + " " + search + " " + str(at(c, "resolution"));
        FitToMapResult got;
        const auto error =
            errorOf([&] { got = fitToMap(PdbSource::text(str(at(models, model))), map, resolution, o); });
        ++cases;
        if (const Value* e = c.find("error")) {
            expectError(error, *e, where);
            continue;
        }
        ASSERT_FALSE(error) << where << ": " << error->second;
        expectFittedText(got.text, str(at(c, "text")), where);
        EXPECT_NEAR(got.cc, *at(c, "cc").get<double>(), 1e-3) << where;
        expectFitNotes(got.notes, strings(at(c, "notes")), where);
    }
    EXPECT_EQ(cases, 14u);
}

namespace {

// tests/test_builder.py _moved_map: a map simulated from the model moved by (rotvec, shift)
struct MovedMap {
    fs::path path;
    std::vector<core::Vec3> moved;
};
MovedMap movedMap(const fs::path& dir, const std::string& text, const core::Vec3& rotvec, const core::Vec3& shift,
                  double resolution = 4.0) {
    const ModelAtoms atoms = modelAtoms(io::readPdbText(text));
    core::RigidTransform t;
    t.rotation = core::rotvecMatrix(rotvec);
    t.translation = shift;
    MovedMap out;
    out.moved = t.apply(atoms.xyz);
    out.path = dir / "sim.mrc";
    io::writeMap(out.path, io::simulateMap(out.moved, io::atomWeights(atoms.elements), resolution));
    return out;
}

double rmsdTo(const std::string& text, const std::vector<core::Vec3>& reference) {
    const ModelAtoms atoms = modelAtoms(io::readPdbText(text));
    EXPECT_EQ(atoms.xyz.size(), reference.size());
    double sum = 0.0;
    for (std::size_t k = 0; k < std::min(atoms.xyz.size(), reference.size()); ++k) {
        const double d = dist(atoms.xyz[k], reference[k]);
        sum += d * d;
    }
    return std::sqrt(sum / static_cast<double>(reference.size()));
}

}  // namespace

// test_global_fit_recovers_the_pose
TEST(FitToMap, GlobalFitRecoversThePose) {
    const test::TempDir dir;
    const std::string model = buildHelix("GAAGAGAAGGC");
    const MovedMap m = movedMap(dir.path(), model, {0.4, -0.9, 0.3}, {30.0, -12.0, 55.0});
    const FitToMapResult r = fitToMap(PdbSource::text(model), m.path, 4.0);
    EXPECT_LT(rmsdTo(r.text, m.moved), 1.0);
    EXPECT_GT(r.cc, 0.9);
    EXPECT_TRUE(std::any_of(r.notes.begin(), r.notes.end(),
                            [](const std::string& n) { return n.rfind("Density fit (global search)", 0) == 0; }));
    // only the coordinates changed
    const std::vector<std::string> a = lines(model);
    const std::vector<std::string> b = lines(r.text);
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t k = 0; k < a.size(); ++k) {
        if (atomLine(a[k])) {
            EXPECT_EQ(a[k].substr(0, 30), b[k].substr(0, 30));
            EXPECT_EQ(a[k].substr(54), b[k].substr(54));
        } else {
            EXPECT_EQ(a[k], b[k]);
        }
    }
}

// test_local_fit_and_anisou_rotation
TEST(FitToMap, LocalFitAndAnisouRotation) {
    const test::TempDir dir;
    const std::string model = buildHelix("GGCAU", helix("none"));
    const MovedMap m = movedMap(dir.path(), model, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    const std::optional<std::string> start =
        io::transformPdbText(model, core::rotvecMatrix({0.03, -0.04, 0.02}), {0.8, -0.6, 0.4});
    ASSERT_TRUE(start);
    FitToMapOptions local;
    local.search = "local";
    const FitToMapResult r = fitToMap(PdbSource::text(*start), m.path, 4.0, local);
    EXPECT_LT(rmsdTo(r.text, m.moved), 0.3);
    EXPECT_GT(r.cc, 0.9);
    const auto error = errorOf([&] { (void)fitToMap(PdbSource::text(model), m.path, 0.0, local); });
    ASSERT_TRUE(error);
    EXPECT_NE(error->second.find("resolution"), std::string::npos);
}

// The damaged-column rewrite: atom lines whose coordinates are shifted out of their columns (a
// ChimeraX residue field " G5A:" and the later columns one to the right) are rewritten in standard columns, ANISOU
// records dropped, CONECT kept
TEST(FitToMap, DamagedColumnsAreRewritten) {
    const test::TempDir dir;
    HelixOptions o;
    o.partner = "none";
    o.scheme = "amber";
    const std::string model = buildHelix("GGCAU", o);
    const MovedMap m = movedMap(dir.path(), model, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    std::string damaged;
    for (const std::string& ln : lines(model)) {
        if (atomLine(ln)) {
            std::string fused = core::pythonStrip(ln.substr(17, 4)) + ln.substr(21, 1) + ":";
            fused.insert(0, 5 - std::min<std::size_t>(5, fused.size()), ' ');
            const std::string shifted = ln.substr(0, 17) + fused + " " + ln.substr(22);  // one column to the right
            damaged += shifted + "\n";
            if (ln.substr(12, 4) == " P  ") {
                damaged += "ANISOU" + ln.substr(6, 22) + "   1000   2000   3000      0      0      0       P\n";
            }
        } else if (ln == "END") {
            damaged += "CONECT    1    2\nEND\n";
        } else {
            damaged += ln + "\n";
        }
    }
    ASSERT_FALSE(io::transformPdbText(damaged, core::rotvecMatrix({0.0, 0.0, 0.0}), {0.0, 0.0, 0.0}));
    FitToMapOptions local;
    local.search = "local";
    const FitToMapResult r = fitToMap(PdbSource::text(damaged), m.path, 4.0, local);
    EXPECT_TRUE(anyContains(r.notes, "Repaired"));
    EXPECT_TRUE(anyContains(r.notes, "Rewrote the file in standard PDB columns (the input had damaged ATOM records); "
                                     "ANISOU records dropped"));
    EXPECT_EQ(r.text.find("ANISOU"), std::string::npos);
    EXPECT_NE(r.text.find("CONECT"), std::string::npos);
    const io::PdbStructure s = io::readPdbText(r.text);
    EXPECT_TRUE(s.notes.empty());  // standard columns now
    EXPECT_EQ(s.residues().size(), 5u);
    EXPECT_EQ(s.residues()[0]->name, io::readPdbText(model).residues()[0]->name);
    EXPECT_LT(rmsdTo(r.text, m.moved), 0.3);
    // the same coordinates as the fit of the undamaged model
    const FitToMapResult clean = fitToMap(PdbSource::text(model), m.path, 4.0, local);
    EXPECT_LT(rmsdTo(r.text, modelAtoms(io::readPdbText(clean.text)).xyz), 0.002);
}

TEST(FitToMap, ModelAtomsAndElements) {
    io::PdbAtom a;
    a.name = "MG";
    EXPECT_EQ(atomElement(a, "MG"), "MG");
    EXPECT_EQ(atomElement(a, "HOH"), "M");
    a.name = "C1'";
    EXPECT_EQ(atomElement(a, "G"), "C");
    a.name = "1H5'";
    EXPECT_EQ(atomElement(a, "G"), "H");
    a.element = "Se";
    EXPECT_EQ(atomElement(a, "G"), "Se");
    a.element.clear();
    a.name = "MG2+";
    EXPECT_EQ(atomElement(a, "MG+2"), "MG");
    HelixOptions h;
    h.hydrogens = true;
    h.partner = "none";
    const io::PdbStructure s = io::readPdbText(buildHelix("GC", h));
    const ModelAtoms heavy = modelAtoms(s);
    const ModelAtoms all = modelAtoms(s, false);
    EXPECT_LT(heavy.xyz.size(), all.xyz.size());
    for (const std::string& e : heavy.elements) {
        EXPECT_NE(e, "H");
    }
    EXPECT_TRUE(errorOf([&] {
        FitToMapOptions o;
        o.search = "local";
        (void)fitToMap(PdbSource::text("ATOM      1  P     G A   1       0.000   0.000   0.000  1.00  0.00           P\n"),
                       "none.mrc", 4.0, o);
    }));
}
