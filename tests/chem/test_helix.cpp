// The structure builder's templates, sequence functions and Watson-Crick helices (src/chem/helix.cpp)
// against what MDATK's mdatk/builder.py computes with MDBridge's improvements applied
// (tests/parity/make_builder_cases.py -> data/chem/builder_cases.json: "sequences", "templates",
// "helices"), and the geometry MDATK's tests/test_builder.py checks (test_ideal_duplex_geometry,
// test_l_helix_is_the_mirror_image_and_left_handed, test_helix_sugars_are_c3_endo,
// test_names_terminals_hydrogens_and_phosphate, test_single_strand_partner_and_validation_messages,
// test_explicit_partner_with_wobble_and_mismatch).
//
// The PDB text is compared byte for byte with glibc's sin / cos (the recorded bits are glibc's);
// with another C library every coordinate within 0.001 Angstrom and everything else equal.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/reference_tables.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/pdb.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::chem;
using core::Array;
using core::Value;
using test::kData;

namespace {

const Value& cases() {
    static const Value data =
        core::parseJson(test::readFile(kData / "chem" / "builder_cases.json"), "builder_cases.json");
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
std::optional<std::string> optStr(const Value& value) {
    if (value.isNull()) {
        return std::nullopt;
    }
    return str(value);
}
std::optional<std::string_view> view(const std::optional<std::string>& value) {
    if (!value) {
        return std::nullopt;
    }
    return std::string_view(*value);
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

bool glibcMath() {
#if defined(__GLIBC__)
    return true;
#else
    return false;
#endif
}

// Python's error type name of an exception thrown by the builder
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

void expectPdbText(const std::string& got, const std::string& expected, const std::string& where) {
    if (glibcMath()) {
        EXPECT_EQ(got, expected) << where;
        return;
    }
    const std::vector<std::string> a = lines(got);
    const std::vector<std::string> b = lines(expected);
    ASSERT_EQ(a.size(), b.size()) << where;
    for (std::size_t k = 0; k < a.size(); ++k) {
        if (a[k].rfind("ATOM", 0) == 0 || a[k].rfind("HETATM", 0) == 0) {
            ASSERT_GE(a[k].size(), 54u) << where;
            EXPECT_EQ(a[k].substr(0, 30), b[k].substr(0, 30)) << where;
            EXPECT_EQ(a[k].substr(54), b[k].substr(54)) << where;
            for (std::size_t c = 30; c < 54; c += 8) {
                EXPECT_NEAR(std::stod(a[k].substr(c, 8)), std::stod(b[k].substr(c, 8)), 0.0011) << where << ": " << a[k];
            }
        } else {
            EXPECT_EQ(a[k], b[k]) << where;
        }
    }
}

HelixOptions helixOptions(const Value& options) {
    HelixOptions out;
    for (const auto& [key, value] : *options.get<core::Table>()) {
        if (key == "partner") {
            out.partner = optStr(value);
        } else if (key == "chirality") {
            out.chirality = str(value);
        } else if (key == "chains") {
            out.chains = str(value);
        } else if (key == "start") {
            out.start = value.isString() ? str(value) : std::to_string(*value.get<std::int64_t>());
        } else if (key == "scheme") {
            out.scheme = optStr(value);
        } else if (key == "terminals") {
            out.terminals = *value.get<bool>();
        } else if (key == "phosphate5") {
            out.phosphate5 = *value.get<bool>();
        } else if (key == "hydrogens") {
            out.hydrogens = *value.get<bool>();
        } else {
            throw std::runtime_error("unknown helix option " + key);
        }
    }
    return out;
}

using Vec = std::array<double, 3>;

Vec sub(const Vec& a, const Vec& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
double dot(const Vec& a, const Vec& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vec crossv(const Vec& a, const Vec& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dist(const Vec& a, const Vec& b) {
    return std::sqrt(dot(sub(a, b), sub(a, b)));
}
double dihedral(const Vec& a, const Vec& b, const Vec& c, const Vec& d) {
    const Vec b0 = sub(a, b);
    Vec b1 = sub(c, b);
    const Vec b2 = sub(d, c);
    const double n = std::sqrt(dot(b1, b1));
    b1 = {b1[0] / n, b1[1] / n, b1[2] / n};
    const double s0 = dot(b0, b1);
    const double s2 = dot(b2, b1);
    const Vec v = {b0[0] - s0 * b1[0], b0[1] - s0 * b1[1], b0[2] - s0 * b1[2]};
    const Vec w = {b2[0] - s2 * b1[0], b2[1] - s2 * b1[1], b2[2] - s2 * b1[2]};
    return std::atan2(dot(crossv(b1, v), w), dot(v, w)) * 180.0 / 3.14159265358979323846;
}
// the signed volume (a - centre) . ((b - centre) x (c - centre)), independent of chem::chiralVolume
double tripleVolume(const Vec& centre, const Vec& a, const Vec& b, const Vec& c) {
    return dot(sub(a, centre), crossv(sub(b, centre), sub(c, centre)));
}

std::map<std::string, Vec> coordsOf(const io::PdbResidue& res) {
    std::map<std::string, Vec> out;
    for (const auto& [name, xyz] : res.coords()) {
        out.emplace(name, xyz);
    }
    return out;
}

std::map<std::string, Vec> templateMap(std::string_view base, bool hydrogens) {
    std::map<std::string, Vec> out;
    for (const TemplateAtom& a : helixTemplate(base, hydrogens)) {
        out[a.name] = a.xyz;
    }
    return out;
}

}  // namespace

// --- recorded cases ----------------------------------------------------------------------------------

TEST(HelixCases, CheckHelixMessages) {
    const Array& rows = list(at(at(cases(), "sequences"), "check_helix"));
    ASSERT_GT(rows.size(), 1000u);
    for (const Value& row : rows) {
        const Array& r = list(row);
        const std::string seq = str(r[0]);
        const auto partner = optStr(r[1]);
        const auto pairing = optStr(r[2]);
        const auto triplex = optStr(r[3]);
        const auto third = optStr(r[4]);
        std::vector<HelixMessage> expected;
        for (const Value& m : list(r[5])) {
            expected.push_back({str(list(m)[0]), str(list(m)[1])});
        }
        const std::string where = "check_helix(" + seq + ", " + partner.value_or("None") + ", " +
                                  pairing.value_or("None") + ", " + triplex.value_or("None") + ", " +
                                  third.value_or("None") + ")";
        EXPECT_EQ(checkHelix(seq, view(partner), view(pairing), view(triplex), view(third)), expected) << where;
    }
}

TEST(HelixCases, ParseSequence) {
    for (const Value& row : list(at(at(cases(), "sequences"), "parse_sequence"))) {
        const std::string text = str(list(row)[0]);
        const Value& result = list(row)[1];
        if (const Value* value = result.find("value")) {
            const ParsedSequence got = parseSequence(text);
            EXPECT_EQ(got.sequence, str(list(*value)[0])) << text;
            EXPECT_EQ(got.notes, strings(list(*value)[1])) << text;
        } else {
            expectError(errorOf([&] { parseSequence(text); }), at(result, "error"), text);
        }
    }
}

TEST(HelixCases, ResolvePartner) {
    for (const Value& row : list(at(at(cases(), "sequences"), "resolve_partner"))) {
        const Array& r = list(row);
        const std::string seq = str(r[0]);
        const auto partner = optStr(r[1]);
        const std::string pairing = str(r[2]);
        const std::string where = "resolve_partner(" + seq + ", " + partner.value_or("None") + ", " + pairing + ")";
        if (const Value* value = r[3].find("value")) {
            const ResolvedPartner got = resolvePartner(seq, view(partner), pairing);
            EXPECT_EQ(got.sequence, str(list(*value)[0])) << where;
            EXPECT_EQ(got.partner, optStr(list(*value)[1])) << where;
            EXPECT_EQ(got.notes, strings(list(*value)[2])) << where;
        } else {
            expectError(errorOf([&] { resolvePartner(seq, view(partner), pairing); }), at(r[3], "error"), where);
        }
    }
}

TEST(HelixCases, ResolveThird) {
    for (const Value& row : list(at(at(cases(), "sequences"), "resolve_third"))) {
        const Array& r = list(row);
        const std::string seq = str(r[0]);
        const auto third = optStr(r[1]);
        const std::string motif = str(r[2]);
        const std::string where = "resolve_third(" + seq + ", " + third.value_or("None") + ", " + motif + ")";
        if (const Value* value = r[3].find("value")) {
            const ResolvedThird got = resolveThird(seq, view(third), motif);
            EXPECT_EQ(got.sequence, str(list(*value)[0])) << where;
            EXPECT_EQ(got.notes, strings(list(*value)[1])) << where;
        } else {
            expectError(errorOf([&] { resolveThird(seq, view(third), motif); }), at(r[3], "error"), where);
        }
    }
}

TEST(HelixCases, HelixKind) {
    for (const Value& row : list(at(at(cases(), "sequences"), "helix_kind"))) {
        const Array& r = list(row);
        const auto pairing = optStr(r[0]);
        const auto triplex = optStr(r[1]);
        const std::string where = "helix_kind(" + pairing.value_or("None") + ", " + triplex.value_or("None") + ")";
        if (const Value* value = r[2].find("value")) {
            EXPECT_EQ(helixKind(view(pairing), view(triplex)), str(*value)) << where;
        } else {
            expectError(errorOf([&] { helixKind(view(pairing), view(triplex)); }), at(r[2], "error"), where);
        }
    }
}

TEST(HelixCases, PairTypesAndKinds) {
    const Value& seqs = at(cases(), "sequences");
    for (const Value& row : list(at(seqs, "pair_types"))) {
        const Array& r = list(row);
        EXPECT_EQ(pairTypes(str(r[0]), str(r[1])), strings(r[2])) << str(r[0]) << " " << str(r[1]);
    }
    const auto kindsOf = [](const Value& value) {
        std::vector<PairKind> out;
        for (const Value& k : list(value)) {
            out.push_back({str(list(k)[0])[0], str(list(k)[1])[0], str(list(k)[2])});
        }
        return out;
    };
    std::size_t n = 0;
    for (const Value& row : list(at(seqs, "pairing_kinds"))) {
        const Array& r = list(row);
        EXPECT_EQ(pairingKinds(str(r[0]), str(r[1]), str(r[2])), kindsOf(r[3])) << str(r[0]) << " " << str(r[2]);
        ++n;
    }
    for (const Value& row : list(at(seqs, "triplex_kinds"))) {
        const Array& r = list(row);
        EXPECT_EQ(triplexKinds(str(r[0]), str(r[1]), str(r[2])), kindsOf(r[3])) << str(r[0]) << " " << str(r[2]);
        ++n;
    }
    EXPECT_EQ(n, 28u);
}

TEST(HelixCases, ParseIds) {
    for (const Value& row : list(at(at(cases(), "sequences"), "parse_ids"))) {
        const Array& r = list(row);
        const std::string value = str(r[0]);
        const int count = static_cast<int>(*r[1].get<std::int64_t>());
        const bool numbers = str(r[2]) == "int";
        const std::string where = "_parse_ids(" + value + ", " + std::to_string(count) + ", " + str(r[2]) + ")";
        if (const Value* expected = r[3].find("value")) {
            if (numbers) {
                std::vector<std::int64_t> want;
                for (const Value& v : list(*expected)) {
                    want.push_back(*v.get<std::int64_t>());
                }
                EXPECT_EQ(parseStartNumbers(value, count), want) << where;
            } else {
                EXPECT_EQ(parseChainIds(value, count), strings(*expected)) << where;
            }
        } else {
            expectError(errorOf([&] {
                            if (numbers) {
                                parseStartNumbers(value, count);
                            } else {
                                parseChainIds(value, count);
                            }
                        }),
                        at(r[3], "error"), where);
        }
    }
}

TEST(HelixCases, Templates) {
    for (const auto& [base, atoms] : *at(cases(), "templates").get<core::Table>()) {
        const std::vector<TemplateAtom> got = helixTemplate(base, true);
        const core::Table& want = *atoms.get<core::Table>();
        ASSERT_EQ(got.size(), want.size()) << base;
        for (std::size_t k = 0; k < got.size(); ++k) {
            EXPECT_EQ(got[k].name, want[k].first) << base;
            const Array& xyz = list(want[k].second);
            for (std::size_t c = 0; c < 3; ++c) {
                if (glibcMath()) {
                    EXPECT_EQ(got[k].xyz[c], hexDouble(xyz[c])) << base << " " << got[k].name;
                } else {
                    EXPECT_NEAR(got[k].xyz[c], hexDouble(xyz[c]), 1e-12) << base << " " << got[k].name;
                }
            }
        }
        EXPECT_EQ(helixTemplate(base, false).size(),
                  reference::arnottBackbone().size() + reference::arnottBase(base[0]).size());
    }
    const auto err = errorOf([] { helixTemplate("T"); });
    ASSERT_TRUE(err);
    EXPECT_EQ(err->second, "no template for base 'T' (A, C, G, U)");
}

TEST(HelixCases, WatsonCrickHelices) {
    const Array& helices = list(at(cases(), "helices"));
    ASSERT_EQ(helices.size(), 80u);
    for (const Value& c : helices) {
        const std::string seq = str(at(c, "sequence"));
        const HelixOptions options = helixOptions(at(c, "options"));
        const std::string where = "build_helix(" + seq + ", " + core::pythonJsonCompact(at(c, "options")) + ")";
        if (const Value* error = c.find("error")) {
            std::vector<std::string> notes;
            expectError(errorOf([&] { buildHelix(seq, options, &notes); }), *error, where);
            EXPECT_TRUE(notes.empty()) << where;
            continue;
        }
        std::vector<std::string> notes;
        const std::string text = buildHelix(seq, options, &notes);
        expectPdbText(text, str(at(c, "text")), where);
        EXPECT_EQ(notes, strings(at(c, "notes"))) << where;
    }
}

TEST(HelixCases, OtherPairingsWaitForTheFittedBuilder) {
    HelixOptions options;
    options.pairing = "hh";
    const auto err = errorOf([&] { buildHelix("AAAA", options); });
    ASSERT_TRUE(err);
    EXPECT_EQ(err->first, "BuildError");
    options.pairing = "WC";
    EXPECT_NO_THROW(buildHelix("AAAA", options));
}

// --- geometry (tests/test_builder.py) ------------------------------------------------------------

namespace {

std::vector<double> o3pDistances(const io::PdbStructure& s) {
    std::vector<double> out;
    for (const auto& seg : s.segments) {
        for (std::size_t k = 0; k + 1 < seg.size(); ++k) {
            out.push_back(dist(coordsOf(seg[k]).at("O3'"), coordsOf(seg[k + 1]).at("P")));
        }
    }
    return out;
}

}  // namespace

class IdealDuplex : public ::testing::TestWithParam<const char*> {};

TEST_P(IdealDuplex, Geometry) {
    const std::string seq = "GGCAUCGAUGCCAG";
    HelixOptions options;
    options.chirality = GetParam();
    std::vector<std::string> notes;
    const io::PdbStructure s = io::readPdbText(buildHelix(seq, options, &notes));
    ASSERT_EQ(s.segments.size(), 2u);
    ASSERT_EQ(s.segments[0].size(), 14u);
    ASSERT_EQ(s.segments[1].size(), 14u);
    EXPECT_EQ(s.chainIds(), (std::vector<std::string>{"A", "B"}));
    // O3'(i)-P(i+1) on both strands: the Arnott links, 1.59 Angstrom
    const std::vector<double> links = o3pDistances(s);
    ASSERT_EQ(links.size(), 26u);
    for (const double d : links) {
        EXPECT_NEAR(d, 1.59, 0.006);
    }
    // the Watson-Crick pairs: C1'-C1' 10.9 Angstrom, donor-acceptor 2.7-3.1 Angstrom
    const std::size_t n = seq.size();
    for (std::size_t k = 0; k < n; ++k) {
        const io::PdbResidue& a = s.segments[0][k];
        const io::PdbResidue& b = s.segments[1][n - 1 - k];
        const auto ca = coordsOf(a);
        const auto cb = coordsOf(b);
        EXPECT_NEAR(dist(ca.at("C1'"), cb.at("C1'")), 10.9, 0.06) << a.label();
        const char b1 = io::parseResidueName(a.name)->base;
        const char b2 = io::parseResidueName(b.name)->base;
        const auto hbonds = pairHbonds("wc", b1, b2);
        ASSERT_FALSE(hbonds.empty()) << a.label();
        for (const PairHbond& h : hbonds) {
            const double d = dist(ca.at(std::string(h.atom1)), cb.at(std::string(h.atom2)));
            EXPECT_GT(d, 2.6) << a.label() << " " << h.atom1;
            EXPECT_LT(d, 3.2) << a.label() << " " << h.atom1;
        }
    }
    // the sugars' hand: the signs of the L-ribose chiral volumes (all D or all L)
    const double sign = options.chirality == "L" ? 1.0 : -1.0;
    for (const auto& seg : s.segments) {
        for (const io::PdbResidue& r : seg) {
            const auto c = coordsOf(r);
            const char base = io::parseResidueName(r.name)->base;
            for (const reference::ChiralCentre& centre : reference::lRiboseChirality()) {
                const auto atom = [&](std::string_view name) {
                    return c.at(name == "N" ? (base == 'A' || base == 'G' ? "N9" : "N1") : std::string(name));
                };
                if (!c.count(std::string(centre.centre)) || !c.count(std::string(centre.a)) ||
                    (centre.b != "N" && !c.count(std::string(centre.b))) || !c.count(std::string(centre.c))) {
                    continue;
                }
                const double v = tripleVolume(atom(centre.centre), atom(centre.a), atom(centre.b), atom(centre.c));
                EXPECT_GT(v * centre.sign * sign, 0.5) << r.label() << " " << centre.centre;
            }
        }
    }
    std::set<std::string> names;
    for (const io::PdbResidue* r : s.residues()) {
        names.insert(r->name);
    }
    EXPECT_EQ(names, (options.chirality == "D" ? std::set<std::string>{"A", "C", "G", "U"}
                                               : std::set<std::string>{"0A", "0C", "0G", "0U"}));
    EXPECT_EQ(s.header[0].rfind("REMARK 951 MDBridge IDEAL A-FORM RNA", 0), 0u);
    EXPECT_EQ(notes.back().rfind("Built an ideal A-form RNA duplex", 0), 0u);
}

INSTANTIATE_TEST_SUITE_P(Hands, IdealDuplex, ::testing::Values("D", "L"));

TEST(HelixGeometry, LHelixIsTheMirrorImageAndLeftHanded) {
    HelixOptions left;
    left.chirality = "L";
    left.hydrogens = true;
    HelixOptions right;
    right.hydrogens = true;
    const io::PdbStructure d = io::readPdbText(buildHelix("GGCAUC", right));
    const io::PdbStructure l = io::readPdbText(buildHelix("GGCAUC", left));
    const auto rd = d.residues();
    const auto rl = l.residues();
    ASSERT_EQ(rd.size(), rl.size());
    for (std::size_t k = 0; k < rd.size(); ++k) {
        const auto cl = coordsOf(*rl[k]);
        for (const auto& [name, x] : coordsOf(*rd[k])) {
            ASSERT_TRUE(cl.count(name)) << name;
            EXPECT_NEAR(-x[0], cl.at(name)[0], 2e-3) << name;
            EXPECT_NEAR(x[1], cl.at(name)[1], 2e-3) << name;
            EXPECT_NEAR(x[2], cl.at(name)[2], 2e-3) << name;
        }
    }
    // strand 1 turns by -32.7 degrees per step about +z going down for D, the opposite for L
    const auto turn = [](const io::PdbStructure& s) {
        const Vec p0 = coordsOf(s.segments[0][0]).at("C1'");
        const Vec p1 = coordsOf(s.segments[0][1]).at("C1'");
        const double pi = 3.14159265358979323846;
        double t = std::atan2(p1[1], p1[0]) - std::atan2(p0[1], p0[0]);
        t = std::fmod(t + pi + 2 * pi, 2 * pi) - pi;
        return std::make_pair(t * 180.0 / pi, p1[2] - p0[2]);
    };
    const auto [td, zd] = turn(d);
    const auto [tl, zl] = turn(l);
    EXPECT_NEAR(td, -kHelixTwist, 0.05);
    EXPECT_NEAR(zd, -kHelixRise, 0.01);
    EXPECT_NEAR(tl, kHelixTwist, 0.05);
    EXPECT_NEAR(zl, -kHelixRise, 0.01);
}

TEST(HelixGeometry, SugarsAreC3Endo) {
    for (const char* base : {"A", "C", "G", "U"}) {
        const auto t = templateMap(base, false);
        const char* n = (base[0] == 'A' || base[0] == 'G') ? "N9" : "N1";
        const char* c = (base[0] == 'A' || base[0] == 'G') ? "C4" : "C2";
        const double delta = dihedral(t.at("C5'"), t.at("C4'"), t.at("C3'"), t.at("O3'"));
        const double nu2 = dihedral(t.at("C1'"), t.at("C2'"), t.at("C3'"), t.at("C4'"));
        const double chi = dihedral(t.at("O4'"), t.at("C1'"), t.at(n), t.at(c));
        EXPECT_GT(delta, 75) << base;
        EXPECT_LT(delta, 90) << base;
        EXPECT_GT(nu2, 30) << base;
        EXPECT_GT(chi, -175) << base;
        EXPECT_LT(chi, -150) << base;
    }
}

TEST(HelixGeometry, NamesTerminalsHydrogensAndPhosphate) {
    const auto names = [](const HelixOptions& options, const std::string& seq) {
        std::vector<std::string> out;
        const io::PdbStructure built = io::readPdbText(buildHelix(seq, options));
        for (const io::PdbResidue* r : built.residues()) {
            out.push_back(r->name);
        }
        return out;
    };
    HelixOptions amber;
    amber.scheme = "amber";
    EXPECT_EQ(names(amber, "GCA"), (std::vector<std::string>{"G5", "C", "A3", "U5", "G", "C3"}));
    amber.chirality = "L";
    EXPECT_EQ(names(amber, "GCA"), (std::vector<std::string>{"LG5", "LC", "LA3", "LU5", "LG", "LC3"}));
    HelixOptions ccd;
    ccd.scheme = "ccd";
    ccd.chirality = "L";
    ccd.phosphate5 = true;
    EXPECT_EQ(names(ccd, "GGAC")[0], "0G");  // a 5' end with its phosphate is the CCD component
    HelixOptions pdb;
    pdb.scheme = "pdb";
    pdb.chirality = "L";
    for (const std::string& line : lines(buildHelix("GCA", pdb))) {
        if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) {
            EXPECT_EQ(line.substr(0, 6), "HETATM");
        }
    }
    HelixOptions phosphate;
    phosphate.phosphate5 = true;
    const io::PdbStructure p = io::readPdbText(buildHelix("GCA", phosphate));
    const auto& first = p.segments[0][0];
    ASSERT_GE(first.atoms.size(), 4u);
    EXPECT_EQ(first.atoms[0].name, "OP3");
    EXPECT_EQ(first.atoms[1].name, "P");
    EXPECT_NEAR(dist(coordsOf(first).at("OP3"), coordsOf(first).at("P")), 1.48, 0.01);
    EXPECT_FALSE(coordsOf(io::readPdbText(buildHelix("GCA")).segments[0][0]).count("P"));  // 5'-OH by default

    HelixOptions withH;
    withH.hydrogens = true;
    const io::PdbStructure s = io::readPdbText(buildHelix("GCAU", withH));
    const auto res = s.residues();
    const auto g5 = coordsOf(*res[0]);
    for (const char* h : {"H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", "H1", "H21", "H22", "H8", "HO5'"}) {
        EXPECT_TRUE(g5.count(h)) << h;
    }
    EXPECT_TRUE(coordsOf(*res[3]).count("HO3'"));
    for (const io::PdbResidue* r : res) {  // every hydrogen 0.9-1.12 Angstrom from its nearest heavy atom
        const auto c = coordsOf(*r);
        for (const auto& [name, x] : c) {
            if (name[0] != 'H') {
                continue;
            }
            double nearest = 1e9;
            for (const auto& [other, y] : c) {
                if (other[0] != 'H') {
                    nearest = std::min(nearest, dist(x, y));
                }
            }
            EXPECT_GT(nearest, 0.9) << r->label() << " " << name;
            EXPECT_LT(nearest, 1.12) << r->label() << " " << name;
        }
    }
    // PDB v3 prochiral names on the D sugar: positive volume C5'(O5', C4', H5')
    EXPECT_GT(tripleVolume(g5.at("C5'"), g5.at("O5'"), g5.at("C4'"), g5.at("H5'")), 0.0);
    // CHD-01: guanine's H21 trans to N1 (N3 side), H22 cis to N1 (Watson-Crick edge)
    EXPECT_LT(dist(g5.at("H22"), g5.at("N1")), dist(g5.at("H21"), g5.at("N1")));
    EXPECT_LT(dist(g5.at("H21"), g5.at("N3")), dist(g5.at("H22"), g5.at("N3")));
    // the H22 of G1 donates to O2 of its partner C (B:4)
    const auto c4 = coordsOf(*res[7]);
    EXPECT_NEAR(dist(g5.at("H22"), c4.at("O2")), 1.9, 0.25);
}

TEST(HelixGeometry, SingleStrandPartnerAndValidationMessages) {
    HelixOptions single;
    single.partner = "none";
    single.chains = "C";
    single.start = "10";
    const io::PdbStructure s = io::readPdbText(buildHelix("GGACU", single));
    ASSERT_EQ(s.segments.size(), 1u);
    std::vector<std::int64_t> numbers;
    for (const io::PdbResidue* r : s.residues()) {
        numbers.push_back(r->resseq);
    }
    EXPECT_EQ(numbers, (std::vector<std::int64_t>{10, 11, 12, 13, 14}));
    EXPECT_EQ(s.chainIds(), (std::vector<std::string>{"C"}));
    EXPECT_EQ(parseSequence("5'-GGT ACU-3'").sequence, "GGUACU");
    const auto msgs = checkHelix("GGCAU", "AUGCC");
    EXPECT_NE(std::find(msgs.begin(), msgs.end(), HelixMessage{"info", "the explicit partner is the Watson-Crick complement"}),
              msgs.end());
    EXPECT_EQ(checkHelix("GGC", "GG")[0].level, "error");
    EXPECT_EQ(checkHelix("GGCAU", "")[0].level, "error");
}

TEST(HelixGeometry, NumbersOutsideThePdbColumns) {
    // CH-01: below -999 an error before anything is written; CH-02: above 9999 a warning
    HelixOptions low;
    low.start = "1,-1000";
    const auto err = errorOf([&] { buildHelix("GGC", low); });
    ASSERT_TRUE(err);
    EXPECT_EQ(err->first, "PDBError");
    EXPECT_EQ(err->second, "residue number -1000 (B) does not fit the PDB columns (-999 to 9999)");
    HelixOptions high;
    high.start = "9997,1";
    std::vector<std::string> notes;
    const io::PdbStructure s = io::readPdbText(buildHelix("GGCA", high, &notes));
    EXPECT_EQ(notes.back(), "WARNING: residue numbers above 9999 wrap in the PDB columns");
    EXPECT_EQ(s.segments[0].back().resseq, 0);
    high.start = "9996,1";
    notes.clear();
    buildHelix("GGCA", high, &notes);
    EXPECT_EQ(notes.back().rfind("Built", 0), 0u);
}

TEST(HelixGeometry, ExplicitPartnerWithWobbleAndMismatch) {
    HelixOptions options;
    options.partner = "AAGCU";
    std::vector<std::string> notes;
    const io::PdbStructure s = io::readPdbText(buildHelix("GGCAU", options, &notes));
    bool wobble = false;
    for (const std::string& n : notes) {
        wobble = wobble || n.find("wobble geometry") != std::string::npos;
    }
    EXPECT_TRUE(wobble);
    for (const double d : o3pDistances(s)) {  // the wobble turns leave the backbone alone
        EXPECT_NEAR(d, 1.59, 0.006);
    }
    // G1·U wobble: G O6 - U N3 and G N1 - U O2 near 2.8-3.1 Angstrom
    const auto g = coordsOf(s.segments[0][0]);
    const auto u = coordsOf(s.segments[1][4]);
    EXPECT_NEAR(dist(g.at("O6"), u.at("N3")), 2.85, 0.15);
    EXPECT_NEAR(dist(g.at("N1"), u.at("O2")), 2.95, 0.2);
    options.partner = "CCAC";
    notes.clear();
    buildHelix("GGGG", options, &notes);
    bool collide = false;
    for (const std::string& n : notes) {
        collide = collide || n.rfind("WARNING: mismatch", 0) == 0;
    }
    EXPECT_TRUE(collide);
}
