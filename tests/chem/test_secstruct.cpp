// Secondary structure (src/chem/secstruct.cpp) against what MDATK's mdatk/secstruct.py computes
// (tests/parity/make_secstruct_cases.py, data/chem, see data/README.md), and MDATK's own tests of
// tests/test_secstruct.py and the builtin cases of tests/test_basepairs.py. Everything exact: the
// detector's scores and deviations bit for bit.

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::chem;
using core::Array;
using core::Value;
using test::kData;

namespace {

const Value& cases() {
    static const Value data =
        core::parseJson(test::readFile(kData / "chem" / "secstruct_cases.json"), "secstruct_cases.json");
    return data;
}

const Array& list(const Value& value) {
    return *value.get<Array>();
}
const Array& items(std::string_view key) {
    return list(*cases().find(key));
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
std::int64_t integer(const Value& value) {
    return *value.get<std::int64_t>();
}

// float.hex() text, "nan", "inf", "-inf"
double hexDouble(const Value& value) {
    const std::string text = str(value);
    if (text == "nan") {
        return std::nan("");
    }
    if (text == "inf") {
        return HUGE_VAL;
    }
    if (text == "-inf") {
        return -HUGE_VAL;
    }
    return std::strtod(text.c_str(), nullptr);
}

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

std::vector<ResiduePair> pairs(const Value& value) {
    std::vector<ResiduePair> out;
    for (const Value& p : list(value)) {
        out.emplace_back(integer(list(p).at(0)), integer(list(p).at(1)));
    }
    return out;
}

std::vector<std::int64_t> ints(const Value& value) {
    std::vector<std::int64_t> out;
    for (const Value& v : list(value)) {
        out.push_back(integer(v));
    }
    return out;
}

std::vector<BasePair> basePairs(const Value& value) {
    std::vector<BasePair> out;
    for (const Value& p : list(value)) {
        const Array& a = list(p);
        out.push_back({integer(a.at(0)), integer(a.at(1)), static_cast<int>(integer(a.at(2)))});
    }
    return out;
}

std::vector<PairCandidate> candidateRows(const Value& value) {
    std::vector<PairCandidate> out;
    for (const Value& r : list(value)) {
        const Array& a = list(r);
        out.push_back({integer(a.at(0)), integer(a.at(1)), static_cast<int>(integer(a.at(2))), hexDouble(a.at(3)),
                       hexDouble(a.at(4))});
    }
    return out;
}

void expectCandidates(const std::vector<PairCandidate>& got, const std::vector<PairCandidate>& expected,
                      const std::string& context) {
    ASSERT_EQ(got.size(), expected.size()) << context;
    for (std::size_t k = 0; k < got.size(); ++k) {
        const std::string where = context + " candidate " + std::to_string(k);
        EXPECT_EQ(got[k].first, expected[k].first) << where;
        EXPECT_EQ(got[k].second, expected[k].second) << where;
        EXPECT_EQ(got[k].kind, expected[k].kind) << where;
        EXPECT_TRUE(sameDouble(got[k].score, expected[k].score))
            << where << ": score " << got[k].score << " vs " << expected[k].score;
        EXPECT_TRUE(sameDouble(got[k].c1Deviation, expected[k].c1Deviation))
            << where << ": dev " << got[k].c1Deviation << " vs " << expected[k].c1Deviation;
    }
}

std::string pairsText(const std::vector<BasePair>& p) {
    std::string out;
    for (const BasePair& b : p) {
        out += "(" + std::to_string(b.i) + "," + std::to_string(b.j) + "," + std::to_string(b.kind) + ")";
    }
    return out;
}

template <class F>
std::string errorOf(F&& f) {
    try {
        f();
    } catch (const core::ValueError& e) {
        return e.what();
    }
    return "<no error>";
}

}  // namespace

// --- data-driven cases ------------------------------------------------------------------------------

TEST(Secstruct, Constants) {
    const Value& c = at(cases(), "constants");
    EXPECT_EQ(str(at(c, "openers")), kOpeners);
    EXPECT_EQ(str(at(c, "closers")), kClosers);
    EXPECT_EQ(str(at(c, "unpaired")), kUnpaired);
    EXPECT_EQ(str(at(c, "chain_separator")), std::string(1, kChainSeparator));
    const Array& rules = list(at(c, "pair_rules"));
    ASSERT_EQ(rules.size(), pairRules().size());
    for (std::size_t k = 0; k < rules.size(); ++k) {
        const PairRule& rule = pairRules()[k];
        EXPECT_EQ(str(at(rules[k], "base")), std::string(1, rule.base));
        std::string partners;
        for (const Value& p : list(at(rules[k], "partners"))) {
            partners += str(p);
        }
        EXPECT_EQ(partners, rule.partners);
        EXPECT_EQ(integer(at(rules[k], "kind")), rule.kind);
        EXPECT_EQ(integer(at(rules[k], "required")), rule.required);
        const Array& atoms = list(at(rules[k], "atoms"));
        ASSERT_EQ(atoms.size(), rule.atoms.size());
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            EXPECT_EQ(str(list(atoms[a]).at(0)), rule.atoms[a].first);
            EXPECT_EQ(str(list(atoms[a]).at(1)), rule.atoms[a].second);
        }
    }
    EXPECT_EQ(hexDouble(list(at(c, "c1_window")).at(0)), kC1WindowLow);
    EXPECT_EQ(hexDouble(list(at(c, "c1_window")).at(1)), kC1WindowHigh);
    EXPECT_EQ(hexDouble(at(c, "c1_ideal")), kC1Ideal);
    EXPECT_EQ(hexDouble(at(c, "score_tie")), kScoreTie);
    std::vector<std::string> atoms;
    for (const Value& a : list(at(c, "base_pair_atoms"))) {
        atoms.push_back(str(a));
    }
    EXPECT_EQ(atoms, basePairAtoms());
    for (const auto& [kind, label] : *at(c, "kind_labels").get<core::Table>()) {
        EXPECT_EQ(kindLabel(std::stoi(kind)), str(label));
    }
    EXPECT_EQ(integer(at(at(c, "kinds"), "WATSON_CRICK")), kWatsonCrick);
    EXPECT_EQ(integer(at(at(c, "kinds"), "WOBBLE")), kWobble);
    EXPECT_EQ(integer(at(at(c, "kinds"), "OTHER")), kOther);
}

TEST(Secstruct, ResidueNames) {
    for (const Value& c : items("residue_names")) {
        const std::string name = str(at(c, "name"));
        const Value& base = at(c, "base");
        const auto got = baseLetter(name);
        if (base.isNull()) {
            EXPECT_FALSE(got) << name;
        } else {
            ASSERT_TRUE(got) << name;
            EXPECT_EQ(std::string(1, *got), str(base)) << name;
        }
        const Value& canonical = at(c, "canonical");
        const auto gotName = canonicalName(name);
        if (canonical.isNull()) {
            EXPECT_FALSE(gotName) << name;
        } else {
            ASSERT_TRUE(gotName) << name;
            EXPECT_EQ(*gotName, str(canonical)) << name;
        }
    }
}

TEST(Secstruct, LwKind) {
    for (const Value& c : items("lw_kind")) {
        EXPECT_EQ(lwKind(str(at(c, "code"))), integer(at(c, "kind"))) << str(at(c, "code"));
    }
}

TEST(Secstruct, ParseCases) {
    for (const Value& c : items("parse")) {
        const std::string text = str(at(c, "text"));
        if (const Value* error = c.find("error")) {
            EXPECT_EQ(errorOf([&] { dotbracketToPairs(text); }), str(*error)) << text;
            continue;
        }
        DotBracket got;
        ASSERT_NO_THROW(got = dotbracketToPairs(text)) << text;
        EXPECT_EQ(got.pairs, pairs(at(c, "pairs"))) << text;
        EXPECT_EQ(got.n, integer(at(c, "n"))) << text;
        EXPECT_EQ(got.chainBreaks, ints(at(c, "breaks"))) << text;
    }
}

TEST(Secstruct, WriteCases) {
    for (const Value& c : items("write")) {
        const auto p = pairs(at(c, "pairs"));
        const std::int64_t n = integer(at(c, "n"));
        const std::vector<std::int64_t> breaks = c.find("breaks") ? ints(at(c, "breaks")) : std::vector<std::int64_t>{};
        const std::string context = core::toJson(c, -1);
        if (const Value* error = c.find("error")) {
            EXPECT_EQ(errorOf([&] { pairsToDotbracket(p, n, breaks); }), str(*error)) << context;
            continue;
        }
        std::string got;
        ASSERT_NO_THROW(got = pairsToDotbracket(p, n, breaks)) << context;
        EXPECT_EQ(got, str(at(c, "text"))) << context;
    }
}

TEST(Secstruct, StemsCrossesUnique) {
    for (const Value& c : items("stems")) {
        const auto got = stems(pairs(at(c, "pairs")));
        const Array& expected = list(at(c, "stems"));
        ASSERT_EQ(got.size(), expected.size()) << core::toJson(c, -1);
        for (std::size_t k = 0; k < got.size(); ++k) {
            EXPECT_EQ(got[k], pairs(expected[k])) << core::toJson(c, -1);
        }
    }
    for (const Value& c : items("crosses")) {
        const auto a = ints(at(c, "a"));
        const auto b = ints(at(c, "b"));
        EXPECT_EQ(crosses({a[0], a[1]}, {b[0], b[1]}), *at(c, "crosses").get<bool>()) << core::toJson(c, -1);
    }
    for (const Value& c : items("unique_pairs")) {
        const auto got = uniquePairs(pairs(at(c, "pairs")));
        EXPECT_EQ(got.keep, pairs(at(c, "keep"))) << core::toJson(c, -1);
        EXPECT_EQ(got.dropped, pairs(at(c, "dropped"))) << core::toJson(c, -1);
    }
}

TEST(Secstruct, DetectorCases) {
    for (const Value& c : items("detector")) {
        const std::string name = str(at(c, "name"));
        std::vector<std::string> resnames;
        for (const Value& r : list(at(c, "resnames"))) {
            resnames.push_back(str(r));
        }
        std::vector<AtomIndex> index;
        for (const Value& residue : list(at(c, "atoms"))) {
            AtomIndex atoms;
            for (const Value& a : list(residue)) {
                atoms[str(list(a).at(0))] = integer(list(a).at(1));
            }
            index.push_back(std::move(atoms));
        }
        const PairDetector det(resnames, index, hexDouble(at(c, "cutoff")), *at(c, "wobble").get<bool>());
        const Array& frames = list(at(c, "frames"));
        for (std::size_t f = 0; f < frames.size(); ++f) {
            const std::string context = name + " frame " + std::to_string(f);
            std::vector<double> xyz;
            for (const Value& v : list(at(frames[f], "xyz"))) {
                xyz.push_back(hexDouble(v));
            }
            expectCandidates(det.candidates(xyz), candidateRows(at(frames[f], "candidates")), context);
            const auto expected = basePairs(at(frames[f], "pairs"));
            const auto got = det.detect(xyz);
            EXPECT_EQ(got, expected) << context << ": " << pairsText(got) << " vs " << pairsText(expected);
        }
    }
}

TEST(Secstruct, SelectCases) {
    for (const Value& c : items("select")) {
        const std::string name = str(at(c, "name"));
        const auto cand = candidateRows(at(c, "candidates"));
        const Array& rounded = list(at(c, "rounded"));
        for (std::size_t k = 0; k < cand.size(); ++k) {
            const double r = std::nearbyint(cand[k].score / kScoreTie) * kScoreTie;
            EXPECT_TRUE(sameDouble(r, hexDouble(rounded[k]))) << name << " row " << k;
        }
        const auto expected = basePairs(at(c, "pairs"));
        const auto got = PairDetector::select(cand, static_cast<std::size_t>(integer(at(c, "n"))));
        EXPECT_EQ(got, expected) << name << ": " << pairsText(got) << " vs " << pairsText(expected);
    }
}

// The formula of make_secstruct_cases.py's lattice_xyz.
TEST(Secstruct, LatticeCases) {
    const std::vector<std::string> names = {"G", "C", "C", "C", "A", "U", "U", "DT", "GUA", "LC5"};
    const auto atomsOf = [](char base) -> std::vector<std::string> {
        switch (base) {
        case 'G':
            return {"N1", "N2", "O6", "C1'"};
        case 'C':
            return {"N3", "O2", "N4", "C1'"};
        case 'A':
            return {"N6", "N1", "C1'"};
        default:
            return {"O4", "N3", "O2", "C1'"};
        }
    };
    for (const Value& c : items("lattice")) {
        const std::int64_t nRes = integer(at(c, "n_residues"));
        std::vector<std::string> resnames;
        std::vector<AtomIndex> index;
        std::vector<double> xyz;
        for (std::int64_t r = 0; r < nRes; ++r) {
            const std::string& name = names[static_cast<std::size_t>(r) % names.size()];
            resnames.push_back(name);
            AtomIndex idx;
            const auto atoms = atomsOf(*baseLetter(name));
            for (std::int64_t a = 0; a < static_cast<std::int64_t>(atoms.size()); ++a) {
                idx[atoms[static_cast<std::size_t>(a)]] = static_cast<std::int64_t>(xyz.size() / 3);
                xyz.push_back(static_cast<double>((37 * r + 11 * a) % 53) / 4);
                xyz.push_back(static_cast<double>((17 * r + 5 * a) % 47) / 4);
                xyz.push_back(static_cast<double>((29 * r + 3 * a) % 41) / 4);
            }
            index.push_back(std::move(idx));
        }
        const PairDetector det(resnames, index, hexDouble(at(c, "cutoff")));
        const auto cand = det.candidates(xyz);
        const std::string context = "lattice " + std::to_string(nRes);
        ASSERT_EQ(static_cast<std::int64_t>(cand.size()), integer(at(c, "n_candidates"))) << context;
        const auto first = candidateRows(at(c, "first_candidates"));
        expectCandidates(std::vector<PairCandidate>(cand.begin(), cand.begin() + static_cast<std::ptrdiff_t>(first.size())),
                         first, context);
        // the candidate array as numpy holds it: float64 rows (i, j, kind, score, dev), little-endian
        std::string bytes;
        const auto put = [&](double v) {
            unsigned char b[8];
            std::uint64_t bits;
            std::memcpy(&bits, &v, 8);
            for (int k = 0; k < 8; ++k) {
                b[k] = static_cast<unsigned char>(bits >> (8 * k));
            }
            bytes.append(reinterpret_cast<const char*>(b), 8);
        };
        for (const PairCandidate& p : cand) {
            put(static_cast<double>(p.first));
            put(static_cast<double>(p.second));
            put(static_cast<double>(p.kind));
            put(p.score);
            put(p.c1Deviation);
        }
        EXPECT_EQ(core::sha256Hex(bytes), str(at(c, "candidates_sha256"))) << context;
        const auto detected = det.detect(xyz);
        ASSERT_EQ(static_cast<std::int64_t>(detected.size()), integer(at(c, "n_pairs"))) << context;
        std::string pairBytes;
        for (const BasePair& p : detected) {
            for (const std::int64_t v : {p.i, p.j, static_cast<std::int64_t>(p.kind)}) {
                const auto u = static_cast<std::uint64_t>(v);
                for (int k = 0; k < 8; ++k) {
                    pairBytes.push_back(static_cast<char>(static_cast<unsigned char>(u >> (8 * k))));
                }
            }
        }
        EXPECT_EQ(core::sha256Hex(pairBytes), str(at(c, "pairs_sha256"))) << context;
    }
}

// --- MDATK's tests/test_secstruct.py -----------------------------------------------------------------

TEST(Secstruct, ParseSimpleHairpin) {
    const DotBracket db = dotbracketToPairs("((..))");
    EXPECT_EQ(db.pairs, (std::vector<ResiduePair>{{0, 5}, {1, 4}}));
    EXPECT_EQ(db.n, 6);
    EXPECT_TRUE(db.chainBreaks.empty());
}

TEST(Secstruct, ParsePseudoknotLevelsAndLetters) {
    const DotBracket db = dotbracketToPairs("((..[[..))..]]{<Aa>}");
    EXPECT_EQ(db.n, 20);
    EXPECT_EQ(db.pairs, (std::vector<ResiduePair>{{0, 9}, {1, 8}, {4, 13}, {5, 12}, {14, 19}, {15, 18}, {16, 17}}));
}

TEST(Secstruct, ChainSeparators) {
    for (const char* text : {"(((&)))", "((( )))", "(((  &  )))", "&(((&)))&"}) {
        const DotBracket db = dotbracketToPairs(text);
        EXPECT_EQ(db.n, 6) << text;
        EXPECT_EQ(db.chainBreaks, std::vector<std::int64_t>{3}) << text;
        EXPECT_EQ(db.pairs, (std::vector<ResiduePair>{{0, 5}, {1, 4}, {2, 3}})) << text;
    }
}

TEST(Secstruct, UnbalancedMessages) {
    EXPECT_EQ(errorOf([] { dotbracketToPairs("..)"); }),
              "unbalanced dot-bracket: ')' at position 3 (residue 3) has no matching '('");
    EXPECT_EQ(errorOf([] { dotbracketToPairs("([.)"); }), "unbalanced dot-bracket: '[' at position 2 is never closed");
    EXPECT_EQ(errorOf([] { dotbracketToPairs("((."); }),
              "unbalanced dot-bracket: '(' at position 1 is never closed (and 1 more)");
    EXPECT_EQ(errorOf([] { dotbracketToPairs("(#)"); }), "invalid dot-bracket character '#' at position 2");
}

TEST(Secstruct, RoundTrip) {
    for (const char* text : {"((((....))))", "((..[[..))..]]", "(((&)))..((..))", "(([[{{<<..))]]}}>>", ".....", "(&)"}) {
        const DotBracket db = dotbracketToPairs(text);
        EXPECT_EQ(pairsToDotbracket(db.pairs, db.n, db.chainBreaks), text);
    }
}

TEST(Secstruct, RoundTripNormalisesLevels) {
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"((..[[..))..{{..]]..}}", "((..[[..))..((..]]..))"}, {"AB..((..ab..))", "[{..((..]}..))"}};
    for (const auto& [text, written] : cases) {
        const DotBracket db = dotbracketToPairs(text);
        EXPECT_EQ(pairsToDotbracket(db.pairs, db.n, db.chainBreaks), written);
        EXPECT_EQ(dotbracketToPairs(written), db);
    }
}

TEST(Secstruct, LongestHelixGetsRoundBrackets) {
    const std::vector<ResiduePair> p = {{0, 7}, {1, 6}, {3, 13}, {4, 12}, {5, 11}, {2, 14}};
    EXPECT_EQ(pairsToDotbracket(p, 15), "[[((((]]...))))");
}

TEST(Secstruct, PairsToDotbracketValidation) {
    EXPECT_EQ(pairsToDotbracket(std::vector<ResiduePair>{{5, 0}}, 6), "(....)");
    EXPECT_EQ(pairsToDotbracket(std::vector<ResiduePair>{{0, 5}, {0, 5}}, 6), "(....)");
    EXPECT_EQ(errorOf([] { pairsToDotbracket(std::vector<ResiduePair>{{0, 5}, {0, 4}}, 6); }),
              "residue 0 is in more than one base pair ((0, 5) and (0, 4))");
    EXPECT_EQ(errorOf([] { pairsToDotbracket(std::vector<ResiduePair>{{0, 6}}, 6); }),
              "invalid base pair (0, 6) for 6 residues");
    EXPECT_EQ(pairsToDotbracket({}, 4, std::vector<std::int64_t>{0, 2, 4}), "..&..");
}

TEST(Secstruct, UniquePairs) {
    const UniquePairs u = uniquePairs(std::vector<ResiduePair>{{0, 5}, {1, 5}, {2, 4}});
    EXPECT_EQ(u.keep, (std::vector<ResiduePair>{{0, 5}, {2, 4}}));
    EXPECT_EQ(u.dropped, (std::vector<ResiduePair>{{1, 5}}));
}

// --- MDATK's tests/test_basepairs.py (builtin detector; the geometries are in the data cases) -------

TEST(Secstruct, ScoreTiesGoToTheIdealC1Distance) {
    const std::vector<PairCandidate> cand = {{0, 1, kWatsonCrick, 5.63, 1.3}, {0, 2, kWatsonCrick, 5.63, 0.0}};
    EXPECT_EQ(PairDetector::select(cand, 3), (std::vector<BasePair>{{0, 2, kWatsonCrick}}));
}

TEST(Secstruct, NonNucleotidesHaveNoRules) {
    // a lanthanum ion named LA has no base atoms; ALA is no nucleotide
    const std::vector<std::string> names = {"A", "LA", "ALA"};
    const std::vector<AtomIndex> index = {{{"N6", 0}, {"N1", 1}}, {{"LA", 2}}, {{"N3", 3}, {"O4", 4}}};
    const PairDetector det(names, index);
    const std::vector<double> xyz(15, 0.0);
    EXPECT_TRUE(det.candidates(xyz).empty());
    EXPECT_TRUE(det.detect(xyz).empty());
}

TEST(Secstruct, RowsOutOfRange) {
    const std::vector<std::string> names = {"G", "C"};
    const std::vector<AtomIndex> index = {{{"N1", 0}, {"N2", 1}, {"O6", 2}}, {{"N3", 3}, {"O2", 4}, {"N4", 9}}};
    const PairDetector det(names, index);
    EXPECT_THROW(det.candidates(std::vector<double>(15, 0.0)), std::out_of_range);
    EXPECT_THROW(det.candidates(std::vector<double>(4, 0.0)), std::invalid_argument);
}
