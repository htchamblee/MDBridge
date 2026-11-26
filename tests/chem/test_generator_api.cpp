// chem/generator_api.hpp, the functions milestone 7's generators call (docs/specs/m6-plan.md 4.4),
// against the generator cases of the builder (tests/data/chem/builder_cases.json: check_helix,
// resolve_partner, pair_types) and on files the builder and the annotation write.
#include <gtest/gtest.h>

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "mdbridge/chem/annotate.hpp"
#include "mdbridge/chem/generator_api.hpp"
#include "mdbridge/chem/helix.hpp"
#include "mdbridge/core/json.hpp"
#include "test_support.hpp"

namespace mdbridge {
namespace {

namespace fs = std::filesystem;
using core::Value;
namespace gen = chem::generator;

const Value& sequences() {
    static const Value data = test::readJson(fs::path(MDBRIDGE_TEST_DATA) / "chem" / "builder_cases.json");
    return *data.find("sequences");
}

std::optional<std::string> optionalText(const Value& v) {
    return v.isNull() ? std::nullopt : std::optional<std::string>(v.stringOr());
}

TEST(GeneratorApi, CheckHelixAsTheBuilderCases) {
    std::size_t n = 0;
    for (const Value& c : *sequences().find("check_helix")->get<core::Array>()) {
        const core::Array& a = *c.get<core::Array>();
        std::vector<gen::HelixMessage> expected;
        for (const Value& m : *a[5].get<core::Array>()) {
            const core::Array& pair = *m.get<core::Array>();
            expected.push_back({pair[0].stringOr(), pair[1].stringOr()});
        }
        const std::optional<std::string> partner = optionalText(a[1]), pairing = optionalText(a[2]),
                                         triplex = optionalText(a[3]), third = optionalText(a[4]);
        const auto view = [](const std::optional<std::string>& s) {
            return s ? std::optional<std::string_view>(*s) : std::nullopt;
        };
        EXPECT_EQ(gen::checkHelix(a[0].stringOr(), view(partner), view(pairing), view(triplex), view(third)), expected)
            << core::toJson(c, 0);
        ++n;
    }
    EXPECT_GT(n, 1000u);
}

TEST(GeneratorApi, ResolvePartnerAndPairTypesAsTheBuilderCases) {
    for (const Value& c : *sequences().find("resolve_partner")->get<core::Array>()) {
        const core::Array& a = *c.get<core::Array>();
        const Value* value = a[3].find("value");
        if (value == nullptr) {
            continue;  // the error cases: chem's own tests
        }
        const std::optional<std::string> partner = optionalText(a[1]), pairing = optionalText(a[2]);
        const gen::ResolvedPartner r =
            gen::resolvePartner(a[0].stringOr(), partner ? std::optional<std::string_view>(*partner) : std::nullopt,
                                pairing ? std::optional<std::string_view>(*pairing) : std::nullopt);
        const core::Array& v = *value->get<core::Array>();
        EXPECT_EQ(r.sequence, v[0].stringOr());
        EXPECT_EQ(r.partner, optionalText(v[1])) << core::toJson(c, 0);
    }
    for (const Value& c : *sequences().find("pair_types")->get<core::Array>()) {
        const core::Array& a = *c.get<core::Array>();
        std::vector<std::string> expected;
        for (const Value& k : *a[2].get<core::Array>()) expected.push_back(k.stringOr());
        EXPECT_EQ(gen::pairTypes(a[0].stringOr(), a[1].stringOr()), expected);
    }
}

TEST(GeneratorApi, ReadModelAssignsChiralityAndNames) {
    test::TempDir dir;
    chem::HelixOptions o;
    o.chirality = "L";
    test::writeFile(dir.path() / "l.pdb", chem::buildHelix("GGCAU", o));
    const gen::Model plain = gen::readModel(dir.path() / "l.pdb", false);
    ASSERT_EQ(plain.residues.size(), 10u);
    EXPECT_EQ(plain.residues[0].name, "0G");
    ASSERT_TRUE(plain.residues[0].nucleotide);
    EXPECT_EQ(plain.residues[0].nucleotide->chirality, 'L');
    EXPECT_TRUE(plain.chiralityNotes.empty());

    const gen::Model ccd = gen::readModel(dir.path() / "l.pdb", true, std::nullopt, std::string("ccd"));
    EXPECT_EQ(ccd.chiralityNotes, (chem::Notes{"Chirality (sugar geometry): 10 L and 0 D nucleotides"}));
    EXPECT_EQ(ccd.residues[0].name, "0G5");
    EXPECT_EQ(ccd.residues[0].chirality, "L");
    EXPECT_EQ(ccd.residues[0].label, "A:1 0G5");
    EXPECT_EQ(ccd.residues[0].atoms.front(), "O5'");
    const std::vector<std::string> variants = gen::codesInModel(ccd);
    for (const std::string& code : variants) {
        ASSERT_NE(chem::findComponent(code), nullptr);
        EXPECT_TRUE(gen::isVariant(*chem::findComponent(code)));
    }
    EXPECT_FALSE(variants.empty());

    // amber.py: the L ranges are checked at plan time
    EXPECT_EQ(gen::parseRanges(std::string("A:1-2, B")).size(), 2u);
    EXPECT_TRUE(gen::parseRanges(std::nullopt).empty());
    EXPECT_THROW(gen::parseRanges(std::string("A:4-1")), io::PdbError);  // CH-15
    EXPECT_THROW(gen::readModel(dir.path() / "missing.pdb"), fs::filesystem_error);
}

TEST(GeneratorApi, ReadAnnotationReadsWhatAnnotateWrites) {
    const chem::AnnotationResult res = chem::annotate(chem::PdbSource::text(chem::buildHelix("GGCAUGCC")));
    const gen::Annotation a = gen::readAnnotation(res.records);
    EXPECT_TRUE(a.annotated);
    ASSERT_EQ(a.pairs.size(), 8u);
    EXPECT_TRUE(a.stacks.empty());
    EXPECT_EQ(a.pairs[0].a, (gen::RemarkResidue{"A", 1, "", "G"}));
    EXPECT_EQ(a.pairs[0].b, (gen::RemarkResidue{"B", 8, "", "C"}));
    EXPECT_EQ(a.pairs[0].code, "WCc");
    EXPECT_NE(gen::lwKind(a.pairs[0].code), gen::lwKind("xx"));
    // unreadable records are skipped
    const std::vector<std::string> broken{"REMARK 950 PAIR  A  xx G   B   8 C   WCc", "REMARK 950 STACK"};
    EXPECT_TRUE(gen::readAnnotation(broken).pairs.empty());
    EXPECT_TRUE(gen::readAnnotation(broken).annotated);
}

TEST(GeneratorApi, DictionaryComponents) {
    EXPECT_EQ(gen::components().size(), 12u);
    EXPECT_EQ(gen::kDefaultGroup, "RNA");
    std::size_t variants = 0;
    for (const gen::Component& c : gen::components()) variants += gen::isVariant(c) ? 1 : 0;
    EXPECT_EQ(variants, 8u);
    EXPECT_NO_THROW(gen::checkSmiles(gen::components().front().code, gen::components().front().smiles));
    EXPECT_THROW(gen::parseSmiles("C1CC"), io::DictionaryError);
}

}  // namespace
}  // namespace mdbridge
