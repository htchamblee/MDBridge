// The contacts and hbonds analyses (mdatk/analysis/contacts.py, hbonds.py) against what MDATK
// computes for every case of tests/parity/make_analysis_cases.py that runs them (data/analysis:
// micro_cases.json, demo_dry_cases.json, demo_dry_variants.json, systems_cases.json; see
// data/README.md). Each case is run as the generator runs it, restricted to its contacts / hbonds
// instances (and unknown analyses, whose error is part of the case): every array (dtype, shape,
// values) exactly, params, meta, output order, cache keys, the log lines of these analyses,
// progress, times and errors.
//
// A count that differs only where the generator's margin to a threshold is below kFlipMargin is
// reported as a possible last-bit flip; anywhere else it fails. Cases with use_pbc and a cell need
// the periodic neighbour search (W1.1): until it exists pairsWithin throws, and those cases are
// skipped ("needs periodic search").
//
// Also: the counterparts of MDATK's tests/test_analyses.py tests of these analyses, and the
// deliberate difference of the hydrogen guess on an exact distance tie (lowest donor index wins;
// the data's h_guess_tie margins show no tie, so the cases above do not depend on it).

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/spatial.hpp"

using namespace mdbridge;
using namespace mdbridge::test::analysis_cases;

namespace {

bool covered(const std::string& base) { return base == "contacts" || base == "hbonds"; }

const std::vector<std::string> kWarningPrefixes{"contacts", "H-bonds"};

struct CaseRef {
    std::string file;
    std::size_t index;
    std::string name;
};

std::vector<CaseRef> coveredCases() {
    std::vector<CaseRef> out;
    for (const std::string file : {"micro_cases.json", "demo_dry_cases.json", "demo_dry_variants.json", "systems_cases.json"}) {
        const Array& cases = list(at(caseFile(file), "cases"));
        for (std::size_t i = 0; i < cases.size(); ++i) {
            const std::vector<std::string> instances = coveredInstances(cases[i], covered);
            const bool any = std::any_of(instances.begin(), instances.end(), [](const std::string& instance) {
                return covered(analysis::parseInstance(instance).base);
            });
            if (any) {
                out.push_back({file, i, str(at(cases[i], "name"))});
            }
        }
    }
    return out;
}

std::string caseName(const ::testing::TestParamInfo<CaseRef>& info) {
    std::string prefix = info.param.file.substr(0, info.param.file.find("_cases"));
    prefix = prefix.substr(0, prefix.find(".json"));
    std::string name = prefix + "__" + info.param.name;
    for (char& ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch))) {
            ch = '_';
        }
    }
    return name;
}

class ContactsHBondsCase : public ::testing::TestWithParam<CaseRef> {};

TEST_P(ContactsHBondsCase, MatchesMdatk) {
    const CaseRef& ref = GetParam();
    const Value& file = caseFile(ref.file);
    const Value& c = list(at(file, "cases"))[ref.index];
    const std::vector<std::string> instances = coveredInstances(c, covered);
    const bool whole = instances.size() == list(at(c, "instances")).size();
    const CaseRun run = runCase(file, c, instances);

    const Value* expectedError = c.find("error");
    if (run.error && run.error->type == "C++" && run.error->message == "periodic neighbour search not implemented yet") {
        GTEST_SKIP() << "needs periodic search (use_pbc with a cell; pairsWithinPeriodic throws until W1.1)";
    }
    if (expectedError) {
        if (!whole) {
            GTEST_SKIP() << "the case's error comes from an analysis this test does not run";
        }
        ASSERT_TRUE(run.error.has_value()) << "expected " << str(at(*expectedError, "type"));
        EXPECT_EQ(run.error->type, str(at(*expectedError, "type")));
        EXPECT_EQ(run.error->message, str(at(*expectedError, "message")));
        EXPECT_EQ(run.logs, expectedLog(c, instances, kWarningPrefixes));
        return;
    }
    ASSERT_FALSE(run.error.has_value()) << run.error->type << ": " << run.error->message;

    EXPECT_EQ(run.logs, expectedLog(c, instances, kWarningPrefixes));
    const Array& progress = list(at(c, "progress"));
    ASSERT_EQ(run.progress.size(), progress.size());
    for (std::size_t n = 0; n < progress.size(); ++n) {
        EXPECT_EQ(run.progress[n].first, hexFloat(list(progress[n])[0]));
        EXPECT_EQ(run.progress[n].second, str(list(progress[n])[1]));
    }
    const Array& frames = list(at(c, "frames"));
    ASSERT_EQ(run.frames.size(), frames.size());
    for (std::size_t n = 0; n < frames.size(); ++n) {
        EXPECT_EQ(run.frames[n], integer(frames[n]));
    }
    const Array& times = list(at(c, "times"));
    ASSERT_EQ(run.run.timePs.size(), times.size());
    for (std::size_t n = 0; n < times.size(); ++n) {
        EXPECT_EQ(run.run.timePs[n], hexFloat(times[n])) << "time " << n;
    }

    ASSERT_EQ(run.run.results.size(), instances.size());
    for (std::size_t r = 0; r < instances.size(); ++r) {
        const core::AnalysisResult& got = run.run.results[r];
        const Value* expected = nullptr;
        for (const Value& e : list(at(c, "results"))) {
            if (str(at(e, "instance")) == instances[r]) {
                expected = &e;
            }
        }
        ASSERT_NE(expected, nullptr) << instances[r];
        SCOPED_TRACE(instances[r]);
        EXPECT_EQ(got.name, instances[r]);
        EXPECT_EQ(run.keys[r], str(at(*expected, "key")));
        expectSameValue(at(*expected, "params"), Value(got.params), "params");
        expectSameValue(at(*expected, "meta"), Value(got.meta), "meta");
        // output order (add_series / add_matrix / add_pairmap order)
        const Value meta(got.meta);
        const Value& outputs = at(meta, "outputs");
        for (const auto& [kind, names] : *at(*expected, "output_order").get<Table>()) {
            std::vector<std::string> order;
            if (const Value* entries = outputs.find(kind)) {
                for (const auto& [out, desc] : *entries->get<Table>()) {
                    order.push_back(out);
                }
            }
            EXPECT_EQ(order, strings(names)) << kind;
        }
        // arrays: the same names in the same order, every value exact
        const Table& arrays = *at(*expected, "arrays").get<Table>();
        std::vector<std::string> gotNames, wantNames;
        for (const auto& [name, a] : got.arrays) {
            gotNames.push_back(name);
        }
        for (const auto& [name, a] : arrays) {
            wantNames.push_back(name);
        }
        std::sort(gotNames.begin(), gotNames.end());
        std::sort(wantNames.begin(), wantNames.end());
        EXPECT_EQ(gotNames, wantNames);
        const double margin = smallestMargin(*expected);
        for (const auto& [name, want] : arrays) {
            const core::NumArray* a = got.array(name);
            if (!a) {
                continue;
            }
            const ArrayDiff diff = compareArray(want, *a);
            EXPECT_TRUE(diff.sameType) << name << ": dtype / shape differ (expected " << str(at(want, "dtype")) << " "
                                       << core::toJson(at(want, "shape"), -1) << ", got " << numpyDtype(*a) << " "
                                       << a->size() << " values)";
            if (diff.differing.empty()) {
                continue;
            }
            if (margin > 0.0 && margin < kFlipMargin) {
                std::cout << "[margin-explained] " << ref.file << " " << ref.name << " " << instances[r] << " " << name
                          << ": " << diff.differing.size() << " value(s) differ, smallest margin " << margin << "\n";
            } else {
                ADD_FAILURE() << name << ": " << diff.differing.size() << " value(s) differ, first at "
                              << diff.differing.front() << " (smallest margin " << margin << ")";
            }
        }
        // the cosine compared with: libm's cos(radians(angle_cutoff)), as the generator records it
        if (const Value* thresholds = at(*expected, "margins").find("thresholds")) {
            if (const Value* angle = thresholds->find("angle_cutoff")) {
                const double deg = hexFloat(at(*angle, "angle"));
                EXPECT_EQ(std::cos(deg * (std::numbers::pi / 180.0)), hexFloat(at(*angle, "libm")));
                if (str(at(*angle, "numpy")) != str(at(*angle, "libm"))) {
                    std::cout << "[numpy-cos] " << ref.name << ": numpy " << str(at(*angle, "numpy")) << " vs libm "
                              << str(at(*angle, "libm")) << "\n";
                }
            }
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Data, ContactsHBondsCase, ::testing::ValuesIn(coveredCases()), caseName);

std::vector<std::int64_t> int64s(const core::NumArray& a) {
    const auto v = a.values<std::int64_t>();
    return {v.begin(), v.end()};
}

// MDATK tests/test_analyses.py::test_contacts_counts_and_native_fraction
TEST(ContactsHBonds, ContactsCountsAndNativeFraction) {
    core::System u = makeMicroSystem({{"C1", 1, "A", "C"}, {"O1", 1, "A", "O"}, {"N1", 2, "U", "N"}, {"H1", 2, "U", "H"}},
                                     {{0, 0, 0, 1, 0, 0, 4, 0, 0, 4.5, 0, 0}, {0, 0, 0, 1, 0, 0, 10, 0, 0, 10.5, 0, 0}});
    core::Project p = microProject({{"Aptamer", "1"}, {"Target", "2"}});
    p.analyses.contacts.cutoff = 3.5;
    const core::AnalysisResult res = runOne("contacts", u, p);
    EXPECT_EQ(int64s(*res.array("atom_pairs")), (std::vector<std::int64_t>{1, 0}));  // O1-N1 only (H ignored)
    const auto q = res.array("q_native")->values<double>();
    EXPECT_EQ(std::vector<double>(q.begin(), q.end()), (std::vector<double>{1.0, 0.0}));
    EXPECT_EQ(res.array("per_res_a")->shape(), (std::vector<std::size_t>{1, 2}));
    ASSERT_EQ(res.meta[0].first, "labels_a");
    EXPECT_EQ(res.meta[0].second, Value(Array{"A1"}));
    ASSERT_EQ(res.meta[1].first, "labels_b");
    EXPECT_EQ(res.meta[1].second, Value(Array{"U2"}));
}

// MDATK tests/test_analyses.py::test_hbond_distance_and_angle (with and without bonds)
TEST(ContactsHBonds, HBondDistanceAndAngle) {
    for (const bool withBonds : {true, false}) {
        SCOPED_TRACE(withBonds ? "with bonds" : "without bonds");
        const std::vector<std::vector<double>> frames{{0, 0, 0, 1.0, 0, 0, 2.9, 0, 0},   // linear
                                                      {0, 0, 0, 1.0, 0, 0, 0, 2.9, 0},   // bent
                                                      {0, 0, 0, 1.0, 0, 0, 4.0, 0, 0}};  // far
        std::vector<std::pair<std::int32_t, std::int32_t>> bonds;
        if (withBonds) {
            bonds.emplace_back(0, 1);
        }
        core::System u = makeMicroSystem({{"N3", 1, "U", "N"}, {"H3", 1, "U", "H"}, {"O6", 2, "G", "O"}}, frames, bonds);
        core::Project p = microProject({{"Aptamer", "1"}, {"Target", "2"}});
        EXPECT_EQ(int64s(*runOne("hbonds", u, p).array("count")), (std::vector<std::int64_t>{1, 0, 0}));
        p.analyses.hbonds.use_angle = false;
        EXPECT_EQ(int64s(*runOne("hbonds", u, p).array("count")), (std::vector<std::int64_t>{1, 1, 0}));
    }
}

// Deliberate difference (hbonds.cpp): without bonds, a hydrogen exactly as far from two candidate
// donors goes to the one with the lower index (MDATK: the first pair in scipy's tree order, which
// MDBridge does not reproduce). H1 lies 1 Å from N1 (index 0) and N2 (index 2), so only N1 donates:
// to O1 when O1 is 2.5 Å from N1, to nothing when O1 is 2.5 Å from N2.
TEST(ContactsHBonds, HydrogenGuessTieGoesToLowestDonorIndex) {
    for (const bool nearN2 : {false, true}) {
        SCOPED_TRACE(nearN2 ? "O1 near N2" : "O1 near N1");
        const double o = nearN2 ? 3.5 : -3.5;
        core::System u = makeMicroSystem(
            {{"N1", 1, "U", "N"}, {"H1", 1, "U", "H"}, {"N2", 1, "U", "N"}, {"O1", 2, "G", "O"}},
            {{-1, 0, 0, 0, 0, 0, 1, 0, 0, o, 0, 0}});
        core::Project p = microProject({{"Aptamer", "1"}, {"Target", "2"}});
        p.analyses.hbonds.use_angle = false;
        EXPECT_EQ(int64s(*runOne("hbonds", u, p).array("count")), (std::vector<std::int64_t>{nearN2 ? 0 : 1}));
    }
}

// The contacts / H-bond part of MDATK tests/test_analyses.py::test_demo_system_signals
TEST(ContactsHBonds, DemoSystemSignals) {
    core::System u = io::loadSystem(test::kData / "topology" / "demo_dry.parm7",
                                    {test::kData / "io" / "trajectories" / "demo_dry.nc"});
    const std::string fileName = "demo_dry_cases.json";
    const Value& file = caseFile(fileName);
    const core::Project p = core::projectFromPlain(untag(at(list(at(file, "cases"))[0], "project")));
    const std::vector<std::int64_t> hb = int64s(*runOne("hbonds", u, p).array("count"));
    ASSERT_EQ(hb.size(), 24u);
    // 6 A-U (2) + 6 G-C (3) Watson-Crick H-bonds, +/- thermal-noise extras
    EXPECT_GE(hb[0], 28);
    EXPECT_LE(hb[0], 34);
    const std::size_t mid = hb.size() / 2;
    EXPECT_LT(hb[mid], hb[0]);
    EXPECT_LT(runOne("contacts", u, p).array("q_native")->values<double>()[mid], 1.0);
}

}  // namespace
