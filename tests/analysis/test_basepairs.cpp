// The basepairs analysis (mdatk/analysis/basepairs.py) against what MDATK computes for every case of
// tests/parity/make_analysis_cases.py that runs it (data/analysis: micro_cases.json — the duplex
// with strides 1 and 2, wobble off / cutoff 3.0, structure residues only, L-RNA names, wobble
// default / off / cutoff 4 —, demo_dry_cases.json — default stride 10 —, demo_dry_variants.json —
// readers, use_pbc, frame slices, stride 1 / 7, no wobble, cutoff 3.0, no ions, Target only —, and
// systems_cases.json — periodic cells open and with use_pbc, PDB / GRO without bonds, heavy atoms
// only, aromatic residues). Each case is run as the generator runs it, restricted to its basepairs
// instances: every array (dtype, shape, values) exactly, params, meta, output order, cache keys,
// the log lines of the analysis, progress, times and errors. Base pairs ignore the cell, so the
// periodic cases need no periodic neighbour search.
//
// A count that differs only where the generator's margin to a threshold (H-bond cutoff, C1'
// window, distance of the score to the 0.05 grid boundary, score-tie C1' deviation) is below
// kFlipMargin is reported as a possible last-bit flip; anywhere else it fails.
//
// Also: the counterparts of the analysis tests of MDATK's tests/test_basepairs.py (builtin engine),
// the error without structure residues, the stride clamp, and the deliberate difference for
// engine = "barnaba" (warning, then the builtin detector).

#include <algorithm>
#include <cctype>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/core/selections.hpp"

using namespace mdbridge;
using namespace mdbridge::test::analysis_cases;

namespace {

bool covered(const std::string& base) { return base == "basepairs"; }

const std::vector<std::string> kWarningPrefixes{"basepairs", "base pairs"};

const std::string kMicro = "micro_cases.json";
const std::string kDemo = "demo_dry_cases.json";
const std::vector<std::string> kFiles{kMicro, kDemo, "demo_dry_variants.json", "systems_cases.json"};

struct CaseRef {
    std::string file;
    std::size_t index;
    std::string name;
};

std::vector<CaseRef> coveredCases() {
    std::vector<CaseRef> out;
    for (const std::string& file : kFiles) {
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

// (names as const char*: no temporary arguments for a returned reference)
const Value& findCase(const std::string& file, const char* name) {
    for (const Value& c : list(at(caseFile(file), "cases"))) {
        if (str(at(c, "name")) == name) {
            return c;
        }
    }
    throw std::runtime_error("no case " + std::string(name));
}

const Value& expectedResult(const Value& c, const char* instance) {
    for (const Value& e : list(at(c, "results"))) {
        if (str(at(e, "instance")) == instance) {
            return e;
        }
    }
    throw std::runtime_error("no result " + std::string(instance));
}

class BasePairsCase : public ::testing::TestWithParam<CaseRef> {};

TEST_P(BasePairsCase, MatchesMdatk) {
    const CaseRef& ref = GetParam();
    const Value& file = caseFile(ref.file);
    const Value& c = list(at(file, "cases"))[ref.index];
    const std::vector<std::string> instances = coveredInstances(c, covered);
    const bool whole = instances.size() == list(at(c, "instances")).size();
    const CaseRun run = runCase(file, c, instances);

    if (const Value* expectedError = c.find("error")) {
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
        const Value& expected = expectedResult(c, instances[r].c_str());
        SCOPED_TRACE(instances[r]);
        EXPECT_EQ(got.name, instances[r]);
        EXPECT_EQ(run.keys[r], str(at(expected, "key")));
        expectSameValue(at(expected, "params"), Value(got.params), "params");
        expectSameValue(at(expected, "meta"), Value(got.meta), "meta");
        const Value meta(got.meta);
        const Value& outputs = at(meta, "outputs");
        for (const auto& [kind, names] : *at(expected, "output_order").get<Table>()) {
            std::vector<std::string> order;
            if (const Value* entries = outputs.find(kind)) {
                for (const auto& [out, desc] : *entries->get<Table>()) {
                    order.push_back(out);
                }
            }
            EXPECT_EQ(order, strings(names)) << kind;
        }
        // arrays: the same names (MDATK's insertion order: map_*, count, resindex), values exact
        const Table& arrays = *at(expected, "arrays").get<Table>();
        std::vector<std::string> gotNames, wantNames;
        for (const auto& [name, a] : got.arrays) {
            gotNames.push_back(name);
        }
        EXPECT_EQ(gotNames, (std::vector<std::string>{"map_offsets", "map_i", "map_j", "map_v", "count", "resindex"}));
        for (const auto& [name, a] : arrays) {
            wantNames.push_back(name);
        }
        std::sort(gotNames.begin(), gotNames.end());
        EXPECT_EQ(gotNames, wantNames);
        const double margin = smallestMargin(expected);
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
    }
}

INSTANTIATE_TEST_SUITE_P(Data, BasePairsCase, ::testing::ValuesIn(coveredCases()), caseName);

// Every group of cases has at least one basepairs case (none silently dropped from the data).
TEST(BasePairs, EveryCaseFileCovered) {
    const std::vector<CaseRef> cases = coveredCases();
    for (const std::string& file : kFiles) {
        const auto n = std::count_if(cases.begin(), cases.end(), [&](const CaseRef& c) { return c.file == file; });
        EXPECT_GT(n, 0) << file;
    }
    EXPECT_EQ(cases.size(), 33u);  // micro 10, demo 1, variants 12, systems 10
}

std::vector<std::int64_t> int64s(const core::NumArray& a) {
    const auto v = a.values<std::int64_t>();
    return {v.begin(), v.end()};
}

// The micro case `name` run with `edit` applied to its project, basepairs only.
struct MicroRun {
    std::vector<std::string> logs;
    core::AnalysisResult result;
    std::string key;
};

template <class Edit>
MicroRun runMicro(const char* name, Edit&& edit) {
    const Value& file = caseFile(kMicro);
    const Value& c = findCase(kMicro, name);
    core::Project project = core::projectFromPlain(untag(at(c, "project")));
    edit(project);
    core::System system = caseSystem(file, str(at(c, "system")), std::nullopt);
    MicroRun out;
    const core::Reporter reporter([&](const std::string& m) { out.logs.push_back(m); });
    analysis::AnalysisContext context(system, project, "", reporter);
    out.key = analysis::basepairsAnalysis().cacheKey(project, untag(at(c, "fingerprint")));
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    analyses.push_back(analysis::basepairsAnalysis().create(context));
    std::vector<std::int64_t> frames;
    for (const Value& f : list(at(c, "frames"))) {
        frames.push_back(integer(f));
    }
    out.result = std::move(analysis::runAnalyses(context, analyses, frames).results.front());
    return out;
}

// MDATK tests/test_basepairs.py::test_analysis_rows_align_with_trace_and_stride
TEST(BasePairs, RowsAlignWithTraceAndStride) {
    const Value& c = findCase(kMicro, "duplex_stride_1");
    const CaseRun run = runCase(caseFile(kMicro), c, {"basepairs", "trace"});
    ASSERT_FALSE(run.error.has_value()) << run.error->message;
    const core::AnalysisResult& bp = run.run.results[0];
    const core::AnalysisResult& trace = run.run.results[1];
    EXPECT_EQ(core::toJson(*Value(bp.meta).find("residue_labels"), -1), core::toJson(*Value(trace.meta).find("residue_labels"), -1));
    EXPECT_EQ(int64s(*bp.array("resindex")), int64s(*trace.array("resindex")));
    EXPECT_EQ(int64s(*bp.array("count")), (std::vector<std::int64_t>{2, 1, 1, 2, 2}));
    // frame 0: G1-C4 and A2-U3, Watson-Crick, i < j only
    const auto offsets = bp.array("map_offsets")->values<std::int64_t>();
    const auto i = bp.array("map_i")->values<std::int32_t>();
    const auto j = bp.array("map_j")->values<std::int32_t>();
    const auto v = bp.array("map_v")->values<float>();
    ASSERT_EQ(offsets[1], 2);
    EXPECT_EQ(std::vector<std::int32_t>(i.begin(), i.begin() + 2), (std::vector<std::int32_t>{0, 1}));
    EXPECT_EQ(std::vector<std::int32_t>(j.begin(), j.begin() + 2), (std::vector<std::int32_t>{3, 2}));
    EXPECT_EQ(v[0], 1.0f);
    EXPECT_EQ(v[1], 1.0f);
    // frame 1: only A2-U3
    ASSERT_EQ(offsets[2], 3);
    EXPECT_EQ(i[2], 1);
    EXPECT_EQ(*Value(bp.meta).find("engine")->get<std::string>(), "builtin");
    // stride 2: frames 1 and 3 reuse frames 0 and 2
    const MicroRun stride2 = runMicro("duplex_stride_1", [](core::Project& p) { p.analyses.basepairs.stride = 2; });
    EXPECT_EQ(int64s(*stride2.result.array("count")), (std::vector<std::int64_t>{2, 2, 1, 1, 2}));
}

// MDATK tests/test_basepairs.py::test_analysis_uses_structure_residues_only and
// test_l_rna_residues_pair (as recorded in the data)
TEST(BasePairs, StructureResiduesOnlyAndLRna) {
    const MicroRun only = runMicro("structure_residues_only", [](core::Project&) {});
    EXPECT_EQ(int64s(*only.result.array("count")), (std::vector<std::int64_t>{1, 1, 1, 1, 1}));
    EXPECT_EQ(core::toJson(*Value(only.result.meta).find("residue_labels"), -1),
              core::toJson(Value(Array{"G51", "A2", "U3", "ALA5"}), -1));
    const MicroRun lrna = runMicro("l_rna", [](core::Project&) {});
    EXPECT_EQ(int64s(*lrna.result.array("count")), (std::vector<std::int64_t>{1, 1}));
}

// A stride below 1 detects on every frame (max(1, stride)); meta keeps the clamped value.
TEST(BasePairs, StrideBelowOneIsOne) {
    const MicroRun one = runMicro("duplex_stride_1", [](core::Project&) {});
    for (const std::int64_t stride : {0, -3}) {
        const MicroRun r = runMicro("duplex_stride_1", [&](core::Project& p) { p.analyses.basepairs.stride = stride; });
        EXPECT_EQ(int64s(*r.result.array("count")), int64s(*one.result.array("count")));
        EXPECT_EQ(*Value(r.result.meta).find("stride")->get<std::int64_t>(), 1);
        EXPECT_EQ(*Value(r.result.params).find("stride")->get<std::int64_t>(), stride);
    }
}

// structure_residues without residues (every atom an ion): MDATK's SelectionError.
TEST(BasePairs, NoStructureResiduesIsAnError) {
    const Value& file = caseFile(kDemo);
    const Value& c = findCase(kDemo, "default");
    core::Project project = core::projectFromPlain(untag(at(c, "project")));
    project.groups.push_back({"Ions", "25-27", "", "#000000"});
    project.visual.structure_groups = {"Ions"};
    core::System system = caseSystem(file, str(at(c, "system")), std::nullopt);
    analysis::AnalysisContext context(system, project);
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    analyses.push_back(analysis::basepairsAnalysis().create(context));
    const std::vector<std::int64_t> frames{0};
    try {
        (void)analysis::runAnalyses(context, analyses, frames);
        ADD_FAILURE() << "expected SelectionError";
    } catch (const core::SelectionError& e) {
        EXPECT_EQ(std::string(e.what()), "base pairs: no residues selected for the structural view");
    }
}

// Deliberate difference: MDBridge does not bundle Barnaba (GPL). engine = "barnaba" warns and runs
// the builtin detector (same arrays as engine = "builtin"; meta engine "builtin", no
// lw_first_frame); params, and so the cache key, keep engine "barnaba" as MDATK's would.
TEST(BasePairs, BarnabaFallsBackToBuiltin) {
    const MicroRun builtin = runMicro("duplex_stride_1", [](core::Project&) {});
    const MicroRun barnaba = runMicro("duplex_stride_1", [](core::Project& p) {
        p.analyses.basepairs.engine = "barnaba";
        p.analyses.basepairs.noncanonical = true;
    });
    EXPECT_EQ(barnaba.logs.front(), "Preparing Base pairs ...");
    EXPECT_NE(std::find(barnaba.logs.begin(), barnaba.logs.end(),
                        "WARNING: basepairs: Barnaba is not available in MDBridge; using the builtin detector"),
              barnaba.logs.end());
    for (const std::string& line : builtin.logs) {
        EXPECT_EQ(line.rfind("WARNING", 0), std::string::npos) << line;
    }
    ASSERT_EQ(barnaba.result.arrays.size(), builtin.result.arrays.size());
    for (std::size_t n = 0; n < builtin.result.arrays.size(); ++n) {
        EXPECT_EQ(barnaba.result.arrays[n].first, builtin.result.arrays[n].first);
        EXPECT_EQ(barnaba.result.arrays[n].second.toDouble(), builtin.result.arrays[n].second.toDouble());
    }
    const Value meta(barnaba.result.meta);
    EXPECT_EQ(*meta.find("engine")->get<std::string>(), "builtin");
    EXPECT_EQ(meta.find("lw_first_frame"), nullptr);
    EXPECT_EQ(*Value(barnaba.result.params).find("engine")->get<std::string>(), "barnaba");
    EXPECT_EQ(*Value(barnaba.result.params).find("noncanonical")->get<bool>(), true);
    EXPECT_NE(barnaba.key, builtin.key);
    EXPECT_EQ(builtin.key, str(at(expectedResult(findCase(kMicro, "duplex_stride_1"), "basepairs"), "key")));
}

}  // namespace
