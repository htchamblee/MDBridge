// The mg and stacking analyses (mdatk/analysis/mg.py, stacking.py) against what MDATK computes for
// the cases of tests/parity/make_analysis_cases.py (data/analysis: micro systems, the demo project
// and its variants, the generated systems; see data/README.md). Every case that runs mg or
// stacking (also as instances) is rerun here with those analyses only, and compared: arrays
// (dtype, shape, values), params, meta, output order, cache key, warnings, times.
//
// Counts are exact, except that a frame whose decision margin (recorded by the generator: the
// smallest distance of a compared quantity from its threshold) is positive but below kFlipMargin
// may differ, as a last-bit flip of numpy / LAPACK against MDBridge's numerics; such frames are
// printed. Deliberate difference: within a frame, MDATK lists the stacking map's entries in
// scipy's cKDTree pair order and MDBridge in (i, j) order, so each frame's entries are compared as
// a sorted list. Cases with use_pbc are skipped while the periodic neighbour search is missing.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/project.hpp"

using namespace mdbridge;
using namespace mdbridge::analysis;
using namespace mdbridge::test::analysis_cases;

namespace {

const char* const kFiles[] = {"micro_cases.json", "demo_dry_cases.json", "demo_dry_variants.json",
                              "systems_cases.json"};

bool isMine(const std::string& instance) {
    const std::string base = parseInstance(instance).base;
    return base == "mg" || base == "stacking";
}

struct CaseRef {
    std::string file;
    std::size_t index = 0;
    std::string name;
};

void PrintTo(const CaseRef& ref, std::ostream* out) { *out << ref.file << " " << ref.name; }

std::vector<CaseRef> caseRefs() {
    std::vector<CaseRef> out;
    for (const std::string file : kFiles) {
        const Value& data = caseFile(file);
        const Array& cases = list(at(data, "cases"));
        for (std::size_t i = 0; i < cases.size(); ++i) {
            const std::vector<std::string> instances = strings(at(cases[i], "instances"));
            if (std::any_of(instances.begin(), instances.end(), isMine)) {
                out.push_back({file, i, str(at(cases[i], "name"))});
            }
        }
    }
    return out;
}

std::string testName(const testing::TestParamInfo<CaseRef>& info) {
    std::string out = info.param.file.substr(0, info.param.file.find('_')) + "_" + info.param.name;
    for (char& c : out) {
        if (!std::isalnum(static_cast<unsigned char>(c))) {
            c = '_';
        }
    }
    return out;
}

// The Reporter's log lines that belong to the analyses with these labels: their "Preparing"
// lines and the warnings logged until the next analysis is prepared.
std::vector<std::string> ownLines(const std::vector<std::string>& lines, const std::vector<std::string>& labels) {
    std::vector<std::string> out;
    bool own = false;
    for (const std::string& line : lines) {
        if (line.rfind("Preparing ", 0) == 0) {
            own = std::any_of(labels.begin(), labels.end(),
                              [&line](const std::string& l) { return line == "Preparing " + l + " ..."; });
        } else if (line.rfind("Analysing ", 0) == 0) {
            own = false;
            continue;
        }
        if (own) {
            out.push_back(line);
        }
    }
    return out;
}

// Records frames that differ from MDATK: allowed (and printed) where the margin is a last-bit flip.
class FrameCheck {
public:
    FrameCheck(const Value& expectedResult, std::string context)
        : m_margins(expectedResult.find("margins")), m_context(std::move(context)) {}

    void differs(std::size_t k, const std::string& what) {
        const double margin = m_margins ? frameMargin(*m_margins, k) : std::numeric_limits<double>::infinity();
        if (margin > 0 && margin < kFlipMargin) {
            std::cout << "[last-bit flip] " << m_context << " " << what << " frame " << k << " (margin " << margin
                      << ")\n";
            return;
        }
        ADD_FAILURE() << m_context << ": " << what << " differs at frame " << k << " (margin " << margin << ")";
    }

private:
    const Value* m_margins;
    std::string m_context;
};

// One frame's stacking map entries (i, j, v), sorted
using Entries = std::vector<std::tuple<std::int64_t, std::int64_t, double>>;

void comparePairmap(const core::AnalysisResult& got, const Value& expected, const std::string& prefix,
                    std::size_t nFrames, FrameCheck& check, int& reordered) {
    const Value& arrays = at(expected, "arrays");
    const auto expectedInts = [&](const std::string& name) {
        std::vector<std::int64_t> out;
        for (const Value& v : list(at(at(arrays, name), "values"))) {
            out.push_back(integer(v));
        }
        return out;
    };
    std::vector<double> expectedV;
    for (const Value& v : list(at(at(arrays, prefix + "_v"), "values"))) {
        expectedV.push_back(hexFloat(v));
    }
    const std::vector<std::int64_t> eo = expectedInts(prefix + "_offsets"), ei = expectedInts(prefix + "_i"),
                                     ej = expectedInts(prefix + "_j");
    const std::vector<std::int64_t> go = integerValues(*got.array(prefix + "_offsets")),
                                     gi = integerValues(*got.array(prefix + "_i")),
                                     gj = integerValues(*got.array(prefix + "_j"));
    const auto gvSpan = got.array(prefix + "_v")->values<float>();
    const std::vector<double> gv(gvSpan.begin(), gvSpan.end());
    ASSERT_EQ(eo.size(), nFrames + 1);
    ASSERT_EQ(go.size(), nFrames + 1);
    EXPECT_EQ(go.front(), 0);
    const auto frameEntries = [](const std::vector<std::int64_t>& o, const std::vector<std::int64_t>& i,
                                 const std::vector<std::int64_t>& j, const std::vector<double>& v, std::size_t k) {
        Entries out;
        for (auto e = static_cast<std::size_t>(o[k]); e < static_cast<std::size_t>(o[k + 1]); ++e) {
            out.emplace_back(i[e], j[e], v[e]);
        }
        return out;
    };
    for (std::size_t k = 0; k < nFrames; ++k) {
        Entries a = frameEntries(go, gi, gj, gv, k);
        Entries b = frameEntries(eo, ei, ej, expectedV, k);
        if (a != b) {
            ++reordered;  // MDATK's order is the cKDTree's
        }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        if (a != b) {
            check.differs(k, prefix + " entries");
        }
    }
}

void compareResult(const core::AnalysisResult& got, const Value& expected, const core::Project& project,
                   const Value& fingerprint, std::size_t nFrames, const std::string& context) {
    SCOPED_TRACE(context);
    EXPECT_EQ(got.name, str(at(expected, "instance")));
    EXPECT_EQ(sortedTables(Value(got.params)), sortedTables(untag(at(expected, "params"))));
    EXPECT_EQ(sortedTables(Value(got.meta)), sortedTables(untag(at(expected, "meta"))))
        << core::toJson(Value(got.meta), -1) << "\nvs\n" << core::toJson(at(expected, "meta"), -1);
    // output order (add_series / add_matrix / add_pairmap order)
    const Value metaValue(got.meta);
    const Value& outputs = at(metaValue, "outputs");
    for (const char* kind : {"series", "matrices", "pairmaps"}) {
        std::vector<std::string> names;
        for (const auto& entry : *outputs.find(kind)->get<Table>()) {
            names.push_back(entry.first);
        }
        EXPECT_EQ(names, strings(at(at(expected, "output_order"), kind))) << kind;
    }
    // the cosine thresholds: MDBridge computes libm's cos(x * (pi / 180)); numpy's may differ
    if (const Value* thresholds = at(expected, "margins").find("thresholds")) {
        for (const auto& [name, t] : *thresholds->get<Table>()) {
            const double libm = std::cos(hexFloat(at(t, "angle")) * (3.141592653589793 / 180.0));
            EXPECT_EQ(libm, hexFloat(at(t, "libm"))) << name;
            if (hexFloat(at(t, "numpy")) != hexFloat(at(t, "libm"))) {
                std::cout << "[numpy cos] " << context << " " << name << ": numpy " << str(at(t, "numpy"))
                          << " libm " << str(at(t, "libm")) << "\n";
            }
        }
    }
    const AnalysisDescriptor& descriptor = analysisFor(got.name);
    EXPECT_EQ(descriptor.cacheKey(project, untag(fingerprint)), str(at(expected, "key")));

    // arrays
    const Table& expectedArrays = *at(expected, "arrays").get<Table>();
    std::vector<std::string> gotNames, expectedNames;
    for (const auto& [name, array] : got.arrays) {
        gotNames.push_back(name);
    }
    for (const auto& [name, array] : expectedArrays) {
        expectedNames.push_back(name);
    }
    std::sort(gotNames.begin(), gotNames.end());
    EXPECT_EQ(gotNames, expectedNames);
    FrameCheck check(expected, context);
    int reordered = 0;
    for (const auto& [name, e] : expectedArrays) {
        const core::NumArray* a = got.array(name);
        if (!a) {
            continue;
        }
        EXPECT_EQ(numpyDtype(*a), str(at(e, "dtype"))) << name;
        std::vector<std::size_t> shape;
        for (const Value& s : list(at(e, "shape"))) {
            shape.push_back(static_cast<std::size_t>(integer(s)));
        }
        EXPECT_EQ(a->shape(), shape) << name;
        if (a->shape() != shape) {
            continue;
        }
        if (name.rfind("map_", 0) == 0) {
            continue;  // below, per frame
        }
        const std::vector<std::int64_t> values = integerValues(*a);
        const Array& ev = list(at(e, "values"));
        ASSERT_EQ(values.size(), ev.size()) << name;
        // per-frame arrays [n], or [rows, n] (occupancy)
        for (std::size_t x = 0; x < values.size(); ++x) {
            if (values[x] != integer(ev[x])) {
                check.differs(x % nFrames, name + "[" + std::to_string(x) + "] " + std::to_string(values[x]) +
                                               " vs " + std::to_string(integer(ev[x])));
            }
        }
    }
    if (got.array("map_offsets")) {
        comparePairmap(got, expected, "map", nFrames, check, reordered);
        if (reordered) {
            std::cout << "[pair order] " << context << ": " << reordered
                      << " frame(s) list the stacking map in another order than MDATK (cKDTree order)\n";
        }
    }
}

class MgStackingParity : public testing::TestWithParam<CaseRef> {};

TEST_P(MgStackingParity, MatchesMdatk) {
    const CaseRef& ref = GetParam();
    const Value& file = caseFile(ref.file);
    const Value& c = list(at(file, "cases"))[ref.index];
    ASSERT_EQ(c.find("error"), nullptr) << "a case of these analyses with an error: extend the test";

    const core::Project project = caseProject(c);
    core::System system = caseSystem(file, str(at(c, "system")), project.inputs.dt_ps);

    std::vector<std::string> logs;
    const core::Reporter reporter([&logs](const std::string& m) { logs.push_back(m); });
    std::vector<std::string> instances;
    for (const std::string& inst : strings(at(c, "instances"))) {
        if (isMine(inst)) {
            instances.push_back(inst);
        }
    }
    std::vector<std::unique_ptr<AnalysisContext>> contexts;
    std::vector<std::unique_ptr<Analysis>> analyses;
    std::vector<std::string> labels;
    for (const std::string& inst : instances) {
        contexts.push_back(std::make_unique<AnalysisContext>(system, instanceProject(project, inst), "", reporter));
        analyses.push_back(analysisFor(inst).create(*contexts.back(), inst));
        labels.push_back(analyses.back()->label());
    }
    std::vector<std::int64_t> frames;
    for (const Value& f : list(at(c, "frames"))) {
        frames.push_back(integer(f));
    }
    // the frame current when MDATK set the first of these analyses up (selections may depend on it)
    system.trajectory->seek(integer(at(at(c, "setup_frame"), instances.front())));

    RunResult run;
    try {
        run = runAnalyses(*contexts.front(), analyses, frames, &reporter);
    } catch (const std::logic_error& e) {
        if (std::string(e.what()).find("periodic neighbour search not implemented yet") != std::string::npos) {
            GTEST_SKIP() << "use_pbc: " << e.what();
        }
        throw;
    }

    // times and the analyses' own log lines (Preparing, warnings)
    const Array& times = list(at(c, "times"));
    ASSERT_EQ(run.timePs.size(), times.size());
    for (std::size_t k = 0; k < times.size(); ++k) {
        EXPECT_EQ(run.timePs[k], hexFloat(times[k])) << "time " << k;
    }
    EXPECT_EQ(ownLines(logs, labels), ownLines(strings(at(c, "log")), labels));

    const Array& expectedResults = list(at(c, "results"));
    ASSERT_EQ(run.results.size(), instances.size());
    for (std::size_t r = 0; r < instances.size(); ++r) {
        const Value* expected = nullptr;
        for (const Value& e : expectedResults) {
            if (str(at(e, "instance")) == instances[r]) {
                expected = &e;
            }
        }
        ASSERT_NE(expected, nullptr) << instances[r];
        compareResult(run.results[r], *expected, contexts[r]->project(), at(c, "fingerprint"), frames.size(),
                      ref.name + " " + instances[r]);
    }
}

INSTANTIATE_TEST_SUITE_P(Cases, MgStackingParity, testing::ValuesIn(caseRefs()), testName);

// ---- the stacking module's helpers ------------------------------------------------------------

// ring_geometry against numpy (spatial_extra_cases.json "rings"): centroids exactly, planarity and
// normals (up to sign) within LAPACK's last bits
TEST(Stacking, RingGeometryMatchesNumpy) {
    const Value data = test::readJson(test::kData / "analysis" / "spatial_extra_cases.json");
    int checked = 0;
    for (const Value& c : list(at(data, "rings"))) {
        const auto k = static_cast<std::size_t>(integer(at(c, "k")));
        const auto points = [](const Value& v) {
            std::vector<core::Vec3> out;
            const Array& a = list(v);
            for (std::size_t i = 0; i + 2 < a.size(); i += 3) {
                out.push_back({hexFloat(a[i]), hexFloat(a[i + 1]), hexFloat(a[i + 2])});
            }
            return out;
        };
        const std::vector<core::Vec3> pos = points(at(c, "positions"));
        const std::vector<core::Vec3> centroids = points(at(c, "centroids"));
        const std::vector<core::Vec3> normals = points(at(c, "normals"));
        const Array& planarity = list(at(c, "planarity"));
        const RingGeometry g = ringGeometry(pos, k);
        ASSERT_EQ(g.centroids.size(), centroids.size());
        for (std::size_t m = 0; m < centroids.size(); ++m) {
            const std::string context = "ring size " + std::to_string(k) + " #" + std::to_string(m);
            for (std::size_t i = 0; i < 3; ++i) {
                EXPECT_EQ(g.centroids[m][i], centroids[m][i]) << context;
            }
            EXPECT_NEAR(g.planarity[m], hexFloat(planarity[m]), 1e-14) << context;
            const core::Vec3& n = g.normals[m];
            const double dot = n[0] * normals[m][0] + n[1] * normals[m][1] + n[2] * normals[m][2];
            const double sign = dot < 0 ? -1.0 : 1.0;
            for (std::size_t i = 0; i < 3; ++i) {
                EXPECT_NEAR(sign * n[i], normals[m][i], 1e-13) << context;
            }
            ++checked;
        }
    }
    EXPECT_GT(checked, 0);
    EXPECT_THROW((void)ringGeometry(std::vector<core::Vec3>(5), 3), std::invalid_argument);
    EXPECT_THROW((void)ringGeometry(std::vector<core::Vec3>(3), 0), std::invalid_argument);
}

// ring_atoms: purine before pyrimidine before the side-chain table (by the upper-cased residue
// name), every atom required, a repeated name refers to its last atom
TEST(Stacking, RingAtoms) {
    auto t = std::make_shared<core::Topology>();
    const auto residue = [&t](const std::string& resname, const std::vector<std::string>& names) {
        const auto r = static_cast<std::int32_t>(t->resids.size());
        t->resids.push_back(r + 1);
        t->resnames.push_back(t->strings.intern(resname));
        t->residueSegment.push_back(0);
        for (const std::string& n : names) {
            t->names.push_back(t->strings.intern(n));
            t->atomResidue.push_back(r);
        }
    };
    residue("G", {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "O6"});  // atoms 0-9
    residue("C", {"N1", "C2", "N3", "C4", "C5", "C6", "C4"});  // 10-16: C4 twice
    residue("phe", {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"});  // 17-22
    residue("HSD", {"CG", "ND1", "CD2", "CE1"});  // 23-26: NE2 missing
    residue("TRP", {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"});  // 27-32
    t->segids.push_back(t->strings.intern("SYS"));
    t->finish();
    using V = std::vector<std::int32_t>;
    EXPECT_EQ(stackingRingAtoms(*t, 0), (V{0, 1, 2, 3, 4, 5, 6, 7, 8}));
    EXPECT_EQ(stackingRingAtoms(*t, 1), (V{10, 11, 12, 16, 14, 15}));
    EXPECT_EQ(stackingRingAtoms(*t, 2), (V{17, 18, 19, 20, 21, 22}));
    EXPECT_EQ(stackingRingAtoms(*t, 3), std::nullopt);
    EXPECT_EQ(stackingRingAtoms(*t, 4), (V{27, 28, 29, 30, 31, 32}));
}

// ---- mg: the octahedral test with linear pairs (the parity data has none outside use_pbc) -----

// One Mg at the origin and five oxygens 2 Å away: four of residue 1 (group A) in the xy plane at
// 0, 180 (150 in frame 1), 90 and 120 degrees, one of residue 2 (group B) on +z. With the default
// tolerance of 25 degrees a pair bridges when its angle is within 65-115 or >= 155 degrees.
TEST(Mg, OctahedralGeometry) {
    auto t = std::make_shared<core::Topology>();
    t->elements.emplace();
    const auto atom = [&t](const std::string& name, const std::string& element, std::int32_t residue) {
        t->names.push_back(t->strings.intern(name));
        t->elements->push_back(t->strings.intern(element));
        t->atomResidue.push_back(residue);
    };
    for (const auto& [resid, resname] : std::vector<std::pair<int, std::string>>{{1, "A"}, {2, "U"}, {3, "MG"}}) {
        t->resids.push_back(resid);
        t->resnames.push_back(t->strings.intern(resname));
        t->residueSegment.push_back(0);
    }
    atom("O1", "O", 0);
    atom("O2", "O", 0);
    atom("O3", "O", 0);
    atom("O4", "O", 0);
    atom("O5", "O", 1);
    atom("MG", "MG", 2);
    t->segids.push_back(t->strings.intern("SYS"));
    t->finish();
    const auto onCircle = [](double degrees) {
        const double r = degrees * (std::acos(-1.0) / 180.0);
        return std::vector<float>{static_cast<float>(2 * std::cos(r)), static_cast<float>(2 * std::sin(r)), 0.0f};
    };
    std::vector<core::Frame> frames(2);
    for (std::size_t f = 0; f < 2; ++f) {
        for (const double deg : {0.0, f == 0 ? 180.0 : 150.0, 90.0, 120.0}) {
            const std::vector<float> p = onCircle(deg);
            frames[f].positions.insert(frames[f].positions.end(), p.begin(), p.end());
        }
        frames[f].positions.insert(frames[f].positions.end(), {0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f});
    }
    core::System system{t, std::make_unique<core::MemoryTrajectory>(std::move(frames), 1.0)};
    core::Project project;
    project.groups = {{"Aptamer", "1", "", "#2a78d6"}, {"Target", "2", "", "#2a78d6"}};
    project.analyses.mg.mg_selection = "resname MG";
    const std::vector<std::int64_t> frameList{0, 1};
    const auto arrays = [&](bool geometry) {
        project.analyses.mg.geometry = geometry;
        AnalysisContext context(system, project);
        std::vector<std::unique_ptr<Analysis>> analyses;
        analyses.push_back(mgAnalysis().create(context));
        return runAnalyses(context, analyses, frameList).results.at(0);
    };
    using V = std::vector<std::int64_t>;
    const core::AnalysisResult with = arrays(true);
    EXPECT_EQ(integerValues(*with.array("contacts")), (V{5, 5}));
    EXPECT_EQ(integerValues(*with.array("bridges")), (V{7, 5}));
    EXPECT_EQ(integerValues(*with.array("inter")), (V{4, 4}));
    EXPECT_EQ(integerValues(*with.array("bridging_ions")), (V{1, 1}));
    EXPECT_EQ(integerValues(*with.array("occupancy")), (V{2, 2, 2, 2}));
    EXPECT_EQ(integerValues(*with.array("inner_residues")), (V{2, 2}));
    const core::AnalysisResult without = arrays(false);
    EXPECT_EQ(integerValues(*without.array("bridges")), (V{10, 10}));
    EXPECT_EQ(integerValues(*without.array("inter")), (V{4, 4}));
}

}  // namespace
