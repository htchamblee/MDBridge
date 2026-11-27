// The rmsd analysis (src/analysis/rmsd.cpp, mdatk/analysis/rmsd.py) against what MDATK computes on
// every recorded case that runs it (data/analysis: micro systems, the demo project and its
// variants, the generated systems): series, meta, outputs, params, cache keys, warnings, errors and
// the frame the trajectory is on at every setup (rmsd seeks its reference frame in setup, and the
// analyses set up after it start from there).
//
// Tolerances (docs/specs/m3-plan.md 4.3): the rms.rmsd path (superposition without a fit
// selection) and the plain path (no superposition) must match bit for bit (NaN equals NaN); the
// fit-selection path goes through kabsch, whose SVD and BLAS products differ from LAPACK / OpenBLAS
// in the last bits, so it is compared within |a - e| <= 1e-12 * max(1, |e|). So are the values of
// the improved RMSD (AN-16 / AN-17: the exact superposition RMSD where the expectation holds numpy's
// SVD value); values below 1e-4 Å may also differ by up to 1e-4 Å (QCP or exact value near the
// AN-17 threshold).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/rmsd_qcp.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

using namespace mdbridge;
using namespace mdbridge::test::analysis_cases;

namespace {

constexpr double kFitTolerance = 1e-12;
constexpr double kNearZero = 1e-4;  // Å

struct Tally {
    int cases = 0;
    int results = 0;
    int exact = 0;
    int withinTolerance = 0;
    double worstRelative = 0.0;
};

bool runsRmsd(const Value& c) {
    for (const std::string& instance : strings(at(c, "instances"))) {
        if (analysis::parseInstance(instance).base == "rmsd") {
            return true;
        }
    }
    return false;
}

void compareRmsdResult(const ProbedRun& run, std::size_t i, const core::AnalysisResult& actual, const Value& expected,
                       const Value& fingerprint, Tally& tally) {
    const std::string& instance = run.instances[i];
    SCOPED_TRACE(instance);
    const core::Project& project = run.projects[i];
    ++tally.results;
    EXPECT_EQ(actual.name, instance);
    EXPECT_EQ(analysis::rmsdAnalysis().cacheKey(project, fingerprint), str(at(expected, "key")));
    EXPECT_EQ(core::pythonJsonCompact(Value(actual.params)), core::pythonJsonCompact(at(expected, "params")));
    EXPECT_EQ(core::pythonJsonCompact(Value(actual.meta)), core::pythonJsonCompact(at(expected, "meta")));
    // output order: the series in the order they were added
    std::vector<std::string> series;
    const Value meta(actual.meta);
    const Value* outputs = meta.path({"outputs", "series"});
    ASSERT_NE(outputs, nullptr);
    for (const auto& entry : *outputs->get<Table>()) {
        series.push_back(entry.first);
    }
    EXPECT_EQ(series, strings(at(at(expected, "output_order"), "series")));

    const auto& targets = project.analyses.rmsd.targets;
    const Table& arrays = *at(expected, "arrays").get<Table>();
    ASSERT_EQ(actual.arrays.size(), arrays.size());
    for (std::size_t n = 0; n < arrays.size(); ++n) {
        const auto& [name, want] = arrays[n];
        SCOPED_TRACE(name);
        // the recorded arrays are sorted by name; ours are in insertion order
        const core::NumArray* got = actual.array(name);
        ASSERT_NE(got, nullptr);
        EXPECT_EQ(str(at(want, "dtype")), "<f8");
        EXPECT_EQ(got->dtypeName(), "float64");
        ASSERT_EQ(got->shape(), (std::vector<std::size_t>{static_cast<std::size_t>(integer(list(at(want, "shape"))[0]))}));
        const std::size_t target = std::stoul(name.substr(5));
        const bool fitPath = project.analyses.rmsd.superposition && target < targets.size() &&
                             !core::pythonStrip(targets[target].fit_selection).empty();
        const auto values = got->values<double>();
        const Array& wantValues = list(at(want, "values"));
        ASSERT_EQ(values.size(), wantValues.size());
        for (std::size_t k = 0; k < values.size(); ++k) {
            const double e = hexFloat(wantValues[k]);
            const double a = values[k];
            if (sameDouble(a, e)) {
                ++tally.exact;
                continue;
            }
            const double deviation = std::fabs(a - e);
            // AN-16 / AN-17 values come from the exact path (numpy's SVD there, a Jacobi
            // eigen-solution here): the same tolerance as kabsch's; below kNearZero the QCP and the
            // exact value may be chosen differently near the AN-17 threshold
            const bool exactPath = !fitPath && project.analyses.rmsd.superposition;
            const bool nearZero = std::min(std::fabs(a), std::fabs(e)) < kNearZero && deviation <= kNearZero;
            if ((fitPath || exactPath) && std::isfinite(e) &&
                (deviation <= kFitTolerance * std::max(1.0, std::fabs(e)) || (exactPath && nearZero))) {
                ++tally.withinTolerance;
                tally.worstRelative = std::max(tally.worstRelative, deviation / std::max(1.0, std::fabs(e)));
                continue;
            }
            ADD_FAILURE() << "frame position " << k << ": expected " << std::hexfloat << e << ", got " << a
                          << std::defaultfloat << (fitPath ? " (fit path)" : "");
        }
    }
}

// Every case of `fileName` that runs rmsd (or fails in its setup).
void checkCaseFile(const std::string& fileName, Tally& tally) {
    const Value& file = caseFile(fileName);
    for (const Value& c : list(at(file, "cases"))) {
        if (!runsRmsd(c)) {
            continue;
        }
        const std::string& name = str(at(c, "name"));
        SCOPED_TRACE(fileName + " " + name);
        ++tally.cases;
        const test::LogCapture quiet;  // the warnings are checked through the reporter
        const ProbedRun run = runCaseWithProbes(file, c, [](const std::string& base) { return base == "rmsd"; });

        // frame state at every setup (rmsd moves it to its reference frame)
        std::map<std::string, std::int64_t> setupFrames;
        for (const auto& [instance, frame] : *at(c, "setup_frame").get<Table>()) {
            setupFrames[instance] = integer(frame);
        }
        EXPECT_EQ(run.setupFrames, setupFrames);
        EXPECT_EQ(run.log, logWithoutOtherWarnings(c, "RMSD"));

        if (const Value* error = c.find("error")) {
            ASSERT_TRUE(run.error.has_value()) << "expected " << str(at(*error, "type"));
            EXPECT_EQ(run.error->type, str(at(*error, "type")));
            EXPECT_EQ(run.error->message, str(at(*error, "message")));
            continue;
        }
        ASSERT_FALSE(run.error.has_value()) << run.error->type << ": " << run.error->message;
        std::vector<std::int64_t> frames;
        for (const Value& f : list(at(c, "frames"))) {
            frames.push_back(integer(f));
        }
        EXPECT_EQ(run.frames, frames);
        const Array& times = list(at(c, "times"));
        ASSERT_EQ(run.run.timePs.size(), times.size());
        for (std::size_t k = 0; k < times.size(); ++k) {
            EXPECT_EQ(run.run.timePs[k], hexFloat(times[k])) << "time " << k;
        }
        ASSERT_EQ(run.run.results.size(), run.instances.size());
        for (std::size_t i = 0; i < run.instances.size(); ++i) {
            if (!run.real[i]) {
                continue;
            }
            const Value* expected = recordedResult(c, run.instances[i]);
            ASSERT_NE(expected, nullptr) << run.instances[i];
            compareRmsdResult(run, i, run.run.results[i], *expected, at(c, "fingerprint"), tally);
        }
    }
    std::printf("%s: %d cases, %d rmsd results, %d values exact, %d within %g only (worst %.3g)\n", fileName.c_str(),
                tally.cases, tally.results, tally.exact, tally.withinTolerance, kFitTolerance, tally.worstRelative);
}

}  // namespace

TEST(RmsdAnalysis, MicroCases) {
    Tally tally;
    checkCaseFile("micro_cases.json", tally);
    EXPECT_EQ(tally.cases, 7);
    EXPECT_EQ(tally.results, 8);
}

TEST(RmsdAnalysis, DemoDryDefault) {
    Tally tally;
    checkCaseFile("demo_dry_cases.json", tally);
    EXPECT_EQ(tally.cases, 1);
    EXPECT_EQ(tally.results, 1);
}

TEST(RmsdAnalysis, DemoDryVariants) {
    Tally tally;
    checkCaseFile("demo_dry_variants.json", tally);
    EXPECT_EQ(tally.cases, 18);
    EXPECT_EQ(tally.results, 18);
}

TEST(RmsdAnalysis, Systems) {
    Tally tally;
    checkCaseFile("systems_cases.json", tally);  // pdb_blown*: NaN for a blown-up fit frame (AN-20)
    EXPECT_EQ(tally.cases, 13);
    EXPECT_EQ(tally.results, 13);
}

// ---------------------------------------------------------------- hand-made checks

namespace {

// 12 atoms in residues 1 and 2 (6 each), frames given as positions
core::System smallSystem(const std::vector<std::vector<float>>& frames, double mass = 12.0) {
    auto t = std::make_shared<core::Topology>();
    for (int r = 0; r < 2; ++r) {
        t->resids.push_back(r + 1);
        t->resnames.push_back(t->strings.intern("A"));
        t->residueSegment.push_back(0);
    }
    for (int k = 0; k < 12; ++k) {
        t->names.push_back(t->strings.intern("C" + std::to_string(k)));
        t->atomResidue.push_back(k / 6);
    }
    t->masses = std::vector<double>(12, mass);
    t->segids.push_back(t->strings.intern("SYS"));
    t->finish();
    std::vector<core::Frame> list;
    for (const auto& positions : frames) {
        core::Frame f;
        f.positions = positions;
        list.push_back(std::move(f));
    }
    return core::System{t, std::make_unique<core::MemoryTrajectory>(std::move(list), 1.0)};
}

std::vector<float> points(double scale, double shift) {
    std::vector<float> out;
    for (int k = 0; k < 12; ++k) {
        out.push_back(static_cast<float>(scale * std::sin(1.3 * k) + shift));
        out.push_back(static_cast<float>(scale * std::cos(0.7 * k * k)));
        out.push_back(static_cast<float>(scale * std::sin(0.4 * k + 1.0) * 2.0));
    }
    return out;
}

// the same points rotated by 90 degrees about z and moved
std::vector<float> rotated(const std::vector<float>& p) {
    std::vector<float> out;
    for (std::size_t i = 0; i < p.size(); i += 3) {
        out.push_back(-p[i + 1] + 5.0f);
        out.push_back(p[i]);
        out.push_back(p[i + 2] - 2.0f);
    }
    return out;
}

// conftest.micro_project with the groups of the two residues and one rmsd target over both
core::Project rmsdProject() {
    core::Project p = microProject({{"Aptamer", "1"}, {"Target", "2"}});
    p.analyses.rmsd.targets = {{"All", "group:Aptamer or group:Target", ""}};
    return p;
}

core::AnalysisResult runRmsd(core::System& system, const core::Project& project) {
    return runOne("rmsd", system, project);
}

}  // namespace

// MDATK's tests/test_analyses.py::test_rmsd_superposition on our own rigid motion
TEST(RmsdAnalysis, SuperpositionRemovesRigidMotion) {
    const std::vector<float> x = points(3.0, 0.0);
    core::System system = smallSystem({x, rotated(x)});
    core::Project p = rmsdProject();
    EXPECT_NEAR(runRmsd(system, p).array("rmsd_0")->values<double>()[1], 0.0, 1e-4);
    p.analyses.rmsd.targets[0].fit_selection = "group:Aptamer";
    EXPECT_NEAR(runRmsd(system, p).array("rmsd_0")->values<double>()[1], 0.0, 1e-4);
    p.analyses.rmsd.mass_weighted = true;
    EXPECT_NEAR(runRmsd(system, p).array("rmsd_0")->values<double>()[1], 0.0, 1e-4);
    p.analyses.rmsd.superposition = false;
    EXPECT_GT(runRmsd(system, p).array("rmsd_0")->values<double>()[1], 1.0);
}

TEST(RmsdAnalysis, ZeroMassesAndErrors) {
    const std::vector<float> x = points(3.0, 0.0);
    core::System system = smallSystem({x, rotated(x)}, 0.0);
    core::Project p = rmsdProject();
    p.analyses.rmsd.mass_weighted = true;
    // rms.rmsd's and np.average's ZeroDivisionError, raised while processing
    EXPECT_THROW((void)runRmsd(system, p), analysis::ZeroDivisionError);
    p.analyses.rmsd.superposition = false;
    try {
        (void)runRmsd(system, p);
        ADD_FAILURE() << "no error";
    } catch (const analysis::ZeroDivisionError& e) {
        EXPECT_EQ(std::string(e.what()), "Weights sum to zero, can't be normalized");
    }
    // a fit selection matching nothing is a selection error (allow_empty only for the target)
    p = rmsdProject();
    p.analyses.rmsd.targets[0].fit_selection = "name ZZ";
    EXPECT_THROW((void)runRmsd(system, p), core::ValueError);
    // a blank fit selection counts as none
    p.analyses.rmsd.targets[0].fit_selection = " \t";
    EXPECT_NO_THROW((void)runRmsd(system, p));
}

// SEL-17 / AN-13 (docs/DIFFERENCES.md): every setup starts on the first analysed frame, also after
// rmsd's setup moved the trajectory to its reference frame (MDATK: the later setups start there).
TEST(RmsdAnalysis, EverySetupStartsOnTheFirstAnalysedFrame) {
    std::vector<std::vector<float>> frames;
    for (int f = 0; f < 6; ++f) {
        frames.push_back(points(3.0, f * 0.1));
    }
    core::System system = smallSystem(frames);
    core::Project p = rmsdProject();
    p.analyses.rmsd.reference_frame = 2;
    analysis::AnalysisContext context(system, p);
    std::map<std::string, std::int64_t> setupFrames;
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    analyses.push_back(std::make_unique<Recorded>(analysis::contactsAnalysis(), context, "contacts", nullptr,
                                                  setupFrames));
    analyses.push_back(std::make_unique<Recorded>(analysis::rmsdAnalysis(), context, "rmsd",
                                                  analysis::rmsdAnalysis().create(context), setupFrames));
    analyses.push_back(std::make_unique<Recorded>(analysis::traceAnalysis(), context, "trace", nullptr, setupFrames));
    const std::vector<std::int64_t> analysed{1, 3, 4, 5};
    const analysis::RunResult run = analysis::runAnalyses(context, analyses, analysed);
    // reference position 2 of the analysed frames is frame 4; the setups all start on frame 1
    EXPECT_EQ(setupFrames, (std::map<std::string, std::int64_t>{{"contacts", 1}, {"rmsd", 1}, {"trace", 1}}));
    const auto values = run.results.at(1).array("rmsd_0")->values<double>();
    EXPECT_EQ(values[2], 0.0);  // AN-17: the reference frame against itself is exactly 0
    EXPECT_GT(values[0], 0.0);
    EXPECT_EQ(*core::Value(run.results.at(1).meta).find("reference_frame_index")->get<std::int64_t>(), 2);
}
