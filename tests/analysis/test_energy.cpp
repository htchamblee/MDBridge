// The energy analysis (mdatk/analysis/energy.py) against what MDATK computes for the cases of
// tests/parity/make_analysis_cases.py (data/analysis, see data/README.md): the micro systems of
// MDATK's tests/test_analyses.py (analytic Coulomb + LJ, bonded exclusions), the demo project and
// its variants, and the generated systems; plus the counterparts of MDATK's
// test_energy_matches_analytic_coulomb_and_lj and test_energy_excludes_bonded_pairs_across_groups.
//
// Tolerances (docs/specs/m3-plan.md 4.3): pair counts, exclusions, meta, params, outputs, logs and
// times exact. MDBridge sums the pair terms in (i, j) order, MDATK in scipy's tree order, and libm's
// pow may differ from numpy's in the last bit, so elec / vdw / total / force are compared within
// |delta| <= kSumUlps * eps * sum(|term|), with the sums of |term| the generator records per frame
// (`margins`: abs_elec, abs_vdw, their sum for total, abs_force). Single-pair micro cases are also
// compared within a relative 1e-14 of the value itself.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/project.hpp"

using namespace mdbridge;
using namespace mdbridge::analysis;
using namespace mdbridge::test::analysis_cases;

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();
// k of |delta| <= k * eps * sum(|term|) for sums of up to ~4e4 terms taken in another order. The
// worst-case bound grows with log2(n) for numpy's pairwise sums (elec, vdw) and with n for the
// sequential net-force sum; the observed ratios (printed by the tests as "[ worst ]") stay below
// 2.1 for the energies and 25 for the force on the committed data.
constexpr double kSumUlps = 64.0;
constexpr double kForceUlps = 256.0;
// The pair search keeps d <= cutoff exactly as scipy does; a count may differ only where a pair
// lies within this relative margin of the cutoff (plan 4.3).
constexpr double kCutoffMargin = 1e-9;

const Value& microCases() { return caseFile("micro_cases.json"); }
const Value& demoCases() { return caseFile("demo_dry_cases.json"); }
const Value& variantCases() { return caseFile("demo_dry_variants.json"); }
const Value& systemCases() { return caseFile("systems_cases.json"); }

// The worst |delta| / (eps * sum(|term|)) seen per case group and quantity.
std::map<std::string, double>& worstRatios() {
    static std::map<std::string, double> ratios;
    return ratios;
}

void printWorst(const std::string& group) {
    for (const auto& [key, ratio] : worstRatios()) {
        if (key.starts_with(group + "/")) {
            const double k = key.ends_with("/force") ? kForceUlps : kSumUlps;
            std::cout << "[ worst ] " << key << ": " << ratio << " (tolerance " << k << ")\n";
        }
    }
}

bool isEnergyInstance(const std::string& instance) { return parseInstance(instance).base == "energy"; }

enum class Outcome { Ran, PeriodicMissing };

std::vector<double> sumOf(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> out(a.size());
    for (std::size_t k = 0; k < a.size(); ++k) {
        out[k] = a[k] + b[k];
    }
    return out;
}

void compareSeries(const std::string& ratioKey, const std::string& name, std::span<const double> ours,
                   const std::vector<double>& expected, const std::vector<double>& absSums, bool singlePair) {
    const double ulps = name == "force" ? kForceUlps : kSumUlps;
    ASSERT_EQ(ours.size(), expected.size()) << name;
    ASSERT_EQ(absSums.size(), expected.size()) << name;
    double& worst = worstRatios()[ratioKey + "/" + name];
    for (std::size_t k = 0; k < ours.size(); ++k) {
        const double delta = std::abs(ours[k] - expected[k]);
        const double scale = kEps * absSums[k];
        if (scale > 0) {
            worst = std::max(worst, delta / scale);
        }
        EXPECT_LE(delta, ulps * scale) << name << " frame " << k << ": " << ours[k] << " vs " << expected[k];
        if (singlePair) {
            EXPECT_LE(delta, 1e-14 * std::abs(expected[k])) << name << " frame " << k;
        }
    }
}

// Runs the energy instances of a case (others are left to their own tests) and compares them.
Outcome runEnergyCase(const Value& c, const Value& file, const std::string& group) {
    const std::string name = str(at(c, "name"));
    SCOPED_TRACE(group + " case " + name);
    const std::vector<std::string> allInstances = strings(at(c, "instances"));
    std::vector<std::string> instances;
    std::copy_if(allInstances.begin(), allInstances.end(), std::back_inserter(instances), isEnergyInstance);
    const bool energyOnly = instances.size() == allInstances.size();

    core::System system = caseSystem(file, str(at(c, "system")), std::nullopt);
    const std::string topologyPath = caseTopologyPath(file, str(at(c, "system")));
    const core::Project base = core::projectFromPlain(at(c, "project"));
    std::vector<std::string> logs;
    std::vector<std::pair<double, std::string>> progress;
    const core::Reporter reporter([&](const std::string& m) { logs.push_back(m); },
                                  [&](double f, const std::string& m) { progress.emplace_back(f, m); });
    std::vector<std::unique_ptr<AnalysisContext>> contexts;
    std::vector<std::unique_ptr<Analysis>> analyses;
    for (const std::string& instance : instances) {
        contexts.push_back(std::make_unique<AnalysisContext>(system, instanceProject(base, instance), topologyPath,
                                                             reporter));
        analyses.push_back(analysisFor(instance).create(*contexts.back(), instance));
    }
    const Value* framesValue = c.find("frames");  // absent for errors
    const std::vector<std::int64_t> frames = framesValue ? integers(*framesValue) : std::vector<std::int64_t>{};
    const Array& expectedLog = list(at(c, "log"));

    const auto checkLogs = [&] {
        if (energyOnly) {
            EXPECT_EQ(logs, strings(at(c, "log")));
        } else {
            std::vector<std::string> expected;
            for (const Value& line : expectedLog) {
                if (str(line).starts_with("WARNING: energy:")) {
                    expected.push_back(str(line));
                }
            }
            std::vector<std::string> ours;
            std::copy_if(logs.begin(), logs.end(), std::back_inserter(ours),
                         [](const std::string& line) { return line.starts_with("WARNING:"); });
            EXPECT_EQ(ours, expected);
        }
    };

    if (const Value* error = c.find("error")) {
        EXPECT_TRUE(energyOnly);
        try {
            (void)runAnalyses(*contexts.front(), analyses, frames);
            ADD_FAILURE() << "expected " << str(at(*error, "type"));
        } catch (const core::ValueError& e) {
            EXPECT_EQ(str(at(*error, "type")), "ValueError");
            EXPECT_EQ(std::string(e.what()), str(at(*error, "message")));
        }
        checkLogs();
        return Outcome::Ran;
    }

    RunResult run;
    try {
        run = runAnalyses(*contexts.front(), analyses, frames);
    } catch (const std::logic_error& e) {
        if (std::string(e.what()) == "periodic neighbour search not implemented yet") {
            return Outcome::PeriodicMissing;
        }
        throw;
    }
    checkLogs();
    if (energyOnly) {
        const Array& expectedProgress = list(at(c, "progress"));
        EXPECT_EQ(progress.size(), expectedProgress.size());
        for (std::size_t i = 0; i < std::min(progress.size(), expectedProgress.size()); ++i) {
            EXPECT_EQ(progress[i].first, hexFloat(list(expectedProgress[i])[0]));
            EXPECT_EQ(progress[i].second, str(list(expectedProgress[i])[1]));
        }
    }
    EXPECT_EQ(run.timePs, hexFloats(at(c, "times")));

    std::vector<const Value*> expectedResults;
    for (const Value& r : list(at(c, "results"))) {
        if (isEnergyInstance(str(at(r, "instance")))) {
            expectedResults.push_back(&r);
        }
    }
    EXPECT_EQ(expectedResults.size(), instances.size());
    for (std::size_t n = 0; n < std::min(expectedResults.size(), run.results.size()); ++n) {
        const Value& expected = *expectedResults[n];
        const core::AnalysisResult& ours = run.results[n];
        SCOPED_TRACE(ours.name);
        EXPECT_EQ(ours.name, str(at(expected, "instance")));
        EXPECT_EQ(ours.key, "");
        // the cache key for the case's fixed fingerprint
        EXPECT_EQ(analysisFor(ours.name).cacheKey(contexts[n]->project(), at(c, "fingerprint")),
                  str(at(expected, "key")));
        expectSameValue(at(expected, "params"), Value(ours.params), "params");
        const Value meta(ours.meta);
        expectSameValue(at(expected, "meta"), meta, "meta");
        // output order
        std::vector<std::string> seriesOrder;
        if (const Value* series = meta.path({"outputs", "series"})) {
            for (const auto& [key, value] : *series->get<Table>()) {
                seriesOrder.push_back(key);
            }
        }
        EXPECT_EQ(seriesOrder, strings(at(at(expected, "output_order"), "series")));

        // arrays: names in MDATK's order, dtypes, shapes, values
        std::vector<std::string> names;
        for (const auto& [key, value] : ours.arrays) {
            names.push_back(key);
        }
        EXPECT_EQ(names, (std::vector<std::string>{"total", "elec", "vdw", "force", "n_pairs"}));
        const Value& arrays = at(expected, "arrays");
        EXPECT_EQ(arrays.size(), ours.arrays.size());
        const Value& margins = at(at(expected, "margins"), "frames");
        const std::vector<double> absElec = hexFloats(at(margins, "abs_elec"));
        const std::vector<double> absVdw = hexFloats(at(margins, "abs_vdw"));
        const std::vector<double> absForce = hexFloats(at(margins, "abs_force"));
        const std::vector<double> cutoffMargin = hexFloats(at(margins, "cutoff"));

        const core::NumArray* pairs = ours.array("n_pairs");
        if (!pairs) {
            ADD_FAILURE() << "no n_pairs";
            continue;
        }
        EXPECT_EQ(pairs->dtypeName(), "int64");
        EXPECT_EQ(str(at(at(arrays, "n_pairs"), "dtype")), "<i8");
        const std::vector<std::int64_t> expectedPairs = integers(at(at(arrays, "n_pairs"), "values"));
        const std::vector<std::int64_t> ourPairs(pairs->values<std::int64_t>().begin(),
                                                 pairs->values<std::int64_t>().end());
        if (ourPairs.size() != expectedPairs.size()) {
            ADD_FAILURE() << "n_pairs: " << ourPairs.size() << " frames vs " << expectedPairs.size();
            continue;
        }
        const double cutoff = contexts[n]->project().analyses.energy.cutoff;
        bool samePairs = true;
        for (std::size_t k = 0; k < ourPairs.size(); ++k) {
            if (ourPairs[k] != expectedPairs[k]) {
                samePairs = false;
                EXPECT_LT(cutoffMargin.at(k), kCutoffMargin * cutoff)
                    << "n_pairs frame " << k << ": " << ourPairs[k] << " vs " << expectedPairs[k];
            }
        }
        if (!samePairs) {
            continue;  // sums over different pair sets are not comparable
        }
        bool singlePair = true;
        for (const std::int64_t p : ourPairs) {
            singlePair = singlePair && p <= 1;
        }
        const std::string ratioKey = group;
        for (const auto& [arrayName, sums] :
             std::vector<std::pair<std::string, std::vector<double>>>{{"total", sumOf(absElec, absVdw)},
                                                                      {"elec", absElec},
                                                                      {"vdw", absVdw},
                                                                      {"force", absForce}}) {
            const core::NumArray* array = ours.array(arrayName);
            if (!array) {
                ADD_FAILURE() << "no " << arrayName;
                continue;
            }
            EXPECT_EQ(array->dtypeName(), "float64");
            const Value& e = at(arrays, arrayName);
            EXPECT_EQ(str(at(e, "dtype")), "<f8");
            EXPECT_EQ(integers(at(e, "shape")), (std::vector<std::int64_t>{static_cast<std::int64_t>(frames.size())}));
            compareSeries(ratioKey, arrayName, array->values<double>(), hexFloats(at(e, "values")), sums, singlePair);
        }
    }
    return Outcome::Ran;
}

bool hasEnergy(const Value& c) {
    for (const std::string& instance : strings(at(c, "instances"))) {
        if (isEnergyInstance(instance)) {
            return true;
        }
    }
    return false;
}

bool usesPbc(const Value& c) {
    const Value* pbc = c.path({"project", "analyses", "use_pbc"});
    return pbc && pbc->get<bool>() && *pbc->get<bool>();
}

// Runs the energy cases of a file accepted by `filter`; returns the names skipped for want of the
// periodic search.
template <class Filter>
std::vector<std::string> runCases(const Value& file, const std::string& group, Filter&& filter) {
    std::vector<std::string> skipped;
    int ran = 0;
    for (const Value& c : list(at(file, "cases"))) {
        if (!hasEnergy(c) || !filter(c)) {
            continue;
        }
        if (runEnergyCase(c, file, group) == Outcome::PeriodicMissing) {
            skipped.push_back(str(at(c, "name")));
        } else {
            ++ran;
        }
    }
    printWorst(group);
    EXPECT_GT(ran + static_cast<int>(skipped.size()), 0);
    return skipped;
}

std::string joined(const std::vector<std::string>& names) {
    std::string out;
    for (const std::string& n : names) {
        out += (out.empty() ? "" : ", ") + n;
    }
    return out;
}

core::AtomGroup groupOf(const core::System& system, const std::vector<std::int64_t>& indices) {
    std::vector<std::int32_t> atoms;
    for (const std::int64_t i : indices) {
        atoms.push_back(static_cast<std::int32_t>(i));
    }
    return core::AtomGroup(system.topology, atoms);
}

}  // namespace

// ---------------------------------------------------------------- bonded exclusions

TEST(Energy, BondedExclusionKeys) {
    const Value& be = at(microCases(), "bonded_exclusions");
    const core::System system = microSystem(at(at(microCases(), "systems"), str(at(be, "system"))));
    const core::AtomGroup a = groupOf(system, integers(at(be, "a")));
    const core::AtomGroup b = groupOf(system, integers(at(be, "b")));
    const Array& cases = list(at(be, "cases"));
    ASSERT_EQ(cases.size(), 6u);
    for (const Value& c : cases) {
        const std::int64_t depth = integer(at(c, "depth"));
        EXPECT_EQ(energyBondedExclusions(a, b, depth), integers(at(c, "keys"))) << "depth " << depth;
    }
    EXPECT_TRUE(energyBondedExclusions(a, b, -1).empty());
}

// MDATK's test_energy_excludes_bonded_pairs_across_groups
TEST(Energy, BondedExclusionPairsAcrossGroups) {
    const core::System system = microSystem(at(at(microCases(), "systems"), "energy_chain"));
    const core::AtomGroup a = groupOf(system, {0, 1, 2});
    const core::AtomGroup b = groupOf(system, {3, 4, 5});
    std::vector<std::pair<std::int64_t, std::int64_t>> pairs;
    for (const std::int64_t k : energyBondedExclusions(a, b, 3)) {
        pairs.emplace_back(k / 3, k % 3);
    }
    std::sort(pairs.begin(), pairs.end());
    EXPECT_EQ(pairs, (std::vector<std::pair<std::int64_t, std::int64_t>>{{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}}));
    EXPECT_EQ(energyBondedExclusions(a, b, 1).size(), 1u);
    EXPECT_TRUE(energyBondedExclusions(a, b, 0).empty());
}

TEST(Energy, BondedExclusionsWithoutBonds) {
    const core::System system = microSystem(at(at(microCases(), "systems"), "energy_ions"));
    EXPECT_TRUE(energyBondedExclusions(groupOf(system, {0}), groupOf(system, {1}), 3).empty());
}

// ---------------------------------------------------------------- micro systems

TEST(Energy, MicroCases) {
    const auto skipped = runCases(microCases(), "micro", [](const Value&) { return true; });
    EXPECT_TRUE(skipped.empty()) << joined(skipped);
}

// MDATK's test_energy_matches_analytic_coulomb_and_lj, with its own formulas
TEST(Energy, AnalyticCoulombAndLennardJones) {
    const Value& spec = at(at(microCases(), "systems"), "energy_ions");
    core::Project project;
    project.groups = {{"Aptamer", "1", "", "#2a78d6"}, {"Target", "2", "", "#2a78d6"}};
    const auto run = [&](const core::Project& p) {
        core::System system = microSystem(spec);
        AnalysisContext context(system, p);
        std::vector<std::unique_ptr<Analysis>> analyses;
        analyses.push_back(energyAnalysis().create(context));
        const std::vector<std::int64_t> frames{0};
        return runAnalyses(context, analyses, frames).results.at(0);
    };
    const double r = 5.0;
    const double rmin = 1.369 + 2.513;  // GENERIC_LJ NA, CL
    const double eps = std::sqrt(0.0874 * 0.0355);
    const auto lj = [&](double x) { return eps * (std::pow(rmin / x, 12) - 2 * std::pow(rmin / x, 6)); };
    const auto value = [](const core::AnalysisResult& res, const char* name) {
        return res.array(name)->values<double>()[0];
    };
    const auto near = [](double ours, double expected, double rel) {
        EXPECT_LE(std::abs(ours - expected), rel * std::abs(expected)) << ours << " vs " << expected;
    };
    core::AnalysisResult res = run(project);
    near(value(res, "elec"), -kEnergyCoulomb / r, 1e-9);
    near(value(res, "vdw"), lj(r), 1e-9);
    const double h = 1e-5;
    const double de = (-kEnergyCoulomb / (r + h) + lj(r + h) - (-kEnergyCoulomb / (r - h) + lj(r - h))) / (2 * h);
    near(value(res, "force"), std::abs(de), 1e-5);
    EXPECT_EQ(*res.meta[0].second.get<std::string>(), "generic");
    project.analyses.energy.dielectric = 4.0;
    project.analyses.energy.distance_dependent = true;
    res = run(project);
    near(value(res, "elec"), -kEnergyCoulomb / (4.0 * r * r), 1e-9);
}

// use_parm7 with a parm7 file whose atom count differs from the system's: MDATK's warning and the
// generic parameters
TEST(Energy, Parm7AtomCountMismatch) {
    core::Project project;
    project.groups = {{"Aptamer", "1", "", "#2a78d6"}, {"Target", "2", "", "#2a78d6"}};
    core::System system = microSystem(at(at(microCases(), "systems"), "energy_ions"));
    std::vector<std::string> logs;
    const core::Reporter reporter([&](const std::string& m) { logs.push_back(m); });
    AnalysisContext context(system, project, (test::kData / "topology" / "demo_dry.parm7").string(), reporter);
    std::vector<std::unique_ptr<Analysis>> analyses;
    analyses.push_back(energyAnalysis().create(context));
    const std::vector<std::int64_t> frames{0};
    const core::AnalysisResult res = runAnalyses(context, analyses, frames).results.at(0);
    EXPECT_EQ(logs.at(1), "WARNING: energy: could not read LJ parameters from parm7 (atom count mismatch); "
                          "using generic element parameters");
    EXPECT_EQ(*res.meta[0].second.get<std::string>(), "generic");
    EXPECT_EQ(*res.meta[1].second.get<std::string>(), "topology");
}

// ---------------------------------------------------------------- demo project

TEST(Energy, DemoDefault) {
    const auto skipped = runCases(demoCases(), "demo", [](const Value&) { return true; });
    EXPECT_TRUE(skipped.empty()) << joined(skipped);
}

TEST(Energy, DemoVariants) {
    const auto skipped = runCases(variantCases(), "variants", [](const Value& c) { return !usesPbc(c); });
    EXPECT_TRUE(skipped.empty()) << joined(skipped);
}

TEST(Energy, DemoVariantsPeriodic) {
    const auto skipped = runCases(variantCases(), "variants_pbc", usesPbc);
    if (!skipped.empty()) {
        GTEST_SKIP() << "periodic neighbour search not implemented yet: " << joined(skipped);
    }
}

// ---------------------------------------------------------------- generated systems

TEST(Energy, Systems) {
    const auto skipped = runCases(systemCases(), "systems", [](const Value& c) { return !usesPbc(c); });
    EXPECT_TRUE(skipped.empty()) << joined(skipped);
}

TEST(Energy, SystemsPeriodic) {
    const auto skipped = runCases(systemCases(), "systems_pbc", usesPbc);
    if (!skipped.empty()) {
        GTEST_SKIP() << "periodic neighbour search not implemented yet: " << joined(skipped);
    }
}
