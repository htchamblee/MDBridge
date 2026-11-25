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
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
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
            std::vector<std::string> ours;
            for (const std::string& line : logs) {
                ours.push_back(withPlaceholders(line));  // AN-06 names the topology file
            }
            EXPECT_EQ(ours, strings(at(c, "log")));
        } else {
            std::vector<std::string> expected;
            for (const Value& line : expectedLog) {
                if (str(line).starts_with("WARNING: energy:")) {
                    expected.push_back(str(line));
                }
            }
            std::vector<std::string> ours;
            for (const std::string& line : logs) {
                if (line.starts_with("WARNING:")) {
                    ours.push_back(withPlaceholders(line));  // AN-06 names the topology file
                }
            }
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
    // AN-25: the one-atom residues NA and CL are the TIP3P ion templates (Joung-Cheatham 2008)
    const double rmin = 1.369 + 2.513;
    const double eps = std::sqrt(0.0874393 * 0.035591);
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
    EXPECT_EQ(*res.meta[0].second.get<std::string>(), "templates");
    project.analyses.energy.dielectric = 4.0;
    project.analyses.energy.distance_dependent = true;
    res = run(project);
    near(value(res, "elec"), -kEnergyCoulomb / (4.0 * r * r), 1e-9);
}

// AN-05 (docs/DIFFERENCES.md): for overlapping atoms (r < 0.5 A) the net force is the force at the
// clamped distance 0.5 along the unit vector; MDATK scaled it by r / 0.5. Open and periodic paths.
TEST(Energy, NetForceAtClampedDistanceHasFullMagnitude) {
    const double r = 0.3;
    const std::vector<MicroAtom> atoms = {{"NA", 1, "NA", "NA", 22.99, 1.0}, {"CL", 2, "CL", "CL", 35.45, -1.0}};
    core::Project project;
    project.groups = {{"Aptamer", "1", "", "#2a78d6"}, {"Target", "2", "", "#2a78d6"}};
    const double rmin = 1.369 + 2.513;  // the TIP3P ion templates NA, CL (AN-25)
    const double eps = std::sqrt(0.0874393 * 0.035591);
    const double d = 0.5;
    const double dE = kEnergyCoulomb / (d * d) - 12 * eps * std::pow(rmin, 12) / std::pow(d, 13) +
                      12 * eps * std::pow(rmin, 6) / std::pow(d, 7);  // qq = -332.0522
    for (const bool periodic : {false, true}) {
        // across the face of a 20 A cell with use_pbc: 19.85 and 0.15 -> 0.3 A apart
        const std::vector<double> xyz = periodic ? std::vector<double>{19.85, 5, 5, 0.15, 5, 5}
                                                 : std::vector<double>{5, 5, 5, 5.3, 5, 5};
        core::System system = makeMicroSystem(atoms, {xyz});
        core::Project p = project;
        if (periodic) {
            core::Frame frame;
            for (const double v : xyz) {
                frame.positions.push_back(static_cast<float>(v));
            }
            frame.box = std::array<float, 6>{20, 20, 20, 90, 90, 90};
            std::vector<core::Frame> frames;
            frames.push_back(std::move(frame));
            system.trajectory = std::make_unique<core::MemoryTrajectory>(std::move(frames));
            p.analyses.use_pbc = true;
        }
        const core::AnalysisResult res = runOne("energy", system, p);
        const double force = res.array("force")->values<double>()[0];
        EXPECT_LE(std::abs(force - std::abs(dE)), 1e-6 * std::abs(dE)) << periodic << ": " << force << " vs " << dE
                                                                        << " (MDATK: " << std::abs(dE) * r / d << ")";
    }
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
                          "Lennard-Jones parameters come from force-field residue templates and Amber atom types");
    EXPECT_EQ(*res.meta[0].second.get<std::string>(), "templates");
    EXPECT_EQ(*res.meta[1].second.get<std::string>(), "topology");
}

// AN-07 / CIT-03: exclusions along guessed bonds are announced, citing the covalent radii the guess
// uses (Cordero et al. 2008)
TEST(Energy, GuessedBondsWarningCitesCovalentRadii) {
    core::Project project;
    project.groups = {{"Aptamer", "1", "", "#2a78d6"}, {"Target", "2", "", "#2a78d6"}};
    core::System system = microSystem(at(at(microCases(), "systems"), "energy_ions"));
    std::vector<std::string> logs;
    const core::Reporter reporter([&](const std::string& m) { logs.push_back(m); });
    AnalysisContext context(system, project, "", reporter);
    std::vector<std::unique_ptr<Analysis>> analyses;
    analyses.push_back(energyAnalysis().create(context));
    const std::vector<std::int64_t> frames{0};
    runAnalyses(context, analyses, frames);
    EXPECT_NE(std::find(logs.begin(), logs.end(),
                        "WARNING: energy: the topology has no bonds; bonded pairs are excluded along bonds guessed "
                        "from distances on the first analysed frame (covalent radii: Cordero et al., Dalton Trans. "
                        "2008, 2832-2838, doi:10.1039/b801115j)"),
              logs.end())
        << joined(logs);
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

// ---------------------------------------------------------------- AN-25 - AN-27, CIT-04

namespace {

// A micro system without topology charges (as a PDB or GRO file has none)
core::System withoutCharges(core::System system) {
    auto t = std::make_shared<core::Topology>(*system.topology);
    t->charges.reset();
    system.topology = t;
    return system;
}

std::int64_t metaInt(const core::AnalysisResult& res, std::string_view table, std::string_view key) {
    for (const auto& [name, value] : res.meta) {
        if (name == table) {
            return *value.find(key)->get<std::int64_t>();
        }
    }
    throw std::logic_error("no meta " + std::string(table));
}

std::string metaString(const core::AnalysisResult& res, std::string_view key) {
    for (const auto& [name, value] : res.meta) {
        if (name == key) {
            return *value.get<std::string>();
        }
    }
    throw std::logic_error("no meta " + std::string(key));
}

void expectLj(const EnergyParameters& p, std::size_t k, double rstar, double epsilon) {
    EXPECT_NEAR(p.rstar[k], rstar, 1e-12) << "atom " << k;
    EXPECT_NEAR(p.epsilon[k], epsilon, 1e-12) << "atom " << k;
}

}  // namespace

// AN-25: every residue of the demo (OL3 nucleotides with every atom, Mg2+) matches a template; its
// atoms take the template's Lennard-Jones parameters and charges (Amber's, converted back from
// OpenMM's sigma / epsilon)
TEST(EnergyParameters, DemoMatchesTheTemplates) {
    core::System system = io::loadSystem(test::kData / "topology" / "demo_dry.parm7",
                                         {test::kData / "io" / "trajectories" / "demo_dry.nc"}, std::nullopt);
    const core::AtomGroup all = system.atoms();
    const EnergyParameters p = energyParameters(all, *system.currentFrame());
    const core::Topology& t = *system.topology;
    ASSERT_EQ(p.rstar.size(), all.size());
    for (std::size_t k = 0; k < all.size(); ++k) {
        EXPECT_EQ(p.layer[k], EnergyParameters::Layer::Template) << k;
        EXPECT_TRUE(p.templateCharge[k]) << k;
        const std::string& name = t.name(static_cast<std::int32_t>(k));
        const std::string& resname = t.resname(t.atomResidue[k]);
        if (resname == "G" && name == "P") {
            expectLj(p, k, 2.1, 0.2);  // P
            EXPECT_DOUBLE_EQ(p.charge[k], 1.1662);
        } else if (name == "C1'") {
            expectLj(p, k, 1.908, 0.1094);  // CT
        } else if (name == "O2'") {
            expectLj(p, k, 1.721, 0.2104);  // OH
        } else if (name == "O4'") {
            expectLj(p, k, 1.6837, 0.17);  // OS
        } else if (name == "H8") {
            expectLj(p, k, 1.359, 0.015);  // H5
        } else if (name == "MG") {
            expectLj(p, k, 1.36, 0.01020237);  // Li-Merz 2013 12-6 CM, TIP3P
            EXPECT_EQ(p.charge[k], 2.0);
        }
    }
    EXPECT_EQ(p.forceFields, (std::vector<std::string>{"OL3", "TIP3P"}));
    // CIT-04: the originators of the charges and Lennard-Jones parameters used (OL3's are ff94's)
    EXPECT_EQ(p.references, (std::vector<std::string>{"cornell1995", "li2013"}));
    EXPECT_EQ(energyShortCitation("li2013"), "Li et al. 2013, doi:10.1021/ct400146w");
}

// AN-25: residues with missing atoms (no hydrogens) take the Lennard-Jones parameters of the template
// for their position in the chain (here an N-terminal and a C-terminal alanine, bonded by the
// topology), but no charges; an ion by its residue name; a residue name of pdbNames.xml (HIS) by the
// templates it stands for (HID / HIE / HIP agree on every heavy atom)
TEST(EnergyParameters, IncompleteResiduesAndTermini) {
    // ALA - HIS - ALA bonded by the topology (heavy atoms only), then a sodium ion
    const std::vector<std::vector<const char*>> residues{
        {"N", "CA", "C", "O", "CB"}, {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"}, {"N", "CA", "C", "O", "CB"}};
    std::vector<MicroAtom> atoms;
    std::vector<double> xyz;
    std::vector<std::pair<std::int32_t, std::int32_t>> bonds;
    std::int32_t previousC = -1;
    for (std::size_t r = 0; r < residues.size(); ++r) {
        const auto base = static_cast<std::int32_t>(atoms.size());
        for (const char* name : residues[r]) {
            atoms.push_back({name, static_cast<std::int64_t>(r + 1), r == 1 ? "HIS" : "ALA", std::string(1, name[0])});
            xyz.insert(xyz.end(), {static_cast<double>(atoms.size()) * 3.0, 0.0, 0.0});
        }
        bonds.insert(bonds.end(), {{base, base + 1}, {base + 1, base + 2}, {base + 2, base + 3}, {base + 1, base + 4}});
        if (previousC >= 0) {
            bonds.emplace_back(previousC, base);  // C(i-1) - N(i)
        }
        previousC = base + 2;
    }
    atoms.push_back({"Na+", 4, "Na+", "NA"});
    xyz.insert(xyz.end(), {100.0, 0.0, 0.0});
    core::System system = withoutCharges(makeMicroSystem(atoms, {xyz}, bonds));
    const EnergyParameters p = energyParameters(system.atoms(), *system.currentFrame());
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        EXPECT_EQ(p.layer[k], EnergyParameters::Layer::Template) << k;
    }
    expectLj(p, 1, 1.908, 0.1094);  // CA: CX
    expectLj(p, 3, 1.6612, 0.21);   // O
    expectLj(p, 11, 1.824, 0.17);   // ND1 of HIS: NA in HID / HIP, NB in HIE, the same values
    expectLj(p, 20, 1.369, 0.0874393);  // Na+ (Joung-Cheatham 2008)
    EXPECT_EQ(p.charge[20], 1.0);
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        EXPECT_EQ(p.templateCharge[k], k == 20) << k;
    }
    EXPECT_EQ(p.chargelessResidues, 3u);
    EXPECT_EQ(p.forceFields, (std::vector<std::string>{"ff14SB", "TIP3P"}));
}

// AN-26: atoms of residues without a template take Amber atom types: CT for a carbon with four
// neighbours, OH / OS / OW for oxygens by their hydrogens, hydrogens by their bonded atom (here by
// distance: no topology bonds), ions alone in a residue by element, carbon's for unknown elements
TEST(EnergyParameters, TypedLayer) {
    // methanol (C bonded to O and three H; O-H), dimethyl ether's O between two carbons, a water
    // named LIG, a sodium named SOD, an unknown element
    const std::vector<MicroAtom> atoms{
        {"C1", 1, "MOH", "C"},  {"O1", 1, "MOH", "O"},  {"H1", 1, "MOH", "H"}, {"H2", 1, "MOH", "H"},
        {"H3", 1, "MOH", "H"},  {"HO", 1, "MOH", "H"},  {"CA", 2, "DME", "C"}, {"OE", 2, "DME", "O"},
        {"CB", 2, "DME", "C"},  {"OW", 3, "LIG", "O"},  {"HW1", 3, "LIG", "H"}, {"HW2", 3, "LIG", "H"},
        {"SOD", 4, "SOD", "NA"}, {"X1", 5, "UNK", "XX"}, {"NZ", 6, "AMN", "N"}, {"HZ", 6, "AMN", "H"}};
    const std::vector<double> xyz{
        0, 0, 0,     1.43, 0, 0,  -0.36, 1.0, 0,  -0.36, -0.5, 0.87, -0.36, -0.5, -0.87, 1.75, 0.9, 0,  // methanol
        20, 0, 0,    21.41, 0, 0, 22.0, 1.3, 0,                                                        // ether
        40, 0, 0,    40.96, 0, 0, 39.76, 0.93, 0,                                                      // water
        60, 0, 0,    80, 0, 0,    100, 0, 0,  101.01, 0, 0};
    core::System system = withoutCharges(makeMicroSystem(atoms, {xyz}));
    const EnergyParameters p = energyParameters(system.atoms(), *system.currentFrame());
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        EXPECT_EQ(p.layer[k], EnergyParameters::Layer::Typed) << k;
        EXPECT_FALSE(p.templateCharge[k]) << k;
    }
    expectLj(p, 0, 1.908, 0.1094);     // CT
    expectLj(p, 1, 1.721, 0.2104);     // OH
    expectLj(p, 2, 1.387, 0.0157);     // H1 (on carbon)
    expectLj(p, 5, 0.0, 0.0);          // HO (on oxygen)
    expectLj(p, 6, 1.908, 0.086);      // C (two neighbours)
    expectLj(p, 7, 1.6837, 0.17);      // OS
    expectLj(p, 9, 1.7683, 0.152);     // OW
    expectLj(p, 12, 1.369, 0.0874393);  // Na+ by element
    expectLj(p, 13, 1.908, 0.086);     // unknown: carbon's
    expectLj(p, 15, 0.6, 0.0157);      // H (on nitrogen)
    EXPECT_EQ(p.defaultElements, std::vector<std::string>{"XX"});
    EXPECT_EQ(p.defaultAtoms, 1u);
    EXPECT_EQ(p.references, (std::vector<std::string>{"cornell1995", "jorgensen1983", "joung2008"}));
}

// AN-27 / CIT-04 in the analysis: template charges only for complete residues, a warning naming the
// rest, the metadata per layer and one log line with the references; topology charges keep
// precedence
TEST(EnergyParameters, AnalysisChargesMetadataAndLog) {
    const test::LogCapture quiet;
    const std::vector<MicroAtom> atoms{{"O", 1, "WAT", "O"}, {"H1", 1, "WAT", "H"}, {"H2", 1, "WAT", "H"},
                                       {"O", 2, "HOH", "O"}, {"OW", 3, "SOL", "O"}, {"HW1", 3, "SOL", "H"},
                                       {"HW2", 3, "SOL", "H"}};
    const std::vector<double> xyz{0, 0, 0, 0.9572, 0, 0, -0.24, 0.927, 0, 3, 0, 0,
                                  6, 0, 0, 6.9572, 0, 0, 5.76, 0.927, 0};
    core::Project project;
    project.groups = {{"Aptamer", "1-2", "", "#2a78d6"}, {"Target", "3", "", "#2a78d6"}};
    project.analyses.energy.exclude_bonded = 0;
    const auto run = [&](core::System& system, std::vector<std::string>& log) {
        const core::Reporter reporter([&log](const std::string& m) { log.push_back(m); });
        AnalysisContext context(system, project, "", reporter);
        std::vector<std::unique_ptr<Analysis>> analyses;
        analyses.push_back(energyAnalysis().create(context));
        return runAnalyses(context, analyses, std::vector<std::int64_t>{0}).results.at(0);
    };
    core::System bare = withoutCharges(makeMicroSystem(atoms, {xyz}));
    std::vector<std::string> log;
    const core::AnalysisResult res = run(bare, log);
    ASSERT_GE(log.size(), 3u);
    EXPECT_EQ(log[1], "WARNING: energy: the topology has no partial charges; 6 of 7 atoms take the charges of the "
                      "force-field residue template they match completely (hydrogens included); the other 1 atom "
                      "(1 residue) have no charge, so their electrostatic terms are zero");
    EXPECT_EQ(log[2], "energy: Lennard-Jones parameters without parm7: 7 atoms from Amber force-field residue "
                      "templates (TIP3P); partial charges of 6 atoms from the templates. Cite: Jorgensen et al. 1983, "
                      "doi:10.1063/1.445869");
    EXPECT_EQ(metaString(res, "lj_source"), "templates");
    // MDATK 0.1.1's key: the table the parameters come from (docs/GENERIC_LJ.md)
    EXPECT_NE(metaString(res, "lj_table").find("Amber force-field residue templates"), std::string::npos);
    EXPECT_NE(metaString(res, "lj_table").find("(docs/GENERIC_LJ.md)"), std::string::npos);
    EXPECT_EQ(metaString(res, "charge_source"), "partial");
    EXPECT_EQ(metaInt(res, "lj_atoms", "template"), 7);
    EXPECT_EQ(metaInt(res, "charge_atoms", "template"), 6);
    EXPECT_EQ(metaInt(res, "charge_atoms", "none"), 1);
    // water - water at 3 and 6 Å: TIP3P charges -0.834 / 0.417
    EXPECT_NE(res.array("elec")->values<double>()[0], 0.0);

    core::System charged = makeMicroSystem(atoms, {xyz});  // topology charges (all 0 here)
    log.clear();
    const core::AnalysisResult topo = run(charged, log);
    EXPECT_EQ(metaString(topo, "charge_source"), "topology");
    EXPECT_EQ(metaInt(topo, "charge_atoms", "topology"), 7);
    EXPECT_EQ(metaInt(topo, "charge_atoms", "template"), 0);
    EXPECT_EQ(topo.array("elec")->values<double>()[0], 0.0);
    EXPECT_EQ(log[1], "energy: Lennard-Jones parameters without parm7: 7 atoms from Amber force-field residue "
                      "templates (TIP3P); partial charges from the topology. Cite: Jorgensen et al. 1983, "
                      "doi:10.1063/1.445869");
}

// AN-26: extra points (virtual sites, lone pairs) have no Lennard-Jones term (Amber's EP) and are
// nobody's neighbour; a residue with one matches a template by its other atoms but takes no charges
TEST(EnergyParameters, ExtraPoints) {
    const std::vector<MicroAtom> atoms{
        {"OW", 1, "SOL", "O"}, {"HW1", 1, "SOL", "H"}, {"HW2", 1, "SOL", "H"}, {"MW", 1, "SOL", ""},    // GROMACS TIP4P
        {"O", 2, "WAT", "O"},  {"H1", 2, "WAT", "H"},  {"H2", 2, "WAT", "H"},  {"EPW", 2, "WAT", ""},   // Amber TIP4P-Ew
        {"OH2", 3, "TP4E", ""}, {"H1", 3, "TP4E", ""}, {"H2", 3, "TP4E", ""}, {"OM", 3, "TP4E", ""},  // CHARMM
        {"C1", 4, "LIG", "C"}, {"OM", 4, "LIG", "O"}};
    std::vector<double> xyz;
    for (int w = 0; w < 3; ++w) {
        const double x = 10.0 * w;
        xyz.insert(xyz.end(), {x, 0, 0, x + 0.9572, 0, 0, x - 0.24, 0.927, 0, x + 0.09, 0.12, 0});
    }
    xyz.insert(xyz.end(), {40, 0, 0, 41.2, 0, 0});
    core::System system = withoutCharges(makeMicroSystem(atoms, {xyz}));
    const EnergyParameters p = energyParameters(system.atoms(), *system.currentFrame());
    for (const std::size_t k : {3u, 7u, 11u}) {
        EXPECT_EQ(p.layer[k], EnergyParameters::Layer::Typed) << k;
        expectLj(p, k, 0.0, 0.0);
    }
    for (const std::size_t k : {0u, 1u, 2u, 4u, 5u, 6u}) {
        EXPECT_EQ(p.layer[k], EnergyParameters::Layer::Template) << k;  // TIP3P's HOH by the other atoms
        EXPECT_FALSE(p.templateCharge[k]) << k;
    }
    expectLj(p, 0, 1.7683, 0.152);
    expectLj(p, 8, 1.7683, 0.152);   // TP4E (no template): OW by its two hydrogens, OM no neighbour
    expectLj(p, 9, 0.0, 0.0);        // H on O
    expectLj(p, 13, 1.6612, 0.21);   // OM of a ligand: an oxygen
    EXPECT_TRUE(p.defaultElements.empty());
    EXPECT_EQ(p.references, (std::vector<std::string>{"jorgensen1983", "cornell1995"}));
}

// AN-25: chain ends without hydrogens. A chain breaks between two residues only on a new segment /
// chain ID, or on a gap in the numbering without a bond; consecutive numbering in one chain links
// residues even far apart (a synthetic backbone), an insertion code counts as consecutive. Seen on
// the middle residue's O3': OS (middle form A) when linked to the next residue, OH (3' form A3)
// when its chain ends there.
TEST(EnergyParameters, ChainEndsFromChainIdsNumberingAndBonds) {
    const std::vector<const char*> heavy{"OP1", "OP2", "P", "O5'", "C5'", "C4'", "O4'", "C1'", "N9", "C8", "N7",
                                         "C5", "C6", "N6", "N1", "C2", "N3", "C4", "C3'", "C2'", "O2'", "O3'"};
    // three A nucleotides 30 Å apart (no bond guessed between them); the first without its phosphate
    const auto build = [&](std::array<std::int64_t, 3> resids, std::array<const char*, 3> chains,
                           std::array<const char*, 3> icodes, bool bond23) {
        std::vector<MicroAtom> atoms;
        std::vector<double> xyz;
        std::vector<std::pair<std::int32_t, std::int32_t>> bonds;
        std::int32_t o3 = -1;
        std::int32_t p3 = -1;
        std::size_t third = 0;  // the first atom of residue 3 (with its phosphate a 5' end has no template)
        for (std::size_t r = 0; r < 3; ++r) {
            if (r == 2) {
                third = atoms.size();
            }
            for (std::size_t k = r == 0 ? 3 : 0; k < heavy.size(); ++k) {
                const std::string name = heavy[k];
                if (r == 1 && name == "O3'") o3 = static_cast<std::int32_t>(atoms.size());
                if (r == 2 && name == "P") p3 = static_cast<std::int32_t>(atoms.size());
                atoms.push_back({name, static_cast<std::int64_t>(r + 1), "A", std::string(1, name[0])});
                xyz.insert(xyz.end(), {30.0 * static_cast<double>(r), 1.5 * static_cast<double>(k), 0.0});
            }
        }
        if (bond23) {
            bonds.emplace_back(o3, p3);
        }
        core::System system = withoutCharges(makeMicroSystem(atoms, {xyz}, bonds));
        auto t = std::make_shared<core::Topology>(*system.topology);
        std::vector<std::uint32_t> chainIds;
        for (const std::int32_t r : t->atomResidue) {
            chainIds.push_back(t->strings.intern(chains[static_cast<std::size_t>(r)]));
        }
        t->chainIDs = std::move(chainIds);
        std::vector<std::uint32_t> codes;
        for (const char* code : icodes) {
            codes.push_back(t->strings.intern(code));
        }
        t->icodes = std::move(codes);
        for (std::size_t r = 0; r < 3; ++r) {  // the numbering under test (micro systems number 1, 2, 3)
            t->resids[r] = resids[r];
            t->resnums[r] = resids[r];
        }
        system.topology = t;
        const EnergyParameters p = energyParameters(system.atoms(), *system.currentFrame());
        for (std::size_t k = 0; k < third; ++k) {
            EXPECT_EQ(p.layer[k], EnergyParameters::Layer::Template) << k;
        }
        return std::pair{p.rstar[static_cast<std::size_t>(o3)], p.epsilon[static_cast<std::size_t>(o3)]};
    };
    const std::pair<double, double> os{1.6837, 0.17};
    const std::pair<double, double> oh{1.721, 0.2104};
    const auto expect = [](std::pair<double, double> got, std::pair<double, double> want, const char* what) {
        EXPECT_NEAR(got.first, want.first, 1e-12) << what;
        EXPECT_NEAR(got.second, want.second, 1e-12) << what;
    };
    expect(build({1, 2, 3}, {"A", "A", "A"}, {"", "", ""}, false), os, "consecutive, one chain: no break");
    expect(build({1, 2, 2}, {"A", "A", "A"}, {"", "", "A"}, false), os, "insertion code: consecutive");
    expect(build({1, 2, 3}, {"A", "A", "B"}, {"", "", ""}, false), oh, "new chain ID: a break");
    expect(build({1, 2, 7}, {"A", "A", "A"}, {"", "", ""}, false), oh, "numbering gap, no bond: a break");
    expect(build({1, 2, 7}, {"A", "A", "A"}, {"", "", ""}, true), os, "numbering gap with a bond: no break");
}

// MDATK 0.1.0's parm7 format lines: a CHAMBER topology ('(a80)' CTITLE, '(i2,a78)' FORCE_FIELD_TYPE)
// or an ff19SB one ('(8(F9.5))' CMAP parameters) gives its own Lennard-Jones parameters (MDATK
// 0.0.4.1: "cannot parse format line", and the energy fell back to other parameters)
TEST(Energy, Parm7WithChamberAndCmapFormatLines) {
    const std::filesystem::path demo = test::kData / "topology" / "demo_dry.parm7";
    const std::filesystem::path trajectory = test::kData / "io" / "trajectories" / "demo_dry.nc";
    test::TempDir dir;
    const std::filesystem::path chamber = dir / "chamber.parm7";
    std::string text = test::readFile(demo);
    const std::size_t at = text.find("%FLAG POINTERS");
    ASSERT_NE(at, std::string::npos);
    text.insert(at, "%FLAG CTITLE\n%FORMAT(a80)\nCHARMM topology\n%FLAG FORCE_FIELD_TYPE\n%FORMAT(i2,a78)\n"
                    " 1 CHARMM 36\n%FLAG CMAP_COUNT\n%FORMAT(2I8)\n       0       1\n%FLAG CMAP_PARAMETER_01\n"
                    "%FORMAT(8(F9.5))\n  0.12500 -1.50000  2.00000\n");
    test::writeFile(chamber, text);
    core::Project project;
    project.groups = {{"Aptamer", "1-13", "", "#2a78d6"}, {"Target", "14-26", "", "#2a78d6"}};
    core::System plain = io::loadSystem(demo, {trajectory}, std::nullopt);
    core::System modified = io::loadSystem(chamber, {trajectory}, std::nullopt);
    const core::AnalysisResult expected = runOne("energy", plain, project, demo.string());
    const core::AnalysisResult got = runOne("energy", modified, project, chamber.string());
    EXPECT_EQ(metaString(got, "lj_source"), "parm7");
    EXPECT_TRUE(std::none_of(got.meta.begin(), got.meta.end(), [](const auto& e) { return e.first == "lj_table"; }));
    EXPECT_EQ(got.array("vdw")->values<double>()[0], expected.array("vdw")->values<double>()[0]);
    EXPECT_EQ(got.array("elec")->values<double>()[0], expected.array("elec")->values<double>()[0]);
}

namespace {

std::vector<std::string> splitOn(const std::string& line, char sep) {
    std::vector<std::string> out(1);
    for (const char c : line) {
        if (c == sep) {
            out.emplace_back();
        } else {
            out.back() += c;
        }
    }
    return out;
}

std::string trimmed(const std::string& text) {
    const std::size_t first = text.find_first_not_of(" \r");
    const std::size_t last = text.find_last_not_of(" \r");
    return first == std::string::npos ? std::string() : text.substr(first, last - first + 1);
}

// The rows of a tab-separated table with a header line: one map per row, by column name.
std::vector<std::map<std::string, std::string>> readTsv(const std::filesystem::path& path) {
    std::vector<std::map<std::string, std::string>> out;
    const std::vector<std::string> lines = splitOn(test::readFile(path), '\n');
    const std::vector<std::string> header = splitOn(trimmed(lines.at(0)), '\t');
    for (std::size_t k = 1; k < lines.size(); ++k) {
        if (trimmed(lines[k]).empty()) {
            continue;
        }
        const std::vector<std::string> cells = splitOn(trimmed(lines[k]), '\t');
        std::map<std::string, std::string> row;
        for (std::size_t c = 0; c < header.size() && c < cells.size(); ++c) {
            row[header[c]] = cells[c];
        }
        out.push_back(std::move(row));
    }
    return out;
}

bool sameNumber(double a, double b) { return std::abs(a - b) <= 1e-9 * std::max(std::abs(a), std::abs(b)) + 1e-15; }

}  // namespace

// docs/GENERIC_LJ.md records every number of both parameter layers: its tables are those of the
// extracts the analysis' tables are generated from (data/sources/forcefields; tools/gen_energy_lj.py)
TEST(EnergyParameters, GenericLjDocument) {
    const std::filesystem::path root = test::kData.parent_path().parent_path();
    const std::filesystem::path sources = root / "data" / "sources" / "forcefields";
    std::set<std::string> typedDoc, typedSource;
    std::map<std::string, std::vector<std::string>> templateDoc;  // class -> {R*, epsilon, files}
    for (const std::string& line : splitOn(test::readFile(root / "docs" / "GENERIC_LJ.md"), '\n')) {
        if (line.rfind("| `", 0) != 0) {
            continue;
        }
        std::vector<std::string> cells = splitOn(trimmed(line), '|');
        cells = std::vector<std::string>(cells.begin() + 1, cells.end() - 1);
        for (std::string& c : cells) {
            c = trimmed(c);
        }
        const std::string name = cells[0].substr(1, cells[0].size() - 2);
        if (cells.size() == 6) {  // layer 2: type, R*, epsilon, file, line, reference
            typedDoc.insert(name + "|" + cells[1] + "|" + cells[2] + "|" + cells[3] + "|" + cells[4]);
        } else if (cells.size() == 5 && cells[3].find(".xml`") != std::string::npos) {  // layer 1
            std::string cls = name;
            if (cells[3] == "`amber14/tip3p.xml`") {
                const bool ion = name.find('+') != std::string::npos || name.find('-') != std::string::npos;
                cls = (ion ? "tip3p_standard-" : "tip3p-") + name;
            }
            EXPECT_TRUE(templateDoc.emplace(cls, std::vector<std::string>{cells[1], cells[2], cells[3]}).second) << cls;
        }
    }
    for (const auto& row : readTsv(sources / "amber_lj_types.tsv")) {
        typedSource.insert(row.at("type") + "|" + row.at("rstar_angstrom") + "|" + row.at("epsilon_kcal_per_mol") +
                           "|`" + row.at("file") + "`|" + row.at("line"));
    }
    EXPECT_EQ(typedDoc, typedSource);

    std::map<std::string, std::pair<double, double>> classes;
    std::map<std::string, std::set<std::string>> files;
    for (const auto& row : readTsv(sources / "openmm_amber14_templates.tsv")) {
        if (row.at("sigma_nm").empty()) {
            continue;
        }
        const double sigma = std::stod(row.at("sigma_nm"));
        const double epsilon = std::stod(row.at("epsilon_kj_per_mol"));
        classes.emplace(row.at("class"), std::pair{sigma * std::pow(2.0, 1.0 / 6.0) / 2.0 * 10.0, epsilon / 4.184});
        std::string file = row.at("file");
        file.replace(0, std::string("amber14_").size(), "amber14/");
        files[row.at("class")].insert("`" + file + "`");
    }
    EXPECT_EQ(templateDoc.size(), classes.size());
    for (const auto& [cls, values] : classes) {
        const auto found = templateDoc.find(cls);
        if (found == templateDoc.end()) {
            ADD_FAILURE() << "docs/GENERIC_LJ.md has no row of class " << cls;
            continue;
        }
        const std::vector<std::string>& doc = found->second;
        EXPECT_TRUE(sameNumber(std::stod(doc[1]), values.second)) << cls << ": " << doc[1];
        if (values.second == 0.0) {
            EXPECT_EQ(doc[0], "-") << cls;  // no Lennard-Jones term
        } else {
            EXPECT_TRUE(sameNumber(std::stod(doc[0]), values.first)) << cls << ": " << doc[0];
        }
        std::string listed;
        for (const std::string& f : files[cls]) {
            listed += (listed.empty() ? "" : ", ") + f;
        }
        EXPECT_EQ(doc[2], listed) << cls;
    }
}
