// Hand-made checks of MDBridge's deliberate improvements of the analyses (docs/DIFFERENCES.md,
// "Analyses"; docs/IMPROVEMENTS.md section 4 and SEL-17) that the recorded MDATK cases do not reach:
// the runner's warnings for unusable cells (AN-19) and non-finite coordinates (AN-20), the
// periodicity of the analyses' selections (SEL-17), bonded exclusions through atoms outside both
// groups (AN-08), the covalent-radius bond guess (AN-14) and translation-only fits (AN-18).

#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"

using namespace mdbridge;
using namespace mdbridge::test::analysis_cases;
using analysis::Analysis;
using analysis::AnalysisContext;
using analysis::AnalysisDescriptor;

namespace {

// An analysis that does nothing: for the runner's own behaviour.
class Idle final : public Analysis {
public:
    using Analysis::Analysis;
    void process(std::size_t, const core::Frame&) override {}
    core::AnalysisResult finalize() override { return newResult(); }
};

class IdleDescriptor final : public AnalysisDescriptor {
public:
    IdleDescriptor() : AnalysisDescriptor({"idle", "Idle", 1, true, 0, {}}) {}
    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<Idle>(*this, context, std::move(instance));
    }
};

const IdleDescriptor& idle() {
    static const IdleDescriptor descriptor;
    return descriptor;
}

std::vector<MicroAtom> twoAtoms() {
    return {{"A", 1, "AAA", "C"}, {"B", 2, "BBB", "C"}};
}

std::vector<std::string> warnings(const std::vector<std::string>& log) {
    std::vector<std::string> out;
    for (const std::string& line : log) {
        if (line.starts_with("WARNING: ")) {
            out.push_back(line);
        }
    }
    return out;
}

std::vector<std::string> runIdle(core::System& system, bool usePbc, const std::vector<std::int64_t>& frames,
                                 std::size_t threads) {
    std::vector<std::string> log;
    core::Reporter reporter([&log](const std::string& m) { log.push_back(m); });
    core::Project project;
    project.analyses.use_pbc = usePbc;
    AnalysisContext context(system, project, "", reporter);
    std::vector<std::unique_ptr<Analysis>> analyses;
    analyses.push_back(idle().create(context, "idle"));
    (void)analysis::runAnalyses(context, analyses, frames, nullptr, analysis::RunOptions{threads});
    return warnings(log);
}

}  // namespace

// AN-19: one warning per run for the first analysed frame without a usable cell, with use_pbc only.
TEST(AnalysisImprovements, UsePbcWithoutUsableCell) {
    const test::LogCapture quiet;
    for (const std::size_t threads : {1u, 4u}) {
        core::System system = makeMicroSystem(twoAtoms(), {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 1, 0, 0}});
        EXPECT_EQ(runIdle(system, true, {1, 2}, threads),
                  (std::vector<std::string>{"WARNING: use_pbc is on but trajectory frame 1 has no usable periodic "
                                            "cell; distances are not periodic"}));
        EXPECT_TRUE(runIdle(system, false, {1, 2}, threads).empty());
    }
}

// AN-20: frames with non-finite coordinates are reported once per run of consecutive analysed
// frames; the open pair search leaves such points out instead of failing.
TEST(AnalysisImprovements, NonFiniteCoordinates) {
    const test::LogCapture quiet;
    const double nan = std::nan("");
    const std::vector<double> good{0, 0, 0, 1, 0, 0};
    const std::vector<double> one{nan, 0, 0, 1, 0, 0};
    const std::vector<double> two{nan, 0, 0, 1, INFINITY, 0};
    for (const std::size_t threads : {1u, 4u}) {
        core::System system = makeMicroSystem(twoAtoms(), {good, one, two, good, one, one});
        EXPECT_EQ(runIdle(system, false, {0, 1, 2, 3, 4, 5}, threads),
                  (std::vector<std::string>{
                      "WARNING: frames 1-2: up to 2 atoms have non-finite coordinates; they are left out of "
                      "distance searches",
                      "WARNING: frames 4-5: up to 1 atom has non-finite coordinates; they are left out of "
                      "distance searches"}));
        EXPECT_EQ(runIdle(system, false, {2}, threads),
                  (std::vector<std::string>{"WARNING: frame 2: 2 atoms have non-finite coordinates; they are left "
                                            "out of distance searches"}));
    }
}

// SEL-17: distance keywords in the analyses' selections use the cell only with use_pbc, and are
// evaluated on the first analysed frame whatever moved the trajectory before.
TEST(AnalysisImprovements, SelectionsHonourUsePbc) {
    core::System system = makeMicroSystem(twoAtoms(), {{0.5, 5, 5, 9.5, 5, 5}, {0.5, 5, 5, 5.0, 5, 5}});
    auto* memory = dynamic_cast<core::MemoryTrajectory*>(system.trajectory.get());
    ASSERT_NE(memory, nullptr);
    // a 10 Å cubic cell on frame 0: A and B are 1 Å apart through the face
    std::vector<core::Frame> frames;
    for (std::int64_t f = 0; f < 2; ++f) {
        core::Frame frame = memory->seek(f);
        frame.box = std::array<float, 6>{10, 10, 10, 90, 90, 90};
        frames.push_back(std::move(frame));
    }
    system.trajectory = std::make_unique<core::MemoryTrajectory>(std::move(frames));
    system.trajectory->seek(0);
    core::Project project;
    project.analyses.use_pbc = false;
    const AnalysisContext open(system, project);
    EXPECT_TRUE(open.select("around 2 name A", true).empty());
    project.analyses.use_pbc = true;
    const AnalysisContext periodic(system, project);
    EXPECT_EQ(periodic.select("around 2 name A", true).indices(), (std::vector<std::int32_t>{1}));
}

// AN-08: a pair joined only through an atom outside both groups is excluded within the depth.
TEST(AnalysisImprovements, ExclusionsThroughOtherAtoms) {
    std::vector<MicroAtom> atoms{{"A", 1, "AAA", "C"}, {"X", 2, "XXX", "C"}, {"B", 3, "BBB", "C"}};
    core::System system = makeMicroSystem(atoms, {{0, 0, 0, 1.5, 0, 0, 3, 0, 0}}, {{0, 1}, {1, 2}});
    const core::AtomGroup a(system.topology, {0});
    const core::AtomGroup b(system.topology, {2});
    EXPECT_TRUE(analysis::energyBondedExclusions(a, b, 1).empty());
    EXPECT_EQ(analysis::energyBondedExclusions(a, b, 2), (std::vector<std::int64_t>{0}));  // MDATK: none
    // explicit bonds (the guessed ones of AN-07) replace the topology's
    const std::vector<core::Bond> none{{0, 1}};
    EXPECT_TRUE(analysis::energyBondedExclusions(a, b, 5, &none).empty());
}

// AN-14: covalent radii (a 2.05 Å disulfide is a bond, MDATK 0.0.4.1's 1.95 Å limit misses it), no
// bonds to one-atom residues; MDATK 0.1.2's residue rule: neighbouring residues whatever their
// chains, two sulfurs across any residues.
TEST(AnalysisImprovements, CovalentRadiusBondGuess) {
    std::vector<MicroAtom> atoms{{"SG", 1, "CYS", "S"}, {"SG", 2, "CYS", "S"}, {"CB", 2, "CYS", "C"},
                                 {"MG", 3, "MG", "MG"}, {"C1", 4, "LIG", "C"}, {"CB", 1, "CYS", "C"}};
    // S-S 2.05 Å; CB 1.8 Å from either S; Mg 2.0 Å from CB; C1 of residue 4 1.5 Å from Mg
    core::System system =
        makeMicroSystem(atoms, {{0, 0, 0, 2.05, 0, 0, 3.85, 0, 0, 5.85, 0, 0, 7.35, 0, 0, -1.8, 0, 0}});
    const core::AtomGroup all = system.atoms();
    EXPECT_EQ(analysis::guessBonds(all, *system.currentFrame()), (std::vector<std::int64_t>{0, 1, 0, 5, 1, 2}));

    // a disulfide between residues 1 and 3; two carbons as close across them are no bond
    std::vector<MicroAtom> bridge{{"SG", 1, "CYS", "S"}, {"CB", 1, "CYS", "C"}, {"CA", 2, "GLY", "C"},
                                  {"C", 2, "GLY", "C"},  {"SG", 3, "CYS", "S"}, {"CB", 3, "CYS", "C"}};
    core::System disulfide = makeMicroSystem(bridge, {{0, 0, 0, -1.8, 0, 0, 20, 0, 0, 21.5, 0, 0, 0, 2.05, 0,
                                                       -1.8, 1.5, 0}});
    EXPECT_EQ(analysis::guessBonds(disulfide.atoms(), *disulfide.currentFrame()),
              (std::vector<std::int64_t>{0, 1, 0, 4, 2, 3, 4, 5}));

    // two chains: neighbouring residues bond within bonding distance (O3'-P across a chain break),
    // not farther apart (MDATK 0.0.4.1's rule; MDBridge before 0.1.2's design: never across chains)
    std::vector<MicroAtom> pair{{"C3'", 1, "AAA", "C"}, {"O3'", 1, "AAA", "O"}, {"P", 2, "BBB", "P"},
                                {"O5'", 2, "BBB", "O"}};
    for (const double gap : {1.6, 3.0}) {
        core::System chains = makeMicroSystem(pair, {{-1.43, 0, 0, 0, 0, 0, gap, 0, 0, gap + 1.6, 0, 0}});
        auto t = std::make_shared<core::Topology>(*chains.topology);
        t->segids.push_back(t->strings.intern("B"));
        t->residueSegment[1] = 1;
        chains.topology = t;
        EXPECT_EQ(analysis::guessBonds(chains.atoms(), *chains.currentFrame()),
                  gap < 2.0 ? (std::vector<std::int64_t>{0, 1, 1, 2, 2, 3}) : (std::vector<std::int64_t>{0, 1, 2, 3}));
    }
}

// CIT-05: the atoms analysis says when it guesses its bonds, citing the covalent radii (once per
// run, only when the topology has no bonds among the structure atoms).
TEST(AnalysisImprovements, AtomsBondGuessWarningCitesCovalentRadii) {
    const test::LogCapture quiet;
    std::vector<MicroAtom> atoms{
        {"C1", 1, "LIG", "C"}, {"O1", 1, "LIG", "O"}, {"C2", 2, "LIG", "C"}, {"C3", 2, "LIG", "C"}};
    const std::vector<std::vector<double>> frames{{0, 0, 0, 1.2, 0, 0, 0, 1.5, 0, 0, 3, 0},
                                                  {0, 0, 0, 1.3, 0, 0, 0, 1.5, 0, 0, 3, 0}};
    core::Project project;
    project.groups = {{"Aptamer", "1-2", "", "#2a78d6"}};
    project.visual.structure_groups = {"Aptamer"};
    const auto run = [&](core::System& system) {
        std::vector<std::string> log;
        core::Reporter reporter([&log](const std::string& m) { log.push_back(m); });
        AnalysisContext context(system, project, "", reporter);
        std::vector<std::unique_ptr<Analysis>> analyses;
        analyses.push_back(analysis::atomsAnalysis().create(context, "atoms"));
        const analysis::RunResult result = analysis::runAnalyses(context, analyses, std::vector<std::int64_t>{0, 1});
        const auto bonds = result.results.at(0).array("bonds")->values<std::int64_t>();
        return std::pair{warnings(log), std::vector<std::int64_t>(bonds.begin(), bonds.end())};
    };
    core::System guessed = makeMicroSystem(atoms, frames);
    const auto [guessedWarnings, guessedBonds] = run(guessed);
    EXPECT_EQ(guessedWarnings,
              (std::vector<std::string>{
                  "WARNING: structural view: the topology has no bonds among the structure atoms; bonds are guessed "
                  "from distances on the first analysed frame (covalent radii: Cordero et al., Dalton Trans. 2008, "
                  "2832-2838, doi:10.1039/b801115j)"}));
    EXPECT_EQ(guessedBonds, (std::vector<std::int64_t>{0, 1, 0, 2, 2, 3}));
    core::System bonded = makeMicroSystem(atoms, frames, {{0, 1}});
    const auto [bondedWarnings, bondedBonds] = run(bonded);
    EXPECT_TRUE(bondedWarnings.empty());
    EXPECT_EQ(bondedBonds, (std::vector<std::int64_t>{0, 1}));
}

// AN-18: a fit selection of fewer than 3 atoms is translation only, with one warning.
TEST(AnalysisImprovements, TranslationOnlyFit) {
    const test::LogCapture quiet;
    std::vector<MicroAtom> atoms;
    for (int k = 0; k < 4; ++k) {
        atoms.push_back({"C" + std::to_string(k), 1, "AAA", "C"});
    }
    const std::vector<double> f0{0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    const std::vector<double> f1{2, 0, 0, 3, 0, 0, 2, 1, 0, 2, 0, 1};  // translated by 2 along x
    core::System system = makeMicroSystem(atoms, {f0, f1});
    core::Project project;
    project.analyses.rmsd.targets = {{"t", "all", "name C0 C1"}};
    std::vector<std::string> log;
    core::Reporter reporter([&log](const std::string& m) { log.push_back(m); });
    AnalysisContext context(system, project, "", reporter);
    std::vector<std::unique_ptr<Analysis>> analyses;
    analyses.push_back(analysis::rmsdAnalysis().create(context, "rmsd"));
    const analysis::RunResult run = analysis::runAnalyses(context, analyses, std::vector<std::int64_t>{0, 1});
    EXPECT_EQ(warnings(log), (std::vector<std::string>{"WARNING: RMSD target 't': fit selection has 2 atoms; "
                                                       "superposition is translation only"}));
    const auto values = run.results.at(0).array("rmsd_0")->values<double>();
    EXPECT_EQ(values[0], 0.0);
    EXPECT_NEAR(values[1], 0.0, 1e-12);
}
