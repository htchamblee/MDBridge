// The structural-view analyses trace, proximity and atoms (mdatk/analysis/structure.py) against what
// MDATK computes for the cases of tests/parity/make_analysis_cases.py (data/analysis, see
// data/README.md): the micro systems, the demo project and its variants, and the new small systems
// (periodic cells, PDB / GRO without bonds, heavy atoms only, aromatic residues).
//
// Each structure instance of a case runs in its own pass over the case's frames, with the
// trajectory at the frame MDATK's setup saw (setup_frame: after the rmsd analysis' setup, its
// reference frame). Compared: params, meta, the order of the outputs, the cache key, every array
// (dtype, shape, values), the log lines of the instance and errors. Tolerances (docs/specs/m3-plan.md
// 4.3): exact, except the coordinates of frames superposed onto the first one (trace / atoms with
// align_frames), which may differ by one float32 ulp (kabsch's SVD and BLAS products differ in the
// last bits; for coordinates near zero: 1e-12 of the frame's largest coordinate), and the atoms'
// distance-guessed bonds, whose rows MDBridge sorts (MDATK: kd-tree order; same set).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <set>
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
#include "mdbridge/core/selections.hpp"

using namespace mdbridge;
using namespace mdbridge::analysis;
using namespace mdbridge::test::analysis_cases;

namespace {

constexpr std::string_view kPeriodicMissing = "periodic neighbour search not implemented yet";

bool isStructure(const std::string& instance) {
    const std::string base = parseInstance(instance).base;
    return base == "trace" || base == "proximity" || base == "atoms";
}

struct Counts {
    int instances = 0;  // instances compared
    int errors = 0;  // expected errors reproduced
    int periodicSkipped = 0;  // instances skipped for want of the periodic search
};

std::size_t product(const std::vector<std::size_t>& shape) {
    std::size_t n = 1;
    for (const std::size_t s : shape) {
        n *= s;
    }
    return n;
}

// Rows [m, 2] of an int64 array, sorted.
std::vector<std::pair<std::int64_t, std::int64_t>> sortedRows(std::span<const std::int64_t> v) {
    std::vector<std::pair<std::int64_t, std::int64_t>> rows;
    for (std::size_t k = 0; k + 1 < v.size(); k += 2) {
        rows.emplace_back(v[k], v[k + 1]);
    }
    std::sort(rows.begin(), rows.end());
    return rows;
}

void compareStructureArray(const std::string& name, const core::NumArray& actual, const Value& expected, bool aligned,
                  bool rowsAsSet) {
    SCOPED_TRACE("array " + name);
    ASSERT_EQ(numpyDtype(actual), str(at(expected, "dtype")));
    std::vector<std::size_t> shape;
    for (const Value& s : list(at(expected, "shape"))) {
        shape.push_back(static_cast<std::size_t>(integer(s)));
    }
    ASSERT_EQ(actual.shape(), shape);
    const Array& values = list(at(expected, "values"));
    ASSERT_EQ(values.size(), product(shape)) << "the generator recorded every value";
    switch (actual.dtype()) {
    case core::NumArray::DType::Bool: {
        const auto v = actual.values<std::uint8_t>();
        for (std::size_t k = 0; k < v.size(); ++k) {
            ASSERT_EQ(v[k] != 0, *values[k].get<bool>()) << "at " << k;
        }
        break;
    }
    case core::NumArray::DType::Int64: {
        const auto v = actual.values<std::int64_t>();
        std::vector<std::int64_t> e;
        for (const Value& x : values) {
            e.push_back(integer(x));
        }
        if (rowsAsSet) {
            const auto rows = sortedRows(v);
            std::vector<std::int64_t> flat;
            for (const auto& [i, j] : rows) {
                flat.push_back(i);
                flat.push_back(j);
            }
            EXPECT_EQ(std::vector<std::int64_t>(v.begin(), v.end()), flat) << "MDBridge's bond rows are sorted";
            EXPECT_EQ(rows, sortedRows(e));
        } else {
            EXPECT_EQ(std::vector<std::int64_t>(v.begin(), v.end()), e);
        }
        break;
    }
    case core::NumArray::DType::Float32: {
        const auto v = actual.values<float>();
        // per frame (first axis) the largest |coordinate|, for the near-zero tolerance
        const std::size_t perFrame = shape.empty() || shape[0] == 0 ? 1 : v.size() / shape[0];
        int differ = 0;
        for (std::size_t k = 0; k < v.size(); ++k) {
            const double e = hexFloat(values[k]);
            const double a = v[k];
            if (a == e || (std::isnan(a) && std::isnan(e))) {
                continue;
            }
            const std::size_t frame = perFrame ? k / perFrame : 0;
            if (aligned && frame > 0) {
                double scale = 0.0;
                for (std::size_t q = frame * perFrame; q < (frame + 1) * perFrame; ++q) {
                    scale = std::max(scale, std::abs(hexFloat(values[q])));
                }
                const double diff = std::abs(a - e);
                if (diff <= ulp32(std::max(std::abs(a), std::abs(e))) || diff <= 1e-12 * scale) {
                    ++differ;
                    continue;
                }
            }
            ADD_FAILURE() << "at " << k << ": " << a << " vs MDATK " << e;
            if (++differ > 20) {
                return;
            }
        }
        break;
    }
    case core::NumArray::DType::String: {
        const auto v = actual.values<std::string>();
        for (std::size_t k = 0; k < v.size(); ++k) {
            EXPECT_EQ(v[k], str(values[k])) << "at " << k;
        }
        break;
    }
    default:
        ADD_FAILURE() << "unexpected dtype " << numpyDtype(actual);
    }
}

void compareResult(const core::AnalysisResult& res, const Value& expected, const core::Project& project,
                   const Value& fingerprint) {
    const std::string& instance = str(at(expected, "instance"));
    SCOPED_TRACE("instance " + instance);
    EXPECT_EQ(res.name, instance);
    const AnalysisDescriptor& d = analysisFor(instance);
    EXPECT_EQ(d.cacheKey(project, fingerprint), str(at(expected, "key")));
    EXPECT_EQ(sortedTables(Value(res.params)), sortedTables(untag(at(expected, "params"))));
    EXPECT_EQ(sortedTables(Value(res.meta)), sortedTables(untag(at(expected, "meta"))));
    // the order of the outputs per kind
    const Value meta(res.meta);
    for (const char* kind : {"series", "matrices", "pairmaps"}) {
        std::vector<std::string> names;
        if (const Value* outputs = meta.find("outputs")) {
            if (const Value* table = outputs->find(kind)) {
                for (const auto& [k, v] : *table->get<Table>()) {
                    names.push_back(k);
                }
            }
        }
        std::vector<std::string> want;
        for (const Value& n : list(at(at(expected, "output_order"), kind))) {
            want.push_back(str(n));
        }
        EXPECT_EQ(names, want) << kind;
    }
    // arrays: the same names in the same order
    const Table& arrays = *at(expected, "arrays").get<Table>();
    std::vector<std::string> actualNames;
    for (const auto& [k, v] : res.arrays) {
        actualNames.push_back(k);
    }
    std::vector<std::string> expectedNames;
    for (const auto& [k, v] : arrays) {
        expectedNames.push_back(k);
    }
    std::vector<std::string> sortedActual = actualNames;
    std::sort(sortedActual.begin(), sortedActual.end());
    EXPECT_EQ(sortedActual, expectedNames);  // the generator's JSON sorts the keys
    const std::string base = parseInstance(instance).base;
    const bool aligned = base != "proximity" && *untag(at(expected, "params")).find("align")->get<bool>();
    for (const auto& [name, values] : arrays) {
        const core::NumArray* actual = res.array(name);
        if (!actual) {
            ADD_FAILURE() << "missing array " << name;
            continue;
        }
        compareStructureArray(name, *actual, values, aligned && name.find("coords") != std::string::npos,
                     base == "atoms" && name == "bonds");
    }
}

// One case's run of the structure analyses.
struct StructureRun {
    std::string file;
    const Value* c = nullptr;
    core::Project project;
    Value fingerprint;
    std::vector<std::int64_t> frames;
    const Value* error = nullptr;  // MDATK's error, if any
    std::set<std::string> expectedLog;
};

// Runs one instance of the case in its own pass and compares its result (errors propagate).
void runInstance(const StructureRun& run, const std::string& instance, Counts& counts) {
    SCOPED_TRACE("instance " + instance);
    const Value& c = *run.c;
    core::System system = caseSystem(caseFile(run.file), str(at(c, "system")), std::nullopt);
    std::vector<std::string> log;
    core::Reporter reporter([&](const std::string& m) { log.push_back(m); });
    const core::Project project = instanceProject(run.project, instance);
    AnalysisContext context(system, project, "", reporter);
    const Value* setupFrame = at(c, "setup_frame").find(instance);
    std::vector<std::unique_ptr<Analysis>> analyses;
    analyses.push_back(std::make_unique<SeekBeforeSetup>(analysisFor(instance).create(context, instance),
                                                         setupFrame ? integer(*setupFrame) : 0));
    RunResult result;
    try {
        result = runAnalyses(context, analyses, run.frames);
    } catch (const std::logic_error& e) {
        if (std::string_view(e.what()) == kPeriodicMissing) {
            ++counts.periodicSkipped;
            return;
        }
        throw;
    }
    ASSERT_FALSE(run.error) << "MDATK failed: " << str(at(*run.error, "message"));
    // log: this instance's lines are MDATK's (warnings included)
    for (const std::string& line : log) {
        if (line.rfind("Analysing ", 0) != 0) {
            EXPECT_TRUE(run.expectedLog.contains(line)) << line;
        }
    }
    EXPECT_NE(std::find(log.begin(), log.end(), "Preparing " + analysisFor(instance).label() + " ..."), log.end());
    // times of the analysed frames
    const Array& times = list(at(c, "times"));
    ASSERT_EQ(result.timePs.size(), times.size());
    for (std::size_t k = 0; k < times.size(); ++k) {
        EXPECT_EQ(result.timePs[k], hexFloat(times[k])) << "frame " << k;
    }
    const Value* expected = nullptr;
    for (const Value& r : list(at(c, "results"))) {
        if (str(at(r, "instance")) == instance) {
            expected = &r;
        }
    }
    ASSERT_NE(expected, nullptr);
    ASSERT_EQ(result.results.size(), 1u);
    compareResult(result.results[0], *expected, project, run.fingerprint);
    ++counts.instances;
}

// Runs the structure instances of the cases of `file` selected by `want`; returns what was compared.
template <class Want>
Counts runCases(const std::string& file, Want&& want) {
    Counts counts;
    for (const Value& c : list(at(caseFile(file), "cases"))) {
        const std::string& name = str(at(c, "name"));
        std::vector<std::string> instances;
        bool all = true;
        for (const Value& i : list(at(c, "instances"))) {
            if (isStructure(str(i))) {
                instances.push_back(str(i));
            } else {
                all = false;
            }
        }
        if (instances.empty() || !want(name)) {
            continue;
        }
        SCOPED_TRACE(file + " case " + name);
        StructureRun run{file, &c, caseProject(c), untag(at(c, "fingerprint")), {}, c.find("error"), {}};
        if (run.error && !all) {
            continue;  // an error of another analysis
        }
        if (run.error) {
            run.frames = {0};
        } else {
            for (const Value& f : list(at(c, "frames"))) {
                run.frames.push_back(integer(f));
            }
        }
        for (const Value& line : list(at(c, "log"))) {
            run.expectedLog.insert(str(line));
        }
        try {
            for (const std::string& instance : instances) {
                runInstance(run, instance, counts);
            }
            if (run.error) {
                ADD_FAILURE() << "expected " << str(at(*run.error, "type")) << ": " << str(at(*run.error, "message"));
            }
        } catch (const core::ValueError& e) {
            if (!run.error) {
                ADD_FAILURE() << e.pythonName() << ": " << e.what();
                continue;
            }
            EXPECT_EQ(e.pythonName(), str(at(*run.error, "type")));
            EXPECT_EQ(e.what(), str(at(*run.error, "message")));
            ++counts.errors;
        }
    }
    return counts;
}

bool any(std::string_view) { return true; }

const Value& findCase(const std::string& file, std::string_view name) {
    const Value& data = caseFile(file);
    for (const Value& c : list(at(data, "cases"))) {
        if (str(at(c, "name")) == name) {
            return c;
        }
    }
    throw std::runtime_error("no case " + std::string(name));
}

// `instances`: the instances of the selected cases that run (compared, or skipped while the
// periodic neighbour search is missing); `errors`: the expected errors reproduced.
void report(const Counts& c, int instances, int errors) {
    EXPECT_EQ(c.instances + c.periodicSkipped, instances);
    EXPECT_EQ(c.errors, errors);
    if (c.periodicSkipped) {
        GTEST_SKIP() << c.periodicSkipped << " instance(s) need the periodic neighbour search: " << kPeriodicMissing;
    }
}

}  // namespace

// ---------------------------------------------------------------- micro systems

TEST(StructureAnalyses, MicroSystems) {
    // contacts_micro_all, duplex_all (3 each), duplex_stride_1, structure_residues_only (2),
    // proximity_unknown_set (error)
    report(runCases("micro_cases.json", any), 9, 1);
}

// ---------------------------------------------------------------- demo

TEST(StructureAnalyses, DemoDefault) { report(runCases("demo_dry_cases.json", any), 3, 0); }

TEST(StructureAnalyses, DemoVariants) {
    // every variant but use_pbc: readers (3 x 3), frames (3 + 1), mg_no_ions (2), trace / atoms
    // (2 + 1 + 2 + 2), structure groups (3 + error), proximity sets (2 + error), atoms_hydrogens,
    // instance_no_groups (error)
    report(runCases("demo_dry_variants.json", [](std::string_view n) { return n != "use_pbc"; }), 28, 3);
}

TEST(StructureAnalyses, DemoUsePbc) {
    // trace and atoms do not search neighbours; proximity does (periodic search)
    report(runCases("demo_dry_variants.json", [](std::string_view n) { return n == "use_pbc"; }), 3, 0);
}

// ---------------------------------------------------------------- new small systems

TEST(StructureAnalyses, PeriodicSystemsOpen) {
    report(runCases("systems_cases.json", [](std::string_view n) { return n.find("_open") != std::string_view::npos; }),
           6, 0);
}

TEST(StructureAnalyses, PeriodicSystemsPbc) {
    report(runCases("systems_cases.json",
                    [](std::string_view n) {
                        return n.rfind("periodic", 0) == 0 && n.find("_pbc") != std::string_view::npos;
                    }),
           8, 0);  // trace and atoms of the *_pbc cases; proximity also in *_pbc_cutoff_12
}

TEST(StructureAnalyses, PdbAndGroWithoutBonds) {
    // elements guessed, chain IDs (PDB), bonds guessed from distances (rows compared as a set)
    report(runCases("systems_cases.json",
                    [](std::string_view n) {
                        return (n.rfind("pdb_nobonds", 0) == 0 || n.rfind("gro_nobonds", 0) == 0) &&
                               n != "pdb_nobonds_dcd_pbc";
                    }),
           12, 0);
}

TEST(StructureAnalyses, PdbWithoutBondsPbc) {
    // atoms only: its bond guess never uses the cell
    report(runCases("systems_cases.json", [](std::string_view n) { return n == "pdb_nobonds_dcd_pbc"; }), 1, 0);
}

TEST(StructureAnalyses, HeavyAtomsAndAromatic) {
    report(runCases("systems_cases.json", [](std::string_view n) { return n == "noh" || n == "aromatic"; }), 6, 0);
}

// ---------------------------------------------------------------- helpers

TEST(StructureAnalyses, GuessElement) {
    EXPECT_EQ(structureGuessElement("C4'"), "C");
    EXPECT_EQ(structureGuessElement(" 1HB "), "H");
    EXPECT_EQ(structureGuessElement("HO2'"), "H");
    EXPECT_EQ(structureGuessElement("ca"), "Ca");
    EXPECT_EQ(structureGuessElement("CA1"), "C");
    EXPECT_EQ(structureGuessElement("Mg"), "Mg");
    EXPECT_EQ(structureGuessElement("ZN"), "Zn");
    EXPECT_EQ(structureGuessElement("na"), "Na");
    EXPECT_EQ(structureGuessElement("cl"), "Cl");
    EXPECT_EQ(structureGuessElement("o1p"), "O");
    EXPECT_EQ(structureGuessElement("123"), "X");
    EXPECT_EQ(structureGuessElement(""), "X");
    EXPECT_EQ(structureGuessElement("  "), "X");
}

TEST(StructureAnalyses, SelectionText) {
    core::Project p;
    EXPECT_EQ(structureSelectionText(p), "group:Aptamer or group:Target");
    p.visual.structure_groups = {"Target"};
    EXPECT_EQ(structureSelectionText(p), "group:Target");
    p.visual.structure_groups = {"Nope"};
    try {
        (void)structureSelectionText(p);
        ADD_FAILURE() << "expected SelectionError";
    } catch (const core::SelectionError& e) {
        EXPECT_EQ(std::string(e.what()), "no residue groups available for the structural view");
    }
}

TEST(StructureAnalyses, StructureResiduesLeaveOutIons) {
    const std::string fileName = "demo_dry_cases.json";
    const Value& c = findCase(fileName, "default");
    core::System system = caseSystem(caseFile(fileName), str(at(c, "system")), std::nullopt);
    core::Project project = caseProject(c);
    project.groups.push_back({"Ions", "25-27", "", "#000000"});
    project.visual.structure_groups = {"Aptamer", "Ions"};
    AnalysisContext context(system, project);
    const StructureResidues s = structureResidues(context);
    ASSERT_EQ(s.residues.size(), 12u);
    EXPECT_EQ(s.residues.front(), 0);
    EXPECT_EQ(s.residues.back(), 11);
    EXPECT_EQ(s.atoms.residues(), s.residues);
    // only ions: no residues, and trace / atoms say so
    project.visual.structure_groups = {"Ions"};
    AnalysisContext ions(system, project);
    EXPECT_TRUE(structureResidues(ions).residues.empty());
    for (const AnalysisDescriptor* d : {&traceAnalysis(), &atomsAnalysis()}) {
        std::vector<std::unique_ptr<Analysis>> analyses;
        analyses.push_back(d->create(ions));
        const std::vector<std::int64_t> frames{0};
        try {
            (void)runAnalyses(ions, analyses, frames);
            ADD_FAILURE() << d->name() << ": expected SelectionError";
        } catch (const core::SelectionError& e) {
            EXPECT_EQ(std::string(e.what()), "structural view: no residues selected");
        }
    }
}
