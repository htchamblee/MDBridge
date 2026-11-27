#pragma once
// Shared helpers of the analysis tests that replay the cases of tests/parity/make_analysis_cases.py
// (data/analysis/micro_cases.json, demo_dry_cases.json, demo_dry_variants.json, systems_cases.json;
// format in data/README.md and the generator's docstring):
//
//   * JSON accessors, float.hex() parsing and the generator's tagged JSON ({"$float": ...});
//   * microSystem(): the core::System of a micro system recorded as atoms + float32 frames (what
//     MDATK's tests/conftest.make_universe builds), caseSystem(): the system of any case (micro, or
//     topology / trajectory files loaded with io::loadSystem);
//   * runCase(): one case run as the generator runs it (instance projects, cache keys, analyses
//     created in order, the frames of inputs.start / stop / step, runAnalyses with one reporter),
//     restricted to the instances a test covers, catching the error as {type, message};
//     runCaseWithProbes(): every instance of a case, those not under test as probes, recording the
//     frame current at each setup;
//   * wrappers around an analysis: SeekBeforeSetup (the trajectory moved to MDATK's setup frame)
//     and Recorded (a probe or a real analysis, recording its setup frame); both pass prepare()
//     and concurrentProcess() through to the wrapped analysis;
//   * comparison of arrays ({dtype, shape, values}), of tagged JSON (tables compared without regard
//     to key order), of log lines restricted to the covered analyses, and the decision margins.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/topology.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/loader.hpp"
#include "test_support.hpp"

namespace mdbridge::test::analysis_cases {

using core::Array;
using core::Table;
using core::Value;

// ---- JSON accessors --------------------------------------------------------------------------------

inline const Value& at(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
inline const Array& list(const Value& value) { return *value.get<Array>(); }
inline const std::string& str(const Value& value) { return *value.get<std::string>(); }
inline std::int64_t integer(const Value& value) { return *value.get<std::int64_t>(); }

// float.hex() text, "inf" / "-inf" / "nan" included
inline double hexText(std::string_view text) {
    if (text == "inf") {
        return std::numeric_limits<double>::infinity();
    }
    if (text == "-inf") {
        return -std::numeric_limits<double>::infinity();
    }
    if (text == "nan") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::strtod(std::string(text).c_str(), nullptr);
}
inline double hexFloat(const Value& value) { return hexText(str(value)); }

inline std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    for (const Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}
inline std::vector<std::int64_t> integers(const Value& value) {
    std::vector<std::int64_t> out;
    for (const Value& v : list(value)) {
        out.push_back(integer(v));
    }
    return out;
}
inline std::vector<double> hexFloats(const Value& value) {
    std::vector<double> out;
    for (const Value& v : list(value)) {
        out.push_back(hexFloat(v));
    }
    return out;
}

// Doubles equal bit for bit up to NaN payloads: NaN equals NaN, 0.0 differs from -0.0.
inline bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

// One float32 ulp of |x| (the spacing above it)
inline double ulp32(double x) {
    const float f = static_cast<float>(std::abs(x));
    return static_cast<double>(std::nextafter(f, std::numeric_limits<float>::infinity()) - f);
}

// The generator's tagged JSON: {"$float": "nan" | "inf" | "-inf"} stands for a non-finite float.
inline Value untag(const Value& value) {
    if (const Table* t = value.get<Table>()) {
        if (t->size() == 1 && t->front().first == "$float") {
            return hexText(str(t->front().second));
        }
        Table out;
        for (const auto& [k, v] : *t) {
            out.emplace_back(k, untag(v));
        }
        return out;
    }
    if (const Array* a = value.get<Array>()) {
        Array out;
        for (const Value& v : *a) {
            out.push_back(untag(v));
        }
        return out;
    }
    return value;
}

// Tables with their entries sorted by key, at every depth (the generator writes sorted keys).
inline Value sortedTables(const Value& value) {
    if (const Table* t = value.get<Table>()) {
        Table out;
        for (const auto& [k, v] : *t) {
            out.emplace_back(k, sortedTables(v));
        }
        std::stable_sort(out.begin(), out.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
        return out;
    }
    if (const Array* a = value.get<Array>()) {
        Array out;
        for (const Value& v : *a) {
            out.push_back(sortedTables(v));
        }
        return out;
    }
    return value;
}

// Tagged JSON from the data against a value computed here, key order ignored.
inline void expectSameValue(const Value& expected, const Value& actual, const std::string& what) {
    EXPECT_EQ(core::toJson(sortedTables(untag(expected)), -1), core::toJson(sortedTables(actual), -1)) << what;
}

// ---- case files ------------------------------------------------------------------------------------

inline const Value& caseFile(const std::string& file) {
    static std::vector<std::pair<std::string, std::unique_ptr<Value>>> cache;
    for (const auto& [name, value] : cache) {
        if (name == file) {
            return *value;
        }
    }
    cache.emplace_back(file, std::make_unique<Value>(readJson(kData / "analysis" / file)));
    return *cache.back().second;
}

// "{data}/..." (tests/data) and "{out}/..." (tests/data/analysis) of the generator's paths
inline std::filesystem::path dataPath(std::string_view path) {
    const auto replace = [&](std::string_view prefix, const std::filesystem::path& dir) -> std::optional<std::filesystem::path> {
        if (path.substr(0, prefix.size()) == prefix) {
            return dir / std::filesystem::path(std::string(path.substr(prefix.size())));
        }
        return std::nullopt;
    };
    if (auto p = replace("{data}/", kData)) {
        return *p;
    }
    if (auto p = replace("{out}/", kData / "analysis")) {
        return *p;
    }
    return std::filesystem::path(std::string(path));
}

// The generator's placeholders put back into a message: "<tests/data/analysis>/x" -> "{out}/x",
// "<tests/data>/x" -> "{data}/x" (either separator after the folder).
inline std::string withPlaceholders(std::string text) {
    const auto put = [&text](const std::filesystem::path& dir, std::string_view placeholder) {
        for (const std::string& prefix : {core::pathToUtf8(dir), dir.generic_string()}) {
            for (const char sep : {'/', '\\'}) {
                const std::string from = prefix + sep;
                for (std::size_t at = text.find(from); at != std::string::npos; at = text.find(from, at)) {
                    text.replace(at, from.size(), std::string(placeholder) + "/");
                }
            }
        }
    };
    put(kData / "analysis", "{out}");
    put(kData, "{data}");
    return text;
}

// ---- systems ---------------------------------------------------------------------------------------

// One atom of a micro system (conftest.make_universe's atom dicts; MDATK's test_analyses.atom).
struct MicroAtom {
    std::string name;
    std::int64_t resid = 1;
    std::string resname;
    std::string element;
    double mass = 12.0;
    double charge = 0.0;
};

// conftest.make_universe(atoms, coords, bonds, dt): residues in order of first appearance of their
// resid (resname of the first atom), one segment `segid`, names / elements / masses / charges per
// atom, no types (""), bonds only when some are given (else none at all: bonds unknown), and an
// in-memory trajectory of `frames` (x, y, z per atom, rounded to float32) without cells, frame k at
// time k * dt.
inline core::System makeMicroSystem(const std::vector<MicroAtom>& atoms, const std::vector<std::vector<double>>& frames,
                                    const std::vector<std::pair<std::int32_t, std::int32_t>>& bonds = {},
                                    double dt = 1.0, const std::string& segid = "SYS") {
    auto t = std::make_shared<core::Topology>();
    std::vector<std::int64_t> resids;
    std::vector<std::uint32_t> elements;
    std::vector<double> masses, charges;
    for (const MicroAtom& atom : atoms) {
        auto it = std::find(resids.begin(), resids.end(), atom.resid);
        if (it == resids.end()) {
            resids.push_back(atom.resid);
            t->resnames.push_back(t->strings.intern(atom.resname));
            it = resids.end() - 1;
        }
        t->names.push_back(t->strings.intern(atom.name));
        t->types.push_back(t->strings.intern(""));
        t->atomResidue.push_back(static_cast<std::int32_t>(it - resids.begin()));
        elements.push_back(t->strings.intern(atom.element));
        masses.push_back(atom.mass);
        charges.push_back(atom.charge);
    }
    t->elements = std::move(elements);
    t->masses = std::move(masses);
    t->charges = std::move(charges);
    t->resids = resids;
    t->resnums = resids;
    t->residueSegment.assign(resids.size(), 0);
    t->segids.push_back(t->strings.intern(segid));
    if (!bonds.empty()) {
        std::vector<core::Bond> pairs;
        for (const auto& [x, y] : bonds) {
            pairs.push_back({std::min(x, y), std::max(x, y)});
        }
        t->bonds = std::move(pairs);
    }
    t->finish();
    std::vector<core::Frame> list_;
    for (const std::vector<double>& coords : frames) {
        core::Frame f;
        for (const double v : coords) {
            f.positions.push_back(static_cast<float>(v));
        }
        list_.push_back(std::move(f));
    }
    return core::System{std::move(t), std::make_unique<core::MemoryTrajectory>(std::move(list_), dt)};
}

// A micro system as micro_cases.json records it ({atoms, bonds, coords, dt, segid}).
inline core::System microSystem(const Value& system) {
    std::vector<MicroAtom> atoms;
    for (const Value& atom : list(at(system, "atoms"))) {
        atoms.push_back({str(at(atom, "name")), integer(at(atom, "resid")), str(at(atom, "resname")),
                         str(at(atom, "element")), hexFloat(at(atom, "mass")), hexFloat(at(atom, "charge"))});
    }
    std::vector<std::pair<std::int32_t, std::int32_t>> bonds;
    for (const Value& bond : list(at(system, "bonds"))) {
        bonds.emplace_back(static_cast<std::int32_t>(integer(list(bond)[0])), static_cast<std::int32_t>(integer(list(bond)[1])));
    }
    std::vector<std::vector<double>> frames;
    for (const Value& coords : list(at(system, "coords"))) {
        std::vector<double> xyz;
        for (const Value& v : list(coords)) {
            xyz.push_back(hexFloat(v));
        }
        frames.push_back(std::move(xyz));
    }
    return makeMicroSystem(atoms, frames, bonds, hexFloat(at(system, "dt")), str(at(system, "segid")));
}

// conftest.micro_project(groups): the default project with residue groups {name: residues}.
inline core::Project microProject(const std::vector<std::pair<std::string, std::string>>& groups) {
    core::Project p;
    p.groups.clear();  // conftest replaces the default groups
    for (const auto& [name, residues] : groups) {
        p.groups.push_back({name, residues, "", "#2a78d6"});
    }
    return p;
}

// conftest.run_one(name, universe, project): one analysis over every frame.
inline core::AnalysisResult runOne(std::string_view name, core::System& system, const core::Project& project,
                                   const std::string& topologyPath = {}) {
    analysis::AnalysisContext context(system, project, topologyPath);
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    analyses.push_back(analysis::analysisFor(name).create(context));
    std::vector<std::int64_t> frames(static_cast<std::size_t>(system.trajectory->frameCount()));
    for (std::size_t k = 0; k < frames.size(); ++k) {
        frames[k] = static_cast<std::int64_t>(k);
    }
    return std::move(analysis::runAnalyses(context, analyses, frames).results.front());
}

// The system named `name` of a case file: a micro system (atoms) or files ({topology,
// trajectories}), loaded with the case's inputs.dt_ps.
inline core::System caseSystem(const Value& file, const std::string& name, std::optional<double> dtPs) {
    const Value& system = at(at(file, "systems"), name);
    if (system.find("atoms")) {
        return microSystem(system);
    }
    std::vector<std::filesystem::path> trajectories;
    for (const Value& t : list(at(system, "trajectories"))) {
        trajectories.push_back(dataPath(str(t)));
    }
    return io::loadSystem(dataPath(str(at(system, "topology"))), trajectories, dtPs);
}

// The topology file of the system `name` of a case file ("" for a micro system: MDATK's run_one
// passes none).
inline std::string caseTopologyPath(const Value& file, const std::string& name) {
    if (const Value* top = at(at(file, "systems"), name).find("topology")) {
        return dataPath(str(*top)).string();
    }
    return {};
}

// The project of a case, its tagged JSON read as a project file's plain data.
inline core::Project caseProject(const Value& c) { return core::projectFromPlain(untag(at(c, "project"))); }

// ---- running a case --------------------------------------------------------------------------------

struct CaseError {
    std::string type;
    std::string message;
};

struct CaseRun {
    std::vector<std::string> instances;  // the instances run (those the predicate kept)
    std::vector<std::string> keys;  // cache key per instance
    std::vector<std::int64_t> frames;
    std::vector<std::string> logs;
    std::vector<std::pair<double, std::string>> progress;
    analysis::RunResult run;
    std::optional<CaseError> error;
};

// {type, message} as the generator records a Python exception
inline CaseError caseError(const std::exception& e) {
    if (const auto* v = dynamic_cast<const core::ValueError*>(&e)) {
        return {v->pythonName(), v->what()};
    }
    if (const auto* n = dynamic_cast<const core::NamedError*>(&e)) {
        return {n->pythonName(), n->what()};
    }
    if (dynamic_cast<const std::out_of_range*>(&e)) {
        return {"KeyError", core::pythonRepr(Value(std::string(e.what())))};  // str(KeyError) is the repr
    }
    return {"C++", e.what()};
}

// The instances of a case a test covers: those whose analysis `covered` accepts, and unknown
// analyses (their error is part of the case).
inline std::vector<std::string> coveredInstances(const Value& c, const std::function<bool(const std::string&)>& covered) {
    std::vector<std::string> out;
    for (const std::string& instance : strings(at(c, "instances"))) {
        const std::string base = analysis::parseInstance(instance).base;
        if (covered(base) || !analysis::findAnalysis(base)) {
            out.push_back(instance);
        }
    }
    return out;
}

// One case as the generator runs it, restricted to `instances` (runAnalyses with `options`).
inline CaseRun runCase(const Value& file, const Value& c, const std::vector<std::string>& instances,
                       const analysis::RunOptions& options = {}) {
    CaseRun out;
    out.instances = instances;
    const core::Project project = caseProject(c);
    const core::Reporter reporter([&](const std::string& m) { out.logs.push_back(m); },
                                  [&](double f, const std::string& m) { out.progress.emplace_back(f, m); });
    const Value fingerprint = untag(at(c, "fingerprint"));
    try {
        core::System system = caseSystem(file, str(at(c, "system")), project.inputs.dt_ps);
        out.frames = io::frameIndices(system.trajectory->frameCount(), project.inputs.start, project.inputs.stop,
                                      project.inputs.step);
        std::vector<core::Project> projects;
        for (const std::string& instance : instances) {
            projects.push_back(analysis::instanceProject(project, instance));
        }
        for (std::size_t n = 0; n < instances.size(); ++n) {
            out.keys.push_back(analysis::analysisFor(instances[n]).cacheKey(projects[n], fingerprint));
        }
        const std::string topologyPath = caseTopologyPath(file, str(at(c, "system")));
        std::vector<std::unique_ptr<analysis::AnalysisContext>> contexts;
        std::vector<std::unique_ptr<analysis::Analysis>> analyses;
        for (std::size_t n = 0; n < instances.size(); ++n) {
            contexts.push_back(std::make_unique<analysis::AnalysisContext>(system, projects[n], topologyPath, reporter));
            analyses.push_back(analysis::analysisFor(instances[n]).create(*contexts.back(), instances[n]));
        }
        if (!analyses.empty()) {
            out.run = analysis::runAnalyses(*contexts.front(), analyses, out.frames, &reporter, options);
        }
    } catch (const std::exception& e) {
        out.error = caseError(e);
    }
    return out;
}

// The case's log lines restricted to the covered instances: "Preparing <label> ..." of their
// labels, warnings whose message starts with one of `warningPrefixes`, and the "Analysing ..."
// line with the covered labels. The whole log when every instance is covered.
inline std::vector<std::string> expectedLog(const Value& c, const std::vector<std::string>& instances,
                                            const std::vector<std::string>& warningPrefixes) {
    const std::vector<std::string> all = strings(at(c, "log"));
    if (instances.size() == list(at(c, "instances")).size()) {
        return all;
    }
    std::vector<std::string> labels;
    for (const std::string& instance : instances) {
        if (const analysis::AnalysisDescriptor* d = analysis::findAnalysis(analysis::parseInstance(instance).base)) {
            labels.push_back(d->label());
        }
    }
    std::vector<std::string> out;
    for (const std::string& line : all) {
        if (line.rfind("Analysing ", 0) == 0) {
            std::string text = line.substr(0, line.find(": ") + 2);
            for (std::size_t n = 0; n < labels.size(); ++n) {
                text += (n ? ", " : "") + labels[n];
            }
            out.push_back(text);
            continue;
        }
        bool keep = std::any_of(labels.begin(), labels.end(), [&](const std::string& l) { return line == "Preparing " + l + " ..."; });
        for (const std::string& prefix : warningPrefixes) {
            keep = keep || line.rfind("WARNING: " + prefix, 0) == 0;
        }
        if (keep) {
            out.push_back(line);
        }
    }
    return out;
}

// ---- wrapped analyses ------------------------------------------------------------------------------

// Runs `inner` with the trajectory moved to `setupFrame` before its setup, as MDATK's run left it
// (the rmsd analysis' setup seeks its reference frame; the atoms analysis guesses bonds from the
// frame current at its setup).
class SeekBeforeSetup final : public analysis::Analysis {
public:
    SeekBeforeSetup(std::unique_ptr<analysis::Analysis> inner, std::int64_t setupFrame)
        : Analysis(inner->descriptor(), inner->context(), inner->instance()), m_inner(std::move(inner)),
          m_setupFrame(setupFrame) {}

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        context().system().trajectory->seek(m_setupFrame);
        m_inner->setup(frames);
    }
    void prepare(const core::Frame& first) override { m_inner->prepare(first); }
    bool concurrentProcess() const override { return m_inner->concurrentProcess(); }
    void process(std::size_t k, const core::Frame& frame) override { m_inner->process(k, frame); }
    core::AnalysisResult finalize() override { return m_inner->finalize(); }

private:
    std::unique_ptr<analysis::Analysis> m_inner;
    std::int64_t m_setupFrame;
};

// An instance of a run: the real analysis `inner`, or a probe (no `inner`: an analysis with the
// descriptor's label and params that computes nothing); either records the trajectory frame
// current at its setup in `setupFrames` (-1 when none).
class Recorded final : public analysis::Analysis {
public:
    Recorded(const analysis::AnalysisDescriptor& d, analysis::AnalysisContext& context, std::string instance,
             std::unique_ptr<analysis::Analysis> inner, std::map<std::string, std::int64_t>& setupFrames)
        : Analysis(d, context, std::move(instance)), m_inner(std::move(inner)), m_setupFrames(&setupFrames) {}

    void setup(std::span<const std::int64_t> frames) override {
        const core::Frame* current = context().system().currentFrame();
        (*m_setupFrames)[instance()] = current ? current->index : -1;
        Analysis::setup(frames);
        if (m_inner) {
            m_inner->setup(frames);
        }
    }
    void prepare(const core::Frame& first) override {
        if (m_inner) {
            m_inner->prepare(first);
        }
    }
    bool concurrentProcess() const override { return !m_inner || m_inner->concurrentProcess(); }
    void process(std::size_t k, const core::Frame& frame) override {
        if (m_inner) {
            m_inner->process(k, frame);
        }
    }
    core::AnalysisResult finalize() override { return m_inner ? m_inner->finalize() : newResult(); }
    bool real() const { return m_inner != nullptr; }

private:
    std::unique_ptr<analysis::Analysis> m_inner;
    std::map<std::string, std::int64_t>* m_setupFrames;
};

// ---- running a case with probes --------------------------------------------------------------------

struct ProbedRun {
    core::Project project;
    std::vector<std::string> instances;  // every instance of the case
    std::vector<core::Project> projects;  // per instance
    std::vector<std::int64_t> frames;
    std::vector<std::string> log;
    std::map<std::string, std::int64_t> setupFrames;
    std::optional<CaseError> error;
    analysis::RunResult run;
    std::vector<bool> real;  // per instance: ran the analysis under test
};

// One run of every instance of a case as MDATK's Pipeline.analyze builds it: the analyses whose
// name `underTest` accepts are created for real, the others run as probes (Recorded), so the frame
// state across setups is that of the whole run. The project is read without untagging and with
// its warnings collected (not logged); a ValueError / NamedError is caught as {type, message}.
inline ProbedRun runCaseWithProbes(const Value& file, const Value& c,
                                   const std::function<bool(const std::string&)>& underTest) {
    ProbedRun out;
    std::vector<std::string> warnings;
    out.project = core::projectFromPlain(at(c, "project"), {}, &warnings);
    const core::InputParams& inputs = out.project.inputs;
    const std::string& systemName = str(at(c, "system"));
    core::System system = caseSystem(file, systemName, inputs.dt_ps);
    const std::string topologyPath = caseTopologyPath(file, systemName);
    const core::Reporter reporter([&out](const std::string& m) { out.log.push_back(m); });
    out.instances = strings(at(c, "instances"));
    std::vector<std::unique_ptr<analysis::AnalysisContext>> contexts;
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    try {
        out.frames = io::frameIndices(system.trajectory->frameCount(), inputs.start, inputs.stop, inputs.step);
        for (const std::string& instance : out.instances) {
            out.projects.push_back(analysis::instanceProject(out.project, instance));
        }
        for (std::size_t i = 0; i < out.instances.size(); ++i) {
            const std::string& instance = out.instances[i];
            const analysis::AnalysisDescriptor& d = analysis::analysisFor(instance);
            contexts.push_back(
                std::make_unique<analysis::AnalysisContext>(system, out.projects[i], topologyPath, reporter));
            std::unique_ptr<analysis::Analysis> inner;
            if (underTest(d.name())) {
                inner = d.create(*contexts.back(), instance);
            }
            out.real.push_back(inner != nullptr);
            analyses.push_back(
                std::make_unique<Recorded>(d, *contexts.back(), instance, std::move(inner), out.setupFrames));
        }
        out.run = analysis::runAnalyses(*contexts.front(), analyses, out.frames, &reporter);
    } catch (const core::ValueError& e) {
        out.error = caseError(e);
    } catch (const core::NamedError& e) {
        out.error = caseError(e);
    }
    return out;
}

// The recorded result of `instance` in a case, or nullptr.
inline const Value* recordedResult(const Value& c, const std::string& instance) {
    const Value* results = c.find("results");
    if (!results) {
        return nullptr;
    }
    for (const Value& r : list(*results)) {
        if (str(at(r, "instance")) == instance) {
            return &r;
        }
    }
    return nullptr;
}

// The log lines of a case that a run with probes reproduces: all but the warnings of the other
// analyses (warnings whose text does not start with `prefix`, and the energy analysis' parameter line).
// The runner's own frame warnings (AN-19, AN-20) stay: every run gives them.
inline std::vector<std::string> logWithoutOtherWarnings(const Value& c, std::string_view prefix) {
    std::vector<std::string> out;
    const auto runners = [](const std::string& line) {
        return line.rfind("WARNING: frame ", 0) == 0 || line.rfind("WARNING: frames ", 0) == 0 ||
               line.rfind("WARNING: use_pbc is on but trajectory frame ", 0) == 0;
    };
    for (const std::string& line : strings(at(c, "log"))) {
        if (line.rfind("WARNING: ", 0) == 0 && line.compare(9, prefix.size(), prefix) != 0 && !runners(line)) {
            continue;
        }
        // the energy analysis' parameter line (CIT-04), the one log line of an analysis that is no warning
        if (prefix != "energy" && line.rfind("energy: Lennard-Jones parameters", 0) == 0) {
            continue;
        }
        out.push_back(line);
    }
    return out;
}

// ---- arrays and margins ----------------------------------------------------------------------------

// numpy's dtype.str of an array
inline std::string numpyDtype(const core::NumArray& a) {
    switch (a.dtype()) {
    case core::NumArray::DType::Bool: return "|b1";
    case core::NumArray::DType::Int8: return "|i1";
    case core::NumArray::DType::Int32: return "<i4";
    case core::NumArray::DType::Int64: return "<i8";
    case core::NumArray::DType::Float32: return "<f4";
    case core::NumArray::DType::Float64: return "<f8";
    case core::NumArray::DType::String: return std::string(a.dtypeName());
    }
    return "?";
}

// The values of an integer array as int64 (empty for other types).
inline std::vector<std::int64_t> integerValues(const core::NumArray& a) {
    std::vector<std::int64_t> out;
    switch (a.dtype()) {
    case core::NumArray::DType::Int8:
        for (const std::int8_t v : a.values<std::int8_t>()) {
            out.push_back(v);
        }
        break;
    case core::NumArray::DType::Int32:
        for (const std::int32_t v : a.values<std::int32_t>()) {
            out.push_back(v);
        }
        break;
    case core::NumArray::DType::Int64:
        for (const std::int64_t v : a.values<std::int64_t>()) {
            out.push_back(v);
        }
        break;
    default:
        break;
    }
    return out;
}

struct ArrayDiff {
    bool sameType = true;  // dtype and shape
    std::vector<std::size_t> differing;  // flat positions whose values differ
};

// An array of the data ({dtype, shape, values}) against one computed here: dtype and shape, and every
// value exactly (floats bit for bit as float.hex() gives them; NaN equals NaN).
inline ArrayDiff compareArray(const Value& expected, const core::NumArray& actual) {
    ArrayDiff diff;
    std::vector<std::size_t> shape;
    for (const Value& s : list(at(expected, "shape"))) {
        shape.push_back(static_cast<std::size_t>(integer(s)));
    }
    if (str(at(expected, "dtype")) != numpyDtype(actual) || shape != actual.shape()) {
        diff.sameType = false;
        return diff;
    }
    const Array& values = list(at(expected, "values"));
    if (actual.dtype() == core::NumArray::DType::String) {
        const auto got = actual.values<std::string>();
        for (std::size_t n = 0; n < values.size(); ++n) {
            if (str(values[n]) != got[n]) {
                diff.differing.push_back(n);
            }
        }
        return diff;
    }
    const std::vector<double> got = actual.toDouble();
    for (std::size_t n = 0; n < values.size(); ++n) {
        const Value& v = values[n];
        double want = 0.0;
        if (const bool* b = v.get<bool>()) {
            want = *b ? 1.0 : 0.0;
        } else if (const std::int64_t* i = v.get<std::int64_t>()) {
            want = static_cast<double>(*i);
        } else {
            want = hexFloat(v);
        }
        if (!sameDouble(want, got[n])) {
            diff.differing.push_back(n);
        }
    }
    return diff;
}

// The smallest recorded margin of a result over its quantities and frames (+inf when none):
// margins.frames.<quantity> is a list per analysed frame of float.hex / "inf" / null.
inline double smallestMargin(const Value& result) {
    double smallest = std::numeric_limits<double>::infinity();
    const Value* frames = at(result, "margins").find("frames");
    if (!frames) {
        return smallest;
    }
    for (const auto& [quantity, values] : *frames->get<Table>()) {
        for (const Value& v : list(values)) {
            if (v.isString()) {
                smallest = std::min(smallest, hexFloat(v));
            }
        }
    }
    return smallest;
}

// The smallest |margin| of analysed frame k over every quantity of a result's `margins` (+inf when
// nothing was recorded for the frame).
inline double frameMargin(const Value& margins, std::size_t k) {
    double out = std::numeric_limits<double>::infinity();
    const Value* frames = margins.find("frames");
    if (!frames) {
        return out;
    }
    for (const auto& [name, values] : *frames->get<Table>()) {
        const Value& v = list(values).at(k);
        if (!v.isNull()) {
            out = std::min(out, std::abs(hexFloat(v)));
        }
    }
    return out;
}

// A count that differs where the smallest margin is below this (and not an exact tie, margin 0) may
// be a last-bit threshold flip (plan 4.3: relative 1e-9 of thresholds of a few Å / cosines).
inline constexpr double kFlipMargin = 1e-8;

}  // namespace mdbridge::test::analysis_cases
