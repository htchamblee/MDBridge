#pragma once
// Running the recorded cases of tests/parity/make_analysis_cases.py (data/analysis/micro_cases.json,
// demo_dry_cases.json, demo_dry_variants.json, systems_cases.json) for the rmsd analysis' tests:
// the systems (micro systems built from their JSON, files loaded with io::loadSystem), the project
// and its instances, and one run of runAnalyses as MDATK's Pipeline.analyze builds it.
//
// Instances whose analysis is not under test run as probes: an analysis with the same descriptor
// (label, params) that computes nothing. Every instance records the trajectory frame current at its
// setup (the cases' `setup_frame`), so the frame state across setups is checked for whole runs.
//
// To be unified with the other analyses' shared support header at merge time.

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"
#include "mdbridge/io/loader.hpp"
#include "test_support.hpp"

namespace mdbridge::test::rmsd_support {

using core::Array;
using core::Table;
using core::Value;

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
inline bool isNull(const Value& value) { return value.isNull(); }

// float.hex() text, "nan", "inf", "-inf"
inline double hexDouble(const Value& value) {
    const std::string& text = str(value);
    if (text == "nan") {
        return std::nan("");
    }
    if (text == "inf") {
        return HUGE_VAL;
    }
    if (text == "-inf") {
        return -HUGE_VAL;
    }
    return std::strtod(text.c_str(), nullptr);
}

inline std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    for (const Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}

// One case file, parsed once.
inline const Value& caseFile(const std::string& name) {
    static std::map<std::string, Value> files;
    auto it = files.find(name);
    if (it == files.end()) {
        it = files.emplace(name, readJson(kData / "analysis" / name)).first;
    }
    return it->second;
}

// "{data}/x" -> tests/data/x, "{out}/x" -> tests/data/analysis/x
inline std::filesystem::path resolvePlaceholder(const std::string& path) {
    if (path.rfind("{data}/", 0) == 0) {
        return kData / path.substr(7);
    }
    if (path.rfind("{out}/", 0) == 0) {
        return kData / "analysis" / path.substr(6);
    }
    throw std::runtime_error("unknown placeholder in " + path);
}

// conftest.make_universe: residues in order of first appearance of their resid, segment SYS, no
// types, the given elements / masses / charges, bonds only when listed, float32 frames in memory.
inline core::System microSystem(const Value& system) {
    auto t = std::make_shared<core::Topology>();
    std::vector<std::int64_t> resids;
    std::vector<double> masses;
    std::vector<double> charges;
    std::vector<std::uint32_t> elements;
    for (const Value& atom : list(at(system, "atoms"))) {
        const std::int64_t resid = integer(at(atom, "resid"));
        std::size_t residue = 0;
        while (residue < resids.size() && resids[residue] != resid) {
            ++residue;
        }
        if (residue == resids.size()) {
            resids.push_back(resid);
            t->resids.push_back(resid);
            t->resnames.push_back(t->strings.intern(str(at(atom, "resname"))));
            t->residueSegment.push_back(0);
        }
        t->names.push_back(t->strings.intern(str(at(atom, "name"))));
        t->atomResidue.push_back(static_cast<std::int32_t>(residue));
        elements.push_back(t->strings.intern(str(at(atom, "element"))));
        masses.push_back(hexDouble(at(atom, "mass")));
        charges.push_back(hexDouble(at(atom, "charge")));
    }
    t->elements = std::move(elements);
    t->masses = std::move(masses);
    t->charges = std::move(charges);
    t->segids.push_back(t->strings.intern(str(at(system, "segid"))));
    const Array& bonds = list(at(system, "bonds"));
    if (!bonds.empty()) {
        std::vector<core::Bond> out;
        for (const Value& b : bonds) {
            out.push_back({static_cast<std::int32_t>(integer(list(b)[0])), static_cast<std::int32_t>(integer(list(b)[1]))});
        }
        t->bonds = std::move(out);
    }
    t->finish();
    std::vector<core::Frame> frames;
    for (const Value& coords : list(at(system, "coords"))) {
        core::Frame f;
        for (const Value& v : list(coords)) {
            f.positions.push_back(static_cast<float>(hexDouble(v)));
        }
        frames.push_back(std::move(f));
    }
    return core::System{t, std::make_unique<core::MemoryTrajectory>(std::move(frames), hexDouble(at(system, "dt")))};
}

struct LoadedSystem {
    core::System system;
    std::string topologyPath;
};

// The system of a case: a micro system, or files loaded as load_universe loads them.
inline LoadedSystem loadCaseSystem(const Value& file, const Value& c, std::optional<double> dtPs) {
    const std::string& name = str(at(c, "system"));
    const Value& system = at(at(file, "systems"), name);
    if (!system.find("topology")) {
        return {microSystem(system), ""};
    }
    const std::filesystem::path topology = resolvePlaceholder(str(at(system, "topology")));
    std::vector<std::filesystem::path> trajectories;
    for (const Value& t : list(at(system, "trajectories"))) {
        trajectories.push_back(resolvePlaceholder(str(t)));
    }
    return {io::loadSystem(topology, trajectories, dtPs), topology.string()};
}

// An instance of the run: the real analysis (when under test) or a probe, recording the frame
// current at setup.
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

struct CaseRun {
    core::Project project;
    std::vector<std::string> instances;
    std::vector<core::Project> projects;  // per instance
    std::vector<std::int64_t> frames;
    std::vector<std::string> log;
    std::map<std::string, std::int64_t> setupFrames;
    std::optional<std::pair<std::string, std::string>> error;  // {type, message}
    analysis::RunResult run;
    std::vector<bool> real;  // per instance: ran the analysis under test
};

// The project of a case (its recorded sections read as a project file's plain data).
inline core::Project caseProject(const Value& c) {
    std::vector<std::string> warnings;
    return core::projectFromPlain(at(c, "project"), {}, &warnings);
}

// One run of a case. `underTest(base)`: the analyses created for real; the others are probes.
inline CaseRun runCase(const Value& file, const Value& c, const std::function<bool(const std::string&)>& underTest) {
    CaseRun out;
    out.project = caseProject(c);
    const core::InputParams& inputs = out.project.inputs;
    LoadedSystem loaded = loadCaseSystem(file, c, inputs.dt_ps);
    const core::Reporter reporter([&out](const std::string& m) { out.log.push_back(m); });
    out.instances = strings(at(c, "instances"));
    std::vector<std::unique_ptr<analysis::AnalysisContext>> contexts;
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    try {
        out.frames = io::frameIndices(loaded.system.trajectory->frameCount(), inputs.start, inputs.stop, inputs.step);
        for (const std::string& instance : out.instances) {
            out.projects.push_back(analysis::instanceProject(out.project, instance));
        }
        for (std::size_t i = 0; i < out.instances.size(); ++i) {
            const std::string& instance = out.instances[i];
            const analysis::AnalysisDescriptor& d = analysis::analysisFor(instance);
            contexts.push_back(std::make_unique<analysis::AnalysisContext>(loaded.system, out.projects[i],
                                                                           loaded.topologyPath, reporter));
            std::unique_ptr<analysis::Analysis> inner;
            if (underTest(d.name())) {
                inner = d.create(*contexts.back(), instance);
            }
            out.real.push_back(inner != nullptr);
            analyses.push_back(std::make_unique<Recorded>(d, *contexts.back(), instance, std::move(inner),
                                                          out.setupFrames));
        }
        out.run = analysis::runAnalyses(*contexts.front(), analyses, out.frames, &reporter);
    } catch (const core::ValueError& e) {
        out.error = {e.pythonName(), e.what()};
    } catch (const core::NamedError& e) {
        out.error = {std::string(e.pythonName()), e.what()};
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
// analyses (warnings whose text does not start with `prefix`).
inline std::vector<std::string> expectedLog(const Value& c, std::string_view prefix) {
    std::vector<std::string> out;
    for (const std::string& line : strings(at(c, "log"))) {
        if (line.rfind("WARNING: ", 0) == 0 && line.compare(9, prefix.size(), prefix) != 0) {
            continue;
        }
        out.push_back(line);
    }
    return out;
}

inline bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

}  // namespace mdbridge::test::rmsd_support
