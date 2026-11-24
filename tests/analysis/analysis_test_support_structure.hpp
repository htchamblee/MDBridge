#pragma once
// Helpers of the structure analyses' tests (to be unified with analysis_test_support.hpp): the
// cases of tests/parity/make_analysis_cases.py (data/analysis, see data/README.md) as systems,
// projects and analysis runs, and the comparison of results with MDATK's.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/loader.hpp"
#include "test_support.hpp"

namespace mdbridge::test::structure {

using core::Array;
using core::Table;
using core::Value;

// One of the case files of data/analysis, read once.
inline const Value& caseFile(std::string_view name) {
    static std::map<std::string, Value, std::less<>> files;
    auto it = files.find(name);
    if (it == files.end()) {
        const std::string key(name);
        it = files.emplace(key, readJson(kData / "analysis" / (key + ".json"))).first;
    }
    return it->second;
}

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
inline double hexFloat(const Value& value) { return std::strtod(str(value).c_str(), nullptr); }

// Tagged JSON of the generator -> plain values: {"$float": "nan" | "inf" | "-inf"} -> double.
inline Value untag(const Value& value) {
    if (const Table* t = value.get<Table>()) {
        if (t->size() == 1 && t->front().first == "$float") {
            const std::string& s = str(t->front().second);
            if (s == "nan") {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return (s[0] == '-' ? -1.0 : 1.0) * std::numeric_limits<double>::infinity();
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

// Tables with their keys sorted (recursively), integral floats as they are: for comparing a
// result's meta / params with the generator's sorted JSON.
inline Value sortedTables(const Value& value) {
    if (const Table* t = value.get<Table>()) {
        Table out;
        for (const auto& [k, v] : *t) {
            out.emplace_back(k, sortedTables(v));
        }
        std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
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

// conftest.make_universe of a micro system: atoms (name, resid, resname, element, mass, charge),
// residues in order of first appearance, one segment, bonds only when there are some, float32
// coordinates in memory, times index * dt, no cell.
inline core::System microSystem(const Value& system) {
    auto t = std::make_shared<core::Topology>();
    std::vector<std::int64_t> resids;
    std::vector<std::uint32_t> elements;
    std::vector<double> masses;
    std::vector<double> charges;
    for (const Value& atom : list(at(system, "atoms"))) {
        const std::int64_t resid = integer(at(atom, "resid"));
        auto it = std::find(resids.begin(), resids.end(), resid);
        if (it == resids.end()) {
            resids.push_back(resid);
            t->resids.push_back(resid);
            t->resnums.push_back(resid);
            t->resnames.push_back(t->strings.intern(str(at(atom, "resname"))));
            t->residueSegment.push_back(0);
            it = resids.end() - 1;
        }
        t->names.push_back(t->strings.intern(str(at(atom, "name"))));
        t->types.push_back(0);
        t->atomResidue.push_back(static_cast<std::int32_t>(it - resids.begin()));
        elements.push_back(t->strings.intern(str(at(atom, "element"))));
        masses.push_back(hexFloat(at(atom, "mass")));
        charges.push_back(hexFloat(at(atom, "charge")));
    }
    t->elements = std::move(elements);
    t->masses = std::move(masses);
    t->charges = std::move(charges);
    t->segids.push_back(t->strings.intern(str(at(system, "segid"))));
    const Array& bonds = list(at(system, "bonds"));
    if (!bonds.empty()) {
        std::vector<core::Bond> out;
        for (const Value& b : bonds) {
            const auto x = static_cast<std::int32_t>(integer(list(b)[0]));
            const auto y = static_cast<std::int32_t>(integer(list(b)[1]));
            out.push_back({std::min(x, y), std::max(x, y)});
        }
        t->bonds = std::move(out);
    }
    t->finish();
    std::vector<core::Frame> frames;
    for (const Value& coords : list(at(system, "coords"))) {
        core::Frame f;
        for (const Value& v : list(coords)) {
            f.positions.push_back(static_cast<float>(hexFloat(v)));
        }
        frames.push_back(std::move(f));
    }
    return core::System{std::move(t),
                        std::make_unique<core::MemoryTrajectory>(std::move(frames), hexFloat(at(system, "dt")))};
}

// "{data}/..." / "{out}/..." of the generator -> a path under data/ (OUT is data/analysis).
inline fs::path dataPath(const std::string& text) {
    for (const auto& [prefix, base] : {std::pair<std::string, fs::path>{"{out}/", kData / "analysis"},
                                       std::pair<std::string, fs::path>{"{data}/", kData}}) {
        if (text.rfind(prefix, 0) == 0) {
            return base / text.substr(prefix.size());
        }
    }
    throw std::runtime_error("unexpected path " + text);
}

// The system of a case: a micro system (micro_cases.json) or files loaded like load_universe.
inline core::System caseSystem(const std::string& file, const Value& c) {
    const Value& systems = at(caseFile(file), "systems");
    const Value& system = at(systems, str(at(c, "system")));
    if (file == "micro_cases") {
        return microSystem(system);
    }
    std::vector<fs::path> trajectories;
    for (const Value& p : list(at(system, "trajectories"))) {
        trajectories.push_back(dataPath(str(p)));
    }
    return io::loadSystem(dataPath(str(at(system, "topology"))), trajectories);
}

inline core::Project caseProject(const Value& c) { return core::projectFromPlain(untag(at(c, "project"))); }

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
    void process(std::size_t k, const core::Frame& frame) override { m_inner->process(k, frame); }
    core::AnalysisResult finalize() override { return m_inner->finalize(); }

private:
    std::unique_ptr<analysis::Analysis> m_inner;
    std::int64_t m_setupFrame;
};

// numpy's dtype.str of an array
inline std::string numpyDtype(const core::NumArray& a) {
    switch (a.dtype()) {
    case core::NumArray::DType::Bool: return "|b1";
    case core::NumArray::DType::Int8: return "|i1";
    case core::NumArray::DType::Int32: return "<i4";
    case core::NumArray::DType::Int64: return "<i8";
    case core::NumArray::DType::Float32: return "<f4";
    case core::NumArray::DType::Float64: return "<f8";
    case core::NumArray::DType::String: return "<U" + std::to_string(a.stringWidth());
    }
    return "?";
}

inline double ulp32(double x) {
    const float f = static_cast<float>(std::abs(x));
    return static_cast<double>(std::nextafter(f, std::numeric_limits<float>::infinity()) - f);
}

}  // namespace mdbridge::test::structure
