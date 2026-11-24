#pragma once
// Helpers for the analysis parity tests of tests/parity/make_analysis_cases.py (data/analysis, see
// data/README.md): the systems of the cases (micro systems rebuilt from their atom tables as
// conftest.make_universe builds them, files loaded as MDATK loads them), the tagged JSON of
// projects / params / meta, numpy arrays as the generator records them, and the decision margins.
// (W2.2's copy; to be unified with analysis_test_support.hpp.)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/loader.hpp"
#include "test_support.hpp"

namespace mdbridge::test::mgsupport {

inline const core::Value& at(const core::Value& value, std::string_view key) {
    const core::Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
inline const core::Array& list(const core::Value& value) { return *value.get<core::Array>(); }
inline const std::string& str(const core::Value& value) { return *value.get<std::string>(); }
inline std::int64_t integer(const core::Value& value) { return *value.get<std::int64_t>(); }
// float.hex(), "inf", "-inf", "nan"
inline double hexFloat(const core::Value& value) { return std::strtod(str(value).c_str(), nullptr); }

inline std::vector<std::string> strings(const core::Value& value) {
    std::vector<std::string> out;
    for (const core::Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}

// The generator's strict JSON back to values: {"$float": "nan" | "inf" | "-inf"} -> double.
inline core::Value untagged(const core::Value& value) {
    if (const core::Table* table = value.get<core::Table>()) {
        if (table->size() == 1 && table->front().first == "$float") {
            return hexFloat(table->front().second);
        }
        core::Table out;
        for (const auto& [k, v] : *table) {
            out.emplace_back(k, untagged(v));
        }
        return out;
    }
    if (const core::Array* array = value.get<core::Array>()) {
        core::Array out;
        for (const core::Value& v : *array) {
            out.push_back(untagged(v));
        }
        return out;
    }
    return value;
}

// Tables with their keys sorted, recursively (the generator writes JSON with sorted keys).
inline core::Value sortedKeys(const core::Value& value) {
    if (const core::Table* table = value.get<core::Table>()) {
        core::Table out;
        for (const auto& [k, v] : *table) {
            out.emplace_back(k, sortedKeys(v));
        }
        std::sort(out.begin(), out.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
        return out;
    }
    if (const core::Array* array = value.get<core::Array>()) {
        core::Array out;
        for (const core::Value& v : *array) {
            out.push_back(sortedKeys(v));
        }
        return out;
    }
    return value;
}

// A micro system of micro_cases.json as conftest.make_universe builds it: residues in order of
// first appearance of their resid (resname of the first atom), one segment "SYS", elements,
// masses, charges, optional bonds, float32 coordinates in memory, no cell.
inline core::System microSystem(const core::Value& system) {
    auto t = std::make_shared<core::Topology>();
    std::vector<std::int64_t> resids;
    t->elements.emplace();
    t->masses.emplace();
    t->charges.emplace();
    for (const core::Value& atom : list(at(system, "atoms"))) {
        const std::int64_t resid = integer(at(atom, "resid"));
        auto it = std::find(resids.begin(), resids.end(), resid);
        if (it == resids.end()) {
            resids.push_back(resid);
            t->resids.push_back(resid);
            t->resnames.push_back(t->strings.intern(str(at(atom, "resname"))));
            t->residueSegment.push_back(0);
            it = std::prev(resids.end());
        }
        t->atomResidue.push_back(static_cast<std::int32_t>(it - resids.begin()));
        t->names.push_back(t->strings.intern(str(at(atom, "name"))));
        t->elements->push_back(t->strings.intern(str(at(atom, "element"))));
        t->masses->push_back(hexFloat(at(atom, "mass")));
        t->charges->push_back(hexFloat(at(atom, "charge")));
    }
    t->segids.push_back(t->strings.intern(str(at(system, "segid"))));
    const core::Array& bonds = list(at(system, "bonds"));
    if (!bonds.empty()) {
        t->bonds.emplace();
        for (const core::Value& bond : bonds) {
            const auto a = static_cast<std::int32_t>(integer(list(bond)[0]));
            const auto b = static_cast<std::int32_t>(integer(list(bond)[1]));
            t->bonds->push_back({std::min(a, b), std::max(a, b)});
        }
    }
    t->finish();
    std::vector<core::Frame> frames;
    for (const core::Value& coords : list(at(system, "coords"))) {
        core::Frame frame;
        for (const core::Value& v : list(coords)) {
            frame.positions.push_back(static_cast<float>(hexFloat(v)));
        }
        frames.push_back(std::move(frame));
    }
    return core::System{t, std::make_unique<core::MemoryTrajectory>(std::move(frames), hexFloat(at(system, "dt")))};
}

// "{data}/..." (tests/data) and "{out}/..." (tests/data/analysis) paths of the case files
inline fs::path dataPath(const std::string& text) {
    if (text.rfind("{data}/", 0) == 0) {
        return kData / text.substr(7);
    }
    if (text.rfind("{out}/", 0) == 0) {
        return kData / "analysis" / text.substr(6);
    }
    throw std::runtime_error("unexpected path " + text);
}

// A system of a case file's "systems" table, loaded as MDATK's Pipeline loads it.
inline core::System fileSystem(const core::Value& system, std::optional<double> dtPs) {
    std::vector<fs::path> trajectories;
    for (const core::Value& t : list(at(system, "trajectories"))) {
        trajectories.push_back(dataPath(str(t)));
    }
    return io::loadSystem(dataPath(str(at(system, "topology"))), trajectories, dtPs);
}

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
            for (const std::int8_t v : a.values<std::int8_t>()) out.push_back(v);
            break;
        case core::NumArray::DType::Int32:
            for (const std::int32_t v : a.values<std::int32_t>()) out.push_back(v);
            break;
        case core::NumArray::DType::Int64:
            for (const std::int64_t v : a.values<std::int64_t>()) out.push_back(v);
            break;
        default:
            break;
    }
    return out;
}

// The smallest recorded margin of frame k over every quantity of a result's margins
// (+inf when nothing was recorded for the frame).
inline double smallestMargin(const core::Value& margins, std::size_t k) {
    double out = std::numeric_limits<double>::infinity();
    const core::Value* frames = margins.find("frames");
    if (!frames) {
        return out;
    }
    for (const auto& [name, values] : *frames->get<core::Table>()) {
        const core::Value& v = list(values).at(k);
        if (!v.isNull()) {
            out = std::min(out, std::abs(hexFloat(v)));
        }
    }
    return out;
}

}  // namespace mdbridge::test::mgsupport
