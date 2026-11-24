#pragma once
// Helpers of the energy tests for the analysis cases of tests/parity/make_analysis_cases.py
// (data/analysis/{micro,demo_dry,demo_dry_variants,systems}_cases.json, see data/README.md): JSON
// access, the micro systems built from their atom tables (MDATK's tests/conftest.make_universe)
// and the file systems of the cases. To be unified with analysis_test_support.hpp.

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/topology.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/loader.hpp"
#include "test_support.hpp"

namespace mdbridge::test::energy_support {

inline const core::Value& at(const core::Value& value, std::string_view key) {
    const core::Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
inline const core::Array& list(const core::Value& value) {
    const core::Array* a = value.get<core::Array>();
    if (!a) {
        throw std::runtime_error("not an array");
    }
    return *a;
}
inline const std::string& str(const core::Value& value) {
    const std::string* s = value.get<std::string>();
    if (!s) {
        throw std::runtime_error("not a string");
    }
    return *s;
}
inline std::int64_t integer(const core::Value& value) { return *value.get<std::int64_t>(); }
// float.hex() text
inline double hexFloat(const core::Value& value) { return std::strtod(str(value).c_str(), nullptr); }

inline std::vector<std::string> strings(const core::Value& value) {
    std::vector<std::string> out;
    for (const core::Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}

inline std::vector<std::int64_t> integers(const core::Value& value) {
    std::vector<std::int64_t> out;
    for (const core::Value& v : list(value)) {
        out.push_back(integer(v));
    }
    return out;
}

inline std::vector<double> hexFloats(const core::Value& value) {
    std::vector<double> out;
    for (const core::Value& v : list(value)) {
        out.push_back(hexFloat(v));
    }
    return out;
}

// A JSON value equal to another with table keys in any order (the generator sorts them).
inline bool sameValue(const core::Value& a, const core::Value& b) {
    if (const core::Table* ta = a.get<core::Table>()) {
        const core::Table* tb = b.get<core::Table>();
        if (!tb || ta->size() != tb->size()) {
            return false;
        }
        for (const auto& [key, value] : *ta) {
            const core::Value* other = b.find(key);
            if (!other || !sameValue(value, *other)) {
                return false;
            }
        }
        return true;
    }
    if (const core::Array* aa = a.get<core::Array>()) {
        const core::Array* ab = b.get<core::Array>();
        if (!ab || aa->size() != ab->size()) {
            return false;
        }
        for (std::size_t i = 0; i < aa->size(); ++i) {
            if (!sameValue((*aa)[i], (*ab)[i])) {
                return false;
            }
        }
        return true;
    }
    return a == b;
}

// The micro system `spec` of micro_cases.json as conftest.make_universe builds it: residues in
// order of first appearance of their resid, segment SYS, names / resnames / elements / masses /
// charges from the atom table, bonds only when listed, float32 coordinates in memory.
inline core::System microSystem(const core::Value& spec) {
    auto t = std::make_shared<core::Topology>();
    const core::Array& atoms = list(at(spec, "atoms"));
    std::vector<std::int64_t> resids;
    std::vector<std::uint32_t> elements;
    std::vector<double> masses, charges;
    for (const core::Value& atom : atoms) {
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
        masses.push_back(hexFloat(at(atom, "mass")));
        charges.push_back(hexFloat(at(atom, "charge")));
    }
    t->elements = std::move(elements);
    t->masses = std::move(masses);
    t->charges = std::move(charges);
    t->segids.push_back(t->strings.intern(str(at(spec, "segid"))));
    const core::Array& bonds = list(at(spec, "bonds"));
    if (!bonds.empty()) {
        std::vector<core::Bond> out;
        for (const core::Value& bond : bonds) {
            const auto pair = integers(bond);
            out.push_back({static_cast<std::int32_t>(pair.at(0)), static_cast<std::int32_t>(pair.at(1))});
        }
        t->bonds = std::move(out);
    }
    t->finish();
    std::vector<core::Frame> frames;
    for (const core::Value& coords : list(at(spec, "coords"))) {
        core::Frame frame;
        for (const double x : hexFloats(coords)) {
            frame.positions.push_back(static_cast<float>(x));
        }
        frames.push_back(std::move(frame));
    }
    return core::System{std::move(t),
                        std::make_unique<core::MemoryTrajectory>(std::move(frames), hexFloat(at(spec, "dt")))};
}

// "{data}/..." (tests/data) and "{out}/..." (tests/data/analysis) paths of the case files.
inline std::filesystem::path casePath(const std::string& text) {
    for (const auto& [prefix, root] : {std::pair<std::string_view, std::filesystem::path>{"{data}/", kData},
                                       std::pair<std::string_view, std::filesystem::path>{"{out}/", kData / "analysis"}}) {
        if (text.starts_with(prefix)) {
            return root / text.substr(prefix.size());
        }
    }
    throw std::runtime_error("unknown path " + text);
}

struct LoadedSystem {
    core::System system;
    std::string topologyPath;  // "" for micro systems (MDATK's run_one passes none)
};

// The system named `name` of a case file: a micro system (micro_cases.json) or {topology,
// trajectories} loaded from files.
inline LoadedSystem caseSystem(const core::Value& systems, const std::string& name) {
    const core::Value& spec = at(systems, name);
    if (spec.find("atoms")) {
        return {microSystem(spec), {}};
    }
    const std::filesystem::path topology = casePath(str(at(spec, "topology")));
    std::vector<std::filesystem::path> trajectories;
    for (const std::string& t : strings(at(spec, "trajectories"))) {
        trajectories.push_back(casePath(t));
    }
    return {io::loadSystem(topology, trajectories), topology.string()};
}

}  // namespace mdbridge::test::energy_support
