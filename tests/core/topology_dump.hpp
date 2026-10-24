#pragma once
// Topologies as MDAnalysis builds them, from the dumps of tests/parity/dump_topologies.py
// (tests/data/topology/*.topology.json): the reference the selection tests run on and the file
// readers are compared with.

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "mdbridge/core/topology.hpp"
#include "test_support.hpp"

namespace mdbridge::test {

namespace detail {

inline const core::Array* array(const core::Value& table, std::string_view key) {
    const core::Value* v = table.find(key);
    return v ? v->get<core::Array>() : nullptr;
}

inline std::vector<std::uint32_t> strings(core::StringPool& pool, const core::Array& values) {
    std::vector<std::uint32_t> out;
    out.reserve(values.size());
    for (const core::Value& v : values) {
        out.push_back(pool.intern(*v.get<std::string>()));
    }
    return out;
}

inline std::optional<std::vector<std::uint32_t>> optionalStrings(core::StringPool& pool, const core::Value& table,
                                                                 std::string_view key) {
    const core::Array* values = array(table, key);
    if (!values) {
        return std::nullopt;
    }
    return strings(pool, *values);
}

inline double number(const core::Value& v) {
    if (const double* d = v.get<double>()) {
        return *d;
    }
    return static_cast<double>(*v.get<std::int64_t>());
}

inline std::optional<std::vector<double>> optionalNumbers(const core::Value& table, std::string_view key) {
    const core::Array* values = array(table, key);
    if (!values) {
        return std::nullopt;
    }
    std::vector<double> out;
    for (const core::Value& v : *values) {
        out.push_back(number(v));
    }
    return out;
}

template <class T>
std::vector<T> integers(const core::Array& values) {
    std::vector<T> out;
    out.reserve(values.size());
    for (const core::Value& v : values) {
        out.push_back(static_cast<T>(*v.get<std::int64_t>()));
    }
    return out;
}

template <class T>
std::optional<std::vector<T>> optionalIntegers(const core::Value& table, std::string_view key) {
    const core::Array* values = array(table, key);
    if (!values) {
        return std::nullopt;
    }
    return integers<T>(*values);
}

}  // namespace detail

struct TopologyDump {
    std::shared_ptr<const core::Topology> topology;
    std::optional<core::Frame> frame;  // the first frame, when the file has coordinates
};

inline TopologyDump readTopologyDump(const fs::path& path) {
    using namespace detail;
    const core::Value dump = readJson(path);
    auto topology = std::make_shared<core::Topology>();
    core::StringPool& pool = topology->strings;

    const core::Value& atoms = *dump.find("atoms");
    topology->names = strings(pool, *array(atoms, "names"));
    if (auto types = optionalStrings(pool, atoms, "types")) {
        topology->types = std::move(*types);
    }
    topology->atomResidue = integers<std::int32_t>(*array(atoms, "resindex"));
    topology->elements = optionalStrings(pool, atoms, "elements");
    topology->chainIDs = optionalStrings(pool, atoms, "chainIDs");
    topology->masses = optionalNumbers(atoms, "masses");
    topology->charges = optionalNumbers(atoms, "charges");
    topology->ids = optionalIntegers<std::int64_t>(atoms, "ids");
    topology->altLocs = optionalStrings(pool, atoms, "altLocs");
    topology->recordTypes = optionalStrings(pool, atoms, "record_types");
    topology->occupancies = optionalNumbers(atoms, "occupancies");
    topology->tempFactors = optionalNumbers(atoms, "tempfactors");
    topology->formalCharges = optionalIntegers<std::int64_t>(atoms, "formalcharges");

    const core::Value& residues = *dump.find("residues");
    topology->resids = integers<std::int64_t>(*array(residues, "resids"));
    if (auto resnums = optionalIntegers<std::int64_t>(residues, "resnums")) {
        topology->resnums = std::move(*resnums);
    }
    topology->resnames = strings(pool, *array(residues, "resnames"));
    topology->icodes = optionalStrings(pool, residues, "icodes");
    topology->residueSegment = integers<std::int32_t>(*array(residues, "segindex"));

    topology->segids = strings(pool, *array(*dump.find("segments"), "segids"));

    if (const core::Array* bonds = array(dump, "bonds")) {
        std::vector<core::Bond> out;
        for (const core::Value& bond : *bonds) {
            const core::Array& pair = *bond.get<core::Array>();
            out.push_back({static_cast<std::int32_t>(*pair[0].get<std::int64_t>()),
                           static_cast<std::int32_t>(*pair[1].get<std::int64_t>())});
        }
        topology->bonds = std::move(out);
    }
    topology->finish();

    TopologyDump out;
    if (const core::Value* frame = dump.find("frame"); frame && frame->isTable()) {
        core::Frame f;
        for (const core::Value& row : *array(*frame, "positions")) {
            for (const core::Value& x : *row.get<core::Array>()) {
                f.positions.push_back(static_cast<float>(number(x)));
            }
        }
        if (const core::Array* box = array(*frame, "dimensions")) {
            std::array<float, 6> cell{};
            for (std::size_t i = 0; i < 6 && i < box->size(); ++i) {
                cell[i] = static_cast<float>(number((*box)[i]));
            }
            f.box = cell;
        }
        f.time = number(*frame->find("time"));
        out.frame = std::move(f);
    }
    out.topology = std::move(topology);
    return out;
}

inline TopologyDump topologyFixture(const std::string& file) {
    return readTopologyDump(kData / "topology" / (file + ".topology.json"));
}

}  // namespace mdbridge::test
