// Loading systems as MDATK does (mdatk/io/loader.py): the PDB / GRO topologies and coordinates
// MDAnalysis 2.11.0 reads, load_universe over a matrix of files, frame_indices, system_summary and
// `mdatk info`, against tests/parity/make_loader_cases.py (data/io/loader, see data/README.md) and
// the topology fixtures of data/topology.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/io/loader.hpp"
#include "topology_dump.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const fs::path kLoader = kData / "io" / "loader";

// Coordinate files as topologies: MDAnalysis builds a minimal topology (atoms only) that MDATK's
// summary then fails on; MDBridge rejects the file as a topology.
const std::set<std::string> kMinimalTopologies = {"{data}/io/trajectories/files/small.xtc"};

const core::Array& list(const core::Value& value) {
    return *value.get<core::Array>();
}

std::string pathText(const fs::path& path) {
    const std::u8string text = path.u8string();
    return std::string(text.begin(), text.end());
}

std::string replaceAll(std::string text, std::string_view from, const std::string& to) {
    for (std::size_t at = text.find(from); at != std::string::npos;
         at = text.find(from, at + to.size())) {
        text.replace(at, from.size(), to);
    }
    return text;
}

std::string withPaths(const std::string& text) {
    return replaceAll(replaceAll(text, "{out}", pathText(kLoader)), "{data}", pathText(kData));
}

std::string hexBits(const float* values, std::size_t n) {
    static const char digits[] = "0123456789abcdef";
    std::string out;
    for (std::size_t i = 0; i < n; ++i) {
        std::uint32_t bits;
        std::memcpy(&bits, &values[i], 4);
        for (int shift = 28; shift >= 0; shift -= 4) {
            out += digits[bits >> shift & 15];
        }
    }
    return out;
}

double hexDouble(const std::string& text) {
    return std::strtod(text.c_str(), nullptr);
}

// Runs fn and compares its error with {"type", "message"}.
template <class Fn>
void expectError(Fn&& fn, const core::Value& error, const std::string& context) {
    const std::string message = withPaths(text(error, "message"));
    try {
        fn();
        ADD_FAILURE() << context << ": no error, expected " << message;
    } catch (const core::ValueError& e) {
        EXPECT_EQ(std::string(e.pythonName()), text(error, "type")) << context;
        EXPECT_EQ(e.what(), message) << context;
    } catch (const fs::filesystem_error& e) {
        EXPECT_EQ(text(error, "type"), "FileNotFoundError") << context;
        EXPECT_EQ(e.what(), message) << context;
    }
}

bool same(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || a == b;
}

template <class T>
bool same(const T& a, const T& b) {
    return a == b;
}

template <class T>
void expectSame(const std::optional<std::vector<T>>& actual,
                const std::optional<std::vector<T>>& expected, const std::string& what,
                const std::string& context) {
    ASSERT_EQ(actual.has_value(), expected.has_value()) << context << " " << what;
    if (actual) {
        ASSERT_EQ(actual->size(), expected->size()) << context << " " << what;
        for (std::size_t i = 0; i < actual->size(); ++i) {
            EXPECT_TRUE(same((*actual)[i], (*expected)[i])) << context << " " << what << " " << i;
        }
    }
}

std::vector<std::string> texts(const core::Topology& top, const std::vector<std::uint32_t>& ids) {
    std::vector<std::string> out;
    for (const std::uint32_t id : ids) {
        out.push_back(top.strings[id]);
    }
    return out;
}

std::optional<std::vector<std::string>> texts(
    const core::Topology& top, const std::optional<std::vector<std::uint32_t>>& ids) {
    if (!ids) {
        return std::nullopt;
    }
    return texts(top, *ids);
}

void expectTopology(const core::Topology& a, const core::Topology& e, const std::string& context) {
    EXPECT_EQ(texts(a, a.names), texts(e, e.names)) << context << " names";
    EXPECT_EQ(texts(a, a.types), texts(e, e.types)) << context << " types";
    EXPECT_EQ(a.atomResidue, e.atomResidue) << context << " resindex";
    expectSame(texts(a, a.elements), texts(e, e.elements), "elements", context);
    expectSame(texts(a, a.chainIDs), texts(e, e.chainIDs), "chainIDs", context);
    expectSame(a.masses, e.masses, "masses", context);
    expectSame(a.charges, e.charges, "charges", context);
    expectSame(a.ids, e.ids, "ids", context);
    expectSame(texts(a, a.altLocs), texts(e, e.altLocs), "altLocs", context);
    expectSame(texts(a, a.recordTypes), texts(e, e.recordTypes), "record_types", context);
    expectSame(a.occupancies, e.occupancies, "occupancies", context);
    expectSame(a.tempFactors, e.tempFactors, "tempfactors", context);
    expectSame(a.formalCharges, e.formalCharges, "formalcharges", context);
    EXPECT_EQ(a.resids, e.resids) << context << " resids";
    EXPECT_EQ(a.resnums, e.resnums) << context << " resnums";
    EXPECT_EQ(texts(a, a.resnames), texts(e, e.resnames)) << context << " resnames";
    expectSame(texts(a, a.icodes), texts(e, e.icodes), "icodes", context);
    EXPECT_EQ(a.residueSegment, e.residueSegment) << context << " segindex";
    EXPECT_EQ(texts(a, a.segids), texts(e, e.segids)) << context << " segids";
    expectSame(a.bonds, e.bonds, "bonds", context);
}

void expectFrame(const core::Frame& frame, const core::Frame& expected,
                 const std::string& context) {
    EXPECT_EQ(hexBits(frame.positions.data(), frame.positions.size()),
              hexBits(expected.positions.data(), expected.positions.size()))
        << context;
    ASSERT_EQ(frame.box.has_value(), expected.box.has_value()) << context;
    if (frame.box) {
        EXPECT_EQ(hexBits(frame.box->data(), 6), hexBits(expected.box->data(), 6)) << context;
    }
    EXPECT_EQ(frame.time, expected.time) << context;
}

void expectReader(core::Trajectory& trajectory, const core::Value& expected,
                  const std::string& context) {
    ASSERT_EQ(trajectory.frameCount(), *expected.find("n_frames")->get<std::int64_t>()) << context;
    EXPECT_EQ(trajectory.dt(), hexDouble(text(expected, "dt"))) << context;
    const core::TrajectoryInfo info = trajectory.info();
    const core::Value& header = *expected.find("header");
    if (header.isNull()) {
        EXPECT_FALSE(info.header.has_value()) << context;
        EXPECT_FALSE(info.compounds.has_value()) << context;
        EXPECT_FALSE(info.remarks.has_value()) << context;
    } else {
        EXPECT_EQ(info.header, *header.get<std::string>()) << context;
        std::vector<std::string> compounds;
        for (const core::Value& v : list(*expected.find("compound"))) {
            compounds.push_back(*v.get<std::string>());
        }
        std::vector<std::string> remarks;
        for (const core::Value& v : list(*expected.find("remarks"))) {
            remarks.push_back(*v.get<std::string>());
        }
        EXPECT_EQ(info.compounds, compounds) << context;
        EXPECT_EQ(info.remarks, remarks) << context;
    }
    for (const core::Value& f : list(*expected.find("frames"))) {
        const std::int64_t k = *f.find("frame")->get<std::int64_t>();
        const std::string where = context + " frame " + std::to_string(k);
        if (const core::Value* error = f.find("error")) {
            expectError([&] { trajectory.seek(k); }, *error, where);
            break;
        }
        const core::Frame& frame = trajectory.seek(k);
        EXPECT_EQ(hexBits(frame.positions.data(), frame.positions.size()), text(f, "positions"))
            << where;
        const core::Value& box = *f.find("dimensions");
        if (box.isNull()) {
            EXPECT_FALSE(frame.box.has_value()) << where;
        } else {
            ASSERT_TRUE(frame.box.has_value()) << where;
            EXPECT_EQ(hexBits(frame.box->data(), 6), *box.get<std::string>()) << where;
        }
        EXPECT_EQ(frame.time, hexDouble(text(f, "time"))) << where;
    }
}

// MDATK's summary with top_resnames in MDBridge's order: numpy's argsort leaves the order of equal
// counts to the CPU's sort kernel, MDBridge orders them by name.
core::Value withStableNames(const core::Value& summary, const core::Value& resnameCounts) {
    std::vector<std::pair<std::string, std::int64_t>> order;
    for (const core::Value& entry : list(resnameCounts)) {
        order.emplace_back(*list(entry)[0].get<std::string>(), *list(entry)[1].get<std::int64_t>());
    }
    std::stable_sort(order.begin(), order.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });
    core::Table names;
    for (std::size_t k = 0; k < order.size() && k < 12; ++k) {
        names.emplace_back(order[k].first, core::Value(order[k].second));
    }
    core::Table out = *summary.get<core::Table>();
    for (auto& [key, value] : out) {
        if (key == "top_resnames") {
            value = core::Value(std::move(names));
            break;
        }
    }
    return core::Value(std::move(out));
}

std::vector<std::int64_t> counts(const core::Value& summary) {
    std::vector<std::int64_t> out;
    for (const auto& [name, count] : *summary.find("top_resnames")->get<core::Table>()) {
        out.push_back(*count.get<std::int64_t>());
    }
    return out;
}

}  // namespace

TEST(Loader, TopologyFixtures) {
    for (const char* name : {"fixture_chains.pdb", "fixture_noelements.pdb", "fixture_segids.pdb",
                             "fixture_chains_ba.pdb"}) {
        const io::TopologyFile file = io::readTopologyFile(kData / "topology" / name);
        const TopologyDump dump = topologyFixture(name);
        expectTopology(*file.topology, *dump.topology, name);
        ASSERT_TRUE(file.coordinates);
        expectFrame(file.coordinates->current(), *dump.frame, name);
    }
    const io::TopologyFile parm7 = io::readTopologyFile(kData / "topology" / "demo_dry.parm7");
    EXPECT_FALSE(parm7.coordinates);
    expectTopology(*parm7.topology, *topologyFixture("demo_dry.parm7").topology, "demo_dry.parm7");
}

TEST(Loader, PdbAndGroFiles) {
    const core::Value cases = readJson(kLoader / "topology_cases.json");
    int checked = 0;
    for (const core::Value& c : list(cases)) {
        const std::string name = text(c, "file");
        const fs::path path = kLoader / "files" / name;
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::readTopologyFile(path); }, *error, name);
            ++checked;
            continue;
        }
        const TopologyDump dump = readTopologyDump(kLoader / "files" / text(c, "dump"));
        io::TopologyFile file;
        try {
            file = io::readTopologyFile(path);
        } catch (const std::exception& e) {
            ADD_FAILURE() << name << ": " << e.what();
            continue;
        }
        expectTopology(*file.topology, *dump.topology, name);
        expectFrame(file.coordinates->current(), *dump.frame, name);
        for (const core::Value& reader : list(*c.find("readers"))) {
            const core::Value& dt = *reader.find("dt_given");
            const std::optional<double> given =
                dt.isNull() ? std::nullopt : std::optional<double>(*dt.get<double>());
            io::TopologyFile again = io::readTopologyFile(path, "", given);
            expectReader(*again.coordinates, reader, name + (given ? " dt" : ""));
        }
        ++checked;
    }
    EXPECT_GT(checked, 50);
}

TEST(Loader, LoadSystemAndSummary) {
    const core::Value data = readJson(kLoader / "loader_cases.json");
    for (const core::Value& c : list(*data.find("summaries"))) {
        const fs::path topology = withPaths(text(c, "topology"));
        std::vector<fs::path> trajectories;
        for (const core::Value& t : list(*c.find("trajectories"))) {
            trajectories.push_back(withPaths(*t.get<std::string>()));
        }
        std::optional<double> dt;
        if (const core::Value* v = c.find("dt"); v && !v->isNull()) {
            dt =
                v->get<double>() ? *v->get<double>() : static_cast<double>(*v->get<std::int64_t>());
        }
        const std::string context = text(c, "topology") + " " +
                                    std::to_string(trajectories.size()) + " " +
                                    (trajectories.empty() ? "" : pathText(trajectories.front()));
        if (const core::Value* error = c.find("error")) {
            if (kMinimalTopologies.count(text(c, "topology")) == 0) {
                expectError([&] { io::loadSystem(topology, trajectories, dt); }, *error, context);
            } else {
                EXPECT_THROW(io::loadSystem(topology, trajectories, dt), core::ValueError)
                    << context;
            }
            continue;
        }
        core::System system;
        try {
            system = io::loadSystem(topology, trajectories, dt);
        } catch (const std::exception& e) {
            ADD_FAILURE() << context << ": " << e.what();
            continue;
        }
        const core::Value summary = io::systemSummary(system);
        const core::Value expected = withStableNames(*c.find("summary"), *c.find("resname_counts"));
        EXPECT_EQ(core::toJson(summary, 2), core::toJson(expected, 2)) << context;
        EXPECT_EQ(counts(summary), counts(*c.find("summary"))) << context;
        if (core::toJson(expected, 2) == core::toJson(*c.find("summary"), 2)) {
            EXPECT_EQ(core::toJson(summary, 2), text(c, "summary_json")) << context;
        }
    }
}

TEST(Loader, NaNTimeStep) {
    // dt_ps = nan is truthy in Python: the readers get it
    const core::System system =
        io::loadSystem(kData / "topology" / "demo_dry.parm7",
                       {kData / "io" / "trajectories" / "demo_dry.nc"}, std::nan(""));
    EXPECT_EQ(system.trajectory->frameCount(), 24);
}

TEST(Loader, FrameIndices) {
    const core::Value data = readJson(kLoader / "loader_cases.json");
    for (const core::Value& c : list(*data.find("frame_indices"))) {
        const std::int64_t n = *c.find("n")->get<std::int64_t>();
        const std::int64_t start = *c.find("start")->get<std::int64_t>();
        const core::Value& stopValue = *c.find("stop");
        const std::optional<std::int64_t> stop =
            stopValue.isNull() ? std::nullopt
                               : std::optional<std::int64_t>(*stopValue.get<std::int64_t>());
        const std::int64_t step = *c.find("step")->get<std::int64_t>();
        const std::string context = std::to_string(n) + " " + std::to_string(start) + " " +
                                    (stop ? std::to_string(*stop) : "None") + " " +
                                    std::to_string(step);
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::frameIndices(n, start, stop, step); }, *error, context);
            continue;
        }
        std::vector<std::int64_t> expected;
        for (const core::Value& v : list(*c.find("frames"))) {
            expected.push_back(*v.get<std::int64_t>());
        }
        EXPECT_EQ(io::frameIndices(n, start, stop, step), expected) << context;
    }
}

TEST(Loader, InfoCommand) {
    const core::Value data = readJson(kLoader / "loader_cases.json");
    for (const core::Value& c : list(*data.find("info"))) {
        std::vector<fs::path> trajectories;
        for (const core::Value& t : list(*c.find("trajectories"))) {
            trajectories.push_back(withPaths(*t.get<std::string>()));
        }
        const core::System system = io::loadSystem(withPaths(text(c, "topology")), trajectories);
        const core::Value printed = core::parseJson(text(c, "stdout"), "stdout");
        core::Table head;
        for (const auto& [key, value] : *printed.get<core::Table>()) {
            if (key != "groups" && key != "problems") {
                head.emplace_back(key, value);
            }
        }
        const core::Value summary = io::systemSummary(system);
        const core::Value expected = withStableNames(core::Value(head), *c.find("resname_counts"));
        EXPECT_EQ(core::toJson(summary, 2), core::toJson(expected, 2));
        EXPECT_EQ(counts(summary), counts(core::Value(head)));
        std::string json = core::toJson(expected, 2);
        json.resize(json.size() - 2);
        if (core::toJson(expected, 2) == core::toJson(core::Value(head), 2)) {
            EXPECT_TRUE(text(c, "stdout").starts_with(json + ",\n  \"groups\": {"))
                << text(c, "stdout");
        }
    }
}

TEST(Loader, MissingFiles) {
    try {
        io::loadSystem(kData / "no-such.parm7", {});
        ADD_FAILURE() << "no error";
    } catch (const fs::filesystem_error& e) {
        EXPECT_EQ(std::string(e.what()),
                  "topology not found: " + pathText(kData / "no-such.parm7"));
        EXPECT_EQ(e.path1(), kData / "no-such.parm7");
    }
    EXPECT_THROW(io::loadSystem("", {}), core::ValueError);
}
