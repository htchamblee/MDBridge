// The system model: topologies as MDAnalysis builds them (dumps of tests/parity/dump_topologies.py),
// atom groups with MDAnalysis' ordering rules, chains as MDATK numbers and names them
// (analysis/structure.py chain_ids / chain_names; data/topology/golden_model.json) and frames.

#include <numeric>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/topology.hpp"
#include "topology_dump.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

std::vector<std::int32_t> allResidues(const Topology& topology) {
    std::vector<std::int32_t> out(topology.residueCount());
    std::iota(out.begin(), out.end(), 0);
    return out;
}

template <class T>
std::vector<T> numbers(const Value& array) {
    std::vector<T> out;
    for (const Value& v : *array.get<Array>()) {
        out.push_back(static_cast<T>(*v.get<std::int64_t>()));
    }
    return out;
}

std::vector<std::string> texts(const Value& array) {
    std::vector<std::string> out;
    for (const Value& v : *array.get<Array>()) {
        out.push_back(*v.get<std::string>());
    }
    return out;
}

// Universe.empty(n, n_residues=n, atom_resindex=range(n)) with one segment
std::shared_ptr<Topology> emptyUniverse(std::size_t n, const std::string& segid) {
    auto topology = std::make_shared<Topology>();
    for (std::size_t i = 0; i < n; ++i) {
        topology->names.push_back(0);
        topology->atomResidue.push_back(static_cast<std::int32_t>(i));
        topology->resids.push_back(static_cast<std::int64_t>(i + 1));
        topology->resnames.push_back(0);
        topology->residueSegment.push_back(0);
    }
    topology->segids.push_back(topology->strings.intern(segid));
    return topology;
}

}  // namespace

TEST(Topology, FixturesAsMdanalysisBuildsThem) {
    const TopologyDump chains = topologyFixture("fixture_chains.pdb");
    const Topology& t = *chains.topology;
    EXPECT_EQ(t.atomCount(), 60u);
    EXPECT_EQ(t.residueCount(), 19u);
    EXPECT_EQ(t.segmentCount(), 5u);
    EXPECT_TRUE(t.icodes && t.chainIDs && t.elements);
    EXPECT_FALSE(t.bonds);
    ASSERT_TRUE(chains.frame);
    EXPECT_EQ(chains.frame->positions.size(), 180u);
    EXPECT_EQ(chains.frame->positions[0], 1.2f);
    // residue 3 and 3A are two residues with their own atoms
    EXPECT_EQ(t.resids[2], 3);
    EXPECT_EQ(t.resids[3], 3);
    EXPECT_EQ(t.strings[(*t.icodes)[3]], "A");
    ASSERT_EQ(t.residueAtoms(3).size(), 2u);
    EXPECT_EQ(t.name(t.residueAtoms(3)[0]), "P");

    const TopologyDump parm = topologyFixture("demo_dry.parm7");
    EXPECT_FALSE(parm.topology->chainIDs);
    EXPECT_FALSE(parm.topology->icodes);
    ASSERT_TRUE(parm.topology->bonds);
    EXPECT_EQ(parm.topology->bonds->size(), 828u);
    EXPECT_FALSE(parm.frame);
}

TEST(Topology, ChainIdsAndNamesMatchMdatk) {
    const Value golden = readJson(kData / "topology" / "golden_model.json");
    for (const Value& c : *golden.find("chain_ids")->get<Array>()) {
        const std::string& file = text(c, "file");
        std::shared_ptr<const Topology> topology;
        if (file == "chains B then A") {
            topology = topologyFixture("fixture_chains_ba.pdb").topology;
        } else if (file.rfind("Universe.empty", 0) == 0) {
            auto t = emptyUniverse(4, "SYSTEM");
            const std::uint32_t a = t->strings.intern("RNA_chain_A");
            const std::uint32_t b = t->strings.intern("RNA_chain_B");
            t->chainIDs = std::vector<std::uint32_t>{a, a, b, b};
            t->finish();
            topology = t;
        } else if (file == "empty no bonds") {
            auto t = emptyUniverse(3, "SYS");
            t->finish();
            topology = t;
        } else {
            topology = topologyFixture(file).topology;
        }
        const std::vector<std::int32_t> residues = allResidues(*topology);
        EXPECT_EQ(chainIds(*topology, residues), numbers<std::int64_t>(*c.find("chain_ids"))) << file;
        EXPECT_EQ(chainNames(*topology, residues), texts(*c.find("chain_names"))) << file;
    }
}

TEST(Topology, AtomGroupsKeepMdanalysisOrder) {
    const auto topology = topologyFixture("fixture_chains.pdb").topology;
    const AtomGroup all(topology);
    EXPECT_EQ(all.size(), 60u);
    const std::vector<std::int64_t> picks{5, 2, 2, -1};
    const AtomGroup picked = all.at(picks);
    EXPECT_EQ(picked.indices(), (std::vector<std::int32_t>{5, 2, 2, 59}));
    EXPECT_EQ(picked.unique().indices(), (std::vector<std::int32_t>{2, 5, 59}));
    EXPECT_EQ(picked.residues(), (std::vector<std::int32_t>{0, 18}));
    const AtomGroup other(topology, {59, 1});
    EXPECT_EQ((picked | other).indices(), (std::vector<std::int32_t>{1, 2, 5, 59}));
    std::vector<bool> mask(picked.size(), false);
    mask[0] = mask[3] = true;
    EXPECT_EQ(picked.masked(mask).indices(), (std::vector<std::int32_t>{5, 59}));
    const std::vector<std::int32_t> residues{3, 0};
    EXPECT_EQ(residueAtoms(topology, residues).indices(),
              (std::vector<std::int32_t>{17, 18, 0, 1, 2, 3, 4, 5, 6, 7}));
    EXPECT_THROW((void)all.bonds(), NoDataError);
    EXPECT_THROW((void)fragmentIndices(*topology), NoDataError);
    const std::vector<std::int64_t> outside{60};
    EXPECT_THROW((void)all.at(outside), std::out_of_range);
}

TEST(Topology, FragmentsAreNumberedByLowestAtom) {
    auto t = emptyUniverse(6, "S");
    t->bonds = std::vector<Bond>{{4, 2}, {5, 0}, {3, 1}};
    t->finish();
    EXPECT_EQ(fragmentIndices(*t), (std::vector<std::int32_t>{0, 1, 2, 1, 2, 0}));
    EXPECT_EQ(*t->bonds, (std::vector<Bond>{{0, 5}, {1, 3}, {2, 4}}));
    const AtomGroup group(t, {3});
    EXPECT_EQ(group.bonds(), (std::vector<Bond>{{1, 3}}));
}

TEST(Topology, FinishRejectsInconsistentData) {
    auto t = emptyUniverse(2, "S");
    t->atomResidue[1] = 5;
    EXPECT_THROW(t->finish(), std::invalid_argument);
    auto u = emptyUniverse(2, "S");
    u->masses = std::vector<double>{1.0};
    EXPECT_THROW(u->finish(), std::invalid_argument);
}

TEST(Topology, MemoryTrajectory) {
    std::vector<Frame> frames(3);
    for (Frame& f : frames) {
        f.positions = {0.0f, 0.0f, 0.0f};
    }
    MemoryTrajectory trajectory(frames, 2.5);
    EXPECT_EQ(trajectory.frameCount(), 3);
    EXPECT_EQ(trajectory.current().index, 0);
    EXPECT_EQ(trajectory.seek(-1).time, 5.0);
    EXPECT_EQ(trajectory.current().index, 2);
    EXPECT_THROW(trajectory.seek(3), std::out_of_range);
}
