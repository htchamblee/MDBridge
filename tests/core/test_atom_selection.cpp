// MDAnalysis' selection language (mdbridge/core/atom_selection.hpp) against MDAnalysis 2.11.0:
// thousands of selection strings on the fixtures of tests/data/topology and the distance keywords
// on frames with and without boxes (tests/data/selections,
// tests/parity/make_atom_selection_cases.py).

#include <chrono>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/atom_selection.hpp"
#include "topology_dump.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

std::vector<std::int32_t> expandRuns(const Value& runs) {
    std::vector<std::int32_t> out;
    for (const Value& run : *runs.get<Array>()) {
        const Array& pair = *run.get<Array>();
        for (std::int64_t i = *pair[0].get<std::int64_t>(); i <= *pair[1].get<std::int64_t>();
             ++i) {
            out.push_back(static_cast<std::int32_t>(i));
        }
    }
    return out;
}

std::string describe(const std::vector<std::int32_t>& atoms) {
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < atoms.size() && i < 40; ++i) {
        out << (i ? " " : "") << atoms[i];
    }
    out << (atoms.size() > 40 ? " ...]" : "]") << " (" << atoms.size() << ")";
    return out.str();
}

// Runs every case; returns a description of each disagreement.
std::vector<std::string> runCases(const Topology& topology, const Frame* frame,
                                  const Array& cases) {
    std::vector<std::string> problems;
    for (const Value& c : cases) {
        const std::string& sel = text(c, "sel");
        const bool unsupported = c.find("unsupported") != nullptr;
        std::string got;
        std::vector<std::int32_t> atoms;
        bool failed = false;
        try {
            atoms = selectAtoms(topology, frame, sel);
        } catch (const AtomSelectionError& e) {
            failed = true;
            got = e.what();
        }
        if (unsupported && failed && got.find("not supported by MDBridge") != std::string::npos) {
            continue;
        }
        if (const Value* error = c.find("error")) {
            if (!failed || got != *error->get<std::string>()) {
                problems.push_back("'" + sel + "': expected error '" + *error->get<std::string>() +
                                   "', got " + (failed ? "error '" + got + "'" : describe(atoms)));
            }
            continue;
        }
        const std::vector<std::int32_t> expected = expandRuns(*c.find("atoms"));
        if (failed || atoms != expected) {
            problems.push_back("'" + sel + "': expected " + describe(expected) + ", got " +
                               (failed ? "error '" + got + "'" : describe(atoms)));
        }
    }
    return problems;
}

void report(const std::vector<std::string>& problems, const std::string& where) {
    std::ostringstream out;
    for (std::size_t i = 0; i < problems.size() && i < 30; ++i) {
        out << "  " << problems[i] << "\n";
    }
    EXPECT_TRUE(problems.empty()) << where << ": " << problems.size() << " disagreements\n"
                                  << out.str();
}

}  // namespace

TEST(AtomSelection, MatchesMdanalysisOnTheFixtures) {
    const Value data = readJson(kData / "selections" / "atom_selection_cases.json");
    const Table& fixtures = *data.find("fixtures")->get<Table>();
    std::size_t total = 0;
    for (const auto& [fixture, cases] : fixtures) {
        const TopologyDump dump = topologyFixture(fixture);
        const Frame* frame = dump.frame ? &*dump.frame : nullptr;
        report(runCases(*dump.topology, frame, *cases.get<Array>()), fixture);
        total += cases.get<Array>()->size();
    }
    EXPECT_GT(total, 10000u);
}

TEST(AtomSelection, DistanceKeywordsMatchMdanalysis) {
    const Value data = readJson(kData / "selections" / "atom_selection_distance_cases.json");
    for (const Value& f : *data.find("frames")->get<Array>()) {
        const TopologyDump dump = topologyFixture(text(f, "fixture"));
        Frame frame;
        for (const Value& row : *f.find("positions")->get<Array>()) {
            for (const Value& x : *row.get<Array>()) {
                frame.positions.push_back(static_cast<float>(detail::number(x)));
            }
        }
        if (const Value* dims = f.find("dimensions"); dims && dims->isArray()) {
            std::array<float, 6> box{};
            for (std::size_t i = 0; i < 6; ++i) {
                box[i] = static_cast<float>(detail::number((*dims->get<Array>())[i]));
            }
            frame.box = box;
        }
        report(runCases(*dump.topology, &frame, *f.find("cases")->get<Array>()),
               text(f, "fixture") + " " + text(f, "label"));
    }
}

TEST(AtomSelection, ErrorsAreSelectionErrors) {
    const TopologyDump dump = topologyFixture("fixture_chains.pdb");
    EXPECT_TRUE(selectAtoms(*dump.topology, nullptr, "").empty());
    try {
        selectAtoms(*dump.topology, nullptr, "foo");
        FAIL() << "no error";
    } catch (const ValueError& e) {
        EXPECT_STREQ(e.what(), "Unknown selection token: 'foo'");
        EXPECT_STREQ(e.pythonName(), "SelectionError");
    }
    // distance keywords without coordinates
    EXPECT_THROW(selectAtoms(*dump.topology, nullptr, "around 3 name P"), AtomSelectionError);
}

// A million atoms: string matching runs once per distinct value, boolean logic on masks.
TEST(AtomSelection, MillionAtoms) {
    auto topology = std::make_shared<Topology>();
    Topology& t = *topology;
    const std::uint32_t names[] = {t.strings.intern("P"), t.strings.intern("C4'"),
                                   t.strings.intern("N1"), t.strings.intern("O"),
                                   t.strings.intern("H1")};
    const std::uint32_t resnames[] = {t.strings.intern("G"), t.strings.intern("ALA"),
                                      t.strings.intern("HOH")};
    constexpr std::size_t kResidues = 200'000;
    for (std::size_t r = 0; r < kResidues; ++r) {
        for (std::size_t k = 0; k < 5; ++k) {
            t.names.push_back(names[k]);
            t.atomResidue.push_back(static_cast<std::int32_t>(r));
        }
        t.resids.push_back(static_cast<std::int64_t>(r + 1));
        t.resnames.push_back(resnames[r % 3]);
        t.residueSegment.push_back(r < kResidues / 2 ? 0 : 1);
    }
    t.segids = {t.strings.intern("A"), t.strings.intern("B")};
    t.masses = std::vector<double>(t.names.size(), 12.011);
    t.finish();
    Frame frame;
    frame.positions.resize(3 * t.atomCount());
    for (std::size_t a = 0; a < t.atomCount(); ++a) {
        frame.positions[3 * a] = static_cast<float>(a % 1000) * 1.5f;
        frame.positions[3 * a + 1] = static_cast<float>((a / 1000) % 100) * 1.5f;
        frame.positions[3 * a + 2] = static_cast<float>(a / 100000) * 1.5f;
    }
    const auto start = std::chrono::steady_clock::now();
    EXPECT_EQ(selectAtoms(t, &frame, "name P").size(), kResidues);
    EXPECT_EQ(selectAtoms(t, &frame, "name C4* and resname G").size(), (kResidues + 2) / 3);
    // (nucleic or protein) and not name H*: and / or are left-associative
    EXPECT_EQ(selectAtoms(t, &frame, "nucleic or protein and not name H*").size(), 133'334u * 4);
    EXPECT_EQ(selectAtoms(t, &frame, "resid 1000:1999 and segid A").size(), 5000u);
    EXPECT_EQ(selectAtoms(t, &frame, "byres index 7").size(), 5u);
    EXPECT_EQ(selectAtoms(t, &frame, "mass 12.011 and bynum 1-10").size(), 10u);
    EXPECT_EQ(selectAtoms(t, &frame, "around 1.6 index 0").size(), 3u);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
#ifdef NDEBUG
    EXPECT_LT(seconds, 1.0);
#else
    EXPECT_LT(seconds, 20.0);
#endif
}
