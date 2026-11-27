// MDAnalysis' selection language (mdbridge/core/atom_selection.hpp) against MDAnalysis 2.11.0:
// thousands of selection strings on the fixtures of tests/data/topology and the distance keywords
// on frames with and without boxes (tests/data/selections,
// tests/parity/make_atom_selection_cases.py).

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/atom_selection.hpp"
#include "time_scale.hpp"
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

// The messages of docs/DIFFERENCES.md (SEL-01 ... SEL-10), on a fixture with coordinates.
TEST(AtomSelection, ImprovedMessages) {
    const TopologyDump dump = topologyFixture("fixture_chains.pdb");
    const Frame* frame = dump.frame ? &*dump.frame : nullptr;
    ASSERT_NE(frame, nullptr);
    auto error = [&](const std::string& sel) -> std::string {
        try {
            selectAtoms(*dump.topology, frame, sel);
        } catch (const AtomSelectionError& e) {
            return e.what();
        }
        return "(no error)";
    };
    EXPECT_EQ(error("name P and"), "Unexpected end of selection after 'and'");
    EXPECT_EQ(error("(name P"), "Unexpected end of selection (expected ')')");
    EXPECT_EQ(error(" \t "), "empty selection");
    EXPECT_EQ(error("around 3"), "around: expected a distance and a selection, e.g. 'around 3.5 resname MG'");
    EXPECT_EQ(error("point 1 2 3"), "point: expected x y z radius");
    EXPECT_EQ(error("atom A 1"), "atom: expected segid resid name");
    EXPECT_EQ(error("atom A 1 and"), "atom: expected segid resid name");
    EXPECT_EQ(error("same residue"), "same: expected 'same <property> as <selection>'");
    EXPECT_EQ(error("prop x <"), "prop: expected 'prop [abs] <property> <op> <value>'");
    EXPECT_EQ(error("prop x ><5"), "invalid operator '><' in prop (use <, >, <=, >=, ==, !=)");
    EXPECT_EQ(error("around -1 name P"), "around: distance must not be negative");
    EXPECT_EQ(error("sphlayer 1 -2 name P"), "sphlayer: radius must not be negative");
    EXPECT_EQ(error("AND"), "Unknown selection token: 'AND' (keywords are lower case: did you mean 'and'?)");
    EXPECT_EQ(error("mass :3"), "Selection failed: 'Failed to parse number: :3'");
    EXPECT_EQ(error("same RESID as name P").substr(0, 52),
              "Unknown 'same' property 'RESID' (choose one of: frag");
    // SEL-01: reversed operands mirror the operator
    EXPECT_EQ(selectAtoms(*dump.topology, frame, "prop 1.2 > x"),
              selectAtoms(*dump.topology, frame, "prop x < 1.2"));
    EXPECT_EQ(selectAtoms(*dump.topology, frame, "prop abs 3 <= z"),
              selectAtoms(*dump.topology, frame, "prop abs z >= 3"));
    // SEL-09 / SEL-10
    EXPECT_EQ(selectAtoms(*dump.topology, frame, "mass -inf:+inf").size(), dump.topology->atomCount());
    EXPECT_TRUE(selectAtoms(*dump.topology, frame, "resid 3:1").empty());
    EXPECT_EQ(selectAtoms(*dump.topology, frame, "atom \\A 1 P"), selectAtoms(*dump.topology, frame, "atom A 1 P"));
}

// SEL-06: a cutoff larger than half the cell (here a strongly skewed one, where MDAnalysis' grid
// search refused it) gets exact minimum-image distances.
TEST(AtomSelection, CutoffsBeyondTheBox) {
    auto topology = std::make_shared<Topology>();
    Topology& t = *topology;
    constexpr std::size_t kAtoms = 24;
    for (std::size_t a = 0; a < kAtoms; ++a) {
        t.names.push_back(t.strings.intern("C"));
        t.atomResidue.push_back(static_cast<std::int32_t>(a));
        t.resids.push_back(static_cast<std::int64_t>(a + 1));
        t.resnames.push_back(t.strings.intern("X"));
        t.residueSegment.push_back(0);
    }
    t.segids = {t.strings.intern("A")};
    t.finish();
    Frame frame;
    for (std::size_t a = 0; a < kAtoms; ++a) {
        frame.positions.push_back(static_cast<float>((a * 7) % 10));
        frame.positions.push_back(static_cast<float>(a % 3) * 0.5f);
        frame.positions.push_back(static_cast<float>((a * 3) % 10));
    }
    frame.box = std::array<float, 6>{10, 10, 10, 90, 90, 170};
    // the cell vectors as the selections build them (b in the xy plane at 170 degrees from a)
    const double gamma = 170.0 * std::acos(-1.0) / 180.0;
    const double cell[3][3] = {{10, 0, 0},
                               {static_cast<float>(10 * std::cos(gamma)), static_cast<float>(10 * std::sin(gamma)), 0},
                               {0, 0, 10}};
    auto nearest = [&](std::size_t i, std::size_t j) {
        double best = 1e300;
        for (int u = -12; u <= 12; ++u) {
            for (int v = -12; v <= 12; ++v) {
                for (int w = -2; w <= 2; ++w) {
                    double d2 = 0;
                    for (std::size_t c = 0; c < 3; ++c) {
                        const double d = static_cast<double>(frame.positions[3 * j + c]) -
                                         static_cast<double>(frame.positions[3 * i + c]) + u * cell[0][c] +
                                         v * cell[1][c] + w * cell[2][c];
                        d2 += d * d;
                    }
                    best = std::min(best, d2);
                }
            }
        }
        return std::sqrt(best);
    };
    for (const double cutoff : {0.4, 0.9, 2.5}) {
        std::vector<std::int32_t> expected;
        for (std::size_t j = 12; j < kAtoms; ++j) {
            for (std::size_t i = 0; i < 12; ++i) {
                if (nearest(i, j) <= cutoff) {
                    expected.push_back(static_cast<std::int32_t>(j));
                    break;
                }
            }
        }
        std::ostringstream sel;
        sel << "around " << cutoff << " index 0:11";
        EXPECT_EQ(selectAtoms(t, &frame, sel.str()), expected) << sel.str();
    }
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
    EXPECT_LT(seconds, 1.0 * mdbridge::test::timeScale());  // MDBRIDGE_TIME_SCALE: e.g. the sanitizer build
#else
    EXPECT_LT(seconds, 20.0 * mdbridge::test::timeScale());
#endif
}
