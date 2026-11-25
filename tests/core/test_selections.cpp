// MDATK's selection layer (mdbridge/core/selections.hpp) against mdatk/selections.py and the
// selection / frame helpers of mdatk/io/pdb_export.py (tests/data/selections/selections_cases.json,
// tests/parity/make_selections_cases.py), plus the examples of MDATK's
// tests/test_config_selections.py.

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/suggest.hpp"
#include "mdbridge/core/text.hpp"
#include "topology_dump.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

const Value& cases() {
    static const Value data = readJson(kData / "selections" / "selections_cases.json");
    return data;
}

const Array& list(const Value& v, std::string_view key) { return *v.find(key)->get<Array>(); }

std::optional<std::string> errorOf(const Value& c) {
    if (const Value* e = c.find("error")) {
        return *e->get<std::string>();
    }
    return std::nullopt;
}

bool hugeInt(const Value& v) { return v.isTable(); }  // {"$int": "..."}

std::vector<std::string> strings(const Value& v) {
    std::vector<std::string> out;
    for (const Value& x : *v.get<Array>()) {
        out.push_back(*x.get<std::string>());
    }
    return out;
}

template <class T>
std::vector<T> ints(const Value& v) {
    std::vector<T> out;
    for (const Value& x : *v.get<Array>()) {
        out.push_back(static_cast<T>(*x.get<std::int64_t>()));
    }
    return out;
}

struct Fixture {
    TopologyDump dump;
    const Topology* topology() const { return dump.topology.get(); }
    const Frame* frame() const { return dump.frame ? &*dump.frame : nullptr; }
};

const Fixture& fixture(const std::string& name) {
    static std::map<std::string, Fixture> loaded;
    auto it = loaded.find(name);
    if (it == loaded.end()) {
        it = loaded.emplace(name, Fixture{topologyFixture(name)}).first;
    }
    return it->second;
}

const Topology* topologyOf(const Value& c) {
    const Value* f = c.find("fixture");
    return f && f->isString() ? fixture(*f->get<std::string>()).topology() : nullptr;
}

std::vector<ResidueGroup> groupsOf(const Value& set) {
    std::vector<ResidueGroup> out;
    for (const Value& g : list(set, "groups")) {
        ResidueGroup group;
        group.name = text(g, "name");
        group.residues = text(g, "residues");
        group.selection = text(g, "selection");
        out.push_back(group);
    }
    return out;
}

// Runs fn and compares its string result (or SelectionError / ValueError text) with `expected`.
template <class F>
void expectText(const Value& expected, F&& fn, const std::string& what) {
    try {
        const std::string got = fn();
        ASSERT_FALSE(errorOf(expected))
            << what << ": got '" << got << "', expected error " << *errorOf(expected);
        EXPECT_EQ(got, text(expected, "result")) << what;
    } catch (const ValueError& e) {
        ASSERT_TRUE(errorOf(expected)) << what << ": unexpected error " << e.what();
        EXPECT_EQ(e.what(), *errorOf(expected)) << what;
    }
}

}  // namespace

TEST(Selections, MdatkExamples) {
    EXPECT_EQ(parseResidueRanges("1-40, 45 ,50 - 60"),
              (std::vector<ResidueRange>{
                  {std::nullopt, 1, 40}, {std::nullopt, 45, 45}, {std::nullopt, 50, 60}}));
    EXPECT_EQ(parseResidueRanges("A:1-4;B:7"),
              (std::vector<ResidueRange>{{"A", 1, 4}, {"B", 7, 7}}));
    EXPECT_EQ(rangesToSelection("1-3, 7"), "(resid 1:3) or (resid 7)");
    EXPECT_THROW(parseResidueRanges("1-x"), SelectionError);
    const std::vector<ResidueGroup> groups = {{"A", "1-5", "", ""},
                                              {"B", "", "group:A and name P", ""},
                                              {"Loop1", "", "group:Loop2", ""},
                                              {"Loop2", "", "group:Loop1", ""}};
    EXPECT_EQ(expandGroups("group:B or resid 9", groups), "(((resid 1:5)) and name P) or resid 9");
    EXPECT_EQ(resolveSelection("A", groups), "((resid 1:5))");
    EXPECT_EQ(resolveSelection("sel:name P", groups), "name P");
    EXPECT_EQ(parseFrameSpec("0, 10, 20-50:10, last, -1", 100),
              (std::vector<std::int64_t>{0, 10, 20, 30, 40, 50, 99}));
    EXPECT_THROW(parseFrameSpec("100", 100), SelectionError);
    EXPECT_EQ(parseTimeSpec("1.5ns, 2000ps, 3", "ns"), (std::vector<double>{1500, 2000, 3000}));
    EXPECT_EQ(guessElement("Na+", "Na+", true), "NA");
}

TEST(Selections, ResidueRanges) {
    for (const Value& c : list(cases(), "parse_residue_ranges")) {
        const std::string& input = text(c, "text");
        if (const auto error = errorOf(c)) {
            try {
                parseResidueRanges(input);
                ADD_FAILURE() << input << ": no error";
            } catch (const SelectionError& e) {
                EXPECT_EQ(e.what(), *error) << input;
            }
            continue;
        }
        bool huge = false;
        std::vector<ResidueRange> expected;
        for (const Value& item : list(c, "result")) {
            const Array& t = *item.get<Array>();
            huge = huge || hugeInt(t[1]) || hugeInt(t[2]);
            if (!huge) {
                expected.push_back({t[0].isNull()
                                        ? std::nullopt
                                        : std::optional<std::string>(*t[0].get<std::string>()),
                                    *t[1].get<std::int64_t>(), *t[2].get<std::int64_t>()});
            }
        }
        if (huge) {  // beyond int64: MDBridge refuses (rangesToSelection keeps them, see below)
            EXPECT_THROW(parseResidueRanges(input), SelectionError) << input;
        } else {
            EXPECT_EQ(parseResidueRanges(input), expected) << input;
        }
    }
    for (const Value& c : list(cases(), "ranges_to_selection")) {
        expectText(
            c, [&] { return rangesToSelection(text(c, "text"), topologyOf(c)); }, text(c, "text"));
    }
}

TEST(Selections, GroupsExpandAndResolve) {
    const Table& sets = *cases().find("groups")->get<Table>();
    for (const auto& [name, set] : sets) {
        const std::vector<ResidueGroup> groups = groupsOf(set);
        for (const Value& c : list(set, "cases")) {
            const Topology* topology = topologyOf(c);
            if (c.find("group")) {
                continue;
            }
            const std::string& input = text(c, "text");
            expectText(
                *c.find("expand"), [&] { return expandGroups(input, groups, topology); },
                name + " expand " + input);
            expectText(
                *c.find("resolve"), [&] { return resolveSelection(input, groups, topology); },
                name + " resolve " + input);
        }
        // group_selection of every group, in order
        std::size_t k = 0;
        for (const Value& c : list(set, "cases")) {
            if (!c.find("group")) {
                continue;
            }
            const ResidueGroup& group = groups[k % groups.size()];
            ++k;
            EXPECT_EQ(group.name, text(c, "group"));
            expectText(
                *c.find("group_selection"), [&] { return groupSelection(group, topologyOf(c)); },
                name + " group_selection " + group.name);
        }
    }
}

TEST(Selections, SelectAndExportSelection) {
    for (const Value& c : list(cases(), "select")) {
        const Fixture& f = fixture(text(c, "fixture"));
        const std::vector<ResidueGroup> groups =
            groupsOf(*cases().path({"groups", text(c, "groups")}));
        const std::string& input = text(c, "text");
        const std::string what = text(c, "fixture") + " " + text(c, "groups") + " '" + input + "'";
        try {
            AtomGroup atoms;
            if (c.find("export")) {
                System system;
                system.topology = f.dump.topology;
                if (f.frame()) {
                    system.trajectory =
                        std::make_unique<MemoryTrajectory>(std::vector<Frame>{*f.frame()});
                }
                atoms = resolveExportSelection(system, input, groups);
            } else {
                atoms = select(f.dump.topology, f.frame(), input, groups,
                               *c.find("allow_empty")->get<bool>());
            }
            ASSERT_FALSE(errorOf(c)) << what << ": no error, expected " << *errorOf(c);
            EXPECT_EQ(atoms.indices(), ints<std::int32_t>(*c.find("result"))) << what;
        } catch (const SelectionError& e) {
            ASSERT_TRUE(errorOf(c)) << what << ": unexpected error " << e.what();
            EXPECT_EQ(e.what(), *errorOf(c)) << what;
        }
    }
}

TEST(Selections, Elements) {
    for (const Value& c : list(cases(), "guess_element")) {
        const Array& args = list(c, "args");
        EXPECT_EQ(guessElement(*args[0].get<std::string>(), *args[1].get<std::string>(),
                               *args[2].get<bool>()),
                  text(c, "result"))
            << *args[0].get<std::string>() << " " << *args[1].get<std::string>();
    }
    for (const Value& c : list(cases(), "atoms")) {
        const Fixture& f = fixture(text(c, "fixture"));
        const AtomGroup atoms(f.dump.topology, ints<std::int32_t>(*c.find("indices")));
        const std::string what = text(c, "fixture") + " " + text(c, "subset");
        for (const char* key : {"elements", "heavy_mask", "heavy", "filter_NO", "filter_MG"}) {
            ASSERT_TRUE(c.path({key, "result"})) << what << " " << key;
        }
        EXPECT_EQ(atomElements(atoms), strings(*c.path({"elements", "result"}))) << what;
        std::vector<bool> mask;
        for (const Value& v : *c.path({"heavy_mask", "result"})->get<Array>()) {
            mask.push_back(*v.get<bool>());
        }
        EXPECT_EQ(heavyMask(atoms), mask) << what;
        EXPECT_EQ(heavy(atoms).indices(), ints<std::int32_t>(*c.path({"heavy", "result"}))) << what;
        EXPECT_EQ(elementFilter(atoms, {"n", " O ", ""}).indices(),
                  ints<std::int32_t>(*c.path({"filter_NO", "result"})))
            << what;
        EXPECT_EQ(elementFilter(atoms, {"Mg", "na", "cl", "K"}).indices(),
                  ints<std::int32_t>(*c.path({"filter_MG", "result"})))
            << what;
        EXPECT_TRUE(elementFilter(atoms, {}).empty()) << what;
    }
}

TEST(Selections, LabelsAndCsv) {
    for (const Value& c : list(cases(), "residue_labels")) {
        const Fixture& f = fixture(text(c, "fixture"));
        const auto residues = ints<std::int32_t>(*c.find("residues"));
        EXPECT_EQ(residueLabels(*f.topology(), residues), strings(*c.find("result")))
            << text(c, "fixture");
    }
    for (const Value& c : list(cases(), "split_csv")) {
        const Value& input = *c.find("text");
        const std::string in = input.isNull() ? "" : *input.get<std::string>();
        EXPECT_EQ(splitCsv(in), strings(*c.find("result"))) << in;
    }
}

TEST(Selections, FrameAndTimeSpecs) {
    for (const Value& c : list(cases(), "parse_frame_spec")) {
        const std::string& input = text(c, "text");
        const std::int64_t n = *c.find("n")->get<std::int64_t>();
        try {
            const auto frames = parseFrameSpec(input, n);
            ASSERT_FALSE(errorOf(c)) << input << ": no error, expected " << *errorOf(c);
            EXPECT_EQ(frames, ints<std::int64_t>(*c.find("result"))) << input;
        } catch (const SelectionError& e) {
            ASSERT_TRUE(errorOf(c)) << input << ": unexpected error " << e.what();
            EXPECT_EQ(e.what(), *errorOf(c)) << input;
        }
    }
    for (const Value& c : list(cases(), "parse_time_spec")) {
        const std::string& input = text(c, "text");
        const std::string& unit = text(c, "unit");
        try {
            const auto times = parseTimeSpec(input, unit);
            ASSERT_FALSE(errorOf(c)) << input << ": no error, expected " << *errorOf(c);
            const Array& expected = list(c, "result");
            ASSERT_EQ(times.size(), expected.size()) << input;
            for (std::size_t i = 0; i < times.size(); ++i) {
                EXPECT_EQ(times[i], detail::number(expected[i])) << input;
            }
        } catch (const ValueError& e) {
            ASSERT_TRUE(errorOf(c)) << input << " " << unit << ": unexpected error " << e.what();
            EXPECT_EQ(e.what(), *errorOf(c)) << input;
        }
    }
    for (const Value& c : list(cases(), "time_unit_factor")) {
        const std::string& unit = text(c, "unit");
        if (const auto error = errorOf(c)) {
            try {
                timeUnitFactor(unit);
                ADD_FAILURE() << unit << ": no error";
            } catch (const ValueError& e) {
                EXPECT_EQ(e.what(), *error);
            }
        } else {
            EXPECT_EQ(timeUnitFactor(unit), detail::number(*c.find("result"))) << unit;
        }
    }
}

TEST(Selections, FramesFromSpecs) {
    for (const Value& c : list(cases(), "frames_from_specs")) {
        std::vector<Frame> frames;
        for (const Value& t : list(c, "times")) {
            Frame f;
            f.time = detail::number(t);
            frames.push_back(f);
        }
        MemoryTrajectory trajectory(std::move(frames), detail::number(*c.find("dt")), true);
        std::optional<double> dtPs;
        if (const Value* dt = c.find("dt_ps"); dt && !dt->isNull()) {
            dtPs = detail::number(*dt);
        }
        const std::string what =
            text(c, "frames") + " | " + text(c, "times_text") + " | " + text(c, "unit");
        LogCapture log;
        try {
            const auto got = framesFromSpecs(trajectory, text(c, "frames"), text(c, "times_text"),
                                             text(c, "unit"), dtPs);
            ASSERT_FALSE(errorOf(c)) << what << ": no error, expected " << *errorOf(c);
            EXPECT_EQ(got, ints<std::int64_t>(*c.find("result"))) << what;
        } catch (const ValueError& e) {
            ASSERT_TRUE(errorOf(c)) << what << ": unexpected error " << e.what();
            EXPECT_EQ(e.what(), *errorOf(c)) << what;
        }
        std::vector<std::string> warnings;
        for (const std::string& w : strings(*c.find("warnings"))) {
            warnings.push_back("WARNING: " + w);
        }
        EXPECT_EQ(log.lines(), warnings) << what;
        EXPECT_EQ(trajectory.current().index, 0) << what;
    }
}

// int() and float() of str (mdbridge/core/python.hpp), which the selection parsers rely on
TEST(Selections, PythonIntAndFloat) {
    for (const Value& c : list(cases(), "python_numbers")) {
        const std::string& input = text(c, "text");
        if (const Value* v = c.find("int")) {
            if (hugeInt(*v)) {
                EXPECT_THROW(pythonInt(input), std::overflow_error) << input;
                EXPECT_FALSE(tryPythonInt(input)) << input;
            } else {
                EXPECT_EQ(pythonInt(input), *v->get<std::int64_t>()) << input;
                EXPECT_EQ(tryPythonInt(input), *v->get<std::int64_t>()) << input;
            }
        } else {
            try {
                pythonInt(input);
                ADD_FAILURE() << input << ": no error";
            } catch (const ValueError& e) {
                EXPECT_EQ(e.what(), text(c, "int_error")) << input;
            }
            EXPECT_FALSE(tryPythonInt(input)) << input;
        }
        if (const Value* v = c.find("float")) {
            EXPECT_EQ(pythonFloatRepr(pythonFloat(input)), *v->get<std::string>()) << input;
            EXPECT_TRUE(tryPythonFloat(input)) << input;
        } else {
            try {
                pythonFloat(input);
                ADD_FAILURE() << input << ": no error";
            } catch (const ValueError& e) {
                EXPECT_EQ(e.what(), text(c, "float_error")) << input;
            }
            EXPECT_FALSE(tryPythonFloat(input)) << input;
        }
    }
}

// SEL-13: an atom without an element is guessed from its own name; the others keep theirs, whatever
// group they are in (MDATK guesses the whole group from names as soon as one element is missing).
TEST(Selections, ElementFallbackPerAtom) {
    auto topology = std::make_shared<Topology>();
    Topology& t = *topology;
    const char* names[] = {"CL1", "C2", "EP1", "MG"};
    const char* elements[] = {"Cl", "C", "", "Mg"};
    for (std::size_t a = 0; a < 4; ++a) {
        t.names.push_back(t.strings.intern(names[a]));
        t.atomResidue.push_back(a < 3 ? 0 : 1);
    }
    t.elements = std::vector<std::uint32_t>{};
    for (const char* e : elements) {
        t.elements->push_back(t.strings.intern(e));
    }
    t.resids = {1, 2};
    t.resnames = {t.strings.intern("LIG"), t.strings.intern("MG")};
    t.residueSegment = {0, 0};
    t.segids = {t.strings.intern("A")};
    t.finish();
    const AtomGroup all(topology, std::vector<std::int32_t>{0, 1, 2, 3});
    EXPECT_EQ(atomElements(all), (std::vector<std::string>{"CL", "C", "E", "MG"}));
    const AtomGroup ligand(topology, std::vector<std::int32_t>{0, 3});
    EXPECT_EQ(atomElements(ligand), (std::vector<std::string>{"CL", "MG"}));
    EXPECT_EQ(elementFilter(all, {"Cl"}).indices(), (std::vector<std::int32_t>{0}));
}

// SEL-19, SEL-20: "did you mean" (MDATK 0.1.x selections.suggest, Python's difflib)
TEST(Selections, SuggestLikeDifflib) {
    for (const Value& c : list(cases(), "similarity")) {
        const std::u32string a = decodeUtf8(text(c, "a")), b = decodeUtf8(text(c, "b"));
        EXPECT_EQ(similarity(a, b), *c.find("ratio")->get<double>()) << text(c, "a") << " / " << text(c, "b");
    }
    for (const Value& c : list(cases(), "suggest")) {
        const std::optional<std::string> got = suggestName(text(c, "name"), strings(*c.find("candidates")));
        const Value& expected = *c.find("result");
        if (expected.isNull()) {
            EXPECT_FALSE(got) << text(c, "name") << " -> " << got.value_or("");
        } else {
            EXPECT_EQ(got.value_or("<none>"), *expected.get<std::string>()) << text(c, "name");
        }
    }
}

// SEL-19 ... SEL-22 on a small example (the cases above cover the fixtures)
TEST(Selections, MessagesNameTheProblem) {
    const Fixture& f = fixture("fixture_chains.pdb");
    const std::vector<ResidueGroup> groups = {{"Aptamer", "A:1-4", "", ""},
                                              {"Target", "B:10-12", "", ""},
                                              {"Far", "900-910", "", ""},
                                              {"A", "", "group:B", ""},
                                              {"B", "", "group:A", ""},
                                              {"C", "", "group:Tagret", ""}};
    const auto message = [&](const std::string& input) {
        try {
            select(f.dump.topology, f.frame(), input, groups);
        } catch (const SelectionError& e) {
            return std::string(e.what());
        }
        return std::string("no error");
    };
    EXPECT_EQ(message("Target and name P"),
              "invalid selection 'Target and name P': Unknown selection token: 'Target' (to use the residue "
              "group write group:Target)");
    EXPECT_EQ(message("resnam G"),
              "invalid selection 'resnam G': Unknown selection token: 'resnam' (did you mean 'resname'?)");
    EXPECT_EQ(message("group:Aptamer and resnam G"),
              "invalid selection 'group:Aptamer and resnam G' (expanded: '((chainID A and resid 1:4)) and resnam "
              "G'): Unknown selection token: 'resnam' (did you mean 'resname'?)");
    EXPECT_EQ(message("group Target"), "invalid selection 'group Target': Selection failed: 'Failed to find "
                                       "group: Target' (write group:Target, with a colon, to use the residue group)");
    EXPECT_EQ(message("group:Tagret"),
              "unknown residue group 'Tagret' (defined: Aptamer, Target, Far, A, B, C); did you mean 'Target'?");
    EXPECT_EQ(message("C"), "unknown residue group 'Tagret' in the selection of group 'C' (defined: Aptamer, "
                            "Target, Far, A, B, C); did you mean 'Target'?");
    EXPECT_EQ(message("A"), "group definitions are nested in a circle: group:A \xe2\x86\x92 group:B "
                            "\xe2\x86\x92 group:A");
    const std::string far = message("Far");
    EXPECT_EQ(far.rfind("selection 'Far' matched no atoms: group 'Far' (residues '900-910') matches no atoms - "
                        "the structure's residues are numbered ",
                        0),
              0u)
        << far;
    EXPECT_EQ(message("name P and resid 999"),
              "selection 'name P and resid 999' matched no atoms: the part 'resid 999' matches no atoms");
    try {
        select(f.dump.topology, f.frame(), "group:Far", groups);
        ADD_FAILURE() << "no error";
    } catch (const SelectionError& e) {
        EXPECT_EQ(e.kind(), SelectionError::Kind::Empty);
        EXPECT_EQ(e.token(), "Far");
    }
}
