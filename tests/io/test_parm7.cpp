// Amber parm7 files: MDATK's section reader / writer and Lennard-Jones table (mdatk/io/parm7.py)
// and the topology MDAnalysis' TOPParser builds, against what they give for the files and values of
// tests/parity/make_parm7_cases.py (data/io/parm7, see data/README.md).

#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/io/parm7.hpp"
#include "mdbridge/io/warnings.hpp"
#include "topology_dump.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const fs::path kParm7 = kData / "io" / "parm7";

const core::Array& items(const core::Value& value, std::string_view key) {
    return *value.find(key)->get<core::Array>();
}

double hexFloat(const core::Value& v) {
    return std::strtod(v.get<std::string>()->c_str(), nullptr);
}

std::string replaceAll(std::string text, std::string_view from, const std::string& to) {
    for (std::size_t at = text.find(from); at != std::string::npos;
         at = text.find(from, at + to.size())) {
        text.replace(at, from.size(), to);
    }
    return text;
}

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

// Runs `fn` and compares what it throws with {"type", "message"} ("{path}" replaced).
template <class Fn>
void expectError(Fn&& fn, const core::Value& error, const std::string& path,
                 const std::string& context) {
    const std::string type = text(error, "type");
    const std::string message = replaceAll(text(error, "message"), "{path}", path);
    try {
        fn();
        ADD_FAILURE() << context << ": no error, expected " << message;
    } catch (const core::ValueError& e) {
        EXPECT_TRUE(type != "OverflowError") << context;
        EXPECT_EQ(e.what(), message) << context;
    } catch (const std::overflow_error& e) {
        EXPECT_EQ(type, "OverflowError") << context;
        EXPECT_EQ(e.what(), message) << context;
    }
}

void expectSections(const io::Parm7Sections& got, const core::Array& expected,
                    const std::string& context) {
    ASSERT_EQ(got.size(), expected.size()) << context;
    for (std::size_t s = 0; s < got.size(); ++s) {
        const core::Array& e = *expected[s].get<core::Array>();
        const std::string& kind = *e[1].get<std::string>();
        const core::Array& values = *e[2].get<core::Array>();
        EXPECT_EQ(got[s].first, *e[0].get<std::string>()) << context;
        if (kind == "a") {
            const auto* v = std::get_if<std::vector<std::string>>(&got[s].second);
            ASSERT_TRUE(v) << context << " " << got[s].first;
            ASSERT_EQ(v->size(), values.size()) << context;
            for (std::size_t k = 0; k < v->size(); ++k) {
                EXPECT_EQ((*v)[k], *values[k].get<std::string>()) << context << " " << k;
            }
        } else if (kind == "i") {
            const auto* v = std::get_if<std::vector<std::int64_t>>(&got[s].second);
            ASSERT_TRUE(v) << context << " " << got[s].first;
            ASSERT_EQ(v->size(), values.size()) << context;
            for (std::size_t k = 0; k < v->size(); ++k) {
                EXPECT_EQ((*v)[k], *values[k].get<std::int64_t>()) << context << " " << k;
            }
        } else {
            const auto* v = std::get_if<std::vector<double>>(&got[s].second);
            ASSERT_TRUE(v) << context << " " << got[s].first;
            ASSERT_EQ(v->size(), values.size()) << context;
            for (std::size_t k = 0; k < v->size(); ++k) {
                EXPECT_TRUE(sameDouble((*v)[k], hexFloat(values[k])))
                    << context << " " << k << ": " << (*v)[k] << " vs "
                    << *values[k].get<std::string>();
            }
        }
    }
}

std::optional<std::set<std::string>> wantedOf(const core::Value& c) {
    const core::Value* wanted = c.find("wanted");
    if (!wanted || wanted->isNull()) {
        return std::nullopt;
    }
    std::set<std::string> out;
    for (const core::Value& w : *wanted->get<core::Array>()) {
        out.insert(*w.get<std::string>());
    }
    return out;
}

io::Parm7Values valuesOf(const core::Value& tagged) {
    if (const core::Value* s = tagged.find("strings")) {
        std::vector<std::string> out;
        for (const core::Value& v : *s->get<core::Array>()) {
            out.push_back(*v.get<std::string>());
        }
        return out;
    }
    if (const core::Value* i = tagged.find("ints")) {
        std::vector<std::int64_t> out;
        for (const core::Value& v : *i->get<core::Array>()) {
            out.push_back(std::stoll(*v.get<std::string>()));
        }
        return out;
    }
    std::vector<double> out;
    for (const core::Value& v : items(tagged, "floats")) {
        out.push_back(hexFloat(v));
    }
    return out;
}

template <class T>
void expectOptional(const std::optional<T>& got, const std::optional<T>& expected, const char* what,
                    const std::string& context) {
    ASSERT_EQ(got.has_value(), expected.has_value()) << context << " " << what;
    if (got) {
        EXPECT_EQ(*got, *expected) << context << " " << what;
    }
}

std::vector<std::string> resolve(const core::StringPool& pool,
                                 const std::vector<std::uint32_t>& ids) {
    std::vector<std::string> out;
    for (const std::uint32_t id : ids) {
        out.push_back(pool[id]);
    }
    return out;
}

std::optional<std::vector<std::string>> resolve(
    const core::StringPool& pool, const std::optional<std::vector<std::uint32_t>>& ids) {
    if (!ids) {
        return std::nullopt;
    }
    return resolve(pool, *ids);
}

void expectSameTopology(const core::Topology& got, const core::Topology& expected,
                        const std::string& context) {
    const core::StringPool& g = got.strings;
    const core::StringPool& e = expected.strings;
    EXPECT_EQ(resolve(g, got.names), resolve(e, expected.names)) << context;
    EXPECT_EQ(resolve(g, got.types), resolve(e, expected.types)) << context;
    EXPECT_EQ(got.atomResidue, expected.atomResidue) << context;
    expectOptional(resolve(g, got.elements), resolve(e, expected.elements), "elements", context);
    expectOptional(resolve(g, got.chainIDs), resolve(e, expected.chainIDs), "chainIDs", context);
    expectOptional(got.masses, expected.masses, "masses", context);
    expectOptional(got.charges, expected.charges, "charges", context);
    expectOptional(got.ids, expected.ids, "ids", context);
    expectOptional(resolve(g, got.altLocs), resolve(e, expected.altLocs), "altLocs", context);
    expectOptional(resolve(g, got.recordTypes), resolve(e, expected.recordTypes), "recordTypes",
                   context);
    expectOptional(got.occupancies, expected.occupancies, "occupancies", context);
    expectOptional(got.tempFactors, expected.tempFactors, "tempFactors", context);
    expectOptional(got.formalCharges, expected.formalCharges, "formalCharges", context);
    EXPECT_EQ(got.resids, expected.resids) << context;
    EXPECT_EQ(got.resnums, expected.resnums) << context;
    EXPECT_EQ(resolve(g, got.resnames), resolve(e, expected.resnames)) << context;
    expectOptional(resolve(g, got.icodes), resolve(e, expected.icodes), "icodes", context);
    EXPECT_EQ(got.residueSegment, expected.residueSegment) << context;
    EXPECT_EQ(resolve(g, got.segids), resolve(e, expected.segids)) << context;
    expectOptional(got.bonds, expected.bonds, "bonds", context);
}

}  // namespace

TEST(Parm7, ReadSections) {
    const core::Value cases = readJson(kParm7 / "read_cases.json");
    TempDir dir;
    for (const core::Value& c : *cases.get<core::Array>()) {
        const std::string name = text(c, "name");
        fs::path path;
        if (const core::Value* file = c.find("file")) {
            path = kData / *file->get<std::string>();
        } else {
            path = dir / (name + ".parm7");
            writeFile(path, fromHex(text(c, "bytes")));
        }
        const auto wanted = wantedOf(c);
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::readParm7Sections(path, wanted); }, *error, path.string(), name);
            continue;
        }
        expectSections(io::readParm7Sections(path, wanted), items(c, "sections"), name);
    }
}

TEST(Parm7, FindSection) {
    const io::Parm7Sections sections{{"A", std::vector<std::int64_t>{1}},
                                     {"B", std::vector<double>{2.0}}};
    EXPECT_EQ(io::findSection(sections, "B"), &sections[1].second);
    EXPECT_EQ(io::findSection(sections, "C"), nullptr);
}

TEST(Parm7, MissingFile) {
    EXPECT_THROW(io::readParm7Sections(kParm7 / "no-such-file.parm7"), fs::filesystem_error);
}

TEST(Parm7, FormatValues) {
    const core::Value cases = readJson(kParm7 / "writer_cases.json");
    for (const core::Value& c : items(cases, "format_values")) {
        const char kind = text(c, "kind")[0];
        std::vector<std::string> expected;
        for (const core::Value& line : items(c, "lines")) {
            expected.push_back(*line.get<std::string>());
        }
        EXPECT_EQ(io::formatParm7Values(kind, valuesOf(*c.find("values"))), expected) << kind;
    }
}

TEST(Parm7, Write) {
    const core::Value cases = readJson(kParm7 / "writer_cases.json");
    TempDir dir;
    for (const core::Value& c : items(cases, "write")) {
        const std::string name = text(c, "name");
        std::vector<io::Parm7Section> sections;
        for (const core::Value& s : items(c, "sections")) {
            sections.push_back({text(s, "flag"), text(s, "kind")[0], valuesOf(*s.find("values"))});
        }
        const fs::path path = dir / (name + ".parm7");
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::writeParm7(path, sections, text(c, "title")); }, *error,
                        path.string(), name);
        } else {
            io::writeParm7(path, sections, text(c, "title"));
        }
        EXPECT_EQ(readFile(path), fromHex(text(c, "bytes"))) << name;
    }
}

TEST(Parm7, WriteReadRoundTrip) {
    TempDir dir;
    const std::vector<io::Parm7Section> sections{
        {"POINTERS", 'i', std::vector<std::int64_t>{3, -2, 1234567}},
        {"ATOM_NAME", 'a', std::vector<std::string>{"N", "CA", "HO2'"}},
        {"CHARGE", 'e', std::vector<double>{0.5, -1.25e-7, 3e300}},
        {"NOTE", 's', std::vector<std::string>{"free text"}},
    };
    io::writeParm7(dir / "x.parm7", sections, "rt");
    const io::Parm7Sections read = io::readParm7Sections(dir / "x.parm7");
    ASSERT_EQ(read.size(), 5u);
    EXPECT_EQ(std::get<std::vector<std::string>>(read[0].second)[0], "rt");
    EXPECT_EQ(std::get<std::vector<std::int64_t>>(read[1].second),
              (std::vector<std::int64_t>{3, -2, 1234567}));
    EXPECT_EQ(std::get<std::vector<std::string>>(read[2].second),
              (std::vector<std::string>{"N", "CA", "HO2'"}));
    EXPECT_EQ(std::get<std::vector<double>>(read[3].second),
              (std::vector<double>{0.5, -1.25e-7, 3e300}));
    EXPECT_EQ(std::get<std::vector<std::string>>(read[4].second),
              (std::vector<std::string>{"free text"}));
}

TEST(Parm7, LjTable) {
    const core::Value cases = readJson(kParm7 / "lj_cases.json");
    const io::LjTable table = io::LjTable::fromParm7(kData / text(cases, "file"));
    EXPECT_EQ(table.ntypes, *cases.find("ntypes")->get<std::int64_t>());
    std::vector<std::int64_t> typeIndex, nbIndex;
    for (const core::Value& v : items(cases, "type_index")) {
        typeIndex.push_back(*v.get<std::int64_t>());
    }
    for (const core::Value& v : items(cases, "nb_index")) {
        nbIndex.push_back(*v.get<std::int64_t>());
    }
    EXPECT_EQ(table.typeIndex, typeIndex);
    EXPECT_EQ(table.nbIndex, nbIndex);
    ASSERT_EQ(table.acoef.size(), items(cases, "acoef").size());
    for (std::size_t k = 0; k < table.acoef.size(); ++k) {
        EXPECT_EQ(table.acoef[k], hexFloat(items(cases, "acoef")[k]));
        EXPECT_EQ(table.bcoef[k], hexFloat(items(cases, "bcoef")[k]));
    }

    const core::Value& pairs = *cases.find("pairs");
    std::vector<std::int64_t> i, j;
    for (std::size_t k = 0; k < items(pairs, "i").size(); ++k) {
        i.push_back(*items(pairs, "i")[k].get<std::int64_t>());
        j.push_back(*items(pairs, "j")[k].get<std::int64_t>());
    }
    const auto [a, b] = table.pairCoefficients(i, j);
    for (std::size_t k = 0; k < i.size(); ++k) {
        EXPECT_EQ(a[k], hexFloat(items(pairs, "a")[k])) << k;
        EXPECT_EQ(b[k], hexFloat(items(pairs, "b")[k])) << k;
    }
    const std::vector<std::int64_t> outside{static_cast<std::int64_t>(table.typeIndex.size())};
    EXPECT_THROW(table.pairCoefficients(outside, outside), std::out_of_range);

    TempDir dir;
    for (const core::Value& c : items(cases, "errors")) {
        const fs::path path = dir / "lj_missing.parm7";
        writeFile(path, fromHex(text(c, "bytes")));
        expectError([&] { io::LjTable::fromParm7(path); }, *c.find("error"), path.string(), "lj");
    }
}

// MDATK's tests/test_analyses.py::test_parm7_lj_table_roundtrip: the demo's coefficients follow
// from the Lorentz-Berthelot rules applied to mdatk.demo.LJ_TYPES.
TEST(Parm7, LjTableDemoRoundTrip) {
    const core::Value cases = readJson(kParm7 / "lj_cases.json");
    const io::LjTable table = io::LjTable::fromParm7(kData / "topology" / "demo_dry.parm7");
    const auto demo = readTopologyDump(kData / "topology" / "demo_dry.parm7.topology.json");
    const core::Topology& topology = *demo.topology;
    const std::vector<std::int64_t> i{0, 5, 30, 100};
    const std::vector<std::int64_t> j{1, 7, 64, 200};
    const auto [a, b] = table.pairCoefficients(i, j);
    const core::Value& ljTypes = *cases.find("lj_types");
    for (std::size_t k = 0; k < i.size(); ++k) {
        const std::string& ti = topology.strings[topology.types[static_cast<std::size_t>(i[k])]];
        const std::string& tj = topology.strings[topology.types[static_cast<std::size_t>(j[k])]];
        const core::Array& pi = *ljTypes.find(ti)->get<core::Array>();
        const core::Array& pj = *ljTypes.find(tj)->get<core::Array>();
        const double rmin = hexFloat(pi[0]) + hexFloat(pj[0]);
        const double eps = std::sqrt(hexFloat(pi[1]) * hexFloat(pj[1]));
        EXPECT_NEAR(a[k], eps * std::pow(rmin, 12), 1e-6 * eps * std::pow(rmin, 12));
        EXPECT_NEAR(b[k], 2 * eps * std::pow(rmin, 6), 1e-6 * 2 * eps * std::pow(rmin, 6));
    }
}

TEST(Parm7, TopologyOfDemo) {
    const auto expected = readTopologyDump(kData / "topology" / "demo_dry.parm7.topology.json");
    const auto got = io::readParm7Topology(kData / "topology" / "demo_dry.parm7");
    expectSameTopology(*got, *expected.topology, "demo_dry");
}

TEST(Parm7, TopologyVariants) {
    // LD-03: the reader warnings of each variant (none for the others). chain_short / chain_long:
    // RESIDUE_CHAINID with 7 / 9 values for the 8 residues of make_parm7_cases.py's system.
    const auto chainCount = [](const fs::path& file, int values) {
        const std::u8string text = file.u8string();  // the path as given (UTF-8)
        return std::string(text.begin(), text.end()) + ": %FLAG RESIDUE_CHAINID has " + std::to_string(values) +
               " values for 8 residues; chain IDs are not read";
    };
    const std::map<std::string, int> chainValues = {{"chain_short.parm7", 7}, {"chain_long.parm7", 9}};
    int compared = 0;
    for (const auto& entry : fs::directory_iterator(kParm7 / "topology")) {
        const fs::path dump = entry.path();
        if (dump.extension() != ".json") {
            continue;
        }
        const fs::path file = dump.parent_path() / dump.stem().stem();  // x.parm7.topology.json
        const auto expected = readTopologyDump(dump);
        try {
            io::ReaderWarningScope warnings;
            const auto got = io::readParm7Topology(file);
            expectSameTopology(*got, *expected.topology, file.filename().string());
            const auto values = chainValues.find(file.filename().string());
            EXPECT_EQ(warnings.take(), values == chainValues.end()
                                           ? std::vector<std::string>{}
                                           : std::vector<std::string>{chainCount(file, values->second)})
                << file.filename().string();
        } catch (const std::exception& e) {
            ADD_FAILURE() << file.filename().string() << ": " << e.what();
        }
        ++compared;
    }
    EXPECT_GE(compared, 25);
}

TEST(Parm7, TopologyErrors) {
    const core::Value cases = readJson(kParm7 / "topology_errors.json");
    for (const core::Value& c : *cases.get<core::Array>()) {
        const fs::path path = kParm7 / "topology" / text(c, "file");
        const core::Value& error = *c.find("error");
        if (text(error, "type") ==
            "StopIteration") {  // Python's message is empty; MDBridge names the file
            EXPECT_THROW(io::readParm7Topology(path), core::ValueError) << text(c, "file");
            continue;
        }
        expectError([&] { io::readParm7Topology(path); }, error, path.string(), text(c, "file"));
    }
}

TEST(Parm7, IsParm7) {
    const core::Value cases = readJson(kParm7 / "is_parm7_cases.json");
    for (const core::Value& c : *cases.get<core::Array>()) {
        const core::Array& pair = *c.get<core::Array>();
        EXPECT_EQ(io::isParm7(*pair[0].get<std::string>()), *pair[1].get<bool>())
            << *pair[0].get<std::string>();
    }
}

namespace {

// small.parm7 with `from` replaced by `to` (once), read as a topology; the error text without the
// path, or "" when it loads.
std::string readVariant(const std::string& from, const std::string& to,
                        std::shared_ptr<core::Topology>* out = nullptr) {
    std::string text = readFile(kParm7 / "topology" / "small.parm7");
    const std::size_t at = text.find(from);
    EXPECT_NE(at, std::string::npos) << from;
    text.replace(at, from.size(), to);
    TempDir dir;
    const fs::path path = dir / "v.parm7";
    writeFile(path, text);
    try {
        auto topology = io::readParm7Topology(path);
        if (out) {
            *out = topology;
        }
        return "";
    } catch (const core::ValueError& e) {
        std::string message = e.what();
        const std::u8string p = path.u8string();
        const std::string ps(p.begin(), p.end());
        if (message.compare(0, ps.size(), ps) == 0) {
            message = "{path}" + message.substr(ps.size());
        }
        return message;
    }
}

}  // namespace

// MDBridge's parm7 rules beyond MDAnalysis (docs/DIFFERENCES.md, RD-P05 ... RD-P10).
TEST(Parm7, ImprovedTopologyRules) {
    // RD-P07: Fortran D exponents; RD-P12: charges divided by 18.2223 in float64
    std::shared_ptr<core::Topology> top;
    EXPECT_EQ(readVariant("  1.82223000E+01 -1.82223000E+01", "  1.82223000D+01 -1.82223000d+01", &top), "");
    ASSERT_TRUE(top && top->charges);
    EXPECT_EQ((*top->charges)[2], 18.2223 / 18.2223);
    EXPECT_EQ((*top->charges)[3], -1.0);
    EXPECT_EQ((*top->charges)[4], -3.70815013 / 18.2223);
    // RD-P10: bond and dihedral indices
    EXPECT_EQ(readVariant("       0       3       1       6", "       0       4       1       6"),
              "{path}: %FLAG BONDS_INC_HYDROGEN entry 1: atom index 4 is not a multiple of 3");
    EXPECT_EQ(readVariant("       0       3       1       6", "       0     120       1       6"),
              "{path}: %FLAG BONDS_INC_HYDROGEN entry 1: atom index 120 refers to atom 41 of 40");
    EXPECT_EQ(readVariant("       0       3       1       6", "      -3       3       1       6"),
              "{path}: %FLAG BONDS_INC_HYDROGEN entry 1: atom index -3 is negative");
    EXPECT_EQ(readVariant("      21       1      24      27       1\n",
                          "      21       1      24      27        \n"),
              "{path}: %FLAG BONDS_INC_HYDROGEN entry 5 has 2 of 3 values");
    // negative third and fourth dihedral indices are Amber flags, not errors
    EXPECT_EQ(readVariant("      30      33     -36      39", "      30      33     -36     -39"), "");
    EXPECT_EQ(readVariant("      30      33     -36      39", "      30     -33      36      39"),
              "{path}: %FLAG DIHEDRALS_WITHOUT_HYDROGEN entry 2: atom index -33 is negative");
    // RD-P05: a blank numeric field inside a line
    EXPECT_EQ(readVariant(" -0.00000000E+00  1.82223000E+01", "                  1.82223000E+01"),
              "{path} line 17: blank field 2 in %FLAG CHARGE");
    // RD-P02: the file ends inside a section
    {
        std::string text = readFile(kParm7 / "topology" / "small.parm7");
        text.resize(text.find("%FLAG MASS") + 30);
        TempDir dir;
        writeFile(dir / "cut.parm7", text);
        try {
            io::readParm7Topology(dir / "cut.parm7");
            ADD_FAILURE();
        } catch (const core::ValueError& e) {
            const std::string message = e.what();
            EXPECT_NE(message.find(": unexpected end of file in %FLAG MASS (expected 8 data lines, "
                                   "found 1)"),
                      std::string::npos)
                << message;
        }
    }
    // RD-P03: POINTERS
    {
        const std::string text = readFile(kParm7 / "topology" / "small.parm7");
        const std::size_t first = text.find("%FLAG POINTERS\n%FORMAT(10I8)\n") + 29;
        const std::string block = text.substr(first, text.find("%FLAG", first) - first);
        EXPECT_EQ(readVariant(block, "      40       3       5      10       2       1       1       3       0       0\n"),
                  "{path}: %FLAG POINTERS needs at least 12 integers (found 10)");
    }
    EXPECT_EQ(readVariant("      40       3       5      10", "      40       x       5      10"),
              "{path} line 7: %FLAG POINTERS value 'x' is not an integer");
}
