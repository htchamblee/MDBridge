// MDATK's PDB tools (mdatk/pdbtools.py read_pdb / write_pdb / parse_residue_name,
// annotate._insert_records, builder.transform_pdb_text / _tail) against what they give for the
// inputs of tests/parity/make_pdb_cases.py (data/io/pdb, see data/README.md), byte for byte.

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/io/pdb.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const core::Value& cases() {
    static const core::Value data = readJson(kData / "io" / "pdb" / "pdb_cases.json");
    return data;
}

const core::Array& list(const core::Value& value) {
    return *value.get<core::Array>();
}

const core::Array& items(const core::Value& value, std::string_view key) {
    return list(*value.find(key));
}

const std::string& str(const core::Value& value) {
    return *value.get<std::string>();
}

std::int64_t integer(const core::Value& value) {
    return *value.get<std::int64_t>();
}

double hexFloat(const core::Value& v) {
    const std::string& text = str(v);
    if (text == "nan") {
        return std::nan("");
    }
    return std::strtod(text.c_str(), nullptr);
}

std::vector<std::string> strings(const core::Value& value) {
    std::vector<std::string> out;
    for (const core::Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}

std::string pythonName(const core::ValueError& e) {
    return e.pythonName();
}

template <class Fn>
void expectError(Fn&& fn, const core::Value& error, const std::string& context) {
    try {
        fn();
        ADD_FAILURE() << context << ": no error, expected " << text(error, "message");
    } catch (const core::ValueError& e) {
        EXPECT_EQ(e.what(), text(error, "message")) << context;
        EXPECT_EQ(pythonName(e), text(error, "type")) << context;
    }
}

void expectStructure(const io::PdbStructure& s, const core::Value& expected,
                     const std::string& context) {
    EXPECT_EQ(s.header, strings(*expected.find("header"))) << context;
    EXPECT_EQ(s.trailer, strings(*expected.find("trailer"))) << context;
    EXPECT_EQ(s.conect, strings(*expected.find("conect"))) << context;
    EXPECT_EQ(s.model, text(expected, "model")) << context;
    EXPECT_EQ(s.notes, strings(*expected.find("notes"))) << context;
    EXPECT_EQ(s.chainIds(), strings(*expected.find("chain_ids"))) << context;
    const core::Array& segments = items(expected, "segments");
    ASSERT_EQ(s.segments.size(), segments.size()) << context;
    for (std::size_t k = 0; k < segments.size(); ++k) {
        const core::Array& residues = list(segments[k]);
        ASSERT_EQ(s.segments[k].size(), residues.size()) << context;
        for (std::size_t r = 0; r < residues.size(); ++r) {
            const io::PdbResidue& res = s.segments[k][r];
            const core::Array& e = list(residues[r]);
            EXPECT_EQ(res.name, str(e[0])) << context;
            EXPECT_EQ(res.chain, str(e[1])) << context;
            EXPECT_EQ(res.resseq, integer(e[2])) << context;
            EXPECT_EQ(res.icode, str(e[3])) << context;
            EXPECT_EQ(res.label(), str(e[4])) << context;
            const auto coords = res.coords();
            const core::Array& expectedCoords = list(e[5]);
            ASSERT_EQ(coords.size(), expectedCoords.size()) << context;
            for (std::size_t c = 0; c < coords.size(); ++c) {
                const core::Array& entry = list(expectedCoords[c]);
                EXPECT_EQ(coords[c].first, str(entry[0])) << context;
                for (std::size_t i = 0; i < 3; ++i) {
                    EXPECT_EQ(coords[c].second[i], hexFloat(entry[i + 1])) << context;
                }
            }
            const core::Array& atoms = list(e[6]);
            ASSERT_EQ(res.atoms.size(), atoms.size()) << context;
            for (std::size_t a = 0; a < atoms.size(); ++a) {
                const io::PdbAtom& atom = res.atoms[a];
                const core::Array& f = list(atoms[a]);
                EXPECT_EQ(atom.record, str(f[0])) << context;
                EXPECT_EQ(atom.name, str(f[1])) << context;
                EXPECT_EQ(atom.alt, str(f[2])) << context;
                EXPECT_EQ(atom.tail, str(f[3])) << context;
                EXPECT_EQ(atom.element, str(f[4])) << context;
                EXPECT_EQ(atom.charge, str(f[5])) << context;
                if (f[6].isNull()) {
                    EXPECT_FALSE(atom.serial.has_value()) << context;
                } else {
                    EXPECT_EQ(atom.serial, integer(f[6])) << context;
                }
                EXPECT_EQ(atom.anisou, str(f[7])) << context;
            }
        }
    }
}

}  // namespace

TEST(Pdb, ReadAndWrite) {
    int n = 0;
    for (const core::Value& c : items(cases(), "read")) {
        const std::string input = text(c, "text");
        const std::string context =
            "read case " + std::to_string(n++) + ": " + input.substr(0, 200);
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::readPdbText(input); }, *error, context);
            continue;
        }
        const io::PdbStructure s = io::readPdbText(input);
        expectStructure(s, *c.find("structure"), context);
        EXPECT_EQ(io::writePdb(s), text(c, "write")) << context;
        const core::Value* conect = c.find("write_conect");
        EXPECT_EQ(io::writePdb(s, false), conect ? str(*conect) : text(c, "write")) << context;
        std::vector<io::SeqresChain> chains;
        for (const core::Value& chain : items(c, "seqres_chains")) {
            chains.push_back({str(list(chain)[0]), strings(list(chain)[1])});
        }
        EXPECT_EQ(io::writePdb(s, false, chains), text(c, "write_seqres")) << context;
    }
}

TEST(Pdb, ReadFile) {
    TempDir dir;
    int n = 0;
    for (const core::Value& c : items(cases(), "files")) {
        const fs::path path = dir / ("case" + std::to_string(n) + ".pdb");
        writeFile(path, fromHex(text(c, "bytes")));
        const std::string context = "file case " + std::to_string(n++);
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::readPdbFile(path); }, *error, context);
            continue;
        }
        expectStructure(io::readPdbFile(path), *c.find("structure"), context);
    }
    EXPECT_THROW(io::readPdbFile(dir / "missing.pdb"), fs::filesystem_error);
}

TEST(Pdb, ParseResidueName) {
    for (const core::Value& c : items(cases(), "residue_names")) {
        const std::string& name = str(list(c)[0]);
        const auto parsed = io::parseResidueName(name);
        const core::Value& expected = list(c)[1];
        if (expected.isNull()) {
            EXPECT_FALSE(parsed.has_value()) << name;
            continue;
        }
        ASSERT_TRUE(parsed.has_value()) << name;
        const core::Array& e = list(expected);
        EXPECT_EQ(std::string(1, parsed->base), str(e[0])) << name;
        EXPECT_EQ(parsed->sugar, str(e[1])) << name;
        EXPECT_EQ(std::string(1, parsed->chirality), str(e[2])) << name;
        EXPECT_EQ(parsed->terminal, str(e[3])) << name;
    }
}

// NA-01 (MDATK 0.1.x pdbtools.has_backbone_atoms / standard_atom_name; data/sources/ccd)
TEST(Pdb, BackboneAtomsAndComponentAtomNames) {
    EXPECT_TRUE(io::hasBackboneAtoms(std::vector<std::string>{"N1", "C1'"}));
    EXPECT_TRUE(io::hasBackboneAtoms(std::vector<std::string>{"C1*"}));
    EXPECT_TRUE(io::hasBackboneAtoms(std::vector<std::string>{" P "}));
    EXPECT_FALSE(io::hasBackboneAtoms(std::vector<std::string>{"LA"}));
    EXPECT_FALSE(io::hasBackboneAtoms(std::vector<std::string>{"N1", "C2", "PA"}));
    EXPECT_EQ(io::standardAtomName("0DA", "N9A"), "N9");
    EXPECT_EQ(io::standardAtomName(" 0da ", "HN61"), "H61");
    EXPECT_EQ(io::standardAtomName("0DA", "H2P"), "HOP2");
    EXPECT_EQ(io::standardAtomName("0DT", "C5M"), "C7");
    EXPECT_EQ(io::standardAtomName("0DA", "N9"), "N9");
    EXPECT_EQ(io::standardAtomName("DA", "N9A"), "N9A");
    EXPECT_EQ(io::standardAtomName("0DC", "N1"), "N1");
}

TEST(Pdb, InsertHeaderRecords) {
    for (const core::Value& c : items(cases(), "insert")) {
        EXPECT_EQ(io::insertHeaderRecords(text(c, "text"), strings(*c.find("records"))),
                  text(c, "out"))
            << text(c, "text");
    }
}

TEST(Pdb, TransformText) {
    int n = 0;
    for (const core::Value& c : items(cases(), "transform")) {
        const std::string context = "transform case " + std::to_string(n++);
        std::array<std::array<double, 3>, 3> rotation{};
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                rotation[i][j] = hexFloat(list(items(c, "rotation")[i])[j]);
            }
        }
        std::array<double, 3> translation{};
        for (std::size_t i = 0; i < 3; ++i) {
            translation[i] = hexFloat(items(c, "translation")[i]);
        }
        const std::string input = text(c, "text");
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::transformPdbText(input, rotation, translation); }, *error,
                        context);
            continue;
        }
        const auto out = io::transformPdbText(input, rotation, translation);
        const core::Value& expected = *c.find("out");
        if (expected.isNull()) {
            EXPECT_FALSE(out.has_value()) << context;
        } else {
            ASSERT_TRUE(out.has_value()) << context;
            EXPECT_EQ(*out, str(expected)) << context;
        }
    }
}

TEST(Pdb, NewAtomTail) {
    for (const core::Value& c : items(cases(), "tail")) {
        const std::array<double, 3> xyz = {hexFloat(items(c, "xyz")[0]),
                                           hexFloat(items(c, "xyz")[1]),
                                           hexFloat(items(c, "xyz")[2])};
        const double occupancy = hexFloat(*c.find("occupancy"));
        const double bfactor = hexFloat(*c.find("bfactor"));
        const std::string context = str(items(c, "xyz")[0]) + " " + str(items(c, "xyz")[1]) + " " +
                                    str(items(c, "xyz")[2]) + " " + text(c, "occupancy");
        if (const core::Value* error = c.find("error")) {
            expectError([&] { io::newAtomTail(xyz, occupancy, bfactor); }, *error, context);
            continue;
        }
        EXPECT_EQ(io::newAtomTail(xyz, occupancy, bfactor), text(c, "out")) << context;
    }
}
