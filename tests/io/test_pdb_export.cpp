// The PDB export of trajectory frames (mdatk/io/pdb_export.py export_frames through MDAnalysis'
// PDBWriter) against the files MDATK writes for the cases of tests/parity/make_export_cases.py
// (data/io/pdb/export, see data/README.md), byte for byte.

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/io/loader.hpp"
#include "mdbridge/io/pdb_export.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const fs::path kExport = kData / "io" / "pdb" / "export";

const core::Array& list(const core::Value& value) {
    return *value.get<core::Array>();
}

std::string pathText(const fs::path& path) {
    const std::u8string text = path.u8string();
    return std::string(text.begin(), text.end());
}

std::string withPaths(std::string text) {
    const std::string data = pathText(kData);
    for (std::size_t at = text.find("{data}"); at != std::string::npos;
         at = text.find("{data}", at + data.size())) {
        text.replace(at, 6, data);
    }
    return text;
}

void expectFile(const fs::path& path, const core::Value& expected, const std::string& context) {
    const std::string actual = readFile(path);
    if (const core::Value* text = expected.find("text")) {
        EXPECT_EQ(actual, *text->get<std::string>()) << context;
        return;
    }
    EXPECT_EQ(static_cast<std::int64_t>(actual.size()), *expected.find("size")->get<std::int64_t>())
        << context;
    if (core::sha256Hex(actual) != text(expected, "sha256")) {
        const std::string& head = text(expected, "head");
        EXPECT_EQ(actual.substr(0, head.size()), head) << context;
        ADD_FAILURE() << context << ": SHA-256 differs";
    }
}

std::vector<std::string> fileNames(const fs::path& dir) {
    std::vector<std::string> names;
    if (fs::is_directory(dir)) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            names.push_back(pathText(entry.path().filename()));
        }
    }
    std::sort(names.begin(), names.end());
    return names;
}

// The 100001-atom system of make_export_cases.big_universe
core::System bigSystem() {
    const std::int32_t atoms = 100001;
    const std::int32_t residues = (atoms + 2) / 3;
    auto top = std::make_shared<core::Topology>();
    const char* names[] = {"C1", "N", "OXT"};
    const char* elements[] = {"C", "N", "O"};
    const char* resnames[] = {"ALA", "GLY", "HOH", "MG"};
    std::vector<std::uint32_t> elementIds;
    for (std::int32_t k = 0; k < atoms; ++k) {
        top->names.push_back(top->strings.intern(names[k % 3]));
        top->types.push_back(0);
        top->atomResidue.push_back(k / 3);
        elementIds.push_back(top->strings.intern(elements[k % 3]));
    }
    top->elements = std::move(elementIds);
    for (std::int32_t r = 0; r < residues; ++r) {
        top->resids.push_back(static_cast<std::int64_t>(r) * 7 - 20000);
        top->resnames.push_back(top->strings.intern(resnames[r % 4]));
        top->residueSegment.push_back(r >= residues / 2 ? 1 : 0);
    }
    top->resnums = top->resids;
    top->segids = {top->strings.intern("PROT"), top->strings.intern("SOLVENT")};
    std::vector<core::Bond> bonds;
    for (std::int32_t k = 0; k < 60; k += 3) {
        bonds.push_back({k, k + 1});
    }
    bonds.push_back({atoms - 2, atoms - 1});
    top->bonds = std::move(bonds);
    top->finish();
    core::Frame frame;
    for (std::int32_t k = 0; k < atoms; ++k) {
        const double x = static_cast<double>(k % 1000) * 0.5 - 250.0;
        const double y = static_cast<double>(k / 1000) * 0.25 - 12.5;
        const double z =
            static_cast<double>((static_cast<std::int64_t>(k) * 37) % 2001 - 1000) * 0.01;
        frame.positions.push_back(static_cast<float>(x));
        frame.positions.push_back(static_cast<float>(y));
        frame.positions.push_back(static_cast<float>(z));
    }
    frame.box = std::array<float, 6>{500.0f, 260.0f, 30.0f, 90.0f, 90.0f, 90.0f};
    core::System system;
    system.topology = std::move(top);
    system.trajectory = std::make_unique<core::MemoryTrajectory>(std::vector<core::Frame>{frame});
    return system;
}

}  // namespace

TEST(PdbExport, MdatkCases) {
    const core::Value data = readJson(kExport / "export_cases.json");
    int n = 0;
    for (const core::Value& c : list(*data.find("cases"))) {
        const std::string context =
            "case " + std::to_string(n++) + " " + text(c, "topology") + " " + text(c, "prefix");
        std::vector<fs::path> trajectories;
        for (const core::Value& t : list(*c.find("trajectories"))) {
            trajectories.push_back(withPaths(*t.get<std::string>()));
        }
        std::optional<double> dt;
        if (const core::Value* v = c.find("dt"); !v->isNull()) {
            dt =
                v->get<double>() ? *v->get<double>() : static_cast<double>(*v->get<std::int64_t>());
        }
        core::System system;
        try {
            system = io::loadSystem(withPaths(text(c, "topology")), trajectories, dt);
        } catch (const std::exception& e) {
            ADD_FAILURE() << context << ": " << e.what();
            continue;
        }
        const std::string selection = text(c, "selection");
        const core::AtomGroup atoms =
            selection.empty() ? system.atoms() : core::select(system, selection);
        std::vector<std::int64_t> frames;
        for (const core::Value& f : list(*c.find("frames"))) {
            frames.push_back(*f.get<std::int64_t>());
        }
        TempDir dir;
        const fs::path out = dir / "out";
        std::vector<std::pair<double, std::string>> calls;
        const auto progress = [&](double fraction, const std::string& message) {
            calls.emplace_back(fraction, message);
        };
        try {
            const auto paths = io::exportFrames(system, frames, atoms, out, text(c, "prefix"),
                                                progress, text(c, "stamp"));
            ASSERT_FALSE(c.find("error")) << context << ": no error";
            const core::Array& expected = list(*c.find("paths"));
            ASSERT_EQ(paths.size(), expected.size()) << context;
            for (std::size_t k = 0; k < paths.size(); ++k) {
                EXPECT_EQ(paths[k], out / *expected[k].get<std::string>()) << context;
            }
        } catch (const core::ValueError& e) {
            const core::Value* error = c.find("error");
            ASSERT_TRUE(error) << context << ": " << e.what();
            EXPECT_EQ(std::string(e.pythonName()), text(*error, "type")) << context;
            EXPECT_EQ(e.what(), text(*error, "message")) << context;
        }
        const core::Array& expectedCalls = list(*c.find("progress"));
        ASSERT_EQ(calls.size(), expectedCalls.size()) << context;
        for (std::size_t k = 0; k < calls.size(); ++k) {
            EXPECT_EQ(calls[k].first,
                      std::strtod(list(expectedCalls[k])[0].get<std::string>()->c_str(), nullptr))
                << context;
            EXPECT_EQ(calls[k].second, *list(expectedCalls[k])[1].get<std::string>()) << context;
        }
        const core::Array& files = list(*c.find("files"));
        std::vector<std::string> expectedNames;
        for (const core::Value& f : files) {
            expectedNames.push_back(text(f, "name"));
        }
        EXPECT_EQ(fileNames(out), expectedNames) << context;
        for (const core::Value& f : files) {
            if (fs::exists(out / text(f, "name"))) {
                expectFile(out / text(f, "name"), f, context + " " + text(f, "name"));
            }
        }
    }
    EXPECT_GT(n, 60);
}

TEST(PdbExport, SerialAndResidueWraps) {
    core::System system = bigSystem();
    const core::Value data = readJson(kExport / "export_cases.json");
    for (const core::Value& c : list(*data.find("big"))) {
        const std::string name = text(c, "selection");
        std::vector<std::int32_t> indices;
        if (name == "all") {
            indices.resize(system.topology->atomCount());
        } else if (name == "first") {
            indices.resize(60);
        }
        for (std::size_t k = 0; k < indices.size(); ++k) {
            indices[k] = static_cast<std::int32_t>(k);
        }
        if (name == "tail") {
            for (std::int32_t k = 99990; k < 100001; ++k) {
                indices.push_back(k);
            }
        }
        const core::AtomGroup atoms(system.topology, indices);
        const std::string text =
            io::pdbFrameText(atoms, system.trajectory->current(), system.trajectory->info());
        EXPECT_EQ(static_cast<std::int64_t>(text.size()), *c.find("size")->get<std::int64_t>())
            << name;
        if (core::sha256Hex(text) != mdbridge::test::text(c, "sha256")) {
            const std::string& head = mdbridge::test::text(c, "head");
            EXPECT_EQ(text.substr(0, head.size()), head) << name;
            const std::string& tail = mdbridge::test::text(c, "tail");
            EXPECT_TRUE(text.ends_with(tail + "\n")) << name << "\n"
                                                     << text.substr(text.size() - 1200);
            ADD_FAILURE() << name << ": SHA-256 differs";
        }
    }
}

TEST(PdbExport, EmptyFrameList) {
    core::System system = bigSystem();
    TempDir dir;
    EXPECT_THROW(io::exportFrames(system, {}, system.atoms(), dir / "x", "p"),
                 core::SelectionError);
    EXPECT_FALSE(fs::exists(dir / "x"));
}
