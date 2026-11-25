// Amber NetCDF: the reader on its own and the writer, whose files must equal MDAnalysis' NCDFWriter
// output (scipy backend) byte for byte (data/io/trajectories/writer, see data/README.md). The
// reader's MDAnalysis parity cases are in test_trajectory.cpp.

#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/io/netcdf.hpp"
#include "mdbridge/io/trajectory.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const fs::path kTrajectories = kData / "io" / "trajectories";

const core::Array& list(const core::Value& value) {
    return *value.get<core::Array>();
}

std::string utf8Path(const fs::path& path) {
    const std::u8string text = path.u8string();
    return std::string(text.begin(), text.end());
}

}  // namespace

TEST(Netcdf, WriterMatchesMDAnalysis) {
    const core::Value data = readJson(kTrajectories / "cases.json");
    TempDir dir;
    int written = 0;
    for (const core::Value& w : list(*data.find("writer"))) {
        const std::string name = text(w, "name");
        const auto atoms = static_cast<std::size_t>(*w.find("atoms")->get<std::int64_t>());
        auto source =
            io::openTrajectory(kTrajectories / text(w, "source"), "", std::nullopt, atoms);
        std::vector<core::Frame> frames;
        for (const core::Value& k : list(*w.find("frames"))) {
            frames.push_back(source->seek(*k.get<std::int64_t>()));
        }
        const fs::path out = dir / (name + ".nc");
        io::writeNetcdfTrajectory(out, frames, atoms);
        EXPECT_TRUE(readFile(out) == readFile(kTrajectories / text(w, "expected"))) << name;
        ++written;
    }
    EXPECT_EQ(written, 6);
}

TEST(Netcdf, DemoTrajectory) {
    // MDATK's dry demo (773 atoms, 24 frames of 10 ps, a box): the file's own atom count is taken
    auto demo = io::openNetcdfTrajectory(kTrajectories / "demo_dry.nc");
    EXPECT_EQ(demo->frameCount(), 24);
    EXPECT_EQ(demo->dt(), 10.0);
    EXPECT_EQ(demo->current().positions.size(), 773u * 3);
    EXPECT_EQ(demo->seek(23).time, 230.0);
    EXPECT_TRUE(demo->current().box.has_value());
    EXPECT_FALSE(demo->info().header.has_value());
    // the given dt changes dt() but not the file's times
    auto given = io::openNetcdfTrajectory(kTrajectories / "demo_dry.nc", 2.5, 773);
    EXPECT_EQ(given->dt(), 2.5);
    EXPECT_EQ(given->seek(3).time, 30.0);
    try {
        io::openNetcdfTrajectory(kTrajectories / "demo_dry.nc", std::nullopt, 10);
        ADD_FAILURE();
    } catch (const core::ValueError& e) {
        // RD-N02: the message of the other readers
        EXPECT_EQ(std::string(e.what()),
                  "The topology and NCDF trajectory files don't have the same number of atoms!\n"
                  "Topology number of atoms 10\nTrajectory: " +
                      utf8Path(kTrajectories / "demo_dry.nc") + " Number of atoms 773");
    }
}

TEST(Netcdf, WriterRoundTrip) {
    TempDir dir;
    std::vector<core::Frame> frames(3);
    for (std::size_t k = 0; k < frames.size(); ++k) {
        frames[k].time = 0.1 * static_cast<double>(k);
        frames[k].positions = {1.5f, -2.25f, 3.0f, static_cast<float>(k), 0.1f, 1e6f};
        frames[k].box = std::array<float, 6>{10.0f, 11.0f, 12.0f, 90.0f, 100.0f, 120.0f};
    }
    const fs::path path = dir / "round.nc";
    io::writeNetcdfTrajectory(path, frames, 2);
    auto t = io::openNetcdfTrajectory(path);
    ASSERT_EQ(t->frameCount(), 3);
    EXPECT_EQ(t->dt(), static_cast<double>(static_cast<float>(0.1)));
    for (std::int64_t k = 0; k < 3; ++k) {
        const core::Frame& f = t->seek(k);
        EXPECT_EQ(f.positions, frames[static_cast<std::size_t>(k)].positions);
        EXPECT_EQ(f.time, static_cast<double>(
                              static_cast<float>(frames[static_cast<std::size_t>(k)].time)));
        EXPECT_EQ(f.box, frames[static_cast<std::size_t>(k)].box);
    }
    // no frames: no file (NCDFWriter opens its file with the first frame)
    io::writeNetcdfTrajectory(dir / "empty.nc", {}, 2);
    EXPECT_FALSE(fs::exists(dir / "empty.nc"));
    EXPECT_THROW(io::writeNetcdfTrajectory(dir / "bad.nc", frames, 3), core::ValueError);
    frames[1].box.reset();
    EXPECT_THROW(io::writeNetcdfTrajectory(dir / "bad.nc", frames, 2), core::ValueError);
}
