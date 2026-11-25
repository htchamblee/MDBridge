// Coordinate files read as MDAnalysis 2.11.0 reads them for MDATK: every case of
// tests/parity/make_trajectory_cases.py (data/io/trajectories, see data/README.md) - frame counts,
// dt, times, boxes and float32 positions bit for bit, chains, and the errors.

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/io/trajectory.hpp"
#include "mdbridge/io/warnings.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const fs::path kTrajectories = kData / "io" / "trajectories";

const core::Array& list(const core::Value& value) {
    return *value.get<core::Array>();
}

double hexDouble(const std::string& text) {
    return std::strtod(text.c_str(), nullptr);
}

std::uint32_t bitsOf(float x) {
    std::uint32_t bits;
    std::memcpy(&bits, &x, 4);
    return bits;
}

std::string hexBits(const std::vector<float>& values) {
    static const char digits[] = "0123456789abcdef";
    std::string out;
    for (const float x : values) {
        const std::uint32_t bits = bitsOf(x);
        for (int shift = 28; shift >= 0; shift -= 4) {
            out += digits[bits >> shift & 15];
        }
    }
    return out;
}

std::string littleEndianBytes(const std::vector<float>& values) {
    std::string out;
    for (const float x : values) {
        const std::uint32_t bits = bitsOf(x);
        for (int shift = 0; shift < 32; shift += 8) {
            out += static_cast<char>(bits >> shift & 0xFF);
        }
    }
    return out;
}

std::string pathText(const fs::path& path) {
    const std::u8string text = path.u8string();
    return std::string(text.begin(), text.end());
}

// A path as Python's repr() shows it inside quotes: a '\\' separator (Windows) is doubled.
std::string reprText(const fs::path& path) {
    std::string out;
    for (const char c : pathText(path)) {
        if (c == '\\') {
            out += c;
        }
        out += c;
    }
    return out;
}

// The paths in a message replaced by "{path<k>}": MDAnalysis names a file as given or, in a list
// of files, by its repr(); the two differ when the path has a '\\' (Windows).
std::string withPlaceholders(std::string message, const std::vector<fs::path>& paths) {
    std::vector<std::size_t> order(paths.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    for (std::string (*const form)(const fs::path&) : {reprText, pathText}) {
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return form(paths[a]).size() > form(paths[b]).size();
        });
        for (const std::size_t k : order) {
            const std::string from = form(paths[k]);
            const std::string to = "{path" + std::to_string(k) + "}";
            for (std::size_t at = message.find(from); at != std::string::npos;
                 at = message.find(from, at + to.size())) {
                message.replace(at, from.size(), to);
            }
        }
    }
    return message;
}

void expectError(const core::ValueError& error, const core::Value& expected,
                 const std::vector<fs::path>& paths, const std::string& context) {
    EXPECT_EQ(error.pythonName(), text(expected, "type")) << context;
    EXPECT_EQ(withPlaceholders(error.what(), paths), text(expected, "message")) << context;
}

// Cases whose MDAnalysis error comes from deep inside scipy / numpy or from a reader MDBridge does
// not have: MDBridge rejects the file with its own message.
const std::set<std::string> kOwnErrors = {
    "crd_guessed",
    "chain_crd_guessed",
};

}  // namespace

TEST(Trajectory, MDAnalysisCases) {
    const core::Value data = readJson(kTrajectories / "cases.json");
    int checked = 0;
    for (const core::Value& c : list(*data.find("cases"))) {
        const std::string name = text(c, "name");
        std::vector<fs::path> paths;
        for (const core::Value& f : list(*c.find("files"))) {
            paths.push_back(kTrajectories / *f.get<std::string>());
        }
        std::optional<double> dt;
        if (!c.find("dt")->isNull()) {
            dt = hexDouble(text(c, "dt"));
        }
        const auto atoms = static_cast<std::size_t>(*c.find("atoms")->get<std::int64_t>());
        std::unique_ptr<core::Trajectory> trajectory;
        io::ReaderWarningScope warnings;  // LD-03: the case's "warnings" (none when absent)
        try {
            trajectory = io::openTrajectories(paths, text(c, "format"), dt, atoms);
        } catch (const core::ValueError& error) {
            const core::Value* expected = c.find("error");
            if (expected == nullptr) {
                ADD_FAILURE() << name << ": " << error.what();
            } else if (kOwnErrors.count(name) == 0) {
                expectError(error, *expected, paths, name);
            }
            ++checked;
            continue;
        }
        const core::Value* result = c.find("result");
        if (result == nullptr) {
            ADD_FAILURE() << name << ": no error, expected " << text(*c.find("error"), "message");
            continue;
        }
        ++checked;
        ASSERT_EQ(trajectory->frameCount(), *result->find("n_frames")->get<std::int64_t>()) << name;
        const double dt0 = hexDouble(text(*result, "dt"));
        EXPECT_EQ(trajectory->dt(), dt0) << name;
        if (const core::Value* header = result->find("header")) {
            EXPECT_EQ(trajectory->info().header.value_or("<none>"), *header->get<std::string>())
                << name;
        } else {
            EXPECT_FALSE(trajectory->info().header.has_value()) << name;
        }
        const core::Array& frames = list(*result->find("frames"));
        for (std::size_t k = 0; k < frames.size(); ++k) {
            const core::Value& expected = frames[k];
            const std::string context = name + " frame " + std::to_string(k);
            if (const core::Value* error = expected.find("error")) {
                try {
                    trajectory->seek(static_cast<std::int64_t>(k));
                    ADD_FAILURE() << context << ": no error";
                } catch (const core::ValueError& e) {
                    expectError(e, *error, paths, context);
                }
                continue;
            }
            const core::Frame& frame = trajectory->seek(static_cast<std::int64_t>(k));
            EXPECT_EQ(frame.index, static_cast<std::int64_t>(k)) << context;
            EXPECT_EQ(frame.time, hexDouble(text(expected, "time"))) << context;
            const core::Value* frameDt = expected.find("dt");
            EXPECT_EQ(trajectory->dt(), frameDt ? hexDouble(*frameDt->get<std::string>()) : dt0)
                << context;
            const core::Value& box = *expected.find("box");
            if (box.isNull()) {
                EXPECT_FALSE(frame.box.has_value()) << context;
            } else if (!frame.box) {
                ADD_FAILURE() << context << ": no box";
            } else {
                std::string bits;
                for (const core::Value& v : list(box)) {
                    bits += *v.get<std::string>();
                }
                EXPECT_EQ(hexBits(std::vector<float>(frame.box->begin(), frame.box->end())), bits)
                    << context;
            }
            if (const core::Value* pos = expected.find("pos")) {
                EXPECT_EQ(hexBits(frame.positions), *pos->get<std::string>()) << context;
            } else {
                EXPECT_EQ(core::sha256Hex(littleEndianBytes(frame.positions)),
                          text(expected, "pos_sha256"))
                    << context;
            }
        }
        std::vector<std::string> expectedWarnings;
        if (const core::Value* w = c.find("warnings")) {
            for (const core::Value& v : list(*w)) {
                expectedWarnings.push_back(*v.get<std::string>());
            }
        }
        std::vector<std::string> issued;
        for (const std::string& w : warnings.take()) {
            issued.push_back(withPlaceholders(w, paths));
        }
        EXPECT_EQ(issued, expectedWarnings) << name;
    }
    EXPECT_GT(checked, 200);
}

TEST(Trajectory, FormatForExtension) {
    const core::Value data = readJson(kTrajectories / "formats.json");
    for (const auto& [ext, format] : *data.find("mdatk_table")->get<core::Table>()) {
        EXPECT_EQ(io::coordinateFormatForExtension("dir/x" + ext), *format.get<std::string>())
            << ext;
        std::string upper = ext;
        for (char& ch : upper) {
            ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        }
        EXPECT_EQ(io::coordinateFormatForExtension("x" + upper), *format.get<std::string>())
            << upper;
    }
    for (const core::Value& g : list(*data.find("guess"))) {
        const std::string path = text(g, "path");
        if (const core::Value* guess = g.find("guess")) {
            const std::string ext =
                path.substr(path.rfind('.') == std::string::npos ? path.size() : path.rfind('.'));
            std::string lower = ext;
            for (char& ch : lower) {
                ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
            }
            const bool inTable = data.find("mdatk_table")->find(lower) != nullptr;
            if (!inTable) {
                EXPECT_EQ(io::coordinateFormatForExtension(path), *guess->get<std::string>())
                    << path;
            }
        }
    }
}

TEST(Trajectory, GuessedFormats) {
    // per-file guesses: the reader MDAnalysis picks, or its "Unknown coordinate trajectory format"
    const core::Value data = readJson(kTrajectories / "formats.json");
    TempDir dir;
    for (const core::Value& g : list(*data.find("guess"))) {
        const core::Value* error = g.find("load_error");
        if (error == nullptr || text(*error, "type") != "TypeError") {
            continue;
        }
        const fs::path path = dir / text(g, "path");
        fs::create_directories(path.parent_path());
        writeFile(path, "");
        try {
            io::openTrajectory(path, "", std::nullopt, 7);
            ADD_FAILURE() << text(g, "path");
        } catch (const core::ValueError& e) {
            EXPECT_STREQ(e.pythonName(), "TypeError");
            EXPECT_EQ(withPlaceholders(e.what(), {path}), text(*error, "message"))
                << text(g, "path");
        }
    }
    for (const core::Value& f : list(*data.find("format_errors"))) {
        const core::Value* error = f.find("error");
        if (error == nullptr || text(*error, "type") != "TypeError") {
            continue;
        }
        const fs::path path = kTrajectories / "files" / "mda_box.nc";
        try {
            io::openTrajectory(path, text(f, "format"), std::nullopt, 7);
            ADD_FAILURE() << text(f, "format");
        } catch (const core::ValueError& e) {
            EXPECT_EQ(e.pythonName(), std::string("TypeError"));
            EXPECT_EQ(withPlaceholders(e.what(), {path}), text(*error, "message"))
                << text(f, "format");
        }
    }
    // formats MDAnalysis reads but MDBridge does not: a clear error naming the format
    for (const char* format : {"CRD", "XYZ", "TNG", "crd"}) {
        try {
            io::openTrajectory(kTrajectories / "files" / "box.crd", format, std::nullopt, 7);
            ADD_FAILURE() << format;
        } catch (const core::ValueError& e) {
            EXPECT_NE(std::string(e.what()).find("MDBridge cannot read"), std::string::npos)
                << e.what();
        }
    }
    EXPECT_THROW(io::openTrajectory(kTrajectories / "files" / "box.crd", "", std::nullopt, 7),
                 core::ValueError);
}

TEST(Trajectory, ChainOfMemoryTrajectories) {
    // parts without stored times (PDB / GRO / memory) count time from the global frame
    std::vector<core::Frame> a(2), b(3);
    for (auto* frames : {&a, &b}) {
        for (auto& f : *frames) {
            f.positions = {1.0f, 2.0f, 3.0f};
        }
    }
    std::vector<std::unique_ptr<core::Trajectory>> parts;
    parts.push_back(std::make_unique<core::MemoryTrajectory>(a, 2.0));
    parts.push_back(std::make_unique<core::MemoryTrajectory>(b, 0.5));
    auto chain = io::chainTrajectories(std::move(parts));
    EXPECT_EQ(chain->frameCount(), 5);
    EXPECT_EQ(chain->dt(), 2.0);
    EXPECT_EQ(chain->seek(1).time, 2.0);
    EXPECT_EQ(chain->seek(-1).time, 2.0);  // 4 * 0.5
    EXPECT_EQ(chain->dt(), 0.5);
    EXPECT_EQ(chain->current().index, 4);
    EXPECT_FALSE(chain->info().header.has_value());
    EXPECT_THROW(chain->seek(5), std::out_of_range);
    EXPECT_THROW(chain->seek(-6), std::out_of_range);
    EXPECT_THROW(io::chainTrajectories({}), core::ValueError);
}

TEST(Trajectory, RandomAccess) {
    // frames are decoded on demand, in any order
    auto t = io::openTrajectory(kTrajectories / "files" / "water.xtc", "XTC", std::nullopt, 90);
    const std::vector<float> last = t->seek(4).positions;
    const std::vector<float> first = t->seek(0).positions;
    EXPECT_EQ(t->seek(-1).positions, last);
    EXPECT_EQ(t->seek(-5).positions, first);
    EXPECT_THROW(t->seek(5), std::out_of_range);
    EXPECT_THROW(
        io::openTrajectory(kTrajectories / "files" / "missing.nc", "NCDF", std::nullopt, 0),
        fs::filesystem_error);
}
