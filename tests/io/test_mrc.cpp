// Density maps (mdatk/io/mrc.py and the map functions of mdatk/density.py) against what MDATK
// gives for the cases of tests/parity/make_map_cases.py (data/io/maps, see data/README.md), and
// ports of MDATK's own I/O, simulation and info tests (tests/test_density.py).

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/mrc.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;
using core::Array;
using core::Value;
using core::Vec3;
using io::DensityMap;

namespace {

const fs::path kMaps = kData / "io" / "maps";

// simulateMap rounds as the C library's exp(); the committed hashes are glibc's (which is what
// numpy's scalar exp is on Linux)
#ifdef __GLIBC__
constexpr bool kGlibcExp = true;
#else
constexpr bool kGlibcExp = false;
#endif

const Value& load(const char* name) {
    static std::map<std::string, Value, std::less<>> files;
    auto it = files.find(name);
    if (it == files.end()) {
        it = files.emplace(name, core::parseJson(readFile(kMaps / name), name)).first;
    }
    return it->second;
}

const Array& list(const Value& value) {
    return *value.get<Array>();
}
const Array& items(const Value& value, std::string_view key) {
    return list(*value.find(key));
}
const std::string& str(const Value& value) {
    return *value.get<std::string>();
}
std::int64_t integer(const Value& value) {
    return *value.get<std::int64_t>();
}

double number(const Value& value) {
    if (const double* d = value.get<double>()) {
        return *d;
    }
    return static_cast<double>(integer(value));
}

Vec3 vec3(const Value& value) {
    const Array& a = list(value);
    return {number(a.at(0)), number(a.at(1)), number(a.at(2))};
}

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

void expectVec(const Vec3& got, const Value& expected, const std::string& context) {
    const Vec3 e = vec3(expected);
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_TRUE(sameDouble(got[k], e[k]))
            << context << " [" << k << "] " << got[k] << " vs " << e[k];
    }
}

std::uint32_t bitsOf(double value) {
    const auto f = static_cast<float>(value);
    std::uint32_t bits;
    std::memcpy(&bits, &f, 4);
    return bits;
}

void expectStats(const io::MapStats& st, const Value& expected, const std::string& context) {
    const Array& e = list(expected);
    const double got[4] = {st.min, st.max, st.mean, st.rms};
    const char* names[4] = {"min", "max", "mean", "rms"};
    for (std::size_t k = 0; k < 4; ++k) {
        const auto want = static_cast<std::uint32_t>(integer(e.at(k)));
        float wantFloat;
        std::memcpy(&wantFloat, &want, 4);
        if (std::isnan(wantFloat)) {
            EXPECT_TRUE(std::isnan(got[k])) << context << " " << names[k];
        } else {
            EXPECT_EQ(bitsOf(got[k]), want)
                << context << " " << names[k] << ": " << got[k] << " vs " << wantFloat;
        }
    }
}

std::string dataSha(const std::vector<float>& data) {
    std::string bytes;
    bytes.reserve(data.size() * 4);
    for (const float x : data) {
        std::uint32_t bits;
        std::memcpy(&bits, &x, 4);
        for (int shift = 0; shift < 32; shift += 8) {
            bytes += static_cast<char>(bits >> shift & 0xFF);
        }
    }
    return core::sha256Hex(bytes);
}

std::string hex(std::string_view bytes) {
    static const char* digits = "0123456789abcdef";
    std::string out;
    for (const char c : bytes) {
        const auto b = static_cast<unsigned char>(c);
        out += digits[b >> 4];
        out += digits[b & 15];
    }
    return out;
}

std::string replaceAll(std::string text, std::string_view from, std::string_view to) {
    for (std::size_t pos = text.find(from); pos != std::string::npos;
         pos = text.find(from, pos + to.size())) {
        text.replace(pos, from.size(), to);
    }
    return text;
}

// ---- the generated arrays of make_map_cases.py ------------------------------------------------

std::uint64_t splitmix(std::uint64_t seed, std::uint64_t i) {
    std::uint64_t z = seed + i * 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

std::vector<float> generate(const Value& spec) {
    const auto n = static_cast<std::size_t>(integer(*spec.find("n")));
    const std::string kind = str(*spec.find("kind"));
    std::vector<float> out(n);
    if (kind == "uniform") {
        const auto seed = static_cast<std::uint64_t>(integer(*spec.find("seed")));
        const double scale = number(*spec.find("scale")), offset = number(*spec.find("offset"));
        for (std::size_t i = 0; i < n; ++i) {
            const double u = static_cast<double>(splitmix(seed, i + 1) >> 40) / 16777216.0;
            out[i] = static_cast<float>((u - 0.5) * scale + offset);
        }
    } else if (kind == "bits") {
        const auto seed = static_cast<std::uint64_t>(integer(*spec.find("seed")));
        const auto emin = static_cast<std::uint64_t>(integer(*spec.find("emin")));
        const auto span = static_cast<std::uint64_t>(integer(*spec.find("emax"))) - emin + 1;
        for (std::size_t i = 0; i < n; ++i) {
            const std::uint64_t z = splitmix(seed, i + 1);
            const auto bits =
                static_cast<std::uint32_t>((z & 0x807FFFFFu) | (emin + (z >> 32) % span) << 23);
            std::memcpy(&out[i], &bits, 4);
        }
    } else if (kind == "constant") {
        std::fill(out.begin(), out.end(), static_cast<float>(number(*spec.find("value"))));
    } else if (kind == "ramp") {
        const double start = number(*spec.find("start")), step = number(*spec.find("step"));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = static_cast<float>(start + static_cast<double>(i) * step);
        }
    } else {
        ADD_FAILURE() << "unknown data kind " << kind;
    }
    return out;
}

std::string mapError(const fs::path& path) {
    try {
        io::readMap(path);
    } catch (const io::MapError& e) {
        return e.what();
    }
    return "";
}

// ---- writing map files by hand (MDATK's tests/test_density.py::_raw_map) ----------------------

void putInt(std::string& raw, std::size_t offset, std::int32_t v, bool big) {
    auto u = static_cast<std::uint32_t>(v);
    for (int k = 0; k < 4; ++k) {
        const int shift = big ? 24 - 8 * k : 8 * k;
        raw[offset + static_cast<std::size_t>(k)] = static_cast<char>(u >> shift & 0xFF);
    }
}

void putFloat(std::string& raw, std::size_t offset, float v, bool big) {
    std::int32_t bits;
    std::memcpy(&bits, &v, 4);
    putInt(raw, offset, bits, big);
}

struct RawMap {
    std::array<int, 3> shapeSrc{};  // ns, nr, nc
    std::vector<double> dataSrc;    // [section][row][column]
    std::array<int, 3> mapcrs{1, 2, 3};
    Vec3 voxel{1, 1, 1};
    Vec3 origin{0, 0, 0};
    std::array<int, 3> nstart{0, 0, 0};
    bool big = false;
    int mode = 2;
    int nsymbt = 0;

    std::string bytes() const {
        const int ns = shapeSrc[0], nr = shapeSrc[1], nc = shapeSrc[2];
        std::string raw(1024, '\0');
        const std::array<int, 3> dims{nc, nr, ns};
        putInt(raw, 0, nc, big);
        putInt(raw, 4, nr, big);
        putInt(raw, 8, ns, big);
        putInt(raw, 12, mode, big);
        std::array<int, 3> nCart{};
        for (std::size_t k = 0; k < 3; ++k) {
            putInt(raw, 16 + 4 * k, nstart[k], big);
            nCart[static_cast<std::size_t>(mapcrs[k] - 1)] = dims[k];
        }
        for (std::size_t k = 0; k < 3; ++k) {
            putInt(raw, 28 + 4 * k, nCart[k], big);
            putFloat(raw, 40 + 4 * k, static_cast<float>(nCart[k] * voxel[k]), big);
            putFloat(raw, 52 + 4 * k, 90.0f, big);
            putInt(raw, 64 + 4 * k, mapcrs[k], big);
            putFloat(raw, 196 + 4 * k, static_cast<float>(origin[k]), big);
        }
        putInt(raw, 92, nsymbt, big);
        putInt(raw, 108, 20140, big);
        raw.replace(208, 4, "MAP ");
        raw[212] = raw[213] = big ? '\x11' : '\x44';
        raw += std::string(static_cast<std::size_t>(nsymbt), '\0');
        for (const double v : dataSrc) {
            std::string item(4, '\0');
            switch (mode) {
            case 0:
                raw += static_cast<char>(static_cast<std::int8_t>(v));
                continue;
            case 1:
            case 6: {
                const auto u = mode == 1 ? static_cast<std::uint16_t>(static_cast<std::int16_t>(v))
                                         : static_cast<std::uint16_t>(v);
                raw += static_cast<char>(big ? u >> 8 : u & 0xFF);
                raw += static_cast<char>(big ? u & 0xFF : u >> 8);
                continue;
            }
            default:
                putFloat(item, 0, static_cast<float>(v), big);
                raw += item;
            }
        }
        return raw;
    }
};

// exp(-|p - c|^2 / (2 w^2)) on a [z][y][x] grid, as float32 (MDATK's _blob)
std::vector<float> blob(std::array<int, 3> shapeZyx, Vec3 centreXyz, double width = 1.5) {
    std::vector<float> out;
    for (int z = 0; z < shapeZyx[0]; ++z) {
        for (int y = 0; y < shapeZyx[1]; ++y) {
            for (int x = 0; x < shapeZyx[2]; ++x) {
                const double d = (x - centreXyz[0]) * (x - centreXyz[0]) +
                                 (y - centreXyz[1]) * (y - centreXyz[1]) +
                                 (z - centreXyz[2]) * (z - centreXyz[2]);
                out.push_back(static_cast<float>(std::exp(-d / (2 * width * width))));
            }
        }
    }
    return out;
}

std::uint32_t crc32(std::string_view data) {
    std::uint32_t crc = 0xFFFFFFFFu;
    for (const char c : data) {
        crc ^= static_cast<unsigned char>(c);
        for (int k = 0; k < 8; ++k) {
            crc = crc & 1 ? (crc >> 1) ^ 0xEDB88320u : crc >> 1;
        }
    }
    return ~crc;
}

// a gzip member of stored deflate blocks
std::string gzipStored(std::string_view data, std::size_t block = 65535) {
    std::string out("\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03", 10);
    std::size_t pos = 0;
    do {
        const std::size_t n = std::min(block, data.size() - pos);
        const bool last = pos + n == data.size();
        out += static_cast<char>(last ? 1 : 0);
        out += static_cast<char>(n & 0xFF);
        out += static_cast<char>(n >> 8);
        out += static_cast<char>(~n & 0xFF);
        out += static_cast<char>((~n >> 8) & 0xFF);
        out.append(data.substr(pos, n));
        pos += n;
    } while (pos < data.size());
    std::string tail(8, '\0');
    putInt(tail, 0, static_cast<std::int32_t>(crc32(data)), false);
    putInt(tail, 4, static_cast<std::int32_t>(data.size()), false);
    return out + tail;
}

}  // namespace

// ---- MDATK's results -------------------------------------------------------------------------

TEST(Mrc, ReadCasesMatchMdatk) {
    int checked = 0;
    for (const Value& c : items(load("read_cases.json"), "cases")) {
        const std::string name = str(*c.find("name"));
        const fs::path path = kMaps / "files" / name;
        const std::string given = core::pathToUtf8(path);
        SCOPED_TRACE(name);
        if (const Value* error = c.find("error")) {
            EXPECT_EQ(replaceAll(mapError(path), given, "{path}"), str(*error));
            ++checked;
            continue;
        }
        const DensityMap map = io::readMap(path);
        const Array& shape = items(c, "shape_xyz");
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_EQ(map.shapeXyz[k], integer(shape.at(k)));
        }
        expectVec(map.origin, *c.find("origin"), "origin");
        expectVec(map.voxel, *c.find("voxel"), "voxel");
        EXPECT_EQ(dataSha(map.data), str(*c.find("sha256")));
        expectStats(map.stats(), *c.find("stats"), "stats");
        const Array& first = items(c, "first");
        for (std::size_t k = 0; k < first.size(); ++k) {
            EXPECT_TRUE(sameDouble(map.data.at(k), number(first[k])));
        }
        EXPECT_EQ(map.path, given);
        ++checked;
    }
    EXPECT_GT(checked, 100);
}

TEST(Mrc, DensityInfoMatchesMdatk) {
    for (const Value& c : items(load("info_cases.json"), "cases")) {
        const fs::path path = kMaps / "files" / str(*c.find("name"));
        Value info = io::densityInfo(path);
        Value* map = info.find("map");
        ASSERT_NE(map, nullptr);
        EXPECT_TRUE(core::posixpath::isabs(str(*map)) ||
                    fs::path(core::pathFromUtf8(str(*map))).is_absolute());
        EXPECT_TRUE(fs::equivalent(core::pathFromUtf8(str(*map)), path));
        *map = Value("{path}");
        EXPECT_EQ(core::toJson(info, 2) + "\n", str(*c.find("output"))) << str(*c.find("name"));
    }
}

TEST(Mrc, StatsMatchNumpy) {
    for (const Value& c : items(load("stats_cases.json"), "cases")) {
        const Value& spec = *c.find("data");
        const std::string context = core::toJson(spec, -1);
        DensityMap map;
        map.data = generate(spec);
        map.shapeXyz = {static_cast<std::int64_t>(map.data.size()), 1, 1};
        map.voxel = {1, 1, 1};
        const io::MapStats st = map.stats();
        expectStats(st, *c.find("stats"), context);
        for (const Value& level : items(c, "levels")) {
            const double sigmas = number(list(level).at(0));
            // mapLevel recomputes the statistics: only on the smaller arrays
            const double got =
                map.data.size() <= 100000 ? io::mapLevel(map, sigmas) : st.mean + sigmas * st.rms;
            EXPECT_TRUE(sameDouble(got, number(list(level).at(1)))) << context;
        }
    }
}

TEST(Mrc, WriteMatchesMdatk) {
    TempDir dir;
    const Value& data = load("write_cases.json");
    for (const Value& c : items(data, "cases")) {
        const std::string context = core::toJson(c.find("data") ? *c.find("data") : Value(), -1);
        DensityMap map;
        map.data = generate(*c.find("data"));
        const Array& shape = items(c, "shape_xyz");
        map.shapeXyz = {integer(shape.at(0)), integer(shape.at(1)), integer(shape.at(2))};
        map.origin = vec3(*c.find("origin"));
        map.voxel = vec3(*c.find("voxel"));
        const fs::path path = dir / "w.mrc";
        io::writeMap(path, map);
        const std::string raw = readFile(path);
        EXPECT_EQ(raw.size(), static_cast<std::size_t>(integer(*c.find("size")))) << context;
        EXPECT_EQ(hex(raw.substr(0, 224)), str(*c.find("header"))) << context;
        EXPECT_EQ(core::sha256Hex(raw), str(*c.find("sha256"))) << context;
        const DensityMap back = io::readMap(path);
        expectVec(back.origin, *c.find("read")->find("origin"), context + " origin");
        expectVec(back.voxel, *c.find("read")->find("voxel"), context + " voxel");
        EXPECT_EQ(dataSha(back.data), str(*c.find("read")->find("sha256"))) << context;
    }
    DensityMap empty;
    empty.shapeXyz = {2, 2, 0};
    try {
        io::writeMap(dir / "empty.mrc", empty);
        ADD_FAILURE() << "an empty map was written";
    } catch (const core::ValueError& e) {
        EXPECT_EQ(e.what(), str(*data.find("empty_error")));
    }
    EXPECT_FALSE(fs::exists(dir / "empty.mrc"));
    DensityMap bad;
    bad.shapeXyz = {2, 2, 2};
    bad.data = {1, 2, 3};
    EXPECT_THROW(io::writeMap(dir / "bad.mrc", bad), std::invalid_argument);
}

TEST(Mrc, SimulateMatchesMdatk) {
    const Value& data = load("simulate_cases.json");
    EXPECT_EQ(io::kSigmaFactor, number(*data.find("sigma_factor")));
    const Value& weights = *data.find("atom_weights");
    std::vector<std::string> elements;
    for (const Value& e : items(weights, "elements")) {
        elements.push_back(str(e));
    }
    const std::vector<double> w = io::atomWeights(elements);
    const Array& expectedWeights = items(weights, "weights");
    ASSERT_EQ(w.size(), expectedWeights.size());
    for (std::size_t k = 0; k < w.size(); ++k) {
        EXPECT_EQ(w[k], number(expectedWeights[k])) << elements[k];
    }
    for (const Value& c : items(data, "cases")) {
        std::vector<Vec3> xyz;
        for (const Value& p : items(c, "xyz")) {
            xyz.push_back(vec3(p));
        }
        std::vector<double> atomWeights;
        if (const Value* els = c.find("elements")) {
            std::vector<std::string> names;
            for (const Value& e : list(*els)) {
                names.push_back(str(e));
            }
            atomWeights = io::atomWeights(names);
        } else {
            for (const Value& v : items(c, "weights")) {
                atomWeights.push_back(number(v));
            }
        }
        std::optional<double> spacing;
        if (const Value* s = c.find("spacing")) {
            spacing = number(*s);
        }
        const double sigmaFactor =
            c.find("sigma_factor") ? number(*c.find("sigma_factor")) : io::kSigmaFactor;
        const double pad = c.find("pad_sigmas") ? number(*c.find("pad_sigmas")) : 5.0;
        const DensityMap map = io::simulateMap(xyz, atomWeights, number(*c.find("resolution")),
                                               spacing, sigmaFactor, pad);
        const std::string context = std::to_string(xyz.size()) + " atoms";
        const Array& shape = items(c, "shape_xyz");
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_EQ(map.shapeXyz[k], integer(shape.at(k))) << context;
        }
        expectVec(map.origin, *c.find("origin"), context + " origin");
        expectVec(map.voxel, *c.find("voxel"), context + " voxel");
        if constexpr (kGlibcExp) {
            EXPECT_EQ(dataSha(map.data), str(*c.find("sha256"))) << context;
            expectStats(map.stats(), *c.find("stats"), context);
        }
    }
    for (const Value& e : items(data, "errors")) {
        const std::vector<Vec3> xyz(static_cast<std::size_t>(integer(*e.find("atoms"))),
                                    Vec3{0, 0, 0});
        const std::vector<double> ones(xyz.size(), 1.0);
        try {
            io::simulateMap(xyz, ones, number(*e.find("resolution")));
            ADD_FAILURE() << "no error";
        } catch (const io::FitError& error) {
            EXPECT_EQ(error.what(), str(*e.find("error")));
            EXPECT_STREQ(error.pythonName(), "FitError");
        }
    }
}

// ---- behaviour beyond the case files ----------------------------------------------------------

TEST(Mrc, LargeStoredGzipAndMissingFiles) {
    TempDir dir;
    DensityMap map;
    map.shapeXyz = {50, 41, 37};
    Value spec = core::Table{{"kind", "uniform"},
                             {"seed", std::int64_t{42}},
                             {"n", std::int64_t{50 * 41 * 37}},
                             {"scale", 2.0},
                             {"offset", 0.0}};
    map.data = generate(spec);
    map.origin = {-3.0, 1.5, 6.0};
    map.voxel = {1.5, 1.5, 1.5};
    io::writeMap(dir / "m.mrc", map);
    const std::string raw = readFile(dir / "m.mrc");
    // stored blocks: more compressed input than one read chunk, several output buffers
    writeFile(dir / "m.mrc.gz", gzipStored(raw, 40000) + gzipStored(raw.substr(0, 10)));
    const DensityMap back = io::readMap(dir / "m.mrc.gz");
    EXPECT_EQ(back.data, map.data);
    EXPECT_EQ(back.origin, map.origin);
    writeFile(dir / "cut.mrc.gz", gzipStored(raw, 40000).substr(0, 90000));
    EXPECT_EQ(mapError(dir / "cut.mrc.gz"),
              core::pathToUtf8(dir / "cut.mrc.gz") + ": file is shorter than its header says");
    EXPECT_THROW(io::readMap(dir / "missing.mrc"), fs::filesystem_error);
    EXPECT_THROW(io::readMap(dir.path()), fs::filesystem_error);
    EXPECT_THROW(io::loadMapCached(dir / "missing.mrc"), fs::filesystem_error);
    EXPECT_EQ(back.at(1, 2, 3), map.data[(3 * 41 + 2) * 50 + 1]);
    EXPECT_THROW((void)back.at(50, 0, 0), std::out_of_range);
}

TEST(Mrc, LoadMapCachedSharesMaps) {
    TempDir dir;
    DensityMap map;
    map.shapeXyz = {3, 2, 2};
    map.data = std::vector<float>(12, 1.0f);
    map.voxel = {1, 1, 1};
    io::writeMap(dir / "a.mrc", map);
    io::writeMap(dir / "b.mrc", map);
    io::writeMap(dir / "c.mrc", map);
    const auto a1 = io::loadMapCached(dir / "a.mrc");
    const auto a2 = io::loadMapCached(dir / "a.mrc");
    EXPECT_EQ(a1, a2);
    const auto b1 = io::loadMapCached(dir / "b.mrc");
    EXPECT_EQ(io::loadMapCached(dir / "a.mrc"), a1);
    io::loadMapCached(dir / "c.mrc");  // two maps are kept: a.mrc, the oldest, goes
    EXPECT_NE(io::loadMapCached(dir / "a.mrc"), a1);
    map.data[0] = 5.0f;
    io::writeMap(dir / "b.mrc", map);
    fs::last_write_time(dir / "b.mrc",
                        fs::last_write_time(dir / "b.mrc") + std::chrono::seconds(5));
    const auto b2 = io::loadMapCached(dir / "b.mrc");
    EXPECT_NE(b2, b1);
    EXPECT_EQ(b2->data[0], 5.0f);
}

TEST(Mrc, MapMethods) {
    DensityMap map;
    map.shapeXyz = {4, 3, 2};
    map.data.resize(24);
    map.origin = {1.0, -2.0, 0.5};
    map.voxel = {0.5, 2.0, 1.5};
    const Vec3 extent = map.extent();
    EXPECT_EQ(extent, (Vec3{1.0 + 3 * 0.5, -2.0 + 2 * 2.0, 0.5 + 1.5}));
    const Vec3 xyz = map.indexToXyz({1.5, 2.0, -1.0});
    EXPECT_EQ(xyz, (Vec3{1.0 + 1.5 * 0.5, -2.0 + 2.0 * 2.0, 0.5 - 1.5}));
    EXPECT_EQ(map.xyzToIndex(xyz), (Vec3{1.5, 2.0, -1.0}));
    DensityMap empty;
    EXPECT_THROW((void)empty.stats(), core::ValueError);
}

// ---- MDATK's tests (tests/test_density.py) ---------------------------------------------------

// test_read_map_axis_orders_and_byte_order
TEST(Mrc, MdatkAxisOrdersAndByteOrder) {
    TempDir dir;
    const std::array<int, 3> shape{14, 12, 10};  // z, y, x
    const std::vector<float> cart = blob(shape, {3, 7, 9});
    const auto cartAt = [&](int x, int y, int z) {
        return cart[static_cast<std::size_t>((z * 12 + y) * 10 + x)];
    };
    const std::array<std::array<int, 3>, 4> orders{{{1, 2, 3}, {3, 1, 2}, {2, 3, 1}, {3, 2, 1}}};
    for (const auto& mapcrs : orders) {
        for (const bool big : {false, true}) {
            const std::array<int, 3> nXyz{10, 12, 14};
            RawMap raw;
            raw.shapeSrc = {nXyz[static_cast<std::size_t>(mapcrs[2] - 1)],
                            nXyz[static_cast<std::size_t>(mapcrs[1] - 1)],
                            nXyz[static_cast<std::size_t>(mapcrs[0] - 1)]};
            for (int s = 0; s < raw.shapeSrc[0]; ++s) {
                for (int r = 0; r < raw.shapeSrc[1]; ++r) {
                    for (int c = 0; c < raw.shapeSrc[2]; ++c) {
                        std::array<int, 3> p{};
                        p[static_cast<std::size_t>(mapcrs[0] - 1)] = c;
                        p[static_cast<std::size_t>(mapcrs[1] - 1)] = r;
                        p[static_cast<std::size_t>(mapcrs[2] - 1)] = s;
                        raw.dataSrc.push_back(cartAt(p[0], p[1], p[2]));
                    }
                }
            }
            raw.mapcrs = mapcrs;
            raw.voxel = {0.8, 1.1, 1.3};
            raw.big = big;
            raw.origin = {10.0, -5.0, 2.0};
            raw.nsymbt = 64;
            writeFile(dir / "m.mrc", raw.bytes());
            const DensityMap m = io::readMap(dir / "m.mrc");
            EXPECT_EQ(m.shapeXyz, (std::array<std::int64_t, 3>{10, 12, 14}));
            EXPECT_EQ(m.data, cart);
            const auto peak = static_cast<std::size_t>(
                std::max_element(m.data.begin(), m.data.end()) - m.data.begin());
            const Vec3 index{double(peak % 10), double(peak / 10 % 12), double(peak / 120)};
            const Vec3 xyz = m.indexToXyz(index);
            EXPECT_NEAR(xyz[0], 10 + 3 * 0.8, 1e-5);
            EXPECT_NEAR(xyz[1], -5 + 7 * 1.1, 1e-5);
            EXPECT_NEAR(xyz[2], 2 + 9 * 1.3, 1e-5);
        }
    }
}

// test_read_map_nstart_origin_modes_and_errors
TEST(Mrc, MdatkNstartOriginModesAndErrors) {
    TempDir dir;
    RawMap raw;
    raw.shapeSrc = {8, 8, 8};
    for (const float v : blob({8, 8, 8}, {4, 4, 4})) {
        raw.dataSrc.push_back(static_cast<double>(static_cast<std::int16_t>(v * 100)));
    }
    raw.voxel = {2.0, 2.0, 2.0};
    raw.nstart = {-4, 2, 10};
    raw.mode = 1;
    writeFile(dir / "s.map", raw.bytes());
    const DensityMap m = io::readMap(dir / "s.map");
    EXPECT_EQ(m.origin, (Vec3{-8.0, 4.0, 20.0}));
    EXPECT_EQ(*std::max_element(m.data.begin(), m.data.end()), 100.0f);
    writeFile(dir / "bad.mrc", std::string(100, '\0'));
    EXPECT_THROW(io::readMap(dir / "bad.mrc"), io::MapError);
    writeFile(dir / "trunc.mrc", raw.bytes().substr(0, 1100));
    EXPECT_NE(mapError(dir / "trunc.mrc").find("shorter"), std::string::npos);
}

// test_write_read_roundtrip
TEST(Mrc, MdatkWriteReadRoundtrip) {
    TempDir dir;
    DensityMap m;
    m.shapeXyz = {7, 6, 5};
    Value spec = core::Table{{"kind", "uniform"},
                             {"seed", std::int64_t{0}},
                             {"n", std::int64_t{210}},
                             {"scale", 3.0},
                             {"offset", 0.0}};
    m.data = generate(spec);
    m.origin = {1.5, -2.0, 3.25};
    m.voxel = {1.1, 1.1, 1.1};
    io::writeMap(dir / "r.mrc", m);
    const DensityMap back = io::readMap(dir / "r.mrc");
    EXPECT_EQ(back.data, m.data);
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_NEAR(back.origin[k], m.origin[k], 1e-6);
        EXPECT_NEAR(back.voxel[k], m.voxel[k], 1e-6 * m.voxel[k]);
    }
}

// test_simulated_map_peak_and_width
TEST(Mrc, MdatkSimulatedPeakAndWidth) {
    const std::vector<Vec3> xyz{{0, 0, 0}, {12, 0, 0}};
    const DensityMap m = io::simulateMap(xyz, std::vector<double>{6.0, 6.0}, 4.0, 0.5);
    const double sigma = io::kSigmaFactor * 4.0;
    for (const Vec3& atom : xyz) {
        const Vec3 f = m.xyzToIndex(atom);
        const auto i = static_cast<std::int64_t>(std::round(f[0])),
                   j = static_cast<std::int64_t>(std::round(f[1])),
                   k = static_cast<std::int64_t>(std::round(f[2]));
        float local = 0;
        for (std::int64_t z = k - 3; z <= k + 3; ++z) {
            for (std::int64_t y = j - 3; y <= j + 3; ++y) {
                for (std::int64_t x = i - 3; x <= i + 3; ++x) {
                    local = std::max(local, m.at(x, y, z));
                }
            }
        }
        EXPECT_EQ(local, m.at(i, j, k));
    }
    const Vec3 f = m.xyzToIndex(xyz[0]);
    const auto j = static_cast<std::int64_t>(std::round(f[1])),
               k = static_cast<std::int64_t>(std::round(f[2]));
    for (std::int64_t x = 0; x < m.shapeXyz[0]; ++x) {
        const double xs = m.origin[0] + m.voxel[0] * static_cast<double>(x);
        if (std::abs(xs) < 4) {
            const double expected = 6.0 * std::exp(-xs * xs / (2 * sigma * sigma));
            EXPECT_NEAR(m.at(x, j, k), expected, 1e-5 + 1e-4 * expected);
        }
    }
}

// test_gzip_imod_bytes_and_nstart
TEST(Mrc, MdatkGzipImodBytesAndNstart) {
    TempDir dir;
    DensityMap cart;
    cart.shapeXyz = {6, 6, 6};
    cart.data = blob({6, 6, 6}, {2, 3, 4});
    cart.origin = {-3.0, 1.5, 6.0};
    cart.voxel = {1.5, 1.5, 1.5};
    io::writeMap(dir / "g.mrc", cart);
    const std::string raw = readFile(dir / "g.mrc");
    EXPECT_EQ(raw[16], static_cast<char>(-2));
    EXPECT_EQ(raw[20], 1);
    EXPECT_EQ(raw[24], 4);
    writeFile(dir / "g.map.gz", gzipStored(raw));
    const DensityMap m = io::readMap(dir / "g.map.gz");
    EXPECT_EQ(m.data, cart.data);
    EXPECT_EQ(m.origin, cart.origin);
    // IMOD mode 0 maps hold unsigned bytes
    RawMap imod;
    imod.shapeSrc = {2, 2, 2};
    imod.dataSrc.assign(8, -56.0);  // 200 as int8
    imod.mode = 0;
    std::string bytes = imod.bytes();
    putInt(bytes, 152, 1146047817, false);
    putInt(bytes, 156, 0, false);
    writeFile(dir / "imod.mrc", bytes);
    const DensityMap u = io::readMap(dir / "imod.mrc");
    EXPECT_EQ(*std::min_element(u.data.begin(), u.data.end()), 200.0f);
    const DensityMap sim =
        io::simulateMap(std::vector<Vec3>{{0.3, 0.7, -1.1}}, std::vector<double>{6.0}, 3.0);
    for (std::size_t k = 0; k < 3; ++k) {
        const double start = sim.origin[k] / sim.voxel[k];
        EXPECT_NEAR(start, std::round(start), 1e-8);
    }
}

// test_cli_density (the info part)
TEST(Mrc, MdatkDensityInfo) {
    TempDir dir;
    DensityMap m;
    m.shapeXyz = {20, 12, 9};
    m.data = blob({9, 12, 20}, {10, 6, 4}, 2.0);
    m.origin = {-12.0, 3.0, 4.5};
    m.voxel = {1.5, 1.5, 1.5};
    io::writeMap(dir / "t.mrc", m);
    const Value info = io::densityInfo(dir / "t.mrc");
    EXPECT_EQ(*info.find("voxel_A"), Value(Array{1.5, 1.5, 1.5}));
    EXPECT_GT(integer(list(*info.find("grid_xyz")).at(0)), 10);
    std::vector<std::string> keys;
    for (const auto& [key, value] : *info.get<core::Table>()) {
        keys.push_back(key);
    }
    EXPECT_EQ(keys, (std::vector<std::string>{"map", "grid_xyz", "voxel_A", "origin_A", "extent_A",
                                              "min", "max", "mean", "rms", "level_2sigma"}));
}
