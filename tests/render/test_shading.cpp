// Residue shading, orientation, the example dataset and numpy's PCG64 (src/render/core shading.cpp,
// orientation.cpp, example.cpp, numpy_random.cpp) against MDATK v0.0.4.1 and numpy 2.5.5
// (tests/data/render/shading, written by tests/parity/make_shading_cases.py with MDBridge's
// improvements RN-21 ... RN-25 from tests/parity/improvements_shading.py). Counterparts of MDATK's
// test_render_cache.py (glow modes, static colours, Euler round trip / drag), test_v002.py
// (example dataset matches the layout) and test_seqpanels.py (one-letter codes, residue rows).
//
// Random numbers, colours, glow, tags and the example arrays are compared bit for bit; numbers that
// go through the C library's sin / cos / exp (rotations, schematic coordinates, example waves) are
// bit for bit with glibc (the generator's numpy uses the same libm) and within a few ulps elsewhere.

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <ranges>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/analysis/cache.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/example.hpp"
#include "mdbridge/render/numpy_random.hpp"
#include "mdbridge/render/orientation.hpp"
#include "mdbridge/render/shading.hpp"

namespace {

using mdbridge::core::Array;
using mdbridge::core::NumArray;
using mdbridge::core::Table;
using mdbridge::core::Value;
using namespace mdbridge::render;

#if defined(__GLIBC__)
constexpr bool kExactLibm = true;  // the generator's numpy calls glibc's sin / cos / exp
#else
constexpr bool kExactLibm = false;
#endif

std::filesystem::path dataDir() { return std::filesystem::path(MDBRIDGE_TEST_DATA) / "render/shading"; }

Value unplain(Value v) {
    if (Table* t = v.get<Table>()) {
        if (t->size() == 1 && (*t)[0].first == "$float") {
            const std::string s = *(*t)[0].second.get<std::string>();
            return Value(s == "nan" ? std::nan("") : (s == "inf" ? INFINITY : -INFINITY));
        }
        for (auto& [k, e] : *t) {
            e = unplain(std::move(e));
        }
    } else if (Array* a = v.get<Array>()) {
        for (Value& e : *a) {
            e = unplain(std::move(e));
        }
    }
    return v;
}

Value load(const std::string& name) {
    const std::filesystem::path path = dataDir() / name;
    std::ifstream in(path, std::ios::binary);
    std::ostringstream text;
    text << in.rdbuf();
    return unplain(mdbridge::core::parseJson(text.str(), path.string()));
}

const Value& cases() {
    static const Value doc = load("shading_cases.json");
    return doc;
}
const Value& examples() {
    static const Value doc = load("examples.json");
    return doc;
}
const std::map<std::string, NumArray>& arrays() {
    static const std::map<std::string, NumArray> all = [] {
        std::map<std::string, NumArray> out;
        for (auto& [k, v] : mdbridge::analysis::readNpz(dataDir() / "arrays.npz")) {
            out.emplace(k, std::move(v));
        }
        return out;
    }();
    return all;
}
const NumArray& npz(const std::string& key) {
    const auto it = arrays().find(key);
    if (it == arrays().end()) {
        throw std::runtime_error("arrays.npz has no " + key);
    }
    return it->second;
}

const Value& at(const Value& v, std::string_view key) {
    const Value* found = v.find(key);
    if (found == nullptr) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *found;
}
const Array& arr(const Value& v) { return *v.get<Array>(); }
const std::string& str(const Value& v) { return *v.get<std::string>(); }
double num(const Value& v) {
    if (const auto* i = v.get<std::int64_t>()) {
        return static_cast<double>(*i);
    }
    return *v.get<double>();
}
std::int64_t integer(const Value& v) { return *v.get<std::int64_t>(); }
std::vector<std::int64_t> ints(const Value& v) {
    std::vector<std::int64_t> out;
    for (const Value& e : arr(v)) {
        out.push_back(integer(e));
    }
    return out;
}
std::vector<double> nums(const Value& v) {
    std::vector<double> out;
    for (const Value& e : arr(v)) {
        out.push_back(num(e));
    }
    return out;
}
std::vector<std::string> strs(const Value& v) {
    std::vector<std::string> out;
    for (const Value& e : arr(v)) {
        out.push_back(str(e));
    }
    return out;
}
std::vector<bool> bools(const Value& v) {
    std::vector<bool> out;
    for (const Value& e : arr(v)) {
        out.push_back(*e.get<bool>());
    }
    return out;
}

bool same(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}
// Equal bit for bit, or (without glibc) within `ulps` units in the last place.
bool close(double a, double b, int ulps = 4) {
    if (same(a, b)) {
        return true;
    }
    if (kExactLibm || std::isnan(a) || std::isnan(b)) {
        return false;
    }
    return std::abs(a - b) <= ulps * std::numeric_limits<double>::epsilon() * std::max({std::abs(a), std::abs(b), 1.0});
}

mdbridge::core::Project projectOf(const Value& v) { return mdbridge::core::projectFromPlain(v); }

std::vector<Vec3d> vecs(const Value& v) {
    std::vector<Vec3d> out;
    for (const Value& p : arr(v)) {
        out.push_back({num(arr(p)[0]), num(arr(p)[1]), num(arr(p)[2])});
    }
    return out;
}
std::vector<std::array<std::int64_t, 2>> pairs(const Value& v) {
    std::vector<std::array<std::int64_t, 2>> out;
    for (const Value& p : arr(v)) {
        out.push_back({integer(arr(p)[0]), integer(arr(p)[1])});
    }
    return out;
}

StructurePreview previewOf(const Value& v) {
    StructurePreview p;
    TraceSource& t = p.trace;
    t.coords = vecs(at(v, "coords"));
    t.ions = vecs(at(v, "ions"));
    t.links = pairs(at(v, "links"));
    t.chain = ints(at(v, "chain"));
    t.groupIndex = ints(at(v, "group_index"));
    t.resid = ints(at(v, "resid"));
    t.resname = strs(at(v, "resname"));
    t.chainName = strs(at(v, "chain_name"));
    std::vector<std::vector<bool>> member;
    for (const Value& row : arr(at(v, "group_member"))) {
        member.push_back(bools(row));
    }
    t.groupMember = member;
    t.labels = strs(at(v, "labels"));
    if (const Value* a = v.find("atoms")) {
        AtomsSource& s = p.atoms.emplace();
        s.coords = vecs(at(*a, "coords"));
        s.names = strs(at(*a, "names"));
        s.elements = strs(at(*a, "elements"));
        s.resnames = strs(at(*a, "resnames"));
        s.row = ints(at(*a, "row"));
        s.bonds = pairs(at(*a, "bonds"));
    }
    return p;
}

// --- arrays ---------------------------------------------------------------------------------------


// dtype, shape and values; floats bit for bit (`waves`: within a few ulps without glibc)
void expectArray(const NumArray& got, const NumArray& want, const std::string& what, bool waves = true) {
    ASSERT_EQ(got.dtypeName(), want.dtypeName()) << what;
    ASSERT_EQ(got.shape(), want.shape()) << what;
    if (want.dtype() == NumArray::DType::String) {
        EXPECT_TRUE(std::ranges::equal(got.values<std::string>(), want.values<std::string>())) << what;
        return;
    }
    const std::vector<double> g = got.toDouble(), w = want.toDouble();
    std::size_t bad = 0;
    for (std::size_t k = 0; k < g.size(); ++k) {
        const bool ok = waves ? close(g[k], w[k]) : same(g[k], w[k]);
        if (!ok && bad++ < 5) {
            ADD_FAILURE() << what << "[" << k << "]: " << g[k] << " != " << w[k];
        }
    }
    EXPECT_EQ(bad, 0U) << what;
}

NumArray floats(const FrameMatrix& m) { return NumArray(m.values, {m.rows, m.frames}); }

// --- crc32 and numpy's generator ------------------------------------------------------------------

TEST(NumpyRandom, Crc32OfKeys) {
    for (const Value& c : arr(at(at(cases(), "rng"), "crc32"))) {
        const std::string key = str(arr(c)[0]);
        EXPECT_EQ(crc32(key.data(), key.size()), static_cast<std::uint32_t>(integer(arr(c)[1]))) << key;
    }
}

TEST(NumpyRandom, Pcg64StreamsAreNumpys) {
    for (const Value& s : arr(at(at(cases(), "rng"), "seeds"))) {
        const std::uint64_t seed = std::stoull(str(at(s, "seed")));
        const std::vector<std::uint32_t> state = seedSequenceState(seed, 8);
        const std::vector<std::int64_t> wantState = ints(at(s, "state"));
        for (std::size_t k = 0; k < 8; ++k) {
            EXPECT_EQ(state[k], static_cast<std::uint32_t>(wantState[k])) << "seed " << seed;
        }
        NumpyGenerator raw(seed);
        for (const Value& r : arr(at(s, "raw"))) {
            EXPECT_EQ(raw.next64(), std::stoull(str(r))) << "seed " << seed;
        }
        NumpyGenerator g(seed);
        for (const Value& step : arr(at(s, "script"))) {
            const std::string op = str(arr(step)[0]);
            const Array& a = arr(arr(step)[1]);
            const Value& result = arr(step)[2];
            const std::string what = "seed " + std::to_string(seed) + " " + op;
            if (op == "uniform") {
                EXPECT_TRUE(same(g.uniform(num(a[0]), num(a[1])), num(result))) << what;
            } else if (op == "random") {
                EXPECT_TRUE(same(g.random(), num(result))) << what;
            } else if (op == "normal") {
                const std::vector<double> got = g.normal(num(a[0]), num(a[1]), static_cast<std::size_t>(integer(a[2])));
                const std::vector<double> want = nums(result);
                ASSERT_EQ(got.size(), want.size());
                for (std::size_t k = 0; k < got.size(); ++k) {
                    EXPECT_TRUE(close(got[k], want[k], 1)) << what << "[" << k << "]";
                }
            } else {
                const bool endpoint = a.size() > 2 && *a[2].get<bool>();
                EXPECT_EQ(g.integers(integer(a[0]), integer(a[1]), endpoint), integer(result)) << what;
            }
        }
    }
    for (const std::uint64_t seed : {std::uint64_t{0}, std::uint64_t{42}, (std::uint64_t{1} << 40) + 7}) {
        NumpyGenerator g(seed);
        const std::string p = "rng." + std::to_string(seed) + ".";
        expectArray(NumArray(g.normal(0.0, 1.0, 8000)), npz(p + "normal"), p + "normal");
        expectArray(NumArray(g.uniform(0.0, 1.0, 500)), npz(p + "uniform"), p + "uniform", false);
        expectArray(NumArray(g.integerArray(0, 1000003, 300)), npz(p + "integers"), p + "integers", false);
    }
}

TEST(NumpyRandom, Errors) {
    NumpyGenerator g(1);
    EXPECT_THROW(g.integers(5, 5), mdbridge::core::ValueError);
    EXPECT_THROW(g.normal(0.0, -1.0), mdbridge::core::ValueError);
    EXPECT_THROW(g.uniform(-1e308, 1e308), std::overflow_error);
    EXPECT_EQ(g.integers(7, 7, true), 7);
}

// --- coloring.py ----------------------------------------------------------------------------------

TEST(Shading, GlowModes) {
    const Value& glow = at(cases(), "glow");
    for (const Value& c : arr(at(glow, "cases"))) {
        const Value& input = at(at(glow, "inputs"), str(at(c, "input")));
        std::vector<double> d;
        const Array& rows = arr(at(input, "d"));
        for (const Value& r : rows) {
            for (double x : nums(r)) {
                d.push_back(x);
            }
        }
        const std::size_t frames = d.size() / rows.size();
        const std::string what = str(at(c, "input")) + " " + str(at(c, "mode")) + " fade " +
                                 std::to_string(integer(at(c, "fade_frames"))) + " cutoff " + std::to_string(num(at(c, "cutoff")));
        auto run = [&] {
            return glowIntensity(d, rows.size(), frames, bools(at(input, "member")), num(at(c, "cutoff")), str(at(c, "mode")),
                                 integer(at(c, "fade_frames")));
        };
        if (const Value* e = c.find("error")) {
            try {
                run();
                ADD_FAILURE() << what << ": no error";
            } catch (const mdbridge::core::ValueError& ex) {
                EXPECT_EQ(std::string(ex.what()), str(*e));
            }
            continue;
        }
        const FrameMatrix got = run();
        std::size_t k = 0;
        for (const Value& r : arr(at(c, "glow"))) {
            for (double x : nums(r)) {
                EXPECT_TRUE(same(static_cast<double>(got.values[k]), x)) << what << " [" << k << "]";
                ++k;
            }
        }
    }
}

TEST(Shading, MdatkGlowTest) {  // test_render_cache.test_glow_modes
    const std::vector<double> d{9, 9, 2, 2, 9, 9, 9, 9, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const std::vector<bool> member{false, true};
    const FrameMatrix fade = glowIntensity(d, 2, 10, member, 4.0, "fade", 3);
    const std::vector<float> want{0, 0, 1, 1, 0.75f, 0.5f, 0.25f, 0, 0, 0};
    for (std::size_t f = 0; f < 10; ++f) {
        EXPECT_FLOAT_EQ(fade.at(0, f), want[f]);
        EXPECT_EQ(fade.at(1, f), 0.0f);  // members never glow
    }
    EXPECT_THROW(glowIntensity(d, 2, 10, member, 4.0, "disco"), mdbridge::core::ValueError);
}

TEST(Shading, StaticAndFrameColors) {
    for (const Value& c : arr(at(cases(), "colors"))) {
        const std::string what = str(at(c, "mode")) + " " + str(at(c, "palette"));
        const std::vector<Rgba> got = staticColors(str(at(c, "mode")), ints(at(c, "chain")), ints(at(c, "group_index")),
                                                   strs(at(c, "group_colors")), str(at(c, "base")), str(at(c, "palette")));
        const Array& want = arr(at(c, "static"));
        ASSERT_EQ(got.size(), want.size()) << what;
        auto check = [&](const std::vector<Rgba>& g, const Array& w, const std::string& label) {
            for (std::size_t k = 0; k < g.size(); ++k) {
                const std::vector<double> e = nums(w[k]);
                EXPECT_TRUE(same(g[k].r, e[0]) && same(g[k].g, e[1]) && same(g[k].b, e[2]) && same(g[k].a, e[3]))
                    << label << " residue " << k;
            }
        };
        check(got, want, what);
        if (const Value* frame = c.find("frame")) {
            std::vector<float> intensity;
            for (double x : nums(at(c, "intensity"))) {
                intensity.push_back(static_cast<float>(x));
            }
            for (const char* mode : {"proximity", "groups"}) {
                check(frameColors(got, mode, "#b5b3ad", "#e34948", intensity, bools(at(c, "member"))), arr(at(*frame, mode)),
                      what + " frame " + mode);
            }
        }
    }
}

TEST(Shading, OneLetterCodes) {  // test_seqpanels.test_one_letter_codes and more
    for (const Value& c : arr(at(cases(), "one_letter"))) {
        EXPECT_EQ(oneLetter(str(arr(c)[0])), str(arr(c)[1])) << str(arr(c)[0]);
    }
}

// --- orientation.py -------------------------------------------------------------------------------

void expectMatrix(const Mat3& got, const Value& want, const std::string& what, double tol = 0.0) {
    for (std::size_t i = 0; i < 3; ++i) {
        const std::vector<double> row = nums(arr(want)[i]);
        for (std::size_t j = 0; j < 3; ++j) {
            if (tol > 0.0) {
                EXPECT_NEAR(got[i][j], row[j], tol) << what;
            } else {
                EXPECT_TRUE(close(got[i][j], row[j])) << what << " [" << i << "][" << j << "] " << got[i][j] << " " << row[j];
            }
        }
    }
}

Vec3d vec(const Value& v) { return {num(arr(v)[0]), num(arr(v)[1]), num(arr(v)[2])}; }

Mat3 matrix(const Value& v) {
    Mat3 m{};
    for (std::size_t i = 0; i < 3; ++i) {
        const std::vector<double> row = nums(arr(v)[i]);
        m[i] = {row[0], row[1], row[2]};
    }
    return m;
}

TEST(Orientation, MatchesMdatk) {
    const Value& o = at(cases(), "orientation");
    for (const Value& c : arr(at(o, "rx"))) {
        const double a = num(arr(c)[0]);
        expectMatrix(rotationX(a), arr(c)[1], "rx");
        expectMatrix(rotationY(a), arr(c)[2], "ry");
        expectMatrix(rotationZ(a), arr(c)[3], "rz");
    }
    // products: numpy's matmul goes through BLAS, whose summation order is its own
    for (const Value& c : arr(at(o, "euler"))) {
        expectMatrix(eulerToMatrix(vec(arr(c)[0])), arr(c)[1], "euler", 1e-15);
    }
    for (const Value& c : arr(at(o, "to_euler"))) {
        const Vec3d got = matrixToEuler(matrix(arr(c)[0]));
        const Vec3d want = vec(arr(c)[1]);
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_NEAR(got[k], want[k], 1e-9) << "matrix_to_euler";
        }
    }
    for (const Value& c : arr(at(o, "wrap"))) {
        EXPECT_TRUE(same(wrapDegrees(num(arr(c)[0])), num(arr(c)[1]))) << num(arr(c)[0]);
    }
    for (const Value& c : arr(at(o, "drag"))) {
        const double dpp = arr(c).size() > 4 ? num(arr(c)[4]) : 0.5;
        const Vec3d got = drag(vec(arr(c)[0]), num(arr(c)[1]), num(arr(c)[2]), dpp);
        const Vec3d want = vec(arr(c)[3]);
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_NEAR(got[k], want[k], 1e-9) << "drag";
        }
    }
    for (const Value& c : arr(at(o, "format"))) {
        EXPECT_EQ(formatAngles(vec(arr(c)[0])), str(arr(c)[1]));
    }
    for (const Value& c : arr(at(o, "project"))) {
        const std::vector<Vec3d> got = project(vecs(arr(c)[1]), vec(arr(c)[0]), vec(arr(c)[2]));
        const std::vector<Vec3d> want = vecs(arr(c)[3]);
        for (std::size_t k = 0; k < got.size(); ++k) {
            for (std::size_t i = 0; i < 3; ++i) {
                EXPECT_NEAR(got[k][i], want[k][i], 1e-12) << "project";
            }
        }
    }
}

TEST(Orientation, EulerRoundTripAndDrag) {  // test_render_cache.test_euler_roundtrip_and_drag
    NumpyGenerator rng(0);
    for (int k = 0; k < 50; ++k) {
        Vec3d ang{rng.uniform(-179, 179), rng.uniform(-179, 179), rng.uniform(-179, 179)};
        ang[1] = rng.uniform(-85, 85);
        const Vec3d back = matrixToEuler(eulerToMatrix(ang));
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_NEAR(back[i], ang[i], 1e-6);
        }
    }
    const Mat3 r = eulerToMatrix(drag({10, 20, 30}, 40, -15));
    const Mat3 rrt = multiply(r, transpose(r));
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(rrt[i][j], i == j ? 1.0 : 0.0, 1e-9);
        }
    }
    const Vec3d still = drag({0, 0, 0}, 0, 0);
    EXPECT_EQ(still, (Vec3d{0, 0, 0}));
    EXPECT_THROW(eulerToMatrix(std::vector<double>{1.0, 2.0}), mdbridge::core::ValueError);
}

TEST(Orientation, AnglesNextToMinus180StayRN22) {
    // MDATK's np.isclose turned -179.999 into +180 (a 0.001 degree jump)
    EXPECT_EQ(wrapDegrees(-179.999), -179.999);
    EXPECT_EQ(wrapDegrees(-180.0), 180.0);
    EXPECT_EQ(wrapDegrees(540.0), 180.0);
    const Vec3d back = matrixToEuler(eulerToMatrix(Vec3d{-179.999, 0.0, 0.0}));
    EXPECT_NEAR(back[0], -179.999, 1e-9);
}

// --- example.py -----------------------------------------------------------------------------------

struct ExampleCase {
    std::string name;
    mdbridge::core::Project project;
    mdbridge::core::Dataset dataset;
};

const std::vector<ExampleCase>& exampleCases() {
    static const std::vector<ExampleCase> all = [] {
        std::vector<ExampleCase> out;
        for (const Value& c : arr(at(examples(), "cases"))) {
            mdbridge::core::Project p = projectOf(at(c, "project"));
            std::optional<StructurePreview> preview;
            if (!at(c, "preview").isNull()) {
                preview = previewOf(at(c, "preview"));
            }
            mdbridge::core::Dataset ds =
                exampleDataset(p, preview ? &*preview : nullptr, static_cast<int>(integer(at(c, "n_frames"))));
            out.push_back({str(at(c, "case")), std::move(p), std::move(ds)});
        }
        return out;
    }();
    return all;
}

TEST(ExampleDataset, MatchesMdatk) {
    const Array& want = arr(at(examples(), "cases"));
    ASSERT_EQ(want.size(), exampleCases().size());
    for (std::size_t c = 0; c < want.size(); ++c) {
        const ExampleCase& ex = exampleCases()[c];
        const Value& d = at(want[c], "dataset");
        SCOPED_TRACE(ex.name);
        EXPECT_EQ(ex.dataset.frames(), ints(at(d, "frames")));
        const std::vector<double> t = nums(at(d, "time_ps"));
        ASSERT_EQ(ex.dataset.timePs().size(), t.size());
        for (std::size_t k = 0; k < t.size(); ++k) {
            EXPECT_TRUE(same(ex.dataset.timePs()[k], t[k]));
        }
        EXPECT_EQ(ex.dataset.timeUnit(), str(at(d, "time_unit")));
        const Array& results = arr(at(d, "results"));
        ASSERT_EQ(ex.dataset.results().size(), results.size());
        for (std::size_t r = 0; r < results.size(); ++r) {
            const mdbridge::core::AnalysisResult& got = ex.dataset.results()[r];
            const Value& w = results[r];
            EXPECT_EQ(got.name, str(at(w, "name")));
            EXPECT_EQ(got.key, str(at(w, "key")));
            EXPECT_EQ(Value(got.params), at(w, "params")) << got.name;
            EXPECT_EQ(Value(got.meta), at(w, "meta")) << got.name << ": " << mdbridge::core::toJson(Value(got.meta), -1);
            const std::vector<std::string> names = strs(at(w, "arrays"));
            ASSERT_EQ(got.arrays.size(), names.size()) << got.name;
            for (std::size_t a = 0; a < names.size(); ++a) {
                EXPECT_EQ(got.arrays[a].first, names[a]);
                expectArray(got.arrays[a].second, npz(ex.name + "." + str(at(w, "name")) + "." + names[a]),
                            ex.name + "." + got.name + "." + names[a]);
            }
        }
    }
}

TEST(ExampleDataset, MdatkV002Layout) {  // test_v002.test_example_dataset_matches_layout
    const ExampleCase& ex = exampleCases()[1];
    ASSERT_EQ(ex.name, "v002");
    const std::vector<std::string> keys = ex.dataset.seriesKeys();
    EXPECT_EQ(std::set<std::string>(keys.begin(), keys.end()),
              (std::set<std::string>{"stacking@Aptamer+Target/count", "rmsd/Complex", "energy/total"}));
    for (double v : ex.dataset.series("stacking@Aptamer+Target/count").values) {
        EXPECT_EQ(v, std::round(v));
    }
    for (double v : ex.dataset.matrix("mg/occupancy").values.toDouble()) {
        EXPECT_LE(v, 2.0);
    }
    EXPECT_EQ(ex.dataset.pairmap("hbonds/map").shape(), (std::pair<std::size_t, std::size_t>{12, 12}));
}

TEST(ExampleDataset, GroupResids) {
    std::map<std::string, const ExampleCase*> byName;
    for (const ExampleCase& ex : exampleCases()) {
        byName[ex.name] = &ex;
    }
    for (const Value& c : arr(at(cases(), "group_resids"))) {
        const ExampleCase& ex = *byName.at(str(arr(c)[0]));
        EXPECT_EQ(groupResids(ex.project, str(arr(c)[1]), static_cast<int>(integer(arr(c)[2]))), ints(arr(c)[3]))
            << ex.name << " " << str(arr(c)[1]);
    }
}

TEST(ExampleDataset, SingleStrandAndLargeGroupsRN21RN25) {
    std::map<std::string, const ExampleCase*> byName;
    for (const ExampleCase& ex : exampleCases()) {
        byName[ex.name] = &ex;
    }
    // RN-21: MDATK raised an IndexError for a structure with one strand
    EXPECT_EQ(byName.at("single_strand")->dataset.pairmap("basepairs/map").i.size(), 0U);
    // RN-25: a range of a billion residues lists 400, the sequence repeats beyond 120 residues
    const auto* tr = byName.at("big_groups")->dataset.result("trace");
    ASSERT_NE(tr, nullptr);
    EXPECT_EQ(tr->array("resid")->size(), 400U + 152U);
    EXPECT_THROW(exampleDataset(byName.at("default")->project, nullptr, 0), mdbridge::core::ValueError);
}

// --- residues.py ----------------------------------------------------------------------------------

struct LogCapture {
    std::vector<std::string> lines;
    mdbridge::core::LogSink previous;
    LogCapture() {
        previous = mdbridge::core::setLogSink([this](mdbridge::core::LogLevel level, std::string_view message) {
            lines.push_back(std::string(mdbridge::core::logLevelName(level)) + " " + std::string(message));
        });
    }
    ~LogCapture() { mdbridge::core::setLogSink(previous); }
};

TEST(Residues, RowsTagsAndShading) {
    std::map<std::string, const ExampleCase*> byName;
    for (const ExampleCase& ex : exampleCases()) {
        byName[ex.name] = &ex;
    }
    for (const Value& rec : arr(at(cases(), "residues"))) {
        const ExampleCase& ex = *byName.at(str(at(rec, "case")));
        SCOPED_TRACE(ex.name);
        const mdbridge::core::AnalysisResult& tr = *ex.dataset.result("trace");
        EXPECT_EQ(chainNames(tr), strs(at(rec, "chain_names")));
        EXPECT_EQ(residueNames(tr), strs(at(rec, "resnames")));
        for (const Value& c : arr(at(rec, "rows"))) {
            const std::string spec = str(arr(c)[0]);
            if (const Value* e = arr(c)[1].find("error")) {
                try {
                    residueRows(tr, spec);
                    ADD_FAILURE() << spec << ": no error";
                } catch (const mdbridge::core::SelectionError& ex2) {
                    EXPECT_EQ(std::string(ex2.what()), str(*e)) << spec;
                }
            } else {
                EXPECT_EQ(residueRows(tr, spec), ints(arr(c)[1])) << spec;
            }
        }
        for (const Value& c : arr(at(rec, "tags"))) {
            mdbridge::core::VisualParams v = ex.project.visual;
            v.label_every = integer(arr(c)[0]);
            v.label_residues = str(arr(c)[1]);
            resetWarnOnce();
            LogCapture log;
            EXPECT_EQ(tagRows(tr, v), ints(arr(c)[2])) << v.label_every << " '" << v.label_residues << "'";
            EXPECT_EQ(log.lines, strs(arr(c)[3]));
        }
        for (const Value& s : arr(at(rec, "shading"))) {
            mdbridge::core::Project p = ex.project;
            for (const auto& [k, val] : *at(s, "visual").get<Table>()) {
                if (k == "color_mode") p.visual.color_mode = str(val);
                if (k == "glow_mode") p.visual.glow_mode = str(val);
                if (k == "proximity_cutoff") p.visual.proximity_cutoff = num(val);
                if (k == "proximity_color") p.visual.proximity_color = str(val);
                if (k == "glow_color") p.visual.glow_color = str(val);
                if (k == "fade_frames") p.visual.fade_frames = integer(val);
                if (k == "chain_palette") p.visual.chain_palette = str(val);
            }
            const Value& need = at(s, "need");
            const std::optional<bool> needProx = need.isNull() ? std::nullopt : std::optional<bool>(*need.get<bool>());
            const Shading sh = residueShading(p, ex.dataset, tr, needProx);
            const std::string key = str(at(s, "key"));
            std::vector<double> flat;
            for (const Rgba& c : sh.staticColors) {
                flat.insert(flat.end(), {c.r, c.g, c.b, c.a});
            }
            expectArray(NumArray(flat, {sh.staticColors.size(), 4}), npz(key + ".static"), key + ".static", false);
            expectArray(NumArray(sh.member), npz(key + ".member"), key + ".member", false);
            expectArray(floats(sh.glow), npz(key + ".glow"), key + ".glow", false);
            expectArray(floats(sh.colorIntensity), npz(key + ".color_i"), key + ".color_i", false);
            const std::vector<double> rgb = nums(at(s, "glow_rgb"));
            EXPECT_TRUE(same(sh.glowRgb[0], rgb[0]) && same(sh.glowRgb[1], rgb[1]) && same(sh.glowRgb[2], rgb[2]));
        }
    }
}

TEST(Residues, ProximityRowsMatchedByResindex) {
    const Value& c = at(cases(), "remap");
    const mdbridge::core::Project p = projectOf(at(c, "project"));
    const mdbridge::core::Dataset base = exampleDataset(p, nullptr, static_cast<int>(integer(at(c, "n_frames"))));
    std::vector<mdbridge::core::AnalysisResult> results = base.results();
    for (auto& r : results) {
        if (r.name == "proximity") {
            r.setArray("resindex", npz("remap.resindex"));
            r.setArray("min_distance", npz("remap.min_distance"));
            r.setArray("member", npz("remap.prox_member"));
        }
    }
    const mdbridge::core::Dataset ds(base.frames(), base.timePs(), results, base.timeUnit());
    const Shading sh = residueShading(p, ds, *ds.result("trace"));
    expectArray(NumArray(sh.member), npz("remap.member"), "member", false);
    expectArray(floats(sh.glow), npz("remap.glow"), "glow", false);
    expectArray(floats(sh.colorIntensity), npz("remap.color_i"), "color_i", false);
    // no proximity result: KeyError('proximity') in MDATK
    const mdbridge::core::Dataset none(base.frames(), base.timePs(), {base.results()[0]}, base.timeUnit());
    EXPECT_THROW(residueShading(p, none, *ds.result("trace"), true), std::out_of_range);
    EXPECT_EQ(str(at(at(c, "missing"), "type")), "KeyError");
}

TEST(Residues, ChainQualifiedTagsRN24) {
    const ExampleCase& ex = exampleCases()[0];
    const mdbridge::core::AnalysisResult& tr = *ex.dataset.result("trace");
    mdbridge::core::VisualParams v = ex.project.visual;
    v.label_residues = "B:13-15";
    EXPECT_EQ(tagRows(tr, v), (std::vector<std::int64_t>{12, 13, 14}));
    v.label_residues = "A:13-15";  // MDATK tagged the B residues anyway
    EXPECT_TRUE(tagRows(tr, v).empty());
}

TEST(Residues, SchematicSourcesMatchMdatk) {
    const Value& s = at(cases(), "schematic");
    const mdbridge::core::Project p = exampleCases()[0].project;
    const TraceSource tr = schematicTrace(p);
    const Value& wt = at(s, "trace");
    const Array& frames = arr(at(wt, "coords"));
    const std::vector<Vec3d> want = vecs(frames[0]);
    ASSERT_EQ(tr.coords.size(), want.size());
    for (std::size_t k = 0; k < want.size(); ++k) {
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_TRUE(close(tr.coords[k][i], want[k][i]));
        }
    }
    EXPECT_EQ(*tr.resname, strs(at(wt, "resname")));
    const StructurePreview st = schematicStructure(p);
    const Value& wa = at(at(s, "structure"), "atoms");
    EXPECT_EQ(st.atoms->names, strs(at(wa, "names")));
    EXPECT_EQ(st.atoms->elements, strs(at(wa, "elements")));
    EXPECT_EQ(st.atoms->bonds, pairs(at(wa, "bonds")));
    const std::vector<Vec3d> xyz = vecs(at(wa, "coords"));
    ASSERT_EQ(st.atoms->coords.size(), xyz.size());
    for (std::size_t k = 0; k < xyz.size(); ++k) {
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_TRUE(close(st.atoms->coords[k][i], xyz[k][i])) << k;
        }
    }
}

}  // namespace
