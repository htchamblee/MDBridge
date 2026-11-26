// The rigid-body density fit (src/chem/density_fit.cpp) against what MDATK's mdatk/density.py
// computes (tests/parity/make_density_fit_cases.py: data/chem/density_fit_cases.json and the maps of
// data/chem/maps, with MDBridge's CH-11 applied), and the counterparts of MDATK's tests of the fit in
// tests/test_density.py. Fitted poses within owner decision D10's tolerance (0.5 degrees, 0.2 A at the
// model centre, CC 1e-3); the deterministic pieces (rotations, grids, envelopes, scores) much closer:
// MDBridge sums in double precision where MDATK sums float32 samples and runs float32 FFTs of the
// coarse map, so scores agree to about 1e-6.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/density_fit.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/io/mrc.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::chem;
using core::Array;
using core::Mat3;
using core::RigidTransform;
using core::Value;
using core::Vec3;
using test::kData;

namespace {

constexpr double kAngleTolerance = 0.5;  // degrees (D10)
constexpr double kShiftTolerance = 0.2;  // A at the model centre
constexpr double kCcTolerance = 1e-3;
constexpr double kScoreTolerance = 1e-5;  // scores of a given pose: double vs MDATK's float32 sums

const Value& cases() {
    static const Value data = core::parseJson(test::readFile(kData / "chem" / "density_fit_cases.json"),
                                              "density_fit_cases.json");
    return data;
}

const Value& at(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
bool has(const Value& value, std::string_view key) {
    return value.find(key) != nullptr && !value.find(key)->isNull();
}
const Array& list(const Value& value) {
    return *value.get<Array>();
}
std::string str(const Value& value) {
    return *value.get<std::string>();
}
std::size_t count(const Value& value) {
    return static_cast<std::size_t>(*value.get<std::int64_t>());
}
double hexDouble(const Value& value) {
    const std::string text = str(value);
    if (text == "nan") {
        return std::nan("");
    }
    if (text == "inf") {
        return HUGE_VAL;
    }
    if (text == "-inf") {
        return -HUGE_VAL;
    }
    return std::strtod(text.c_str(), nullptr);
}
Vec3 vec(const Value& value) {
    const Array& a = list(value);
    return {hexDouble(a.at(0)), hexDouble(a.at(1)), hexDouble(a.at(2))};
}
Mat3 mat(const Value& value) {
    const Array& a = list(value);
    return {vec(a.at(0)), vec(a.at(1)), vec(a.at(2))};
}
RigidTransform transform(const Value& value) {
    RigidTransform t;
    t.rotation = mat(at(value, "rotation"));
    t.translation = vec(at(value, "translation"));
    return t;
}

struct Model {
    std::vector<Vec3> xyz;
    std::vector<std::string> elements;
};
const Model& model(const std::string& name) {
    static std::map<std::string, Model> models;
    static std::mutex lock;
    const std::lock_guard<std::mutex> guard(lock);
    auto found = models.find(name);
    if (found == models.end()) {
        Model m;
        const Value& v = at(at(cases(), "models"), name);
        for (const Value& p : list(at(v, "xyz"))) {
            m.xyz.push_back(vec(p));
        }
        for (const Value& e : list(at(v, "elements"))) {
            m.elements.push_back(str(e));
        }
        found = models.emplace(name, std::move(m)).first;
    }
    return found->second;
}

double maxAbs(const Mat3& a, const Mat3& b) {
    double out = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out = std::max(out, std::abs(a[i][j] - b[i][j]));
        }
    }
    return out;
}
double maxAbs(const Vec3& a, const Vec3& b) {
    return std::max({std::abs(a[0] - b[0]), std::abs(a[1] - b[1]), std::abs(a[2] - b[2])});
}
double distance(const Vec3& a, const Vec3& b) {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}
// rotation angle of a b^T (degrees)
double angleBetween(const Mat3& a, const Mat3& b) {
    Mat3 rel{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            rel[i][j] = a[i][0] * b[j][0] + a[i][1] * b[j][1] + a[i][2] * b[j][2];
        }
    }
    return core::rotationAngle(rel);
}
Vec3 meanOf(const std::vector<Vec3>& xyz) {
    Vec3 c{0, 0, 0};
    for (const Vec3& p : xyz) {
        for (std::size_t k = 0; k < 3; ++k) {
            c[k] += p[k];
        }
    }
    for (std::size_t k = 0; k < 3; ++k) {
        c[k] /= static_cast<double>(xyz.size());
    }
    return c;
}

// the pose within D10's tolerance: angle and the model centre's position
void expectPose(const RigidTransform& got, const RigidTransform& want, const Vec3& centre, const std::string& what) {
    EXPECT_LT(angleBetween(got.rotation, want.rotation), kAngleTolerance) << what;
    EXPECT_LT(distance(got.apply(centre), want.apply(centre)), kShiftTolerance) << what;
}

std::string errorMessage(const Value& error) {
    return str(at(error, "message"));
}

}  // namespace

// rotation_samples: counts, matrices (to the last bits of sin / cos), and the step errors
TEST(DensityFit, RotationSamplesAsMdatk) {
    for (const Value& c : list(at(cases(), "rotations"))) {
        const double step = hexDouble(at(c, "step"));
        SCOPED_TRACE(str(at(c, "step")));
        if (has(c, "error")) {
            try {
                (void)rotationSamples(step);
                ADD_FAILURE() << "no error";
            } catch (const FitError& e) {
                EXPECT_EQ(e.what(), errorMessage(at(c, "error")));
            }
            continue;
        }
        const std::vector<Mat3> r = rotationSamples(step);
        ASSERT_EQ(r.size(), count(at(c, "count")));
        EXPECT_LT(maxAbs(r.front(), mat(at(c, "first"))), 1e-14);
        EXPECT_LT(maxAbs(r[1], mat(at(c, "second"))), 1e-14);
        EXPECT_LT(maxAbs(r.back(), mat(at(c, "last"))), 1e-12);
        if (has(c, "all")) {
            const Array& all = list(at(c, "all"));
            ASSERT_EQ(all.size(), r.size());
            double worst = 0.0;
            for (std::size_t k = 0; k < r.size(); ++k) {
                worst = std::max(worst, maxAbs(r[k], mat(all[k])));
            }
            EXPECT_LT(worst, 1e-13);
        }
    }
}

TEST(DensityFit, FftSizeAsMdatk) {
    for (const Value& c : list(at(cases(), "fft_size"))) {
        const Array& pair = list(c);
        EXPECT_EQ(fftSize(count(pair.at(0))), count(pair.at(1))) << count(pair.at(0));
    }
    EXPECT_EQ(fftSize(0), 0u);
}

TEST(DensityFit, ErrorsAsMdatk) {
    const Model& m = model("helix");
    io::DensityMap small;
    small.shapeXyz = {3, 3, 3};
    small.data.assign(27, 1.0f);
    small.voxel = {1, 1, 1};
    const std::vector<double> w = io::atomWeights(m.elements);
    for (const Value& c : list(at(cases(), "errors"))) {
        const std::string label = str(list(c).at(0));
        const Value& want = list(c).at(1);
        SCOPED_TRACE(label);
        std::optional<FitResult> result;
        std::string error;
        try {
            if (label == "step below 3") {
                (void)rotationSamples(2.0);
            } else if (label == "no resolution") {
                FitProblem p(small, m.xyz, w, 0.0);
            } else if (label == "negative resolution") {
                FitProblem p(small, m.xyz, w, -1.0);
            } else if (label == "two atoms") {
                FitProblem p(small, std::span(m.xyz).first(2), std::span(w).first(2), 6.0);
            } else if (label == "no atoms to simulate") {
                (void)io::simulateMap({}, {}, 6.0);
            } else if (label == "tiny map global") {
                result = fitModel(small, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Global, 45.0);
            } else if (label == "tiny map manual") {
                result = fitModel(small, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Manual);
            } else {
                FAIL() << "unknown case";
            }
        } catch (const FitError& e) {
            error = e.what();
        }
        if (has(want, "error")) {
            EXPECT_EQ(error, errorMessage(at(want, "error")));
        } else if (has(at(want, "value"), "cc")) {
            ASSERT_TRUE(result) << error;
            const Value& v = at(want, "value");
            EXPECT_NEAR(result->cc, hexDouble(at(v, "cc")), kCcTolerance);
            EXPECT_EQ(result->message, str(at(v, "message")));
            EXPECT_EQ(result->start, str(at(v, "start")));
        } else {
            EXPECT_TRUE(error.empty()) << error;
        }
    }
}

// CH-10: a non-finite resolution is refused before any work
TEST(DensityFit, NonFiniteResolutionIsRefused) {
    const Model& m = model("helix");
    const io::DensityMap map = io::readMap(kData / "chem" / "maps" / "fit_helix.mrc");
    const std::vector<double> w = io::atomWeights(m.elements);
    for (const double res : {std::nan(""), HUGE_VAL, -HUGE_VAL}) {
        try {
            FitProblem p(map, m.xyz, w, res);
            ADD_FAILURE() << res;
        } catch (const FitError& e) {
            EXPECT_EQ(std::string(e.what()), "the resolution must be a positive number of Å");
        }
    }
}

// Every problem case: the FitProblem, scores, flips, the coarse search, refinements and fits
class DensityFitProblem : public ::testing::TestWithParam<std::size_t> {};

TEST_P(DensityFitProblem, AsMdatk) {
    const Value& c = list(at(cases(), "problems")).at(GetParam());
    const Model& m = model(str(at(c, "model")));
    const io::DensityMap map = io::readMap(kData / "chem" / str(at(c, "map")));
    const double res = hexDouble(at(c, "resolution"));
    SCOPED_TRACE(str(at(c, "map")));
    ASSERT_FALSE(has(c, "error"));
    const FitProblem prob(map, m.xyz, io::atomWeights(m.elements), res);
    const Value& p = at(c, "problem");
    EXPECT_LT(maxAbs(prob.centre(), vec(at(p, "centre"))), 1e-12);
    EXPECT_NEAR(prob.radius(), hexDouble(at(p, "radius")), 1e-12);
    EXPECT_NEAR(prob.rg(), hexDouble(at(p, "rg")), 1e-12);
    ASSERT_EQ(prob.points().size(), count(at(p, "points")));
    for (std::size_t i = 0; i < list(at(p, "first_points")).size(); ++i) {
        EXPECT_LT(maxAbs(prob.points()[i], vec(list(at(p, "first_points"))[i])), 1e-12);
        EXPECT_EQ(prob.values()[i], hexDouble(list(at(p, "first_values"))[i]));
    }
    double total = 0.0;
    for (const double v : prob.values()) {
        total += v;
    }
    EXPECT_NEAR(total, hexDouble(at(p, "sum")), 1e-9 * std::abs(total));
    EXPECT_EQ(prob.nonFiniteVoxels(), str(at(c, "map")) == "maps/fit_nan.mrc" ? 1u : 0u);

    // scores of given poses, at the resolution and against the 8 A copy
    std::vector<RigidTransform> poses;
    const std::optional<FitProblem> coarse8 =
        res < kCoarseResolution ? std::optional<FitProblem>(prob.atResolution(kCoarseResolution)) : std::nullopt;
    const FitProblem& at8 = coarse8 ? *coarse8 : prob;
    for (std::size_t i = 0; i < list(at(c, "cc")).size(); ++i) {
        const Array& row = list(list(at(c, "cc"))[i]);
        const RigidTransform t = transform(row.at(0));
        poses.push_back(t);
        EXPECT_NEAR(prob.cc(t), hexDouble(row.at(1)), kScoreTolerance) << "pose " << i;
        EXPECT_NEAR(prob.insideFraction(t), hexDouble(row.at(2)), 1e-12) << "pose " << i;
        EXPECT_NEAR(at8.cc(t), hexDouble(list(at(c, "cc_at_8"))[i]), kScoreTolerance) << "pose " << i;
    }
    const std::vector<RigidTransform> flips = prob.flips();
    for (std::size_t k = 0; k < 3; ++k) {
        const RigidTransform want = transform(list(at(c, "flips"))[k]);
        EXPECT_LT(maxAbs(flips[k].rotation, want.rotation), 1e-9) << "flip " << k;
        EXPECT_LT(maxAbs(flips[k].translation, want.translation), 1e-8) << "flip " << k;
    }

    // the coarse search: grid and the peaks of three rotations
    const Value& cs = at(c, "coarse");
    const double step = str(at(c, "model")) == "strand" ? 45.0 : 30.0;
    const CoarseSearch search(prob, step);
    EXPECT_EQ(search.spacing(), hexDouble(at(cs, "h")));
    EXPECT_LT(maxAbs(search.lowCorner(), vec(at(cs, "lo"))), 1e-12);
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_EQ(search.gridXyz()[k], count(list(at(cs, "n"))[k]));
    }
    EXPECT_NEAR(search.varFloor(), hexDouble(at(cs, "var_floor")), 1e-5 * search.varFloor());
    const std::vector<Mat3> rots = rotationSamples(step);
    for (const Value& entry : list(at(cs, "peaks"))) {
        const std::size_t k = count(list(entry).at(0));
        const Array& want = list(list(entry).at(1));
        const std::vector<CoarseSearch::Peak> got = search.peaks(rots[k], 2);
        ASSERT_EQ(got.size(), want.size()) << "rotation " << k;
        for (std::size_t j = 0; j < got.size(); ++j) {
            EXPECT_NEAR(got[j].score, hexDouble(list(want[j]).at(0)), kScoreTolerance) << "rotation " << k;
            EXPECT_LT(maxAbs(got[j].centre, vec(list(want[j]).at(1))), 1e-9) << "rotation " << k;
        }
    }

    // local refinements from the recorded starts
    for (const Value& r : list(at(c, "refine"))) {
        const auto [fit, value] = prob.refine(transform(at(r, "start")));
        expectPose(fit, transform(at(r, "transform")), prob.centre(), "refine");
        EXPECT_NEAR(value, hexDouble(at(r, "cc")), kCcTolerance);
    }

    // fit_model: global, manual from the first start, from the identity
    for (const Value& f : list(at(c, "fit_model"))) {
        const bool global = str(at(f, "search")) == "global";
        SCOPED_TRACE(str(at(f, "search")));
        ASSERT_FALSE(has(f, "error"));
        const std::optional<RigidTransform> start =
            has(f, "start") ? std::optional<RigidTransform>(transform(at(f, "start"))) : std::nullopt;
        const FitResult got = fitModel(map, m.xyz, m.elements, res, start,
                                       global ? FitSearch::Global : FitSearch::Manual, step);
        expectPose(got.transform, transform(at(f, "transform")), prob.centre(), "fit_model");
        EXPECT_NEAR(got.cc, hexDouble(at(f, "cc")), kCcTolerance);
        EXPECT_EQ(got.start, str(at(f, "start_kind")));
        EXPECT_EQ(got.candidates, count(at(f, "candidates")));
        EXPECT_EQ(got.message, str(at(f, "message")));
        // Powell's path follows every score: the count agrees only roughly
        const auto want = static_cast<double>(count(at(f, "evaluations")));
        EXPECT_NEAR(static_cast<double>(got.evaluations), want, 0.5 * want);
    }
}

INSTANTIATE_TEST_SUITE_P(Cases, DensityFitProblem,
                         ::testing::Range<std::size_t>(0, list(at(cases(), "problems")).size()));

// The global search gives the same result for any thread count; progress in MDATK's steps
TEST(DensityFit, SameResultForAnyThreadCountAndProgress) {
    const Model& m = model("helix");
    const io::DensityMap map = io::readMap(kData / "chem" / "maps" / "fit_helix.mrc");
    std::vector<std::pair<double, std::string>> events;
    const core::Reporter reporter({}, [&](double fraction, const std::string& message) {
        events.emplace_back(fraction, message);
    });
    const FitResult one = fitModel(map, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Global, 30.0, {&reporter, 1});
    const FitResult three = fitModel(map, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Global, 30.0, {nullptr, 3});
    const FitResult many = fitModel(map, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Global, 30.0, {nullptr, 8});
    for (const FitResult* r : {&three, &many}) {
        EXPECT_EQ(r->transform.rotation, one.transform.rotation);
        EXPECT_EQ(r->transform.translation, one.transform.translation);
        EXPECT_EQ(r->cc, one.cc);
        EXPECT_EQ(r->evaluations, one.evaluations);
        EXPECT_EQ(r->message, one.message);
    }
    ASSERT_FALSE(events.empty());
    EXPECT_EQ(events.front().second, "global search: rotation 1/551");
    EXPECT_EQ(events[1].second, "global search: rotation 17/551");
    EXPECT_EQ(events.back().second, "global search done");
    EXPECT_EQ(events.back().first, 1.0);
    std::size_t rotations = 0, poses = 0, finals = 0;
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (i > 0) {
            EXPECT_GE(events[i].first, events[i - 1].first) << events[i].second;
        }
        if (events[i].second.rfind("global search: rotation", 0) == 0) {
            ++rotations;
        } else if (events[i].second.rfind("refining pose", 0) == 0) {
            ++poses;
        } else if (events[i].second == "final refinement") {
            ++finals;
        }
    }
    EXPECT_EQ(rotations, 36u);  // 1, 17, ..., 545 and 551
    EXPECT_EQ(poses, 12u);
    EXPECT_EQ(finals, 9u);  // the best 3 and the 3 flips of the best 2
}

TEST(DensityFit, CancelStopsTheSearch) {
    const Model& m = model("helix");
    const io::DensityMap map = io::readMap(kData / "chem" / "maps" / "fit_helix.mrc");
    auto flag = std::make_shared<std::atomic<bool>>(false);
    std::size_t calls = 0;
    const core::Reporter reporter(
        {},
        [&](double, const std::string&) {
            if (++calls == 3) {
                flag->store(true);
            }
        },
        flag);
    EXPECT_THROW(fitModel(map, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Global, 30.0, {&reporter, 2}),
                 core::Cancelled);
    EXPECT_EQ(calls, 3u);
    reporter.cancel();
    EXPECT_THROW(fitModel(map, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Manual, 30.0, {&reporter, 1}),
                 core::Cancelled);
}

// CH-11: a NaN voxel reads as no density: the fit equals the fit of the map with a 0 there, plus a note
TEST(DensityFit, NonFiniteVoxelsCarryNoDensity) {
    const Model& m = model("helix");
    const io::DensityMap nan = io::readMap(kData / "chem" / "maps" / "fit_nan.mrc");
    io::DensityMap zero = nan;
    std::size_t bad = 0;
    for (float& v : zero.data) {
        if (!std::isfinite(v)) {
            v = 0.0f;
            ++bad;
        }
    }
    ASSERT_EQ(bad, 1u);
    const Value& nanCase = list(at(cases(), "problems")).at(3);
    ASSERT_EQ(str(at(nanCase, "map")), "maps/fit_nan.mrc");
    const RigidTransform start = transform(at(list(at(nanCase, "refine")).at(0), "start"));
    const FitResult a = fitModel(nan, m.xyz, m.elements, 6.0, start, FitSearch::Manual);
    const FitResult b = fitModel(zero, m.xyz, m.elements, 6.0, start, FitSearch::Manual);
    EXPECT_EQ(a.transform.rotation, b.transform.rotation);
    EXPECT_EQ(a.cc, b.cc);
    EXPECT_EQ(a.message, b.message + "; 1 non-finite map voxel(s) read as no density");
    // the density centre of a local start uses the finite voxels
    const RigidTransform centred = centreTransform(nan, io::mapLevel(zero, 1.0), {0, 0, 0});
    EXPECT_TRUE(std::isfinite(centred.translation[0]) && std::isfinite(centred.translation[1]) &&
                std::isfinite(centred.translation[2]));
}

TEST(DensityFit, CentreTransformPutsTheModelOnTheDensity) {
    io::DensityMap map;
    map.shapeXyz = {5, 4, 3};
    map.data.assign(60, 0.0f);
    map.voxel = {2, 2, 2};
    map.origin = {10, 20, 30};
    map.data[(1 * 4 + 2) * 5 + 3] = 5.0f;  // x 3, y 2, z 1
    const RigidTransform t = centreTransform(map, 1.0, {1, 1, 1});
    EXPECT_EQ(t.rotation, RigidTransform::identity().rotation);
    EXPECT_EQ(t.translation, (Vec3{10 + 6 - 1, 20 + 4 - 1, 30 + 2 - 1}));
}

// tests/test_density.py: test_rotation_samples_cover_so3
TEST(DensityFit, RotationSamplesCoverSo3) {
    std::mt19937_64 rng(1);
    std::normal_distribution<double> normal;
    for (const double step : {20.0, 30.0}) {
        const std::vector<Mat3> rots = rotationSamples(step);
        for (const Mat3& r : rots) {
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    const double d = r[i][0] * r[j][0] + r[i][1] * r[j][1] + r[i][2] * r[j][2];
                    ASSERT_NEAR(d, i == j ? 1.0 : 0.0, 1e-9);
                }
            }
        }
        double worst = 0.0;
        for (int n = 0; n < 200; ++n) {
            const Mat3 r = core::quaternionMatrix({normal(rng), normal(rng), normal(rng), normal(rng)});
            double best = 180.0;
            for (const Mat3& s : rots) {
                best = std::min(best, angleBetween(s, r));
            }
            worst = std::max(worst, best);
        }
        EXPECT_LT(worst, 1.05 * step) << step;
    }
}

namespace {

// tests/test_density.py's _target_map: the model simulated at the true pose, padded, with noise (here
// from std::mt19937_64: a property test, not a parity case)
io::DensityMap targetMap(std::span<const Vec3> xyz, std::span<const double> weights, double res, double spacing,
                         std::size_t pad, double noise, std::uint64_t seed) {
    const io::DensityMap sim = io::simulateMap(xyz, weights, res, spacing);
    io::DensityMap out;
    for (std::size_t k = 0; k < 3; ++k) {
        out.shapeXyz[k] = sim.shapeXyz[k] + static_cast<std::int64_t>(2 * pad);
        out.origin[k] = sim.origin[k] - static_cast<double>(pad) * spacing;
    }
    out.voxel = sim.voxel;
    const auto nx = static_cast<std::size_t>(out.shapeXyz[0]), ny = static_cast<std::size_t>(out.shapeXyz[1]),
               nz = static_cast<std::size_t>(out.shapeXyz[2]);
    out.data.assign(nx * ny * nz, 0.0f);
    float vmax = 0.0f;
    for (const float v : sim.data) {
        vmax = std::max(vmax, v);
    }
    const auto sx = static_cast<std::size_t>(sim.shapeXyz[0]), sy = static_cast<std::size_t>(sim.shapeXyz[1]),
               sz = static_cast<std::size_t>(sim.shapeXyz[2]);
    for (std::size_t z = 0; z < sz; ++z) {
        for (std::size_t y = 0; y < sy; ++y) {
            for (std::size_t x = 0; x < sx; ++x) {
                out.data[((z + pad) * ny + y + pad) * nx + x + pad] = sim.data[(z * sy + y) * sx + x];
            }
        }
    }
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, noise * static_cast<double>(vmax));
    for (float& v : out.data) {
        v += static_cast<float>(normal(rng));
    }
    return out;
}

}  // namespace

// tests/test_density.py: test_global_search_recovers_the_pose (the demo model, pad 15, 1.5 A voxels)
TEST(DensityFit, GlobalSearchRecoversThePose) {
    const Model& m = model("demo");
    RigidTransform truth;
    truth.rotation = core::quaternionMatrix({0.3, -0.8, 0.4, 0.2});
    truth.translation = {25.0, -40.0, 12.0};
    const std::vector<double> w = io::atomWeights(m.elements);
    const io::DensityMap map = targetMap(truth.apply(m.xyz), w, 6.0, 1.5, 15, 0.1, 3);
    const FitResult res = fitModel(map, m.xyz, m.elements, 6.0, std::nullopt, FitSearch::Global, 30.0);
    const Vec3 c = meanOf(m.xyz);
    EXPECT_GT(res.cc, 0.9);
    EXPECT_LT(angleBetween(res.transform.rotation, truth.rotation), 3.0);
    EXPECT_LT(distance(res.transform.apply(c), truth.apply(c)), 1.0);
    // local refinement from a nearby placement converges to the same pose
    RigidTransform nearby;
    nearby.rotation = core::rotvecMatrix({0.15, -0.1, 0.1});
    nearby.translation = {2.0, -1.5, 1.0};
    const RigidTransform start = nearby.about(truth.apply(c)).compose(truth);
    const FitResult loc = fitModel(map, m.xyz, m.elements, 6.0, start, FitSearch::Manual);
    EXPECT_EQ(loc.start, "manual");
    EXPECT_GT(loc.cc, 0.9);
    EXPECT_LT(angleBetween(loc.transform.rotation, truth.rotation), 3.0);
}

// tests/test_density.py: test_global_search_ignores_unrelated_density_and_validates_steps
TEST(DensityFit, GlobalSearchIgnoresUnrelatedDensityAndValidatesSteps) {
    const Model& m = model("demo");
    RigidTransform truth;
    truth.rotation = core::quaternionMatrix({0.6, -0.2, 0.7, 0.3});
    truth.translation = {60.0, 50.0, 55.0};
    std::vector<Vec3> pts = truth.apply(m.xyz);
    const Vec3 placedCentre = meanOf(pts);
    std::vector<double> w = io::atomWeights(m.elements);
    std::mt19937_64 rng(5);
    std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (int i = 0; i < 2500; ++i) {  // a protein-like blob 42 A along x
        const Vec3 d{normal(rng), normal(rng), normal(rng)};
        const double len = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        const double r = 22.0 * std::cbrt(uniform(rng));
        pts.push_back({d[0] / len * r + placedCentre[0] + 42.0, d[1] / len * r + placedCentre[1],
                       d[2] / len * r + placedCentre[2]});
        w.push_back(6.7);
    }
    io::DensityMap map = io::simulateMap(pts, w, 8.0, 1.6);
    float vmax = 0.0f;
    for (const float v : map.data) {
        vmax = std::max(vmax, v);
    }
    std::normal_distribution<double> noise(0.0, 0.05 * static_cast<double>(vmax));
    for (float& v : map.data) {
        v += static_cast<float>(noise(rng));
    }
    const FitResult res = fitModel(map, m.xyz, m.elements, 8.0, std::nullopt, FitSearch::Global, 30.0);
    const Vec3 c = meanOf(m.xyz);
    EXPECT_LT(distance(res.transform.apply(c), truth.apply(c)), 1.5);
    EXPECT_LT(angleBetween(res.transform.rotation, truth.rotation), 5.0);
    EXPECT_THROW((void)rotationSamples(0), FitError);
    // tiny maps do not crash
    io::DensityMap tiny;
    tiny.shapeXyz = {30, 30, 2};
    tiny.voxel = {1, 1, 1};
    for (int i = 0; i < 1800; ++i) {
        tiny.data.push_back(static_cast<float>(uniform(rng)));
    }
    try {
        (void)fitModel(tiny, m.xyz, m.elements, 20.0, std::nullopt, FitSearch::Global);
    } catch (const FitError&) {
    }
}
