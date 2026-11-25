// The optimal-superposition RMSD (src/analysis/rmsd_qcp.cpp) against MDAnalysis 2.11.0's
// rms.rmsd(a, b, weights=w, center=True, superposition=True) (docs/specs/m3-rmsd.md), on the
// probes of data/analysis/rmsd_cases.json: every case and the demo series compared exactly (NaN
// equals NaN). The spec's fallback tolerance |r_MDA^2 - r^2| <= 1024 eps E0/n is accepted only for
// non-degenerate cases and reported.
//
// Improvements (docs/DIFFERENCES.md, AN-16 / AN-17): where the expectation is marked "improved"
// (degenerate sets, and E0 - lambda lost to cancellation), it is the exact superposition RMSD
// (numpy's SVD, residual after rotating), which MDBridge reaches with a Jacobi eigen-solution:
// compared within kExactTolerance * sqrt(E0 / n) (far below float32 coordinate precision), and exactly
// where the exact value is 0 (identical inputs, one atom). "near" AN-17 decisions also accept
// MDAnalysis' value bit for bit.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/analysis/rmsd_qcp.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/io/trajectory.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using analysis::ShapedArray;
using analysis::superposedRmsd;
using core::Array;
using core::Value;
using test::kData;

namespace {

const Value& data() {
    static const Value value =
        core::parseJson(test::readFile(kData / "analysis" / "rmsd_cases.json"), "rmsd_cases.json");
    return value;
}

const Value& at(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
const Array& list(const Value& value) {
    return *value.get<Array>();
}
std::string str(const Value& value) {
    return *value.get<std::string>();
}

// float.hex() text, "nan", "inf", "-inf"
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
std::vector<double> doubles(const Value& value) {
    std::vector<double> out;
    for (const Value& v : list(value)) {
        out.push_back(hexDouble(v));
    }
    return out;
}
std::vector<std::size_t> shape(const Value& value) {
    std::vector<std::size_t> out;
    for (const Value& v : list(value)) {
        out.push_back(static_cast<std::size_t>(*v.get<std::int64_t>()));
    }
    return out;
}

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

// within the spec's fallback tolerance (§7)
bool withinTolerance(double expected, double actual, double scale) {
    const double eps = std::ldexp(1.0, -52);
    return std::isfinite(expected) && std::isfinite(actual) &&
           std::fabs(expected * expected - actual * actual) <= 1024.0 * eps * scale;
}

struct Tally {
    int exact = 0;
    int tolerance = 0;
    int mismatch = 0;
    int improved = 0;
};

// The exact value of a nearly collinear set is ill-conditioned: numpy's SVD itself moves by 1e-12
// of the coordinates' size between equivalent inputs (a and b swapped)
constexpr double kExactTolerance = 1e-10;

// An improved expectation (exact superposition RMSD): exact for 0, else within kExactTolerance of
// the coordinates' size; `mdanalysis`: the QCP value also accepted (near the AN-17 threshold).
void compareImproved(Tally& tally, const std::string& name, double expected, double actual, double scale,
                     std::optional<double> mdanalysis) {
    ++tally.improved;
    const bool ok = expected == 0.0 ? actual == 0.0
                                    : std::fabs(actual - expected) <= kExactTolerance * std::sqrt(scale) + 1e-300;
    if (ok || (mdanalysis && sameDouble(*mdanalysis, actual))) {
        ++tally.exact;
        return;
    }
    ++tally.mismatch;
    ADD_FAILURE() << name << ": expected the exact value " << std::hexfloat << expected << ", got " << actual
                  << std::defaultfloat << " (difference " << actual - expected << ", scale " << scale << ")";
}

// E0 / n of centred a and b (unweighted): the size the improved values are compared at
double meanSquare(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = a.size() / 3;
    double total = 0.0;
    for (const std::vector<double>* x : {&a, &b}) {
        for (std::size_t k = 0; k < 3; ++k) {
            double centre = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                centre += (*x)[3 * i + k];
            }
            centre /= static_cast<double>(n);
            for (std::size_t i = 0; i < n; ++i) {
                const double d = (*x)[3 * i + k] - centre;
                total += d * d;
            }
        }
    }
    return total / static_cast<double>(2 * n);
}

void compare(Tally& tally, const std::string& name, double expected, double actual, bool degenerate,
             std::optional<double> scale) {
    if (sameDouble(expected, actual)) {
        ++tally.exact;
    } else if (!degenerate && scale && withinTolerance(expected, actual, *scale)) {
        ++tally.tolerance;
        std::printf("within tolerance only: %s expected %a got %a\n", name.c_str(), expected, actual);
    } else {
        ++tally.mismatch;
        ADD_FAILURE() << name << ": expected " << expected << " (" << std::hexfloat << expected << "), got "
                      << actual << " (" << actual << ")" << std::defaultfloat
                      << (degenerate ? " [degenerate]" : "");
    }
}

}  // namespace

TEST(RmsdQcp, ProbeCases) {
    Tally tally;
    int errors = 0;
    int degenerateExact = 0;
    for (const Value& c : list(at(data(), "cases"))) {
        const std::string name = str(at(c, "name"));
        const std::vector<double> a = doubles(at(c, "a"));
        const std::vector<double> b = doubles(at(c, "b"));
        const std::vector<std::size_t> aShape = shape(at(c, "a_shape"));
        const std::vector<std::size_t> bShape = shape(at(c, "b_shape"));
        std::vector<double> w;
        std::vector<std::size_t> wShape;
        std::optional<ShapedArray> weights;
        if (at(c, "weights").get<Array>()) {
            w = doubles(at(c, "weights"));
            wShape = shape(at(c, "weights_shape"));
            weights = ShapedArray{w, wShape};
        }
        const ShapedArray sa{a, aShape};
        const ShapedArray sb{b, bShape};
        const Value* error = c.find("error");
        if (error && error->find("type")) {
            const std::string type = str(at(*error, "type"));
            const std::string message = str(at(*error, "message"));
            try {
                superposedRmsd(sa, sb, weights);
                ADD_FAILURE() << name << ": no error, expected " << type;
            } catch (const analysis::ZeroDivisionError& e) {
                EXPECT_EQ(type, "ZeroDivisionError") << name;
                EXPECT_EQ(e.pythonName(), "ZeroDivisionError") << name;
                EXPECT_EQ(std::string(e.what()), message) << name;
            } catch (const core::ValueError& e) {
                EXPECT_EQ(type, "ValueError") << name;
                EXPECT_STREQ(e.pythonName(), "ValueError") << name;
                EXPECT_EQ(std::string(e.what()), message) << name;
            }
            ++errors;
            continue;
        }
        const bool degenerate = *at(c, "degenerate").get<bool>();
        std::optional<double> scale;
        if (at(c, "scale").get<std::string>()) {
            scale = hexDouble(at(c, "scale"));
        }
        const double expected = hexDouble(at(c, "result"));
        const double actual = superposedRmsd(sa, sb, weights);
        if (c.find("improved")) {
            std::optional<double> mdanalysis;
            if (*at(c, "near").get<bool>()) {
                mdanalysis = hexDouble(at(c, "mdanalysis"));
            }
            compareImproved(tally, name, expected, actual, std::max(*scale, meanSquare(a, b)), mdanalysis);
            continue;
        }
        const int before = tally.exact;
        compare(tally, name, expected, actual, degenerate, scale);
        if (degenerate && tally.exact > before) {
            ++degenerateExact;
        }
    }
    std::printf("probe cases: %d exact or improved (%d improved, %d degenerate exact), %d within tolerance only, "
                "%d mismatches, %d errors\n",
                tally.exact, tally.improved, degenerateExact, tally.tolerance, tally.mismatch, errors);
    EXPECT_EQ(tally.exact + tally.tolerance + tally.mismatch, 233);
    EXPECT_EQ(tally.improved, 136);
    EXPECT_EQ(errors, 10);
    EXPECT_EQ(tally.mismatch, 0);
}

TEST(RmsdQcp, DemoSeries) {
    // MDATK's loop: the reference frame's positions, then rmsd(frame, reference) for every frame
    const Value& demo = at(data(), "demo");
    auto trajectory = io::openTrajectory(kData / str(at(demo, "trajectory")), "", std::nullopt, 0);
    const auto frames = *at(demo, "n_frames").get<std::int64_t>();
    ASSERT_EQ(trajectory->frameCount(), frames);
    Tally tally;
    int series = 0;
    for (const Value& target : list(at(demo, "targets"))) {
        std::vector<std::size_t> atoms;
        for (const Value& run : list(at(target, "atoms"))) {
            const auto first = *list(run).at(0).get<std::int64_t>();
            const auto last = *list(run).at(1).get<std::int64_t>();
            for (auto i = first; i <= last; ++i) {
                atoms.push_back(static_cast<std::size_t>(i));
            }
        }
        std::optional<std::vector<double>> weights;
        if (at(target, "weights").get<Array>()) {
            weights = doubles(at(target, "weights"));
            ASSERT_EQ(weights->size(), atoms.size());
        }
        const auto positions = [&](std::int64_t frame) {
            const core::Frame& f = trajectory->seek(frame);
            std::vector<double> out;
            for (std::size_t i : atoms) {
                for (std::size_t k = 0; k < 3; ++k) {
                    out.push_back(static_cast<double>(f.positions.at(3 * i + k)));
                }
            }
            return out;
        };
        const std::vector<double> reference = positions(*at(target, "reference_frame").get<std::int64_t>());
        const std::vector<double> expected = doubles(at(target, "values"));
        ASSERT_EQ(static_cast<std::int64_t>(expected.size()), frames);
        const std::string label = str(at(target, "label")) + " ref " +
                                  std::to_string(*at(target, "reference_frame").get<std::int64_t>()) +
                                  (weights ? " mass" : "");
        // one or two atoms are degenerate (no tolerance); larger sets of the demo are not collinear
        const bool degenerate = atoms.size() <= 2;
        std::set<std::int64_t> improvedFrames, nearFrames;
        for (const Value& f : list(at(target, "improved"))) {
            improvedFrames.insert(*f.get<std::int64_t>());
        }
        for (const Value& f : list(at(target, "near"))) {
            nearFrames.insert(*f.get<std::int64_t>());
        }
        const std::vector<double> mdanalysis = doubles(at(target, "mdanalysis"));
        for (std::int64_t frame = 0; frame < frames; ++frame) {
            const std::vector<double> current = positions(frame);
            const double actual =
                weights ? superposedRmsd(current, reference, std::span<const double>(*weights))
                        : superposedRmsd(current, reference);
            const auto f = static_cast<std::size_t>(frame);
            const std::string where = label + " frame " + std::to_string(frame);
            if (improvedFrames.contains(frame)) {
                compareImproved(tally, where, expected[f], actual, meanSquare(current, reference),
                                nearFrames.contains(frame) ? std::optional<double>(mdanalysis[f]) : std::nullopt);
                continue;
            }
            compare(tally, where, expected[f], actual, degenerate, std::nullopt);
        }
        ++series;
    }
    std::printf("demo: %d series, %d values exact, %d within tolerance only, %d mismatches\n", series,
                tally.exact, tally.tolerance, tally.mismatch);
    EXPECT_EQ(series, 32);
    EXPECT_EQ(tally.mismatch, 0);
}

TEST(RmsdQcp, ErrorsInOrder) {
    const std::vector<double> five(15, 1.0);
    const std::vector<double> four(12, 1.0);
    const std::vector<double> zeros(5, 0.0);
    // the shape check comes before the weight checks
    EXPECT_THROW(superposedRmsd(five, four, std::span<const double>(zeros)), core::ValueError);
    EXPECT_THROW(superposedRmsd(five, five, std::span<const double>(zeros).first(4)), core::ValueError);
    EXPECT_THROW(superposedRmsd(five, five, std::span<const double>(zeros)), analysis::ZeroDivisionError);
    // not n x 3
    EXPECT_THROW(superposedRmsd(std::vector<double>(4, 0.0), std::vector<double>(4, 0.0)), std::invalid_argument);
}

TEST(RmsdQcp, SimpleValues) {
    // a square against itself and against its mirror image (a proper rotation by pi superposes it)
    const std::vector<double> square{0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<double> moved{5, 5, 5, 5, 6, 5, 4, 6, 5, 4, 5, 5};  // rotated by 90 degrees, shifted
    EXPECT_LT(superposedRmsd(square, moved), 1e-7);
    // stretched by 2 along x: each centred atom is 0.5 off in x, the best rotation is the identity
    const std::vector<double> stretched{0, 0, 0, 2, 0, 0, 2, 1, 0, 0, 1, 0};
    EXPECT_NEAR(superposedRmsd(square, stretched), 0.5, 1e-12);
}
