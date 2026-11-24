// Rigid transforms and rotation helpers against what MDATK computes (tests/parity/
// make_rigid_cases.py, data/rigid, see data/README.md), and MDATK's own transform test.

#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

const Value& cases() {
    static const Value data =
        parseJson(readFile(kData / "rigid" / "rigid_cases.json"), "rigid_cases.json");
    return data;
}

const Array& list(const Value& value) {
    return *value.get<Array>();
}
const Array& items(const Value& value, std::string_view key) {
    return list(*value.find(key));
}

double number(const Value& value) {
    if (const double* d = value.get<double>()) {
        return *d;
    }
    return static_cast<double>(*value.get<std::int64_t>());
}

std::vector<double> numbers(const Value& value) {
    std::vector<double> out;
    for (const Value& v : list(value)) {
        out.push_back(number(v));
    }
    return out;
}

Vec3 vec3(const Value& value) {
    const std::vector<double> v = numbers(value);
    return {v.at(0), v.at(1), v.at(2)};
}

Mat3 mat3(const Value& value) {
    const std::vector<double> v = numbers(value);
    Mat3 m{};
    for (std::size_t i = 0; i < 9; ++i) {
        m[i / 3][i % 3] = v.at(i);
    }
    return m;
}

RigidTransform transform(const Value& value) {
    return {mat3(*value.find("rotation")), vec3(*value.find("translation"))};
}

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

// numpy multiplies with BLAS (fused, other order) and decomposes with LAPACK
void expectNear(double got, double expected, double tolerance, const std::string& context) {
    if (std::isnan(expected)) {
        EXPECT_TRUE(std::isnan(got)) << context;
        return;
    }
    EXPECT_LE(std::abs(got - expected), tolerance * std::max(1.0, std::abs(expected)))
        << context << ": " << got << " vs " << expected;
}

void expectMat(const Mat3& got, const Value& expected, double tolerance,
               const std::string& context) {
    const Mat3 e = mat3(expected);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            expectNear(got[i][j], e[i][j], tolerance,
                       context + " [" + std::to_string(i) + "][" + std::to_string(j) + "]");
        }
    }
}

void expectVec(const Vec3& got, const Value& expected, double tolerance,
               const std::string& context) {
    const Vec3 e = vec3(expected);
    for (std::size_t k = 0; k < 3; ++k) {
        expectNear(got[k], e[k], tolerance, context + " [" + std::to_string(k) + "]");
    }
}

}  // namespace

TEST(RigidTransform, FromListMatchesMdatk) {
    int index = 0;
    for (const Value& c : items(cases(), "from_list")) {
        const std::string context = "case " + std::to_string(index++);
        const std::vector<double> input = numbers(*c.find("input"));
        if (const Value* error = c.find("error")) {
            try {
                RigidTransform::fromList(input);
                ADD_FAILURE() << context << ": no error, expected " << error->stringOr();
            } catch (const ValueError& e) {
                EXPECT_EQ(e.what(), error->stringOr()) << context;
            }
            continue;
        }
        const RigidTransform t = RigidTransform::fromList(input);
        // numpy's LAPACK factor is itself up to 3e-15 from the exact one for coarse inputs
        expectMat(t.rotation, *c.find("rotation"), 4e-15, context);
        const Vec3 translation = vec3(*c.find("translation"));
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_TRUE(sameDouble(t.translation[k], translation[k])) << context;
        }
    }
    EXPECT_GT(index, 100);
}

TEST(RigidTransform, ToListMatchesMdatk) {
    for (const Value& c : items(cases(), "to_list")) {
        const RigidTransform t = transform(c);
        const int digits = static_cast<int>(*c.find("digits")->get<std::int64_t>());
        const std::vector<double> got = t.toList(digits);
        const std::vector<double> expected = numbers(*c.find("list"));
        ASSERT_EQ(got.size(), expected.size());
        for (std::size_t k = 0; k < got.size(); ++k) {
            EXPECT_TRUE(sameDouble(got[k], expected[k]))
                << "digits " << digits << " [" << k << "] " << got[k] << " vs " << expected[k];
        }
    }
}

TEST(RigidTransform, AlgebraMatchesMdatk) {
    int index = 0;
    for (const Value& c : items(cases(), "algebra")) {
        const std::string context = "case " + std::to_string(index++);
        const RigidTransform a = transform(*c.find("a"));
        const RigidTransform b = transform(*c.find("b"));
        std::vector<Vec3> points;
        for (const Value& p : items(c, "points")) {
            points.push_back(vec3(p));
        }
        const std::vector<Vec3> applied = a.apply(points);
        const Array& expected = items(c, "apply");
        ASSERT_EQ(applied.size(), expected.size());
        for (std::size_t k = 0; k < applied.size(); ++k) {
            expectVec(applied[k], expected[k], 4e-15, context + " apply");
            expectVec(a.apply(points[k]), expected[k], 4e-15, context + " apply one");
        }
        const RigidTransform inv = a.inverse();
        expectMat(inv.rotation, *c.find("inverse")->find("rotation"), 0, context + " inverse");
        expectVec(inv.translation, *c.find("inverse")->find("translation"), 4e-15,
                  context + " inverse");
        const RigidTransform comp = a.compose(b);
        expectMat(comp.rotation, *c.find("compose")->find("rotation"), 4e-16, context + " compose");
        expectVec(comp.translation, *c.find("compose")->find("translation"), 4e-15,
                  context + " compose");
        const RigidTransform about = a.about(vec3(*c.find("centre")));
        expectMat(about.rotation, *c.find("about")->find("rotation"), 0, context + " about");
        expectVec(about.translation, *c.find("about")->find("translation"), 4e-15,
                  context + " about");
    }
}

TEST(RigidTransform, RotationHelpersMatchMdatk) {
    for (const Value& c : items(cases(), "rotvec")) {
        expectMat(rotvecMatrix(vec3(*c.find("w"))), *c.find("matrix"), 1e-15, "rotvec");
    }
    for (const Value& c : items(cases(), "quaternion")) {
        const std::vector<double> q = numbers(*c.find("q"));
        expectMat(quaternionMatrix({q.at(0), q.at(1), q.at(2), q.at(3)}), *c.find("matrix"), 1e-15,
                  "quaternion");
    }
    for (const Value& c : items(cases(), "angle")) {
        expectNear(rotationAngle(mat3(*c.find("matrix"))), number(*c.find("degrees")), 1e-13,
                   "angle");
    }
    for (const Value& c : items(cases(), "anchor_points")) {
        std::vector<Vec3> points;
        for (const Value& p : items(c, "points")) {
            points.push_back(vec3(p));
        }
        const std::vector<Vec3> got = anchorPoints(points);
        const Array& expected = items(c, "anchor");
        ASSERT_EQ(got.size(), expected.size());
        for (std::size_t k = 0; k < got.size(); ++k) {
            const Vec3 e = vec3(expected[k]);
            for (std::size_t j = 0; j < 3; ++j) {
                EXPECT_TRUE(sameDouble(got[k][j], e[j])) << got[k][j] << " vs " << e[j];
            }
        }
    }
}

TEST(RigidTransform, IdentityAndMessages) {
    const RigidTransform id = RigidTransform::identity();
    EXPECT_EQ(id.toList(), (std::vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}));
    const auto message = [](std::vector<double> v) -> std::string {
        try {
            RigidTransform::fromList(v);
        } catch (const ValueError& e) {
            return e.what();
        }
        return "";
    };
    EXPECT_EQ(message({1, 2, 3}),
              "a rigid transform needs 12 numbers (3x3 rotation + translation)");
    EXPECT_EQ(message({1, 0, 0, 0, 1, 0, 0, 0, 1, std::numeric_limits<double>::infinity(), 0, 0}),
              "the transform contains non-finite numbers");
    EXPECT_EQ(message({2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}),
              "the rotation part of the transform is not a rotation");
    EXPECT_EQ(message({1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0}),
              "the rotation part of the transform is a reflection");
    EXPECT_EQ(message({1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}), "");
}

// test_density.py::test_rigid_transform_algebra
TEST(RigidTransform, MdatkAlgebraTest) {
    const RigidTransform t{rotvecMatrix({0.3, -0.2, 0.5}), {1.0, 2.0, 3.0}};
    std::mt19937 gen(2);
    std::normal_distribution<double> normal;
    std::vector<Vec3> p(5);
    for (Vec3& x : p) {
        x = {normal(gen), normal(gen), normal(gen)};
    }
    const std::vector<Vec3> back = t.inverse().apply(t.apply(p));
    const RigidTransform restored = RigidTransform::fromList(t.toList(12));
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Vec3 a = t.apply(p[k]), b = restored.apply(p[k]);
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(back[k][j], p[k][j], 1e-12);
            EXPECT_NEAR(a[j], b[j], 1e-9);
        }
    }
    EXPECT_THROW(RigidTransform::fromList(std::vector<double>{1, 2, 3}), ValueError);
    // a rotation applied about a point keeps that point fixed
    const Vec3 c{4.0, -2.0, 7.5};
    const RigidTransform r{rotvecMatrix({0.1, 0.4, -0.3}), {0, 0, 0}};
    const Vec3 fixed = r.about(c).apply(c);
    for (std::size_t j = 0; j < 3; ++j) {
        EXPECT_NEAR(fixed[j], c[j], 1e-12);
    }
    EXPECT_NEAR(rotationAngle(rotvecMatrix({0, 0, 0.5})), 0.5 * 180 / 3.141592653589793, 1e-12);
}
