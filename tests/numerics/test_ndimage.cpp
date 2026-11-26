// The scipy.ndimage pieces of src/numerics/ndimage.cpp that the density fit adds (per-axis Gaussian
// in modes "reflect" and "constant", order-1 map_coordinates in mode "constant" with a cval) against
// SciPy 1.18.1 (data/numerics/numerics_cases.json "gaussian", "map_coords"), and density_centroid.
// The surface's and the shadow's paths are tested in tests/render.

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/io/mrc.hpp"
#include "mdbridge/numerics/density.hpp"
#include "mdbridge/numerics/ndimage.hpp"
#include "numerics_cases.hpp"

using namespace mdbridge;
using namespace mdbridge::test::numerics_cases;

TEST(Ndimage, GaussianFilterPerAxis) {
    for (const core::Value& c : cases("gaussian")) {
        const std::vector<std::size_t> s = counts(at(c, "shape"));
        const std::vector<double> sigma = nums(at(c, "sigma"));
        const std::string mode = str(at(c, "mode"));
        std::vector<float> d(s[0] * s[1] * s[2]);
        for (std::size_t i = 0; i < d.size(); ++i) {
            d[i] = static_cast<float>(value(i));
        }
        numerics::gaussianFilter(d, {s[0], s[1], s[2]}, {sigma[0], sigma[1], sigma[2]},
                                 mode == "reflect" ? numerics::ExtendMode::Reflect : numerics::ExtendMode::Constant);
        const std::vector<double> want = nums(at(c, "out"));
        ASSERT_EQ(d.size(), want.size());
        std::size_t exact = 0;
        for (std::size_t i = 0; i < d.size(); ++i) {
            // the kernel's exp may differ from numpy's in the last bit: float32 results within one ulp
            const auto w = static_cast<float>(want[i]);
            EXPECT_LE(std::fabs(d[i] - w), std::nextafter(std::fabs(w), 1.0f) - std::fabs(w) + 1e-12f)
                << mode << " " << i;
            if (d[i] == w) {
                ++exact;
            }
        }
        EXPECT_GE(exact * 100, d.size() * 95) << "most values are scipy's bit for bit";
    }
}

TEST(Ndimage, MapCoordinatesConstant) {
    for (const core::Value& c : cases("map_coords")) {
        const std::vector<std::size_t> s = counts(at(c, "shape"));
        const std::vector<std::size_t> nanAt = counts(at(c, "nan_at"));
        std::vector<float> d(s[0] * s[1] * s[2]);
        for (std::size_t i = 0; i < d.size(); ++i) {
            d[i] = static_cast<float>(value(i));
        }
        d[(nanAt[0] * s[1] + nanAt[1]) * s[2] + nanAt[2]] = std::numeric_limits<float>::quiet_NaN();
        const double cval = num(at(c, "cval"));
        std::vector<double> zyx;
        for (const core::Value& p : list(at(c, "points"))) {
            for (double v : nums(p)) {
                zyx.push_back(v);
            }
        }
        const std::vector<float> got = numerics::mapCoordinatesLinearConstant(d, {s[0], s[1], s[2]}, zyx, cval);
        const std::vector<double> want = nums(at(c, "out"));
        ASSERT_EQ(got.size(), want.size());
        for (std::size_t i = 0; i < got.size(); ++i) {
            const auto w = static_cast<float>(want[i]);
            if (std::isnan(w)) {
                EXPECT_TRUE(std::isnan(got[i])) << i;
            } else {
                EXPECT_EQ(got[i], w) << "point " << i << " cval " << cval;
            }
        }
    }
    EXPECT_THROW(numerics::mapCoordinatesLinearConstant(std::vector<float>(5), {1, 2, 3}, std::vector<double>{0, 0, 0}),
                 std::invalid_argument);
}

TEST(Ndimage, ReflectRepeatsForLongKernels) {
    // a constant line stays constant under "reflect" whatever the kernel length (the reflection
    // repeats), and loses mass under "constant" with cval 0
    std::vector<float> d(3, 2.0f);
    numerics::gaussianFilter(d, {1, 1, 3}, {0.0, 0.0, 5.0}, numerics::ExtendMode::Reflect);
    for (float v : d) {
        EXPECT_NEAR(v, 2.0f, 1e-6f);
    }
    std::vector<float> e(3, 2.0f);
    numerics::gaussianFilter(e, {1, 1, 3}, {0.0, 0.0, 5.0}, numerics::ExtendMode::Constant);
    EXPECT_LT(e[1], 1.0f);
}

TEST(Density, CentroidOfTheDensityAboveTheLevel) {
    io::DensityMap map;
    map.shapeXyz = {4, 3, 2};
    map.origin = {10.0, 20.0, 30.0};
    map.voxel = {1.0, 2.0, 0.5};
    map.data.assign(24, 0.0f);
    map.data[(1 * 3 + 2) * 4 + 3] = 3.0f;  // x 3, y 2, z 1
    map.data[0] = std::numeric_limits<float>::quiet_NaN();  // carries no weight
    core::Vec3 c = numerics::densityCentroid(map, 1.0);
    EXPECT_DOUBLE_EQ(c[0], 13.0);
    EXPECT_DOUBLE_EQ(c[1], 24.0);
    EXPECT_DOUBLE_EQ(c[2], 30.5);
    // nothing above the level: the centre of the box
    c = numerics::densityCentroid(map, 5.0);
    EXPECT_DOUBLE_EQ(c[0], 11.5);
    EXPECT_DOUBLE_EQ(c[1], 22.0);
    EXPECT_DOUBLE_EQ(c[2], 30.25);
}
