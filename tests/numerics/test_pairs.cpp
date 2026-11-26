// queryPairs (src/numerics/pairs.cpp) against cKDTree.query_pairs of SciPy 1.18.1: the same set
// (data/numerics/numerics_cases.json "pairs"), sorted.

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/numerics/pairs.hpp"
#include "numerics_cases.hpp"

using namespace mdbridge;
using namespace mdbridge::test::numerics_cases;

TEST(Pairs, SameSetAsCKDTree) {
    for (const core::Value& c : cases("pairs")) {
        const std::vector<double> points = nums(at(c, "points"));
        const double r = num(at(c, "r"));
        std::vector<std::pair<std::size_t, std::size_t>> want;
        for (const core::Value& p : list(at(c, "pairs"))) {
            const std::vector<std::size_t> ij = counts(p);
            want.emplace_back(ij[0], ij[1]);
        }
        EXPECT_EQ(numerics::queryPairs(points, r), want) << points.size() / 3 << " points, r = " << r;
    }
}

TEST(Pairs, BruteForceAgreesAndOddInputs) {
    // a spread the cell grid refuses (far more cells than points) goes pair by pair: same answer
    std::vector<double> pts;
    for (int i = 0; i < 50; ++i) {
        pts.insert(pts.end(), {std::pow(1.7, i), 0.0, 0.5 * (i % 3)});
    }
    const auto pairs = numerics::queryPairs(pts, 1.0);
    for (const auto& [i, j] : pairs) {
        EXPECT_LT(i, j);
    }
    std::vector<std::pair<std::size_t, std::size_t>> brute;
    for (std::size_t i = 0; i < 50; ++i) {
        for (std::size_t j = i + 1; j < 50; ++j) {
            const double dx = pts[3 * i] - pts[3 * j], dy = pts[3 * i + 1] - pts[3 * j + 1],
                         dz = pts[3 * i + 2] - pts[3 * j + 2];
            if ((dx * dx + dy * dy) + dz * dz <= 1.0) {
                brute.emplace_back(i, j);
            }
        }
    }
    EXPECT_EQ(pairs, brute);
    // a non-finite point is an error, as cKDTree's; a negative (NU-01) or NaN radius gives nothing
    const std::vector<double> odd = {0, 0, 0, std::numeric_limits<double>::quiet_NaN(), 0, 0, 0.5, 0, 0};
    EXPECT_THROW(numerics::queryPairs(odd, 1.0), std::invalid_argument);
    const std::vector<double> two = {0, 0, 0, 0.5, 0, 0};
    EXPECT_EQ(numerics::queryPairs(two, 1.0), (std::vector<std::pair<std::size_t, std::size_t>>{{0, 1}}));
    EXPECT_TRUE(numerics::queryPairs(two, -1.0).empty());  // cKDTree: {(0, 1)}
    EXPECT_TRUE(numerics::queryPairs(two, std::nan("")).empty());
    EXPECT_THROW(numerics::queryPairs(std::vector<double>{1.0, 2.0}, 1.0), std::invalid_argument);
}
