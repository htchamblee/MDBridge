// MDATK's spatial helpers (mdatk/spatial.py) against what MDATK computes (tests/parity/
// make_spatial_cases.py, data/analysis/spatial_extra_cases.json, see data/README.md), ports of MDATK's
// test_pairs_within_* / test_kabsch_* (tests/test_analyses.py), and the shared numerics under them
// (core/linalg.hpp, core/detail/numpy_sum.hpp, core/detail/neighbors.hpp).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/linalg.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::analysis;
using core::Array;
using core::Value;
using test::kData;
using test::readFile;

namespace {

const Value& cases() {
    static const Value data = core::parseJson(
        readFile(kData / "analysis" / "spatial_extra_cases.json"), "spatial_extra_cases.json");
    return data;
}

const Array& list(const Value& value) {
    return *value.get<Array>();
}
const Array& items(const Value& value, std::string_view key) {
    return list(*value.find(key));
}
const std::string& text(const Value& value, std::string_view key) {
    return *value.find(key)->get<std::string>();
}

double hexValue(const Value& value) {
    const std::string& s = *value.get<std::string>();
    return std::strtod(s.c_str(), nullptr);  // float.hex(), "inf", "-inf", "nan"
}
double hexValue(const Value& value, std::string_view key) {
    return hexValue(*value.find(key));
}
std::vector<double> hexValues(const Value& value) {
    std::vector<double> out;
    for (const Value& v : list(value)) {
        out.push_back(hexValue(v));
    }
    return out;
}
std::vector<double> hexValues(const Value& value, std::string_view key) {
    return hexValues(*value.find(key));
}
std::vector<std::int64_t> integers(const Value& value, std::string_view key) {
    std::vector<std::int64_t> out;
    for (const Value& v : items(value, key)) {
        out.push_back(*v.get<std::int64_t>());
    }
    return out;
}
std::vector<Vec3> points(const std::vector<double>& xyz) {
    std::vector<Vec3> out;
    for (std::size_t k = 0; k + 2 < xyz.size(); k += 3) {
        out.push_back({xyz[k], xyz[k + 1], xyz[k + 2]});
    }
    return out;
}

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

void expectSameValues(const std::vector<double>& got, const std::vector<double>& expected,
                      const std::string& context) {
    ASSERT_EQ(got.size(), expected.size()) << context;
    for (std::size_t k = 0; k < got.size(); ++k) {
        EXPECT_TRUE(sameDouble(got[k], expected[k]))
            << context << " [" << k << "]: " << got[k] << " vs " << expected[k];
    }
}

void expectNear(double got, double expected, double tolerance, const std::string& context) {
    EXPECT_LE(std::abs(got - expected), tolerance * std::max(1.0, std::abs(expected)))
        << context << ": " << got << " vs " << expected;
}

using Triple = std::tuple<std::int64_t, std::int64_t, double>;

std::vector<Triple> triples(const PairList& p) {
    std::vector<Triple> out;
    for (std::size_t k = 0; k < p.size(); ++k) {
        out.emplace_back(p.i[k], p.j[k], p.d[k]);
    }
    return out;
}

// the pairs MDBridge returns: sorted by (i, j), each pair once
void expectCanonicalOrder(const PairList& p, const std::string& context) {
    ASSERT_EQ(p.i.size(), p.j.size()) << context;
    ASSERT_EQ(p.i.size(), p.d.size()) << context;
    for (std::size_t k = 1; k < p.size(); ++k) {
        EXPECT_TRUE(std::pair(p.i[k - 1], p.j[k - 1]) < std::pair(p.i[k], p.j[k])) << context;
    }
}

bool float32Exact(const std::vector<double>& xyz) {
    return std::all_of(xyz.begin(), xyz.end(),
                       [](double v) { return static_cast<double>(static_cast<float>(v)) == v; });
}

std::vector<float> asFloat(const std::vector<double>& xyz) {
    std::vector<float> out;
    for (const double v : xyz) {
        out.push_back(static_cast<float>(v));
    }
    return out;
}

// what scipy keeps: ((dx dx + dy dy) + dz dz) <= r * r, distance its square root
PairList bruteForce(const std::vector<double>& a, const std::vector<double>& b, double r) {
    PairList out;
    for (std::size_t i = 0; 3 * i < a.size(); ++i) {
        for (std::size_t j = 0; 3 * j < b.size(); ++j) {
            double s = 0;
            for (std::size_t k = 0; k < 3; ++k) {
                const double d = a[3 * i + k] - b[3 * j + k];
                s += d * d;
            }
            if (s <= r * r) {
                out.i.push_back(static_cast<std::int64_t>(i));
                out.j.push_back(static_cast<std::int64_t>(j));
                out.d.push_back(std::sqrt(s));
            }
        }
    }
    return out;
}

}  // namespace

// ---- pairs_within ------------------------------------------------------------------------------

TEST(PairsWithin, MatchesMdatkExactly) {
    const Array& all = items(cases(), "pairs");
    ASSERT_GT(all.size(), 70u);
    for (const Value& c : all) {
        const std::string name = text(c, "name");
        const std::vector<double> a = hexValues(c, "a"), b = hexValues(c, "b");
        const double cutoff = hexValue(c, "cutoff");
        if (const Value* error = c.find("error")) {
            try {
                (void)pairsWithin(a, b, cutoff);
                ADD_FAILURE() << name << ": no error";
            } catch (const core::ValueError& e) {
                EXPECT_EQ(e.what(), text(*error, "message")) << name;
                EXPECT_EQ(text(*error, "type"), "ValueError") << name;
            }
            continue;
        }
        std::vector<Triple> expected;
        const std::vector<std::int64_t> ei = integers(c, "i"), ej = integers(c, "j");
        const std::vector<double> ed = hexValues(c, "d");
        for (std::size_t k = 0; k < ei.size(); ++k) {
            expected.emplace_back(ei[k], ej[k], ed[k]);
        }
        std::sort(expected.begin(), expected.end());
        const PairList got = pairsWithin(a, b, cutoff);
        expectCanonicalOrder(got, name);
        ASSERT_EQ(got.size(), expected.size()) << name;
        const std::vector<Triple> mine = triples(got);
        for (std::size_t k = 0; k < mine.size(); ++k) {
            EXPECT_EQ(std::get<0>(mine[k]), std::get<0>(expected[k])) << name;
            EXPECT_EQ(std::get<1>(mine[k]), std::get<1>(expected[k])) << name;
            EXPECT_TRUE(sameDouble(std::get<2>(mine[k]), std::get<2>(expected[k])))
                << name << ": " << std::get<2>(mine[k]) << " vs " << std::get<2>(expected[k]);
        }
        // float32 positions (a Frame's) give the same pairs: MDATK widens them to float64
        if (float32Exact(a) && float32Exact(b)) {
            const PairList narrow = pairsWithin(asFloat(a), asFloat(b), cutoff);
            EXPECT_EQ(narrow.i, got.i) << name;
            EXPECT_EQ(narrow.j, got.j) << name;
            expectSameValues(narrow.d, got.d, name);
        }
    }
}

TEST(PairsWithin, BoundaryCasesIncludeSqrtSquareDisagreements) {
    // the generator crafted pairs where sqrt(d2) <= r and d2 <= r * r disagree: scipy (and
    // MDBridge) follow the squares
    int disagreements = 0;
    for (const Value& c : items(cases(), "pairs")) {
        if (c.find("error")) {
            continue;
        }
        const std::vector<double> a = hexValues(c, "a"), b = hexValues(c, "b");
        const double r = hexValue(c, "cutoff");
        for (std::size_t i = 0; 3 * i < a.size(); ++i) {
            for (std::size_t j = 0; 3 * j < b.size(); ++j) {
                double s = 0;
                for (std::size_t k = 0; k < 3; ++k) {
                    const double d = a[3 * i + k] - b[3 * j + k];
                    s += d * d;
                }
                disagreements += (std::sqrt(s) <= r) != (s <= r * r);
            }
        }
    }
    EXPECT_GT(disagreements, 0);
}

TEST(PairsWithin, CellGridMatchesBruteForce) {
    std::mt19937_64 random(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    struct Setup {
        std::size_t na, nb;
        double size, cutoff;
    };
    for (const Setup s : {Setup{1, 1, 3, 2}, Setup{50, 70, 10, 2.5}, Setup{900, 1100, 40, 4.0},
                          Setup{2000, 1500, 25, 3.2}, Setup{300, 300, 5, 0.3},
                          Setup{400, 500, 3000, 3.0}, Setup{200, 200, 1e5, 50.0}}) {
        std::vector<double> a(3 * s.na), b(3 * s.nb);
        for (double& v : a) {
            v = unit(random) * s.size - s.size / 3;
        }
        for (double& v : b) {
            v = unit(random) * s.size;
        }
        // clustered copies: shared and nearby coordinates
        for (std::size_t k = 0; k < std::min(s.na, s.nb) / 4; ++k) {
            for (std::size_t c = 0; c < 3; ++c) {
                b[3 * k + c] = a[3 * k + c] + (k % 3 == 0 ? 0.0 : unit(random) * s.cutoff);
            }
        }
        for (const double cutoff : {s.cutoff, -s.cutoff, s.cutoff * 0.5}) {
            const PairList got = pairsWithin(a, b, cutoff);
            const PairList expected = bruteForce(a, b, cutoff);
            const std::string context =
                std::to_string(s.na) + "x" + std::to_string(s.nb) + " r=" + std::to_string(cutoff);
            expectCanonicalOrder(got, context);
            EXPECT_EQ(got.i, expected.i) << context;
            EXPECT_EQ(got.j, expected.j) << context;
            expectSameValues(got.d, expected.d, context);
        }
    }
}

TEST(PairsWithin, SpecialInputs) {
    const std::vector<double> a = {0, 0, 0}, b = {1, 0, 0, 3, 0, 0};
    EXPECT_TRUE(pairsWithin(std::vector<double>{}, b, 2.0).empty());
    EXPECT_TRUE(pairsWithin(a, std::vector<double>{}, 2.0).empty());
    // empty input wins over non-finite points, as in MDATK
    EXPECT_TRUE(pairsWithin(std::vector<double>{NAN, 0, 0}, std::vector<double>{}, 2.0).empty());
    EXPECT_THROW((void)pairsWithin(std::vector<double>{NAN, 0, 0}, b, 2.0), core::ValueError);
    EXPECT_EQ(pairsWithin(a, b, -2.0).size(), 1u);
    EXPECT_EQ(pairsWithin(a, b, INFINITY).size(), 2u);
    EXPECT_TRUE(pairsWithin(a, b, NAN).empty());
    EXPECT_THROW((void)pairsWithin(std::vector<double>{0, 0}, b, 2.0), std::invalid_argument);
    // a box that is not usable takes the open path
    const PairList open = pairsWithin(a, b, 2.0, Box{0, 10, 10, 90, 90, 90});
    EXPECT_EQ(open.size(), 1u);
}

TEST(PairsWithin, PeriodicPathIsAHookForW1_1) {
    const std::vector<double> a = {0.2, 5, 5}, b = {9.2, 5, 5};
    EXPECT_THROW((void)pairsWithin(a, b, 2.0, Box{10, 10, 10, 90, 90, 90}), std::logic_error);
    EXPECT_THROW((void)pairsWithin(asFloat(a), asFloat(b), 2.0, Box{10, 10, 10, 90, 90, 90}),
                 std::logic_error);
    // empty inputs return before the box is looked at
    EXPECT_TRUE(pairsWithin(std::vector<double>{}, b, 2.0, Box{10, 10, 10, 90, 90, 90}).empty());
}

// tests/test_analyses.py::test_pairs_within_matches_brute_force_and_pbc, open part
TEST(PairsWithin, MdatkTestMatchesBruteForce) {
    const Value& t = *cases().find("tests")->find("pairs_within_matches_brute_force");
    const std::vector<double> a = hexValues(t, "a"), b = hexValues(t, "b");
    const PairList got = pairsWithin(a, b, hexValue(t, "cutoff"));
    std::set<std::pair<std::int64_t, std::int64_t>> mine, brute;
    for (std::size_t k = 0; k < got.size(); ++k) {
        mine.emplace(got.i[k], got.j[k]);
    }
    for (std::size_t i = 0; i < 60; ++i) {
        for (std::size_t j = 0; j < 50; ++j) {
            const double dx = a[3 * i] - b[3 * j], dy = a[3 * i + 1] - b[3 * j + 1],
                         dz = a[3 * i + 2] - b[3 * j + 2];
            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);  // np.linalg.norm
            if (dist <= 4.0) {
                brute.emplace(i, j);
            }
        }
    }
    EXPECT_EQ(mine, brute);
    EXPECT_FALSE(mine.empty());
    for (std::size_t k = 0; k < got.size(); ++k) {
        const auto i = static_cast<std::size_t>(got.i[k]), j = static_cast<std::size_t>(got.j[k]);
        const double dx = a[3 * i] - b[3 * j], dy = a[3 * i + 1] - b[3 * j + 1],
                     dz = a[3 * i + 2] - b[3 * j + 2];
        const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        EXPECT_LE(std::abs(got.d[k] - dist), 1e-7 * std::abs(dist));  // assert_allclose rtol
    }
}

// the periodic part of the same MDATK test: enabled by W1.1 (periodic neighbour search)
TEST(PairsWithin, DISABLED_MdatkTestPeriodic) {
    const PairList got =
        pairsWithin(std::vector<double>{0.2, 5, 5}, std::vector<double>{9.2, 5, 5}, 2.0,
                    Box{10, 10, 10, 90, 90, 90});
    ASSERT_EQ(got.size(), 1u);
    EXPECT_NEAR(got.d[0], 1.0, 1e-4);
}

// ---- boxes and minimum images ------------------------------------------------------------------

TEST(UsableBox, MatchesMdatk) {
    for (const Value& c : items(cases(), "usable_box")) {
        const std::vector<double> box = hexValues(c, "box");
        const std::optional<Box> got = usableBox(box);
        const Value& out = *c.find("out");
        if (out.isNull()) {
            EXPECT_FALSE(got) << box.size();
        } else {
            ASSERT_TRUE(got) << box.size();
            expectSameValues({got->begin(), got->end()}, hexValues(out), "usable_box");
        }
    }
    EXPECT_FALSE(usableBox(std::optional<Box>{}));
    EXPECT_FALSE(usableBox(std::optional<std::array<float, 6>>{}));
    const std::optional<Box> widened =
        usableBox(std::optional<std::array<float, 6>>{{23.7f, 31.3f, 27.1f, 90.f, 90.f, 90.f}});
    ASSERT_TRUE(widened);
    EXPECT_EQ((*widened)[0], static_cast<double>(23.7f));
}

TEST(MinimumImageVectors, MatchesMdatk) {
    for (const Value& c : items(cases(), "min_image")) {
        const std::vector<double> vectors = hexValues(c, "vectors");
        std::optional<Box> box;
        const Value& given = *c.find("box");
        std::string context = "no box";
        if (!given.isNull()) {
            box = usableBox(hexValues(given));  // MDATK's _usable_box (shape 6, lengths > 0)
            context = "box of " + std::to_string(list(given).size());
            for (const double v : hexValues(given)) {
                context += " " + std::to_string(v);
            }
        }
        expectSameValues(minimumImageVectors(vectors, box), hexValues(c, "out"), context);
    }
    const Vec3 v = minimumImageVector({15.0, -5.0, 4.0}, Box{10, 10, 10, 90, 90, 90});
    EXPECT_EQ(v, (Vec3{-5.0, -5.0, 4.0}));  // 1.5 -> 2, -0.5 -> -0: half to even
}

// ---- kabsch ------------------------------------------------------------------------------------

TEST(Kabsch, MatchesMdatk) {
    for (const Value& c : items(cases(), "kabsch")) {
        const std::string name = text(c, "name");
        const std::vector<Vec3> mobile = points(hexValues(c, "mobile"));
        const std::vector<Vec3> reference = points(hexValues(c, "reference"));
        std::vector<double> weights;
        const bool weighted = !c.find("weights")->isNull();
        if (weighted) {
            weights = hexValues(c, "weights");
        }
        const KabschResult got =
            weighted ? kabsch(mobile, reference, std::span<const double>(weights))
                     : kabsch(mobile, reference);
        const std::vector<double> r = hexValues(c, "r"), mc = hexValues(c, "mc"),
                                  rc = hexValues(c, "rc");
        for (std::size_t k = 0; k < 9; ++k) {
            expectNear(got.rotation[k / 3][k % 3], r[k], 1e-14, name + " R");
        }
        for (std::size_t k = 0; k < 3; ++k) {
            expectNear(got.mobileCenter[k], mc[k], 1e-14, name + " mobile centre");
            expectNear(got.referenceCenter[k], rc[k], 1e-14, name + " reference centre");
        }
        const double det = core::det3(got.rotation);
        if (mobile.size() >= 3) {
            EXPECT_NEAR(det, 1.0, 1e-14) << name;
        }
    }
}

// tests/test_analyses.py::test_kabsch_recovers_rotation
TEST(Kabsch, MdatkTestRecoversRotation) {
    const Value& t = *cases().find("tests")->find("kabsch_recovers_rotation");
    const std::vector<Vec3> x = points(hexValues(t, "x")), y = points(hexValues(t, "y"));
    const KabschResult k = kabsch(x, y);
    for (std::size_t p = 0; p < x.size(); ++p) {
        for (std::size_t i = 0; i < 3; ++i) {
            double v = k.referenceCenter[i];
            for (std::size_t j = 0; j < 3; ++j) {
                v += (x[p][j] - k.mobileCenter[j]) * k.rotation[i][j];
            }
            EXPECT_NEAR(v, y[p][i], 1e-8);
        }
    }
}

TEST(Kabsch, Errors) {
    const std::vector<Vec3> three = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const std::vector<Vec3> two = {{0, 0, 0}, {1, 0, 0}};
    EXPECT_THROW((void)kabsch(three, two), std::invalid_argument);
    const std::vector<double> w = {1, 2};
    EXPECT_THROW((void)kabsch(three, three, std::span<const double>(w)), std::invalid_argument);
    const std::vector<Vec3> bad = {{0, 0, 0}, {1, 0, NAN}, {0, 1, 0}};
    try {
        (void)kabsch(bad, three);
        ADD_FAILURE() << "no error";
    } catch (const LinAlgError& e) {
        EXPECT_STREQ(e.what(), "SVD did not converge");
        EXPECT_STREQ(e.pythonName(), "LinAlgError");
    }
}

// ---- core/linalg: the SVDs behind kabsch and the ring planes ----------------------------------

TEST(Linalg, RingGeometrySvdMatchesNumpy) {
    for (const Value& c : items(cases(), "rings")) {
        const auto k = static_cast<std::size_t>(*c.find("k")->get<std::int64_t>());
        const std::vector<Vec3> pos = points(hexValues(c, "positions"));
        const std::vector<Vec3> centroids = points(hexValues(c, "centroids"));
        const std::vector<Vec3> normals = points(hexValues(c, "normals"));
        const std::vector<double> planarity = hexValues(c, "planarity");
        std::vector<Vec3> centred = pos;
        for (std::size_t p = 0; p < pos.size(); ++p) {
            for (std::size_t i = 0; i < 3; ++i) {
                centred[p][i] = pos[p][i] - centroids[p / k][i];
            }
        }
        const std::vector<core::SvdRows3> svd = core::svdRows3Batch(centred, k);
        ASSERT_EQ(svd.size(), centroids.size());
        for (std::size_t m = 0; m < svd.size(); ++m) {
            const std::string context = "ring size " + std::to_string(k) + " #" + std::to_string(m);
            EXPECT_NEAR(svd[m].s[2] / std::sqrt(static_cast<double>(k)), planarity[m], 1e-14)
                << context;
            // the normal (last row of vt) up to LAPACK's sign choice
            const Vec3& n = svd[m].vt[2];
            const double dot = n[0] * normals[m][0] + n[1] * normals[m][1] + n[2] * normals[m][2];
            const double sign = dot < 0 ? -1.0 : 1.0;
            for (std::size_t i = 0; i < 3; ++i) {
                EXPECT_NEAR(sign * n[i], normals[m][i], 1e-13) << context;
            }
        }
    }
    EXPECT_THROW((void)core::svdRows3Batch(std::vector<Vec3>(5), 3), std::invalid_argument);
}

TEST(Linalg, Svd3x3Reconstructs) {
    std::mt19937_64 random(11);
    std::normal_distribution<double> normal;
    std::vector<core::Mat3> inputs;
    for (int k = 0; k < 200; ++k) {
        core::Mat3 a{};
        for (auto& row : a) {
            for (double& v : row) {
                v = normal(random);
            }
        }
        if (k % 4 == 1) {  // rank 2
            for (std::size_t i = 0; i < 3; ++i) {
                a[i][2] = a[i][0] * 0.5 - a[i][1] * 2.0;
            }
        } else if (k % 4 == 2) {  // rank 1
            for (std::size_t i = 0; i < 3; ++i) {
                a[i][1] = a[i][0] * 3.0;
                a[i][2] = -a[i][0];
            }
        }
        inputs.push_back(a);
    }
    inputs.push_back({});  // zero
    inputs.push_back({{{2, 0, 0}, {0, 2, 0}, {0, 0, 2}}});  // repeated singular values
    for (const core::Mat3& a : inputs) {
        const core::Svd3 svd = core::svd3x3(a);
        EXPECT_GE(svd.s[0], svd.s[1]);
        EXPECT_GE(svd.s[1], svd.s[2]);
        EXPECT_GE(svd.s[2], 0.0);
        double scale = 1;
        for (const auto& row : a) {
            for (const double v : row) {
                scale = std::max(scale, std::abs(v));
            }
        }
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                double usv = 0, uu = 0, vv = 0;
                for (std::size_t m = 0; m < 3; ++m) {
                    usv += svd.u[i][m] * svd.s[m] * svd.vt[m][j];
                    uu += svd.u[m][i] * svd.u[m][j];
                    vv += svd.vt[i][m] * svd.vt[j][m];
                }
                EXPECT_NEAR(usv, a[i][j], 1e-14 * scale);
                EXPECT_NEAR(uu, i == j ? 1.0 : 0.0, 1e-14);
                EXPECT_NEAR(vv, i == j ? 1.0 : 0.0, 1e-14);
            }
        }
    }
}

// ---- core/detail/numpy_sum: numpy's pairwise float64 summation ---------------------------------

TEST(NumpySum, MatchesNumpyFloat64) {
    const auto value = [](std::int64_t k) {
        return static_cast<double>((k * 7919) % 1000 - 500) / 37.0 *
               (1.0 + static_cast<double>(k) / 3.0);
    };
    const Array& all = items(cases(), "sums");
    ASSERT_GT(all.size(), 300u);
    for (const Value& c : all) {
        const std::int64_t n = *c.find("n")->get<std::int64_t>();
        std::vector<double> x;
        for (std::int64_t k = 0; k < n; ++k) {
            x.push_back(value(k));
        }
        const std::string context = "n = " + std::to_string(n);
        EXPECT_TRUE(sameDouble(core::detail::numpySum(x), hexValue(c, "sum"))) << context;
        EXPECT_TRUE(sameDouble(core::detail::numpyMean(x), hexValue(c, "mean"))) << context;
    }
    // np.sum of negative zeros is +0: the reduction starts from the identity 0
    const std::vector<double> zeros(5, -0.0);
    EXPECT_TRUE(sameDouble(core::detail::numpySum(zeros), 0.0));
}
