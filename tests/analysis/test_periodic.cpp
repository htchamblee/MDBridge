// The periodic pair search (pairsWithinPeriodic) on every case of data/analysis/distance_cases.json
// (data/README.md): MDBridge's float64 minimum-image search (docs/DIFFERENCES.md AN-01 to AN-03)
// against the case's `exact` column, computed independently with numpy
// (tests/parity/improvements_periodic.py). Pairs and order (sorted by (i, j)) must match exactly
// and so must the distance bits; for triclinic cells, whose box vectors come from libm's cos / sin,
// distances may differ by a few float64 ulps on another C library and a pair within that of the
// cutoff may be in or out. The MDAnalysis columns (capped_distance, the old reference behaviour,
// docs/specs/m3-capped-distance.md) are only compared to count where MDBridge now differs.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::analysis;
using core::Array;
using core::Value;
using test::kData;
using test::readFile;

namespace {

const Value& distanceCases() {
    static const Value data = core::parseJson(
        readFile(kData / "analysis" / "distance_cases.json"), "distance_cases.json");
    return data;
}

const Array& list(const Value& value) {
    return *value.get<Array>();
}

double number(const Value& value) {
    if (const double* d = value.get<double>()) {
        return *d;
    }
    return static_cast<double>(*value.get<std::int64_t>());
}

double hexNumber(const Value& value) {
    return std::strtod(value.get<std::string>()->c_str(), nullptr);  // float.hex(), "nan", "inf"
}

std::vector<float> points(const Value& value) {
    std::vector<float> out;
    for (const Value& p : list(value)) {
        for (const Value& x : list(p)) {
            out.push_back(static_cast<float>(number(x)));
        }
    }
    return out;
}

// a lattice case's points: origin + spacing * (ix, iy, iz), x fastest and z slowest
std::vector<float> lattice(const Value& spec) {
    const Array& origin = list(*spec.find("origin"));
    const Array& counts = list(*spec.find("counts"));
    const double spacing = number(*spec.find("spacing"));
    std::vector<float> out;
    const auto count = [&](std::size_t k) { return *counts[k].get<std::int64_t>(); };
    for (std::int64_t iz = 0; iz < count(2); ++iz) {
        for (std::int64_t iy = 0; iy < count(1); ++iy) {
            for (std::int64_t ix = 0; ix < count(0); ++ix) {
                const std::int64_t idx[3] = {ix, iy, iz};
                for (std::size_t k = 0; k < 3; ++k) {
                    out.push_back(static_cast<float>(number(origin[k]) +
                                                     spacing * static_cast<double>(idx[k])));
                }
            }
        }
    }
    return out;
}

template <class T>
std::string sha256Of(const std::vector<T>& values) {
    std::string bytes(values.size() * sizeof(T), '\0');  // little-endian hosts only (x86, ARM)
    std::memcpy(bytes.data(), values.data(), bytes.size());
    return core::sha256Hex(bytes);
}

using Pair = std::pair<std::int64_t, std::int64_t>;

}  // namespace

TEST(PairsWithinPeriodic, MatchesExactMinimumImage) {
    const Array& all = list(*distanceCases().find("cases"));
    ASSERT_EQ(all.size(), 495u);
    int exact = 0, nearly = 0, failed = 0, differFromMdanalysis = 0;
    for (const Value& c : all) {
        const std::string name = *c.find("name")->get<std::string>();
        SCOPED_TRACE(name);
        const bool isLattice = c.find("lattice_a") != nullptr;
        const std::vector<float> a = isLattice ? lattice(*c.find("lattice_a")) : points(*c.find("a"));
        const std::vector<float> b = isLattice ? lattice(*c.find("lattice_b")) : points(*c.find("b"));
        Box box{};
        const Array& boxValues = list(*c.find("box"));
        for (std::size_t k = 0; k < 6; ++k) {
            box[k] = number(boxValues[k]);
        }
        const double cutoff = number(*c.find("cutoff"));
        const bool orthogonal = static_cast<float>(box[3]) == 90.0f &&
                                static_cast<float>(box[4]) == 90.0f &&
                                static_cast<float>(box[5]) == 90.0f;
        const Value& expected = *c.find("exact");

        PairList got;
        ASSERT_NO_THROW(got = pairsWithinPeriodic(a, b, cutoff, box));  // never an error (AN-02)
        ASSERT_EQ(got.i.size(), got.j.size());
        ASSERT_EQ(got.i.size(), got.d.size());
        for (std::size_t k = 1; k < got.size(); ++k) {
            EXPECT_TRUE(std::pair(got.i[k - 1], got.j[k - 1]) < std::pair(got.i[k], got.j[k]));
        }
        if (isLattice) {
            std::vector<std::int64_t> flat;
            for (std::size_t k = 0; k < got.size(); ++k) {
                flat.push_back(got.i[k]);
                flat.push_back(got.j[k]);
            }
            EXPECT_EQ(got.size(), static_cast<std::size_t>(*expected.find("count")->get<std::int64_t>()));
            EXPECT_EQ(sha256Of(flat), *expected.find("pairs_sha256")->get<std::string>());
            const bool bits = sha256Of(got.d) == *expected.find("distances_sha256")->get<std::string>();
            if (orthogonal) {
                EXPECT_TRUE(bits);
            }
            ++(bits ? exact : nearly);
            if (c.find("error") != nullptr || *c.find("count")->get<std::int64_t>() !=
                                                  *expected.find("count")->get<std::int64_t>() ||
                *c.find("pairs_sha256")->get<std::string>() != *expected.find("pairs_sha256")->get<std::string>() ||
                *c.find("distances_sha256")->get<std::string>() != *expected.find("distances_sha256")->get<std::string>()) {
                ++differFromMdanalysis;
            }
            continue;
        }

        std::vector<Pair> expectedPairs, gotPairs;
        for (const Value& p : list(*expected.find("pairs"))) {
            expectedPairs.emplace_back(*list(p)[0].get<std::int64_t>(), *list(p)[1].get<std::int64_t>());
        }
        std::vector<double> expectedDistances;
        for (const Value& d : list(*expected.find("distances"))) {
            expectedDistances.push_back(hexNumber(d));
        }
        for (std::size_t k = 0; k < got.size(); ++k) {
            gotPairs.emplace_back(got.i[k], got.j[k]);
        }
        bool same = gotPairs == expectedPairs;
        for (std::size_t k = 0; same && k < got.size(); ++k) {
            same = std::memcmp(&got.d[k], &expectedDistances[k], sizeof(double)) == 0;
        }
        if (const Value* mda = c.find("pairs"); mda == nullptr || list(*mda).size() != expectedPairs.size()) {
            ++differFromMdanalysis;
        } else {
            for (std::size_t k = 0; k < expectedPairs.size(); ++k) {
                const Array& p = list(list(*mda)[k]);
                const double md = hexNumber(list(*c.find("distances"))[k]);
                if (Pair(*p[0].get<std::int64_t>(), *p[1].get<std::int64_t>()) != expectedPairs[k] ||
                    std::memcmp(&md, &expectedDistances[k], sizeof(double)) != 0) {
                    ++differFromMdanalysis;
                    break;
                }
            }
        }
        if (same) {
            ++exact;
            continue;
        }
        if (orthogonal) {
            ADD_FAILURE() << "orthogonal cell: not an exact match";
            ++failed;
            continue;
        }
        // triclinic: distances within a few ulps; pairs that close to the cutoff may differ
        std::map<Pair, double> mine, theirs;
        for (std::size_t k = 0; k < got.size(); ++k) {
            mine.emplace(gotPairs[k], got.d[k]);
        }
        for (std::size_t k = 0; k < expectedPairs.size(); ++k) {
            theirs.emplace(expectedPairs[k], expectedDistances[k]);
        }
        const double band = 1e-13 * std::max(1.0, std::abs(cutoff));
        bool ok = true;
        for (const auto& [pair, d] : theirs) {
            const auto found = mine.find(pair);
            if (found == mine.end()) {
                ok = ok && std::abs(d - std::abs(cutoff)) <= band;
            } else {
                ok = ok && std::abs(found->second - d) <= 1e-13 * std::max(1.0, d);
            }
        }
        for (const auto& [pair, d] : mine) {
            if (theirs.count(pair) == 0) {
                ok = ok && std::abs(d - std::abs(cutoff)) <= band;
            }
        }
        EXPECT_TRUE(ok);
        ++(ok ? nearly : failed);
    }
    std::printf("periodic cases: %d exact, %d within float64 rounding, %d failed; %d differ from "
                "MDAnalysis' capped_distance\n",
                exact, nearly, failed, differFromMdanalysis);
}

// The behaviours MDATK / MDAnalysis got wrong (docs/IMPROVEMENTS.md AN-01, AN-02), directly.
TEST(PairsWithinPeriodic, ExactBoundaryForEveryBoxSize) {
    // a pair at exactly 3.5 A is kept whatever the box length (MDAnalysis' brute force dropped it
    // for L = 30, 18.3, 40, 60)
    for (const double length : {25.0, 30.0, 31.7, 18.3, 40.0, 60.0, 77.1}) {
        const Box box{length, length, length, 90, 90, 90};
        const PairList got = pairsWithinPeriodic(std::vector<float>{1, 1, 1}, std::vector<float>{4.5f, 1, 1}, 3.5, box);
        ASSERT_EQ(got.size(), 1u) << length;
        EXPECT_EQ(got.d[0], 3.5) << length;
        // across the face: from L - 0.5 to 3.0 is 3.5 too
        const auto x = static_cast<float>(length - 0.5);
        const PairList across =
            pairsWithinPeriodic(std::vector<float>{x, 1, 1}, std::vector<float>{3.0f, 1, 1}, 3.5, box);
        ASSERT_EQ(across.size(), 1u) << length;
        EXPECT_NEAR(across.d[0], 3.5, 1e-12) << length;
    }
}

TEST(PairsWithinPeriodic, FarOutsidePointsAreExact) {
    const Box box{25, 25, 25, 90, 90, 90};
    const PairList got = pairsWithinPeriodic(std::vector<float>{0, 0, 0},
                                             std::vector<float>{1e4f, 0, 0, 1e10f + 0.0f, 0, 0, 1e4f + 3, 0, 0}, 5.0, box);
    ASSERT_EQ(got.size(), 3u);
    EXPECT_EQ(got.d[0], 0.0);  // MDAnalysis: 2.2e-4
    EXPECT_EQ(got.d[2], 3.0);
}

TEST(PairsWithinPeriodic, NoSizeDependentErrorsAndSameResultForAnySize) {
    // cutoffs MDAnalysis' grid refused ("too large for box", "must be positive") on large inputs
    std::vector<float> a, b;
    for (int k = 0; k < 40; ++k) {
        a.insert(a.end(), {static_cast<float>(k * 0.61), static_cast<float>(k * 1.37), static_cast<float>(k * 0.29)});
        b.insert(b.end(), {static_cast<float>(24.9 - k * 0.53), static_cast<float>(k * 0.71), static_cast<float>(k * 1.13)});
    }
    for (const Box& box : {Box{25, 25, 25, 90, 90, 90}, Box{25, 27, 29, 60, 70, 80}}) {
        for (const double cutoff : {13.0, 30.0, -3.0, 0.0, 3.0}) {
            PairList all;
            ASSERT_NO_THROW(all = pairsWithinPeriodic(a, b, cutoff, box));
            // the first 5 points alone give the same pairs among them
            const PairList few = pairsWithinPeriodic(std::span(a).first(15), std::span(b).first(15), cutoff, box);
            std::size_t k = 0;
            for (std::size_t p = 0; p < all.size(); ++p) {
                if (all.i[p] < 5 && all.j[p] < 5) {
                    ASSERT_LT(k, few.size());
                    EXPECT_EQ(all.i[p], few.i[k]);
                    EXPECT_EQ(all.j[p], few.j[k]);
                    EXPECT_EQ(all.d[p], few.d[k]);
                    ++k;
                }
            }
            EXPECT_EQ(k, few.size());
        }
    }
}

// The cell list only prunes: the same pairs as comparing every pair, for many cells and cutoffs.
TEST(PairsWithinPeriodic, CellListMatchesAllPairs) {
    std::uint64_t state = 12345;
    const auto next = [&state] {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(state >> 11) * 0x1p-53;
    };
    for (const Box& box : {Box{25, 25, 25, 90, 90, 90}, Box{40, 33, 51, 90, 90, 90}, Box{25, 27, 29, 60, 70, 80},
                           Box{40, 40, 40, 109.4712, 109.4712, 109.4712}, Box{60, 60, 60, 90, 90, 120},
                           Box{30, 1e30, 30, 90, 90, 90}}) {
        std::vector<float> a, b;
        for (int k = 0; k < 300; ++k) {
            for (int c = 0; c < 3; ++c) {
                a.push_back(static_cast<float>(next() * 90 - 30));
                b.push_back(static_cast<float>(next() * 90 - 30));
            }
        }
        b.insert(b.end(), {2e5f, 1, 1});  // far outside: compared with everything
        const std::optional<PeriodicCell> cell = PeriodicCell::from(box);
        ASSERT_TRUE(cell);
        for (const double cutoff : {0.5, 2.5, 4.0, 7.9}) {
            const PairList got = pairsWithinPeriodic(a, b, cutoff, box);
            PairList brute;
            for (std::size_t i = 0; i < a.size() / 3; ++i) {
                for (std::size_t j = 0; j < b.size() / 3; ++j) {
                    const Vec3 v = cell->minimumImage({static_cast<double>(a[3 * i]) - b[3 * j],
                                                       static_cast<double>(a[3 * i + 1]) - b[3 * j + 1],
                                                       static_cast<double>(a[3 * i + 2]) - b[3 * j + 2]});
                    const double s = (v[0] * v[0] + v[1] * v[1]) + v[2] * v[2];
                    if (s <= cutoff * cutoff) {
                        brute.i.push_back(static_cast<std::int64_t>(i));
                        brute.j.push_back(static_cast<std::int64_t>(j));
                        brute.d.push_back(std::sqrt(s));
                    }
                }
            }
            EXPECT_EQ(got.i, brute.i) << cutoff;
            EXPECT_EQ(got.j, brute.j) << cutoff;
            EXPECT_EQ(got.d, brute.d) << cutoff;
        }
    }
}

TEST(PeriodicCell, OneOrthogonalityTestAndTriclinicImages) {
    // 90.0005 is not 90 in float32: a triclinic cell for pairs and vectors alike (MDATK's vectors
    // used np.allclose and were left uncorrected for the pair search's triclinic image)
    const std::optional<PeriodicCell> slanted = PeriodicCell::from(Box{10, 10, 10, 90.0005, 90, 90});
    ASSERT_TRUE(slanted);
    EXPECT_FALSE(slanted->orthogonal());
    EXPECT_TRUE(PeriodicCell::from(Box{10, 10, 10, 90.000001, 90, 90})->orthogonal());
    // a triclinic cell: the image across the slanted face is the short one
    const Box tric{25, 25, 25, 60, 70, 80};
    const Vec3 v = minimumImageVector({24.0, 0.0, 0.0}, tric);
    EXPECT_NEAR(std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]), 1.0, 1e-12);
    const Vec3 w = minimumImageVector({3.0, 4.0, 0.0}, tric);
    EXPECT_EQ(w, (Vec3{3.0, 4.0, 0.0}));  // within the cell: unchanged
    // no periodic geometry: unchanged vectors, plain distances
    EXPECT_FALSE(PeriodicCell::from(Box{25, 25, 25, 60, 60, 150}));
    EXPECT_EQ(minimumImageVector({24.0, 0.0, 0.0}, Box{25, 25, 25, 60, 60, 150}), (Vec3{24.0, 0.0, 0.0}));
    const PairList plain = pairsWithinPeriodic(std::vector<float>{0, 0, 0}, std::vector<float>{24, 0, 0, NAN, 0, 0},
                                               30.0, Box{25, 25, 25, 60, 60, 150});
    ASSERT_EQ(plain.size(), 1u);
    EXPECT_EQ(plain.d[0], 24.0);
}
