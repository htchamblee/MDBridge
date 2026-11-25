// The periodic pair search (pairsWithinPeriodic) against what MDAnalysis 2.11.0's capped_distance
// returns for every case of data/analysis/distance_cases.json (docs/specs/m3-capped-distance.md
// §7, data/README.md). Orthogonal cells must match exactly: method, pairs, order and distance
// bits. Triclinic and invalid cells within the spec's tolerance (§6): distances within
// 8 ulp32(M), M the largest of the pair's absolute coordinates and the cell lengths; the pair set
// exact except for the pairs the case lists in `band`; the order is reported.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
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

double ulp32(double x) {
    const double m = std::abs(x);
    if (!std::isfinite(m)) {
        return m;
    }
    const int exponent = std::ilogb(std::max(m, 0x1p-126));  // m = f * 2^exponent, 1 <= f < 2
    return std::ldexp(1.0, exponent - 23);
}

struct Counts {
    int exact = 0;      // pairs, order and every distance bit
    int tolerance = 0;  // same pairs, distances within the tolerance (order may differ)
    int band = 0;       // pair sets differ only by listed band pairs
    int failed = 0;
    int orderDiffers = 0;  // of tolerance / band cases: common pairs in another order
    int deviations = 0;    // documented deviations from MDAnalysis
    std::vector<std::string> failures;
};

}  // namespace

TEST(PairsWithinPeriodic, MatchesCappedDistance) {
    const Array& all = list(*distanceCases().find("cases"));
    ASSERT_EQ(all.size(), 495u);
    Counts counts;
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

        const bool grid = *c.find("method")->get<std::string>() == "grid";
        EXPECT_EQ(periodicSearchUsesGrid(a.size() / 3, b.size() / 3, cutoff, box), grid);

        if (const Value* error = c.find("error")) {
            const std::string expected = *error->find("message")->get<std::string>();
            try {
                (void)pairsWithinPeriodic(a, b, cutoff, box);
                ADD_FAILURE() << "no error, expected " << expected;
                ++counts.failed;
                counts.failures.push_back(name);
            } catch (const core::ValueError& e) {
                EXPECT_EQ(std::string(e.what()), expected);
                if (std::string(e.what()) == expected) {
                    ++counts.exact;
                } else {
                    ++counts.failed;
                    counts.failures.push_back(name);
                }
            }
            continue;
        }

        const PairList got = pairsWithinPeriodic(a, b, cutoff, box);
        ASSERT_EQ(got.i.size(), got.j.size());
        ASSERT_EQ(got.i.size(), got.d.size());
        if (name == "lengths [inf, 25.0, 25.0] triclinic") {
            // Deviation: MDAnalysis reports 126 of the 144 pairs, all at distance 0 (meaningless:
            // the points are 4 Å or more apart, cutoff 3). MDBridge treats the infinite vector as
            // not periodic, as for an orthogonal cell: each point only with itself.
            ASSERT_EQ(got.size(), a.size() / 3);
            for (std::size_t k = 0; k < got.size(); ++k) {
                EXPECT_EQ(got.i[k], static_cast<std::int64_t>(k));
                EXPECT_EQ(got.j[k], static_cast<std::int64_t>(k));
                EXPECT_EQ(got.d[k], 0.0);
            }
            ++counts.deviations;
            continue;
        }
        if (isLattice) {
            std::vector<std::int64_t> flat;
            for (std::size_t k = 0; k < got.size(); ++k) {
                flat.push_back(got.i[k]);
                flat.push_back(got.j[k]);
            }
            const auto expectedCount = static_cast<std::size_t>(*c.find("count")->get<std::int64_t>());
            const bool same = got.size() == expectedCount &&
                              sha256Of(flat) == *c.find("pairs_sha256")->get<std::string>() &&
                              sha256Of(got.d) == *c.find("distances_sha256")->get<std::string>();
            EXPECT_EQ(got.size(), expectedCount);
            EXPECT_EQ(sha256Of(flat), *c.find("pairs_sha256")->get<std::string>());
            if (orthogonal) {
                EXPECT_EQ(sha256Of(got.d), *c.find("distances_sha256")->get<std::string>());
            }
            if (same) {
                ++counts.exact;
            } else if (!orthogonal && got.size() == expectedCount &&
                       sha256Of(flat) == *c.find("pairs_sha256")->get<std::string>()) {
                ++counts.tolerance;  // pairs and order exact; distances not hashed per pair
            } else {
                ++counts.failed;
                counts.failures.push_back(name);
            }
            continue;
        }

        std::vector<Pair> expectedPairs;
        for (const Value& p : list(*c.find("pairs"))) {
            expectedPairs.emplace_back(*list(p)[0].get<std::int64_t>(),
                                       *list(p)[1].get<std::int64_t>());
        }
        std::vector<double> expectedDistances;
        for (const Value& d : list(*c.find("distances"))) {
            expectedDistances.push_back(hexNumber(d));
        }
        std::vector<Pair> gotPairs;
        for (std::size_t k = 0; k < got.size(); ++k) {
            gotPairs.emplace_back(got.i[k], got.j[k]);
        }
        bool exact = gotPairs == expectedPairs;
        for (std::size_t k = 0; exact && k < got.size(); ++k) {
            exact = std::memcmp(&got.d[k], &expectedDistances[k], sizeof(double)) == 0;
        }
        if (exact) {
            ++counts.exact;
            continue;
        }
        if (orthogonal) {
            ADD_FAILURE() << "orthogonal cell: not an exact match";
            ++counts.failed;
            counts.failures.push_back(name);
            continue;
        }

        // triclinic / invalid cell: the tolerance of §6
        std::set<Pair> band;
        if (const Value* listed = c.find("band")) {
            for (const Value& p : list(*listed)) {
                band.emplace(*list(p)[0].get<std::int64_t>(), *list(p)[1].get<std::int64_t>());
            }
        }
        std::map<Pair, double> mine, theirs;
        for (std::size_t k = 0; k < got.size(); ++k) {
            mine.emplace(gotPairs[k], got.d[k]);
        }
        for (std::size_t k = 0; k < expectedPairs.size(); ++k) {
            theirs.emplace(expectedPairs[k], expectedDistances[k]);
        }
        double lengths = 0;
        for (std::size_t k = 0; k < 3; ++k) {
            lengths = std::max(lengths, std::abs(static_cast<double>(static_cast<float>(box[k]))));
        }
        bool ok = true, bandOnly = false;
        for (const auto& [pair, d] : theirs) {
            const auto found = mine.find(pair);
            if (found == mine.end()) {
                bandOnly = true;
                ok = ok && band.count(pair) > 0;
                EXPECT_TRUE(band.count(pair) > 0) << "missing " << pair.first << ", " << pair.second;
                continue;
            }
            double m = lengths;
            for (std::size_t k = 0; k < 3; ++k) {
                m = std::max({m, std::abs(static_cast<double>(a[3 * static_cast<std::size_t>(pair.first) + k])),
                              std::abs(static_cast<double>(b[3 * static_cast<std::size_t>(pair.second) + k]))});
            }
            const double allowed = 8 * ulp32(m);
            const bool near = std::abs(found->second - d) <= allowed;
            ok = ok && near;
            EXPECT_TRUE(near) << pair.first << ", " << pair.second << ": " << found->second
                              << " vs " << d;
        }
        for (const auto& [pair, d] : mine) {
            if (theirs.count(pair) == 0) {
                bandOnly = true;
                ok = ok && band.count(pair) > 0;
                EXPECT_TRUE(band.count(pair) > 0) << "extra " << pair.first << ", " << pair.second;
            }
        }
        if (!ok) {
            ++counts.failed;
            counts.failures.push_back(name);
            continue;
        }
        // the order of the pairs both lists hold
        std::vector<Pair> commonMine, commonTheirs;
        for (const Pair& p : gotPairs) {
            if (theirs.count(p) != 0) {
                commonMine.push_back(p);
            }
        }
        for (const Pair& p : expectedPairs) {
            if (mine.count(p) != 0) {
                commonTheirs.push_back(p);
            }
        }
        if (commonMine != commonTheirs) {
            ++counts.orderDiffers;
            std::printf("  order differs: %s\n", name.c_str());
        }
        ++(bandOnly ? counts.band : counts.tolerance);
    }
    std::printf("periodic cases: %d exact, %d within tolerance, %d band differences, %d documented "
                "deviation(s), %d failed (%d with another order)\n",
                counts.exact, counts.tolerance, counts.band, counts.deviations, counts.failed,
                counts.orderDiffers);
    for (const std::string& f : counts.failures) {
        std::printf("  failed: %s\n", f.c_str());
    }
}
