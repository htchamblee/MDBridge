#pragma once
// The cases tests/parity/make_numerics_cases.py recorded from SciPy 1.18.1
// (data/numerics/numerics_cases.json) and small helpers to read them.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/value.hpp"
#include "test_support.hpp"

namespace mdbridge::test::numerics_cases {

inline const core::Value& all() {
    static const core::Value data = core::parseJson(
        readFile(kData / "numerics" / "numerics_cases.json"), "numerics_cases.json");
    return data;
}

inline const core::Value& at(const core::Value& v, std::string_view key) {
    const core::Value* found = v.find(key);
    if (found == nullptr) {
        throw std::runtime_error("numerics_cases.json: no key " + std::string(key));
    }
    return *found;
}

inline const core::Array& list(const core::Value& v) {
    return *v.get<core::Array>();
}

inline const core::Array& cases(std::string_view key) {
    return list(at(all(), key));
}

inline bool has(const core::Value& v, std::string_view key) {
    return v.find(key) != nullptr && !v.find(key)->isNull();
}

inline double num(const core::Value& v) {
    if (const double* d = v.get<double>()) {
        return *d;
    }
    return static_cast<double>(*v.get<std::int64_t>());
}

inline std::size_t count(const core::Value& v) {
    return static_cast<std::size_t>(*v.get<std::int64_t>());
}

inline std::vector<double> nums(const core::Value& v) {
    std::vector<double> out;
    for (const core::Value& e : list(v)) {
        out.push_back(num(e));
    }
    return out;
}

inline std::vector<std::size_t> counts(const core::Value& v) {
    std::vector<std::size_t> out;
    for (const core::Value& e : list(v)) {
        out.push_back(count(e));
    }
    return out;
}

inline std::string str(const core::Value& v) {
    return v.stringOr();
}

// value(i) of the generator: ((i * 2654435761) mod 2^32) / 2^32 - 0.5, exact in double
inline double value(std::uint64_t i) {
    const std::uint64_t k = (i * 2654435761ULL) % 4294967296ULL;
    return static_cast<double>(k) / 4294967296.0 - 0.5;
}

// |a - b| <= tol * max(1, |b|), or both NaN, or the same infinity
inline bool close(double a, double b, double tol) {
    if (std::isnan(a) || std::isnan(b)) {
        return std::isnan(a) && std::isnan(b);
    }
    if (std::isinf(a) || std::isinf(b)) {
        return a == b;
    }
    return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

inline void expectClose(std::span<const double> got, std::span<const double> want, double tol,
                        const std::string& what) {
    ASSERT_EQ(got.size(), want.size()) << what;
    for (std::size_t i = 0; i < got.size(); ++i) {
        EXPECT_TRUE(close(got[i], want[i], tol)) << what << "[" << i << "]: " << got[i] << " vs " << want[i];
    }
}

}  // namespace mdbridge::test::numerics_cases
