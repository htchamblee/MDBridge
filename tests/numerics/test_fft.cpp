// rfftn / irfftn (src/numerics/fft.cpp, the ducc0 FFT SciPy vendors) against scipy.fft 1.18.1:
// bit for bit, for any thread count (data/numerics/numerics_cases.json "fft", "next_fast_len").

#include <complex>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/hash.hpp"
#include "mdbridge/numerics/fft.hpp"
#include "numerics_cases.hpp"

using namespace mdbridge;
using namespace mdbridge::test::numerics_cases;

namespace {

template <class T>
std::string bytesSha(const std::vector<T>& v) {  // the little-endian bytes numpy hashes
    return core::sha256Hex(std::string_view(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T)));
}

std::vector<double> asDoubles(const std::vector<std::complex<double>>& c) {
    std::vector<double> out;
    out.reserve(2 * c.size());
    for (const std::complex<double>& z : c) {
        out.push_back(z.real());
        out.push_back(z.imag());
    }
    return out;
}

bool sameBits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

}  // namespace

TEST(Fft, MatchesScipyBitForBit) {
    for (const core::Value& c : cases("fft")) {
        const std::vector<std::size_t> s = counts(at(c, "shape"));
        const numerics::FftShape3 shape{s[0], s[1], s[2]};
        const std::size_t n = s[0] * s[1] * s[2];
        std::vector<double> a(n), b(n);
        for (std::size_t i = 0; i < n; ++i) {
            a[i] = value(i);
            b[i] = value(i + 7);
        }
        const std::string label = std::to_string(s[0]) + "x" + std::to_string(s[1]) + "x" + std::to_string(s[2]);
        const auto spec = numerics::rfftn(a, shape);
        EXPECT_EQ(bytesSha(spec), str(at(c, "rfftn_sha256"))) << label;
        const auto back = numerics::irfftn(spec, shape);
        EXPECT_EQ(bytesSha(back), str(at(c, "irfftn_sha256"))) << label;
        // MDATK's correlation product: irfftn(conj(rfftn(b)) * rfftn(a)), numpy's complex product
        const auto specB = numerics::rfftn(b, shape);
        std::vector<std::complex<double>> prod(spec.size());
        for (std::size_t k = 0; k < spec.size(); ++k) {
            const double ar = specB[k].real(), ai = -specB[k].imag(), br = spec[k].real(), bi = spec[k].imag();
            prod[k] = {ar * br - ai * bi, ar * bi + ai * br};
        }
        EXPECT_EQ(bytesSha(numerics::irfftn(prod, shape)), str(at(c, "product_sha256"))) << label;
        std::vector<std::complex<double>> half(spec.size());
        for (std::size_t k = 0; k < half.size(); ++k) {
            half[k] = {value(k + 3), value(k + 11)};
        }
        const auto raw = numerics::irfftn(half, shape);
        EXPECT_EQ(bytesSha(raw), str(at(c, "raw_sha256"))) << label;
        if (has(c, "rfftn")) {
            EXPECT_TRUE(sameBits(asDoubles(spec), nums(at(c, "rfftn")))) << label;
            EXPECT_TRUE(sameBits(back, nums(at(c, "irfftn")))) << label;
            EXPECT_TRUE(sameBits(raw, nums(at(c, "raw")))) << label;
        }
        // threads share out whole lines: the same bits for any count
        for (const std::size_t threads : {std::size_t{2}, std::size_t{3}, std::size_t{8}}) {
            EXPECT_TRUE(sameBits(asDoubles(numerics::rfftn(a, shape, threads)), asDoubles(spec))) << label;
            EXPECT_TRUE(sameBits(numerics::irfftn(prod, shape, threads), numerics::irfftn(prod, shape))) << label;
        }
    }
}

TEST(Fft, NextFastLen) {
    for (const core::Value& c : cases("next_fast_len")) {
        EXPECT_EQ(numerics::nextFastLenReal(count(at(c, "n"))), count(at(c, "out"))) << count(at(c, "n"));
    }
}

TEST(Fft, RejectsMismatchedShapes) {
    std::vector<double> a(10);
    EXPECT_THROW(numerics::rfftn(a, {2, 2, 2}), std::invalid_argument);
    EXPECT_THROW(numerics::rfftn(std::vector<double>{}, {0, 1, 1}), std::invalid_argument);
    std::vector<std::complex<double>> c(5);
    EXPECT_THROW(numerics::irfftn(c, {2, 2, 2}), std::invalid_argument);
}
