// rfftn / irfftn (src/numerics/fft.cpp, the ducc0 FFT SciPy vendors) against scipy.fft 1.18.1:
// bit for bit with glibc, for any thread count (data/numerics/numerics_cases.json "fft",
// "next_fast_len"). ducc0 computes its twiddle factors with the C library's sin / cos
// (ducc0/math/unity_roots.h), so the bits are those of the generator's libm (glibc): with another C
// library (MSVC's: CI run 45882807857, every case with a radix above 5 or a long line) they move in
// the last places. Everywhere, the transforms are also compared with a direct (separable) DFT within
// a relative 1e-12 of the array's largest value, which is what stands in for the bits without glibc.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstring>
#include <numbers>
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

#if defined(__GLIBC__)
constexpr bool kGlibcLibm = true;  // the generator's SciPy ran on glibc's sin / cos
#else
constexpr bool kGlibcLibm = false;
#endif

constexpr double kDirectTol = 1e-12;

// exp(sign * 2 pi i k / n), k reduced modulo n first (accurate for any k)
std::complex<double> root(std::size_t k, std::size_t n, double sign) {
    const double angle = sign * 2.0 * std::numbers::pi * static_cast<double>(k % n) / static_cast<double>(n);
    return {std::cos(angle), std::sin(angle)};
}

// out[j] = sum_k in[k] exp(sign 2 pi i j k / n) along one axis of a C-order array of `dims`
void directAxis(std::vector<std::complex<double>>& v, const std::array<std::size_t, 3>& dims, int axis,
                double sign) {
    const std::size_t n = dims[static_cast<std::size_t>(axis)];
    const std::size_t stride = axis == 2 ? 1 : axis == 1 ? dims[2] : dims[1] * dims[2];
    std::vector<std::complex<double>> line(n), w(n);
    for (std::size_t k = 0; k < n; ++k) {
        w[k] = root(k, n, sign);
    }
    for (std::size_t i = 0; i < v.size(); ++i) {
        if ((i / stride) % n != 0) {
            continue;  // i is the first element of a line
        }
        for (std::size_t j = 0; j < n; ++j) {
            std::complex<double> sum = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                sum += v[i + k * stride] * w[(j * k) % n];
            }
            line[j] = sum;
        }
        for (std::size_t j = 0; j < n; ++j) {
            v[i + j * stride] = line[j];
        }
    }
}

// rfftn by the definition: the full complex DFT over the three axes, the half kept
std::vector<std::complex<double>> directRfftn(const std::vector<double>& a, const std::vector<std::size_t>& s) {
    const std::array<std::size_t, 3> dims{s[0], s[1], s[2]};
    std::vector<std::complex<double>> v(a.begin(), a.end());
    for (int axis = 0; axis < 3; ++axis) {
        directAxis(v, dims, axis, -1.0);
    }
    const std::size_t h = s[2] / 2 + 1;
    std::vector<std::complex<double>> out;
    out.reserve(s[0] * s[1] * h);
    for (std::size_t i = 0; i < s[0] * s[1]; ++i) {
        out.insert(out.end(), v.begin() + static_cast<std::ptrdiff_t>(i * s[2]),
                   v.begin() + static_cast<std::ptrdiff_t>(i * s[2] + h));
    }
    return out;
}

// irfftn(c, s) by the definition, as pocketfft / ducc0 read a half spectrum: inverse complex DFTs
// over the first two axes, then the real inverse along the last (the imaginary parts of bin 0 and,
// for even n2, bin n2 / 2 ignored), scaled by 1 / (n0 n1 n2)
std::vector<double> directIrfftn(const std::vector<std::complex<double>>& c, const std::vector<std::size_t>& s) {
    const std::size_t h = s[2] / 2 + 1;
    const std::array<std::size_t, 3> dims{s[0], s[1], h};
    std::vector<std::complex<double>> v(c);
    directAxis(v, dims, 0, 1.0);
    directAxis(v, dims, 1, 1.0);
    const std::size_t n = s[2];
    const double scale = 1.0 / static_cast<double>(s[0] * s[1] * s[2]);
    std::vector<std::complex<double>> w(n);
    for (std::size_t k = 0; k < n; ++k) {
        w[k] = root(k, n, 1.0);
    }
    std::vector<double> out(s[0] * s[1] * n);
    for (std::size_t i = 0; i < s[0] * s[1]; ++i) {
        const std::complex<double>* x = v.data() + i * h;
        for (std::size_t j = 0; j < n; ++j) {
            double sum = x[0].real();
            for (std::size_t k = 1; k < h; ++k) {
                const double term = (x[k] * w[(j * k) % n]).real();
                sum += (2 * k == n) ? term : 2.0 * term;
            }
            out[i * n + j] = sum * scale;
        }
    }
    return out;
}

double largest(const std::vector<double>& v) {
    double m = 0.0;
    for (const double x : v) {
        m = std::max(m, std::fabs(x));
    }
    return m;
}

// every value within kDirectTol of the reference's largest magnitude (at least 1)
bool closeScaled(const std::vector<double>& got, const std::vector<double>& want) {
    if (got.size() != want.size()) {
        return false;
    }
    const double tol = kDirectTol * std::max(1.0, largest(want));
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (!(std::fabs(got[i] - want[i]) <= tol)) {
            return false;
        }
    }
    return true;
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
        const auto back = numerics::irfftn(spec, shape);
        // MDATK's correlation product: irfftn(conj(rfftn(b)) * rfftn(a)), numpy's complex product
        const auto specB = numerics::rfftn(b, shape);
        std::vector<std::complex<double>> prod(spec.size());
        for (std::size_t k = 0; k < spec.size(); ++k) {
            const double ar = specB[k].real(), ai = -specB[k].imag(), br = spec[k].real(), bi = spec[k].imag();
            prod[k] = {ar * br - ai * bi, ar * bi + ai * br};
        }
        const auto product = numerics::irfftn(prod, shape);
        std::vector<std::complex<double>> half(spec.size());
        for (std::size_t k = 0; k < half.size(); ++k) {
            half[k] = {value(k + 3), value(k + 11)};
        }
        const auto raw = numerics::irfftn(half, shape);
        if constexpr (kGlibcLibm) {
            EXPECT_EQ(bytesSha(spec), str(at(c, "rfftn_sha256"))) << label;
            EXPECT_EQ(bytesSha(back), str(at(c, "irfftn_sha256"))) << label;
            EXPECT_EQ(bytesSha(product), str(at(c, "product_sha256"))) << label;
            EXPECT_EQ(bytesSha(raw), str(at(c, "raw_sha256"))) << label;
        }
        if (has(c, "rfftn")) {
            if constexpr (kGlibcLibm) {
                EXPECT_TRUE(sameBits(asDoubles(spec), nums(at(c, "rfftn")))) << label;
                EXPECT_TRUE(sameBits(back, nums(at(c, "irfftn")))) << label;
                EXPECT_TRUE(sameBits(raw, nums(at(c, "raw")))) << label;
            } else {
                EXPECT_TRUE(closeScaled(asDoubles(spec), nums(at(c, "rfftn")))) << label;
                EXPECT_TRUE(closeScaled(back, nums(at(c, "irfftn")))) << label;
                EXPECT_TRUE(closeScaled(raw, nums(at(c, "raw")))) << label;
            }
        }
        // the definition, on every C library
        EXPECT_TRUE(closeScaled(asDoubles(spec), asDoubles(directRfftn(a, s)))) << label;
        EXPECT_TRUE(closeScaled(back, a)) << label;
        EXPECT_TRUE(closeScaled(product, directIrfftn(prod, s))) << label;
        EXPECT_TRUE(closeScaled(raw, directIrfftn(half, s))) << label;
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
