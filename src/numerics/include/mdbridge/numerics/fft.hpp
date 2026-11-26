#pragma once
// scipy.fft.rfftn / irfftn of 3-D float64 arrays in C order ([z, y, x], shape {n0, n1, n2}) with the
// default ("backward") normalisation, computed by the FFT scipy.fft itself runs: Martin Reinecke's ducc0
// FFT 0.40.0 as SciPy 1.18.1 vendors it (third_party/ducc0, BSD-3-Clause; owner decision D3). The
// results equal SciPy's bit for bit on the same CPU family (ducc0 vectorises over lines with the
// compiler's SIMD types where SciPy's build does, never inside a line), and do not depend on the
// thread count: threads only share out whole lines (MDATK passes workers=-1).

#include <array>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace mdbridge::numerics {

using FftShape3 = std::array<std::size_t, 3>;

// rfftn(a) (axes (0, 1, 2)): the half spectrum, shape {n0, n1, n2 / 2 + 1}, C order. `a` holds
// n0 * n1 * n2 values (std::invalid_argument otherwise); threads >= 1 (0 counts as 1).
std::vector<std::complex<double>> rfftn(std::span<const double> a, const FftShape3& shape, std::size_t threads = 1);

// irfftn(a, s=shape, axes=(0, 1, 2)) of a half spectrum of shape {n0, n1, n2 / 2 + 1}: the real
// array of `shape`, scaled by 1 / (n0 n1 n2) as scipy scales it (the factor rounded from long double,
// as SciPy's C++ computes it). The imaginary parts of the bins a real signal cannot have (0 and, for
// even n2, n2 / 2 along the last axis) are ignored as pocketfft / ducc0 ignore them.
std::vector<double> irfftn(std::span<const std::complex<double>> a, const FftShape3& shape, std::size_t threads = 1);

// scipy.fft.next_fast_len(n, real=True): n itself up to 6, else the smallest 2^a 3^b 5^c >= n
// (ducc0's good_size_real).
std::size_t nextFastLenReal(std::size_t n);

}  // namespace mdbridge::numerics
