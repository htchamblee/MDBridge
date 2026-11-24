#pragma once
// numpy's sums of float32 and float64 arrays, rounded as numpy rounds them: np.add.reduce of a
// contiguous 1-D array (np.sum, the sum of np.mean) is the identity 0 plus numpy's pairwise sum.
// Results match numpy bit for bit only when the caller's translation unit is compiled without
// floating-point contraction (-ffp-contract=off, MSVC /fp:precise), as the analysis and io code is.
//
// pairwiseSum is a port of numpy's pairwise summation (@TYPE@_pairwise_sum,
// numpy/_core/src/umath/loops_utils.h.src, numpy 2.5.5), used under numpy's licence:
//
// Copyright (c) 2005-2025, NumPy Developers.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright notice, this list of
//       conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright notice, this list of
//       conditions and the following disclaimer in the documentation and/or other materials
//       provided with the distribution.
//
//     * Neither the name of the NumPy Developers nor the names of any contributors may be used to
//       endorse or promote products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <cstddef>
#include <span>

namespace mdbridge::core::detail {

// numpy's pairwise summation of n values of type T (float or double) starting at `first`;
// value(i) gives element i. Blocks of up to 128 values are summed with 8 accumulators.
template <class T, class F>
T pairwiseSum(const F& value, std::size_t first, std::size_t n) {
    if (n < 8) {
        T res = T(-0.0);  // keeps the sign of a sum of negative zeros
        for (std::size_t i = 0; i < n; ++i) {
            res += value(first + i);
        }
        return res;
    }
    if (n <= 128) {
        T r[8];
        for (std::size_t j = 0; j < 8; ++j) {
            r[j] = value(first + j);
        }
        std::size_t i = 8;
        for (; i < n - n % 8; i += 8) {
            for (std::size_t j = 0; j < 8; ++j) {
                r[j] += value(first + i + j);
            }
        }
        T res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) {
            res += value(first + i);
        }
        return res;
    }
    std::size_t n2 = n / 2;
    n2 -= n2 % 8;
    return pairwiseSum<T>(value, first, n2) + pairwiseSum<T>(value, first + n2, n - n2);
}

// np.add.reduce of a contiguous array of n values: the identity 0 plus the pairwise sum
template <class T, class F>
T numpySum(const F& value, std::size_t n) {
    return T(0) + pairwiseSum<T>(value, 0, n);
}

inline double numpySum(std::span<const double> values) {
    const double* d = values.data();
    return numpySum<double>([d](std::size_t i) { return d[i]; }, values.size());
}

inline float numpySum(std::span<const float> values) {
    const float* d = values.data();
    return numpySum<float>([d](std::size_t i) { return d[i]; }, values.size());
}

// np.mean of a contiguous float64 array: the sum divided by the count (NaN when empty, where numpy
// also warns "Mean of empty slice")
inline double numpyMean(std::span<const double> values) {
    return numpySum(values) / static_cast<double>(values.size());
}

}  // namespace mdbridge::core::detail
