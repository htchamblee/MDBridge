// Gaussian filtering (mode "constant") and order-1 interpolation (mode "nearest") as
// scipy.ndimage computes them, ported from SciPy 1.18.1: _gaussian_kernel1d / gaussian_filter1d /
// gaussian_filter (scipy/ndimage/_filters.py), NI_Correlate1D and its line buffers
// (scipy/ndimage/src/ni_filters.c, ni_support.c), NI_GeometricTransform's interpolation (its
// matrix path for affine_transform, scipy/ndimage/_interpolation.py) and map_coordinate
// (src/ni_interpolation.c) and get_spline_interpolation_weights (src/ni_splines.c),
// used under SciPy's licence and the licence of its ndimage package:
//
// Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
// All rights reserved.
//
// Copyright (C) 2003-2005 Peter J. Verveer
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors (for the ndimage
//    package: the name of the author) may be used to endorse or promote products derived from this
//    software without specific prior written permission.
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
//
// Only the paths MDATK uses are ported (order-0 kernels, which are symmetric; order-1
// interpolation in mode "nearest" for map_coordinates and mode "constant" for affine_transform). The
// kernel's normalising sum is numpy's pairwise sum (core/detail/numpy_sum.hpp). Built without
// floating-point contraction, so that every product and sum rounds as in scipy.

#include "mdbridge/render/molecule/ndimage.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "mdbridge/core/detail/numpy_sum.hpp"

namespace mdbridge::render::molecule {

std::vector<double> gaussianKernel1d(double sigma, std::size_t radius) {
    const double sigma2 = sigma * sigma;
    const double factor = -0.5 / sigma2;
    const std::size_t n = 2 * radius + 1;
    std::vector<double> phi(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto x = static_cast<std::int64_t>(i) - static_cast<std::int64_t>(radius);
        phi[i] = std::exp(factor * static_cast<double>(x * x));
    }
    const double total = core::detail::numpySum(std::span<const double>(phi));
    for (double& v : phi) {
        v = v / total;
    }
    return phi;
}

void gaussianFilter1d(std::vector<float>& data, const Shape3& shape, std::size_t axis, double sigma, double truncate,
                      double cval) {
    if (axis > 2) {
        throw std::invalid_argument("gaussianFilter1d: axis must be 0, 1 or 2");
    }
    if (data.size() != shape[0] * shape[1] * shape[2]) {
        throw std::invalid_argument("gaussianFilter1d: data does not match the shape");
    }
    const double lw_real = truncate * sigma + 0.5;
    if (!(lw_real >= 0.0) || !std::isfinite(lw_real)) {
        throw std::invalid_argument("Radius must be a nonnegative integer.");
    }
    const auto radius = static_cast<std::size_t>(lw_real);  // int(): truncation toward zero
    // correlate with the reversed kernel; an order-0 kernel is symmetric
    const std::vector<double> weights = gaussianKernel1d(sigma, radius);
    const double* fw = weights.data() + radius;  // fw[-radius .. radius]
    const std::size_t length = shape[axis];
    const std::size_t stride = axis == 0 ? shape[1] * shape[2] : axis == 1 ? shape[2] : 1;
    const std::size_t lines = length == 0 ? 0 : data.size() / length;
    std::vector<double> line(length + 2 * radius, cval);
    const auto size1 = static_cast<std::ptrdiff_t>(radius);
    for (std::size_t l = 0; l < lines; ++l) {
        // the l-th line along `axis`, lines in C order of the other axes
        std::size_t start = 0;
        if (axis == 0) {
            start = l;
        } else if (axis == 1) {
            start = (l / shape[2]) * shape[1] * shape[2] + l % shape[2];
        } else {
            start = l * shape[2];
        }
        for (std::size_t i = 0; i < length; ++i) {
            line[radius + i] = static_cast<double>(data[start + i * stride]);
        }
        const double* in = line.data() + radius;
        for (std::size_t i = 0; i < length; ++i) {
            const double* p = in + i;
            double o = p[0] * fw[0];
            for (std::ptrdiff_t jj = -size1; jj < 0; ++jj) {
                o += (p[jj] + p[-jj]) * fw[jj];
            }
            data[start + i * stride] = static_cast<float>(o);
        }
    }
}

void gaussianFilter(std::vector<float>& data, const Shape3& shape, double sigma, double truncate, double cval) {
    if (!(sigma > 1e-15)) {
        return;
    }
    for (std::size_t axis = 0; axis < 3; ++axis) {
        gaussianFilter1d(data, shape, axis, sigma, truncate, cval);
    }
}

namespace detail {

LinearAxis linearAxis(double coordinate, std::size_t length) {
    const double fl = std::floor(coordinate);
    const auto start = static_cast<std::ptrdiff_t>(fl);
    LinearAxis out{};
    const auto last = static_cast<std::ptrdiff_t>(length) - 1;
    for (int k = 0; k < 2; ++k) {
        std::ptrdiff_t idx = start + k;
        if (idx < 0) {  // mode "nearest": the border value continues
            idx = 0;
        } else if (idx > last) {
            idx = last;
        }
        out.index[k] = idx;
    }
    const double x = coordinate - fl;
    out.weight[0] = 1.0 - x;
    out.weight[1] = 1.0 - out.weight[0];  // all weights add to 1
    return out;
}

bool linearAxisConstant(double coordinate, std::size_t length, LinearAxis& out) {
    // map_coordinate(): mode "constant" maps a coordinate outside [0, length - 1] to -1, which (as
    // NaN does) fails the `cc > -1.0` test: the point is cval
    double cc = coordinate;
    const auto len = static_cast<double>(length);
    if (cc < 0) {
        cc = -1;
    } else if (cc > len - 1) {
        cc = -1;
    }
    if (!(cc > -1.0)) {
        return false;
    }
    const double fl = std::floor(cc);
    const auto start = static_cast<std::ptrdiff_t>(fl);
    const auto n = static_cast<std::ptrdiff_t>(length);
    for (int k = 0; k < 2; ++k) {
        std::ptrdiff_t idx = start + k;
        if (start < 0 || start + 1 >= n) {
            // the spline boundary of mode "constant" is "mirror" (only idx == length is reached)
            if (n <= 1) {
                idx = 0;
            } else if (idx > n - 1) {
                const std::ptrdiff_t sz2 = 2 * n - 2;
                idx -= sz2 * (idx / sz2);
                if (idx >= n) {
                    idx = sz2 - idx;
                }
            }
        }
        out.index[k] = idx;
    }
    const double x = cc - fl;
    out.weight[0] = 1.0 - x;
    out.weight[1] = 1.0 - out.weight[0];
    return true;
}

}  // namespace detail

double affineSampleLinearConstant(std::span<const float> data, const Shape3& shape, const Matrix3& matrix,
                                  const std::array<double, 3>& offset, std::size_t i, std::size_t j, std::size_t k,
                                  double cval) {
    const double o[3] = {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)};
    detail::LinearAxis axes[3];
    for (std::size_t h = 0; h < 3; ++h) {
        double tmp = offset[h];
        tmp += o[0] * matrix[h][0];
        tmp += o[1] * matrix[h][1];
        tmp += o[2] * matrix[h][2];
        if (!detail::linearAxisConstant(tmp, shape[h], axes[h])) {
            return cval;
        }
    }
    const auto ny = static_cast<std::ptrdiff_t>(shape[1]);
    const auto nx = static_cast<std::ptrdiff_t>(shape[2]);
    double t = 0.0;
    for (int fz = 0; fz < 2; ++fz) {  // the last axis varies fastest, as in scipy's filter footprint
        for (int fy = 0; fy < 2; ++fy) {
            for (int fx = 0; fx < 2; ++fx) {
                const std::ptrdiff_t idx = (axes[0].index[fz] * ny + axes[1].index[fy]) * nx + axes[2].index[fx];
                double coeff = static_cast<double>(data[static_cast<std::size_t>(idx)]);
                coeff *= axes[0].weight[fz];
                coeff *= axes[1].weight[fy];
                coeff *= axes[2].weight[fx];
                t += coeff;
            }
        }
    }
    return t;
}

std::vector<float> affineTransformLinearConstant(std::span<const float> data, const Shape3& shape, const Matrix3& matrix,
                                                 const std::array<double, 3>& offset, const Shape3& outputShape,
                                                 double cval) {
    if (data.size() != shape[0] * shape[1] * shape[2]) {
        throw std::invalid_argument("affineTransformLinearConstant: data does not match the shape");
    }
    std::vector<float> out(outputShape[0] * outputShape[1] * outputShape[2]);
    std::size_t n = 0;
    for (std::size_t i = 0; i < outputShape[0]; ++i) {
        for (std::size_t j = 0; j < outputShape[1]; ++j) {
            for (std::size_t k = 0; k < outputShape[2]; ++k) {
                out[n++] = static_cast<float>(affineSampleLinearConstant(data, shape, matrix, offset, i, j, k, cval));
            }
        }
    }
    return out;
}

}  // namespace mdbridge::render::molecule
