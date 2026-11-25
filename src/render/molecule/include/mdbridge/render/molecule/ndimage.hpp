#pragma once
// The pieces of scipy.ndimage that MDATK's molecular surface uses, on 3-D float32 arrays in C order
// ([z, y, x], shape {nz, ny, nx}), rounded as scipy rounds them: gaussian_filter1d /
// gaussian_filter with mode "constant" (separable correlation in double precision, each pass
// stored as float32) and map_coordinates with order 1 and mode "nearest". Ported from SciPy 1.18.1
// (scipy/ndimage/_filters.py, src/ni_filters.c, src/ni_interpolation.c, src/ni_splines.c; BSD
// licence, see ndimage.cpp and THIRD_PARTY_NOTICES.md).

#include <array>
#include <cstddef>
#include <vector>

namespace mdbridge::render::molecule {

using Shape3 = std::array<std::size_t, 3>;

// scipy.ndimage._filters._gaussian_kernel1d(sigma, 0, radius): exp(-x^2 / (2 sigma^2)) for x in
// -radius..radius, normalised by its (numpy pairwise) sum.
std::vector<double> gaussianKernel1d(double sigma, std::size_t radius);

// gaussian_filter1d(data, sigma, axis, mode="constant", cval, truncate) in place: radius
// int(truncate * sigma + 0.5).
void gaussianFilter1d(std::vector<float>& data, const Shape3& shape, std::size_t axis, double sigma,
                      double truncate = 4.0, double cval = 0.0);

// gaussian_filter(data, sigma, mode="constant", cval, truncate) in place: the 1-D filter along
// axes 0, 1, 2 (none when sigma <= 1e-15).
void gaussianFilter(std::vector<float>& data, const Shape3& shape, double sigma, double truncate = 4.0,
                    double cval = 0.0);

// map_coordinates(data, [[z], [y], [x]], order=1, mode="nearest") at one point (finite coordinates).
// `value(i)` gives the array value at flat index i, so the array need not exist (e.g. gradients).
template <class F>
float mapCoordinatesLinearNearest(const F& value, const Shape3& shape, double z, double y, double x);

// Implementation --------------------------------------------------------------------------------------

namespace detail {
// start index and the two spline weights of order 1 along one axis
struct LinearAxis {
    std::ptrdiff_t index[2];
    double weight[2];
};
LinearAxis linearAxis(double coordinate, std::size_t length);
}  // namespace detail

template <class F>
float mapCoordinatesLinearNearest(const F& value, const Shape3& shape, double z, double y, double x) {
    const detail::LinearAxis az = detail::linearAxis(z, shape[0]);
    const detail::LinearAxis ay = detail::linearAxis(y, shape[1]);
    const detail::LinearAxis ax = detail::linearAxis(x, shape[2]);
    const auto ny = static_cast<std::ptrdiff_t>(shape[1]);
    const auto nx = static_cast<std::ptrdiff_t>(shape[2]);
    double t = 0.0;
    for (int fz = 0; fz < 2; ++fz) {  // the last axis varies fastest, as in scipy's filter footprint
        for (int fy = 0; fy < 2; ++fy) {
            for (int fx = 0; fx < 2; ++fx) {
                const std::ptrdiff_t idx = (az.index[fz] * ny + ay.index[fy]) * nx + ax.index[fx];
                double coeff = static_cast<double>(value(static_cast<std::size_t>(idx)));
                coeff *= az.weight[fz];
                coeff *= ay.weight[fy];
                coeff *= ax.weight[fx];
                t += coeff;
            }
        }
    }
    return static_cast<float>(t);
}

}  // namespace mdbridge::render::molecule
