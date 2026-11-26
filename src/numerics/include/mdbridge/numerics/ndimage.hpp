#pragma once
// The pieces of scipy.ndimage that MDATK uses, on 3-D float32 arrays in C order ([z, y, x], shape
// {nz, ny, nx}), rounded as scipy rounds them: gaussian_filter1d / gaussian_filter with mode
// "constant" (the molecular surface) or "reflect" (scipy's default: the density fit's low-pass
// filters), one sigma or one per axis (separable correlation in double precision, each pass stored as
// float32), uniform_filter1d with mode "constant" (the software renderer's edge anti-aliasing),
// map_coordinates with order 1 and mode "nearest" (surface normals) or mode "constant" with a cval
// (mdatk/density.py's sample_map), and affine_transform with order 1 and mode "constant" (the density
// shadow). Moved here from render/molecule (milestone 6) so that chem's density fit can use it.
// Ported from SciPy 1.18.1
// (scipy/ndimage/_filters.py, src/ni_filters.c, src/ni_interpolation.c, src/ni_splines.c; BSD
// licence, see ndimage.cpp and THIRD_PARTY_NOTICES.md).

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace mdbridge::numerics {

using Shape3 = std::array<std::size_t, 3>;

// scipy.ndimage._filters._gaussian_kernel1d(sigma, 0, radius): exp(-x^2 / (2 sigma^2)) for x in
// -radius..radius, normalised by its (numpy pairwise) sum.
std::vector<double> gaussianKernel1d(double sigma, std::size_t radius);

// How a filter extends a line beyond its ends (scipy.ndimage's `mode`): "constant" (cval) or
// "reflect" (d c b a | a b c d | d c b a; the reflection repeats for kernels longer than the line).
enum class ExtendMode { Constant, Reflect };

// gaussian_filter1d(data, sigma, axis, mode, cval, truncate) in place: radius
// int(truncate * sigma + 0.5).
void gaussianFilter1d(std::vector<float>& data, const Shape3& shape, std::size_t axis, double sigma,
                      double truncate = 4.0, double cval = 0.0, ExtendMode mode = ExtendMode::Constant);

// gaussian_filter(data, sigma, mode="constant", cval, truncate) in place: the 1-D filter along
// axes 0, 1, 2 (none when sigma <= 1e-15).
void gaussianFilter(std::vector<float>& data, const Shape3& shape, double sigma, double truncate = 4.0,
                    double cval = 0.0);

// gaussian_filter(data, sigma=[s0, s1, s2], mode, cval, truncate) in place: the 1-D filter along each
// axis whose sigma is above 1e-15, axes in order (MDATK: sigma per axis = a length over the voxel
// size of that axis, mode "reflect").
void gaussianFilter(std::vector<float>& data, const Shape3& shape, const std::array<double, 3>& sigma,
                    ExtendMode mode, double truncate = 4.0, double cval = 0.0);

// uniform_filter1d(line, size, mode="constant", cval) of one line of `length` values `stride` apart:
// NI_UniformFilter1D's running sum in double, each value stored as float32 into `out` (same
// stride). size >= 1.
void uniformFilter1dLine(const float* in, std::size_t stride, std::size_t length, std::size_t size, float* out,
                         double cval = 0.0);

// map_coordinates(data, [[z], [y], [x]], order=1, mode="nearest") at one point (finite coordinates).
// `value(i)` gives the array value at flat index i, so the array need not exist (e.g. gradients).
template <class F>
float mapCoordinatesLinearNearest(const F& value, const Shape3& shape, double z, double y, double x);

// map_coordinates(data, [[z], [y], [x]], order=1, mode="constant", cval, prefilter=False) at one
// point, before the float32 rounding of a float32 input: cval when a coordinate is below 0, above
// length - 1 or NaN along any axis (checked in axis order), else the order-1 interpolation in double
// precision (NaN data propagate). MDATK's sample_map passes cval = NaN to mark points outside the map.
double mapCoordinatesLinearConstant(std::span<const float> data, const Shape3& shape, double z, double y, double x,
                                    double cval = 0.0);

// map_coordinates as above for n points (coordinates [z0, y0, x0, z1, ...]), each result stored as
// float32 as scipy stores them for a float32 input.
std::vector<float> mapCoordinatesLinearConstant(std::span<const float> data, const Shape3& shape,
                                                std::span<const double> zyx, double cval = 0.0);

// affine_transform(data, matrix, offset, output_shape, order=1, mode="constant", cval,
// prefilter=False) for a 3 x 3 matrix (row-major): output[o] = the order-1 interpolation of data at
// matrix @ o + offset (output index o as doubles, products added in order to the offset). A
// coordinate below 0 or above length - 1 along any axis (or NaN) gives cval; the result of each point
// is stored as float32, as for a float32 input with no output array.
using Matrix3 = std::array<std::array<double, 3>, 3>;
std::vector<float> affineTransformLinearConstant(std::span<const float> data, const Shape3& shape, const Matrix3& matrix,
                                                 const std::array<double, 3>& offset, const Shape3& outputShape,
                                                 double cval = 0.0);

// One output point of affineTransformLinearConstant (output index (i, j, k)), before the float32
// rounding: the interpolated value in double precision, or cval.
double affineSampleLinearConstant(std::span<const float> data, const Shape3& shape, const Matrix3& matrix,
                                  const std::array<double, 3>& offset, std::size_t i, std::size_t j, std::size_t k,
                                  double cval = 0.0);

// Implementation --------------------------------------------------------------------------------------

namespace detail {
// start index and the two spline weights of order 1 along one axis
struct LinearAxis {
    std::ptrdiff_t index[2];
    double weight[2];
};
LinearAxis linearAxis(double coordinate, std::size_t length);
// mode "constant": false (cval) when the coordinate is below 0, above length - 1 or NaN; taps beyond
// the last grid line (coordinate == length - 1) are mirrored, with weight 0
bool linearAxisConstant(double coordinate, std::size_t length, LinearAxis& out);
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

}  // namespace mdbridge::numerics
