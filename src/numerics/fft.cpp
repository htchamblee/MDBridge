// rfftn / irfftn on the ducc0 FFT that SciPy 1.18.1 vendors (third_party/ducc0, BSD-3-Clause; see
// third_party/ducc0/LICENSE.md and THIRD_PARTY_NOTICES.md), called as scipy/fft/_duccfft/pyduccfft.cxx
// calls it (r2c_internal / c2r_internal, norm_fct).

#include "mdbridge/numerics/fft.hpp"

#include <stdexcept>

// ducc0 is third-party code: its headers are on a system include path (no warnings with GCC) and,
// for MSVC, warnings are off around the include only (the project's own warning level is unchanged)
#if defined(_MSC_VER)
#pragma warning(push, 0)
#endif
#include "ducc0/fft/fft.h"
#include "ducc0/fft/fftnd_impl.h"
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

namespace mdbridge::numerics {

namespace {

std::size_t checkedThreads(std::size_t threads) {
    return threads == 0 ? 1 : threads;
}

}  // namespace

std::vector<std::complex<double>> rfftn(std::span<const double> a, const FftShape3& shape, std::size_t threads) {
    if (a.size() != shape[0] * shape[1] * shape[2]) {
        throw std::invalid_argument("rfftn: the data do not match the shape");
    }
    if (shape[0] == 0 || shape[1] == 0 || shape[2] == 0) {
        throw std::invalid_argument("rfftn: every axis needs at least one point");
    }
    const ducc0::fmav_info::shape_t in{shape[0], shape[1], shape[2]};
    const ducc0::fmav_info::shape_t outShape{shape[0], shape[1], shape[2] / 2 + 1};
    std::vector<std::complex<double>> out(outShape[0] * outShape[1] * outShape[2]);
    const ducc0::cfmav<double> src(a.data(), in);
    const ducc0::vfmav<std::complex<double>> dst(out.data(), outShape);
    ducc0::r2c(src, dst, ducc0::fmav_info::shape_t{0, 1, 2}, true, 1.0, checkedThreads(threads));
    return out;
}

std::vector<double> irfftn(std::span<const std::complex<double>> a, const FftShape3& shape, std::size_t threads) {
    if (shape[0] == 0 || shape[1] == 0 || shape[2] == 0) {
        throw std::invalid_argument("irfftn: every axis needs at least one point");
    }
    const ducc0::fmav_info::shape_t inShape{shape[0], shape[1], shape[2] / 2 + 1};
    if (a.size() != inShape[0] * inShape[1] * inShape[2]) {
        throw std::invalid_argument("irfftn: the half spectrum does not match the shape");
    }
    // norm_fct(inorm=2, N): T(1 / ldbl_t(N))
    const long double n = static_cast<long double>(shape[0]) * static_cast<long double>(shape[1]) *
                          static_cast<long double>(shape[2]);
    const auto fct = static_cast<double>(1.0L / n);
    std::vector<double> out(shape[0] * shape[1] * shape[2]);
    const ducc0::cfmav<std::complex<double>> src(a.data(), inShape);
    const ducc0::vfmav<double> dst(out.data(), ducc0::fmav_info::shape_t{shape[0], shape[1], shape[2]});
    ducc0::c2r(src, dst, ducc0::fmav_info::shape_t{0, 1, 2}, false, fct, checkedThreads(threads));
    return out;
}

std::size_t nextFastLenReal(std::size_t n) {
    return ducc0::good_size_real(n);
}

}  // namespace mdbridge::numerics
