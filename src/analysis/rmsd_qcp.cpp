#include "mdbridge/analysis/rmsd_qcp.hpp"

// Written from the behaviour specification docs/specs/m3-rmsd.md (clean room) and the published
// QCP mathematics; every operation is rounded to f64 in the spec's order. This file must be
// compiled without floating-point contraction (-ffp-contract=off, /fp:precise; CMakeLists.txt).
//
// Sums of weights use numpy's pairwise summation (core/detail/numpy_sum.hpp).

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "mdbridge/core/detail/numpy_sum.hpp"

namespace mdbridge::analysis {

namespace {

constexpr std::size_t kMaxNewtonSteps = 50;
constexpr double kNewtonTolerance = 1e-14;

// The coordinates minus their (weighted) centre, §3: sums in atom order, one division.
std::vector<double> centred(std::span<const double> x, std::size_t n, const double* w, double weightSum) {
    std::array<double, 3> centre{};
    for (std::size_t k = 0; k < 3; ++k) {
        double sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double term = w ? x[3 * i + k] * w[i] : x[3 * i + k];
            sum = i == 0 ? term : sum + term;
        }
        centre[k] = sum / (w ? weightSum : static_cast<double>(n));
    }
    std::vector<double> out(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            out[3 * i + k] = x[3 * i + k] - centre[k];
        }
    }
    return out;
}

// §4: QCP on centred coordinates a (moved) and b (reference), weights normalised to mean 1.
double qcpRmsd(const std::vector<double>& a, const std::vector<double>& b, std::size_t n, const double* w) {
    double gA = 0.0;
    double gB = 0.0;
    double sxx = 0.0, sxy = 0.0, sxz = 0.0;
    double syx = 0.0, syy = 0.0, syz = 0.0;
    double szx = 0.0, szy = 0.0, szz = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double ax = a[3 * i], ay = a[3 * i + 1], az = a[3 * i + 2];
        const double bx = b[3 * i], by = b[3 * i + 1], bz = b[3 * i + 2];
        double ux = bx, uy = by, uz = bz;
        if (w) {
            ux = w[i] * bx;
            uy = w[i] * by;
            uz = w[i] * bz;
            gA += w[i] * ((ax * ax + ay * ay) + az * az);
        } else {
            gA += (ax * ax + ay * ay) + az * az;
        }
        gB += (ux * bx + uy * by) + uz * bz;
        sxx += ux * ax;
        sxy += ux * ay;
        sxz += ux * az;
        syx += uy * ax;
        syy += uy * ay;
        syz += uy * az;
        szx += uz * ax;
        szy += uz * ay;
        szz += uz * az;
    }
    const double e0 = (gB + gA) * 0.5;

    const double sxx2 = sxx * sxx, syy2 = syy * syy, szz2 = szz * szz;
    const double sxy2 = sxy * sxy, syz2 = syz * syz, sxz2 = sxz * sxz;
    const double syx2 = syx * syx, szy2 = szy * szy, szx2 = szx * szx;

    const double c2 = -2.0 * (sxx2 + syy2 + szz2 + sxy2 + syx2 + sxz2 + szx2 + syz2 + szy2);
    const double c1 = 8.0 * (sxx * syz * szy + syy * szx * sxz + szz * sxy * syx - sxx * syy * szz -
                             syz * szx * sxy - szy * syx * sxz);

    const double d1 = sxy2 + sxz2 - syx2 - szx2;
    const double e = syy2 + szz2 - sxx2 + syz2 + szy2;
    const double f = 2.0 * (syz * szy - syy * szz);
    const double pxz = sxz + szx, pyz = syz + szy, pxy = sxy + syx;
    const double myz = syz - szy, mxz = sxz - szx, mxy = sxy - syx;
    const double pxx = sxx + syy, mxx = sxx - syy;

    const double t1 = d1 * d1;
    const double t2 = (e + f) * (e - f);
    const double t3 = ((-pxz) * myz + mxy * (mxx - szz)) * ((-mxz) * pyz + mxy * (mxx + szz));
    const double t4 = ((-pxz) * pyz - pxy * (pxx - szz)) * ((-mxz) * myz - pxy * (pxx + szz));
    const double t5 = (pxy * pyz + pxz * (mxx + szz)) * ((-mxy) * myz + pxz * (pxx + szz));
    const double t6 = (pxy * myz + mxz * (mxx - szz)) * ((-mxy) * pyz + mxz * (pxx - szz));
    const double c0 = t1 + t2 + t3 + t4 + t5 + t6;

    double lambda = e0;
    for (std::size_t step = 0; step < kMaxNewtonSteps; ++step) {
        const double old = lambda;
        const double x2 = lambda * lambda;
        const double beta = (x2 + c2) * lambda;
        const double alpha = beta + c1;
        const double delta = (alpha * lambda + c0) / (2.0 * x2 * lambda + beta + alpha);
        lambda -= delta;
        if (std::fabs(lambda - old) < std::fabs(kNewtonTolerance * lambda)) {
            break;
        }
    }
    return std::sqrt(std::fabs(2.0 * (e0 - lambda) / static_cast<double>(n)));
}

}  // namespace

double superposedRmsd(const ShapedArray& a, const ShapedArray& b, const std::optional<ShapedArray>& weights) {
    if (!std::equal(a.shape.begin(), a.shape.end(), b.shape.begin(), b.shape.end())) {
        throw core::ValueError("a and b must have same shape");
    }
    if (a.shape.size() != 2 || a.shape[1] != 3 || a.values.size() != 3 * a.shape[0] ||
        b.values.size() != a.values.size()) {
        throw std::invalid_argument("superposedRmsd: coordinates must be n x 3 arrays");
    }
    const std::size_t n = a.shape[0];
    if (!weights) {
        return qcpRmsd(centred(a.values, n, nullptr, 0.0), centred(b.values, n, nullptr, 0.0), n, nullptr);
    }
    if (weights->shape.size() != 1 || weights->shape[0] != n || weights->values.size() != n) {
        throw core::ValueError("Shape of weights must be consistent with shape of a along specified axis.");
    }
    const std::span<const double> w = weights->values;
    const double sum = core::detail::numpySum(w);
    if (sum == 0.0) {
        throw ZeroDivisionError("Weights sum to zero, can't be normalized");
    }
    std::vector<double> ac = centred(a.values, n, w.data(), sum);
    std::vector<double> bc = centred(b.values, n, w.data(), sum);
    const double mean = sum / static_cast<double>(n);
    std::vector<double> normalised(n);
    for (std::size_t i = 0; i < n; ++i) {
        normalised[i] = w[i] / mean;
    }
    return qcpRmsd(ac, bc, n, normalised.data());
}

double superposedRmsd(std::span<const double> a, std::span<const double> b,
                      std::optional<std::span<const double>> weights) {
    if (a.size() % 3 != 0 || b.size() % 3 != 0) {
        throw std::invalid_argument("superposedRmsd: coordinates must be n x 3 arrays");
    }
    const std::array<std::size_t, 2> aShape{a.size() / 3, 3};
    const std::array<std::size_t, 2> bShape{b.size() / 3, 3};
    std::optional<ShapedArray> w;
    std::array<std::size_t, 1> wShape{};
    if (weights) {
        wShape[0] = weights->size();
        w = ShapedArray{*weights, wShape};
    }
    return superposedRmsd(ShapedArray{a, aShape}, ShapedArray{b, bShape}, w);
}

}  // namespace mdbridge::analysis
