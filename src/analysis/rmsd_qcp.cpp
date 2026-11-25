#include "mdbridge/analysis/rmsd_qcp.hpp"

// Written from the behaviour specification docs/specs/m3-rmsd.md (clean room) and the published
// QCP mathematics; every operation is rounded to f64 in the spec's order. This file must be
// compiled without floating-point contraction (-ffp-contract=off, /fp:precise; CMakeLists.txt).
//
// Sums of weights use numpy's pairwise summation (core/detail/numpy_sum.hpp).

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "mdbridge/core/detail/numpy_sum.hpp"

namespace mdbridge::analysis {

namespace {

constexpr std::size_t kMaxNewtonSteps = 50;
constexpr double kNewtonTolerance = 1e-14;
// AN-16: a centred set whose weighted inner-product matrix has q <= kCollinear * p^2 (p its trace,
// q the sum of its principal 2 x 2 minors; equivalently the second-largest eigenvalue is within
// about 1e-10 of the largest) has no unique optimal rotation.
constexpr double kCollinear = 1e-10;
// AN-17: E0 - lambda below kCancellation * eps * E0 is dominated by rounding: the iteration stops
// at a relative step of 1e-14 (about 45 eps) and E0 - lambda cancels; identical structures were
// seen at 250 eps (the register's 64 eps catches too few of them).
constexpr double kCancellation = 4096.0;

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

// AN-16: whether the centred set x (weights w, none: 1) is collinear or a single point: the
// weighted inner-product matrix M = sum w x x^T has (almost) one non-zero eigenvalue.
bool degenerate(const std::vector<double>& x, std::size_t n, const double* w) {
    if (n < 3) {
        return true;
    }
    double m[3][3] = {};
    for (std::size_t i = 0; i < n; ++i) {
        const double wi = w ? w[i] : 1.0;
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                m[r][c] += wi * (x[3 * i + r] * x[3 * i + c]);
            }
        }
    }
    const double p = m[0][0] + m[1][1] + m[2][2];
    const double q = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) + (m[0][0] * m[2][2] - m[0][2] * m[2][0]) +
                     (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
    return q <= kCollinear * p * p;
}

// The eigenvector of the largest eigenvalue of the symmetric 4 x 4 matrix k (cyclic Jacobi
// rotations; ties go to the lowest index). Rows of k already diagonal are never touched, so an
// exactly diagonal first row gives exactly (1, 0, 0, 0).
std::array<double, 4> topEigenvector(std::array<std::array<double, 4>, 4> k) {
    std::array<std::array<double, 4>, 4> v{};
    for (std::size_t i = 0; i < 4; ++i) {
        v[i][i] = 1.0;
    }
    for (int sweep = 0; sweep < 64; ++sweep) {
        double off = 0.0;
        for (std::size_t p = 0; p < 4; ++p) {
            for (std::size_t q = p + 1; q < 4; ++q) {
                off += k[p][q] * k[p][q];
            }
        }
        if (off == 0.0) {
            break;
        }
        for (std::size_t p = 0; p < 4; ++p) {
            for (std::size_t q = p + 1; q < 4; ++q) {
                if (k[p][q] == 0.0) {
                    continue;
                }
                const double theta = (k[q][q] - k[p][p]) / (2.0 * k[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (std::size_t r = 0; r < 4; ++r) {  // k <- J^T k J
                    const double krp = k[r][p];
                    const double krq = k[r][q];
                    k[r][p] = c * krp - s * krq;
                    k[r][q] = s * krp + c * krq;
                }
                for (std::size_t r = 0; r < 4; ++r) {
                    const double kpr = k[p][r];
                    const double kqr = k[q][r];
                    k[p][r] = c * kpr - s * kqr;
                    k[q][r] = s * kpr + c * kqr;
                }
                k[p][q] = 0.0;
                k[q][p] = 0.0;
                for (std::size_t r = 0; r < 4; ++r) {  // v <- v J
                    const double vrp = v[r][p];
                    const double vrq = v[r][q];
                    v[r][p] = c * vrp - s * vrq;
                    v[r][q] = s * vrp + c * vrq;
                }
            }
        }
    }
    std::size_t best = 0;
    for (std::size_t i = 1; i < 4; ++i) {
        if (k[i][i] > k[best][best]) {
            best = i;
        }
    }
    std::array<double, 4> q{v[0][best], v[1][best], v[2][best], v[3][best]};
    const double norm = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    for (double& x : q) {
        x /= norm;
    }
    return q;
}

// AN-16 / AN-17: the exact superposition RMSD of centred a (moved) onto centred b: the optimal
// rotation from the top eigenvector of Horn's key matrix, then the residual summed directly after
// rotating (never E0 - lambda). One atom and identical centred sets give exactly 0.
double exactRmsd(const std::vector<double>& a, const std::vector<double>& b, std::size_t n, const double* w) {
    if (n == 1 || a == b) {
        return 0.0;  // one atom, identical structures: exactly 0
    }
    double m[3][3] = {};  // sum w a b^T
    for (std::size_t i = 0; i < n; ++i) {
        const double wi = w ? w[i] : 1.0;
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                m[r][c] += wi * (a[3 * i + r] * b[3 * i + c]);  // symmetric for a == b
            }
        }
    }
    const std::array<std::array<double, 4>, 4> k{{
        {m[0][0] + m[1][1] + m[2][2], m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0]},
        {m[1][2] - m[2][1], m[0][0] - m[1][1] - m[2][2], m[0][1] + m[1][0], m[2][0] + m[0][2]},
        {m[2][0] - m[0][2], m[0][1] + m[1][0], -m[0][0] + m[1][1] - m[2][2], m[1][2] + m[2][1]},
        {m[0][1] - m[1][0], m[2][0] + m[0][2], m[1][2] + m[2][1], -m[0][0] - m[1][1] + m[2][2]},
    }};
    const auto [q0, q1, q2, q3] = topEigenvector(k);
    const double r[3][3] = {
        {q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q1 * q3 + q0 * q2)},
        {2.0 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2.0 * (q2 * q3 - q0 * q1)},
        {2.0 * (q1 * q3 - q0 * q2), 2.0 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3},
    };
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double d2 = 0.0;
        for (std::size_t c = 0; c < 3; ++c) {
            const double moved = (r[c][0] * a[3 * i] + r[c][1] * a[3 * i + 1]) + r[c][2] * a[3 * i + 2];
            const double d = moved - b[3 * i + c];
            d2 += d * d;
        }
        sum += (w ? w[i] : 1.0) * d2;
    }
    return std::sqrt(sum / static_cast<double>(n));
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
    bool converged = false;
    for (std::size_t step = 0; step < kMaxNewtonSteps; ++step) {
        const double old = lambda;
        const double x2 = lambda * lambda;
        const double beta = (x2 + c2) * lambda;
        const double alpha = beta + c1;
        const double delta = (alpha * lambda + c0) / (2.0 * x2 * lambda + beta + alpha);
        lambda -= delta;
        if (std::fabs(lambda - old) < std::fabs(kNewtonTolerance * lambda)) {
            converged = true;
            break;
        }
    }
    const double qcp = std::sqrt(std::fabs(2.0 * (e0 - lambda) / static_cast<double>(n)));
    if (n == 0) {
        return qcp;  // NaN
    }
    for (std::size_t i = 0; i < 3 * n; ++i) {
        if (!std::isfinite(a[i]) || !std::isfinite(b[i])) {
            return qcp;  // non-finite input: NaN as before
        }
    }
    if (w) {
        for (std::size_t i = 0; i < n; ++i) {
            if (!(w[i] >= 0.0) || !std::isfinite(w[i])) {
                return qcp;  // negative or non-finite weights: no minimum to speak of; MDAnalysis' value
            }
        }
    }
    const bool cancelled = e0 - lambda < kCancellation * std::numeric_limits<double>::epsilon() * e0;
    if (std::isfinite(qcp) && converged && !cancelled && !degenerate(a, n, w) && !degenerate(b, n, w)) {
        return qcp;  // the QCP value, bit for bit
    }
    return exactRmsd(a, b, n, w);
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
