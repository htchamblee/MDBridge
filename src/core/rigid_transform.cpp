#include "mdbridge/core/rigid_transform.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "mdbridge/core/python.hpp"

// Compiled without floating-point contraction (see CMakeLists.txt): every product and sum is
// rounded as numpy rounds it.

namespace mdbridge::core {

namespace {

// Singular values of a 3x3 matrix: one-sided Jacobi (Hestenes) rotations of column pairs until
// the columns are orthogonal to working precision; the column norms are then the singular values.
Vec3 singularValues(const Mat3& a) {
    Mat3 b = a;
    constexpr double eps = 1e-17;
    for (int sweep = 0; sweep < 64; ++sweep) {
        bool rotated = false;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                double alpha = 0, beta = 0, gamma = 0;
                for (int i = 0; i < 3; ++i) {
                    alpha += b[i][p] * b[i][p];
                    beta += b[i][q] * b[i][q];
                    gamma += b[i][p] * b[i][q];
                }
                if (gamma == 0.0 || std::abs(gamma) <= eps * std::sqrt(alpha * beta)) {
                    continue;
                }
                rotated = true;
                const double zeta = (beta - alpha) / (2.0 * gamma);
                const double t =
                    (zeta >= 0 ? 1.0 : -1.0) / (std::abs(zeta) + std::hypot(1.0, zeta));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = c * t;
                for (int i = 0; i < 3; ++i) {
                    const double bp = b[i][p], bq = b[i][q];
                    b[i][p] = c * bp - s * bq;
                    b[i][q] = s * bp + c * bq;
                }
            }
        }
        if (!rotated) {
            break;
        }
    }
    Vec3 out{};
    for (int k = 0; k < 3; ++k) {
        out[static_cast<std::size_t>(k)] =
            std::sqrt(b[0][k] * b[0][k] + b[1][k] * b[1][k] + b[2][k] * b[2][k]);
    }
    return out;
}

// The orthogonal polar factor U V^T of a = U S V^T (what MDATK takes from numpy's SVD), by Newton's
// iteration x <- (x + x^-T) / 2: for singular values within 1 % of 1 it converges in four or five
// steps to within an ulp of the exact factor (numpy's LAPACK result is within a few 1e-15 of it).
Mat3 polarFactor(const Mat3& a) {
    Mat3 x = a;
    for (int step = 0; step < 12; ++step) {
        Mat3 cof{};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                const int i1 = (i + 1) % 3, i2 = (i + 2) % 3, j1 = (j + 1) % 3, j2 = (j + 2) % 3;
                cof[i][j] = x[i1][j1] * x[i2][j2] - x[i1][j2] * x[i2][j1];
            }
        }
        const double det = x[0][0] * cof[0][0] + x[0][1] * cof[0][1] + x[0][2] * cof[0][2];
        Mat3 next{};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                next[i][j] = (x[i][j] + cof[i][j] / det) * 0.5;  // x^-T = cofactors / det
            }
        }
        if (next == x) {
            break;
        }
        x = next;
    }
    return x;
}

double det3(const Mat3& r) {
    return r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) -
           r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) +
           r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
}

Mat3 matmul(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    return out;
}

Vec3 matvec(const Mat3& a, const Vec3& x) {
    return {a[0][0] * x[0] + a[0][1] * x[1] + a[0][2] * x[2],
            a[1][0] * x[0] + a[1][1] * x[1] + a[1][2] * x[2],
            a[2][0] * x[0] + a[2][1] * x[1] + a[2][2] * x[2]};
}

Mat3 transposed(const Mat3& a) {
    Mat3 out{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out[i][j] = a[j][i];
        }
    }
    return out;
}

}  // namespace

RigidTransform RigidTransform::fromList(std::span<const double> numbers) {
    if (numbers.size() != 12) {
        throw ValueError("a rigid transform needs 12 numbers (3x3 rotation + translation)");
    }
    if (!std::all_of(numbers.begin(), numbers.end(), [](double x) { return std::isfinite(x); })) {
        throw ValueError("the transform contains non-finite numbers");
    }
    Mat3 r{};
    for (std::size_t i = 0; i < 9; ++i) {
        r[i / 3][i % 3] = numbers[i];
    }
    const char* notRotation = "the rotation part of the transform is not a rotation";
    // every entry is bounded by the largest singular value; this also keeps the SVD finite
    for (const auto& row : r) {
        for (const double x : row) {
            if (std::abs(x) > 1.01) {
                throw ValueError(notRotation);
            }
        }
    }
    for (const double s : singularValues(r)) {
        if (std::abs(s - 1.0) > 0.01) {
            throw ValueError(notRotation);
        }
    }
    RigidTransform t;
    t.rotation = polarFactor(r);  // re-orthonormalise rounded values
    if (det3(t.rotation) < 0) {
        throw ValueError("the rotation part of the transform is a reflection");
    }
    t.translation = {numbers[9], numbers[10], numbers[11]};
    return t;
}

std::vector<double> RigidTransform::toList(int digits) const {
    std::vector<double> out;
    out.reserve(12);
    for (const auto& row : rotation) {
        for (const double x : row) {
            out.push_back(pythonRound(x, digits));
        }
    }
    for (const double x : translation) {
        out.push_back(pythonRound(x, digits));
    }
    return out;
}

Vec3 RigidTransform::apply(const Vec3& x) const {
    const Vec3 r = matvec(rotation, x);
    return {r[0] + translation[0], r[1] + translation[1], r[2] + translation[2]};
}

std::vector<Vec3> RigidTransform::apply(std::span<const Vec3> points) const {
    std::vector<Vec3> out;
    out.reserve(points.size());
    for (const Vec3& p : points) {
        out.push_back(apply(p));
    }
    return out;
}

RigidTransform RigidTransform::inverse() const {
    RigidTransform out;
    out.rotation = transposed(rotation);
    Mat3 negated = out.rotation;
    for (auto& row : negated) {
        for (double& x : row) {
            x = -x;
        }
    }
    out.translation = matvec(negated, translation);
    return out;
}

RigidTransform RigidTransform::compose(const RigidTransform& other) const {
    RigidTransform out;
    out.rotation = matmul(rotation, other.rotation);
    const Vec3 r = matvec(rotation, other.translation);
    out.translation = {r[0] + translation[0], r[1] + translation[1], r[2] + translation[2]};
    return out;
}

RigidTransform RigidTransform::about(const Vec3& centre) const {
    RigidTransform out;
    out.rotation = rotation;
    const Vec3 rc = matvec(rotation, centre);
    for (std::size_t k = 0; k < 3; ++k) {
        out.translation[k] = centre[k] - rc[k] + translation[k];
    }
    return out;
}

Mat3 rotvecMatrix(const Vec3& w) {
    const double theta = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    Mat3 out{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    if (theta < 1e-12) {
        return out;
    }
    const Vec3 k{w[0] / theta, w[1] / theta, w[2] / theta};
    const Mat3 kx{{{0, -k[2], k[1]}, {k[2], 0, -k[0]}, {-k[1], k[0], 0}}};
    const Mat3 kx2 = matmul(kx, kx);
    const double s = std::sin(theta), c = 1 - std::cos(theta);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out[i][j] = out[i][j] + s * kx[i][j] + c * kx2[i][j];
        }
    }
    return out;
}

Mat3 quaternionMatrix(const std::array<double, 4>& q) {
    const double norm = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    const double a = q[0] / norm, b = q[1] / norm, c = q[2] / norm, d = q[3] / norm;
    return {{{a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)},
             {2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)},
             {2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d}}};
}

double rotationAngle(const Mat3& r) {
    const double c = (r[0][0] + r[1][1] + r[2][2] - 1.0) / 2.0;
    // Python's min(1.0, max(-1.0, c)): NaN gives -1
    const double lower = c > -1.0 ? c : -1.0;
    const double clamped = lower < 1.0 ? lower : 1.0;
    return std::acos(clamped) * (180.0 / std::numbers::pi);
}

std::vector<Vec3> anchorPoints(std::span<const Vec3> points) {
    std::vector<Vec3> out;
    out.reserve(points.size());
    for (const Vec3& p : points) {
        out.push_back({pythonRound(p[0], 3), pythonRound(p[1], 3), pythonRound(p[2], 3)});
    }
    return out;
}

}  // namespace mdbridge::core
