#include "mdbridge/core/linalg.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// Compiled without floating-point contraction (see CMakeLists.txt).

namespace mdbridge::core {

namespace {

// One-sided Jacobi (Hestenes): rotate column pairs of the m x 3 matrix b (row-major) until the
// columns are orthogonal to working precision, accumulating the rotations in v (b_in v = b_out).
// The column norms are then the singular values and v's columns the right singular vectors.
void jacobi(std::vector<double>& b, std::size_t m, Mat3& v) {
    v = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    constexpr double eps = 1e-17;
    for (int sweep = 0; sweep < 64; ++sweep) {
        bool rotated = false;
        for (std::size_t p = 0; p < 2; ++p) {
            for (std::size_t q = p + 1; q < 3; ++q) {
                double alpha = 0, beta = 0, gamma = 0;
                for (std::size_t i = 0; i < m; ++i) {
                    const double bp = b[3 * i + p], bq = b[3 * i + q];
                    alpha += bp * bp;
                    beta += bq * bq;
                    gamma += bp * bq;
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
                for (std::size_t i = 0; i < m; ++i) {
                    const double bp = b[3 * i + p], bq = b[3 * i + q];
                    b[3 * i + p] = c * bp - s * bq;
                    b[3 * i + q] = s * bp + c * bq;
                }
                for (std::size_t i = 0; i < 3; ++i) {
                    const double vp = v[i][p], vq = v[i][q];
                    v[i][p] = c * vp - s * vq;
                    v[i][q] = s * vp + c * vq;
                }
            }
        }
        if (!rotated) {
            break;
        }
    }
}

Vec3 columnNorms(const std::vector<double>& b, std::size_t m) {
    Vec3 out{};
    for (std::size_t k = 0; k < 3; ++k) {
        double sum = 0;
        for (std::size_t i = 0; i < m; ++i) {
            sum += b[3 * i + k] * b[3 * i + k];
        }
        out[k] = std::sqrt(sum);
    }
    return out;
}

// column indices by descending norm (stable: equal norms keep column order)
std::array<std::size_t, 3> descending(const Vec3& norms) {
    std::array<std::size_t, 3> order{0, 1, 2};
    std::stable_sort(order.begin(), order.end(),
                     [&norms](std::size_t x, std::size_t y) { return norms[x] > norms[y]; });
    return order;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

Vec3 unit(const Vec3& a) {
    const double n = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    return {a[0] / n, a[1] / n, a[2] / n};
}

SvdRows3 sortedRight(const Vec3& norms, const Mat3& v, std::array<std::size_t, 3>& order) {
    order = descending(norms);
    SvdRows3 out;
    for (std::size_t k = 0; k < 3; ++k) {
        out.s[k] = norms[order[k]];
        for (std::size_t i = 0; i < 3; ++i) {
            out.vt[k][i] = v[i][order[k]];
        }
    }
    return out;
}

}  // namespace

Vec3 jacobiSingularValues(const Mat3& a) {
    std::vector<double> b(9);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            b[3 * i + j] = a[i][j];
        }
    }
    Mat3 v{};
    jacobi(b, 3, v);
    return columnNorms(b, 3);
}

Svd3 svd3x3(const Mat3& a) {
    std::vector<double> b(9);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            b[3 * i + j] = a[i][j];
        }
    }
    Mat3 v{};
    jacobi(b, 3, v);
    std::array<std::size_t, 3> order{};
    const SvdRows3 right = sortedRight(columnNorms(b, 3), v, order);
    Svd3 out;
    out.s = right.s;
    out.vt = right.vt;
    // left singular vectors: the rotated columns divided by their norms; columns of (numerically)
    // zero singular values complete an orthonormal basis
    const double tiny = out.s[0] * 1e-13;
    std::array<Vec3, 3> u{};
    std::size_t rank = 0;
    for (std::size_t k = 0; k < 3; ++k) {
        const double sk = out.s[k];
        if (!(sk > tiny) || sk == 0.0) {
            break;
        }
        const std::size_t c = order[k];
        u[k] = {b[c] / sk, b[3 + c] / sk, b[6 + c] / sk};
        rank = k + 1;
    }
    if (rank == 0) {
        u[0] = {1, 0, 0};
        rank = 1;
    }
    if (rank == 1) {
        // the coordinate axis least aligned with u0, made orthogonal to it
        std::size_t axis = 0;
        for (std::size_t i = 1; i < 3; ++i) {
            if (std::abs(u[0][i]) < std::abs(u[0][axis])) {
                axis = i;
            }
        }
        Vec3 e{0, 0, 0};
        e[axis] = 1;
        u[1] = unit(cross(cross(u[0], e), u[0]));
        rank = 2;
    }
    if (rank == 2) {
        u[2] = cross(u[0], u[1]);
    }
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            out.u[i][k] = u[k][i];
        }
    }
    return out;
}

SvdRows3 svdRows3(std::span<const Vec3> rows) {
    const std::size_t m = rows.size();
    std::vector<double> b(3 * m);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            b[3 * i + j] = rows[i][j];
        }
    }
    Mat3 v{};
    jacobi(b, m, v);
    std::array<std::size_t, 3> order{};
    return sortedRight(columnNorms(b, m), v, order);
}

std::vector<SvdRows3> svdRows3Batch(std::span<const Vec3> rows, std::size_t k) {
    if (k == 0 || rows.size() % k != 0) {
        throw std::invalid_argument("svdRows3Batch: the rows do not form k x 3 matrices");
    }
    std::vector<SvdRows3> out;
    out.reserve(rows.size() / k);
    for (std::size_t first = 0; first < rows.size(); first += k) {
        out.push_back(svdRows3(rows.subspan(first, k)));
    }
    return out;
}

double det3(const Mat3& r) {
    return r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) -
           r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) +
           r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
}

}  // namespace mdbridge::core
