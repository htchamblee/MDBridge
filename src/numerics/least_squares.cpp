// The trust region reflective method of scipy.optimize.least_squares (no bounds, linear loss,
// 2-point finite-difference Jacobian, dense or sparse), ported from SciPy 1.18.1:
// scipy/optimize/_lsq/least_squares.py (least_squares, check_x_scale, check_tolerance,
// check_jac_sparsity), _lsq/trf.py (trf_no_bounds), _lsq/common.py (solve_lsq_trust_region,
// solve_trust_region_2d, update_tr_radius, build_quadratic_1d, minimize_quadratic_1d,
// evaluate_quadratic, compute_grad, compute_jac_scale, right_multiplied_operator,
// check_termination), scipy/optimize/_numdiff.py (approx_derivative's 2-point scheme,
// _compute_absolute_step, _dense_difference, _sparse_difference, group_columns) and
// scipy/optimize/_group_columns.pyx (group_sparse), used under SciPy's licence:
//
// Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
// All rights reserved.
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
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to
//    endorse or promote products derived from this software without specific prior written
//    permission.
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
// numpy.random.RandomState(0).permutation (the column order of group_columns) follows NumPy 2.5.5
// (numpy/random/mtrand.pyx shuffle, distributions.c random_interval, the MT19937 generator of
// Matsumoto and Nishimura as numpy/random/src/mt19937 implements it; BSD-3-Clause, see
// THIRD_PARTY_NOTICES.md: Copyright (c) 2005-2025, NumPy Developers; MT19937: Copyright (C) 1997 -
// 2002, Makoto Matsumoto and Takuji Nishimura).
//
// SciPy's LAPACK calls are replaced by MDBridge's own linear algebra (least_squares.hpp): the thin
// SVD (one-sided Jacobi), the 2-column QR (Householder, LAPACK's sign convention) and the real roots
// of np.roots' quartic. Built without floating-point contraction.

#include "mdbridge/numerics/least_squares.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>

#include "mdbridge/core/detail/numpy_sum.hpp"
#include "numerics_internal.hpp"

namespace mdbridge::numerics {

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

using detail::dot;
using detail::norm2;

double normInf(std::span<const double> a) {  // norm(g, ord=np.inf): max |g| (NaN propagates)
    double m = 0.0;
    for (double v : a) {
        const double x = std::fabs(v);
        if (std::isnan(x)) {
            return x;
        }
        m = std::max(m, x);
    }
    return m;
}

double npSum(std::span<const double> a) {  // np.sum: numpy's pairwise sum
    return core::detail::numpySum(a);
}

// ---- the Jacobian ---------------------------------------------------------------------------------

struct Jacobian {
    std::size_t m = 0, n = 0;
    bool sparse = false;
    std::vector<double> dense;  // m x n, row-major
    CsrMatrix csr;

    // J.dot(x)
    std::vector<double> matvec(std::span<const double> x) const {
        std::vector<double> y(m, 0.0);
        if (sparse) {
            for (std::size_t i = 0; i < m; ++i) {
                double s = 0.0;
                for (std::size_t k = csr.indptr[i]; k < csr.indptr[i + 1]; ++k) {
                    s += csr.data[k] * x[csr.indices[k]];
                }
                y[i] = s;
            }
        } else {
            for (std::size_t i = 0; i < m; ++i) {
                double s = 0.0;
                for (std::size_t j = 0; j < n; ++j) {
                    s += dense[i * n + j] * x[j];
                }
                y[i] = s;
            }
        }
        return y;
    }

    // J.T.dot(y)
    std::vector<double> rmatvec(std::span<const double> y) const {
        std::vector<double> g(n, 0.0);
        if (sparse) {
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t k = csr.indptr[i]; k < csr.indptr[i + 1]; ++k) {
                    g[csr.indices[k]] += csr.data[k] * y[i];
                }
            }
        } else {
            for (std::size_t j = 0; j < n; ++j) {
                double s = 0.0;
                for (std::size_t i = 0; i < m; ++i) {
                    s += dense[i * n + j] * y[i];
                }
                g[j] = s;
            }
        }
        return g;
    }

    // np.sum(J**2, axis=0) ** 0.5 (dense: rows added in order), J.power(2).sum(axis=0) (sparse)
    std::vector<double> columnNorms() const {
        std::vector<double> s(n, 0.0);
        if (sparse) {
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t k = csr.indptr[i]; k < csr.indptr[i + 1]; ++k) {
                    s[csr.indices[k]] += csr.data[k] * csr.data[k];
                }
            }
        } else {
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    s[j] += dense[i * n + j] * dense[i * n + j];
                }
            }
        }
        for (double& v : s) {
            v = std::pow(v, 0.5);
        }
        return s;
    }
};

// _compute_absolute_step(rel_step, x0, f0, "2-point")
std::vector<double> absoluteStep(std::span<const double> x0, const std::optional<double>& relStep) {
    const double rstep = std::pow(kEps, 0.5);
    std::vector<double> h(x0.size());
    for (std::size_t i = 0; i < x0.size(); ++i) {
        const double sign = (x0[i] >= 0 ? 1.0 : 0.0) * 2 - 1;
        const double defaultStep = rstep * sign * std::max(1.0, std::fabs(x0[i]));
        if (!relStep) {
            h[i] = defaultStep;
        } else {
            const double step = *relStep * sign * std::fabs(x0[i]);
            const double dx = (x0[i] + step) - x0[i];
            h[i] = dx == 0 ? defaultStep : step;
        }
    }
    return h;
}

class FiniteDifferences {
public:
    FiniteDifferences(const ResidualFunction& fun, std::size_t m, std::size_t n, const LeastSquaresOptions& options)
        : m_fun(fun), m_m(m), m_n(n), m_relStep(options.diffStep), m_threads(options.jacobianThreads) {
        if (options.jacSparsity) {
            m_sparse = true;
            m_densify = options.denseFromSparsity;
            m_structure = *options.jacSparsity;
            m_groups = groupColumns(m_structure);
            m_groupCount = 0;
            for (std::size_t g : m_groups) {
                m_groupCount = std::max(m_groupCount, g + 1);
            }
            // the structure by column (csc), rows ascending within a column
            m_colStart.assign(n + 1, 0);
            for (std::size_t k = 0; k < m_structure.indices.size(); ++k) {
                ++m_colStart[m_structure.indices[k] + 1];
            }
            for (std::size_t j = 0; j < n; ++j) {
                m_colStart[j + 1] += m_colStart[j];
            }
            m_colRows.resize(m_structure.indices.size());
            std::vector<std::size_t> fill(m_colStart.begin(), m_colStart.end() - 1);
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t k = m_structure.indptr[i]; k < m_structure.indptr[i + 1]; ++k) {
                    m_colRows[fill[m_structure.indices[k]]++] = i;
                }
            }
            // find(structure[:, cols]): the nonzero entries only (stored zeros group, but get no value)
            m_colNonzero.assign(m_colRows.size(), true);
            if (!m_structure.data.empty()) {
                std::vector<std::size_t> at(m_colStart.begin(), m_colStart.end() - 1);
                for (std::size_t i = 0; i < m; ++i) {
                    for (std::size_t k = m_structure.indptr[i]; k < m_structure.indptr[i + 1]; ++k) {
                        m_colNonzero[at[m_structure.indices[k]]++] = m_structure.data[k] != 0;
                    }
                }
            }
        }
    }

    bool sparse() const { return m_sparse && !m_densify; }

    Jacobian operator()(std::span<const double> x0, std::span<const double> f0) const {
        const std::vector<double> h = absoluteStep(x0, m_relStep);
        Jacobian jac;
        jac.m = m_m;
        jac.n = m_n;
        jac.sparse = m_sparse && !m_densify;
        const std::size_t evaluations = m_sparse ? m_groupCount : m_n;
        std::vector<std::vector<double>> points(evaluations, std::vector<double>(x0.begin(), x0.end()));
        for (std::size_t e = 0; e < evaluations; ++e) {
            if (m_sparse) {  // x0 + h * (groups == e)
                for (std::size_t k = 0; k < m_n; ++k) {
                    points[e][k] = x0[k] + (m_groups[k] == e ? h[k] : 0.0);
                }
            } else {
                points[e][e] = x0[e] + h[e];
            }
        }
        std::vector<std::vector<double>> values(evaluations);
        detail::parallelFor(evaluations, m_threads, [&](std::size_t e) {
            values[e] = m_fun(std::span<const double>(points[e]));
            if (values[e].size() != m_m) {
                throw std::invalid_argument("the residual function returned a different number of values");
            }
        });
        if (!m_sparse) {
            jac.dense.assign(m_m * m_n, 0.0);
            for (std::size_t j = 0; j < m_n; ++j) {
                const double dx = (x0[j] + h[j]) - x0[j];
                for (std::size_t i = 0; i < m_m; ++i) {
                    jac.dense[i * m_n + j] = (values[j][i] - f0[i]) / dx;
                }
            }
            return jac;
        }
        if (m_densify) {  // _sparse_difference(...).toarray()
            jac.dense.assign(m_m * m_n, 0.0);
            for (std::size_t j = 0; j < m_n; ++j) {
                const std::size_t g = m_groups[j];
                const double dx = points[g][j] - x0[j];
                for (std::size_t k = m_colStart[j]; k < m_colStart[j + 1]; ++k) {
                    if (m_colNonzero[k]) {
                        const std::size_t i = m_colRows[k];
                        jac.dense[i * m_n + j] = (values[g][i] - f0[i]) / dx;
                    }
                }
            }
            return jac;
        }
        // _sparse_difference: entry (i, j) of each group's columns, then csr_array (sorted)
        std::vector<std::vector<std::pair<std::size_t, double>>> rows(m_m);
        for (std::size_t g = 0; g < m_groupCount; ++g) {
            for (std::size_t j = 0; j < m_n; ++j) {
                if (m_groups[j] != g) {
                    continue;
                }
                const double dx = points[g][j] - x0[j];
                for (std::size_t k = m_colStart[j]; k < m_colStart[j + 1]; ++k) {
                    if (!m_colNonzero[k]) {
                        continue;
                    }
                    const std::size_t i = m_colRows[k];
                    const double df = values[g][i] - f0[i];
                    rows[i].emplace_back(j, df / dx);
                }
            }
        }
        CsrMatrix& c = jac.csr;
        c.rows = m_m;
        c.cols = m_n;
        c.indptr.assign(m_m + 1, 0);
        for (std::size_t i = 0; i < m_m; ++i) {
            std::sort(rows[i].begin(), rows[i].end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            for (const auto& [j, v] : rows[i]) {
                c.indices.push_back(j);
                c.data.push_back(v);
            }
            c.indptr[i + 1] = c.indices.size();
        }
        return jac;
    }

private:
    const ResidualFunction& m_fun;
    std::size_t m_m, m_n;
    std::optional<double> m_relStep;
    std::size_t m_threads;
    bool m_sparse = false;
    bool m_densify = false;
    CsrMatrix m_structure;
    std::vector<std::size_t> m_groups;
    std::size_t m_groupCount = 0;
    std::vector<std::size_t> m_colStart, m_colRows;
    std::vector<bool> m_colNonzero;
};

// ---- dense linear algebra ---------------------------------------------------------------------

// svd(a, full_matrices=False) of an m x n matrix (row-major): u (m x k), s (k, descending), v (n x k,
// the right singular vectors as columns, i.e. Vt.T), k = min(m, n). One-sided Jacobi (Hestenes)
// rotations on the columns of a (m >= n) or of a.T (m < n). Columns of u for a zero singular value
// are zero (the callers multiply them by s).
struct ThinSvd {
    std::size_t k = 0;
    std::vector<double> u, s, v;
};

ThinSvd thinSvd(const std::vector<double>& a, std::size_t m, std::size_t n) {
    const bool transposed = m < n;
    const std::size_t rows = transposed ? n : m;  // the matrix whose columns are rotated: rows x cols
    const std::size_t cols = transposed ? m : n;
    std::vector<double> w(rows * cols);  // column-major: column j at w[j * rows]
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (transposed) {
                w[i * rows + j] = a[i * n + j];  // a.T: column i is row i of a
            } else {
                w[j * rows + i] = a[i * n + j];
            }
        }
    }
    std::vector<double> rot(cols * cols, 0.0);  // column-major, starts as the identity
    for (std::size_t j = 0; j < cols; ++j) {
        rot[j * cols + j] = 1.0;
    }
    for (int sweep = 0; sweep < 80; ++sweep) {
        bool rotated = false;
        for (std::size_t p = 0; p + 1 < cols; ++p) {
            for (std::size_t q = p + 1; q < cols; ++q) {
                double* wp = &w[p * rows];
                double* wq = &w[q * rows];
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (std::size_t i = 0; i < rows; ++i) {
                    alpha += wp[i] * wp[i];
                    beta += wq[i] * wq[i];
                    gamma += wp[i] * wq[i];
                }
                if (gamma == 0.0 || std::fabs(gamma) <= kEps * std::sqrt(alpha * beta)) {
                    continue;
                }
                rotated = true;
                const double zeta = (beta - alpha) / (2.0 * gamma);
                const double t = (zeta >= 0 ? 1.0 : -1.0) / (std::fabs(zeta) + std::sqrt(1.0 + zeta * zeta));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = c * t;
                for (std::size_t i = 0; i < rows; ++i) {
                    const double x = wp[i], y = wq[i];
                    wp[i] = c * x - s * y;
                    wq[i] = s * x + c * y;
                }
                double* rp = &rot[p * cols];
                double* rq = &rot[q * cols];
                for (std::size_t i = 0; i < cols; ++i) {
                    const double x = rp[i], y = rq[i];
                    rp[i] = c * x - s * y;
                    rq[i] = s * x + c * y;
                }
            }
        }
        if (!rotated) {
            break;
        }
    }
    std::vector<double> sv(cols);
    for (std::size_t j = 0; j < cols; ++j) {
        double ss = 0.0;
        for (std::size_t i = 0; i < rows; ++i) {
            ss += w[j * rows + i] * w[j * rows + i];
        }
        sv[j] = std::sqrt(ss);
    }
    std::vector<std::size_t> order(cols);
    for (std::size_t j = 0; j < cols; ++j) {
        order[j] = j;
    }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return sv[x] > sv[y]; });
    ThinSvd out;
    out.k = cols;
    out.s.resize(cols);
    // left vectors of the rotated matrix: w_j / s_j (rows x cols); right vectors: rot (cols x cols)
    std::vector<double> left(rows * cols, 0.0), right(cols * cols, 0.0);  // row-major
    for (std::size_t r = 0; r < cols; ++r) {
        const std::size_t j = order[r];
        out.s[r] = sv[j];
        for (std::size_t i = 0; i < rows; ++i) {
            left[i * cols + r] = sv[j] > 0 ? w[j * rows + i] / sv[j] : 0.0;
        }
        for (std::size_t i = 0; i < cols; ++i) {
            right[i * cols + r] = rot[j * cols + i];
        }
    }
    if (transposed) {  // a.T = L S R.T  ->  a = R S L.T: u = R (m x m), v = L (n x m)
        out.u = std::move(right);
        out.v = std::move(left);
    } else {
        out.u = std::move(left);
        out.v = std::move(right);
    }
    return out;
}

// qr(s, mode="economic") of an n x c matrix (c = 1 or 2 columns, column-major), Q only: Householder
// reflections as LAPACK's dgeqr2 / dorg2r make them (beta = -sign(alpha) * norm)
std::vector<std::vector<double>> qrEconomic(const std::vector<std::vector<double>>& columns, std::size_t n) {
    const std::size_t c = std::min(columns.size(), n);
    std::vector<std::vector<double>> a = columns;
    std::vector<std::vector<double>> vs(c);
    std::vector<double> taus(c, 0.0);
    for (std::size_t j = 0; j < c; ++j) {
        // dlarfg on a[j][j:]
        const double alpha = a[j][j];
        double xnorm = 0.0;
        {
            double ss = 0.0;
            for (std::size_t i = j + 1; i < n; ++i) {
                ss += a[j][i] * a[j][i];
            }
            xnorm = std::sqrt(ss);
        }
        std::vector<double> v(n, 0.0);
        v[j] = 1.0;
        double tau = 0.0;
        if (xnorm != 0.0) {
            const double beta = -std::copysign(std::hypot(alpha, xnorm), alpha);
            tau = (beta - alpha) / beta;
            const double scale = 1.0 / (alpha - beta);
            for (std::size_t i = j + 1; i < n; ++i) {
                v[i] = a[j][i] * scale;
            }
            a[j][j] = beta;
            for (std::size_t i = j + 1; i < n; ++i) {
                a[j][i] = 0.0;
            }
        }
        // apply H = I - tau v v^T to the remaining columns
        for (std::size_t q = j + 1; q < columns.size(); ++q) {
            double w = 0.0;
            for (std::size_t i = j; i < n; ++i) {
                w += v[i] * a[q][i];
            }
            for (std::size_t i = j; i < n; ++i) {
                a[q][i] -= tau * w * v[i];
            }
        }
        vs[j] = std::move(v);
        taus[j] = tau;
    }
    // Q = H_0 H_1 ... applied to the first c columns of the identity (dorg2r: from the last reflector)
    std::vector<std::vector<double>> q(c, std::vector<double>(n, 0.0));
    for (std::size_t j = 0; j < c; ++j) {
        q[j][j] = 1.0;
    }
    for (std::size_t r = c; r-- > 0;) {
        for (std::size_t col = 0; col < c; ++col) {
            double w = 0.0;
            for (std::size_t i = r; i < n; ++i) {
                w += vs[r][i] * q[col][i];
            }
            for (std::size_t i = r; i < n; ++i) {
                q[col][i] -= taus[r] * w * vs[r][i];
            }
        }
    }
    return q;
}

// The SVD of a square matrix a (n x n, row-major) as the exact subproblem needs it: singular values
// (descending), right singular vectors and U.T y, without forming U. Householder bidiagonalisation
// (y and V take the reflections), then implicit-shift QR sweeps on the bidiagonal (Golub and Kahan's
// step with Wilkinson's shift; Golub and Van Loan, Matrix Computations, 8.6). nullopt when a diagonal
// entry of the bidiagonal vanishes or the sweeps do not converge: the caller then uses the Jacobi SVD.
std::optional<ThinSvd> bidiagonalSvd(std::vector<double> a, std::size_t n, std::vector<double>& y) {
    std::vector<double> v(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        v[i * n + i] = 1.0;
    }
    std::vector<double> h(n);
    for (std::size_t k = 0; k < n; ++k) {
        // left reflection: zeros below a[k][k]; applied to a's columns k.. and to y
        double norm = 0.0;
        for (std::size_t i = k; i < n; ++i) {
            norm += a[i * n + k] * a[i * n + k];
        }
        norm = std::sqrt(norm);
        if (norm != 0.0) {
            const double alpha = a[k * n + k] >= 0 ? -norm : norm;
            for (std::size_t i = k; i < n; ++i) {
                h[i] = a[i * n + k];
            }
            h[k] -= alpha;
            double hh = 0.0;
            for (std::size_t i = k; i < n; ++i) {
                hh += h[i] * h[i];
            }
            for (std::size_t j = k + 1; j < n; ++j) {
                double t = 0.0;
                for (std::size_t i = k; i < n; ++i) {
                    t += h[i] * a[i * n + j];
                }
                t = 2.0 * t / hh;
                for (std::size_t i = k; i < n; ++i) {
                    a[i * n + j] -= t * h[i];
                }
            }
            double t = 0.0;
            for (std::size_t i = k; i < n; ++i) {
                t += h[i] * y[i];
            }
            t = 2.0 * t / hh;
            for (std::size_t i = k; i < n; ++i) {
                y[i] -= t * h[i];
            }
            a[k * n + k] = alpha;
            for (std::size_t i = k + 1; i < n; ++i) {
                a[i * n + k] = 0.0;
            }
        }
        if (k + 2 > n) {
            continue;
        }
        // right reflection: zeros right of a[k][k + 1]; applied to a's rows k + 1.. and to V
        norm = 0.0;
        for (std::size_t j = k + 1; j < n; ++j) {
            norm += a[k * n + j] * a[k * n + j];
        }
        norm = std::sqrt(norm);
        if (norm == 0.0) {
            continue;
        }
        const double alpha = a[k * n + k + 1] >= 0 ? -norm : norm;
        for (std::size_t j = k + 1; j < n; ++j) {
            h[j] = a[k * n + j];
        }
        h[k + 1] -= alpha;
        double hh = 0.0;
        for (std::size_t j = k + 1; j < n; ++j) {
            hh += h[j] * h[j];
        }
        const auto reflectRow = [&](double* row) {
            double t = 0.0;
            for (std::size_t j = k + 1; j < n; ++j) {
                t += row[j] * h[j];
            }
            t = 2.0 * t / hh;
            for (std::size_t j = k + 1; j < n; ++j) {
                row[j] -= t * h[j];
            }
        };
        for (std::size_t i = k + 1; i < n; ++i) {
            reflectRow(&a[i * n]);
        }
        for (std::size_t i = 0; i < n; ++i) {
            reflectRow(&v[i * n]);
        }
        a[k * n + k + 1] = alpha;
        for (std::size_t j = k + 2; j < n; ++j) {
            a[k * n + j] = 0.0;
        }
    }
    // a is upper bidiagonal: QR sweeps on it (dense storage, every rotation touches a few entries)
    const auto B = [&](std::size_t i, std::size_t j) -> double& { return a[i * n + j]; };
    const auto givens = [](double f, double g, double& c, double& sn) {
        const double r = std::hypot(f, g);
        if (r == 0.0) {
            c = 1.0;
            sn = 0.0;
        } else {
            c = f / r;
            sn = g / r;
        }
    };
    std::size_t hi = n - 1;
    std::size_t sweeps = 0;
    while (hi > 0) {
        for (std::size_t i = 0; i < hi; ++i) {
            if (std::fabs(B(i, i + 1)) <= kEps * (std::fabs(B(i, i)) + std::fabs(B(i + 1, i + 1)))) {
                B(i, i + 1) = 0.0;
            }
        }
        if (B(hi - 1, hi) == 0.0) {
            --hi;
            continue;
        }
        std::size_t lo = hi - 1;
        while (lo > 0 && B(lo - 1, lo) != 0.0) {
            --lo;
        }
        for (std::size_t i = lo; i <= hi; ++i) {
            if (B(i, i) == 0.0) {
                return std::nullopt;
            }
        }
        if (++sweeps > 100 * n) {
            return std::nullopt;
        }
        const double dm = B(hi - 1, hi - 1), dn = B(hi, hi), em = B(hi - 1, hi);
        const double el = hi - 1 > lo ? B(hi - 2, hi - 1) : 0.0;
        const double t11 = dm * dm + el * el, t12 = dm * em, t22 = dn * dn + em * em;
        const double delta = (t11 - t22) / 2.0;
        const double denom = delta + (delta >= 0 ? 1.0 : -1.0) * std::hypot(delta, t12);
        const double mu = denom == 0.0 ? t22 : t22 - t12 * t12 / denom;
        double f = B(lo, lo) * B(lo, lo) - mu;
        double g = B(lo, lo) * B(lo, lo + 1);
        for (std::size_t k = lo; k < hi; ++k) {
            double c = 1.0, sn = 0.0;
            givens(f, g, c, sn);  // columns k, k + 1: [f g] -> [r 0]
            for (std::size_t r = (k > lo ? k - 1 : lo); r <= k + 1; ++r) {
                const double x = B(r, k), z = B(r, k + 1);
                B(r, k) = c * x + sn * z;
                B(r, k + 1) = -sn * x + c * z;
            }
            for (std::size_t r = 0; r < n; ++r) {
                const double x = v[r * n + k], z = v[r * n + k + 1];
                v[r * n + k] = c * x + sn * z;
                v[r * n + k + 1] = -sn * x + c * z;
            }
            if (k > lo) {
                B(k - 1, k + 1) = 0.0;
            }
            f = B(k, k);
            g = B(k + 1, k);
            givens(f, g, c, sn);  // rows k, k + 1: the bulge below the diagonal
            for (std::size_t col = k; col <= std::min(hi, k + 2); ++col) {
                const double x = B(k, col), z = B(k + 1, col);
                B(k, col) = c * x + sn * z;
                B(k + 1, col) = -sn * x + c * z;
            }
            B(k + 1, k) = 0.0;
            const double yk = y[k], yk1 = y[k + 1];
            y[k] = c * yk + sn * yk1;
            y[k + 1] = -sn * yk + c * yk1;
            if (k + 1 < hi) {
                f = B(k, k + 1);
                g = B(k, k + 2);
            }
        }
    }
    std::vector<double> sv(n);
    for (std::size_t i = 0; i < n; ++i) {
        sv[i] = B(i, i);
        if (sv[i] < 0) {
            sv[i] = -sv[i];
            for (std::size_t r = 0; r < n; ++r) {
                v[r * n + i] = -v[r * n + i];
            }
        }
    }
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t z) { return sv[x] > sv[z]; });
    ThinSvd out;
    out.k = n;
    out.s.resize(n);
    out.v.resize(n * n);
    std::vector<double> uf(n);
    for (std::size_t r = 0; r < n; ++r) {
        out.s[r] = sv[order[r]];
        uf[r] = y[order[r]];
        for (std::size_t i = 0; i < n; ++i) {
            out.v[i * n + r] = v[i * n + order[r]];
        }
    }
    y = std::move(uf);
    return out;
}

// svd(J_h, full_matrices=False) and uf = U.T f for the exact subproblem. Rows of J_h that are all
// zero (inactive one-sided terms) add nothing to U.T f (their row of U is zero), so they are left out
// when at least n rows remain; a matrix that is then taller than wide is first reduced by
// Householder reflections (J_h = Q R, R n x n upper triangular, Q.T applied to f), and U = Q U_R
// gives uf = U_R.T (Q.T f)[:n]. The same numbers up to rounding as the SVD of the whole matrix,
// much faster for the tall Jacobians of many contact terms.
ThinSvd exactSvd(const std::vector<double>& a, std::size_t m, std::size_t n, const std::vector<double>& f,
                 std::vector<double>& uf) {
    std::vector<std::size_t> live;
    live.reserve(m);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (a[i * n + j] != 0.0) {
                live.push_back(i);
                break;
            }
        }
    }
    if (live.size() < n) {
        live.resize(m);
        for (std::size_t i = 0; i < m; ++i) {
            live[i] = i;
        }
    }
    const std::size_t rows = live.size();
    if (rows < 2 * n) {  // not much taller than wide: the SVD of the rows themselves
        std::vector<double> b(rows * n);
        std::vector<double> fb(rows);
        for (std::size_t r = 0; r < rows; ++r) {
            std::copy(a.begin() + static_cast<std::ptrdiff_t>(live[r] * n),
                      a.begin() + static_cast<std::ptrdiff_t>((live[r] + 1) * n), b.begin() + static_cast<std::ptrdiff_t>(r * n));
            fb[r] = f[live[r]];
        }
        ThinSvd svd = thinSvd(b, rows, n);
        uf.assign(svd.k, 0.0);
        for (std::size_t r = 0; r < svd.k; ++r) {
            double acc = 0.0;
            for (std::size_t i = 0; i < rows; ++i) {
                acc += svd.u[i * svd.k + r] * fb[i];
            }
            uf[r] = acc;
        }
        return svd;
    }
    // Householder QR, column-major copy of the live rows; qf = Q.T f
    std::vector<double> w(rows * n);
    std::vector<double> qf(rows);
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t j = 0; j < n; ++j) {
            w[j * rows + r] = a[live[r] * n + j];
        }
        qf[r] = f[live[r]];
    }
    for (std::size_t j = 0; j < n; ++j) {
        double* col = &w[j * rows];
        double norm = 0.0;
        for (std::size_t i = j; i < rows; ++i) {
            norm += col[i] * col[i];
        }
        norm = std::sqrt(norm);
        if (norm == 0.0) {
            continue;
        }
        const double alpha = col[j] >= 0 ? -norm : norm;
        // v = x - alpha e1 (stored in col[j..]), H = I - 2 v v.T / (v.T v)
        col[j] -= alpha;
        double vv = 0.0;
        for (std::size_t i = j; i < rows; ++i) {
            vv += col[i] * col[i];
        }
        const auto reflect = [&](double* y) {
            double s = 0.0;
            for (std::size_t i = j; i < rows; ++i) {
                s += col[i] * y[i];
            }
            s = 2.0 * s / vv;
            for (std::size_t i = j; i < rows; ++i) {
                y[i] -= s * col[i];
            }
        };
        for (std::size_t c = j + 1; c < n; ++c) {
            reflect(&w[c * rows]);
        }
        reflect(qf.data());
        col[j] = alpha;
        for (std::size_t i = j + 1; i < rows; ++i) {
            col[i] = 0.0;
        }
    }
    std::vector<double> r(n * n, 0.0);  // row-major upper triangle
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            r[i * n + j] = w[j * rows + i];
        }
    }
    if (n >= 16) {  // large: the bidiagonal QR (much faster than Jacobi rotations)
        std::vector<double> y(qf.begin(), qf.begin() + static_cast<std::ptrdiff_t>(n));
        if (std::optional<ThinSvd> fast = bidiagonalSvd(r, n, y)) {
            uf = std::move(y);
            return std::move(*fast);
        }
    }
    ThinSvd svd = thinSvd(r, n, n);
    uf.assign(svd.k, 0.0);
    for (std::size_t k = 0; k < svd.k; ++k) {
        double acc = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            acc += svd.u[i * svd.k + k] * qf[i];
        }
        uf[k] = acc;
    }
    return svd;
}

// ---- the trust-region subproblems -------------------------------------------------------------

struct LsqStep {
    std::vector<double> p;
    double alpha;
};

// solve_lsq_trust_region(n, m, uf, s, V, Delta, initial_alpha, rtol=0.01, max_iter=10)
LsqStep solveLsqTrustRegion(std::size_t n, std::size_t m, const std::vector<double>& uf, const ThinSvd& svd,
                            double Delta, double initialAlpha) {
    const std::size_t k = svd.k;
    const std::vector<double>& s = svd.s;
    std::vector<double> suf(k);
    for (std::size_t i = 0; i < k; ++i) {
        suf[i] = s[i] * uf[i];
    }
    auto vdot = [&](const std::vector<double>& c) {  // V.dot(c): n values
        std::vector<double> p(n, 0.0);
        for (std::size_t r = 0; r < n; ++r) {
            double acc = 0.0;
            for (std::size_t i = 0; i < k; ++i) {
                acc += svd.v[r * k + i] * c[i];
            }
            p[r] = acc;
        }
        return p;
    };
    auto phiAndDerivative = [&](double alpha) {
        std::vector<double> ratio(k), terms(k);
        for (std::size_t i = 0; i < k; ++i) {
            const double denom = s[i] * s[i] + alpha;
            ratio[i] = suf[i] / denom;
            terms[i] = suf[i] * suf[i] / std::pow(denom, 3.0);
        }
        const double pNorm = norm2(ratio);
        const double phi = pNorm - Delta;
        const double phiPrime = -npSum(terms) / pNorm;
        return std::pair<double, double>(phi, phiPrime);
    };
    bool fullRank = false;
    if (m >= n) {
        const double threshold = kEps * static_cast<double>(m) * s[0];
        fullRank = s[k - 1] > threshold;
    }
    if (fullRank) {
        std::vector<double> c(k);
        for (std::size_t i = 0; i < k; ++i) {
            c[i] = uf[i] / s[i];
        }
        std::vector<double> p = vdot(c);
        for (double& v : p) {
            v = -v;
        }
        if (norm2(p) <= Delta) {
            return {std::move(p), 0.0};
        }
    }
    double alphaUpper = norm2(suf) / Delta;
    double alphaLower = 0.0;
    if (fullRank) {
        const auto [phi, phiPrime] = phiAndDerivative(0.0);
        alphaLower = -phi / phiPrime;
    }
    double alpha = 0.0;
    if (!fullRank && initialAlpha == 0) {
        alpha = std::max(0.001 * alphaUpper, std::pow(alphaLower * alphaUpper, 0.5));
    } else {
        alpha = initialAlpha;
    }
    constexpr double rtol = 0.01;
    for (int it = 0; it < 10; ++it) {
        if (alpha < alphaLower || alpha > alphaUpper) {
            alpha = std::max(0.001 * alphaUpper, std::pow(alphaLower * alphaUpper, 0.5));
        }
        const auto [phi, phiPrime] = phiAndDerivative(alpha);
        if (phi < 0) {
            alphaUpper = alpha;
        }
        const double ratio = phi / phiPrime;
        alphaLower = std::max(alphaLower, alpha - ratio);
        alpha -= (phi + Delta) * ratio / Delta;
        if (std::fabs(phi) < rtol * Delta) {
            break;
        }
    }
    std::vector<double> c(k);
    for (std::size_t i = 0; i < k; ++i) {
        c[i] = suf[i] / (s[i] * s[i] + alpha);
    }
    std::vector<double> p = vdot(c);
    for (double& v : p) {
        v = -v;
    }
    const double scale = Delta / norm2(p);
    for (double& v : p) {
        v *= scale;
    }
    return {std::move(p), alpha};
}

// solve_trust_region_2d(B, g, Delta) for a 2 x 2 (or 1 x 1) B
std::vector<double> solveTrustRegion2d(const std::array<double, 4>& B, const std::array<double, 2>& g, std::size_t dim,
                                       double Delta) {
    for (std::size_t i = 0; i < dim * dim; ++i) {
        const double v = dim == 1 ? B[0] : B[i];
        if (!std::isfinite(v)) {
            throw std::invalid_argument("array must not contain infs or NaNs");
        }
    }
    // cho_factor (upper) and cho_solve; a non-positive pivot is LinAlgError (skipped)
    if (dim == 1) {
        if (B[0] > 0) {
            const double r = std::sqrt(B[0]);
            const double p = -((g[0] / r) / r);
            if (p * p <= Delta * Delta) {
                return {p};
            }
        }
        throw std::runtime_error("solveTrustRegion2d: a one-dimensional subspace whose step leaves the region");
    }
    const double b00 = B[0], b01 = B[1], b11 = B[3];
    if (b00 > 0) {
        const double r11 = std::sqrt(b00);
        const double r12 = b01 / r11;
        const double a22 = b11 - r12 * r12;
        if (a22 > 0) {
            const double r22 = std::sqrt(a22);
            const double y0 = g[0] / r11;
            const double y1 = (g[1] - r12 * y0) / r22;
            const double p1 = y1 / r22;
            const double p0 = (y0 - r12 * p1) / r11;
            const double q0 = -p0, q1 = -p1;
            if (q0 * q0 + q1 * q1 <= Delta * Delta) {
                return {q0, q1};
            }
        }
    }
    const double Delta2 = Delta * Delta;
    const double a = B[0] * Delta2;
    const double b = B[1] * Delta2;
    const double c = B[3] * Delta2;
    const double d = g[0] * Delta;
    const double f = g[1] * Delta;
    const std::array<double, 5> coeffs = {-b + d, 2 * (a - c + f), 6 * b, 2 * (-a + c + f), -b - d};
    const std::vector<double> ts = realPolynomialRoots(coeffs);
    if (ts.empty()) {
        throw std::runtime_error("solveTrustRegion2d: the quartic has no real root");
    }
    double best = 0.0;
    std::array<double, 2> bestP{};
    for (std::size_t i = 0; i < ts.size(); ++i) {
        const double t = ts[i];
        const double p0 = Delta * (2 * t / (1 + t * t));
        const double p1 = Delta * ((1 - t * t) / (1 + t * t));
        const double bp0 = B[0] * p0 + B[1] * p1;
        const double bp1 = B[2] * p0 + B[3] * p1;
        const double value = 0.5 * (p0 * bp0 + p1 * bp1) + (g[0] * p0 + g[1] * p1);
        if (i == 0 || value < best) {
            best = value;
            bestP = {p0, p1};
        }
    }
    return {bestP[0], bestP[1]};
}

// ---- MT19937 (numpy.random.RandomState) ---------------------------------------------------------

class Mt19937 {
public:
    explicit Mt19937(std::uint32_t seed) {
        std::uint32_t s = seed;
        for (std::size_t pos = 0; pos < kN; ++pos) {
            m_key[pos] = s;
            s = static_cast<std::uint32_t>(1812433253UL * (s ^ (s >> 30)) + pos + 1);
        }
        m_pos = kN;
    }

    std::uint32_t next32() {
        if (m_pos >= kN) {
            generate();
        }
        std::uint32_t y = m_key[m_pos++];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);
        return y;
    }

    std::uint64_t next64() {
        const std::uint64_t hi = next32();
        const std::uint64_t lo = next32();
        return (hi << 32) | lo;
    }

private:
    static constexpr std::size_t kN = 624;
    static constexpr std::size_t kM = 397;

    void generate() {
        constexpr std::uint32_t upper = 0x80000000UL, lower = 0x7fffffffUL, matrixA = 0x9908b0dfUL;
        std::size_t i = 0;
        for (; i < kN - kM; ++i) {
            const std::uint32_t y = (m_key[i] & upper) | (m_key[i + 1] & lower);
            m_key[i] = m_key[i + kM] ^ (y >> 1) ^ ((y & 1U) ? matrixA : 0U);
        }
        for (; i < kN - 1; ++i) {
            const std::uint32_t y = (m_key[i] & upper) | (m_key[i + 1] & lower);
            m_key[i] = m_key[i + kM - kN] ^ (y >> 1) ^ ((y & 1U) ? matrixA : 0U);
        }
        const std::uint32_t y = (m_key[kN - 1] & upper) | (m_key[0] & lower);
        m_key[kN - 1] = m_key[kM - 1] ^ (y >> 1) ^ ((y & 1U) ? matrixA : 0U);
        m_pos = 0;
    }

    std::array<std::uint32_t, kN> m_key{};
    std::size_t m_pos = 0;
};

// random_interval(bitgen, max): uniform in [0, max] by masked rejection
std::uint64_t randomInterval(Mt19937& rng, std::uint64_t max) {
    if (max == 0) {
        return 0;
    }
    std::uint64_t mask = max;
    mask |= mask >> 1;
    mask |= mask >> 2;
    mask |= mask >> 4;
    mask |= mask >> 8;
    mask |= mask >> 16;
    mask |= mask >> 32;
    std::uint64_t value = 0;
    if (max <= 0xffffffffULL) {
        while ((value = (rng.next32() & mask)) > max) {
        }
    } else {
        while ((value = (rng.next64() & mask)) > max) {
        }
    }
    return value;
}

// the real roots of c[0] t^d + ... + c[d] (no leading or trailing zero, c[0] != 0), ascending
std::vector<double> realRootsTrimmed(const std::vector<double>& c) {
    const std::size_t d = c.size() - 1;
    if (d == 0) {
        return {};
    }
    if (d == 1) {
        return {-c[1] / c[0]};
    }
    auto value = [&](double t) {  // Horner
        double v = c[0];
        for (std::size_t i = 1; i <= d; ++i) {
            v = v * t + c[i];
        }
        return v;
    };
    auto slope = [&](double t) {
        double v = static_cast<double>(d) * c[0];
        for (std::size_t i = 1; i < d; ++i) {
            v = v * t + static_cast<double>(d - i) * c[i];
        }
        return v;
    };
    // the stationary points split the line into monotone pieces
    std::vector<double> dc(d);
    for (std::size_t i = 0; i < d; ++i) {
        dc[i] = static_cast<double>(d - i) * c[i];
    }
    std::vector<double> critical = realRootsTrimmed(dc);
    // Cauchy's bound: every root lies within 1 + max |c[i] / c[0]|
    double bound = 0.0;
    for (std::size_t i = 1; i <= d; ++i) {
        bound = std::max(bound, std::fabs(c[i] / c[0]));
    }
    bound = 1.0 + bound;
    std::vector<double> edges;
    edges.push_back(-bound);
    for (double t : critical) {
        if (t > -bound && t < bound && (edges.empty() || t > edges.back())) {
            edges.push_back(t);
        }
    }
    edges.push_back(bound);
    std::vector<double> roots;
    for (std::size_t e = 0; e + 1 < edges.size(); ++e) {
        double lo = edges[e], hi = edges[e + 1];
        double flo = value(lo), fhi = value(hi);
        if (flo == 0.0) {
            if (roots.empty() || roots.back() != lo) {
                roots.push_back(lo);
            }
            continue;
        }
        if (fhi == 0.0) {
            if (e + 2 == edges.size()) {
                roots.push_back(hi);
            }
            continue;  // the next piece starts with it
        }
        if ((flo < 0) == (fhi < 0)) {
            continue;
        }
        // safeguarded Newton steps inside the bracket, down to adjacent doubles
        double t = 0.5 * (lo + hi);
        for (int it = 0; it < 200; ++it) {
            const double ft = value(t);
            if (ft == 0.0) {
                break;
            }
            if ((ft < 0) == (flo < 0)) {
                lo = t;
                flo = ft;
            } else {
                hi = t;
            }
            const double mid = 0.5 * (lo + hi);
            if (!(mid > lo && mid < hi)) {
                t = std::fabs(value(lo)) <= std::fabs(value(hi)) ? lo : hi;
                break;
            }
            const double dt = slope(t);
            double next = dt != 0.0 ? t - ft / dt : mid;
            if (!(next > lo && next < hi)) {
                next = mid;
            }
            t = next;
        }
        roots.push_back(t);
    }
    return roots;
}

}  // namespace

// ---- public pieces -------------------------------------------------------------------------------

CsrMatrix CsrMatrix::pattern(std::size_t rows, std::size_t cols, std::span<const std::size_t> rowOf,
                             std::span<const std::size_t> colOf) {
    if (rowOf.size() != colOf.size()) {
        throw std::invalid_argument("CsrMatrix::pattern: rows and columns differ in length");
    }
    std::vector<std::vector<std::size_t>> byRow(rows);
    for (std::size_t k = 0; k < rowOf.size(); ++k) {
        if (rowOf[k] >= rows || colOf[k] >= cols) {
            throw std::invalid_argument("CsrMatrix::pattern: an entry lies outside the shape");
        }
        byRow[rowOf[k]].push_back(colOf[k]);
    }
    CsrMatrix out;
    out.rows = rows;
    out.cols = cols;
    out.indptr.assign(rows + 1, 0);
    for (std::size_t i = 0; i < rows; ++i) {
        std::vector<std::size_t>& r = byRow[i];
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
        out.indices.insert(out.indices.end(), r.begin(), r.end());
        out.indptr[i + 1] = out.indices.size();
    }
    out.data.assign(out.indices.size(), 1.0);
    return out;
}

std::vector<std::size_t> randomStatePermutation(std::uint32_t seed, std::size_t n) {
    std::vector<std::size_t> a(n);
    for (std::size_t i = 0; i < n; ++i) {
        a[i] = i;
    }
    Mt19937 rng(seed);
    for (std::size_t i = n; i-- > 1;) {  // reversed(range(1, n))
        const auto j = static_cast<std::size_t>(randomInterval(rng, i));
        std::swap(a[i], a[j]);
    }
    return a;
}

std::vector<std::size_t> groupColumns(const CsrMatrix& structure) {
    const std::size_t m = structure.rows, n = structure.cols;
    // csc of A[:, order]: the rows of each permuted column
    const std::vector<std::size_t> order = randomStatePermutation(0, n);
    std::vector<std::vector<std::size_t>> colRows(n);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t k = structure.indptr[i]; k < structure.indptr[i + 1]; ++k) {
            colRows[structure.indices[k]].push_back(i);  // every stored entry, as csc_array(A) keeps them
        }
    }
    std::vector<std::vector<std::size_t>> permuted(n);
    for (std::size_t k = 0; k < n; ++k) {
        permuted[k] = colRows[order[k]];
    }
    // group_sparse
    std::vector<std::ptrdiff_t> groups(n, -1);
    std::ptrdiff_t current = 0;
    std::vector<int> unionRows(m, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (groups[i] >= 0) {
            continue;
        }
        groups[i] = current;
        bool allGrouped = true;
        std::fill(unionRows.begin(), unionRows.end(), 0);
        for (std::size_t r : permuted[i]) {
            unionRows[r] = 1;
        }
        for (std::size_t j = 0; j < n; ++j) {
            if (groups[j] < 0) {
                allGrouped = false;
            } else {
                continue;
            }
            bool intersect = false;
            for (std::size_t r : permuted[j]) {
                if (unionRows[r] == 1) {
                    intersect = true;
                    break;
                }
            }
            if (!intersect) {
                for (std::size_t r : permuted[j]) {
                    unionRows[r] = 1;
                }
                groups[j] = current;
            }
        }
        if (allGrouped) {
            break;
        }
        ++current;
    }
    std::vector<std::size_t> out(n);  // groups[order] = groups.copy()
    for (std::size_t k = 0; k < n; ++k) {
        out[order[k]] = static_cast<std::size_t>(groups[k]);
    }
    return out;
}

std::vector<double> realPolynomialRoots(std::span<const double> c) {
    std::size_t first = 0, last = c.size();
    while (first < last && c[first] == 0.0) {
        ++first;
    }
    if (first == last) {
        return {};
    }
    std::size_t trailing = 0;
    while (c[last - 1] == 0.0) {
        --last;
        ++trailing;
    }
    std::vector<double> trimmed(c.begin() + static_cast<std::ptrdiff_t>(first),
                                c.begin() + static_cast<std::ptrdiff_t>(last));
    if (trimmed.size() > 5) {
        throw std::invalid_argument("realPolynomialRoots: degree above 4");
    }
    for (double v : trimmed) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("realPolynomialRoots: non-finite coefficient");
        }
    }
    std::vector<double> roots = realRootsTrimmed(trimmed);
    roots.insert(roots.end(), trailing, 0.0);
    std::sort(roots.begin(), roots.end());
    return roots;
}

// ---- least_squares --------------------------------------------------------------------------------

LeastSquaresResult leastSquares(const ResidualFunction& fun, std::span<const double> x0span,
                                const LeastSquaresOptions& options, const LeastSquaresCallback& callback) {
    const std::size_t n = x0span.size();
    if (n == 0) {
        throw std::invalid_argument("`x0` must have at least one variable.");
    }
    if (options.maxNfev && *options.maxNfev == 0) {
        throw std::invalid_argument("`max_nfev` must be None or positive integer.");
    }
    std::vector<double> x0(x0span.begin(), x0span.end());
    // check_x_scale
    std::vector<double> xScale = options.xScale;
    if (options.xScaleJac) {
        if (!xScale.empty()) {
            throw std::invalid_argument("`x_scale` must be 'jac' or array_like with positive numbers.");
        }
    } else {
        if (xScale.empty()) {
            xScale.assign(n, 1.0);
        }
        for (double v : xScale) {
            if (!(std::isfinite(v) && v > 0)) {
                throw std::invalid_argument("`x_scale` must be 'jac' or array_like with positive numbers.");
            }
        }
        if (xScale.size() == 1 && n != 1) {
            xScale.assign(n, xScale[0]);
        }
        if (xScale.size() != n) {
            throw std::invalid_argument("Inconsistent shapes between `x_scale` and `x0`.");
        }
    }
    // check_tolerance
    const double ftol = options.ftol, xtol = options.xtol, gtol = options.gtol;
    if (ftol < kEps && xtol < kEps && gtol < kEps) {
        throw std::invalid_argument("At least one of the tolerances must be higher than machine epsilon (2.22e-16).");
    }
    TrustRegionSolver trSolver = options.trSolver;
    if (options.jacSparsity && trSolver == TrustRegionSolver::Exact && !options.denseFromSparsity) {
        throw std::invalid_argument("tr_solver='exact' is incompatible with `jac_sparsity`.");
    }
    if (options.denseFromSparsity && trSolver == TrustRegionSolver::Lsmr) {
        throw std::invalid_argument("a Jacobian made dense from `jac_sparsity` is solved with tr_solver='exact'.");
    }

    std::vector<double> f = fun(std::span<const double>(x0));
    const std::size_t m = f.size();
    for (double v : f) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("Residuals are not finite in the initial point.");
        }
    }
    if (options.jacSparsity && (options.jacSparsity->rows != m || options.jacSparsity->cols != n ||
                                options.jacSparsity->indptr.size() != m + 1)) {
        throw std::invalid_argument("`jac_sparsity` has wrong shape.");
    }
    auto checkedFun = [&](std::span<const double> x) {
        std::vector<double> v = fun(x);
        if (v.size() != m) {
            throw std::invalid_argument("the residual function returned a different number of values");
        }
        return v;
    };
    const ResidualFunction wrapped = checkedFun;
    const FiniteDifferences jacobian(wrapped, m, n, options);
    if (trSolver == TrustRegionSolver::Auto) {
        trSolver = jacobian.sparse() ? TrustRegionSolver::Lsmr : TrustRegionSolver::Exact;
    }

    // trf_no_bounds
    std::vector<double> x = x0;
    std::size_t nfev = 1;
    Jacobian J = jacobian(x, f);
    std::size_t njev = 1;
    double cost = 0.5 * dot(f, f);
    std::vector<double> g = J.rmatvec(f);
    std::vector<double> scale(n), scaleInv(n);
    if (options.xScaleJac) {
        scaleInv = J.columnNorms();
        for (double& v : scaleInv) {
            if (v == 0) {
                v = 1;
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            scale[i] = 1 / scaleInv[i];
        }
    } else {
        scale = xScale;
        for (std::size_t i = 0; i < n; ++i) {
            scaleInv[i] = 1 / xScale[i];
        }
    }
    double Delta = 0.0;
    {
        std::vector<double> t(n);
        for (std::size_t i = 0; i < n; ++i) {
            t[i] = x0[i] * scaleInv[i];
        }
        Delta = norm2(t);
    }
    if (Delta == 0) {
        Delta = 1.0;
    }
    double regTerm = 0;
    const double damp = options.lsmr.damp;
    const bool regularize = options.lsmr.regularize;
    const std::size_t maxNfev = options.maxNfev ? *options.maxNfev : n * 100;
    double alpha = 0.0;  // the Levenberg-Marquardt parameter
    int terminationStatus = 0;
    bool terminated = false;
    std::size_t iteration = 0;
    double gNorm = 0.0;
    std::vector<double> xNew, fNew;
    double costNew = 0.0;

    while (true) {
        gNorm = normInf(g);
        if (gNorm < gtol) {
            terminationStatus = 1;
            terminated = true;
        }
        if (terminated || nfev == maxNfev) {
            break;
        }
        const std::vector<double>& d = scale;
        std::vector<double> gH(n);
        for (std::size_t i = 0; i < n; ++i) {
            gH[i] = d[i] * g[i];
        }
        // J_h = J * d: products with it
        auto jhDot = [&](std::span<const double> v) {
            std::vector<double> t(n);
            for (std::size_t i = 0; i < n; ++i) {
                t[i] = v[i] * d[i];
            }
            return J.matvec(t);
        };
        auto jhTDot = [&](std::span<const double> y) {
            std::vector<double> t = J.rmatvec(y);
            for (std::size_t i = 0; i < n; ++i) {
                t[i] = d[i] * t[i];
            }
            return t;
        };
        ThinSvd svd;
        std::vector<double> uf;
        std::vector<double> jhDense;
        std::vector<std::vector<double>> S;
        std::array<double, 4> BS{};
        std::array<double, 2> gS{};
        std::size_t dimS = 0;
        if (trSolver == TrustRegionSolver::Exact) {
            jhDense = J.dense;
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    jhDense[i * n + j] *= d[j];
                }
            }
            svd = exactSvd(jhDense, m, n, f, uf);
        } else {
            if (regularize) {
                std::vector<double> sDir(n);
                for (std::size_t i = 0; i < n; ++i) {
                    sDir[i] = -gH[i];
                }
                const std::vector<double> v = jhDot(sDir);
                double a = dot(v, v);
                a *= 0.5;
                const double b = dot(gH, sDir);
                const double toTr = Delta / norm2(gH);
                // minimize_quadratic_1d(a, b, 0, to_tr)[1]
                std::vector<double> ts = {0.0, toTr};
                if (a != 0) {
                    const double extremum = -0.5 * b / a;
                    if (0.0 < extremum && extremum < toTr) {
                        ts.push_back(extremum);
                    }
                }
                double agValue = 0.0;
                for (std::size_t i = 0; i < ts.size(); ++i) {
                    const double y = ts[i] * (a * ts[i] + b) + 0;
                    if (i == 0 || y < agValue) {
                        agValue = y;
                    }
                }
                regTerm = -agValue / (Delta * Delta);
            }
            const double dampFull = std::pow(damp * damp + regTerm, 0.5);
            const LsmrResult gn = lsmr(m, n, jhDot, jhTDot, f, dampFull, options.lsmr.atol, options.lsmr.btol,
                                       options.lsmr.conlim, options.lsmr.maxiter);
            S = qrEconomic({gH, gn.x}, n);
            dimS = S.size();
            std::vector<std::vector<double>> JS(dimS);
            for (std::size_t c = 0; c < dimS; ++c) {
                JS[c] = jhDot(S[c]);
            }
            for (std::size_t r = 0; r < dimS; ++r) {
                for (std::size_t c = 0; c < dimS; ++c) {
                    BS[r * 2 + c] = dot(JS[r], JS[c]);
                }
                gS[r] = dot(S[r], gH);
            }
            if (dimS == 1) {
                BS[1] = BS[2] = BS[3] = 0.0;
            }
        }

        double actualReduction = -1;
        double stepNorm = 0.0;
        while (actualReduction <= 0 && nfev < maxNfev) {
            std::vector<double> stepH;
            if (trSolver == TrustRegionSolver::Exact) {
                LsqStep st = solveLsqTrustRegion(n, m, uf, svd, Delta, alpha);
                stepH = std::move(st.p);
                alpha = st.alpha;
            } else {
                const std::vector<double> pS = solveTrustRegion2d(BS, gS, dimS, Delta);
                stepH.assign(n, 0.0);
                for (std::size_t i = 0; i < n; ++i) {
                    double acc = 0.0;
                    for (std::size_t c = 0; c < dimS; ++c) {
                        acc += S[c][i] * pS[c];
                    }
                    stepH[i] = acc;
                }
            }
            // evaluate_quadratic(J_h, g_h, step_h)
            std::vector<double> Js;
            if (trSolver == TrustRegionSolver::Exact) {
                Js.assign(m, 0.0);
                for (std::size_t i = 0; i < m; ++i) {
                    double acc = 0.0;
                    for (std::size_t j = 0; j < n; ++j) {
                        acc += jhDense[i * n + j] * stepH[j];
                    }
                    Js[i] = acc;
                }
            } else {
                Js = jhDot(stepH);
            }
            const double predictedReduction = -(0.5 * dot(Js, Js) + dot(stepH, gH));
            std::vector<double> step(n);
            xNew.assign(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                step[i] = d[i] * stepH[i];
                xNew[i] = x[i] + step[i];
            }
            fNew = wrapped(std::span<const double>(xNew));
            ++nfev;
            const double stepHNorm = norm2(stepH);
            bool finite = true;
            for (double v : fNew) {
                finite = finite && std::isfinite(v);
            }
            if (!finite) {
                Delta = 0.25 * stepHNorm;
                continue;
            }
            costNew = 0.5 * dot(fNew, fNew);
            actualReduction = cost - costNew;
            // update_tr_radius
            double ratio = 0.0;
            if (predictedReduction > 0) {
                ratio = actualReduction / predictedReduction;
            } else if (predictedReduction == actualReduction && actualReduction == 0) {
                ratio = 1;
            } else {
                ratio = 0;
            }
            double DeltaNew = Delta;
            if (ratio < 0.25) {
                DeltaNew = 0.25 * stepHNorm;
            } else if (ratio > 0.75 && stepHNorm > 0.95 * Delta) {
                DeltaNew *= 2.0;
            }
            stepNorm = norm2(step);
            // check_termination
            const bool ftolSatisfied = actualReduction < ftol * cost && ratio > 0.25;
            const bool xtolSatisfied = stepNorm < xtol * (xtol + norm2(x));
            if (ftolSatisfied && xtolSatisfied) {
                terminationStatus = 4;
            } else if (ftolSatisfied) {
                terminationStatus = 2;
            } else if (xtolSatisfied) {
                terminationStatus = 3;
            }
            if (ftolSatisfied || xtolSatisfied) {
                terminated = true;
                break;
            }
            alpha *= Delta / DeltaNew;
            Delta = DeltaNew;
        }

        if (actualReduction > 0) {
            x = xNew;
            f = fNew;
            cost = costNew;
            J = jacobian(x, f);
            ++njev;
            g = J.rmatvec(f);
            if (options.xScaleJac) {  // compute_jac_scale(J, scale_inv)
                const std::vector<double> norms = J.columnNorms();
                for (std::size_t i = 0; i < n; ++i) {
                    scaleInv[i] = std::max(norms[i], scaleInv[i]);
                    scale[i] = 1 / scaleInv[i];
                }
            }
        }
        ++iteration;
        if (callback && callback(std::span<const double>(x), cost, iteration, nfev)) {
            terminationStatus = -2;
            terminated = true;
            break;
        }
    }
    if (!terminated) {
        terminationStatus = 0;
    }

    static const char* const messages[] = {
        "Stopped because `callback` function raised `StopIteration` or returned `True`",
        "Improper input parameters status returned from `leastsq`",
        "The maximum number of function evaluations is exceeded.",
        "`gtol` termination condition is satisfied.",
        "`ftol` termination condition is satisfied.",
        "`xtol` termination condition is satisfied.",
        "Both `ftol` and `xtol` termination conditions are satisfied."};
    LeastSquaresResult out;
    out.x = std::move(x);
    out.cost = cost;
    out.fun = std::move(f);
    out.grad = std::move(g);
    out.optimality = gNorm;
    out.nfev = nfev;
    out.njev = njev;
    out.nit = iteration;
    out.status = terminationStatus;
    out.success = terminationStatus > 0;
    out.message = messages[terminationStatus + 2];
    return out;
}

}  // namespace mdbridge::numerics
