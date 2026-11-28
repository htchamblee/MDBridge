// Powell / Brent (src/numerics/powell.cpp) and least_squares' TRF with LSMR (least_squares.cpp,
// lsmr.cpp) against SciPy 1.18.1 (data/numerics/numerics_cases.json "powell", "brent",
// "least_squares", "groups", "permutations", "lsmr", "roots"). The objectives and residuals are the
// generator's, written with the same operations in the same order.

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/numerics/least_squares.hpp"
#include "mdbridge/numerics/powell.hpp"
#include "numerics_cases.hpp"

using namespace mdbridge;
using namespace mdbridge::test::numerics_cases;

namespace {

constexpr double kIterateTol = 1e-9;

// gauss6 and gauss6_shift go through the C library's exp: the recorded iterates are glibc's (the
// generator's math.exp). Another libm rounds some exp() an ulp away, and Brent's line search then
// stops at other points of the flat valley near the minimum: MSVC's moved the iterates by up to
// 1.6e-8 and nfev by 8 of 311 and 14 of 599, with the same nit (CI run 45882807857). Without glibc
// those two are compared within kOtherLibmTol (about 100 sqrt(eps), the scale of a line-search
// minimum moved by an ulp of f) and nfev within 5 %; rosen and quad3 use no libm and stay exact.
#if defined(__GLIBC__)
constexpr bool kGlibcLibm = true;
#else
constexpr bool kGlibcLibm = false;
#endif
constexpr double kOtherLibmTol = 1e-6;

bool usesLibm(const std::string& objective) {
    return objective == "gauss6" || objective == "gauss6_shift";
}

double sq(double v) {
    return v * v;  // Python's v ** 2 (pow(v, 2) is exact rounding of v * v)
}

// ---- the Powell objectives (make_numerics_cases.py) ----

double rosen(std::span<const double> x) {
    double s = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        s += 100.0 * sq(x[i + 1] - x[i] * x[i]) + sq(1.0 - x[i]);
    }
    return s;
}

const double kQuadA[3][3] = {{4.0, 1.0, 0.5}, {1.0, 3.0, -0.25}, {0.5, -0.25, 2.0}};
const double kQuadC[3] = {0.5, -1.25, 2.0};

double quad3(std::span<const double> x) {
    double s = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            s += (x[i] - kQuadC[i]) * kQuadA[i][j] * (x[j] - kQuadC[j]);
        }
    }
    return s + 0.125;
}

const double kGaussC[6] = {0.4, -0.3, 0.2, 1.1, -0.7, 0.5};
const double kGaussW[6] = {0.9, 1.4, 0.6, 0.3, 0.5, 0.8};

double gauss6(std::span<const double> x) {
    double s = 0.0;
    for (int i = 0; i < 6; ++i) {
        s += kGaussW[i] * sq(x[i] - kGaussC[i]);
    }
    s += 0.1 * x[0] * x[3];
    return -std::exp(-s);
}

double gauss6Shift(std::span<const double> x) {
    if (std::sqrt(x[3] * x[3] + x[4] * x[4] + x[5] * x[5]) > 1.2) {
        return 1.0;
    }
    return gauss6(x);
}

numerics::ScalarObjective powellObjective(const std::string& name) {
    if (name == "rosen") return rosen;
    if (name == "quad3") return quad3;
    if (name == "gauss6") return gauss6;
    if (name == "gauss6_shift") return gauss6Shift;
    throw std::runtime_error("unknown objective " + name);
}

numerics::ScalarFunction1d brentObjective(const std::string& name) {
    if (name == "parabola") return [](double t) { return sq(t - 2.0) + 1.0; };
    if (name == "quartic") return [](double t) { return std::pow(t, 4.0) - 3.0 * std::pow(t, 3.0) + 2.0; };
    if (name == "double_well") return [](double t) { return sq(t - 0.3) * sq(t + 1.0) + 0.1 * t; };
    if (name == "decreasing") return [](double t) { return -t; };
    if (name == "flat_far") return [](double t) { return sq(t - 250.0); };
    throw std::runtime_error("unknown objective " + name);
}

// ---- the least-squares residuals ----

std::vector<double> rRosen(std::span<const double> x) {
    return {10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]};
}

std::vector<double> rExpfit(std::span<const double> x) {
    std::vector<double> out;
    for (int i = 0; i < 12; ++i) {
        const double t = i / 4.0;
        const double y = 2.5 * std::exp(-1.3 * t) + 0.3 + 0.01 * (((i * 37) % 11) - 5) / 5.0;
        out.push_back(x[0] * std::exp(-x[1] * t) + x[2] - y);
    }
    return out;
}

std::vector<double> rUnder(std::span<const double> x) {
    return {x[0] + 2.0 * x[1] - 3.0 * x[2] - 1.0, x[0] * x[1] - 0.5};
}

std::array<double, 3> rodrigues(const double* v, const std::array<double, 3>& p) {
    const double th = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (th == 0.0) {
        return p;
    }
    const double k[3] = {v[0] / th, v[1] / th, v[2] / th};
    const double c = std::cos(th), s = std::sin(th);
    const double kxp[3] = {k[1] * p[2] - k[2] * p[1], k[2] * p[0] - k[0] * p[2], k[0] * p[1] - k[1] * p[0]};
    const double kp = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
    std::array<double, 3> out{};
    for (int i = 0; i < 3; ++i) {
        out[i] = p[i] * c + kxp[i] * s + k[i] * kp * (1.0 - c);
    }
    return out;
}

const std::array<std::array<double, 3>, 5> kRotP = {
    {{1.0, 0.0, 0.0}, {0.0, 1.5, 0.0}, {0.0, 0.0, 2.0}, {1.0, 1.0, -0.5}, {-0.5, 0.25, 1.0}}};

std::vector<double> rRigid(std::span<const double> x) {
    static const std::vector<std::array<double, 3>> q = [] {
        const double truth[3] = {0.3, -0.2, 0.5};
        const double shift[3] = {0.4, -1.0, 0.25};
        std::vector<std::array<double, 3>> out;
        for (const auto& p : kRotP) {
            std::array<double, 3> r = rodrigues(truth, p);
            for (int i = 0; i < 3; ++i) {
                r[i] = r[i] + shift[i];
            }
            out.push_back(r);
        }
        return out;
    }();
    std::vector<double> res;
    for (std::size_t k = 0; k < kRotP.size(); ++k) {
        const std::array<double, 3> r = rodrigues(x.data(), kRotP[k]);
        for (int i = 0; i < 3; ++i) {
            res.push_back(r[i] + x[3 + i] - q[k][i]);
        }
    }
    for (int i = 0; i < 6; ++i) {
        res.push_back(0.05 * x[i]);
    }
    return res;
}

constexpr std::size_t kChainN = 10;

std::vector<double> rChain(std::span<const double> x) {
    std::vector<double> res;
    for (std::size_t i = 0; i + 1 < kChainN; ++i) {
        res.push_back(x[i] * x[i] + x[i + 1] - (0.5 + 0.1 * static_cast<double>(i)));
    }
    for (std::size_t i = 0; i + 2 < kChainN; ++i) {
        res.push_back(0.5 * (x[i] - x[i + 1] * x[i + 2]));
    }
    for (std::size_t i = 0; i < kChainN; ++i) {
        res.push_back(0.3 * (x[i] - 1.0));
    }
    return res;
}

numerics::ResidualFunction residuals(const std::string& name) {
    if (name == "rosen") return rRosen;
    if (name == "expfit") return rExpfit;
    if (name == "under") return rUnder;
    if (name == "rigid") return rRigid;
    if (name == "chain") return rChain;
    throw std::runtime_error("unknown residuals " + name);
}

numerics::CsrMatrix patternOf(const core::Value& c) {
    const std::vector<std::size_t> shape = counts(at(c, "shape"));
    return numerics::CsrMatrix::pattern(shape[0], shape[1], counts(at(c, "rows")), counts(at(c, "cols")));
}

}  // namespace

TEST(Powell, IteratesMatchScipy) {
    for (const core::Value& c : cases("powell")) {
        const std::string name = str(at(c, "objective"));
        const core::Value& o = at(c, "options");
        numerics::PowellOptions options;
        if (has(o, "xtol")) options.xtol = num(at(o, "xtol"));
        if (has(o, "ftol")) options.ftol = num(at(o, "ftol"));
        if (has(o, "maxfev")) options.maxfev = count(at(o, "maxfev"));
        if (has(o, "maxiter")) options.maxiter = count(at(o, "maxiter"));
        if (has(o, "direc")) {
            for (const core::Value& row : list(at(o, "direc"))) {
                options.direc.push_back(nums(row));
            }
        }
        std::vector<std::vector<double>> xs;
        std::vector<double> funs;
        const std::vector<double> x0 = nums(at(c, "x0"));
        const numerics::PowellResult r = numerics::minimizePowell(
            powellObjective(name), x0, options, [&](std::span<const double> x, double f) {
                xs.emplace_back(x.begin(), x.end());
                funs.push_back(f);
                return false;
            });
        const core::Array& iterates = list(at(c, "iterates"));
        const std::string label = name + " from " + std::to_string(x0[0]);
        const bool exact = kGlibcLibm || !usesLibm(name);
        const double tol = exact ? kIterateTol : kOtherLibmTol;
        ASSERT_EQ(xs.size(), iterates.size()) << label;
        for (std::size_t k = 0; k < xs.size(); ++k) {
            expectClose(xs[k], nums(at(iterates[k], "x")), tol, label + " iterate " + std::to_string(k));
            EXPECT_TRUE(close(funs[k], num(at(iterates[k], "fun")), tol)) << label << " " << k;
        }
        expectClose(r.x, nums(at(c, "x")), tol, label + " x");
        EXPECT_TRUE(close(r.fun, num(at(c, "fun")), tol)) << label;
        const core::Array& direc = list(at(c, "direc"));
        ASSERT_EQ(r.direc.size(), direc.size());
        for (std::size_t k = 0; k < direc.size(); ++k) {
            expectClose(r.direc[k], nums(direc[k]), tol, label + " direc");
        }
        EXPECT_EQ(r.nit, count(at(c, "nit"))) << label;
        if (exact) {
            EXPECT_EQ(r.nfev, count(at(c, "nfev"))) << label;
        } else {
            const std::size_t want = count(at(c, "nfev"));
            EXPECT_LE((r.nfev > want ? r.nfev - want : want - r.nfev) * 20, want) << label << ": nfev " << r.nfev;
        }
        EXPECT_EQ(r.status, static_cast<int>(count(at(c, "status")))) << label;
        EXPECT_EQ(r.message, str(at(c, "message"))) << label;
    }
}

TEST(Powell, CallbackStopsAndInputIsChecked) {
    std::size_t calls = 0;
    const std::vector<double> x0 = {-1.2, 1.0};
    const auto r = numerics::minimizePowell(rosen, x0, {}, [&](std::span<const double>, double) {
        return ++calls == 2;
    });
    EXPECT_EQ(r.nit, 2u);
    EXPECT_THROW(numerics::minimizePowell(rosen, std::vector<double>{}), std::invalid_argument);
    numerics::PowellOptions bad;
    bad.direc = {{1.0, 0.0}};
    EXPECT_THROW(numerics::minimizePowell(rosen, x0, bad), std::invalid_argument);
}

TEST(Brent, MatchesScipy) {
    for (const core::Value& c : cases("brent")) {
        const std::string name = str(at(c, "objective"));
        const numerics::BrentResult r = numerics::minimizeScalarBrent(brentObjective(name), num(at(c, "xtol")));
        EXPECT_TRUE(close(r.x, num(at(c, "x")), 1e-12)) << name << " " << r.x;
        EXPECT_TRUE(close(r.fun, num(at(c, "fun")), 1e-12)) << name;
        EXPECT_EQ(r.nit, count(at(c, "nit"))) << name;
        EXPECT_EQ(r.nfev, count(at(c, "nfev"))) << name;
        EXPECT_EQ(r.success, *at(c, "success").get<bool>()) << name;
    }
}

TEST(LeastSquares, IteratesMatchScipy) {
    const numerics::CsrMatrix chain = patternOf(at(all(), "chain_sparsity"));
    for (const core::Value& c : cases("least_squares")) {
        const std::string name = str(at(c, "residuals"));
        const core::Value& o = at(c, "options");
        numerics::LeastSquaresOptions options;
        std::string label = name;
        if (has(o, "x_scale")) {
            if (at(o, "x_scale").isString()) {
                options.xScaleJac = true;
            } else {
                options.xScale = nums(at(o, "x_scale"));
            }
            label += " x_scale";
        }
        if (has(o, "ftol")) options.ftol = num(at(o, "ftol"));
        if (has(o, "diff_step")) options.diffStep = num(at(o, "diff_step"));
        if (has(o, "max_nfev")) options.maxNfev = count(at(o, "max_nfev"));
        if (has(o, "tr_solver")) {
            options.trSolver = str(at(o, "tr_solver")) == "lsmr" ? numerics::TrustRegionSolver::Lsmr
                                                                  : numerics::TrustRegionSolver::Exact;
            label += " " + str(at(o, "tr_solver"));
        }
        if (has(o, "jac_sparsity")) {
            options.jacSparsity = chain;
            label += " sparse";
        }
        if (has(o, "tr_options")) {
            const core::Value& t = at(o, "tr_options");
            if (has(t, "atol")) options.lsmr.atol = num(at(t, "atol"));
            if (has(t, "btol")) options.lsmr.btol = num(at(t, "btol"));
            if (has(t, "regularize")) options.lsmr.regularize = *at(t, "regularize").get<bool>();
        }
        // SciPy's iterates move by ulp_spread when its residuals move by an ulp (the finite
        // differences amplify rounding by 1 / sqrt(eps)): the port, whose SVD, QR and sums round
        // differently from LAPACK and BLAS, is held to that, and to 1e-9 at least
        const double spread = num(at(c, "ulp_spread"));
        const bool stable = *at(c, "ulp_same_count").get<bool>();
        const double tol = std::max(kIterateTol, 4.0 * spread);
        for (const std::size_t threads : {std::size_t{1}, std::size_t{3}}) {
            options.jacobianThreads = threads;
            std::vector<std::vector<double>> xs;
            std::vector<double> costs;
            std::vector<std::size_t> nfevs;
            const numerics::LeastSquaresResult r = numerics::leastSquares(
                residuals(name), nums(at(c, "x0")), options,
                [&](std::span<const double> x, double cost, std::size_t, std::size_t nfev) {
                    xs.emplace_back(x.begin(), x.end());
                    costs.push_back(cost);
                    nfevs.push_back(nfev);
                    return false;
                });
            if (!stable) {
                // an ulp changes SciPy's own path (LSMR with loose tolerances): the same minimum
                EXPECT_TRUE(close(r.cost, num(at(c, "cost")), 1e-6)) << label << ": " << r.cost;
                EXPECT_TRUE(r.success) << label;
                continue;
            }
            const core::Array& iterates = list(at(c, "iterates"));
            ASSERT_EQ(xs.size(), iterates.size()) << label;
            for (std::size_t k = 0; k < xs.size(); ++k) {
                expectClose(xs[k], nums(at(iterates[k], "x")), tol, label + " iterate " + std::to_string(k));
                EXPECT_TRUE(close(costs[k], num(at(iterates[k], "cost")), tol)) << label << " " << k;
                EXPECT_EQ(nfevs[k], count(at(iterates[k], "nfev"))) << label << " " << k;
            }
            expectClose(r.x, nums(at(c, "x")), tol, label + " x");
            EXPECT_TRUE(close(r.cost, num(at(c, "cost")), tol)) << label;
            expectClose(r.fun, nums(at(c, "fun")), tol, label + " fun");
            expectClose(r.grad, nums(at(c, "grad")), std::max(1e-7, 4.0 * spread), label + " grad");
            EXPECT_TRUE(close(r.optimality, num(at(c, "optimality")), std::max(1e-7, 4.0 * spread))) << label;
            EXPECT_EQ(r.nfev, count(at(c, "nfev"))) << label;
            EXPECT_EQ(r.njev, count(at(c, "njev"))) << label;
            EXPECT_EQ(r.status, static_cast<int>(count(at(c, "status")))) << label;
            EXPECT_EQ(r.message, str(at(c, "message"))) << label;
        }
    }
}

TEST(LeastSquares, ChecksItsInput) {
    const std::vector<double> x0 = {1.0, 1.0, 0.0};
    numerics::LeastSquaresOptions o;
    o.xScale = {1.0, -1.0, 1.0};
    EXPECT_THROW(numerics::leastSquares(rExpfit, x0, o), std::invalid_argument);
    o.xScale = {1.0, 1.0};
    EXPECT_THROW(numerics::leastSquares(rExpfit, x0, o), std::invalid_argument);
    o = {};
    o.ftol = o.xtol = o.gtol = 0.0;
    EXPECT_THROW(numerics::leastSquares(rExpfit, x0, o), std::invalid_argument);
    o = {};
    o.maxNfev = 0;
    EXPECT_THROW(numerics::leastSquares(rExpfit, x0, o), std::invalid_argument);
    o = {};
    const auto nanAtStart = [](std::span<const double>) { return std::vector<double>{std::nan("")}; };
    EXPECT_THROW(numerics::leastSquares(nanAtStart, x0, o), std::invalid_argument);
    o.jacSparsity = numerics::CsrMatrix::pattern(12, 3, std::vector<std::size_t>{0}, std::vector<std::size_t>{0});
    o.trSolver = numerics::TrustRegionSolver::Exact;
    EXPECT_THROW(numerics::leastSquares(rExpfit, x0, o), std::invalid_argument);
    o.trSolver = numerics::TrustRegionSolver::Auto;
    o.jacSparsity = numerics::CsrMatrix::pattern(5, 3, std::vector<std::size_t>{0}, std::vector<std::size_t>{0});
    EXPECT_THROW(numerics::leastSquares(rExpfit, x0, o), std::invalid_argument);
}

TEST(LeastSquares, CallbackStops) {
    const std::vector<double> x0 = {1.0, 1.0, 0.0};
    const auto r = numerics::leastSquares(rExpfit, x0, {}, [](std::span<const double>, double, std::size_t nit,
                                                              std::size_t) { return nit == 2; });
    EXPECT_EQ(r.status, -2);
    EXPECT_EQ(r.nit, 2u);
    EXPECT_FALSE(r.success);
}

// MDBridge's extension (docs/DIFFERENCES.md NU-03): with denseFromSparsity the structure groups the
// finite differences and the dense Jacobian goes to the exact subproblem, whose SVD of a tall matrix
// leaves out the zero rows and reduces the rest by Householder QR (40 variables: the bidiagonal SVD).
// A banded linear least-squares problem with one-sided (clipped) rows reaches the solution of its
// normal equations, solved here by Cholesky (to the differences' rounding), and the same point as the
// dense differences.
// 40 variables: the Householder QR of the tall matrix; 240: the row-wise Givens QR (200 or more, NU-03)
void denseFromSparsityCase(std::size_t n) {
    constexpr std::size_t bands = 3;
    std::vector<std::size_t> rowOf, colOf;
    std::vector<std::vector<std::pair<std::size_t, double>>> rows;
    for (std::size_t i = 0; i < n; ++i) {  // band rows: x_i .. x_{i+2}
        std::vector<std::pair<std::size_t, double>> row;
        for (std::size_t k = 0; k < bands && i + k < n; ++k) {
            row.emplace_back(i + k, 1.0 + 0.3 * std::sin(static_cast<double>(7 * i + 3 * k)));
        }
        rows.push_back(row);
    }
    for (std::size_t i = 0; i < n; ++i) {  // plain rows
        rows.push_back({{i, 0.5 + 0.1 * std::cos(static_cast<double>(i))}});
    }
    std::vector<double> b(rows.size());
    for (std::size_t r = 0; r < rows.size(); ++r) {
        b[r] = std::sin(static_cast<double>(r) * 0.7) * 2.0;
        for (const auto& [c, v] : rows[r]) {
            rowOf.push_back(r);
            colOf.push_back(c);
            (void)v;
        }
    }
    const std::size_t linear = rows.size();
    for (std::size_t i = 0; i < n; ++i) {  // one-sided rows, inactive near the solution (zero rows)
        rowOf.push_back(linear + i);
        colOf.push_back(i);
    }
    const auto fun = [&](std::span<const double> x) {
        std::vector<double> f(linear + n);
        for (std::size_t r = 0; r < linear; ++r) {
            double acc = -b[r];
            for (const auto& [c, v] : rows[r]) acc += v * x[c];
            f[r] = acc;
        }
        for (std::size_t i = 0; i < n; ++i) f[linear + i] = std::max(0.0, x[i] - 50.0);
        return f;
    };
    // the normal equations A^T A x = A^T b of the linear rows, by Cholesky
    std::vector<double> ata(n * n, 0.0), atb(n, 0.0);
    for (std::size_t r = 0; r < linear; ++r) {
        for (const auto& [c, v] : rows[r]) {
            atb[c] += v * b[r];
            for (const auto& [d, w] : rows[r]) ata[c * n + d] += v * w;
        }
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t k = 0; k < j; ++k) ata[j * n + j] -= ata[j * n + k] * ata[j * n + k];
        ata[j * n + j] = std::sqrt(ata[j * n + j]);
        for (std::size_t i = j + 1; i < n; ++i) {
            for (std::size_t k = 0; k < j; ++k) ata[i * n + j] -= ata[i * n + k] * ata[j * n + k];
            ata[i * n + j] /= ata[j * n + j];
        }
    }
    std::vector<double> y(n), want(n);
    for (std::size_t i = 0; i < n; ++i) {
        double acc = atb[i];
        for (std::size_t k = 0; k < i; ++k) acc -= ata[i * n + k] * y[k];
        y[i] = acc / ata[i * n + i];
    }
    for (std::size_t i = n; i-- > 0;) {
        double acc = y[i];
        for (std::size_t k = i + 1; k < n; ++k) acc -= ata[k * n + i] * want[k];
        want[i] = acc / ata[i * n + i];
    }
    const std::vector<double> x0(n, 1.0);
    numerics::LeastSquaresOptions grouped;
    grouped.ftol = grouped.xtol = grouped.gtol = 1e-12;
    grouped.jacSparsity = numerics::CsrMatrix::pattern(linear + n, n, rowOf, colOf);
    grouped.denseFromSparsity = true;
    grouped.jacobianThreads = 3;
    const auto a = numerics::leastSquares(fun, x0, grouped);
    numerics::LeastSquaresOptions dense;
    dense.ftol = dense.xtol = dense.gtol = 1e-12;
    const auto d = numerics::leastSquares(fun, x0, dense);
    EXPECT_TRUE(a.success) << a.message;
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(a.x[i], want[i], 1e-6) << i;  // the differences' rounding: about 3e-8
        EXPECT_NEAR(d.x[i], a.x[i], 1e-12) << i;
    }
    grouped.trSolver = numerics::TrustRegionSolver::Lsmr;
    EXPECT_THROW(numerics::leastSquares(fun, x0, grouped), std::invalid_argument);
}

TEST(LeastSquares, DenseFromSparsityAndTheTallExactSubproblem) {
    denseFromSparsityCase(40);
    denseFromSparsityCase(240);
}

TEST(LeastSquares, ColumnGroupsAndPermutations) {
    for (const core::Value& c : cases("permutations")) {
        EXPECT_EQ(numerics::randomStatePermutation(0, count(at(c, "n"))), counts(at(c, "permutation")))
            << count(at(c, "n"));
    }
    for (const core::Value& c : cases("groups")) {
        EXPECT_EQ(numerics::groupColumns(patternOf(c)), counts(at(c, "groups")));
    }
}

TEST(LeastSquares, LsmrMatchesScipy) {
    for (const core::Value& c : cases("lsmr")) {
        std::vector<std::vector<double>> a;
        for (const core::Value& row : list(at(c, "a"))) {
            a.push_back(nums(row));
        }
        const std::size_t m = a.size(), n = a[0].size();
        const auto matvec = [&](std::span<const double> x) {
            std::vector<double> y(m, 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    y[i] += a[i][j] * x[j];
                }
            }
            return y;
        };
        const auto rmatvec = [&](std::span<const double> y) {
            std::vector<double> x(n, 0.0);
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t i = 0; i < m; ++i) {
                    x[j] += a[i][j] * y[i];
                }
            }
            return x;
        };
        std::optional<std::size_t> maxiter;
        if (has(c, "maxiter")) {
            maxiter = count(at(c, "maxiter"));
        }
        const numerics::LsmrResult r = numerics::lsmr(m, n, matvec, rmatvec, nums(at(c, "b")), num(at(c, "damp")),
                                                      num(at(c, "atol")), num(at(c, "btol")), 1e8, maxiter);
        // the Krylov iterations amplify rounding (BLAS sums differently): within 4x SciPy's own
        // movement when b moves by an ulp, and 1e-12 at least
        const double tol = std::max(1e-12, 4.0 * num(at(c, "ulp_spread")));
        expectClose(r.x, nums(at(c, "x")), tol, "lsmr x");
        EXPECT_EQ(r.istop, static_cast<int>(count(at(c, "istop"))));
        EXPECT_EQ(r.itn, count(at(c, "itn")));
        for (const char* key : {"normr", "normar", "normA", "condA", "normx"}) {
            const double want = num(at(c, key));
            const double got = std::string(key) == "normr"    ? r.normr
                               : std::string(key) == "normar" ? r.normar
                               : std::string(key) == "normA"  ? r.normA
                               : std::string(key) == "condA"  ? r.condA
                                                              : r.normx;
            // the residual norms of a consistent system are rounding noise near 0: 1e-8 absolute
            const bool residual = std::string(key) == "normr" || std::string(key) == "normar";
            EXPECT_TRUE(close(got, want, residual ? 1e-8 : std::max(1e-10, 4.0 * num(at(c, "ulp_spread")))))
                << key << ": " << got << " vs " << want;
        }
    }
}

TEST(LeastSquares, RealRootsOfTheSubspaceQuartic) {
    for (const core::Value& c : cases("roots")) {
        const std::vector<double> coeffs = nums(at(c, "c"));
        const std::vector<double> got = numerics::realPolynomialRoots(coeffs);
        const std::vector<double> want = nums(at(c, "real"));
        if (*at(c, "multiple").get<bool>()) {
            continue;  // LAPACK splits a multiple root by about sqrt(eps); see the generator
        }
        expectClose(got, want, 1e-12, "roots");
    }
}
