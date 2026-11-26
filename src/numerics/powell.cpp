// Powell's method and Brent's line search as scipy.optimize computes them, ported from SciPy 1.18.1
// (scipy/optimize/_optimize.py: _minimize_powell, _linesearch_powell, _minimize_scalar_brent, Brent,
// bracket, _recover_from_bracket_error, _wrap_scalar_function_maxfun_validation), used under SciPy's
// licence and the notice of its optimize module:
//
// ******NOTICE***************
// optimize.py module by Travis E. Oliphant
//
// You may copy and use this module as you see fit with no
// guarantee implied provided you keep this notice in all copies.
// *****END NOTICE************
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
// Only the unbounded paths are ported (MDATK passes no bounds). Built without floating-point
// contraction, so that every product and sum rounds as in Python.

#include "mdbridge/numerics/powell.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mdbridge::numerics {

namespace {

constexpr std::size_t kUnlimited = std::numeric_limits<std::size_t>::max();

// _MaxFuncCallError: the objective was called maxfun times already
struct MaxFunctionCalls {};

// BracketError: bracket() ended without a valid bracket; the points it had
struct BracketFailure {
    double xa, xb, xc, fa, fb, fc;
    std::size_t funcalls;
};

struct Bracket {
    double xa, xb, xc, fa, fb, fc;
    std::size_t funcalls;
};

// scipy.optimize.bracket(func, xa=0.0, xb=1.0, grow_limit=110.0, maxiter=1000)
template <class F>
Bracket bracket(F& func) {
    constexpr double gold = 1.618034;
    constexpr double verysmall = 1e-21;
    constexpr double growLimit = 110.0;
    constexpr std::size_t maxiter = 1000;
    double xa = 0.0, xb = 1.0;
    double fa = func(xa);
    double fb = func(xb);
    if (fa < fb) {
        std::swap(xa, xb);
        std::swap(fa, fb);
    }
    double xc = xb + gold * (xb - xa);
    double fc = func(xc);
    std::size_t funcalls = 3;
    std::size_t iter = 0;
    while (fc < fb) {
        const double tmp1 = (xb - xa) * (fb - fc);
        const double tmp2 = (xb - xc) * (fb - fa);
        const double val = tmp2 - tmp1;
        double denom = 0.0;
        if (std::fabs(val) < verysmall) {
            denom = 2.0 * verysmall;
        } else {
            denom = 2.0 * val;
        }
        double w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
        const double wlim = xb + growLimit * (xc - xb);
        if (iter > maxiter) {
            throw std::runtime_error(
                "No valid bracket was found before the iteration limit was reached. Consider trying different "
                "initial points or increasing `maxiter`.");
        }
        ++iter;
        double fw = 0.0;
        if ((w - xc) * (xb - w) > 0.0) {
            fw = func(w);
            ++funcalls;
            if (fw < fc) {
                xa = xb;
                xb = w;
                fa = fb;
                fb = fw;
                break;
            }
            if (fw > fb) {
                xc = w;
                fc = fw;
                break;
            }
            w = xc + gold * (xc - xb);
            fw = func(w);
            ++funcalls;
        } else if ((w - wlim) * (wlim - xc) >= 0.0) {
            w = wlim;
            fw = func(w);
            ++funcalls;
        } else if ((w - wlim) * (xc - w) > 0.0) {
            fw = func(w);
            ++funcalls;
            if (fw < fc) {
                xb = xc;
                xc = w;
                w = xc + gold * (xc - xb);
                fb = fc;
                fc = fw;
                fw = func(w);
                ++funcalls;
            }
        } else {
            w = xc + gold * (xc - xb);
            fw = func(w);
            ++funcalls;
        }
        xa = xb;
        xb = xc;
        xc = w;
        fa = fb;
        fb = fc;
        fc = fw;
    }
    const bool cond1 = (fb < fc && fb <= fa) || (fb < fa && fb <= fc);
    const bool cond2 = (xa < xb && xb < xc) || (xc < xb && xb < xa);
    const bool cond3 = std::isfinite(xa) && std::isfinite(xb) && std::isfinite(xc);
    if (!(cond1 && cond2 && cond3)) {
        throw BracketFailure{xa, xb, xc, fa, fb, fc, funcalls};
    }
    return {xa, xb, xc, fa, fb, fc, funcalls};
}

// Brent.optimize() on bracket() from (0, 1), then _minimize_scalar_brent's success flag;
// BracketFailure handled as _recover_from_bracket_error does
template <class F>
BrentResult brent(F& func, double tol, std::size_t maxiter) {
    Bracket br{};
    try {
        br = bracket(func);
    } catch (const BracketFailure& e) {
        BrentResult out;
        const double xs[3] = {e.xa, e.xb, e.xc};
        const double fs[3] = {e.fa, e.fb, e.fc};
        bool anyNan = false;
        for (int k = 0; k < 3; ++k) {
            anyNan = anyNan || std::isnan(xs[k]) || std::isnan(fs[k]);
        }
        if (anyNan) {
            out.x = out.fun = std::numeric_limits<double>::quiet_NaN();
        } else {
            int imin = 0;  // np.argmin: the first minimum
            for (int k = 1; k < 3; ++k) {
                if (fs[k] < fs[imin]) {
                    imin = k;
                }
            }
            out.x = xs[imin];
            out.fun = fs[imin];
        }
        out.nfev = e.funcalls;
        out.nit = 0;
        out.success = false;
        return out;
    }
    constexpr double mintol = 1.0e-11;
    constexpr double cg = 0.3819660;
    double x = br.xb, w = br.xb, v = br.xb;
    double fw = br.fb, fv = br.fb, fx = br.fb;
    double a = 0.0, b = 0.0;
    if (br.xa < br.xc) {
        a = br.xa;
        b = br.xc;
    } else {
        a = br.xc;
        b = br.xa;
    }
    double deltax = 0.0;
    double rat = 0.0;
    std::size_t funcalls = br.funcalls;
    std::size_t iter = 0;
    while (iter < maxiter) {
        const double tol1 = tol * std::fabs(x) + mintol;
        const double tol2 = 2.0 * tol1;
        const double xmid = 0.5 * (a + b);
        if (std::fabs(x - xmid) < (tol2 - 0.5 * (b - a))) {
            break;
        }
        if (std::fabs(deltax) <= tol1) {
            if (x >= xmid) {
                deltax = a - x;  // a golden section step
            } else {
                deltax = b - x;
            }
            rat = cg * deltax;
        } else {  // a parabolic step
            const double tmp1 = (x - w) * (fx - fv);
            double tmp2 = (x - v) * (fx - fw);
            double p = (x - v) * tmp2 - (x - w) * tmp1;
            tmp2 = 2.0 * (tmp2 - tmp1);
            if (tmp2 > 0.0) {
                p = -p;
            }
            tmp2 = std::fabs(tmp2);
            const double dxTemp = deltax;
            deltax = rat;
            if ((p > tmp2 * (a - x)) && (p < tmp2 * (b - x)) && (std::fabs(p) < std::fabs(0.5 * tmp2 * dxTemp))) {
                rat = p * 1.0 / tmp2;  // the parabolic step is useful
                const double u = x + rat;
                if ((u - a) < tol2 || (b - u) < tol2) {
                    if (xmid - x >= 0) {
                        rat = tol1;
                    } else {
                        rat = -tol1;
                    }
                }
            } else {
                if (x >= xmid) {
                    deltax = a - x;  // not useful: a golden section step
                } else {
                    deltax = b - x;
                }
                rat = cg * deltax;
            }
        }
        double u = 0.0;
        if (std::fabs(rat) < tol1) {  // move by at least tol1
            if (rat >= 0) {
                u = x + tol1;
            } else {
                u = x - tol1;
            }
        } else {
            u = x + rat;
        }
        const double fu = func(u);
        ++funcalls;
        if (fu > fx) {
            if (u < x) {
                a = u;
            } else {
                b = u;
            }
            if ((fu <= fw) || (w == x)) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if ((fu <= fv) || (v == x) || (v == w)) {
                v = u;
                fv = fu;
            }
        } else {
            if (u >= x) {
                a = x;
            } else {
                b = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        ++iter;
    }
    BrentResult out;
    out.x = x;
    out.fun = fx;
    out.nit = iter;
    out.nfev = funcalls;
    out.success = iter < maxiter && !(std::isnan(x) || std::isnan(fx));
    return out;
}

}  // namespace

BrentResult minimizeScalarBrent(const ScalarFunction1d& f, double xtol, std::size_t maxiter) {
    if (xtol < 0) {
        throw std::invalid_argument("tolerance should be >= 0");
    }
    auto func = [&f](double t) { return f(t); };
    return brent(func, xtol, maxiter);
}

PowellResult minimizePowell(const ScalarObjective& f, std::span<const double> x0, const PowellOptions& options,
                            const PowellCallback& callback) {
    const std::size_t n = x0.size();
    if (n == 0) {
        throw std::invalid_argument("minimizePowell: x0 is empty");
    }
    std::size_t maxiter = 0, maxfun = 0;
    if (!options.maxiter && !options.maxfev) {
        maxiter = maxfun = n * 1000;
    } else if (!options.maxiter) {
        maxfun = *options.maxfev;
        maxiter = kUnlimited;
    } else if (!options.maxfev) {
        maxiter = *options.maxiter;
        maxfun = kUnlimited;
    } else {
        maxiter = *options.maxiter;
        maxfun = *options.maxfev;
    }

    std::size_t fcalls = 0;
    std::vector<double> scratch(n);
    auto func = [&](const std::vector<double>& x) {  // _wrap_scalar_function_maxfun_validation
        if (fcalls >= maxfun) {
            throw MaxFunctionCalls{};
        }
        ++fcalls;
        scratch = x;  // the objective gets a copy
        return f(std::span<const double>(scratch));
    };

    std::vector<std::vector<double>> direc;
    if (options.direc.empty()) {
        direc.assign(n, std::vector<double>(n, 0.0));
        for (std::size_t i = 0; i < n; ++i) {
            direc[i][i] = 1.0;
        }
    } else {
        if (options.direc.size() != n) {
            throw std::invalid_argument("minimizePowell: direc needs one row per variable");
        }
        for (const auto& row : options.direc) {
            if (row.size() != n) {
                throw std::invalid_argument("minimizePowell: every row of direc needs one number per variable");
            }
        }
        direc = options.direc;
    }

    // _linesearch_powell(func, p, xi, tol=xtol * 100, fval=fval) without bounds: (fret, p + xi, xi)
    // with xi scaled to the step taken
    std::vector<double> point(n);
    auto linesearch = [&](double& fval, std::vector<double>& p, std::vector<double>& xi) {
        bool any = false;
        for (double v : xi) {
            any = any || v != 0.0;
        }
        if (!any) {
            return;  // a zero direction: nothing to optimise
        }
        auto along = [&](double alpha) {  // func(p + alpha * xi)
            for (std::size_t k = 0; k < n; ++k) {
                point[k] = p[k] + alpha * xi[k];
            }
            return func(point);
        };
        const BrentResult res = brent(along, options.xtol * 100, 500);
        for (std::size_t k = 0; k < n; ++k) {
            xi[k] = res.x * xi[k];
        }
        for (std::size_t k = 0; k < n; ++k) {
            p[k] = p[k] + xi[k];
        }
        fval = res.fun;
    };

    std::vector<double> x(x0.begin(), x0.end());
    PowellResult result;
    double fval = func(x);
    std::vector<double> x1 = x;
    std::size_t iter = 0;
    try {
        while (true) {
            const double fx = fval;
            std::size_t bigind = 0;
            double delta = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                std::vector<double> direc1 = direc[i];
                const double fx2 = fval;
                // the line search updates fval and x only when it completes
                double fnew = fval;
                std::vector<double> xnew = x;
                linesearch(fnew, xnew, direc1);
                fval = fnew;
                x = std::move(xnew);
                if ((fx2 - fval) > delta) {
                    delta = fx2 - fval;
                    bigind = i;
                }
            }
            ++iter;
            if (callback && callback(std::span<const double>(x), fval)) {
                break;
            }
            const double bnd = options.ftol * (std::fabs(fx) + std::fabs(fval)) + 1e-20;
            if (2.0 * (fx - fval) <= bnd) {
                break;
            }
            if (fcalls >= maxfun) {
                break;
            }
            if (iter >= maxiter) {
                break;
            }
            if (std::isnan(fx) && std::isnan(fval)) {
                break;  // ended up in a NaN region
            }
            // the extrapolated point
            std::vector<double> direc1(n);
            for (std::size_t k = 0; k < n; ++k) {
                direc1[k] = x[k] - x1[k];
            }
            x1 = x;
            std::vector<double> x2(n);
            for (std::size_t k = 0; k < n; ++k) {
                x2[k] = x[k] + 1.0 * direc1[k];  // min(lmax, 1) with lmax = 1
            }
            const double fx2 = func(x2);
            if (fx > fx2) {
                double t = 2.0 * (fx + fx2 - 2.0 * fval);
                double temp = (fx - fval - delta);
                t *= temp * temp;
                temp = fx - fx2;
                t -= delta * temp * temp;
                if (t < 0.0) {
                    double fnew = fval;
                    std::vector<double> xnew = x;
                    linesearch(fnew, xnew, direc1);
                    fval = fnew;
                    x = std::move(xnew);
                    bool any = false;
                    for (double v : direc1) {
                        any = any || v != 0.0;
                    }
                    if (any) {
                        direc[bigind] = direc[n - 1];
                        direc[n - 1] = direc1;
                    }
                }
            }
        }
    } catch (const MaxFunctionCalls&) {
        // the point and value before the interrupted step stand
    }

    result.status = 0;
    result.message = "Optimization terminated successfully.";
    bool nanX = false;
    for (double v : x) {
        nanX = nanX || std::isnan(v);
    }
    if (fcalls >= maxfun) {
        result.status = 1;
        result.message = "Maximum number of function evaluations has been exceeded.";
    } else if (iter >= maxiter) {
        result.status = 2;
        result.message = "Maximum number of iterations has been exceeded.";
    } else if (std::isnan(fval) || nanX) {
        result.status = 3;
        result.message = "NaN result encountered.";
    }
    result.success = result.status == 0;
    result.x = std::move(x);
    result.fun = fval;
    result.direc = std::move(direc);
    result.nit = iter;
    result.nfev = fcalls;
    return result;
}

}  // namespace mdbridge::numerics
