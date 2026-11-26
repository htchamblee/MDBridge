// LSMR as scipy.sparse.linalg.lsmr computes it, ported from SciPy 1.18.1
// (scipy/sparse/linalg/_isolve/lsmr.py, and _sym_ortho from lsqr.py). LSMR is David Fong and
// Michael Saunders' method and code (Copyright (C) 2010 David Fong and Michael Saunders; LSQR's
// _sym_ortho: Copyright (c) 2006, Systems Optimization Laboratory), distributed with SciPy under its
// BSD licence and LSQR's:
//
// Copyright (c) 2006, Systems Optimization Laboratory
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright notice, this list of
//       conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright notice, this list of
//       conditions and the following disclaimer in the documentation and/or other materials
//       provided with the distribution.
//
//     * Neither the name of Stanford University nor the names of its contributors may be used to
//       endorse or promote products derived from this software without specific prior written
//       permission.
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
// (SciPy's own licence: Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers, BSD-3-Clause,
// as in least_squares.cpp.) Real double precision only, without x0 and without printing. Built
// without floating-point contraction.

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "mdbridge/numerics/least_squares.hpp"
#include "numerics_internal.hpp"

namespace mdbridge::numerics {

namespace {

double sign(double a) {  // np.sign
    if (std::isnan(a)) {
        return a;
    }
    return a > 0 ? 1.0 : a < 0 ? -1.0 : 0.0;
}

struct Rotation {
    double c, s, r;
};

// _sym_ortho(a, b): a stable Givens rotation
Rotation symOrtho(double a, double b) {
    if (b == 0) {
        return {sign(a), 0, std::fabs(a)};
    }
    if (a == 0) {
        return {0, sign(b), std::fabs(b)};
    }
    if (std::fabs(b) > std::fabs(a)) {
        const double tau = a / b;
        const double s = sign(b) / std::sqrt(1 + tau * tau);
        const double c = s * tau;
        const double r = b / s;
        return {c, s, r};
    }
    const double tau = b / a;
    const double c = sign(a) / std::sqrt(1 + tau * tau);
    const double s = c * tau;
    const double r = a / c;
    return {c, s, r};
}

}  // namespace

LsmrResult lsmr(std::size_t m, std::size_t n, const MatVec& matvec, const MatVec& rmatvec, std::span<const double> b,
                double damp, double atol, double btol, double conlim, std::optional<std::size_t> maxiterOpt) {
    if (b.size() != m) {
        throw std::invalid_argument("lsmr: b does not match the operator");
    }
    const std::size_t maxiter = maxiterOpt ? *maxiterOpt : std::min(m, n);
    std::vector<double> u(b.begin(), b.end());
    const double normb = detail::norm2(b);
    LsmrResult out;
    std::vector<double>& x = out.x;
    x.assign(n, 0.0);
    double beta = normb;
    std::vector<double> v;
    double alpha = 0;
    if (beta > 0) {
        for (double& e : u) {
            e = (1 / beta) * e;
        }
        v = rmatvec(u);
        alpha = detail::norm2(v);
    } else {
        v.assign(n, 0.0);
        alpha = 0;
    }
    if (alpha > 0) {
        for (double& e : v) {
            e = (1 / alpha) * e;
        }
    }
    std::size_t itn = 0;
    double zetabar = alpha * beta;
    double alphabar = alpha;
    double rho = 1;
    double rhobar = 1;
    double cbar = 1;
    double sbar = 0;
    std::vector<double> h = v;
    std::vector<double> hbar(n, 0.0);
    double betadd = beta;
    double betad = 0;
    double rhodold = 1;
    double tautildeold = 0;
    double thetatilde = 0;
    double zeta = 0;
    double d = 0;
    double normA2 = alpha * alpha;
    double maxrbar = 0;
    double minrbar = 1e+100;
    double normA = std::sqrt(normA2);
    double condA = 1;
    double normx = 0;
    int istop = 0;
    double ctol = 0;
    if (conlim > 0) {
        ctol = 1 / conlim;
    }
    double normr = beta;
    double normar = alpha * beta;
    auto finish = [&]() {
        out.istop = istop;
        out.itn = itn;
        out.normr = normr;
        out.normar = normar;
        out.normA = normA;
        out.condA = condA;
        out.normx = normx;
        return out;
    };
    if (normar == 0) {
        return finish();
    }
    if (normb == 0) {
        std::fill(x.begin(), x.end(), 0.0);
        return finish();
    }
    while (itn < maxiter) {
        itn = itn + 1;
        // the next step of the bidiagonalisation: beta u = A v - alpha u, alpha v = A' u - beta v
        for (double& e : u) {
            e *= -alpha;
        }
        const std::vector<double> av = matvec(v);
        for (std::size_t i = 0; i < m; ++i) {
            u[i] += av[i];
        }
        beta = detail::norm2(u);
        if (beta > 0) {
            for (double& e : u) {
                e *= (1 / beta);
            }
            for (double& e : v) {
                e *= -beta;
            }
            const std::vector<double> atu = rmatvec(u);
            for (std::size_t i = 0; i < n; ++i) {
                v[i] += atu[i];
            }
            alpha = detail::norm2(v);
            if (alpha > 0) {
                for (double& e : v) {
                    e *= (1 / alpha);
                }
            }
        }
        // the rotation Qhat_{k, 2k+1}
        const Rotation hat = symOrtho(alphabar, damp);
        const double chat = hat.c, shat = hat.s, alphahat = hat.r;
        // the rotation Q_{k, k+1}
        const double rhoold = rho;
        const Rotation q = symOrtho(alphahat, beta);
        const double c = q.c, s = q.s;
        rho = q.r;
        const double thetanew = s * alpha;
        alphabar = c * alpha;
        // the rotation Qbar_{k, k+1}
        const double rhobarold = rhobar;
        const double zetaold = zeta;
        const double thetabar = sbar * rho;
        const double rhotemp = cbar * rho;
        const Rotation qbar = symOrtho(cbar * rho, thetanew);
        cbar = qbar.c;
        sbar = qbar.s;
        rhobar = qbar.r;
        zeta = cbar * zetabar;
        zetabar = -sbar * zetabar;
        // h, hbar, x
        const double hbarFactor = -(thetabar * rho / (rhoold * rhobarold));
        for (std::size_t i = 0; i < n; ++i) {
            hbar[i] *= hbarFactor;
            hbar[i] += h[i];
        }
        const double xFactor = (zeta / (rho * rhobar));
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += xFactor * hbar[i];
        }
        const double hFactor = -(thetanew / rho);
        for (std::size_t i = 0; i < n; ++i) {
            h[i] *= hFactor;
            h[i] += v[i];
        }
        // estimate of ||r||: apply the rotation Qhat_{k, 2k+1}, then Q_{k, k+1}
        const double betaacute = chat * betadd;
        const double betacheck = -shat * betadd;
        const double betahat = c * betaacute;
        betadd = -s * betaacute;
        // the rotation Qtilde_{k-1}
        const double thetatildeold = thetatilde;
        const Rotation qt = symOrtho(rhodold, thetabar);
        const double ctildeold = qt.c, stildeold = qt.s, rhotildeold = qt.r;
        thetatilde = stildeold * rhobar;
        rhodold = ctildeold * rhobar;
        betad = -stildeold * betad + ctildeold * betahat;
        tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold;
        const double taud = (zeta - thetatilde * tautildeold) / rhodold;
        d = d + betacheck * betacheck;
        normr = std::sqrt(d + (betad - taud) * (betad - taud) + betadd * betadd);
        // estimate of ||A||
        normA2 = normA2 + beta * beta;
        normA = std::sqrt(normA2);
        normA2 = normA2 + alpha * alpha;
        // estimate of cond(A)
        maxrbar = std::max(maxrbar, rhobarold);
        if (itn > 1) {
            minrbar = std::min(minrbar, rhobarold);
        }
        condA = std::max(maxrbar, rhotemp) / std::min(minrbar, rhotemp);
        // the stopping tests
        normar = std::fabs(zetabar);
        normx = detail::norm2(x);
        const double test1 = normr / normb;
        double test2 = 0;
        if ((normA * normr) != 0) {
            test2 = normar / (normA * normr);
        } else {
            test2 = std::numeric_limits<double>::infinity();
        }
        const double test3 = 1 / condA;
        const double t1 = test1 / (1 + normA * normx / normb);
        const double rtol = btol + atol * normA * normx / normb;
        if (itn >= maxiter) {
            istop = 7;
        }
        if (1 + test3 <= 1) {
            istop = 6;
        }
        if (1 + test2 <= 1) {
            istop = 5;
        }
        if (1 + t1 <= 1) {
            istop = 4;
        }
        if (test3 <= ctol) {
            istop = 3;
        }
        if (test2 <= atol) {
            istop = 2;
        }
        if (test1 <= rtol) {
            istop = 1;
        }
        if (istop > 0) {
            break;
        }
    }
    return finish();
}

}  // namespace mdbridge::numerics
