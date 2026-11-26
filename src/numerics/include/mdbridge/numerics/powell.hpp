#pragma once
// scipy.optimize.minimize(method="Powell") without bounds, with its line search (Brent's method on a
// bracket found by scipy.optimize.bracket), ported from SciPy 1.18.1 (scipy/optimize/_optimize.py:
// _minimize_powell, _linesearch_powell, Brent, bracket, _recover_from_bracket_error; BSD licence and
// Travis E. Oliphant's notice, see powell.cpp). Every operation is done in SciPy's order in double
// precision, so the iterates are SciPy's for an objective that rounds as the Python one does.
// MDATK's density fit refines poses with it (FitProblem.refine: xtol, ftol 1e-7, maxfev 4000, direc).

#include <cstddef>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace mdbridge::numerics {

// f(x) of the minimiser; it gets a copy of the current point
using ScalarObjective = std::function<double(std::span<const double>)>;
// called after each Powell iteration with the point and its value; returning true stops the
// minimisation (SciPy: a callback raising StopIteration)
using PowellCallback = std::function<bool(std::span<const double>, double)>;

struct PowellOptions {
    double xtol = 1e-4;  // the line searches' relative tolerance is 100 xtol (Brent's tol)
    double ftol = 1e-4;
    // SciPy's defaults: both unset -> 1000 n each; one unset -> unlimited, unless the other is
    // unlimited (std::nullopt here never means "infinite" for both)
    std::optional<std::size_t> maxiter;
    std::optional<std::size_t> maxfev;
    // initial directions, one row per direction (n rows of n numbers); empty: the identity
    std::vector<std::vector<double>> direc;
};

struct PowellResult {
    std::vector<double> x;
    double fun = 0.0;
    std::vector<std::vector<double>> direc;  // the final direction set
    std::size_t nit = 0;
    std::size_t nfev = 0;
    int status = 0;  // 0 success, 1 maxfev, 2 maxiter, 3 NaN
    bool success = true;
    std::string message;  // SciPy's _status_message text
};

// minimize(f, x0, method="Powell", options=...): throws std::invalid_argument for an empty x0 or
// directions of the wrong shape, std::runtime_error with SciPy's message when bracket() runs out of
// iterations ("No valid bracket was found before the iteration limit was reached. ...").
PowellResult minimizePowell(const ScalarObjective& f, std::span<const double> x0, const PowellOptions& options = {},
                            const PowellCallback& callback = {});

// minimize_scalar(f, method="brent", options={"xtol": xtol, "maxiter": maxiter}) with the bracket
// search from (0, 1): x, fun, nit, nfev; success false with x / fun from the best bracket point when
// no valid bracket is found (_recover_from_bracket_error).
using ScalarFunction1d = std::function<double(double)>;
struct BrentResult {
    double x = 0.0;
    double fun = 0.0;
    std::size_t nit = 0;
    std::size_t nfev = 0;
    bool success = true;
};
BrentResult minimizeScalarBrent(const ScalarFunction1d& f, double xtol = 1.48e-8, std::size_t maxiter = 500);

}  // namespace mdbridge::numerics
