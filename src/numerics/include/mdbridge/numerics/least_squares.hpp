#pragma once
// scipy.optimize.least_squares(method="trf") without bounds and with the linear loss, the Jacobian
// estimated by 2-point finite differences (dense, or with a sparsity structure and SciPy's column
// groups), ported from SciPy 1.18.1 (scipy/optimize/_lsq/least_squares.py, trf.py (trf_no_bounds),
// common.py, scipy/optimize/_numdiff.py, _group_columns.pyx; BSD licence, see least_squares.cpp) with
// LSMR (scipy/sparse/linalg/_isolve/lsmr.py, Fong and Saunders' BSD code, lsmr.cpp) for
// tr_solver="lsmr". Owner decision D4 (docs/specs/m6-plan.md): MDATK's helix builder solves its pair
// fits, helix closures and level closures with it (mdatk/builder.py), so MDBridge follows the same
// iterates.
//
// Linear algebra that SciPy hands to LAPACK is MDBridge's own: the SVD of tr_solver="exact"
// (one-sided Jacobi rotations; for a tall Jacobian of its rows that are not all zero, after a
// Householder QR, as LAPACK's dgesdd also reduces a tall matrix first), the QR of the 2-D subspace (Householder reflections with LAPACK's sign
// convention) and np.roots of the subspace's quartic (its real roots by bracketing and Newton steps,
// not the eigenvalues of the companion matrix). They agree with LAPACK to rounding, and BLAS sums
// (J^T f, norms) are taken in index order, so the iterates agree with SciPy's to rounding, not bit for
// bit: the finite differences amplify an ulp in a residual to about 1e-8 relative in the next
// Jacobian (steps of sqrt(eps)), and SciPy's own iterates move by that much when its residuals move by
// an ulp. tests/numerics/test_optimize.cpp holds the port to 1e-9, or to four times that spread of
// SciPy's own where it is larger (at most 2.4e-8 relative on its smooth cases).

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace mdbridge::numerics {

// residuals f(x) (m values; m fixed by the first call)
using ResidualFunction = std::function<std::vector<double>(std::span<const double>)>;

// A sparse matrix in compressed sparse row form (scipy.sparse.csr_array with sorted column indices).
struct CsrMatrix {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<std::size_t> indptr;   // rows + 1 offsets into indices / data
    std::vector<std::size_t> indices;  // column of each stored value, ascending within a row
    std::vector<double> data;

    // csr_matrix((ones, (rowOf, colOf)), shape=(rows, cols)): the structure of the given entries
    // (duplicates stored once), every value 1; std::invalid_argument for an entry outside the shape
    static CsrMatrix pattern(std::size_t rows, std::size_t cols, std::span<const std::size_t> rowOf,
                             std::span<const std::size_t> colOf);
};

// tr_solver: "exact" (dense SVD) or "lsmr"; Auto as SciPy chooses (exact for a dense Jacobian, lsmr
// with a sparsity structure)
enum class TrustRegionSolver { Auto, Exact, Lsmr };

struct LsmrOptions {  // tr_options of tr_solver="lsmr" (lsmr's keywords and TRF's own)
    double atol = 1e-6;
    double btol = 1e-6;
    double conlim = 1e8;
    std::optional<std::size_t> maxiter;  // lsmr's default: min(m, n)
    double damp = 0.0;
    bool regularize = true;
};

struct LeastSquaresOptions {
    double ftol = 1e-8;
    double xtol = 1e-8;
    double gtol = 1e-8;
    // x_scale: empty = 1 for every variable; one number per variable otherwise (positive, finite)
    std::vector<double> xScale;
    bool xScaleJac = false;  // x_scale="jac" (then xScale must be empty)
    std::optional<std::size_t> maxNfev;  // default 100 n
    std::optional<double> diffStep;      // diff_step (relative), default sqrt(eps)
    TrustRegionSolver trSolver = TrustRegionSolver::Auto;
    LsmrOptions lsmr;
    // jac_sparsity: the structure of the Jacobian (rows = residuals, cols = variables), or none (dense)
    std::optional<CsrMatrix> jacSparsity;
    // With jacSparsity: the structure only groups the finite differences (SciPy's column groups, one
    // evaluation per group) and the Jacobian is then dense, solved with tr_solver "exact" (trSolver
    // Auto or Exact). SciPy refuses tr_solver="exact" with jac_sparsity; this is its
    // least_squares(jac=lambda x: approx_derivative(fun, x, sparsity=...).toarray(), tr_solver="exact").
    // MDBridge's extension, for the helix builder's strand closure (docs/DIFFERENCES.md CH-16).
    bool denseFromSparsity = false;
    // finite-difference evaluations run on up to this many threads (each column or group of columns
    // is one evaluation; the Jacobian does not depend on the count). The residual function must
    // then be safe to call concurrently. 1: in the caller's thread.
    std::size_t jacobianThreads = 1;
};

// called after each iteration (the point, its cost, iteration count, nfev); returning true stops
// with status -2
using LeastSquaresCallback = std::function<bool(std::span<const double> x, double cost, std::size_t nit,
                                                std::size_t nfev)>;

struct LeastSquaresResult {
    std::vector<double> x;
    double cost = 0.0;          // 0.5 * sum(f^2)
    std::vector<double> fun;    // f(x)
    std::vector<double> grad;   // J^T f
    double optimality = 0.0;    // max |grad| at the last check
    std::size_t nfev = 0;       // evaluations of f by TRF (not those of the finite differences)
    std::size_t njev = 0;
    std::size_t nit = 0;
    int status = 0;             // -2 callback, 0 max_nfev, 1 gtol, 2 ftol, 3 xtol, 4 ftol and xtol
    bool success = false;       // status > 0
    std::string message;        // SciPy's TERMINATION_MESSAGES text
};

// least_squares(fun, x0, method="trf", jac="2-point", ...): std::invalid_argument with SciPy's
// messages for bad input ("Residuals are not finite in the initial point.", "`x_scale` must be 'jac'
// or array_like with positive numbers.", "tr_solver='exact' is incompatible with `jac_sparsity`.", ...).
LeastSquaresResult leastSquares(const ResidualFunction& fun, std::span<const double> x0,
                                const LeastSquaresOptions& options = {}, const LeastSquaresCallback& callback = {});

// The pieces, exposed for tests --------------------------------------------------------------------

// scipy.optimize._numdiff.group_columns(structure, order=0): the group of every column (columns of a
// group share no row), the columns visited in numpy.random.RandomState(0).permutation(n) order
std::vector<std::size_t> groupColumns(const CsrMatrix& structure);

// numpy.random.RandomState(seed).permutation(n) (MT19937, legacy seeding and shuffle)
std::vector<std::size_t> randomStatePermutation(std::uint32_t seed, std::size_t n);

// lsmr(A, b, damp, atol, btol, conlim, maxiter) for A given by its products: x, istop, itn
struct LsmrResult {
    std::vector<double> x;
    int istop = 0;
    std::size_t itn = 0;
    double normr = 0.0;
    double normar = 0.0;
    double normA = 0.0;
    double condA = 0.0;
    double normx = 0.0;
};
using MatVec = std::function<std::vector<double>(std::span<const double>)>;
LsmrResult lsmr(std::size_t m, std::size_t n, const MatVec& matvec, const MatVec& rmatvec, std::span<const double> b,
                double damp = 0.0, double atol = 1e-6, double btol = 1e-6, double conlim = 1e8,
                std::optional<std::size_t> maxiter = std::nullopt);

// the real roots of the polynomial c[0] x^d + ... + c[d] (degree at most 4 after leading zeros are
// dropped), ascending, a root at 0 for every trailing zero (np.roots(c)[np.isreal(...)] up to order)
std::vector<double> realPolynomialRoots(std::span<const double> c);

}  // namespace mdbridge::numerics
