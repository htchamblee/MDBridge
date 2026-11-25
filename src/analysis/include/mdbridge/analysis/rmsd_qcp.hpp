#pragma once
// The RMSD after optimal superposition, as MDATK's rmsd analysis gets it from MDAnalysis 2.11.0
// (rms.rmsd(a, b, weights=w, center=True, superposition=True); docs/specs/m3-rmsd.md). Both
// structures are centred (with the raw weights), the weights normalised to mean 1, and the minimum
// over proper rotations of sqrt(sum_i w_i |R a_i - b_i|^2 / n) found with the quaternion
// characteristic polynomial (QCP) method of Theobald (Acta Cryst. A 61 (2005) 478) and Liu,
// Agrafiotis & Theobald (J. Comput. Chem. 31 (2010) 1561), evaluated in f64 in exactly the spec's
// order of operations without fused multiply-adds, so that results match MDAnalysis bit for bit
// for non-degenerate inputs.
//
// Deliberate differences (docs/DIFFERENCES.md): for degenerate inputs (fewer than 3 atoms, a
// collinear or coincident set by its weighted inner-product matrix, or no convergence in 50 Newton
// steps; AN-16) and where E0 - lambda < 64 eps E0 is lost to cancellation (AN-17), the value is the
// exact superposition RMSD: the optimal rotation from a Jacobi eigen-solution of the 4 x 4 key
// matrix, then the residual summed directly after rotating; one atom and identical structures give
// exactly 0. MDAnalysis gives NaN or values off by up to ~10 Å there (5e-6 Å for identical sets).
// Non-finite inputs, negative weights and empty sets keep MDAnalysis' values.
//
// Not reproduced: numpy's RuntimeWarnings for non-finite or empty inputs (MDATK ignores them).

#include <cstddef>
#include <optional>
#include <span>
#include <string>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::analysis {

// Python's ZeroDivisionError (not a ValueError): the weights sum to zero.
class ZeroDivisionError : public core::NamedError {
public:
    explicit ZeroDivisionError(const std::string& message) : core::NamedError("ZeroDivisionError", message) {}
};

// A row-major f64 array with the shape numpy would report (for the checks of §2 of the spec).
struct ShapedArray {
    std::span<const double> values;
    std::span<const std::size_t> shape;
};

// rmsd(a, b, weights, center=True, superposition=True) of n x 3 arrays (row-major, a and b
// typically the current and the reference positions). Errors, in this order:
// * a and b of different shapes: core::ValueError "a and b must have same shape";
// * weights whose shape is not (n,): core::ValueError
//   "Shape of weights must be consistent with shape of a along specified axis.";
// * weights summing to exactly 0: ZeroDivisionError "Weights sum to zero, can't be normalized".
// Equal shapes other than (n, 3) throw std::invalid_argument (not specified). n = 0 without
// weights gives NaN; non-finite inputs give NaN.
double superposedRmsd(const ShapedArray& a, const ShapedArray& b, const std::optional<ShapedArray>& weights);

// The same for flat n x 3 arrays (sizes multiples of 3, else std::invalid_argument; different sizes
// are the shape error) and 1-D weights.
double superposedRmsd(std::span<const double> a, std::span<const double> b,
                      std::optional<std::span<const double>> weights = std::nullopt);

}  // namespace mdbridge::analysis
