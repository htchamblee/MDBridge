#pragma once
// MDATK's spatial helpers (mdatk/spatial.py): the neighbour search every geometric analysis uses
// (pairs_within), the minimum-image convention for difference vectors (minimum_image_vectors),
// the box filter (_usable_box) and the optimal superposition (kabsch).
//
// Points are x, y, z interleaved (n * 3 values), float32 as a Frame holds them or float64 as the
// analyses compute them; sizes that are not a multiple of 3 throw std::invalid_argument.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/rigid_transform.hpp"  // Vec3, Mat3

namespace mdbridge::analysis {

using core::Mat3;
using core::Vec3;

// A periodic cell as MDAnalysis' ts.dimensions gives it: lengths a, b, c (Å), then the angles
// alpha, beta, gamma (degrees), widened to float64.
using Box = std::array<double, 6>;

// The pairs (i[k], j[k]) with |a[i] - b[j]| <= cutoff and their distances d[k] (int64, int64,
// float64 as MDATK returns them).
struct PairList {
    std::vector<std::int64_t> i;
    std::vector<std::int64_t> j;
    std::vector<double> d;

    std::size_t size() const { return i.size(); }
    bool empty() const { return i.empty(); }
};

// _usable_box: the box when its three lengths are positive, else std::nullopt (a missing box, or
// lengths <= 0). As in MDATK a NaN length passes (NaN <= 0 is false) and angles are not checked.
// The span form also returns std::nullopt unless it holds exactly 6 values (MDATK's shape check).
std::optional<Box> usableBox(std::span<const double> box);
std::optional<Box> usableBox(const std::optional<Box>& box);
std::optional<Box> usableBox(const std::optional<std::array<float, 6>>& dims);  // ts.dimensions

// pairs_within(a, b, cutoff, box): every pair within the cutoff.
//
// Without a usable box (MDATK's default, use_pbc = false): scipy's cKDTree path. The points are
// widened to float64 (float32 positions exactly), the squared distance is
// ((dx * dx + dy * dy) + dz * dz) with d = a[i] - b[j] in float64, a pair is kept when it is
// <= cutoff * cutoff (scipy compares squares, not sqrt(d2) <= cutoff: they disagree for some pairs
// within an ulp of the cutoff), and the distance is sqrt of that sum. Hence a negative cutoff acts
// as its absolute value, an infinite one keeps every pair and a NaN one none. Empty a or b give no
// pairs; otherwise a non-finite coordinate in a or b throws core::ValueError("data must be finite,
// check for nan or inf values") (scipy's message).
//
// Order (open path only): MDATK's order is cKDTree's traversal order, which depends on the tree's
// construction (std::nth_element of the C++ library scipy was built with), so MDBridge does not
// reproduce it. MDBridge returns the pairs sorted by (i, j), i.e. by a index, then by b index. The pair set and
// every distance equal MDATK's exactly; consumers must not depend on the order beyond what MDATK's
// analyses need. In MDATK the order matters only to: float sums over pairs (energy elec / vdw
// sums, net force: equal within rounding), the H-bond hydrogen guess (an exact tie between two
// donors goes to the first pair: here the lowest donor index), the stacking pair map (CSR entries
// within a frame in pair order) and atom_bonds' distance fallback (bond rows in pair order); see
// docs/specs/m3-plan.md and the W0.7 notes.
//
// With a usable box (use_pbc): MDAnalysis' capped_distance on float32 coordinates (float64 points
// are rounded to float32 first), which pairsWithinPeriodic implements, in MDAnalysis' order (not
// sorted; see there).
PairList pairsWithin(std::span<const double> a, std::span<const double> b, double cutoff,
                     const std::optional<Box>& box = std::nullopt);
PairList pairsWithin(std::span<const float> a, std::span<const float> b, double cutoff,
                     const std::optional<Box>& box = std::nullopt);

// The open-boundary search alone (the path above without a box).
PairList pairsWithinOpen(std::span<const double> a, std::span<const double> b, double cutoff);

// The periodic search: pairs within the cutoff under periodic boundaries, as MDAnalysis'
// capped_distance(a, b, max_cutoff=cutoff, box=box, return_distances=True) returns them, written
// clean-room from docs/specs/m3-capped-distance.md (section numbers below refer to it). box is a
// usable box; the cell is rounded to float32 and is orthogonal iff its three float32 angles are
// exactly 90, otherwise it becomes float32 lower-triangular box vectors (all zero, "invalid", for
// angles <= 0 or >= 180, NaN, lengths <= 0 or no third vector: plain distances).
//
// The method is chosen as MDAnalysis chooses it (periodicSearchUsesGrid): brute force or a cell
// grid. The two differ in rounding (brute force on an orthogonal cell scales by the float32
// reciprocal of each length, so even pairs that need no image are off by up to 2^-24 relative),
// in the boundary test (brute force keeps sqrt(d2) <= cutoff, the grid d2 <= cutoff * cutoff) and
// in their errors. Each (i, j) occurs at most once, with its nearest-image distance, even for
// cutoffs beyond half the cell. A point with a NaN or infinite coordinate is in no pair (MDAnalysis
// is only consistent about this for orthogonal cells). Errors, grid method only:
// core::ValueError("Cutoff must be positive") for a negative cutoff, "Cutoff <c> too large for box
// (max <M>)" beyond the grid's limit, "Any of the box dimensions cannot be 0" for a length that
// is 0 in float32. A NaN cutoff gives no pairs.
//
// Order: MDAnalysis' order, NOT the (i, j) order of the open path, because the analyses' outputs
// that depend on pair order (energy sums, the H-bond hydrogen tie, stacking CSR entries, atom-bond
// rows) must follow MDATK's. Brute force: i ascending, then j ascending. Grid: i ascending; for
// each i the 27 grid cells around its cell (x offset outermost, z innermost, -1 / 0 / +1, the grid
// wrapping periodically), and within a cell j DESCENDING.
//
// Precision: orthogonal cells reproduce MDAnalysis bit for bit (pairs, order, distances).
// Triclinic cells are wrapped in single precision as the spec describes; every probe case of
// tests/data/analysis/distance_cases.json matches bit for bit too (test_periodic.cpp), but the
// spec only guarantees that wrap to within 8 ulp32 of the largest of the pair's coordinates and
// the cell lengths: where it differs, pairs within that band of the cutoff may be in or out and a
// point moved to a neighbouring grid cell changes the order of its reference point's pairs.
// Deliberate deviations where MDAnalysis' result is undefined or meaningless: a point more than
// 2^31 cells away along a vector is not moved along it (MDAnalysis overflows a 32-bit integer);
// a triclinic box vector with an infinite component is not periodic (MDAnalysis reports distinct
// points at distance 0).
PairList pairsWithinPeriodic(std::span<const float> a, std::span<const float> b, double cutoff,
                             const Box& box);

// The method rule of pairsWithinPeriodic for n reference and m configuration points: brute force
// (false) when n < 10 or m < 10; the grid (true) when n * m >= 1e8; otherwise brute force when
// float32(cutoff) > float32(0.3) * S on any axis (float32 product, strictly greater), S the three
// float32 lengths of an orthogonal cell or, per Cartesian axis, the extent of the three box
// vectors' components (0 for an invalid cell), and the grid otherwise.
bool periodicSearchUsesGrid(std::size_t n, std::size_t m, double cutoff, const Box& box);

// minimum_image_vectors(vectors, box): with a usable box whose three angles are all close to 90
// degrees (np.allclose: |angle - 90| <= 1e-8 + 1e-5 * 90), v - L * round(v / L) per component
// with round half to even (np.round); otherwise (no box, unusable, any other angle: triclinic
// cells are not handled, as in MDATK) the vectors unchanged.
std::vector<double> minimumImageVectors(std::span<const double> vectors,
                                        const std::optional<Box>& box);
Vec3 minimumImageVector(const Vec3& vector, const std::optional<Box>& box);

// np.linalg.LinAlgError (a ValueError in Python)
class LinAlgError : public core::ValueError {
public:
    using core::ValueError::ValueError;
    const char* pythonName() const noexcept override { return "LinAlgError"; }
};

// kabsch(mobile, reference, weights): the optimal (optionally weighted) rotation superposing
// mobile onto reference; (x - mobileCenter) R^T + referenceCenter maps mobile coordinates onto the
// reference. Weights are normalised by their sum (numpy's pairwise sum); without weights every
// point weighs 1 / max(1, n). Centres are the weighted means; fewer than 3 points give the
// identity. R = V diag(1, 1, d) U^T from the SVD H = U S V^T of the weighted covariance H, with
// d = sign(det(V U^T)) (1 when 0): never a reflection. Throws std::invalid_argument when the sizes
// differ, LinAlgError("SVD did not converge") when H is not finite (as numpy's SVD does). Equals
// numpy's result within a few ulp (LAPACK / BLAS differ in the last bits) except for degenerate
// sets (collinear points, rank 1), whose rotation is not unique.
struct KabschResult {
    Mat3 rotation{};
    Vec3 mobileCenter{};
    Vec3 referenceCenter{};
};
KabschResult kabsch(std::span<const Vec3> mobile, std::span<const Vec3> reference,
                    std::optional<std::span<const double>> weights = std::nullopt);

}  // namespace mdbridge::analysis
