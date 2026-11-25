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
// pairs. A point with a non-finite coordinate is in no pair (AN-20; MDATK: scipy's ValueError
// "data must be finite, check for nan or inf values").
//
// Order (open path only): MDATK's order is cKDTree's traversal order, which depends on the tree's
// construction (std::nth_element of the C++ library scipy was built with), so MDBridge does not
// reproduce it. MDBridge returns the pairs sorted by (i, j), i.e. by a index, then by b index. The pair set and
// every distance equal MDATK's exactly; consumers must not depend on the order beyond what MDATK's
// analyses need. In MDATK the order matters only to: float sums over pairs (energy elec / vdw
// sums, net force: equal within rounding), the H-bond hydrogen guess (an exact tie between two
// donors goes to the first pair: here the lowest donor index), the stacking pair map (CSR entries
// within a frame in pair order) and atom_bonds' distance fallback (bond rows in pair order); see
// docs/specs/m3-plan.md (risks: pair order) and docs/DIFFERENCES.md.
//
// With a usable box (use_pbc): pairsWithinPeriodic on the points rounded to float32 (as MDATK
// passes them to its periodic search), also sorted by (i, j).
PairList pairsWithin(std::span<const double> a, std::span<const double> b, double cutoff,
                     const std::optional<Box>& box = std::nullopt);
PairList pairsWithin(std::span<const float> a, std::span<const float> b, double cutoff,
                     const std::optional<Box>& box = std::nullopt);

// The open-boundary search alone (the path above without a box).
PairList pairsWithinOpen(std::span<const double> a, std::span<const double> b, double cutoff);

// A periodic cell for minimum-image geometry (MDBridge's own; docs/DIFFERENCES.md AN-01 / AN-04).
// The six values are rounded to float32 (as trajectories store them) and widened again. The cell
// is orthogonal iff its three float32 angles are exactly 90 (the one orthogonality test for pairs
// and vectors alike); then each axis whose length is finite and > 0 is periodic (the others are
// not). Otherwise its float64 box vectors are a = (a, 0, 0), b = (b cos g, b sin g, 0),
// c = (c cos be, c (cos al - cos be cos g) / sin g, sqrt(c^2 - cx^2 - cy^2)) (cos 90 = 0 and
// sin 90 = 1 exactly), and from() gives std::nullopt (no periodic geometry: plain distances) when
// a length is not finite or not > 0, an angle is not in (0, 180) or c_z is not a positive number.
class PeriodicCell {
public:
    static std::optional<PeriodicCell> from(const Box& box);
    // usableBox first (std::nullopt without a usable box)
    static std::optional<PeriodicCell> from(const std::optional<Box>& box);
    static std::optional<PeriodicCell> from(const std::optional<std::array<float, 6>>& dims);

    // The minimum image of the difference vector d, in float64: orthogonal cells per periodic
    // axis v - L * nearbyint(v / L) (half to even), repeated while |v| > L / 2 (so a vector within
    // half a cell is returned unchanged, bit for bit); triclinic cells reduce along c (while
    // |v_z| > c_z / 2), then b (|v_y| > b_y / 2), then a (|v_x| > a_x / 2) the same way, and then
    // take the shortest of the 27 neighbouring images v + i a + j b + k c (i, j, k in -1, 0, 1;
    // x = ((v_x + i a_x) + j b_x) + k c_x, y = (v_y + j b_y) + k c_y, z = v_z + k c_z; squared
    // length (x x + y y) + z z; v itself unless another is strictly shorter, a outermost).
    Vec3 minimumImage(const Vec3& d) const;

    bool orthogonal() const { return m_orthogonal; }
    const std::array<bool, 3>& periodic() const { return m_periodic; }
    const std::array<Vec3, 3>& vectors() const { return m_vectors; }
    // the distances between opposite faces (orthogonal: the lengths)
    std::array<double, 3> widths() const;
    // fractional coordinates of a point along a, b, c (orthogonal: 0 along a non-periodic axis)
    std::array<double, 3> fractions(const Vec3& x) const;

private:
    bool m_orthogonal = false;
    std::array<bool, 3> m_periodic{};
    std::array<Vec3, 3> m_vectors{};
};

// The periodic search (MDBridge's own; replaces MDATK's call of MDAnalysis' capped_distance, whose
// behaviour docs/specs/m3-capped-distance.md records): every pair (i, j) whose minimum-image
// distance is within the cutoff, once, sorted by (i, j) like the open path. The difference is
// a[i] - b[j] in float64 (the float32 points widened), PeriodicCell::minimumImage of it, the
// squared distance s = (x x + y y) + z z; a pair is kept iff s <= cutoff * cutoff (the open
// path's test: a negative cutoff acts as its absolute value, an infinite one keeps every pair, a
// NaN one none) and its distance is sqrt(s). A cell without periodic geometry (PeriodicCell::from
// gives std::nullopt) gives plain distances. A point with a NaN or infinite coordinate is in no
// pair. The result does not depend on the number of points or the cutoff's size relative to the
// cell (a cell list only prunes candidates), and nothing throws for any cutoff or cell.
PairList pairsWithinPeriodic(std::span<const float> a, std::span<const float> b, double cutoff,
                             const Box& box);
PairList pairsWithinPeriodic(std::span<const double> a, std::span<const double> b, double cutoff,
                             const Box& box);

// minimum_image_vectors(vectors, box): with a usable box and PeriodicCell geometry, each vector's
// PeriodicCell::minimumImage (orthogonal and triclinic cells: the same image the pair search
// uses); otherwise (no box, unusable, no periodic geometry) the vectors unchanged.
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
