#pragma once
// Small dense linear algebra where MDATK calls np.linalg: singular value decompositions of 3x3
// matrices (kabsch, rigid transforms) and of k x 3 point sets (ring planes). One-sided Jacobi
// (Hestenes) rotations, accurate to a few ulp; LAPACK (numpy) differs in the last bits and in the
// signs of singular vector pairs, which every MDATK caller tolerates (see each function).

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "mdbridge/core/rigid_transform.hpp"  // Vec3, Mat3

namespace mdbridge::core {

// a = U diag(s) Vt as np.linalg.svd gives it: singular values in descending order, U and Vt
// orthogonal (rows of Vt are the right singular vectors). Each pair (U column k, Vt row k) is
// unique up to a common sign when s is simple; LAPACK may choose the other sign. For a singular
// value 0 (rank < 3) the matching U columns complete an orthonormal basis (the last one is
// u0 x u1).
struct Svd3 {
    Mat3 u{};
    Vec3 s{};
    Mat3 vt{};
};
Svd3 svd3x3(const Mat3& a);

// Singular values (descending) and right singular vectors (rows of vt) of a k x 3 matrix whose rows
// are `rows`: np.linalg.svd(rows)[1:] (U, k x k, is not formed). Rows of vt carry LAPACK's sign
// ambiguity.
struct SvdRows3 {
    Vec3 s{};
    Mat3 vt{};
};
SvdRows3 svdRows3(std::span<const Vec3> rows);

// np.linalg.svd of a stack of (n, k, 3) matrices: rows holds n * k points, matrix m being rows
// [m k, (m + 1) k). k > 0.
std::vector<SvdRows3> svdRows3Batch(std::span<const Vec3> rows, std::size_t k);

// The column norms of a after the Jacobi rotations, in column order (unsorted): the singular
// values as RigidTransform::fromList checks them.
Vec3 jacobiSingularValues(const Mat3& a);

double det3(const Mat3& a);

}  // namespace mdbridge::core
