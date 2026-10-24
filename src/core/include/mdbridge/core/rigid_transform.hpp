#pragma once
// Rigid placements (rotation + translation) as MDATK stores and composes them (mdatk/density.py
// RigidTransform and the rotation helpers): density-map placements in project files, fits, and the
// map following the molecule.

#include <array>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::core {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;  // row-major

struct RigidTransform {
    Mat3 rotation{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    Vec3 translation{0, 0, 0};

    static RigidTransform identity() { return {}; }
    // 12 numbers: rotation row-major, then translation; the rotation re-orthonormalised (SVD) as
    // MDATK does. Throws ValueError with MDATK's four messages ("a rigid transform needs 12
    // numbers (3x3 rotation + translation)", "the transform contains non-finite numbers", "the
    // rotation part of the transform is not a rotation", "... is a reflection").
    static RigidTransform fromList(std::span<const double> numbers);
    // to_list(digits): the 12 numbers rounded to `digits` decimals (Python round), -0.0 as MDATK
    // gives it.
    std::vector<double> toList(int digits = 6) const;

    Vec3 apply(const Vec3& x) const;
    std::vector<Vec3> apply(std::span<const Vec3> points) const;
    RigidTransform inverse() const;
    // this ∘ other: first other, then this
    RigidTransform compose(const RigidTransform& other) const;
    // rotation `r` about the point `centre`
    static RigidTransform about(const Mat3& r, const Vec3& centre);
};

// The rotation of a rotation vector (axis * angle in radians), of a unit quaternion (w, x, y, z),
// and the angle of a rotation in degrees, as MDATK computes them.
Mat3 rotvecMatrix(const Vec3& w);
Mat3 quaternionMatrix(const std::array<double, 4>& q);
double rotationAngle(const Mat3& r);

}  // namespace mdbridge::core
