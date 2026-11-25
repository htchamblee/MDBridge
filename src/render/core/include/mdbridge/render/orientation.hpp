#pragma once
// Model orientation of the structure view (MDATK's mdatk/render/orientation.py): three angles in
// degrees about X, Y and Z composed as R = Rz Ry Rx and applied about the reference frame's
// centroid; the same matrix for every frame of the movie. Qt-free.
//
// Differences from MDATK (docs/DIFFERENCES.md): wrapDegrees folds only -180 itself (within 1e-9)
// to +180, where MDATK's np.isclose folded everything within 0.0018 degrees of it (RN-22).

#include <array>
#include <span>
#include <string>
#include <vector>

namespace mdbridge::render {

using Vec3d = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;  // rows

Mat3 rotationX(double degrees);
Mat3 rotationY(double degrees);
Mat3 rotationZ(double degrees);
// a @ b (each element summed k = 0, 1, 2 in order)
Mat3 multiply(const Mat3& a, const Mat3& b);
Mat3 transpose(const Mat3& m);

// euler_to_matrix(angles) = Rz(az) Ry(ay) Rx(ax).
Mat3 eulerToMatrix(const Vec3d& angles);
// The same for visual.rotation; throws core::ValueError("expected 3 angles, got N") for another
// count (MDATK: Python's unpacking error).
Mat3 eulerToMatrix(std::span<const double> angles);

// wrap_degrees: into (-180, 180] ((a + 180) mod 360 - 180 with Python's float modulo; -180 -> 180).
double wrapDegrees(double degrees);
// matrix_to_euler: the inverse of eulerToMatrix, each angle in (-180, 180]; at gimbal lock
// (|cos(beta)| <= 1e-8) Z is folded into X.
Vec3d matrixToEuler(const Mat3& r);
// drag(angles, dx, dy): new angles after a mouse drag of (dx, dy) pixels (screen-space trackball:
// Ry(dx * dpp) Rx(dy * dpp) R).
Vec3d drag(const Vec3d& angles, double dx, double dy, double degreesPerPixel = 0.5);
// project(points, angles, center): (p - center) @ R^T: view coordinates (x right, y up, z to the
// viewer).
std::vector<Vec3d> project(std::span<const Vec3d> points, const Vec3d& angles, const Vec3d& center);
// format_angles: "X/Y/Z = +10.0°, -20.0°, +0.0°  (+0.1745, -0.3491, +0.0000 rad)".
std::string formatAngles(const Vec3d& angles);

}  // namespace mdbridge::render
