// Port of MDATK's mdatk/render/orientation.py (the owner's code). numpy's radians / degrees are
// x * (pi / 180) and x * (180 / pi); its sin, cos, arcsin and arctan2 are the C library's.

#include "mdbridge/render/orientation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::render {

namespace {

constexpr double kDegToRad = std::numbers::pi / 180.0;
constexpr double kRadToDeg = 180.0 / std::numbers::pi;

// Python's float modulo (result has the sign of the divisor)
double pyMod(double x, double y) {
    double mod = std::fmod(x, y);
    if (mod != 0.0) {
        if ((y < 0.0) != (mod < 0.0)) {
            mod += y;
        }
    } else {
        mod = std::copysign(0.0, y);
    }
    return mod;
}

// Python's format(value, "+.Nf")
std::string signedFixed(double value, int digits) {
    if (std::isnan(value)) {
        return "+nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "+inf" : "-inf";
    }
    char buf[64];
    std::snprintf(buf, sizeof buf, "%+.*f", digits, value);
    return buf;
}

}  // namespace

Mat3 rotationX(double degrees) {
    const double a = degrees * kDegToRad;
    const double c = std::cos(a), s = std::sin(a);
    return {{{1.0, 0.0, 0.0}, {0.0, c, -s}, {0.0, s, c}}};
}

Mat3 rotationY(double degrees) {
    const double a = degrees * kDegToRad;
    const double c = std::cos(a), s = std::sin(a);
    return {{{c, 0.0, s}, {0.0, 1.0, 0.0}, {-s, 0.0, c}}};
}

Mat3 rotationZ(double degrees) {
    const double a = degrees * kDegToRad;
    const double c = std::cos(a), s = std::sin(a);
    return {{{c, -s, 0.0}, {s, c, 0.0}, {0.0, 0.0, 1.0}}};
}

Mat3 multiply(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    return out;
}

Mat3 transpose(const Mat3& m) {
    Mat3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out[i][j] = m[j][i];
        }
    }
    return out;
}

Mat3 eulerToMatrix(const Vec3d& angles) {
    return multiply(multiply(rotationZ(angles[2]), rotationY(angles[1])), rotationX(angles[0]));
}

Mat3 eulerToMatrix(std::span<const double> angles) {
    if (angles.size() != 3) {
        throw core::ValueError("expected 3 angles, got " + std::to_string(angles.size()));
    }
    return eulerToMatrix(Vec3d{angles[0], angles[1], angles[2]});
}

double wrapDegrees(double degrees) {
    const double a = pyMod(degrees + 180.0, 360.0) - 180.0;
    // RN-22: MDATK used np.isclose(a, -180.0) (relative 1e-5: 0.0018 degrees)
    return std::abs(a + 180.0) <= 1e-9 ? 180.0 : a;
}

Vec3d matrixToEuler(const Mat3& r) {
    double sy = -r[2][0];
    sy = std::min(1.0, std::max(-1.0, sy));
    const double beta = std::asin(sy);
    double alpha = 0.0, gamma = 0.0;
    if (std::abs(std::cos(beta)) > 1e-8) {
        alpha = std::atan2(r[2][1], r[2][2]);
        gamma = std::atan2(r[1][0], r[0][0]);
    } else {  // gimbal lock: fold Z into X
        gamma = 0.0;
        alpha = std::atan2(-r[1][2], r[1][1]);
    }
    return {wrapDegrees(alpha * kRadToDeg), wrapDegrees(beta * kRadToDeg), wrapDegrees(gamma * kRadToDeg)};
}

Vec3d drag(const Vec3d& angles, double dx, double dy, double degreesPerPixel) {
    const Mat3 r = eulerToMatrix(angles);
    const Mat3 turned = multiply(multiply(rotationY(dx * degreesPerPixel), rotationX(dy * degreesPerPixel)), r);
    return matrixToEuler(turned);
}

std::vector<Vec3d> project(std::span<const Vec3d> points, const Vec3d& angles, const Vec3d& center) {
    const Mat3 r = eulerToMatrix(angles);
    std::vector<Vec3d> out;
    out.reserve(points.size());
    for (const Vec3d& p : points) {
        const double x = p[0] - center[0], y = p[1] - center[1], z = p[2] - center[2];
        out.push_back({x * r[0][0] + y * r[0][1] + z * r[0][2], x * r[1][0] + y * r[1][1] + z * r[1][2],
                       x * r[2][0] + y * r[2][1] + z * r[2][2]});
    }
    return out;
}

std::string formatAngles(const Vec3d& angles) {
    std::string deg, rad;
    for (std::size_t k = 0; k < 3; ++k) {
        if (k != 0) {
            deg += ", ";
            rad += ", ";
        }
        deg += signedFixed(angles[k], 1) + "°";
        rad += signedFixed(angles[k] * kDegToRad, 4);
    }
    return "X/Y/Z = " + deg + "  (" + rad + " rad)";
}

}  // namespace mdbridge::render
