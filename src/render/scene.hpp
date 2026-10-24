#pragma once
// The placeholder scene of milestone 1, shared by the Vulkan and the software view: a triangle in
// the residue-group colours turning on a dark background. The Vulkan view's shader
// (shaders/placeholder.vert) holds the same corners and colours.

#include <QColor>
#include <QPointF>

namespace mdbridge::render::scene {

inline constexpr float kBackground[3] = {0.102f, 0.118f, 0.153f};  // #1a1e27
inline constexpr double kTurnsPerSecond = 0.1;

inline constexpr QPointF kCorners[3] = {{0.0, -0.62}, {0.54, 0.31}, {-0.54, 0.31}};  // clip space, y down
inline const QColor kColours[3] = {QColor(0x2a, 0x78, 0xd6), QColor(0xeb, 0x68, 0x34), QColor(0x4c, 0xaf, 0x50)};

// Rotation after `milliseconds` of animation.
inline double angleAt(qint64 milliseconds) {
    return 2.0 * 3.14159265358979323846 * kTurnsPerSecond * static_cast<double>(milliseconds) / 1000.0;
}

}  // namespace mdbridge::render::scene
