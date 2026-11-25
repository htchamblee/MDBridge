#pragma once
// The projected density map ("shadow") of the structure view, as MDATK computes it
// (mdatk/density.py: project_density, _above_level, shadow_rgba, anchored): the map above a
// contour level integrated along the viewing direction, normalised to its peak, and drawn in one
// colour with alpha ~ density. Qt-free; the structure panel (W2-B) places and draws the image.
//
// Conventions (MDATK's): a placement is a RigidTransform taking world (superposed model)
// coordinates to map coordinates; view coordinates are (x_world - centre) @ view^T (x right, y up,
// z toward the viewer); the image's rows run bottom -> top (v), its columns left -> right (u), and
// the extent (u0, u1, v0, v1) is in view coordinates (Å) at the pixels' outer edges.
//
// Improvements over MDATK v0.0.4.1 (docs/DIFFERENCES.md):
// - RN-30: non-finite voxels (NaN, +-inf) carry no density; a non-finite level gives no shadow
//   (MDATK: one NaN voxel makes the whole image NaN; an infinite one makes it NaN or blank).
// - RN-31: a non-finite placement, view, centre or map grid raises ValueError("density map
//   placement, view or grid is not finite") (MDATK: "cannot convert float NaN to integer" or an
//   OverflowError from int(ceil(...))).
// - RN-32: anchored() with an anchor that is not n x 3 numbers or with non-finite anchor /
//   reference points uses the stored placement and says why (MDATK: numpy's ValueError or
//   LinAlgError "SVD did not converge", and the map is not drawn).

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/io/mrc.hpp"

namespace mdbridge::render::molecule {

struct Shadow {
    std::size_t rows = 0;  // v, bottom -> top
    std::size_t cols = 0;  // u
    std::vector<float> image;  // rows x cols, row-major, values in [0, 1] (peak 1)
    std::array<double, 4> extent{};  // u0, u1, v0, v1 (Å, view coordinates)
};

// _above_level(map, level): max(map - float32(level), 0) as float32, in the map's [z][y][x] order;
// non-finite differences are 0 (RN-30).
std::vector<float> aboveLevel(const io::DensityMap& map, double level);

// The same, cached for the two most recently used (map, level) pairs (MDATK's _LEVEL_CACHE); the
// cache holds the maps, so a key never names a freed map. Thread-safe.
std::shared_ptr<const std::vector<float>> aboveLevelCached(const std::shared_ptr<const io::DensityMap>& map,
                                                           double level);

// The visible window (u0, u1, v0, v1) in view coordinates; NaN entries are ignored (Python's
// max / min keep the map's own bound).
using ShadowWindow = std::array<double, 4>;

// project_density(map, world_to_map, view, centre, level, window, max_pixels, max_depth): the map
// above `level` integrated along view z over the map box (clipped to the window in u and v), on a
// square grid of step max(smallest voxel, longest side / max_pixels, depth / max_depth), sampled
// trilinearly (scipy affine_transform, order 1, mode "constant") and summed in float32. nullopt when
// the box misses the window, when nothing is above the level, or when the level is not finite
// (RN-30). Throws core::ValueError for non-finite inputs (RN-31) and std::invalid_argument for
// max_pixels / max_depth < 1. The shared_ptr overload uses aboveLevelCached.
std::optional<Shadow> projectDensity(const io::DensityMap& map, const core::RigidTransform& worldToMap,
                                     const core::Mat3& view, const core::Vec3& centre, double level,
                                     const std::optional<ShadowWindow>& window = std::nullopt, int maxPixels = 220,
                                     int maxDepth = 220);
std::optional<Shadow> projectDensity(const std::shared_ptr<const io::DensityMap>& map,
                                     const core::RigidTransform& worldToMap, const core::Mat3& view,
                                     const core::Vec3& centre, double level,
                                     const std::optional<ShadowWindow>& window = std::nullopt, int maxPixels = 220,
                                     int maxDepth = 220);

// shadow_rgba(shadow, rgb, opacity, gamma): float32 RGBA (rows x cols x 4, rows bottom -> top):
// the colour (as float32) everywhere, alpha clip(image, 0, 1)^gamma * clip(opacity, 0, 1).
std::vector<float> shadowRgba(const Shadow& shadow, const std::array<double, 3>& rgb, double opacity,
                              double gamma = 0.7);

// anchored(values, anchor, reference): the stored placement (12 numbers, RigidTransform::fromList)
// expressed for the current model frame: when the anchor (the reference frame's points when the map
// was placed) and the current reference points have the same n x 3 shape, n >= 3, and differ
// (np.allclose, atol 1e-3), the reference is superposed onto the anchor first (kabsch) and the
// placement composed with that superposition. Throws core::ValueError for bad values (fromList's
// messages). `note` is empty, or says why an unusable anchor was ignored (RN-32), for a warning.
struct AnchoredPlacement {
    core::RigidTransform transform;
    std::string note;
};
AnchoredPlacement anchored(std::span<const double> values, std::span<const std::vector<double>> anchor,
                           std::span<const core::Vec3> reference);

}  // namespace mdbridge::render::molecule
