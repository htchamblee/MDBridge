#pragma once
// The built-in 3D renderer (MDATK's mdatk/render/molecule/software.py): a CPU z-buffer rasteriser
// with no GPU and no Qt, so it works the same on a laptop, in the GUI and on a headless cluster
// node, and it is the fallback of the Vulkan renderer.
//
// Primitives are rasterised exactly in an orthographic view matched to the panel's axes:
//   * triangles by scanline spans with plane-linear depth and interpolated normals,
//   * spheres as per-pixel impostors (sub-pixel phase in half pixels),
//   * sticks by analytic ray / cylinder intersection.
// Every fragment carries a 64-bit key (quantised depth | fragment order), so the maximum per pixel
// keeps the front-most fragment (ties: the later one, as MDATK's np.maximum.at on fragment numbers);
// only those pixels are shaded (ambient + Lambert + Phong highlight from a camera light). The
// compositor then applies the depth cue (fog toward the panel colour), the glow of residues near the
// target (a tint plus a soft halo at a quarter of the resolution), ink outlines at depth jumps and
// silhouettes, and alpha anti-aliasing of the edges, and returns a straight-alpha RGBA image (row 0
// at the top). Representations with opacity < 1 are drawn in a second pass: the front-most
// translucent layer in front of the opaque surface is blended over it once. Spheres hidden behind
// others are culled per screen tile (the image is the same with and without culling).
//
// Arithmetic follows MDATK's numpy code operation by operation in its dtypes (float32 where numpy
// works in float32, double elsewhere; built without floating-point contraction), so the images are
// MDATK's to the last bit except where a product of 3 x 3 matrices or a libm function rounds
// differently (docs/specs/m4-plan.md section 3: >= 99.5 % of pixels equal, all within 2/255).
// Rendering runs on several threads (bands of rows); the image does not depend on their number.
//
// Improvements over MDATK v0.0.4.1 (docs/DIFFERENCES.md, "Rendering: software 3D renderer"):
// - RN-26: the edge anti-aliasing sees the real neighbours of its region (MDATK filtered the
//   molecule's bounding box with zeros outside it, which drew a darker rectangle into a glow halo).
// - RN-27: "fast" quality at an odd width or height covers the whole image at the requested scale
//   (MDATK rounded the half size down: a transparent last column / row, the view 1 px too small).
// - RN-28: a view with empty, reversed or non-finite limits, or a non-finite rotation, depth range
//   or size below 1 px, throws std::invalid_argument("3D view ...") (MDATK: ZeroDivisionError or
//   garbage); primitives whose screen position is not finite or is beyond 2^40 px are skipped like
//   non-finite atoms; non-finite glow counts as none; a sphere far larger than the image is
//   rasterised only where it is on the image (MDATK: MemoryError for its sprite).
// - RN-29: the glow halo's normalisation is computed for the exact halo size (MDATK cached it under
//   the size rounded to 0.01 px, so a halo's strength depended on the renders before it).

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/topology.hpp"

namespace mdbridge::render::molecule {

// Orthographic view matching the axes: view = xyz @ rotation^T (x right, y up, z toward the
// viewer), the axes' data limits (Å) and their size in pixels, and the panel's depth range (fog).
struct Viewport {
    std::array<Vec3, 3> rotation{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};  // rows
    std::array<double, 2> xlim{-1.0, 1.0};
    std::array<double, 2> ylim{-1.0, 1.0};
    int width = 2;
    int height = 2;
    std::array<double, 2> zrange{-1.0, 1.0};
};

enum class RenderQuality {
    High,  // full resolution with outlines and anti-aliasing
    Fast,  // half resolution, no outlines or anti-aliasing (GUI drags), scaled back to the size
};

// Straight-alpha RGBA, row 0 at the top, rows of width * 4 bytes.
struct RgbaImage {
    std::size_t width = 0;
    std::size_t height = 0;
    std::vector<std::uint8_t> rgba;

    std::uint8_t* pixel(std::size_t x, std::size_t y) { return rgba.data() + 4 * (y * width + x); }
    const std::uint8_t* pixel(std::size_t x, std::size_t y) const { return rgba.data() + 4 * (y * width + x); }
};

struct SoftwareRenderOptions {
    unsigned threads = 0;        // worker threads; 0 = the hardware's (at most 16)
    std::size_t cull_min = 64;   // spheres per list from which hidden ones are culled
};

// The front-most fragment of every covered pixel (pixels in ascending order, index y * width + x):
// lit colour (clipped to 0..1), view depth (Å), owner residue row (-1 none) and opacity.
struct ResolvedPixels {
    std::size_t width = 0;
    std::size_t height = 0;
    double scale = 0.0;  // pixels per Å
    std::vector<std::int64_t> pixel;
    std::vector<std::array<float, 3>> rgb;
    std::vector<float> depth;
    std::vector<std::int64_t> owner;
    std::vector<float> alpha;
};

// Stateless apart from its options: one renderer may serve several figures and threads.
class SoftwareRenderer {
public:
    explicit SoftwareRenderer(SoftwareRenderOptions options = {}) : options_(options) {}

    // The items as an image of the viewport's size. Throws std::invalid_argument for an unusable
    // viewport (RN-28) or inconsistent items (colour or owner counts, triangle indices).
    RgbaImage render(const std::vector<Colored>& items, const Viewport& viewport, const Shading& shading,
                     RenderQuality quality = RenderQuality::High) const;

    // The z-buffer of all items (whatever their opacity) at full resolution, resolved: MDATK's
    // _ZPass.resolve() (tests, other backends).
    ResolvedPixels resolve(const std::vector<Colored>& items, const Viewport& viewport) const;

    // Spheres of `spheres` the occlusion culling skips in this view (depth range from `items`).
    std::size_t culledSpheres(const Spheres& spheres, const std::vector<Colored>& items,
                              const Viewport& viewport) const;

    const SoftwareRenderOptions& options() const { return options_; }

private:
    SoftwareRenderOptions options_;
};

// Items with glowing residues tinted toward the glow colour (for renderers that draw the molecule
// themselves; MDATK's tinted).
std::vector<Colored> tinted(const std::vector<Colored>& items, const Shading& shading);

// The glow halo alone (straight-alpha RGBA), nullopt when nothing glows (MDATK's glow_underlay).
std::optional<RgbaImage> glowUnderlay(const std::vector<Colored>& items, const Viewport& viewport,
                                      const Shading& shading, const SoftwareRenderOptions& options = {});

// Straight-alpha `top` over `under` (same size; MDATK's over). Throws std::invalid_argument for
// different sizes.
RgbaImage over(const RgbaImage& top, const RgbaImage& under);

}  // namespace mdbridge::render::molecule
