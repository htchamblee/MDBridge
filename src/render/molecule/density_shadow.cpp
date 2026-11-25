// The density shadow: a port of MDATK's project_density / _above_level / shadow_rgba / anchored
// (mdatk/density.py, the owner's code), with the resampling of scipy.ndimage.affine_transform
// (order 1, mode "constant") from ndimage.cpp. Built without floating-point contraction, so that the
// index arithmetic rounds as numpy's does; the float32 sums and the normalisation follow numpy's
// float32 arithmetic (vol.sum(axis=0) adds the depth slices in order).

#include "mdbridge/render/molecule/density_shadow.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <tuple>

#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/render/molecule/ndimage.hpp"

namespace mdbridge::render::molecule {

namespace {

constexpr const char* kNotFinite = "density map placement, view or grid is not finite";

// Python's max(a, b) / min(a, b): the first argument unless the second compares greater / less
// (a NaN second argument is ignored, a NaN first one kept).
double pyMax(double a, double b) { return b > a ? b : a; }
double pyMin(double a, double b) { return b < a ? b : a; }

// np.min / np.max of a few doubles: NaN propagates
double npMin(std::span<const double> v) {
    double m = v[0];
    for (double x : v) {
        if (std::isnan(x)) {
            return x;
        }
        m = x < m ? x : m;
    }
    return m;
}
double npMax(std::span<const double> v) {
    double m = v[0];
    for (double x : v) {
        if (std::isnan(x)) {
            return x;
        }
        m = x > m ? x : m;
    }
    return m;
}

// np.float32(level): out-of-range values round to +-inf (as numpy casts them)
float toFloat32(double x) {
    if (x > static_cast<double>(FLT_MAX) || x < -static_cast<double>(FLT_MAX)) {
        // beyond the largest float32: inf unless round-to-nearest still gives FLT_MAX
        const double limit = static_cast<double>(FLT_MAX) + std::ldexp(1.0, 103);  // half an ulp above FLT_MAX
        if (std::isnan(x)) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        if (std::fabs(x) >= limit) {
            return x > 0 ? std::numeric_limits<float>::infinity() : -std::numeric_limits<float>::infinity();
        }
        return x > 0 ? FLT_MAX : -FLT_MAX;
    }
    return static_cast<float>(x);
}

std::optional<Shadow> project(const io::DensityMap& map, std::span<const float> data,
                              const core::RigidTransform& toMap, const core::Mat3& view, const core::Vec3& centre,
                              const std::optional<ShadowWindow>& window, int maxPixels, int maxDepth) {
    // the map box's corners in view coordinates
    const core::Vec3 lo = map.origin;
    const core::Vec3 hi = map.extent();
    const core::RigidTransform toWorld = toMap.inverse();
    std::array<double, 8> cu{}, cv{}, cw{};
    int n = 0;
    for (double x : {lo[0], hi[0]}) {
        for (double y : {lo[1], hi[1]}) {
            for (double z : {lo[2], hi[2]}) {
                const core::Vec3 p = toWorld.apply(core::Vec3{x, y, z});
                const double d[3] = {p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]};
                double out[3];
                for (int j = 0; j < 3; ++j) {
                    out[j] = d[0] * view[j][0] + d[1] * view[j][1] + d[2] * view[j][2];
                }
                cu[n] = out[0];
                cv[n] = out[1];
                cw[n] = out[2];
                ++n;
            }
        }
    }
    double u0 = npMin(cu), v0 = npMin(cv), w0 = npMin(cw);
    double u1 = npMax(cu), v1 = npMax(cv), w1 = npMax(cw);
    if (window) {
        u0 = pyMax(u0, (*window)[0]);
        u1 = pyMin(u1, (*window)[1]);
        v0 = pyMax(v0, (*window)[2]);
        v1 = pyMin(v1, (*window)[3]);
    }
    if (u1 <= u0 || v1 <= v0) {
        return std::nullopt;
    }
    const double minVoxel = npMin(map.voxel);
    const double step = pyMax(pyMax(minVoxel, pyMax(u1 - u0, v1 - v0) / maxPixels), (w1 - w0) / maxDepth);
    const double qu = (u1 - u0) / step, qv = (v1 - v0) / step, qw = (w1 - w0) / step;
    if (!std::isfinite(qu) || !std::isfinite(qv) || !std::isfinite(qw)) {
        throw core::ValueError(kNotFinite);  // RN-31 (MDATK: int(ceil(nan / inf)) raises)
    }
    // bounded by max_pixels / max_depth (+ 1 for rounding), unless the voxel term dominates (then fewer)
    const auto count = [](double q) { return std::max<std::int64_t>(2, static_cast<std::int64_t>(std::ceil(q)) + 1); };
    const std::int64_t nu = count(qu), nv = count(qv), nw = count(qw);

    // view (u, v, w) -> map index (x, y, z): idx = A @ (u, v, w) + b
    core::Mat3 r{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {  // rotation @ view^T
            r[i][j] = toMap.rotation[i][0] * view[j][0] + toMap.rotation[i][1] * view[j][1] +
                      toMap.rotation[i][2] * view[j][2];
        }
    }
    const core::Vec3 off = toMap.apply(centre);
    core::Mat3 aIdx{};
    core::Vec3 bIdx{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            aIdx[i][j] = r[i][j] / map.voxel[i];
        }
        bIdx[i] = (off[i] - map.origin[i]) / map.voxel[i];
    }
    // output axes (w, v, u) with spacing `step`; the map array's axes are (z, y, x)
    Matrix3 m{};
    for (int outAx = 0; outAx < 3; ++outAx) {
        for (int inAx = 0; inAx < 3; ++inAx) {
            m[inAx][outAx] = aIdx[2 - inAx][2 - outAx] * step;
        }
    }
    const double viewStart[3] = {u0, v0, w0};
    core::Vec3 idx0{};
    for (int i = 0; i < 3; ++i) {
        idx0[i] = aIdx[i][0] * viewStart[0] + aIdx[i][1] * viewStart[1] + aIdx[i][2] * viewStart[2] + bIdx[i];
    }
    const std::array<double, 3> offset = {idx0[2], idx0[1], idx0[0]};
    const Shape3 shape = {static_cast<std::size_t>(map.shapeXyz[2]), static_cast<std::size_t>(map.shapeXyz[1]),
                          static_cast<std::size_t>(map.shapeXyz[0])};

    // vol = affine_transform(...) as float32, img = vol.sum(axis=0): depth slices added in order
    Shadow out;
    out.rows = static_cast<std::size_t>(nv);
    out.cols = static_cast<std::size_t>(nu);
    out.image.assign(out.rows * out.cols, 0.0f);
    for (std::size_t v = 0; v < out.rows; ++v) {
        for (std::size_t u = 0; u < out.cols; ++u) {
            float acc = static_cast<float>(affineSampleLinearConstant(data, shape, m, offset, 0, v, u));
            for (std::size_t w = 1; w < static_cast<std::size_t>(nw); ++w) {
                acc += static_cast<float>(affineSampleLinearConstant(data, shape, m, offset, w, v, u));
            }
            out.image[v * out.cols + u] = acc;
        }
    }
    // peak = float(img.max()) (NaN propagates); img / peak in float32
    float peak = out.image[0];
    for (float x : out.image) {
        if (std::isnan(x)) {
            peak = x;
            break;
        }
        peak = x > peak ? x : peak;
    }
    if (peak <= 0.0f) {
        return std::nullopt;
    }
    for (float& x : out.image) {
        x = x / peak;
    }
    out.extent = {u0 - 0.5 * step, u0 + (static_cast<double>(nu) - 0.5) * step, v0 - 0.5 * step,
                  v0 + (static_cast<double>(nv) - 0.5) * step};
    return out;
}

void checkLimits(int maxPixels, int maxDepth) {
    if (maxPixels < 1 || maxDepth < 1) {
        throw std::invalid_argument("projectDensity: max_pixels and max_depth must be positive");
    }
}

}  // namespace

std::vector<float> aboveLevel(const io::DensityMap& map, double level) {
    const float l = toFloat32(level);
    std::vector<float> out(map.data.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        const float d = map.data[i] - l;
        // np.maximum(d, 0) keeps d when d >= 0; non-finite differences carry no density (RN-30)
        out[i] = std::isfinite(d) ? (d >= 0.0f ? d : 0.0f) : 0.0f;
    }
    return out;
}

std::shared_ptr<const std::vector<float>> aboveLevelCached(const std::shared_ptr<const io::DensityMap>& map,
                                                           double level) {
    using Entry = std::tuple<std::shared_ptr<const io::DensityMap>, double, std::shared_ptr<const std::vector<float>>>;
    static std::mutex mutex;
    static std::deque<Entry> cache;
    {
        const std::lock_guard<std::mutex> lock(mutex);
        for (const auto& [m, l, values] : cache) {
            if (m == map && l == level) {
                return values;
            }
        }
    }
    auto values = std::make_shared<const std::vector<float>>(aboveLevel(*map, level));
    const std::lock_guard<std::mutex> lock(mutex);
    if (cache.size() >= 2) {
        cache.pop_front();
    }
    cache.emplace_back(map, level, values);
    return values;
}

std::optional<Shadow> projectDensity(const io::DensityMap& map, const core::RigidTransform& worldToMap,
                                     const core::Mat3& view, const core::Vec3& centre, double level,
                                     const std::optional<ShadowWindow>& window, int maxPixels, int maxDepth) {
    checkLimits(maxPixels, maxDepth);
    if (!std::isfinite(level)) {
        return std::nullopt;  // RN-30: nothing is above a NaN level (and everything above -inf is infinite)
    }
    const std::vector<float> data = aboveLevel(map, level);
    return project(map, data, worldToMap, view, centre, window, maxPixels, maxDepth);
}

std::optional<Shadow> projectDensity(const std::shared_ptr<const io::DensityMap>& map,
                                     const core::RigidTransform& worldToMap, const core::Mat3& view,
                                     const core::Vec3& centre, double level,
                                     const std::optional<ShadowWindow>& window, int maxPixels, int maxDepth) {
    if (!map) {
        throw std::invalid_argument("projectDensity: no map");
    }
    checkLimits(maxPixels, maxDepth);
    if (!std::isfinite(level)) {
        return std::nullopt;
    }
    const auto data = aboveLevelCached(map, level);
    return project(*map, *data, worldToMap, view, centre, window, maxPixels, maxDepth);
}

std::vector<float> shadowRgba(const Shadow& shadow, const std::array<double, 3>& rgb, double opacity, double gamma) {
    // np.clip(opacity, 0, 1) of a Python float (NaN stays NaN), then float32 arithmetic
    const double op = std::isnan(opacity) ? opacity : std::clamp(opacity, 0.0, 1.0);
    const float opf = static_cast<float>(op);
    const float g = static_cast<float>(gamma);
    const float color[3] = {static_cast<float>(rgb[0]), static_cast<float>(rgb[1]), static_cast<float>(rgb[2])};
    std::vector<float> out(shadow.image.size() * 4);
    for (std::size_t i = 0; i < shadow.image.size(); ++i) {
        const float x = shadow.image[i];
        const float c = std::isnan(x) ? x : std::clamp(x, 0.0f, 1.0f);
        out[4 * i + 0] = color[0];
        out[4 * i + 1] = color[1];
        out[4 * i + 2] = color[2];
        out[4 * i + 3] = std::pow(c, g) * opf;
    }
    return out;
}

AnchoredPlacement anchored(std::span<const double> values, std::span<const std::vector<double>> anchor,
                           std::span<const core::Vec3> reference) {
    AnchoredPlacement out{core::RigidTransform::fromList(values), {}};
    if (anchor.empty()) {
        return out;  // np.asarray([]) has one dimension
    }
    const std::size_t width = anchor[0].size();
    for (const auto& row : anchor) {
        if (row.size() != width) {
            out.note = "density anchor ignored: its points are not all x, y, z";  // RN-32
            return out;
        }
    }
    if (width != 3 || anchor.size() != reference.size() || anchor.size() < 3) {
        return out;
    }
    std::vector<core::Vec3> a(anchor.size());
    for (std::size_t i = 0; i < anchor.size(); ++i) {
        a[i] = {anchor[i][0], anchor[i][1], anchor[i][2]};
    }
    const auto finite = [](std::span<const core::Vec3> pts) {
        return std::all_of(pts.begin(), pts.end(), [](const core::Vec3& p) {
            return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
        });
    };
    if (!finite(a)) {
        out.note = "density anchor ignored: it holds non-finite numbers";  // RN-32
        return out;
    }
    if (!finite(reference)) {
        out.note = "density anchor ignored: the reference frame holds non-finite coordinates";  // RN-32
        return out;
    }
    bool close = true;  // np.allclose(a, ref, atol=1e-3): |a - ref| <= 1e-3 + 1e-5 |ref|
    for (std::size_t i = 0; i < a.size() && close; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (!(std::fabs(a[i][k] - reference[i][k]) <= 1e-3 + 1e-5 * std::fabs(reference[i][k]))) {
                close = false;
                break;
            }
        }
    }
    if (close) {
        return out;
    }
    const analysis::KabschResult k = analysis::kabsch(reference, a);
    core::RigidTransform superpose;
    superpose.rotation = k.rotation;
    for (int i = 0; i < 3; ++i) {
        const double rm = k.rotation[i][0] * k.mobileCenter[0] + k.rotation[i][1] * k.mobileCenter[1] +
                          k.rotation[i][2] * k.mobileCenter[2];
        superpose.translation[i] = k.referenceCenter[i] - rm;
    }
    out.transform = out.transform.compose(superpose);
    return out;
}

}  // namespace mdbridge::render::molecule
