// The built-in 3D renderer, ported from MDATK's mdatk/render/molecule/software.py (MDATK v0.0.4.1,
// the same owner's code). See software.hpp for the method and the differences from MDATK.
//
// Port notes. numpy decides the pixels, so every expression keeps numpy's operation order and dtype:
// float32 arrays combined with Python floats stay float32 (NEP 50: the Python scalar is cast to
// float32 first), float32 combined with float64 / int64 arrays becomes float64, and casts to integers
// truncate. The library is built without floating-point contraction (CMakeLists.txt), so no fused
// multiply-add changes a rounding. Where MDATK concatenates fragments and keeps the maximum key per
// pixel (np.maximum.at), this port rasterises bands of rows on several threads and keys each fragment
// by its quantised depth and its primitive's rank in MDATK's fragment order (later fragments have
// larger numbers there, and one primitive covers a pixel once), so the front-most fragment and the
// tie-breaking are MDATK's; the winner's attributes are recomputed from its primitive when the
// buffer is resolved. The scipy.ndimage pieces (gaussian_filter, uniform_filter) are ndimage.cpp's
// ports of SciPy's line-buffer arithmetic.

#include "mdbridge/render/molecule/software.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <variant>

#include "mdbridge/render/molecule/ndimage.hpp"

namespace mdbridge::render::molecule {
namespace {

constexpr int kIdBits = 32;  // fragment order bits (the depth gets the 31 above them, as in MDATK)
constexpr float kAmbient = 0.36f, kDiffuse = 0.64f, kSpecular = 0.22f;
constexpr float kInk = 0.12f;  // INK, every channel
// glow: halo radius / disc (Å), halo opacity, tint toward the glow colour
constexpr double kGlowRadius = 2.0, kGlowDisc = 2.5, kGlowAlpha = 0.75, kGlowTint = 0.25;
constexpr double kHuge = 1099511627776.0;  // 2^40 px: beyond it a primitive is skipped (RN-28)
// LIGHT = [-0.35, 0.55, 0.76] / norm, as float32 (the shader's constants)
constexpr float kLight[3] = {-0x1.65ef5cp-2f, 0x1.193c1p-1f, 0x1.849d7ap-1f};

float f32(double v) { return static_cast<float>(v); }

// np.clip for one value (NaN stays NaN)
template <class T>
T clip(T x, T lo, T hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}
// np.maximum / np.minimum (NaN propagates)
template <class T>
T npMax(T a, T b) {
    if (a != a) {
        return a;
    }
    if (b != b) {
        return b;
    }
    return a > b ? a : b;
}
template <class T>
T npMin(T a, T b) {
    if (a != a) {
        return a;
    }
    if (b != b) {
        return b;
    }
    return a < b ? a : b;
}
// a float that is already an integer value, as an index (clamped far outside any image)
std::int64_t toIndex(double v) {
    constexpr double lim = 4611686018427387904.0;  // 2^62
    if (!(v == v)) {
        return 0;
    }
    return static_cast<std::int64_t>(std::clamp(v, -lim, lim));
}
// astype(np.uint8) of a value clipped to 0..255
std::uint8_t toByte(float v) {
    if (!(v >= 0.0f)) {
        return 0;
    }
    return static_cast<std::uint8_t>(static_cast<int>(std::min(v, 255.0f)));
}

using Rgb = std::array<float, 3>;

// ---- threads -----------------------------------------------------------------------------------

unsigned workerCount(const SoftwareRenderOptions& options, std::size_t pixels) {
    if (pixels < 16384) {
        return 1;
    }
    unsigned n = options.threads;
    if (n == 0) {
        n = std::thread::hardware_concurrency();
        n = n == 0 ? 1 : std::min(n, 16u);
    }
    return std::max(n, 1u);
}

std::size_t bandCount(std::int64_t rows, unsigned threads) {
    if (threads <= 1 || rows <= 1) {
        return 1;
    }
    return static_cast<std::size_t>(std::min<std::int64_t>(rows, static_cast<std::int64_t>(threads) * 4));
}

// fn(band, row0, row1) for the bandCount(rows, threads) bands of rows covering 0..rows, spread over
// `threads` threads. Every band writes its own rows, so the result does not depend on the thread
// count.
template <class F>
void forBands(std::int64_t rows, unsigned threads, F&& fn) {
    const std::size_t nb = bandCount(rows, threads);
    const auto n = static_cast<std::int64_t>(nb);
    auto run = [&](std::size_t b) {
        const auto k = static_cast<std::int64_t>(b);
        fn(b, rows * k / n, rows * (k + 1) / n);
    };
    if (nb == 1) {
        run(0);
        return;
    }
    std::atomic<std::size_t> next{0};
    std::exception_ptr error;
    std::mutex error_mutex;
    auto work = [&] {
        for (;;) {
            const std::size_t b = next.fetch_add(1);
            if (b >= nb) {
                return;
            }
            try {
                run(b);
            } catch (...) {
                const std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) {
                    error = std::current_exception();
                }
            }
        }
    };
    const auto workers = static_cast<unsigned>(std::min<std::size_t>(threads, nb));
    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    for (unsigned t = 1; t < workers; ++t) {
        pool.emplace_back(work);
    }
    work();
    for (std::thread& t : pool) {
        t.join();
    }
    if (error) {
        std::rethrow_exception(error);
    }
}

// ---- camera ------------------------------------------------------------------------------------

// Pixel mapping of the panel's axes (half the size for "fast") and the depth range of the drawn
// geometry (keys) and of the panel (fog).
struct Cam {
    std::array<Vec3, 3> R{};
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    std::int64_t W = 0, H = 0;
    double s = 0, sy = 0;  // pixels per Å
    double zmin = 0, zmax = 0;
    double fog_lo = 0, fog_hi = 0;
    // float32 constants of the depth keys
    float zmin_f = 0, zspan_f = 0;
};

bool finite3(const Vec3& v) { return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]); }

void checkViewport(const Viewport& vp) {
    bool ok = std::isfinite(vp.xlim[0]) && std::isfinite(vp.xlim[1]) && std::isfinite(vp.ylim[0]) &&
              std::isfinite(vp.ylim[1]) && vp.xlim[1] > vp.xlim[0] && vp.ylim[1] > vp.ylim[0];
    if (!ok) {
        throw std::invalid_argument("3D view limits are empty, reversed or not finite");
    }
    for (const Vec3& row : vp.rotation) {
        ok = ok && finite3(row);
    }
    if (!ok || !std::isfinite(vp.zrange[0]) || !std::isfinite(vp.zrange[1])) {
        throw std::invalid_argument("3D view rotation or depth range is not finite");
    }
    if (vp.width < 1 || vp.height < 1) {
        throw std::invalid_argument("3D view size must be at least 1 x 1 pixels");
    }
}

Cam makeCam(const Viewport& vp, const std::vector<Colored>& items, bool fast) {
    checkViewport(vp);
    Cam c;
    c.R = vp.rotation;
    c.xlo = vp.xlim[0];
    c.xhi = vp.xlim[1];
    c.ylo = vp.ylim[0];
    c.yhi = vp.ylim[1];
    if (!fast) {
        c.W = std::max<std::int64_t>(2, vp.width);
        c.H = std::max<std::int64_t>(2, vp.height);
    } else {
        // RN-27: half of the size rounded up, the limits widened to keep two pixels per half pixel
        // (MDATK: int(round(size * 0.5)), which loses the last column / row of odd sizes)
        c.W = std::max<std::int64_t>(2, (vp.width + 1) / 2);
        c.H = std::max<std::int64_t>(2, (vp.height + 1) / 2);
        if (2 * c.W != vp.width) {
            const double sx = static_cast<double>(vp.width) / (vp.xlim[1] - vp.xlim[0]);
            c.xhi = vp.xlim[0] + static_cast<double>(2 * c.W) / sx;
        }
        if (2 * c.H != vp.height) {
            const double sy = static_cast<double>(vp.height) / (vp.ylim[1] - vp.ylim[0]);
            c.ylo = vp.ylim[1] - static_cast<double>(2 * c.H) / sy;
        }
    }
    c.s = static_cast<double>(c.W) / (c.xhi - c.xlo);
    c.sy = static_cast<double>(c.H) / (c.yhi - c.ylo);
    std::tie(c.zmin, c.zmax) = depthRange(items, c.R[2]);
    c.fog_lo = vp.zrange[0];
    c.fog_hi = vp.zrange[1];
    c.zmin_f = f32(c.zmin);
    c.zspan_f = f32(c.zmax - c.zmin);
    return c;
}

struct ScreenPoint {
    float x, y, z;
};

// (xyz @ R^T): pixel x, pixel y (row, downward), view z (Å), each rounded to float32
ScreenPoint screen(const Cam& c, const Vec3& p) {
    const double v0 = p[0] * c.R[0][0] + p[1] * c.R[0][1] + p[2] * c.R[0][2];
    const double v1 = p[0] * c.R[1][0] + p[1] * c.R[1][1] + p[2] * c.R[1][2];
    const double v2 = p[0] * c.R[2][0] + p[1] * c.R[2][1] + p[2] * c.R[2][2];
    return {f32((v0 - c.xlo) * c.s), f32((c.yhi - v1) * c.sy), f32(v2)};
}

// RN-28: a screen position the renderer can place (finite in float32, within 2^40 px)
bool placeable(const ScreenPoint& p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) && std::fabs(p.x) < kHuge &&
           std::fabs(p.y) < kHuge;
}

Vec3 rotate(const Cam& c, const Vec3& n) {
    return {n[0] * c.R[0][0] + n[1] * c.R[0][1] + n[2] * c.R[0][2],
            n[0] * c.R[1][0] + n[1] * c.R[1][1] + n[2] * c.R[1][2],
            n[0] * c.R[2][0] + n[1] * c.R[2][1] + n[2] * c.R[2][2]};
}

Rgb rgb32(const Vec3& v) { return {f32(v[0]), f32(v[1]), f32(v[2])}; }

// ---- spans -------------------------------------------------------------------------------------

// One row of a convex region given as half-planes a*x + b*y + c >= 0, tested at pixel centres
// (_spans): per half-plane its kind (+1: a > 1e-9, -1: a < -1e-9, 0: horizontal) and b / a, c / a
// (horizontal ones: b, c). Columns x0..x1 (empty when x1 < x0).
template <class T, std::size_t K>
struct HalfPlanes {
    std::array<std::int8_t, K> kind{};
    std::array<T, K> b{}, c{};
};

template <class T, std::size_t K>
HalfPlanes<T, K> halfPlanes(const std::array<std::array<T, 3>, K>& G) {
    HalfPlanes<T, K> h;
    for (std::size_t k = 0; k < K; ++k) {
        const T a = G[k][0];
        const bool pos = a > static_cast<T>(1e-9);
        const bool neg = a < static_cast<T>(-1e-9);
        if (pos || neg) {
            const T ia = static_cast<T>(1) / a;
            h.kind[k] = pos ? 1 : -1;
            h.b[k] = G[k][1] * ia;
            h.c[k] = G[k][2] * ia;
        } else {
            h.kind[k] = 0;
            h.b[k] = G[k][1];
            h.c[k] = G[k][2];
        }
    }
    return h;
}

template <class T, std::size_t K>
bool spanRow(const HalfPlanes<T, K>& h, T yc, std::int64_t W, std::int64_t& x0, std::int64_t& x1) {
    const T big = static_cast<T>(1e9);
    T lo = -big;
    T hi = big;
    for (std::size_t k = 0; k < K; ++k) {
        if (h.kind[k] == 0) {  // horizontal constraint: the whole row is in or out
            if (h.b[k] * yc + h.c[k] < 0) {
                return false;
            }
            continue;
        }
        const T bound = -(h.b[k] * yc + h.c[k]);
        if (h.kind[k] > 0) {
            lo = npMax(lo, bound);
        } else {
            hi = npMin(hi, bound);
        }
    }
    const T half = static_cast<T>(0.5);
    const T l = npMax(static_cast<T>(std::ceil(lo - half)), static_cast<T>(0));
    const T r = npMin(static_cast<T>(std::floor(hi - half)), static_cast<T>(W - 1));
    if (!(l == l) || !(r == r)) {
        return false;
    }
    x0 = toIndex(static_cast<double>(std::min(l, static_cast<T>(W))));
    x1 = toIndex(static_cast<double>(std::max(r, static_cast<T>(-1))));
    return x0 <= x1;
}

// ---- primitives --------------------------------------------------------------------------------

struct TriSetup {
    std::array<std::int32_t, 3> v{};
    double a0 = 0, b0 = 0, c0 = 0, a1 = 0, b1 = 0, c1 = 0;
    HalfPlanes<double, 3> h;
    std::int64_t by0 = 0, by1 = 0, bx0 = 0, bx1 = 0;
    double za = 0, zb = 0, zc = 0, zlo = 0, zhi = 0;
};

struct MeshPass {
    std::vector<float> nx, ny, nz;  // view-space normals per vertex (float32)
    std::vector<Rgb> col;
    std::vector<std::int64_t> own;
    std::vector<TriSetup> tris;  // rank order
};

struct Sphere {
    std::int64_t cx = 0, cy = 0, n = 0;
    double rp = 0, rp2 = 0, sx = 0, sy = 0;
    float Z = 0, rpf = 0, sf = 0;
    std::size_t index = 0;  // into col / own
};

struct SpherePass {
    std::vector<Rgb> col;
    std::vector<std::int64_t> own;
    std::vector<Sphere> spheres;  // rank order
};

struct Stick {
    HalfPlanes<float, 4> h;
    std::int64_t by0 = 0, by1 = 0;
    float Ax = 0, Ay = 0, Az = 0, ux = 0, uy = 0, uz = 0, a = 0, R = 0, L = 0;
    std::size_t index = 0;
};

struct StickPass {
    std::vector<Rgb> col;
    std::vector<std::int64_t> own;
    std::vector<Stick> sticks;  // rank order
    float H_half = 0;           // float32(H - 0.5)
    float s = 0;
};

struct Pass {
    std::variant<MeshPass, SpherePass, StickPass> data;
    std::uint64_t base = 0;  // rank of its first primitive in the buffer's fragment order
    std::size_t count = 0;
    float opacity = 1.0f;
};

// Depth key of a fragment: quantised depth (float32, as numpy) above the primitive's rank.
std::int64_t keyOf(const Cam& c, float z, std::uint64_t id) {
    float t = clip((z - c.zmin_f) / c.zspan_f, 0.0f, 1.0f);
    const float q = t * 2147483648.0f;  // float32((1 << 31) - 1)
    const std::int64_t zq = q == q ? static_cast<std::int64_t>(q) : 0;
    // as numpy's int64 shift: a depth key of 2^31 wraps negative (never wins), as in MDATK
    return static_cast<std::int64_t>((static_cast<std::uint64_t>(zq) << kIdBits) | id);
}

struct MeshFrag {
    float l0, l1, z;
};
MeshFrag meshFrag(const TriSetup& t, std::int64_t x, std::int64_t y) {
    const double fx = static_cast<double>(x) + 0.5;
    const double fy = static_cast<double>(y) + 0.5;
    const float l0 = f32(t.a0 * fx + t.b0 * fy + t.c0);
    const float l1 = f32(t.a1 * fx + t.b1 * fy + t.c1);
    double z = t.za * fx + t.zb * fy + t.zc;
    // pixel centres just outside a sliver's edges extrapolate far: stay within the triangle
    z = npMin(npMax(z, t.zlo), t.zhi);
    return {l0, l1, f32(z)};
}

struct SphereFrag {
    bool in;
    float nx, ny, nz, z;
};
SphereFrag sphereFrag(const Sphere& s, std::int64_t x, std::int64_t y) {
    const double dx = (static_cast<double>(x - s.cx) + 0.5) - s.sx;
    const double dy = (static_cast<double>(y - s.cy) + 0.5) - s.sy;
    const double d2 = dx * dx + dy * dy;
    if (!(d2 < s.rp2)) {
        return {false, 0, 0, 0, 0};
    }
    const float dz = f32(std::sqrt(s.rp2 - d2));
    return {true, f32(dx / s.rp), f32(-dy / s.rp), dz / s.rpf, s.Z + dz / s.sf};
}

struct StickFrag {
    bool in;
    float nx, ny, nz, z;
};
StickFrag stickFrag(const Stick& k, float H_half, float s, std::int64_t x, std::int64_t y) {
    const float p0x = (static_cast<float>(x) + 0.5f) - k.Ax;
    const float p0y = (H_half - static_cast<float>(y)) - k.Ay;
    const float m = p0x * k.ux + p0y * k.uy;
    const float pu = m - k.Az * k.uz;
    const float b = -k.uz * m - k.Az * k.a;
    const float c = p0x * p0x + p0y * p0y + k.Az * k.Az - pu * pu - k.R * k.R;
    const float disc = b * b - k.a * c;
    const float tt = (-b + std::sqrt(npMax(disc, 0.0f))) / k.a;
    const float sp = pu + tt * k.uz;
    if (!(disc >= 0.0f && sp >= 0.0f && sp <= k.L)) {
        return {false, 0, 0, 0, 0};
    }
    return {true, p0x - sp * k.ux, p0y - sp * k.uy, -k.Az + tt - sp * k.uz, tt / s};
}

// ---- building the passes -----------------------------------------------------------------------

void checkSize(std::size_t got, std::size_t want, const char* what) {
    if (got != want) {
        throw std::invalid_argument(std::string("3D item: ") + what + " count does not match");
    }
}

std::vector<Rgb> colours(const std::vector<Vec3>& rgb) {
    std::vector<Rgb> out(rgb.size());
    for (std::size_t i = 0; i < rgb.size(); ++i) {
        out[i] = rgb32(rgb[i]);
    }
    return out;
}

std::vector<std::int64_t> owners(const std::vector<std::int32_t>& owner) {
    return std::vector<std::int64_t>(owner.begin(), owner.end());
}

std::optional<MeshPass> meshPass(const Cam& c, const Mesh& m, const std::vector<Vec3>& rgb) {
    const std::size_t nv = m.verts.size();
    checkSize(m.normals.size(), nv, "normal");
    checkSize(m.owner.size(), nv, "owner");
    checkSize(rgb.size(), nv, "colour");
    MeshPass p;
    p.nx.resize(nv);
    p.ny.resize(nv);
    p.nz.resize(nv);
    std::vector<ScreenPoint> sp(nv);
    std::vector<bool> bad(nv);
    for (std::size_t i = 0; i < nv; ++i) {
        sp[i] = screen(c, m.verts[i]);
        bad[i] = !finite3(m.verts[i]) || !finite3(m.normals[i]) || !placeable(sp[i]);
        const Vec3 n = rotate(c, m.normals[i]);
        p.nx[i] = f32(n[0]);
        p.ny[i] = f32(n[1]);
        p.nz[i] = f32(n[2]);
    }
    const double W1 = static_cast<double>(c.W - 1);
    const double H1 = static_cast<double>(c.H - 1);
    for (const molecule::Tri& t : m.tris) {
        for (std::int32_t v : t) {
            if (v < 0 || static_cast<std::size_t>(v) >= nv) {
                throw std::invalid_argument("3D item: mesh triangle index out of range");
            }
        }
        const auto i0 = static_cast<std::size_t>(t[0]);
        const auto i1 = static_cast<std::size_t>(t[1]);
        const auto i2 = static_cast<std::size_t>(t[2]);
        if (bad[i0] || bad[i1] || bad[i2]) {  // a blown-up frame draws what it can
            continue;
        }
        const float back = -0.05f;  // back faces: all vertex normals facing away
        if (!(p.nz[i0] > back || p.nz[i1] > back || p.nz[i2] > back)) {
            continue;
        }
        const double x0 = sp[i0].x, x1 = sp[i1].x, x2 = sp[i2].x;
        const double y0 = sp[i0].y, y1 = sp[i1].y, y2 = sp[i2].y;
        const double A = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        const double bx0 = clip(std::floor(npMin(npMin(x0, x1), x2) - 0.5), 0.0, W1);
        const double bx1 = clip(std::ceil(npMax(npMax(x0, x1), x2) - 0.5), 0.0, W1);
        const double by0 = clip(std::floor(npMin(npMin(y0, y1), y2) - 0.5), 0.0, H1);
        const double by1 = clip(std::ceil(npMax(npMax(y0, y1), y2) - 0.5), 0.0, H1);
        if (!(std::fabs(A) > 1e-4 && bx1 >= bx0 && by1 >= by0)) {
            continue;
        }
        TriSetup r;
        r.v = t;
        const double inv = 1.0 / A;
        r.a0 = (y1 - y2) * inv;
        r.b0 = (x2 - x1) * inv;
        r.c0 = (x1 * y2 - x2 * y1) * inv;
        r.a1 = (y2 - y0) * inv;
        r.b1 = (x0 - x2) * inv;
        r.c1 = (x2 * y0 - x0 * y2) * inv;
        std::array<std::array<double, 3>, 3> G{{{r.a0, r.b0, r.c0},
                                                {r.a1, r.b1, r.c1},
                                                {-r.a0 - r.a1, -r.b0 - r.b1, 1 - r.c0 - r.c1}}};
        for (auto& g : G) {  // slightly conservative: no cracks
            g[2] += 0.02 * std::hypot(g[0], g[1]);
        }
        r.h = halfPlanes(G);
        r.by0 = toIndex(by0);
        r.by1 = toIndex(by1);
        r.bx0 = toIndex(bx0);
        r.bx1 = toIndex(bx1);
        const double z0 = sp[i0].z, z1 = sp[i1].z, z2 = sp[i2].z;  // depth is plane-linear too
        r.za = r.a0 * (z0 - z2) + r.a1 * (z1 - z2);
        r.zb = r.b0 * (z0 - z2) + r.b1 * (z1 - z2);
        r.zc = z2 + r.c0 * (z0 - z2) + r.c1 * (z1 - z2);
        r.zlo = npMin(npMin(z0, z1), z2);
        r.zhi = npMax(npMax(z0, z1), z2);
        p.tris.push_back(r);
    }
    if (p.tris.empty()) {
        return std::nullopt;
    }
    p.col = colours(rgb);
    p.own = owners(m.owner);
    return p;
}

// Spheres that may be visible (conservative tile occlusion culling, MDATK's _unoccluded): for every
// screen tile, the depth each of its pixels is guaranteed to reach is the largest depth any sphere
// has over the whole tile (spheres covering the tile, shrunk by half a pixel for the sub-pixel
// placement); a sphere whose front point lies behind that depth on every tile it touches cannot
// show. The image is the same with and without culling.
std::vector<bool> unoccluded(const Cam& c, const std::vector<ScreenPoint>& P, const std::vector<float>& rad,
                             std::size_t cull_min) {
    const std::size_t n = P.size();
    if (n < cull_min) {
        return std::vector<bool>(n, true);
    }
    std::vector<double> rp(n);
    double rmax = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n; ++i) {
        rp[i] = static_cast<double>(rad[i]) * c.s;
        rmax = npMax(rmax, rp[i]);
    }
    const std::int64_t T = std::max<std::int64_t>(8, toIndex(std::ceil(rmax / 4)));  // ~10 tiles per sphere
    const std::int64_t tw = (c.W + T - 1) / T;
    const std::int64_t th = (c.H + T - 1) / T;
    const double Td = static_cast<double>(T);
    struct Range {
        std::int64_t tx0, tx1, ty0, ty1;
    };
    std::vector<Range> rg(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double X = P[i].x, Y = P[i].y;
        rg[i] = {toIndex(std::floor((X - rp[i]) / Td)), toIndex(std::floor((X + rp[i]) / Td)),
                 toIndex(std::floor((Y - rp[i]) / Td)), toIndex(std::floor((Y + rp[i]) / Td))};
    }
    std::vector<double> occ(static_cast<std::size_t>(tw * th), -std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < n; ++i) {
        const double X = P[i].x, Y = P[i].y, Z = P[i].z;
        const double ro = npMax(rp[i] - 0.5, 0.0);
        const double ro2 = ro * ro;
        for (std::int64_t ty = std::max<std::int64_t>(rg[i].ty0, 0); ty <= std::min(rg[i].ty1, th - 1); ++ty) {
            const double ey = npMax(std::fabs(static_cast<double>(ty * T) - Y),
                                    std::fabs(static_cast<double>((ty + 1) * T) - Y));  // farthest tile edge
            for (std::int64_t tx = std::max<std::int64_t>(rg[i].tx0, 0); tx <= std::min(rg[i].tx1, tw - 1); ++tx) {
                const double ex = npMax(std::fabs(static_cast<double>(tx * T) - X),
                                        std::fabs(static_cast<double>((tx + 1) * T) - X));
                const double r2 = ex * ex + ey * ey;
                if (r2 < ro2) {
                    double& o = occ[static_cast<std::size_t>(ty * tw + tx)];
                    o = npMax(o, Z + std::sqrt(npMax(ro2 - r2, 0.0)) / c.s);
                }
            }
        }
    }
    std::vector<bool> keep(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double front = static_cast<double>(P[i].z) + static_cast<double>(rad[i]);
        bool hidden = true;
        for (std::int64_t ty = std::max<std::int64_t>(rg[i].ty0, 0); hidden && ty <= std::min(rg[i].ty1, th - 1);
             ++ty) {
            for (std::int64_t tx = std::max<std::int64_t>(rg[i].tx0, 0); tx <= std::min(rg[i].tx1, tw - 1); ++tx) {
                if (!(occ[static_cast<std::size_t>(ty * tw + tx)] > front)) {  // tiles off the image draw nothing
                    hidden = false;
                    break;
                }
            }
        }
        keep[i] = !hidden;
    }
    return keep;
}

struct SphereInput {
    std::vector<ScreenPoint> P;
    std::vector<float> rad;
    std::vector<std::size_t> index;  // into the item's arrays
};

SphereInput sphereInput(const Cam& c, const Spheres& s) {
    SphereInput in;
    for (std::size_t i = 0; i < s.xyz.size(); ++i) {
        if (!finite3(s.xyz[i]) || !std::isfinite(s.radius[i])) {
            continue;
        }
        const ScreenPoint p = screen(c, s.xyz[i]);
        const float r = f32(s.radius[i]);
        if (!placeable(p) || !(std::fabs(static_cast<double>(r) * c.s) < kHuge)) {
            continue;
        }
        in.P.push_back(p);
        in.rad.push_back(r);
        in.index.push_back(i);
    }
    return in;
}

void checkSpheres(const Spheres& s, std::size_t n_rgb) {
    const std::size_t n = s.xyz.size();
    checkSize(s.radius.size(), n, "radius");
    checkSize(s.owner.size(), n, "owner");
    checkSize(n_rgb, n, "colour");
}

std::optional<SpherePass> spherePass(const Cam& c, const Spheres& s, const std::vector<Vec3>& rgb,
                                     std::size_t cull_min) {
    checkSpheres(s, rgb.size());
    const SphereInput in = sphereInput(c, s);
    const std::vector<bool> keep = unoccluded(c, in.P, in.rad, cull_min);
    struct Entry {
        float rad;
        int q;
        std::size_t k;
    };
    std::vector<Entry> order;
    for (std::size_t k = 0; k < in.P.size(); ++k) {
        if (!keep[k]) {
            continue;
        }
        const float X = in.P[k].x, Y = in.P[k].y;
        const double cx = std::floor(static_cast<double>(X));
        const double cy = std::floor(static_cast<double>(Y));
        const int qx = std::min(static_cast<int>((static_cast<double>(X) - cx) * 2), 1);
        const int qy = std::min(static_cast<int>((static_cast<double>(Y) - cy) * 2), 1);
        order.push_back({in.rad[k], qx * 2 + qy, k});
    }
    // MDATK's fragment order: by radius (np.unique), sub-pixel phase, then sphere
    std::sort(order.begin(), order.end(), [](const Entry& a, const Entry& b) {
        if (a.rad != b.rad) {
            return a.rad < b.rad;
        }
        if (a.q != b.q) {
            return a.q < b.q;
        }
        return a.k < b.k;
    });
    SpherePass p;
    for (const Entry& e : order) {
        const ScreenPoint& P = in.P[e.k];
        Sphere sp;
        sp.cx = toIndex(std::floor(static_cast<double>(P.x)));
        sp.cy = toIndex(std::floor(static_cast<double>(P.y)));
        sp.rp = static_cast<double>(e.rad) * c.s;
        sp.rp2 = sp.rp * sp.rp;
        sp.n = toIndex(std::ceil(sp.rp)) + 1;
        if (sp.n < 0) {
            continue;  // a negative radius has no pixels
        }
        sp.sx = e.q < 2 ? 0.25 : 0.75;
        sp.sy = e.q % 2 == 0 ? 0.25 : 0.75;
        sp.Z = P.z;
        sp.rpf = f32(sp.rp);
        sp.sf = f32(c.s);
        sp.index = in.index[e.k];
        p.spheres.push_back(sp);
    }
    if (p.spheres.empty()) {
        return std::nullopt;
    }
    p.col = colours(rgb);
    p.own = owners(s.owner);
    return p;
}

std::optional<StickPass> stickPass(const Cam& c, const Sticks& st, const std::vector<Vec3>& rgb) {
    const std::size_t n = st.a.size();
    checkSize(st.b.size(), n, "stick end");
    checkSize(st.radius.size(), n, "radius");
    checkSize(st.owner.size(), n, "owner");
    checkSize(rgb.size(), n, "colour");
    const float s = f32(c.s);
    const float Hf = static_cast<float>(c.H);
    const float H1 = static_cast<float>(c.H - 1);
    StickPass p;
    p.s = s;
    p.H_half = f32(static_cast<double>(c.H) - 0.5);
    for (std::size_t i = 0; i < n; ++i) {
        if (!finite3(st.a[i]) || !finite3(st.b[i]) || !std::isfinite(st.radius[i])) {
            continue;
        }
        const ScreenPoint A = screen(c, st.a[i]);
        const ScreenPoint B = screen(c, st.b[i]);
        if (!placeable(A) || !placeable(B)) {
            continue;
        }
        const float ax = A.x, ay = A.y, az = A.z, bx = B.x, by = B.y, bz = B.z;
        Stick k;
        k.Ax = ax;  // pixel units, y up
        k.Ay = Hf - ay;
        k.Az = az * s;
        const float Dx = bx - ax, Dy = ay - by, Dz = (bz - az) * s;
        k.L = std::sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
        k.R = f32(st.radius[i]) * s;
        // end-on sticks are covered by their balls
        const bool ok = (k.L > 1e-3f) && (k.L * k.L - Dz * Dz > 1e-4f);
        const float inv = 1.0f / npMax(k.L, 1e-3f);
        k.ux = Dx * inv;
        k.uy = Dy * inv;
        k.uz = Dz * inv;
        k.a = npMax(1.0f - k.uz * k.uz, 1e-6f);
        const float fy0 = clip(std::floor(npMin(ay, by) - k.R), 0.0f, H1);
        const float fy1 = clip(std::ceil(npMax(ay, by) + k.R), 0.0f, H1);
        if (!(ok && fy1 >= fy0)) {
            continue;
        }
        k.by0 = toIndex(static_cast<double>(fy0));
        k.by1 = toIndex(static_cast<double>(fy1));
        const float ddx = bx - ax, ddy = by - ay;
        const float L2 = std::sqrt(ddx * ddx + ddy * ddy);
        const float dx2 = ddx / npMax(L2, 1e-6f), dy2 = ddy / npMax(L2, 1e-6f);
        const float nx2 = -dy2, ny2 = dx2;
        const float Rs = k.R + 0.5f;
        const std::array<std::array<float, 3>, 4> G{{  // screen-space capsule bound
            {nx2, ny2, Rs - nx2 * ax - ny2 * ay},
            {-nx2, -ny2, Rs + nx2 * ax + ny2 * ay},
            {dx2, dy2, Rs - dx2 * ax - dy2 * ay},
            {-dx2, -dy2, L2 + Rs + dx2 * ax + dy2 * ay}}};
        k.h = halfPlanes(G);
        k.index = i;
        p.sticks.push_back(k);
    }
    if (p.sticks.empty()) {
        return std::nullopt;
    }
    p.col = colours(rgb);
    p.own = owners(st.owner);
    return p;
}

// ---- the z-buffer ------------------------------------------------------------------------------

struct Fragments {  // the resolved front-most fragments, pixels ascending
    std::vector<std::int64_t> pix;
    std::vector<Rgb> rgb;  // lit, clipped to 0..1
    std::vector<float> depth;
    std::vector<std::int64_t> owner;
    std::vector<float> alpha;
};

// Ambient + Lambert + Phong highlight (shininess 24) from a light near the camera (MDATK's shade).
Rgb shade(float nx, float ny, float nz, const Rgb& rgb) {
    const float inv = 1.0f / npMax(std::sqrt(nx * nx + ny * ny + nz * nz), 1e-6f);
    const float lx = kLight[0], ly = kLight[1], lz = kLight[2];
    const float ndl_raw = (nx * lx + ny * ly + nz * lz) * inv;
    const float ndl = clip(ndl_raw, 0.0f, 1.0f);
    const float rv = clip(2.0f * ndl_raw * nz * inv - lz, 0.0f, 1.0f) * (ndl_raw > 0.0f ? 1.0f : 0.0f);
    const float x2 = rv * rv;
    const float x4 = x2 * x2;
    const float x8 = x4 * x4;
    const float diffuse = kAmbient + kDiffuse * ndl;
    const float spec = kSpecular * x8 * x8 * x8;
    Rgb out;
    for (int k = 0; k < 3; ++k) {
        out[k] = clip(rgb[k] * diffuse + spec, 0.0f, 1.0f);
    }
    return out;
}

class ZPass {
public:
    ZPass(const Cam& cam, unsigned threads) : cam_(cam), threads_(threads) {}

    void add(const Colored& it, std::size_t cull_min) {
        Pass p;
        p.opacity = f32(it.opacity);
        if (const auto* m = std::get_if<Mesh>(&it.prim)) {
            if (auto d = meshPass(cam_, *m, it.rgb)) {
                p.count = d->tris.size();
                p.data = std::move(*d);
            }
        } else if (const auto* s = std::get_if<Spheres>(&it.prim)) {
            if (auto d = spherePass(cam_, *s, it.rgb, cull_min)) {
                p.count = d->spheres.size();
                p.data = std::move(*d);
            }
        } else if (const auto* k = std::get_if<Sticks>(&it.prim)) {
            if (auto d = stickPass(cam_, *k, it.rgb)) {
                p.count = d->sticks.size();
                p.data = std::move(*d);
            }
        }
        if (p.count == 0) {
            return;
        }
        if (n_prim_ + p.count >= (std::uint64_t{1} << kIdBits)) {
            throw std::runtime_error("scene too large for the built-in renderer");
        }
        p.base = n_prim_;
        n_prim_ += p.count;
        passes_.push_back(std::move(p));
    }

    // Rasterise every pass and resolve the front-most fragment per pixel.
    Fragments resolve() const {
        const std::int64_t W = cam_.W, H = cam_.H;
        std::vector<std::int64_t> zbuf(static_cast<std::size_t>(W * H), -1);
        if (!passes_.empty()) {
            forBands(H, threads_, [&](std::size_t, std::int64_t y0, std::int64_t y1) { raster(zbuf, y0, y1); });
        }
        std::vector<Fragments> parts(bandCount(H, threads_));
        forBands(H, threads_, [&](std::size_t b, std::int64_t y0, std::int64_t y1) {
            Fragments& f = parts[b];
            for (std::int64_t pix = y0 * W; pix < y1 * W; ++pix) {
                const std::int64_t key = zbuf[static_cast<std::size_t>(pix)];
                if (key >= 0) {
                    resolvePixel(f, pix, static_cast<std::uint64_t>(key) & ((std::uint64_t{1} << kIdBits) - 1));
                }
            }
        });
        Fragments out;
        for (Fragments& f : parts) {
            out.pix.insert(out.pix.end(), f.pix.begin(), f.pix.end());
            out.rgb.insert(out.rgb.end(), f.rgb.begin(), f.rgb.end());
            out.depth.insert(out.depth.end(), f.depth.begin(), f.depth.end());
            out.owner.insert(out.owner.end(), f.owner.begin(), f.owner.end());
            out.alpha.insert(out.alpha.end(), f.alpha.begin(), f.alpha.end());
        }
        return out;
    }

private:
    void raster(std::vector<std::int64_t>& zbuf, std::int64_t ya, std::int64_t yb) const {
        const std::int64_t W = cam_.W, H = cam_.H;
        auto write = [&](std::int64_t x, std::int64_t y, float z, std::uint64_t id) {
            std::int64_t& k = zbuf[static_cast<std::size_t>(y * W + x)];
            k = std::max(k, keyOf(cam_, z, id));
        };
        for (const Pass& p : passes_) {
            if (const auto* m = std::get_if<MeshPass>(&p.data)) {
                for (std::size_t r = 0; r < m->tris.size(); ++r) {
                    const TriSetup& t = m->tris[r];
                    const std::int64_t y0 = std::max(t.by0, ya), y1 = std::min(t.by1, yb - 1);
                    for (std::int64_t y = y0; y <= y1; ++y) {
                        std::int64_t x0 = 0, x1 = 0;
                        if (!spanRow(t.h, static_cast<double>(y) + 0.5, W, x0, x1)) {
                            continue;
                        }
                        // widened edges of a sliver meet far outside it: keep to its box
                        x0 = std::max(x0, t.bx0);
                        x1 = std::min(x1, t.bx1);
                        for (std::int64_t x = x0; x <= x1; ++x) {
                            write(x, y, meshFrag(t, x, y).z, p.base + r);
                        }
                    }
                }
            } else if (const auto* s = std::get_if<SpherePass>(&p.data)) {
                for (std::size_t r = 0; r < s->spheres.size(); ++r) {
                    const Sphere& sp = s->spheres[r];
                    const std::int64_t y0 = std::max({sp.cy - sp.n, ya, std::int64_t{0}});
                    const std::int64_t y1 = std::min({sp.cy + sp.n, yb - 1, H - 1});
                    const std::int64_t x0 = std::max(sp.cx - sp.n, std::int64_t{0});
                    const std::int64_t x1 = std::min(sp.cx + sp.n, W - 1);
                    for (std::int64_t y = y0; y <= y1; ++y) {
                        for (std::int64_t x = x0; x <= x1; ++x) {
                            const SphereFrag f = sphereFrag(sp, x, y);
                            if (f.in) {
                                write(x, y, f.z, p.base + r);
                            }
                        }
                    }
                }
            } else if (const auto* k = std::get_if<StickPass>(&p.data)) {
                for (std::size_t r = 0; r < k->sticks.size(); ++r) {
                    const Stick& st = k->sticks[r];
                    const std::int64_t y0 = std::max(st.by0, ya), y1 = std::min(st.by1, yb - 1);
                    for (std::int64_t y = y0; y <= y1; ++y) {
                        std::int64_t x0 = 0, x1 = 0;
                        if (!spanRow(st.h, static_cast<float>(y) + 0.5f, W, x0, x1)) {
                            continue;
                        }
                        for (std::int64_t x = x0; x <= x1; ++x) {
                            const StickFrag f = stickFrag(st, k->H_half, k->s, x, y);
                            if (f.in) {
                                write(x, y, f.z, p.base + r);
                            }
                        }
                    }
                }
            }
        }
    }

    void resolvePixel(Fragments& out, std::int64_t pix, std::uint64_t id) const {
        const auto it = std::upper_bound(passes_.begin(), passes_.end(), id,
                                         [](std::uint64_t v, const Pass& p) { return v < p.base; });
        const Pass& p = *(it - 1);
        const std::size_t r = static_cast<std::size_t>(id - p.base);
        const std::int64_t x = pix % cam_.W, y = pix / cam_.W;
        float nx = 0, ny = 0, nz = 0, z = 0;
        Rgb rgb{};
        std::int64_t owner = -1;
        if (const auto* m = std::get_if<MeshPass>(&p.data)) {
            const TriSetup& t = m->tris[r];
            const MeshFrag f = meshFrag(t, x, y);
            float a = clip(f.l0, 0.0f, 1.0f), b = clip(f.l1, 0.0f, 1.0f);
            const float s = npMax(a + b, 1.0f);
            a = a / s;
            b = b / s;
            const float c = 1.0f - a - b;
            const auto v0 = static_cast<std::size_t>(t.v[0]), v1 = static_cast<std::size_t>(t.v[1]),
                       v2 = static_cast<std::size_t>(t.v[2]);
            nx = m->nx[v0] * a + m->nx[v1] * b + m->nx[v2] * c;
            ny = m->ny[v0] * a + m->ny[v1] * b + m->ny[v2] * c;
            nz = m->nz[v0] * a + m->nz[v1] * b + m->nz[v2] * c;
            // colour of the residue with the largest share of the barycentric weight: crisp,
            // straight boundaries between residues (no zig-zag along diagonals)
            const std::int64_t o0 = m->own[v0], o1 = m->own[v1], o2 = m->own[v2];
            auto eq = [](std::int64_t u, std::int64_t v) { return u == v ? 1.0f : 0.0f; };
            const float w0 = a + b * eq(o1, o0) + c * eq(o2, o0);
            const float w1 = b + a * eq(o0, o1) + c * eq(o2, o1);
            const float w2 = c + a * eq(o0, o2) + b * eq(o1, o2);
            const std::int64_t best = w0 >= w1 ? (w0 >= w2 ? o0 : o2) : (w1 >= w2 ? o1 : o2);
            // within that residue, the nearest vertex (per-atom colours, e.g. surfaces)
            const float s0 = o0 == best ? a : -1.0f, s1 = o1 == best ? b : -1.0f, s2 = o2 == best ? c : -1.0f;
            const std::size_t vert = s0 >= s1 ? (s0 >= s2 ? v0 : v2) : (s1 >= s2 ? v1 : v2);
            rgb = m->col[vert];
            owner = m->own[vert];
            z = f.z;
        } else if (const auto* s = std::get_if<SpherePass>(&p.data)) {
            const Sphere& sp = s->spheres[r];
            const SphereFrag f = sphereFrag(sp, x, y);
            nx = f.nx;
            ny = f.ny;
            nz = f.nz;
            z = f.z;
            rgb = s->col[sp.index];
            owner = s->own[sp.index];
        } else if (const auto* k = std::get_if<StickPass>(&p.data)) {
            const Stick& st = k->sticks[r];
            const StickFrag f = stickFrag(st, k->H_half, k->s, x, y);
            nx = f.nx;
            ny = f.ny;
            nz = f.nz;
            z = f.z;
            rgb = k->col[st.index];
            owner = k->own[st.index];
        }
        out.pix.push_back(pix);
        out.rgb.push_back(shade(nx, ny, nz, rgb));
        out.depth.push_back(z);
        out.owner.push_back(owner);
        out.alpha.push_back(p.opacity);
    }

    const Cam& cam_;
    unsigned threads_;
    std::vector<Pass> passes_;
    std::uint64_t n_prim_ = 0;
};

Fragments draw(const Cam& cam, const std::vector<const Colored*>& items, const SoftwareRenderOptions& options) {
    ZPass zp(cam, workerCount(options, static_cast<std::size_t>(cam.W * cam.H)));
    for (const Colored* it : items) {
        zp.add(*it, options.cull_min);
    }
    return zp.resolve();
}

// ---- compositing -------------------------------------------------------------------------------

std::vector<float> glowValues(const Shading& shading) {
    std::vector<float> g;
    if (shading.glow) {
        g.reserve(shading.glow->size());
        for (double v : *shading.glow) {
            g.push_back(std::isfinite(v) ? f32(v) : 0.0f);  // RN-28: non-finite glow counts as none
        }
    }
    return g;
}

float glowOf(const std::vector<float>& g, std::int64_t owner) {
    return owner >= 0 && owner < static_cast<std::int64_t>(g.size()) ? g[static_cast<std::size_t>(owner)] : 0.0f;
}

// Glowing residues lighten toward the glow colour (screen blend).
Rgb tint(const Rgb& rgb, float gi, const Rgb& glow_rgb) {
    const float t = static_cast<float>(kGlowTint) * gi;
    Rgb out;
    for (int k = 0; k < 3; ++k) {
        out[k] = 1.0f - (1.0f - rgb[k]) * (1.0f - t * glow_rgb[k]);
    }
    return out;
}

// Depth cue toward the panel colour, over the panel's depth range (as the 2D trace).
Rgb fog(const Cam& c, const Rgb& rgb, float depth, const Shading& shading) {
    const float span = f32(std::max(c.fog_hi - c.fog_lo, 1e-6));
    const float f = f32(shading.fog_strength) * (1.0f - clip((depth - f32(c.fog_lo)) / span, 0.0f, 1.0f));
    const Rgb fr = rgb32(shading.fog_rgb);
    Rgb out;
    for (int k = 0; k < 3; ++k) {
        out[k] = rgb[k] * (1.0f - f) + fr[k] * f;
    }
    return out;
}

// 1 / peak of a blurred disc of radius `rad` (a glowing residue reaches full alpha). RN-29: for the
// exact sizes (MDATK cached it under the sizes rounded to 0.01).
double discNorm(double sig, double rad) {
    const std::int64_t n = toIndex(std::ceil(3 * sig + rad)) + 2;
    const auto side = static_cast<std::size_t>(2 * n + 1);
    const double r = std::max(rad, 0.5);
    const double r2 = r * r;
    std::vector<float> d(side * side);
    for (std::int64_t yy = -n; yy <= n; ++yy) {
        for (std::int64_t xx = -n; xx <= n; ++xx) {
            d[static_cast<std::size_t>((yy + n) * (2 * n + 1) + (xx + n))] =
                static_cast<double>(xx * xx + yy * yy) <= r2 ? 1.0f : 0.0f;
        }
    }
    const Shape3 shape{1, side, side};
    gaussianFilter1d(d, shape, 1, sig, 2.5);
    gaussianFilter1d(d, shape, 2, sig, 2.5);
    float peak = -std::numeric_limits<float>::infinity();
    for (float v : d) {
        peak = npMax(peak, v);
    }
    return 1.0 / std::max(static_cast<double>(peak), 1e-6);
}

// Bilinear upsampling of a coarse grid (values at cell centres) by f along its rows (axis 0),
// edges clamped (MDATK's _upsample, one axis).
std::vector<float> upsampleRows(const std::vector<float>& x, std::size_t rows, std::size_t cols, int f) {
    std::vector<float> out(rows * static_cast<std::size_t>(f) * cols);
    for (int j = 0; j < f; ++j) {
        const double o = (j + 0.5) / f - 0.5;  // offset of the fine pixel from its coarse centre
        const float wa = o < 0 ? f32(-o) : f32(1 - o);
        const float wb = o < 0 ? f32(1 + o) : f32(o);
        for (std::size_t i = 0; i < rows; ++i) {
            const std::size_t ia = o < 0 ? (i == 0 ? 0 : i - 1) : i;
            const std::size_t ib = o < 0 ? i : std::min(i + 1, rows - 1);
            float* dst = out.data() + (i * static_cast<std::size_t>(f) + static_cast<std::size_t>(j)) * cols;
            const float* a = x.data() + ia * cols;
            const float* b = x.data() + ib * cols;
            for (std::size_t c = 0; c < cols; ++c) {
                dst[c] = a[c] * wa + b[c] * wb;
            }
        }
    }
    return out;
}

std::vector<float> transposed(const std::vector<float>& x, std::size_t rows, std::size_t cols) {
    std::vector<float> out(x.size());
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            out[c * rows + r] = x[r * cols + c];
        }
    }
    return out;
}

RgbaImage blank(std::int64_t W, std::int64_t H) {
    RgbaImage img;
    img.width = static_cast<std::size_t>(W);
    img.height = static_cast<std::size_t>(H);
    img.rgba.assign(img.width * img.height * 4, 0);
    return img;
}

// scipy's gaussian_filter (mode "constant", truncate 2.5) of a rows x cols float32 image: axis 0,
// then axis 1, each line on its own (blocks of lines on several threads).
void gaussian2d(std::vector<float>& data, std::size_t rows, std::size_t cols, double sigma, unsigned threads) {
    auto alongRows = [&](std::vector<float>& d, std::size_t nr, std::size_t nc) {  // each row filtered
        forBands(static_cast<std::int64_t>(nr), threads, [&](std::size_t, std::int64_t r0, std::int64_t r1) {
            const auto a = static_cast<std::size_t>(r0), b = static_cast<std::size_t>(r1);
            if (a == b) {
                return;
            }
            std::vector<float> block(d.begin() + static_cast<std::ptrdiff_t>(a * nc),
                                     d.begin() + static_cast<std::ptrdiff_t>(b * nc));
            gaussianFilter1d(block, Shape3{1, b - a, nc}, 2, sigma, 2.5);
            std::copy(block.begin(), block.end(), d.begin() + static_cast<std::ptrdiff_t>(a * nc));
        });
    };
    std::vector<float> t = transposed(data, rows, cols);  // axis 0: the columns
    alongRows(t, cols, rows);
    data = transposed(t, cols, rows);
    alongRows(data, rows, cols);
}

// A soft halo of the glow colour around glowing pixels (under the molecule).
void glowHalo(RgbaImage& out, const Cam& c, const std::vector<std::int64_t>& pix, const std::vector<float>& gi,
              const Shading& shading, unsigned threads) {
    const std::int64_t W = c.W, H = c.H;
    std::int64_t pymin = std::numeric_limits<std::int64_t>::max(), pymax = -1;
    std::int64_t pxmin = std::numeric_limits<std::int64_t>::max(), pxmax = -1;
    for (std::int64_t p : pix) {
        pymin = std::min(pymin, p / W);
        pymax = std::max(pymax, p / W);
        pxmin = std::min(pxmin, p % W);
        pxmax = std::max(pxmax, p % W);
    }
    const std::int64_t y0 = std::max<std::int64_t>(pymin - 2, 0), y1 = std::min(pymax + 3, H);
    const std::int64_t x0 = std::max<std::int64_t>(pxmin - 2, 0), x1 = std::min(pxmax + 3, W);
    const std::int64_t f = 4;
    const std::int64_t pad = toIndex(std::ceil(2.5 * kGlowRadius * c.s));
    const std::int64_t gy0 = (std::max<std::int64_t>(y0 - pad, 0) / f) * f;
    const std::int64_t gx0 = (std::max<std::int64_t>(x0 - pad, 0) / f) * f;
    const std::int64_t gy1 = std::min(y1 + pad, H), gx1 = std::min(x1 + pad, W);
    const std::int64_t hh = (gy1 - gy0 + f - 1) / f, ww = (gx1 - gx0 + f - 1) / f;
    std::vector<float> low(static_cast<std::size_t>(hh * ww), 0.0f);
    for (std::size_t i = 0; i < pix.size(); ++i) {
        const std::int64_t py = pix[i] / W, px = pix[i] % W;
        float& v = low[static_cast<std::size_t>(((py - gy0) / f) * ww + (px - gx0) / f)];
        v = npMax(v, gi[i]);
    }
    const double sig = kGlowRadius * c.s / static_cast<double>(f);
    gaussian2d(low, static_cast<std::size_t>(hh), static_cast<std::size_t>(ww), sig, threads);
    const float norm = f32(discNorm(sig, kGlowDisc * c.s / static_cast<double>(f)));
    for (float& v : low) {
        v = clip(v * norm, 0.0f, 1.0f);
    }
    const auto uh = static_cast<std::size_t>(hh), uw = static_cast<std::size_t>(ww);
    std::vector<float> a = upsampleRows(low, uh, uw, static_cast<int>(f));       // (hh f) x ww
    a = upsampleRows(transposed(a, uh * 4, uw), uw, uh * 4, static_cast<int>(f));  // (ww f) x (hh f)
    const Rgb gr = rgb32(shading.glow_rgb);
    std::uint8_t g8[3];
    for (int k = 0; k < 3; ++k) {
        g8[k] = toByte(clip(gr[k] * 255.0f + 0.5f, 0.0f, 255.0f));
    }
    const float scale = static_cast<float>(kGlowAlpha * 255);
    const std::size_t rows = uh * 4;  // a[c][r]
    for (std::int64_t y = gy0; y < gy1; ++y) {
        for (std::int64_t x = gx0; x < gx1; ++x) {
            const float v = a[static_cast<std::size_t>(x - gx0) * rows + static_cast<std::size_t>(y - gy0)];
            std::uint8_t* q = out.pixel(static_cast<std::size_t>(x), static_cast<std::size_t>(y));
            q[0] = g8[0];
            q[1] = g8[1];
            q[2] = g8[2];
            q[3] = toByte(clip(v * scale + 0.5f, 0.0f, 255.0f));
        }
    }
}

// Alpha anti-aliasing of the edges in rows y0..y1, columns x0..x1: where a pixel's alpha differs from
// its 3 x 3 mean by more than 1, the average of the two. RN-26: the mean sees the image around the
// region (MDATK filtered the region alone, with zeros outside it).
void antialias(RgbaImage& out, std::int64_t x0, std::int64_t x1, std::int64_t y0, std::int64_t y1, unsigned threads) {
    const auto W = static_cast<std::int64_t>(out.width), H = static_cast<std::int64_t>(out.height);
    const std::int64_t ey0 = std::max<std::int64_t>(y0 - 1, 0), ey1 = std::min(y1 + 1, H);
    const std::int64_t ex0 = std::max<std::int64_t>(x0 - 1, 0), ex1 = std::min(x1 + 1, W);
    const auto h = static_cast<std::size_t>(ey1 - ey0), w = static_cast<std::size_t>(ex1 - ex0);
    if (h == 0 || w == 0) {
        return;
    }
    std::vector<float> a(h * w), t(h * w), a3(h * w);
    for (std::size_t r = 0; r < h; ++r) {
        for (std::size_t c = 0; c < w; ++c) {
            a[r * w + c] = static_cast<float>(out.pixel(static_cast<std::size_t>(ex0) + c, static_cast<std::size_t>(ey0) + r)[3]);
        }
    }
    auto lines = [&](std::size_t count, auto&& fn) {
        forBands(static_cast<std::int64_t>(count), threads, [&](std::size_t, std::int64_t i0, std::int64_t i1) {
            for (auto i = static_cast<std::size_t>(i0); i < static_cast<std::size_t>(i1); ++i) {
                fn(i);
            }
        });
    };
    // axis 0, then axis 1 (lines are independent)
    lines(w, [&](std::size_t c) { uniformFilter1dLine(a.data() + c, w, h, 3, t.data() + c); });
    lines(h, [&](std::size_t r) { uniformFilter1dLine(t.data() + r * w, 1, w, 3, a3.data() + r * w); });
    lines(static_cast<std::size_t>(y1 - y0), [&](std::size_t row) {
        const std::int64_t y = y0 + static_cast<std::int64_t>(row);
        for (std::int64_t x = x0; x < x1; ++x) {
            const std::size_t i = static_cast<std::size_t>(y - ey0) * w + static_cast<std::size_t>(x - ex0);
            const float v = a[i];
            const float m = a3[i];
            const float r = std::fabs(m - v) > 1.0f ? 0.5f * (v + m) : v;
            out.pixel(static_cast<std::size_t>(x), static_cast<std::size_t>(y))[3] = toByte(r);
        }
    });
}

struct Front {  // the front-most translucent layer
    std::vector<std::int64_t> pix;
    std::vector<Rgb> rgb;
    std::vector<float> depth;
    std::vector<std::int64_t> owner;
    std::vector<float> alpha;
};

// Covered pixels -> straight-alpha RGBA: fog, glow (tint + halo), ink outlines at depth jumps and
// silhouettes, alpha anti-aliasing of the edges. `extra` (translucent pixels drawn later) adds to
// the halo.
RgbaImage composite(const Cam& c, Fragments frag, const Shading& shading, float outline, bool aa,
                    const Front* extra, unsigned threads) {
    // per pixel (or per row) work in parallel: every element is computed on its own
    auto each = [&](std::size_t count, auto&& fn) {
        forBands(static_cast<std::int64_t>(count), threads, [&](std::size_t, std::int64_t i0, std::int64_t i1) {
            for (auto i = static_cast<std::size_t>(i0); i < static_cast<std::size_t>(i1); ++i) {
                fn(i);
            }
        });
    };
    const std::int64_t W = c.W, H = c.H;
    RgbaImage out = blank(W, H);
    const std::vector<float> g = glowValues(shading);
    std::vector<float> gi(frag.pix.size());
    float gmax = 0.0f;
    for (std::size_t i = 0; i < gi.size(); ++i) {
        gi[i] = glowOf(g, frag.owner[i]);
        gmax = npMax(gmax, gi[i]);
    }
    {
        std::vector<std::int64_t> hpix = frag.pix;
        std::vector<float> hgi = gi;
        float hmax = gmax;
        if (extra) {
            for (std::size_t i = 0; i < extra->pix.size(); ++i) {
                hpix.push_back(extra->pix[i]);
                hgi.push_back(glowOf(g, extra->owner[i]));
                hmax = npMax(hmax, hgi.back());
            }
        }
        if (hmax > 0.01f) {
            glowHalo(out, c, hpix, hgi, shading, threads);
        }
    }
    if (frag.pix.empty()) {
        return out;
    }
    const std::size_t n = frag.pix.size();
    if (shading.fog_strength > 0) {
        each(n, [&](std::size_t i) { frag.rgb[i] = fog(c, frag.rgb[i], frag.depth[i], shading); });
    }
    std::int64_t pymin = std::numeric_limits<std::int64_t>::max(), pymax = -1;
    std::int64_t pxmin = std::numeric_limits<std::int64_t>::max(), pxmax = -1;
    for (std::int64_t p : frag.pix) {
        pymin = std::min(pymin, p / W);
        pymax = std::max(pymax, p / W);
        pxmin = std::min(pxmin, p % W);
        pxmax = std::max(pxmax, p % W);
    }
    const std::int64_t y0 = std::max<std::int64_t>(pymin - 2, 0), y1 = std::min(pymax + 3, H);
    const std::int64_t x0 = std::max<std::int64_t>(pxmin - 2, 0), x1 = std::min(pxmax + 3, W);
    if (gmax > 0.01f) {
        const Rgb gr = rgb32(shading.glow_rgb);
        each(n, [&](std::size_t i) { frag.rgb[i] = tint(frag.rgb[i], gi[i], gr); });
    }
    if (outline > 0) {
        const std::int64_t h = y1 - y0, w = x1 - x0;
        auto li = [&](std::size_t i) { return (frag.pix[i] / W - y0) * w + (frag.pix[i] % W - x0); };
        std::vector<float> dz(static_cast<std::size_t>(h * w), f32(c.zmin - 1e3));
        std::vector<std::uint8_t> cover(dz.size(), 0);
        for (std::size_t i = 0; i < n; ++i) {
            dz[static_cast<std::size_t>(li(i))] = frag.depth[i];
            cover[static_cast<std::size_t>(li(i))] = 1;
        }
        const float t = 1.6f;
        std::vector<std::uint8_t> edge(dz.size(), 0);
        each(static_cast<std::size_t>(h), [&](std::size_t row) {
            const auto r = static_cast<std::int64_t>(row);
            for (std::int64_t q = 0; q < w; ++q) {
                const float d = dz[static_cast<std::size_t>(r * w + q)];
                bool e = false;
                if (r > 0) {
                    e = e || dz[static_cast<std::size_t>((r - 1) * w + q)] - d > t;
                }
                if (r + 1 < h) {
                    e = e || dz[static_cast<std::size_t>((r + 1) * w + q)] - d > t;
                }
                if (q > 0) {
                    e = e || dz[static_cast<std::size_t>(r * w + q - 1)] - d > t;
                }
                if (q + 1 < w) {
                    e = e || dz[static_cast<std::size_t>(r * w + q + 1)] - d > t;
                }
                edge[static_cast<std::size_t>(r * w + q)] = e ? 1 : 0;
            }
        });
        const float keep = f32(1 - static_cast<double>(outline));
        const float ink = kInk * outline;
        each(n, [&](std::size_t i) {
            if (edge[static_cast<std::size_t>(li(i))]) {
                for (int k = 0; k < 3; ++k) {
                    frag.rgb[i][k] = frag.rgb[i][k] * keep + ink;
                }
            }
        });
        const std::uint8_t ink8 = toByte(clip(kInk * 255.0f + 0.5f, 0.0f, 255.0f));
        const auto sil_alpha = static_cast<std::uint8_t>(static_cast<int>(static_cast<double>(outline) * 0.8 * 255));
        for (std::int64_t r = 0; r < h; ++r) {
            for (std::int64_t q = 0; q < w; ++q) {
                const auto j = static_cast<std::size_t>(r * w + q);
                if (edge[j] && !cover[j]) {
                    std::uint8_t* px = out.pixel(static_cast<std::size_t>(q + x0), static_cast<std::size_t>(r + y0));
                    px[0] = px[1] = px[2] = ink8;
                    px[3] = std::max(px[3], sil_alpha);
                }
            }
        }
    }
    each(n, [&](std::size_t i) {
        std::uint8_t* px = out.rgba.data() + 4 * static_cast<std::size_t>(frag.pix[i]);
        for (int k = 0; k < 3; ++k) {
            px[k] = toByte(clip(frag.rgb[i][k] * 255.0f + 0.5f, 0.0f, 255.0f));
        }
        px[3] = 255;
    });
    if (aa) {
        antialias(out, x0, x1, y0, y1, threads);
    }
    return out;
}

// Translucent pixels (straight alpha) over the opaque image.
void blend(const Cam& c, RgbaImage& out, const Front& front, const Shading& shading) {
    const std::vector<float> g = glowValues(shading);
    const Rgb gr = rgb32(shading.glow_rgb);
    bool any = false;
    for (std::int64_t o : front.owner) {
        any = any || glowOf(g, o) != 0.0f;
    }
    for (std::size_t i = 0; i < front.pix.size(); ++i) {
        Rgb rgb = shading.fog_strength > 0 ? fog(c, front.rgb[i], front.depth[i], shading) : front.rgb[i];
        if (any) {
            rgb = tint(rgb, glowOf(g, front.owner[i]), gr);
        }
        const float a = front.alpha[i];
        std::uint8_t* px = out.rgba.data() + 4 * static_cast<std::size_t>(front.pix[i]);
        const float dst_a = static_cast<float>(px[3]) / 255.0f;
        const float out_a = a + dst_a * (1.0f - a);
        for (int k = 0; k < 3; ++k) {
            const float dst = static_cast<float>(px[k]) / 255.0f;
            const float col = (rgb[k] * a + dst * dst_a * (1.0f - a)) / npMax(out_a, 1e-6f);
            px[k] = toByte(clip(col * 255.0f + 0.5f, 0.0f, 255.0f));
        }
        px[3] = toByte(clip(out_a * 255.0f + 0.5f, 0.0f, 255.0f));
    }
}

}  // namespace

// ---- public ------------------------------------------------------------------------------------

RgbaImage SoftwareRenderer::render(const std::vector<Colored>& items, const Viewport& vp, const Shading& shading,
                                   RenderQuality quality) const {
    const bool fast = quality == RenderQuality::Fast;
    const Cam cam = makeCam(vp, items, fast);
    std::vector<const Colored*> opaque, glass;
    for (const Colored& it : items) {
        (it.opacity >= 0.999 ? opaque : glass).push_back(&it);
    }
    Fragments solid = draw(cam, opaque, options_);
    std::optional<Front> front;
    if (!glass.empty()) {  // the front-most translucent surface in front of the opaque one, per pixel
        const Fragments g = draw(cam, glass, options_);
        std::vector<float> opaque_z(static_cast<std::size_t>(cam.W * cam.H), -std::numeric_limits<float>::infinity());
        for (std::size_t i = 0; i < solid.pix.size(); ++i) {
            opaque_z[static_cast<std::size_t>(solid.pix[i])] = solid.depth[i];
        }
        front.emplace();
        for (std::size_t i = 0; i < g.pix.size(); ++i) {
            if (g.depth[i] > opaque_z[static_cast<std::size_t>(g.pix[i])] + 1e-3f) {
                front->pix.push_back(g.pix[i]);
                front->rgb.push_back(g.rgb[i]);
                front->depth.push_back(g.depth[i]);
                front->owner.push_back(g.owner[i]);
                front->alpha.push_back(g.alpha[i]);
            }
        }
    }
    RgbaImage out = composite(cam, std::move(solid), shading, fast ? 0.0f : 0.5f, !fast, front ? &*front : nullptr,
                              workerCount(options_, static_cast<std::size_t>(cam.W * cam.H)));
    if (front && !front->pix.empty()) {
        blend(cam, out, *front, shading);
    }
    if (fast) {  // back to the requested size
        RgbaImage full = blank(vp.width, vp.height);
        for (std::size_t y = 0; y < full.height && y / 2 < out.height; ++y) {
            for (std::size_t x = 0; x < full.width && x / 2 < out.width; ++x) {
                std::copy_n(out.pixel(x / 2, y / 2), 4, full.pixel(x, y));
            }
        }
        out = std::move(full);
    }
    return out;
}

ResolvedPixels SoftwareRenderer::resolve(const std::vector<Colored>& items, const Viewport& vp) const {
    const Cam cam = makeCam(vp, items, false);
    std::vector<const Colored*> all;
    for (const Colored& it : items) {
        all.push_back(&it);
    }
    Fragments f = draw(cam, all, options_);
    ResolvedPixels out;
    out.width = static_cast<std::size_t>(cam.W);
    out.height = static_cast<std::size_t>(cam.H);
    out.scale = cam.s;
    out.pixel = std::move(f.pix);
    out.rgb = std::move(f.rgb);
    out.depth = std::move(f.depth);
    out.owner = std::move(f.owner);
    out.alpha = std::move(f.alpha);
    return out;
}

std::size_t SoftwareRenderer::culledSpheres(const Spheres& spheres, const std::vector<Colored>& items,
                                            const Viewport& vp) const {
    const Cam cam = makeCam(vp, items, false);
    checkSpheres(spheres, spheres.xyz.size());
    const SphereInput in = sphereInput(cam, spheres);
    const std::vector<bool> keep = unoccluded(cam, in.P, in.rad, options_.cull_min);
    return static_cast<std::size_t>(std::count(keep.begin(), keep.end(), false));
}

std::vector<Colored> tinted(const std::vector<Colored>& items, const Shading& shading) {
    const std::vector<float> g = glowValues(shading);
    const Rgb gr = rgb32(shading.glow_rgb);
    std::vector<Colored> out;
    out.reserve(items.size());
    for (const Colored& it : items) {
        const std::vector<std::int32_t>* owner_ptr =
            std::visit([](const auto& p) { return &p.owner; }, it.prim);
        const std::vector<std::int32_t>& owner = *owner_ptr;
        float gmax = 0.0f;
        for (std::int32_t o : owner) {
            gmax = npMax(gmax, glowOf(g, o));
        }
        Colored c = it;
        if (gmax > 0.01f) {
            checkSize(it.rgb.size(), owner.size(), "colour");
            for (std::size_t i = 0; i < c.rgb.size(); ++i) {
                const Rgb t = tint(rgb32(it.rgb[i]), glowOf(g, owner[i]), gr);
                c.rgb[i] = {t[0], t[1], t[2]};
            }
        }
        out.push_back(std::move(c));
    }
    return out;
}

std::optional<RgbaImage> glowUnderlay(const std::vector<Colored>& items, const Viewport& vp, const Shading& shading,
                                      const SoftwareRenderOptions& options) {
    const std::vector<float> g = glowValues(shading);
    auto hot = [&](std::int32_t o) { return glowOf(g, o) > 0.01f; };
    std::vector<Colored> glowing;
    for (const Colored& it : items) {
        if (const auto* m = std::get_if<Mesh>(&it.prim)) {
            Mesh sub = *m;
            sub.tris.clear();
            for (const molecule::Tri& t : m->tris) {
                bool any = false;
                for (std::int32_t v : t) {
                    any = any || (v >= 0 && static_cast<std::size_t>(v) < m->owner.size() &&
                                  hot(m->owner[static_cast<std::size_t>(v)]));
                }
                if (any) {
                    sub.tris.push_back(t);
                }
            }
            bool any_vertex = false;
            for (std::int32_t o : m->owner) {
                any_vertex = any_vertex || hot(o);
            }
            if (any_vertex) {
                glowing.push_back({sub, it.rgb, it.opacity, it.representation});
            }
            continue;
        }
        auto subset = [&](const auto& p, auto make) {
            std::vector<std::size_t> keep;
            for (std::size_t i = 0; i < p.owner.size(); ++i) {
                if (hot(p.owner[i])) {
                    keep.push_back(i);
                }
            }
            if (keep.empty()) {
                return;
            }
            checkSize(it.rgb.size(), p.owner.size(), "colour");
            std::vector<Vec3> rgb;
            for (std::size_t i : keep) {
                rgb.push_back(it.rgb[i]);
            }
            glowing.push_back({make(p, keep), std::move(rgb), it.opacity, it.representation});
        };
        if (const auto* s = std::get_if<Spheres>(&it.prim)) {
            checkSpheres(*s, it.rgb.size());
            subset(*s, [](const Spheres& p, const std::vector<std::size_t>& keep) {
                Spheres q;
                for (std::size_t i : keep) {
                    q.xyz.push_back(p.xyz[i]);
                    q.radius.push_back(p.radius[i]);
                    q.owner.push_back(p.owner[i]);
                    q.atom.push_back(i < p.atom.size() ? p.atom[i] : -1);
                }
                return q;
            });
        } else if (const auto* k = std::get_if<Sticks>(&it.prim)) {
            subset(*k, [](const Sticks& p, const std::vector<std::size_t>& keep) {
                Sticks q;
                for (std::size_t i : keep) {
                    q.a.push_back(p.a[i]);
                    q.b.push_back(p.b[i]);
                    q.radius.push_back(p.radius[i]);
                    q.owner.push_back(p.owner[i]);
                    q.atom.push_back(i < p.atom.size() ? p.atom[i] : -1);
                }
                return q;
            });
        }
    }
    if (glowing.empty()) {
        return std::nullopt;
    }
    const Cam cam = makeCam(vp, items, false);
    std::vector<const Colored*> ptrs;
    for (const Colored& it : glowing) {
        ptrs.push_back(&it);
    }
    const Fragments f = draw(cam, ptrs, options);
    RgbaImage out = blank(cam.W, cam.H);
    if (!f.pix.empty()) {
        std::vector<float> gi(f.pix.size());
        for (std::size_t i = 0; i < gi.size(); ++i) {
            gi[i] = glowOf(g, f.owner[i]);
        }
        glowHalo(out, cam, f.pix, gi, shading, workerCount(options, static_cast<std::size_t>(cam.W * cam.H)));
    }
    return out;
}

RgbaImage over(const RgbaImage& top, const RgbaImage& under) {
    if (top.width != under.width || top.height != under.height || top.rgba.size() != under.rgba.size()) {
        throw std::invalid_argument("over: the images differ in size");
    }
    RgbaImage out = top;
    for (std::size_t i = 0; i + 3 < top.rgba.size(); i += 4) {
        const float ta = static_cast<float>(top.rgba[i + 3]) / 255.0f;
        const float ua = static_cast<float>(under.rgba[i + 3]) / 255.0f;
        const float a = ta + ua * (1.0f - ta);
        for (std::size_t k = 0; k < 3; ++k) {
            const float rgb = (static_cast<float>(top.rgba[i + k]) * ta +
                               static_cast<float>(under.rgba[i + k]) * ua * (1.0f - ta)) /
                              npMax(a, 1e-6f);
            out.rgba[i + k] = toByte(clip(rgb + 0.5f, 0.0f, 255.0f));
        }
        out.rgba[i + 3] = toByte(clip(a * 255.0f + 0.5f, 0.0f, 255.0f));
    }
    return out;
}

}  // namespace mdbridge::render::molecule
