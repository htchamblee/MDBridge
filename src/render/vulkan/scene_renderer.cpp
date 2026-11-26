// The GPU renderer of the molecule scene (see scene_renderer.hpp). The CPU sets every primitive up
// as the software renderer does (src/render/molecule/software.cpp: makeCam, meshPass, spherePass,
// stickPass; MDATK's own code, ported there) and writes one record per primitive; the shaders in
// src/render/shaders (molecule_*) do the per-pixel work in float32.
//
// Records (32-bit words in the primitive buffer, binding 0; every record starts with its pixel box
// x0, y0, x1, y1, inclusive, which molecule_quad.vert covers with a quad):
//   sphere (16): box, cx, cy, sx, sy, rp, rp^2, Z, r, g, b, owner, 0
//   stick  (32): box, 4 kinds, 4 b, 4 c (capsule bound), Ax, Ay, Az, ux, uy, uz, a, R, L, r, g, b,
//                owner, 0, 0, 0
//   mesh triangle (28): box, 3 kinds, 3 b, 3 c, a0, b0, c0, a1, b1, c1, za, zb, zc, zlo, zhi, v0, v1,
//                v2, 0 (planes relative to the box's corner); vertex (8): normal, owner, r, g, b, 0
// Spheres are written in the software renderer's fragment order (radius, sub-pixel phase, index),
// items in their order: with the depth test GREATER_OR_EQUAL, a tie goes to the later primitive as
// the software renderer's depth keys break it. The occlusion culling of the software renderer is
// not needed (the depth test hides what it would skip).

#include "mdbridge/render/vulkan/scene_renderer.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <variant>

#include "mdbridge/numerics/ndimage.hpp"
#include "mdbridge/render/vulkan/pipeline.hpp"
#include "shaders/spirv.hpp"

namespace mdbridge::render::vulkan {

using molecule::Colored;
using molecule::Mesh;
using molecule::Spheres;
using molecule::Sticks;
using molecule::Vec3;

namespace {

// the software renderer's constants (software.cpp)
constexpr double kGlowRadius = 2.0, kGlowDisc = 2.5;
constexpr double kHuge = 1099511627776.0;  // 2^40 px: beyond it a primitive is skipped (RN-28)

constexpr std::uint32_t kSphereWords = 16, kStickWords = 32, kTriWords = 28, kVertexWords = 8;
constexpr std::uint32_t kFlagFog = 1, kFlagAa = 2;
constexpr std::uint32_t kStateWords = 16;

float f32(double v) { return static_cast<float>(v); }

template <class T>
T clip(T x, T lo, T hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}
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
std::int64_t toIndex(double v) {
    constexpr double lim = 4611686018427387904.0;  // 2^62
    if (!(v == v)) {
        return 0;
    }
    return static_cast<std::int64_t>(std::clamp(v, -lim, lim));
}
bool finite3(const Vec3& v) { return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]); }

// ---- camera (software.cpp Cam / makeCam) --------------------------------------------------------

struct Cam {
    std::array<Vec3, 3> R{};
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    std::int64_t W = 0, H = 0;
    double s = 0, sy = 0;
    double zmin = 0, zmax = 0;
    double fog_lo = 0, fog_hi = 0;
};

void checkViewport(const molecule::Viewport& vp) {
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

Cam makeCam(const molecule::Viewport& vp, const std::vector<Colored>& items, bool fast) {
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
    } else {  // RN-27: half the size rounded up, the limits widened
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
    std::tie(c.zmin, c.zmax) = molecule::depthRange(items, c.R[2]);
    c.fog_lo = vp.zrange[0];
    c.fog_hi = vp.zrange[1];
    return c;
}

struct ScreenPoint {
    float x, y, z;
};

ScreenPoint screen(const Cam& c, const Vec3& p) {
    const double v0 = p[0] * c.R[0][0] + p[1] * c.R[0][1] + p[2] * c.R[0][2];
    const double v1 = p[0] * c.R[1][0] + p[1] * c.R[1][1] + p[2] * c.R[1][2];
    const double v2 = p[0] * c.R[2][0] + p[1] * c.R[2][1] + p[2] * c.R[2][2];
    return {f32((v0 - c.xlo) * c.s), f32((c.yhi - v1) * c.sy), f32(v2)};
}

bool placeable(const ScreenPoint& p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) && std::fabs(p.x) < kHuge &&
           std::fabs(p.y) < kHuge;
}

Vec3 rotate(const Cam& c, const Vec3& n) {
    return {n[0] * c.R[0][0] + n[1] * c.R[0][1] + n[2] * c.R[0][2],
            n[0] * c.R[1][0] + n[1] * c.R[1][1] + n[2] * c.R[1][2],
            n[0] * c.R[2][0] + n[1] * c.R[2][1] + n[2] * c.R[2][2]};
}

// half-planes a*x + b*y + c >= 0 as the software renderer stores them (halfPlanes)
template <class T, std::size_t K>
struct HalfPlanes {
    std::array<std::int32_t, K> kind{};
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

void checkSize(std::size_t got, std::size_t want, const char* what) {
    if (got != want) {
        throw std::invalid_argument(std::string("3D item: ") + what + " count does not match");
    }
}

// ---- the frame's records ------------------------------------------------------------------------

enum class Kind { Mesh, Sphere, Stick };

struct Draw {
    Kind kind;
    std::uint32_t base = 0;   // first record (words)
    std::uint32_t vbase = 0;  // first vertex record (meshes)
    std::uint32_t count = 0;  // records
    float opacity = 1.0f;
};

// push constants of every pipeline (molecule_common.glsl, block Push)
struct Push {
    std::uint32_t W = 0, H = 0, base = 0, vbase = 0;
    float zmin = 0, zspan = 0, opacity = 1, hHalf = 0;
    float s = 0, fogLo = 0, fogSpan = 1, fogStrength = 0;
    float fogR = 0, fogG = 0, fogB = 0;
    float glowR = 0, glowG = 0, glowB = 0;
    std::uint32_t nGlow = 0;
    std::int32_t pad = 0;
    std::int32_t radius = 0;
    std::uint32_t flags = 0;
    float norm = 1, outline = 0, zFar = 0;
    std::uint32_t cw = 1, ch = 1;
    std::uint32_t stride = 0;
};
static_assert(sizeof(Push) == 112, "the push constant block of molecule_common.glsl");

struct Plan {
    bool fast = false;
    int outWidth = 0, outHeight = 0;  // the viewport's size
    std::uint32_t W = 0, H = 0;       // the drawn image's size
    std::vector<std::uint32_t> words;
    std::vector<Draw> opaque, glass;
    std::vector<float> glow;
    std::vector<float> weights;
    Push push;
};

class Writer {
public:
    explicit Writer(std::vector<std::uint32_t>& w) : m_w(w) {}
    void u(std::int64_t v) { m_w.push_back(static_cast<std::uint32_t>(static_cast<std::int32_t>(v))); }
    void f(float v) {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &v, 4);
        m_w.push_back(bits);
    }
    void zeros(std::size_t n) { m_w.insert(m_w.end(), n, 0u); }
    std::uint32_t size() const { return static_cast<std::uint32_t>(m_w.size()); }

private:
    std::vector<std::uint32_t>& m_w;
};

void writeRgb(Writer& w, const Vec3& rgb) {
    w.f(f32(rgb[0]));
    w.f(f32(rgb[1]));
    w.f(f32(rgb[2]));
}

// software.cpp spherePass (without the culling)
std::uint32_t addSpheres(Writer& w, const Cam& c, const Spheres& s, const std::vector<Vec3>& rgb) {
    const std::size_t n = s.xyz.size();
    checkSize(s.radius.size(), n, "radius");
    checkSize(s.owner.size(), n, "owner");
    checkSize(rgb.size(), n, "colour");
    struct Entry {
        ScreenPoint P;
        float rad;
        int q;
        std::size_t i;
    };
    std::vector<Entry> order;
    for (std::size_t i = 0; i < n; ++i) {
        if (!finite3(s.xyz[i]) || !std::isfinite(s.radius[i])) {
            continue;
        }
        const ScreenPoint p = screen(c, s.xyz[i]);
        const float r = f32(s.radius[i]);
        if (!placeable(p) || !(std::fabs(static_cast<double>(r) * c.s) < kHuge)) {
            continue;
        }
        const double cx = std::floor(static_cast<double>(p.x));
        const double cy = std::floor(static_cast<double>(p.y));
        const int qx = std::min(static_cast<int>((static_cast<double>(p.x) - cx) * 2), 1);
        const int qy = std::min(static_cast<int>((static_cast<double>(p.y) - cy) * 2), 1);
        order.push_back({p, r, qx * 2 + qy, i});
    }
    std::sort(order.begin(), order.end(), [](const Entry& a, const Entry& b) {
        if (a.rad != b.rad) {
            return a.rad < b.rad;
        }
        if (a.q != b.q) {
            return a.q < b.q;
        }
        return a.i < b.i;
    });
    std::uint32_t count = 0;
    for (const Entry& e : order) {
        const std::int64_t cx = toIndex(std::floor(static_cast<double>(e.P.x)));
        const std::int64_t cy = toIndex(std::floor(static_cast<double>(e.P.y)));
        const double rp = static_cast<double>(e.rad) * c.s;
        const std::int64_t nn = toIndex(std::ceil(rp)) + 1;
        if (nn < 0) {
            continue;
        }
        const std::int64_t x0 = std::max<std::int64_t>(cx - nn, 0), x1 = std::min(cx + nn, c.W - 1);
        const std::int64_t y0 = std::max<std::int64_t>(cy - nn, 0), y1 = std::min(cy + nn, c.H - 1);
        if (x1 < x0 || y1 < y0) {
            continue;
        }
        w.u(x0);
        w.u(y0);
        w.u(x1);
        w.u(y1);
        w.u(cx);
        w.u(cy);
        w.f(e.q < 2 ? 0.25f : 0.75f);
        w.f(e.q % 2 == 0 ? 0.25f : 0.75f);
        w.f(f32(rp));
        w.f(f32(rp * rp));
        w.f(e.P.z);
        writeRgb(w, rgb[e.i]);
        w.u(s.owner[e.i]);
        w.zeros(1);
        ++count;
    }
    return count;
}

// software.cpp stickPass
std::uint32_t addSticks(Writer& w, const Cam& c, const Sticks& st, const std::vector<Vec3>& rgb) {
    const std::size_t n = st.a.size();
    checkSize(st.b.size(), n, "stick end");
    checkSize(st.radius.size(), n, "radius");
    checkSize(st.owner.size(), n, "owner");
    checkSize(rgb.size(), n, "colour");
    const float s = f32(c.s);
    const float Hf = static_cast<float>(c.H);
    const float H1 = static_cast<float>(c.H - 1);
    const float W1 = static_cast<float>(c.W - 1);
    std::uint32_t count = 0;
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
        const float Ax = ax, Ay = Hf - ay, Az = az * s;
        const float Dx = bx - ax, Dy = ay - by, Dz = (bz - az) * s;
        const float L = std::sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
        const float R = f32(st.radius[i]) * s;
        const bool ok = (L > 1e-3f) && (L * L - Dz * Dz > 1e-4f);
        const float inv = 1.0f / npMax(L, 1e-3f);
        const float ux = Dx * inv, uy = Dy * inv, uz = Dz * inv;
        const float a = npMax(1.0f - uz * uz, 1e-6f);
        const float fy0 = clip(std::floor(npMin(ay, by) - R), 0.0f, H1);
        const float fy1 = clip(std::ceil(npMax(ay, by) + R), 0.0f, H1);
        if (!(ok && fy1 >= fy0)) {
            continue;
        }
        const float ddx = bx - ax, ddy = by - ay;
        const float L2 = std::sqrt(ddx * ddx + ddy * ddy);
        const float dx2 = ddx / npMax(L2, 1e-6f), dy2 = ddy / npMax(L2, 1e-6f);
        const float nx2 = -dy2, ny2 = dx2;
        const float Rs = R + 0.5f;
        const std::array<std::array<float, 3>, 4> G{{{nx2, ny2, Rs - nx2 * ax - ny2 * ay},
                                                     {-nx2, -ny2, Rs + nx2 * ax + ny2 * ay},
                                                     {dx2, dy2, Rs - dx2 * ax - dy2 * ay},
                                                     {-dx2, -dy2, L2 + Rs + dx2 * ax + dy2 * ay}}};
        const HalfPlanes<float, 4> h = halfPlanes(G);
        // columns: the cylinder stays within R of its axis on the screen
        const float fx0 = clip(std::floor(npMin(ax, bx) - Rs) - 1.0f, 0.0f, W1);
        const float fx1 = clip(std::ceil(npMax(ax, bx) + Rs) + 1.0f, 0.0f, W1);
        if (!(fx1 >= fx0)) {
            continue;
        }
        w.u(toIndex(fx0));
        w.u(toIndex(fy0));
        w.u(toIndex(fx1));
        w.u(toIndex(fy1));
        for (std::int32_t k : h.kind) {
            w.u(k);
        }
        for (float v : h.b) {
            w.f(v);
        }
        for (float v : h.c) {
            w.f(v);
        }
        for (float v : {Ax, Ay, Az, ux, uy, uz, a, R, L}) {
            w.f(v);
        }
        writeRgb(w, rgb[i]);
        w.u(st.owner[i]);
        w.zeros(3);
        ++count;
    }
    return count;
}

// software.cpp meshPass: vertex records at `vbase`, then the triangles
std::uint32_t addMesh(Writer& w, const Cam& c, const Mesh& m, const std::vector<Vec3>& rgb, std::uint32_t* vbase) {
    const std::size_t nv = m.verts.size();
    checkSize(m.normals.size(), nv, "normal");
    checkSize(m.owner.size(), nv, "owner");
    checkSize(rgb.size(), nv, "colour");
    for (const molecule::Tri& t : m.tris) {
        for (std::int32_t v : t) {
            if (v < 0 || static_cast<std::size_t>(v) >= nv) {
                throw std::invalid_argument("3D item: mesh triangle index out of range");
            }
        }
    }
    std::vector<ScreenPoint> sp(nv);
    std::vector<float> nz(nv);
    std::vector<bool> bad(nv);
    *vbase = w.size();
    for (std::size_t i = 0; i < nv; ++i) {
        sp[i] = screen(c, m.verts[i]);
        bad[i] = !finite3(m.verts[i]) || !finite3(m.normals[i]) || !placeable(sp[i]);
        const Vec3 n = rotate(c, m.normals[i]);
        nz[i] = f32(n[2]);
        w.f(f32(n[0]));
        w.f(f32(n[1]));
        w.f(f32(n[2]));
        w.u(m.owner[i]);
        writeRgb(w, rgb[i]);
        w.zeros(1);
    }
    const double W1 = static_cast<double>(c.W - 1);
    const double H1 = static_cast<double>(c.H - 1);
    std::uint32_t count = 0;
    for (const molecule::Tri& t : m.tris) {
        const auto i0 = static_cast<std::size_t>(t[0]);
        const auto i1 = static_cast<std::size_t>(t[1]);
        const auto i2 = static_cast<std::size_t>(t[2]);
        if (bad[i0] || bad[i1] || bad[i2]) {
            continue;
        }
        const float back = -0.05f;
        if (!(nz[i0] > back || nz[i1] > back || nz[i2] > back)) {
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
        const double inv = 1.0 / A;
        const double a0 = (y1 - y2) * inv, b0 = (x2 - x1) * inv, c0 = (x1 * y2 - x2 * y1) * inv;
        const double a1 = (y2 - y0) * inv, b1 = (x0 - x2) * inv, c1 = (x2 * y0 - x0 * y2) * inv;
        std::array<std::array<double, 3>, 3> G{{{a0, b0, c0}, {a1, b1, c1}, {-a0 - a1, -b0 - b1, 1 - c0 - c1}}};
        for (auto& g : G) {
            g[2] += 0.02 * std::hypot(g[0], g[1]);  // slightly conservative: no cracks
            g[2] += g[0] * bx0 + g[1] * by0;        // relative to the box's corner
        }
        const HalfPlanes<double, 3> h = halfPlanes(G);
        const double z0 = sp[i0].z, z1 = sp[i1].z, z2 = sp[i2].z;
        const double za = a0 * (z0 - z2) + a1 * (z1 - z2);
        const double zb = b0 * (z0 - z2) + b1 * (z1 - z2);
        const double zc = z2 + c0 * (z0 - z2) + c1 * (z1 - z2);
        w.u(toIndex(bx0));
        w.u(toIndex(by0));
        w.u(toIndex(bx1));
        w.u(toIndex(by1));
        for (std::int32_t k : h.kind) {
            w.u(k);
        }
        for (double v : h.b) {
            w.f(f32(v));
        }
        for (double v : h.c) {
            w.f(f32(v));
        }
        w.f(f32(a0));
        w.f(f32(b0));
        w.f(f32(c0 + a0 * bx0 + b0 * by0));
        w.f(f32(a1));
        w.f(f32(b1));
        w.f(f32(c1 + a1 * bx0 + b1 * by0));
        w.f(f32(za));
        w.f(f32(zb));
        w.f(f32(zc + za * bx0 + zb * by0));
        w.f(f32(npMin(npMin(z0, z1), z2)));
        w.f(f32(npMax(npMax(z0, z1), z2)));
        w.u(t[0]);
        w.u(t[1]);
        w.u(t[2]);
        w.zeros(1);
        ++count;
    }
    return count;
}

std::optional<Draw> addItem(std::vector<std::uint32_t>& words, const Cam& c, const Colored& it) {
    Writer w(words);
    Draw d{};
    d.opacity = f32(it.opacity);
    const std::size_t start = words.size();
    if (const auto* m = std::get_if<Mesh>(&it.prim)) {
        d.kind = Kind::Mesh;
        std::uint32_t vbase = 0;
        d.count = addMesh(w, c, *m, it.rgb, &vbase);  // the triangles follow the vertex block
        d.vbase = vbase;
        d.base = vbase + static_cast<std::uint32_t>(m->verts.size()) * kVertexWords;
    } else if (const auto* s = std::get_if<Spheres>(&it.prim)) {
        d.kind = Kind::Sphere;
        d.base = w.size();
        d.count = addSpheres(w, c, *s, it.rgb);
    } else if (const auto* k = std::get_if<Sticks>(&it.prim)) {
        d.kind = Kind::Stick;
        d.base = w.size();
        d.count = addSticks(w, c, *k, it.rgb);
    }
    if (d.count == 0) {
        words.resize(start);
        return std::nullopt;
    }
    return d;
}

std::vector<float> glowValues(const molecule::Shading& shading) {
    std::vector<float> g;
    if (shading.glow) {
        g.reserve(shading.glow->size());
        for (double v : *shading.glow) {
            g.push_back(std::isfinite(v) ? f32(v) : 0.0f);  // RN-28: non-finite glow counts as none
        }
    }
    return g;
}

// software.cpp discNorm: 1 / peak of a blurred disc of radius `rad` (RN-29: for the exact sizes)
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
    const numerics::Shape3 shape{1, side, side};
    numerics::gaussianFilter1d(d, shape, 1, sig, 2.5);
    numerics::gaussianFilter1d(d, shape, 2, sig, 2.5);
    float peak = -std::numeric_limits<float>::infinity();
    for (float v : d) {
        peak = npMax(peak, v);
    }
    return 1.0 / std::max(static_cast<double>(peak), 1e-6);
}

// `all`: one pass with every item (resolve); else opaque and translucent items apart
Plan makePlan(const std::vector<Colored>& items, const molecule::Viewport& vp, const molecule::Shading& shading,
              bool fast, bool all) {
    const Cam c = makeCam(vp, items, fast);
    if (c.W > 65536 || c.H > 65536) {
        throw std::runtime_error("3D view of " + std::to_string(c.W) + " x " + std::to_string(c.H) +
                                 " px is too large for the Vulkan renderer");
    }
    Plan p;
    p.fast = fast;
    p.outWidth = vp.width;
    p.outHeight = vp.height;
    p.W = static_cast<std::uint32_t>(c.W);
    p.H = static_cast<std::uint32_t>(c.H);
    for (const Colored& it : items) {
        if (std::optional<Draw> d = addItem(p.words, c, it)) {
            (all || it.opacity >= 0.999 ? p.opaque : p.glass).push_back(*d);
        }
    }
    if (p.words.empty()) {
        p.words.push_back(0);  // buffers are never empty
    }
    p.glow = glowValues(shading);
    const bool glows = std::any_of(p.glow.begin(), p.glow.end(), [](float g) { return g > 0.01f; });
    Push& u = p.push;
    u.W = p.W;
    u.H = p.H;
    u.zmin = f32(c.zmin);
    u.zspan = f32(c.zmax - c.zmin);
    u.hHalf = f32(static_cast<double>(c.H) - 0.5);
    u.s = f32(c.s);
    u.fogLo = f32(c.fog_lo);
    u.fogSpan = f32(std::max(c.fog_hi - c.fog_lo, 1e-6));
    u.fogStrength = f32(shading.fog_strength);
    u.fogR = f32(shading.fog_rgb[0]);
    u.fogG = f32(shading.fog_rgb[1]);
    u.fogB = f32(shading.fog_rgb[2]);
    u.glowR = f32(shading.glow_rgb[0]);
    u.glowG = f32(shading.glow_rgb[1]);
    u.glowB = f32(shading.glow_rgb[2]);
    u.nGlow = static_cast<std::uint32_t>(p.glow.size());
    u.flags = (shading.fog_strength > 0 ? kFlagFog : 0u) | (fast ? 0u : kFlagAa);
    u.outline = fast ? 0.0f : 0.5f;
    u.zFar = f32(c.zmin - 1e3);
    u.cw = (p.W + 3) / 4 + 1;
    u.ch = (p.H + 3) / 4 + 1;
    u.pad = static_cast<std::int32_t>(std::min<std::int64_t>(toIndex(std::ceil(2.5 * kGlowRadius * c.s)), 1 << 20));
    if (glows) {  // the halo's Gaussian (software.cpp glowHalo)
        const double sig = kGlowRadius * c.s / 4.0;
        const std::size_t radius = static_cast<std::size_t>(std::max(0.0, std::floor(2.5 * sig + 0.5)));
        if (radius > 4096) {
            throw std::runtime_error("3D view zoomed in too far for the Vulkan renderer's glow");
        }
        for (double v : numerics::gaussianKernel1d(sig, radius)) {
            p.weights.push_back(f32(v));
        }
        u.radius = static_cast<std::int32_t>(radius);
        u.norm = f32(discNorm(sig, kGlowDisc * c.s / 4.0));
    } else {
        p.weights.push_back(1.0f);
    }
    if (p.glow.empty()) {
        p.glow.push_back(0.0f);  // buffers are never empty (nGlow stays 0)
    }
    return p;
}

// ---- Vulkan resources -----------------------------------------------------------------------------

[[noreturn]] void fail(const std::string& what, VkResult r) {
    throw std::runtime_error(what + " (" + resultText(r) + ")");
}

struct Buffer {
    VkBuffer buffer = VK_NULL_HANDLE;
    VkDeviceMemory memory = VK_NULL_HANDLE;
    void* mapped = nullptr;
    VkDeviceSize size = 0;
    bool coherent = true;
};

struct ImageRes {
    VkImage image = VK_NULL_HANDLE;
    VkDeviceMemory memory = VK_NULL_HANDLE;
    VkImageView view = VK_NULL_HANDLE;
};

VkImageMemoryBarrier2 imageBarrier(VkImage image, VkImageAspectFlags aspect, VkImageLayout from, VkImageLayout to,
                                   VkPipelineStageFlags2 srcStage, VkAccessFlags2 srcAccess,
                                   VkPipelineStageFlags2 dstStage, VkAccessFlags2 dstAccess) {
    VkImageMemoryBarrier2 b{};
    b.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
    b.srcStageMask = srcStage;
    b.srcAccessMask = srcAccess;
    b.dstStageMask = dstStage;
    b.dstAccessMask = dstAccess;
    b.oldLayout = from;
    b.newLayout = to;
    b.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.image = image;
    b.subresourceRange = {aspect, 0, 1, 0, 1};
    return b;
}

std::uint32_t groups(std::uint32_t n, std::uint32_t size) { return (n + size - 1) / size; }

}  // namespace

// ---- the renderer's state -----------------------------------------------------------------------

struct SceneRenderer::Impl {
    Context& ctx;
    const DeviceFunctions& vk;
    VkDevice device;
    VkDescriptorSetLayout setLayout = VK_NULL_HANDLE;
    VkDescriptorPool pool = VK_NULL_HANDLE;
    std::unique_ptr<Pipeline> mesh, sphere, stick, present;
    std::unique_ptr<ComputePipeline> bounds, coarse, blurY, blurX, composite, finish;
    VkFormat depthFormat = VK_FORMAT_UNDEFINED;

    struct Slot {
        Buffer prims, glow, weights;
        VkDescriptorSet set = VK_NULL_HANDLE;
    };
    std::vector<Slot> slots;

    // size-dependent
    std::uint32_t W = 0, H = 0, cw = 0, ch = 0;
    ImageRes gOpaque, gGlass, depth;
    Buffer coarseA, coarseB, comp, out, state, readback;

    // the last recorded frame (recordPresent)
    std::uint32_t lastSlot = 0, lastW = 0, lastH = 0;
    bool lastFast = false, haveFrame = false;

    explicit Impl(Context& c) : ctx(c), vk(c.vk), device(c.device()) {}
    ~Impl();

    void createBuffer(Buffer& b, VkDeviceSize size, VkBufferUsageFlags usage, bool host);
    void destroyBuffer(Buffer& b);
    void createImage(ImageRes& im, VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect);
    void destroyImage(ImageRes& im);
    void init(std::uint32_t frameSlots);
    void ensureSize(std::uint32_t width, std::uint32_t height, bool waitIdle);
    void ensureReadback(VkDeviceSize bytes);
    void upload(Slot& slot, const Plan& plan);
    void writeDescriptors(Slot& slot);
    void record(VkCommandBuffer cmd, Slot& slot, const Plan& plan, bool resolveOnly);
    void geometryPass(VkCommandBuffer cmd, const ImageRes& target, const std::vector<Draw>& draws, Push push);
    void barrier(VkCommandBuffer cmd, VkPipelineStageFlags2 src, VkAccessFlags2 srcAccess, VkPipelineStageFlags2 dst,
                 VkAccessFlags2 dstAccess) const;
    void dispatch(VkCommandBuffer cmd, const ComputePipeline& p, const Push& push, std::uint32_t gx,
                  std::uint32_t gy) const;
    void invalidate(const Buffer& b) const;
};

void SceneRenderer::Impl::createBuffer(Buffer& b, VkDeviceSize size, VkBufferUsageFlags usage, bool host) {
    destroyBuffer(b);
    VkBufferCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    info.size = size;
    info.usage = usage;
    info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    if (const VkResult r = vk.vkCreateBuffer(device, &info, nullptr, &b.buffer); r != VK_SUCCESS) {
        b.buffer = VK_NULL_HANDLE;
        fail("could not create a Vulkan buffer of " + std::to_string(size) + " bytes", r);
    }
    VkMemoryRequirements needs{};
    vk.vkGetBufferMemoryRequirements(device, b.buffer, &needs);
    const std::optional<std::uint32_t> type =
        host ? ctx.memoryType(needs.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT,
                              VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)
             : ctx.memoryType(needs.memoryTypeBits, 0, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    if (!type) {
        throw std::runtime_error("no Vulkan memory type for the scene renderer's buffers");
    }
    VkPhysicalDeviceMemoryProperties memory{};
    ctx.vki.vkGetPhysicalDeviceMemoryProperties(ctx.physicalDevice(), &memory);
    b.coherent = (memory.memoryTypes[*type].propertyFlags & VK_MEMORY_PROPERTY_HOST_COHERENT_BIT) != 0;
    VkMemoryAllocateInfo allocate{};
    allocate.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocate.allocationSize = needs.size;
    allocate.memoryTypeIndex = *type;
    if (const VkResult r = vk.vkAllocateMemory(device, &allocate, nullptr, &b.memory); r != VK_SUCCESS) {
        b.memory = VK_NULL_HANDLE;
        fail("could not allocate " + std::to_string(needs.size) + " bytes of Vulkan memory", r);
    }
    vk.vkBindBufferMemory(device, b.buffer, b.memory, 0);
    if (host) {
        if (const VkResult r = vk.vkMapMemory(device, b.memory, 0, VK_WHOLE_SIZE, 0, &b.mapped); r != VK_SUCCESS) {
            b.mapped = nullptr;
            fail("could not map Vulkan memory", r);
        }
    }
    b.size = size;
}

void SceneRenderer::Impl::destroyBuffer(Buffer& b) {
    if (b.mapped) {
        vk.vkUnmapMemory(device, b.memory);
    }
    if (b.buffer) {
        vk.vkDestroyBuffer(device, b.buffer, nullptr);
    }
    if (b.memory) {
        vk.vkFreeMemory(device, b.memory, nullptr);
    }
    b = Buffer{};
}

void SceneRenderer::Impl::createImage(ImageRes& im, VkFormat format, VkImageUsageFlags usage,
                                      VkImageAspectFlags aspect) {
    destroyImage(im);
    VkImageCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    info.imageType = VK_IMAGE_TYPE_2D;
    info.format = format;
    info.extent = {W, H, 1};
    info.mipLevels = 1;
    info.arrayLayers = 1;
    info.samples = VK_SAMPLE_COUNT_1_BIT;
    info.tiling = VK_IMAGE_TILING_OPTIMAL;
    info.usage = usage;
    info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    if (const VkResult r = vk.vkCreateImage(device, &info, nullptr, &im.image); r != VK_SUCCESS) {
        im.image = VK_NULL_HANDLE;
        fail("could not create a " + std::to_string(W) + " x " + std::to_string(H) + " Vulkan image", r);
    }
    VkMemoryRequirements needs{};
    vk.vkGetImageMemoryRequirements(device, im.image, &needs);
    const std::optional<std::uint32_t> type =
        ctx.memoryType(needs.memoryTypeBits, 0, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    if (!type) {
        throw std::runtime_error("no Vulkan memory type for the scene renderer's images");
    }
    VkMemoryAllocateInfo allocate{};
    allocate.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocate.allocationSize = needs.size;
    allocate.memoryTypeIndex = *type;
    if (const VkResult r = vk.vkAllocateMemory(device, &allocate, nullptr, &im.memory); r != VK_SUCCESS) {
        im.memory = VK_NULL_HANDLE;
        fail("could not allocate " + std::to_string(needs.size) + " bytes for a Vulkan image", r);
    }
    vk.vkBindImageMemory(device, im.image, im.memory, 0);
    VkImageViewCreateInfo view{};
    view.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    view.image = im.image;
    view.viewType = VK_IMAGE_VIEW_TYPE_2D;
    view.format = format;
    view.subresourceRange = {aspect, 0, 1, 0, 1};
    if (const VkResult r = vk.vkCreateImageView(device, &view, nullptr, &im.view); r != VK_SUCCESS) {
        im.view = VK_NULL_HANDLE;
        fail("could not create a Vulkan image view", r);
    }
}

void SceneRenderer::Impl::destroyImage(ImageRes& im) {
    if (im.view) {
        vk.vkDestroyImageView(device, im.view, nullptr);
    }
    if (im.image) {
        vk.vkDestroyImage(device, im.image, nullptr);
    }
    if (im.memory) {
        vk.vkFreeMemory(device, im.memory, nullptr);
    }
    im = ImageRes{};
}

SceneRenderer::Impl::~Impl() {
    if (device) {
        vk.vkDeviceWaitIdle(device);
    }
    for (Slot& s : slots) {
        destroyBuffer(s.prims);
        destroyBuffer(s.glow);
        destroyBuffer(s.weights);
    }
    for (ImageRes* im : {&gOpaque, &gGlass, &depth}) {
        destroyImage(*im);
    }
    for (Buffer* b : {&coarseA, &coarseB, &comp, &out, &state, &readback}) {
        destroyBuffer(*b);
    }
    mesh.reset();
    sphere.reset();
    stick.reset();
    present.reset();
    bounds.reset();
    coarse.reset();
    blurY.reset();
    blurX.reset();
    composite.reset();
    finish.reset();
    if (pool) {
        vk.vkDestroyDescriptorPool(device, pool, nullptr);
    }
    if (setLayout) {
        vk.vkDestroyDescriptorSetLayout(device, setLayout, nullptr);
    }
}

void SceneRenderer::Impl::init(std::uint32_t frameSlots) {
    // one descriptor set layout for every pass: 8 storage buffers, 2 storage images
    std::array<VkDescriptorSetLayoutBinding, 10> bindings{};
    for (std::uint32_t i = 0; i < bindings.size(); ++i) {
        bindings[i].binding = i;
        bindings[i].descriptorType = i < 8 ? VK_DESCRIPTOR_TYPE_STORAGE_BUFFER : VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT | VK_SHADER_STAGE_COMPUTE_BIT;
    }
    VkDescriptorSetLayoutCreateInfo layoutInfo{};
    layoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    layoutInfo.bindingCount = static_cast<std::uint32_t>(bindings.size());
    layoutInfo.pBindings = bindings.data();
    if (const VkResult r = vk.vkCreateDescriptorSetLayout(device, &layoutInfo, nullptr, &setLayout); r != VK_SUCCESS) {
        setLayout = VK_NULL_HANDLE;
        fail("could not create a Vulkan descriptor set layout", r);
    }
    const std::uint32_t n = std::max<std::uint32_t>(frameSlots, 1);
    const VkDescriptorPoolSize sizes[2] = {{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 8 * n},
                                           {VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 2 * n}};
    VkDescriptorPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    poolInfo.maxSets = n;
    poolInfo.poolSizeCount = 2;
    poolInfo.pPoolSizes = sizes;
    if (const VkResult r = vk.vkCreateDescriptorPool(device, &poolInfo, nullptr, &pool); r != VK_SUCCESS) {
        pool = VK_NULL_HANDLE;
        fail("could not create a Vulkan descriptor pool", r);
    }
    slots.resize(n);
    for (Slot& s : slots) {
        VkDescriptorSetAllocateInfo allocate{};
        allocate.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        allocate.descriptorPool = pool;
        allocate.descriptorSetCount = 1;
        allocate.pSetLayouts = &setLayout;
        if (const VkResult r = vk.vkAllocateDescriptorSets(device, &allocate, &s.set); r != VK_SUCCESS) {
            fail("could not allocate a Vulkan descriptor set", r);
        }
    }

    for (VkFormat format : {VK_FORMAT_D32_SFLOAT, VK_FORMAT_X8_D24_UNORM_PACK32, VK_FORMAT_D16_UNORM}) {
        if (ctx.supportsFormat(format, VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT)) {
            depthFormat = format;
            break;
        }
    }
    if (depthFormat == VK_FORMAT_UNDEFINED) {
        throw std::runtime_error("the Vulkan device has no depth format");
    }
    const VkFormatFeatureFlags gbufferFeatures = VK_FORMAT_FEATURE_COLOR_ATTACHMENT_BIT |
                                                 VK_FORMAT_FEATURE_STORAGE_IMAGE_BIT | VK_FORMAT_FEATURE_TRANSFER_SRC_BIT;
    if (!ctx.supportsFormat(VK_FORMAT_R32G32B32A32_UINT, gbufferFeatures)) {
        throw std::runtime_error("the Vulkan device cannot render into RGBA32UI storage images");
    }

    std::string problem;
    auto graphics = [&](std::span<const std::uint32_t> frag) {
        GraphicsPipelineDesc d;
        d.stages = {{VK_SHADER_STAGE_VERTEX_BIT, spirv::kMoleculeQuadVertSpirv},
                    {VK_SHADER_STAGE_FRAGMENT_BIT, frag}};
        d.colorFormat = VK_FORMAT_R32G32B32A32_UINT;
        d.depthFormat = depthFormat;
        d.depthTest = true;
        d.depthWrite = true;
        d.depthCompare = VK_COMPARE_OP_GREATER_OR_EQUAL;
        d.pushConstantBytes = sizeof(Push);
        d.setLayouts = {setLayout};
        std::unique_ptr<Pipeline> p = Pipeline::create(ctx, d, &problem);
        if (!p) {
            throw std::runtime_error(problem);
        }
        return p;
    };
    mesh = graphics(spirv::kMoleculeMeshFragSpirv);
    sphere = graphics(spirv::kMoleculeSphereFragSpirv);
    stick = graphics(spirv::kMoleculeStickFragSpirv);
    auto compute = [&](std::span<const std::uint32_t> code) {
        ComputePipelineDesc d;
        d.spirv = code;
        d.pushConstantBytes = sizeof(Push);
        d.setLayouts = {setLayout};
        std::unique_ptr<ComputePipeline> p = ComputePipeline::create(ctx, d, &problem);
        if (!p) {
            throw std::runtime_error(problem);
        }
        return p;
    };
    bounds = compute(spirv::kMoleculeBoundsCompSpirv);
    coarse = compute(spirv::kMoleculeCoarseCompSpirv);
    blurY = compute(spirv::kMoleculeBlurYCompSpirv);
    blurX = compute(spirv::kMoleculeBlurXCompSpirv);
    composite = compute(spirv::kMoleculeCompositeCompSpirv);
    finish = compute(spirv::kMoleculeFinishCompSpirv);
    createBuffer(state, kStateWords * 4, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT, false);
    for (Slot& s : slots) {
        createBuffer(s.prims, 4096, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, true);
        createBuffer(s.glow, 4096, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, true);
        createBuffer(s.weights, 4096, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, true);
    }
    ensureSize(2, 2, false);
}

void SceneRenderer::Impl::ensureSize(std::uint32_t width, std::uint32_t height, bool waitIdle) {
    if (width == W && height == H) {
        return;
    }
    const std::uint32_t limit = ctx.properties().limits.maxImageDimension2D;
    if (width > limit || height > limit) {
        throw std::runtime_error("3D view of " + std::to_string(width) + " x " + std::to_string(height) +
                                 " px is larger than the Vulkan device's limit of " + std::to_string(limit) + " px");
    }
    const VkDeviceSize pixels = VkDeviceSize{width} * height;
    if (pixels * 4 > ctx.properties().limits.maxStorageBufferRange) {
        throw std::runtime_error("3D view of " + std::to_string(width) + " x " + std::to_string(height) +
                                 " px is too large for the Vulkan device's storage buffers");
    }
    if (waitIdle) {
        vk.vkDeviceWaitIdle(device);
    }
    W = width;
    H = height;
    cw = (W + 3) / 4 + 1;
    ch = (H + 3) / 4 + 1;
    const VkImageUsageFlags gbuffer =
        VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT;
    createImage(gOpaque, VK_FORMAT_R32G32B32A32_UINT, gbuffer, VK_IMAGE_ASPECT_COLOR_BIT);
    createImage(gGlass, VK_FORMAT_R32G32B32A32_UINT, gbuffer, VK_IMAGE_ASPECT_COLOR_BIT);
    createImage(depth, depthFormat, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, VK_IMAGE_ASPECT_DEPTH_BIT);
    const VkDeviceSize coarseBytes = VkDeviceSize{cw} * ch * 4;
    createBuffer(coarseA, coarseBytes, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, false);
    createBuffer(coarseB, coarseBytes, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, false);
    createBuffer(comp, pixels * 4, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, false);
    createBuffer(out, pixels * 4, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT, false);
    haveFrame = false;
    for (Slot& s : slots) {
        writeDescriptors(s);
    }
}

void SceneRenderer::Impl::ensureReadback(VkDeviceSize bytes) {
    if (readback.size < bytes) {
        createBuffer(readback, bytes, VK_BUFFER_USAGE_TRANSFER_DST_BIT, true);
    }
}

void SceneRenderer::Impl::writeDescriptors(Slot& s) {
    const Buffer* buffers[8] = {&s.prims, &s.glow, &s.weights, &state, &coarseA, &coarseB, &comp, &out};
    VkDescriptorBufferInfo bufferInfo[8]{};
    VkDescriptorImageInfo imageInfo[2]{};
    VkWriteDescriptorSet writes[10]{};
    for (std::uint32_t i = 0; i < 10; ++i) {
        writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[i].dstSet = s.set;
        writes[i].dstBinding = i;
        writes[i].descriptorCount = 1;
        if (i < 8) {
            bufferInfo[i] = {buffers[i]->buffer, 0, VK_WHOLE_SIZE};
            writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            writes[i].pBufferInfo = &bufferInfo[i];
        } else {
            imageInfo[i - 8] = {VK_NULL_HANDLE, (i == 8 ? gOpaque : gGlass).view, VK_IMAGE_LAYOUT_GENERAL};
            writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
            writes[i].pImageInfo = &imageInfo[i - 8];
        }
    }
    vk.vkUpdateDescriptorSets(device, 10, writes, 0, nullptr);
}

void SceneRenderer::Impl::upload(Slot& s, const Plan& plan) {
    const VkDeviceSize limit = ctx.properties().limits.maxStorageBufferRange;
    bool rewrite = false;
    auto put = [&](Buffer& b, const void* data, VkDeviceSize bytes) {
        if (bytes > limit) {
            throw std::runtime_error("scene too large for the Vulkan renderer (" + std::to_string(bytes) +
                                     " bytes of primitives; the device allows " + std::to_string(limit) + ")");
        }
        if (b.size < bytes) {
            createBuffer(b, std::max<VkDeviceSize>(bytes, b.size * 2), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, true);
            rewrite = true;
        }
        std::memcpy(b.mapped, data, static_cast<std::size_t>(bytes));
        if (!b.coherent) {
            VkMappedMemoryRange range{};
            range.sType = VK_STRUCTURE_TYPE_MAPPED_MEMORY_RANGE;
            range.memory = b.memory;
            range.size = VK_WHOLE_SIZE;
            vk.vkFlushMappedMemoryRanges(device, 1, &range);
        }
    };
    put(s.prims, plan.words.data(), plan.words.size() * 4);
    put(s.glow, plan.glow.data(), plan.glow.size() * 4);
    put(s.weights, plan.weights.data(), plan.weights.size() * 4);
    if (rewrite) {
        writeDescriptors(s);
    }
}

void SceneRenderer::Impl::invalidate(const Buffer& b) const {
    if (!b.coherent) {
        VkMappedMemoryRange range{};
        range.sType = VK_STRUCTURE_TYPE_MAPPED_MEMORY_RANGE;
        range.memory = b.memory;
        range.size = VK_WHOLE_SIZE;
        vk.vkInvalidateMappedMemoryRanges(device, 1, &range);
    }
}

void SceneRenderer::Impl::barrier(VkCommandBuffer cmd, VkPipelineStageFlags2 src, VkAccessFlags2 srcAccess,
                                  VkPipelineStageFlags2 dst, VkAccessFlags2 dstAccess) const {
    VkMemoryBarrier2 b{};
    b.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER_2;
    b.srcStageMask = src;
    b.srcAccessMask = srcAccess;
    b.dstStageMask = dst;
    b.dstAccessMask = dstAccess;
    VkDependencyInfo dep{};
    dep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
    dep.memoryBarrierCount = 1;
    dep.pMemoryBarriers = &b;
    vk.vkCmdPipelineBarrier2(cmd, &dep);
}

void SceneRenderer::Impl::dispatch(VkCommandBuffer cmd, const ComputePipeline& p, const Push& push, std::uint32_t gx,
                                   std::uint32_t gy) const {
    p.bind(cmd);
    p.push(cmd, &push, sizeof(Push));
    vk.vkCmdDispatch(cmd, gx, gy, 1);
    barrier(cmd, VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT,
            VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT,
            VK_ACCESS_2_SHADER_STORAGE_READ_BIT | VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT);
}

void SceneRenderer::Impl::geometryPass(VkCommandBuffer cmd, const ImageRes& target, const std::vector<Draw>& draws,
                                       Push push) {
    VkRenderingAttachmentInfo color{};
    color.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
    color.imageView = target.view;
    color.imageLayout = VK_IMAGE_LAYOUT_GENERAL;
    color.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    color.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    color.clearValue.color.uint32[0] = 0;
    VkRenderingAttachmentInfo depthInfo{};
    depthInfo.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
    depthInfo.imageView = depth.view;
    depthInfo.imageLayout = VK_IMAGE_LAYOUT_DEPTH_ATTACHMENT_OPTIMAL;
    depthInfo.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    depthInfo.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    depthInfo.clearValue.depthStencil = {0.0f, 0};
    VkRenderingInfo rendering{};
    rendering.sType = VK_STRUCTURE_TYPE_RENDERING_INFO;
    rendering.renderArea = {{0, 0}, {W, H}};
    rendering.layerCount = 1;
    rendering.colorAttachmentCount = 1;
    rendering.pColorAttachments = &color;
    rendering.pDepthAttachment = &depthInfo;
    vk.vkCmdBeginRendering(cmd, &rendering);
    const VkViewport viewport{0.0f, 0.0f, static_cast<float>(W), static_cast<float>(H), 0.0f, 1.0f};
    const VkRect2D scissor{{0, 0}, {W, H}};
    vk.vkCmdSetViewport(cmd, 0, 1, &viewport);
    vk.vkCmdSetScissor(cmd, 0, 1, &scissor);
    for (const Draw& d : draws) {
        const Pipeline& p = d.kind == Kind::Mesh ? *mesh : d.kind == Kind::Sphere ? *sphere : *stick;
        push.base = d.base;
        push.vbase = d.vbase;
        push.opacity = d.opacity;
        push.stride = d.kind == Kind::Mesh ? kTriWords : d.kind == Kind::Sphere ? kSphereWords : kStickWords;
        p.bind(cmd);
        p.push(cmd, &push, sizeof(Push));
        vk.vkCmdDraw(cmd, 6, d.count, 0, 0);
    }
    vk.vkCmdEndRendering(cmd);
}

void SceneRenderer::Impl::record(VkCommandBuffer cmd, Slot& s, const Plan& plan, bool resolveOnly) {
    // earlier submissions (an earlier frame of a window) are done with every shared resource, and the
    // G-buffers start undefined
    VkMemoryBarrier2 all{};
    all.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER_2;
    all.srcStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    all.srcAccessMask = VK_ACCESS_2_MEMORY_WRITE_BIT;
    all.dstStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    all.dstAccessMask = VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT;
    const VkPipelineStageFlags2 frag = VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT |
                                       VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT;
    VkImageMemoryBarrier2 layouts[3] = {
        imageBarrier(gOpaque.image, VK_IMAGE_ASPECT_COLOR_BIT, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL,
                     VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0, VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT,
                     VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT),
        imageBarrier(gGlass.image, VK_IMAGE_ASPECT_COLOR_BIT, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL,
                     VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0, VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT,
                     VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT),
        imageBarrier(depth.image, VK_IMAGE_ASPECT_DEPTH_BIT, VK_IMAGE_LAYOUT_UNDEFINED,
                     VK_IMAGE_LAYOUT_DEPTH_ATTACHMENT_OPTIMAL, VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0, frag,
                     VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_READ_BIT | VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT),
    };
    VkDependencyInfo dep{};
    dep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
    dep.memoryBarrierCount = 1;
    dep.pMemoryBarriers = &all;
    dep.imageMemoryBarrierCount = 3;
    dep.pImageMemoryBarriers = layouts;
    vk.vkCmdPipelineBarrier2(cmd, &dep);

    const VkPipelineBindPoint points[2] = {VK_PIPELINE_BIND_POINT_GRAPHICS, VK_PIPELINE_BIND_POINT_COMPUTE};
    const VkPipelineLayout layoutOf[2] = {sphere->layout(), bounds->layout()};
    for (int k = 0; k < 2; ++k) {
        vk.vkCmdBindDescriptorSets(cmd, points[k], layoutOf[k], 0, 1, &s.set, 0, nullptr);
    }

    Push push = plan.push;
    push.cw = cw;
    push.ch = ch;
    geometryPass(cmd, gOpaque, plan.opaque, push);
    if (resolveOnly) {
        barrier(cmd, VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT, VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT,
                VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);
        VkBufferImageCopy region{};
        region.imageSubresource = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1};
        region.imageExtent = {W, H, 1};
        vk.vkCmdCopyImageToBuffer(cmd, gOpaque.image, VK_IMAGE_LAYOUT_GENERAL, readback.buffer, 1, &region);
        barrier(cmd, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT, VK_PIPELINE_STAGE_2_HOST_BIT,
                VK_ACCESS_2_HOST_READ_BIT);
        return;
    }
    barrier(cmd, frag, VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT, frag,
            VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_READ_BIT | VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT);
    geometryPass(cmd, gGlass, plan.glass, push);

    // the compositing state: empty boxes, no glow
    std::int32_t init[kStateWords] = {};
    for (int k = 0; k < 8; ++k) {
        init[k] = (k % 4) < 2 ? std::numeric_limits<std::int32_t>::max() : -1;
    }
    vk.vkCmdUpdateBuffer(cmd, state.buffer, 0, sizeof init, init);
    VkMemoryBarrier2 toCompute[2]{};
    for (VkMemoryBarrier2& b : toCompute) {
        b.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER_2;
        b.dstStageMask = VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT;
        b.dstAccessMask = VK_ACCESS_2_SHADER_STORAGE_READ_BIT | VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT;
    }
    toCompute[0].srcStageMask = VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
    toCompute[0].srcAccessMask = VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT;
    toCompute[1].srcStageMask = VK_PIPELINE_STAGE_2_COPY_BIT;
    toCompute[1].srcAccessMask = VK_ACCESS_2_TRANSFER_WRITE_BIT;
    VkDependencyInfo computeDep{};
    computeDep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
    computeDep.memoryBarrierCount = 2;
    computeDep.pMemoryBarriers = toCompute;
    vk.vkCmdPipelineBarrier2(cmd, &computeDep);

    const std::uint32_t px = groups(W, 16), py = groups(H, 16);
    const std::uint32_t cx = groups(cw, 8), cy = groups(ch, 8);
    dispatch(cmd, *bounds, push, groups(groups(W, 16), 64), H);
    dispatch(cmd, *coarse, push, cx, cy);
    dispatch(cmd, *blurY, push, cx, cy);
    dispatch(cmd, *blurX, push, cx, cy);
    dispatch(cmd, *composite, push, px, py);
    finish->bind(cmd);
    finish->push(cmd, &push, sizeof(Push));
    vk.vkCmdDispatch(cmd, px, py, 1);
}

// ---- public -------------------------------------------------------------------------------------

std::unique_ptr<SceneRenderer> SceneRenderer::create(Context& context, std::string* problem,
                                                     const SceneRendererOptions& options) {
    std::unique_ptr<SceneRenderer> renderer(new SceneRenderer(context));
    renderer->m_impl = std::make_unique<Impl>(context);
    try {
        renderer->m_impl->init(options.frameSlots);
    } catch (const std::exception& e) {
        *problem = e.what();
        return nullptr;
    }
    problem->clear();
    return renderer;
}

SceneRenderer::SceneRenderer(Context& context) : m_context(context) {}

SceneRenderer::~SceneRenderer() = default;

namespace {

molecule::RgbaImage toImage(const std::uint32_t* pixels, std::uint32_t W, std::uint32_t H) {
    molecule::RgbaImage img;
    img.width = W;
    img.height = H;
    img.rgba.resize(std::size_t{W} * H * 4);
    std::memcpy(img.rgba.data(), pixels, img.rgba.size());  // RGBA8 per uint, little-endian order
    return img;
}

// "fast": each pixel 2 x 2, cut to the requested size (software.cpp render)
molecule::RgbaImage upscaled(const molecule::RgbaImage& half, int width, int height) {
    molecule::RgbaImage full;
    full.width = static_cast<std::size_t>(width);
    full.height = static_cast<std::size_t>(height);
    full.rgba.assign(full.width * full.height * 4, 0);
    for (std::size_t y = 0; y < full.height && y / 2 < half.height; ++y) {
        for (std::size_t x = 0; x < full.width && x / 2 < half.width; ++x) {
            std::copy_n(half.pixel(x / 2, y / 2), 4, full.pixel(x, y));
        }
    }
    return full;
}

}  // namespace

molecule::RgbaImage SceneRenderer::render(const std::vector<Colored>& items, const molecule::Viewport& viewport,
                                          const molecule::Shading& shading, molecule::RenderQuality quality) {
    const bool fast = quality == molecule::RenderQuality::Fast;
    const Plan plan = makePlan(items, viewport, shading, fast, false);
    molecule::RgbaImage img;
    {
        const std::lock_guard<std::mutex> lock(m_mutex);
        Impl& im = *m_impl;
        im.ensureSize(plan.W, plan.H, true);
        Impl::Slot& slot = im.slots[0];
        im.upload(slot, plan);
        im.ensureReadback(VkDeviceSize{plan.W} * plan.H * 4);
        const VkResult r = m_context.submitAndWait([&](VkCommandBuffer cmd) {
            im.record(cmd, slot, plan, false);
            im.barrier(cmd, VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT,
                       VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);
            const VkBufferCopy region{0, 0, VkDeviceSize{plan.W} * plan.H * 4};
            m_context.vk.vkCmdCopyBuffer(cmd, im.out.buffer, im.readback.buffer, 1, &region);
            im.barrier(cmd, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT, VK_PIPELINE_STAGE_2_HOST_BIT,
                       VK_ACCESS_2_HOST_READ_BIT);
        });
        if (r != VK_SUCCESS) {
            fail("Vulkan rendering of the 3D view failed", r);
        }
        im.invalidate(im.readback);
        img = toImage(static_cast<const std::uint32_t*>(im.readback.mapped), plan.W, plan.H);
        im.haveFrame = false;  // the output buffer is the read-back's, not a window's frame
    }
    if (fast) {
        img = upscaled(img, viewport.width, viewport.height);
    }
    return img;
}

molecule::ResolvedPixels SceneRenderer::resolve(const std::vector<Colored>& items, const molecule::Viewport& viewport) {
    const Plan plan = makePlan(items, viewport, {}, false, true);
    std::vector<std::uint32_t> texels;
    {
        const std::lock_guard<std::mutex> lock(m_mutex);
        Impl& im = *m_impl;
        im.ensureSize(plan.W, plan.H, true);
        Impl::Slot& slot = im.slots[0];
        im.upload(slot, plan);
        const VkDeviceSize bytes = VkDeviceSize{plan.W} * plan.H * 16;
        im.ensureReadback(bytes);
        const VkResult r = m_context.submitAndWait([&](VkCommandBuffer cmd) { im.record(cmd, slot, plan, true); });
        if (r != VK_SUCCESS) {
            fail("Vulkan rendering of the 3D view failed", r);
        }
        im.invalidate(im.readback);
        texels.resize(static_cast<std::size_t>(bytes / 4));
        std::memcpy(texels.data(), im.readback.mapped, static_cast<std::size_t>(bytes));
    }
    const Cam cam = makeCam(viewport, items, false);
    molecule::ResolvedPixels out;
    out.width = plan.W;
    out.height = plan.H;
    out.scale = cam.s;
    for (std::size_t i = 0; i < std::size_t{plan.W} * plan.H; ++i) {
        const std::uint32_t* t = texels.data() + 4 * i;
        if (t[3] == 0) {
            continue;
        }
        auto unorm = [](std::uint32_t v, int half) {
            return static_cast<float>((v >> (16 * half)) & 0xffffu) / 65535.0f;
        };
        float depth = 0;
        std::memcpy(&depth, &t[2], 4);
        out.pixel.push_back(static_cast<std::int64_t>(i));
        out.rgb.push_back({unorm(t[0], 0), unorm(t[0], 1), unorm(t[1], 0)});
        out.alpha.push_back(unorm(t[1], 1));
        out.depth.push_back(depth);
        out.owner.push_back(static_cast<std::int64_t>(t[3]) - 2);
    }
    return out;
}

void SceneRenderer::recordFrame(VkCommandBuffer cmd, std::uint32_t slot, const std::vector<Colored>& items,
                                const molecule::Viewport& viewport, const molecule::Shading& shading,
                                molecule::RenderQuality quality) {
    const bool fast = quality == molecule::RenderQuality::Fast;
    const Plan plan = makePlan(items, viewport, shading, fast, false);
    const std::lock_guard<std::mutex> lock(m_mutex);
    Impl& im = *m_impl;
    if (slot >= im.slots.size()) {
        throw std::invalid_argument("SceneRenderer::recordFrame: no upload slot " + std::to_string(slot));
    }
    im.ensureSize(plan.W, plan.H, true);
    im.upload(im.slots[slot], plan);
    im.record(cmd, im.slots[slot], plan, false);
    im.barrier(cmd, VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT,
               VK_PIPELINE_STAGE_2_FRAGMENT_SHADER_BIT, VK_ACCESS_2_SHADER_STORAGE_READ_BIT);
    im.lastSlot = slot;
    im.lastW = plan.W;
    im.lastH = plan.H;
    im.lastFast = fast;
    im.haveFrame = true;
}

bool SceneRenderer::createPresentPipeline(VkRenderPass renderPass, VkSampleCountFlagBits samples,
                                          std::string* problem) {
    const std::lock_guard<std::mutex> lock(m_mutex);
    GraphicsPipelineDesc d;
    d.stages = {{VK_SHADER_STAGE_VERTEX_BIT, spirv::kFullscreenVertSpirv},
                {VK_SHADER_STAGE_FRAGMENT_BIT, spirv::kMoleculePresentFragSpirv}};
    d.renderPass = renderPass;
    d.samples = samples;
    d.pushConstantBytes = sizeof(Push);
    d.setLayouts = {m_impl->setLayout};
    m_impl->present = Pipeline::create(m_context, d, problem);
    return m_impl->present != nullptr;
}

void SceneRenderer::recordPresent(VkCommandBuffer cmd, std::uint32_t targetWidth, std::uint32_t targetHeight,
                                  const std::array<float, 3>& background) const {
    const Impl& im = *m_impl;
    if (!im.present) {
        return;
    }
    Push push;
    push.W = im.haveFrame ? im.lastW : 0;
    push.H = im.haveFrame ? im.lastH : 0;
    push.base = im.lastFast ? 1 : 0;
    push.fogR = background[0];
    push.fogG = background[1];
    push.fogB = background[2];
    const VkViewport viewport{0.0f, 0.0f, static_cast<float>(targetWidth), static_cast<float>(targetHeight), 0.0f, 1.0f};
    const VkRect2D scissor{{0, 0}, {targetWidth, targetHeight}};
    m_context.vk.vkCmdSetViewport(cmd, 0, 1, &viewport);
    m_context.vk.vkCmdSetScissor(cmd, 0, 1, &scissor);
    im.present->bind(cmd);
    m_context.vk.vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, im.present->layout(), 0, 1,
                                         &im.slots[im.lastSlot].set, 0, nullptr);
    im.present->push(cmd, &push, sizeof(Push));
    m_context.vk.vkCmdDraw(cmd, 3, 1, 0, 0);
}

}  // namespace mdbridge::render::vulkan
