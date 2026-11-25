// The built-in software 3D renderer (src/render/molecule/software.cpp): counterparts of MDATK's
// tests/test_molecule_render.py (exact rasterisation, depth, compositing, threads, culling,
// translucency, non-finite input), MDBridge's improvements (RN-26 ... RN-29), and the RGBA of
// MDATK's renderer for the same items and viewports (data/render/molecule/renders.npz,
// tests/parity/make_molecule_cases.py): >= 99.5 % of pixels equal, all within 2/255, coverage IoU
// >= 0.999 (docs/specs/m4-plan.md, section 3).

#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <filesystem>
#include <fstream>
#include <map>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <cstdio>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/molecule/geometry.hpp"
#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/software.hpp"
#include "mdbridge/render/molecule/surface.hpp"
#include "mdbridge/render/orientation.hpp"

namespace mdbridge::render::molecule {
namespace {

using core::Array;
using core::Value;

// ---- data ------------------------------------------------------------------------------------------

std::filesystem::path dataDir() { return std::filesystem::path(MDBRIDGE_TEST_DATA) / "render" / "molecule"; }

const Value& cases() {
    static const Value doc = [] {
        std::ifstream in(dataDir() / "molecule_cases.json", std::ios::binary);
        if (!in) {
            throw std::runtime_error("cannot read molecule_cases.json");
        }
        std::stringstream text;
        text << in.rdbuf();
        return core::parseJson(text.str(), "molecule_cases.json");
    }();
    return doc;
}

// An array of a .npz archive written by render_common.write_npz (ZIP, members stored), widened to
// double: '<f4', '<f8', '<i4', '<i8', '|i1', '|u1', '|b1' in C order.
struct Arr {
    std::vector<std::size_t> shape;
    std::vector<double> values;
    std::vector<double> toDouble() const { return values; }
};

std::uint32_t le32(const unsigned char* p) {
    return std::uint32_t{p[0]} | (std::uint32_t{p[1]} << 8) | (std::uint32_t{p[2]} << 16) | (std::uint32_t{p[3]} << 24);
}
std::uint16_t le16(const unsigned char* p) { return static_cast<std::uint16_t>(p[0] | (p[1] << 8)); }

Arr parseNpy(const unsigned char* d, std::size_t n, const std::string& name) {
    if (n < 10 || d[0] != 0x93 || std::string(reinterpret_cast<const char*>(d) + 1, 5) != "NUMPY") {
        throw std::runtime_error(name + ": not a .npy");
    }
    const std::size_t hl = d[6] == 1 ? le16(d + 8) : le32(d + 8);
    const std::size_t off = d[6] == 1 ? 10 : 12;
    const std::string header(reinterpret_cast<const char*>(d) + off, hl);
    auto field = [&](const std::string& key) {
        const std::size_t k = header.find("'" + key + "'");
        if (k == std::string::npos) {
            throw std::runtime_error(name + ": no " + key);
        }
        return header.substr(header.find(':', k) + 1);
    };
    std::string descr = field("descr");
    descr = descr.substr(descr.find('\'') + 1);
    descr = descr.substr(0, descr.find('\''));
    if (field("fortran_order").find("False") != 1 && field("fortran_order").find("False") != 0) {
        throw std::runtime_error(name + ": Fortran order");
    }
    std::string sh = field("shape");
    sh = sh.substr(sh.find('(') + 1, sh.find(')') - sh.find('(') - 1);
    Arr a;
    std::size_t count = 1;
    std::stringstream ss(sh);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.find_first_not_of(' ') != std::string::npos) {
            a.shape.push_back(static_cast<std::size_t>(std::stoull(tok)));
            count *= a.shape.back();
        }
    }
    const unsigned char* p = d + off + hl;
    a.values.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        if (descr == "<f4") {
            float f;
            std::memcpy(&f, p + 4 * i, 4);
            a.values[i] = f;
        } else if (descr == "<f8") {
            std::memcpy(&a.values[i], p + 8 * i, 8);
        } else if (descr == "<i4") {
            std::int32_t v;
            std::memcpy(&v, p + 4 * i, 4);
            a.values[i] = v;
        } else if (descr == "<i8") {
            std::int64_t v;
            std::memcpy(&v, p + 8 * i, 8);
            a.values[i] = static_cast<double>(v);
        } else if (descr == "|i1") {
            a.values[i] = static_cast<std::int8_t>(p[i]);
        } else if (descr == "|u1" || descr == "|b1") {
            a.values[i] = p[i];
        } else {
            throw std::runtime_error(name + ": dtype " + descr);
        }
    }
    return a;
}

const std::map<std::string, Arr>& npzFile(const std::string& file) {
    static std::map<std::string, std::map<std::string, Arr>> cache;
    auto it = cache.find(file);
    if (it == cache.end()) {
        std::ifstream in(dataDir() / file, std::ios::binary);
        const std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::map<std::string, Arr> arrays;
        std::size_t pos = 0;
        while (pos + 30 <= bytes.size() && le32(bytes.data() + pos) == 0x04034b50) {
            const unsigned char* h = bytes.data() + pos;
            if (le16(h + 8) != 0) {
                throw std::runtime_error(file + ": compressed member");
            }
            const std::size_t size = le32(h + 18), nl = le16(h + 26), el = le16(h + 28);
            std::string name(reinterpret_cast<const char*>(h) + 30, nl);
            const unsigned char* data = h + 30 + nl + el;
            if (name.size() > 4 && name.substr(name.size() - 4) == ".npy") {
                name.resize(name.size() - 4);
            }
            arrays.emplace(name, parseNpy(data, size, file + "/" + name));
            pos += 30 + nl + el + size;
        }
        it = cache.emplace(file, std::move(arrays)).first;
    }
    return it->second;
}

const Arr& npz(const std::string& file, const std::string& key) {
    const auto& all = npzFile(file);
    const auto it = all.find(key);
    if (it == all.end()) {
        throw std::runtime_error(file + " has no " + key);
    }
    return it->second;
}

const Value& at(const Value& v, std::string_view key) {
    const Value* p = v.find(key);
    if (!p) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *p;
}
double number(const Value& v) {
    if (const auto* d = v.get<double>()) {
        return *d;
    }
    return static_cast<double>(*v.get<std::int64_t>());
}
std::vector<double> doubles(const Value& v) {
    std::vector<double> out;
    for (const Value& x : *v.get<Array>()) {
        out.push_back(number(x));
    }
    return out;
}
Vec3 vec3(const Value& v) {
    const std::vector<double> d = doubles(v);
    return {d.at(0), d.at(1), d.at(2)};
}

std::vector<Vec3> vec3s(const Arr& a) {
    const std::vector<double> d = a.toDouble();
    std::vector<Vec3> out(d.size() / 3);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
    }
    return out;
}
std::vector<std::int32_t> ints(const Arr& a) {
    std::vector<std::int32_t> out;
    for (double v : a.values) {
        out.push_back(static_cast<std::int32_t>(v));
    }
    return out;
}

// The items as make_molecule_cases.py stored them ("<prefix>/<k>/<field>", float32 widened to double).
std::vector<Colored> storedItems(const std::string& file, const std::string& prefix, const Value& records) {
    std::vector<Colored> items;
    std::size_t k = 0;
    for (const Value& rec : *records.get<Array>()) {
        const std::string p = prefix + "/" + std::to_string(k++) + "/";
        const std::string& kind = *at(rec, "kind").get<std::string>();
        Colored it;
        it.opacity = number(at(rec, "opacity"));
        it.rgb = vec3s(npz(file, p + "rgb"));
        if (kind == "mesh") {
            Mesh m;
            m.verts = vec3s(npz(file, p + "verts"));
            m.normals = vec3s(npz(file, p + "normals"));
            const std::vector<std::int32_t> t = ints(npz(file, p + "tris"));
            for (std::size_t i = 0; i + 2 < t.size(); i += 3) {
                m.tris.push_back({t[i], t[i + 1], t[i + 2]});
            }
            m.owner = ints(npz(file, p + "owner"));
            m.role.assign(m.verts.size(), 0);
            it.prim = std::move(m);
        } else if (kind == "spheres") {
            Spheres s;
            s.xyz = vec3s(npz(file, p + "xyz"));
            s.radius = npz(file, p + "radius").toDouble();
            s.owner = ints(npz(file, p + "owner"));
            s.atom = ints(npz(file, p + "atom"));
            it.prim = std::move(s);
        } else {
            Sticks s;
            s.a = vec3s(npz(file, p + "a"));
            s.b = vec3s(npz(file, p + "b"));
            s.radius = npz(file, p + "radius").toDouble();
            s.owner = ints(npz(file, p + "owner"));
            s.atom = ints(npz(file, p + "atom"));
            it.prim = std::move(s);
        }
        items.push_back(std::move(it));
    }
    return items;
}

struct RenderCase {
    std::string name;
    std::vector<Colored> items;
    Viewport vp;
    Shading shading;
    RenderQuality quality = RenderQuality::High;
    std::vector<std::uint8_t> expected;
    std::size_t width = 0, height = 0;
};

RenderCase renderCase(const std::string& name) {
    const Value& rec = at(at(cases(), "renders"), name);
    RenderCase c;
    c.name = name;
    const Value& src = at(rec, "source");
    if (const Value* scene = src.find("scene")) {
        const std::string& s = *scene->get<std::string>();
        c.items = storedItems("scenes.npz", s, at(at(at(cases(), "scenes"), s), "items"));
    } else {
        c.items = storedItems("renders.npz", name, at(rec, "items"));
    }
    const Value& vp = at(rec, "viewport");
    const Array* rot = at(vp, "rotation").get<Array>();
    for (std::size_t i = 0; i < 3; ++i) {
        c.vp.rotation[i] = vec3((*rot)[i]);
    }
    const std::vector<double> xl = doubles(at(vp, "xlim")), yl = doubles(at(vp, "ylim")), zr = doubles(at(vp, "zrange"));
    c.vp.xlim = {xl[0], xl[1]};
    c.vp.ylim = {yl[0], yl[1]};
    c.vp.zrange = {zr[0], zr[1]};
    c.vp.width = static_cast<int>(number(at(vp, "width")));
    c.vp.height = static_cast<int>(number(at(vp, "height")));
    const Value& sh = at(rec, "shading");
    c.shading.fog_rgb = vec3(at(sh, "fog_rgb"));
    c.shading.fog_strength = number(at(sh, "fog_strength"));
    if (!at(sh, "glow").isNull()) {
        c.shading.glow = doubles(at(sh, "glow"));
    }
    c.shading.glow_rgb = vec3(at(sh, "glow_rgb"));
    c.quality = *at(rec, "quality").get<std::string>() == "fast" ? RenderQuality::Fast : RenderQuality::High;
    const std::string key = name + "/rgba";
    const Arr img = npz("renders.npz", key);
    c.height = img.shape.at(0);
    c.width = img.shape.at(1);
    for (double v : img.values) {
        c.expected.push_back(static_cast<std::uint8_t>(v));
    }
    return c;
}

std::vector<std::string> renderNames() {
    std::vector<std::string> out;
    for (const auto& [k, v] : *at(cases(), "renders").get<core::Table>()) {
        out.push_back(k);
    }
    return out;
}

struct Compare {
    double equal = 0;      // fraction of pixels with all four channels equal
    int max_diff = 0;      // largest channel difference
    double iou = 0;        // coverage (alpha > 0) intersection over union
};

Compare compare(const std::vector<std::uint8_t>& got, const std::vector<std::uint8_t>& want) {
    Compare c;
    std::size_t same = 0, inter = 0, uni = 0;
    const std::size_t n = want.size() / 4;
    for (std::size_t i = 0; i < n; ++i) {
        bool eq = true;
        for (std::size_t k = 0; k < 4; ++k) {
            const int d = std::abs(static_cast<int>(got[4 * i + k]) - static_cast<int>(want[4 * i + k]));
            c.max_diff = std::max(c.max_diff, d);
            eq = eq && d == 0;
        }
        same += eq ? 1 : 0;
        const bool a = got[4 * i + 3] > 0, b = want[4 * i + 3] > 0;
        inter += (a && b) ? 1 : 0;
        uni += (a || b) ? 1 : 0;
    }
    c.equal = n ? static_cast<double>(same) / static_cast<double>(n) : 1.0;
    c.iou = uni ? static_cast<double>(inter) / static_cast<double>(uni) : 1.0;
    return c;
}

// ---- synthetic items (as MDATK's tests build them) -----------------------------------------------------

Viewport view(const std::array<Vec3, 3>& R = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}, double lim = 10.0, int size = 100,
              std::array<double, 2> zrange = {-10.0, 10.0}, int height = 0) {
    Viewport v;
    v.rotation = R;
    v.xlim = {-lim, lim};
    v.ylim = {-lim, lim};
    v.width = size;
    v.height = height ? height : size;
    v.zrange = zrange;
    return v;
}

std::array<Vec3, 3> rotation(double ax, double ay, double az) {
    const Mat3 m = eulerToMatrix(Vec3d{ax, ay, az});
    return {Vec3{m[0][0], m[0][1], m[0][2]}, Vec3{m[1][0], m[1][1], m[1][2]}, Vec3{m[2][0], m[2][1], m[2][2]}};
}

Colored spheres(const std::vector<Vec3>& xyz, double r, Vec3 rgb = {0.8, 0.2, 0.2},
                std::vector<std::int32_t> owner = {}, double opacity = 1.0) {
    Spheres s;
    s.xyz = xyz;
    s.radius.assign(xyz.size(), r);
    if (owner.empty()) {
        for (std::size_t i = 0; i < xyz.size(); ++i) {
            owner.push_back(static_cast<std::int32_t>(i));
        }
    }
    s.owner = owner;
    s.atom.assign(xyz.size(), -1);
    Colored c;
    c.prim = std::move(s);
    c.rgb.assign(xyz.size(), rgb);
    c.opacity = opacity;
    return c;
}

Colored stick(Vec3 a, Vec3 b, double r, Vec3 rgb, double opacity = 1.0) {
    Sticks s;
    s.a = {a};
    s.b = {b};
    s.radius = {r};
    s.owner = {0};
    s.atom = {-1};
    Colored c;
    c.prim = std::move(s);
    c.rgb = {rgb};
    c.opacity = opacity;
    return c;
}

Colored mesh(std::vector<Vec3> verts, std::vector<Tri> tris, std::vector<std::int32_t> owner = {},
             std::vector<Vec3> rgb = {}, std::vector<Vec3> normals = {}) {
    Mesh m;
    const std::size_t n = verts.size();
    m.verts = std::move(verts);
    m.normals = normals.empty() ? std::vector<Vec3>(n, Vec3{0, 0, 1}) : std::move(normals);
    m.tris = std::move(tris);
    m.owner = owner.empty() ? std::vector<std::int32_t>(n, 0) : std::move(owner);
    m.role.assign(n, 0);
    Colored c;
    c.prim = std::move(m);
    c.rgb = rgb.empty() ? std::vector<Vec3>(n, Vec3{0.5, 0.5, 0.5}) : std::move(rgb);
    return c;
}

const Shading kPlain{{1.0, 1.0, 1.0}, 0.0, std::nullopt, {1.0, 0.69, 0.0}};

Shading glowShading(std::vector<double> glow, Vec3 rgb = {1.0, 0.7, 0.0}) {
    Shading s{{1.0, 1.0, 1.0}, 0.0, std::move(glow), rgb};
    return s;
}

const std::uint8_t* px(const RgbaImage& img, std::size_t x, std::size_t y) { return img.pixel(x, y); }

std::vector<Vec3> uniformPoints(std::mt19937_64& rng, std::size_t n, double lo, double hi) {
    std::uniform_real_distribution<double> u(lo, hi);
    std::vector<Vec3> out(n);
    for (Vec3& p : out) {
        p = {u(rng), u(rng), u(rng)};
    }
    return out;
}
std::vector<Vec3> normalPoints(std::mt19937_64& rng, std::size_t n, double scale) {
    std::normal_distribution<double> g(0.0, scale);
    std::vector<Vec3> out(n);
    for (Vec3& p : out) {
        p = {g(rng), g(rng), g(rng)};
    }
    return out;
}

// ---- parity with MDATK's renderer ------------------------------------------------------------------

TEST(MoleculeRender, MatchesMdatkRenders) {
    const std::vector<std::string> names = renderNames();
    ASSERT_GE(names.size(), 37u);
    const SoftwareRenderer renderer;
    for (const std::string& name : names) {
        SCOPED_TRACE(name);
        const RenderCase c = renderCase(name);
        const RgbaImage img = renderer.render(c.items, c.vp, c.shading, c.quality);
        ASSERT_EQ(img.width, c.width);
        ASSERT_EQ(img.height, c.height);
        const Compare cmp = compare(img.rgba, c.expected);
        std::printf("%-26s equal %.5f  max diff %d  IoU %.5f\n", name.c_str(), cmp.equal, cmp.max_diff, cmp.iou);
        EXPECT_GE(cmp.equal, 0.995);
        EXPECT_LE(cmp.max_diff, 2);
        EXPECT_GE(cmp.iou, 0.999);
        const std::size_t cover = static_cast<std::size_t>(number(at(at(at(cases(), "renders"), name), "coverage")));
        std::size_t got = 0;
        for (std::size_t i = 3; i < img.rgba.size(); i += 4) {
            got += img.rgba[i] > 0 ? 1 : 0;
        }
        EXPECT_NEAR(static_cast<double>(got), static_cast<double>(cover), 0.001 * static_cast<double>(cover) + 1);
    }
}

TEST(MoleculeRender, ThreadCountDoesNotChangeTheImage) {
    for (const std::string& name : {std::string("demo_surf_residue"), std::string("demo_mixed"),
                                    std::string("culling_40_-70_10"), std::string("glow_on")}) {
        SCOPED_TRACE(name);
        const RenderCase c = renderCase(name);
        Viewport big = c.vp;  // large enough for several threads
        big.width = 300;
        big.height = 260;
        const RgbaImage one = SoftwareRenderer({1, 64}).render(c.items, big, c.shading, c.quality);
        for (unsigned t : {2u, 3u, 8u}) {
            EXPECT_EQ(SoftwareRenderer({t, 64}).render(c.items, big, c.shading, c.quality).rgba, one.rgba) << t;
        }
    }
}

// ---- counterparts of MDATK's tests/test_molecule_render.py ---------------------------------------------

TEST(MoleculeRender, SphereAreaAndDepthAreExact) {
    const ResolvedPixels r = SoftwareRenderer().resolve({spheres({{0, 0, 0}}, 5.0)}, view());
    const double rp = 5.0 * r.scale;  // 25 px
    EXPECT_NEAR(static_cast<double>(r.pixel.size()), std::numbers::pi * rp * rp, 0.02 * std::numbers::pi * rp * rp);
    const auto centre = static_cast<std::int64_t>((r.height / 2) * r.width + r.width / 2);
    const auto it = std::find(r.pixel.begin(), r.pixel.end(), centre);
    ASSERT_NE(it, r.pixel.end());
    EXPECT_NEAR(r.depth[static_cast<std::size_t>(it - r.pixel.begin())], 5.0, 0.02);  // the front, toward the viewer
    EXPECT_EQ(std::set<std::int64_t>(r.owner.begin(), r.owner.end()), std::set<std::int64_t>{0});
}

TEST(MoleculeRender, StickWidthAndEnds) {
    const ResolvedPixels r = SoftwareRenderer().resolve({stick({-6, 0, 0}, {6, 0, 0}, 1.0, {0.2, 0.4, 0.8})}, view());
    std::int64_t ymin = 1 << 30, ymax = -1, xmin = 1 << 30, xmax = -1;
    const auto W = static_cast<std::int64_t>(r.width);
    for (std::int64_t p : r.pixel) {
        if (p % W == W / 2) {
            ymin = std::min(ymin, p / W);
            ymax = std::max(ymax, p / W);
        }
        xmin = std::min(xmin, p % W);
        xmax = std::max(xmax, p % W);
    }
    EXPECT_NEAR(static_cast<double>(ymax - ymin + 1), 2 * 1.0 * r.scale, 1.0);
    EXPECT_GE(static_cast<double>(xmin), (4 - 0.1) * r.scale - 1);  // open, unextended ends
    EXPECT_LE(static_cast<double>(xmax), (16 + 0.1) * r.scale + 1);
}

// cover with its holes filled (scipy's binary_fill_holes: background not connected to the border)
std::vector<bool> fillHoles(const std::vector<bool>& cover, std::size_t W, std::size_t H) {
    std::vector<bool> outside(cover.size(), false);
    std::vector<std::size_t> stack;
    for (std::size_t y = 0; y < H; ++y) {
        for (std::size_t x = 0; x < W; ++x) {
            if ((x == 0 || y == 0 || x + 1 == W || y + 1 == H) && !cover[y * W + x]) {
                stack.push_back(y * W + x);
                outside[y * W + x] = true;
            }
        }
    }
    while (!stack.empty()) {
        const std::size_t p = stack.back();
        stack.pop_back();
        const std::size_t x = p % W, y = p / W;
        const std::size_t nb[4] = {x > 0 ? p - 1 : p, x + 1 < W ? p + 1 : p, y > 0 ? p - W : p, y + 1 < H ? p + W : p};
        for (std::size_t q : nb) {
            if (!outside[q] && !cover[q]) {
                outside[q] = true;
                stack.push_back(q);
            }
        }
    }
    std::vector<bool> out(cover.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = !outside[i];
    }
    return out;
}

// MDATK's test asserts binary_fill_holes(cover) == cover; a Gaussian surface of random atoms can
// have a tunnel, whose projection is a real hole (seed 3 below: one pixel, in MDATK's renderer too),
// so here an enclosed uncovered pixel must also lie outside every triangle's projection.
TEST(MoleculeRender, ClosedMeshHasNoPinholes) {
    for (unsigned seed = 0; seed < 8; ++seed) {
        SCOPED_TRACE(seed);
        std::mt19937_64 rng(seed);
        const std::vector<Vec3> atoms = normalPoints(rng, 12, 2.5);
        const GaussianSurface s = gaussianSurface(atoms);
        Colored m = mesh(s.verts, s.tris, {}, {}, s.normals);
        std::uniform_real_distribution<double> ang(-180, 180);
        const double a = ang(rng), b = ang(rng), c = ang(rng);
        const Viewport v = view(rotation(a, b, c), 12, 160);
        const ResolvedPixels r = SoftwareRenderer().resolve({m}, v);
        std::vector<bool> cover(r.width * r.height, false);
        for (std::int64_t p : r.pixel) {
            cover[static_cast<std::size_t>(p)] = true;
        }
        std::vector<std::array<double, 2>> xy;
        for (const Vec3& p : s.verts) {
            const auto& R = v.rotation;
            xy.push_back({(p[0] * R[0][0] + p[1] * R[0][1] + p[2] * R[0][2] + 12) * r.scale,
                          (12 - (p[0] * R[1][0] + p[1] * R[1][1] + p[2] * R[1][2])) * r.scale});
        }
        const std::vector<bool> filled = fillHoles(cover, r.width, r.height);
        for (std::size_t i = 0; i < cover.size(); ++i) {
            if (filled[i] == cover[i]) {
                continue;
            }
            const double px0 = static_cast<double>(i % r.width) + 0.5, py0 = static_cast<double>(i / r.width) + 0.5;
            for (const Tri& t : s.tris) {
                const auto& p0 = xy[static_cast<std::size_t>(t[0])];
                const auto& p1 = xy[static_cast<std::size_t>(t[1])];
                const auto& p2 = xy[static_cast<std::size_t>(t[2])];
                auto cross = [&](const std::array<double, 2>& u, const std::array<double, 2>& w) {
                    return (w[0] - u[0]) * (py0 - u[1]) - (px0 - u[0]) * (w[1] - u[1]);
                };
                const double A = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
                const double l0 = cross(p1, p2) / A, l1 = cross(p2, p0) / A, l2 = cross(p0, p1) / A;
                EXPECT_FALSE(std::fabs(A) > 1e-3 && l0 >= 0 && l1 >= 0 && l2 >= 0)
                    << "pinhole at " << i % r.width << ", " << i / r.width;
            }
        }
    }
}

TEST(MoleculeRender, PrimitiveOrderDoesNotChangeTheImage) {
    std::mt19937_64 rng(3);
    const std::vector<Vec3> xyz = uniformPoints(rng, 40, -6, 6);
    const std::vector<Vec3> rgb = uniformPoints(rng, 40, 0.2, 0.9);
    Colored a = spheres(xyz, 1.6);
    a.rgb = rgb;
    std::vector<std::size_t> perm(40);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<Vec3> sx, sr;
    std::vector<std::int32_t> so;
    for (std::size_t i : perm) {
        sx.push_back(xyz[i]);
        sr.push_back(rgb[i]);
        so.push_back(static_cast<std::int32_t>(i));
    }
    Colored b = spheres(sx, 1.6, {}, so);
    b.rgb = sr;
    const SoftwareRenderer r;
    EXPECT_EQ(r.render({a}, view(), kPlain).rgba, r.render({b}, view(), kPlain).rgba);
}

TEST(MoleculeRender, GeometryBeyondThePanelDepthRangeKeepsItsOrder) {
    const Colored front = spheres({{0, 0, 60.0}}, 3.0, {1, 0, 0});
    const Colored back = spheres({{1.0, 0, 57.0}}, 3.0, {0, 0, 1});
    const RgbaImage img = SoftwareRenderer().render({back, front}, view(view().rotation, 10, 100, {-10, 10}), kPlain);
    const std::uint8_t* c = px(img, 50, 50);
    EXPECT_GT(c[0], 2 * c[2]);  // the front (red) sphere wins at the centre
}

TEST(MoleculeRender, ImageKeepsHandedness) {
    // rotation about z: model +x goes to view +y (right-handed); rows run downward
    const RgbaImage img = SoftwareRenderer().render({spheres({{4.0, 0, 0}}, 1.2, {1, 0, 0})},
                                                    view(rotation(0, 0, 90), 10, 300), kPlain);
    double sx = 0, sy = 0, n = 0;
    for (std::size_t y = 0; y < img.height; ++y) {
        for (std::size_t x = 0; x < img.width; ++x) {
            if (px(img, x, y)[3] == 255) {
                sx += static_cast<double>(x) + 0.5;
                sy += static_cast<double>(y) + 0.5;
                n += 1;
            }
        }
    }
    ASSERT_GT(n, 0);
    EXPECT_NEAR(sx / n, 150.0, 1.5);                 // view x = 0
    EXPECT_NEAR(sy / n, (10.0 - 4.0) * 15.0, 1.5);   // view y = +4: above the centre
}

TEST(MoleculeRender, FogFadesDistantPartsTowardThePanelColour) {
    const Colored near = spheres({{-5, 0, 8.0}}, 2.0, {0.9, 0.1, 0.1});
    const Colored far = spheres({{5, 0, -8.0}}, 2.0, {0.9, 0.1, 0.1});
    const Shading fogged{{1.0, 1.0, 1.0}, 0.45, std::nullopt, {1.0, 0.69, 0.0}};
    const RgbaImage img = SoftwareRenderer().render({near, far}, view(), fogged);
    EXPECT_GT(px(img, 75, 50)[1], px(img, 25, 50)[1] + 40);  // the far sphere is paler (closer to white)
}

TEST(MoleculeRender, GlowHaloOnlyWithGlow) {
    const std::vector<Colored> items{spheres({{0, 0, 0}}, 2.0, {0.3, 0.3, 0.8}, {0})};
    const SoftwareRenderer r;
    const RgbaImage plain = r.render(items, view(), glowShading({0.0}));
    const RgbaImage glow = r.render(items, view(), glowShading({1.0}));
    const std::size_t x = static_cast<std::size_t>(std::lround((10 + 3.2) * 5));  // beyond the sphere and outline
    EXPECT_EQ(px(plain, x, 50)[3], 0);
    EXPECT_GT(px(glow, x, 50)[3], 60);
    EXPECT_EQ(px(glow, x, 50)[0], 255);
    EXPECT_EQ(px(glow, x, 50)[1], 179);
    EXPECT_EQ(px(glow, x, 50)[2], 0);
    EXPECT_GT(px(glow, 50, 50)[0], px(plain, 50, 50)[0]);  // tinted toward the glow colour
}

TEST(MoleculeRender, TranslucentRepresentationBlendsOverTheOpaqueOne) {
    const Colored opaque = spheres({{0, 0, 0}}, 3.0, {1, 0, 0});
    const Colored glass = spheres({{0, 0, 4.0}}, 3.0, {0, 0, 1}, {}, 0.5);
    const RgbaImage img = SoftwareRenderer().render({opaque, glass}, view(), kPlain);
    const std::uint8_t* c = px(img, 50, 50);
    EXPECT_EQ(c[3], 255);
    EXPECT_GT(c[0], 40);
    EXPECT_GT(c[2], 40);  // both show
    const RgbaImage alone = SoftwareRenderer().render({glass}, view(), kPlain);
    EXPECT_GT(px(alone, 50, 50)[3], 100);
    EXPECT_LT(px(alone, 50, 50)[3], 160);  // straight alpha of the translucent surface
}

TEST(MoleculeRender, FastQualityReturnsTheRequestedSize) {
    const RgbaImage img = SoftwareRenderer().render({spheres({{0, 0, 0}}, 4.0)}, view(view().rotation, 10, 101, {-10, 10}, 77),
                                                    kPlain, RenderQuality::Fast);
    EXPECT_EQ(img.width, 101u);
    EXPECT_EQ(img.height, 77u);
    bool any = false;
    for (std::size_t i = 3; i < img.rgba.size(); i += 4) {
        any = any || img.rgba[i] > 0;
    }
    EXPECT_TRUE(any);
}

TEST(MoleculeRender, ThreadsRenderLikeSerial) {
    std::mt19937_64 rng(5);
    std::vector<std::vector<Colored>> scenes;
    for (int k = 0; k < 4; ++k) {
        const std::vector<Vec3> xyz = uniformPoints(rng, 30, -6, 6);
        const Vec3 rgb = uniformPoints(rng, 1, 0, 1)[0];
        scenes.push_back({spheres(xyz, 1.5, rgb)});
    }
    const SoftwareRenderer r;
    std::vector<RgbaImage> serial, out(4);
    for (const auto& s : scenes) {
        serial.push_back(r.render(s, view(view().rotation, 10, 150), kPlain));
    }
    std::vector<std::thread> threads;
    for (std::size_t k = 0; k < 4; ++k) {
        threads.emplace_back([&, k] {
            for (int i = 0; i < 3; ++i) {
                out[k] = r.render(scenes[k], view(view().rotation, 10, 150), kPlain);
            }
        });
    }
    for (std::thread& t : threads) {
        t.join();
    }
    for (std::size_t k = 0; k < 4; ++k) {
        EXPECT_EQ(out[k].rgba, serial[k].rgba);
    }
}

TEST(MoleculeRender, SphereCullingLeavesTheImageUnchanged) {
    std::mt19937_64 rng(11);
    const std::vector<Vec3> xyz = normalPoints(rng, 400, 4.0);
    Colored items = spheres(xyz, 1.7);
    const double radii[3] = {1.52, 1.7, 1.8};
    std::uniform_int_distribution<int> pick(0, 2);
    for (double& r : std::get<Spheres>(items.prim).radius) {
        r = radii[pick(rng)];
    }
    items.rgb = uniformPoints(rng, 400, 0.2, 0.9);
    const std::vector<Colored> list{items};
    const SoftwareRenderer culling;
    const SoftwareRenderer all({0, std::size_t{1} << 30});
    for (const auto& a : {Vec3{0, 0, 0}, Vec3{40, -70, 10}, Vec3{120, 30, -60}}) {
        const Viewport v = view(rotation(a[0], a[1], a[2]), 14, 240);
        EXPECT_GT(culling.culledSpheres(std::get<Spheres>(items.prim), list, v), 0u);  // some culled
        EXPECT_EQ(all.culledSpheres(std::get<Spheres>(items.prim), list, v), 0u);
        EXPECT_EQ(culling.render(list, v, kPlain).rgba, all.render(list, v, kPlain).rgba);
    }
}

TEST(MoleculeRender, ManyItemsRender) {
    std::mt19937_64 rng(2);
    std::vector<Colored> items;
    for (int i = 0; i < 80; ++i) {
        items.push_back(spheres(uniformPoints(rng, 1, -8, 8), 0.8, uniformPoints(rng, 1, 0.2, 0.9)[0]));
    }
    const RgbaImage img = SoftwareRenderer().render(items, view(view().rotation, 10, 120), kPlain);
    EXPECT_EQ(img.width, 120u);
    EXPECT_FALSE(SoftwareRenderer().resolve(items, view(view().rotation, 10, 120)).pixel.empty());
    bool any = false;
    for (std::size_t i = 3; i < img.rgba.size(); i += 4) {
        any = any || img.rgba[i] > 0;
    }
    EXPECT_TRUE(any);
}

TEST(MoleculeRender, TranslucentSurfacesResolveFrontMostOnce) {
    const Colored front = spheres({{0, 0, 2.0}}, 3.0, {1, 0, 0}, {}, 0.5);
    const Colored back = spheres({{0, 0, -2.0}}, 3.0, {0, 0, 1}, {}, 0.5);
    const SoftwareRenderer r;
    const RgbaImage a = r.render({front, back}, view(), kPlain);
    const RgbaImage b = r.render({back, front}, view(), kPlain);
    EXPECT_EQ(a.rgba, b.rgba);
    EXPECT_GT(px(a, 50, 50)[0], 150);
    EXPECT_LT(px(a, 50, 50)[2], 30);
    EXPECT_EQ(px(a, 50, 50)[3], 128);
    const Colored st = stick({-5, 0, 0}, {5, 0, 0}, 0.6, {0.2, 0.4, 0.9}, 0.5);
    const Colored ball = spheres({{-5, 0, 0}, {5, 0, 0}}, 0.6, {0.2, 0.4, 0.9}, {0, 0}, 0.5);
    const RgbaImage img = r.render({st, ball}, view(), kPlain);
    for (std::size_t x = 0; x < img.width; ++x) {  // one layer everywhere
        const std::uint8_t v = px(img, x, 50)[3];
        EXPECT_TRUE(v == 0 || v == 127 || v == 128) << x << ": " << int(v);
    }
}

TEST(MoleculeRender, TranslucentGlowingResidueKeepsItsHalo) {
    const Colored item = spheres({{0, 0, 0}}, 2.0, {0.3, 0.3, 0.8}, {0}, 0.6);
    const RgbaImage img = SoftwareRenderer().render({item}, view(), glowShading({1.0}));
    EXPECT_GT(px(img, static_cast<std::size_t>(std::lround((10 + 3.2) * 5)), 50)[3], 60);
}

TEST(MoleculeRender, SliverBehindASurfaceStaysHidden) {
    const Colored plane = mesh({{-9, -9, 0}, {9, -9, 0}, {9, 9, 0}, {-9, 9, 0}}, {{0, 1, 2}, {0, 2, 3}}, {},
                               std::vector<Vec3>(4, Vec3{0, 0, 1}));
    std::mt19937_64 rng(0);
    std::uniform_real_distribution<double> u(-5, 5);
    for (int i = 0; i < 40; ++i) {
        const double y0 = u(rng);  // nearly edge-on: projected altitude ~1e-4 px
        const Colored sl = mesh({{-6, y0, -0.5}, {6, y0, -0.5}, {0, y0 + 2e-5, -3.0}}, {{0, 1, 2}}, {},
                                std::vector<Vec3>(3, Vec3{1, 0, 0}));
        const RgbaImage img = SoftwareRenderer().render({plane, sl}, view(view().rotation, 10, 200), kPlain);
        for (std::size_t p = 0; p < img.rgba.size(); p += 4) {
            ASSERT_FALSE(img.rgba[p] > 150 && img.rgba[p + 2] < 100) << y0;
        }
    }
}

TEST(MoleculeRender, NonFiniteAtomsDoNotBreakDepthOrSurfaces) {
    std::mt19937_64 rng(1);
    std::vector<Vec3> xyz = normalPoints(rng, 60, 4);
    std::vector<Vec3> bad = xyz;
    bad[0] = {std::nan(""), std::nan(""), std::nan("")};
    const std::vector<Vec3> rest(xyz.begin() + 1, xyz.end());
    const auto a = depthRange({spheres(bad, 1.5)}, {0, 0, 1});
    const auto b = depthRange({spheres(rest, 1.5)}, {0, 0, 1});
    EXPECT_DOUBLE_EQ(a.first, b.first);
    EXPECT_DOUBLE_EQ(a.second, b.second);
    const RgbaImage img = SoftwareRenderer().render({spheres(bad, 1.5)}, view(view().rotation, 14), kPlain);
    bool any = false;
    for (std::size_t i = 3; i < img.rgba.size(); i += 4) {
        any = any || img.rgba[i] > 0;
    }
    EXPECT_TRUE(any);
}

TEST(MoleculeRender, ColourBoundaryBetweenResiduesIsStraight) {
    std::vector<Vec3> verts;
    for (double x : {-6.0, 0.0, 6.0}) {
        for (double y : {-4.0, 4.0}) {
            verts.push_back({x, y, 0});
        }
    }
    const std::vector<std::int32_t> owner{0, 0, 0, 0, 1, 1};  // the middle ring belongs to residue 0
    std::vector<Vec3> rgb;
    for (std::int32_t o : owner) {
        rgb.push_back(o == 0 ? Vec3{1, 0, 0} : Vec3{0, 0, 1});
    }
    const Colored m = mesh(verts, {{0, 2, 3}, {0, 3, 1}, {2, 4, 5}, {2, 5, 3}}, owner, rgb);
    const RgbaImage img = SoftwareRenderer().render({m}, view(view().rotation, 10, 200), kPlain);
    std::size_t lo = 1000, hi = 0;
    for (std::size_t y = 70; y < 130; ++y) {
        std::size_t x = 100;
        while (x < 200 && px(img, x, y)[0] > px(img, x, y)[2]) {
            ++x;
        }
        lo = std::min(lo, x);
        hi = std::max(hi, x);
    }
    EXPECT_LE(hi - lo, 1u);  // a vertical seam, not a diagonal staircase
}

TEST(MoleculeRender, GlowHelpersForOtherBackends) {
    const std::vector<Colored> items{spheres({{-4, 0, 0}, {4, 0, 0}}, 1.5, {0.2, 0.2, 0.8}, {0, 1})};
    const Shading sh = glowShading({1.0, 0.0});
    const std::vector<Colored> t = tinted(items, sh);
    EXPECT_GT(t[0].rgb[0][0], 0.2);
    EXPECT_NEAR(t[0].rgb[1][0], 0.2, 1e-6);
    EXPECT_NEAR(t[0].rgb[1][2], 0.8, 1e-6);
    const std::optional<RgbaImage> halo = glowUnderlay(items, view(), sh);
    ASSERT_TRUE(halo.has_value());
    const std::size_t left = static_cast<std::size_t>(std::lround(6 - 3.0) * 5);
    const std::size_t right = static_cast<std::size_t>(std::lround((14 + 3.0) * 5));
    EXPECT_GT(px(*halo, left, 50)[3], 0);
    EXPECT_EQ(px(*halo, right, 50)[3], 0);
    EXPECT_FALSE(glowUnderlay(items, view(), kPlain).has_value());
    RgbaImage clear = *halo;
    std::fill(clear.rgba.begin(), clear.rgba.end(), 0);
    const RgbaImage comp = over(clear, *halo);  // a transparent image over the halo leaves the halo
    for (std::size_t i = 0; i < comp.rgba.size(); i += 4) {
        EXPECT_EQ(comp.rgba[i + 3], halo->rgba[i + 3]);
        if (halo->rgba[i + 3] > 0) {
            for (std::size_t k = 0; k < 3; ++k) {
                EXPECT_LE(std::abs(int(comp.rgba[i + k]) - int(halo->rgba[i + k])), 1);
            }
        }
    }
}

// ---- MDBridge's improvements ------------------------------------------------------------------------

// RN-26: the anti-aliasing of the molecule's box no longer draws a darker rectangle into the halo.
TEST(MoleculeRender, HaloHasNoSeamAtTheMoleculeBox) {
    const std::vector<Colored> items{spheres({{0, 0, 0}}, 2.0, {0.3, 0.3, 0.8}, {0})};
    const RgbaImage img = SoftwareRenderer().render(items, view(), glowShading({1.0}));
    // the sphere covers rows / columns 40..59, so the anti-aliased box is 38..61. Toward the sphere the
    // halo rises smoothly across the box edge: along column 40 (rows 20..41) and along row 38
    // (columns 30..50). MDATK: row 38 dropped from 61 to 48 at column 38 and stayed ~17 % darker.
    for (std::size_t y = 21; y <= 41; ++y) {
        EXPECT_GE(px(img, 40, y)[3], px(img, 40, y - 1)[3]) << y;
    }
    for (std::size_t x = 31; x <= 50; ++x) {
        EXPECT_GE(px(img, x, 38)[3], px(img, x - 1, 38)[3]) << x;
    }
    EXPECT_EQ(px(img, 40, 38)[3], 77);  // MDATK's value without the seam: 77 (with it: 65)
}

// RN-27: fast quality at odd sizes covers the last column / row, and keeps the scale.
TEST(MoleculeRender, FastQualityAtOddSizesCoversTheWholeImage) {
    Mesh m;
    m.verts = {{-10, -10, 0}, {10, -10, 0}, {10, 10, 0}, {-10, 10, 0}};
    m.normals.assign(4, Vec3{0, 0, 1});
    m.tris = {{0, 1, 2}, {0, 2, 3}};
    m.owner.assign(4, 0);
    m.role.assign(4, 0);
    Colored plane;
    plane.prim = m;
    plane.rgb.assign(4, Vec3{0.5, 0.5, 0.5});
    for (const auto& [w, h] : {std::pair{101, 77}, std::pair{103, 51}, std::pair{3, 5}}) {
        SCOPED_TRACE(std::to_string(w) + "x" + std::to_string(h));
        const RgbaImage img = SoftwareRenderer().render({plane}, view(view().rotation, 10, w, {-10, 10}, h), kPlain,
                                                        RenderQuality::Fast);
        ASSERT_EQ(img.width, static_cast<std::size_t>(w));
        ASSERT_EQ(img.height, static_cast<std::size_t>(h));
        for (std::size_t i = 3; i < img.rgba.size(); i += 4) {
            ASSERT_EQ(img.rgba[i], 255) << i / 4;
        }
    }
    // a sphere stays where the high-quality image draws it
    const Colored s = spheres({{6, -5, 0}}, 2.0);
    const Viewport v = view(view().rotation, 10, 101, {-10, 10}, 77);
    auto centroid = [](const RgbaImage& img) {
        double sx = 0, sy = 0, n = 0;
        for (std::size_t y = 0; y < img.height; ++y) {
            for (std::size_t x = 0; x < img.width; ++x) {
                if (img.pixel(x, y)[3] == 255) {
                    sx += static_cast<double>(x);
                    sy += static_cast<double>(y);
                    n += 1;
                }
            }
        }
        return std::pair{sx / n, sy / n};
    };
    const auto hq = centroid(SoftwareRenderer().render({s}, v, kPlain));
    const auto fq = centroid(SoftwareRenderer().render({s}, v, kPlain, RenderQuality::Fast));
    EXPECT_NEAR(hq.first, fq.first, 1.0);
    EXPECT_NEAR(hq.second, fq.second, 1.0);
}

// RN-28: unusable views raise a clear error; unplaceable primitives are skipped; huge spheres draw.
TEST(MoleculeRender, UnusableViewsAndPrimitives) {
    const std::vector<Colored> items{spheres({{0, 0, 0}}, 1.0)};
    const SoftwareRenderer r;
    Viewport v = view();
    v.xlim = {3.0, 3.0};
    EXPECT_THROW(r.render(items, v, kPlain), std::invalid_argument);
    v = view();
    v.ylim = {1.0, std::nan("")};
    EXPECT_THROW(r.render(items, v, kPlain), std::invalid_argument);
    v = view();
    v.rotation[1][1] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(r.render(items, v, kPlain), std::invalid_argument);
    v = view();
    v.width = 0;
    EXPECT_THROW(r.render(items, v, kPlain), std::invalid_argument);
    try {
        v = view();
        v.xlim = {1.0, -1.0};
        r.render(items, v, kPlain);
        FAIL();
    } catch (const std::invalid_argument& e) {
        EXPECT_EQ(std::string(e.what()), "3D view limits are empty, reversed or not finite");
    }
    // a sphere 1e30 Å away is skipped (its screen position is beyond 2^40 px)
    std::vector<Colored> far = items;
    far.push_back(spheres({{1e30, 0, 0}}, 1.0, {0, 0, 1}));
    const RgbaImage a = r.render(far, view(), kPlain);
    EXPECT_EQ(px(a, 50, 50)[3], 255);
    // a sphere far larger than the image covers it (MDATK: a sprite of 10^10 entries); the image's
    // border is softened by the anti-aliasing
    const RgbaImage big = r.render({spheres({{0, 0, -1e4}}, 2e4)}, view(), kPlain);
    for (std::size_t y = 1; y + 1 < big.height; ++y) {
        for (std::size_t x = 1; x + 1 < big.width; ++x) {
            ASSERT_EQ(px(big, x, y)[3], 255) << x << ", " << y;
        }
    }
    // non-finite glow counts as none
    const RgbaImage g = r.render(items, view(), glowShading({std::nan("")}));
    EXPECT_EQ(g.rgba, r.render(items, view(), kPlain).rgba);
}

// RN-29: the halo's strength depends only on its own size, not on earlier renders.
TEST(MoleculeRender, HaloDoesNotDependOnEarlierRenders) {
    const std::vector<Colored> items{spheres({{0, 0, 0}}, 2.0, {0.3, 0.3, 0.8}, {0})};
    const SoftwareRenderer r;
    const Viewport v1 = view(view().rotation, 10.0, 100), v2 = view(view().rotation, 10.0004, 100);
    const RgbaImage first = r.render(items, v2, glowShading({1.0}));
    (void)r.render(items, v1, glowShading({1.0}));
    EXPECT_EQ(r.render(items, v2, glowShading({1.0})).rgba, first.rgba);
}

// Rendering speed on the largest committed scene (printed; the threads must not be slower).
TEST(MoleculeRender, ThreadsSpeedUpLargeImages) {
    const RenderCase c = renderCase("demo_surf_residue");
    Viewport big = c.vp;
    big.width = 1280;
    big.height = 1280;
    auto time = [&](unsigned threads) {
        const SoftwareRenderer r({threads, 64});
        const auto t0 = std::chrono::steady_clock::now();
        const RgbaImage img = r.render(c.items, big, c.shading, c.quality);
        const auto t1 = std::chrono::steady_clock::now();
        EXPECT_EQ(img.width, 1280u);
        return std::chrono::duration<double>(t1 - t0).count();
    };
    const double one = time(1), many = time(0);
    std::printf("1280 x 1280 surface: %.3f s on 1 thread, %.3f s on %u\n", one, many,
                std::max(1u, std::thread::hardware_concurrency()));
    SUCCEED();
}

}  // namespace
}  // namespace mdbridge::render::molecule
