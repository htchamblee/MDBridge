// The 3D representations' geometry (src/render/molecule): counterparts of MDATK's
// tests/test_molecule_geometry.py, and the topology roles, DrawLists, surfaces, scenes and scipy
// filters MDATK computes (data/render/molecule_geometry, tests/parity/make_molecule_geometry_cases.py):
// vertices within 1e-5 Å, indices exact.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/molecule/geometry.hpp"
#include "mdbridge/render/molecule/ndimage.hpp"
#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/surface.hpp"
#include "mdbridge/render/molecule/topology.hpp"

namespace mdbridge::render::molecule {
namespace {

using core::Array;
using core::Value;

constexpr double kTol = 1e-5;  // Å (docs/specs/m4-plan.md, section 3)

// ---- data ------------------------------------------------------------------------------------------

std::filesystem::path dataDir() { return std::filesystem::path(MDBRIDGE_TEST_DATA) / "render" / "molecule_geometry"; }

const Value& loadJson(std::string_view file) {
    static std::map<std::string, Value> cache;
    const std::string name(file);
    auto it = cache.find(name);
    if (it == cache.end()) {
        std::ifstream in(dataDir() / name, std::ios::binary);
        if (!in) {
            throw std::runtime_error("cannot read " + (dataDir() / name).string());
        }
        std::stringstream text;
        text << in.rdbuf();
        it = cache.emplace(name, core::parseJson(text.str(), name)).first;
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
const Array& list(const Value& v) { return *v.get<Array>(); }
double number(const Value& v) {
    if (const auto* d = v.get<double>()) {
        return *d;
    }
    return static_cast<double>(*v.get<std::int64_t>());
}
std::int64_t integer(const Value& v) {
    if (const auto* b = v.get<bool>()) {
        return *b ? 1 : 0;
    }
    return *v.get<std::int64_t>();
}
const std::string& str(const Value& v) { return *v.get<std::string>(); }

std::vector<double> doubles(const Value& v) {
    std::vector<double> out;
    for (const Value& x : list(v)) {
        out.push_back(number(x));
    }
    return out;
}
std::vector<std::int64_t> ints(const Value& v) {
    std::vector<std::int64_t> out;
    for (const Value& x : list(v)) {
        out.push_back(integer(x));
    }
    return out;
}
std::vector<std::int32_t> ints32(const Value& v) {
    std::vector<std::int32_t> out;
    for (const Value& x : list(v)) {
        out.push_back(static_cast<std::int32_t>(integer(x)));
    }
    return out;
}
std::vector<bool> bools(const Value& v) {
    std::vector<bool> out;
    for (const Value& x : list(v)) {
        out.push_back(integer(x) != 0);
    }
    return out;
}
std::vector<std::string> strings(const Value& v) {
    std::vector<std::string> out;
    for (const Value& x : list(v)) {
        out.push_back(str(x));
    }
    return out;
}
std::vector<Vec3> vec3s(const Value& v) {
    const std::vector<double> d = doubles(v);
    std::vector<Vec3> out(d.size() / 3);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
    }
    return out;
}

struct SystemData {
    TopologyInput input;
    std::vector<std::string> group_names;
    std::vector<std::vector<bool>> member;
    std::map<std::int64_t, std::vector<Vec3>> frames;
};

const SystemData& system(std::string_view system_name) {
    static std::map<std::string, SystemData> cache;
    const std::string name(system_name);
    auto it = cache.find(name);
    if (it != cache.end()) {
        return it->second;
    }
    const Value& s = at(loadJson("systems.json"), name);
    SystemData out;
    out.input.names = strings(at(s, "names"));
    out.input.elements = strings(at(s, "elements"));
    out.input.resnames = strings(at(s, "resnames"));
    out.input.row = ints32(at(s, "row"));
    out.input.n_res = static_cast<std::int32_t>(integer(at(s, "n_res")));
    const std::vector<std::int32_t> b = ints32(at(s, "bonds"));
    for (std::size_t i = 0; i + 1 < b.size(); i += 2) {
        out.input.bonds.push_back({b[i], b[i + 1]});
    }
    out.input.chain = ints(at(s, "chain"));
    out.input.resindex = ints(at(s, "resindex"));
    for (const Value& g : list(at(s, "groups"))) {
        out.group_names.push_back(str(at(g, "name")));
        out.member.push_back(bools(at(g, "member")));
    }
    for (const auto& [key, frame] : *at(s, "frames").get<core::Table>()) {
        out.frames[std::stoll(key)] = vec3s(frame);
    }
    out.input.ref_xyz = out.frames.at(0);
    return cache.emplace(name, std::move(out)).first->second;
}

std::shared_ptr<const Topology> topology(std::string_view name, std::string_view backbone = "P") {
    TopologyInput input = system(name).input;
    input.backbone_atom = std::string(backbone);
    return std::make_shared<const Topology>(input);
}

// ---- comparisons -----------------------------------------------------------------------------------

bool close(double a, double b, double tol) {
    if (std::isnan(a) || std::isnan(b)) {
        return std::isnan(a) && std::isnan(b);
    }
    return std::abs(a - b) <= tol;
}

void expectVec3s(const std::vector<Vec3>& got, const Value& want, const std::string& what, double tol = kTol) {
    const std::vector<double> w = doubles(want);
    ASSERT_EQ(got.size() * 3, w.size()) << what;
    int reported = 0;
    for (std::size_t i = 0; i < got.size(); ++i) {
        for (std::size_t c = 0; c < 3; ++c) {
            if (!close(got[i][c], w[3 * i + c], tol) && reported++ < 5) {
                ADD_FAILURE() << what << "[" << i << "][" << c << "]: " << got[i][c] << " vs " << w[3 * i + c];
            }
        }
    }
}

void expectDoubles(const std::vector<double>& got, const Value& want, const std::string& what, double tol = kTol) {
    const std::vector<double> w = doubles(want);
    ASSERT_EQ(got.size(), w.size()) << what;
    int reported = 0;
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (!close(got[i], w[i], tol) && reported++ < 5) {
            ADD_FAILURE() << what << "[" << i << "]: " << got[i] << " vs " << w[i];
        }
    }
}

template <class T>
void expectInts(const std::vector<T>& got, const Value& want, const std::string& what) {
    const std::vector<std::int64_t> w = ints(want);
    ASSERT_EQ(got.size(), w.size()) << what;
    for (std::size_t i = 0; i < got.size(); ++i) {
        ASSERT_EQ(static_cast<std::int64_t>(got[i]), w[i]) << what << "[" << i << "]";
    }
}

void expectTris(const std::vector<Tri>& got, const Value& want, const std::string& what) {
    const std::vector<std::int64_t> w = ints(want);
    ASSERT_EQ(got.size() * 3, w.size()) << what;
    for (std::size_t i = 0; i < got.size(); ++i) {
        for (std::size_t c = 0; c < 3; ++c) {
            ASSERT_EQ(got[i][c], w[3 * i + c]) << what << "[" << i << "]";
        }
    }
}

void expectBools(const std::vector<bool>& got, const Value& want, const std::string& what) {
    const std::vector<bool> w = bools(want);
    ASSERT_EQ(got, w) << what;
}

void expectMesh(const Mesh& m, const Value& want, const std::string& what) {
    expectVec3s(m.verts, at(want, "verts"), what + ".verts");
    expectVec3s(m.normals, at(want, "normals"), what + ".normals");
    expectTris(m.tris, at(want, "tris"), what + ".tris");
    expectInts(m.owner, at(want, "owner"), what + ".owner");
    expectInts(m.role, at(want, "role"), what + ".role");
    const Value* atom = want.find("atom");
    ASSERT_EQ(m.has_atom, atom != nullptr) << what;
    if (atom) {
        expectInts(m.atom, *atom, what + ".atom");
    }
}

void expectSpheres(const Spheres& s, const Value& want, const std::string& what) {
    expectVec3s(s.xyz, at(want, "xyz"), what + ".xyz");
    expectDoubles(s.radius, at(want, "radius"), what + ".radius", 1e-12);
    expectInts(s.owner, at(want, "owner"), what + ".owner");
    expectInts(s.atom, at(want, "atom"), what + ".atom");
}

void expectSticks(const Sticks& s, const Value& want, const std::string& what) {
    expectVec3s(s.a, at(want, "a"), what + ".a");
    expectVec3s(s.b, at(want, "b"), what + ".b");
    expectDoubles(s.radius, at(want, "radius"), what + ".radius", 1e-12);
    expectInts(s.owner, at(want, "owner"), what + ".owner");
    expectInts(s.atom, at(want, "atom"), what + ".atom");
}

void expectDrawList(const DrawList& dl, const Value& want, const std::string& what) {
    const Array& meshes = list(at(want, "meshes"));
    const Array& spheres = list(at(want, "spheres"));
    const Array& sticks = list(at(want, "sticks"));
    ASSERT_EQ(dl.meshes.size(), meshes.size()) << what;
    ASSERT_EQ(dl.spheres.size(), spheres.size()) << what;
    ASSERT_EQ(dl.sticks.size(), sticks.size()) << what;
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        expectMesh(dl.meshes[i], meshes[i], what + ".mesh" + std::to_string(i));
    }
    for (std::size_t i = 0; i < spheres.size(); ++i) {
        expectSpheres(dl.spheres[i], spheres[i], what + ".spheres" + std::to_string(i));
    }
    for (std::size_t i = 0; i < sticks.size(); ++i) {
        expectSticks(dl.sticks[i], sticks[i], what + ".sticks" + std::to_string(i));
    }
}

std::vector<bool> maskOf(std::size_t n, std::initializer_list<int> rows) {
    std::vector<bool> m(n, false);
    for (int r : rows) {
        m[static_cast<std::size_t>(r)] = true;
    }
    return m;
}

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// Face normals agree with the vertex normals (non-degenerate faces).
void expectOutward(const Mesh& m, const std::string& what) {
    for (const Tri& t : m.tris) {
        const Vec3& v0 = m.verts[static_cast<std::size_t>(t[0])];
        const Vec3 face = cross(sub(m.verts[static_cast<std::size_t>(t[1])], v0), sub(m.verts[static_cast<std::size_t>(t[2])], v0));
        if (norm(face) <= 1e-9) {
            continue;
        }
        Vec3 n{0, 0, 0};
        for (std::int32_t v : t) {
            for (std::size_t c = 0; c < 3; ++c) {
                n[c] += m.normals[static_cast<std::size_t>(v)][c];
            }
        }
        ASSERT_GT(dot(face, n), 0.0) << what;
    }
}

// ---- MDATK's test_molecule_geometry.py -------------------------------------------------------------

TEST(MoleculeGeometry, NamesAndLetters) {
    EXPECT_EQ(normaliseName("C1*"), "C1'");
    EXPECT_EQ(normaliseName("O1P"), "OP1");
    std::string letters;
    for (const char* r : {"G", "RA5", "C3", "DT", "PSU", "U5", "MG", "LYS"}) {
        letters += baseLetter(r);
    }
    EXPECT_EQ(letters, "GACTUU??");
    const Value& data = loadJson("topology.json");
    for (const Value& pair : list(at(data, "normalise_name"))) {
        EXPECT_EQ(normaliseName(str(list(pair)[0])), str(list(pair)[1])) << str(list(pair)[0]);
    }
    for (const Value& pair : list(at(data, "base_letter"))) {
        EXPECT_EQ(std::string(1, baseLetter(str(list(pair)[0]))), str(list(pair)[1])) << str(list(pair)[0]);
    }
}

std::vector<Vec3> testPoints(std::size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-2.0, 2.0);
    std::vector<Vec3> p(n);
    for (Vec3& v : p) {
        v = {u(rng), u(rng), u(rng)};
    }
    return p;
}

std::vector<Vec3> applyDense(const std::vector<double>& a, std::size_t samples, const std::vector<Vec3>& p) {
    const std::size_t m = p.size();
    std::vector<Vec3> out(samples, Vec3{0, 0, 0});
    for (std::size_t s = 0; s < samples; ++s) {
        for (std::size_t j = 0; j < m; ++j) {
            for (std::size_t c = 0; c < 3; ++c) {
                out[s][c] += a[s * m + j] * p[j][c];
            }
        }
    }
    return out;
}

TEST(MoleculeGeometry, CatmullRomInterpolatesAndDifferentiates) {
    const std::vector<Vec3> p = testPoints(6, 0);
    const DenseSpline sp = catmullRom(p.size(), 8);
    const std::vector<Vec3> s = applyDense(sp.a, sp.samples, p);
    for (std::size_t i = 0; i < p.size(); ++i) {  // passes through the control points
        for (std::size_t c = 0; c < 3; ++c) {
            EXPECT_NEAR(s[8 * i][c], p[i][c], 1e-12);
        }
    }
    // the derivative matches finite differences of a densely sampled spline
    const DenseSpline d2 = catmullRom(p.size(), 200);
    const std::vector<Vec3> s2 = applyDense(d2.a, d2.samples, p);
    const std::vector<Vec3> ds2 = applyDense(d2.d, d2.samples, p);
    for (std::size_t i = 1; i + 1 < s2.size(); ++i) {
        if ((i - 1) % 200 == 199) {
            continue;  // the kinks at control points (C1 only)
        }
        for (std::size_t c = 0; c < 3; ++c) {
            const double fd = (s2[i + 1][c] - s2[i - 1][c]) / (2.0 / 200.0);
            EXPECT_NEAR(ds2[i][c], fd, 1e-3 + 1e-3 * std::abs(fd));
        }
    }
    EXPECT_EQ(sp.owner.front(), 0);
    EXPECT_EQ(sp.owner.back(), 5);
    const std::vector<std::int32_t> mid(sp.owner.begin() + 8, sp.owner.begin() + 16);
    EXPECT_EQ(mid, (std::vector<std::int32_t>{1, 1, 1, 1, 2, 2, 2, 2}));
}

TEST(MoleculeGeometry, SweepFacesCloseAndFramesAreOrthonormal) {
    const std::vector<Tri> tris = sweepFaces(5, 6);
    std::int32_t mx = 0;
    for (const Tri& t : tris) {
        mx = std::max({mx, t[0], t[1], t[2]});
    }
    EXPECT_EQ(mx, 5 * 6 + 1);
    EXPECT_EQ(tris.size(), 2u * 4 * 6 + 2 * 6);
    const std::size_t n = 50;
    std::vector<Vec3> pos(n), tan(n), side(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double t = 4 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n - 1);
        pos[i] = {std::cos(t), std::sin(t), 0.2 * t};
        side[i] = {-pos[i][0], -pos[i][1], 0.0};
    }
    for (std::size_t i = 0; i < n; ++i) {  // np.gradient
        const std::size_t a = i == 0 ? 0 : i - 1;
        const std::size_t b = i + 1 == n ? i : i + 1;
        const double h = static_cast<double>(b - a);
        tan[i] = {(pos[b][0] - pos[a][0]) / h, (pos[b][1] - pos[a][1]) / h, (pos[b][2] - pos[a][2]) / h};
    }
    const std::vector<Vec3>* sides[] = {nullptr, &side};
    for (const std::vector<Vec3>* s : sides) {
        const Frames f = framesAlong(pos, tan, s);
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(dot(f.t[i], f.n[i]), 0.0, 1e-9);
            EXPECT_NEAR(dot(f.t[i], f.b[i]), 0.0, 1e-9);
            EXPECT_NEAR(dot(f.n[i], f.b[i]), 0.0, 1e-9);
            EXPECT_NEAR(norm(f.n[i]), 1.0, 1e-9);
        }
    }
}

TEST(MoleculeGeometry, RingSlabIsClosedWithOutwardNormals) {
    std::vector<Vec3> ring(6);
    for (std::size_t j = 0; j < 6; ++j) {
        const double ang = 2 * std::numbers::pi * static_cast<double>(j) / 6.0;
        ring[j] = {1.4 * std::cos(ang), 1.4 * std::sin(ang), 0.0};
    }
    const RingSlabs rs = ringSlabs(ring, 6, 0.45);
    Mesh m;
    m.verts = rs.verts;
    m.normals = rs.normals;
    m.tris = rs.tris;
    for (const Tri& t : m.tris) {  // every face, none degenerate here
        const Vec3& v0 = m.verts[static_cast<std::size_t>(t[0])];
        const Vec3 face = cross(sub(m.verts[static_cast<std::size_t>(t[1])], v0), sub(m.verts[static_cast<std::size_t>(t[2])], v0));
        Vec3 n{0, 0, 0};
        for (std::int32_t v : t) {
            for (std::size_t c = 0; c < 3; ++c) {
                n[c] += m.normals[static_cast<std::size_t>(v)][c];
            }
        }
        EXPECT_GT(dot(face, n), 0.0);
    }
    // closed surface: the signed volume equals area * thickness
    double vol = 0.0;
    for (const Tri& t : m.tris) {
        vol += dot(m.verts[static_cast<std::size_t>(t[0])],
                   cross(m.verts[static_cast<std::size_t>(t[1])], m.verts[static_cast<std::size_t>(t[2])]));
    }
    vol /= 6.0;
    const double area = 3 * std::sqrt(3.0) / 2 * 1.4 * 1.4;
    EXPECT_NEAR(vol, area * 0.9, 1e-6 * area * 0.9);
    for (std::int32_t s : rs.slab) {
        EXPECT_EQ(s, 0);
    }
}

TEST(MoleculeGeometry, TopologyRoles) {
    const auto topo = topology("demo_dry");
    EXPECT_EQ(topo->n_res, 24);
    ASSERT_EQ(topo->segments.size(), 2u);  // two strands
    EXPECT_EQ(topo->segments[0].size(), 12u);
    EXPECT_EQ(topo->segments[1].size(), 12u);
    EXPECT_EQ(topo->purine_res.size() + topo->pyrimidine_res.size(), 24u);
    EXPECT_EQ(topo->sugar_res.size(), 24u);
    for (std::size_t q = 0; q < topo->purine_res.size(); ++q) {
        for (std::int32_t a : topo->purine_ring[q]) {
            EXPECT_EQ(topo->row[static_cast<std::size_t>(a)], topo->purine_res[q]);
        }
        EXPECT_EQ(topo->names[static_cast<std::size_t>(topo->purine_gly[q])], "N9");
    }
    for (std::size_t q = 0; q < topo->pyrimidine_res.size(); ++q) {
        for (std::int32_t a : topo->pyrimidine_ring[q]) {
            EXPECT_EQ(topo->row[static_cast<std::size_t>(a)], topo->pyrimidine_res[q]);
        }
        EXPECT_EQ(topo->names[static_cast<std::size_t>(topo->pyrimidine_gly[q])], "N1");
    }
    for (std::int32_t c : topo->control) {
        if (c >= 0) {
            EXPECT_EQ(topo->names[static_cast<std::size_t>(c)], "P");
        }
    }
}

std::vector<std::int32_t> ownersOf(const DrawList& dl) {
    std::vector<std::int32_t> out;
    for (const Mesh& m : dl.meshes) {
        out.insert(out.end(), m.owner.begin(), m.owner.end());
    }
    for (const Spheres& s : dl.spheres) {
        out.insert(out.end(), s.owner.begin(), s.owner.end());
    }
    for (const Sticks& s : dl.sticks) {
        out.insert(out.end(), s.owner.begin(), s.owner.end());
    }
    return out;
}

TEST(MoleculeGeometry, StylesBuildOutwardMeshesForSelectedResidues) {
    const auto topo = topology("demo_dry");
    const auto& xyz = system("demo_dry").frames.at(3);
    std::vector<bool> mask(24, false);
    std::fill(mask.begin(), mask.begin() + 12, true);  // first strand only
    for (const std::string& name : styleNames()) {
        SCOPED_TRACE(name);
        StyleOptions options;
        options.sugars = true;
        const auto style = makeStyle(name, topo, mask, options);
        ASSERT_TRUE(style);
        EXPECT_EQ(style->name(), name);
        const DrawList dl = style->build(xyz);
        const std::vector<std::int32_t> owners = ownersOf(dl);
        ASSERT_FALSE(owners.empty());
        for (std::int32_t o : owners) {
            EXPECT_TRUE(o >= 0 && o < 12);
        }
        for (const Mesh& m : dl.meshes) {
            expectOutward(m, name);
        }
        if (name == "cartoon") {
            EXPECT_TRUE(std::any_of(dl.meshes.begin(), dl.meshes.end(), [](const Mesh& m) {
                return !m.role.empty() && std::all_of(m.role.begin(), m.role.end(), [](std::int8_t r) { return r == 1; });
            }));
            ASSERT_FALSE(dl.sticks.empty());
            EXPECT_GE(dl.sticks[0].a.size(), 12u);
        }
        if (name == "vdw" || name == "licorice") {
            const auto n_atoms = static_cast<std::size_t>(
                std::count_if(topo->row.begin(), topo->row.end(), [](std::int32_t r) { return r < 12; }));
            EXPECT_EQ(dl.spheres[0].xyz.size(), n_atoms);
        }
    }
    EXPECT_EQ(makeStyle("bogus", topo, mask), nullptr);
}

TEST(MoleculeGeometry, BandedSplineMatchesDense) {
    const std::vector<Vec3> p = testPoints(7, 4);
    const DenseSpline dense = catmullRom(7, 8);
    const BandedSpline banded = catmullRomBanded(7, 8);
    const std::vector<Vec3> ext = extendControls(p);
    const std::vector<Vec3> a = applyDense(dense.a, dense.samples, p);
    const std::vector<Vec3> d = applyDense(dense.d, dense.samples, p);
    ASSERT_EQ(banded.idx.size(), dense.samples);
    for (std::size_t s = 0; s < dense.samples; ++s) {
        for (std::size_t c = 0; c < 3; ++c) {
            double pa = 0.0, pd = 0.0;
            for (std::size_t k = 0; k < 4; ++k) {
                pa += banded.w[s][k] * ext[static_cast<std::size_t>(banded.idx[s][k])][c];
                pd += banded.d[s][k] * ext[static_cast<std::size_t>(banded.idx[s][k])][c];
            }
            EXPECT_NEAR(pa, a[s][c], 1e-12);
            EXPECT_NEAR(pd, d[s][c], 1e-12);
        }
    }
    EXPECT_EQ(dense.owner, banded.owner);
}

TEST(MoleculeGeometry, CharmmNamesAndBasesToldByTheirAtoms) {
    std::string letters;
    for (const char* r : {"ADE", "GUA", "CYT", "URA", "THY", "URI"}) {
        letters += baseLetter(r);
    }
    EXPECT_EQ(letters, "AGCUTU");
    const auto topo = topology("charmm");
    EXPECT_EQ(std::string(topo->letter.begin(), topo->letter.end()), "GACUTT?");
}

TEST(MoleculeGeometry, UndrawableResiduesFallBack) {
    // residues no representation actually draws get the trace; residues without a backbone atom
    // (ligands, which backbone styles cannot draw) get licorice
    auto topo = topology("demo_dry");
    core::RepSpec tube;
    tube.style = "tube";
    MoleculeScene sc(topo, {tube}, {maskOf(24, {0, 2, 3})});
    const std::vector<bool> d = sc.drawn();
    EXPECT_TRUE(d[0] && d[2] && d[3]);
    EXPECT_FALSE(d[1]);
    sc.fillRest();
    ASSERT_EQ(sc.size(), 2u);
    EXPECT_EQ(sc.specs()[0].style, "tube");
    EXPECT_EQ(sc.specs()[1].style, "trace");
    EXPECT_TRUE(sc.style(1)->mask()[1]);
    EXPECT_FALSE(sc.style(1)->mask()[0]);

    auto edited = std::make_shared<Topology>(*topo);
    edited->control[5] = -1;  // pretend residue 5 has no backbone atom
    core::RepSpec cartoon;
    std::vector<bool> first(24, false);
    std::fill(first.begin(), first.begin() + 3, true);
    MoleculeScene sc2(edited, {cartoon}, {first});
    sc2.fillRest();
    const Style* lic = nullptr;
    for (std::size_t k = 0; k < sc2.size(); ++k) {
        if (sc2.specs()[k].style == "licorice") {
            lic = sc2.style(k);
        }
    }
    ASSERT_NE(lic, nullptr);
    EXPECT_TRUE(lic->mask()[5]);
    EXPECT_EQ(std::count(lic->mask().begin(), lic->mask().end(), true), 1);
}

TEST(MoleculeGeometry, ProteinSegmentsGetARoundTube) {
    Topology t;
    t.n_res = 6;
    t.segments = {{0, 1, 2, 3, 4, 5}};
    t.control = {0, 1, 2, 3, 4, 5};
    t.side = {6, 7, 8, 9, 10, 11};
    std::vector<Vec3> xyz(12, Vec3{0, 0, 0});
    for (std::size_t i = 0; i < 6; ++i) {
        xyz[i][0] = 3.8 * static_cast<double>(i);
        xyz[6 + i][0] = xyz[i][0];
        xyz[6 + i][1] = (i % 2 == 0 ? 1.0 : -1.0) * 1.5;  // side atoms alternate (CB on a strand)
    }
    const Backbone bb(t, std::vector<bool>(6, true), "ellipse", 1.6, 0.45, 8, 12, 0.6);
    const auto [meshes, anchor] = bb.build(xyz);
    ASSERT_EQ(meshes.size(), 1u);
    const auto& v = meshes[0].verts;
    ASSERT_EQ((v.size() - 2) % 12, 0u);
    for (std::size_t ring = 0; ring < (v.size() - 2) / 12; ++ring) {
        Vec3 c{0, 0, 0};
        for (std::size_t j = 0; j < 12; ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                c[k] += v[ring * 12 + j][k] / 12.0;
            }
        }
        for (std::size_t j = 0; j < 12; ++j) {
            EXPECT_NEAR(norm(sub(v[ring * 12 + j], c)), 0.6, 1e-6);
        }
    }
    EXPECT_EQ(anchor.size(), 6u);
}

TEST(MoleculeGeometry, CartoonRungsAreClosedByBalls) {
    const auto topo = topology("demo_dry");
    const auto& xyz = system("demo_dry").frames.at(0);
    for (const auto& [sugars, per_res] : {std::pair{false, 2}, std::pair{true, 3}}) {
        StyleOptions options;
        options.sugars = sugars;
        const DrawList dl = makeStyle("cartoon", topo, maskOf(24, {0, 1, 2, 3}), options)->build(xyz);
        const auto& balls = dl.spheres.at(0).xyz;
        EXPECT_EQ(balls.size(), static_cast<std::size_t>(4 * per_res));
        std::vector<Vec3> ends = dl.sticks.at(0).a;
        ends.insert(ends.end(), dl.sticks[0].b.begin(), dl.sticks[0].b.end());
        std::size_t inside = 0;  // every stick end except the backbone anchor sits in a ball
        for (const Vec3& e : ends) {
            double best = 1e300;
            for (const Vec3& b : balls) {
                best = std::min(best, norm(sub(e, b)));
            }
            inside += best < 1e-9 ? 1 : 0;
        }
        EXPECT_GE(inside, ends.size() - 4);
    }
}

// ---- against MDATK's numbers ----------------------------------------------------------------------

TEST(MoleculeGeometryParity, TopologyRoles) {
    for (const Value& c : list(at(loadJson("topology.json"), "topologies"))) {
        const std::string name = str(at(c, "system"));
        const std::string backbone = str(at(c, "backbone_atom"));
        SCOPED_TRACE(name + " " + backbone);
        const auto topo = topology(name, backbone);
        const Value& t = at(c, "topology");
        EXPECT_EQ(topo->names, strings(at(t, "names")));
        EXPECT_EQ(topo->elements, strings(at(t, "elements")));
        EXPECT_EQ(topo->resname, strings(at(t, "resname")));
        std::vector<std::string> letters;
        for (char l : topo->letter) {
            letters.emplace_back(1, l);
        }
        EXPECT_EQ(letters, strings(at(t, "letter")));
        expectDoubles(topo->radius, at(t, "radius"), "radius", 0.0);
        expectVec3s(topo->element_rgb, at(t, "element_rgb"), "element_rgb", 0.0);
        expectBools(topo->is_carbon, at(t, "is_carbon"), "is_carbon");
        expectBools(topo->is_hydrogen, at(t, "is_hydrogen"), "is_hydrogen");
        expectInts(topo->control, at(t, "control"), "control");
        expectInts(topo->anchor_atom, at(t, "anchor_atom"), "anchor_atom");
        expectInts(topo->side, at(t, "side"), "side");
        expectInts(topo->c1, at(t, "c1"), "c1");
        expectInts(topo->c4, at(t, "c4"), "c4");
        expectInts(topo->purine_res, at(t, "purine_res"), "purine_res");
        expectInts(topo->purine_gly, at(t, "purine_gly"), "purine_gly");
        expectInts(topo->pyrimidine_res, at(t, "pyrimidine_res"), "pyrimidine_res");
        expectInts(topo->pyrimidine_gly, at(t, "pyrimidine_gly"), "pyrimidine_gly");
        expectInts(topo->sugar_res, at(t, "sugar_res"), "sugar_res");
        auto flatRings = [](const auto& rings) {
            std::vector<std::int32_t> out;
            for (const auto& r : rings) {
                out.insert(out.end(), r.begin(), r.end());
            }
            return out;
        };
        expectInts(flatRings(topo->purine_ring), at(t, "purine_ring"), "purine_ring");
        expectInts(flatRings(topo->pyrimidine_ring), at(t, "pyrimidine_ring"), "pyrimidine_ring");
        expectInts(flatRings(topo->sugar_ring), at(t, "sugar_ring"), "sugar_ring");
        const Array& segs = list(at(t, "segments"));
        ASSERT_EQ(topo->segments.size(), segs.size());
        for (std::size_t s = 0; s < segs.size(); ++s) {
            expectInts(topo->segments[s], segs[s], "segment");
        }
        expectInts(topo->nucleic(), at(t, "nucleic"), "nucleic");
    }
}

TEST(MoleculeGeometryParity, StylesMatchMdatk) {
    const Array& cases = list(at(loadJson("styles.json"), "cases"));
    ASSERT_GT(cases.size(), 40u);
    std::size_t k = 0;
    for (const Value& c : cases) {
        const std::string name = str(at(c, "system"));
        const std::string style_name = str(at(c, "style"));
        const std::string what = "case " + std::to_string(k++) + " (" + name + " " + style_name + " frame " +
                                 std::to_string(integer(at(c, "frame"))) + " " + str(at(c, "backbone_atom")) + ")";
        SCOPED_TRACE(what);
        const auto topo = topology(name, str(at(c, "backbone_atom")));
        StyleOptions options;
        options.size = number(at(c, "size"));
        options.sugars = integer(at(c, "sugars")) != 0;
        options.samples = static_cast<std::size_t>(integer(at(c, "samples")));
        const auto style = makeStyle(style_name, topo, bools(at(c, "mask")), options);
        ASSERT_TRUE(style);
        expectBools(style->drawn(), at(c, "drawn"), what + " drawn");
        const DrawList dl = style->build(system(name).frames.at(integer(at(c, "frame"))));
        expectDrawList(dl, at(c, "drawlist"), what);
    }
}

TEST(MoleculeGeometryParity, SurfacesMatchMdatk) {
    for (const Value& c : list(at(loadJson("surface.json"), "cases"))) {
        const std::string name = str(at(c, "name"));
        SCOPED_TRACE(name);
        const std::vector<Vec3> xyz = vec3s(at(c, "xyz"));
        const double spacing = number(at(c, "spacing"));
        const double sigma = number(at(c, "sigma"));
        const double level = number(at(c, "level"));
        const DensityGrid grid = densityGrid(xyz, spacing, sigma);
        const std::vector<std::int64_t> shape = ints(at(c, "grid_shape"));
        ASSERT_EQ(shape, (std::vector<std::int64_t>{static_cast<std::int64_t>(grid.shape[0]),
                                                    static_cast<std::int64_t>(grid.shape[1]),
                                                    static_cast<std::int64_t>(grid.shape[2])}));
        expectDoubles({grid.origin[0], grid.origin[1], grid.origin[2]}, at(c, "grid_origin"), "origin", 0.0);
        if (const Value* values = c.find("grid")) {  // float32, rounded as scipy rounds it
            std::vector<double> got(grid.values.begin(), grid.values.end());
            expectDoubles(got, *values, "grid", 0.0);
        }
        const SurfaceNets nets = surfaceNets(grid, level);
        const Value& want_nets = at(c, "nets");
        expectVec3s(nets.verts, at(want_nets, "verts"), "nets.verts");
        expectTris(nets.tris, at(want_nets, "tris"), "nets.tris");
        const std::vector<double> coords = doubles(at(want_nets, "coords"));  // (3, n): gz, gy, gx
        ASSERT_EQ(coords.size(), 3 * nets.coords.size());
        for (std::size_t i = 0; i < nets.coords.size(); ++i) {
            for (std::size_t a = 0; a < 3; ++a) {
                EXPECT_NEAR(nets.coords[i][a], coords[a * nets.coords.size() + i], 1e-9);
            }
        }
        const GaussianSurface surf = gaussianSurface(xyz, spacing, sigma, level);
        const Value& want = at(c, "surface");
        expectVec3s(surf.verts, at(want, "verts"), "verts");
        expectVec3s(surf.normals, at(want, "normals"), "normals");
        expectTris(surf.tris, at(want, "tris"), "tris");
        // the nearest atom: cKDTree's choice, or an atom exactly as near
        const std::vector<std::int64_t> near = ints(at(want, "nearest"));
        const std::vector<double> dist = doubles(at(want, "nearest_distance"));
        ASSERT_EQ(surf.nearest.size(), near.size());
        for (std::size_t i = 0; i < near.size(); ++i) {
            if (surf.nearest[i] != near[i]) {
                EXPECT_NEAR(norm(sub(surf.verts[i], xyz[static_cast<std::size_t>(surf.nearest[i])])), dist[i], 1e-12) << i;
            }
        }
    }
}

TEST(MoleculeGeometryParity, ScenesMatchMdatk) {
    const Array& cases = list(at(loadJson("scenes.json"), "cases"));
    std::size_t k = 0;
    for (const Value& c : cases) {
        const std::string name = str(at(c, "system"));
        SCOPED_TRACE("scene case " + std::to_string(k++) + " (" + name + ")");
        const SystemData& sys = system(name);
        std::vector<core::RepSpec> specs;
        for (const Value& s : list(at(c, "specs"))) {
            core::RepSpec spec;
            spec.style = str(at(s, "style"));
            spec.groups = strings(at(s, "groups"));
            spec.color = str(at(s, "color"));
            spec.size = number(at(s, "size"));
            spec.opacity = number(at(s, "opacity"));
            spec.sugars = integer(at(s, "sugars")) != 0;
            specs.push_back(spec);
        }
        auto scene = buildScene(sys.input, sys.member, sys.group_names, specs, str(at(c, "rest")));
        const Array& want_specs = list(at(c, "scene_specs"));
        if (specs.empty()) {  // MDATK's build_scene: no scene without representations
            EXPECT_EQ(scene, nullptr);
            continue;
        }
        ASSERT_TRUE(scene);
        ASSERT_EQ(scene->size(), want_specs.size());
        const Array& masks = list(at(c, "scene_masks"));
        for (std::size_t i = 0; i < want_specs.size(); ++i) {
            EXPECT_EQ(scene->specs()[i].style, str(at(want_specs[i], "style")));
            EXPECT_EQ(scene->specs()[i].color, str(at(want_specs[i], "color")));
            EXPECT_EQ(scene->specs()[i].size, number(at(want_specs[i], "size")));
            if (masks[i].isNull()) {
                EXPECT_EQ(scene->style(i), nullptr) << i;
            } else {
                ASSERT_NE(scene->style(i), nullptr) << i;
                expectBools(scene->style(i)->mask(), masks[i], "mask");
            }
        }
        expectBools(scene->drawn(), at(c, "drawn"), "drawn");
        const std::vector<Vec3> residue_rgb = vec3s(at(c, "residue_rgb"));
        const std::vector<Colored> items = scene->frame(sys.frames.at(integer(at(c, "frame"))), residue_rgb);
        const Array& want_items = list(at(c, "items"));
        ASSERT_EQ(items.size(), want_items.size());
        for (std::size_t i = 0; i < items.size(); ++i) {
            const Value& w = want_items[i];
            const std::string what = "item " + std::to_string(i);
            const std::string kind = str(at(w, "kind"));
            EXPECT_EQ(items[i].opacity, number(at(w, "opacity"))) << what;
            expectVec3s(items[i].rgb, at(w, "rgb"), what + " rgb", 1e-12);
            if (const auto* m = std::get_if<Mesh>(&items[i].prim)) {
                ASSERT_EQ(kind, "mesh") << what;
                expectInts(m->owner, at(w, "owner"), what + " owner");
                expectInts(m->role, at(w, "role"), what + " role");
                EXPECT_EQ(static_cast<std::int64_t>(m->tris.size()), integer(at(w, "n_tris"))) << what;
                if (const Value* atom = w.find("atom")) {
                    expectInts(m->atom, *atom, what + " atom");
                }
            } else if (const auto* s = std::get_if<Spheres>(&items[i].prim)) {
                ASSERT_EQ(kind, "spheres") << what;
                expectInts(s->owner, at(w, "owner"), what + " owner");
                expectInts(s->atom, at(w, "atom"), what + " atom");
                expectDoubles(s->radius, at(w, "radius"), what + " radius", 1e-12);
            } else {
                const auto& st = std::get<Sticks>(items[i].prim);
                ASSERT_EQ(kind, "sticks") << what;
                expectInts(st.owner, at(w, "owner"), what + " owner");
                expectInts(st.atom, at(w, "atom"), what + " atom");
                expectDoubles(st.radius, at(w, "radius"), what + " radius", 1e-12);
            }
        }
        const std::vector<double> axis = doubles(at(c, "depth_axis"));
        const auto [lo, hi] = depthRange(items, {axis[0], axis[1], axis[2]});
        expectDoubles({lo, hi}, at(c, "depth_range"), "depth_range");
    }
}

TEST(MoleculeGeometryParity, NdimageMatchesScipy) {
    const Value& data = loadJson("ndimage.json");
    for (const Value& c : list(at(data, "filter1d"))) {
        const std::vector<double> in = doubles(at(c, "input"));
        std::vector<float> x(in.begin(), in.end());
        const double sigma = number(at(c, "sigma"));
        const double truncate = number(at(c, "truncate"));
        const auto radius = static_cast<std::size_t>(truncate * sigma + 0.5);
        expectDoubles(gaussianKernel1d(sigma, radius), at(c, "kernel"), "kernel", 1e-16);
        gaussianFilter1d(x, {1, 1, x.size()}, 2, sigma, truncate);
        expectDoubles(std::vector<double>(x.begin(), x.end()), at(c, "output"), "filter1d", 0.0);
    }
    for (const Value& c : list(at(data, "filter3d"))) {
        const std::vector<std::int64_t> s = ints(at(c, "shape"));
        const std::vector<double> in = doubles(at(c, "input"));
        std::vector<float> x(in.begin(), in.end());
        gaussianFilter(x, {static_cast<std::size_t>(s[0]), static_cast<std::size_t>(s[1]), static_cast<std::size_t>(s[2])},
                       number(at(c, "sigma")), 3.0);
        expectDoubles(std::vector<double>(x.begin(), x.end()), at(c, "output"), "filter3d", 0.0);
    }
    for (const Value& c : list(at(data, "map_coordinates"))) {
        const std::vector<std::int64_t> s = ints(at(c, "shape"));
        const std::vector<double> in = doubles(at(c, "input"));
        const std::vector<float> x(in.begin(), in.end());
        const std::vector<Vec3> pts = vec3s(at(c, "points"));
        std::vector<double> got;
        const Shape3 shape{static_cast<std::size_t>(s[0]), static_cast<std::size_t>(s[1]), static_cast<std::size_t>(s[2])};
        for (const Vec3& p : pts) {
            got.push_back(mapCoordinatesLinearNearest([&](std::size_t i) { return x[i]; }, shape, p[0], p[1], p[2]));
        }
        expectDoubles(got, at(c, "output"), "map_coordinates", 0.0);
    }
}

// ---- MDBridge's own --------------------------------------------------------------------------------

TEST(MoleculeGeometry, NearestPointsAgreeWithBruteForce) {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-15.0, 15.0);
    std::vector<Vec3> pts(400), queries(500);
    for (Vec3& p : pts) {
        p = {u(rng), u(rng), u(rng)};
    }
    pts[17] = pts[3];  // duplicates: the lower index wins
    for (Vec3& q : queries) {
        q = {u(rng) * 1.3, u(rng), u(rng) * 0.5};
    }
    queries[0] = pts[3];
    const std::vector<std::int32_t> near = nearestPoints(pts, queries);
    for (std::size_t q = 0; q < queries.size(); ++q) {
        std::size_t best = 0;
        for (std::size_t i = 1; i < pts.size(); ++i) {
            if (norm(sub(pts[i], queries[q])) < norm(sub(pts[best], queries[q]))) {
                best = i;
            }
        }
        EXPECT_EQ(near[q], static_cast<std::int32_t>(best)) << q;
    }
    EXPECT_EQ(near[0], 3);
}

TEST(MoleculeGeometry, SurfaceOfAnIsolatedAtomIsASphere) {
    const std::vector<Vec3> xyz{{0.3, -1.2, 2.0}};
    const GaussianSurface s = gaussianSurface(xyz);
    ASSERT_FALSE(s.tris.empty());
    for (std::size_t i = 0; i < s.verts.size(); ++i) {
        const double r = norm(sub(s.verts[i], xyz[0]));
        EXPECT_NEAR(r, 1.55, 0.15);  // the level-0.5 shell of a unit Gaussian (1.65 Å) on a 1 Å grid
        EXPECT_GT(dot(s.normals[i], sub(s.verts[i], xyz[0])), 0.0);
        EXPECT_EQ(s.nearest[i], 0);
    }
}

// RN-19 (docs/DIFFERENCES.md): a surface grid of more than 2^26 points (a blown-up frame) is
// refused before anything is allocated; the scene leaves the surface out of that frame with one
// warning, and draws the other representations.
TEST(MoleculeGeometry, OversizedSurfaceGridIsRefused) {
    const std::vector<Vec3> far{{0.0, 0.0, 0.0}, {5000.0, 0.0, 0.0}, {0.0, 3000.0, 1000.0}};
    try {
        (void)densityGrid(far);
        FAIL() << "no error";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("points is too large (at most 67108864 points)"), std::string::npos)
            << e.what();
    }
    auto topo = topology("demo_dry");
    std::vector<Vec3> xyz = system("demo_dry").frames.at(0);
    xyz[7] = {1e6, -1e6, 1e6};
    core::RepSpec surf, vdw;
    surf.style = "surf";
    vdw.style = "vdw";
    MoleculeScene scene(topo, {surf, vdw}, {std::vector<bool>(24, true), std::vector<bool>(24, true)});
    std::vector<std::string> warnings;
    const core::LogSink previous = core::setLogSink([&](core::LogLevel level, std::string_view message) {
        if (level == core::LogLevel::Warning) {
            warnings.emplace_back(message);
        }
    });
    const std::vector<Vec3> rgb(24, Vec3{0.5, 0.5, 0.5});
    const std::vector<Colored> first = scene.frame(xyz, rgb);
    const std::vector<Colored> second = scene.frame(xyz, rgb);
    core::setLogSink(previous);
    ASSERT_EQ(first.size(), 1u);  // the vdw spheres only
    EXPECT_TRUE(std::holds_alternative<Spheres>(first[0].prim));
    EXPECT_EQ(second.size(), 1u);
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_EQ(warnings[0].rfind("3D representation 1 (surf) not drawn on some frames: surface grid of ", 0), 0u)
        << warnings[0];
    // the normal frame still has its surface
    EXPECT_EQ(scene.frame(system("demo_dry").frames.at(0), rgb).size(), 2u);
}

}  // namespace
}  // namespace mdbridge::render::molecule
