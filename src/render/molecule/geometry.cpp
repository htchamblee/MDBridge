// Port of MDATK's mdatk/render/molecule/geometry.py (MDATK v0.0.4.1, the owner's code).

#include "mdbridge/render/molecule/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "mdbridge/render/molecule/surface.hpp"

namespace mdbridge::render::molecule {
namespace {

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 scale(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Vec3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Vec3 neg(const Vec3& a) { return {-a[0], -a[1], -a[2]}; }

std::size_t at(std::int32_t i) { return static_cast<std::size_t>(i); }

// The Catmull-Rom basis rows of one interval at t = j / n_sub: position and derivative weights.
void basisWeights(std::size_t n_sub, std::vector<std::array<double, 4>>& w, std::vector<std::array<double, 4>>& dw,
                  std::vector<double>& ts) {
    static constexpr double basis[4][4] = {{0.0, 1.0, 0.0, 0.0},
                                           {-0.5, 0.0, 0.5, 0.0},
                                           {1.0, -2.5, 2.0, -0.5},
                                           {-0.5, 1.5, -1.5, 0.5}};
    w.assign(n_sub, {});
    dw.assign(n_sub, {});
    ts.assign(n_sub, 0.0);
    for (std::size_t s = 0; s < n_sub; ++s) {
        const double t = static_cast<double>(s) / static_cast<double>(n_sub);
        ts[s] = t;
        const double pw[4] = {1.0, t, t * t, t * t * t};
        const double pd[4] = {0.0, 1.0, 2.0 * t, 3.0 * t * t};
        for (std::size_t j = 0; j < 4; ++j) {
            double a = 0.0;
            double d = 0.0;
            for (std::size_t i = 0; i < 4; ++i) {
                a += pw[i] * basis[i][j];
                d += pd[i] * basis[i][j];
            }
            w[s][j] = a;
            dw[s][j] = d;
        }
    }
}

}  // namespace

void DrawList::extend(DrawList other) {
    for (auto& m : other.meshes) {
        meshes.push_back(std::move(m));
    }
    for (auto& s : other.spheres) {
        spheres.push_back(std::move(s));
    }
    for (auto& s : other.sticks) {
        sticks.push_back(std::move(s));
    }
}

// ---- splines and swept profiles -------------------------------------------------------------------

DenseSpline catmullRom(std::size_t m, std::size_t n_sub) {
    if (m < 2 || n_sub < 1) {
        throw std::invalid_argument("catmullRom needs at least 2 control points and 1 sample per interval");
    }
    std::vector<std::array<double, 4>> w, dw;
    std::vector<double> ts;
    basisWeights(n_sub, w, dw, ts);
    const std::size_t ns = (m - 1) * n_sub + 1;
    const std::size_t me = m + 2;
    std::vector<double> a(ns * me, 0.0), d(ns * me, 0.0);
    DenseSpline out;
    out.samples = ns;
    out.m = m;
    out.owner.assign(ns, 0);
    for (std::size_t i = 0; i + 1 < m; ++i) {
        for (std::size_t s = 0; s < n_sub; ++s) {
            const std::size_t r = i * n_sub + s;
            for (std::size_t j = 0; j < 4; ++j) {
                a[r * me + i + j] = w[s][j];
                d[r * me + i + j] = dw[s][j];
            }
            out.owner[r] = static_cast<std::int32_t>(ts[s] < 0.5 ? i : i + 1);
        }
    }
    a[(ns - 1) * me + m] = 1.0;
    d[(ns - 1) * me + m - 1] = -0.5;
    d[(ns - 1) * me + m + 1] = 0.5;
    out.owner[ns - 1] = static_cast<std::int32_t>(m - 1);
    // phantom end points: p[-1] = 2 p0 - p1, p[m] = 2 p[m-1] - p[m-2]
    auto fold = [&](const std::vector<double>& x, std::vector<double>& y) {
        y.assign(ns * m, 0.0);
        for (std::size_t r = 0; r < ns; ++r) {
            const double* row = &x[r * me];
            double* dst = &y[r * m];
            for (std::size_t c = 0; c < m; ++c) {
                dst[c] += row[c + 1];
            }
            dst[0] += 2.0 * row[0];
            dst[1] -= row[0];
            dst[m - 1] += 2.0 * row[m + 1];
            dst[m - 2] -= row[m + 1];
        }
    };
    fold(a, out.a);
    fold(d, out.d);
    return out;
}

BandedSpline catmullRomBanded(std::size_t m, std::size_t n_sub) {
    if (m < 2 || n_sub < 1) {
        throw std::invalid_argument("catmullRomBanded needs at least 2 control points and 1 sample per interval");
    }
    std::vector<std::array<double, 4>> w, dw;
    std::vector<double> ts;
    basisWeights(n_sub, w, dw, ts);
    BandedSpline out;
    const std::size_t ns = (m - 1) * n_sub + 1;
    out.idx.reserve(ns);
    out.w.reserve(ns);
    out.d.reserve(ns);
    out.owner.reserve(ns);
    for (std::size_t i = 0; i + 1 < m; ++i) {
        const auto ii = static_cast<std::int32_t>(i);
        for (std::size_t s = 0; s < n_sub; ++s) {
            out.idx.push_back({ii, ii + 1, ii + 2, ii + 3});
            out.w.push_back(w[s]);
            out.d.push_back(dw[s]);
            out.owner.push_back(ts[s] < 0.5 ? ii : ii + 1);
        }
    }
    const auto mm = static_cast<std::int32_t>(m);
    out.idx.push_back({mm - 1, mm, mm + 1, mm + 1});
    out.w.push_back({0.0, 1.0, 0.0, 0.0});
    out.d.push_back({-0.5, 0.0, 0.5, 0.0});
    out.owner.push_back(mm - 1);
    return out;
}

std::vector<Vec3> extendControls(std::span<const Vec3> p) {
    const std::size_t m = p.size();
    std::vector<Vec3> out;
    out.reserve(m + 2);
    out.push_back(sub(scale(2.0, p[0]), p[1]));
    out.insert(out.end(), p.begin(), p.end());
    out.push_back(sub(scale(2.0, p[m - 1]), p[m - 2]));
    return out;
}

Profile profile(std::string_view kind, double half_width, double half_height, std::size_t k) {
    const double w = half_width;
    const double h = half_height;
    Profile out;
    if (kind == "ellipse") {
        const double step = 2.0 * std::numbers::pi / static_cast<double>(k);
        for (std::size_t i = 0; i < k; ++i) {
            const double ang = static_cast<double>(i) * step;
            const double c = std::cos(ang);
            const double s = std::sin(ang);
            out.pts.push_back({w * c, h * s});
            const double nx = h * c;
            const double ny = w * s;
            const double len = std::sqrt(nx * nx + ny * ny);
            out.normals.push_back({nx / len, ny / len});
        }
        return out;
    }
    if (kind == "band") {  // flat ribbon with duplicated corners for crisp edges
        out.pts = {{-w, -h}, {w, -h}, {w, -h}, {w, h}, {w, h}, {-w, h}, {-w, h}, {-w, -h}};
        out.normals = {{0, -1}, {0, -1}, {1, 0}, {1, 0}, {0, 1}, {0, 1}, {-1, 0}, {-1, 0}};
        return out;
    }
    throw std::invalid_argument("unknown profile '" + std::string(kind) + "'");
}

std::vector<Tri> sweepFaces(std::size_t ns, std::size_t k, bool caps) {
    std::vector<Tri> tris;
    const std::size_t rings = ns > 0 ? ns - 1 : 0;
    tris.reserve(2 * rings * k + (caps ? 2 * k : 0));
    auto v = [](std::size_t x) { return static_cast<std::int32_t>(x); };
    for (int pass = 0; pass < 2; ++pass) {
        for (std::size_t i = 0; i < rings; ++i) {
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t jn = (j + 1) % k;
                const std::size_t a = i * k + j;
                const std::size_t b = i * k + jn;
                const std::size_t c = (i + 1) * k + jn;
                const std::size_t d = (i + 1) * k + j;
                tris.push_back(pass == 0 ? Tri{v(a), v(c), v(b)} : Tri{v(a), v(d), v(c)});
            }
        }
    }
    if (caps) {
        for (std::size_t j = 0; j < k; ++j) {
            tris.push_back({v(ns * k), v(j), v((j + 1) % k)});
        }
        const std::size_t last = (ns - 1) * k;
        for (std::size_t j = 0; j < k; ++j) {
            tris.push_back({v(ns * k + 1), v(last + (j + 1) % k), v(last + j)});
        }
    }
    return tris;
}

Frames framesAlong(std::span<const Vec3> pos, std::span<const Vec3> tan, const std::vector<Vec3>* side) {
    const std::size_t ns = pos.size();
    Frames f;
    f.t.resize(ns);
    for (std::size_t i = 0; i < ns; ++i) {
        const double len = norm(tan[i]) + 1e-12;
        f.t[i] = {tan[i][0] / len, tan[i][1] / len, tan[i][2] / len};
    }
    if (side != nullptr) {
        std::vector<Vec3> n(ns);
        bool bad = false;
        for (std::size_t i = 0; i < ns; ++i) {
            const Vec3& s = (*side)[i];
            n[i] = sub(s, scale(dot(s, f.t[i]), f.t[i]));
            if (norm(n[i]) < 1e-6) {
                bad = true;
            }
        }
        if (!bad) {
            for (std::size_t i = 0; i < ns; ++i) {
                const double len = norm(n[i]);
                n[i] = {n[i][0] / len, n[i][1] / len, n[i][2] / len};
            }
            f.n = std::move(n);
            f.b.resize(ns);
            for (std::size_t i = 0; i < ns; ++i) {
                f.b[i] = cross(f.t[i], f.n[i]);
            }
            return f;
        }
    }
    // rotation-minimising frame
    f.n.assign(ns, {0.0, 0.0, 0.0});
    if (ns > 0) {
        const Vec3 ref = std::abs(f.t[0][2]) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
        Vec3 n0 = sub(ref, scale(dot(ref, f.t[0]), f.t[0]));
        const double len = norm(n0);
        f.n[0] = {n0[0] / len, n0[1] / len, n0[2] / len};
        for (std::size_t i = 1; i < ns; ++i) {
            const Vec3 v = sub(f.n[i - 1], scale(dot(f.n[i - 1], f.t[i]), f.t[i]));
            const double nv = norm(v);
            f.n[i] = nv > 1e-9 ? Vec3{v[0] / nv, v[1] / nv, v[2] / nv} : f.n[i - 1];
        }
    }
    f.b.resize(ns);
    for (std::size_t i = 0; i < ns; ++i) {
        f.b[i] = cross(f.t[i], f.n[i]);
    }
    return f;
}

RingSlabs ringSlabs(std::span<const Vec3> ring_points, std::size_t k, double half) {
    RingSlabs out;
    if (k == 0 || ring_points.empty()) {
        return out;
    }
    const std::size_t n = ring_points.size() / k;
    const std::size_t m = 6 * k + 2;
    out.verts.reserve(n * m);
    out.normals.reserve(n * m);
    out.slab.reserve(n * m);
    out.tris.reserve(n * 4 * k);
    for (std::size_t s = 0; s < n; ++s) {
        const Vec3* ring = ring_points.data() + s * k;
        Vec3 centre{0.0, 0.0, 0.0};
        if (k == 9) {  // purine: fan from the C4-C5 bond midpoint, which sees every perimeter edge
            centre = scale(0.5, add(ring[3], ring[8]));
        } else {
            for (std::size_t j = 0; j < k; ++j) {
                centre = add(centre, ring[j]);
            }
            centre = {centre[0] / static_cast<double>(k), centre[1] / static_cast<double>(k),
                      centre[2] / static_cast<double>(k)};
        }
        Vec3 nrm{0.0, 0.0, 0.0};  // Newell normal
        for (std::size_t j = 0; j < k; ++j) {
            nrm = add(nrm, cross(ring[j], ring[(j + 1) % k]));
        }
        const double len = norm(nrm) + 1e-12;
        nrm = {nrm[0] / len, nrm[1] / len, nrm[2] / len};
        const Vec3 h = scale(half, nrm);
        std::vector<Vec3> top(k), bot(k), side_n(k);
        for (std::size_t j = 0; j < k; ++j) {
            top[j] = add(ring[j], h);
            bot[j] = sub(ring[j], h);
            const Vec3 sn = cross(sub(ring[(j + 1) % k], ring[j]), nrm);
            const double sl = norm(sn) + 1e-12;
            side_n[j] = {sn[0] / sl, sn[1] / sl, sn[2] / sl};
        }
        out.verts.push_back(add(centre, h));
        out.verts.insert(out.verts.end(), top.begin(), top.end());
        out.verts.push_back(sub(centre, h));
        out.verts.insert(out.verts.end(), bot.begin(), bot.end());
        for (std::size_t j = 0; j < k; ++j) {
            const std::size_t jn = (j + 1) % k;
            out.verts.insert(out.verts.end(), {top[j], top[jn], bot[jn], bot[j]});
        }
        for (std::size_t j = 0; j < k + 1; ++j) {
            out.normals.push_back(nrm);
        }
        for (std::size_t j = 0; j < k + 1; ++j) {
            out.normals.push_back(neg(nrm));
        }
        for (std::size_t j = 0; j < k; ++j) {
            out.normals.insert(out.normals.end(), {side_n[j], side_n[j], side_n[j], side_n[j]});
        }
        const auto o = static_cast<std::int32_t>(s * m);
        const auto kk = static_cast<std::int32_t>(k);
        for (std::int32_t j = 0; j < kk; ++j) {
            out.tris.push_back({o, o + 1 + j, o + 1 + (j + 1) % kk});
        }
        for (std::int32_t j = 0; j < kk; ++j) {
            out.tris.push_back({o + kk + 1, o + kk + 2 + (j + 1) % kk, o + kk + 2 + j});
        }
        for (std::int32_t j = 0; j < kk; ++j) {
            const std::int32_t base = o + 2 * kk + 2 + 4 * j;
            out.tris.push_back({base, base + 2, base + 1});
        }
        for (std::int32_t j = 0; j < kk; ++j) {
            const std::int32_t base = o + 2 * kk + 2 + 4 * j;
            out.tris.push_back({base, base + 3, base + 2});
        }
        out.slab.insert(out.slab.end(), m, static_cast<std::int32_t>(s));
    }
    return out;
}

// ---- Backbone -------------------------------------------------------------------------------------

Backbone::Backbone(const Topology& topo, const std::vector<bool>& res_mask, std::string_view kind, double half_width,
                   double half_height, std::size_t n_sub, std::size_t k, std::optional<double> round_radius)
    : topo_(&topo), flat_(profile(kind, half_width, half_height, k)) {
    const bool has_round = round_radius.has_value();
    if (has_round) {
        round_ = profile("ellipse", *round_radius, *round_radius, k);
    }
    const auto nr = static_cast<std::size_t>(topo.n_res);
    drawn_.assign(nr, false);
    std::vector<bool> has_base(nr, false);
    for (std::int32_t r : topo.nucleic()) {
        has_base[at(r)] = true;
    }
    auto selected = [&](std::int32_t r) { return at(r) < res_mask.size() && res_mask[at(r)]; };
    for (const auto& seg : topo.segments) {
        if (std::none_of(seg.begin(), seg.end(), selected)) {
            continue;
        }
        const bool rnd = has_round && std::none_of(seg.begin(), seg.end(), [&](std::int32_t r) { return has_base[at(r)]; });
        const BandedSpline spline = catmullRomBanded(seg.size(), n_sub);
        const std::size_t ns = spline.owner.size();
        std::vector<std::size_t> keep;
        for (std::size_t s = 0; s < ns; ++s) {
            if (selected(seg[at(spline.owner[s])])) {
                keep.push_back(s);
            }
        }
        std::size_t start = 0;
        while (start < keep.size()) {
            std::size_t stop = start + 1;
            while (stop < keep.size() && keep[stop] - keep[stop - 1] == 1) {
                ++stop;
            }
            std::vector<std::size_t> run(keep.begin() + static_cast<std::ptrdiff_t>(start),
                                         keep.begin() + static_cast<std::ptrdiff_t>(stop));
            start = stop;
            std::vector<std::int32_t> own;
            own.reserve(run.size() + 1);
            for (std::size_t s : run) {
                own.push_back(seg[at(spline.owner[s])]);
            }
            if (run.back() + 1 < ns) {  // meet the next stretch without a gap
                run.push_back(run.back() + 1);
                own.push_back(own.back());
            }
            if (run.size() < 2) {
                continue;
            }
            Part part;
            part.seg = seg;
            part.round = rnd;
            for (std::size_t s : run) {
                part.idx.push_back(spline.idx[s]);
                part.w.push_back(spline.w[s]);
                part.d.push_back(spline.d[s]);
            }
            for (std::size_t q = 0; q < run.size(); ++q) {
                if (run[q] % n_sub == 0) {
                    const std::int32_t r = seg[run[q] / n_sub];
                    if (selected(r)) {  // samples at the selected control atoms
                        part.at_res.push_back(r);
                        part.at.push_back(q);
                    }
                }
            }
            for (std::int32_t r : own) {
                drawn_[at(r)] = true;
            }
            part.own = std::move(own);
            parts_.push_back(std::move(part));
        }
    }
}

std::pair<std::vector<Mesh>, std::map<std::int32_t, Vec3>> Backbone::build(std::span<const Vec3> xyz) const {
    const Topology& topo = *topo_;
    std::vector<Mesh> meshes;
    std::map<std::int32_t, Vec3> anchor;
    for (const Part& part : parts_) {
        const std::size_t m = part.seg.size();
        std::vector<Vec3> p(m);
        for (std::size_t i = 0; i < m; ++i) {
            p[i] = xyz[at(topo.control[at(part.seg[i])])];
        }
        const std::vector<Vec3> ext = extendControls(p);
        const std::size_t ns = part.idx.size();
        std::vector<Vec3> pos(ns), tan(ns);
        auto combine = [](const std::array<double, 4>& w, const std::vector<Vec3>& e,
                          const std::array<std::int32_t, 4>& idx) {
            Vec3 out{0.0, 0.0, 0.0};
            for (std::size_t c = 0; c < 3; ++c) {
                double acc = 0.0;
                for (std::size_t j = 0; j < 4; ++j) {
                    acc += w[j] * e[at(idx[j])][c];
                }
                out[c] = acc;
            }
            return out;
        };
        for (std::size_t s = 0; s < ns; ++s) {
            pos[s] = combine(part.w[s], ext, part.idx[s]);
            tan[s] = combine(part.d[s], ext, part.idx[s]);
        }
        std::vector<Vec3> side;
        bool use_side = !part.round;
        if (use_side) {
            for (std::int32_t r : part.seg) {
                if (topo.side[at(r)] < 0) {
                    use_side = false;
                    break;
                }
            }
        }
        if (use_side) {
            std::vector<Vec3> rel(m);
            for (std::size_t i = 0; i < m; ++i) {
                rel[i] = sub(xyz[at(topo.side[at(part.seg[i])])], p[i]);
            }
            const std::vector<Vec3> ext_side = extendControls(rel);
            side.resize(ns);
            for (std::size_t s = 0; s < ns; ++s) {
                side[s] = combine(part.w[s], ext_side, part.idx[s]);
            }
        }
        const Frames f = framesAlong(pos, tan, use_side ? &side : nullptr);
        const Profile& prof = part.round ? round_ : flat_;
        const std::size_t k = prof.pts.size();
        Mesh mesh;
        mesh.verts.reserve(ns * k + 2);
        mesh.normals.reserve(ns * k + 2);
        mesh.owner.reserve(ns * k + 2);
        for (std::size_t s = 0; s < ns; ++s) {
            for (std::size_t j = 0; j < k; ++j) {
                mesh.verts.push_back(add(add(pos[s], scale(prof.pts[j][0], f.b[s])), scale(prof.pts[j][1], f.n[s])));
                mesh.normals.push_back(add(scale(prof.normals[j][0], f.b[s]), scale(prof.normals[j][1], f.n[s])));
                mesh.owner.push_back(part.own[s]);
            }
        }
        mesh.verts.push_back(pos.front());
        mesh.verts.push_back(pos.back());
        mesh.normals.push_back(neg(f.t.front()));
        mesh.normals.push_back(f.t.back());
        mesh.owner.push_back(part.own.front());
        mesh.owner.push_back(part.own.back());
        mesh.tris = sweepFaces(ns, k);
        mesh.role.assign(mesh.verts.size(), kRoleResidue);
        meshes.push_back(std::move(mesh));
        for (std::size_t q = 0; q < part.at.size(); ++q) {
            anchor[part.at_res[q]] = pos[part.at[q]];
        }
    }
    return {std::move(meshes), std::move(anchor)};
}

// ---- styles ---------------------------------------------------------------------------------------

Style::Style(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options)
    : topo_(std::move(topo)), mask_(std::move(mask)), options_(options) {
    mask_.resize(static_cast<std::size_t>(topo_->n_res), false);
}

std::vector<std::int32_t> Style::selectedAtoms() const {
    std::vector<std::int32_t> out;
    for (std::size_t i = 0; i < topo_->row.size(); ++i) {
        const std::int32_t r = topo_->row[i];
        if (r >= 0 && r < topo_->n_res && mask_[at(r)]) {
            out.push_back(static_cast<std::int32_t>(i));
        }
    }
    return out;
}

std::vector<bool> Style::drawn() const {
    std::vector<bool> out(mask_.size(), false);
    for (std::int32_t a : selectedAtoms()) {
        out[at(topo_->row[at(a)])] = true;
    }
    return out;
}

namespace {

class TubeStyle : public Style {
public:
    TubeStyle(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options,
              std::string_view kind = "ellipse", double half_width = 0.45, double half_height = 0.45)
        : Style(std::move(topo), std::move(mask), options),
          backbone_(*topo_, mask_, kind, half_width * options.size, half_height * options.size, options.samples, 12,
                    0.45 * options.size) {}

    std::string_view name() const override { return "tube"; }
    DrawList build(std::span<const Vec3> xyz) const override {
        DrawList out;
        out.meshes = backbone_.build(xyz).first;
        return out;
    }
    std::vector<bool> drawn() const override {
        std::vector<bool> out = backbone_.drawn();
        for (std::size_t r = 0; r < out.size(); ++r) {
            out[r] = out[r] && mask_[r];
        }
        return out;
    }

private:
    Backbone backbone_;
};

class RibbonStyle : public TubeStyle {
public:
    RibbonStyle(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options)
        : TubeStyle(std::move(topo), std::move(mask), options, "band", 1.2, 0.12) {}
    std::string_view name() const override { return "ribbon"; }
};

// New Cartoon for nucleic acids: backbone + rungs + base slabs (+ optional sugar slabs).
class CartoonStyle : public Style {
public:
    CartoonStyle(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options)
        : Style(std::move(topo), std::move(mask), options),
          backbone_(*topo_, mask_, "ellipse", 1.6 * options.size, 0.45 * options.size, options.samples, 12,
                    0.6 * options.size) {
        const Topology& t = *topo_;
        for (std::size_t q = 0; q < t.purine_res.size(); ++q) {
            if (mask_[at(t.purine_res[q])]) {
                pur_.push_back(q);
            }
        }
        for (std::size_t q = 0; q < t.pyrimidine_res.size(); ++q) {
            if (mask_[at(t.pyrimidine_res[q])]) {
                pyr_.push_back(q);
            }
        }
        if (options.sugars) {
            for (std::size_t q = 0; q < t.sugar_res.size(); ++q) {
                if (mask_[at(t.sugar_res[q])]) {
                    sug_.push_back(q);
                }
            }
        }
        std::vector<bool> in_seg(mask_.size(), false);
        for (const auto& seg : t.segments) {
            for (std::int32_t r : seg) {
                in_seg[at(r)] = true;
            }
        }
        for (std::size_t r = 0; r < mask_.size(); ++r) {
            if (mask_[r] && !in_seg[r] && t.control[r] >= 0) {
                lone_.push_back(static_cast<std::int32_t>(r));
            }
        }
    }

    std::string_view name() const override { return "cartoon"; }

    std::vector<bool> drawn() const override {
        const Topology& t = *topo_;
        std::vector<bool> out = backbone_.drawn();
        for (std::int32_t r : lone_) {
            out[at(r)] = true;
        }
        for (std::size_t q : pur_) {
            out[at(t.purine_res[q])] = true;
        }
        for (std::size_t q : pyr_) {
            out[at(t.pyrimidine_res[q])] = true;
        }
        for (std::size_t q : sug_) {
            out[at(t.sugar_res[q])] = true;
        }
        for (std::size_t r = 0; r < out.size(); ++r) {
            out[r] = out[r] && mask_[r];
        }
        return out;
    }

    DrawList build(std::span<const Vec3> xyz) const override {
        const Topology& t = *topo_;
        const double s = options_.size;
        auto [meshes, anchor] = backbone_.build(xyz);
        DrawList out;
        out.meshes = std::move(meshes);
        Spheres joints;
        for (std::int32_t r : lone_) {  // no backbone neighbours: a ball at the control atom
            anchor[r] = xyz[at(t.control[at(r)])];
            joints.xyz.push_back(anchor[r]);
            joints.owner.push_back(r);
            joints.radius.push_back(0.9 * s);
        }
        // base slabs, coloured by nucleotide
        auto slabs = [&](const std::vector<std::size_t>& sel, const auto& rows, const auto& rings, double half,
                         Role role) {
            if (sel.empty()) {
                return;
            }
            const std::size_t k = rings.empty() ? 0 : rings.front().size();
            std::vector<Vec3> pts;
            pts.reserve(sel.size() * k);
            for (std::size_t q : sel) {
                for (std::int32_t a : rings[q]) {
                    pts.push_back(xyz[at(a)]);
                }
            }
            RingSlabs rs = ringSlabs(pts, k, half);
            Mesh mesh;
            mesh.owner.reserve(rs.slab.size());
            for (std::int32_t slab : rs.slab) {
                mesh.owner.push_back(rows[sel[at(slab)]]);
            }
            mesh.verts = std::move(rs.verts);
            mesh.normals = std::move(rs.normals);
            mesh.tris = std::move(rs.tris);
            mesh.role.assign(mesh.verts.size(), role);
            out.meshes.push_back(std::move(mesh));
        };
        slabs(pur_, t.purine_res, t.purine_ring, 0.45 * s, kRoleNucleotide);
        slabs(pyr_, t.pyrimidine_res, t.pyrimidine_ring, 0.45 * s, kRoleNucleotide);
        slabs(sug_, t.sugar_res, t.sugar_ring, 0.25 * s, kRoleSugar);
        // rungs: backbone spline point -> C1' (or C4' with sugar slabs) -> glycosidic atom
        Sticks rungs;
        auto stick = [&](const Vec3& a, const Vec3& b, std::int32_t r) {
            rungs.a.push_back(a);
            rungs.b.push_back(b);
            rungs.owner.push_back(r);
        };
        auto rung = [&](std::int32_t r, std::int32_t g) {
            const std::int32_t knee = options_.sugars ? t.c4[at(r)] : t.c1[at(r)];
            const auto it = anchor.find(r);
            if (it == anchor.end() || knee < 0) {
                return;
            }
            std::vector<Vec3> ends;
            // sticks are open cylinders: a ball closes every rung end and bend
            if (options_.sugars) {
                stick(it->second, xyz[at(knee)], r);
                ends.push_back(xyz[at(knee)]);
                const std::int32_t c1 = t.c1[at(r)];
                if (c1 >= 0) {
                    stick(xyz[at(c1)], xyz[at(g)], r);
                    ends.push_back(xyz[at(c1)]);
                    ends.push_back(xyz[at(g)]);
                }
            } else {
                stick(it->second, xyz[at(knee)], r);
                stick(xyz[at(knee)], xyz[at(g)], r);
                ends.push_back(xyz[at(knee)]);
                ends.push_back(xyz[at(g)]);
            }
            for (const Vec3& e : ends) {
                joints.xyz.push_back(e);
                joints.owner.push_back(r);
                joints.radius.push_back(0.42 * s);
            }
        };
        for (std::size_t q : pur_) {
            rung(t.purine_res[q], t.purine_gly[q]);
        }
        for (std::size_t q : pyr_) {
            rung(t.pyrimidine_res[q], t.pyrimidine_gly[q]);
        }
        if (!rungs.a.empty()) {
            rungs.radius.assign(rungs.a.size(), 0.42 * s);
            rungs.atom.assign(rungs.a.size(), -1);
            out.sticks.push_back(std::move(rungs));
        }
        if (!joints.xyz.empty()) {
            joints.atom.assign(joints.xyz.size(), -1);
            out.spheres.push_back(std::move(joints));
        }
        return out;
    }

private:
    Backbone backbone_;
    std::vector<std::size_t> pur_, pyr_, sug_;
    std::vector<std::int32_t> lone_;
};

class LicoriceStyle : public Style {
public:
    LicoriceStyle(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options)
        : Style(std::move(topo), std::move(mask), options), atoms_(selectedAtoms()), sel_(topo_->row.size(), false) {
        for (std::int32_t a : atoms_) {
            sel_[at(a)] = true;
        }
        // bonds with at least one selected atom are drawn from the selected side, so licorice joins
        // a neighbouring style (e.g. the cartoon at the phosphate)
        for (const auto& b : topo_->bonds) {
            if (sel_[at(b[0])] || sel_[at(b[1])]) {
                bonds_.push_back(b);
            }
        }
    }

    std::string_view name() const override { return "licorice"; }

    DrawList build(std::span<const Vec3> xyz) const override {
        const Topology& t = *topo_;
        const double r = 0.3 * options_.size;
        Sticks sticks;
        // half-bonds inside the selection; a bond leaving it runs to the neighbouring atom
        for (int pass = 0; pass < 2; ++pass) {
            for (const auto& bond : bonds_) {
                const std::int32_t i = pass == 0 ? bond[0] : bond[1];
                const std::int32_t j = pass == 0 ? bond[1] : bond[0];
                if (!sel_[at(i)]) {
                    continue;
                }
                const Vec3 mid = scale(0.5, add(xyz[at(bond[0])], xyz[at(bond[1])]));
                sticks.a.push_back(xyz[at(i)]);
                sticks.b.push_back(sel_[at(j)] ? mid : xyz[at(j)]);
                sticks.atom.push_back(i);
                sticks.owner.push_back(t.row[at(i)]);
            }
        }
        DrawList out;
        if (!sticks.a.empty()) {
            sticks.radius.assign(sticks.a.size(), r);
            out.sticks.push_back(std::move(sticks));
        }
        Spheres balls;
        for (std::int32_t a : atoms_) {
            balls.xyz.push_back(xyz[at(a)]);
            balls.owner.push_back(t.row[at(a)]);
            balls.atom.push_back(a);
        }
        balls.radius.assign(atoms_.size(), r);
        out.spheres.push_back(std::move(balls));
        return out;
    }

private:
    std::vector<std::int32_t> atoms_;
    std::vector<bool> sel_;
    std::vector<std::array<std::int32_t, 2>> bonds_;
};

class VdwStyle : public Style {
public:
    VdwStyle(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options)
        : Style(std::move(topo), std::move(mask), options), atoms_(selectedAtoms()) {}

    std::string_view name() const override { return "vdw"; }

    DrawList build(std::span<const Vec3> xyz) const override {
        const Topology& t = *topo_;
        Spheres balls;
        for (std::int32_t a : atoms_) {
            balls.xyz.push_back(xyz[at(a)]);
            balls.radius.push_back(t.radius[at(a)] * options_.size);
            balls.owner.push_back(t.row[at(a)]);
            balls.atom.push_back(a);
        }
        DrawList out;
        out.spheres.push_back(std::move(balls));
        return out;
    }

private:
    std::vector<std::int32_t> atoms_;
};

// Smooth molecular surface (like VMD QuickSurf): iso-surface of summed atomic Gaussians (surface
// nets, surface.hpp), coloured by the residue of the nearest atom.
class SurfStyle : public Style {
public:
    static constexpr double kSigma = 1.4;  // an isolated atom's surface lies ~1.65 Å from its centre
    static constexpr double kIso = 0.5;
    static constexpr std::size_t kLarge = 5000;  // atoms above which the grid coarsens to 1.5 Å

    SurfStyle(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options)
        : Style(std::move(topo), std::move(mask), options), atoms_(selectedAtoms()) {
        const double fallback = atoms_.size() > kLarge ? 1.5 : 1.0;
        spacing_ = options.spacing && *options.spacing != 0.0 ? *options.spacing : fallback;
    }

    std::string_view name() const override { return "surf"; }

    DrawList build(std::span<const Vec3> xyz) const override {
        std::vector<std::int32_t> finite;  // a blown-up frame must not break the video
        std::vector<Vec3> pts;
        for (std::int32_t a : atoms_) {
            const Vec3& p = xyz[at(a)];
            if (std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2])) {
                finite.push_back(a);
                pts.push_back(p);
            }
        }
        if (finite.empty()) {
            return {};
        }
        GaussianSurface surf = gaussianSurface(pts, spacing_, kSigma * options_.size, kIso);
        if (surf.tris.empty()) {
            return {};
        }
        Mesh mesh;
        mesh.verts = std::move(surf.verts);
        mesh.normals = std::move(surf.normals);
        mesh.tris = std::move(surf.tris);
        mesh.has_atom = true;
        mesh.atom.reserve(surf.nearest.size());
        mesh.owner.reserve(surf.nearest.size());
        for (std::int32_t q : surf.nearest) {
            const std::int32_t a = finite[at(q)];
            mesh.atom.push_back(a);
            mesh.owner.push_back(topo_->row[at(a)]);
        }
        mesh.role.assign(mesh.verts.size(), kRoleResidue);
        DrawList out;
        out.meshes.push_back(std::move(mesh));
        return out;
    }

private:
    std::vector<std::int32_t> atoms_;
    double spacing_ = 1.0;
};

// The 2D trace in the scene: a ball at each residue's backbone point and half-sticks to its
// neighbours along the backbone segments (so it hides, and is hidden by, the other
// representations by depth).
class TraceStyle : public Style {
public:
    TraceStyle(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options)
        : Style(std::move(topo), std::move(mask), options) {
        const Topology& t = *topo_;
        for (std::size_t r = 0; r < mask_.size(); ++r) {
            if (mask_[r] && t.anchor_atom[r] >= 0) {
                rows_.push_back(static_cast<std::int32_t>(r));
            }
        }
        for (const auto& seg : t.segments) {
            for (std::size_t q = 0; q + 1 < seg.size(); ++q) {
                if (mask_[at(seg[q])] || mask_[at(seg[q + 1])]) {
                    pairs_.push_back({seg[q], seg[q + 1]});
                }
            }
        }
    }

    std::string_view name() const override { return "trace"; }

    std::vector<bool> drawn() const override {
        std::vector<bool> out(mask_.size(), false);
        for (std::int32_t r : rows_) {
            out[at(r)] = true;
        }
        return out;
    }

    DrawList build(std::span<const Vec3> xyz) const override {
        const Topology& t = *topo_;
        const double s = options_.size;
        DrawList out;
        if (!pairs_.empty()) {
            Sticks sticks;
            for (int pass = 0; pass < 2; ++pass) {
                for (const auto& [i, j] : pairs_) {
                    const std::int32_t r = pass == 0 ? i : j;
                    if (!mask_[at(r)]) {
                        continue;
                    }
                    const Vec3& a = xyz[at(t.control[at(i)])];
                    const Vec3& b = xyz[at(t.control[at(j)])];
                    sticks.a.push_back(pass == 0 ? a : b);
                    sticks.b.push_back(scale(0.5, add(a, b)));
                    sticks.owner.push_back(r);
                }
            }
            sticks.radius.assign(sticks.a.size(), 0.38 * s);
            sticks.atom.assign(sticks.a.size(), -1);
            out.sticks.push_back(std::move(sticks));
        }
        Spheres balls;
        for (std::int32_t r : rows_) {
            balls.xyz.push_back(xyz[at(t.anchor_atom[at(r)])]);
            balls.owner.push_back(r);
        }
        balls.radius.assign(rows_.size(), 0.85 * s);
        balls.atom.assign(rows_.size(), -1);
        out.spheres.push_back(std::move(balls));
        return out;
    }

private:
    std::vector<std::int32_t> rows_;
    std::vector<std::pair<std::int32_t, std::int32_t>> pairs_;
};

}  // namespace

std::unique_ptr<Style> makeStyle(std::string_view name, std::shared_ptr<const Topology> topo, std::vector<bool> mask,
                                 const StyleOptions& options) {
    if (name == "cartoon") {
        return std::make_unique<CartoonStyle>(std::move(topo), std::move(mask), options);
    }
    if (name == "ribbon") {
        return std::make_unique<RibbonStyle>(std::move(topo), std::move(mask), options);
    }
    if (name == "tube") {
        return std::make_unique<TubeStyle>(std::move(topo), std::move(mask), options);
    }
    if (name == "licorice") {
        return std::make_unique<LicoriceStyle>(std::move(topo), std::move(mask), options);
    }
    if (name == "vdw") {
        return std::make_unique<VdwStyle>(std::move(topo), std::move(mask), options);
    }
    if (name == "surf") {
        return std::make_unique<SurfStyle>(std::move(topo), std::move(mask), options);
    }
    if (name == "trace") {
        return std::make_unique<TraceStyle>(std::move(topo), std::move(mask), options);
    }
    return nullptr;
}

const std::vector<std::string>& styleNames() {
    static const std::vector<std::string> names = {"cartoon", "ribbon", "tube", "licorice", "vdw", "surf", "trace"};
    return names;
}

}  // namespace mdbridge::render::molecule
