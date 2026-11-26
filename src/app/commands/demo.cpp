// The synthetic demonstration system: a port of MDATK's mdatk/demo.py (build_system,
// _bridge_site, write_topology, write_trajectory, demo_project, build_demo). MDATK's own code; the
// arithmetic follows numpy's evaluation order term by term (no fused multiply-adds: this file is
// compiled with -ffp-contract=off / /fp:precise) so that the float32 coordinates and the parm7
// numbers are MDATK's bit for bit. scipy's cKDTree queries are exact searches here (brute force
// over a few thousand points), which give the same sets and minimum distances.
#include "demo.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <string_view>
#include <utility>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/process.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/netcdf.hpp"
#include "mdbridge/io/parm7.hpp"
#include "mdbridge/render/numpy_random.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;

namespace {

using Vec = std::array<double, 3>;
using Mat = std::array<Vec, 3>;
using Vec2 = std::array<double, 2>;

constexpr double kPi = std::numbers::pi;  // numpy.pi
// np.radians(x): x * (pi / 180)
double radians(double deg) { return deg * (kPi / 180.0); }

constexpr double kRise = 3.0;
const double kTwist = radians(32.7);
constexpr double kXPyr = 4.046;
constexpr double kXPur = 1.634;
constexpr double kRingR = 1.39;
constexpr double kAmberChargeFactor = 18.2223;

struct LjType {
    std::string_view name;
    double rstar, epsilon;
};
// LJ_TYPES: (r*, epsilon) per atom type, in MDATK's order
constexpr LjType kLjTypes[] = {
    {"P", 2.1, 0.2},      {"O2", 1.6612, 0.21}, {"OS", 1.6837, 0.17},   {"OH", 1.721, 0.2104},
    {"O", 1.6612, 0.21},  {"CT", 1.908, 0.1094}, {"CA", 1.908, 0.086},  {"N", 1.824, 0.17},
    {"H1", 1.387, 0.0157}, {"HA", 1.359, 0.015}, {"H", 0.6, 0.0157},    {"HO", 0.2, 0.03},
    {"MG", 0.7926, 0.8947}, {"Na+", 1.369, 0.0874}, {"OW", 1.7683, 0.152}, {"HW", 0.0, 0.0},
};

std::int64_t elementZ(std::string_view el) {
    static const std::map<std::string_view, std::int64_t> z{{"H", 1}, {"C", 6}, {"N", 7}, {"O", 8},
                                                             {"P", 15}, {"MG", 12}, {"NA", 11}};
    return z.at(el);
}
double elementMass(std::string_view el) {
    static const std::map<std::string_view, double> m{{"H", 1.008}, {"C", 12.01}, {"N", 14.01}, {"O", 16.0},
                                                       {"P", 30.97}, {"MG", 24.305}, {"NA", 22.99}};
    return m.at(el);
}

struct Offset {
    std::string_view name;
    Vec xyz;
    std::string_view type;
    double charge;
};
// SUGAR offsets from C1' (pair frame, strand A side), SUGAR_TYPES and SUGAR_CHARGES
constexpr Offset kSugar[] = {
    {"O4'", {-0.78, 1.08, 0.35}, "OS", -0.355},  {"C2'", {-0.58, -1.07, 0.85}, "CT", 0.067},
    {"O2'", {-0.13, -2.37, 0.55}, "OH", -0.614}, {"C3'", {-1.93, -0.52, 1.25}, "CT", 0.202},
    {"C4'", {-2.08, 0.93, 0.85}, "CT", 0.107},   {"O3'", {-2.83, -1.37, 2.00}, "OS", -0.525},
    {"C5'", {-3.28, 1.93, 0.40}, "CT", 0.056},   {"O5'", {-4.18, 0.93, 0.95}, "OS", -0.499},
    {"P", {-4.48, 2.03, -0.35}, "P", 1.166},     {"OP1", {-5.53, 2.78, -0.95}, "O2", -0.776},
    {"OP2", {-5.18, 0.98, -1.15}, "O2", -0.776},
};
constexpr double kC1Charge = 0.039;

struct SugarH {
    std::string_view name, parent;
    Vec hint;
};
constexpr SugarH kSugarH[] = {
    {"H1'", "C1'", {0.0, -0.3, -1.0}}, {"H2'", "C2'", {0.2, -0.4, 1.0}}, {"HO2'", "O2'", {0.8, -0.5, -0.3}},
    {"H3'", "C3'", {-0.2, 0.0, 1.0}},  {"H4'", "C4'", {0.0, 0.6, 1.0}},  {"H5'", "C5'", {-0.3, 1.0, -0.6}},
    {"H5''", "C5'", {-0.6, 0.6, 1.0}},
};

struct Atom {
    std::string name;
    Vec xyz{};
    std::string el;
    std::string type;
    double q = 0.0;
};

// ---- small vector algebra in numpy's order ----------------------------------------------------

Vec add(const Vec& a, const Vec& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec sub(const Vec& a, const Vec& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec scale(double s, const Vec& a) { return {s * a[0], s * a[1], s * a[2]}; }
// x.dot(x) of three values, summed in order
double dot(const Vec& a, const Vec& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
// np.linalg.norm: sqrt(x.dot(x))
double norm(const Vec& a) { return std::sqrt(dot(a, a)); }
// _unit(v) = v / norm(v)
Vec unit(const Vec& v) {
    const double n = norm(v);
    return {v[0] / n, v[1] / n, v[2] / n};
}
// np.cross of 3-vectors
Vec cross(const Vec& a, const Vec& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
// m @ v
Vec mul(const Mat& m, const Vec& v) { return {dot(m[0], v), dot(m[1], v), dot(m[2], v)}; }

Vec2 polar(const Vec2& center, double radius, double deg) {
    const double a = radians(deg);
    return {center[0] + radius * std::cos(a), center[1] + radius * std::sin(a)};
}
Vec2 rot2(const Vec2& v, double deg) {
    const double a = radians(deg);
    const double c = std::cos(a), s = std::sin(a);
    return {c * v[0] - s * v[1], s * v[0] + c * v[1]};
}
Vec2 add2(const Vec2& a, const Vec2& b) { return {a[0] + b[0], a[1] + b[1]}; }
Vec2 scale2(double s, const Vec2& a) { return {s * a[0], s * a[1]}; }

double smoothStep(double x) {
    x = std::clamp(x, 0.0, 1.0);  // np.clip
    return x * x * (3 - 2 * x);
}

// _level_frame(k): rotation about z by k * TWIST, shift k * RISE along z
Mat levelRotation(std::int64_t k) {
    const double t = static_cast<double>(k) * kTwist;
    return Mat{Vec{std::cos(t), -std::sin(t), 0.0}, Vec{std::sin(t), std::cos(t), 0.0}, Vec{0.0, 0.0, 1.0}};
}

// _rotation(axis, angle)
Mat rotation(const Vec& axisIn, double angle) {
    const Vec axis = unit(axisIn);
    const double x = axis[0], y = axis[1], z = axis[2];
    const double c = std::cos(angle), s = std::sin(angle);
    const double C = 1 - c;
    return Mat{Vec{c + x * x * C, x * y * C - z * s, x * z * C + y * s},
               Vec{y * x * C + z * s, c + y * y * C, y * z * C - x * s},
               Vec{z * x * C - y * s, z * y * C + x * s, c + z * z * C}};
}

// ---- nucleotides ----------------------------------------------------------------------------------

struct BaseAtom {
    std::string name;
    Vec2 xy;
    std::string el, type;
    double q;
};

// base_template(base): base atoms (strand A side) and the C1' position
std::pair<std::vector<BaseAtom>, Vec2> baseTemplate(char base) {
    std::vector<BaseAtom> atoms;
    const auto ringAtom = [](const std::string& n) {
        return std::pair<std::string, double>(n[0] == 'N' ? "N" : "CA", n[0] == 'N' ? -0.45 : 0.15);
    };
    Vec2 c1{};
    if (base == 'C' || base == 'U') {
        const Vec2 c{-kXPyr, 0.0};
        const std::map<std::string, double> ring{{"N3", 0}, {"C4", 60}, {"C5", 120}, {"C6", 180}, {"N1", 240}, {"C2", 300}};
        std::map<std::string, Vec2> pos;
        for (const auto& [n, d] : ring) {
            pos[n] = polar(c, kRingR, d);
        }
        for (const char* n : {"N1", "C2", "N3", "C4", "C5", "C6"}) {
            const auto [type, q] = ringAtom(n);
            atoms.push_back({n, pos[n], std::string(1, n[0]), type, q});
        }
        atoms.push_back({"O2", polar(c, kRingR + 1.23, 300), "O", "O", -0.58});
        atoms.push_back({"H5", polar(c, kRingR + 1.08, 120), "H", "HA", 0.15});
        atoms.push_back({"H6", polar(c, kRingR + 1.08, 180), "H", "HA", 0.15});
        if (base == 'U') {
            atoms.push_back({"O4", polar(c, kRingR + 1.23, 60), "O", "O", -0.58});
            atoms.push_back({"H3", add2(pos["N3"], {1.01, 0.0}), "H", "H", 0.40});
        } else {
            const Vec2 n4 = polar(c, kRingR + 1.34, 60);
            atoms.push_back({"N4", n4, "N", "N", -0.85});
            atoms.push_back({"H41", add2(n4, {1.01, 0.0}), "H", "H", 0.42});
            atoms.push_back({"H42", add2(n4, scale2(1.01, rot2({1.0, 0.0}, 120))), "H", "H", 0.42});
        }
        c1 = polar(c, kRingR + 1.47, 240);
    } else {
        const Vec2 c{-kXPur, 0.0};
        const std::map<std::string, double> ring{{"N1", 0}, {"C6", 60}, {"C5", 120}, {"C4", 180}, {"N3", 240}, {"C2", 300}};
        std::map<std::string, Vec2> pos;
        for (const auto& [n, d] : ring) {
            pos[n] = polar(c, kRingR, d);
        }
        const Vec2 c5 = add2(c, {-1.871, 1.080});
        const double r5 = 1.1824;
        pos["N7"] = polar(c5, r5, 78);
        pos["C8"] = polar(c5, r5, 150);
        pos["N9"] = polar(c5, r5, 222);
        for (const char* n : {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"}) {
            const auto [type, q] = ringAtom(n);
            atoms.push_back({n, pos[n], std::string(1, n[0]), type, q});
        }
        atoms.push_back({"H8", polar(c5, r5 + 1.08, 150), "H", "HA", 0.15});
        if (base == 'A') {
            const Vec2 n6 = polar(c, kRingR + 1.34, 60);
            atoms.push_back({"N6", n6, "N", "N", -0.85});
            atoms.push_back({"H61", add2(n6, {1.01, 0.0}), "H", "H", 0.42});
            atoms.push_back({"H62", add2(n6, scale2(1.01, rot2({1.0, 0.0}, 120))), "H", "H", 0.42});
            atoms.push_back({"H2", polar(c, kRingR + 1.08, 300), "H", "HA", 0.15});
        } else {  // G
            atoms.push_back({"O6", polar(c, kRingR + 1.23, 60), "O", "O", -0.58});
            atoms.push_back({"H1", add2(pos["N1"], {1.01, 0.0}), "H", "H", 0.40});
            const Vec2 n2 = polar(c, kRingR + 1.34, 300);
            atoms.push_back({"N2", n2, "N", "N", -0.85});
            atoms.push_back({"H21", add2(n2, {1.01, 0.0}), "H", "H", 0.42});
            atoms.push_back({"H22", add2(n2, scale2(1.01, rot2({1.0, 0.0}, -120))), "H", "H", 0.42});
        }
        c1 = polar(c5, r5 + 1.47, 222);
    }
    return {atoms, c1};
}

using Bond = std::pair<std::size_t, std::size_t>;

// nucleotide(base): the atoms (pair frame, strand A side) and bonds of one nucleotide
std::pair<std::vector<Atom>, std::vector<Bond>> nucleotide(char base) {
    const auto [baseAtoms, c1] = baseTemplate(base);
    std::vector<Atom> atoms;
    for (const BaseAtom& b : baseAtoms) {
        atoms.push_back({b.name, Vec{b.xy[0], b.xy[1], 0.0}, b.el, b.type, b.q});
    }
    const Vec c1p{c1[0], c1[1], 0.0};
    atoms.push_back({"C1'", c1p, "C", "CT", kC1Charge});
    for (const Offset& o : kSugar) {
        atoms.push_back({std::string(o.name), add(c1p, o.xyz), std::string(1, o.name[0]), std::string(o.type), o.charge});
    }
    const auto indexOf = [&atoms](std::string_view name) -> std::optional<std::size_t> {
        for (std::size_t k = 0; k < atoms.size(); ++k) {
            if (atoms[k].name == name) {
                return k;
            }
        }
        return std::nullopt;
    };
    for (const SugarH& h : kSugarH) {
        const Vec p = atoms[*indexOf(h.parent)].xyz;
        const bool hydroxyl = h.parent[0] == 'O';
        atoms.push_back({std::string(h.name), add(p, scale(hydroxyl ? 0.96 : 1.09, unit(h.hint))), "H",
                         hydroxyl ? "HO" : "H1", hydroxyl ? 0.42 : 0.08});
    }
    // neutralise to a net charge of -1 per nucleotide (Python's sum: in order, from 0)
    double total = 0.0;
    for (const Atom& a : atoms) {
        total += a.q;
    }
    const double corr = (-1.0 - total) / static_cast<double>(atoms.size());
    for (Atom& a : atoms) {
        a.q += corr;
    }
    std::vector<Bond> bonds;
    const auto bond = [&](std::string_view x, std::string_view y) {
        const auto ix = indexOf(x), iy = indexOf(y);
        if (ix && iy) {
            bonds.emplace_back(*ix, *iy);
        }
    };
    constexpr std::pair<std::string_view, std::string_view> kRingBonds[] = {
        {"N1", "C2"}, {"C2", "N3"}, {"N3", "C4"}, {"C4", "C5"}, {"C5", "C6"}, {"C6", "N1"}, {"C5", "N7"},
        {"N7", "C8"}, {"C8", "N9"}, {"N9", "C4"}, {"C2", "O2"}, {"C4", "O4"}, {"C4", "N4"}, {"N4", "H41"},
        {"N4", "H42"}, {"N3", "H3"}, {"C5", "H5"}, {"C6", "H6"}, {"C8", "H8"}, {"C6", "N6"}, {"N6", "H61"},
        {"N6", "H62"}, {"C2", "H2"}, {"C6", "O6"}, {"N1", "H1"}, {"C2", "N2"}, {"N2", "H21"}, {"N2", "H22"}};
    for (const auto& [x, y] : kRingBonds) {
        bond(x, y);
    }
    bond(base == 'A' || base == 'G' ? "N9" : "N1", "C1'");
    constexpr std::pair<std::string_view, std::string_view> kSugarBonds[] = {
        {"C1'", "O4'"}, {"C1'", "C2'"}, {"C2'", "O2'"}, {"C2'", "C3'"}, {"C3'", "C4'"}, {"C4'", "O4'"},
        {"C3'", "O3'"}, {"C4'", "C5'"}, {"C5'", "O5'"}, {"O5'", "P"},   {"P", "OP1"},   {"P", "OP2"},
        {"C1'", "H1'"}, {"C2'", "H2'"}, {"O2'", "HO2'"}, {"C3'", "H3'"}, {"C4'", "H4'"}, {"C5'", "H5'"},
        {"C5'", "H5''"}};
    for (const auto& [x, y] : kSugarBonds) {
        bond(x, y);
    }
    return {atoms, bonds};
}

char complement(char b) {
    switch (b) {
    case 'A':
        return 'U';
    case 'U':
        return 'A';
    case 'G':
        return 'C';
    default:
        return 'G';
    }
}

// ---- the system ---------------------------------------------------------------------------------

struct Residue {
    std::string name;
    std::size_t start;
};

struct System {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Residue> residues;
    std::int64_t soluteResidues = 0;  // n_solute_res
    std::optional<std::size_t> waterStart;
    std::vector<float> coords;  // frames x atoms x 3
    std::array<double, 6> box{};
};

// _bridge_site: a position for an Mg2+ bridging an aptamer and a target oxygen (Mg-O = 2.08 A)
Vec bridgeSite(const std::vector<Atom>& atoms, const std::vector<std::int64_t>& strand,
               const std::vector<std::int64_t>& level, std::int64_t nBp) {
    std::vector<std::size_t> aIdx, bIdx;
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        if (atoms[k].el == "O" && level[k] >= nBp - 4) {
            if (strand[k] == 0) {
                aIdx.push_back(k);
            } else if (strand[k] == 1) {
                bIdx.push_back(k);
            }
        }
    }
    const double twoPi = 2 * kPi;
    std::optional<std::pair<double, Vec>> best;
    for (std::size_t i : aIdx) {
        for (std::size_t j : bIdx) {
            const Vec& xi = atoms[i].xyz;
            const Vec& xj = atoms[j].xyz;
            const double d = norm(sub(xi, xj));
            if (!(2.6 <= d && d <= 4.0)) {
                continue;
            }
            const Vec axis = unit(sub(xj, xi));
            const Vec perp = unit(std::abs(axis[2]) < 0.9 ? cross(axis, {0.0, 0.0, 1.0}) : cross(axis, {1.0, 0.0, 0.0}));
            const Vec perp2 = cross(axis, perp);
            const Vec mid = scale(0.5, add(xi, xj));
            const double half = d / 2;
            const double lift = std::sqrt(std::max(0.0, 2.08 * 2.08 - half * half));
            const double step = twoPi / 24;  // np.linspace(0, 2 pi, 24, endpoint=False)
            for (int n = 0; n < 24; ++n) {
                const double ang = static_cast<double>(n) * step + 0.0;
                const double c = std::cos(ang), s = std::sin(ang);
                const Vec dir{c * perp[0] + s * perp2[0], c * perp[1] + s * perp2[1], c * perp[2] + s * perp2[2]};
                const Vec pos = add(mid, scale(lift, dir));
                // tree.query_ball_point(pos, 2.2) without i and j, then the smallest distance
                double clearance = 2.2;
                bool any = false;
                for (std::size_t k = 0; k < atoms.size(); ++k) {
                    if (k == i || k == j) {
                        continue;
                    }
                    const double dk = norm(sub(atoms[k].xyz, pos));
                    if (dk <= 2.2) {
                        clearance = any ? std::min(clearance, dk) : dk;
                        any = true;
                    }
                }
                if (!best || clearance > best->first) {
                    best = std::pair<double, Vec>(clearance, pos);
                }
            }
        }
    }
    if (!best) {  // next to the first aptamer oxygen
        return add(atoms[aIdx.at(0)].xyz, {2.08, 0.0, 0.0});
    }
    return best->second;
}

// np.arange(start, stop, step): start, start + delta, then start + i * delta with delta the
// difference of the first two values (numpy's fill)
std::vector<double> arange(double start, double stop, double step) {
    const double count = std::ceil((stop - start) / step);
    std::vector<double> out;
    if (!(count > 0)) {
        return out;
    }
    const auto n = static_cast<std::size_t>(count);
    out.reserve(n);
    out.push_back(start);
    if (n > 1) {
        out.push_back(start + step);
    }
    const double delta = n > 1 ? out[1] - out[0] : step;
    for (std::size_t i = 2; i < n; ++i) {
        out.push_back(start + static_cast<double>(i) * delta);
    }
    return out;
}

System buildSystem(const DemoOptions& options, const core::Reporter& reporter) {
    render::NumpyGenerator rng(options.seed);
    const std::string& sequence = options.sequence;
    const auto nBp = static_cast<std::int64_t>(sequence.size());
    std::string target(sequence.rbegin(), sequence.rend());
    for (char& c : target) {
        c = complement(c);
    }
    System sys;
    std::vector<std::int64_t> levelOf, strandOf;

    const auto addResidue = [&](char base, std::int64_t strand, std::int64_t level, bool dyad) {
        auto [nt, nb] = nucleotide(base);
        const Mat rot = levelRotation(level);
        const Vec shift{0.0, 0.0, static_cast<double>(level) * kRise};
        const std::size_t start = sys.atoms.size();
        sys.residues.push_back({std::string(1, base), start});
        for (Atom& a : nt) {
            Vec v = a.xyz;
            if (dyad) {  // DYAD = diag(-1, 1, -1)
                v = Vec{-1.0 * v[0] + 0.0 * v[1] + 0.0 * v[2], 0.0 * v[0] + 1.0 * v[1] + 0.0 * v[2],
                        0.0 * v[0] + 0.0 * v[1] + -1.0 * v[2]};
            }
            a.xyz = add(mul(rot, v), shift);
            sys.atoms.push_back(a);
            levelOf.push_back(level);
            strandOf.push_back(strand);
        }
        for (const auto& [x, y] : nb) {
            sys.bonds.emplace_back(start + x, start + y);
        }
        std::size_t p = 0, o3 = 0;
        for (std::size_t k = 0; k < nt.size(); ++k) {
            if (nt[k].name == "P") {
                p = start + k;
            }
            if (nt[k].name == "O3'") {
                o3 = start + k;
            }
        }
        return std::pair<std::size_t, std::size_t>(p, o3);
    };
    for (int strandIndex = 0; strandIndex < 2; ++strandIndex) {
        std::optional<std::size_t> prevO3;
        const std::string& bases = strandIndex == 0 ? sequence : target;
        for (std::int64_t m = 0; m < nBp; ++m) {
            const std::int64_t k = strandIndex == 0 ? m : nBp - 1 - m;
            const auto [p, o3] = addResidue(bases[static_cast<std::size_t>(m)], strandIndex, k, strandIndex == 1);
            if (prevO3) {
                sys.bonds.emplace_back(*prevO3, p);
            }
            prevO3 = o3;
        }
    }
    const std::size_t nSoluteNt = sys.atoms.size();

    const auto atomIndex = [&sys, nBp](std::int64_t strand, std::int64_t resInStrand, std::string_view name) {
        const auto ri = static_cast<std::size_t>(resInStrand + (strand == 1 ? nBp : 0));
        const std::size_t end = ri + 1 < sys.residues.size() ? sys.residues[ri + 1].start : sys.atoms.size();
        for (std::size_t idx = sys.residues[ri].start; idx < end; ++idx) {
            if (sys.atoms[idx].name == name) {
                return idx;
            }
        }
        throw std::out_of_range(std::string(name));
    };

    // Mg ions: bound to an aptamer phosphate, bridging aptamer and target, visiting a target phosphate
    const Vec bridge = bridgeSite(sys.atoms, strandOf, levelOf, nBp);
    struct MgSpec {
        int kind;  // 0 bound, 1 bridge, 2 transient
        std::int64_t strand, residue;
        std::string_view oxygen, phosphorus;
    };
    const MgSpec mgSpecs[] = {{0, 0, 2, "OP1", "P"}, {1, 0, 0, "", ""}, {2, 1, 3, "OP2", "P"}};
    std::vector<std::size_t> mgAtoms;
    for (const MgSpec& spec : mgSpecs) {
        const std::size_t idx = sys.atoms.size();
        sys.residues.push_back({"MG", idx});
        sys.atoms.push_back({"MG", spec.kind == 1 ? bridge : Vec{0.0, 0.0, 0.0}, "MG", "MG", 2.0});
        levelOf.push_back(-1);
        strandOf.push_back(2);
        mgAtoms.push_back(idx);
    }

    // the solute's centre and bounding box (numpy's mean along axis 0: summed row by row)
    Vec center{0.0, 0.0, 0.0}, lo{}, hi{};
    for (std::size_t k = 0; k < nSoluteNt; ++k) {
        const Vec& x = sys.atoms[k].xyz;
        for (int d = 0; d < 3; ++d) {
            center[static_cast<std::size_t>(d)] += x[static_cast<std::size_t>(d)];
            lo[static_cast<std::size_t>(d)] = k == 0 ? x[static_cast<std::size_t>(d)] : std::min(lo[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]);
            hi[static_cast<std::size_t>(d)] = k == 0 ? x[static_cast<std::size_t>(d)] : std::max(hi[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]);
        }
    }
    for (double& c : center) {
        c /= static_cast<double>(nSoluteNt);
    }
    for (int d = 0; d < 3; ++d) {
        lo[static_cast<std::size_t>(d)] -= 12.0;
        hi[static_cast<std::size_t>(d)] += 12.0;
    }
    const Vec boxLen = sub(hi, lo);
    const Vec shiftAll{-lo[0], -lo[1], -lo[2]};
    for (Atom& a : sys.atoms) {
        a.xyz = add(a.xyz, shiftAll);
    }

    const std::int64_t nNa = 2 * nBp - 2 * 3;
    std::vector<std::size_t> naAtoms;
    if (options.solvent) {
        const std::vector<double> grid = arange(1.6, std::max({boxLen[0], boxLen[1], boxLen[2]}), 3.1);
        std::vector<double> pts;
        for (double x : grid) {
            if (!(x < boxLen[0] - 1.0)) {
                continue;
            }
            for (double y : grid) {
                if (!(y < boxLen[1] - 1.0)) {
                    continue;
                }
                for (double z : grid) {
                    if (z < boxLen[2] - 1.0) {
                        pts.insert(pts.end(), {x, y, z});
                    }
                }
            }
        }
        // cKDTree(solute).query(pts): the distance to the nearest nucleotide atom, kept above 3 A
        std::vector<double> kept;
        for (std::size_t p = 0; p + 2 < pts.size(); p += 3) {
            const Vec point{pts[p], pts[p + 1], pts[p + 2]};
            double nearest = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k < nSoluteNt; ++k) {
                const Vec diff = sub(sys.atoms[k].xyz, point);
                nearest = std::min(nearest, diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
            }
            if (std::sqrt(nearest) > 3.0) {
                kept.insert(kept.end(), point.begin(), point.end());
            }
        }
        pts = std::move(kept);
        rng.shuffleRows(pts, 3);
        const std::size_t nPts = pts.size() / 3;
        const std::size_t naCount = nNa > 0 ? std::min(nPts, static_cast<std::size_t>(nNa)) : 0;
        for (std::size_t n = 0; n < naCount; ++n) {
            const std::size_t idx = sys.atoms.size();
            sys.residues.push_back({"Na+", idx});
            sys.atoms.push_back({"Na+", Vec{pts[3 * n], pts[3 * n + 1], pts[3 * n + 2]}, "NA", "Na+", 1.0});
            naAtoms.push_back(idx);
            levelOf.push_back(-1);
            strandOf.push_back(3);
        }
        sys.waterStart = sys.atoms.size();
        const double ang = radians(104.52) / 2;
        const double ca = std::cos(ang), sa = std::sin(ang);
        for (std::size_t n = naCount; n < nPts; ++n) {
            const Vec p{pts[3 * n], pts[3 * n + 1], pts[3 * n + 2]};
            Vec r1{};
            for (double& v : r1) {
                v = rng.normal(0.0, 1.0);
            }
            const Vec u = unit(r1);
            Vec r2{};
            for (double& v : r2) {
                v = rng.normal(0.0, 1.0);
            }
            const Vec v = unit(cross(u, r2));
            const Vec h1 = add(p, scale(0.9572, Vec{ca * u[0] + sa * v[0], ca * u[1] + sa * v[1], ca * u[2] + sa * v[2]}));
            const Vec h2 = add(p, scale(0.9572, Vec{ca * u[0] - sa * v[0], ca * u[1] - sa * v[1], ca * u[2] - sa * v[2]}));
            const std::size_t idx = sys.atoms.size();
            sys.residues.push_back({"WAT", idx});
            sys.atoms.push_back({"O", p, "O", "OW", -0.834});
            sys.atoms.push_back({"H1", h1, "H", "HW", 0.417});
            sys.atoms.push_back({"H2", h2, "H", "HW", 0.417});
            sys.bonds.emplace_back(idx, idx + 1);
            sys.bonds.emplace_back(idx, idx + 2);
            for (int k = 0; k < 3; ++k) {
                levelOf.push_back(-1);
                strandOf.push_back(4);
            }
        }
    }

    const std::size_t nAtoms = sys.atoms.size();
    std::vector<char> rigid(nAtoms, 0);
    for (std::size_t k = 0; k < nSoluteNt; ++k) {
        rigid[k] = 1;
    }
    for (std::size_t idx : mgAtoms) {
        rigid[idx] = 1;
    }
    const Vec centre = add(center, shiftAll);
    // fray_w = clip((level - (n_bp - 5)) / 4, 0, 1) * (strand == 1): one weight per level
    std::vector<double> frayWeight(static_cast<std::size_t>(nBp), 0.0);
    for (std::int64_t k = 0; k < nBp; ++k) {
        frayWeight[static_cast<std::size_t>(k)] =
            std::clamp(static_cast<double>(k - (nBp - 5)) / 4.0, 0.0, 1.0) * 1.0;
    }
    std::vector<Mat> levelRot;
    for (std::int64_t k = 0; k < nBp; ++k) {
        levelRot.push_back(levelRotation(k));
    }

    const auto nFrames = static_cast<std::size_t>(options.frames);
    sys.coords.assign(nFrames * nAtoms * 3, 0.0f);
    std::vector<Vec> x(nAtoms);
    for (std::size_t t = 0; t < nFrames; ++t) {
        reporter.checkCancel();
        for (std::size_t k = 0; k < nAtoms; ++k) {
            x[k] = sys.atoms[k].xyz;
        }
        const double s = static_cast<double>(t) / static_cast<double>(std::max<std::int64_t>(1, options.frames - 1));
        const double fray = smoothStep((s - 0.30) / 0.15) * (1.0 - smoothStep((s - 0.65) / 0.15));
        // fray: target nucleotides near the top move away from the aptamer
        for (std::int64_t k = 0; k < nBp; ++k) {
            const double w = frayWeight[static_cast<std::size_t>(k)];
            if (w == 0) {
                continue;
            }
            const Vec away = mul(levelRot[static_cast<std::size_t>(k)], Vec{1.0, 0.0, 0.0});
            const Vec move = scale(6.5 * fray * w, away);
            for (std::size_t a = 0; a < nSoluteNt; ++a) {
                if (strandOf[a] == 1 && levelOf[a] == k) {
                    x[a] = add(x[a], move);
                }
            }
        }
        for (std::size_t a = 0; a < nSoluteNt; ++a) {
            for (double& v : x[a]) {
                v += rng.normal(0.0, 0.12);
            }
        }
        for (std::size_t m = 0; m < 3; ++m) {
            const MgSpec& spec = mgSpecs[m];
            const std::size_t idx = mgAtoms[m];
            if (spec.kind != 1) {  // the bridging ion stays with the aptamer
                const Vec o = x[atomIndex(spec.strand, spec.residue, spec.oxygen)];
                const Vec p = x[atomIndex(spec.strand, spec.residue, spec.phosphorus)];
                double dist = 2.08;
                if (spec.kind == 2) {  // inner -> outer sphere -> unbound -> back
                    const double outer = smoothStep((s - 0.35) / 0.1) * (1 - smoothStep((s - 0.85) / 0.1));
                    const double away = smoothStep((s - 0.55) / 0.1) * (1 - smoothStep((s - 0.8) / 0.1));
                    dist += 2.0 * outer + 5.0 * away;
                }
                x[idx] = add(o, scale(dist, unit(sub(o, p))));
            }
            for (double& v : x[idx]) {
                v += rng.normal(0.0, 0.05);
            }
        }
        // rigid-body drift of the solute and its ions
        const Mat rot = rotation({1.0, 1.0, 0.2}, radians(25.0) * std::sin(2 * kPi * s));
        const Vec drift{2.0 * s, -1.0 * s, 1.5 * s};
        for (std::size_t a = 0; a < nAtoms; ++a) {
            if (rigid[a]) {
                const Vec d = sub(x[a], centre);
                x[a] = add(add(mul(rot, d), centre), drift);
            }
        }
        if (options.solvent) {
            for (std::size_t a = 0; a < nAtoms; ++a) {
                if (!rigid[a]) {
                    for (double& v : x[a]) {
                        v += rng.normal(0.0, 0.25);
                    }
                }
            }
        }
        float* out = sys.coords.data() + t * nAtoms * 3;
        for (std::size_t a = 0; a < nAtoms; ++a) {
            for (std::size_t d = 0; d < 3; ++d) {
                out[3 * a + d] = static_cast<float>(x[a][d]);
            }
        }
        reporter.progress(static_cast<double>(t + 1) / static_cast<double>(nFrames), "demo trajectory");
    }
    sys.box = {boxLen[0], boxLen[1], boxLen[2], 90.0, 90.0, 90.0};
    sys.soluteResidues = 2 * nBp + static_cast<std::int64_t>(mgAtoms.size() + naAtoms.size());
    return sys;
}

// write_topology
void writeTopology(const fs::path& path, const System& sys) {
    using io::Parm7Section;
    const std::vector<Atom>& atoms = sys.atoms;
    const std::size_t nAtoms = atoms.size();
    std::vector<std::size_t> typeNames;  // indices into kLjTypes, in its order, those used
    for (std::size_t t = 0; t < std::size(kLjTypes); ++t) {
        if (std::any_of(atoms.begin(), atoms.end(), [&](const Atom& a) { return a.type == kLjTypes[t].name; })) {
            typeNames.push_back(t);
        }
    }
    const auto typeIndex = [&typeNames](const std::string& type) {
        for (std::size_t k = 0; k < typeNames.size(); ++k) {
            if (kLjTypes[typeNames[k]].name == type) {
                return static_cast<std::int64_t>(k + 1);
            }
        }
        throw std::out_of_range(type);
    };
    const std::size_t ntypes = typeNames.size();
    std::vector<double> acoef, bcoef;
    for (std::size_t i = 0; i < ntypes; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const LjType& ti = kLjTypes[typeNames[i]];
            const LjType& tj = kLjTypes[typeNames[j]];
            const double rmin = ti.rstar + tj.rstar;
            const double eps = std::sqrt(ti.epsilon * tj.epsilon);
            acoef.push_back(eps * std::pow(rmin, 12.0));
            bcoef.push_back(2.0 * eps * std::pow(rmin, 6.0));
        }
    }
    std::vector<std::int64_t> nbIndex;
    for (std::int64_t i = 1; i <= static_cast<std::int64_t>(ntypes); ++i) {
        for (std::int64_t j = 1; j <= static_cast<std::int64_t>(ntypes); ++j) {
            const std::int64_t h = std::max(i, j), l = std::min(i, j);
            nbIndex.push_back(h * (h - 1) / 2 + l);
        }
    }
    std::vector<std::int64_t> bondsH, bondsA;
    for (const auto& [a, b] : sys.bonds) {
        const bool h = atoms[a].el == "H" || atoms[b].el == "H";
        std::vector<std::int64_t>& list = h ? bondsH : bondsA;
        list.insert(list.end(), {static_cast<std::int64_t>(3 * a), static_cast<std::int64_t>(3 * b), h ? 1 : 2});
    }
    std::vector<std::int64_t> resPtr;
    std::int64_t largest = 0;
    for (std::size_t r = 0; r < sys.residues.size(); ++r) {
        resPtr.push_back(static_cast<std::int64_t>(sys.residues[r].start) + 1);
        const std::size_t end = r + 1 < sys.residues.size() ? sys.residues[r + 1].start : nAtoms;
        largest = std::max(largest, static_cast<std::int64_t>(end - sys.residues[r].start));
    }
    // molecules: runs of atoms with the same union-find root (MDATK's grouping)
    std::vector<std::size_t> parent(nAtoms);
    for (std::size_t a = 0; a < nAtoms; ++a) {
        parent[a] = a;
    }
    const auto find = [&parent](std::size_t a) {
        while (parent[a] != a) {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        return a;
    };
    for (const auto& [a, b] : sys.bonds) {
        const std::size_t ra = find(a), rb = find(b);
        if (ra != rb) {
            parent[std::max(ra, rb)] = std::min(ra, rb);
        }
    }
    std::vector<std::int64_t> molSizes;
    std::optional<std::size_t> last;
    for (std::size_t a = 0; a < nAtoms; ++a) {
        const std::size_t r = find(a);
        if (!last || r != *last) {
            molSizes.push_back(0);
            last = r;
        }
        ++molSizes.back();
    }
    const auto n = static_cast<std::int64_t>(nAtoms);
    const std::int64_t ifbox = 1;
    const std::vector<std::int64_t> pointers{n, static_cast<std::int64_t>(ntypes), static_cast<std::int64_t>(bondsH.size() / 3),
                                             static_cast<std::int64_t>(bondsA.size() / 3), 0, 0, 0, 0, 0, 0, n,
                                             static_cast<std::int64_t>(sys.residues.size()),
                                             static_cast<std::int64_t>(bondsA.size() / 3), 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0,
                                             0, 0, 0, ifbox, largest, 0, 0};
    std::vector<std::string> names, types, resLabels;
    std::vector<double> charges, masses, radii;
    std::vector<std::int64_t> numbers, typeIdx;
    for (const Atom& a : atoms) {
        names.push_back(a.name);
        types.push_back(a.type);
        charges.push_back(a.q * kAmberChargeFactor);
        masses.push_back(elementMass(a.el));
        radii.push_back(a.el == "H" ? 1.2 : 1.5);
        numbers.push_back(elementZ(a.el));
        typeIdx.push_back(typeIndex(a.type));
    }
    for (const Residue& r : sys.residues) {
        resLabels.push_back(r.name);
    }
    const auto ints = [](std::vector<std::int64_t> v) { return io::Parm7Values(std::move(v)); };
    const auto reals = [](std::vector<double> v) { return io::Parm7Values(std::move(v)); };
    const auto texts = [](std::vector<std::string> v) { return io::Parm7Values(std::move(v)); };
    std::vector<Parm7Section> sections{
        {"POINTERS", 'i', ints(pointers)},
        {"ATOM_NAME", 'a', texts(names)},
        {"CHARGE", 'e', reals(charges)},
        {"ATOMIC_NUMBER", 'i', ints(numbers)},
        {"MASS", 'e', reals(masses)},
        {"ATOM_TYPE_INDEX", 'i', ints(typeIdx)},
        {"NUMBER_EXCLUDED_ATOMS", 'i', ints(std::vector<std::int64_t>(nAtoms, 1))},
        {"NONBONDED_PARM_INDEX", 'i', ints(nbIndex)},
        {"RESIDUE_LABEL", 'a', texts(resLabels)},
        {"RESIDUE_POINTER", 'i', ints(resPtr)},
        {"BOND_FORCE_CONSTANT", 'e', reals({340.0, 300.0})},
        {"BOND_EQUIL_VALUE", 'e', reals({1.09, 1.50})},
        {"ANGLE_FORCE_CONSTANT", 'e', reals({})},
        {"ANGLE_EQUIL_VALUE", 'e', reals({})},
        {"DIHEDRAL_FORCE_CONSTANT", 'e', reals({})},
        {"DIHEDRAL_PERIODICITY", 'e', reals({})},
        {"DIHEDRAL_PHASE", 'e', reals({})},
        {"SCEE_SCALE_FACTOR", 'e', reals({})},
        {"SCNB_SCALE_FACTOR", 'e', reals({})},
        {"SOLTY", 'e', reals({0.0})},
        {"LENNARD_JONES_ACOEF", 'e', reals(acoef)},
        {"LENNARD_JONES_BCOEF", 'e', reals(bcoef)},
        {"BONDS_INC_HYDROGEN", 'i', ints(bondsH)},
        {"BONDS_WITHOUT_HYDROGEN", 'i', ints(bondsA)},
        {"ANGLES_INC_HYDROGEN", 'i', ints({})},
        {"ANGLES_WITHOUT_HYDROGEN", 'i', ints({})},
        {"DIHEDRALS_INC_HYDROGEN", 'i', ints({})},
        {"DIHEDRALS_WITHOUT_HYDROGEN", 'i', ints({})},
        {"EXCLUDED_ATOMS_LIST", 'i', ints(std::vector<std::int64_t>(nAtoms, 0))},
        {"HBOND_ACOEF", 'e', reals({})},
        {"HBOND_BCOEF", 'e', reals({})},
        {"HBCUT", 'e', reals({})},
        {"AMBER_ATOM_TYPE", 'a', texts(types)},
        {"TREE_CHAIN_CLASSIFICATION", 'a', texts(std::vector<std::string>(nAtoms, "BLA"))},
        {"JOIN_ARRAY", 'i', ints(std::vector<std::int64_t>(nAtoms, 0))},
        {"IROTAT", 'i', ints(std::vector<std::int64_t>(nAtoms, 0))},
    };
    auto firstSolventMol = static_cast<std::int64_t>(molSizes.size()) + 1;
    if (sys.waterStart) {
        std::int64_t start = 0, before = 0;
        for (std::int64_t size : molSizes) {
            if (start < static_cast<std::int64_t>(*sys.waterStart)) {
                ++before;
            }
            start += size;
        }
        firstSolventMol = 1 + before;
    }
    sections.push_back({"SOLVENT_POINTERS", 'i',
                        ints({sys.soluteResidues, static_cast<std::int64_t>(molSizes.size()), firstSolventMol})});
    sections.push_back({"ATOMS_PER_MOLECULE", 'i', ints(molSizes)});
    sections.push_back({"BOX_DIMENSIONS", 'e', reals({90.0, sys.box[0], sys.box[1], sys.box[2]})});
    sections.push_back({"RADIUS_SET", 's', texts({"modified Bondi radii (mbondi2)"})});
    sections.push_back({"RADII", 'e', reals(radii)});
    sections.push_back({"SCREEN", 'e', reals(std::vector<double>(nAtoms, 0.8))});
    sections.push_back({"IPOL", 'i', ints({0})});
    io::writeParm7(path, sections, "MDATK synthetic RNA aptamer/target demo");
}

// write_trajectory: MDAnalysis' NCDFWriter over a MemoryReader (time = frame x dt, the box in float32)
void writeTrajectory(const fs::path& path, const System& sys, double dtPs) {
    const std::size_t nAtoms = sys.atoms.size();
    const std::size_t nFrames = nAtoms == 0 ? 0 : sys.coords.size() / (nAtoms * 3);
    std::vector<core::Frame> frames(nFrames);
    std::array<float, 6> box{};
    for (std::size_t k = 0; k < 6; ++k) {
        box[k] = static_cast<float>(sys.box[k]);
    }
    for (std::size_t t = 0; t < nFrames; ++t) {
        core::Frame& f = frames[t];
        f.index = static_cast<std::int64_t>(t);
        f.time = static_cast<double>(t) * dtPs;
        f.positions.assign(sys.coords.begin() + static_cast<std::ptrdiff_t>(t * nAtoms * 3),
                           sys.coords.begin() + static_cast<std::ptrdiff_t>((t + 1) * nAtoms * 3));
        f.box = box;
    }
    io::writeNetcdfTrajectory(path, frames, nAtoms);
}

// os.path.join(out_dir, name)
std::string joinText(const std::string& dir, const std::string& name) {
#ifdef _WIN32
    if (dir.empty()) {
        return name;
    }
    const char lastChar = dir.back();
    return dir + (lastChar == '/' || lastChar == '\\' || lastChar == ':' ? "" : "\\") + name;
#else
    return core::posixpath::join(dir, name);
#endif
}

// f"{x:g}"
std::string formatG(double x) {
    char text[64];
    std::snprintf(text, sizeof text, "%g", x);
    return text;
}

}  // namespace

core::Project demoProject(const std::string& topology, const std::string& trajectory, const std::string& outDir,
                          std::int64_t basePairs, bool cleaning, double midTimePs) {
    core::Project p;
    p.inputs.topology = topology;
    p.inputs.trajectories = {trajectory};
    p.inputs.time_unit = "ns";
    p.groups[0].residues = "1-" + std::to_string(basePairs);
    p.groups[1].residues = std::to_string(basePairs + 1) + "-" + std::to_string(2 * basePairs);
    p.cleaning.enabled = cleaning;
    // screened electrostatics (epsilon = 4r) suit a stripped, solvent-free complex
    p.analyses.energy.dielectric = 4.0;
    p.analyses.energy.distance_dependent = true;
    p.analyses.energy.target_value = -60.0;
    p.visual.rotation = {-90.0, 0.0, 0.0};  // helix axis vertical, viewed from the side
    p.output.out_dir = outDir;
    p.output.pdb.enabled = true;
    p.output.pdb.frames = "0, last";
    p.output.pdb.times = formatG(midTimePs) + "ps";
    p.output.pdb.selection = "Aptamer, Target";
    p.layout.preset = "large_plus_4";
    p.layout.panels.push_back({"heatmap", {"mg/occupancy"}, "Mg occupancy", {}});
    return p;
}

DemoFiles buildDemo(const std::string& outDirText, const DemoOptions& options, const core::Reporter& reporter) {
    if (options.frames < 1) {
        throw core::ValueError("the demo needs at least one frame (--n-frames " + std::to_string(options.frames) + ")");
    }
    if (!(std::isfinite(options.dtPs) && options.dtPs > 0)) {
        throw core::ValueError("the time between frames must be a positive number of ps (--dt " +
                               core::pythonRepr(core::Value(options.dtPs)) + ")");
    }
    if (options.sequence.size() < 4 ||
        options.sequence.find_first_not_of("ACGU") != std::string::npos) {
        throw core::ValueError("demo sequence " + core::pythonRepr(core::Value(options.sequence)) +
                               ": at least four of the letters A, C, G and U");
    }
    const fs::path outDir = core::pathFromUtf8(outDirText);
    fs::create_directories(outDir);
    const System sys = buildSystem(options, reporter);
    DemoFiles files;
    files.topology = joinText(outDirText, "demo.parm7");
    files.trajectory = joinText(outDirText, "demo.nc");
    files.project = joinText(outDirText, "demo.mdatk.toml");
    writeTopology(outDir / "demo.parm7", sys);
    writeTrajectory(outDir / "demo.nc", sys, options.dtPs);
    const auto nBp = static_cast<std::int64_t>(options.sequence.size());
    core::Project project = demoProject("demo.parm7", "demo.nc", "MDATK_output", nBp,
                                        options.solvent && core::process::which("cpptraj").has_value(),
                                        static_cast<double>(options.frames / 2) * options.dtPs);
    project.notes =
        "Synthetic 12-bp RNA aptamer/target duplex written by 'mdatk demo' to try every feature; its geometry is "
        "schematic.";
    core::saveProject(outDir / "demo.mdatk.toml", project);
    return files;
}

}  // namespace mdbridge::app
