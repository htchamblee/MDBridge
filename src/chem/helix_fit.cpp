// MDATK's structure builder, part 2 (mdatk/builder.py v0.0.4.1, MDATK's own code, ported): the other
// pairings (reverse-wc, hoogsteen, reverse-hoogsteen, hh) and the triplexes. Every nucleotide is the
// 4JRD nucleotide (in the fitted 4JRD screw's frame, CHD-02) carrying the Arnott base; partner bases
// are placed on the H-bond edges (_pair_fit), the helix is the screw that repeats the reference level
// with closed O3'-P links and no close contacts (_helix_fit), mismatches are moved off the bases they
// would hit (_level_pose), and a helix whose levels differ is closed jointly (_close_levels) - by
// continuation (docs/DIFFERENCES.md CH-16, owner decision D17) instead of MDATK's four truncated
// passes. Then the notes: links, H-bonds, mismatches, close contacts (hydrogens included), backbone
// torsions (_geometry_notes, _close_contacts, _torsion_note), and the terminal hydroxyls turned away
// from the other atoms (_turn_terminal_hydroxyls).
//
// The fits are numpy / SciPy code in MDATK: least_squares (TRF, numerics::leastSquares), cKDTree's
// query_pairs (numerics::queryPairs), np.linalg.svd of 3 x 3 and k x 3 matrices (core::svd3x3,
// core::svdRows3). They are converged here (tolerances 1e-14, the closure's steps 1e-10), so their
// results agree with the patched MDATK's within docs/specs/m6-plan.md D10, not bit for bit.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "helix_internal.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/linalg.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/numerics/least_squares.hpp"
#include "mdbridge/numerics/pairs.hpp"
#include "mdbridge/chem/reference_tables.hpp"

namespace mdbridge::chem::helix_detail {

namespace {

// --- constants (builder.py) ------------------------------------------------------------------------

constexpr double kPi = 3.141592653589793;
constexpr double kHbond = 2.9;    // Å, donor-acceptor distance of the built pairs
constexpr double kLink = 1.607;   // Å, O3'-P
constexpr double kLinkLow = 1.5;  // LINK_RANGE
constexpr double kLinkHigh = 1.7;
constexpr double kClashNb = 2.4;    // Å, heavy atoms three bonds apart
constexpr double kJrdChi = -170.0;  // chi of the nucleosides of _core (dchi 0)
constexpr double kChiLow = -185.0;  // CHI_RANGE
constexpr double kChiHigh = -145.0;
// _GEOM
constexpr double kC5O5 = 1.443;
constexpr double kC4C5O5 = 109.9;
constexpr double kO5P = 1.592;
constexpr double kC5O5P = 119.6;
constexpr double kC3O3P = 119.7;
constexpr double kO3PO5 = 104.0;
constexpr double kPOp = 1.485;
constexpr double kOpPOp = 119.6;
// _TORSIONS
constexpr double kAlpha = -62.0;
constexpr double kBeta = 174.0;
constexpr double kGamma = 58.0;
constexpr double kEpsilon = -146.0;
constexpr double kZeta = -66.0;
constexpr double kBondCH2 = 1.09;
constexpr double kLevelRotation = 10.0;  // LEVEL_MOVE: degrees, Å
constexpr double kLevelShift = 0.8;
constexpr double kFitContactWeight = 0.03;

// CONTACT_FIT / CONTACT_REPORT
struct ContactTable {
    double heavy, hbond, hh, polar, hx;
};
constexpr ContactTable kContactFit = {3.1, 2.7, 2.1, 1.8, 2.4};
constexpr ContactTable kContactReport = {3.0, 2.6, 2.0, 1.6, 1.9};

// The closure (CH-16): continuation in kSteps steps, a step halved while its solve stops on its
// evaluation limit or its contact pairs need many refreshes (tests/parity/improvements_chem.py
// _patch_closure: the same constants)
constexpr int kSteps = 10;
constexpr int kMaxHalvings = 5;
constexpr double kClosureWeight = 0.01;
constexpr double kClosureTol = 1e-10;
constexpr std::size_t kClosureMaxNfev = 1000;
constexpr std::size_t kClosureFinalNfev = 20000;
constexpr double kClosureCutoff = 4.5;
constexpr int kRoundsOk = 4;
constexpr int kMaxRounds = 50;
constexpr double kSmallFitTol = 1e-14;  // _pair_fit, _level_pose, _helix_fit
constexpr std::size_t kSmallFitMaxNfev = 5000;

constexpr std::array<std::string_view, 7> kSugar = {"C2'", "O2'", "C3'", "O3'", "C4'", "O4'", "C5'"};
constexpr std::array<std::string_view, 5> kSugarH = {"H1'", "H2'", "H3'", "H4'", "HO2'"};
// _RIGID = _SUGAR + _SUGAR_H
constexpr std::array<std::string_view, 12> kRigid = {"C2'", "O2'", "C3'", "O3'", "C4'", "O4'",
                                                     "C5'", "H1'", "H2'", "H3'", "H4'", "HO2'"};
constexpr std::array<std::string_view, 19> kStrandNames = {"P",   "OP1", "OP2",  "O5'", "C5'", "C4'", "O4'",
                                                           "C3'", "O3'", "C2'",  "O2'", "C1'", "H1'", "H2'",
                                                           "H3'", "H4'", "HO2'", "H5'", "H5''"};
constexpr std::size_t kBaseSlots = 15;
constexpr std::size_t kSlots = kStrandNames.size() + kBaseSlots;  // 34
// slot of each _RIGID atom in kStrandNames
constexpr std::array<std::size_t, 12> kRigidSlot = {9, 10, 7, 8, 5, 6, 4, 12, 13, 14, 15, 16};
constexpr std::size_t kSlotP = 0, kSlotOP1 = 1, kSlotOP2 = 2, kSlotO5 = 3, kSlotC5 = 4, kSlotC4 = 5, kSlotC3 = 7,
                      kSlotO3 = 8, kSlotC1 = 11, kSlotH5 = 17, kSlotH5b = 18;

// --- small vector helpers ------------------------------------------------------------------------

double dot(const Vec3& a, const Vec3& b) {
    double s = a[0] * b[0];
    s += a[1] * b[1];
    s += a[2] * b[2];
    return s;
}

double degrees(double radiansValue) {
    return radiansValue * (180.0 / kPi);
}

// numpy's (x + 180) % 360 - 180 (Python's modulo: the sign of the divisor)
double wrap(double x) {
    double r = std::fmod(x + 180.0, 360.0);
    if (r != 0.0 && r < 0.0) {
        r += 360.0;
    }
    return r - 180.0;
}

double clip0(double x) {
    return x > 0.0 ? x : 0.0;
}

Mat3 matMul(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            double s = a[i][0] * b[0][j];
            s += a[i][1] * b[1][j];
            s += a[i][2] * b[2][j];
            out[i][j] = s;
        }
    }
    return out;
}

Mat3 identity() {
    return {Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
}

// _dihedral (degrees)
double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const Vec3 b0 = a - b;
    Vec3 b1 = c - b;
    const Vec3 b2 = d - c;
    b1 = unit(b1);
    const Vec3 v = b0 - dot(b0, b1) * b1;
    const Vec3 w = b2 - dot(b2, b1) * b1;
    return degrees(std::atan2(dot(cross(b1, v), w), dot(v, w)));
}

// _angle (degrees)
double angle(const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 u = a - b;
    const Vec3 v = c - b;
    const double cosine = dot(u, v) / (norm(u) * norm(v));
    return degrees(std::acos(std::clamp(cosine, -1.0, 1.0)));
}

// _rotvecs: Rodrigues of a rotation vector in radians
Mat3 rodrigues(const Vec3& v) {
    const double th = norm(v);
    const double d = th > 1e-12 ? th : 1.0;
    const Vec3 k = {v[0] / d, v[1] / d, v[2] / d};
    const Mat3 kx = {Vec3{0.0, -k[2], k[1]}, Vec3{k[2], 0.0, -k[0]}, Vec3{-k[1], k[0], 0.0}};
    const Mat3 kk = matMul(kx, kx);
    const double s = std::sin(th);
    const double c = 1.0 - std::cos(th);
    Mat3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out[i][j] = ((i == j ? 1.0 : 0.0) + s * kx[i][j]) + c * kk[i][j];
        }
    }
    return out;
}

// _rotvec (radians)
Mat3 rotvec(const Vec3& v) {
    const double th = norm(v);
    if (th < 1e-12) {
        return identity();
    }
    return turnMatrix({v[0] / th, v[1] / th, v[2] / th}, degrees(th));
}

// the rotation vector (radians) of a rotation matrix (angle below pi)
Vec3 rotationVector(const Mat3& r) {
    const double c = ((r[0][0] + r[1][1] + r[2][2]) - 1.0) / 2.0;
    const double th = std::acos(std::min(1.0, std::max(-1.0, c)));
    if (th < 1e-12) {
        return {0.0, 0.0, 0.0};
    }
    const double s = 2.0 * std::sin(th);
    return {(r[2][1] - r[1][2]) / s * th, (r[0][2] - r[2][0]) / s * th, (r[1][0] - r[0][1]) / s * th};
}

// x @ rot.T + shift
Vec3 transform(const Mat3& rot, const Vec3& shift, const Vec3& x) {
    return matVec(rot, x) + shift;
}

Vec3 mean(std::span<const Vec3> points) {
    Vec3 s{0.0, 0.0, 0.0};
    for (const Vec3& p : points) {
        s = s + p;
    }
    const double n = static_cast<double>(points.size());
    return {s[0] / n, s[1] / n, s[2] / n};
}

// mdatk.spatial.kabsch: (R, mobile centre, reference centre), (x - mc) @ R.T + rc maps mobile on
// reference
struct Kabsch {
    Mat3 rotation = identity();
    Vec3 mobileCentre{};
    Vec3 referenceCentre{};
};

Kabsch kabsch(std::span<const Vec3> mobile, std::span<const Vec3> reference, std::span<const double> weights) {
    const std::size_t n = mobile.size();
    double total = 0.0;
    for (const double w : weights) {
        total += w;
    }
    std::vector<double> w(n);
    for (std::size_t k = 0; k < n; ++k) {
        w[k] = weights[k] / total;
    }
    Kabsch out;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t c = 0; c < 3; ++c) {
            out.mobileCentre[c] += w[k] * mobile[k][c];
            out.referenceCentre[c] += w[k] * reference[k][c];
        }
    }
    if (n < 3) {
        return out;
    }
    core::Mat3 h{};
    for (std::size_t k = 0; k < n; ++k) {
        const Vec3 x = (w[k] * (mobile[k] - out.mobileCentre));
        const Vec3 y = reference[k] - out.referenceCentre;
        for (std::size_t p = 0; p < 3; ++p) {
            for (std::size_t q = 0; q < 3; ++q) {
                h[p][q] += x[p] * y[q];
            }
        }
    }
    const core::Svd3 svd = core::svd3x3(h);
    core::Mat3 vu{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            vu[i][j] = svd.vt[0][i] * svd.u[j][0] + svd.vt[1][i] * svd.u[j][1] + svd.vt[2][i] * svd.u[j][2];
        }
    }
    const double det = core::det3(vu);
    const double d = det < 0 ? -1.0 : 1.0;
    for (std::size_t i = 0; i < 3; ++i) {
        const Vec3 v{svd.vt[0][i], svd.vt[1][i], svd.vt[2][i] * d};
        for (std::size_t j = 0; j < 3; ++j) {
            out.rotation[i][j] = v[0] * svd.u[j][0] + v[1] * svd.u[j][1] + v[2] * svd.u[j][2];
        }
    }
    return out;
}

Vec3 kabschMove(const Kabsch& k, const Vec3& x) {
    return matVec(k.rotation, x - k.mobileCentre) + k.referenceCentre;
}

// _place_many / _place (the same arithmetic)
Vec3 placeAtom(const Vec3& a, const Vec3& b, const Vec3& c, double bond, double angleDeg, double torsion) {
    const Vec3 bc = unit(c - b);
    const Vec3 n = unit(cross(b - a, bc));
    const Vec3 m = cross(n, bc);
    const double th = radians(angleDeg);
    const double ph = radians(torsion);
    return c + bond * ((-std::cos(th) * bc + (std::sin(th) * std::cos(ph)) * m) + (std::sin(th) * std::sin(ph)) * n);
}

// _phosphate_oxygens: OP1, OP2 on the bisector plane of O3'-P-O5'
std::pair<Vec3, Vec3> phosphateOxygens(const Vec3& p, const Vec3& o3, const Vec3& o5) {
    const Vec3 u1 = unit(o3 - p);
    const Vec3 u2 = unit(o5 - p);
    const Vec3 bis = -unit(u1 + u2);
    const Vec3 perp = unit(cross(u1, u2));
    const double half = radians(kOpPOp / 2);
    const Vec3 a = p + kPOp * (std::cos(half) * bis + std::sin(half) * perp);
    const Vec3 b = p + kPOp * (std::cos(half) * bis - std::sin(half) * perp);
    if (helixChiralVolume(p, o5, a, b) > 0) {
        return {b, a};
    }
    return {a, b};
}

// _methylene_h of C5' (bonded to O5' and C4')
std::pair<Vec3, Vec3> methyleneH(const Vec3& c, const Vec3& a, const Vec3& b) {
    const Vec3 u1 = unit(a - c);
    const Vec3 u2 = unit(b - c);
    const Vec3 bis = -unit(u1 + u2);
    const Vec3 perp = unit(cross(u1, u2));
    const double half = radians(109.5 / 2);
    const Vec3 h1 = c + kBondCH2 * (std::cos(half) * bis + std::sin(half) * perp);
    const Vec3 h2 = c + kBondCH2 * (std::cos(half) * bis - std::sin(half) * perp);
    if (dot(a - c, cross(b - c, h1 - c)) < 0) {
        return {h2, h1};
    }
    return {h1, h2};
}

bool contains(std::span<const std::string_view> names, std::string_view name) {
    return std::find(names.begin(), names.end(), name) != names.end();
}

// --- the helices ---------------------------------------------------------------------------------

struct StrandSpec {
    char base;
    std::string_view rule;  // "" for strand 1
    int step;               // +1 with strand 1, -1 against it
};

std::vector<StrandSpec> helices(std::string_view kind) {
    if (kind == "reverse-wc") return {{'A', "", 1}, {'U', "reverse-wc", 1}};
    if (kind == "hoogsteen") return {{'A', "", 1}, {'U', "hoogsteen", 1}};
    if (kind == "reverse-hoogsteen") return {{'A', "", 1}, {'U', "reverse-hoogsteen", -1}};
    if (kind == "hh") return {{'A', "", 1}, {'A', "hh", 1}};
    if (kind == "triplex-hoogsteen") return {{'A', "", 1}, {'U', "wc", -1}, {'U', "hoogsteen", 1}};
    if (kind == "triplex-reverse-hoogsteen") return {{'A', "", 1}, {'U', "wc", -1}, {'A', "reverse-hoogsteen", -1}};
    throw std::logic_error("unknown helix kind " + std::string(kind));
}

std::vector<double> helixStart(std::string_view kind) {
    if (kind == "hh") return {1.00, 53.39, 172.39};
    if (kind == "reverse-wc")
        return {-10.46, -4.59, 54.04, 27.20, 22.08, 4.18, 27.75, 73.75, 121.99, 25.10, 99.98, 163.60};
    if (kind == "hoogsteen")
        return {3.60, -4.93, -26.01, 4.09, 36.63, 3.73, -5.87, 132.12, 165.56, 18.69, 117.52, 201.05};
    if (kind == "reverse-hoogsteen")
        return {-4.56, 0.76, 15.30, -16.81, 48.08, 4.39, -15.06, 162.79, 152.69, 20.07, 105.73, 215.23};
    if (kind == "triplex-hoogsteen")
        return {0.95,   -5.35, -23.66, 0.26,   25.61, 4.85,   3.15,  142.15,
                173.97, -5.92, 55.06,  170.43, 21.23, 135.34, 170.28};
    return {-0.44, -4.42, -26.00, 5.62, 21.95, 5.06, 8.96, 155.86, 157.34, -4.39, 52.09, 168.87, 18.85, 98.76, 116.35};
}

std::pair<double, double> screwPrior(std::string_view kind) {
    if (kind == "reverse-wc") return {20.66, 3.77};
    if (kind == "hoogsteen") return {37.98, 3.48};
    if (kind == "reverse-hoogsteen") return {45.7, 4.21};
    if (kind == "triplex-hoogsteen") return {36.19, 2.9};
    return {27.27, 4.18};
}

// _Screw: turn by -twist degrees about the axis through the origin and move -rise along it per level
struct Screw {
    Vec3 origin{};
    Vec3 axis{0.0, 0.0, 1.0};
    double twist = 0.0;
    double rise = 0.0;
};

// Screw.apply for one level, the matrix made once
struct ScrewLevel {
    Mat3 rot{};
    Vec3 origin{};
    Vec3 offset{};  // rise * level * u
    Vec3 operator()(const Vec3& x) const {
        return (matVec(rot, x - origin) + origin) - offset;
    }
};

ScrewLevel screwAt(const Screw& s, double level) {
    const Vec3 u = unit(s.axis);
    return {turnMatrix(u, -s.twist * level), s.origin, (s.rise * level) * u};
}

// --- nucleosides ---------------------------------------------------------------------------------

std::span<const std::string_view> baseSlotNames(char base) {
    static const std::map<char, std::vector<std::string_view>> slots = [] {
        std::map<char, std::vector<std::string_view>> out;
        for (const char b : std::string_view("ACGU")) {
            std::vector<std::string_view> names(baseAtomNames(b).begin(), baseAtomNames(b).end());
            for (const std::string_view h : baseHydrogenNames(b)) {
                names.push_back(h);
            }
            if (b == 'C') {
                names.emplace_back("H3+");
            }
            out[b] = std::move(names);
        }
        return out;
    }();
    return slots.at(base);
}

int baseCode(char base) {
    return base == 'A' ? 0 : base == 'C' ? 1 : base == 'G' ? 2 : 3;
}

// _core_template: the 4JRD sugar and backbone with the Arnott base superposed on the 4JRD glycosidic
// frame; the base hydrogens and the ring hydrogens of the sugar (_SUGAR_H); C: also H3+
const Atoms& coreTemplate(char base) {
    static const std::map<char, Atoms> cache = [] {
        std::map<char, Atoms> out;
        Atoms jrd;
        for (const reference::CylindricalAtom& a : reference::jrdNucleotide()) {
            jrd.set(a.name, cartesian(a.r, a.phi, a.z));
        }
        for (const char b : std::string_view("ACGU")) {
            const Atoms arn = templateAtoms(b, false);
            const Glycosidic ga = glycosidic('A');
            const Glycosidic gb = glycosidic(b);
            const std::array<Vec3, 4> mob = {arn.at("C1'"), arn.at(gb.n), arn.at(gb.a), arn.at(gb.b)};
            const std::array<Vec3, 4> ref = {jrd.at("C1'"), jrd.at(ga.n), jrd.at(ga.a), jrd.at(ga.b)};
            const std::array<double, 4> weights = {2.0, 2.0, 1.0, 1.0};
            const Kabsch k = kabsch(mob, ref, weights);
            Atoms heavy;
            for (const TemplateAtom& a : jrd.items) {
                if (!contains(baseAtomNames('A'), a.name)) {
                    heavy.items.push_back(a);
                }
            }
            for (const std::string_view name : baseAtomNames(b)) {
                heavy.set(name, kabschMove(k, arn.at(name)));
            }
            Atoms core = heavy;
            for (const TemplateAtom& h : hydrogens(heavy, b).items) {
                if (contains(baseHydrogenNames(b), h.name) || contains(kSugarH, h.name)) {
                    core.items.push_back(h);
                }
            }
            if (b == 'C') {
                const Vec3 n3 = heavy.at("N3");
                core.set("H3+", n3 - kBondNH * unit(unit(heavy.at("C2") - n3) + unit(heavy.at("C4") - n3)));
            }
            out[b] = std::move(core);
        }
        return out;
    }();
    return cache.at(base);
}

Atoms movedAtoms(const Atoms& coords, const Mat3& rot, const Vec3& shift) {
    Atoms out = coords;
    for (TemplateAtom& a : out.items) {
        a.xyz = transform(rot, shift, a.xyz);
    }
    return out;
}

Atoms screwed(const Atoms& coords, const ScrewLevel& op) {
    Atoms out = coords;
    for (TemplateAtom& a : out.items) {
        a.xyz = op(a.xyz);
    }
    return out;
}

// _base_plane: centre and normal (last right singular vector) of the base heavy atoms
std::pair<Vec3, Vec3> basePlane(const Atoms& coords, char base) {
    std::vector<Vec3> p;
    for (const std::string_view name : baseAtomNames(base)) {
        p.push_back(coords.at(name));
    }
    const Vec3 c = mean(p);
    for (Vec3& v : p) {
        v = v - c;
    }
    const core::SvdRows3 svd = core::svdRows3(p);
    return {c, svd.vt[2]};
}

// _base_normal
Vec3 baseNormal(const Atoms& coords, char base) {
    const Glycosidic g = glycosidic(base);
    return unit(cross(coords.at(g.a) - coords.at(g.n), coords.at(g.b) - coords.at(g.n)));
}

numerics::LeastSquaresOptions smallFitOptions() {
    numerics::LeastSquaresOptions o;
    o.ftol = o.xtol = o.gtol = kSmallFitTol;
    o.maxNfev = kSmallFitMaxNfev;
    return o;
}

struct Pose {
    Mat3 rot = identity();
    Vec3 shift{};
};

// _pair_fit
Pose pairFitUncached(std::string_view rule, char b1, char b2) {
    if (rule == "hh" && b1 == 'A' && b2 == 'A') {
        return {Mat3{Vec3{-1.0, 0.0, 0.0}, Vec3{0.0, -1.0, 0.0}, Vec3{0.0, 0.0, 1.0}}, Vec3{0.0, 0.0, 0.0}};
    }
    const std::span<const PairHbond> bonds = pairHbonds(rule, b1, b2);
    const Atoms& t1 = coreTemplate(b1);
    const Atoms& t2 = coreTemplate(b2);
    const auto [ring1, n1] = basePlane(t1, b1);
    const auto [ring2, n2] = basePlane(t2, b2);
    std::vector<Vec3> a, b;
    for (const PairHbond& h : bonds) {
        a.push_back(t1.at(h.atom1));
        b.push_back(t2.at(h.atom2));
    }
    Vec3 out = unit(cross(n1, a.back() - a.front()));
    if (dot(out, mean(a) - ring1) < 0) {
        out = -out;
    }
    std::vector<Vec3> targets;
    for (const Vec3& v : a) {
        targets.push_back(v + kHbond * out);
    }
    const Vec3 centre = mean(targets) + norm(mean(b) - ring2) * out;
    std::vector<std::string_view> names1(baseAtomNames(b1).begin(), baseAtomNames(b1).end());
    std::vector<std::string_view> names2(baseAtomNames(b2).begin(), baseAtomNames(b2).end());
    names1.emplace_back("C1'");
    names2.emplace_back("C1'");
    for (const std::string_view s : kSugar) {
        names1.push_back(s);
        names2.push_back(s);
    }
    std::vector<Vec3> fixedXyz, mobile, planar;
    for (const std::string_view n : names1) fixedXyz.push_back(t1.at(n));
    for (const std::string_view n : names2) mobile.push_back(t2.at(n));
    for (const std::string_view n : baseAtomNames(b2)) planar.push_back(t2.at(n));
    std::vector<std::uint8_t> freePair(names1.size() * names2.size(), 1);
    for (std::size_t i = 0; i < names1.size(); ++i) {
        for (std::size_t j = 0; j < names2.size(); ++j) {
            for (const PairHbond& h : bonds) {
                if (h.atom1 == names1[i] && h.atom2 == names2[j]) {
                    freePair[i * names2.size() + j] = 0;
                }
            }
        }
    }
    struct Bond {
        Vec3 fixedAtom, mobileAtom, hydrogen;
        bool mobileHydrogen;  // the donor is the partner
        int side;
    };
    std::vector<Bond> hb;
    for (const PairHbond& h : bonds) {
        hb.push_back({t1.at(h.atom1), t2.at(h.atom2), h.donor == 1 ? t1.at(h.hydrogen) : t2.at(h.hydrogen),
                      h.donor != 1, h.donor});
    }
    std::optional<std::pair<double, Pose>> best;
    for (const double sign : {1.0, -1.0}) {
        std::vector<Vec3> from = b, to = targets;
        from.push_back(ring2);
        from.push_back(ring2 + n2);
        to.push_back(centre);
        to.push_back(centre + sign * n1);
        std::vector<double> weights(b.size(), 10.0);
        weights.push_back(1.0);
        weights.push_back(1.0);
        const Kabsch k0 = kabsch(from, to, weights);
        const auto residuals = [&](std::span<const double> p) {
            const Mat3 rot = matMul(rotvec({p[0], p[1], p[2]}), k0.rotation);
            const Vec3 shift = {p[3], p[4], p[5]};
            const auto move = [&](const Vec3& x) {
                return (matVec(rot, x - k0.mobileCentre) + k0.referenceCentre) + shift;
            };
            std::vector<double> res;
            for (const Bond& bond : hb) {
                const Vec3 moved = move(bond.mobileAtom);
                const Vec3 acc = bond.side == 1 ? moved : bond.fixedAtom;
                const Vec3 don = bond.side == 1 ? bond.fixedAtom : moved;
                const Vec3 hyd = bond.mobileHydrogen ? move(bond.hydrogen) : bond.hydrogen;
                const double d = norm(don - acc);
                res.push_back((d - kHbond) * 4.0);
                res.push_back((d - norm(hyd - acc) - kBondNH) * 3.0);
            }
            for (const Vec3& x : planar) {
                res.push_back(dot(move(x) - ring1, n1));
            }
            std::vector<Vec3> mob(mobile.size());
            for (std::size_t j = 0; j < mobile.size(); ++j) {
                mob[j] = move(mobile[j]);
            }
            for (std::size_t i = 0; i < fixedXyz.size(); ++i) {
                for (std::size_t j = 0; j < mob.size(); ++j) {
                    const double dist = norm(fixedXyz[i] - mob[j]);
                    res.push_back((clip0(3.1 - dist) * freePair[i * mob.size() + j]) * 3.0);
                }
            }
            return res;
        };
        const std::vector<double> x0(6, 0.0);
        const numerics::LeastSquaresResult sol = numerics::leastSquares(residuals, x0, smallFitOptions());
        if (!best || sol.cost < best->first) {
            const Mat3 rot = matMul(rotvec({sol.x[0], sol.x[1], sol.x[2]}), k0.rotation);
            const Vec3 shift = (k0.referenceCentre + Vec3{sol.x[3], sol.x[4], sol.x[5]}) - matVec(rot, k0.mobileCentre);
            best = std::make_pair(sol.cost, Pose{rot, shift});
        }
    }
    return best->second;
}

std::mutex& cacheMutex() {
    static std::mutex m;
    return m;
}

Pose pairFit(std::string_view rule, char b1, char b2) {
    static std::map<std::string, Pose> cache;
    const std::string key = std::string(rule) + ":" + b1 + b2;
    {
        const std::lock_guard<std::mutex> lock(cacheMutex());
        const auto it = cache.find(key);
        if (it != cache.end()) {
            return it->second;
        }
    }
    const Pose pose = pairFitUncached(rule, b1, b2);
    const std::lock_guard<std::mutex> lock(cacheMutex());
    cache.emplace(key, pose);
    return pose;
}

// _nucleotides for one nucleoside: the sugar (_SUGAR, 7 atoms, or _RIGID, 12) turned by dchi about the
// glycosidic bond, O5' and P from gamma and beta, H5' / H5'' with the ring hydrogens. Names: the sugar
// atoms, C1', O5', P (, H5', H5'').
Atoms nucleotide(const Vec3& c1, const Vec3& glyco, std::span<const Vec3> sugar, double dchi, double gamma,
                 double beta) {
    const Vec3 k = unit(c1 - glyco);
    const double t = radians(dchi);
    const double ct = std::cos(t), st = std::sin(t);
    Atoms out;
    out.items.reserve(sugar.size() + 5);
    for (std::size_t i = 0; i < sugar.size(); ++i) {
        const Vec3 p = sugar[i] - c1;
        const Vec3 v = ((ct * p + st * cross(k, p)) + (dot(k, p) * (1 - ct)) * k) + c1;
        out.items.push_back({std::string(kRigid[i]), v});
    }
    out.items.push_back({"C1'", c1});
    const Vec3 o5 = placeAtom(out.at("C3'"), out.at("C4'"), out.at("C5'"), kC5O5, kC4C5O5, gamma);
    out.items.push_back({"O5'", o5});
    const Vec3 p = placeAtom(out.at("C4'"), out.at("C5'"), o5, kO5P, kC5O5P, beta);
    out.items.push_back({"P", p});
    if (sugar.size() > kSugar.size()) {
        const auto [h1, h2] = methyleneH(out.at("C5'"), o5, out.at("C4'"));
        out.items.push_back({"H5'", h1});
        out.items.push_back({"H5''", h2});
    }
    return out;
}

// _core_arrays of one nucleoside: C1', the glycosidic N, the sugar atoms (_RIGID or _SUGAR)
struct CoreArrays {
    Vec3 c1{};
    Vec3 glyco{};
    std::vector<Vec3> sugar;
};
CoreArrays coreArrays(const Atoms& placed, char base, bool withHydrogens) {
    CoreArrays out{placed.at("C1'"), placed.at(glycosidic(base).n), {}};
    const std::size_t count = withHydrogens ? kRigid.size() : kSugar.size();
    for (std::size_t i = 0; i < count; ++i) {
        out.sugar.push_back(placed.at(kRigid[i]));
    }
    return out;
}

Atoms nucleotideOf(const Atoms& placed, char base, bool withHydrogens, const Vec3& q) {
    const CoreArrays a = coreArrays(placed, base, withHydrogens);
    return nucleotide(a.c1, a.glyco, a.sugar, q[0], q[1], q[2]);
}

// The tilt of a link torsion delta degrees from its reference. MDATK's wrap(delta) / 120 jumps by 3
// where delta crosses +-180 degrees (a discontinuous cost); CH-17's is smooth and periodic, with
// MDATK's slope at the reference (1/120 per degree) and MDATK's values at the other staggered
// positions (+-1 at +-120 degrees): kTiltA sin(delta) + kTiltB sin(2 delta).
enum class TorsionTerm { Continuous, Mdatk };

const double kTiltB = (1.5 / kPi - 2.0 / std::sqrt(3.0)) / 3.0;
const double kTiltA = 2.0 / std::sqrt(3.0) + kTiltB;

double torsionTilt(double delta, TorsionTerm term) {
    if (term == TorsionTerm::Mdatk) {
        return wrap(delta) / 120.0;
    }
    const double r = radians(delta);
    return kTiltA * std::sin(r) + kTiltB * std::sin(2.0 * r);
}

// _link_residuals of one O3'(a) - P(b) link, the six terms
std::array<double, 6> linkResiduals(const Vec3& aC4, const Vec3& aC3, const Vec3& aO3, const Vec3& bP, const Vec3& bO5,
                                    const Vec3& bC5, TorsionTerm term) {
    const auto torsion = [term](double theta, double ref) {
        return (1 + std::cos(radians(3 * theta))) * 1.5 + torsionTilt(theta - ref, term);
    };
    const double eps = dihedral(aC4, aC3, aO3, bP);
    return {(norm(aO3 - bP) - kLink) / 0.01,
            (angle(aC3, aO3, bP) - kC3O3P) / 1.5,
            (angle(aO3, bP, bO5) - kO3PO5) / 1.5,
            torsion(eps, kEpsilon) + clip0(std::fabs(wrap(eps + 120.0)) - 60.0) / 10.0,
            torsion(dihedral(aC3, aO3, bP, bO5), kZeta),
            torsion(dihedral(aO3, bP, bO5, bC5), kAlpha)};
}

// _torsion_priors of one nucleotide: dchi, chi range, gamma, beta
std::array<double, 4> torsionPriors(double dchi, double gamma, double beta) {
    const double lo = kChiLow - kJrdChi;
    const double hi = kChiHigh - kJrdChi;
    return {dchi / 20.0, (clip0(lo - dchi) + clip0(dchi - hi)) / 1.0, wrap(gamma - kGamma) / 25.0,
            wrap(beta - kBeta) / 25.0};
}

// --- contacts ------------------------------------------------------------------------------------

bool polarH(std::string_view name) {
    static constexpr std::array<std::string_view, 11> polar = {"H1",  "H3",  "H21",  "H22",  "H41", "H42",
                                                               "H61", "H62", "HO2'", "HO5'", "HO3'"};
    return contains(polar, name);
}

// _atom_class
int atomClass(char base, std::string_view name) {
    if (name.empty() || name == "H3+") {
        return 6;
    }
    if (name[0] == 'H') {
        return polarH(name) ? 4 : 5;
    }
    if (name == "O2'" || name == "HO5'" || name == "HO3'") {
        return 3;
    }
    static const std::map<char, std::vector<std::string_view>> acceptors = {
        {'A', {"N1", "N3", "N7"}}, {'G', {"N3", "N7"}}, {'C', {"N3"}}, {'U', {}}};
    static const std::map<char, std::vector<std::string_view>> donors = {
        {'A', {"N6"}}, {'G', {"N1", "N2"}}, {'C', {"N4"}}, {'U', {"N3"}}};
    if (name[0] == 'O' || contains(acceptors.at(base), name)) {
        return 2;
    }
    return contains(donors.at(base), name) ? 1 : 0;
}

using LimitMatrix = std::array<std::array<double, 8>, 8>;

// _limit_matrix
LimitMatrix limitMatrix(const ContactTable& table) {
    LimitMatrix out{};
    const auto hyd = [](int c) { return c == 4 || c == 5; };
    for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
            double& v = out[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
            if (a == 6 || b == 6) {
                v = 0.0;
            } else if (hyd(a) && hyd(b)) {
                v = table.hh;
            } else if (hyd(a) || hyd(b)) {
                const int h = hyd(a) ? a : b;
                const int x = hyd(a) ? b : a;
                v = h == 4 && (x == 2 || x == 3 || x == 7) ? table.polar : table.hx;
            } else {
                const bool polarPair =
                    a >= 1 && a <= 3 && b >= 1 && b <= 3 && !(a == 1 && b == 1) && !(a == 2 && b == 2);
                v = a == 7 || b == 7 || polarPair ? table.hbond : table.heavy;
            }
        }
    }
    return out;
}

const LimitMatrix& fitLimits() {
    static const LimitMatrix m = limitMatrix(kContactFit);
    return m;
}
const LimitMatrix& reportLimits() {
    static const LimitMatrix m = limitMatrix(kContactReport);
    return m;
}
double limitOf(const LimitMatrix& m, int a, int b) {
    return m[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
}

// _residue_graph: atom index and bond-path lengths (99: none, at most 15) of the Arnott nucleotide
// with its hydrogens
struct ResidueGraph {
    std::vector<std::string> names;
    std::vector<int> dist;  // n x n
    int index(std::string_view name) const {
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (names[i] == name) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }
};

const ResidueGraph& residueGraph(char base) {
    static const std::map<char, ResidueGraph> graphs = [] {
        std::map<char, ResidueGraph> out;
        for (const char b : std::string_view("ACGU")) {
            const Atoms t = templateAtoms(b, true);
            ResidueGraph g;
            const std::size_t n = t.items.size();
            for (const TemplateAtom& a : t.items) {
                g.names.push_back(a.name);
            }
            std::vector<std::vector<std::size_t>> bonded(n);
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    if (i == j) continue;
                    const bool hi = g.names[i][0] == 'H', hj = g.names[j][0] == 'H';
                    const double d = norm(t.items[i].xyz - t.items[j].xyz);
                    if (d < 1.75 && !(hi && hj) && !((hi || hj) && d > 1.25)) {
                        bonded[i].push_back(j);
                    }
                }
            }
            g.dist.assign(n * n, 99);
            for (std::size_t s = 0; s < n; ++s) {  // breadth first, at most 15 bonds
                std::vector<std::size_t> frontier = {s};
                g.dist[s * n + s] = 0;
                for (int k = 1; k < 16 && !frontier.empty(); ++k) {
                    std::vector<std::size_t> next;
                    for (const std::size_t u : frontier) {
                        for (const std::size_t v : bonded[u]) {
                            if (g.dist[s * n + v] == 99) {
                                g.dist[s * n + v] = k;
                                next.push_back(v);
                            }
                        }
                    }
                    frontier = std::move(next);
                }
            }
            out[b] = std::move(g);
        }
        return out;
    }();
    return graphs.at(base);
}

int bondPath(char base, std::string_view a, std::string_view b) {
    const ResidueGraph& g = residueGraph(base);
    const int i = g.index(a), j = g.index(b);
    if (i < 0 || j < 0) {
        return 99;
    }
    return g.dist[static_cast<std::size_t>(i) * g.names.size() + static_cast<std::size_t>(j)];
}

// _link_path: from atom a of a nucleotide to atom b of the next one (O3'-P link)
int linkPath(char base1, std::string_view a, char base2, std::string_view b) {
    return bondPath(base1, a, "O3'") + 1 + bondPath(base2, "P", b);
}

// --- the reference helix (_helix_fit) ------------------------------------------------------------------
//
// The reference helix of a kind (its screw and every strand's chi change, gamma and beta) is the one
// input the builds take from it: the uniform helices are it, and every closure (CH-16) starts on it.
// It is separable from the closure: what defines it is helixStart() (the fit's start), screwPrior()
// and the torsion term of its fit (kReferenceTerm: CH-17's continuous term, or MDATK's, which gives
// MDATK's reference helices); changing those changes the reference geometry without touching the
// closure. Whether MDATK's or the refitted reference helices are the right ones is open (docs/
// DIFFERENCES.md CH-17).
constexpr TorsionTerm kReferenceTerm = TorsionTerm::Continuous;

struct ReferenceLevel {
    std::vector<Atoms> placed;
    Vec3 centre{}, n{}, e1{}, e2{};
};

// _reference_level
ReferenceLevel referenceLevel(std::string_view kind) {
    const std::vector<StrandSpec> strands = helices(kind);
    ReferenceLevel out;
    out.placed.push_back(coreTemplate(strands[0].base));
    for (std::size_t s = 1; s < strands.size(); ++s) {
        const Pose p = pairFit(strands[s].rule, strands[0].base, strands[s].base);
        out.placed.push_back(movedAtoms(coreTemplate(strands[s].base), p.rot, p.shift));
    }
    out.n = -baseNormal(out.placed[0], strands[0].base);
    out.e1 = unit(cross(out.n, {1.0, 0.0, 0.0}));
    std::vector<Vec3> c1;
    for (const Atoms& a : out.placed) {
        c1.push_back(a.at("C1'"));
    }
    out.centre = mean(c1);
    out.e2 = cross(out.n, out.e1);
    return out;
}

// _level_atoms: names, owner strand and xyz of one level of the uniform helix
struct LevelAtoms {
    std::vector<std::string> names;
    std::vector<std::size_t> owner;
    std::vector<Vec3> xyz;
};
LevelAtoms levelAtoms(const std::vector<Atoms>& placed, const std::vector<StrandSpec>& strands,
                      const std::vector<Vec3>& q, const Screw& screw) {
    LevelAtoms out;
    for (std::size_t s = 0; s < strands.size(); ++s) {
        Atoms atoms = nucleotideOf(placed[s], strands[s].base, true, q[s]);
        const Vec3 prevO3 = screwAt(screw, -strands[s].step)(atoms.at("O3'"));
        const auto [op1, op2] = phosphateOxygens(atoms.at("P"), prevO3, atoms.at("O5'"));
        atoms.set("OP1", op1);
        atoms.set("OP2", op2);
        for (const std::string_view name : baseSlotNames(strands[s].base)) {
            if (name != "H3+") {
                atoms.set(name, placed[s].at(name));
            }
        }
        for (const TemplateAtom& a : atoms.items) {
            out.names.push_back(a.name);
            out.owner.push_back(s);
            out.xyz.push_back(a.xyz);
        }
    }
    return out;
}

struct ContactPairs {  // _contact_pairs
    std::vector<std::size_t> i, j;
    std::vector<int> level;
    std::vector<double> limit;
    std::vector<std::uint8_t> heavy, hh;
};

struct UnpackedFit {
    Screw screw;
    std::vector<Vec3> q;
};

// one reference helix kind: its level, its contact pairs, its fit
struct HelixKind {
    std::string kind;
    std::vector<StrandSpec> strands;
    ReferenceLevel ref;
    ContactPairs pairs;
};

ContactPairs contactPairs(const HelixKind& h) {
    const std::vector<Vec3> q(h.strands.size(), Vec3{0.0, 58.0, 174.0});
    Screw dummy;
    dummy.twist = 30.0;
    dummy.rise = 3.0;
    const LevelAtoms la = levelAtoms(h.ref.placed, h.strands, q, dummy);
    const std::size_t count = la.names.size();
    std::vector<char> bases(count);
    std::vector<int> cls(count);
    for (std::size_t i = 0; i < count; ++i) {
        bases[i] = h.strands[la.owner[i]].base;
        cls[i] = atomClass(bases[i], la.names[i]);
    }
    std::set<std::tuple<std::size_t, std::string, std::size_t, std::string>> hbonded;
    for (std::size_t s = 1; s < h.strands.size(); ++s) {
        for (const PairHbond& hb : pairHbonds(h.strands[s].rule, h.strands[0].base, h.strands[s].base)) {
            hbonded.insert({0, std::string(hb.atom1), s, std::string(hb.atom2)});
            hbonded.insert({s, std::string(hb.atom2), 0, std::string(hb.atom1)});
        }
    }
    ContactPairs out;
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = 0; j < count; ++j) {
            for (int level = 0; level < 3; ++level) {
                if (level == 0 && j <= i) {
                    continue;
                }
                const std::size_t s = la.owner[i], t = la.owner[j];
                int path = 99;
                if (s == t && level == 0) {
                    path = bondPath(bases[i], la.names[i], la.names[j]);
                } else if (s == t && level == 1) {
                    path = h.strands[s].step > 0 ? linkPath(bases[i], la.names[i], bases[j], la.names[j])
                                                 : linkPath(bases[j], la.names[j], bases[i], la.names[i]);
                }
                if (path < 4 || (level == 0 && hbonded.count({s, la.names[i], t, la.names[j]}))) {
                    continue;
                }
                const double limit = limitOf(fitLimits(), cls[i], cls[j]);
                if (limit > 0) {
                    out.i.push_back(i);
                    out.j.push_back(j);
                    out.level.push_back(level);
                    out.limit.push_back(limit);
                    out.heavy.push_back(cls[i] < 4 && cls[j] < 4);
                    out.hh.push_back((cls[i] == 4 || cls[i] == 5) && (cls[j] == 4 || cls[j] == 5));
                }
            }
        }
    }
    return out;
}

// _unpack_fit
UnpackedFit unpackFit(const HelixKind& h, std::span<const double> p) {
    UnpackedFit out;
    if (h.kind == "hh") {
        out.screw = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, reference::jrdTwist(), reference::jrdRise()};
        out.q = {Vec3{p[0], p[1], p[2]}, Vec3{p[0], p[1], p[2]}};
        return out;
    }
    const ReferenceLevel& r = h.ref;
    const Vec3 axis = matVec(matMul(turnMatrix(r.e1, p[2]), turnMatrix(r.e2, p[3])), r.n);
    out.screw = {(r.centre + p[0] * r.e1) + p[1] * r.e2, axis, p[4], p[5]};
    for (std::size_t s = 0; s < h.strands.size(); ++s) {
        out.q.push_back({p[6 + 3 * s], p[7 + 3 * s], p[8 + 3 * s]});
    }
    return out;
}

// _level_distances
std::vector<double> levelDistances(const HelixKind& h, const UnpackedFit& u) {
    const LevelAtoms la = levelAtoms(h.ref.placed, h.strands, u.q, u.screw);
    const ScrewLevel one = screwAt(u.screw, 1), two = screwAt(u.screw, 2);
    std::vector<std::array<Vec3, 3>> levels(la.xyz.size());
    for (std::size_t i = 0; i < la.xyz.size(); ++i) {
        levels[i] = {la.xyz[i], one(la.xyz[i]), two(la.xyz[i])};
    }
    std::vector<double> d(h.pairs.i.size());
    for (std::size_t k = 0; k < d.size(); ++k) {
        d[k] = norm(la.xyz[h.pairs.i[k]] - levels[h.pairs.j[k]][static_cast<std::size_t>(h.pairs.level[k])]);
    }
    return d;
}

// _fit_residuals
std::vector<double> fitResiduals(const HelixKind& h, std::span<const double> p) {
    const UnpackedFit u = unpackFit(h, p);
    std::vector<double> out;
    for (std::size_t s = 0; s < h.strands.size(); ++s) {
        const Atoms here = nucleotideOf(h.ref.placed[s], h.strands[s].base, false, u.q[s]);
        const ScrewLevel op = screwAt(u.screw, h.strands[s].step);
        const std::array<double, 6> link =
            linkResiduals(here.at("C4'"), here.at("C3'"), here.at("O3'"), op(here.at("P")), op(here.at("O5'")),
                          op(here.at("C5'")), kReferenceTerm);
        out.insert(out.end(), link.begin(), link.end());
        const std::array<double, 4> prior = torsionPriors(u.q[s][0], u.q[s][1], u.q[s][2]);
        out.insert(out.end(), prior.begin(), prior.end());
    }
    if (h.kind != "hh") {
        const auto [twist, rise] = screwPrior(h.kind);
        out.push_back(p[2] / 25.0);
        out.push_back(p[3] / 25.0);
        out.push_back((p[4] - twist) / 12.0);
        out.push_back((p[5] - rise) / 0.3);
    }
    const std::vector<double> d = levelDistances(h, u);
    for (std::size_t k = 0; k < d.size(); ++k) {
        out.push_back(clip0(h.pairs.limit[k] - d[k]) / kFitContactWeight);
    }
    return out;
}

struct HelixFitResult {
    Screw screw;
    std::vector<Vec3> torsions;  // per strand: chi change, gamma, beta
    double cost = 0.0;
    double contact = 0.0;
    double contactH = 0.0;
};

HelixFitResult helixFitUncached(const std::string& kind, const core::Reporter* reporter) {
    HelixKind h;
    h.kind = kind;
    h.strands = helices(kind);
    h.ref = referenceLevel(kind);
    h.pairs = contactPairs(h);
    const std::vector<double> start = helixStart(kind);
    numerics::LeastSquaresOptions o = smallFitOptions();
    if (kind == "hh") {
        o.xScale = {10.0, 10.0, 10.0};
    } else {
        o.xScale = {1.0, 1.0, 5.0, 5.0, 5.0, 0.3};
        for (std::size_t k = 0; k < 3 * h.strands.size(); ++k) {
            o.xScale.push_back(10.0);
        }
    }
    const numerics::LeastSquaresCallback callback = [reporter](std::span<const double>, double, std::size_t,
                                                               std::size_t) {
        if (reporter) {
            reporter->checkCancel();
        }
        return false;
    };
    const numerics::LeastSquaresResult sol =
        numerics::leastSquares([&h](std::span<const double> p) { return fitResiduals(h, p); }, start, o, callback);
    const UnpackedFit u = unpackFit(h, sol.x);
    HelixFitResult out{u.screw, u.q, sol.cost, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity()};
    const std::vector<double> d = levelDistances(h, u);
    for (std::size_t k = 0; k < d.size(); ++k) {
        if (h.pairs.heavy[k]) out.contact = std::min(out.contact, d[k]);
        if (h.pairs.hh[k]) out.contactH = std::min(out.contactH, d[k]);
    }
    return out;
}

HelixFitResult helixFit(const std::string& kind, const core::Reporter* reporter) {
    static std::map<std::string, HelixFitResult> cache;
    static std::mutex fitMutex;  // one fit of a kind at a time (it is cached)
    const std::lock_guard<std::mutex> fitLock(fitMutex);
    {
        const std::lock_guard<std::mutex> lock(cacheMutex());
        const auto it = cache.find(kind);
        if (it != cache.end()) {
            return it->second;
        }
    }
    HelixFitResult fit = helixFitUncached(kind, reporter);
    const std::lock_guard<std::mutex> lock(cacheMutex());
    cache.emplace(kind, fit);
    return fit;
}

// --- level poses (_level_pose) -------------------------------------------------------------------

Pose levelPose(const std::string& kind, std::size_t strand, const std::string& level);

Pose levelPoseUncached(const std::string& kind, std::size_t strand, const std::string& level) {
    const std::vector<StrandSpec> spec = helices(kind);
    const std::string_view rule = spec[strand].rule;
    const char a = level[0], b = level[strand];
    if (canonical(rule, a, b)) {
        return pairFit(rule, a, b);
    }
    const Pose p0 = pairFit(rule, spec[0].base, spec[strand].base);
    std::vector<Atoms> others = {coreTemplate(a)};
    for (std::size_t t = 1; t < strand; ++t) {
        const Pose pt = levelPose(kind, t, level.substr(0, t + 1));
        others.push_back(movedAtoms(coreTemplate(level[t]), pt.rot, pt.shift));
    }
    const auto selected = [](const Atoms& coords, char base, std::vector<Vec3>& xyz, std::vector<int>& cls) {
        for (const TemplateAtom& x : coords.items) {
            const bool keep = contains(kRigid, x.name) || x.name == "C1'" ||
                              (contains(baseSlotNames(base), x.name) && x.name != "H3+");
            if (keep) {
                xyz.push_back(x.xyz);
                cls.push_back(atomClass(base, x.name));
            }
        }
    };
    std::vector<Vec3> fixedXyz;
    std::vector<int> fixedCls;
    for (std::size_t t = 0; t < others.size(); ++t) {
        selected(others[t], level[t], fixedXyz, fixedCls);
    }
    const Atoms start = movedAtoms(coreTemplate(b), p0.rot, p0.shift);
    std::vector<Vec3> mobile;
    std::vector<int> mobileCls;
    selected(start, b, mobile, mobileCls);
    std::vector<double> limit(mobile.size() * fixedXyz.size());
    for (std::size_t i = 0; i < mobile.size(); ++i) {
        for (std::size_t j = 0; j < fixedXyz.size(); ++j) {
            limit[i * fixedXyz.size() + j] = limitOf(fitLimits(), mobileCls[i], fixedCls[j]);
        }
    }
    std::vector<Vec3> planar;
    for (const std::string_view name : baseAtomNames(b)) {
        planar.push_back(start.at(name));
    }
    const auto [ring1, n1] = basePlane(others[0], a);
    const Vec3 centre = start.at("C1'");
    const auto residuals = [&](std::span<const double> p) {
        const Mat3 r = rotvec({radians(p[0]), radians(p[1]), radians(p[2])});
        const Vec3 shift = {p[3], p[4], p[5]};
        const auto move = [&](const Vec3& x) { return (matVec(r, x - centre) + centre) + shift; };
        std::vector<double> res = {p[0] / 10.0, p[1] / 10.0, p[2] / 10.0, p[3] / 0.5, p[4] / 0.5, p[5] / 0.5};
        for (const Vec3& x : planar) {
            res.push_back(dot(move(x) - ring1, n1) / 0.25);
        }
        for (std::size_t i = 0; i < mobile.size(); ++i) {
            const Vec3 m = move(mobile[i]);
            for (std::size_t j = 0; j < fixedXyz.size(); ++j) {
                res.push_back(clip0(limit[i * fixedXyz.size() + j] - norm(m - fixedXyz[j])) / 0.05);
            }
        }
        return res;
    };
    std::optional<numerics::LeastSquaresResult> best;
    for (const double turn : {0.0, 25.0, -25.0, 50.0, -50.0}) {  // in-plane turns about C1' as starts
        const double rt = radians(turn);
        const std::vector<double> x0 = {degrees(n1[0] * rt), degrees(n1[1] * rt), degrees(n1[2] * rt), 0.0, 0.0, 0.0};
        numerics::LeastSquaresResult sol = numerics::leastSquares(residuals, x0, smallFitOptions());
        if (!best || sol.cost < best->cost - 1e-9) {
            best = std::move(sol);
        }
    }
    const Mat3 r = rotvec({radians(best->x[0]), radians(best->x[1]), radians(best->x[2])});
    return {matMul(r, p0.rot), (matVec(r, p0.shift - centre) + centre) + Vec3{best->x[3], best->x[4], best->x[5]}};
}

Pose levelPose(const std::string& kind, std::size_t strand, const std::string& level) {
    static std::map<std::string, Pose> cache;
    const std::string key = kind + ":" + std::to_string(strand) + ":" + level;
    {
        const std::lock_guard<std::mutex> lock(cacheMutex());
        const auto it = cache.find(key);
        if (it != cache.end()) {
            return it->second;
        }
    }
    const Pose pose = levelPoseUncached(kind, strand, level);
    const std::lock_guard<std::mutex> lock(cacheMutex());
    cache.emplace(key, pose);
    return pose;
}

// --- the strand closure (_close_levels, CH-16) ---------------------------------------------------------

// _slot_table / _contact_tables: slot classes, slot pairs scored within a nucleotide, across a link
struct ContactTables {
    std::array<std::array<int, kSlots>, 4> cls{};
    std::array<std::vector<std::uint8_t>, 4> inside;               // kSlots x kSlots
    std::array<std::array<std::vector<std::uint8_t>, 4>, 4> link;  // [5'][3'] kSlots x kSlots
};

std::array<std::string_view, kSlots> slotNames(char base) {
    std::array<std::string_view, kSlots> names{};
    std::size_t k = 0;
    for (const std::string_view n : kStrandNames) names[k++] = n;
    for (const std::string_view n : baseSlotNames(base)) names[k++] = n;
    return names;
}

const ContactTables& contactTables() {
    static const ContactTables tables = [] {
        ContactTables t;
        const std::string_view order = "ACGU";
        std::array<std::array<std::string_view, kSlots>, 4> names{};
        for (std::size_t x = 0; x < 4; ++x) {
            const char b = order[x];
            names[x] = slotNames(b);
            for (std::size_t i = 0; i < kSlots; ++i) {
                t.cls[x][i] = atomClass(b, names[x][i]);
            }
            t.inside[x].assign(kSlots * kSlots, 0);
            for (std::size_t i = 0; i < kSlots; ++i) {
                for (std::size_t k = i + 1; k < kSlots; ++k) {
                    const std::string_view a = names[x][i], c = names[x][k];
                    const bool rigidA = contains(kRigid, a) || a == "C1'";
                    const bool rigidC = contains(kRigid, c) || c == "C1'";
                    const bool ringA = contains(baseSlotNames(b), a), ringC = contains(baseSlotNames(b), c);
                    if (t.cls[x][i] == 6 || t.cls[x][k] == 6 || (rigidA && rigidC) || (ringA && ringC)) {
                        continue;
                    }
                    const std::uint8_t v = bondPath(b, a, c) >= 4;
                    t.inside[x][i * kSlots + k] = t.inside[x][k * kSlots + i] = v;
                }
            }
        }
        for (std::size_t x = 0; x < 4; ++x) {
            for (std::size_t y = 0; y < 4; ++y) {
                std::vector<std::uint8_t>& l = t.link[x][y];
                l.assign(kSlots * kSlots, 0);
                for (std::size_t i = 0; i < kSlots; ++i) {
                    for (std::size_t k = 0; k < kSlots; ++k) {
                        l[i * kSlots + k] = t.cls[x][i] != 6 && t.cls[y][k] != 6 &&
                                            linkPath(order[x], names[x][i], order[y], names[y][k]) >= 4;
                    }
                }
            }
        }
        return t;
    }();
    return tables;
}

struct ClosurePair {
    std::size_t p, q;
    double limit;
    std::size_t lo, hi;
};

class Closure {
public:
    Closure(const std::string& kind, std::vector<std::vector<Atoms>> placed, const std::vector<std::string>& bases,
            const std::vector<int>& steps, const std::vector<Vec3>& start, const Screw& screw,
            const HelixOptions& options)
        : m_placed(std::move(placed)), m_bases(bases), m_steps(steps), m_options(options) {
        m_ns = m_placed.size();
        m_n = m_placed[0].size();
        m_width = 6 + 3 * m_ns;
        const std::vector<StrandSpec> spec = helices(kind);
        // per partner nucleoside: the reference pose and the rigid move to its _level_pose pose
        m_moves.resize(m_ns);
        for (std::size_t s = 1; s < m_ns; ++s) {
            const Pose ref = pairFit(spec[s].rule, spec[0].base, spec[s].base);
            for (std::size_t k = 0; k < m_n; ++k) {
                Move mv;
                mv.ref = screwed(movedAtoms(coreTemplate(m_bases[s][k]), ref.rot, ref.shift),
                                 screwAt(screw, static_cast<double>(k)));
                std::vector<Vec3> P, Q;
                for (const TemplateAtom& a : m_placed[s][k].items) {
                    if (a.name != "H3+") {
                        if (const Vec3* r = mv.ref.find(a.name)) {
                            P.push_back(*r);
                            Q.push_back(a.xyz);
                        }
                    }
                }
                const std::vector<double> ones(P.size(), 1.0);
                const Kabsch kb = kabsch(P, Q, ones);
                mv.c = mv.ref.at("C1'");
                mv.w = rotationVector(kb.rotation);
                mv.d = (matVec(kb.rotation, mv.c - kb.mobileCentre) + kb.referenceCentre) - mv.c;
                m_moves[s].push_back(std::move(mv));
            }
        }
        const ContactTables& tables = contactTables();
        m_codes.resize(m_ns * m_n);
        m_cls.resize(m_ns * m_n * kSlots);
        for (std::size_t s = 0; s < m_ns; ++s) {
            for (std::size_t k = 0; k < m_n; ++k) {
                const int code = baseCode(m_bases[s][k]);
                m_codes[s * m_n + k] = code;
                for (std::size_t i = 0; i < kSlots; ++i) {
                    m_cls[(s * m_n + k) * kSlots + i] = tables.cls[static_cast<std::size_t>(code)][i];
                }
            }
        }
        for (std::size_t s = 0; s < m_ns; ++s) {
            std::vector<std::size_t> o(m_n);
            for (std::size_t k = 0; k < m_n; ++k) {
                o[k] = m_steps[s] > 0 ? k : m_n - 1 - k;
            }
            m_order.push_back(std::move(o));
        }
        m_x.assign(m_n * m_width, 0.0);
        for (std::size_t k = 0; k < m_n; ++k) {
            for (std::size_t s = 0; s < m_ns; ++s) {
                for (std::size_t c = 0; c < 3; ++c) {
                    m_x[k * m_width + 6 + 3 * s + c] = start[s][c];
                }
            }
        }
        m_scale.resize(m_n * m_width);
        for (std::size_t k = 0; k < m_n; ++k) {
            for (std::size_t c = 0; c < m_width; ++c) {
                m_scale[k * m_width + c] = c < 3 ? 3.0 : c < 6 ? 0.2 : 10.0;
            }
        }
        if (m_options.reporter) {
            m_reporter = *m_options.reporter;
        }
    }

    // the continuation: the t = 0 solve, then kSteps steps (halved where they have trouble)
    void run() {
        m_x = step(m_x, 0.0, kClosureFinalNfev, false).first;
        for (int i = 0; i < kSteps; ++i) {
            m_x = advance(m_x, static_cast<double>(i) / kSteps, static_cast<double>(i + 1) / kSteps, 0);
        }
        setGeometry(1.0);
    }

    // torsions of strand s by level
    std::vector<Vec3> torsions(std::size_t s) const {
        std::vector<Vec3> out(m_n);
        for (std::size_t k = 0; k < m_n; ++k) {
            for (std::size_t c = 0; c < 3; ++c) {
                out[k][c] = m_x[k * m_width + 6 + 3 * s + c];
            }
        }
        return out;
    }

    // the placed nucleosides moved with their level
    std::vector<std::vector<Atoms>> moved() const {
        std::vector<std::vector<Atoms>> out = m_placed;
        for (std::size_t k = 0; k < m_n; ++k) {
            const Mat3 rot = levelRotation(m_x, k);
            const Vec3 shift = {m_x[k * m_width + 3], m_x[k * m_width + 4], m_x[k * m_width + 5]};
            for (std::size_t s = 0; s < m_ns; ++s) {
                for (TemplateAtom& a : out[s][k].items) {
                    a.xyz = (matVec(rot, a.xyz - m_centre[k]) + m_centre[k]) + shift;
                }
            }
        }
        return out;
    }

private:
    struct Move {
        Atoms ref;
        Vec3 w{}, c{}, d{};
    };
    struct Geometry {  // _core_arrays and the base slots of every (strand, level) at the current t
        std::vector<Vec3> c1, glyco;
        std::vector<std::array<Vec3, 12>> sugar;
        std::vector<std::array<Vec3, kBaseSlots>> base;
    };

    Mat3 levelRotation(const std::vector<double>& x, std::size_t k) const {
        return rodrigues({radians(x[k * m_width]), radians(x[k * m_width + 1]), radians(x[k * m_width + 2])});
    }

    Atoms placedAt(std::size_t s, std::size_t k, double t) const {
        if (s == 0 || t == 1.0) {
            return m_placed[s][k];
        }
        const Move& mv = m_moves[s][k];
        const Mat3 r = rotvec(t * mv.w);
        Atoms out = mv.ref;
        for (TemplateAtom& a : out.items) {
            a.xyz = (matVec(r, a.xyz - mv.c) + mv.c) + t * mv.d;
        }
        return out;
    }

    void setGeometry(double t) {
        Geometry g;
        const std::size_t count = m_ns * m_n;
        g.c1.resize(count);
        g.glyco.resize(count);
        g.sugar.resize(count);
        g.base.resize(count);
        for (std::size_t s = 0; s < m_ns; ++s) {
            for (std::size_t k = 0; k < m_n; ++k) {
                const Atoms a = placedAt(s, k, t);
                const char b = m_bases[s][k];
                const std::size_t at = s * m_n + k;
                g.c1[at] = a.at("C1'");
                g.glyco[at] = a.at(glycosidic(b).n);
                for (std::size_t i = 0; i < 12; ++i) {
                    g.sugar[at][i] = a.at(kRigid[i]);
                }
                const std::span<const std::string_view> slots = baseSlotNames(b);
                for (std::size_t i = 0; i < kBaseSlots; ++i) {
                    g.base[at][i] = i < slots.size() ? a.at(slots[i]) : Vec3{1e4, 1e4, 1e4};
                }
            }
        }
        m_centre.assign(m_n, Vec3{0.0, 0.0, 0.0});
        for (std::size_t k = 0; k < m_n; ++k) {
            Vec3 sum{0.0, 0.0, 0.0};
            for (std::size_t s = 0; s < m_ns; ++s) {
                sum = sum + g.c1[s * m_n + k];
            }
            const double ns = static_cast<double>(m_ns);
            m_centre[k] = {sum[0] / ns, sum[1] / ns, sum[2] / ns};
        }
        m_geometry = std::move(g);
    }

    // build: every atom slot of every nucleotide (strand, level, slot) for the variables x
    std::vector<Vec3> build(std::span<const double> x) const {
        std::vector<Vec3> atoms(m_ns * m_n * kSlots);
        std::vector<Mat3> rot(m_n);
        for (std::size_t k = 0; k < m_n; ++k) {
            rot[k] = rodrigues({radians(x[k * m_width]), radians(x[k * m_width + 1]), radians(x[k * m_width + 2])});
        }
        for (std::size_t s = 0; s < m_ns; ++s) {
            for (std::size_t k = 0; k < m_n; ++k) {
                const std::size_t at = s * m_n + k;
                const Vec3 shift = {x[k * m_width + 3], x[k * m_width + 4], x[k * m_width + 5]};
                const auto move = [&](const Vec3& y) {
                    return (matVec(rot[k], y - m_centre[k]) + m_centre[k]) + shift;
                };
                std::array<Vec3, 12> sugar{};
                for (std::size_t i = 0; i < 12; ++i) {
                    sugar[i] = move(m_geometry.sugar[at][i]);
                }
                const Vec3 c1 = move(m_geometry.c1[at]);
                const Vec3 glyco = move(m_geometry.glyco[at]);
                const std::size_t v = k * m_width + 6 + 3 * s;
                // _nucleotides
                const Vec3 axis = unit(c1 - glyco);
                const double t = radians(x[v]);
                const double ct = std::cos(t), st = std::sin(t);
                Vec3* slot = &atoms[at * kSlots];
                for (std::size_t i = 0; i < 12; ++i) {
                    const Vec3 p = sugar[i] - c1;
                    slot[kRigidSlot[i]] = ((ct * p + st * cross(axis, p)) + (dot(axis, p) * (1 - ct)) * axis) + c1;
                }
                slot[kSlotC1] = c1;
                slot[kSlotO5] = placeAtom(slot[kSlotC3], slot[kSlotC4], slot[kSlotC5], kC5O5, kC4C5O5, x[v + 1]);
                slot[kSlotP] = placeAtom(slot[kSlotC4], slot[kSlotC5], slot[kSlotO5], kO5P, kC5O5P, x[v + 2]);
                const auto [h1, h2] = methyleneH(slot[kSlotC5], slot[kSlotO5], slot[kSlotC4]);
                slot[kSlotH5] = h1;
                slot[kSlotH5b] = h2;
                for (std::size_t i = 0; i < kBaseSlots; ++i) {
                    slot[kStrandNames.size() + i] = move(m_geometry.base[at][i]);
                }
            }
            const std::vector<std::size_t>& o = m_order[s];
            for (std::size_t j = 0; j < m_n; ++j) {
                Vec3* slot = &atoms[(s * m_n + o[j]) * kSlots];
                const Vec3 prev = j > 0 ? atoms[(s * m_n + o[j - 1]) * kSlots + kSlotO3]
                                        : placeAtom(slot[kSlotC5], slot[kSlotO5], slot[kSlotP], kLink, kO3PO5, kAlpha);
                const auto [op1, op2] = phosphateOxygens(slot[kSlotP], prev, slot[kSlotO5]);
                slot[kSlotOP1] = op1;
                slot[kSlotOP2] = op2;
            }
        }
        return atoms;
    }

    std::vector<double> residuals(std::span<const double> x, const std::vector<ClosurePair>& pairs) const {
        const std::vector<Vec3> atoms = build(x);
        std::vector<double> out;
        out.reserve(m_ns * (10 * m_n) + 6 * m_n + pairs.size());
        for (std::size_t s = 0; s < m_ns; ++s) {
            const std::vector<std::size_t>& o = m_order[s];
            const std::size_t links = m_n - 1;
            std::vector<std::array<double, 6>> terms(links);
            for (std::size_t j = 0; j < links; ++j) {
                const Vec3* a = &atoms[(s * m_n + o[j]) * kSlots];
                const Vec3* b = &atoms[(s * m_n + o[j + 1]) * kSlots];
                terms[j] = linkResiduals(a[kSlotC4], a[kSlotC3], a[kSlotO3], b[kSlotP], b[kSlotO5], b[kSlotC5],
                                         TorsionTerm::Continuous);
            }
            for (std::size_t c = 0; c < 6; ++c) {
                for (std::size_t j = 0; j < links; ++j) {
                    out.push_back(terms[j][c]);
                }
            }
        }
        for (std::size_t s = 0; s < m_ns; ++s) {
            std::vector<std::array<double, 4>> priors(m_n);
            for (std::size_t k = 0; k < m_n; ++k) {
                const std::size_t v = k * m_width + 6 + 3 * s;
                priors[k] = torsionPriors(x[v], x[v + 1], x[v + 2]);
            }
            for (std::size_t c = 0; c < 4; ++c) {
                for (std::size_t k = 0; k < m_n; ++k) {
                    out.push_back(priors[k][c]);
                }
            }
        }
        for (std::size_t k = 0; k < m_n; ++k) {
            for (std::size_t c = 0; c < 3; ++c) {
                out.push_back(x[k * m_width + c] / kLevelRotation);
            }
        }
        for (std::size_t k = 0; k < m_n; ++k) {
            for (std::size_t c = 3; c < 6; ++c) {
                out.push_back(x[k * m_width + c] / kLevelShift);
            }
        }
        for (const ClosurePair& pr : pairs) {
            out.push_back(clip0(pr.limit - norm(atoms[pr.p] - atoms[pr.q])) / kClosureWeight);
        }
        return out;
    }

    // the level span (lo, hi) of every residual: MDATK's sparsity structure
    numerics::CsrMatrix structure(const std::vector<ClosurePair>& pairs) const {
        std::vector<std::size_t> lo, hi;
        for (std::size_t s = 0; s < m_ns; ++s) {
            const std::vector<std::size_t>& o = m_order[s];
            for (int c = 0; c < 6; ++c) {
                for (std::size_t j = 0; j + 1 < m_n; ++j) {
                    lo.push_back(std::min(o[j], o[j + 1]));
                    hi.push_back(std::max(o[j], o[j + 1]));
                }
            }
        }
        for (std::size_t r = 0; r < 4 * m_ns; ++r) {
            for (std::size_t k = 0; k < m_n; ++k) {
                lo.push_back(k);
                hi.push_back(k);
            }
        }
        for (int r = 0; r < 2; ++r) {
            for (std::size_t k = 0; k < m_n; ++k) {
                for (int c = 0; c < 3; ++c) {
                    lo.push_back(k);
                    hi.push_back(k);
                }
            }
        }
        for (const ClosurePair& p : pairs) {
            lo.push_back(p.lo);
            hi.push_back(p.hi);
        }
        numerics::CsrMatrix csr;
        csr.rows = lo.size();
        csr.cols = m_n * m_width;
        csr.indptr.assign(csr.rows + 1, 0);
        for (std::size_t r = 0; r < lo.size(); ++r) {
            for (std::size_t c = lo[r] * m_width; c < (hi[r] + 1) * m_width; ++c) {
                csr.indices.push_back(c);
                csr.data.push_back(1.0);
            }
            csr.indptr[r + 1] = csr.indices.size();
        }
        return csr;
    }

    // the pairs within the cutoff (MDATK's filters) added to those scored so far, (p, q) order;
    // whether the set grew
    bool refresh(std::span<const double> x) {
        const std::vector<Vec3> atoms = build(x);
        std::vector<std::size_t> real;
        std::vector<double> xyz;
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            if (m_cls[a] != 6) {
                real.push_back(a);
                xyz.insert(xyz.end(), atoms[a].begin(), atoms[a].end());
            }
        }
        const ContactTables& tables = contactTables();
        const std::size_t nStrand = kStrandNames.size();
        std::vector<ClosurePair> found;
        for (const auto& [fi, fj] : numerics::queryPairs(xyz, kClosureCutoff)) {
            const std::size_t p = real[fi], q = real[fj];
            const std::size_t sp = p / (m_n * kSlots), kp = (p / kSlots) % m_n, ip = p % kSlots;
            const std::size_t sq = q / (m_n * kSlots), kq = (q / kSlots) % m_n, iq = q % kSlots;
            const auto code = [&](std::size_t s, std::size_t k) {
                return static_cast<std::size_t>(m_codes[s * m_n + k]);
            };
            bool keep = kp != kq || (sp != sq && (ip < nStrand || iq < nStrand));
            const long long stepP = m_steps[sp], stepQ = m_steps[sq];
            const long long dk = static_cast<long long>(kq) - static_cast<long long>(kp);
            if (sp == sq && kp == kq) {
                keep = tables.inside[code(sp, kp)][ip * kSlots + iq];
            }
            if (sp == sq && dk == stepP) {  // q follows p on its strand
                keep = tables.link[code(sp, kp)][code(sq, kq)][ip * kSlots + iq];
            }
            if (sp == sq && -dk == stepP) {  // p follows q
                keep = tables.link[code(sq, kq)][code(sp, kp)][iq * kSlots + ip];
            }
            const double limit = limitOf(fitLimits(), m_cls[p], m_cls[q]);
            keep = keep && limit > 0;
            if (!keep) {
                continue;
            }
            const auto extended = [&](std::size_t k, std::size_t i, long long st) {
                long long e = static_cast<long long>(k);
                if (i == kSlotOP1 || i == kSlotOP2) {
                    e -= st;
                }
                return static_cast<std::size_t>(std::clamp<long long>(e, 0, static_cast<long long>(m_n) - 1));
            };
            const std::size_t ep = extended(kp, ip, stepP), eq = extended(kq, iq, stepQ);
            found.push_back({p, q, limit, std::min(std::min(kp, kq), std::min(ep, eq)),
                             std::max(std::max(kp, kq), std::max(ep, eq))});
        }
        const std::size_t before = m_pairs.size();
        std::vector<ClosurePair> merged;
        merged.reserve(m_pairs.size() + found.size());
        std::merge(m_pairs.begin(), m_pairs.end(), found.begin(), found.end(), std::back_inserter(merged),
                   [](const ClosurePair& a, const ClosurePair& b) { return std::tie(a.p, a.q) < std::tie(b.p, b.q); });
        // np.unique(..., return_index=True): the first of equal keys (the pair scored earlier)
        std::vector<ClosurePair> unique;
        unique.reserve(merged.size());
        for (const ClosurePair& c : merged) {
            if (unique.empty() || unique.back().p != c.p || unique.back().q != c.q) {
                unique.push_back(c);
            }
        }
        m_pairs = std::move(unique);
        return m_pairs.size() != before;
    }

    // one solve from x with the current pairs: (x, whether it stopped on max_nfev). The variables are
    // the move from x, so TRF's first trust radius is 1 (in x_scale's units), not the norm of x.
    std::pair<std::vector<double>, bool> solve(const std::vector<double>& x, std::size_t maxNfev) {
        numerics::LeastSquaresOptions o;
        o.ftol = o.xtol = o.gtol = kClosureTol;
        o.xScale = m_scale;
        o.maxNfev = maxNfev;
        o.jacSparsity = structure(m_pairs);
        o.denseFromSparsity = true;
        o.trSolver = numerics::TrustRegionSolver::Exact;
        o.jacobianThreads =
            m_options.threads ? m_options.threads : std::max<std::size_t>(1, std::thread::hardware_concurrency());
        const std::vector<ClosurePair>& pairs = m_pairs;
        const numerics::LeastSquaresCallback callback = [this](std::span<const double>, double, std::size_t,
                                                               std::size_t) {
            if (m_reporter) {
                m_reporter->checkCancel();
            }
            return false;
        };
        const auto moved = [this, &pairs, &x](std::span<const double> y) {
            std::vector<double> v(x.size());
            for (std::size_t i = 0; i < v.size(); ++i) {
                v[i] = x[i] + y[i];
            }
            return residuals(v, pairs);
        };
        const std::vector<double> y0(x.size(), 0.0);
        const numerics::LeastSquaresResult sol = numerics::leastSquares(moved, y0, o, callback);
        std::vector<double> out(x.size());
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = x[i] + sol.x[i];
        }
        return {std::move(out), sol.status == 0};
    }

    // solve at t from x with the pair refreshes: (x, whether the step had trouble)
    // (a step that can still be halved stops at its first trouble)
    std::pair<std::vector<double>, bool> step(std::vector<double> x, double t, std::size_t maxNfev, bool halvable) {
        setGeometry(t);
        bool trouble = false;
        for (int rounds = 0; rounds < kMaxRounds; ++rounds) {
            if (!refresh(x) && rounds > 0) {
                break;
            }
            auto [y, hit] = solve(x, maxNfev);
            x = std::move(y);
            trouble = trouble || hit || rounds >= kRoundsOk;
            if (trouble && halvable) {
                break;
            }
        }
        if (m_reporter) {
            m_reporter->progress(t, "closing the strands");
        }
        return {std::move(x), trouble};
    }

    std::vector<double> advance(const std::vector<double>& x, double t0, double t1, int halvings) {
        const bool last = halvings >= kMaxHalvings;
        const std::vector<ClosurePair> saved = m_pairs;
        auto [y, trouble] = step(x, t1, last ? kClosureFinalNfev : kClosureMaxNfev, !last);
        if (!trouble || last) {
            return y;
        }
        m_pairs = saved;
        const double mid = 0.5 * (t0 + t1);
        return advance(advance(x, t0, mid, halvings + 1), mid, t1, halvings + 1);
    }

    std::vector<std::vector<Atoms>> m_placed;
    std::vector<std::string> m_bases;
    std::vector<int> m_steps;
    const HelixOptions& m_options;
    std::optional<core::Reporter> m_reporter;
    std::size_t m_ns = 0, m_n = 0, m_width = 0;
    std::vector<std::vector<Move>> m_moves;
    std::vector<int> m_codes;
    std::vector<int> m_cls;
    std::vector<std::vector<std::size_t>> m_order;
    std::vector<double> m_x, m_scale;
    Geometry m_geometry;
    std::vector<Vec3> m_centre;
    std::vector<ClosurePair> m_pairs;
};

// --- strands and residues ------------------------------------------------------------------------

// _strand_heavy: atoms of a strand (5'->3', hydrogens included) with "_prev", the virtual O3' before
// the 5' residue
std::vector<Atoms> strandHeavy(const std::vector<Atoms>& placed, const std::string& bases, const std::vector<Vec3>& q) {
    std::vector<Atoms> out;
    for (std::size_t j = 0; j < placed.size(); ++j) {
        const char base = bases[j];
        Atoms heavy;
        for (const std::string_view name : baseSlotNames(base)) {
            if (name != "H3+") {
                heavy.set(name, placed[j].at(name));
            }
        }
        for (const TemplateAtom& a : nucleotideOf(placed[j], base, true, q[j]).items) {
            heavy.set(a.name, a.xyz);
        }
        const Vec3 prev = out.empty()
                              ? placeAtom(heavy.at("C5'"), heavy.at("O5'"), heavy.at("P"), kLink, kO3PO5, kAlpha)
                              : out.back().at("O3'");
        const auto [op1, op2] = phosphateOxygens(heavy.at("P"), prev, heavy.at("O5'"));
        heavy.set("OP1", op1);
        heavy.set("OP2", op2);
        heavy.set("_prev", prev);
        out.push_back(std::move(heavy));
    }
    return out;
}

// _residue_atoms: one rebuilt nucleotide in the order of _strand_residue
FittedResidue residueAtoms(const Atoms& heavy, char base, bool five, bool three, bool phosphate5, bool withHydrogens) {
    Atoms t = heavy;
    std::vector<FittedAtom> first, last;
    if (five) {
        if (phosphate5) {
            first.push_back({"OP3", t.at("P") + kPOp3 * unit(t.at("_prev") - t.at("P"))});
        } else {
            if (withHydrogens) {
                last.push_back({"HO5'", place(t.at("C4'"), t.at("C5'"), t.at("O5'"), kBondOH, 109.5, 180.0)});
            }
            for (const std::string_view name : {"P", "OP1", "OP2"}) {
                t.erase(name);
            }
        }
    }
    if (three && withHydrogens) {
        last.push_back({"HO3'", place(t.at("C4'"), t.at("C3'"), t.at("O3'"), kBondOH, 109.5, 180.0)});
    }
    FittedResidue out;
    out.base = base;
    out.atoms = first;
    for (const std::string_view name : backboneNames()) {
        if (const Vec3* v = t.find(name)) {
            out.atoms.push_back({std::string(name), *v});
        }
    }
    for (const std::string_view name : baseAtomNames(base)) {
        if (const Vec3* v = t.find(name)) {
            out.atoms.push_back({std::string(name), *v});
        }
    }
    if (withHydrogens) {
        const Atoms hyd = hydrogens(heavy, base);
        for (const std::string_view name : sugarHydrogenNames()) {
            out.atoms.push_back({std::string(name), hyd.at(name)});
        }
        for (const std::string_view name : baseHydrogenNames(base)) {
            out.atoms.push_back({std::string(name), hyd.at(name)});
        }
    }
    out.atoms.insert(out.atoms.end(), last.begin(), last.end());
    return out;
}

// _turn_terminal_hydroxyls: HO5' / HO3' of the strand ends on the staggered position (C-C-O-H 180,
// 60 or -60 degrees) farthest from the other atoms
void turnTerminalHydroxyls(std::vector<std::vector<FittedResidue>>& built) {
    struct Owner {
        std::size_t s, r;
        std::string name;
    };
    std::vector<Owner> owners;
    std::vector<Vec3> rows;
    for (std::size_t si = 0; si < built.size(); ++si) {
        for (std::size_t ri = 0; ri < built[si].size(); ++ri) {
            for (const FittedAtom& a : built[si][ri].atoms) {
                owners.push_back({si, ri, a.name});
                rows.push_back(a.xyz);
            }
        }
    }
    for (std::size_t si = 0; si < built.size(); ++si) {
        for (std::size_t ri = 0; ri < built[si].size(); ++ri) {
            std::vector<FittedAtom>& atoms = built[si][ri].atoms;
            const auto find = [&](std::string_view name) -> const Vec3* {
                const Vec3* found = nullptr;
                for (const FittedAtom& a : atoms) {
                    if (a.name == name) found = &a.xyz;  // dict(atoms): the last of a name
                }
                return found;
            };
            const struct {
                std::string_view h, a, b, o;
            } hydroxyls[] = {{"HO5'", "C4'", "C5'", "O5'"}, {"HO3'", "C4'", "C3'", "O3'"}};
            for (const auto& hx : hydroxyls) {
                std::size_t index = atoms.size();
                for (std::size_t k = 0; k < atoms.size(); ++k) {
                    if (atoms[k].name == hx.h) index = k;
                }
                if (index == atoms.size()) {
                    continue;
                }
                const Vec3 a = *find(hx.a), b = *find(hx.b), o = *find(hx.o);
                std::optional<std::pair<double, Vec3>> best;
                for (const double torsion : {180.0, 60.0, -60.0}) {
                    const Vec3 pos = place(a, b, o, kBondOH, 109.5, torsion);
                    double clear = 9.0;
                    for (std::size_t i = 0; i < rows.size(); ++i) {
                        const Owner& w = owners[i];
                        if (w.s == si && w.r == ri && (w.name == hx.h || w.name == hx.o || w.name == hx.b)) {
                            continue;
                        }
                        // query_ball_point(pos, 3.0): squared distance within r * r
                        const Vec3 dv = pos - rows[i];
                        const double d2 = (dv[0] * dv[0] + dv[1] * dv[1]) + dv[2] * dv[2];
                        if (d2 <= 9.0) {
                            clear = std::min(clear, norm(dv));
                        }
                    }
                    if (!best || clear > best->first + 1e-6) {
                        best = std::make_pair(clear, pos);
                    }
                }
                atoms[index].xyz = best->second;
            }
        }
    }
}

// _axis_frame: the screw axis to z through the origin (direction kept)
std::pair<Mat3, Vec3> axisFrame(const Screw& screw) {
    const Vec3 u = unit(screw.axis);
    const Vec3 z = {0.0, 0.0, 1.0};
    const double c = dot(u, z);
    const Vec3 v = cross(u, z);
    const double s = norm(v);
    const Mat3 rot = s < 1e-12 ? identity() : turnMatrix(v, degrees(std::atan2(s, c)));
    return {rot, -matVec(rot, screw.origin)};
}

}  // namespace

// _helix_strands with _residue_atoms and _turn_terminal_hydroxyls (build_helix's rebuilt part)
FittedHelix fittedHelix(const std::string& kind, const std::vector<std::string>& sequences, bool phosphate5,
                        bool withHydrogens, const HelixOptions& options) {
    const core::Reporter* reporter = options.reporter;
    if (reporter) {
        reporter->progress(0.0, "fitting the reference helix");
    }
    const HelixFitResult fit = helixFit(kind, reporter);
    const std::vector<StrandSpec> spec = helices(kind);
    const Screw& screw = fit.screw;
    const std::string& seq = sequences[0];
    const std::size_t n = seq.size();
    std::vector<std::vector<Atoms>> placed(spec.size());
    for (std::size_t k = 0; k < n; ++k) {
        placed[0].push_back(screwed(coreTemplate(seq[k]), screwAt(screw, static_cast<double>(k))));
    }
    std::vector<std::string> levelBases = {seq};
    std::vector<std::string> atLevel(n);
    for (std::size_t k = 0; k < n; ++k) {
        atLevel[k] = std::string(1, seq[k]);
    }
    FittedHelix out;
    bool uniform = true;
    for (std::size_t s = 1; s < spec.size(); ++s) {
        const StrandSpec& st = spec[s];
        const std::string& other = sequences[s];
        std::vector<FitLevel> lvls;
        std::string bases;
        for (std::size_t k = 0; k < n; ++k) {
            const char a = seq[k];
            const char b = st.step > 0 ? other[k] : other[n - 1 - k];
            atLevel[k] += b;
            const bool isCanonical = canonical(st.rule, a, b);
            const Pose pose = levelPose(kind, s, atLevel[k]);
            uniform = uniform && a == spec[0].base && b == st.base;
            placed[s].push_back(
                screwed(movedAtoms(coreTemplate(b), pose.rot, pose.shift), screwAt(screw, static_cast<double>(k))));
            FitLevel lv;
            lv.level = k;
            lv.bases = std::string(1, a) + b;
            lv.canonical = isCanonical;
            if (isCanonical) {
                lv.hbonds = pairHbonds(st.rule, a, b);
            }
            lvls.push_back(lv);
            bases += b;
        }
        levelBases.push_back(bases);
        out.levels.push_back(std::move(lvls));
    }
    std::vector<std::vector<Vec3>> torsions;
    for (std::size_t s = 0; s < spec.size(); ++s) {
        torsions.emplace_back(n, fit.torsions[s]);
    }
    if (!uniform) {  // also a single level: the torsions of the other bases
        std::vector<int> steps;
        for (const StrandSpec& st : spec) {
            steps.push_back(st.step);
        }
        HelixOptions o = options;
        std::optional<core::Reporter> sub;
        if (reporter) {
            sub = reporter->sub(0.05, 1.0);
            o.reporter = &*sub;
        }
        Closure closure(kind, placed, levelBases, steps, fit.torsions, screw, o);
        closure.run();
        for (std::size_t s = 0; s < spec.size(); ++s) {
            torsions[s] = closure.torsions(s);
        }
        placed = closure.moved();
    }
    for (std::size_t s = 0; s < spec.size(); ++s) {
        std::vector<Atoms> nucleosides = placed[s];
        std::string bases = levelBases[s];
        std::vector<Vec3> q = torsions[s];
        if (spec[s].step < 0) {
            std::reverse(nucleosides.begin(), nucleosides.end());
            std::reverse(bases.begin(), bases.end());
            std::reverse(q.begin(), q.end());
        }
        const std::vector<Atoms> heavy = strandHeavy(nucleosides, bases, q);
        std::vector<FittedResidue> residues;
        for (std::size_t j = 0; j < heavy.size(); ++j) {
            residues.push_back(residueAtoms(heavy[j], bases[j], j == 0, j + 1 == n, phosphate5, withHydrogens));
        }
        out.strands.push_back(std::move(residues));
    }
    if (withHydrogens) {
        turnTerminalHydroxyls(out.strands);
    }
    out.twist = screw.twist;
    out.rise = screw.rise;
    out.cost = fit.cost;
    out.contact = fit.contact;
    out.contactH = fit.contactH;
    std::tie(out.axisRotation, out.axisShift) = axisFrame(screw);
    if (reporter) {
        reporter->progress(1.0, "built");
    }
    return out;
}

// --- notes -----------------------------------------------------------------------------------------

namespace {

bool isHydrogenName(std::string_view name) {
    const std::size_t k = name.find_first_not_of("0123456789");
    return k != std::string_view::npos && name[k] == 'H';
}

std::string stripDigits(std::string_view name) {
    const std::size_t k = name.find_first_not_of("0123456789");
    return k == std::string_view::npos ? std::string() : std::string(name.substr(k));
}

bool allBaseName(std::string_view name) {  // _ALL_BASE_NAMES
    for (const char b : std::string_view("ACGU")) {
        if (contains(baseAtomNames(b), name) || contains(baseHydrogenNames(b), name)) {
            return true;
        }
    }
    static constexpr std::array<std::string_view, 7> extra = {"H62", "H22", "H42", "C7", "H71", "H72", "H73"};
    return contains(extra, name);
}

// the last value of a name (Python dict of a list of pairs)
const Vec3* lastOf(const ResAtoms& atoms, std::string_view name) {
    const Vec3* found = nullptr;
    for (const auto& [n, x] : atoms) {
        if (n == name) found = &x;
    }
    return found;
}

}  // namespace

// _close_contacts
std::vector<Contact> closeContacts(
    const std::vector<ResAtoms>& resAtoms, const std::vector<std::vector<std::size_t>>& segments,
    const std::vector<char>& bases,
    const std::vector<std::tuple<std::size_t, std::string, std::size_t, std::string>>& hbonds) {
    // _bonded_within
    std::vector<Vec3> xyz;
    std::vector<std::pair<std::size_t, std::string>> owner;
    std::vector<std::pair<std::size_t, std::size_t>> bonds;
    std::map<std::pair<std::size_t, std::string>, std::size_t> index;
    for (std::size_t r = 0; r < resAtoms.size(); ++r) {
        const std::size_t start = xyz.size();
        std::vector<bool> hyd;
        for (const auto& [name, x] : resAtoms[r]) {
            index[{r, name}] = xyz.size();
            xyz.push_back(x);
            owner.emplace_back(r, name);
            hyd.push_back(isHydrogenName(name));
        }
        const std::size_t count = xyz.size() - start;
        for (std::size_t i = 0; i < count; ++i) {
            for (std::size_t j = i + 1; j < count; ++j) {
                if (hyd[i] && hyd[j]) continue;
                const double cut = hyd[i] || hyd[j] ? 1.25 : 1.75;
                if (norm(xyz[start + i] - xyz[start + j]) < cut) {
                    bonds.emplace_back(start + i, start + j);
                }
            }
        }
    }
    for (const std::vector<std::size_t>& seg : segments) {
        for (std::size_t k = 0; k + 1 < seg.size(); ++k) {
            const auto o3 = index.find({seg[k], "O3'"});
            const auto p = index.find({seg[k + 1], "P"});
            if (o3 != index.end() && p != index.end()) {
                bonds.emplace_back(o3->second, p->second);
            }
        }
    }
    std::vector<Contact> out;
    if (xyz.empty()) {
        return out;
    }
    std::vector<int> cls(xyz.size());
    for (std::size_t a = 0; a < xyz.size(); ++a) {
        cls[a] = atomClass(bases[owner[a].first], stripDigits(owner[a].second));
    }
    std::vector<std::set<std::size_t>> neighbours(xyz.size());
    for (const auto& [a, b] : bonds) {
        neighbours[a].insert(b);
        neighbours[b].insert(a);
    }
    std::map<std::size_t, std::pair<std::set<std::size_t>, std::set<std::size_t>>> shells;
    const auto within = [&](std::size_t a) -> const std::pair<std::set<std::size_t>, std::set<std::size_t>>& {
        auto it = shells.find(a);
        if (it == shells.end()) {
            std::set<std::size_t> two, three;
            for (const std::size_t c : neighbours[a]) two.insert(neighbours[c].begin(), neighbours[c].end());
            for (const std::size_t c : two) three.insert(neighbours[c].begin(), neighbours[c].end());
            it = shells.emplace(a, std::make_pair(std::move(two), std::move(three))).first;
        }
        return it->second;
    };
    std::set<std::pair<std::size_t, std::size_t>> skip;
    for (const auto& [i, x, j, y] : hbonds) {
        const auto ia = index.find({i, x});
        const auto jb = index.find({j, y});
        if (ia != index.end() && jb != index.end()) {
            skip.insert({std::min(ia->second, jb->second), std::max(ia->second, jb->second)});
        }
    }
    std::vector<double> flat;
    for (const Vec3& v : xyz) flat.insert(flat.end(), v.begin(), v.end());
    const double cutoff = std::max(
        {kContactReport.heavy, kContactReport.hbond, kContactReport.hh, kContactReport.polar, kContactReport.hx});
    for (const auto& [a, b] : numerics::queryPairs(flat, cutoff)) {
        if (neighbours[a].count(b)) continue;
        const auto& [two, three] = within(a);
        if (two.count(b)) continue;
        const double d = norm(xyz[a] - xyz[b]);
        const bool hyd = cls[a] == 4 || cls[a] == 5 || cls[b] == 4 || cls[b] == 5;
        if (three.count(b)) {  // 1-4 pairs: heavy atoms only
            if (!hyd && d < kClashNb) {
                out.push_back({d, "heavy, 3 bonds, < " + fixed(kClashNb, 1), owner[a], owner[b]});
            }
            continue;
        }
        if (skip.count({a, b})) continue;
        const double limit = limitOf(reportLimits(), cls[a], cls[b]);
        if (d < limit) {
            std::string what;
            if ((cls[a] == 4 || cls[a] == 5) && (cls[b] == 4 || cls[b] == 5)) {
                what = "H-H";
            } else if (hyd) {
                what = limit == kContactReport.polar ? "H-bond H" : "H-heavy";
            } else {
                what = limit == kContactReport.hbond ? "N/O" : "heavy";
            }
            out.push_back({d, what + " < " + fixed(limit, 1), owner[a], owner[b]});
        }
    }
    std::stable_sort(out.begin(), out.end(), [](const Contact& x, const Contact& y) { return x.d < y.d; });
    return out;
}

// _geometry_notes
std::vector<std::string> geometryNotes(const std::vector<std::string>& labels, const std::vector<ResAtoms>& resAtoms,
                                       const std::vector<std::vector<std::size_t>>& segments,
                                       const std::vector<NotePair>& pairs, const std::vector<char>& bases) {
    std::vector<std::string> out;
    struct Link {
        std::size_t a, b;
        double d;
    };
    std::vector<Link> links;
    for (const std::vector<std::size_t>& seg : segments) {
        for (std::size_t k = 0; k + 1 < seg.size(); ++k) {
            const Vec3* o3 = lastOf(resAtoms[seg[k]], "O3'");
            const Vec3* p = lastOf(resAtoms[seg[k + 1]], "P");
            if (o3 && p) {
                links.push_back({seg[k], seg[k + 1], norm(*o3 - *p)});
            }
        }
    }
    if (!links.empty()) {
        double lo = links[0].d, hi = links[0].d;
        for (const Link& l : links) {
            lo = std::min(lo, l.d);
            hi = std::max(hi, l.d);
        }
        out.push_back("O3'-P links: " + fixed(lo, 2) + "-" + fixed(hi, 2) + " Å over " + std::to_string(links.size()) +
                      " links (ideal " + fixed(kLink, 2) + " Å)");
        std::vector<std::string> gaps;
        for (const Link& l : links) {
            if (!(kLinkLow <= l.d && l.d <= kLinkHigh)) {
                gaps.push_back(labels[l.a] + "-" + labels[l.b] + " " + fixed(l.d, 2) + " Å");
            }
        }
        if (!gaps.empty()) {
            std::string text;
            for (std::size_t k = 0; k < std::min<std::size_t>(8, gaps.size()); ++k) {
                text += (k ? ", " : "") + gaps[k];
            }
            out.push_back("WARNING: " + std::to_string(gaps.size()) +
                          " O3'-P link(s) outside 1.5-1.7 Å (the backbone "
                          "could not be closed there; minimise before use): " +
                          text + (gaps.size() > 8 ? " ..." : ""));
        }
    }
    std::vector<double> hb;
    std::size_t canonicalCount = 0;
    for (const NotePair& p : pairs) {
        if (!p.bonds.empty()) ++canonicalCount;
        for (const PairHbond& h : p.bonds) {
            hb.push_back(norm(*lastOf(resAtoms[p.i], h.atom1) - *lastOf(resAtoms[p.j], h.atom2)));
        }
    }
    if (!hb.empty()) {
        out.push_back("H-bonds of the " + std::to_string(canonicalCount) + " canonical pair(s): donor-acceptor " +
                      fixed(*std::min_element(hb.begin(), hb.end()), 2) + "-" +
                      fixed(*std::max_element(hb.begin(), hb.end()), 2) + " Å");
    }
    for (const NotePair& p : pairs) {
        if (!p.bonds.empty()) continue;
        double d = std::numeric_limits<double>::infinity();
        // coords[i]: a dict (the last atom of a name)
        const auto baseAtoms = [](const ResAtoms& atoms) {
            std::vector<Vec3> out2;
            std::vector<std::string> seen;
            for (const auto& [n, x] : atoms) {
                if (allBaseName(n) && n[0] != 'H' && std::find(seen.begin(), seen.end(), n) == seen.end()) {
                    seen.push_back(n);
                    out2.push_back(*lastOf(atoms, n));
                }
            }
            return out2;
        };
        const std::vector<Vec3> bi = baseAtoms(resAtoms[p.i]), bj = baseAtoms(resAtoms[p.j]);
        for (const Vec3& x : bi) {
            for (const Vec3& y : bj) {
                d = std::min(d, norm(x - y));
            }
        }
        out.push_back("WARNING: mismatch " + labels[p.i] + " · " + labels[p.j] + ": " + p.what +
                      " (closest base heavy atoms " + fixed(d, 1) + " Å" + (d < kClash ? ", the bases collide" : "") +
                      ")");
    }
    std::vector<std::tuple<std::size_t, std::string, std::size_t, std::string>> hbondAtoms;
    for (const NotePair& p : pairs) {
        for (const PairHbond& h : p.bonds) {
            hbondAtoms.emplace_back(p.i, std::string(h.atom1), p.j, std::string(h.atom2));
        }
    }
    const std::vector<Contact> contacts = closeContacts(resAtoms, segments, bases, hbondAtoms);
    const std::string limits = "heavy atoms < " + fixed(kContactReport.heavy, 1) + " Å (N/O pairs that can H-bond < " +
                               fixed(kContactReport.hbond, 1) + " Å; three bonds apart < " + fixed(kClashNb, 1) +
                               " Å), H-H < " + fixed(kContactReport.hh, 1) + " Å, H-heavy < " +
                               fixed(kContactReport.hx, 1) + " Å (H-bonded H < " + fixed(kContactReport.polar, 1) +
                               " Å); atoms four or more bonds apart, hydrogens included";
    if (!contacts.empty()) {
        std::string text;
        for (std::size_t k = 0; k < std::min<std::size_t>(6, contacts.size()); ++k) {
            const Contact& c = contacts[k];
            text += (k ? ", " : "") + labels[c.a.first] + " " + c.a.second + " - " + labels[c.b.first] + " " +
                    c.b.second + " " + fixed(c.d, 2) + " Å (" + c.what + ")";
        }
        out.push_back("WARNING: " + std::to_string(contacts.size()) + " close contact(s) (" + limits +
                      "; minimise before use): " + text + (contacts.size() > 6 ? " ..." : ""));
    } else {
        out.push_back("No close contacts (" + limits + ")");
    }
    return out;
}

namespace {

// _torsion_note: circular means of the backbone torsions per strand
std::string torsionNote(const io::PdbStructure& s) {
    std::vector<std::string> parts;
    for (const std::vector<io::PdbResidue>& seg : s.segments) {
        std::vector<ResAtoms> c;
        for (const io::PdbResidue& r : seg) {
            c.push_back(r.coords());
        }
        static constexpr std::array<std::string_view, 7> keys = {"alpha",   "beta", "gamma", "delta",
                                                                 "epsilon", "zeta", "chi"};
        std::array<std::vector<double>, 7> values;
        for (std::size_t j = 0; j < c.size(); ++j) {
            const ResAtoms& r = c[j];
            const auto at = [](const ResAtoms& atoms, std::string_view name) { return *lastOf(atoms, name); };
            const std::optional<io::ResidueName> info = io::parseResidueName(seg[j].name);
            const Glycosidic g = glycosidic(info ? info->base : 'A');
            values[6].push_back(dihedral(at(r, "O4'"), at(r, "C1'"), at(r, g.n), at(r, g.a)));
            values[3].push_back(dihedral(at(r, "C5'"), at(r, "C4'"), at(r, "C3'"), at(r, "O3'")));
            if (lastOf(r, "P")) {
                values[1].push_back(dihedral(at(r, "P"), at(r, "O5'"), at(r, "C5'"), at(r, "C4'")));
            }
            values[2].push_back(dihedral(at(r, "O5'"), at(r, "C5'"), at(r, "C4'"), at(r, "C3'")));
            if (j) {
                const ResAtoms& p = c[j - 1];
                values[0].push_back(dihedral(at(p, "O3'"), at(r, "P"), at(r, "O5'"), at(r, "C5'")));
                values[4].push_back(dihedral(at(p, "C4'"), at(p, "C3'"), at(p, "O3'"), at(r, "P")));
                values[5].push_back(dihedral(at(p, "C3'"), at(p, "O3'"), at(r, "P"), at(r, "O5'")));
            }
        }
        std::string means;
        for (std::size_t k = 0; k < keys.size(); ++k) {
            if (values[k].empty()) continue;
            std::vector<double> sines, cosines;
            for (const double v : values[k]) {
                sines.push_back(std::sin(radians(v)));
                cosines.push_back(std::cos(radians(v)));
            }
            means += (means.empty() ? "" : " ") + std::string(keys[k]) + " " +
                     fixed(degrees(std::atan2(core::detail::numpyMean(sines), core::detail::numpyMean(cosines))), 0);
        }
        parts.push_back((seg[0].chain.empty() ? std::string("_") : seg[0].chain) + ": " + means);
    }
    std::string text;
    for (std::size_t k = 0; k < parts.size(); ++k) {
        text += (k ? "; " : "") + parts[k];
    }
    return "Backbone torsions (means, degrees; PDB 4JRD A-form: alpha -62 beta 174 gamma 58 delta 80 epsilon -146 "
           "zeta -66 chi -175): " +
           text;
}

}  // namespace

// _model_notes
std::vector<std::string> modelNotes(const io::PdbStructure& s, const std::vector<std::vector<FitLevel>>& levels,
                                    const std::string& kind) {
    std::vector<std::string> labels;
    std::vector<char> bases;
    std::vector<ResAtoms> resAtoms;
    for (const io::PdbResidue* r : s.residues()) {
        labels.push_back(r->chain + ":" + std::to_string(r->resseq));
        const std::optional<io::ResidueName> info = io::parseResidueName(r->name);
        const char base = info ? info->base : 'A';
        bases.push_back(base);
        ResAtoms atoms;
        bool anyH = false;
        for (const io::PdbAtom& a : r->atoms) {
            atoms.emplace_back(a.name, a.xyz());
            anyH = anyH || a.name[0] == 'H';
        }
        if (!anyH) {  // a heavy-atom model: hydrogens for the check
            Atoms heavy;
            for (const auto& [n, x] : atoms) heavy.set(n, x);
            for (const TemplateAtom& h : hydrogens(heavy, base).items) {
                atoms.emplace_back(h.name, h.xyz);
            }
        }
        resAtoms.push_back(std::move(atoms));
    }
    std::vector<std::size_t> offsets;
    std::vector<std::vector<std::size_t>> segments;
    std::size_t k = 0;
    for (const std::vector<io::PdbResidue>& seg : s.segments) {
        offsets.push_back(k);
        std::vector<std::size_t> idx;
        for (std::size_t j = 0; j < seg.size(); ++j) idx.push_back(k + j);
        segments.push_back(std::move(idx));
        k += seg.size();
    }
    const std::size_t n = s.segments[0].size();
    const std::vector<StrandSpec> spec = helices(kind);
    std::vector<NotePair> pairs;
    for (std::size_t strand = 1; strand < spec.size(); ++strand) {
        const StrandSpec& st = spec[strand];
        const std::string ref = std::string(1, spec[0].base) + "·" + std::string(1, st.base);
        const bool third = (st.rule == "hoogsteen" || st.rule == "reverse-hoogsteen") && kind.rfind("triplex", 0) == 0;
        for (const FitLevel& lv : levels[strand - 1]) {
            const std::size_t j = st.step > 0 ? lv.level : n - 1 - lv.level;
            const bool hoog = (st.rule == "hoogsteen" || st.rule == "reverse-hoogsteen") &&
                              (lv.bases[0] == 'C' || lv.bases[0] == 'U');
            const std::string edge = hoog ? " (a strand-1 pyrimidine has no purine Hoogsteen edge)" : "";
            std::string note;
            if (third) {
                note = "no " + std::string(st.rule) + " third-strand pair of " + lv.bases[1] + " on " + lv.bases[0] +
                       edge + "; placed from the " + ref + " reference triple";
            } else if (st.rule == "wc") {
                note = "no Watson-Crick pair" + edge + "; placed from the " + ref + " reference pair";
            } else {
                note = "no " + std::string(pairingLw(st.rule)) + " pair" + edge + "; placed from the " + ref +
                       " reference pair";
            }
            pairs.push_back({offsets[0] + lv.level, offsets[strand] + j, lv.hbonds, note});
        }
    }
    std::vector<std::string> out = geometryNotes(labels, resAtoms, segments, pairs, bases);
    out.push_back(torsionNote(s));
    std::size_t plus = 0;
    for (std::size_t strand = 1; strand < spec.size(); ++strand) {
        if (spec[strand].rule != "hoogsteen") continue;
        for (const FitLevel& lv : levels[strand - 1]) {
            if (lv.canonical && (lv.bases == "GC" || lv.bases == "CG")) ++plus;
        }
    }
    if (plus) {
        out.push_back(std::to_string(plus) +
                      " G·C+ Hoogsteen pair(s): the cytosine is written as plain C (N3 "
                      "unprotonated); protonate N3 (C+) for the force field or restrain the pair (at neutral pH it "
                      "is weak)");
    }
    return out;
}

}  // namespace mdbridge::chem::helix_detail
