#pragma once
// Internals of the structure builder shared by helix.cpp (templates, sequences, Watson-Crick helices)
// and helix_fit.cpp (the other pairings and triplexes, milestone 6 W2-A). Not installed.
//
// Vector arithmetic in NumPy's order of operations; both sources are built without floating-point
// contraction (CMakeLists.txt).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "mdbridge/chem/helix.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::core {
class Reporter;
}

namespace mdbridge::chem::helix_detail {

// --- vectors in NumPy's order of operations ----------------------------------------------------

inline constexpr double kDegToRad = 3.141592653589793 / 180.0;  // math.radians: x * (pi / 180)

inline double radians(double degrees) {
    return degrees * kDegToRad;
}

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline Vec3 operator-(const Vec3& a) {
    return {-a[0], -a[1], -a[2]};
}
inline Vec3 operator*(double s, const Vec3& a) {
    return {s * a[0], s * a[1], s * a[2]};
}

// np.linalg.norm: sqrt(x . x)
inline double norm(const Vec3& v) {
    double sum = v[0] * v[0];
    sum += v[1] * v[1];
    sum += v[2] * v[2];
    return std::sqrt(sum);
}

inline Vec3 unit(const Vec3& v) {
    const double n = norm(v);
    return {v[0] / n, v[1] / n, v[2] / n};
}

// np.cross
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

using Mat3 = std::array<Vec3, 3>;

// v @ m.T
inline Vec3 matVec(const Mat3& m, const Vec3& v) {
    Vec3 out{};
    for (std::size_t j = 0; j < 3; ++j) {
        double sum = v[0] * m[j][0];
        sum += v[1] * m[j][1];
        sum += v[2] * m[j][2];
        out[j] = sum;
    }
    return out;
}

inline Vec3 cartesian(double r, double phi, double z) {
    const double p = radians(phi);
    return {r * std::cos(p), r * std::sin(p), z};
}

// _place: atom d bonded to c with angle b-c-d and torsion a-b-c-d (degrees)
inline Vec3 place(const Vec3& a, const Vec3& b, const Vec3& c, double bond, double angle, double torsion) {
    const Vec3 bc = unit(c - b);
    const Vec3 n = unit(cross(b - a, bc));
    const Vec3 m = cross(n, bc);
    const double th = radians(angle);
    const double ph = radians(torsion);
    const double d0 = -bond * std::cos(th);
    const double d1 = bond * std::sin(th) * std::cos(ph);
    const double d2 = bond * std::sin(th) * std::sin(ph);
    return c + d0 * bc + d1 * m + d2 * n;
}

// _turn_matrix: rotation by `degrees` about `axis` (Rodrigues)
inline Mat3 turnMatrix(const Vec3& axis, double degrees) {
    const Vec3 k = unit(axis);
    const Mat3 kx = {Vec3{0.0, -k[2], k[1]}, Vec3{k[2], 0.0, -k[0]}, Vec3{-k[1], k[0], 0.0}};
    Mat3 kk{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            double sum = kx[i][0] * kx[0][j];
            sum += kx[i][1] * kx[1][j];
            sum += kx[i][2] * kx[2][j];
            kk[i][j] = sum;
        }
    }
    const double t = radians(degrees);
    const double s = std::sin(t);
    const double c = 1.0 - std::cos(t);
    Mat3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out[i][j] = ((i == j ? 1.0 : 0.0) + s * kx[i][j]) + c * kk[i][j];
        }
    }
    return out;
}

// _screw: level k is turned by -k * TWIST about z and moved by -k * RISE
inline Vec3 screw(const Vec3& xyz, double level) {
    const double t = radians(-kHelixTwist * level);
    const Mat3 rot = {Vec3{std::cos(t), -std::sin(t), 0.0}, Vec3{std::sin(t), std::cos(t), 0.0}, Vec3{0.0, 0.0, 1.0}};
    return matVec(rot, xyz) + Vec3{0.0, 0.0, -kHelixRise * level};
}

// _dyad: the fibre two-fold about x
inline Vec3 dyad(const Vec3& xyz) {
    return {xyz[0] * 1.0, xyz[1] * -1.0, xyz[2] * -1.0};
}

// pdbtools.chiral_volume
inline double helixChiralVolume(const Vec3& centre, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 va = a - centre;
    const Vec3 vb = b - centre;
    const Vec3 vc = c - centre;
    const Vec3 cr = {vb[1] * vc[2] - vb[2] * vc[1], vb[2] * vc[0] - vb[0] * vc[2], vb[0] * vc[1] - vb[1] * vc[0]};
    return va[0] * cr[0] + va[1] * cr[1] + va[2] * cr[2];
}

// An ordered name -> coordinates map (Python dict)
class Atoms {
public:
    std::vector<TemplateAtom> items;

    const Vec3* find(std::string_view name) const {
        for (const TemplateAtom& a : items) {
            if (a.name == name) {
                return &a.xyz;
            }
        }
        return nullptr;
    }
    Vec3* find(std::string_view name) {
        for (TemplateAtom& a : items) {
            if (a.name == name) {
                return &a.xyz;
            }
        }
        return nullptr;
    }
    const Vec3& at(std::string_view name) const {
        const Vec3* v = find(name);
        if (!v) {
            throw std::logic_error("helix template: no atom " + std::string(name));
        }
        return *v;
    }
    void set(std::string_view name, const Vec3& xyz) {
        if (Vec3* v = find(name)) {
            *v = xyz;
        } else {
            items.push_back({std::string(name), xyz});
        }
    }
    void erase(std::string_view name) {
        items.erase(std::remove_if(items.begin(), items.end(), [&](const TemplateAtom& a) { return a.name == name; }),
                    items.end());
    }
};

// --- helix.cpp -------------------------------------------------------------------------------------

constexpr double kBondOH = 0.96;  // BOND_H["oh"]
constexpr double kBondNH = 1.01;  // BOND_H["nh"]
constexpr double kPOp3 = 1.48;    // P-OP3 of a 5'-phosphate
constexpr double kClash = 2.5;    // closest allowed heavy-atom distance of a mismatched pair's bases

struct Glycosidic {
    std::string_view n;
    std::string_view a;
    std::string_view b;
};
Glycosidic glycosidic(char base);
// _hydrogens: ideal hydrogens of a D-nucleotide (sugar rules, then the base's)
Atoms hydrogens(const Atoms& heavy, char base);
// template(base, hydrogens) of a checked base
Atoms templateAtoms(char base, bool withHydrogens);
// BACKBONE_NAMES
std::span<const std::string_view> backboneNames();
std::string fixed(double x, int precision);
bool parallel(std::string_view pairing);
std::string_view pairingLw(std::string_view pairing);
std::string_view pairingLabel(std::string_view pairing);
std::string_view triplexLabel(std::string_view motif);
bool canonical(std::string_view rule, char a, char b);
std::string elementOf(std::string_view name);

// --- helix_fit.cpp ---------------------------------------------------------------------------------

// _Level: one pair of a partner strand with strand 1 (strand-1 level, the two bases, canonical?, the
// canonical pair's H-bonds)
struct FitLevel {
    std::size_t level = 0;
    std::string bases;
    bool canonical = false;
    std::span<const PairHbond> hbonds;
};

struct FittedAtom {
    std::string name;
    Vec3 xyz;
};
// The rotation vector (radians, angle 0 ... pi) of a rotation matrix; near pi the axis comes from the
// diagonal (the off-diagonal differences vanish there).
Vec3 rotationVector(const Mat3& r);
// The restrained reference fit of an experimental kind (hoogsteen, triplex-hoogsteen) computed afresh
// (not the builder's cached one) with `threads` finite-difference threads: its parameters.
std::vector<double> referenceFitFresh(const std::string& kind, std::size_t threads);

struct FittedResidue {
    char base = 'A';
    bool protonated = false;  // D19 / CH-21: the C+ of a G·C+ Hoogsteen pair (the CCD's CH; with hydrogens HN3)
    std::vector<FittedAtom> atoms;
};

// The fitted helix of a kind (reverse-wc, hoogsteen, reverse-hoogsteen, hh, triplex-hoogsteen,
// triplex-reverse-hoogsteen) for its sequences (5'->3', strand 1 first): every strand's residues in
// the order _strand_residue writes them (5' end, 3' end, phosphate5, hydrogens as build_helix asks,
// terminal hydroxyls turned), the levels of every partner strand, and the reference helix.
struct FittedHelix {
    std::vector<std::vector<FittedResidue>> strands;
    std::vector<std::vector<FitLevel>> levels;
    double twist = 0.0;
    double rise = 0.0;
    double cost = 0.0;
    double contact = 0.0;
    double contactH = 0.0;
    Mat3 axisRotation{};  // _axis_frame: xyz @ rot.T + shift puts the screw axis on z
    Vec3 axisShift{};
};
FittedHelix fittedHelix(const std::string& kind, const std::vector<std::string>& sequences, bool phosphate5,
                        bool withHydrogens, const HelixOptions& options);

// atoms of one residue (name, xyz) in order, as the notes take them
using ResAtoms = std::vector<std::pair<std::string, Vec3>>;
struct Contact {
    double d;
    std::string what;
    std::pair<std::size_t, std::string> a, b;  // (residue, atom)
};
// _close_contacts(res_atoms, segments, bases, hbonds): sorted by distance
std::vector<Contact> closeContacts(
    const std::vector<ResAtoms>& resAtoms, const std::vector<std::vector<std::size_t>>& segments,
    const std::vector<char>& bases,
    const std::vector<std::tuple<std::size_t, std::string, std::size_t, std::string>>& hbonds);
struct NotePair {
    std::size_t i, j;
    std::span<const PairHbond> bonds;
    std::string what;
};
// _geometry_notes(labels, res_atoms, segments, pairs, bases)
std::vector<std::string> geometryNotes(const std::vector<std::string>& labels, const std::vector<ResAtoms>& resAtoms,
                                       const std::vector<std::vector<std::size_t>>& segments,
                                       const std::vector<NotePair>& pairs, const std::vector<char>& bases);

// _model_notes of the written structure (rounded coordinates, residue names A C G U)
std::vector<std::string> modelNotes(const io::PdbStructure& s, const std::vector<std::vector<FitLevel>>& levels,
                                    const std::string& kind);

}  // namespace mdbridge::chem::helix_detail
