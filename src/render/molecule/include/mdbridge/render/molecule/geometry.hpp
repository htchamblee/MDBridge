#pragma once
// Representation geometry, renderer independent (MDATK's mdatk/render/molecule/geometry.py).
//
// Every style turns one frame's coordinates into a few primitive lists:
//   * Mesh: triangles with per-vertex normals,
//   * Spheres: centres and radii,
//   * Sticks: cylinders between two points.
// Each primitive carries `owner` (the residue row it belongs to) and `atom` (an atom index or -1),
// so colours can be assigned per frame from the residue colours (groups, proximity, glow) or from
// elements / nucleotides without rebuilding the geometry. Sizes are in Å.
//
// Styles (`size` scales the widths and radii):
//   cartoon   elliptical backbone through the control atoms, a rung from the backbone through
//             C1' to the glycosidic N, and the base ring as a slab extruded ±0.45 Å (purines as
//             the fused 9-atom perimeter); protein segments get a round tube (0.6 Å);
//   ribbon    a flat band through the control atoms facing the bases;
//   tube      a round tube through the control atoms;
//   licorice  half-bond sticks and atom balls (radius 0.3 Å);
//   vdw       van der Waals spheres (Bondi radii);
//   surf      a smooth Gaussian molecular surface (surface.hpp);
//   trace     balls at the backbone points joined by thin sticks (the 2D trace in 3D).

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/render/molecule/topology.hpp"

namespace mdbridge::render::molecule {

using Tri = std::array<std::int32_t, 3>;

enum Role : std::int8_t { kRoleResidue = 0, kRoleNucleotide = 1, kRoleSugar = 2 };

struct Mesh {
    std::vector<Vec3> verts;
    std::vector<Vec3> normals;
    std::vector<Tri> tris;
    std::vector<std::int32_t> owner;  // residue row per vertex
    std::vector<std::int8_t> role;    // colour role per vertex (Role)
    std::vector<std::int32_t> atom;   // nearest atom per vertex (surfaces); empty otherwise
    bool has_atom = false;
};

struct Spheres {
    std::vector<Vec3> xyz;
    std::vector<double> radius;
    std::vector<std::int32_t> owner;
    std::vector<std::int32_t> atom;  // atom index (element colour) or -1
};

struct Sticks {
    std::vector<Vec3> a;
    std::vector<Vec3> b;
    std::vector<double> radius;
    std::vector<std::int32_t> owner;
    std::vector<std::int32_t> atom;
};

struct DrawList {
    std::vector<Mesh> meshes;
    std::vector<Spheres> spheres;
    std::vector<Sticks> sticks;

    void extend(DrawList other);
};

// ---- splines and swept profiles -------------------------------------------------------------------

// Dense linear maps from m control points to (m-1)*n_sub + 1 samples of a uniform Catmull-Rom
// spline through them (phantom end points 2 p0 - p1 and 2 p[m-1] - p[m-2]): positions A @ P,
// derivatives D @ P (row-major samples x m), and the control point each sample belongs to (colours
// switch half-way between residues). m >= 2, n_sub >= 1.
struct DenseSpline {
    std::size_t samples = 0, m = 0;
    std::vector<double> a, d;
    std::vector<std::int32_t> owner;
};
DenseSpline catmullRom(std::size_t m, std::size_t n_sub);

// catmullRom in banded form: per sample the 4 indices into the control points extended by the
// two phantom end points (extendControls), position and derivative weights, and the owner.
struct BandedSpline {
    std::vector<std::array<std::int32_t, 4>> idx;
    std::vector<std::array<double, 4>> w, d;
    std::vector<std::int32_t> owner;
};
BandedSpline catmullRomBanded(std::size_t m, std::size_t n_sub);

// Control points with the phantom ends p[-1] = 2 p0 - p1 and p[m] = 2 p[m-1] - p[m-2] (m >= 2).
std::vector<Vec3> extendControls(std::span<const Vec3> p);

// Cross-section points (Å) and outward normals in the (binormal, normal) plane: "ellipse" (k
// points) or "band" (8 points, duplicated corners); other kinds throw std::invalid_argument.
struct Profile {
    std::vector<std::array<double, 2>> pts, normals;
};
Profile profile(std::string_view kind, double half_width, double half_height, std::size_t k = 12);

// Triangles of ns rings of k vertices (outward winding), plus fan caps whose centres are the two
// extra vertices ns*k (start) and ns*k + 1 (end).
std::vector<Tri> sweepFaces(std::size_t ns, std::size_t k, bool caps = true);

// Unit tangent, normal (toward `side` when given and never parallel to the tangent, else
// parallel-transported: a rotation-minimising frame) and binormal along a curve.
struct Frames {
    std::vector<Vec3> t, n, b;
};
Frames framesAlong(std::span<const Vec3> pos, std::span<const Vec3> tan, const std::vector<Vec3>* side);

// Closed slabs from ordered ring atoms (n rings of k points, ring after ring) extruded ±half
// along each ring's Newell normal: verts, normals, triangles and the slab index of each vertex.
struct RingSlabs {
    std::vector<Vec3> verts, normals;
    std::vector<Tri> tris;
    std::vector<std::int32_t> slab;
};
RingSlabs ringSlabs(std::span<const Vec3> ring_points, std::size_t k, double half);

// Spline through the control atoms of whole backbone segments, swept with a profile over the
// stretches owned by the selected residues (so a selection does not change the shape, and one
// selected residue still gets its piece). Segments without bases (proteins) get a round tube of
// radius `round_radius` when given.
class Backbone {
public:
    Backbone(const Topology& topo, const std::vector<bool>& res_mask, std::string_view kind, double half_width,
             double half_height, std::size_t n_sub = 8, std::size_t k = 12,
             std::optional<double> round_radius = std::nullopt);

    // Meshes, and the spline point at each selected residue's control atom (for rungs).
    std::pair<std::vector<Mesh>, std::map<std::int32_t, Vec3>> build(std::span<const Vec3> xyz) const;

    // Residues the backbone draws.
    const std::vector<bool>& drawn() const { return drawn_; }

private:
    struct Part {
        std::vector<std::int32_t> seg;
        std::vector<std::array<std::int32_t, 4>> idx;
        std::vector<std::array<double, 4>> w, d;
        std::vector<std::int32_t> own;
        std::vector<std::int32_t> at_res;
        std::vector<std::size_t> at;
        bool round = false;
    };
    const Topology* topo_;
    Profile flat_, round_;
    std::vector<Part> parts_;
    std::vector<bool> drawn_;
};

// ---- styles ---------------------------------------------------------------------------------------

struct StyleOptions {
    double size = 1.0;
    bool sugars = false;              // cartoon: ribose slabs
    std::size_t samples = 8;          // spline samples per residue (cartoon, ribbon, tube)
    std::optional<double> spacing;    // surf: grid spacing (default 1 Å, 1.5 Å above 5000 atoms)
};

// Builds one style's primitives for the residues in a mask (per residue row).
class Style {
public:
    Style(std::shared_ptr<const Topology> topo, std::vector<bool> mask, const StyleOptions& options);
    virtual ~Style() = default;
    Style(const Style&) = delete;
    Style& operator=(const Style&) = delete;

    virtual std::string_view name() const = 0;
    // The primitives of one frame (per-atom coordinates).
    virtual DrawList build(std::span<const Vec3> xyz) const = 0;
    // Residues this style actually draws (others of the mask produce no geometry).
    virtual std::vector<bool> drawn() const;

    const std::vector<bool>& mask() const { return mask_; }
    const Topology& topology() const { return *topo_; }

protected:
    // Atoms of the selected residues.
    std::vector<std::int32_t> selectedAtoms() const;

    std::shared_ptr<const Topology> topo_;
    std::vector<bool> mask_;
    StyleOptions options_;
};

// The style called `name` (cartoon, ribbon, tube, licorice, vdw, surf, trace), or nullptr.
std::unique_ptr<Style> makeStyle(std::string_view name, std::shared_ptr<const Topology> topo, std::vector<bool> mask,
                                 const StyleOptions& options = {});
// The style names in MDATK's order.
const std::vector<std::string>& styleNames();

}  // namespace mdbridge::render::molecule
