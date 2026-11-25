#pragma once
// Representation scene (MDATK's mdatk/render/molecule/scene.py): which residues each style draws,
// per-frame geometry and colours. Renderer independent: MoleculeScene::frame returns the
// primitives of a frame with per-primitive (or per-vertex) RGB colours, which a renderer turns into
// pixels. No Qt.

#include <cstddef>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "mdbridge/core/project.hpp"
#include "mdbridge/render/molecule/geometry.hpp"
#include "mdbridge/render/molecule/topology.hpp"

namespace mdbridge::render::molecule {

// A primitive with its colours (RGB 0..1: per vertex for a Mesh, per sphere / stick otherwise) and
// opacity; `representation` is the index of the scene's representation that made it.
struct Colored {
    std::variant<Mesh, Spheres, Sticks> prim;
    std::vector<Vec3> rgb;
    double opacity = 1.0;
    std::size_t representation = 0;
};

// Per-render look shared by every renderer: depth cue and glow.
struct Shading {
    Vec3 fog_rgb{1.0, 1.0, 1.0};    // colour distant parts fade toward (the panel surface)
    double fog_strength = 0.0;      // 0 = no depth cue
    std::optional<std::vector<double>> glow;  // per residue row, 0..1
    Vec3 glow_rgb{1.0, 0.69, 0.0};
};

// View-depth extent (Å along `axis`) of every primitive, radii included, plus 1 Å on each side;
// non-finite points are skipped; (-1, 1) when nothing is finite.
std::pair<double, double> depthRange(const std::vector<Colored>& items, const Vec3& axis);

// Residue rows drawn by each representation: its groups by name (unknown names select nothing),
// every residue when it names none. member[g][r]: residue r is in group g.
std::vector<std::vector<bool>> residueMasks(const std::vector<core::RepSpec>& specs,
                                            const std::vector<std::vector<bool>>& member,
                                            const std::vector<std::string>& group_names, std::int32_t n_res);

class MoleculeScene {
public:
    MoleculeScene(std::shared_ptr<const Topology> topo, const std::vector<core::RepSpec>& specs,
                  const std::vector<std::vector<bool>>& masks);

    // One more representation; an unknown style or an empty mask draws nothing (style nullptr).
    void add(const core::RepSpec& spec, std::vector<bool> mask);
    // Residues that some representation actually draws.
    std::vector<bool> drawn() const;
    // Draw the residues no representation draws: the chain as a trace (colour "residue"), ligands
    // (no backbone atom) as licorice (colour "auto").
    void fillRest();
    // Primitives of one frame with colours; residue_rgb is per residue row. A style that fails on a
    // frame (e.g. a blown-up structure) is left out of that frame, with one warning per style
    // ("3D representation K (style) not drawn on some frames: <reason>").
    std::vector<Colored> frame(std::span<const Vec3> xyz, std::span<const Vec3> residue_rgb);

    const Topology& topology() const { return *topo_; }
    const std::vector<core::RepSpec>& specs() const { return specs_; }
    // The style of each representation (nullptr when it draws nothing).
    const Style* style(std::size_t k) const { return styles_[k].get(); }
    std::size_t size() const { return specs_.size(); }
    // Nucleotide colour per residue row.
    const std::vector<Vec3>& nucleotideColours() const { return nuc_rgb_; }

private:
    std::vector<Vec3> meshRgb(const Mesh& m, const core::RepSpec& spec, std::span<const Vec3> res_rgb) const;
    template <class P>
    std::vector<Vec3> atomRgb(const P& p, const core::RepSpec& spec, std::span<const Vec3> res_rgb) const;

    std::shared_ptr<const Topology> topo_;
    std::vector<core::RepSpec> specs_;
    std::vector<std::unique_ptr<Style>> styles_;
    std::vector<Vec3> nuc_rgb_;
    std::set<std::size_t> failed_;
};

// The scene of the cached atoms and trace results (nullptr without representations): the topology
// from the atoms (its frame 0 as reference), masks from the trace's group membership; residues no
// representation draws get a trace (ligands: licorice) unless `rest` is "hide".
std::unique_ptr<MoleculeScene> buildScene(const TopologyInput& atoms, const std::vector<std::vector<bool>>& member,
                                          const std::vector<std::string>& group_names,
                                          const std::vector<core::RepSpec>& specs, std::string_view rest = "trace");

}  // namespace mdbridge::render::molecule
