// Port of MDATK's mdatk/render/molecule/scene.py (MDATK v0.0.4.1, the owner's code).

#include "mdbridge/render/molecule/scene.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <string>

#include "mdbridge/core/log.hpp"

namespace mdbridge::render::molecule {
namespace {

std::size_t at(std::int32_t i) { return static_cast<std::size_t>(i); }

bool any(const std::vector<bool>& v) { return std::find(v.begin(), v.end(), true) != v.end(); }

}  // namespace

std::pair<double, double> depthRange(const std::vector<Colored>& items, const Vec3& axis) {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    auto points = [&](const std::vector<Vec3>& q, double rad) {
        for (const Vec3& p : q) {
            const double z = p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2];
            if (std::isfinite(z)) {  // non-finite points are skipped by the renderers
                lo = std::min(lo, z - rad);
                hi = std::max(hi, z + rad);
            }
        }
    };
    auto maxRadius = [](const std::vector<double>& r) { return r.empty() ? 0.0 : *std::max_element(r.begin(), r.end()); };
    for (const Colored& it : items) {
        if (const auto* m = std::get_if<Mesh>(&it.prim)) {
            points(m->verts, 0.0);
        } else if (const auto* s = std::get_if<Spheres>(&it.prim)) {
            points(s->xyz, maxRadius(s->radius));
        } else if (const auto* k = std::get_if<Sticks>(&it.prim)) {
            const double rad = maxRadius(k->radius);
            points(k->a, rad);
            points(k->b, rad);
        }
    }
    if (!std::isfinite(lo)) {
        return {-1.0, 1.0};
    }
    return {lo - 1.0, hi + 1.0};
}

std::vector<std::vector<bool>> residueMasks(const std::vector<core::RepSpec>& specs,
                                            const std::vector<std::vector<bool>>& member,
                                            const std::vector<std::string>& group_names, std::int32_t n_res) {
    const auto nr = static_cast<std::size_t>(std::max(n_res, 0));
    std::vector<std::vector<bool>> out;
    for (const core::RepSpec& spec : specs) {
        if (spec.groups.empty()) {
            out.emplace_back(nr, true);
            continue;
        }
        std::vector<bool> mask(nr, false);
        for (const std::string& name : spec.groups) {
            const auto it = std::find(group_names.begin(), group_names.end(), name);
            if (it == group_names.end()) {
                continue;
            }
            const auto g = static_cast<std::size_t>(it - group_names.begin());
            if (g < member.size()) {
                for (std::size_t r = 0; r < nr && r < member[g].size(); ++r) {
                    if (member[g][r]) {
                        mask[r] = true;
                    }
                }
            }
        }
        out.push_back(std::move(mask));
    }
    return out;
}

MoleculeScene::MoleculeScene(std::shared_ptr<const Topology> topo, const std::vector<core::RepSpec>& specs,
                             const std::vector<std::vector<bool>>& masks)
    : topo_(std::move(topo)) {
    for (std::size_t k = 0; k < specs.size() && k < masks.size(); ++k) {
        add(specs[k], masks[k]);
    }
    nuc_rgb_.reserve(topo_->letter.size());
    for (char letter : topo_->letter) {
        nuc_rgb_.push_back(nucleotideRgb(letter));
    }
}

void MoleculeScene::add(const core::RepSpec& spec, std::vector<bool> mask) {
    specs_.push_back(spec);
    std::unique_ptr<Style> style;
    mask.resize(static_cast<std::size_t>(topo_->n_res), false);
    if (any(mask)) {
        StyleOptions options;
        options.size = spec.size;
        options.sugars = spec.sugars;
        style = makeStyle(spec.style, topo_, std::move(mask), options);
    }
    styles_.push_back(std::move(style));
}

std::vector<bool> MoleculeScene::drawn() const {
    std::vector<bool> out(static_cast<std::size_t>(topo_->n_res), false);
    for (const auto& style : styles_) {
        if (style) {
            const std::vector<bool> d = style->drawn();
            for (std::size_t r = 0; r < out.size() && r < d.size(); ++r) {
                out[r] = out[r] || d[r];
            }
        }
    }
    return out;
}

void MoleculeScene::fillRest() {
    const Topology& topo = *topo_;
    const auto nr = static_cast<std::size_t>(topo.n_res);
    std::vector<bool> has_atoms(nr, false);  // residues with atoms
    for (std::int32_t r : topo.row) {
        if (r >= 0 && at(r) < nr) {
            has_atoms[at(r)] = true;
        }
    }
    const std::vector<bool> d = drawn();
    std::vector<bool> chain(nr, false), ligands(nr, false);
    for (std::size_t r = 0; r < nr; ++r) {
        const bool rest = !d[r] && has_atoms[r];
        chain[r] = rest && topo.control[r] >= 0;
        ligands[r] = rest && topo.control[r] < 0;
    }
    if (any(chain)) {
        core::RepSpec spec;
        spec.style = "trace";
        spec.color = "residue";
        add(spec, std::move(chain));
    }
    if (any(ligands)) {
        core::RepSpec spec;
        spec.style = "licorice";
        spec.color = "auto";
        add(spec, std::move(ligands));
    }
}

std::vector<Colored> MoleculeScene::frame(std::span<const Vec3> xyz, std::span<const Vec3> residue_rgb) {
    std::vector<Colored> out;
    for (std::size_t k = 0; k < specs_.size(); ++k) {
        const Style* style = styles_[k].get();
        if (style == nullptr) {
            continue;
        }
        const core::RepSpec& spec = specs_[k];
        DrawList dl;
        try {
            dl = style->build(xyz);
        } catch (const std::exception& exc) {
            if (failed_.insert(k).second) {
                core::logWarning("3D representation " + std::to_string(k + 1) + " (" + spec.style +
                                 ") not drawn on some frames: " + exc.what());
            }
            continue;
        }
        for (Mesh& m : dl.meshes) {
            std::vector<Vec3> rgb = meshRgb(m, spec, residue_rgb);
            out.push_back({std::move(m), std::move(rgb), spec.opacity, k});
        }
        for (Sticks& s : dl.sticks) {
            std::vector<Vec3> rgb = atomRgb(s, spec, residue_rgb);
            out.push_back({std::move(s), std::move(rgb), spec.opacity, k});
        }
        for (Spheres& s : dl.spheres) {
            std::vector<Vec3> rgb = atomRgb(s, spec, residue_rgb);
            out.push_back({std::move(s), std::move(rgb), spec.opacity, k});
        }
    }
    return out;
}

std::vector<Vec3> MoleculeScene::meshRgb(const Mesh& m, const core::RepSpec& spec,
                                         std::span<const Vec3> res_rgb) const {
    const std::size_t n = m.owner.size();
    std::vector<Vec3> rgb(n);
    if (spec.color == "nucleotide") {
        for (std::size_t v = 0; v < n; ++v) {
            rgb[v] = nuc_rgb_[at(m.owner[v])];
        }
        return rgb;
    }
    if (spec.color == "element" && m.has_atom) {  // surfaces: the nearest atom's element
        for (std::size_t v = 0; v < n; ++v) {
            rgb[v] = topo_->element_rgb[at(m.atom[v])];
        }
        return rgb;
    }
    const bool by_base = spec.color == "auto" || spec.color == "element";
    for (std::size_t v = 0; v < n; ++v) {
        const std::int8_t role = v < m.role.size() ? m.role[v] : std::int8_t{kRoleResidue};
        Vec3 c = by_base && role == kRoleNucleotide ? nuc_rgb_[at(m.owner[v])] : res_rgb[at(m.owner[v])];
        if (role == kRoleSugar) {  // lighter tint of the residue colour
            c = {0.55 * c[0] + 0.45, 0.55 * c[1] + 0.45, 0.55 * c[2] + 0.45};
        }
        rgb[v] = c;
    }
    return rgb;
}

template <class P>
std::vector<Vec3> MoleculeScene::atomRgb(const P& p, const core::RepSpec& spec, std::span<const Vec3> res_rgb) const {
    const std::size_t n = p.owner.size();
    std::vector<Vec3> rgb(n);
    for (std::size_t i = 0; i < n; ++i) {
        rgb[i] = res_rgb[at(p.owner[i])];
    }
    if (spec.color == "residue") {
        return rgb;
    }
    if (spec.color == "nucleotide") {
        for (std::size_t i = 0; i < n; ++i) {
            rgb[i] = nuc_rgb_[at(p.owner[i])];
        }
        return rgb;
    }
    const bool keep_carbons = spec.color == "auto";  // carbons keep the residue colour (groups, glow)
    for (std::size_t i = 0; i < n; ++i) {
        const std::int32_t a = p.atom[i];
        if (a < 0) {
            continue;
        }
        if (keep_carbons && topo_->is_carbon[at(a)]) {
            continue;
        }
        rgb[i] = topo_->element_rgb[at(a)];
    }
    return rgb;
}

std::unique_ptr<MoleculeScene> buildScene(const TopologyInput& atoms, const std::vector<std::vector<bool>>& member,
                                          const std::vector<std::string>& group_names,
                                          const std::vector<core::RepSpec>& specs, std::string_view rest) {
    if (specs.empty()) {
        return nullptr;
    }
    auto topo = std::make_shared<const Topology>(atoms);
    const auto masks = residueMasks(specs, member, group_names, topo->n_res);
    auto scene = std::make_unique<MoleculeScene>(topo, specs, masks);
    if (rest == "trace") {
        scene->fillRest();
    }
    return scene;
}

}  // namespace mdbridge::render::molecule
